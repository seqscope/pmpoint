#include "pmpoint.h"
#include "qgenlib/tsv_reader.h"
#include "qgenlib/qgen_error.h"

#include <vector>
#include <string>
#include <cstring>
#include <climits>
#include <map>
#include <mutex>
#include <thread>
#include <atomic>
#include <queue>
#include <condition_variable>

#include "pmt_pts.h"
#include "pmt_utils.h"
#include "polygon.h"
#include "mvt_pts.h"
#include "htslib/hts.h"
#include "ext/nlohmann/json.hpp"

#define MAX_BITS 12
class tile_density_summary {
public:
    int32_t maxbits;
    std::vector< std::map<uint64_t, uint64_t> > res2pts2bins;
    
    tile_density_summary(int32_t _maxbits = MAX_BITS) {
        maxbits = _maxbits;
        res2pts2bins.resize(maxbits);
    }

    void merge_summary(const tile_density_summary& tds) {
        for(int32_t i=0; i<maxbits; ++i) {
            const auto& pts2bins = tds.res2pts2bins[i];
            for(auto it=pts2bins.begin(); it!=pts2bins.end(); ++it) {
                uint64_t pts = it->first;
                uint64_t nbins = it->second;
                res2pts2bins[i][pts] += nbins;
            }
        }
    }
};

// calculate density statistics for each tile
// for each resolution bit, 1,2,4,8,...,4096
// keep track of the number of "spots" and the number of "points" per spot
class tile_density_data {
public:
    int32_t maxbits;
    std::vector< std::map<uint64_t, uint64_t> > res2bin2pts;
    tile_density_data(int32_t _maxbits = MAX_BITS) {
        maxbits = _maxbits;
        res2bin2pts.resize(maxbits);
    }

    void add_xycnt(int32_t x, int32_t y, int32_t cnt) {
        if ( cnt > 0 ) {
            for(int32_t i=0; i< maxbits; ++i) {
                int32_t bsize = 1 << i;
                uint64_t bin = ((((uint64_t)x / bsize) << 32) | ((uint64_t)y / bsize));
                res2bin2pts[i][bin] += cnt;
            }
        }
    }

    // summary stats per tile
    // std::map<uint64_t, uint64_t> pts2nbins;
    // void compute_pts2nbins(int32_t res) {
    //     pts2nbins.clear();
    //     std::map<uint64_t, uint64_t> &bin2pts = res2bin2pts[res];
    //     for(auto it=bin2pts.begin(); it!=bin2pts.end(); ++it) {
    //         uint64_t pts = it->second;
    //         ++pts2nbins[pts];
    //     }
    // }
    void compute_summary(tile_density_summary& tds) {
        for(int32_t i=0; i<maxbits; ++i) {
            std::map<uint64_t, uint64_t> &bin2pts = res2bin2pts[i];
            for(auto it=bin2pts.begin(); it!=bin2pts.end(); ++it) {
                uint64_t pts = it->second;
                tds.res2pts2bins[i][pts] += 1;
            }
        }
    }
};

// Thread-safe result aggregator
class ThreadSafeResults {
private:
    std::mutex mutex;
    tile_density_summary tds_all;    
    std::atomic<uint64_t> n_total{0};
    htsFile *tsv_wh;
    bool compact;

public:
    ThreadSafeResults(htsFile* _tsv_wh, bool _compact) : tsv_wh(_tsv_wh), compact(_compact) {}

    void add_result(const tile_density_summary& tds, int32_t z, int32_t x, int32_t y, int32_t n_pts) {
        std::lock_guard<std::mutex> lock(mutex);
        
        n_total += n_pts;
        
        // Write individual tile results if not in compact mode
        if (!compact) {
            for(int32_t i=0; i<MAX_BITS; ++i) {
                const auto& pts2bins = tds.res2pts2bins[i];
                for(auto it=pts2bins.begin(); it!=pts2bins.end(); ++it) {
                    uint64_t pts = it->first;
                    uint64_t nbins = it->second;
                    hprintf(tsv_wh, "%d\t%d\t%d\t%d\t%llu\t%llu\n", z, x, y, 1<<i, pts, nbins);
                }
            }
        }
        
        // Aggregate results
        tds_all.merge_summary(tds);
    }
    
    uint64_t get_total_points() const {
        return n_total;
    }
    
    void write_summary(int32_t zoom) {
        std::lock_guard<std::mutex> lock(mutex);
        for(int32_t i=0; i<MAX_BITS; ++i) {
            const auto& pts2bins = tds_all.res2pts2bins[i];
            for(auto it=pts2bins.begin(); it!=pts2bins.end(); ++it) {
                uint64_t pts = it->first;
                uint64_t nbins = it->second;
                hprintf(tsv_wh, "%d\tALL\tALL\t%d\t%llu\t%llu\n", zoom, 1<<i, pts, nbins);
            }
        }
    }
};

// Thread-safe work queue for tile processing
class TileQueue {
private:
    std::queue<pmtiles::entry_zxy> work_queue;
    std::mutex mutex;
    std::condition_variable cv;
    bool done = false;
    
public:
    void add_tile(const pmtiles::entry_zxy& entry) {
        std::lock_guard<std::mutex> lock(mutex);
        work_queue.push(entry);
        cv.notify_one();
    }
    
    bool get_tile(pmtiles::entry_zxy& entry) {
        std::unique_lock<std::mutex> lock(mutex);
        cv.wait(lock, [this] { return !work_queue.empty() || done; });
        
        if (work_queue.empty() && done) {
            return false;
        }
        
        entry = work_queue.front();
        work_queue.pop();
        return true;
    }
    
    void finish() {
        std::lock_guard<std::mutex> lock(mutex);
        done = true;
        cv.notify_all();
    }
    
    bool is_empty() {
        std::lock_guard<std::mutex> lock(mutex);
        return work_queue.empty();
    }
};

// Worker thread function
void process_tiles(TileQueue& queue, pmt_pts& pmt, ThreadSafeResults& results, const std::string& count_field, std::atomic<int32_t>& processed_count) {
    pmtiles::entry_zxy entry(0,0,0,0,0);
    mvt_pts mvt;
    std::string buffer;
    
    while (queue.get_tile(entry)) {
        notice("Processing tile %d/%d/%d", entry.z, entry.x, entry.y);
        pmt_pts& local_pmt = pmt; // Create a thread-local copy for thread safety
        local_pmt.fetch_tile_to_buffer(entry.z, entry.x, entry.y, buffer);
        
        std::vector<int32_t> xs, ys, cnts;
        //int32_t n_pts = mvt.decode_points_xycnt(local_pmt.tile_data_str, count_field, xs, ys, cnts);
        int32_t n_pts = mvt.decode_points_xycnt(buffer, count_field, xs, ys, cnts);
        
        tile_density_data tdd;
        tile_density_summary tds;
        for (int32_t j = 0; j < n_pts; ++j) {
            tdd.add_xycnt(xs[j], ys[j], cnts[j]);
        }
        
        // Compute statistics
        tdd.compute_summary(tds);
        
        // Add result to the aggregator
        results.add_result(tds, entry.z, entry.x, entry.y, n_pts);
        
        // Update progress
        int32_t current = ++processed_count;
        if (current % 100 == 0) {
            notice("Processed %d tiles", current);
        }
    }
}

/////////////////////////////////////////////////////////////////////////
// tile-density-stats : Calculate statistics of tile densities
////////////////////////////////////////////////////////////////////////
int32_t cmd_tile_density_stats_mt(int32_t argc, char **argv)
{
    std::string pmtilesf;
    int32_t zoom = -1;             // -1 represents the max zoom level
    int32_t verbose_freq = 100000; // not a parameter
    int32_t num_threads = 0;       // Default: use number of hardware threads

    // output format
    std::string out_tsvf;
    std::string count_field("gn");
    bool compact = false;

    paramList pl;

    BEGIN_LONG_PARAMS(longParameters)
    LONG_PARAM_GROUP("Input options", NULL)
    LONG_STRING_PARAM("in", &pmtilesf, "Input PMTiles file")
    LONG_STRING_PARAM("count", &count_field, "Field name for transcript counts")

    LONG_PARAM_GROUP("Output options", NULL)
    LONG_PARAM("compact", &compact, "Skip writing each tile")
    LONG_STRING_PARAM("out", &out_tsvf, "Output TSV file")

    LONG_PARAM_GROUP("Filtering options", NULL)
    LONG_INT_PARAM("zoom", &zoom, "Zoom level (default: -1 -- maximum zoom level)")
    
    LONG_PARAM_GROUP("Performance options", NULL)
    LONG_INT_PARAM("threads", &num_threads, "Number of threads (default: hardware concurrency)")
    END_LONG_PARAMS();

    pl.Add(new longParams("Available Options", longParameters));
    pl.Read(argc, argv);
    pl.Status();

    notice("Analysis started");

    if (pmtilesf.empty())
    {
        error("Missing required options --in");
    }
    if (out_tsvf.empty())
    {
        error("Missing required options --out");
    }
    
    // Set number of threads if not specified
    if (num_threads <= 0) {
        num_threads = std::thread::hardware_concurrency();
        if (num_threads == 0) {
            num_threads = 4; // Fallback if hardware_concurrency returns 0
        }
    }
    notice("Using %d threads", num_threads);

    // Open a PMTiles file
    pmt_pts pmt(pmtilesf.c_str());

    // Open the header and tile entries
    notice("Reading header and tile entries...");
    if (!pmt.read_header_meta_entries())
    {
        error("This pmtiles file is malformed or incompatible with pmpoints, which requires collection of points in MVT format");
    }

    // Identify tiles that intersect with the region
    if (zoom == -1)
    {
        // select the maximum zoom level
        zoom = pmt.hdr.max_zoom;
        notice("Setting the zoom level to the maximum zoom level: %d", zoom);
    }

    if (zoom < pmt.hdr.min_zoom || zoom > pmt.hdr.max_zoom)
    {
        error("Zoom level %d is unavailable in %s", zoom, pmtilesf.c_str());
    }

    // create/open the output files
    htsFile *tsv_wh = NULL;
    if (out_tsvf.compare(out_tsvf.size() - 3, 3, ".gz", 3) == 0)
    {
        tsv_wh = hts_open(out_tsvf.c_str(), "wz");
    }
    else
    {
        tsv_wh = hts_open(out_tsvf.c_str(), "w");
    }
    hprintf(tsv_wh, "zoom\ttile_x\ttile_y\twidth\tnum_pts\tnum_grids\n");

    // Create the tile queue and results aggregator
    TileQueue tile_queue;
    ThreadSafeResults results(tsv_wh, compact);
    
    // Fill the queue with tiles at the specified zoom level
    int32_t num_tiles = 0;
    for (int32_t i = 0; i < pmt.tile_entries.size(); ++i)
    {
        pmtiles::entry_zxy &entry = pmt.tile_entries[i];
        // skip if the zoom level is not the same
        if (entry.z != zoom)
        {
            continue;
        }
        tile_queue.add_tile(entry);
        num_tiles++;
    }
    notice("Found %d tiles at zoom level %d", num_tiles, zoom);
    
    // Signal completion of queue filling
    tile_queue.finish();
    
    // Create and launch worker threads
    std::vector<std::thread> threads;
    std::atomic<int32_t> processed_count(0);
    
    for (int32_t i = 0; i < num_threads; ++i) {
        threads.emplace_back(process_tiles, std::ref(tile_queue), std::ref(pmt), 
                             std::ref(results), std::ref(count_field), std::ref(processed_count));
    }
    
    // Wait for all threads to complete
    for (auto& thread : threads) {
        thread.join();
    }
    
    // Write summary statistics
    results.write_summary(zoom);
    
    hts_close(tsv_wh);

    notice("Finished writing %llu points in total", results.get_total_points());
    notice("Analysis Finished");

    return 0;
}