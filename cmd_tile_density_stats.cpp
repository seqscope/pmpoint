#include "pmpoint.h"
#include "qgenlib/tsv_reader.h"
#include "qgenlib/qgen_error.h"

#include <vector>
#include <string>
#include <cstring>
#include <climits>
#include <map>

#include "pmt_pts.h"
#include "pmt_utils.h"
#include "polygon.h"
#include "mvt_pts.h"
#include "htslib/hts.h"
#include "ext/nlohmann/json.hpp"

#define MAX_BITS 20
// calculate density statistics for each tile
// for each resolution bit, 1,2,4,8,...,4096
// keep track of the number of "spots" and the number of "points" per spot
class tile_density_stats {
public:
    int32_t maxbits;
    std::vector< std::map<uint64_t, uint64_t> > res2bin2pts;
    tile_density_stats(int32_t _maxbits = MAX_BITS) {
        maxbits = _maxbits;
        res2bin2pts.resize(maxbits);
    }

    void add_xycnt(int32_t x, int32_t y, int32_t cnt) {
        for(int32_t i=0; i< maxbits; ++i) {
            int32_t bsize = 1 << i;
            uint64_t bin = ((((uint64_t)x / bsize) << 32) | ((uint64_t)y / bsize));
            res2bin2pts[i][bin] += cnt;
        }
    }

    // summary stats per tile
    std::map<uint64_t, uint64_t> pts2nbins;
    void compute_pts2nbins(int32_t res) {
        pts2nbins.clear();
        std::map<uint64_t, uint64_t> &bin2pts = res2bin2pts[res];
        for(auto it=bin2pts.begin(); it!=bin2pts.end(); ++it) {
            uint64_t pts = it->second;
            ++pts2nbins[pts];
        }
    }
};

/////////////////////////////////////////////////////////////////////////
// tile-density-stats : Calculate statistics of tile densities
////////////////////////////////////////////////////////////////////////
int32_t cmd_tile_density_stats(int32_t argc, char **argv)
{
    std::string pmtilesf;
    int32_t zoom = -1;             // -1 represents the max zoom level
    int32_t verbose_freq = 100000; // not a parameter

    // output format
    std::string out_tsvf;
    std::string count_field("gn");
    std::string feature_field("gene");
    bool compact = false;

    paramList pl;

    BEGIN_LONG_PARAMS(longParameters)
    LONG_PARAM_GROUP("Input options", NULL)
    LONG_STRING_PARAM("in", &pmtilesf, "Input PMTiles file")
    LONG_STRING_PARAM("count", &count_field, "Field name for transcript counts")
    LONG_STRING_PARAM("feature", &feature_field, "Field name for feature name")

    LONG_PARAM_GROUP("Output options", NULL)
    LONG_PARAM("compact", &compact, "Skip writing each tile")
    LONG_STRING_PARAM("out", &out_tsvf, "Output TSV file")

    LONG_PARAM_GROUP("Filtering options", NULL)
    LONG_INT_PARAM("zoom", &zoom, "Zoom level (default: -1 -- maximum zoom level)")
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

    mvt_pts mvt;

    uint64_t n_total = 0;
    std::vector< std::map<uint64_t, uint64_t> > res2pts2nbins(MAX_BITS);
    std::string tile_buffer;
    for (int32_t i = 0; i < pmt.tile_entries.size(); ++i)
    {
        pmtiles::entry_zxy &entry = pmt.tile_entries[i];

        // skip if the zoom level is not the same
        if (entry.z != zoom)
        {
            continue;
        }
        // notice("Fetching tile %d/%d/%d that intersects with the region", entry.z, entry.x, entry.y);
        //pmt.fetch_tile(entry.z, entry.x, entry.y);
        pmt.fetch_tile_to_buffer(entry.z, entry.x, entry.y, tile_buffer);
        std::vector<int32_t> xs, ys, cnts;
        std::vector<std::string> features;
        //int32_t n_pts = mvt.decode_points_xycnt(pmt.tile_data_str, count_field, xs, ys, cnts);
        int32_t n_pts = mvt.decode_points_xycnt_feature(tile_buffer, count_field, feature_field, xs, ys, cnts, features);

        tile_density_stats tds;
        for (int32_t j = 0; j < n_pts; ++j)
        {
            tds.add_xycnt(xs[j], ys[j], cnts[j]);
        }

        // compute statistics
        for(int32_t i=0; i<MAX_BITS; ++i) {
            tds.compute_pts2nbins(i);
            std::map<uint64_t, uint64_t> &pts2nbins = tds.pts2nbins;
            for(auto it=pts2nbins.begin(); it!=pts2nbins.end(); ++it) {
                uint64_t pts = it->first;
                uint64_t nbins = it->second;
                if ( !compact ) {
                    hprintf(tsv_wh, "%d\t%d\t%d\t%d\t%llu\t%llu\n", entry.z, entry.x, entry.y, 1<<i, pts, nbins);
                }
                res2pts2nbins[i][pts] += nbins;
            }
        }

        notice("Tile #%d/%zu %d/%d/%d , %d points", i, pmt.tile_entries.size(), entry.z, entry.x, entry.y, n_pts);
        n_total += n_pts;
    }

    for(int32_t i=0; i<MAX_BITS; ++i) {
        std::map<uint64_t, uint64_t> &pts2nbins = res2pts2nbins[i];
        for(auto it=pts2nbins.begin(); it!=pts2nbins.end(); ++it) {
            uint64_t pts = it->first;
            uint64_t nbins = it->second;
            hprintf(tsv_wh, "%d\tALL\tALL\t%d\t%llu\t%llu\n", zoom, 1<<i, pts, nbins);
        }
    }

    hts_close(tsv_wh);

    notice("Finished writing %llu points in total", n_total);

    notice("Analysis Finished");

    return 0;
}
