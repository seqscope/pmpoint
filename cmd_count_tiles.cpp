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
#include "ext/PMTiles/pmtiles.hpp"

/////////////////////////////////////////////////////////////////////////
// count-tiles : Count number of points in each tiles from a PMTiles file
////////////////////////////////////////////////////////////////////////
int32_t cmd_count_tiles(int32_t argc, char **argv)
{
    std::string pmtilesf;
    int32_t zoom = -1;             // -1 represents the max zoom level available
    int32_t verbose_freq = 100000; // not a parameter

    // output format
    std::string out_tsvf;
    std::string out_jsonf;

    paramList pl;

    BEGIN_LONG_PARAMS(longParameters)
    LONG_PARAM_GROUP("Input options", NULL)
    LONG_STRING_PARAM("in", &pmtilesf, "Input PMTiles file")

    LONG_PARAM_GROUP("Output options", NULL)
    LONG_STRING_PARAM("out-tsv", &out_tsvf, "Output TSV file")
//    LONG_STRING_PARAM("out-json", &out_jsonf, "Output JSON file")

    LONG_PARAM_GROUP("Filtering options", NULL)
    LONG_INT_PARAM("zoom", &zoom, "Zoom level (default: -1 -- all zoom levels)")

    END_LONG_PARAMS();

    pl.Add(new longParams("Available Options", longParameters));
    pl.Read(argc, argv);
    pl.Status();

    notice("Analysis started");

    if (pmtilesf.empty())
    {
        error("Missing required option --in");
    }
    if (out_tsvf.empty())
    {
        error("Missing required option --out-tsv");
    }
    // if (out_tsvf.empty() && out_jsonf.empty())
    // {
    //     error("Missing required options --out-tsv or --out-json (at least 1 required)");
    // }

    // Open a PMTiles file
    pmt_pts pmt(pmtilesf.c_str());

    // Open the header and tile entries
    notice("Reading header and tile entries...");
    if (!pmt.read_header_meta_entries())
    {
        error("This pmtiles file is malformed or incompatible with pmpoints, which requires collection of points in MVT format");
    }

    //create/open the output files
    htsFile *tsv_wh = NULL;
    htsFile *json_wh = NULL;
    if (!out_tsvf.empty())
    {
        if (out_tsvf.compare(out_tsvf.size() - 3, 3, ".gz", 3) == 0)
        {
            tsv_wh = hts_open(out_tsvf.c_str(), "wz");
        }
        else
        {
            tsv_wh = hts_open(out_tsvf.c_str(), "w");
        }
    }
    // if (!out_jsonf.empty())
    // {
    //     if (out_jsonf.compare(out_jsonf.size() - 3, 3, ".gz", 3) == 0)
    //     {
    //         json_wh = hts_open(out_jsonf.c_str(), "wz");
    //     }
    //     else
    //     {
    //         json_wh = hts_open(out_jsonf.c_str(), "w");
    //     }
    //     hprintf(json_wh, "{\n");
    // }

    // for each tile
    // check the following
    // (a) whether the tile is inside the bounding box defined by xmin/xmax/ymin/ymax
    //     [E] No boundary - pass
    //     [Y] Pass
    //     [N] Skip; Do not consider including the tile
    // (b) (aEY-only) whether the tile overlaps with the bounding box of any of the polygon
    //     [E] No polygon (empty)
    //     [Y] Overlaps with some polygon(s)
    //     [N] Skip; Do not consider including the tile
    // (c) (aEY-bY only) identify overlapping polygons and examine individual points to include
    // (d) (aEY-bE only) pass
    mvt_pts mvt;

    // a vector containing [tile id] -> [count] for each zoom level
    std::vector<std::map<uint64_t, uint64_t> > zoom2tile2count(pmt.hdr.max_zoom + 1);
    std::vector<std::map<uint64_t, uint64_t> > zoom2tile2sum(pmt.hdr.max_zoom + 1);

    
    for (int32_t i = 0; i < pmt.tile_entries.size(); ++i)
    {
        pmtiles::entry_zxy &entry = pmt.tile_entries[i];

        if ( zoom >= 0 && entry.z != zoom )
        {
            continue;
        }

        if ( i % 100 == 0 ) {
            notice("Processing %d-th tile %d/%d/%d", i+1, entry.z, entry.x, entry.y);
        }

        // notice("Fetching tile %d/%d/%d that intersects with the region", entry.z, entry.x, entry.y);
        pmt.fetch_tile(entry.z, entry.x, entry.y);
        uint64_t n_points = mvt.count_points(pmt.tile_data_str, entry.z, entry.x, entry.y);
        //notice("Tile %d/%d/%d has %llu points", entry.z, entry.x, entry.y, n_points);

        uint64_t xy = (((uint64_t)entry.x) << 32 | (uint64_t)entry.y);
        zoom2tile2count[entry.z][xy] = n_points;
    }

    // build hierarchical counts
    if ( zoom < 0 ) {
        for(auto it=zoom2tile2count[pmt.hdr.max_zoom].begin(); it!=zoom2tile2count[pmt.hdr.max_zoom].end(); ++it) {
            uint64_t xy = it->first;
            uint64_t x = (xy >> 32) & 0xFFFFFFFF;
            uint64_t y = xy & 0xFFFFFFFF;
            uint64_t count = it->second;

            // add the count to all parent tiles
            for(int i=pmt.hdr.max_zoom; i >= pmt.hdr.min_zoom; --i) {
                uint64_t newxy = (((uint64_t)x) << 32 | (uint64_t)y);
                zoom2tile2sum[i][newxy] += count;
                x >>= 1;
                y >>= 1;
            }
        }
    }
    else {
        zoom2tile2sum = zoom2tile2count;
    }

    // print the count information
    if (!out_tsvf.empty())
    {
        if ( zoom < 0 ) {
            hprintf(tsv_wh, "zoom\tx\ty\ttile_id\ttotal_count\ttile_count\tfrac_in_tile\n");
        }
        else {
            hprintf(tsv_wh, "zoom\tx\ty\ttile_id\ttile_count\n");
        }
    }
    for(int32_t i=pmt.hdr.min_zoom; i <= pmt.hdr.max_zoom; ++i) {
        for(auto it=zoom2tile2sum[i].begin(); it!=zoom2tile2sum[i].end(); ++it) {
            uint64_t xy = it->first;
            uint64_t x = (xy >> 32) & 0xFFFFFFFF;
            uint64_t y = xy & 0xFFFFFFFF;
            uint64_t sum = it->second;
            uint64_t count = zoom2tile2count[i][xy];
            uint64_t tile_id = pmtiles::zxy_to_tileid(i, (uint32_t)x, (uint32_t)y);

            if (!out_tsvf.empty())
            {
                if ( zoom < 0 ) {
                    hprintf(tsv_wh, "%d\t%llu\t%llu\t%llu\t%llu\t%llu\t%.5g\n", i, x, y, tile_id, sum, count, (double)count/(double)sum);
                }
                else {
                    hprintf(tsv_wh, "%d\t%llu\t%llu\t%llu\t%llu\n", i, x, y, tile_id, count);
                }
            }

            //notice("%d\t%d\t%d\t%llu\t%llu\t%.5f", i, x, y, sum, count, (double)count/(double)sum); 
        }
    }

    // if (json_wh != NULL)
    // {
    //     hprintf(json_wh, "}\n");
    //     hts_close(json_wh);
    // }
    if (tsv_wh != NULL)
    {
        hts_close(tsv_wh);
    }

    notice("Analysis Finished");

    return 0;
}
