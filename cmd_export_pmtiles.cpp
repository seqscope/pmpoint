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

/////////////////////////////////////////////////////////////////////////
// extract : Export points from a PMTiles file to a TSV file
////////////////////////////////////////////////////////////////////////
int32_t cmd_export_pmtiles(int32_t argc, char **argv)
{
    std::string pmtilesf;
    int32_t zoom = -1;             // -1 represents the max zoom level available
    int32_t verbose_freq = 100000; // not a parameter

    // parameter for region-based filtering
    double xmin = -std::numeric_limits<double>::infinity();
    double xmax = std::numeric_limits<double>::infinity();
    double ymin = -std::numeric_limits<double>::infinity();
    double ymax = std::numeric_limits<double>::infinity();

    // parameter for geojson-based filtering
    std::string geojsonf;

    // output format
    std::string out_tsvf;
    std::string out_jsonf;

    int32_t precision = 3; // precision of the output

    paramList pl;

    BEGIN_LONG_PARAMS(longParameters)
    LONG_PARAM_GROUP("Input options", NULL)
    LONG_STRING_PARAM("in", &pmtilesf, "Input PMTiles file")

    LONG_PARAM_GROUP("Output options", NULL)
    LONG_STRING_PARAM("out-tsv", &out_tsvf, "Output TSV file")
    LONG_STRING_PARAM("out-json", &out_jsonf, "Output JSON file")

    LONG_PARAM_GROUP("Filtering options", NULL)
    LONG_INT_PARAM("zoom", &zoom, "Zoom level (default: -1 -- maximum zoom level)")
    LONG_DOUBLE_PARAM("xmin", &xmin, "Minimum x-axis value")
    LONG_DOUBLE_PARAM("xmax", &xmax, "Maximum x-axis value")
    LONG_DOUBLE_PARAM("ymin", &ymin, "Minimum y-axis value")
    LONG_DOUBLE_PARAM("ymax", &ymax, "Maximum y-axis value")
    LONG_STRING_PARAM("polygon", &geojsonf, "GeoJSON file (in EPSG:3857) for polygon-based filtering")

    LONG_PARAM_GROUP("Additional options", NULL)
    LONG_INT_PARAM("precision", &precision, "Precision of the output of X/Y coordinates (default: 3)")
    END_LONG_PARAMS();

    pl.Add(new longParams("Available Options", longParameters));
    pl.Read(argc, argv);
    pl.Status();

    notice("Analysis started");

    if (pmtilesf.empty())
    {
        error("Missing required options --in");
    }
    if (out_tsvf.empty() && out_jsonf.empty())
    {
        error("Missing required options --out-tsv or --out-json (at least 1 required)");
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

    // convert the input coordinates to tile space
    pmt_utils::pmt_pt_t min_pt(zoom, xmin, ymin);
    pmt_utils::pmt_pt_t max_pt(zoom, xmax, ymax);
    pmt_utils::pmt_pt_t g_min_pt(zoom, min_pt.tile_x, min_pt.tile_y, min_pt.local_x, min_pt.local_y);
    pmt_utils::pmt_pt_t g_max_pt(zoom, max_pt.tile_x, max_pt.tile_y, max_pt.local_x, max_pt.local_y);

    notice("xmin = %.3f, ymin = %.3f", xmin, ymin);
    notice("xmax = %.3f, ymax = %.3f", xmax, ymax);
    notice("min_pt: %u %.3lf %.3lf -- %lu %lu %.3lf %.3lf", min_pt.zoom, min_pt.global_x, min_pt.global_y, min_pt.tile_x, min_pt.tile_y, min_pt.local_x, min_pt.local_y);
    notice("max_pt: %u %.3lf %.3lf -- %lu %lu %.3lf %.3lf", max_pt.zoom, max_pt.global_x, max_pt.global_y, max_pt.tile_x, max_pt.tile_y, max_pt.local_x, max_pt.local_y);
    notice("g_min_pt: %u %.3lf %.3lf -- %lu %lu %.3lf %.3lf", g_min_pt.zoom, g_min_pt.global_x, g_min_pt.global_y, g_min_pt.tile_x, g_min_pt.tile_y, g_min_pt.local_x, g_min_pt.local_y);
    notice("g_max_pt: %u %.3lf %.3lf -- %lu %lu %.3lf %.3lf", g_max_pt.zoom, g_max_pt.global_x, g_max_pt.global_y, g_max_pt.tile_x, g_max_pt.tile_y, g_max_pt.local_x, g_max_pt.local_y);
    // exit(-1);

    // load geojson
    std::vector<Polygon> polygons;
    if (!geojsonf.empty())
    {
        int32_t npolygons = load_polygons_from_geojson(geojsonf.c_str(), polygons);
    }
    // and compute bounding boxes for each polygon

    std::vector<Rectangle> bounding_boxes;
    for (int32_t i = 0; i < polygons.size(); ++i)
    {
        bounding_boxes.push_back(polygons[i].get_bounding_box());
        Rectangle& r = bounding_boxes.back();
        notice("BBox %d = ll(%lf, %lf) - ur(%lf,%lf)", i, r.p_min.x, r.p_min.y, r.p_max.x, r.p_max.y);
    }

    // create/open the output files
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
    if (!out_jsonf.empty())
    {
        if (out_jsonf.compare(out_jsonf.size() - 3, 3, ".gz", 3) == 0)
        {
            json_wh = hts_open(out_jsonf.c_str(), "wz");
        }
        else
        {
            json_wh = hts_open(out_jsonf.c_str(), "w");
        }
        hprintf(json_wh, "{\n");
    }

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
    pt_dataframe df;
    mvt_pts_filt mvtfilt(&df);

    int xmin_class = std::fpclassify(xmin);
    int xmax_class = std::fpclassify(xmax);
    int ymin_class = std::fpclassify(ymin);
    int ymax_class = std::fpclassify(ymax);

    bool has_boundary = !((xmin_class == FP_INFINITE || xmin_class == FP_NAN) &&
                          (xmax_class == FP_INFINITE || xmax_class == FP_NAN) &&
                          (ymin_class == FP_INFINITE || ymin_class == FP_NAN) &&
                          (ymax_class == FP_INFINITE || ymax_class == FP_NAN));

    bool tsv_hdr_written = false;
    uint64_t n_written = 0;
    std::vector<Polygon *> tile_polygons;
    std::string tile_buffer;
    for (int32_t i = 0; i < pmt.tile_entries.size(); ++i)
    {
        pmtiles::entry_zxy &entry = pmt.tile_entries[i];

        // skip if the zoom level is not the same
        if (entry.z != zoom)
        {
            continue;
        }

        // get the global coordinates of the tile. Note that the y-axis is inverted
        point_t tile_min_pt(0,0), tile_max_pt(0,0);
        pmt_utils::tiletoepsg3857(entry.x, entry.y, entry.z, &tile_min_pt.x, &tile_max_pt.y);
        pmt_utils::tiletoepsg3857(entry.x+1, entry.y+1, entry.z, &tile_max_pt.x, &tile_min_pt.y);
        Rectangle tile_bbox(tile_min_pt.x, tile_min_pt.y, tile_max_pt.x, tile_max_pt.y);

        notice("Checking tile %d/%d/%d, bbox = [(%.5lf, %.5lf) (%.5lf, %.5lf)]", pmt.tile_entries[i].z, pmt.tile_entries[i].x, pmt.tile_entries[i].y, tile_min_pt.x, tile_min_pt.y, tile_max_pt.x, tile_max_pt.y);

        // check if the boundary was set
        // note that the y-axis is inverted, so min/max is swapped in y when comparing the tiles
        //notice("has_boundary = %d", has_boundary);
        if (has_boundary)
        {
            if (entry.x < min_pt.tile_x || entry.x > max_pt.tile_x || entry.y < max_pt.tile_y || entry.y > min_pt.tile_y)
            {
                notice("Skipping (%lu, %lu) as it is outside the rectangle defined by (%lu, %lu) -- (%lu, %lu)",
                       entry.x, entry.y, min_pt.tile_x, min_pt.tile_y, max_pt.tile_x, max_pt.tile_y);
                continue;
            }
            else
            {
                notice("Considering (%lu, %lu) as it is inside the rectangle defined by (%lu, %lu) -- (%lu, %lu)",
                       entry.x, entry.y, min_pt.tile_x, min_pt.tile_y, max_pt.tile_x, max_pt.tile_y);
            }
            // if the min/max point is located at the tile, boundary, then we need to check the points
            if (entry.x == min_pt.tile_x || entry.y == min_pt.tile_y)
            {
                mvtfilt.set_min_filt(&min_pt);
            }
            else
            {
                // mvtfilt.set_min_filt(&min_pt);
                mvtfilt.set_min_filt(NULL);
            }
            if (entry.x == max_pt.tile_x || entry.y == max_pt.tile_y)
            {
                mvtfilt.set_max_filt(&max_pt);
            }
            else
            {
                // mvtfilt.set_max_filt(&max_pt);
                mvtfilt.set_max_filt(NULL);
            }
        }

        // if polygons exists
        if (polygons.size() > 0)
        {
            // check if the polygon is contained in the bounding box
            // check if the bounding 
            tile_polygons.clear();

            //notice("global_min_pt = (%lg, %lg)", min_pt.global_x, min_pt.global_y);
            //notice("global_max_pt = (%lg, %lg)", max_pt.global_x, max_pt.global_y);
            //notice("tile_min_pt = (%lg, %lg)", tile_bbox.p_min.x, tile_bbox.p_min.y);
            //notice("tile_max_pt = (%lg, %lg)", tile_bbox.p_max.x, tile_bbox.p_max.y);

            for (int32_t j = 0; j < bounding_boxes.size(); ++j)
            {
                // if any of the bounding boxes of the polygon intersects with the tile,
                // then we need to include the polygon
                if ( bounding_boxes[j].intersects_rectangle(tile_bbox) )
                {
                    tile_polygons.push_back(&polygons[j]);
                }
            }
            if (tile_polygons.size() == 0)
            {
                notice("Skipping (%lu, %lu) as it does not intersect with any of the polygons", entry.x, entry.y);
                continue;
            }
            mvtfilt.set_polygon_filt(tile_polygons);
        }
        else
        {
            // do not set any filter, automatic pass
            // mvtfilt.set_polygon_filt(tile_polygons);
        }

        // notice("Fetching tile %d/%d/%d that intersects with the region", entry.z, entry.x, entry.y);
        //pmt.fetch_tile(entry.z, entry.x, entry.y);
        //mvtfilt.decode_points(pmt.tile_data_str, entry.z, entry.x, entry.y);
        pmt.fetch_tile_to_buffer(entry.z, entry.x, entry.y, tile_buffer);
        mvtfilt.decode_points_df(tile_buffer, entry.z, entry.x, entry.y, df);

        if (tsv_wh != NULL)
        {
            if (!tsv_hdr_written)
            {
                if (df.points.size() > 0)
                {
                    if (df.points.size() > 0)
                    {
                        hprintf(tsv_wh, "X\tY");
                        for (int32_t i = 0; i < df.feature_names.size(); ++i)
                        {
                            hprintf(tsv_wh, "\t%s", df.feature_names[i].c_str());
                        }
                        hprintf(tsv_wh, "\n");
                        tsv_hdr_written = true;
                    }
                }
            }
            for (int32_t i = 0; i < df.points.size(); ++i)
            {
                hprintf(tsv_wh, "%.*f\t%.*f", precision, df.points[i].global_x, precision, df.points[i].global_y);
                for (int32_t j = 0; j < df.feature_matrix.size(); ++j)
                {
                    hprintf(tsv_wh, "\t%s", df.feature_matrix[j][i].c_str());
                }
                hprintf(tsv_wh, "\n");
                if ((n_written + i + 1) % verbose_freq == 0)
                {
                    notice("Writing %llu points to %s", n_written + i + 1, out_tsvf.c_str());
                }
            }
        }
        if (json_wh != NULL)
        {
            if (df.points.size() > 0)
            {
                for (int32_t i = 0; i < df.points.size(); ++i)
                {
                    hprintf(json_wh, "{\"type\":\"Feature\",\"properties\": {");
                    for (int32_t j = 0; j < df.feature_matrix.size(); ++j)
                    {
                        hprintf(json_wh, "\"%s\":\"%s\"", df.feature_names[j].c_str(), df.feature_matrix[j][i].c_str());
                        if (j < df.feature_matrix.size() - 1)
                        {
                            hprintf(json_wh, ",");
                        }
                    }
                    hprintf(json_wh, "},\"geometry\":{\"type\":\"Point\",\"coordinates\":[%.*f,%.*f]}}\n", precision, df.points[i].global_x, precision, df.points[i].global_y);
                    if ((n_written + i + 1) % verbose_freq == 0)
                    {
                        notice("Writing %llu points to %s", n_written + i + 1, out_jsonf.c_str());
                    }
                }
            }
        }
        n_written += df.points.size();
        notice("Finished writing %zu additional points -- %llu total", df.points.size(), n_written);
        df.clear_values();
    }

    if (json_wh != NULL)
    {
        hprintf(json_wh, "}\n");
        hts_close(json_wh);
    }
    if (tsv_wh != NULL)
    {
        hts_close(tsv_wh);
    }

    notice("Finished writing %llu points in total", n_written);

    notice("Analysis Finished");

    return 0;
}
