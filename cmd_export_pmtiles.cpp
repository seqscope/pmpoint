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

// ---- MLT tile decoding helpers ----


static std::vector<bool> mlt_export_decode_bool_rle(const uint8_t* data, size_t len, size_t count) {
    std::vector<bool> result;
    result.reserve(count);
    size_t i = 0;
    while (i < len && result.size() < count) {
        uint8_t header = data[i++];
        if (header >= 128) {
            size_t run_len = 256 - header;
            for (size_t j = 0; j < run_len && i < len && result.size() < count; ++j, ++i) {
                uint8_t byte = data[i];
                for (int b = 0; b < 8 && result.size() < count; ++b)
                    result.push_back((byte >> b) & 1);
            }
        } else {
            size_t run_len = header + 3;
            if (i < len) {
                uint8_t byte = data[i++];
                for (size_t j = 0; j < run_len && result.size() < count; ++j)
                    for (int b = 0; b < 8 && result.size() < count; ++b)
                        result.push_back((byte >> b) & 1);
            }
        }
    }
    while (result.size() < count) result.push_back(true);
    return result;
}

// Decode an MLT tile and populate pt_dataframe, applying the same
// bounding-box and polygon filters used by the MVT path.
// NOTE: fetch_tile_to_buffer already decompresses, so `tile_buf` is raw MLT bytes.
static void decode_mlt_tile_to_df(const std::string& tile_buf, uint8_t zoom,
                                   int64_t tile_x, int64_t tile_y, pt_dataframe& df,
                                   pmt_utils::pmt_pt_t* p_min_pt,
                                   pmt_utils::pmt_pt_t* p_max_pt,
                                   const std::vector<Polygon*>& polygons) {
    const std::string& buf = tile_buf;
    if (buf.empty()) return;

    double scale_factor = pmt_utils::epsg3857_scale_factor(zoom);
    double offset_x, offset_y;
    pmt_utils::tiletoepsg3857(tile_x, tile_y, zoom, &offset_x, &offset_y);

    const uint8_t* ptr = (const uint8_t*)buf.data();
    const uint8_t* end = ptr + buf.size();
    auto rv = [&]() -> uint64_t {
        uint64_t val = 0; int shift = 0;
        while (ptr < end) {
            uint8_t b = *ptr++;
            val |= (uint64_t)(b & 0x7F) << shift;
            if ((b & 0x80) == 0) break;
            shift += 7;
        }
        return val;
    };

    while (ptr < end) {
        uint64_t layer_len = rv();
        if (layer_len == 0 || ptr >= end) break;
        uint8_t tag = *ptr++;
        if (tag != 1) { ptr += layer_len - 1; continue; }

        // Layer header: name, extent, num_columns
        uint64_t name_len = rv();
        ptr += name_len;
        rv(); // extent (unused here)
        uint64_t num_columns = rv();

        // Read all column metadata first (metadata section)
        struct ColMeta { uint64_t typeCode; std::string name; };
        std::vector<ColMeta> col_metas(num_columns);
        for (uint64_t c = 0; c < num_columns; ++c) {
            col_metas[c].typeCode = rv();
            if (col_metas[c].typeCode >= 10) {
                uint64_t cname_len = rv();
                col_metas[c].name = std::string((char*)ptr, cname_len);
                ptr += cname_len;
            }
        }

        size_t num_attr = num_columns > 0 ? num_columns - 1 : 0;
        // col_types: 2=INT, 1=FLOAT, 0=STRING
        std::vector<int>  col_types(num_attr);
        std::vector<bool> col_nullable(num_attr);
        for (size_t c = 0; c < num_attr; ++c) {
            uint64_t tc = col_metas[c + 1].typeCode;
            col_nullable[c] = (tc % 2 == 1);
            uint64_t base = tc - (tc % 2);
            if      (base >= 20 && base <= 23) col_types[c] = 2;
            else if (base >= 24 && base <= 27) col_types[c] = 1;
            else                               col_types[c] = 0;
        }

        // Data section — GEOMETRY first
        uint64_t geom_num_streams = rv();
        size_t num_features = 0;
        std::vector<double> feat_gx, feat_gy;

        for (uint64_t s = 0; s < geom_num_streams; ++s) {
            if (ptr + 2 > end) break;
            uint8_t h0 = *ptr++;
            uint8_t h1 = *ptr++; (void)h1;
            uint64_t num_vals = rv();
            uint64_t byte_len = rv();
            const uint8_t* sd = ptr;
            ptr += byte_len;
            uint8_t phys = (h0 >> 4) & 0x0F;
            uint8_t dict = h0 & 0x0F;
            if (phys == 1 && dict == 3) { // VERTEX stream
                num_features = (size_t)(num_vals / 2);
                feat_gx.resize(num_features);
                feat_gy.resize(num_features);
                const uint8_t* vp = sd;
                for (size_t i = 0; i < num_features; ++i) {
                    uint64_t zx=0; int sh=0;
                    while(vp<sd+byte_len){uint8_t b=*vp++;zx|=(uint64_t)(b&0x7F)<<sh;sh+=7;if(!(b&0x80))break;}
                    uint64_t zy=0; sh=0;
                    while(vp<sd+byte_len){uint8_t b=*vp++;zy|=(uint64_t)(b&0x7F)<<sh;sh+=7;if(!(b&0x80))break;}
                    int32_t px=(int32_t)((zx>>1)^-(int64_t)(zx&1));
                    int32_t py=(int32_t)((zy>>1)^-(int64_t)(zy&1));
                    feat_gx[i] = offset_x + scale_factor * px;
                    feat_gy[i] = offset_y - scale_factor * py;
                }
            }
        }

        // Decode attribute columns into per-column string arrays
        // Missing values stay as empty strings ("NA" for nullable columns is natural output)
        std::vector<std::vector<std::string>> attr_vals(num_attr,
            std::vector<std::string>(num_features));

        for (size_t c = 0; c < num_attr; ++c) {
            bool nullable = col_nullable[c];
            int ctype = col_types[c];
            bool is_str = (ctype == 0);
            std::vector<bool> present(num_features, true);
            std::vector<uint64_t> str_lens;
            const uint8_t* str_data = nullptr;
            uint64_t str_data_len = 0;

            uint64_t ns = is_str ? rv() : (nullable ? 2 : 1);
            for (uint64_t s = 0; s < ns; ++s) {
                if (ptr + 2 > end) break;
                uint8_t h0 = *ptr++;
                uint8_t h1 = *ptr++; (void)h1;
                uint64_t nv = rv();
                uint64_t bl = rv();
                const uint8_t* sd = ptr;
                ptr += bl;
                uint8_t phys = (h0 >> 4) & 0x0F;

                if (phys == 0) { // PRESENT
                    present = mlt_export_decode_bool_rle(sd, bl, num_features);
                } else if (phys == 1) { // DATA
                    if (ctype == 2) { // INT
                        const uint8_t* dp = sd;
                        size_t fi = 0;
                        for (uint64_t vi = 0; vi < nv; ++vi) {
                            uint64_t zig=0; int sh=0;
                            while(dp<sd+bl){uint8_t b=*dp++;zig|=(uint64_t)(b&0x7F)<<sh;sh+=7;if(!(b&0x80))break;}
                            int64_t val=(int64_t)((zig>>1)^-(int64_t)(zig&1));
                            while(fi<num_features&&!present[fi])++fi;
                            if(fi<num_features) attr_vals[c][fi++]=std::to_string(val);
                        }
                    } else if (ctype == 1) { // FLOAT
                        const uint8_t* dp = sd;
                        size_t fi = 0;
                        for (uint64_t vi = 0; vi < nv; ++vi) {
                            uint32_t bits=(uint32_t)dp[0]|((uint32_t)dp[1]<<8)|
                                          ((uint32_t)dp[2]<<16)|((uint32_t)dp[3]<<24);
                            dp+=4;
                            float fval; memcpy(&fval,&bits,4);
                            while(fi<num_features&&!present[fi])++fi;
                            if(fi<num_features){
                                char tmp[32]; snprintf(tmp,sizeof(tmp),"%.9g",(double)fval);
                                attr_vals[c][fi++]=tmp;
                            }
                        }
                    } else { // STRING DATA
                        str_data=sd; str_data_len=bl; (void)nv;
                    }
                } else if (phys == 3) { // LENGTH (string lengths)
                    const uint8_t* dp = sd;
                    str_lens.reserve(nv);
                    for (uint64_t vi = 0; vi < nv; ++vi) {
                        uint64_t len=0; int sh=0;
                        while(dp<sd+bl){uint8_t b=*dp++;len|=(uint64_t)(b&0x7F)<<sh;sh+=7;if(!(b&0x80))break;}
                        str_lens.push_back(len);
                    }
                }
            }

            if (is_str && str_data && !str_lens.empty()) {
                const uint8_t* dp = str_data;
                size_t fi = 0;
                for (size_t li = 0; li < str_lens.size(); ++li) {
                    while (fi < num_features && !present[fi]) ++fi;
                    if (fi < num_features) {
                        attr_vals[c][fi++] = std::string((char*)dp, str_lens[li]);
                        dp += str_lens[li];
                    }
                }
            }
            (void)str_data_len;
        }

        // Apply filters and add passing features to df
        for (size_t i = 0; i < num_features; ++i) {
            double gx = feat_gx[i], gy = feat_gy[i];
            if (p_min_pt && (gx < p_min_pt->global_x || gy < p_min_pt->global_y)) continue;
            if (p_max_pt && (gx > p_max_pt->global_x || gy > p_max_pt->global_y)) continue;
            if (!polygons.empty()) {
                bool found = false;
                for (auto* p : polygons)
                    if (p->contains_point(gx, gy)) { found = true; break; }
                if (!found) continue;
            }
            pmt_utils::pmt_pt_t pt(zoom, gx, gy);
            df.points.push_back(pt);
            for (size_t c = 0; c < num_attr; ++c) {
                const std::string& v = attr_vals[c][i];
                df.add_feature((int32_t)c, col_metas[c+1].name, v.empty() ? "NA" : v);
            }
        }

        break; // only first layer
    }
}

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

    int xmin_class = std::fpclassify(xmin);
    int xmax_class = std::fpclassify(xmax);
    int ymin_class = std::fpclassify(ymin);
    int ymax_class = std::fpclassify(ymax);

    bool has_boundary = !((xmin_class == FP_INFINITE || xmin_class == FP_NAN) &&
                          (xmax_class == FP_INFINITE || xmax_class == FP_NAN) &&
                          (ymin_class == FP_INFINITE || ymin_class == FP_NAN) &&
                          (ymax_class == FP_INFINITE || ymax_class == FP_NAN));

    if ( has_boundary ) {
        notice("Bounding Box: [(%.3f, %.3f), (%.3f, %.3f)]", xmin, ymin, xmax, ymax);
        notice("min_pt: %u %.3lf %.3lf -- %lu %lu %.3lf %.3lf", min_pt.zoom, min_pt.global_x, min_pt.global_y, min_pt.tile_x, min_pt.tile_y, min_pt.local_x, min_pt.local_y);
        notice("max_pt: %u %.3lf %.3lf -- %lu %lu %.3lf %.3lf", max_pt.zoom, max_pt.global_x, max_pt.global_y, max_pt.tile_x, max_pt.tile_y, max_pt.local_x, max_pt.local_y);
        notice("g_min_pt: %u %.3lf %.3lf -- %lu %lu %.3lf %.3lf", g_min_pt.zoom, g_min_pt.global_x, g_min_pt.global_y, g_min_pt.tile_x, g_min_pt.tile_y, g_min_pt.local_x, g_min_pt.local_y);
        notice("g_max_pt: %u %.3lf %.3lf -- %lu %lu %.3lf %.3lf", g_max_pt.zoom, g_max_pt.global_x, g_max_pt.global_y, g_max_pt.tile_x, g_max_pt.tile_y, g_max_pt.local_x, g_max_pt.local_y);
        // exit(-1);
    }
    else {
        notice("No bounding box is set; all tiles at zoom level %d will be considered", zoom);
    }

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

    bool tsv_hdr_written = false;
    uint64_t n_written = 0;
    std::vector<Polygon *> tile_polygons;
    std::string tile_buffer;
    uint64_t n_skipped_tiles = 0;
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

        //notice("Checking tile %d/%d/%d, bbox = [(%.5lf, %.5lf) (%.5lf, %.5lf)]", pmt.tile_entries[i].z, pmt.tile_entries[i].x, pmt.tile_entries[i].y, tile_min_pt.x, tile_min_pt.y, tile_max_pt.x, tile_max_pt.y);

        // check if the boundary was set
        // note that the y-axis is inverted, so min/max is swapped in y when comparing the tiles
        //notice("has_boundary = %d", has_boundary);
        if (has_boundary)
        {
            if (entry.x < min_pt.tile_x || entry.x > max_pt.tile_x || entry.y < max_pt.tile_y || entry.y > min_pt.tile_y)
            {
                //notice("Skipping (%lu, %lu) as it is outside the rectangle defined by (%lu, %lu) -- (%lu, %lu)",
                //       entry.x, entry.y, min_pt.tile_x, min_pt.tile_y, max_pt.tile_x, max_pt.tile_y);
                n_skipped_tiles++;
                if ( n_skipped_tiles % 100 == 1 ) {
                    notice("Skipped %lu/%d tiles of total %zu tiles...", n_skipped_tiles, i+1, pmt.tile_entries.size());
                }
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
                //notice("Skipping (%lu, %lu) as it does not intersect with any of the polygons", entry.x, entry.y);
                n_skipped_tiles++;
                if ( n_skipped_tiles % 100 == 1 ) {
                    notice("Skipped %lu/%d tiles of total %zu tiles...", n_skipped_tiles, i+1, pmt.tile_entries.size());
                }
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
        pmt.fetch_tile_to_buffer(entry.z, entry.x, entry.y, tile_buffer);
        if (pmt.hdr.tile_type == 0x06) {
            decode_mlt_tile_to_df(tile_buffer, entry.z, entry.x, entry.y, df,
                                  mvtfilt.p_min_pt, mvtfilt.p_max_pt, mvtfilt.polygons);
        } else {
            mvtfilt.decode_points_df(tile_buffer, entry.z, entry.x, entry.y, df);
        }

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
