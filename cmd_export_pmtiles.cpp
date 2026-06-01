#include "pmpoint.h"
#include "qgenlib/tsv_reader.h"
#include "qgenlib/qgen_error.h"

#include <vector>
#include <string>
#include <cstring>
#include <climits>
#include <map>
#include <set>
#include <unordered_map>
#include <algorithm>
#include <atomic>
#include <mutex>
#include <thread>
#include <memory>

#include "pmt_pts.h"
#include "pmt_utils.h"
#include "polygon.h"
#include "mvt_pts.h"
#include "flex_io.h"
#include "htslib/hts.h"
#include "ext/nlohmann/json.hpp"

// ===========================================================================
// Overview
// ---------------------------------------------------------------------------
// `pmpoint export` extracts points from a PMTiles file into TSV (and/or JSON).
// This implementation parallelizes the expensive work (network fetch + decode
// + formatting) across worker threads, each with its own I/O reader so S3
// range requests proceed concurrently, and serializes only the final write.
//
// Output columns: when the PMTiles metadata declares the attribute schema
// (tilestats.layers[].attributes or vector_layers[].fields), it is treated as
// authoritative and used as a fixed rectangular schema (missing -> NA). When
// no schema is available, points are streamed to a temporary binary file while
// the column set is discovered, then expanded to a consistent TSV in a single
// pass over the temp file (no re-fetch / re-decode of the source).
// ===========================================================================

namespace {

// ---- Per-tile work item (filters resolved up front, single-threaded) ------
struct TileTask {
    pmtiles::entry_zxy entry{0, 0, 0, 0, 0};
    pmt_utils::pmt_pt_t* min_filt = nullptr;          // bbox lower bound, or null
    pmt_utils::pmt_pt_t* max_filt = nullptr;          // bbox upper bound, or null
    std::vector<Polygon*> polygons;                   // polygons intersecting this tile
};

// ---- Compact per-tile decode result (names interned at tile scope) --------
// CSR-style attribute storage: for point i, its (col,value) entries are
// vals[starts[i] .. starts[i+1]) with local column index cols[k] into names.
struct TileBatch {
    std::vector<double> xs, ys;             // accepted point coordinates
    std::vector<std::string> names;         // distinct attribute names in this tile
    std::vector<int32_t> starts{0};         // size npoints+1
    std::vector<int32_t> cols;              // local name index per attribute entry
    std::vector<std::string> vals;          // value per attribute entry

    size_t npoints() const { return xs.size(); }
    void clear() {
        xs.clear(); ys.clear(); names.clear();
        starts.assign(1, 0); cols.clear(); vals.clear();
    }
};

// ---- Output schema / formatting context (read-only across threads) --------
struct ExportCtx {
    bool want_tsv = false;
    bool want_json = false;
    bool schema_known = false;                        // fixed-column TSV path
    int32_t precision = 3;
    std::string missing_value;
    std::vector<std::string> col_names;               // display order (TSV header)
    std::unordered_map<std::string, int32_t> name_to_global;  // name -> column position
    size_t ncols() const { return col_names.size(); }
};

// Collect declared attribute names from metadata, tolerating either shape.
bool extract_schema(const nlohmann::json& meta, std::set<std::string>& out) {
    // shape 1: tilestats.layers[*].attributes[*].attribute
    if (meta.contains("tilestats") && meta["tilestats"].is_object()) {
        const auto& ts = meta["tilestats"];
        if (ts.contains("layers") && ts["layers"].is_array()) {
            for (const auto& layer : ts["layers"]) {
                if (layer.contains("attributes") && layer["attributes"].is_array()) {
                    for (const auto& a : layer["attributes"]) {
                        if (a.contains("attribute") && a["attribute"].is_string())
                            out.insert(a["attribute"].get<std::string>());
                    }
                }
            }
        }
    }
    // shape 2: vector_layers[*].fields (object keyed by attribute name)
    if (meta.contains("vector_layers") && meta["vector_layers"].is_array()) {
        for (const auto& vl : meta["vector_layers"]) {
            if (vl.contains("fields") && vl["fields"].is_object()) {
                for (auto it = vl["fields"].begin(); it != vl["fields"].end(); ++it)
                    out.insert(it.key());
            }
        }
    }
    return !out.empty();
}

// Build display column order: priority features first (in given order, if present),
// then the remaining names sorted alphabetically. Mirrors the historical behavior.
void build_column_order(const std::set<std::string>& names,
                        const std::vector<std::string>& priority,
                        std::vector<std::string>& ordered,
                        std::unordered_map<std::string, int32_t>& name_to_idx) {
    ordered.clear();
    name_to_idx.clear();
    std::set<std::string> used;
    for (const auto& p : priority) {
        if (names.count(p) && !used.count(p)) {
            name_to_idx[p] = (int32_t)ordered.size();
            ordered.push_back(p);
            used.insert(p);
        }
    }
    std::vector<std::string> rest;
    for (const auto& n : names)
        if (!used.count(n)) rest.push_back(n);
    std::sort(rest.begin(), rest.end());
    for (const auto& n : rest) {
        name_to_idx[n] = (int32_t)ordered.size();
        ordered.push_back(n);
    }
}

// ---- MLT (MapLibre Tile) helpers ------------------------------------------
std::vector<bool> mlt_decode_bool_rle(const uint8_t* data, size_t len, size_t count) {
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

inline bool point_passes(double gx, double gy, const TileTask& task) {
    if (task.min_filt && (gx < task.min_filt->global_x || gy < task.min_filt->global_y)) return false;
    if (task.max_filt && (gx > task.max_filt->global_x || gy > task.max_filt->global_y)) return false;
    if (!task.polygons.empty()) {
        for (auto* p : task.polygons)
            if (p->contains_point(gx, gy)) return true;
        return false;
    }
    return true;
}

// Decode an MVT tile into a TileBatch, applying the task filters.
void decode_mvt_batch(const std::string& buf, const TileTask& task, TileBatch& b) {
    b.clear();
    mapbox::vector_tile::buffer tile(buf);
    double scale = pmt_utils::epsg3857_scale_factor(task.entry.z);
    double off_x, off_y;
    pmt_utils::tiletoepsg3857(task.entry.x, task.entry.y, task.entry.z, &off_x, &off_y);

    std::unordered_map<std::string, int32_t> local_idx;
    print_value pv;

    for (auto const& name : tile.layerNames()) {
        const mapbox::vector_tile::layer layer = tile.getLayer(name);
        std::size_t fc = layer.featureCount();
        for (std::size_t i = 0; i < fc; ++i) {
            auto const feature = mapbox::vector_tile::feature(layer.getFeature(i), layer);
            if (int(feature.getType()) != 1)
                error("Only points are supported in export");
            auto geom = feature.getGeometries<mapbox::vector_tile::points_arrays_type>(1.0);
            if (geom.size() != 1)
                error("Only single point per feature is supported in export");

            double gx = off_x + scale * geom[0][0].x;
            double gy = off_y - scale * geom[0][0].y;
            if (!point_passes(gx, gy, task)) continue;

            b.xs.push_back(gx);
            b.ys.push_back(gy);
            auto props = feature.getProperties();
            for (auto const& prop : props) {
                int32_t li;
                auto it = local_idx.find(prop.first);
                if (it == local_idx.end()) {
                    li = (int32_t)b.names.size();
                    b.names.push_back(prop.first);
                    local_idx.emplace(prop.first, li);
                } else {
                    li = it->second;
                }
                b.cols.push_back(li);
                b.vals.push_back(mapbox::util::apply_visitor(pv, prop.second));
            }
            b.starts.push_back((int32_t)b.cols.size());
        }
    }
}

// Decode an MLT tile into a TileBatch, applying the task filters.
// NOTE: the tile buffer is already decompressed (raw MLT bytes).
void decode_mlt_batch(const std::string& buf, const TileTask& task, TileBatch& b) {
    b.clear();
    if (buf.empty()) return;

    double scale = pmt_utils::epsg3857_scale_factor(task.entry.z);
    double off_x, off_y;
    pmt_utils::tiletoepsg3857(task.entry.x, task.entry.y, task.entry.z, &off_x, &off_y);

    const uint8_t* ptr = (const uint8_t*)buf.data();
    const uint8_t* end = ptr + buf.size();
    auto rv = [&]() -> uint64_t {
        uint64_t val = 0; int shift = 0;
        while (ptr < end) {
            uint8_t bb = *ptr++;
            val |= (uint64_t)(bb & 0x7F) << shift;
            if ((bb & 0x80) == 0) break;
            shift += 7;
        }
        return val;
    };

    while (ptr < end) {
        uint64_t layer_len = rv();
        if (layer_len == 0 || ptr >= end) break;
        uint8_t tag = *ptr++;
        if (tag != 1) { ptr += layer_len - 1; continue; }

        uint64_t name_len = rv();
        ptr += name_len;
        rv(); // extent (unused)
        uint64_t num_columns = rv();

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
        std::vector<int>  col_types(num_attr);
        std::vector<bool> col_nullable(num_attr);
        for (size_t c = 0; c < num_attr; ++c) {
            uint64_t tc = col_metas[c + 1].typeCode;
            col_nullable[c] = (tc % 2 == 1);
            uint64_t base = tc - (tc % 2);
            if      (base >= 20 && base <= 23) col_types[c] = 2; // INT
            else if (base >= 24 && base <= 27) col_types[c] = 1; // FLOAT
            else                               col_types[c] = 0; // STRING
        }

        // GEOMETRY section
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
            if (phys == 1 && dict == 3) { // VERTEX
                num_features = (size_t)(num_vals / 2);
                feat_gx.resize(num_features);
                feat_gy.resize(num_features);
                const uint8_t* vp = sd;
                for (size_t i = 0; i < num_features; ++i) {
                    uint64_t zx = 0; int sh = 0;
                    while (vp < sd + byte_len) { uint8_t bb = *vp++; zx |= (uint64_t)(bb & 0x7F) << sh; sh += 7; if (!(bb & 0x80)) break; }
                    uint64_t zy = 0; sh = 0;
                    while (vp < sd + byte_len) { uint8_t bb = *vp++; zy |= (uint64_t)(bb & 0x7F) << sh; sh += 7; if (!(bb & 0x80)) break; }
                    int32_t px = (int32_t)((zx >> 1) ^ -(int64_t)(zx & 1));
                    int32_t py = (int32_t)((zy >> 1) ^ -(int64_t)(zy & 1));
                    feat_gx[i] = off_x + scale * px;
                    feat_gy[i] = off_y - scale * py;
                }
            }
        }

        std::vector<std::vector<std::string>> attr_vals(num_attr, std::vector<std::string>(num_features));
        for (size_t c = 0; c < num_attr; ++c) {
            bool nullable = col_nullable[c];
            int ctype = col_types[c];
            bool is_str = (ctype == 0);
            std::vector<bool> present(num_features, true);
            std::vector<uint64_t> str_lens;
            const uint8_t* str_data = nullptr;

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
                if (phys == 0) {        // PRESENT
                    present = mlt_decode_bool_rle(sd, bl, num_features);
                } else if (phys == 1) { // DATA
                    if (ctype == 2) {   // INT
                        const uint8_t* dp = sd; size_t fi = 0;
                        for (uint64_t vi = 0; vi < nv; ++vi) {
                            uint64_t zig = 0; int sh = 0;
                            while (dp < sd + bl) { uint8_t bb = *dp++; zig |= (uint64_t)(bb & 0x7F) << sh; sh += 7; if (!(bb & 0x80)) break; }
                            int64_t val = (int64_t)((zig >> 1) ^ -(int64_t)(zig & 1));
                            while (fi < num_features && !present[fi]) ++fi;
                            if (fi < num_features) attr_vals[c][fi++] = std::to_string(val);
                        }
                    } else if (ctype == 1) { // FLOAT
                        const uint8_t* dp = sd; size_t fi = 0;
                        for (uint64_t vi = 0; vi < nv; ++vi) {
                            uint32_t bits = (uint32_t)dp[0] | ((uint32_t)dp[1] << 8) |
                                            ((uint32_t)dp[2] << 16) | ((uint32_t)dp[3] << 24);
                            dp += 4;
                            float fval; memcpy(&fval, &bits, 4);
                            while (fi < num_features && !present[fi]) ++fi;
                            if (fi < num_features) {
                                char tmp[32]; snprintf(tmp, sizeof(tmp), "%.9g", (double)fval);
                                attr_vals[c][fi++] = tmp;
                            }
                        }
                    } else {            // STRING DATA
                        str_data = sd; (void)nv;
                    }
                } else if (phys == 3) { // LENGTH
                    const uint8_t* dp = sd;
                    str_lens.reserve(nv);
                    for (uint64_t vi = 0; vi < nv; ++vi) {
                        uint64_t len = 0; int sh = 0;
                        while (dp < sd + bl) { uint8_t bb = *dp++; len |= (uint64_t)(bb & 0x7F) << sh; sh += 7; if (!(bb & 0x80)) break; }
                        str_lens.push_back(len);
                    }
                }
            }

            if (is_str && str_data && !str_lens.empty()) {
                const uint8_t* dp = str_data; size_t fi = 0;
                for (size_t li = 0; li < str_lens.size(); ++li) {
                    while (fi < num_features && !present[fi]) ++fi;
                    if (fi < num_features) {
                        attr_vals[c][fi++] = std::string((char*)dp, str_lens[li]);
                        dp += str_lens[li];
                    }
                }
            }
        }

        // Column names are stable within the tile: use local index = c.
        b.names.reserve(num_attr);
        for (size_t c = 0; c < num_attr; ++c)
            b.names.push_back(col_metas[c + 1].name);

        for (size_t i = 0; i < num_features; ++i) {
            double gx = feat_gx[i], gy = feat_gy[i];
            if (!point_passes(gx, gy, task)) continue;
            b.xs.push_back(gx);
            b.ys.push_back(gy);
            for (size_t c = 0; c < num_attr; ++c) {
                const std::string& v = attr_vals[c][i];
                if (!v.empty()) {
                    b.cols.push_back((int32_t)c);
                    b.vals.push_back(v);
                }
            }
            b.starts.push_back((int32_t)b.cols.size());
        }
        break; // only the first layer
    }
}

// ---- JSON line formatting (always sparse; only present attributes) --------
void append_json_rows(const TileBatch& b, const ExportCtx& ctx, std::string& out) {
    char coord[80];
    for (size_t i = 0; i < b.npoints(); ++i) {
        out += "{\"type\":\"Feature\",\"properties\": {";
        bool first = true;
        for (int32_t e = b.starts[i]; e < b.starts[i + 1]; ++e) {
            if (b.vals[e].empty()) continue;
            if (!first) out += ',';
            out += '"'; out += b.names[b.cols[e]]; out += "\":\"";
            out += b.vals[e]; out += '"';
            first = false;
        }
        int n = snprintf(coord, sizeof(coord),
                         "},\"geometry\":{\"type\":\"Point\",\"coordinates\":[%.*f,%.*f]}}\n",
                         ctx.precision, b.xs[i], ctx.precision, b.ys[i]);
        out.append(coord, n);
    }
}

// ---- Dense TSV formatting (fixed schema) ----------------------------------
void append_tsv_rows_fixed(const TileBatch& b, const ExportCtx& ctx, std::string& out) {
    const size_t ncols = ctx.ncols();
    // map this tile's local names -> global columns once
    std::vector<int32_t> l2g(b.names.size(), -1);
    for (size_t k = 0; k < b.names.size(); ++k) {
        auto it = ctx.name_to_global.find(b.names[k]);
        if (it != ctx.name_to_global.end()) l2g[k] = it->second;
    }
    std::vector<const std::string*> row(ncols, nullptr);
    char coord[80];
    for (size_t i = 0; i < b.npoints(); ++i) {
        std::fill(row.begin(), row.end(), nullptr);
        for (int32_t e = b.starts[i]; e < b.starts[i + 1]; ++e) {
            int32_t g = l2g[b.cols[e]];
            if (g >= 0) row[g] = &b.vals[e];
        }
        int n = snprintf(coord, sizeof(coord), "%.*f\t%.*f",
                         ctx.precision, b.xs[i], ctx.precision, b.ys[i]);
        out.append(coord, n);
        for (size_t c = 0; c < ncols; ++c) {
            out += '\t';
            if (row[c] && !row[c]->empty()) out += *row[c];
            else out += ctx.missing_value;
        }
        out += '\n';
    }
}

// ---- Sparse binary temp-file records (schema-unknown fallback) ------------
inline void put_u32(std::string& s, uint32_t v) { s.append((const char*)&v, sizeof(v)); }
inline void put_dbl(std::string& s, double v)    { s.append((const char*)&v, sizeof(v)); }

// Map this tile's names to a shared, growing registry (under lock), then emit
// sparse binary records into `out` using global column indices.
void append_sparse_records(const TileBatch& b, std::string& out,
                           std::vector<std::string>& reg_names,
                           std::unordered_map<std::string, int32_t>& reg_map,
                           std::mutex& reg_mtx) {
    std::vector<int32_t> l2g(b.names.size());
    {
        std::lock_guard<std::mutex> lock(reg_mtx);
        for (size_t k = 0; k < b.names.size(); ++k) {
            auto it = reg_map.find(b.names[k]);
            if (it == reg_map.end()) {
                int32_t gi = (int32_t)reg_names.size();
                reg_names.push_back(b.names[k]);
                reg_map.emplace(b.names[k], gi);
                l2g[k] = gi;
            } else {
                l2g[k] = it->second;
            }
        }
    }
    for (size_t i = 0; i < b.npoints(); ++i) {
        put_dbl(out, b.xs[i]);
        put_dbl(out, b.ys[i]);
        uint32_t np = 0;
        for (int32_t e = b.starts[i]; e < b.starts[i + 1]; ++e)
            if (!b.vals[e].empty()) ++np;
        put_u32(out, np);
        for (int32_t e = b.starts[i]; e < b.starts[i + 1]; ++e) {
            if (b.vals[e].empty()) continue;
            put_u32(out, (uint32_t)l2g[b.cols[e]]);
            put_u32(out, (uint32_t)b.vals[e].size());
            out += b.vals[e];
        }
    }
}

} // namespace

/////////////////////////////////////////////////////////////////////////
// export : Export points from a PMTiles file to TSV/JSON
////////////////////////////////////////////////////////////////////////
int32_t cmd_export_pmtiles(int32_t argc, char **argv)
{
    std::string pmtilesf;
    int32_t zoom = -1;             // -1 represents the max zoom level available
    int32_t num_threads = 0;       // 0 -> hardware concurrency

    // region-based filtering
    double xmin = -std::numeric_limits<double>::infinity();
    double xmax = std::numeric_limits<double>::infinity();
    double ymin = -std::numeric_limits<double>::infinity();
    double ymax = std::numeric_limits<double>::infinity();

    std::string geojsonf;          // polygon-based filtering

    std::string out_tsvf;
    std::string out_jsonf;

    int32_t precision = 3;
    std::string missing_value_str = "NA";
    std::vector<std::string> priority_feature_names;
    bool skip_priority_feature = false;
    bool ignore_metadata = false;  // force column discovery instead of using the metadata schema

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
    LONG_STRING_PARAM("missing-value", &missing_value_str, "String to represent missing values in the output (default: NA)")
    LONG_MULTI_STRING_PARAM("priority-feature", &priority_feature_names, "Feature name to prioritize in the output (default: gene, count)")
    LONG_PARAM("skip-priority-feature", &skip_priority_feature, "Skip priority features in the output if they are missing in the input (default: false)")
    LONG_PARAM("ignore-metadata", &ignore_metadata, "Ignore the attribute schema in the metadata and discover columns from the tiles instead (for comparison/debugging)")

    LONG_PARAM_GROUP("Performance options", NULL)
    LONG_INT_PARAM("threads", &num_threads, "Number of worker threads (default: hardware concurrency)")
    END_LONG_PARAMS();

    pl.Add(new longParams("Available Options", longParameters));
    pl.Read(argc, argv);
    pl.Status();

    notice("Analysis started");

    if (pmtilesf.empty())
        error("Missing required options --in");
    if (out_tsvf.empty() && out_jsonf.empty())
        error("Missing required options --out-tsv or --out-json (at least 1 required)");

    if (priority_feature_names.empty() && !skip_priority_feature) {
        priority_feature_names.push_back("gene");
        priority_feature_names.push_back("count");
    }

    if (num_threads <= 0) {
        num_threads = (int32_t)std::thread::hardware_concurrency();
        if (num_threads <= 0) num_threads = 4;
    }
    notice("Using %d worker threads", num_threads);

    // Open a PMTiles file and read header/metadata/entries
    pmt_pts pmt(pmtilesf.c_str());
    notice("Reading header and tile entries...");
    if (!pmt.read_header_meta_entries())
        error("This pmtiles file is malformed or incompatible with pmpoints, which requires collection of points in MVT format");

    if (zoom == -1) {
        zoom = pmt.hdr.max_zoom;
        notice("Setting the zoom level to the maximum zoom level: %d", zoom);
    }
    if (zoom < pmt.hdr.min_zoom || zoom > pmt.hdr.max_zoom)
        error("Zoom level %d is unavailable in %s", zoom, pmtilesf.c_str());

    // bounding box in tile space
    pmt_utils::pmt_pt_t min_pt(zoom, xmin, ymin);
    pmt_utils::pmt_pt_t max_pt(zoom, xmax, ymax);

    int xmin_class = std::fpclassify(xmin);
    int xmax_class = std::fpclassify(xmax);
    int ymin_class = std::fpclassify(ymin);
    int ymax_class = std::fpclassify(ymax);
    bool has_boundary = !((xmin_class == FP_INFINITE || xmin_class == FP_NAN) &&
                          (xmax_class == FP_INFINITE || xmax_class == FP_NAN) &&
                          (ymin_class == FP_INFINITE || ymin_class == FP_NAN) &&
                          (ymax_class == FP_INFINITE || ymax_class == FP_NAN));
    if (has_boundary)
        notice("Bounding Box: [(%.3f, %.3f), (%.3f, %.3f)]", xmin, ymin, xmax, ymax);
    else
        notice("No bounding box is set; all tiles at zoom level %d will be considered", zoom);

    // load polygons + bounding boxes
    std::vector<Polygon> polygons;
    if (!geojsonf.empty())
        load_polygons_from_geojson(geojsonf.c_str(), polygons);
    std::vector<Rectangle> bounding_boxes;
    for (size_t i = 0; i < polygons.size(); ++i) {
        bounding_boxes.push_back(polygons[i].get_bounding_box());
        Rectangle& r = bounding_boxes.back();
        notice("BBox %zu = ll(%lf, %lf) - ur(%lf,%lf)", i, r.p_min.x, r.p_min.y, r.p_max.x, r.p_max.y);
    }

    // ---- Build the per-tile work list (single-threaded prefilter) ----------
    std::vector<TileTask> tasks;
    uint64_t n_skipped_tiles = 0;
    for (size_t i = 0; i < pmt.tile_entries.size(); ++i) {
        pmtiles::entry_zxy& entry = pmt.tile_entries[i];
        if (entry.z != zoom) continue;

        point_t tile_min_pt(0, 0), tile_max_pt(0, 0);
        pmt_utils::tiletoepsg3857(entry.x, entry.y, entry.z, &tile_min_pt.x, &tile_max_pt.y);
        pmt_utils::tiletoepsg3857(entry.x + 1, entry.y + 1, entry.z, &tile_max_pt.x, &tile_min_pt.y);
        Rectangle tile_bbox(tile_min_pt.x, tile_min_pt.y, tile_max_pt.x, tile_max_pt.y);

        TileTask task;
        task.entry = entry;

        if (has_boundary) {
            // y-axis is inverted, so min/max swap in y when comparing tiles
            if (entry.x < min_pt.tile_x || entry.x > max_pt.tile_x ||
                entry.y < max_pt.tile_y || entry.y > min_pt.tile_y) {
                ++n_skipped_tiles;
                continue;
            }
            if (entry.x == min_pt.tile_x || entry.y == min_pt.tile_y) task.min_filt = &min_pt;
            if (entry.x == max_pt.tile_x || entry.y == max_pt.tile_y) task.max_filt = &max_pt;
        }

        if (!polygons.empty()) {
            for (size_t j = 0; j < bounding_boxes.size(); ++j)
                if (bounding_boxes[j].intersects_rectangle(tile_bbox))
                    task.polygons.push_back(&polygons[j]);
            if (task.polygons.empty()) {
                ++n_skipped_tiles;
                continue;
            }
        }
        tasks.push_back(std::move(task));
    }
    notice("Selected %zu tiles at zoom level %d (%llu skipped by filters)",
           tasks.size(), zoom, (unsigned long long)n_skipped_tiles);

    const bool is_mlt = (pmt.hdr.tile_type == 0x06);

    // ---- Determine output schema -------------------------------------------
    ExportCtx ctx;
    ctx.want_tsv = !out_tsvf.empty();
    ctx.want_json = !out_jsonf.empty();
    ctx.precision = precision;
    ctx.missing_value = missing_value_str;

    std::set<std::string> schema_names;
    bool have_schema = !ignore_metadata && extract_schema(pmt.jmeta, schema_names);
    if (ignore_metadata)
        notice("Ignoring metadata schema (--ignore-metadata); columns will be discovered from tiles");
    if (have_schema) {
        build_column_order(schema_names, priority_feature_names, ctx.col_names, ctx.name_to_global);
        ctx.schema_known = true;
        notice("Using attribute schema from metadata: %zu columns", ctx.col_names.size());
    } else if (ctx.want_tsv && !ignore_metadata) {
        notice("No attribute schema found in metadata; columns will be discovered while exporting");
    }

    // ---- Open outputs ------------------------------------------------------
    auto is_gz = [](const std::string& path) -> bool {
        return path.size() >= 3 && path.compare(path.size() - 3, 3, ".gz") == 0;
    };
    auto open_hts = [&](const std::string& path) -> htsFile* {
        return hts_open(path.c_str(), is_gz(path) ? "wz" : "w");
    };

    htsFile* tsv_wh = nullptr;
    htsFile* json_wh = nullptr;
    if (ctx.want_json) {
        json_wh = open_hts(out_jsonf);
        if (is_gz(out_jsonf)) hts_set_threads(json_wh, num_threads); // background compression
        hprintf(json_wh, "{\n");
    }
    // In the schema-known path we write the final TSV directly; in the fallback
    // path we stream to a temp file first and write the final TSV at the end.
    const bool tsv_direct = ctx.want_tsv && ctx.schema_known;
    if (tsv_direct) {
        tsv_wh = open_hts(out_tsvf);
        if (is_gz(out_tsvf)) hts_set_threads(tsv_wh, num_threads);
        hprintf(tsv_wh, "X\tY");
        for (const auto& c : ctx.col_names)
            hprintf(tsv_wh, "\t%s", c.c_str());
        hprintf(tsv_wh, "\n");
    }

    // Fallback temp file + shared column registry
    std::string tmp_path = out_tsvf + ".pmpoint_tmp";
    FILE* tmp_fp = nullptr;
    std::vector<std::string> reg_names;
    std::unordered_map<std::string, int32_t> reg_map;
    std::mutex reg_mtx;
    if (ctx.want_tsv && !ctx.schema_known) {
        tmp_fp = std::fopen(tmp_path.c_str(), "w+b");
        if (!tmp_fp) error("Failed to create temporary file %s", tmp_path.c_str());
    }

    // ---- Parallel fetch + decode + write -----------------------------------
    std::mutex write_mtx;
    std::atomic<size_t> next_task{0};
    std::atomic<uint64_t> n_written{0};
    std::atomic<size_t> tiles_done{0};

    // pre-create one reader per worker (serial open/HEAD), then run in parallel
    std::vector<std::unique_ptr<FlexReader>> readers(num_threads);
    for (int32_t t = 0; t < num_threads; ++t) {
        readers[t] = pmt.clone_reader();
        if (!readers[t]) error("Failed to create reader #%d for parallel export", t);
    }

    auto worker = [&](int32_t tid) {
        FlexReader* reader = readers[tid].get();
        std::string tile_buf;
        TileBatch batch;
        std::string tsv_block, json_block, tmp_block;
        for (;;) {
            size_t ti = next_task.fetch_add(1);
            if (ti >= tasks.size()) break;
            const TileTask& task = tasks[ti];

            pmt.fetch_tile_decompress(reader, task.entry, tile_buf);
            if (is_mlt) decode_mlt_batch(tile_buf, task, batch);
            else        decode_mvt_batch(tile_buf, task, batch);

            size_t npts = batch.npoints();
            if (npts > 0) {
                tsv_block.clear(); json_block.clear(); tmp_block.clear();
                if (tsv_direct)               append_tsv_rows_fixed(batch, ctx, tsv_block);
                if (ctx.want_json)            append_json_rows(batch, ctx, json_block);
                if (ctx.want_tsv && !ctx.schema_known)
                    append_sparse_records(batch, tmp_block, reg_names, reg_map, reg_mtx);

                std::lock_guard<std::mutex> lock(write_mtx);
                if (tsv_direct && !tsv_block.empty())
                    hprint_str(tsv_wh, tsv_block);
                if (ctx.want_json && !json_block.empty())
                    hprint_str(json_wh, json_block);
                if (tmp_fp && !tmp_block.empty())
                    std::fwrite(tmp_block.data(), 1, tmp_block.size(), tmp_fp);
            }
            n_written.fetch_add(npts);
            size_t done = tiles_done.fetch_add(1) + 1;
            if (done % 200 == 0)
                notice("Processed %zu / %zu tiles -- %llu points so far",
                       done, tasks.size(), (unsigned long long)n_written.load());
        }
    };

    std::vector<std::thread> threads;
    for (int32_t t = 0; t < num_threads; ++t)
        threads.emplace_back(worker, t);
    for (auto& th : threads) th.join();

    // ---- Finalize ----------------------------------------------------------
    if (json_wh) {
        hprintf(json_wh, "}\n");
        hts_close(json_wh);
    }
    if (tsv_direct)
        hts_close(tsv_wh);

    // Fallback: expand sparse temp records into a consistent TSV
    if (tmp_fp) {
        std::vector<std::string> ordered;
        std::unordered_map<std::string, int32_t> name_to_idx;
        std::set<std::string> reg_set(reg_names.begin(), reg_names.end());
        build_column_order(reg_set, priority_feature_names, ordered, name_to_idx);
        // reg global index -> display column position
        std::vector<int32_t> reg_to_disp(reg_names.size(), -1);
        for (size_t g = 0; g < reg_names.size(); ++g)
            reg_to_disp[g] = name_to_idx[reg_names[g]];

        tsv_wh = open_hts(out_tsvf);
        if (is_gz(out_tsvf)) hts_set_threads(tsv_wh, num_threads);
        hprintf(tsv_wh, "X\tY");
        for (const auto& c : ordered) hprintf(tsv_wh, "\t%s", c.c_str());
        hprintf(tsv_wh, "\n");

        std::fflush(tmp_fp);
        std::rewind(tmp_fp);
        const size_t ncols = ordered.size();
        std::vector<std::string> row(ncols);
        std::vector<char> present(ncols);
        std::string line, value;
        char coord[80];
        double x, y; uint32_t np, gidx, vlen;
        while (std::fread(&x, sizeof(x), 1, tmp_fp) == 1) {
            if (std::fread(&y, sizeof(y), 1, tmp_fp) != 1) break;
            if (std::fread(&np, sizeof(np), 1, tmp_fp) != 1) break;
            std::fill(present.begin(), present.end(), 0);
            for (uint32_t k = 0; k < np; ++k) {
                if (std::fread(&gidx, sizeof(gidx), 1, tmp_fp) != 1) { np = 0; break; }
                if (std::fread(&vlen, sizeof(vlen), 1, tmp_fp) != 1) { np = 0; break; }
                value.resize(vlen);
                if (vlen && std::fread(&value[0], 1, vlen, tmp_fp) != vlen) { np = 0; break; }
                int32_t disp = (gidx < reg_to_disp.size()) ? reg_to_disp[gidx] : -1;
                if (disp >= 0) { row[disp] = value; present[disp] = 1; }
            }
            int n = snprintf(coord, sizeof(coord), "%.*f\t%.*f", precision, x, precision, y);
            line.assign(coord, n);
            for (size_t c = 0; c < ncols; ++c) {
                line += '\t';
                line += present[c] ? row[c] : missing_value_str;
            }
            line += '\n';
            hprint_str(tsv_wh, line);
        }
        hts_close(tsv_wh);
        std::fclose(tmp_fp);
        std::remove(tmp_path.c_str());
        notice("Finalized TSV with %zu columns from discovered schema", ncols);
    }

    notice("Finished writing %llu points in total", (unsigned long long)n_written.load());
    notice("Analysis Finished");
    return 0;
}
