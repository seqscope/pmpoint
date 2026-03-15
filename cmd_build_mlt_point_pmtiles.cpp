#include "ext/nlohmann/json.hpp"
#include "pmpoint.h"
#include "qgenlib/params.h"
#include "qgenlib/qgen_error.h"
#include "qgenlib/tsv_reader.h"
#include "pmt_utils.h"
#include "ext/PMTiles/pmtiles.hpp"

#include <vector>
#include <string>
#include <map>
#include <unordered_map>
#include <list>
#include <iostream>
#include <fcntl.h>
#include <unistd.h>
#include <sys/stat.h>
#include <fstream>
#include <sstream>
#include <cmath>
#include <limits>
#include <cstring>
#include <zlib.h>
#include <thread>
#include <mutex>
#include <atomic>
#include <algorithm>

std::string mlt_gzip_compress(const std::string& data) {
    z_stream zs;
    memset(&zs, 0, sizeof(zs));
    if (deflateInit2(&zs, Z_DEFAULT_COMPRESSION, Z_DEFLATED, 15 | 16, 8, Z_DEFAULT_STRATEGY) != Z_OK) {
        error("deflateInit2 failed");
    }
    zs.next_in = (Bytef*)data.data();
    zs.avail_in = data.size();

    int ret;
    char outbuffer[32768];
    std::string outstring;

    do {
        zs.next_out = reinterpret_cast<Bytef*>(outbuffer);
        zs.avail_out = sizeof(outbuffer);
        ret = deflate(&zs, Z_FINISH);
        if (outstring.size() < zs.total_out) {
            outstring.append(outbuffer, zs.total_out - outstring.size());
        }
    } while (ret == Z_OK);

    deflateEnd(&zs);
    return outstring;
}

void epsg3857_to_wgs84(double x, double y, double& lon, double& lat) {
    const double R = 6378137.0;
    lon = (x / R) * (180.0 / M_PI);
    lat = (2.0 * std::atan(std::exp(y / R)) - M_PI / 2.0) * (180.0 / M_PI);
}

// Column type constants (ordered by specificity: INT > FLOAT > STRING)
static const int COL_TYPE_STRING = 0;
static const int COL_TYPE_FLOAT  = 1;
static const int COL_TYPE_INT    = 2;

// Returns true if the value is considered missing (absent)
static bool is_missing(const std::string& v) {
    return v.empty() || v == "NA";
}

struct PointFeature {
    int32_t x, y;
    std::vector<std::string> attrs; // one entry per non-coordinate column
};

// Encode a boolean array as byte RLE (ORC-style).
// The PRESENT stream stores 1 bit per feature (1=present, 0=absent), packed into
// bytes (LSB first), then encoded as byte RLE literal runs of up to 128 bytes each.
static std::string encode_bool_rle(const std::vector<bool>& present) {
    size_t n = present.size();
    size_t numBytes = (n + 7) / 8;

    // Pack bits into bytes
    std::vector<uint8_t> packed(numBytes, 0);
    for (size_t i = 0; i < n; ++i) {
        if (present[i]) packed[i / 8] |= static_cast<uint8_t>(1 << (i % 8));
    }

    // Emit as byte-RLE literal runs (up to 128 packed bytes per run)
    std::string result;
    size_t offset = 0;
    while (offset < numBytes) {
        size_t chunk = std::min(static_cast<size_t>(128), numBytes - offset);
        // Literal header: 256 - chunk  (unsigned; decoder reads 256-header literals)
        result.push_back(static_cast<char>(static_cast<uint8_t>(256 - chunk)));
        for (size_t i = 0; i < chunk; ++i) {
            result.push_back(static_cast<char>(packed[offset + i]));
        }
        offset += chunk;
    }
    return result;
}

class MLTPointEncoder {
    std::string out;

    void write_varint(std::uint64_t value) {
        do {
            uint8_t byte = value & 0x7F;
            value >>= 7;
            if (value > 0) byte |= 0x80;
            out.push_back(byte);
        } while (value > 0);
    }

    uint32_t encode_zigzag32(int32_t v) {
        return (v << 1) ^ (v >> 31);
    }

public:
    // col_types:   0=STRING, 1=FLOAT, 2=INT_64
    // col_nullable: true  → nullable typeCode (21/25/29) + PRESENT stream; DATA holds non-null values only
    //               false → non-nullable typeCode (20/24/28); DATA holds all values
    std::string encode(uint32_t extent, const std::string& layer_name,
                       const std::vector<PointFeature>& features,
                       const std::vector<std::string>& col_names,
                       const std::vector<int>& col_types,
                       const std::vector<bool>& col_nullable) {
        out.clear();

        std::string tmp_ft;

        auto append_varint = [&](std::uint64_t value) {
            do {
                uint8_t byte = value & 0x7F;
                value >>= 7;
                if (value > 0) byte |= 0x80;
                tmp_ft.push_back(byte);
            } while (value > 0);
        };

        // === METADATA SECTION ===
        // Layer name
        append_varint(layer_name.size());
        tmp_ft.append(layer_name);
        // Extent
        append_varint(extent);
        // Number of columns: 1 geometry + N attribute columns
        append_varint(1 + col_names.size());

        // Column 0: GEOMETRY (typeCode=4, no name per spec)
        append_varint(4);

        // Attribute columns metadata: typeCode + name
        // Non-nullable: INT_64=20, FLOAT=24, STRING=28
        // Nullable:     INT_64=21, FLOAT=25, STRING=29
        for (size_t c = 0; c < col_names.size(); ++c) {
            int typeCode;
            if (col_types[c] == COL_TYPE_INT)
                typeCode = col_nullable[c] ? 21 : 20;
            else if (col_types[c] == COL_TYPE_FLOAT)
                typeCode = col_nullable[c] ? 25 : 24;
            else
                typeCode = col_nullable[c] ? 29 : 28;
            append_varint(typeCode);
            append_varint(col_names[c].size());
            tmp_ft.append(col_names[c]);
        }

        // === DATA SECTION ===

        // --- GEOMETRY column data ---
        // numStreams = 2 (GeometryType stream + Vertex stream)
        append_varint(2);

        // Stream 0: GeometryType
        // Header byte 1: (DATA=1 << 4) | (NONE=0) = 0x10
        // Header byte 2: VARINT technique = 0x02
        tmp_ft.push_back(0x10);
        tmp_ft.push_back(0x02);
        append_varint(features.size()); // numValues = numFeatures

        std::string geom_type_data;
        for (size_t i = 0; i < features.size(); ++i) {
            geom_type_data.push_back(0); // 0 = POINT
        }
        append_varint(geom_type_data.size());
        tmp_ft.append(geom_type_data);

        // Stream 1: Vertices
        // Header byte 1: (DATA=1 << 4) | (VERTEX=3) = 0x13
        // Header byte 2: VARINT technique = 0x02
        tmp_ft.push_back(0x13);
        tmp_ft.push_back(0x02);
        append_varint(features.size() * 2); // numValues = X and Y per point

        std::string vertex_data;
        auto append_v_varint = [&](std::uint64_t value) {
            do {
                uint8_t byte = value & 0x7F;
                value >>= 7;
                if (value > 0) byte |= 0x80;
                vertex_data.push_back(byte);
            } while (value > 0);
        };
        for (size_t i = 0; i < features.size(); ++i) {
            append_v_varint(encode_zigzag32(features[i].x));
            append_v_varint(encode_zigzag32(features[i].y));
        }
        append_varint(vertex_data.size());
        tmp_ft.append(vertex_data);

        // --- Attribute column data ---
        for (size_t c = 0; c < col_names.size(); ++c) {
            bool nullable = col_nullable[c];

            // Build the presence bitmap and collect non-null raw values
            std::vector<bool> present(features.size(), true);
            if (nullable) {
                for (size_t i = 0; i < features.size(); ++i) {
                    const std::string& s = (c < features[i].attrs.size()) ? features[i].attrs[c] : "";
                    present[i] = !is_missing(s);
                }
            }
            size_t non_null_count = 0;
            for (size_t i = 0; i < features.size(); ++i) {
                if (present[i]) ++non_null_count;
            }

            if (col_types[c] == COL_TYPE_INT) {
                // INT_64: non-nullable typeCode=20, nullable typeCode=21
                // hasStreamCount=false: NO numStreams written.
                // Streams: [PRESENT (if nullable)] + DATA(non-null values only)

                // Encode non-null zigzag+varint int64 values
                std::string int_data;
                auto append_int_varint = [&](uint64_t value) {
                    do {
                        uint8_t byte = value & 0x7F;
                        value >>= 7;
                        if (value > 0) byte |= 0x80;
                        int_data.push_back(byte);
                    } while (value > 0);
                };
                for (size_t i = 0; i < features.size(); ++i) {
                    if (!present[i]) continue;
                    const std::string& s = (c < features[i].attrs.size()) ? features[i].attrs[c] : "";
                    int64_t val = 0;
                    try { val = std::stoll(s); } catch (...) {}
                    uint64_t zigzag = ((uint64_t)val << 1) ^ (uint64_t)(val >> 63);
                    append_int_varint(zigzag);
                }

                if (nullable) {
                    // PRESENT stream:
                    //   byte 1: (PRESENT=0 << 4) | 0 = 0x00
                    //   byte 2: VARINT technique = 0x02
                    //   numValues = total features; only non-null entries in DATA
                    std::string rle = encode_bool_rle(present);
                    tmp_ft.push_back(0x00);
                    tmp_ft.push_back(0x02);
                    append_varint(features.size()); // numValues = total features
                    append_varint(rle.size());
                    tmp_ft.append(rle);
                }

                // DATA stream
                tmp_ft.push_back(0x10); // DATA, NONE dict
                tmp_ft.push_back(0x02); // VARINT technique
                append_varint(non_null_count); // numValues = non-null count
                append_varint(int_data.size());
                tmp_ft.append(int_data);

            } else if (col_types[c] == COL_TYPE_FLOAT) {
                // FLOAT: non-nullable typeCode=24, nullable typeCode=25
                // hasStreamCount=false: NO numStreams written.
                // Streams: [PRESENT (if nullable)] + DATA(non-null LE float32 values)

                std::string float_data;
                float_data.reserve(non_null_count * 4);
                for (size_t i = 0; i < features.size(); ++i) {
                    if (!present[i]) continue;
                    const std::string& s = (c < features[i].attrs.size()) ? features[i].attrs[c] : "";
                    float val = 0.0f;
                    try { val = std::stof(s); } catch (...) {}
                    uint32_t bits;
                    memcpy(&bits, &val, sizeof(bits));
                    float_data.push_back(static_cast<char>(bits & 0xFF));
                    float_data.push_back(static_cast<char>((bits >> 8) & 0xFF));
                    float_data.push_back(static_cast<char>((bits >> 16) & 0xFF));
                    float_data.push_back(static_cast<char>((bits >> 24) & 0xFF));
                }

                if (nullable) {
                    std::string rle = encode_bool_rle(present);
                    tmp_ft.push_back(0x00); // PRESENT
                    tmp_ft.push_back(0x02); // VARINT
                    append_varint(features.size());
                    append_varint(rle.size());
                    tmp_ft.append(rle);
                }

                // DATA stream:
                //   byte 1: (DATA=1 << 4) | (NONE=0) = 0x10
                //   byte 2: NONE technique = 0x00 (raw LE float32)
                tmp_ft.push_back(0x10);
                tmp_ft.push_back(0x00);
                append_varint(non_null_count); // numValues = non-null count
                append_varint(float_data.size());
                tmp_ft.append(float_data);

            } else {
                // STRING: non-nullable typeCode=28, nullable typeCode=29
                // hasStreamCount=true: numStreams is written.
                // Streams: [PRESENT (if nullable)] + LENGTH + DATA
                // LENGTH and DATA cover only non-null strings.
                int numStreams = nullable ? 3 : 2;
                append_varint(numStreams);

                // Build lengths and bytes for non-null strings only
                std::string lengths_data;
                std::string strings_data;
                auto append_len_varint = [&](uint64_t v) {
                    do {
                        uint8_t b = v & 0x7F;
                        v >>= 7;
                        if (v > 0) b |= 0x80;
                        lengths_data.push_back(b);
                    } while (v > 0);
                };
                for (size_t i = 0; i < features.size(); ++i) {
                    if (!present[i]) continue;
                    const std::string& s = (c < features[i].attrs.size()) ? features[i].attrs[c] : "";
                    append_len_varint(s.size());
                    strings_data.append(s);
                }

                if (nullable) {
                    std::string rle = encode_bool_rle(present);
                    tmp_ft.push_back(0x00); // PRESENT
                    tmp_ft.push_back(0x02); // VARINT
                    append_varint(features.size());
                    append_varint(rle.size());
                    tmp_ft.append(rle);
                }

                // LENGTH stream:
                //   byte 1: (LENGTH=3 << 4) | (VAR_BINARY=0) = 0x30
                //   byte 2: VARINT technique = 0x02
                tmp_ft.push_back(0x30);
                tmp_ft.push_back(0x02);
                append_varint(non_null_count); // numValues = non-null string count
                append_varint(lengths_data.size());
                tmp_ft.append(lengths_data);

                // DATA stream:
                //   byte 1: (DATA=1 << 4) | (NONE=0) = 0x10
                //   byte 2: NONE technique = 0x00 (raw UTF-8 bytes)
                //   numValues = 0 (raw byte stream; byteLength is sufficient)
                tmp_ft.push_back(0x10);
                tmp_ft.push_back(0x00);
                append_varint(0); // numValues = 0
                append_varint(strings_data.size());
                tmp_ft.append(strings_data);
            }
        }

        std::string layer_data = tmp_ft;

        // Layer wrapper: varint(totalLength) | varint(tag=1) | layerData
        // totalLength includes the 1-byte tag varint
        write_varint(1 + layer_data.size());
        out.push_back(1); // Layer Tag = 1
        out.append(layer_data);

        return out;
    }
};

int32_t cmd_build_mlt_point_pmtiles(int32_t argc, char **argv)
{
    std::string in_csv;
    std::string out_pmtiles;
    std::string tmp_dir = ".";
    int32_t zoom = 0;
    std::string colname_x = "X";
    std::string colname_y = "Y";
    std::string delim = ",";
    int32_t n_threads = 1;

    paramList pl;
    BEGIN_LONG_PARAMS(longParameters)
    LONG_STRING_PARAM("in", &in_csv, "Input CSV file")
    LONG_STRING_PARAM("out", &out_pmtiles, "Output pyramidal PMTiles file")
    LONG_INT_PARAM("zoom", &zoom, "Zoom level to build")
    LONG_STRING_PARAM("colname-x", &colname_x, "Column name for X coordinate (EPSG:3857)")
    LONG_STRING_PARAM("colname-y", &colname_y, "Column name for Y coordinate (EPSG:3857)")
    LONG_STRING_PARAM("delim", &delim, "Delimiter for input file")
    LONG_STRING_PARAM("tmp-dir", &tmp_dir, "Temporary directory")
    LONG_INT_PARAM("threads", &n_threads, "Number of threads for encoding [1]")
    END_LONG_PARAMS();

    pl.Add(new longParams("Available Options", longParameters));
    pl.Read(argc, argv);
    pl.Status();

    if (in_csv.empty() || out_pmtiles.empty()) error("Missing --in or --out");
    if (delim.size() != 1) error("Delimiter must be a single character");
    if (n_threads < 1) n_threads = 1;

    std::string base_name = in_csv;
    {
        size_t p = base_name.find_last_of("/\\");
        if (p != std::string::npos) base_name = base_name.substr(p + 1);
        size_t d = base_name.find_last_of(".");
        if (d != std::string::npos) base_name = base_name.substr(0, d);
    }

    tsv_reader tr(in_csv.c_str());
    tr.delimiter = (int)delim[0];

    std::vector<std::string> headers;
    if (tr.read_line()) {
        for (int i = 0; i < tr.nfields; ++i)
            headers.push_back(tr.str_field_at(i));
    } else {
        error("Empty file or failed to read header");
    }

    int col_x = -1, col_y = -1;
    for (size_t i = 0; i < headers.size(); ++i) {
        if (headers[i] == colname_x) col_x = (int)i;
        if (headers[i] == colname_y) col_y = (int)i;
    }
    if (col_x == -1 || col_y == -1) error("Could not find X or Y columns in header.");

    std::vector<int>         attr_col_indices;
    std::vector<std::string> attr_col_names;
    for (size_t i = 0; i < headers.size(); ++i) {
        if ((int)i != col_x && (int)i != col_y) {
            attr_col_indices.push_back((int)i);
            attr_col_names.push_back(headers[i]);
        }
    }
    size_t n_attrs = attr_col_names.size();

    // =========================================================
    // Phase 1: Read CSV → per-tile binary temp files
    //   Memory: O(max_open_tile_buffers) instead of O(all_points)
    // =========================================================
    notice("Phase 1: Reading %s and writing per-tile temp files...", in_csv.c_str());

    // Column type/nullability detection (updated inline, O(n_cols) memory)
    std::vector<int>  attr_col_types(n_attrs, COL_TYPE_INT);
    std::vector<bool> attr_col_nullable(n_attrs, false);

    // Per-tile state
    struct TileInfo {
        std::string path;
        FILE*       fp           = nullptr;
        uint64_t    point_count  = 0;
    };
    std::unordered_map<uint64_t, TileInfo> tile_infos;

    // LRU file-handle cache to avoid hitting OS open-file limits
    const size_t MAX_OPEN_FDS = 512;
    std::list<uint64_t> lru_list;   // front = most recently used
    std::unordered_map<uint64_t, std::list<uint64_t>::iterator> lru_pos;
    size_t open_fd_count = 0;
    std::string pid_str = std::to_string(getpid());

    // Returns an open (append-mode) FILE* for the given tile, managing the LRU.
    auto get_tile_fp = [&](uint64_t tile_id) -> FILE* {
        TileInfo& ti = tile_infos[tile_id];
        if (ti.path.empty())
            ti.path = tmp_dir + "/pmpoint_mlt_" + pid_str + "_" + std::to_string(tile_id) + ".bin";

        auto pos_it = lru_pos.find(tile_id);
        if (pos_it != lru_pos.end()) {
            // Already open: move to front of LRU
            lru_list.erase(pos_it->second);
            lru_list.push_front(tile_id);
            lru_pos[tile_id] = lru_list.begin();
        } else {
            // Not open: evict LRU tail if at limit, then open
            if (open_fd_count >= MAX_OPEN_FDS) {
                uint64_t evict_id = lru_list.back();
                lru_list.pop_back();
                lru_pos.erase(evict_id);
                fclose(tile_infos[evict_id].fp);
                tile_infos[evict_id].fp = nullptr;
                --open_fd_count;
            }
            lru_list.push_front(tile_id);
            lru_pos[tile_id] = lru_list.begin();
            ti.fp = fopen(ti.path.c_str(), "ab");
            if (!ti.fp) error("Cannot create tile temp file: %s", ti.path.c_str());
            ++open_fd_count;
        }
        return ti.fp;
    };

    double min_x = std::numeric_limits<double>::max();
    double min_y = std::numeric_limits<double>::max();
    double max_x = std::numeric_limits<double>::lowest();
    double max_y = std::numeric_limits<double>::lowest();
    uint64_t point_count = 0;
    uint64_t nlines      = 0;
    const uint32_t extent = 4096;

    while (tr.read_line()) {
        if (tr.nfields <= std::max(col_x, col_y)) continue;

        double cx, cy;
        try {
            cx = std::stod(tr.str_field_at(col_x));
            cy = std::stod(tr.str_field_at(col_y));
        } catch (...) {
            continue;
        }

        int64_t tx, ty;
        double lx, ly;
        pmt_utils::epsg3857totilecoord(cx, cy, zoom, &tx, &ty, &lx, &ly);
        int32_t px = (int32_t)std::round(lx * extent / 256.0);
        int32_t py = (int32_t)std::round(ly * extent / 256.0);
        uint64_t tile_id = pmtiles::zxy_to_tileid(zoom, (uint32_t)tx, (uint32_t)ty);

        // Collect attr values + update type detection inline
        for (size_t ai = 0; ai < n_attrs; ++ai) {
            int ci = attr_col_indices[ai];
            const char* raw = (ci < tr.nfields) ? tr.str_field_at(ci) : "";
            std::string v(raw);

            if (is_missing(v)) {
                attr_col_nullable[ai] = true;
            } else if (attr_col_types[ai] != COL_TYPE_STRING) {
                if (attr_col_types[ai] == COL_TYPE_INT) {
                    try {
                        size_t pos;
                        std::stoll(v, &pos);
                        if (pos != v.size()) throw std::invalid_argument("trailing");
                    } catch (...) {
                        try { std::stod(v); attr_col_types[ai] = COL_TYPE_FLOAT; }
                        catch (...) { attr_col_types[ai] = COL_TYPE_STRING; }
                    }
                } else {
                    try { std::stod(v); }
                    catch (...) { attr_col_types[ai] = COL_TYPE_STRING; }
                }
            }
        }

        // Write binary point record to tile temp file:
        //   int32_t px | int32_t py | for each attr: uint32_t len | char[len]
        FILE* fp = get_tile_fp(tile_id);
        fwrite(&px, sizeof(px), 1, fp);
        fwrite(&py, sizeof(py), 1, fp);
        for (size_t ai = 0; ai < n_attrs; ++ai) {
            int ci = attr_col_indices[ai];
            const char* raw = (ci < tr.nfields) ? tr.str_field_at(ci) : "";
            uint32_t len = (uint32_t)strlen(raw);
            fwrite(&len, sizeof(len), 1, fp);
            if (len > 0) fwrite(raw, 1, len, fp);
        }
        tile_infos[tile_id].point_count++;

        min_x = std::min(min_x, cx);
        min_y = std::min(min_y, cy);
        max_x = std::max(max_x, cx);
        max_y = std::max(max_y, cy);
        ++point_count;
        ++nlines;
        if (nlines % 1000000 == 0)
            notice("Read %llu lines, %llu valid points, %zu tiles so far...",
                   nlines, point_count, tile_infos.size());
    }
    tr.close();

    // Close all open tile files
    for (auto& kv : tile_infos) {
        if (kv.second.fp) { fclose(kv.second.fp); kv.second.fp = nullptr; }
    }

    notice("Read %llu valid points into %zu tiles.", point_count, tile_infos.size());

    for (size_t c = 0; c < n_attrs; ++c) {
        const char* ts = attr_col_types[c] == COL_TYPE_INT   ? "int"
                       : attr_col_types[c] == COL_TYPE_FLOAT ? "float" : "string";
        notice("Column '%s': %s%s", attr_col_names[c].c_str(), ts,
               attr_col_nullable[c] ? " (nullable)" : "");
    }

    // Find center tile (tile with the most points)
    uint64_t center_tile_id = 0;
    {
        uint64_t max_pts = 0;
        for (auto& kv : tile_infos) {
            if (kv.second.point_count > max_pts) {
                max_pts = kv.second.point_count;
                center_tile_id = kv.first;
            }
        }
    }

    // Build sorted tile list for PMTiles directory ordering
    std::vector<uint64_t> sorted_tile_ids;
    sorted_tile_ids.reserve(tile_infos.size());
    for (auto& kv : tile_infos) sorted_tile_ids.push_back(kv.first);
    std::sort(sorted_tile_ids.begin(), sorted_tile_ids.end());

    // =========================================================
    // Phase 2: Encode tiles from temp files → compressed output
    //   Memory: O(n_threads × largest_tile) at a time
    // =========================================================
    notice("Phase 2: Encoding %zu tiles with %d thread(s)...",
           sorted_tile_ids.size(), n_threads);

    // Helper: read one tile's points from its temp file (deletes the file after)
    auto read_tile_points = [&](uint64_t tile_id) -> std::vector<PointFeature> {
        TileInfo& ti = tile_infos[tile_id];
        FILE* f = fopen(ti.path.c_str(), "rb");
        if (!f) error("Cannot open tile temp file for reading: %s", ti.path.c_str());
        std::vector<PointFeature> features;
        features.reserve(ti.point_count);
        int32_t px, py;
        while (fread(&px, sizeof(px), 1, f) == 1) {
            if (fread(&py, sizeof(py), 1, f) != 1) break;
            PointFeature pf;
            pf.x = px; pf.y = py;
            pf.attrs.resize(n_attrs);
            for (size_t ai = 0; ai < n_attrs; ++ai) {
                uint32_t len = 0;
                if (fread(&len, sizeof(len), 1, f) != 1) break;
                if (len > 0) {
                    pf.attrs[ai].resize(len);
                    if (fread(&pf.attrs[ai][0], 1, len, f) != len) break;
                }
            }
            features.push_back(std::move(pf));
        }
        fclose(f);
        unlink(ti.path.c_str());
        return features;
    };

    std::string tmp_file = tmp_dir + "/pmpoint_mlt_" + pid_str + ".tmp";
    int out_fd = open(tmp_file.c_str(), O_RDWR | O_CREAT | O_TRUNC, 0666);
    if (out_fd < 0) error("Failed to open temporary output file: %s", tmp_file.c_str());

    std::vector<pmtiles::entryv3> final_entries;
    final_entries.reserve(sorted_tile_ids.size());
    uint64_t current_out_offset = 0;

    // Process in batches of n_threads tiles: load → encode in parallel → write
    size_t n_tiles     = sorted_tile_ids.size();
    size_t batch_size  = (size_t)n_threads;
    size_t tiles_done  = 0;

    for (size_t batch_start = 0; batch_start < n_tiles; batch_start += batch_size) {
        size_t batch_end = std::min(batch_start + batch_size, n_tiles);
        size_t bs        = batch_end - batch_start;

        // Load this batch's tile data from temp files
        std::vector<std::vector<PointFeature>> batch_data(bs);
        for (size_t i = 0; i < bs; ++i)
            batch_data[i] = read_tile_points(sorted_tile_ids[batch_start + i]);

        // Encode + compress in parallel
        std::vector<std::string> batch_compressed(bs);
        if (n_threads > 1) {
            std::vector<std::thread> threads;
            for (size_t i = 0; i < bs; ++i) {
                threads.emplace_back([&, i]() {
                    MLTPointEncoder enc;
                    batch_compressed[i] = mlt_gzip_compress(
                        enc.encode(extent, base_name, batch_data[i],
                                   attr_col_names, attr_col_types, attr_col_nullable));
                    std::vector<PointFeature>().swap(batch_data[i]);
                });
            }
            for (auto& t : threads) t.join();
        } else {
            MLTPointEncoder enc;
            for (size_t i = 0; i < bs; ++i) {
                batch_compressed[i] = mlt_gzip_compress(
                    enc.encode(extent, base_name, batch_data[i],
                               attr_col_names, attr_col_types, attr_col_nullable));
                std::vector<PointFeature>().swap(batch_data[i]);
            }
        }

        // Write batch to the output temp file
        for (size_t i = 0; i < bs; ++i) {
            uint64_t tid = sorted_tile_ids[batch_start + i];
            pwrite(out_fd, batch_compressed[i].data(), batch_compressed[i].size(),
                   current_out_offset);
            final_entries.emplace_back(tid, current_out_offset,
                                       (uint32_t)batch_compressed[i].size(), 1);
            current_out_offset += batch_compressed[i].size();
        }

        tiles_done += bs;
        if (tiles_done % 10000 == 0 || tiles_done == n_tiles)
            notice("Encoded %zu / %zu tiles...", tiles_done, n_tiles);
    }

    notice("Sorting directory...");
    std::sort(final_entries.begin(), final_entries.end(), [](const pmtiles::entryv3& a, const pmtiles::entryv3& b){
        return a.tile_id < b.tile_id;
    });

    notice("Generating PMTiles metadata...");
    std::function<std::string(const std::string&, uint8_t)> compress_fn = [](const std::string& data, uint8_t comp) {
        return mlt_gzip_compress(data);
    };

    auto root_leaves = pmtiles::make_root_leaves(compress_fn, pmtiles::COMPRESSION_GZIP, final_entries);
    std::string root_bytes = std::get<0>(root_leaves);
    std::string leaves_bytes = std::get<1>(root_leaves);

    nlohmann::json jmeta;
    jmeta["name"] = base_name;
    jmeta["type"] = "overlay";
    jmeta["version"] = "2";
    jmeta["format"] = "pbf"; // To mimic tippecanoe
    jmeta["description"] = "Generated PMTiles by pmpoint for MLT points";
    jmeta["generator"] = "pmpoint build-mlt-point-pmtiles";

    nlohmann::json vlayers = nlohmann::json::array();
    nlohmann::json layer1;
    layer1["id"] = base_name;
    nlohmann::json fields = nlohmann::json::object();
    for (size_t c = 0; c < attr_col_names.size(); ++c) {
        fields[attr_col_names[c]] = (attr_col_types[c] != COL_TYPE_STRING) ? "Number" : "String";
    }
    layer1["fields"] = fields;
    layer1["minzoom"] = zoom;
    layer1["maxzoom"] = zoom;
    vlayers.push_back(layer1);

    jmeta["vector_layers"] = vlayers;

    nlohmann::json tstats;
    tstats["layerCount"] = 1;
    nlohmann::json tlayers = nlohmann::json::array();
    nlohmann::json tlayer1;
    tlayer1["layer"] = base_name;
    tlayer1["count"] = point_count;
    tlayer1["geometry"] = "Point";
    tlayer1["attributeCount"] = (int)attr_col_names.size();
    nlohmann::json attributes = nlohmann::json::array();
    for (size_t c = 0; c < attr_col_names.size(); ++c) {
        nlohmann::json attr;
        attr["attribute"] = attr_col_names[c];
        attr["type"] = (attr_col_types[c] != COL_TYPE_STRING) ? "number" : "string";
        attributes.push_back(attr);
    }
    tlayer1["attributes"] = attributes;
    tlayers.push_back(tlayer1);
    tstats["layers"] = tlayers;

    jmeta["tilestats"] = tstats;

    std::string json_metadata = jmeta.dump();
    std::string compressed_json = mlt_gzip_compress(json_metadata);

    pmtiles::headerv3 header;
    header.root_dir_offset = 127;
    header.root_dir_bytes = root_bytes.size();
    header.json_metadata_offset = header.root_dir_offset + header.root_dir_bytes;
    header.json_metadata_bytes = compressed_json.size();
    header.leaf_dirs_offset = header.json_metadata_offset + header.json_metadata_bytes;
    header.leaf_dirs_bytes = leaves_bytes.size();
    header.tile_data_offset = header.leaf_dirs_offset + header.leaf_dirs_bytes;
    header.tile_data_bytes = current_out_offset;
    header.addressed_tiles_count = final_entries.size();
    header.tile_entries_count = final_entries.size();
    header.tile_contents_count = final_entries.size();
    header.clustered = true;
    header.internal_compression = pmtiles::COMPRESSION_GZIP;
    header.tile_compression = pmtiles::COMPRESSION_GZIP;
    header.tile_type = 0x06; // 0x06 = MapLibre Vector Tile (MLT)
    header.min_zoom = zoom;
    header.max_zoom = zoom;
    header.center_zoom = zoom;

    if (point_count > 0) {
        double min_lon, min_lat, max_lon, max_lat;
        epsg3857_to_wgs84(min_x, min_y, min_lon, min_lat);
        epsg3857_to_wgs84(max_x, max_y, max_lon, max_lat);

        pmtiles::zxy center_zxy = pmtiles::tileid_to_zxy(center_tile_id);
        double center_x, center_y;
        pmt_utils::tilecoordtoespg3857(center_zxy.x, center_zxy.y, 0.5, 0.5, center_zxy.z, &center_x, &center_y);
        double center_lon, center_lat;
        epsg3857_to_wgs84(center_x, center_y, center_lon, center_lat);

        header.min_lon_e7 = static_cast<int32_t>(min_lon * 10000000.0);
        header.min_lat_e7 = static_cast<int32_t>(min_lat * 10000000.0);
        header.max_lon_e7 = static_cast<int32_t>(max_lon * 10000000.0);
        header.max_lat_e7 = static_cast<int32_t>(max_lat * 10000000.0);
        header.center_lon_e7 = static_cast<int32_t>(center_lon * 10000000.0);
        header.center_lat_e7 = static_cast<int32_t>(center_lat * 10000000.0);
    } else {
        header.min_lon_e7 = 0;
        header.min_lat_e7 = 0;
        header.max_lon_e7 = 0;
        header.max_lat_e7 = 0;
        header.center_lon_e7 = 0;
        header.center_lat_e7 = 0;
    }

    notice("Writing to %s...", out_pmtiles.c_str());
    int out_final = open(out_pmtiles.c_str(), O_WRONLY | O_CREAT | O_TRUNC, 0666);
    if (out_final < 0) error("Failed to create output file");

    std::string h_str = header.serialize();
    write(out_final, h_str.data(), h_str.size());
    write(out_final, root_bytes.data(), root_bytes.size());
    write(out_final, compressed_json.data(), compressed_json.size());
    if (leaves_bytes.size() > 0) {
        write(out_final, leaves_bytes.data(), leaves_bytes.size());
    }

    // Copy tmp file
    lseek(out_fd, 0, SEEK_SET);
    char cat_buffer[4096 * 16];
    ssize_t rb;
    while ((rb = read(out_fd, cat_buffer, sizeof(cat_buffer))) > 0) {
        write(out_final, cat_buffer, rb);
    }

    close(out_fd);
    close(out_final);

    unlink(tmp_file.c_str());
    notice("Done!");
    return 0;
}
