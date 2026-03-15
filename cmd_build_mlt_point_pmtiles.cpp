#include "ext/nlohmann/json.hpp"
#include "pmpoint.h"
#include "qgenlib/params.h"
#include "qgenlib/qgen_error.h"
#include "pmt_utils.h"
#include "ext/PMTiles/pmtiles.hpp"

#include <vector>
#include <string>
#include <map>
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

    char delimiter = delim.empty() ? ',' : delim[0];

    std::string base_name = in_csv;
    size_t slash_pos = base_name.find_last_of("/\\");
    if (slash_pos != std::string::npos) base_name = base_name.substr(slash_pos + 1);
    size_t dot_pos = base_name.find_last_of(".");
    if (dot_pos != std::string::npos) base_name = base_name.substr(0, dot_pos);

    std::ifstream file(in_csv);
    if (!file.is_open()) error("Could not open input CSV %s", in_csv.c_str());

    std::string line;
    if (!std::getline(file, line)) error("Empty file");

    std::vector<std::string> headers;
    std::stringstream ss(line);
    std::string token;
    while(std::getline(ss, token, delimiter)) {
        headers.push_back(token);
    }

    int col_x = -1;
    int col_y = -1;
    for(size_t i=0; i<headers.size(); ++i) {
        if (headers[i] == colname_x) col_x = i;
        if (headers[i] == colname_y) col_y = i;
    }

    if (col_x == -1 || col_y == -1) error("Could not find X or Y columns in header.");

    // Identify attribute columns (all columns except X and Y)
    std::vector<int> attr_col_indices;
    std::vector<std::string> attr_col_names;
    for (size_t i = 0; i < headers.size(); ++i) {
        if ((int)i != col_x && (int)i != col_y) {
            attr_col_indices.push_back(i);
            attr_col_names.push_back(headers[i]);
        }
    }

    std::map<uint64_t, std::vector<PointFeature>> tile_points;
    uint64_t point_count = 0;

    double min_x = std::numeric_limits<double>::max();
    double min_y = std::numeric_limits<double>::max();
    double max_x = std::numeric_limits<double>::lowest();
    double max_y = std::numeric_limits<double>::lowest();

    notice("Reading CSV...");
    const uint32_t extent = 4096;
    while(std::getline(file, line)) {
        std::stringstream css(line);
        std::vector<std::string> tokens;
        while(std::getline(css, token, delimiter)) {
            tokens.push_back(token);
        }
        if (tokens.size() <= (size_t)std::max(col_x, col_y)) continue;

        try {
            double x = std::stod(tokens[col_x]);
            double y = std::stod(tokens[col_y]);

            int64_t tx, ty;
            double lx, ly;
            pmt_utils::epsg3857totilecoord(x, y, zoom, &tx, &ty, &lx, &ly);

            // lx, ly are in pixel coordinates [0, 256], convert to tile extent [0, extent]
            int32_t px = std::round(lx * extent / 256.0);
            int32_t py = std::round(ly * extent / 256.0);

            uint64_t tile_id = pmtiles::zxy_to_tileid(zoom, tx, ty);

            PointFeature pf;
            pf.x = px;
            pf.y = py;
            pf.attrs.resize(attr_col_indices.size());
            for (size_t ai = 0; ai < attr_col_indices.size(); ++ai) {
                int ci = attr_col_indices[ai];
                if (ci < (int)tokens.size()) {
                    pf.attrs[ai] = tokens[ci];
                }
            }
            tile_points[tile_id].push_back(std::move(pf));
            point_count++;

            min_x = std::min(min_x, x);
            min_y = std::min(min_y, y);
            max_x = std::max(max_x, x);
            max_y = std::max(max_y, y);

        } catch(...) {
            continue;
        }
    }

    notice("Read %llu valid points into %zu tiles.", point_count, tile_points.size());

    // Detect column types and nullability.
    // Type detection: start as INT, downgrade to FLOAT if non-integer, then STRING if non-numeric.
    // Nullability: a column is nullable if any value is missing (empty or "NA").
    // Missing values do not influence type detection.
    std::vector<int>  attr_col_types   (attr_col_names.size(), COL_TYPE_INT);
    std::vector<bool> attr_col_nullable(attr_col_names.size(), false);

    for (auto& kv : tile_points) {
        for (auto& feature : kv.second) {
            for (size_t c = 0; c < feature.attrs.size(); ++c) {
                const std::string& v = feature.attrs[c];
                if (is_missing(v)) {
                    attr_col_nullable[c] = true;
                    continue;
                }
                if (attr_col_types[c] == COL_TYPE_STRING) continue;
                if (attr_col_types[c] == COL_TYPE_INT) {
                    try {
                        size_t pos;
                        std::stoll(v, &pos);
                        if (pos != v.size()) throw std::invalid_argument("trailing chars");
                    } catch (...) {
                        try {
                            std::stod(v);
                            attr_col_types[c] = COL_TYPE_FLOAT;
                        } catch (...) {
                            attr_col_types[c] = COL_TYPE_STRING;
                        }
                    }
                } else { // COL_TYPE_FLOAT
                    try {
                        std::stod(v);
                    } catch (...) {
                        attr_col_types[c] = COL_TYPE_STRING;
                    }
                }
            }
        }
    }

    for (size_t c = 0; c < attr_col_names.size(); ++c) {
        const char* type_str = attr_col_types[c] == COL_TYPE_INT   ? "int"
                             : attr_col_types[c] == COL_TYPE_FLOAT ? "float"
                                                                    : "string";
        notice("Column '%s': %s%s", attr_col_names[c].c_str(), type_str,
               attr_col_nullable[c] ? " (nullable)" : "");
    }

    uint64_t max_points = 0;
    uint64_t center_tile_id = 0;
    for (const auto& kv : tile_points) {
        if (kv.second.size() > max_points) {
            max_points = kv.second.size();
            center_tile_id = kv.first;
        }
    }

    if (n_threads < 1) n_threads = 1;

    // Flatten tile_points map into a vector for indexed parallel access
    std::vector<std::pair<uint64_t, std::vector<PointFeature>>> tile_vec(
        tile_points.begin(), tile_points.end());
    tile_points.clear(); // free memory

    // Generate PMTiles archive
    std::string tmp_file = tmp_dir + "/pmpoint_mlt_" + std::to_string(getpid()) + ".tmp";
    int out_fd = open(tmp_file.c_str(), O_RDWR | O_CREAT | O_TRUNC, 0666);
    if (out_fd < 0) error("Failed to open temporary file: %s", tmp_file.c_str());

    // Per-tile compressed output
    std::vector<std::string> compressed_tiles(tile_vec.size());

    notice("Encoding MLT tiles with %d thread(s)...", n_threads);

    if (n_threads <= 1) {
        MLTPointEncoder encoder;
        for (size_t i = 0; i < tile_vec.size(); ++i) {
            std::string uncompressed_mlt = encoder.encode(
                extent, base_name, tile_vec[i].second,
                attr_col_names, attr_col_types, attr_col_nullable);
            compressed_tiles[i] = mlt_gzip_compress(uncompressed_mlt);
        }
    } else {
        std::atomic<size_t> next_tile(0);
        std::vector<std::thread> threads;
        for (int t = 0; t < n_threads; ++t) {
            threads.emplace_back([&]() {
                MLTPointEncoder encoder;
                while (true) {
                    size_t i = next_tile.fetch_add(1);
                    if (i >= tile_vec.size()) break;
                    std::string uncompressed_mlt = encoder.encode(
                        extent, base_name, tile_vec[i].second,
                        attr_col_names, attr_col_types, attr_col_nullable);
                    compressed_tiles[i] = mlt_gzip_compress(uncompressed_mlt);
                }
            });
        }
        for (auto& th : threads) th.join();
    }

    // Write compressed tiles sequentially and build directory entries
    std::vector<pmtiles::entryv3> final_entries;
    uint64_t current_out_offset = 0;
    for (size_t i = 0; i < tile_vec.size(); ++i) {
        pwrite(out_fd, compressed_tiles[i].data(), compressed_tiles[i].size(), current_out_offset);
        final_entries.emplace_back(tile_vec[i].first, current_out_offset, compressed_tiles[i].size(), 1);
        current_out_offset += compressed_tiles[i].size();
        compressed_tiles[i].clear();
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
