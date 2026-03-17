#include "ext/nlohmann/json.hpp"
#include "pmpoint.h"
#include "qgenlib/params.h"
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
#include <random>
#include <fcntl.h>
#include <unistd.h>
#include <zlib.h>
#include <sys/stat.h>

#include "pmt_pts.h"
#include "pmt_utils.h"
#include "polygon.h"
#include "mvt_pts.h"
#include <cmath>
#include "htslib/hts.h"
#include "ext/protozero/pbf_writer.hpp"
#include "ext/protozero/pbf_reader.hpp"
#include "ext/protozero/varint.hpp"
#include "ext/PMTiles/pmtiles.hpp"

// Define a safe compress wrapper
std::string gzip_compress(const std::string& data) {
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
    if (ret != Z_STREAM_END) {
        error("deflate failed");
    }
    return outstring;
}

std::string gzip_decompress(const std::string& data) {
    z_stream zs;
    memset(&zs, 0, sizeof(zs));
    if (inflateInit2(&zs, 15 | 32) != Z_OK) {
        error("inflateInit2 failed");
    }
    zs.next_in = (Bytef*)data.data();
    zs.avail_in = data.size();
    
    int ret;
    char outbuffer[32768];
    std::string outstring;
    
    do {
        zs.next_out = reinterpret_cast<Bytef*>(outbuffer);
        zs.avail_out = sizeof(outbuffer);
        ret = inflate(&zs, 0);
        if (outstring.size() < zs.total_out) {
            outstring.append(outbuffer, zs.total_out - outstring.size());
        }
    } while (ret == Z_OK);
    
    inflateEnd(&zs);
    if (ret != Z_STREAM_END && ret != Z_OK) {
        error("inflate failed");
    }
    return outstring;
}

struct fast_feature {
    uint64_t id;
    pmt_utils::pmt_pt_t pt;
    std::vector<uint32_t> tags;
    
    fast_feature(pmt_utils::pmt_pt_t p) : id(0), pt(p) {}
};

struct fast_mvt {
    std::vector<std::string> keys;
    std::vector<std::string> values;
    std::vector<fast_feature> features;
};

void decode_mvt_raw(const std::string& buffer, uint8_t z, uint32_t x, uint32_t y, fast_mvt& out) {
    protozero::pbf_reader tile_reader(buffer);
    double scale_factor = pmt_utils::epsg3857_scale_factor(z);
    double offset_x, offset_y;
    pmt_utils::tiletoepsg3857(x, y, z, &offset_x, &offset_y);

    while (tile_reader.next(3)) { // Layer
        protozero::pbf_reader layer_reader = tile_reader.get_message();
        
        while (layer_reader.next()) {
            switch(layer_reader.tag()) {
                case 3: // keys
                    out.keys.push_back(std::string(layer_reader.get_view()));
                    break;
                case 4: // values
                    out.values.push_back(std::string(layer_reader.get_view()));
                    break;
                case 2: { // feature
                    protozero::pbf_reader feature_reader = layer_reader.get_message();
                    fast_feature f(pmt_utils::pmt_pt_t(z, 0, 0));
                    uint32_t geom_type = 0;
                    
                    int32_t cx = 0, cy = 0;
                    
                    while (feature_reader.next()) {
                        switch(feature_reader.tag()) {
                            case 1: f.id = feature_reader.get_uint64(); break;
                            case 2: { // tags
                                auto pi = feature_reader.get_packed_uint32();
                                for (auto it = pi.begin(); it != pi.end(); ++it) {
                                    f.tags.push_back(*it);
                                }
                                break;
                            }
                            case 3: geom_type = feature_reader.get_enum(); break;
                            case 4: { // geometry
                                auto pi = feature_reader.get_packed_uint32();
                                auto it = pi.begin();
                                if (it != pi.end()) {
                                    uint32_t cmd_len = *it++;
                                    if ((cmd_len & 7) == 1 && it != pi.end()) { // MoveTo
                                        cx += protozero::decode_zigzag32(*it++);
                                        if (it != pi.end()) cy += protozero::decode_zigzag32(*it++);
                                    }
                                }
                                break;
                            }
                            default: feature_reader.skip();
                        }
                    }
                    if (geom_type == 1) { // Point
                        f.pt.global_x = offset_x + scale_factor * cx;
                        f.pt.global_y = offset_y - scale_factor * cy;
                        out.features.push_back(std::move(f));
                    }
                    break;
                }
                default: layer_reader.skip();
            }
        }
    }
}

// ============================================================
// MLT-specific types and helpers for pyramid building
// ============================================================

static const int MLT_TYPE_STRING = 0;
static const int MLT_TYPE_FLOAT  = 1;
static const int MLT_TYPE_INT    = 2;

struct mlt_feature_pyr {
    double global_x, global_y;
    std::vector<std::string> attrs;
};

struct mlt_layer_pyr {
    std::string name;
    uint32_t extent;
    std::vector<std::string> col_names;
    std::vector<int>  col_types;
    std::vector<bool> col_nullable;
    std::vector<mlt_feature_pyr> features;
};

// ORC-style boolean RLE decoder (PRESENT stream)
static std::vector<bool> decode_bool_rle_pyr(const uint8_t* data, size_t len, size_t count) {
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

// ORC-style boolean RLE encoder (PRESENT stream)
static std::string encode_bool_rle_pyr(const std::vector<bool>& present) {
    size_t n = present.size();
    size_t numBytes = (n + 7) / 8;
    std::vector<uint8_t> packed(numBytes, 0);
    for (size_t i = 0; i < n; ++i)
        if (present[i]) packed[i/8] |= static_cast<uint8_t>(1 << (i%8));
    std::string result;
    size_t offset = 0;
    while (offset < numBytes) {
        size_t chunk = std::min(static_cast<size_t>(128), numBytes - offset);
        result.push_back(static_cast<char>(static_cast<uint8_t>(256 - chunk)));
        for (size_t i = 0; i < chunk; ++i)
            result.push_back(static_cast<char>(packed[offset + i]));
        offset += chunk;
    }
    return result;
}

// Quick feature count from compressed MLT tile (decompress + parse geometry header only)
static size_t count_mlt_features_quick(const std::string& uncompressed) {
    if (uncompressed.empty()) return 0;
    const uint8_t* ptr = (const uint8_t*)uncompressed.data();
    const uint8_t* end = ptr + uncompressed.size();
    auto read_varint = [&]() -> uint64_t {
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
        uint64_t layer_len = read_varint();
        if (layer_len == 0 || ptr >= end) break;
        uint8_t tag = *ptr++;
        if (tag != 1) { ptr += layer_len - 1; continue; }
        uint64_t name_len = read_varint();
        ptr += name_len;
        read_varint(); // extent
        uint64_t num_columns = read_varint();
        for (uint64_t c = 0; c < num_columns; ++c) {
            uint64_t tc = read_varint();
            if (tc >= 10) { uint64_t cl = read_varint(); ptr += cl; }
        }
        uint64_t geom_num_streams = read_varint();
        for (uint64_t s = 0; s < geom_num_streams; ++s) {
            if (ptr + 2 > end) return 0;
            uint8_t h0 = *ptr++; ptr++;
            uint64_t num_vals = read_varint();
            uint64_t byte_len = read_varint();
            ptr += byte_len;
            if (((h0 >> 4) & 0x0F) == 1 && (h0 & 0x0F) == 3)
                return (size_t)(num_vals / 2);
        }
        return 0;
    }
    return 0;
}

static size_t count_mvt_features_quick(const std::string& uncompressed) {
    size_t count = 0;
    protozero::pbf_reader tile_reader(uncompressed);
    while (tile_reader.next(3)) {
        protozero::pbf_reader layer_reader = tile_reader.get_message();
        while (layer_reader.next()) {
            if (layer_reader.tag() == 2) ++count;
            layer_reader.skip();
        }
    }
    return count;
}

// Decode an MLT tile (new multi-column format) into mlt_layer_pyr
void decode_mlt_layer(const std::string& buffer, uint8_t z, uint32_t x, uint32_t y, mlt_layer_pyr& layer) {
    if (buffer.empty()) return;
    double scale_factor = pmt_utils::epsg3857_scale_factor(z);
    double offset_x, offset_y;
    pmt_utils::tiletoepsg3857(x, y, z, &offset_x, &offset_y);

    const uint8_t* ptr = (const uint8_t*)buffer.data();
    const uint8_t* end = ptr + buffer.size();
    auto read_varint = [&]() -> uint64_t {
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
        uint64_t layer_len = read_varint();
        if (layer_len == 0 || ptr >= end) break;
        uint8_t tag = *ptr++;
        if (tag != 1) { ptr += layer_len - 1; continue; }

        // --- Layer header ---
        uint64_t name_len = read_varint();
        layer.name = std::string((char*)ptr, name_len);
        ptr += name_len;
        layer.extent = (uint32_t)read_varint();
        uint64_t num_columns = read_varint();

        // --- METADATA SECTION: read all column typeCodes + names ---
        struct ColMeta { uint64_t typeCode; std::string name; };
        std::vector<ColMeta> col_metas(num_columns);
        for (uint64_t c = 0; c < num_columns; ++c) {
            col_metas[c].typeCode = read_varint();
            if (col_metas[c].typeCode >= 10) {
                uint64_t cname_len = read_varint();
                col_metas[c].name = std::string((char*)ptr, cname_len);
                ptr += cname_len;
            }
        }

        // Set up schema for attribute columns (column 0 is always GEOMETRY)
        size_t num_attr = num_columns > 0 ? num_columns - 1 : 0;
        layer.col_names.resize(num_attr);
        layer.col_types.resize(num_attr);
        layer.col_nullable.resize(num_attr);
        for (size_t c = 0; c < num_attr; ++c) {
            uint64_t tc = col_metas[c + 1].typeCode;
            layer.col_names[c] = col_metas[c + 1].name;
            layer.col_nullable[c] = (tc % 2 == 1);
            uint64_t base = tc - (tc % 2);
            if      (base >= 16 && base <= 23) layer.col_types[c] = MLT_TYPE_INT;
            else if (base >= 24 && base <= 27) layer.col_types[c] = MLT_TYPE_FLOAT;
            else                               layer.col_types[c] = MLT_TYPE_STRING;
        }

        // --- DATA SECTION ---
        // GEOMETRY column: has numStreams prefix
        uint64_t geom_num_streams = read_varint();
        size_t num_features = 0;
        std::vector<int32_t> vx, vy;

        for (uint64_t s = 0; s < geom_num_streams; ++s) {
            if (ptr + 2 > end) break;
            uint8_t h0 = *ptr++;
            uint8_t h1 = *ptr++; (void)h1;
            uint64_t num_vals = read_varint();
            uint64_t byte_len = read_varint();
            const uint8_t* sd = ptr;
            ptr += byte_len;
            uint8_t phys = (h0 >> 4) & 0x0F;
            uint8_t dict = h0 & 0x0F;
            if (phys == 1 && dict == 3) { // VERTEX stream
                num_features = (size_t)(num_vals / 2);
                vx.resize(num_features);
                vy.resize(num_features);
                const uint8_t* vp = sd;
                for (size_t i = 0; i < num_features; ++i) {
                    uint64_t zx = 0; int sh = 0;
                    while (vp < sd+byte_len) { uint8_t b=*vp++; zx|=(uint64_t)(b&0x7F)<<sh; sh+=7; if(!(b&0x80)) break; }
                    uint64_t zy = 0; sh = 0;
                    while (vp < sd+byte_len) { uint8_t b=*vp++; zy|=(uint64_t)(b&0x7F)<<sh; sh+=7; if(!(b&0x80)) break; }
                    vx[i] = (int32_t)((zx>>1) ^ -(int64_t)(zx&1));
                    vy[i] = (int32_t)((zy>>1) ^ -(int64_t)(zy&1));
                }
            }
        }

        layer.features.resize(num_features);
        for (size_t i = 0; i < num_features; ++i) {
            layer.features[i].global_x = offset_x + scale_factor * vx[i];
            layer.features[i].global_y = offset_y - scale_factor * vy[i];
            layer.features[i].attrs.assign(num_attr, std::string());
        }

        // Attribute columns
        for (size_t c = 0; c < num_attr; ++c) {
            bool nullable = layer.col_nullable[c];
            int ctype = layer.col_types[c];
            bool is_str = (ctype == MLT_TYPE_STRING);
            std::vector<bool> present(num_features, true);
            std::vector<uint64_t> str_lens;
            const uint8_t* str_data = nullptr;
            uint64_t str_data_len = 0;

            uint64_t ns = is_str ? read_varint() : (nullable ? 2 : 1);
            for (uint64_t s = 0; s < ns; ++s) {
                if (ptr + 2 > end) break;
                uint8_t h0 = *ptr++;
                uint8_t h1 = *ptr++; (void)h1;
                uint64_t nv = read_varint();
                uint64_t bl = read_varint();
                const uint8_t* sd = ptr;
                ptr += bl;
                uint8_t phys = (h0 >> 4) & 0x0F;

                if (phys == 0) { // PRESENT
                    present = decode_bool_rle_pyr(sd, bl, num_features);
                } else if (phys == 1) { // DATA
                    if (ctype == MLT_TYPE_INT) {
                        const uint8_t* dp = sd;
                        size_t fi = 0;
                        for (uint64_t vi = 0; vi < nv; ++vi) {
                            uint64_t zig = 0; int sh = 0;
                            while (dp < sd+bl) { uint8_t b=*dp++; zig|=(uint64_t)(b&0x7F)<<sh; sh+=7; if(!(b&0x80)) break; }
                            int64_t val = (int64_t)((zig>>1) ^ -(int64_t)(zig&1));
                            while (fi < num_features && !present[fi]) ++fi;
                            if (fi < num_features) layer.features[fi++].attrs[c] = std::to_string(val);
                        }
                    } else if (ctype == MLT_TYPE_FLOAT) {
                        const uint8_t* dp = sd;
                        size_t fi = 0;
                        for (uint64_t vi = 0; vi < nv; ++vi) {
                            uint32_t bits = (uint32_t)dp[0]|((uint32_t)dp[1]<<8)|((uint32_t)dp[2]<<16)|((uint32_t)dp[3]<<24);
                            dp += 4;
                            float fval; memcpy(&fval, &bits, 4);
                            while (fi < num_features && !present[fi]) ++fi;
                            if (fi < num_features) {
                                char buf[32]; snprintf(buf, sizeof(buf), "%.9g", (double)fval);
                                layer.features[fi++].attrs[c] = buf;
                            }
                        }
                    } else { // STRING DATA
                        str_data = sd; str_data_len = bl; (void)nv;
                    }
                } else if (phys == 3) { // LENGTH (string lengths)
                    const uint8_t* dp = sd;
                    str_lens.reserve(nv);
                    for (uint64_t vi = 0; vi < nv; ++vi) {
                        uint64_t len = 0; int sh = 0;
                        while (dp < sd+bl) { uint8_t b=*dp++; len|=(uint64_t)(b&0x7F)<<sh; sh+=7; if(!(b&0x80)) break; }
                        str_lens.push_back(len);
                    }
                }
            }

            // Decode strings after all streams consumed
            if (is_str && str_data && !str_lens.empty()) {
                const uint8_t* dp = str_data;
                size_t fi = 0;
                for (size_t li = 0; li < str_lens.size(); ++li) {
                    while (fi < num_features && !present[fi]) ++fi;
                    if (fi < num_features) {
                        layer.features[fi++].attrs[c] = std::string((char*)dp, str_lens[li]);
                        dp += str_lens[li];
                    }
                }
            }
            (void)str_data_len;
        }

        break; // only first layer
    }
}

// Encode an mlt_layer_pyr into MLT tile bytes
std::string encode_mlt_layer(const mlt_layer_pyr& layer, const std::vector<int>& indices,
                              uint8_t z, uint32_t x, uint32_t y) {
    double scale_factor = pmt_utils::epsg3857_scale_factor(z);
    double offset_x, offset_y;
    pmt_utils::tiletoepsg3857(x, y, z, &offset_x, &offset_y);

    size_t n = indices.size();
    size_t num_attr = layer.col_names.size();
    const std::string& lname = layer.name;
    uint32_t extent = layer.extent > 0 ? layer.extent : 4096;

    std::string tmp;
    auto av = [&](uint64_t v) {
        do { uint8_t b=v&0x7F; v>>=7; if(v) b|=0x80; tmp.push_back(b); } while(v);
    };

    // METADATA
    av(lname.size()); tmp.append(lname);
    av(extent);
    av(1 + num_attr);
    av(4); // GEOMETRY, no name
    for (size_t c = 0; c < num_attr; ++c) {
        int tc;
        if      (layer.col_types[c] == MLT_TYPE_INT)   tc = layer.col_nullable[c] ? 17 : 16;
        else if (layer.col_types[c] == MLT_TYPE_FLOAT)  tc = layer.col_nullable[c] ? 25 : 24;
        else                                             tc = layer.col_nullable[c] ? 29 : 28;
        av(tc);
        av(layer.col_names[c].size()); tmp.append(layer.col_names[c]);
    }

    // DATA: GEOMETRY
    av(2); // numStreams
    // GeomType stream: 0x10=DATA/NONE, 0x02=VARINT, numValues=n, byteLen=n (each 0=varint POINT)
    tmp.push_back(0x10); tmp.push_back(0x02);
    av(n); av(n);
    for (size_t i = 0; i < n; ++i) tmp.push_back(0);
    // Vertex stream
    std::string vdata;
    auto vv = [&](uint64_t v) {
        do { uint8_t b=v&0x7F; v>>=7; if(v) b|=0x80; vdata.push_back(b); } while(v);
    };
    for (size_t i = 0; i < n; ++i) {
        int fi = indices[i];
        int32_t px = (int32_t)std::round((layer.features[fi].global_x - offset_x) / scale_factor);
        int32_t py = (int32_t)std::round((offset_y - layer.features[fi].global_y) / scale_factor);
        vv(((uint32_t)px << 1) ^ (uint32_t)(px >> 31));
        vv(((uint32_t)py << 1) ^ (uint32_t)(py >> 31));
    }
    tmp.push_back(0x13); tmp.push_back(0x02);
    av(n * 2); av(vdata.size()); tmp.append(vdata);

    // DATA: attribute columns
    for (size_t c = 0; c < num_attr; ++c) {
        bool nullable = layer.col_nullable[c];
        int ctype = layer.col_types[c];
        std::vector<bool> present(n, true);
        if (nullable) {
            for (size_t i = 0; i < n; ++i) {
                int fi = indices[i];
                const std::string& v = (c < layer.features[fi].attrs.size()) ? layer.features[fi].attrs[c] : "";
                present[i] = !v.empty();
            }
        }
        size_t nn = 0;
        for (size_t i = 0; i < n; ++i) if (present[i]) ++nn;

        if (ctype == MLT_TYPE_INT) {
            std::string idata;
            auto iv = [&](uint64_t v) {
                do { uint8_t b=v&0x7F; v>>=7; if(v) b|=0x80; idata.push_back(b); } while(v);
            };
            for (size_t i = 0; i < n; ++i) {
                if (!present[i]) continue;
                int fi = indices[i];
                const std::string& s = (c < layer.features[fi].attrs.size()) ? layer.features[fi].attrs[c] : "0";
                int32_t val = 0; try { val = std::stoi(s); } catch (...) {}
                iv(((uint32_t)val << 1) ^ (uint32_t)(val >> 31));
            }
            if (nullable) {
                std::string rle = encode_bool_rle_pyr(present);
                tmp.push_back(0x00); tmp.push_back(0x02); av(n); av(rle.size()); tmp.append(rle);
            }
            tmp.push_back(0x10); tmp.push_back(0x02); av(nn); av(idata.size()); tmp.append(idata);
        } else if (ctype == MLT_TYPE_FLOAT) {
            std::string fdata;
            for (size_t i = 0; i < n; ++i) {
                if (!present[i]) continue;
                int fi = indices[i];
                const std::string& s = (c < layer.features[fi].attrs.size()) ? layer.features[fi].attrs[c] : "0";
                float val = 0.0f; try { val = std::stof(s); } catch (...) {}
                uint32_t bits; memcpy(&bits, &val, 4);
                fdata.push_back(bits&0xFF); fdata.push_back((bits>>8)&0xFF);
                fdata.push_back((bits>>16)&0xFF); fdata.push_back((bits>>24)&0xFF);
            }
            if (nullable) {
                std::string rle = encode_bool_rle_pyr(present);
                tmp.push_back(0x00); tmp.push_back(0x02); av(n); av(rle.size()); tmp.append(rle);
            }
            tmp.push_back(0x10); tmp.push_back(0x00); av(nn); av(fdata.size()); tmp.append(fdata);
        } else { // STRING
            std::string ldata, sdata;
            auto lv = [&](uint64_t v) {
                do { uint8_t b=v&0x7F; v>>=7; if(v) b|=0x80; ldata.push_back(b); } while(v);
            };
            for (size_t i = 0; i < n; ++i) {
                if (!present[i]) continue;
                int fi = indices[i];
                const std::string& s = (c < layer.features[fi].attrs.size()) ? layer.features[fi].attrs[c] : "";
                lv(s.size()); sdata.append(s);
            }
            av(nullable ? 3 : 2);
            if (nullable) {
                std::string rle = encode_bool_rle_pyr(present);
                tmp.push_back(0x00); tmp.push_back(0x02); av(n); av(rle.size()); tmp.append(rle);
            }
            tmp.push_back(0x30); tmp.push_back(0x02); av(nn); av(ldata.size()); tmp.append(ldata);
            tmp.push_back(0x10); tmp.push_back(0x00); av(0); av(sdata.size()); tmp.append(sdata);
        }
    }

    // Layer wrapper
    std::string out;
    auto ov = [&](uint64_t v) {
        do { uint8_t b=v&0x7F; v>>=7; if(v) b|=0x80; out.push_back(b); } while(v);
    };
    ov(1 + tmp.size());
    out.push_back(1);
    out.append(tmp);
    return out;
}

std::string encode_mvt(const fast_mvt& df, const std::vector<int>& indices, uint8_t z, uint32_t x, uint32_t y, const std::string& layer_name) {
    std::string tile_data;
    protozero::pbf_writer tile_writer(tile_data);
    
    std::string layer_data;
    protozero::pbf_writer layer_writer(layer_data);
    
    layer_writer.add_uint32(15, 2); // version
    layer_writer.add_string(1, layer_name);
    layer_writer.add_uint32(5, 4096); // extent
    
    // Build a map of used keys and values to prune unused ones
    std::vector<uint32_t> key_remap(df.keys.size(), -1);
    std::vector<uint32_t> val_remap(df.values.size(), -1);
    
    std::vector<std::string> final_keys;
    std::vector<std::string> final_values;
    
    for(size_t idx=0; idx<indices.size(); ++idx) {
        size_t i = indices[idx];
        for(size_t t=0; t<df.features[i].tags.size(); t+=2) {
            uint32_t k = df.features[i].tags[t];
            uint32_t v = df.features[i].tags[t+1];
            if (key_remap[k] == (uint32_t)-1) {
                key_remap[k] = final_keys.size();
                final_keys.push_back(df.keys[k]);
            }
            if (val_remap[v] == (uint32_t)-1) {
                val_remap[v] = final_values.size();
                final_values.push_back(df.values[v]);
            }
        }
    }
    
    // Write pruned keys
    for(size_t i=0; i<final_keys.size(); ++i) {
        layer_writer.add_string(3, final_keys[i]);
    }
    
    // Write pruned values directly by copying the raw proto message view
    for(size_t i=0; i<final_values.size(); ++i) {
        layer_writer.add_message(4, final_values[i]);
    }
    
    // features
    double scale_factor = pmt_utils::epsg3857_scale_factor(z);
    double offset_x, offset_y;
    pmt_utils::tiletoepsg3857(x, y, z, &offset_x, &offset_y);
    
    for(size_t idx=0; idx<indices.size(); ++idx) {
        size_t i = indices[idx];
        std::string feature_data;
        protozero::pbf_writer feature_writer(feature_data);
        
        feature_writer.add_uint64(1, idx+1); // id
        feature_writer.add_enum(3, 1); // Point geometry
        
        std::vector<uint32_t> tags;
        for(size_t t=0; t<df.features[i].tags.size(); t+=2) {
            tags.push_back(key_remap[df.features[i].tags[t]]);
            tags.push_back(val_remap[df.features[i].tags[t+1]]);
        }
        if (!tags.empty()) {
            feature_writer.add_packed_uint32(2, tags.begin(), tags.end());
        }
        
        // Geometry
        std::vector<uint32_t> geom;
        geom.push_back((1 << 3) | 1); // MoveTo length 1
        
        int32_t px = std::round((df.features[i].pt.global_x - offset_x) / scale_factor);
        int32_t py = std::round((offset_y - df.features[i].pt.global_y) / scale_factor);
        
        geom.push_back(protozero::encode_zigzag32(px));
        geom.push_back(protozero::encode_zigzag32(py));
        
        feature_writer.add_packed_uint32(4, geom.begin(), geom.end());
        
        layer_writer.add_message(2, feature_data);
    }
    
    tile_writer.add_message(3, layer_data);
    return tile_data;
}

size_t estimate_uncompressed_mvt_size(const fast_mvt& df, const std::vector<int>& indices, const std::string& layer_name) {
    size_t size = 50 + layer_name.size();
    
    std::unordered_map<uint32_t, bool> key_seen;
    std::unordered_map<uint32_t, bool> val_seen;
    
    size_t tag_total = 0;
    
    for(size_t idx=0; idx<indices.size(); ++idx) {
        size_t i = indices[idx];
        tag_total += df.features[i].tags.size();
        for(size_t t=0; t<df.features[i].tags.size(); t+=2) {
            key_seen[df.features[i].tags[t]] = true;
            val_seen[df.features[i].tags[t+1]] = true;
        }
    }
    
    for(const auto& kv : key_seen) {
        size += df.keys[kv.first].size() + 2;
    }
    for(const auto& kv : val_seen) {
        size += df.values[kv.first].size() + 2;
    }
    
    size += indices.size() * 10 + tag_total;
    return size;
}

void get_subset_indices(const fast_mvt& df, uint64_t max_features, uint32_t seed, std::vector<int>& indices) {
    size_t total_points = df.features.size();
    indices.clear();
    if (total_points <= max_features) {
        indices.resize(total_points);
        for (size_t i = 0; i < total_points; ++i) indices[i] = i;
        return;
    }
    std::mt19937 rng(seed);
    indices.resize(total_points);
    for (size_t i = 0; i < total_points; ++i) indices[i] = i;
    std::shuffle(indices.begin(), indices.end(), rng);
    indices.resize(max_features);
    std::sort(indices.begin(), indices.end());
}

void get_subset_indices_mlt(const mlt_layer_pyr& layer, uint64_t max_features, uint32_t seed, std::vector<int>& indices) {
    size_t total_points = layer.features.size();
    indices.clear();
    if (total_points <= max_features) {
        indices.resize(total_points);
        for (size_t i = 0; i < total_points; ++i) indices[i] = i;
        return;
    }
    std::mt19937 rng(seed);
    indices.resize(total_points);
    for (size_t i = 0; i < total_points; ++i) indices[i] = i;
    std::shuffle(indices.begin(), indices.end(), rng);
    indices.resize(max_features);
    std::sort(indices.begin(), indices.end());
}

// Global state for build process
std::mutex out_mutex;
int out_fd = -1;
uint64_t current_out_offset = 0;
std::vector<pmtiles::entryv3> final_entries;
std::string layer_name_global = "data";
std::map<uint64_t, size_t> tile_feature_counts; // track feature counts for uniform subsampling

class PyramidBuilderQueue {
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
        if (work_queue.empty() && done) return false;
        entry = work_queue.front();
        work_queue.pop();
        return true;
    }
    void finish() {
        std::lock_guard<std::mutex> lock(mutex);
        done = true;
        cv.notify_all();
    }
    size_t size() {
        std::lock_guard<std::mutex> lock(mutex);
        return work_queue.size();
    }
};

// Holds decoded parent tile data between pass 1 and pass 2
struct parent_tile_data {
    uint32_t z, x, y;
    mlt_layer_pyr mlt_combined;
    fast_mvt mvt_combined;
    size_t post_density_count;       // raw combined feature count
    size_t estimated_uncompressed;   // estimated uncompressed tile size in bytes
};

// Pass 1: decode child tiles, combine, record feature counts
void builder_worker_pass1(PyramidBuilderQueue& queue,
                          const std::map<uint64_t, pmtiles::entryv3>& next_level_entries,
                          uint8_t tile_type,
                          std::vector<parent_tile_data>& results, std::mutex& results_mutex) {

    pmtiles::entry_zxy entry(0,0,0,0,0);
    while (queue.get_tile(entry)) {
        uint32_t z = entry.z;
        uint32_t x = entry.x;
        uint32_t y = entry.y;

        parent_tile_data ptd;
        ptd.z = z; ptd.x = x; ptd.y = y;
        ptd.post_density_count = 0;

        if (tile_type == 0x06) {
            bool schema_set = false;
            for (int dy = 0; dy < 2; ++dy) {
                for (int dx = 0; dx < 2; ++dx) {
                    uint32_t cz = z + 1, cx = 2*x+dx, cy = 2*y+dy;
                    uint64_t c_id = pmtiles::zxy_to_tileid(cz, cx, cy);
                    auto it = next_level_entries.find(c_id);
                    if (it == next_level_entries.end()) continue;

                    std::string compressed(it->second.length, '\0');
                    pread(out_fd, &compressed[0], it->second.length, it->second.offset);
                    std::string uncompressed = gzip_decompress(compressed);

                    mlt_layer_pyr child;
                    decode_mlt_layer(uncompressed, cz, cx, cy, child);

                    if (!schema_set && !child.features.empty()) {
                        ptd.mlt_combined.name = child.name;
                        ptd.mlt_combined.extent = child.extent;
                        ptd.mlt_combined.col_names = child.col_names;
                        ptd.mlt_combined.col_types = child.col_types;
                        ptd.mlt_combined.col_nullable = child.col_nullable;
                        schema_set = true;
                    }
                    for (auto& f : child.features)
                        ptd.mlt_combined.features.push_back(std::move(f));
                }
            }
            if (ptd.mlt_combined.features.empty()) continue;

            // Record raw combined feature count and estimated size
            ptd.post_density_count = ptd.mlt_combined.features.size();
            ptd.estimated_uncompressed = ptd.post_density_count * (20 + ptd.mlt_combined.col_names.size() * 8);

        } else {
            std::map<std::string, uint32_t> global_key_map;
            std::map<std::string, uint32_t> global_val_map;
            for (int dy = 0; dy < 2; ++dy) {
                for (int dx = 0; dx < 2; ++dx) {
                    uint32_t cz = z + 1, cx = 2*x+dx, cy = 2*y+dy;
                    uint64_t c_id = pmtiles::zxy_to_tileid(cz, cx, cy);
                    auto it = next_level_entries.find(c_id);
                    if (it == next_level_entries.end()) continue;

                    std::string compressed(it->second.length, '\0');
                    pread(out_fd, &compressed[0], it->second.length, it->second.offset);
                    std::string uncompressed = gzip_decompress(compressed);

                    fast_mvt child_mvt;
                    decode_mvt_raw(uncompressed, cz, cx, cy, child_mvt);

                    std::vector<uint32_t> k_remap(child_mvt.keys.size());
                    for (size_t i = 0; i < child_mvt.keys.size(); ++i) {
                        auto gk = global_key_map.find(child_mvt.keys[i]);
                        if (gk == global_key_map.end()) {
                            k_remap[i] = ptd.mvt_combined.keys.size();
                            global_key_map[child_mvt.keys[i]] = k_remap[i];
                            ptd.mvt_combined.keys.push_back(child_mvt.keys[i]);
                        } else {
                            k_remap[i] = gk->second;
                        }
                    }
                    std::vector<uint32_t> v_remap(child_mvt.values.size());
                    for (size_t i = 0; i < child_mvt.values.size(); ++i) {
                        auto gv = global_val_map.find(child_mvt.values[i]);
                        if (gv == global_val_map.end()) {
                            v_remap[i] = ptd.mvt_combined.values.size();
                            global_val_map[child_mvt.values[i]] = v_remap[i];
                            ptd.mvt_combined.values.push_back(child_mvt.values[i]);
                        } else {
                            v_remap[i] = gv->second;
                        }
                    }
                    for (auto& f : child_mvt.features) {
                        for (size_t t = 0; t < f.tags.size(); t += 2) {
                            f.tags[t]   = k_remap[f.tags[t]];
                            f.tags[t+1] = v_remap[f.tags[t+1]];
                        }
                        ptd.mvt_combined.features.push_back(std::move(f));
                    }
                }
            }
            if (ptd.mvt_combined.features.empty()) continue;

            // Record raw combined feature count and estimated size
            ptd.post_density_count = ptd.mvt_combined.features.size();
            // Estimate: each feature ~20 bytes geometry + tags overhead
            size_t avg_tags = 0;
            if (!ptd.mvt_combined.features.empty()) {
                for (size_t i = 0; i < std::min((size_t)100, ptd.mvt_combined.features.size()); ++i)
                    avg_tags += ptd.mvt_combined.features[i].tags.size();
                avg_tags = avg_tags / std::min((size_t)100, ptd.mvt_combined.features.size());
            }
            ptd.estimated_uncompressed = ptd.post_density_count * (20 + avg_tags * 4);
        }

        {
            std::lock_guard<std::mutex> lock(results_mutex);
            results.push_back(std::move(ptd));
        }
    }
}


int32_t cmd_build_pyramid_pmtiles(int32_t argc, char **argv)
{
    std::string in_pmtiles;
    std::string out_pmtiles;
    std::string tmp_dir = ".";
    int32_t min_zoom = 0;
    int32_t max_tile_bytes = 500000;
    int32_t max_tile_features = 50000;
    int32_t num_threads = 0;

    paramList pl;
    BEGIN_LONG_PARAMS(longParameters)
    LONG_STRING_PARAM("in", &in_pmtiles, "Input PMTiles file (finest level only)")
    LONG_STRING_PARAM("out", &out_pmtiles, "Output pyramidal PMTiles file")
    LONG_INT_PARAM("min-zoom", &min_zoom, "Minimum zoom level to construct")
    LONG_INT_PARAM("max-tile-bytes", &max_tile_bytes, "Maximum compressed tile bytes")
    LONG_INT_PARAM("max-tile-features", &max_tile_features, "Maximum features per tile")
    LONG_STRING_PARAM("tmp-dir", &tmp_dir, "Temporary directory for tile data")
    LONG_INT_PARAM("threads", &num_threads, "Number of threads")
    END_LONG_PARAMS();

    pl.Add(new longParams("Available Options", longParameters));
    pl.Read(argc, argv);
    pl.Status();

    if (in_pmtiles.empty() || out_pmtiles.empty()) error("Missing --in or --out");

    if (num_threads <= 0) {
        num_threads = std::thread::hardware_concurrency();
        if (num_threads == 0) num_threads = 4;
    }

    pmt_pts pmt(in_pmtiles.c_str());
    if (!pmt.read_header_meta_entries()) error("Invalid PMTiles file");

    uint8_t tile_type = pmt.hdr.tile_type;
    if (tile_type == pmtiles::TILETYPE_MVT) {
        notice("Input format: MVT-based PMTiles");
    } else if (tile_type == 0x06) {
        notice("Input format: MLT-based PMTiles");
    } else {
        error("Unsupported PMTiles format: tile_type=%d", tile_type);
    }

    uint8_t z_max = pmt.hdr.max_zoom;
    if (min_zoom > z_max) min_zoom = z_max;
    notice("Constructing pyramid from zoom %d down to %d", z_max, min_zoom);

    // Extract layer name from first tile
    if (pmt.tile_entries.size() > 0) {
        std::string buffer;
        pmt.fetch_tile_to_buffer(pmt.tile_entries[0].z, pmt.tile_entries[0].x, pmt.tile_entries[0].y, buffer);
        if (tile_type == 0x06) {
            const char* ptr = buffer.data();
            const char* end = buffer.data() + buffer.size();
            auto read_varint = [&]() -> uint64_t {
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
                uint64_t layer_len = read_varint();
                if (layer_len == 0 || ptr >= end) break;
                uint8_t tag = *ptr++;
                if (tag == 1) {
                    uint64_t name_len = read_varint();
                    layer_name_global = std::string(ptr, name_len);
                    break;
                }
                ptr += layer_len - 1;
            }
        } else {
            mapbox::vector_tile::buffer* p_tile = new mapbox::vector_tile::buffer(buffer);
            if (p_tile && p_tile->layerNames().size() > 0) {
                layer_name_global = p_tile->layerNames()[0];
            }
            delete p_tile;
        }
    }
    notice("Using layer name: %s", layer_name_global.c_str());

    std::string tmp_file = tmp_dir + "/pmpoint_pyramid_" + std::to_string(getpid()) + ".tmp";
    out_fd = open(tmp_file.c_str(), O_RDWR | O_CREAT | O_TRUNC, 0666);
    if (out_fd < 0) error("Failed to open temporary file: %s", tmp_file.c_str());

    std::map<uint64_t, pmtiles::entryv3> level_entries;

    size_t zmax_total_features = 0;
    notice("Copying z_max tiles...");
    for(size_t i=0; i<pmt.tile_entries.size(); ++i) {
        const auto& entry = pmt.tile_entries[i];
        if (entry.z == z_max) {
            std::string buffer;
            pmt.flex_reader_ptr->read_at(entry.offset, entry.length, buffer);

            // Count features for uniform subsampling tracking
            std::string uncompressed = gzip_decompress(buffer);
            size_t nf = 0;
            if (tile_type == 0x06)
                nf = count_mlt_features_quick(uncompressed);
            else
                nf = count_mvt_features_quick(uncompressed);
            zmax_total_features += nf;

            pwrite(out_fd, buffer.data(), buffer.size(), current_out_offset);

            uint64_t tile_id = pmtiles::zxy_to_tileid(entry.z, entry.x, entry.y);
            pmtiles::entryv3 e(tile_id, current_out_offset, buffer.size(), 1);
            final_entries.push_back(e);
            level_entries[tile_id] = e;
            tile_feature_counts[tile_id] = nf;
            current_out_offset += buffer.size();
        }
    }
    notice("  z%d: %zu tiles, %zu total features", z_max, level_entries.size(), zmax_total_features);

    for (int32_t z = z_max - 1; z >= min_zoom; --z) {
        notice("Building zoom level %d...", z);

        // Determine parent tiles from current level entries
        std::map<uint64_t, pmtiles::entry_zxy> parent_tiles;
        for (const auto& kv : level_entries) {
            pmtiles::zxy parent_zxy = pmtiles::tileid_to_zxy(kv.first);
            if (parent_zxy.z > 0) {
                uint32_t px = parent_zxy.x / 2;
                uint32_t py = parent_zxy.y / 2;
                uint64_t pid = pmtiles::zxy_to_tileid(z, px, py);
                if (parent_tiles.find(pid) == parent_tiles.end()) {
                    parent_tiles.insert({pid, pmtiles::entry_zxy(z, px, py, 0, 0)});
                }
            }
        }

        // === Pass 1: Decode all parent tiles, combine children, record raw feature counts ===
        std::vector<parent_tile_data> pass1_results;
        std::mutex pass1_mutex;
        {
            PyramidBuilderQueue queue;
            for (const auto& kv : parent_tiles) {
                queue.add_tile(kv.second);
            }
            queue.finish();

            std::vector<std::thread> threads;
            for (int32_t i = 0; i < num_threads; ++i) {
                threads.emplace_back(builder_worker_pass1, std::ref(queue), std::cref(level_entries),
                                     tile_type,
                                     std::ref(pass1_results), std::ref(pass1_mutex));
            }
            for (auto& t : threads) t.join();
        }

        // Compute per-level uniform subsampling ratio from BOTH feature count and byte size
        size_t max_raw_count = 0;
        size_t max_estimated_uncompressed = 0;
        for (const auto& ptd : pass1_results) {
            max_raw_count = std::max(max_raw_count, ptd.post_density_count);
            max_estimated_uncompressed = std::max(max_estimated_uncompressed, ptd.estimated_uncompressed);
        }

        // Ratio from feature count constraint
        double ratio_features = 1.0;
        if (max_raw_count > (size_t)max_tile_features) {
            ratio_features = (double)max_tile_features / max_raw_count;
        }

        // Ratio from byte size constraint (estimate ~4x gzip compression)
        double ratio_bytes = 1.0;
        if (max_estimated_uncompressed > (size_t)max_tile_bytes * 4) {
            ratio_bytes = (double)((size_t)max_tile_bytes * 4) / max_estimated_uncompressed;
        }

        // Use the more restrictive ratio
        double level_ratio = std::min(ratio_features, ratio_bytes);
        notice("  z%d: pass1 done, %zu tiles, max_raw=%zu, max_est_bytes=%zu, ratio_feat=%.6f, ratio_bytes=%.6f, level_ratio=%.6f",
               z, pass1_results.size(), max_raw_count, max_estimated_uncompressed,
               ratio_features, ratio_bytes, level_ratio);

        // === Pass 2: encode ALL tiles with level_ratio applied uniformly ===
        // Every tile is subsampled: keep = level_ratio * tile_feature_count
        // This ensures uniform subsampling across all tiles at this level
        size_t level_total_features = 0;
        {
            std::vector<std::thread> threads;
            size_t n = pass1_results.size();
            size_t per_thread = (n + num_threads - 1) / num_threads;
            for (int32_t i = 0; i < num_threads; ++i) {
                size_t s = i * per_thread;
                size_t e = std::min(s + per_thread, n);
                if (s >= n) break;
                threads.emplace_back([&pass1_results, s, e, max_tile_features, max_tile_bytes,
                                      tile_type, level_ratio]() {
                    for (size_t ti = s; ti < e; ++ti) {
                        parent_tile_data& ptd = pass1_results[ti];
                        uint32_t tz = ptd.z, tx = ptd.x, ty = ptd.y;
                        uint32_t seed = tz * 1000000 + tx * 1000 + ty;
                        std::vector<int> indices;
                        uint64_t tile_id = pmtiles::zxy_to_tileid(tz, tx, ty);

                        // Every tile keeps level_ratio fraction of its features
                        uint64_t tile_cap = (uint64_t)(level_ratio * ptd.post_density_count);
                        if (tile_cap < 1 && ptd.post_density_count > 0) tile_cap = 1;

                        if (tile_type == 0x06) {
                            if (ptd.mlt_combined.features.empty()) continue;
                            uint64_t current_max = tile_cap;
                            if (current_max > ptd.mlt_combined.features.size())
                                current_max = ptd.mlt_combined.features.size();
                            while (current_max > 0) {
                                get_subset_indices_mlt(ptd.mlt_combined, current_max, seed, indices);
                                size_t est = indices.size() * (20 + ptd.mlt_combined.col_names.size() * 8);
                                if (est <= (size_t)max_tile_bytes * 4) break;
                                double ratio = (double)((size_t)max_tile_bytes * 4) / est;
                                uint64_t next_max = (uint64_t)(current_max * ratio * 0.95);
                                if (next_max >= current_max) current_max--;
                                else current_max = next_max;
                            }
                            if (current_max > 0 && !indices.empty()) {
                                std::string encoded = encode_mlt_layer(ptd.mlt_combined, indices, tz, tx, ty);
                                std::string compressed = gzip_compress(encoded);
                                std::lock_guard<std::mutex> lock(out_mutex);
                                pwrite(out_fd, compressed.data(), compressed.size(), current_out_offset);
                                final_entries.emplace_back(tile_id, current_out_offset, compressed.size(), 1);
                                tile_feature_counts[tile_id] = indices.size();
                                current_out_offset += compressed.size();
                            }
                        } else {
                            if (ptd.mvt_combined.features.empty()) continue;
                            uint64_t current_max_features = tile_cap;
                            if (current_max_features > ptd.mvt_combined.features.size())
                                current_max_features = ptd.mvt_combined.features.size();
                            while (current_max_features > 0) {
                                get_subset_indices(ptd.mvt_combined, current_max_features, seed, indices);
                                size_t estimated_size = estimate_uncompressed_mvt_size(ptd.mvt_combined, indices, layer_name_global);
                                if (estimated_size <= (size_t)max_tile_bytes * 4) break;
                                double ratio = (double)((size_t)max_tile_bytes * 4) / estimated_size;
                                uint64_t next_max = (uint64_t)(current_max_features * ratio * 0.95);
                                if (next_max >= current_max_features) current_max_features--;
                                else current_max_features = next_max;
                            }
                            if (current_max_features > 0 && !indices.empty()) {
                                std::string encoded = encode_mvt(ptd.mvt_combined, indices, tz, tx, ty, layer_name_global);
                                std::string final_compressed = gzip_compress(encoded);
                                std::lock_guard<std::mutex> lock(out_mutex);
                                pwrite(out_fd, final_compressed.data(), final_compressed.size(), current_out_offset);
                                final_entries.emplace_back(tile_id, current_out_offset, final_compressed.size(), 1);
                                tile_feature_counts[tile_id] = indices.size();
                                current_out_offset += final_compressed.size();
                            }
                        }
                    }
                });
            }
            for (auto& t : threads) t.join();
        }
        // Tally level features
        for (const auto& ptd : pass1_results) {
            uint64_t tid = pmtiles::zxy_to_tileid(ptd.z, ptd.x, ptd.y);
            auto it = tile_feature_counts.find(tid);
            if (it != tile_feature_counts.end()) level_total_features += it->second;
        }
        notice("  z%d: %zu tiles, %zu total features (level_ratio=%.6f, max_raw=%zu, vs z%d: %.4f)",
               z, pass1_results.size(), level_total_features, level_ratio, max_raw_count, z_max,
               zmax_total_features > 0 ? (double)level_total_features / zmax_total_features : 0.0);

        // update level_entries
        level_entries.clear();
        for (const auto& e : final_entries) {
            pmtiles::zxy co = pmtiles::tileid_to_zxy(e.tile_id);
            if (co.z == z) {
                level_entries[e.tile_id] = e;
            }
        }
    }

    notice("Sorting directory...");
    std::sort(final_entries.begin(), final_entries.end(), [](const pmtiles::entryv3& a, const pmtiles::entryv3& b){
        return a.tile_id < b.tile_id;
    });

    notice("Generating PMTiles metadata...");
    std::function<std::string(const std::string&, uint8_t)> compress_fn = [](const std::string& data, uint8_t comp) {
        return gzip_compress(data);
    };

    auto root_leaves = pmtiles::make_root_leaves(compress_fn, pmtiles::COMPRESSION_GZIP, final_entries);
    std::string root_bytes = std::get<0>(root_leaves);
    std::string leaves_bytes = std::get<1>(root_leaves);
    
    // Update minzoom in vector_layers metadata
    if (pmt.jmeta.contains("vector_layers") && pmt.jmeta["vector_layers"].is_array()) {
        for (auto& vlayer : pmt.jmeta["vector_layers"]) {
            if (vlayer.contains("minzoom")) {
                vlayer["minzoom"] = min_zoom;
            }
        }
    }
    if (pmt.jmeta.contains("tilestats") && pmt.jmeta["tilestats"].contains("layers") &&
        pmt.jmeta["tilestats"]["layers"].is_array()) {
        for (auto& tlayer : pmt.jmeta["tilestats"]["layers"]) {
            if (tlayer.contains("minzoom")) {
                tlayer["minzoom"] = min_zoom;
            }
        }
    }
    std::string json_metadata = pmt.jmeta.dump();
    std::string compressed_json = gzip_compress(json_metadata);

    pmtiles::headerv3 header = pmt.hdr;
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
    header.clustered = false; 
    header.internal_compression = pmtiles::COMPRESSION_GZIP;
    header.tile_compression = pmtiles::COMPRESSION_GZIP;
    header.tile_type = tile_type;
    header.min_zoom = min_zoom;
    
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
