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
    uint32_t encode_zigzag32(int32_t v) { return (v << 1) ^ (v >> 31); }
public:
    std::string encode(uint32_t extent, const std::string& layer_name, const fast_mvt& df, const std::vector<int>& indices, double offset_x, double offset_y, double scale_factor) {
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
        append_varint(layer_name.size());
        tmp_ft.append(layer_name);
        append_varint(extent);
        append_varint(1); // num columns
        append_varint(4); // typeCode=GEOMETRY
        append_varint(2); // numStreams

        // geometry type
        tmp_ft.push_back(0x10);
        tmp_ft.push_back(0x02);
        append_varint(indices.size());
        std::string geom_type_data(indices.size(), 0);
        append_varint(geom_type_data.size());
        tmp_ft.append(geom_type_data);

        // vertices
        tmp_ft.push_back(0x13);
        tmp_ft.push_back(0x02);
        append_varint(indices.size() * 2);
        std::string vertex_data;
        auto append_v_varint = [&](std::uint64_t value) {
            do {
                uint8_t byte = value & 0x7F;
                value >>= 7;
                if (value > 0) byte |= 0x80;
                vertex_data.push_back(byte);
            } while (value > 0);
        };
        for(size_t idx=0; idx<indices.size(); ++idx) {
            size_t i = indices[idx];
            int32_t px = std::round((df.features[i].pt.global_x - offset_x) / scale_factor);
            int32_t py = std::round((offset_y - df.features[i].pt.global_y) / scale_factor);
            append_v_varint(encode_zigzag32(px));
            append_v_varint(encode_zigzag32(py));
        }
        append_varint(vertex_data.size());
        tmp_ft.append(vertex_data);

        write_varint(1 + tmp_ft.size());
        out.push_back(1);
        out.append(tmp_ft);
        return out;
    }
};

std::string encode_mlt_raw(const fast_mvt& df, const std::vector<int>& indices, uint8_t z, uint32_t x, uint32_t y, const std::string& layer_name) {
    double scale_factor = pmt_utils::epsg3857_scale_factor(z);
    double offset_x, offset_y;
    pmt_utils::tiletoepsg3857(x, y, z, &offset_x, &offset_y);
    MLTPointEncoder encoder;
    return encoder.encode(4096, layer_name, df, indices, offset_x, offset_y, scale_factor);
}

void decode_mlt_raw(const std::string& buffer, uint8_t z, uint32_t x, uint32_t y, fast_mvt& out) {
    if (buffer.empty()) return;
    double scale_factor = pmt_utils::epsg3857_scale_factor(z);
    double offset_x, offset_y;
    pmt_utils::tiletoepsg3857(x, y, z, &offset_x, &offset_y);
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
        if (tag != 1) { ptr += layer_len - 1; continue; }
        uint64_t name_len = read_varint();
        std::string layer_name(ptr, name_len);
        ptr += name_len;
        uint64_t extent = read_varint();
        uint64_t num_columns = read_varint();
        for (uint64_t col=0; col<num_columns; ++col) {
            uint64_t typeCode = read_varint();
            if (typeCode == 4) { // GEOMETRY
                uint64_t numStreams = read_varint();
                uint8_t s0_phys = *ptr++; (void)s0_phys;
                uint8_t s0_tech = *ptr++; (void)s0_tech;
                uint64_t s0_num_vals = read_varint(); (void)s0_num_vals;
                uint64_t s0_len = read_varint();
                ptr += s0_len;
                uint8_t s1_phys = *ptr++; (void)s1_phys;
                uint8_t s1_tech = *ptr++; (void)s1_tech;
                uint64_t s1_num_vals = read_varint(); (void)s1_num_vals;
                uint64_t s1_len = read_varint();
                const char* v_end = ptr + s1_len;
                while (ptr < v_end) {
                    uint64_t zx = read_varint();
                    uint64_t zy = read_varint();
                    int32_t px = protozero::decode_zigzag32(zx);
                    int32_t py = protozero::decode_zigzag32(zy);
                    fast_feature f(pmt_utils::pmt_pt_t(z, 0, 0));
                    f.pt.global_x = offset_x + scale_factor * px;
                    f.pt.global_y = offset_y - scale_factor * py;
                    out.features.push_back(std::move(f));
                }
            }
        }
    }
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

void get_subset_indices(const fast_mvt& df, uint64_t max_features, uint32_t point_density_threshold, uint8_t z, uint32_t x, uint32_t y, uint32_t seed, std::vector<int>& indices) {
    size_t total_points = df.features.size();
    indices.clear();
    indices.reserve(total_points);
    
    std::mt19937 rng(seed);
    
    if (point_density_threshold > 0) {
        // Density aware grid dropping
        double scale_factor = pmt_utils::epsg3857_scale_factor(z);
        double offset_x, offset_y;
        pmt_utils::tiletoepsg3857(x, y, z, &offset_x, &offset_y);
        
        std::unordered_map<uint64_t, std::vector<int>> cell_map;
        
        for(size_t i=0; i<total_points; ++i) {
            int32_t px = std::round((df.features[i].pt.global_x - offset_x) / scale_factor);
            int32_t py = std::round((offset_y - df.features[i].pt.global_y) / scale_factor);
            
            // Limit to tile bounds ideally
            int32_t cell_x = std::max(0, std::min(255, (int)(px / 16.0))); // Assume 4096 extent -> 256 cells
            int32_t cell_y = std::max(0, std::min(255, (int)(py / 16.0)));
            
            uint64_t cell_id = ((uint64_t)cell_x << 32) | (uint32_t)cell_y;
            cell_map[cell_id].push_back(i);
        }
        
        for (auto& kv : cell_map) {
            if (kv.second.size() <= point_density_threshold) {
                for (int idx : kv.second) indices.push_back(idx);
            } else {
                std::vector<int> subset = kv.second;
                std::shuffle(subset.begin(), subset.end(), rng);
                for (int i=0; i < point_density_threshold; ++i) indices.push_back(subset[i]);
            }
        }
        
        // At this point indices are naturally trimmed based on absolute density per 16x16 visual space
        if (indices.size() <= max_features) {
            std::sort(indices.begin(), indices.end());
            return;
        }
    } else {
        // Uniform uniform fast-path array
        indices.resize(total_points);
        for(size_t i=0; i<total_points; ++i) indices[i] = i;
    }
    
    // Final check for max limit constraints
    if (indices.size() > max_features) {
        std::shuffle(indices.begin(), indices.end(), rng);
        indices.resize(max_features);
    }
    std::sort(indices.begin(), indices.end());
}

// Global state for build process
std::mutex out_mutex;
int out_fd = -1;
uint64_t current_out_offset = 0;
std::vector<pmtiles::entryv3> final_entries;
std::string layer_name_global = "data";

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

void builder_worker(PyramidBuilderQueue& queue, 
                    const std::map<uint64_t, pmtiles::entryv3>& next_level_entries,
                    uint64_t max_features_limit, uint64_t max_bytes_limit, uint32_t point_density_threshold, uint8_t tile_type) {
    
    pmtiles::entry_zxy entry(0,0,0,0,0);
    while (queue.get_tile(entry)) {
        fast_mvt combined_df;
        uint32_t z = entry.z;
        uint32_t x = entry.x;
        uint32_t y = entry.y;
        
        std::map<std::string, uint32_t> global_key_map;
        std::map<std::string, uint32_t> global_val_map;
        
        // Fetch 4 children
        for(int dy=0; dy<2; ++dy) {
            for(int dx=0; dx<2; ++dx) {
                uint32_t cz = z + 1;
                uint32_t cx = 2 * x + dx;
                uint32_t cy = 2 * y + dy;
                uint64_t c_id = pmtiles::zxy_to_tileid(cz, cx, cy);
                
                auto it = next_level_entries.find(c_id);
                if (it != next_level_entries.end()) {
                    std::string compressed(it->second.length, '\0');
                    pread(out_fd, &compressed[0], it->second.length, it->second.offset);
                    
                    std::string uncompressed = gzip_decompress(compressed);
                    fast_mvt child_mvt;
                    if (tile_type == 0x06) decode_mlt_raw(uncompressed, cz, cx, cy, child_mvt);
                    else decode_mvt_raw(uncompressed, cz, cx, cy, child_mvt);
                    
                    std::vector<uint32_t> k_remap(child_mvt.keys.size());
                    for(size_t i=0; i<child_mvt.keys.size(); ++i) {
                        auto gk = global_key_map.find(child_mvt.keys[i]);
                        if (gk == global_key_map.end()) {
                            k_remap[i] = combined_df.keys.size();
                            global_key_map[child_mvt.keys[i]] = k_remap[i];
                            combined_df.keys.push_back(child_mvt.keys[i]);
                        } else {
                            k_remap[i] = gk->second;
                        }
                    }
                    
                    std::vector<uint32_t> v_remap(child_mvt.values.size());
                    for(size_t i=0; i<child_mvt.values.size(); ++i) {
                        auto gv = global_val_map.find(child_mvt.values[i]);
                        if (gv == global_val_map.end()) {
                            v_remap[i] = combined_df.values.size();
                            global_val_map[child_mvt.values[i]] = v_remap[i];
                            combined_df.values.push_back(child_mvt.values[i]);
                        } else {
                            v_remap[i] = gv->second;
                        }
                    }
                    
                    for(auto& f : child_mvt.features) {
                        for(size_t t=0; t<f.tags.size(); t+=2) {
                            f.tags[t] = k_remap[f.tags[t]];
                            f.tags[t+1] = v_remap[f.tags[t+1]];
                        }
                        combined_df.features.push_back(std::move(f));
                    }
                }
            }
        }
        
        if (combined_df.features.empty()) continue;
        
        // iteratively sample down based on estimated uncompressed size
        uint64_t current_max_features = combined_df.features.size();
        if (current_max_features > max_features_limit) {
            current_max_features = max_features_limit;
        }
        
        uint32_t seed = z * 1000000 + x * 1000 + y;
        std::vector<int> indices;
        
        while (current_max_features > 0) {
            get_subset_indices(combined_df, current_max_features, point_density_threshold, z, x, y, seed, indices);
            
            size_t estimated_size = estimate_uncompressed_mvt_size(combined_df, indices, layer_name_global);
            
            if (estimated_size <= max_bytes_limit) {
                break;
            }
            
            double ratio = (double)max_bytes_limit / estimated_size;
            uint64_t next_max = (uint64_t)(current_max_features * ratio * 0.95);
            if (next_max >= current_max_features) {
                current_max_features--;
            } else {
                current_max_features = next_max;
            }
        }
        
        if (current_max_features > 0) {
            std::string encoded;
            if (tile_type == 0x06) encoded = encode_mlt_raw(combined_df, indices, z, x, y, layer_name_global);
            else encoded = encode_mvt(combined_df, indices, z, x, y, layer_name_global);

            std::string final_compressed = gzip_compress(encoded);
            
            std::lock_guard<std::mutex> lock(out_mutex);
            pwrite(out_fd, final_compressed.data(), final_compressed.size(), current_out_offset);
            uint64_t tile_id = pmtiles::zxy_to_tileid(z, x, y);
            final_entries.emplace_back(tile_id, current_out_offset, final_compressed.size(), 1);
            current_out_offset += final_compressed.size();
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
    int32_t point_density_threshold = 0;
    int32_t num_threads = 0;

    paramList pl;
    BEGIN_LONG_PARAMS(longParameters)
    LONG_STRING_PARAM("in", &in_pmtiles, "Input PMTiles file (finest level only)")
    LONG_STRING_PARAM("out", &out_pmtiles, "Output pyramidal PMTiles file")
    LONG_INT_PARAM("min-zoom", &min_zoom, "Minimum zoom level to construct")
    LONG_INT_PARAM("max-tile-bytes", &max_tile_bytes, "Maximum compressed tile bytes")
    LONG_INT_PARAM("max-tile-features", &max_tile_features, "Maximum features per tile")
    LONG_INT_PARAM("point-density-threshold", &point_density_threshold, "Maximum points allowed per simulated 16x16 visual area unit (Tippecanoe density limiter) 0=disable")
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

    notice("Copying z_max tiles...");
    for(size_t i=0; i<pmt.tile_entries.size(); ++i) {
        const auto& entry = pmt.tile_entries[i];
        if (entry.z == z_max) {
            std::string buffer;
            pmt.flex_reader_ptr->read_at(entry.offset, entry.length, buffer);
            
            pwrite(out_fd, buffer.data(), buffer.size(), current_out_offset);
            
            uint64_t tile_id = pmtiles::zxy_to_tileid(entry.z, entry.x, entry.y);
            pmtiles::entryv3 e(tile_id, current_out_offset, buffer.size(), 1);
            final_entries.push_back(e);
            level_entries[tile_id] = e;
            current_out_offset += buffer.size();
        }
    }

    for (int32_t z = z_max - 1; z >= min_zoom; --z) {
        notice("Building zoom level %d...", z);
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
        
        PyramidBuilderQueue queue;
        for (const auto& kv : parent_tiles) {
            queue.add_tile(kv.second);
        }
        queue.finish();
        
        std::vector<std::thread> threads;
        for (int32_t i = 0; i < num_threads; ++i) {
            threads.emplace_back(builder_worker, std::ref(queue), std::cref(level_entries), max_tile_features, max_tile_bytes, point_density_threshold, tile_type);
        }
        for (auto& t : threads) t.join();
        
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
