#include "pmpoint.h"
#include "qgenlib/tsv_reader.h"
#include "qgenlib/qgen_error.h"

#include <vector>
#include <string>
#include <cstring>
#include <climits>
#include <cmath>
#include <cstdio>
#include <cstdint>
#include <map>
#include <set>
#include <unordered_map>
#include <unordered_set>
#include <limits>
#include <algorithm>
#include <zlib.h>

#include "pmt_pts.h"
#include "pmt_utils.h"
#include "polygon.h"
#include "mvt_pts.h"
#include "htslib/hts.h"
#include "ext/nlohmann/json.hpp"

// =====================================================================
// MLT tile decoding helpers (duplicated from cmd_export_pmtiles.cpp)
// =====================================================================

static std::vector<bool> mlt_qt_decode_bool_rle(const uint8_t* data, size_t len, size_t count) {
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

static void decode_mlt_tile_to_df(const std::string& tile_buf, uint8_t zoom,
                                   int64_t tile_x, int64_t tile_y, pt_dataframe& df) {
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

        uint64_t name_len = rv();
        ptr += name_len;
        rv();
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
            if      (base >= 20 && base <= 23) col_types[c] = 2;
            else if (base >= 24 && base <= 27) col_types[c] = 1;
            else                               col_types[c] = 0;
        }

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
            if (phys == 1 && dict == 3) {
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

                if (phys == 0) {
                    present = mlt_qt_decode_bool_rle(sd, bl, num_features);
                } else if (phys == 1) {
                    if (ctype == 2) {
                        const uint8_t* dp = sd;
                        size_t fi = 0;
                        for (uint64_t vi = 0; vi < nv; ++vi) {
                            uint64_t zig=0; int sh=0;
                            while(dp<sd+bl){uint8_t b=*dp++;zig|=(uint64_t)(b&0x7F)<<sh;sh+=7;if(!(b&0x80))break;}
                            int64_t val=(int64_t)((zig>>1)^-(int64_t)(zig&1));
                            while(fi<num_features&&!present[fi])++fi;
                            if(fi<num_features) attr_vals[c][fi++]=std::to_string(val);
                        }
                    } else if (ctype == 1) {
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
                    } else {
                        str_data=sd; str_data_len=bl; (void)nv;
                    }
                } else if (phys == 3) {
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

        for (size_t i = 0; i < num_features; ++i) {
            pmt_utils::pmt_pt_t pt(zoom, feat_gx[i], feat_gy[i]);
            df.points.push_back(pt);
            for (size_t c = 0; c < num_attr; ++c) {
                const std::string& v = attr_vals[c][i];
                df.add_feature(col_metas[c+1].name, v.empty() ? "NA" : v);
            }
        }

        break;
    }
}

// =====================================================================
// Quadtree statistical helpers
// =====================================================================

// Build a bottom-up hierarchical sum array. out[d] is a row-major (1<<d)
// x (1<<d) grid; out[max_depth] = leaf_grid; lower levels are sums of 4
// children.
static void build_hierarchy(const std::vector<uint64_t>& leaf_grid,
                            int max_depth,
                            std::vector<std::vector<uint64_t>>& out)
{
    out.clear();
    out.resize(max_depth + 1);
    out[max_depth] = leaf_grid;
    for (int d = max_depth - 1; d >= 0; --d) {
        int n  = 1 << d;
        int n2 = 1 << (d + 1);
        out[d].assign((size_t)n * n, 0);
        for (int r = 0; r < n; ++r) {
            for (int c = 0; c < n; ++c) {
                uint64_t s = 0;
                for (int dr = 0; dr < 2; ++dr)
                    for (int dc = 0; dc < 2; ++dc)
                        s += out[d + 1][(size_t)(2*r + dr) * n2 + (2*c + dc)];
                out[d][(size_t)r * n + c] = s;
            }
        }
    }
}

// Stats from a 4x2 contingency table: 4 children x {gene, background}.
struct NodeStat {
    bool   valid          = false;
    double chi2           = 0.0;
    int    df             = 3;          // (4-1)*(2-1)
    double cramers_v      = 0.0;        // sqrt(chi2/N) since min(k1,k2)-1 = 1
    double g_stat         = 0.0;        // likelihood-ratio G
    double bic_penalized  = 0.0;        // G - df * log(N)
    double N              = 0.0;
    double gene_count     = 0.0;
    double bg_count       = 0.0;
    double residuals[4]   = {0,0,0,0};  // signed Pearson residuals for the gene column
    double lfc[4]         = {0,0,0,0};  // log2 fold-change of (gene fraction in child) / (gene fraction in parent), pseudocount 0.5
};

static NodeStat compute_node_stat(const double child_gene[4], const double child_bg[4])
{
    NodeStat s;
    double Ng = 0.0, Nb = 0.0;
    for (int i = 0; i < 4; ++i) { Ng += child_gene[i]; Nb += child_bg[i]; }
    double N = Ng + Nb;
    s.N = N; s.gene_count = Ng; s.bg_count = Nb;
    if (N <= 0.0 || Ng <= 0.0 || Nb <= 0.0) return s;

    double gene_frac_parent = Ng / N;
    double chi2 = 0.0, g_stat = 0.0;
    for (int i = 0; i < 4; ++i) {
        double Ri = child_gene[i] + child_bg[i];
        if (Ri <= 0.0) continue;

        double e_g = Ng * Ri / N;
        double e_b = Nb * Ri / N;

        if (e_g > 0.0) {
            double d = child_gene[i] - e_g;
            chi2 += d * d / e_g;
            s.residuals[i] = d / std::sqrt(e_g);
            if (child_gene[i] > 0.0)
                g_stat += 2.0 * child_gene[i] * std::log(child_gene[i] / e_g);
        }
        if (e_b > 0.0) {
            double d = child_bg[i] - e_b;
            chi2 += d * d / e_b;
            if (child_bg[i] > 0.0)
                g_stat += 2.0 * child_bg[i] * std::log(child_bg[i] / e_b);
        }

        double gene_frac_child = (child_gene[i] + 0.5) / (Ri + 1.0);
        s.lfc[i] = std::log2(gene_frac_child / (gene_frac_parent <= 0.0 ? 1e-12 : gene_frac_parent));
    }
    s.chi2          = chi2;
    s.g_stat        = g_stat;
    s.df            = 3;
    s.cramers_v     = std::sqrt(chi2 / N);
    s.bic_penalized = g_stat - (double)s.df * std::log(N > 1.0 ? N : 2.0);
    s.valid         = true;
    return s;
}

struct GeneSummary {
    std::vector<double> energy_per_depth;
    int    n_tested_nodes      = 0;
    int    max_depth_used      = 0;
    double total_energy        = 0.0;
    double characteristic_scale= 0.0;
    double scale_entropy       = 0.0;
    // Best-by-BIC node (used for the visualization heading + ranking).
    bool   best_valid          = false;
    double best_bic            = -std::numeric_limits<double>::infinity();
    int    best_depth          = -1;
    int    best_r              = -1;
    int    best_c              = -1;
    double best_chi2           = 0.0;
    double best_cramers_v      = 0.0;
};

static void accumulate_node_to_summary(const NodeStat& s, int depth, int r, int c, GeneSummary& sum)
{
    if (!s.valid) return;
    if ((int)sum.energy_per_depth.size() <= depth) sum.energy_per_depth.resize(depth + 1, 0.0);
    sum.energy_per_depth[depth] += s.chi2;
    sum.n_tested_nodes++;
    if (depth > sum.max_depth_used) sum.max_depth_used = depth;
    if (s.bic_penalized > sum.best_bic) {
        sum.best_bic       = s.bic_penalized;
        sum.best_depth     = depth;
        sum.best_r         = r;
        sum.best_c         = c;
        sum.best_chi2      = s.chi2;
        sum.best_cramers_v = s.cramers_v;
        sum.best_valid     = true;
    }
}

static void finalize_summary(GeneSummary& sum)
{
    double tot = 0.0;
    for (double e : sum.energy_per_depth) tot += e;
    sum.total_energy = tot;
    if (tot > 0.0) {
        double weighted = 0.0, h = 0.0;
        for (size_t d = 0; d < sum.energy_per_depth.size(); ++d) {
            double w = sum.energy_per_depth[d] / tot;
            weighted += (double)d * w;
            if (w > 0.0) h -= w * std::log(w);
        }
        sum.characteristic_scale = weighted;
        sum.scale_entropy        = h;
    }
}

// Recursively emit a JSON tree for a single gene + accumulate summary stats.
// (depth, row, col) are the node coordinates in the local 2^depth grid;
// (x0, y0, x1, y1) is its EPSG:3857 bounding box.
static nlohmann::json build_tree_json(int d, int r, int c,
                                      int max_depth,
                                      double x0, double y0, double x1, double y1,
                                      const std::vector<std::vector<uint64_t>>& gene_h,
                                      const std::vector<std::vector<uint64_t>>& total_h,
                                      double min_node_count,
                                      GeneSummary& summary)
{
    int n = 1 << d;
    double gene_c  = (double)gene_h[d][(size_t)r * n + c];
    double total_c = (double)total_h[d][(size_t)r * n + c];
    double bg_c    = total_c - gene_c;

    nlohmann::json node;
    node["depth"]       = d;
    node["row"]         = r;
    node["col"]         = c;
    node["bbox"]        = {x0, y0, x1, y1};
    node["gene_count"]  = gene_c;
    node["total_count"] = total_c;
    node["bg_count"]    = bg_c;

    if (d >= max_depth) return node;

    int n2 = 1 << (d + 1);
    double mx = 0.5 * (x0 + x1);
    double my = 0.5 * (y0 + y1);

    double child_gene[4], child_bg[4], child_total[4];
    for (int i = 0; i < 4; ++i) {
        int dr = i / 2, dc = i % 2;
        int rr = 2*r + dr, cc = 2*c + dc;
        double g = (double)gene_h[d + 1][(size_t)rr * n2 + cc];
        double t = (double)total_h[d + 1][(size_t)rr * n2 + cc];
        child_gene[i]  = g;
        child_total[i] = t;
        child_bg[i]    = t - g;
    }

    if (total_c >= min_node_count) {
        NodeStat s = compute_node_stat(child_gene, child_bg);
        if (s.valid) {
            nlohmann::json stat;
            stat["chi2"]              = s.chi2;
            stat["df"]                = s.df;
            stat["g_stat"]            = s.g_stat;
            stat["cramers_v"]         = s.cramers_v;
            stat["bic_penalized"]     = s.bic_penalized;
            stat["pearson_residuals"] = {s.residuals[0], s.residuals[1], s.residuals[2], s.residuals[3]};
            stat["log2_fold_change"]  = {s.lfc[0], s.lfc[1], s.lfc[2], s.lfc[3]};
            node["test"] = stat;
            accumulate_node_to_summary(s, d, r, c, summary);
        }
    }

    bool any_child_testable = false;
    for (int i = 0; i < 4; ++i)
        if (child_total[i] >= min_node_count) { any_child_testable = true; break; }
    if (d + 1 <= max_depth && any_child_testable) {
        nlohmann::json children = nlohmann::json::array();
        for (int i = 0; i < 4; ++i) {
            int dr = i / 2, dc = i % 2;
            double cx0 = (dc == 0) ? x0 : mx;
            double cx1 = (dc == 0) ? mx : x1;
            double cy0 = (dr == 0) ? y0 : my;
            double cy1 = (dr == 0) ? my : y1;
            children.push_back(build_tree_json(d + 1, 2*r + dr, 2*c + dc,
                                               max_depth, cx0, cy0, cx1, cy1,
                                               gene_h, total_h, min_node_count, summary));
        }
        node["children"] = children;
    }
    return node;
}

static int find_feature_col(const pt_dataframe& df, const std::string& name) {
    for (int i = 0; i < (int)df.feature_names.size(); ++i)
        if (df.feature_names[i] == name) return i;
    return -1;
}

// =====================================================================
// Minimal PNG writer (RGBA, single IDAT, zlib deflate)
// =====================================================================
namespace mini_png {

static void put_be32(uint8_t* b, uint32_t v) {
    b[0] = (v >> 24) & 0xFF; b[1] = (v >> 16) & 0xFF;
    b[2] = (v >>  8) & 0xFF; b[3] = (v      ) & 0xFF;
}

static void write_chunk(FILE* fp, const char tag[4], const uint8_t* data, uint32_t len) {
    uint8_t lenb[4]; put_be32(lenb, len);
    fwrite(lenb, 1, 4, fp);
    fwrite(tag,  1, 4, fp);
    if (len > 0) fwrite(data, 1, len, fp);
    uint32_t crc = crc32(0L, (const Bytef*)tag, 4);
    if (len > 0) crc = crc32(crc, (const Bytef*)data, len);
    uint8_t crcb[4]; put_be32(crcb, crc);
    fwrite(crcb, 1, 4, fp);
}

static bool write_rgba(const char* path, const uint8_t* px, int w, int h) {
    FILE* fp = fopen(path, "wb");
    if (!fp) return false;

    const uint8_t sig[8] = {137,80,78,71,13,10,26,10};
    fwrite(sig, 1, 8, fp);

    uint8_t ihdr[13];
    put_be32(&ihdr[0], (uint32_t)w);
    put_be32(&ihdr[4], (uint32_t)h);
    ihdr[8]  = 8;   // bit depth
    ihdr[9]  = 6;   // RGBA
    ihdr[10] = 0;
    ihdr[11] = 0;
    ihdr[12] = 0;
    write_chunk(fp, "IHDR", ihdr, 13);

    size_t row_bytes = (size_t)w * 4;
    std::vector<uint8_t> raw((row_bytes + 1) * (size_t)h);
    for (int y = 0; y < h; ++y) {
        raw[(row_bytes + 1) * (size_t)y] = 0; // None filter
        memcpy(&raw[(row_bytes + 1) * (size_t)y + 1], &px[row_bytes * y], row_bytes);
    }

    uLongf comp_len = compressBound((uLong)raw.size());
    std::vector<uint8_t> comp(comp_len);
    if (compress2(comp.data(), &comp_len, raw.data(), (uLong)raw.size(), 6) != Z_OK) {
        fclose(fp);
        return false;
    }
    write_chunk(fp, "IDAT", comp.data(), (uint32_t)comp_len);
    write_chunk(fp, "IEND", nullptr, 0);
    fclose(fp);
    return true;
}

} // namespace mini_png

// =====================================================================
// Visualization: adaptive-partition heatmap for one gene
// =====================================================================

// Diverging blue->white->red color from log2 fold-change.
static void log2fc_to_rgb(double log2fc, double max_log2fc, uint8_t out[3]) {
    if (!std::isfinite(log2fc)) { out[0]=200; out[1]=200; out[2]=200; return; }
    double t = log2fc / max_log2fc;
    if (t > 1.0)  t =  1.0;
    if (t < -1.0) t = -1.0;
    if (t >= 0) {
        out[0] = 255;
        out[1] = (uint8_t)(255.0 * (1.0 - t));
        out[2] = (uint8_t)(255.0 * (1.0 - t));
    } else {
        out[0] = (uint8_t)(255.0 * (1.0 + t));
        out[1] = (uint8_t)(255.0 * (1.0 + t));
        out[2] = 255;
    }
}

struct VizCtx {
    int    img_w = 0, img_h = 0;
    std::vector<uint8_t> img;     // RGBA
    int    max_depth = 0;
    double global_gene_frac = 0.0;
    double max_log2fc = 4.0;
    double min_opacity = 0.05;
    double min_node_count = 20.0;
    double max_density = 1.0;     // count per leaf cell, used to scale opacity
};

static void fill_rect(VizCtx& ctx, int x0, int y0, int x1, int y1,
                      uint8_t rr, uint8_t gg, uint8_t bb, uint8_t aa, bool border) {
    if (x0 < 0) x0 = 0; if (y0 < 0) y0 = 0;
    if (x1 > ctx.img_w) x1 = ctx.img_w;
    if (y1 > ctx.img_h) y1 = ctx.img_h;
    for (int y = y0; y < y1; ++y) {
        for (int x = x0; x < x1; ++x) {
            size_t i = ((size_t)y * ctx.img_w + x) * 4;
            ctx.img[i]   = rr;
            ctx.img[i+1] = gg;
            ctx.img[i+2] = bb;
            ctx.img[i+3] = aa;
        }
    }
    if (border && x1 - x0 > 2 && y1 - y0 > 2) {
        for (int x = x0; x < x1; ++x) {
            size_t i0 = ((size_t)y0     * ctx.img_w + x) * 4;
            size_t i1 = ((size_t)(y1-1) * ctx.img_w + x) * 4;
            ctx.img[i0] = ctx.img[i1] = 60;
            ctx.img[i0+1] = ctx.img[i1+1] = 60;
            ctx.img[i0+2] = ctx.img[i1+2] = 60;
        }
        for (int y = y0; y < y1; ++y) {
            size_t il = ((size_t)y * ctx.img_w + x0    ) * 4;
            size_t ir = ((size_t)y * ctx.img_w + (x1-1)) * 4;
            ctx.img[il] = ctx.img[ir] = 60;
            ctx.img[il+1] = ctx.img[ir+1] = 60;
            ctx.img[il+2] = ctx.img[ir+2] = 60;
        }
    }
}

static void render_leaf(VizCtx& ctx, int d, int r, int c,
                        double gene_c, double total_c) {
    int N = 1 << d;
    int x0 = (int)((double)c       * ctx.img_w / N);
    int x1 = (int)((double)(c + 1) * ctx.img_w / N);
    int y0 = (int)((double)r       * ctx.img_h / N);
    int y1 = (int)((double)(r + 1) * ctx.img_h / N);
    if (total_c <= 0.0) {
        fill_rect(ctx, x0, y0, x1, y1, 240, 240, 240, 255, true);
        return;
    }
    double frac    = gene_c / total_c;
    double log2fc  = std::log2((frac + 1e-9) / (ctx.global_gene_frac + 1e-9));
    uint8_t rgb[3]; log2fc_to_rgb(log2fc, ctx.max_log2fc, rgb);

    // Density-based opacity so cells of different sizes are comparable.
    double leaf_area = std::pow(4.0, ctx.max_depth - d);
    double density  = total_c / leaf_area;
    double op       = density / (ctx.max_density > 0 ? ctx.max_density : 1.0);
    if (op > 1.0) op = 1.0;
    if (op < ctx.min_opacity) op = ctx.min_opacity;

    // Alpha-blend onto white background so output is opaque RGBA=255.
    uint8_t br = (uint8_t)std::round(rgb[0] * op + 255.0 * (1.0 - op));
    uint8_t bg = (uint8_t)std::round(rgb[1] * op + 255.0 * (1.0 - op));
    uint8_t bb = (uint8_t)std::round(rgb[2] * op + 255.0 * (1.0 - op));
    fill_rect(ctx, x0, y0, x1, y1, br, bg, bb, 255, true);
}

// Walk the tree top-down, splitting only where the local 4-way test passes a
// BIC test. Render the leaves of that pruned tree.
static void render_node(VizCtx& ctx, int d, int r, int c,
                        const std::vector<std::vector<uint64_t>>& gene_h,
                        const std::vector<std::vector<uint64_t>>& total_h)
{
    int n = 1 << d;
    double gene_c  = (double)gene_h[d][(size_t)r * n + c];
    double total_c = (double)total_h[d][(size_t)r * n + c];

    bool split = false;
    if (d < ctx.max_depth && total_c >= ctx.min_node_count) {
        int n2 = 1 << (d + 1);
        double cg[4], cb[4];
        for (int i = 0; i < 4; ++i) {
            int dr = i / 2, dc = i % 2;
            double g = (double)gene_h[d + 1][(size_t)(2*r + dr) * n2 + (2*c + dc)];
            double t = (double)total_h[d + 1][(size_t)(2*r + dr) * n2 + (2*c + dc)];
            cg[i] = g;
            cb[i] = t - g;
        }
        NodeStat s = compute_node_stat(cg, cb);
        if (s.valid && s.bic_penalized > 0.0) split = true;
    }
    if (split) {
        for (int i = 0; i < 4; ++i) {
            int dr = i / 2, dc = i % 2;
            render_node(ctx, d + 1, 2*r + dr, 2*c + dc, gene_h, total_h);
        }
    } else {
        render_leaf(ctx, d, r, c, gene_c, total_c);
    }
}

// Pre-compute the maximum per-area density that would be drawn for this gene.
// Mirrors the descent logic in render_node so opacity scales sensibly.
static double compute_max_density(int d, int r, int c,
                                  const std::vector<std::vector<uint64_t>>& gene_h,
                                  const std::vector<std::vector<uint64_t>>& total_h,
                                  int max_depth, double min_node_count)
{
    int n = 1 << d;
    double total_c = (double)total_h[d][(size_t)r * n + c];

    bool split = false;
    if (d < max_depth && total_c >= min_node_count) {
        int n2 = 1 << (d + 1);
        double cg[4], cb[4];
        for (int i = 0; i < 4; ++i) {
            int dr = i / 2, dc = i % 2;
            double g = (double)gene_h[d + 1][(size_t)(2*r + dr) * n2 + (2*c + dc)];
            double t = (double)total_h[d + 1][(size_t)(2*r + dr) * n2 + (2*c + dc)];
            cg[i] = g;
            cb[i] = t - g;
        }
        NodeStat s = compute_node_stat(cg, cb);
        if (s.valid && s.bic_penalized > 0.0) split = true;
    }
    if (split) {
        double best = 0.0;
        for (int i = 0; i < 4; ++i) {
            int dr = i / 2, dc = i % 2;
            double m = compute_max_density(d + 1, 2*r + dr, 2*c + dc,
                                            gene_h, total_h, max_depth, min_node_count);
            if (m > best) best = m;
        }
        return best;
    }
    double leaf_area = std::pow(4.0, max_depth - d);
    return total_c / leaf_area;
}

static bool render_gene_png(const std::string& outpng,
                            const std::vector<std::vector<uint64_t>>& gene_h,
                            const std::vector<std::vector<uint64_t>>& total_h,
                            int max_depth, int img_w, int img_h,
                            double min_node_count, double max_log2fc, double min_opacity)
{
    double gene_total  = (double)gene_h[0][0];
    double total_total = (double)total_h[0][0];
    if (total_total <= 0.0) {
        notice("Skipping %s: empty root", outpng.c_str());
        return false;
    }
    VizCtx ctx;
    ctx.img_w           = img_w;
    ctx.img_h           = img_h;
    ctx.max_depth       = max_depth;
    ctx.global_gene_frac= gene_total / total_total;
    ctx.max_log2fc      = max_log2fc;
    ctx.min_opacity     = min_opacity;
    ctx.min_node_count  = min_node_count;
    ctx.img.assign((size_t)img_w * img_h * 4, 255);

    ctx.max_density = compute_max_density(0, 0, 0, gene_h, total_h, max_depth, min_node_count);
    if (ctx.max_density <= 0.0) ctx.max_density = 1.0;

    render_node(ctx, 0, 0, 0, gene_h, total_h);
    return mini_png::write_rgba(outpng.c_str(), ctx.img.data(), img_w, img_h);
}

// =====================================================================
// Small utilities
// =====================================================================

static std::vector<std::string> parse_csv(const std::string& s) {
    std::vector<std::string> out;
    std::string cur;
    for (char ch : s) {
        if (ch == ',') { if (!cur.empty()) out.push_back(cur); cur.clear(); }
        else cur.push_back(ch);
    }
    if (!cur.empty()) out.push_back(cur);
    return out;
}

// Filename-safe version of a gene name for PNG output paths.
static std::string sanitize_filename(const std::string& s) {
    std::string out; out.reserve(s.size());
    for (char ch : s) {
        if ((ch >= 'A' && ch <= 'Z') ||
            (ch >= 'a' && ch <= 'z') ||
            (ch >= '0' && ch <= '9') ||
            ch == '-' || ch == '_' || ch == '.')
            out.push_back(ch);
        else
            out.push_back('_');
    }
    if (out.empty()) out = "gene";
    return out;
}

// =====================================================================
// Command: quadtree-test
// =====================================================================

int32_t cmd_quadtree_test_pmtiles(int32_t argc, char **argv)
{
    std::string pmtilesf;

    // Column names
    std::string gene_col("gene");
    std::string count_col("count");

    // Zooms (-1 means auto / use PMTiles default).
    int32_t min_zoom  = -1;
    int32_t max_zoom  = -1;
    int32_t data_zoom = -1;

    // Outputs
    std::string out_jsonf;       // --out-tree (required)
    std::string out_summaryf;    // --out-summary (required by spec)
    std::string out_viz_prefix;  // --out-viz-prefix (optional)
    std::string viz_genes_str;   // --viz-genes (required if out_viz_prefix is set)

    // Optional tuning
    std::string gene_listf;
    double  min_node_count = 20.0;
    int32_t viz_width      = 1024;
    int32_t viz_height     = 1024;
    double  max_log2fc     = 4.0;
    double  min_opacity    = 0.05;

    paramList pl;
    BEGIN_LONG_PARAMS(longParameters)
    LONG_PARAM_GROUP("Input options", NULL)
    LONG_STRING_PARAM("in",         &pmtilesf,   "Input PMTiles file (with per-point gene/count attributes)")
    LONG_STRING_PARAM("gene-col",   &gene_col,   "Column name for gene/feature identity (default: gene)")
    LONG_STRING_PARAM("count-col",  &count_col,  "Column name for transcript counts; empty -> each point counts as 1 (default: count)")
    LONG_STRING_PARAM("gene-list",  &gene_listf, "Optional file with one gene per line; restricts the test (and outputs) to listed genes")

    LONG_PARAM_GROUP("Output options", NULL)
    LONG_STRING_PARAM("out-tree",       &out_jsonf,       "Output JSON file with the per-gene quadtrees and node statistics")
    LONG_STRING_PARAM("out-summary",    &out_summaryf,    "Output TSV file with per-gene summary statistics")
    LONG_STRING_PARAM("out-viz-prefix", &out_viz_prefix,  "Output prefix for PNG visualizations (requires --viz-genes)")
    LONG_STRING_PARAM("viz-genes",      &viz_genes_str,   "Comma-separated list of genes to visualize as PNG heatmaps")

    LONG_PARAM_GROUP("Quadtree options", NULL)
    LONG_INT_PARAM   ("min-zoom",       &min_zoom,       "Minimum zoom level (quadtree root). Default: auto = deepest zoom whose tile encloses all data")
    LONG_INT_PARAM   ("max-zoom",       &max_zoom,       "Maximum zoom level (quadtree leaves). May exceed the PMTiles max zoom. Default: PMTiles max zoom")
    LONG_INT_PARAM   ("data-zoom",      &data_zoom,      "Zoom level to read points from. Default: PMTiles max zoom")
    LONG_DOUBLE_PARAM("min-node-count", &min_node_count, "Minimum total count required to test or split a node (default: 20)")

    LONG_PARAM_GROUP("Visualization options", NULL)
    LONG_INT_PARAM   ("viz-width",   &viz_width,   "Width of PNG output in pixels (default: 1024)")
    LONG_INT_PARAM   ("viz-height",  &viz_height,  "Height of PNG output in pixels (default: 1024)")
    LONG_DOUBLE_PARAM("max-log2fc",  &max_log2fc,  "Color scale saturation in absolute log2 fold-change (default: 4)")
    LONG_DOUBLE_PARAM("min-opacity", &min_opacity, "Minimum cell opacity so sparse regions stay visible (default: 0.05)")
    END_LONG_PARAMS();

    pl.Add(new longParams("Available Options", longParameters));
    pl.Read(argc, argv);
    pl.Status();

    notice("Analysis started");

    if (pmtilesf.empty())     error("Missing required option --in");
    if (out_jsonf.empty())    error("Missing required option --out-tree");
    if (out_summaryf.empty()) error("Missing required option --out-summary");
    if (!out_viz_prefix.empty() && viz_genes_str.empty())
        error("--out-viz-prefix was set but --viz-genes was not; specify a comma-separated list of genes to visualize");
    if (viz_width  < 16 || viz_height < 16)
        error("--viz-width and --viz-height must each be >= 16");

    // Parse optional gene-list filter
    std::set<std::string> gene_filter;
    const bool restrict_genes = !gene_listf.empty();
    if (restrict_genes) {
        tsv_reader rd(gene_listf.c_str());
        while (rd.read_line() > 0) {
            if (rd.nfields > 0) {
                const char* g = rd.str_field_at(0);
                if (g && g[0] != '\0' && g[0] != '#') gene_filter.insert(std::string(g));
            }
        }
        notice("Loaded %zu gene(s) from %s", gene_filter.size(), gene_listf.c_str());
        if (gene_filter.empty()) error("--gene-list %s contained no usable gene names", gene_listf.c_str());
    }

    // Parse --viz-genes (don't require them to be in --gene-list; we'll add them
    // to the test set so trees are computed for them too).
    std::vector<std::string> viz_genes = parse_csv(viz_genes_str);
    std::unordered_set<std::string> viz_gene_set(viz_genes.begin(), viz_genes.end());
    if (restrict_genes) for (auto& g : viz_genes) gene_filter.insert(g);

    // -----------------------------------------------------------------
    // Open PMTiles and resolve zoom defaults.
    // -----------------------------------------------------------------
    pmt_pts pmt(pmtilesf.c_str());
    notice("Reading header and tile entries...");
    if (!pmt.read_header_meta_entries())
        error("This pmtiles file is malformed or incompatible with pmpoints (needs MVT or MLT points)");

    int pmt_min_z = pmt.hdr.min_zoom;
    int pmt_max_z = pmt.hdr.max_zoom;
    if (data_zoom < 0) data_zoom = pmt_max_z;
    if (max_zoom  < 0) max_zoom  = pmt_max_z;

    if (data_zoom < pmt_min_z || data_zoom > pmt_max_z)
        error("--data-zoom %d is outside the PMTiles range [%d, %d]", data_zoom, pmt_min_z, pmt_max_z);
    if (max_zoom < 1) error("--max-zoom must be >= 1");

    // Collect tile entries at data_zoom and derive the smallest enclosing tile.
    std::vector<int32_t> data_tile_idxs;
    int64_t any_tx = -1, any_ty = -1;
    int64_t tx_min = LLONG_MAX, tx_max = LLONG_MIN, ty_min = LLONG_MAX, ty_max = LLONG_MIN;
    for (int32_t i = 0; i < (int32_t)pmt.tile_entries.size(); ++i) {
        pmtiles::entry_zxy& e = pmt.tile_entries[i];
        if (e.z != data_zoom) continue;
        data_tile_idxs.push_back(i);
        if (any_tx < 0) { any_tx = e.x; any_ty = e.y; }
        if ((int64_t)e.x < tx_min) tx_min = e.x;
        if ((int64_t)e.x > tx_max) tx_max = e.x;
        if ((int64_t)e.y < ty_min) ty_min = e.y;
        if ((int64_t)e.y > ty_max) ty_max = e.y;
    }
    if (data_tile_idxs.empty()) error("No tile entries found at zoom %d", data_zoom);

    int auto_min_zoom = 0;
    for (int z = data_zoom; z >= 0; --z) {
        int shift = data_zoom - z;
        int64_t rx = tx_min >> shift;
        int64_t ry = ty_min >> shift;
        if ((tx_max >> shift) == rx && (ty_max >> shift) == ry) { auto_min_zoom = z; break; }
    }
    if (min_zoom < 0) {
        min_zoom = auto_min_zoom;
        notice("Auto-detected --min-zoom = %d (deepest zoom whose tile encloses all data)", min_zoom);
    } else if (min_zoom > auto_min_zoom) {
        error("--min-zoom %d is too deep; data spans multiple tiles at that zoom. Try <= %d", min_zoom, auto_min_zoom);
    }
    if (max_zoom < min_zoom)
        error("--max-zoom (%d) must be >= --min-zoom (%d)", max_zoom, min_zoom);

    int depth = max_zoom - min_zoom;
    if (depth > 14) error("max_zoom - min_zoom = %d would create a %d x %d leaf grid, which is too large for this draft", depth, 1 << depth, 1 << depth);

    int64_t root_tx = tx_min >> (data_zoom - min_zoom);
    int64_t root_ty = ty_min >> (data_zoom - min_zoom);

    double r_x0, r_y0, r_x1, r_y1;
    pmt_utils::tiletoepsg3857(root_tx,     root_ty,     min_zoom, &r_x0, &r_y1); // top-left  -> (min_x, max_y)
    pmt_utils::tiletoepsg3857(root_tx + 1, root_ty + 1, min_zoom, &r_x1, &r_y0); // bot-right -> (max_x, min_y)

    notice("Zooms: data=%d, min=%d, max=%d  (tree depth=%d, leaf grid=%d x %d)",
           data_zoom, min_zoom, max_zoom, depth, 1 << depth, 1 << depth);
    notice("Root tile (z=%d, x=%lld, y=%lld) bbox EPSG:3857 = x[%.3f, %.3f] y[%.3f, %.3f]",
           min_zoom, (long long)root_tx, (long long)root_ty, r_x0, r_x1, r_y0, r_y1);

    // -----------------------------------------------------------------
    // Pass: iterate tiles at data_zoom, bin each point to a leaf cell.
    // -----------------------------------------------------------------
    const int N = 1 << depth;
    const size_t n_leaves = (size_t)N * (size_t)N;
    const int64_t leaf_x_base = root_tx << depth;  // root_tx * N
    const int64_t leaf_y_base = root_ty << depth;
    const bool no_count_col = count_col.empty();

    std::vector<uint64_t> bg_leaf(n_leaves, 0);
    std::unordered_map<std::string, std::unordered_map<uint64_t, uint64_t>> gene2leaf;

    mvt_pts mvt;                 // for MVT decode (no filter required)
    pt_dataframe df;
    std::string tile_buffer;
    uint64_t n_points_seen = 0, n_points_binned = 0;
    int n_done = 0;
    const int total_sel = (int)data_tile_idxs.size();

    for (int32_t ix : data_tile_idxs) {
        pmtiles::entry_zxy& entry = pmt.tile_entries[ix];

        pmt.fetch_tile_to_buffer(entry.z, entry.x, entry.y, tile_buffer);
        df.clear_values();
        if (pmt.hdr.tile_type == 0x06) {
            decode_mlt_tile_to_df(tile_buffer, entry.z, entry.x, entry.y, df);
        } else {
            mvt.decode_points_df(tile_buffer, entry.z, entry.x, entry.y, df);
        }

        int gene_col_idx  = find_feature_col(df, gene_col);
        int count_col_idx = no_count_col ? -1 : find_feature_col(df, count_col);
        if (gene_col_idx < 0)
            error("Gene column '%s' not found in tile %d/%d/%d (available columns: %zu)",
                  gene_col.c_str(), entry.z, entry.x, entry.y, df.feature_names.size());

        for (size_t i = 0; i < df.points.size(); ++i) {
            ++n_points_seen;
            double gx = df.points[i].global_x;
            double gy = df.points[i].global_y;

            int64_t tx, ty;
            pmt_utils::epsg3857totile(gx, gy, (uint8_t)max_zoom, &tx, &ty);
            int64_t lc = tx - leaf_x_base;
            int64_t lr = ty - leaf_y_base;
            if (lc < 0 || lc >= N || lr < 0 || lr >= N) continue;  // outside root tile
            uint64_t leaf_idx = (uint64_t)lr * (uint64_t)N + (uint64_t)lc;

            const std::string& gene_val = df.feature_matrix[gene_col_idx][i];
            uint64_t cnt = 1;
            if (!no_count_col) {
                const std::string& s = df.feature_matrix[count_col_idx][i];
                if (!s.empty() && s != "NA") {
                    try { cnt = (uint64_t)std::stoll(s); }
                    catch (...) { cnt = 1; }
                }
            }
            if (cnt == 0) continue;

            bg_leaf[leaf_idx] += cnt;
            if (restrict_genes && gene_filter.find(gene_val) == gene_filter.end()) continue;
            gene2leaf[gene_val][leaf_idx] += cnt;
            ++n_points_binned;
        }

        ++n_done;
        if (n_done % 50 == 0 || n_done == total_sel)
            notice("Binned %d / %d tiles, %llu pts seen, %llu pts in test set (%zu gene(s) so far)",
                   n_done, total_sel,
                   (unsigned long long)n_points_seen, (unsigned long long)n_points_binned,
                   gene2leaf.size());
    }

    if (gene2leaf.empty()) error("No genes matched the criteria; nothing to test");

    // Warn about genes requested for visualization but not seen in the data.
    if (!viz_genes.empty()) {
        std::vector<std::string> missing;
        for (auto& g : viz_genes) if (gene2leaf.find(g) == gene2leaf.end()) missing.push_back(g);
        for (auto& g : missing)
            notice("WARNING: --viz-genes '%s' has no counts in the data; PNG will be skipped", g.c_str());
    }

    // -----------------------------------------------------------------
    // Background hierarchy (shared across all per-gene tests).
    // -----------------------------------------------------------------
    notice("Building background hierarchy...");
    std::vector<std::vector<uint64_t>> bg_h;
    build_hierarchy(bg_leaf, depth, bg_h);

    // -----------------------------------------------------------------
    // Per-gene tree + summary + (optional) PNG
    // -----------------------------------------------------------------
    auto open_w = [](const std::string& path) -> htsFile* {
        if (path.size() >= 3 && path.compare(path.size() - 3, 3, ".gz") == 0)
            return hts_open(path.c_str(), "wz");
        return hts_open(path.c_str(), "w");
    };

    htsFile* json_wh = open_w(out_jsonf);
    if (!json_wh) error("Failed to open --out-tree %s", out_jsonf.c_str());

    htsFile* summary_wh = open_w(out_summaryf);
    if (!summary_wh) error("Failed to open --out-summary %s", out_summaryf.c_str());
    hprintf(summary_wh,
            "gene\ttotal_count\tn_tested_nodes\tmax_depth_used\ttotal_energy\t"
            "characteristic_scale\tscale_entropy\tbest_bic\tbest_depth\tbest_zoom\t"
            "best_chi2\tbest_cramers_v\tenergy_per_depth\n");

    // Stream JSON manually so we don't hold every per-gene tree in memory.
    hprintf(json_wh, "{\n");
    hprintf(json_wh, "  \"metadata\": {\n");
    hprintf(json_wh, "    \"input\": \"%s\",\n", pmtilesf.c_str());
    hprintf(json_wh, "    \"data_zoom\": %d,\n", data_zoom);
    hprintf(json_wh, "    \"min_zoom\": %d,\n",  min_zoom);
    hprintf(json_wh, "    \"max_zoom\": %d,\n",  max_zoom);
    hprintf(json_wh, "    \"depth\": %d,\n",     depth);
    hprintf(json_wh, "    \"root_tile\": [%d, %lld, %lld],\n", min_zoom,
            (long long)root_tx, (long long)root_ty);
    hprintf(json_wh, "    \"bbox\": [%.6f, %.6f, %.6f, %.6f],\n", r_x0, r_y0, r_x1, r_y1);
    hprintf(json_wh, "    \"min_node_count\": %.6g,\n", min_node_count);
    hprintf(json_wh, "    \"n_genes\": %zu,\n", gene2leaf.size());
    hprintf(json_wh, "    \"gene_col\": \"%s\",\n", gene_col.c_str());
    hprintf(json_wh, "    \"count_col\": \"%s\"\n", count_col.c_str());
    hprintf(json_wh, "  },\n");

    // Background tree (gene == total). Useful as a denominator reference and a
    // sanity check that the binning is correct.
    {
        GeneSummary dummy;
        nlohmann::json bg_tree = build_tree_json(0, 0, 0, depth,
                                                  r_x0, r_y0, r_x1, r_y1,
                                                  bg_h, bg_h, min_node_count, dummy);
        std::string bg_dump = bg_tree.dump();
        hprintf(json_wh, "  \"background\": %s,\n", bg_dump.c_str());
    }

    hprintf(json_wh, "  \"genes\": {\n");

    size_t gi = 0;
    const size_t total_genes = gene2leaf.size();
    int n_viz_written = 0;
    for (auto& kv : gene2leaf) {
        const std::string& gene_name = kv.first;
        const auto& sparse = kv.second;

        std::vector<uint64_t> gene_leaf(n_leaves, 0);
        uint64_t gene_total = 0;
        for (auto& sk : sparse) { gene_leaf[sk.first] = sk.second; gene_total += sk.second; }

        std::vector<std::vector<uint64_t>> gene_h;
        build_hierarchy(gene_leaf, depth, gene_h);

        GeneSummary summary;
        nlohmann::json tree = build_tree_json(0, 0, 0, depth,
                                              r_x0, r_y0, r_x1, r_y1,
                                              gene_h, bg_h, min_node_count, summary);
        finalize_summary(summary);

        // JSON entry
        std::string esc_name = nlohmann::json(gene_name).dump();
        std::string dump = tree.dump();
        hprintf(json_wh, "    %s: %s%s\n", esc_name.c_str(), dump.c_str(),
                (gi + 1 < total_genes) ? "," : "");

        // Summary row
        std::string eps;
        for (size_t d = 0; d < summary.energy_per_depth.size(); ++d) {
            if (d > 0) eps += ",";
            char buf[64]; snprintf(buf, sizeof(buf), "%.6g", summary.energy_per_depth[d]);
            eps += buf;
        }
        if (eps.empty()) eps = "NA";
        int best_zoom = summary.best_valid ? (min_zoom + summary.best_depth) : -1;
        hprintf(summary_wh,
                "%s\t%llu\t%d\t%d\t%.6g\t%.6g\t%.6g\t%.6g\t%d\t%d\t%.6g\t%.6g\t%s\n",
                gene_name.c_str(),
                (unsigned long long)gene_total,
                summary.n_tested_nodes,
                summary.max_depth_used,
                summary.total_energy,
                summary.characteristic_scale,
                summary.scale_entropy,
                summary.best_valid ? summary.best_bic : 0.0,
                summary.best_valid ? summary.best_depth : -1,
                best_zoom,
                summary.best_chi2,
                summary.best_cramers_v,
                eps.c_str());

        // Optional PNG visualization
        if (!out_viz_prefix.empty() && viz_gene_set.count(gene_name)) {
            std::string outpng = out_viz_prefix + "." + sanitize_filename(gene_name) + ".png";
            if (render_gene_png(outpng, gene_h, bg_h, depth, viz_width, viz_height,
                                min_node_count, max_log2fc, min_opacity)) {
                notice("Wrote visualization %s", outpng.c_str());
                ++n_viz_written;
            } else {
                notice("WARNING: failed to write %s", outpng.c_str());
            }
        }

        ++gi;
        if (gi % 100 == 0 || gi == total_genes)
            notice("Computed quadtree for %zu / %zu genes", gi, total_genes);
    }

    hprintf(json_wh, "  }\n");
    hprintf(json_wh, "}\n");
    hts_close(json_wh);
    hts_close(summary_wh);

    notice("Wrote %zu gene tree(s) to %s", total_genes, out_jsonf.c_str());
    notice("Wrote per-gene summary to %s",  out_summaryf.c_str());
    if (!out_viz_prefix.empty()) notice("Wrote %d PNG visualization(s)", n_viz_written);
    notice("Analysis Finished");
    return 0;
}
