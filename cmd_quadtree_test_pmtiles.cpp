#include "pmpoint.h"
#include "qgenlib/tsv_reader.h"
#include "qgenlib/qgen_error.h"

#include <vector>
#include <string>
#include <cstring>
#include <climits>
#include <cmath>
#include <map>
#include <set>
#include <unordered_map>
#include <limits>

#include "pmt_pts.h"
#include "pmt_utils.h"
#include "polygon.h"
#include "mvt_pts.h"
#include "htslib/hts.h"
#include "ext/nlohmann/json.hpp"

// ---- MLT tile decoding helpers (duplicated from cmd_export_pmtiles.cpp) ----

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

        break;
    }
}

// ---- Quadtree statistical machinery ----

// Build a bottom-up hierarchical sum array.
// out[d] is a row-major (1<<d) x (1<<d) grid of summed values.
// out[max_depth] is taken directly from leaf_grid (which must be that size).
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

// Result of testing a single node's 4-way split: gene-of-interest vs.
// background counts arranged as a 4x2 contingency table.
struct NodeStat {
    bool   valid          = false;
    double chi2           = 0.0;
    int    df             = 3;          // (4-1)*(2-1)
    double cramers_v      = 0.0;        // sqrt(chi2 / N) since min(k1,k2)-1 = 1
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

        // log2 fold-change of gene's local fraction vs. parent fraction (Haldane pseudocount).
        double gene_frac_child = (child_gene[i] + 0.5) / (Ri + 1.0);
        s.lfc[i] = std::log2(gene_frac_child / (gene_frac_parent <= 0.0 ? 1e-12 : gene_frac_parent));
    }
    s.chi2          = chi2;
    s.g_stat        = g_stat;
    s.df            = 3;
    s.cramers_v     = std::sqrt(chi2 / N);     // 4x2 table => min(k)-1 = 1
    s.bic_penalized = g_stat - (double)s.df * std::log(N > 1.0 ? N : 2.0);
    s.valid         = true;
    return s;
}

// Per-gene scalogram accumulator: sum of chi2 contributions at each depth.
struct GeneSummary {
    std::vector<double>  energy_per_depth;
    int                  n_tested_nodes  = 0;
    int                  max_depth_used  = 0;
    double               total_energy    = 0.0;
    double               characteristic_scale = 0.0;
    double               scale_entropy   = 0.0;
};

static void accumulate_node_to_summary(const NodeStat& s, int depth, GeneSummary& sum)
{
    if (!s.valid) return;
    if ((int)sum.energy_per_depth.size() <= depth) sum.energy_per_depth.resize(depth + 1, 0.0);
    sum.energy_per_depth[depth] += s.chi2;
    sum.n_tested_nodes++;
    if (depth > sum.max_depth_used) sum.max_depth_used = depth;
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

// Recursively emit a JSON tree for a single gene.
// Children are indexed 0..3 in (dr, dc) = (top-bottom, left-right) order using row-major y-increasing coords.
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

    bool testable = (total_c >= min_node_count);
    if (testable) {
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
            accumulate_node_to_summary(s, d, summary);
        }
    }

    // Descend into children that still have enough mass to test something below.
    bool any_child_testable = false;
    for (int i = 0; i < 4; ++i) {
        if (child_total[i] >= min_node_count) { any_child_testable = true; break; }
    }
    if (d + 1 < max_depth && any_child_testable) {
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

// Convenience: locate the index of a feature column in a pt_dataframe.
static int find_feature_col(const pt_dataframe& df, const std::string& name) {
    for (int i = 0; i < (int)df.feature_names.size(); ++i) {
        if (df.feature_names[i] == name) return i;
    }
    return -1;
}

/////////////////////////////////////////////////////////////////////////
// quadtree-test : Quadtree-based spatial variability test for genes
/////////////////////////////////////////////////////////////////////////
int32_t cmd_quadtree_test_pmtiles(int32_t argc, char **argv)
{
    std::string pmtilesf;
    int32_t zoom = -1;

    // Region-based filtering
    double xmin = -std::numeric_limits<double>::infinity();
    double xmax =  std::numeric_limits<double>::infinity();
    double ymin = -std::numeric_limits<double>::infinity();
    double ymax =  std::numeric_limits<double>::infinity();
    std::string geojsonf;

    // Gene/count column names and gene-list filter
    std::string gene_col("gene");
    std::string count_col("count");
    std::string gene_listf;

    // Quadtree configuration
    int32_t max_depth      = 7;
    double  min_node_count = 20.0;

    // Outputs
    std::string out_jsonf;
    std::string out_summaryf;

    paramList pl;

    BEGIN_LONG_PARAMS(longParameters)
    LONG_PARAM_GROUP("Input options", NULL)
    LONG_STRING_PARAM("in",         &pmtilesf,   "Input PMTiles file")
    LONG_STRING_PARAM("gene-col",   &gene_col,   "Column name for gene/feature identity (default: gene)")
    LONG_STRING_PARAM("count-col",  &count_col,  "Column name for transcript counts; empty -> each point counts as 1 (default: count)")
    LONG_STRING_PARAM("gene-list",  &gene_listf, "Optional file with one gene per line to restrict the test to")

    LONG_PARAM_GROUP("Output options", NULL)
    LONG_STRING_PARAM("out-tree",    &out_jsonf,    "Output JSON file with per-gene quadtree trees (required)")
    LONG_STRING_PARAM("out-summary", &out_summaryf, "Output TSV file with per-gene summary statistics (optional)")

    LONG_PARAM_GROUP("Quadtree options", NULL)
    LONG_INT_PARAM   ("max-depth",      &max_depth,      "Maximum quadtree depth, root = 0 (default: 7)")
    LONG_DOUBLE_PARAM("min-node-count", &min_node_count, "Minimum total count in a node to run a test (default: 20)")

    LONG_PARAM_GROUP("Filtering options", NULL)
    LONG_INT_PARAM   ("zoom",    &zoom,     "Zoom level to read points from (default: max)")
    LONG_DOUBLE_PARAM("xmin",    &xmin,     "Minimum x-axis value")
    LONG_DOUBLE_PARAM("xmax",    &xmax,     "Maximum x-axis value")
    LONG_DOUBLE_PARAM("ymin",    &ymin,     "Minimum y-axis value")
    LONG_DOUBLE_PARAM("ymax",    &ymax,     "Maximum y-axis value")
    LONG_STRING_PARAM("polygon", &geojsonf, "GeoJSON file (EPSG:3857) for polygon-based filtering")
    END_LONG_PARAMS();

    pl.Add(new longParams("Available Options", longParameters));
    pl.Read(argc, argv);
    pl.Status();

    notice("Analysis started");

    if (pmtilesf.empty())    error("Missing required option --in");
    if (out_jsonf.empty())   error("Missing required option --out-tree");
    if (max_depth < 1)       error("--max-depth must be >= 1");
    if (max_depth > 12)      error("--max-depth > 12 (>=4096^2 leaves) is unreasonably large for a draft implementation");

    // Optional gene-list filter
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
        notice("Loaded %zu genes from %s", gene_filter.size(), gene_listf.c_str());
        if (gene_filter.empty()) error("--gene-list file %s contained no usable gene names", gene_listf.c_str());
    }

    // Open the PMTiles file
    pmt_pts pmt(pmtilesf.c_str());
    notice("Reading header and tile entries...");
    if (!pmt.read_header_meta_entries())
        error("This pmtiles file is malformed or incompatible with pmpoints (needs MVT or MLT points)");

    if (zoom == -1) {
        zoom = pmt.hdr.max_zoom;
        notice("Using maximum zoom level: %d", zoom);
    }
    if (zoom < pmt.hdr.min_zoom || zoom > pmt.hdr.max_zoom)
        error("Zoom level %d is unavailable in %s", zoom, pmtilesf.c_str());

    // Set up filtering boundary in tile space
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
    if (has_boundary) {
        notice("Bounding Box: [(%.3f, %.3f), (%.3f, %.3f)]", xmin, ymin, xmax, ymax);
    } else {
        notice("No bounding box; using all tiles at zoom %d", zoom);
    }

    // Load polygons + per-polygon bounding boxes
    std::vector<Polygon>   polygons;
    std::vector<Rectangle> polygon_bboxes;
    if (!geojsonf.empty()) {
        load_polygons_from_geojson(geojsonf.c_str(), polygons);
        for (auto& p : polygons) polygon_bboxes.push_back(p.get_bounding_box());
    }

    // --------------------------------------------------------------
    // Pass 1: collect bounding box of selected tiles at zoom
    // (defines the quadtree extent unless user-overridden via --xmin/--xmax/...)
    // --------------------------------------------------------------
    double bbox_x0 =  std::numeric_limits<double>::infinity();
    double bbox_y0 =  std::numeric_limits<double>::infinity();
    double bbox_x1 = -std::numeric_limits<double>::infinity();
    double bbox_y1 = -std::numeric_limits<double>::infinity();

    std::vector<int32_t> selected_tile_idxs;
    selected_tile_idxs.reserve(pmt.tile_entries.size());
    for (int32_t i = 0; i < (int32_t)pmt.tile_entries.size(); ++i) {
        pmtiles::entry_zxy& entry = pmt.tile_entries[i];
        if (entry.z != zoom) continue;

        point_t tile_min_pt(0,0), tile_max_pt(0,0);
        pmt_utils::tiletoepsg3857(entry.x,     entry.y,     entry.z, &tile_min_pt.x, &tile_max_pt.y);
        pmt_utils::tiletoepsg3857(entry.x + 1, entry.y + 1, entry.z, &tile_max_pt.x, &tile_min_pt.y);

        if (has_boundary) {
            // Skip tiles entirely outside the user-supplied filter rectangle.
            // (y-axis is inverted: tile-y increases downward while EPSG:3857-y increases upward.)
            if (entry.x < min_pt.tile_x || entry.x > max_pt.tile_x ||
                entry.y < max_pt.tile_y || entry.y > min_pt.tile_y) continue;
        }

        if (!polygon_bboxes.empty()) {
            Rectangle tile_bbox(tile_min_pt.x, tile_min_pt.y, tile_max_pt.x, tile_max_pt.y);
            bool overlaps = false;
            for (auto& r : polygon_bboxes) {
                if (r.intersects_rectangle(tile_bbox)) { overlaps = true; break; }
            }
            if (!overlaps) continue;
        }

        selected_tile_idxs.push_back(i);
        if (tile_min_pt.x < bbox_x0) bbox_x0 = tile_min_pt.x;
        if (tile_min_pt.y < bbox_y0) bbox_y0 = tile_min_pt.y;
        if (tile_max_pt.x > bbox_x1) bbox_x1 = tile_max_pt.x;
        if (tile_max_pt.y > bbox_y1) bbox_y1 = tile_max_pt.y;
    }

    if (selected_tile_idxs.empty()) error("No tiles selected at zoom %d after filters", zoom);

    // If user gave an explicit bbox, intersect it with the tile-union bbox.
    if (has_boundary) {
        if (std::isfinite(xmin) && xmin > bbox_x0) bbox_x0 = xmin;
        if (std::isfinite(ymin) && ymin > bbox_y0) bbox_y0 = ymin;
        if (std::isfinite(xmax) && xmax < bbox_x1) bbox_x1 = xmax;
        if (std::isfinite(ymax) && ymax < bbox_y1) bbox_y1 = ymax;
    }
    if (!(bbox_x0 < bbox_x1 && bbox_y0 < bbox_y1))
        error("Quadtree bounding box is empty: x=[%g,%g], y=[%g,%g]", bbox_x0, bbox_x1, bbox_y0, bbox_y1);

    notice("Quadtree extent: x=[%.3f, %.3f], y=[%.3f, %.3f]  (%d tile(s) selected)",
           bbox_x0, bbox_x1, bbox_y0, bbox_y1, (int)selected_tile_idxs.size());

    const int N         = 1 << max_depth;
    const size_t leaves = (size_t)N * (size_t)N;
    notice("Leaf grid: %d x %d (%zu cells), max_depth=%d", N, N, leaves, max_depth);

    // --------------------------------------------------------------
    // Pass 2: decode points and bin per (gene, leaf cell)
    //   - bg_leaf: total counts per leaf cell (across all genes)
    //   - gene2leaf: per-gene sparse map (leaf_idx -> count)
    // --------------------------------------------------------------
    std::vector<uint64_t> bg_leaf(leaves, 0);
    std::unordered_map<std::string, std::unordered_map<uint64_t, uint64_t>> gene2leaf;

    const bool no_count_col = count_col.empty();
    const double width  = bbox_x1 - bbox_x0;
    const double height = bbox_y1 - bbox_y0;

    pt_dataframe df;
    mvt_pts_filt mvtfilt(&df);
    std::vector<Polygon*> tile_polygons;
    std::string tile_buffer;

    uint64_t n_points_seen = 0, n_points_binned = 0;
    int      n_tiles_done  = 0;
    const int total_sel    = (int)selected_tile_idxs.size();

    for (int32_t ix : selected_tile_idxs) {
        pmtiles::entry_zxy& entry = pmt.tile_entries[ix];

        point_t tile_min_pt(0,0), tile_max_pt(0,0);
        pmt_utils::tiletoepsg3857(entry.x,     entry.y,     entry.z, &tile_min_pt.x, &tile_max_pt.y);
        pmt_utils::tiletoepsg3857(entry.x + 1, entry.y + 1, entry.z, &tile_max_pt.x, &tile_min_pt.y);
        Rectangle tile_bbox(tile_min_pt.x, tile_min_pt.y, tile_max_pt.x, tile_max_pt.y);

        if (has_boundary) {
            if (entry.x == min_pt.tile_x || entry.y == min_pt.tile_y) mvtfilt.set_min_filt(&min_pt);
            else                                                       mvtfilt.set_min_filt(NULL);
            if (entry.x == max_pt.tile_x || entry.y == max_pt.tile_y) mvtfilt.set_max_filt(&max_pt);
            else                                                       mvtfilt.set_max_filt(NULL);
        }

        if (!polygons.empty()) {
            tile_polygons.clear();
            for (size_t j = 0; j < polygon_bboxes.size(); ++j)
                if (polygon_bboxes[j].intersects_rectangle(tile_bbox))
                    tile_polygons.push_back(&polygons[j]);
            if (tile_polygons.empty()) continue;
            mvtfilt.set_polygon_filt(tile_polygons);
        }

        pmt.fetch_tile_to_buffer(entry.z, entry.x, entry.y, tile_buffer);
        df.clear_values();
        if (pmt.hdr.tile_type == 0x06) {
            decode_mlt_tile_to_df(tile_buffer, entry.z, entry.x, entry.y, df,
                                  mvtfilt.p_min_pt, mvtfilt.p_max_pt, mvtfilt.polygons);
        } else {
            mvtfilt.decode_points_df(tile_buffer, entry.z, entry.x, entry.y, df);
        }

        int gene_col_idx  = find_feature_col(df, gene_col);
        int count_col_idx = no_count_col ? -1 : find_feature_col(df, count_col);
        if (gene_col_idx < 0) {
            // Without gene labels there's nothing to test; bail loudly.
            error("Gene column '%s' not found in tile %d/%d/%d (available columns: %zu)",
                  gene_col.c_str(), entry.z, entry.x, entry.y, df.feature_names.size());
        }

        for (size_t i = 0; i < df.points.size(); ++i) {
            ++n_points_seen;
            double gx = df.points[i].global_x;
            double gy = df.points[i].global_y;
            if (gx < bbox_x0 || gx >= bbox_x1 || gy < bbox_y0 || gy >= bbox_y1) continue;

            int bx = (int)((gx - bbox_x0) / width  * (double)N);
            int by = (int)((gy - bbox_y0) / height * (double)N);
            if (bx < 0) bx = 0; if (bx >= N) bx = N - 1;
            if (by < 0) by = 0; if (by >= N) by = N - 1;
            uint64_t leaf_idx = (uint64_t)by * (uint64_t)N + (uint64_t)bx;

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

        ++n_tiles_done;
        if (n_tiles_done % 50 == 0 || n_tiles_done == total_sel)
            notice("Binned %d / %d tiles, %llu pts seen, %llu pts in test set (%zu gene(s) so far)",
                   n_tiles_done, total_sel,
                   (unsigned long long)n_points_seen, (unsigned long long)n_points_binned,
                   gene2leaf.size());
    }

    if (gene2leaf.empty())
        error("No genes matched the criteria — nothing to test");

    // --------------------------------------------------------------
    // Background hierarchy is shared across all per-gene tests.
    // --------------------------------------------------------------
    notice("Building background hierarchy...");
    std::vector<std::vector<uint64_t>> bg_h;
    build_hierarchy(bg_leaf, max_depth, bg_h);

    // --------------------------------------------------------------
    // Per-gene tree + summary
    // --------------------------------------------------------------
    htsFile* json_wh = NULL;
    if (out_jsonf.compare(out_jsonf.size() >= 3 ? out_jsonf.size() - 3 : 0, 3, ".gz") == 0)
        json_wh = hts_open(out_jsonf.c_str(), "wz");
    else
        json_wh = hts_open(out_jsonf.c_str(), "w");
    if (!json_wh) error("Failed to open output JSON file %s", out_jsonf.c_str());

    htsFile* summary_wh = NULL;
    if (!out_summaryf.empty()) {
        if (out_summaryf.compare(out_summaryf.size() >= 3 ? out_summaryf.size() - 3 : 0, 3, ".gz") == 0)
            summary_wh = hts_open(out_summaryf.c_str(), "wz");
        else
            summary_wh = hts_open(out_summaryf.c_str(), "w");
        if (!summary_wh) error("Failed to open output summary file %s", out_summaryf.c_str());
        hprintf(summary_wh,
                "gene\ttotal_count\tn_tested_nodes\tmax_depth_used\ttotal_energy\t"
                "characteristic_scale\tscale_entropy\tenergy_per_depth\n");
    }

    // Stream the top-level JSON manually so we don't have to hold all per-gene
    // trees in memory at once.
    hprintf(json_wh, "{\n");
    hprintf(json_wh, "  \"metadata\": {\n");
    hprintf(json_wh, "    \"input\": \"%s\",\n", pmtilesf.c_str());
    hprintf(json_wh, "    \"zoom\": %d,\n", zoom);
    hprintf(json_wh, "    \"max_depth\": %d,\n", max_depth);
    hprintf(json_wh, "    \"min_node_count\": %.6g,\n", min_node_count);
    hprintf(json_wh, "    \"bbox\": [%.6f, %.6f, %.6f, %.6f],\n", bbox_x0, bbox_y0, bbox_x1, bbox_y1);
    hprintf(json_wh, "    \"n_genes\": %zu,\n", gene2leaf.size());
    hprintf(json_wh, "    \"gene_col\": \"%s\",\n", gene_col.c_str());
    hprintf(json_wh, "    \"count_col\": \"%s\"\n", count_col.c_str());
    hprintf(json_wh, "  },\n");

    // Emit the background (total counts) tree once so visualizations can reference it.
    {
        std::vector<uint64_t> ones_leaf(leaves, 0);
        // The "background" tree just shows the total counts at every node; we
        // achieve this by passing bg_h as both gene_h and total_h (gene == total).
        GeneSummary dummy;
        nlohmann::json bg_tree = build_tree_json(0, 0, 0, max_depth,
                                                  bbox_x0, bbox_y0, bbox_x1, bbox_y1,
                                                  bg_h, bg_h, min_node_count, dummy);
        std::string bg_dump = bg_tree.dump();
        hprintf(json_wh, "  \"background\": %s,\n", bg_dump.c_str());
    }

    hprintf(json_wh, "  \"genes\": {\n");
    size_t gi = 0;
    const size_t total_genes = gene2leaf.size();
    for (auto& kv : gene2leaf) {
        const std::string& gene_name        = kv.first;
        const std::unordered_map<uint64_t, uint64_t>& sparse = kv.second;

        // Expand sparse gene counts into a dense leaf grid, then build hierarchy.
        std::vector<uint64_t> gene_leaf(leaves, 0);
        uint64_t gene_total = 0;
        for (auto& sk : sparse) { gene_leaf[sk.first] = sk.second; gene_total += sk.second; }
        std::vector<std::vector<uint64_t>> gene_h;
        build_hierarchy(gene_leaf, max_depth, gene_h);

        GeneSummary summary;
        nlohmann::json tree = build_tree_json(0, 0, 0, max_depth,
                                              bbox_x0, bbox_y0, bbox_x1, bbox_y1,
                                              gene_h, bg_h, min_node_count, summary);
        finalize_summary(summary);

        // Write JSON entry for the gene
        std::string esc_name = nlohmann::json(gene_name).dump(); // properly escaped quoted string
        std::string dump = tree.dump();
        hprintf(json_wh, "    %s: %s%s\n", esc_name.c_str(), dump.c_str(),
                (gi + 1 < total_genes) ? "," : "");

        if (summary_wh) {
            std::string eps;
            for (size_t d = 0; d < summary.energy_per_depth.size(); ++d) {
                if (d > 0) eps += ",";
                char buf[64]; snprintf(buf, sizeof(buf), "%.6g", summary.energy_per_depth[d]);
                eps += buf;
            }
            if (eps.empty()) eps = "NA";
            hprintf(summary_wh, "%s\t%llu\t%d\t%d\t%.6g\t%.6g\t%.6g\t%s\n",
                    gene_name.c_str(),
                    (unsigned long long)gene_total,
                    summary.n_tested_nodes,
                    summary.max_depth_used,
                    summary.total_energy,
                    summary.characteristic_scale,
                    summary.scale_entropy,
                    eps.c_str());
        }

        ++gi;
        if (gi % 100 == 0 || gi == total_genes)
            notice("Computed quadtree for %zu / %zu genes", gi, total_genes);
    }
    hprintf(json_wh, "  }\n");
    hprintf(json_wh, "}\n");
    hts_close(json_wh);
    if (summary_wh) hts_close(summary_wh);

    notice("Finished writing %zu gene trees to %s", total_genes, out_jsonf.c_str());
    if (!out_summaryf.empty())
        notice("Per-gene summary written to %s", out_summaryf.c_str());

    notice("Analysis Finished");
    return 0;
}
