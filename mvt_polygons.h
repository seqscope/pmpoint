#ifndef __MVT_POLYGONS_H
#define __MVT_POLYGONS_H

// A class containing MVTiles object that contains a collection of points
// Performs decoding of MVTiles and coordinate transformation
// This implementation relies on the original mapbox implementation of vector_tile

#include "mapbox/vector_tile.hpp"
#include "qgenlib/qgen_error.h"
#include "polygon.h"
#include "pmt_utils.h"
#include <string>
#include <map>
#include <vector>

// MVTile points with filters
class polygon_dataframe
{
public:
    std::vector<std::string> feature_columns;
    std::map<std::string, int32_t> feature_col2idx;
    //std::vector<std::map<int32_t, std::string> > feature_values;
    std::vector< std::vector<std::string> > feature_matrix;
    std::vector<pmt_utils::pmt_polygon_t> polygons;

    inline void set_columns(const std::vector<std::string>& cols) {
        feature_columns = cols;
        feature_col2idx.clear();
        for (int32_t i = 0; i < feature_columns.size(); ++i) {
            feature_col2idx[feature_columns[i]] = i;
        }
    }

    inline void add_column_value(int32_t rowidx, const std::string &colname, const std::string &value)
    {
        std::map<std::string, int32_t>::iterator it = feature_col2idx.find(colname);
        if (it == feature_col2idx.end()) { // do nothing if column not found
            return;
        }
        if ( feature_matrix.size() <= rowidx ) {
            size_t cur_size = feature_matrix.size();
            feature_matrix.resize(rowidx + 1);
            for (int32_t i = cur_size; i < rowidx + 1; ++i) {
                feature_matrix[i].resize(feature_columns.size());
            }
        }
        int32_t colidx = it->second;
        feature_matrix[rowidx][colidx] = value;
        // int32_t colidx = -1;
        // if (it == feature_col2idx.end())
        // {
        //     feature_col2idx[colname] = feature_columns.size();
        //     feature_columns.push_back(colname);
        //     feature_values.push_back(std::map<int32_t, std::string>());
        //     colidx = feature_columns.size() - 1;
        // }
        // else
        // {
        //     colidx = it->second;
        // }

        // // add the value
        // if (feature_values.size() <= rowidx)
        // {
        //     feature_values.resize(rowidx + 1);
        // }
        // feature_values[rowidx][colidx] = value;
    }

    inline void add_column_value(const std::string &colname, const std::string &value) {
        int32_t rowidx = polygons.size();
        add_column_value(rowidx-1, colname, value);
    }

    inline void clear_values() {
        polygons.clear();
        // feature_columns.clear();
        // feature_col2idx.clear();
        feature_matrix.clear();
    }

    // inline void clear_values()
    // {
    //     polygons.clear();
    //     for (int32_t i = 0; i < feature_matrix.size(); ++i)
    //     {
    //         feature_matrix[i].clear();
    //     }
    // }

    // inline void add_feature(int32_t idx, const std::string &name, const std::string &value)
    // {
    //     if (idx >= feature_names.size())
    //     {
    //         feature_names.push_back(name);
    //         feature_matrix.resize(feature_names.size());
    //     }
    //     else if (feature_names[idx] != name)
    //     {
    //         error("Incompatible feature names. %s != %s", feature_names[idx].c_str(), name.c_str());
    //     }
    //     feature_matrix[idx].push_back(value);
    // }
};

class mvt_polygons_filt
{
public:
    std::vector<Polygon*> polygons;
    pmt_utils::pmt_pt_t *p_min_pt = NULL;
    pmt_utils::pmt_pt_t *p_max_pt = NULL;
    //mapbox::vector_tile::buffer *p_tile = NULL;
    //pt_dataframe *p_df;

    mvt_polygons_filt(polygon_dataframe *p = NULL) {} //: p_df(p) {}
    ~mvt_polygons_filt()
    {
        //if (p_tile)
        //    delete p_tile;
    }

    //inline void set_df(pt_dataframe *p) { p_df = p; }
    inline void set_min_filt(pmt_utils::pmt_pt_t *_min_pt) { p_min_pt = _min_pt; }
    inline void set_max_filt(pmt_utils::pmt_pt_t *_max_pt) { p_max_pt = _max_pt; }
    inline void set_polygon_filt(std::vector<Polygon*> &_polygons) { polygons = _polygons; }
    // inline void set_geojson_polygon_filt(const char *jsonf)
    // {
    //     int32_t npolygons = load_polygons_from_geojson(jsonf, polygons);
    //     if (npolygons == 0)
    //     {
    //         error("No polygons are loaded from %s", jsonf);
    //     }
    //     notice("Finished loading %zu polygons", polygons.size());
    // }

    bool decode_polygons_df(const std::string &_buffer, uint8_t zoom, int64_t tile_x, int64_t tile_y, polygon_dataframe& df);
};

class mvt_polygons
{
public:
    //mapbox::vector_tile::buffer *p_tile = NULL;
    mvt_polygons() {}
    ~mvt_polygons()
    {
        //if (p_tile)
        //    delete p_tile;
    }

    bool decode_polygons_df(const std::string &_buffer, uint8_t zoom, int64_t tile_x, int64_t tile_y, polygon_dataframe& df);
    uint64_t count_polygons(const std::string &_buffer); //, double x_offset, double y_offset, uint8_t z);
    //int decode_points_localxy(const std::string &_buffer, std::vector<int32_t>& xs, std::vector<int32_t>& ys);
    //int32_t decode_points_xycnt(const std::string &_buffer, const std::string& colname_cnt, std::vector<int32_t>& xs, std::vector<int32_t>& ys, std::vector<int32_t>& cnts);
    //int32_t decode_polygons_xycnt_feature(const std::string &_buffer, const std::string& colname_cnt, const std::string& colname_feature, std::vector<int32_t>& xs, std::vector<int32_t>& ys, std::vector<int32_t>& cnts, std::vector<std::string>& features);
};
#endif // __MVT_PTS_H