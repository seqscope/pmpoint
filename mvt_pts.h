#ifndef __MVT_PTS_H
#define __MVT_PTS_H

// A class containing MVTiles object that contains a collection of points
// Performs decoding of MVTiles and coordinate transformation
// This implementation relies on the original mapbox implementation of vector_tile

#include "mapbox/vector_tile.hpp"
#include "qgenlib/qgen_error.h"
#include "polygon.h"
#include "pmt_utils.h"
#include <string>

// MVTile points with filters
class pt_dataframe
{
public:
    std::vector<std::string> feature_names;
    std::vector<std::vector<std::string>> feature_matrix;
    std::vector<pmt_utils::pmt_pt_t> points;

    inline void clear_values()
    {
        points.clear();
        for (int32_t i = 0; i < feature_matrix.size(); ++i)
        {
            feature_matrix[i].clear();
        }
    }

    inline void add_feature(int32_t idx, const std::string &name, const std::string &value)
    {
        if (idx >= feature_names.size())
        {
            feature_names.push_back(name);
            feature_matrix.resize(feature_names.size());
        }
        else if (feature_names[idx] != name)
        {
            error("Incompatible feature names. %s != %s", feature_names[idx].c_str(), name.c_str());
        }
        feature_matrix[idx].push_back(value);
    }
};

class mvt_pts_filt
{
public:
    std::vector<Polygon*> polygons;
    pmt_utils::pmt_pt_t *p_min_pt = NULL;
    pmt_utils::pmt_pt_t *p_max_pt = NULL;
    //mapbox::vector_tile::buffer *p_tile = NULL;
    //pt_dataframe *p_df;

    mvt_pts_filt(pt_dataframe *p = NULL) {} //: p_df(p) {}
    ~mvt_pts_filt()
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

    bool decode_points_df(const std::string &_buffer, uint8_t zoom, int64_t tile_x, int64_t tile_y, pt_dataframe& df);
};

class mvt_pts
{
public:
    //mapbox::vector_tile::buffer *p_tile = NULL;
    mvt_pts() {}
    ~mvt_pts()
    {
        //if (p_tile)
        //    delete p_tile;
    }

    bool decode_points_df(const std::string &_buffer, uint8_t zoom, int64_t tile_x, int64_t tile_y, pt_dataframe& df);
    uint64_t count_points(const std::string &_buffer); //, double x_offset, double y_offset, uint8_t z);
    //int decode_points_localxy(const std::string &_buffer, std::vector<int32_t>& xs, std::vector<int32_t>& ys);
    //int32_t decode_points_xycnt(const std::string &_buffer, const std::string& colname_cnt, std::vector<int32_t>& xs, std::vector<int32_t>& ys, std::vector<int32_t>& cnts);
    int32_t decode_points_xycnt_feature(const std::string &_buffer, const std::string& colname_cnt, const std::string& colname_feature, std::vector<int32_t>& xs, std::vector<int32_t>& ys, std::vector<int32_t>& cnts, std::vector<std::string>& features);
};

class print_value
{
public:
    std::string operator()(std::vector<mapbox::feature::value> val)
    {
        return "vector";
    }

    std::string operator()(std::unordered_map<std::string, mapbox::feature::value> val)
    {
        return "unordered_map";
    }

    std::string operator()(mapbox::feature::null_value_t val)
    {
        return "null";
    }

    std::string operator()(std::nullptr_t val)
    {
        return "nullptr";
    }

    std::string operator()(uint64_t val)
    {
        return std::to_string(val);
    }
    std::string operator()(int64_t val)
    {
        return std::to_string(val);
    }
    std::string operator()(double val)
    {
        return std::to_string(val);
    }
    std::string operator()(std::string const &val)
    {
        return val;
    }

    std::string operator()(bool val)
    {
        return val ? "true" : "false";
    }
};
#endif // __MVT_PTS_H