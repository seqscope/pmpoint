#include <cstdio>
#include <string>
#include <zlib.h>

#include "mvt_pts.h"
#include "pmt_utils.h"
#include <fstream>
#include <stdexcept>
#include <iostream>

#include "qgenlib/qgen_error.h"

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

int32_t mvt_pts::decode_points_xycnt(const std::string &_buffer, const std::string& colname_cnt, std::vector<int32_t>& xs, std::vector<int32_t>& ys, std::vector<int32_t>& cnts)
{
    // p_tile = new mapbox::vector_tile::buffer(_buffer);
    // if (p_tile == NULL)
    // {
    //     return false;
    // }

    mapbox::vector_tile::buffer* local_p_tile = new mapbox::vector_tile::buffer(_buffer);
    if (local_p_tile == NULL)
    {
        return false;
    }

    // assuming that all objects are points with rectangular structure, try to decode the tile
    std::vector<std::string> colnames;
    colnames.push_back("X");
    colnames.push_back("Y");
    std::vector<std::vector<std::string>> columns(2);
    int32_t n_points = 0;

    for (auto const &name : local_p_tile->layerNames())
    {
        const mapbox::vector_tile::layer layer = local_p_tile->getLayer(name);
        std::size_t feature_count = layer.featureCount();
        if (feature_count > 0)
        {
            for (std::size_t i = 0; i < feature_count; ++i)
            {
                auto const feature = mapbox::vector_tile::feature(layer.getFeature(i), layer);

                // make sure that the types are points
                if (int(feature.getType()) != 1)
                {
                    error("Only points are supported in decode_points()");
                }

                // obtain vertex
                mapbox::vector_tile::points_arrays_type geom = feature.getGeometries<mapbox::vector_tile::points_arrays_type>(1.0);
                if (geom.size() != 1)
                {
                    error("Only single point per feature is supported in decode_points()");
                }
                xs.push_back(geom[0][0].x);
                ys.push_back(geom[0][0].y);

                // obtain properties;
                auto props = feature.getProperties();
                int32_t j = 2;
                for (auto const &prop : props)
                {
                    if ( colname_cnt.compare(prop.first) == 0 ) {
                        print_value printvisitor;
                        std::string value = mapbox::util::apply_visitor(printvisitor, prop.second);
                        try {
                            cnts.push_back(std::stoi(value));
                        } catch (std::invalid_argument& e) {
                            notice("Invalid count value %s observed at %d/%d, considering as zero count", value.c_str(), geom[0][0].x, geom[0][0].y);
                            cnts.push_back(0);
                        }
                    }
                }                
            }
        }
        n_points += feature_count;
    }
    //delete p_tile;
    //p_tile = NULL;
    delete local_p_tile;

   return n_points;
}

bool mvt_pts_filt::decode_points(const std::string &_buffer, uint8_t zoom, int64_t tile_x, int64_t tile_y)
{
    p_tile = new mapbox::vector_tile::buffer(_buffer);
    if (p_tile == NULL)
    {
        return false;
    }

    uint64_t npass = 0, nskip = 0;
    double scale_factor = pmt_utils::epsg3857_scale_factor(zoom);
    double offset_x, offset_y;
    pmt_utils::tiletoepsg3857(tile_x, tile_y, zoom, &offset_x, &offset_y);
    for (auto const &name : p_tile->layerNames())
    {
        const mapbox::vector_tile::layer layer = p_tile->getLayer(name);
        std::size_t feature_count = layer.featureCount();
        if (feature_count > 0)
        {
            for (std::size_t i = 0; i < feature_count; ++i)
            {
                auto const feature = mapbox::vector_tile::feature(layer.getFeature(i), layer);

                // make sure that the types are points
                if (int(feature.getType()) != 1)
                {
                    error("Only points are supported in decode_points()");
                }

                // obtain vertex
                mapbox::vector_tile::points_arrays_type geom = feature.getGeometries<mapbox::vector_tile::points_arrays_type>(1.0);
                if (geom.size() != 1)
                {
                    error("Only single point per feature is supported in decode_points()");
                }

                pmt_utils::pmt_pt_t pt(zoom, offset_x + scale_factor * geom[0][0].x, offset_y - scale_factor * geom[0][0].y);
                //pmt_utils::pmt_pt_t pt(zoom, tile_x, tile_y, geom[0][0].x/scale_factor, geom[0][0].y/scale_factor);
                //pmt_utils::pmt_pt_t pt(zoom, tile_x, tile_y, 0, 0);
                //notice("Processing point (%.3f %.3f) at z=%u, (%llu, %llu), (%.3f, %.3f)", pt.global_x, pt.global_y, pt.zoom, pt.tile_x, pt.tile_y, pt.local_x, pt.local_y);

                // check the filtering criteria
                if (p_min_pt != NULL)
                {
                    //notice("(%.3f %.3f) > (%.3f %.3f)", pt.global_x, pt.global_y, p_min_pt->global_x, p_min_pt->global_y);
                    if (pt.global_x < p_min_pt->global_x || pt.global_y < p_min_pt->global_y)
                    {
                        ++nskip;
                        continue;
                    }
                }
                if (p_max_pt != NULL)
                {
                    //notice("(%.3f %.3f) < (%.3f %.3f)", pt.global_x, pt.global_y, p_max_pt->global_x, p_max_pt->global_y);
                    if (pt.global_x > p_max_pt->global_x || pt.global_y > p_max_pt->global_y)
                    {
                        ++nskip;
                        continue;
                    }
                }
                if (polygons.size() > 0)
                {
                    bool found = false;
                    for (auto &p_polygon : polygons)
                    {
                        if (p_polygon->contains_point(pt.global_x, pt.global_y))
                        {
                            found = true;
                            break;
                        }
                    }
                    if (!found)
                    {
                        ++nskip;
                        continue;
                    }
                }

                ++npass;

                p_df->points.push_back(pt);

                // obtain properties;
                auto props = feature.getProperties();
                int32_t j = 0;
                for (auto const &prop : props)
                {
                    print_value printvisitor;
                    std::string value = mapbox::util::apply_visitor(printvisitor, prop.second);
                    p_df->add_feature(j, prop.first, value);
                    ++j;
                }
            }
            //notice("npass = %llu, nskip = %llu", npass, nskip);
        }
    }

    delete p_tile;
    p_tile = NULL;
    return true;
}

bool mvt_pts::decode_points(const std::string &_buffer, double x_offset, double y_offset, uint8_t z)
{
    p_tile = new mapbox::vector_tile::buffer(_buffer);
    if (p_tile == NULL)
    {
        return false;
    }

    // assuming that all objects are points with rectangular structure, try to decode the tile
    std::vector<std::string> colnames;
    colnames.push_back("X");
    colnames.push_back("Y");
    std::vector<std::vector<std::string>> columns(2);
    int32_t n_points = 0;
    char bufx[255], bufy[255];

    // scale factor for EPSG:3857 to original coordinates
    // double scale_factor = 2 * pmt_utils::EPSG_3857_bound / ( 1 << (z + 12) );
    double scale_factor = pmt_utils::epsg3857_scale_factor(z);

    for (auto const &name : p_tile->layerNames())
    {
        const mapbox::vector_tile::layer layer = p_tile->getLayer(name);
        std::size_t feature_count = layer.featureCount();
        if (feature_count > 0)
        {
            for (std::size_t i = 0; i < feature_count; ++i)
            {
                auto const feature = mapbox::vector_tile::feature(layer.getFeature(i), layer);

                // make sure that the types are points
                if (int(feature.getType()) != 1)
                {
                    error("Only points are supported in decode_points()");
                }

                // obtain vertex
                mapbox::vector_tile::points_arrays_type geom = feature.getGeometries<mapbox::vector_tile::points_arrays_type>(1.0);
                if (geom.size() != 1)
                {
                    error("Only single point per feature is supported in decode_points()");
                }

                snprintf(bufx, sizeof(bufx), "%.3f", x_offset + scale_factor * geom[0][0].x);
                snprintf(bufy, sizeof(bufy), "%.3f", y_offset - scale_factor * geom[0][0].y);
                columns[0].push_back(bufx);
                columns[1].push_back(bufy);

                // obtain properties;
                auto props = feature.getProperties();
                int32_t j = 2;
                for (auto const &prop : props)
                {
                    print_value printvisitor;
                    std::string value = mapbox::util::apply_visitor(printvisitor, prop.second);
                    if (i == 0)
                    { // for first columns fill in column names
                        colnames.push_back(prop.first);
                        columns.resize(colnames.size());
                    }
                    else
                    {
                        if (colnames[j].compare(prop.first) != 0)
                        {
                            error("Mismatched column names - %s vs %s at index = %d, row = %d in decode_points()", colnames[i].c_str(), prop.first.c_str(), i, j);
                        }
                    }
                    columns[j].push_back(value);
                    ++j;
                }
            }
        }
        n_points += feature_count;
    }

    // // print the decoded points
    // printf("%s", colnames[0].c_str());
    // for (int32_t i = 1; i < colnames.size(); ++i)
    // {
    //     printf("\t%s", colnames[i].c_str());
    // }
    // printf("\n");

    // // print the values
    // for (int32_t i = 0; i < n_points; ++i)
    // {
    //     printf("%s", columns[0][i].c_str());
    //     for (int32_t j = 1; j < colnames.size(); ++j)
    //     {
    //         printf("\t%s", columns[j][i].c_str());
    //     }
    //     printf("\n");
    // }

    delete p_tile;
    p_tile = NULL;

   return true;
}

uint64_t mvt_pts::count_points(const std::string &_buffer, double x_offset, double y_offset, uint8_t z)
{
    p_tile = new mapbox::vector_tile::buffer(_buffer);
    if (p_tile == NULL)
    {
        return false;
    }

    uint64_t n_points = 0;
    for (auto const &name : p_tile->layerNames())
    {
        const mapbox::vector_tile::layer layer = p_tile->getLayer(name);
        std::size_t feature_count = layer.featureCount();
        n_points += feature_count;
    }

    return n_points;
}