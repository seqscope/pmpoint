#include <cstdio>
#include <string>
#include <zlib.h>

#include "mvt_pts.h"
#include "mvt_polygons.h"
#include "pmt_utils.h"
#include <fstream>
#include <stdexcept>
#include <iostream>

#include "qgenlib/qgen_error.h"

// int32_t mvt_pts::decode_points_xycnt_feature(const std::string &_buffer, const std::string& colname_cnt, const std::string& colname_feature, std::vector<int32_t>& xs, std::vector<int32_t>& ys, std::vector<int32_t>& cnts, std::vector<std::string>& features)
// {
//     // p_tile = new mapbox::vector_tile::buffer(_buffer);
//     // if (p_tile == NULL)
//     // {
//     //     return false;
//     // }

//     mapbox::vector_tile::buffer* local_p_tile = new mapbox::vector_tile::buffer(_buffer);
//     if (local_p_tile == NULL)
//     {
//         return false;
//     }

//     // assuming that all objects are points with rectangular structure, try to decode the tile
//     std::vector<std::string> colnames;
//     colnames.push_back("X");
//     colnames.push_back("Y");
//     std::vector<std::vector<std::string>> columns(2);
//     int32_t n_points = 0;

//     for (auto const &name : local_p_tile->layerNames())
//     {
//         const mapbox::vector_tile::layer layer = local_p_tile->getLayer(name);
//         std::size_t feature_count = layer.featureCount();
//         if (feature_count > 0)
//         {
//             for (std::size_t i = 0; i < feature_count; ++i)
//             {
//                 auto const feature = mapbox::vector_tile::feature(layer.getFeature(i), layer);

//                 // make sure that the types are points
//                 if (int(feature.getType()) != 1)
//                 {
//                     error("Only points are supported in decode_points()");
//                 }

//                 // obtain vertex
//                 mapbox::vector_tile::points_arrays_type geom = feature.getGeometries<mapbox::vector_tile::points_arrays_type>(1.0);
//                 if (geom.size() != 1)
//                 {
//                     error("Only single point per feature is supported in decode_points()");
//                 }
//                 xs.push_back(geom[0][0].x);
//                 ys.push_back(geom[0][0].y);

//                 // obtain properties;
//                 auto props = feature.getProperties();
//                 int32_t j = 2;
//                 for (auto const &prop : props)
//                 {
//                     if ( colname_cnt.compare(prop.first) == 0 ) {
//                         print_value printvisitor;
//                         std::string value = mapbox::util::apply_visitor(printvisitor, prop.second);
//                         try {
//                             cnts.push_back(std::stoi(value));
//                         } catch (std::invalid_argument& e) {
//                             notice("Invalid count value %s observed at %d/%d, considering as zero count", value.c_str(), geom[0][0].x, geom[0][0].y);
//                             cnts.push_back(0);
//                         }
//                     }
//                     else if ( colname_feature.compare(prop.first) == 0 ) {
//                         print_value printvisitor;
//                         std::string value = mapbox::util::apply_visitor(printvisitor, prop.second);
//                         features.push_back(value);
//                     }
//                 }
//                 if ( xs.size() != ys.size() || xs.size() != cnts.size() || xs.size() != features.size() ) {
//                     error("Inconsistent sizes of xs, ys, cnts, and features. xs=%zu, ys=%zu, cnts=%zu, features=%zu", xs.size(), ys.size(), cnts.size(), features.size());
//                 }
//             }
//         }
//         n_points += feature_count;
//     }
//     //delete p_tile;
//     //p_tile = NULL;
//     delete local_p_tile;

//    return n_points;
// }

bool mvt_polygons_filt::decode_polygons_df(const std::string &_buffer, uint8_t zoom, int64_t tile_x, int64_t tile_y, polygon_dataframe& df)
{
    mapbox::vector_tile::buffer* p_tile = new mapbox::vector_tile::buffer(_buffer);
    if (p_tile == NULL)
    {
        return false;
    }

    uint64_t npass = 0, nskip = 0;
    double scale_factor = pmt_utils::epsg3857_scale_factor(zoom);
    double offset_x, offset_y;
    pmt_utils::tiletoepsg3857(tile_x, tile_y, zoom, &offset_x, &offset_y);
    notice("tile_x = %ld, tile_y = %ld, zoom = %d, offset_x = %.5f, offset_y = %.5f", tile_x, tile_y, zoom, offset_x, offset_y);
    for (auto const &name : p_tile->layerNames())
    {
        const mapbox::vector_tile::layer layer = p_tile->getLayer(name);
        std::size_t feature_count = layer.featureCount();
        if (feature_count > 0)
        {
            for (std::size_t i = 0; i < feature_count; ++i)
            {
                //notice("Processing polygon %zu", i);
                auto const feature = mapbox::vector_tile::feature(layer.getFeature(i), layer);

                // make sure that the types are polygons
                if (int(feature.getType()) != mapbox::vector_tile::GeomType::POLYGON)
                {
                    error("Only polygons are supported in decode_polygons()");
                }

                // obtain vertex
                mapbox::vector_tile::points_arrays_type geom = feature.getGeometries<mapbox::vector_tile::points_arrays_type>(1.0);

                // make sure that the polygons only have one ring
                if ( geom.size() != 1 ) {
                    error("Only simple polygons (single outer ring) are supported, but found features with %zu rings", geom.size());
                }

                if ( geom[0].empty() ) {
                    error("Empty polygon found");
                }

                // obtain the polygons
                pmt_utils::pmt_polygon_t poly(zoom);
                // pmt_utils::pmt_pt_t ul, lr;
                for (auto const &pt : geom[0]) {
                    //  pmt_utils::pmt_pt_t pt(zoom, offset_x + scale_factor * geom[0][0].x, offset_y - scale_factor * geom[0][0].y);
                    poly.add_global_coord(offset_x + scale_factor * pt.x, offset_y - scale_factor * pt.y);
                    //poly.add_point(pmt_utils::pmt_pt_t(zoom, offset_x + scale_factor * pt.x, offset_y - scale_factor * pt.y));

                    //notice("Adding point: %.5f, %.5f", poly.points.back().global_x, poly.points.back().global_y);
                    //notice("Adding point: %d, %d, %.5f, %.5f, %.5f, (%.5f, %.5f)", pt.x, pt.y, offset_x, offset_y, scale_factor, poly.points.back().global_x, poly.points.back().global_y);

                    // // update the bounding box 
                    // if ( poly.points.size() == 1 ) {
                    //     ul = poly.points[0];
                    //     lr = poly.points[0];
                    // }
                    // else {
                    //     if ( poly.points.back().global_x < ul.global_x ) {
                    //         ul.global_x = poly.points.back().global_x;
                    //     }
                    //     if ( poly.points.back().global_y < ul.global_y ) {
                    //         ul.global_y = poly.points.back().global_y;
                    //     }
                    //     if ( poly.points.back().global_x > lr.global_x ) {
                    //         lr.global_x = poly.points.back().global_x;
                    //     }
                    //     if ( poly.points.back().global_y > lr.global_y ) {
                    //         lr.global_y = poly.points.back().global_y;
                    //     }
                    // }
                }

                // check the filtering criteria
                if (p_min_pt != NULL)
                {
                    if (poly.bbox.lr.global_x < p_min_pt->global_x || poly.bbox.lr.global_y < p_min_pt->global_y)
                    {
                        ++nskip;
                        continue;
                    }
                }
                if (p_max_pt != NULL)
                {
                    if (poly.bbox.ul.global_x > p_max_pt->global_x || poly.bbox.ul.global_y > p_max_pt->global_y)
                    {
                        ++nskip;
                        continue;
                    }
                }

                // Polygon filtering to be implemented later. 
                // check overlaps based on the bounding box first, then check the actual polygon
                if (polygons.size() > 0)
                {
                    //error("Polygon filtering is not implemented yet.");
                    bool found = false;
                    for (auto &p_polygon : polygons)
                    {
                        bool ul_in = p_polygon->contains_point(poly.bbox.ul.global_x, poly.bbox.ul.global_y);
                        bool lr_in = p_polygon->contains_point(poly.bbox.lr.global_x, poly.bbox.lr.global_y);
                        if (ul_in && lr_in)
                        {
                            found = true;
                            break;
                        }
                        // if any point in the polygon is inside the bounding box, then it is a match 
                        for (auto &pt : poly.points)
                        {
                            if (p_polygon->contains_point(pt.global_x, pt.global_y))
                            {
                                found = true;
                                break;
                            }
                        } 
                        if (found)
                        {
                            break;
                        }
                    }
                    if (!found) // no match - not quite exact, but close enough
                    {
                        ++nskip;
                        continue;
                    }
                }

                ++npass;

                //df.points.push_back(pt);
                // for (auto &pt : poly.points) {
                //     notice("Adding point: %.5f, %.5f", pt.global_x, pt.global_y);
                // }
                df.polygons.push_back(poly);

                // obtain properties;
                auto props = feature.getProperties();
                int32_t j = 0;
                for (auto const &prop : props)
                {
                    print_value printvisitor;
                    std::string value = mapbox::util::apply_visitor(printvisitor, prop.second);
                    df.add_column_value(prop.first, value);
                    ++j;
                }
            }
            //notice("npass = %llu, nskip = %llu", npass, nskip);
        }
    }

    delete p_tile;
    //p_tile = NULL;
    return true;
}

// bool mvt_polygons::decode_polygons_df(const std::string &_buffer, uint8_t zoom, int64_t tile_x, int64_t tile_y, pt_dataframe& df)
// {
//     mapbox::vector_tile::buffer* p_tile = new mapbox::vector_tile::buffer(_buffer);
//     if (p_tile == NULL)
//     {
//         return false;
//     }

//     // assuming that all objects are points with rectangular structure, try to decode the tile
//     // std::vector<std::string> colnames;
//     // colnames.push_back("X");
//     // colnames.push_back("Y");
//     // std::vector<std::vector<std::string>> columns(2);
//     int32_t n_points = 0;
//     // char bufx[255], bufy[255];

//     // scale factor for EPSG:3857 to original coordinates
//     // double scale_factor = 2 * pmt_utils::EPSG_3857_bound / ( 1 << (z + 12) );
//     double scale_factor = pmt_utils::epsg3857_scale_factor(zoom);
//     double offset_x, offset_y;
//     pmt_utils::tiletoepsg3857(tile_x, tile_y, zoom, &offset_x, &offset_y);

//     for (auto const &name : p_tile->layerNames())
//     {
//         const mapbox::vector_tile::layer layer = p_tile->getLayer(name);
//         std::size_t feature_count = layer.featureCount();
//         if (feature_count > 0)
//         {
//             for (std::size_t i = 0; i < feature_count; ++i)
//             {
//                 auto const feature = mapbox::vector_tile::feature(layer.getFeature(i), layer);

//                 // make sure that the types are points
//                 if (int(feature.getType()) != 1)
//                 {
//                     error("Only points are supported in decode_points()");
//                 }

//                 // obtain vertex
//                 mapbox::vector_tile::points_arrays_type geom = feature.getGeometries<mapbox::vector_tile::points_arrays_type>(1.0);
//                 if (geom.size() != 1)
//                 {
//                     error("Only single point per feature is supported in decode_points()");
//                 }
//                 pmt_utils::pmt_pt_t pt(zoom, offset_x + scale_factor * geom[0][0].x, offset_y - scale_factor * geom[0][0].y);

//                 df.points.push_back(pt);

//                 // obtain properties;
//                 auto props = feature.getProperties();
//                 int32_t j = 0;
//                 for (auto const &prop : props)
//                 {
//                     print_value printvisitor;
//                     std::string value = mapbox::util::apply_visitor(printvisitor, prop.second);
//                     df.add_feature(j, prop.first, value);
//                     ++j;
//                 }

//                 // snprintf(bufx, sizeof(bufx), "%.3f", x_offset + scale_factor * geom[0][0].x);
//                 // snprintf(bufy, sizeof(bufy), "%.3f", y_offset - scale_factor * geom[0][0].y);
//                 // columns[0].push_back(bufx);
//                 // columns[1].push_back(bufy);

//                 // obtain properties;
//                 // auto props = feature.getProperties();
//                 // int32_t j = 2;
//                 // for (auto const &prop : props)
//                 // {
//                 //     print_value printvisitor;
//                 //     std::string value = mapbox::util::apply_visitor(printvisitor, prop.second);
//                 //     if (i == 0)
//                 //     { // for first columns fill in column names
//                 //         colnames.push_back(prop.first);
//                 //         columns.resize(colnames.size());
//                 //     }
//                 //     else
//                 //     {
//                 //         if (colnames[j].compare(prop.first) != 0)
//                 //         {
//                 //             error("Mismatched column names - %s vs %s at index = %d, row = %d in decode_points()", colnames[i].c_str(), prop.first.c_str(), i, j);
//                 //         }
//                 //     }
//                 //     columns[j].push_back(value);
//                 //     ++j;
//                 // }
//             }
//         }
//         n_points += feature_count;
//     }

//     // // print the decoded points
//     // printf("%s", colnames[0].c_str());
//     // for (int32_t i = 1; i < colnames.size(); ++i)
//     // {
//     //     printf("\t%s", colnames[i].c_str());
//     // }
//     // printf("\n");

//     // // print the values
//     // for (int32_t i = 0; i < n_points; ++i)
//     // {
//     //     printf("%s", columns[0][i].c_str());
//     //     for (int32_t j = 1; j < colnames.size(); ++j)
//     //     {
//     //         printf("\t%s", columns[j][i].c_str());
//     //     }
//     //     printf("\n");
//     // }

//     delete p_tile;
//     //p_tile = NULL;

//    return true;
// }

// uint64_t mvt_pts::count_points(const std::string &_buffer)
// {
//     mapbox::vector_tile::buffer* p_tile = new mapbox::vector_tile::buffer(_buffer);
//     if (p_tile == NULL)
//     {
//         return false;
//     }

//     uint64_t n_points = 0;
//     for (auto const &name : p_tile->layerNames())
//     {
//         const mapbox::vector_tile::layer layer = p_tile->getLayer(name);
//         std::size_t feature_count = layer.featureCount();
//         n_points += feature_count;
//     }

//     delete p_tile;
//     return n_points;
// }