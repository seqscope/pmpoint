#ifndef __POLYGON__H__
#define __POLYGON__H__

#include <climits>
#include <cstdint>
#include <vector>
#include <fstream>
//#include <iostream>
#include "ext/nlohmann/json.hpp"
#include "qgenlib/qgen_error.h"

struct _point_t
{
    double x, y;
    _point_t(double x, double y) : x(x), y(y) {}
};
typedef struct _point_t point_t;

class Rectangle
{
public:
    point_t p_min, p_max;

    Rectangle(double xmin, double ymin, double xmax, double ymax) : p_min(xmin, ymin), p_max(xmax, ymax) {}

    inline void add_point(double x, double y) {
        if (x < p_min.x) p_min.x = x;
        if (x > p_max.x) p_max.x = x;
        if (y < p_min.y) p_min.y = y;
        if (y > p_max.y) p_max.y = y;
    }

    inline bool contains_point(double x, double y) const
    {
        return (x >= p_min.x && x <= p_max.x && y >= p_min.y && y <= p_max.y);
    }

    inline bool contains_point(const point_t &p) const
    {
        return contains_point(p.x, p.y);
    }

    inline bool contains_rectangle(const Rectangle &r) const
    {
        return (contains_point(r.p_min) && contains_point(r.p_max));
    }

    inline bool intersects_rectangle(const Rectangle &r) const
    {
        return (contains_point(r.p_min) || contains_point(r.p_max) || r.contains_point(p_min) || r.contains_point(p_max));
    }
};

class Polygon
{
public:
    std::vector<point_t> vertices;
    inline void add_offset(double x, double y)
    {
        for (size_t i = 0; i < vertices.size(); ++i)
        {
            vertices[i].x += x;
            vertices[i].y += y;
        }
    }
    inline bool contains_point(double x, double y)
    {
        size_t i, j;
        bool result = false;
        for (i = 0, j = vertices.size() - 1; i < vertices.size(); j = i++)
        {
            if ((vertices[i].y > y) != (vertices[j].y > y) &&
                (x < (vertices[j].x - vertices[i].x) * (y - vertices[i].y) / (vertices[j].y - vertices[i].y) + vertices[i].x))
            {
                result = !result;
            }
        }
        return result;
    }

    inline bool contains_point(const point_t &p)
    {
        return contains_point(p.x, p.y);
    }

    Rectangle get_bounding_box()
    {
        double xmin = std::numeric_limits<double>::max();
        double ymin = std::numeric_limits<double>::max();
        double xmax = std::numeric_limits<double>::min();
        double ymax = std::numeric_limits<double>::min();
        for (auto &vertex : vertices)
        {
            notice("Vertex: %lf, %lf, (%lg %lg %lg %lg)", vertex.x, vertex.y, xmin, ymin, xmax, ymax);
            if (vertex.x < xmin)
                xmin = vertex.x;
            if (vertex.x > xmax)
                xmax = vertex.x;
            if (vertex.y < ymin)
                ymin = vertex.y;
            if (vertex.y > ymax)
                ymax = vertex.y;
        }
        return Rectangle(xmin, ymin, xmax, ymax);
    }
};

inline int32_t add_feature_to_polygons(const nlohmann::json &feature, std::vector<Polygon> &polygons)
{
    if (!feature.contains("geometry")) {
        error("Invalid Feature: missing 'geometry' field in the GeoJSON");
    }

    auto &geometry = feature["geometry"];
    
    if (!geometry.contains("type") || !geometry.contains("coordinates")) {
        error("Invalid geometry: missing 'type' or 'coordinates' field");
    }

    if (geometry["type"] != "Polygon") {
        error("Geometry type is not Polygon");
    }

    const auto& coordinates = geometry["coordinates"];

    // Iterate through each ring (exterior and holes)
    Polygon polygon;
    for (size_t ringIndex = 0; ringIndex < coordinates.size(); ++ringIndex) {
        if (ringIndex > 0) {
            error("Holes are not supported");
        }
            
        // Iterate through points in the ring
        for (size_t pointIndex = 0; pointIndex < coordinates[ringIndex].size(); ++pointIndex) {
            const auto& point = coordinates[ringIndex][pointIndex];
            polygon.vertices.push_back(point_t(point[0], point[1]));
        }
    }
    polygons.push_back(polygon);
    return (int32_t)polygons.size();
}

inline int32_t load_polygons_from_geojson(const char *filename, std::vector<Polygon> &polygons)
{
    try {
        std::ifstream file(filename);
        if (!file.is_open())
        {
            throw std::runtime_error("Could not open file");
        }

        nlohmann::json json;
        file >> json;

        // sanity check on 'type'
        if (!json.contains("type")) {
            error("Invalid GeoJSON %s: missing 'type' field", filename);
        }

        if (json["type"] == "FeatureCollection") {
            for (const auto& feature : json["features"]) {
                add_feature_to_polygons(feature, polygons);
            }
        } else if (json["type"] == "Feature") {
            add_feature_to_polygons(json, polygons);
        } else {
            error("Invalid GeoJSON %s: Cannot find 'Feature' of 'FeatureCollection' type", filename);
        }
    }
    catch (const std::exception &e) {
        error("Error loading GeoJSON %s: %s", filename, e.what());
    }
    return (int32_t)polygons.size();
}


inline bool polygons_contain_point(std::vector<Polygon> &polygons, double x, double y)
{
    for (auto &polygon : polygons)
    {
        if (polygon.contains_point(x, y))
        {
            return true;
        }
    }
    return false;
}

#endif // __POLYGON__H__