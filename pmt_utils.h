#ifndef __PMT_UTILS_H
#define __PMT_UTILS_H

#include <cstdint>
#include <vector>
#include "ext/PMTiles/pmtiles.hpp"

// functions useful for projection
namespace pmt_utils
{
    class pmt_pt_t
    {
    public:
        double global_x;
        double global_y;
        uint8_t zoom;
        int64_t tile_x;
        int64_t tile_y;
        double local_x;
        double local_y;

        pmt_pt_t() : global_x(0), global_y(0), zoom(0), tile_x(0), tile_y(0), local_x(0), local_y(0) {}

        void set_global_coord(uint8_t _zoom, double _gx, double _gy);
        void set_tile_coord(uint8_t _zoom, int64_t _tx, int64_t _ty, double _lx = 0, double _ly = 0);

        pmt_pt_t(uint8_t _zoom, double _gx, double _gy)
        {
            set_global_coord(_zoom, _gx, _gy);
        }
        pmt_pt_t(uint8_t _zoom, int64_t _tx, int64_t _ty, double _lx, double _ly)
        {
            set_tile_coord(_zoom, _tx, _ty, _lx, _ly);
        }
    };

    class pmt_rect_t {
    public:
        pmt_pt_t ul;
        pmt_pt_t lr;

        pmt_rect_t() {}
        pmt_rect_t(const pmt_pt_t& _ul, const pmt_pt_t& _lr) : ul(_ul), lr(_lr) {}

        void set_global_coord(uint8_t _zoom, double _ulx, double _uly, double _lrx, double _lry) {
            ul.set_global_coord(_zoom, _ulx, _uly);
            lr.set_global_coord(_zoom, _lrx, _lry);
        }
        void set_tile_coord(uint8_t _zoom, int64_t _tx, int64_t _ty, double _ulx, double _uly, double _lrx, double _lry) {
            ul.set_tile_coord(_zoom, _tx, _ty, _ulx, _uly);
            lr.set_tile_coord(_zoom, _tx, _ty, _lrx, _lry);
        }

        bool intersects(const pmt_rect_t& other) const {
            return ul.global_x <= other.lr.global_x && lr.global_x >= other.ul.global_x &&
                   ul.global_y <= other.lr.global_y && lr.global_y >= other.ul.global_y;
        }

        bool contains(const pmt_rect_t& other) const {
            return ul.global_x <= other.ul.global_x && lr.global_x >= other.lr.global_x &&
                   ul.global_y <= other.ul.global_y && lr.global_y >= other.lr.global_y;
        }

        bool contains_point(const pmt_pt_t& pt) const {
            return pt.global_x >= ul.global_x && pt.global_x <= lr.global_x &&
                   pt.global_y >= ul.global_y && pt.global_y <= lr.global_y;
        }

        bool is_valid() const {
            return ul.global_x <= lr.global_x && ul.global_y <= lr.global_y;
        }
    };

    class pmt_polygon_t {
    public:
        uint8_t zoom;
        std::vector<pmt_pt_t> points;
        pmt_rect_t bbox;

        void set_zoom(uint8_t _zoom) { zoom = _zoom; }

        pmt_polygon_t(uint8_t _zoom = 0) {
            set_zoom(_zoom);
        }

        void update_bbox() {
            double gx = points.back().global_x;
            double gy = points.back().global_y;
            if (points.size() == 1) {
                bbox.set_global_coord(zoom, gx, gy, gx, gy);
            } else {
                double ulx = gx < bbox.ul.global_x ? gx : bbox.ul.global_x;
                double uly = gy < bbox.ul.global_y ? gy : bbox.ul.global_y;
                double lrx = gx > bbox.lr.global_x ? gx : bbox.lr.global_x;
                double lry = gy > bbox.lr.global_y ? gy : bbox.lr.global_y;
                bbox.set_global_coord(zoom, ulx, uly, lrx, lry);
            }
        }

        void add_global_coord(double _gx, double _gy) {
            points.push_back(pmt_pt_t(zoom, _gx, _gy));
            update_bbox();
        }
        void add_tile_coord(int64_t _tx, int64_t _ty, int32_t _lx = 0, double _ly = 0) {
            points.push_back(pmt_pt_t(zoom, _tx, _ty, _lx, _ly));
            update_bbox();
        }
        void add_point(const pmt_pt_t& pt) {
            points.push_back(pt);
            update_bbox();
        }
    };

    // radius of EPSG:3857 in meters
    const double EPSG_3857_radius = 6378137.0;

    // bound of EPSG:3857 in meters
    const double EPSG_3857_bound = 20037508.3428; // = EPSG_3857_radius * M_PI

    // scale factor to convert EPSG:3857 to original coordinates
    inline double epsg3857_scale_factor(uint8_t z)
    {
        return 2 * EPSG_3857_bound / (1 << (z + 12));
    }
    void epsg3857totile(double ix, double iy, uint8_t zoom, int64_t *x, int64_t *y);
    void tiletoepsg3857(int64_t ix, int64_t iy, uint8_t zoom, double *ox, double *oy);
    void epsg3857totilecoord(double ix, double iy, uint8_t zoom, int64_t *tx, int64_t *ty, double *lx, double *ly);
    void tilecoordtoespg3857(int64_t tx, int64_t ty, double lx, double ly, uint8_t zoom, double *ox, double* oy);
};

#endif // __PMT_UTILS_H