#ifndef __PMT_UTILS_H
#define __PMT_UTILS_H

#include <cstdint>
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

        void set_global_coord(double _gx, double _gy, uint8_t _zoom);
        void set_tile_coord(uint8_t _zoom, int64_t _tx, int64_t _ty, int32_t _lx = 0, double _ly = 0);

        pmt_pt_t(uint8_t _zoom, double _gx, double _gy)
        {
            set_global_coord(_gx, _gy, _zoom);
        }
        pmt_pt_t(uint8_t _zoom, int64_t _tx, int64_t _ty, double _lx = 0, double _ly = 0)
        {
            set_tile_coord(_zoom, _tx, _ty, _lx, _ly);
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
};

#endif // __PMT_UTILS_H