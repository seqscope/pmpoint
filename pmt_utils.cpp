#include "pmt_utils.h"
#include <cmath>

void pmt_utils::pmt_pt_t::set_global_coord(double _gx, double _gy, uint8_t _zoom)
{
    global_x = _gx;
    global_y = _gy;
    zoom = _zoom;

    double scale_factor = epsg3857_scale_factor(zoom);

    // get the tile info first
    epsg3857totile(global_x, global_y, zoom, &tile_x, &tile_y);
    // get the offsets for the tile
    double offset_x, offset_y;
    tiletoepsg3857(tile_x, tile_y, zoom, &offset_x, &offset_y);
    // get the local coordinates
    local_x = (global_x - offset_x) / scale_factor;
    local_y = (-global_y + offset_y) / scale_factor;
}

void pmt_utils::pmt_pt_t::set_tile_coord(uint8_t _zoom, int64_t _tx, int64_t _ty, int32_t _lx, double _ly)
{
    zoom = _zoom;
    tile_x = _tx;
    tile_y = _ty;
    local_x = _lx;
    local_y = _ly;

    double scale_factor = epsg3857_scale_factor(zoom);
    // get the offsets for the tile
    double offset_x, offset_y;
    tiletoepsg3857(tile_x, tile_y, zoom, &offset_x, &offset_y);
    // get the global coordinates
    global_x = offset_x + scale_factor * local_x;
    global_y = offset_y - scale_factor * local_y;
}

void pmt_utils::epsg3857totile(double ix, double iy, uint8_t zoom, int64_t *x, int64_t *y)
{
    // Place infinite and NaN coordinates off the edge of the Mercator plane
    int iy_class = std::fpclassify(iy);
    int ix_class = std::fpclassify(ix);

    if (iy_class == FP_INFINITE || iy_class == FP_NAN)
    {
        iy = 40000000.0;
    }
    if (ix_class == FP_INFINITE || ix_class == FP_NAN)
    {
        ix = 40000000.0;
    }

    *x = std::round(ix * (1LL << 31) / 6378137.0 / M_PI + (1LL << 31));
    *y = std::round(((1LL << 32) - 1) - (iy * (1LL << 31) / 6378137.0 / M_PI + (1LL << 31)));

    if (zoom != 0)
    {
        *x = std::round((double)*x / (1LL << (32 - zoom)));
        *y = std::round((double)*y / (1LL << (32 - zoom)));
    }
}

void pmt_utils::tiletoepsg3857(int64_t ix, int64_t iy, uint8_t zoom, double *ox, double *oy)
{
    if (zoom != 0)
    {
        ix <<= (32 - zoom);
        iy <<= (32 - zoom);
    }

    *ox = (ix - (1LL << 31)) * M_PI * 6378137.0 / (1LL << 31);
    *oy = ((1LL << 32) - 1 - iy - (1LL << 31)) * M_PI * 6378137.0 / (1LL << 31);
}