#include "pmt_utils.h"
#include <cmath>
#include "qgenlib/qgen_error.h"

void pmt_utils::pmt_pt_t::set_global_coord(uint8_t _zoom, double _gx, double _gy)
{
    global_x = _gx;
    global_y = _gy;
    zoom = _zoom;

    double scale_factor = epsg3857_scale_factor(zoom);

    //double gx_scaled = global_x / scale_factor;
    //double gy_scaled = global_y / scale_factor;

    // get the tile info first
    epsg3857totilecoord(global_x, global_y, zoom, &tile_x, &tile_y, &local_x, &local_y);
    // epsg3857totile(global_x, global_y, zoom, &tile_x, &tile_y);
    // // get the offsets for the tile
    // double offset_x, offset_y;
    // tiletoepsg3857(tile_x, tile_y, zoom, &offset_x, &offset_y);
    // // get the local coordinates
    // local_x = (global_x - offset_x) / scale_factor;
    // local_y = (-global_y + offset_y) / scale_factor;

    //if ( ( local_x < 0 ) || ( local_y < 0) ) {
    //notice("zoom = %u, global_x = %.3f, global_y = %.3f, tile_x = %lld, tile_y = %lld, local_x = %.3f, local_y = %.3f, scale_factor = %.3f", zoom, global_x, global_y, tile_x, tile_y, local_x, local_y, scale_factor);
    //}
}

void pmt_utils::pmt_pt_t::set_tile_coord(uint8_t _zoom, int64_t _tx, int64_t _ty, int32_t _lx, double _ly)
{
    zoom = _zoom;
    tile_x = _tx;
    tile_y = _ty;
    local_x = _lx;
    local_y = _ly;

    tilecoordtoespg3857(tile_x, tile_y, local_x, local_y, zoom, &global_x, &global_y);
}

void pmt_utils::epsg3857totilecoord(double ix, double iy, uint8_t zoom, int64_t *tx, int64_t *ty, double *lx, double *ly)
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

    const double tileSize = 256;
    // Calculate the number of tiles at the given zoom level
    uint64_t numTiles = 1 << zoom;

    // Convert the global position to tile coordinates
    *tx = static_cast<int64_t>((ix + EPSG_3857_bound) / (2.0 * EPSG_3857_bound / numTiles));
    *ty = static_cast<int64_t>((EPSG_3857_bound - iy) / (2.0 * EPSG_3857_bound / numTiles));

    // Calculate the relative position within the tile
    double tileOriginX = (*tx) * (2.0 * EPSG_3857_bound / numTiles) - EPSG_3857_bound;
    double tileOriginY = EPSG_3857_bound - (*ty) * (EPSG_3857_bound * 2.0 / numTiles);

    *lx = (ix - tileOriginX) / (2. * EPSG_3857_bound / (numTiles * tileSize));
    *ly = (tileOriginY - iy) / (2. * EPSG_3857_bound / (numTiles * tileSize));
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

void pmt_utils::tilecoordtoespg3857(int64_t tx, int64_t ty, double lx, double ly, uint8_t zoom, double *ox, double* oy) {
    const double tileSize = 256;

    // Calculate the number of tiles at the given zoom level
    uint64_t numTiles = (1 << zoom);

    // Calculate the global coordinates of the tile origin
    double tileOriginX = tx * (2 * EPSG_3857_bound / numTiles) - EPSG_3857_bound;
    double tileOriginY = EPSG_3857_bound - ty * (2 * EPSG_3857_bound / numTiles);

    // Calculate the global coordinates from the tile coordinates and relative position
    *ox = tileOriginX + lx * (2 * EPSG_3857_bound / (numTiles * tileSize));
    *oy = tileOriginY - ly * (2 * EPSG_3857_bound / (numTiles * tileSize));
}