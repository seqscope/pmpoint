# Tools in pmpoint

`pmpoint` is a collection of C++ tools that are designed to help process, visualize, and
analyze sub-micron resolution spatial transcriptomics data stored in PMTiles format. 
The categories of software tools provided includes:

* PMTiles Tools
  * [`pmpoint export`](tools/export.md): Extract pixel-level point data from a PMTiles file with spatial filtering options.
  * [`pmpoint summarize`](tools/summarize.md): Summarize the contents of a PMTiles file, including header information, metadata, and tile statistics.
  * [`pmpoint count-tiles`](tools/count_tiles.md): Count the number of points in each tile from a PMTiles file, optionally filtering by zoom level.
  * [`pmpoint tile-density-stats` and `pmpoint tile-density-stats-mt`](tools/tile_density_stats.md): Compute the 2D spatial density statistics in square grids for each tile in a PMTiles file

!!! note
    This documentation contains only a subset of primary tools implemented `pmpoint`.
    Other tools not documented here are considered non-primary tools that have limited support, and they are NOT supported officially.