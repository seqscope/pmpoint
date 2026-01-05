# pmpoint tile-density-stats

## Summary 

`pmpoint tile-density-stats` and `pmpoint tile-density-stats-mt` computes the 2D spatial density statistics in square grids for each tile in a PMTiles file. The latter is a multi-threaded version that can process multiple tiles in parallel.

An example command is given below:

```bash
## Example command to count tiles a remote PMTiles file at zoom level 14
pmpoint tile-density-stats-mt --in https://cartostore.s3.us-east-1.amazonaws.com/data/batch=2025_12/mouse-brain-test-collection/subdata_seqscope_n17t89b_c1d2f/genes_all.pmtiles --count count --out out.tsv.gz --zoom 15 --threads 10
```

## Required options

* `--in`: Input PMTiles file. The file can be either local file or a URL to a remote file (supports HTTP/HTTPS).
* `--out-tsv`: Output TSV file to store the tile counts.

## Additional Options

* `--count`: Field name for transcript counts in the PMTiles file. Default is `count`.
* `--feature`: Field name for feature name in the PMTiles file. Default is `gene`.
* `--compact`: If set, skips writing each tile, and only report aggregated density metrics across all zoom levels.
* `--zoom`: Zoom level to count tiles. Default is -1, which counts density metrics for the highest zoom level. If a specific zoom level is provided, only tiles from that zoom level will be counted.

## Expected Output

The output TSV file contains the following columns:
* `zoom`: Zoom level of the tile.
* `tile_x`: X (integer) coordinate of the tile according to the PMTiles specification.
* `tile_y`: Y (integer) coordinate of the tile according to the PMTiles specification.
* `width` : The width (in the lowest pixel resolution in the zoom level) of square grids 
* `num_pts`: Number of points contained in the square grid (considered as the key of histogram bin).
* `num_grids`: Number of square grids that contain `num_pts` points. 

## Full Usage 

The full usage of `pmpoint count-tiles` can be viewed with the `--help` option:

```
$ ./pmpoint tile-density-stats-mt --help 
[bin/pmpoint tile-density-stats-mt] -- Calculate statistics of tile densities

 Copyright (c) 2022-2025 by Hyun Min Kang
 Licensed under the Apache License v2.0 http://www.apache.org/licenses/

Detailed instructions of parameters are available. Ones with "[]" are in effect:

Available Options:

== Input options ==
   --in      [STR: ]             : Input PMTiles file
   --count   [STR: gn]           : Field name for transcript counts
   --feature [STR: gene]         : Field name for feature name

== Output options ==
   --compact [FLG: OFF]          : Skip writing each tile
   --out     [STR: ]             : Output TSV file

== Filtering options ==
   --zoom    [INT: -1]           : Zoom level (default: -1 -- maximum zoom level)

== Performance options ==
   --threads [INT: 0]            : Number of threads (default: hardware concurrency)


NOTES:
When --help was included in the argument. The program prints the help message but do not actually run
```