# pmpoint count-tiles

## Summary 

`pmpoint count-tiles` counts the number of points in each tile from a PMTiles file, optionally filtering by zoom level.

An example command is given below:

```bash
## Example command to count tiles a remote PMTiles file at zoom level 14
pmpoint count-tiles --in https://cartostore.s3.us-east-1.amazonaws.com/data/batch=2025_12/mouse-brain-test-collection/subdata_seqscope_n17t89b_c1d2f/genes_all.pmtiles --out-tsv out.count_tiles.z15.tsv.gz --zoom 15
```

## Required options

* `--in`: Input PMTiles file. The file can be either local file or a URL to a remote file (supports HTTP/HTTPS).
* `--out-tsv`: Output TSV file to store the tile counts.

## Additional Options

* `--zoom`: Zoom level to count tiles. Default is -1, which counts tiles from all zoom levels. If a specific zoom level is provided, only tiles from that zoom level will be counted.

## Expected Output

The output TSV file contains the following columns:
* `zoom`: Zoom level of the tile.
* `x`: X (integer) coordinate of the tile according to the PMTiles specification.
* `y`: Y (integer) coordinate of the tile according to the PMTiles specification.
* `tile_id` : A 64-bit unique identifier for the tile, according to the PMTiles specification.
* `tile_count`: Number of points in the tile.
* `total_count`: (Available when counting all zoom levels) Total number of points from all its descendant tiles without subsampling.
* `frac_in_tile`: (Available when counting all zoom levels) Fraction of points in the tile compared to the `total_count`.

## Full Usage 

The full usage of `pmpoint count-tiles` can be viewed with the `--help` option:

```
$ ./pmpoint count-tiles --help
[bin/pmpoint count-tiles] -- Count number of points in each tiles from a PMTiles file

 Copyright (c) 2022-2025 by Hyun Min Kang
 Licensed under the Apache License v2.0 http://www.apache.org/licenses/

Detailed instructions of parameters are available. Ones with "[]" are in effect:

Available Options:

== Input options ==
   --in      [STR: ]             : Input PMTiles file

== Output options ==
   --out-tsv [STR: ]             : Output TSV file

== Filtering options ==
   --zoom    [INT: -1]           : Zoom level (default: -1 -- all zoom levels)


NOTES:
When --help was included in the argument. The program prints the help message but do not actually run
```