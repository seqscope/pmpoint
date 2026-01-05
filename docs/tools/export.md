# pmpoint export

## Summary 

`pmpoint export` extracts pixel-level point data stored from a PMTiles files, allowing to spatially filter the data by bounding box or polygon.

An example command is given below:

```bash
pmpoint export --in https://cartostore.s3.us-east-1.amazonaws.com/data/batch=2025_12/mouse-brain-test-collection/subdata_seqscope_n17t89b_c1d2f/genes_all.pmtiles --out-tsv out.tsv.gz --zoom 15 --xmin 500 --xmax 600 --ymin 500 --ymax 600
```

## Required options

* `--in`: Input PMTiles file. The file can be either local file or a URL to a remote file (supports HTTP/HTTPS).
* `--out-tsv`: Output TSV file to store the extracted points.
* `--out-json`: Output JSON file to store the extracted points.

Either `--out-tsv` or `--out-json` should be provided.

## Additional Options

* `--zoom`: Zoom level to extract points. Default is -1, which extracts points from the highest zoom level.
* `--xmin`: Minimum x-axis value for filtering points. Default is -inf (no filtering).
* `--xmax`: Maximum x-axis value for filtering points. Default is inf (no filtering).
* `--ymin`: Minimum y-axis value for filtering points. Default is -inf (no filtering).
* `--ymax`: Maximum y-axis value for filtering points. Default is inf (no filtering).
* `--polygon`: GeoJSON file (in EPSG:3857) for polygon-based filtering. Points within the polygon will be extracted.
* `--precision`: Precision of the output of X/Y coordinates below decimal points (default: 3).

## Expected Output

The output TSV file contains all the extracted points in tabular format with columns for X, Y, gene/feature name, and count, and additional fields such as factor names (typically `[factor_name]_K1`) and the corresponding posterior probabilities (typically `[factor_name]_P1`) if available.

## Full Usage 

The full usage of `pmpoint export` can be viewed with the `--help` option:

```
$ .//pmpoint export --help
[bin/pmpoint export] -- Extract points from a PMTIles file

 Copyright (c) 2022-2025 by Hyun Min Kang
 Licensed under the Apache License v2.0 http://www.apache.org/licenses/

Detailed instructions of parameters are available. Ones with "[]" are in effect:

Available Options:

== Input options ==
   --in        [STR: ]             : Input PMTiles file

== Output options ==
   --out-tsv   [STR: ]             : Output TSV file
   --out-json  [STR: ]             : Output JSON file

== Filtering options ==
   --zoom      [INT: -1]           : Zoom level (default: -1 -- maximum zoom level)
   --xmin      [FLT: -inf]         : Minimum x-axis value
   --xmax      [FLT: inf]          : Maximum x-axis value
   --ymin      [FLT: -inf]         : Minimum y-axis value
   --ymax      [FLT: inf]          : Maximum y-axis value
   --polygon   [STR: ]             : GeoJSON file (in EPSG:3857) for polygon-based filtering

== Additional options ==
   --precision [INT: 3]            : Precision of the output of X/Y coordinates (default: 3)


NOTES:
When --help was included in the argument. The program prints the help message but do not actually run
```