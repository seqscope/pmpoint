# pmpoint summarize

## Summary 

`pmpoint summarize` provides a summary of the contents of a PMTiles file, including header information, metadata, and tile statistics.

An example command is given below:

```bash
## Example command to summarize a PMTiles file
pmpoint summarize --in https://cartostore.s3.us-east-1.amazonaws.com/data/batch=2025_12/mouse-brain-test-collection/subdata_seqscope_n17t89b_c1d2f/genes_all.pmtiles --all
```

## Required options

* `--in`: Input PMTiles file. This can be a local file or a URL to a remote file (supports HTTP/HTTPS).

## Additional Options

* `--all`: Show all information (equivalent to `--hdr`, `--meta`, and `--tile`).
* `--hdr`: Show header information.
* `--meta`: Show metadata information.
* `--tile`: Show tile information.

## Expected Output

With `--hdr` option, it displays the header information of the offsets of various sections in the PMTiles file, number of tiles, zoom levels, and the coordinate range covered.

With `--meta` option, it displays the metadata stored in the PMTiles file, which is typically in JSON format.

With `--tile` option, it provides basic statistics about the tiles in the PMTiles file, such as the offsets to the tile and the size of each tile.

## Full Usage 

The full usage of `pmpoint summerize` can be viewed with the `--help` option:

```
$ ./pmpoint summarize --help
[bin/pmpoint summarize] -- Summarize a PMTIles file

 Copyright (c) 2022-2025 by Hyun Min Kang
 Licensed under the Apache License v2.0 http://www.apache.org/licenses/

Detailed instructions of parameters are available. Ones with "[]" are in effect:

Available Options:

== Input options ==
   --in   [STR: ]             : Input PMTiles file

== Items to show ==
   --all  [FLG: OFF]          : Show all information (equivalent to --hdr --meta --tile)
   --hdr  [FLG: OFF]          : Show header information
   --meta [FLG: OFF]          : Show metadata information
   --tile [FLG: OFF]          : Show tile information


NOTES:
When --help was included in the argument. The program prints the help message but do not actually run
```