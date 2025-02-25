# pmpoint : A C++ toolkit for PMTiles point

## Overview 

`pmpoint` is a collection of C++ tools that are designed to help analyze sub-micron resolution
spatial transcriptomics data stored in PMTiles format.
These tools are under active development, so they may change frequently. 

## Installation

You can install `pmpoint` by following the instructions below:

```bash
## clone the repository
git clone --recursive https://github.com/seqscope/pmpoint.git
cd pmpoint

## build the submodules
cd submodules
sh -x build.sh
cd ..

## build spatula
mkdir build
cd build
cmake ..
make

## list available package
../bin/pmpoint --help
```
