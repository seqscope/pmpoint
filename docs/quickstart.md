# Quickstart for pmpoint

## Installing pmpoint

Please follow the instruction below to install `pmpoint`

```sh
git clone --recursive https://github.com/seqscope/pmpoint.git
cd pmpoint
cd submodules
sh -x build.sh
cd ..
mkdir build
cd build
cmake ..
make
```

## List available tools

To list the available tools, run the following command:

```sh
../bin/pmpoint --help
```

If you encounter any difficulties, see [Install](install.md) for more details.