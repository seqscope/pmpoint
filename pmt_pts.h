#ifndef __PMT_PTS__H
#define __PMT_PTS__H

#include "ext/PMTiles/pmtiles.hpp"
#include "ext/nlohmann/json.hpp"

// A class for handling a PMTiles file that contains a collection of points in MVT format
// Each point may have multiple attributes, so the data can be reformatted in rectangular dataframe
// This implementation relies on the PMTiles v3 specification

class pmt_pts
{
public:
    FILE *fp = NULL;                              // file pointer
    pmtiles::headerv3 hdr;                        // PMTiles v3 header
    nlohmann::json jmeta;                         // metadata as a JSON object
    uint64_t cur_pos = 0;                         // current offset of the file
    std::vector<pmtiles::entry_zxy> tile_entries; // list of tile entries
    std::map<uint64_t, uint32_t> tileid2idx;      // dictionary of the file entries based on the tile ID
    std::function<std::string(const std::string &, uint8_t)> decompress_func;
    std::string tile_data_str; // string to store uncompressed tile data

    bool hdr_read = false;  // flag to indicate if the header is read
    bool meta_read = false; // flag to indicate if the metadata is read

    void init();                     // initialize default parameters
    bool open(const char *fname);    // open a PMTiles file
    void close();                    // close a PMTiles file
    bool read_header_meta_entries(); // read the header of a PMTiles file
    bool read_metadata();            // read the metadata of a PMTiles file
    // bool get_tile_entries();      // get the tile entries of a PMTiles file
    size_t fetch_tile(uint8_t z, uint32_t x, uint32_t y); // fetch the uncompressed tile data of a PMTiles file
    size_t fetch_tile_to_buffer(uint8_t z, uint32_t x, uint32_t y, std::string& buffer); // fetch the uncompressed tile data of a PMTiles file
    // uint32_t parse_fetched_tile_as_mvt();                          // parse the fetched tile data as an MVTile object
    void print_header_info(FILE *fp);
    void print_metadata(FILE *fp);
    void print_tile_entries(FILE *fp, uint8_t min_zoom = 0, uint8_t max_zoom = 255);

    pmt_pts(const char *fname = NULL)
    {
        init();
        if (fname != NULL)
            open(fname);
    }
    ~pmt_pts() { close(); }
};
#endif