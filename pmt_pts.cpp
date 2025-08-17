#include <cstdio>
#include <string>
#include <zlib.h>

#include "pmt_utils.h"
#include "pmt_pts.h"
#include "mvt_pts.h"
#include "ext/PMTiles/pmtiles.hpp"
#include "ext/nlohmann/json.hpp"
#include "qgenlib/qgen_error.h"

void pmt_pts::init()
{
  // define a lambda function to decompress the data, currently, only gzip is supported
  decompress_func = [](const std::string &compressed_data, uint8_t compression_method) -> std::string
  {
    if (compression_method == pmtiles::COMPRESSION_GZIP)
    {
      z_stream inflate_s;
      std::string decompressed_data;
      inflate_s.zalloc = Z_NULL;
      inflate_s.zfree = Z_NULL;
      inflate_s.opaque = Z_NULL;
      inflate_s.avail_in = 0;
      inflate_s.next_in = Z_NULL;
      if (inflateInit2(&inflate_s, 32 + 15) != Z_OK)
      {
        fprintf(stderr, "Decompression error: %s\n", inflate_s.msg);
      }
      inflate_s.next_in = (Bytef *)compressed_data.data();
      inflate_s.avail_in = compressed_data.size();
      inflate_s.next_out = (Bytef *)decompressed_data.data();
      inflate_s.avail_out = decompressed_data.size();

      while (true)
      {
        size_t existing_output = inflate_s.next_out - (Bytef *)decompressed_data.data();

        decompressed_data.resize(existing_output + 2 * inflate_s.avail_in + 100);
        inflate_s.next_out = (Bytef *)decompressed_data.data() + existing_output;
        inflate_s.avail_out = decompressed_data.size() - existing_output;

        int ret = inflate(&inflate_s, 0);
        if (ret < 0)
        {
          fprintf(stderr, "Decompression error: ");
          if (ret == Z_DATA_ERROR)
          {
            fprintf(stderr, "data error");
          }
          if (ret == Z_STREAM_ERROR)
          {
            fprintf(stderr, "stream error");
          }
          if (ret == Z_MEM_ERROR)
          {
            fprintf(stderr, "out of memory");
          }
          if (ret == Z_BUF_ERROR)
          {
            fprintf(stderr, "no data in buffer");
          }
          fprintf(stderr, "\n");
          return "";
        }

        if (ret == Z_STREAM_END)
        {
          break;
        }

        // ret must be Z_OK or Z_NEED_DICT;
        // continue decompresing
      }

      decompressed_data.resize(inflate_s.next_out - (Bytef *)decompressed_data.data());
      inflateEnd(&inflate_s);
      return decompressed_data;
    }
    else
    {
      error("Unsupported compression method");
      return "";
    }
  };
}

bool pmt_pts::open(const char *fname)
{
  std::lock_guard<std::mutex> lock(mtx);
  fp = fopen(fname, "rb");
  if (fp == NULL)
  {
    return false;
  }
  return true;
}

void pmt_pts::close()
{
  std::lock_guard<std::mutex> lock(mtx);
  if (fp != NULL)
  {
    fclose(fp);
    fp = NULL;
  }
}

void pmt_pts::print_header_info(FILE *fp)
{
  const char *comp_strs[5] = {"unknown", "none", "gzip", "brotli", "zstd"};
  const char *ttype_strs[6] = {"unknown", "MVT", "PNG", "JPEG", "WebP", "AVIF"};
  fprintf(fp, "-----------------------------------------\n");
  fprintf(fp, "Root Directory Offset : %llu\n", hdr.root_dir_offset);
  fprintf(fp, "Root Directory Length : %llu\n", hdr.root_dir_bytes);
  fprintf(fp, "Metadata Offset : %llu\n", hdr.json_metadata_offset);
  fprintf(fp, "Metadata Length : %llu\n", hdr.json_metadata_bytes);
  fprintf(fp, "Leaf Directory Offset: %llu\n", hdr.leaf_dirs_offset);
  fprintf(fp, "Leaf Directory Length: %llu\n", hdr.leaf_dirs_bytes);
  fprintf(fp, "Tile Data Offset: %llu\n", hdr.tile_data_offset);
  fprintf(fp, "Tile Data Length: %llu\n", hdr.tile_data_bytes);
  fprintf(fp, "Num of Addressed Tiles: %llu\n", hdr.addressed_tiles_count);
  fprintf(fp, "Number of Tile Entries: %llu\n", hdr.tile_entries_count);
  fprintf(fp, "Number of Tile Contents: %llu\n", hdr.tile_contents_count);
  fprintf(fp, "Clustered: %s\n", hdr.clustered ? "TRUE" : "FALSE");
  fprintf(fp, "Internal Compression: %s\n", comp_strs[hdr.internal_compression]);
  fprintf(fp, "Tile Compression: %s\n", comp_strs[hdr.tile_compression]);
  fprintf(fp, "Tile Type: %s\n", ttype_strs[hdr.tile_type]);
  fprintf(fp, "Min Zoom: %u\n", hdr.min_zoom);
  fprintf(fp, "Max Zoom: %u\n", hdr.max_zoom);
  fprintf(fp, "Min X (x10M) : %d\n", hdr.min_lon_e7);
  fprintf(fp, "Min Y (x10M) : %d\n", hdr.min_lat_e7);
  fprintf(fp, "Max X (x10M) : %d\n", hdr.max_lon_e7);
  fprintf(fp, "Max Y (x10M) : %d\n", hdr.max_lat_e7);
  fprintf(fp, "Center Zoom : %u\n", hdr.center_zoom);
  fprintf(fp, "Center X (x10M) : %d\n", hdr.center_lon_e7);
  fprintf(fp, "Center Y (x10M) : %d\n", hdr.center_lat_e7);
  fprintf(fp, "-----------------------------------------\n");
}

void pmt_pts::print_metadata(FILE *fp)
{
  if (jmeta.empty())
  {
    return;
  }
  fprintf(fp, "Metadata:\n");
  fprintf(fp, "-----------------------------------------\n");
  fprintf(fp, "%s\n", jmeta.dump(2).c_str());
  fprintf(fp, "-----------------------------------------\n");
}

void pmt_pts::print_tile_entries(FILE *fp, uint8_t min_zoom, uint8_t max_zoom)
{
  for (int32_t i = 0; i < tile_entries.size(); i++)
  {
    const pmtiles::entry_zxy &e = tile_entries[i];
    if (e.z < min_zoom || e.z > max_zoom)
    {
      continue;
    }
    uint64_t tile_id = pmtiles::zxy_to_tileid(e.z, e.x, e.y);
    fprintf(fp, "Index: %d, TileID: %llu (Z: %u, X: %u, Y: %u) -- Offset: %llu, Length: %u\n", i, tile_id, e.z, e.x, e.y, e.offset, e.length);
  }
}

bool pmt_pts::read_header_meta_entries()
{
  std::lock_guard<std::mutex> lock(mtx);
  if (fp == NULL)
  {
    return false;
  }
  if (hdr_read)
  {
    return true;
  }
  // the header is 127 bytes field that is at the start of the archive.
  if (cur_pos != 0)
  {
    // reset the current position to the start of the file
    fseek(fp, 0, SEEK_SET);
    cur_pos = 0;
  }
  char hdr_bytes[128];
  if (fread(hdr_bytes, 1, 127, fp) != 127)
  {
    return false;
  }
  cur_pos += 127;
  // create a string
  std::string hdr_str(hdr_bytes, 127);

  hdr = pmtiles::deserialize_header(hdr_str);

  if (hdr.tile_type != pmtiles::TILETYPE_MVT)
  {
    //error("Unsupported tile type: %d -- Only MVTiles allowed", hdr.tile_type);
    notice("Unsupported tile type: %d", hdr.tile_type);
  }

  // print_header_info(stdout);

  char *buffer_before_tiles = new char[hdr.tile_data_offset];
  memcpy(buffer_before_tiles, hdr_bytes, 127);
  if (fread(buffer_before_tiles + 127, 1, hdr.tile_data_offset - 127, fp) != hdr.tile_data_offset - 127)
  {
    return false;
  }
  cur_pos = hdr.tile_data_offset;

  tile_entries = pmtiles::entries_tms(decompress_func, buffer_before_tiles);

  // read metadata information
  std::string meta_compressed(buffer_before_tiles + hdr.json_metadata_offset, hdr.json_metadata_bytes);
  std::string meta_decompressed = decompress_func(meta_compressed, hdr.internal_compression);

  // load the metadata into a json object
  jmeta = nlohmann::json::parse(meta_decompressed);

  delete[] buffer_before_tiles;

  // build a map of tile entries
  notice("A total of %zu tile entries are found", tile_entries.size());
  for (int32_t i = 0; i < tile_entries.size(); i++)
  {
    const pmtiles::entry_zxy &e = tile_entries[i];

    uint64_t tile_id = pmtiles::zxy_to_tileid(e.z, e.x, e.y);
    tileid2idx[tile_id] = i;
  }
  hdr_read = true;
  return hdr.tile_type == pmtiles::TILETYPE_MVT;
}

bool pmt_pts::read_metadata()
{
  std::lock_guard<std::mutex> lock(mtx);
  if (fp == NULL)
  {
    return false;
  }
  if (meta_read)
  {
    return true;
  }
  if (hdr.json_metadata_bytes == 0)
  {
    return false;
  }
  // move to the metadata offset
  if (cur_pos != hdr.json_metadata_offset)
  {
    fseek(fp, hdr.json_metadata_offset, SEEK_SET);
    cur_pos = hdr.json_metadata_offset;
  }
  // read the metadata
  char *p = new char[hdr.json_metadata_bytes];
  if (fread(p, 1, hdr.json_metadata_bytes, fp) != hdr.json_metadata_bytes)
  {
    delete[] p;
    return false;
  }
  // p[hdr.json_metadata_bytes] = '\0';
  cur_pos += hdr.json_metadata_bytes;

  std::string meta_str(p, hdr.json_metadata_bytes);
  std::string decomp_meta_str = decompress_func(meta_str, hdr.tile_compression);

  // load the metadata into a json object
  jmeta = nlohmann::json::parse(decomp_meta_str);

  delete[] p;

  meta_read = true;
  return true;
}

// size_t pmt_pts::fetch_tile(uint8_t z, uint32_t x, uint32_t y)
// {
//   uint64_t tile_id = pmtiles::zxy_to_tileid(z, x, y);
//   if (tileid2idx.find(tile_id) == tileid2idx.end())
//   {
//     error("Tile %u/%lu/%lu not found", z, x, y);
//   }
//   const pmtiles::entry_zxy &e = tile_entries[tileid2idx[tile_id]];
//   // allocate memory for the tile data
//   char *tile_data = new char[e.length];
//   // move to the tile data offset
//   if (cur_pos != e.offset)
//   {
//     fseek(fp, e.offset, SEEK_SET);
//     cur_pos = e.offset;
//   }
//   // read the tile data
//   if (fread(tile_data, 1, e.length, fp) != e.length)
//   {
//     delete[] tile_data;
//     error("Failed to read tile data %u/%lu/%lu with length %", z, x, y);
//   }
//   // uncompress the tile data
//   std::string tile_comp_data_str(tile_data, e.length);
//   delete[] tile_data;
//   tile_data_str = decompress_func(tile_comp_data_str, hdr.tile_compression);

//   return tile_data_str.size();
// }

size_t pmt_pts::fetch_tile_to_buffer(uint8_t z, uint32_t x, uint32_t y, std::string& buffer)
{
  uint64_t tile_id = pmtiles::zxy_to_tileid(z, x, y);
  if (tileid2idx.find(tile_id) == tileid2idx.end())
  {
    error("Tile %u/%lu/%lu not found", z, x, y);
  }
  const pmtiles::entry_zxy &e = tile_entries[tileid2idx[tile_id]];
  // allocate memory for the tile data
  char *tile_data = new char[e.length];
  // move to the tile data offset
  {
    std::lock_guard<std::mutex> lock(mtx);
    if (cur_pos != e.offset)
    {
      fseek(fp, e.offset, SEEK_SET);
      cur_pos = e.offset;
    }
    // read the tile data
    if (fread(tile_data, 1, e.length, fp) != e.length)
    {
      delete[] tile_data;
      error("Failed to read tile data %u/%lu/%lu with length %", z, x, y);
    }
  }
  // uncompress the tile data
  std::string tile_comp_data_str(tile_data, e.length);
  delete[] tile_data;
  buffer = decompress_func(tile_comp_data_str, hdr.tile_compression);

  //return tile_data_str.size();
  return buffer.size();
}