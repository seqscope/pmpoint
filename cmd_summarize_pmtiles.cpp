#include "pmpoint.h"
#include "qgenlib/tsv_reader.h"
#include "qgenlib/qgen_error.h"

#include <vector>
#include <string>
#include <cstring>
#include <climits>
#include <map>

#include "pmt_pts.h"

/////////////////////////////////////////////////////////////////////////
// summarize : Summarize a PMTiles file
////////////////////////////////////////////////////////////////////////
int32_t cmd_summarize_pmtiles(int32_t argc, char** argv) {
  std::string pmtilesf;
  int32_t z = -1, x = -1, y = -1;
  bool show_hdr = false;
  bool show_meta = false;
  bool show_tile = false;
  bool show_all = false;

  paramList pl;

  BEGIN_LONG_PARAMS(longParameters)
    LONG_PARAM_GROUP("Input options", NULL)
    LONG_STRING_PARAM("in", &pmtilesf, "Input PMTiles file")

    LONG_PARAM_GROUP("Items to show", NULL)
    LONG_PARAM("all", &show_all, "Show all information (equivalent to --hdr --meta --tile)")
    LONG_PARAM("hdr", &show_hdr, "Show header information")
    LONG_PARAM("meta", &show_meta, "Show metadata information")
    LONG_PARAM("tile", &show_tile, "Show tile information")
  END_LONG_PARAMS();

  pl.Add(new longParams("Available Options", longParameters));
  pl.Read(argc, argv);
  pl.Status();

  notice("Analysis started");

  if ( pmtilesf.empty() ) {
    error("Missing required options --in");
  }  

  if ( show_all ) {
    show_hdr = show_meta = show_tile = true;
  }
  if ( !show_hdr && !show_meta && !show_tile ) {
    error("No items to show. Please specify --hdr, --meta, --tile, (or --all)");
  }

  // open a PMTiles file
  pmt_pts pmt(pmtilesf.c_str());

  notice("Reading header and tile entries...");
  pmt.read_header_meta_entries();
  
  if ( show_hdr ) 
    pmt.print_header_info(stdout);
  if ( show_meta ) {
    pmt.print_metadata(stdout);
  }
  if ( show_tile ) 
    pmt.print_tile_entries(stdout);

  notice("Analysis Finished");
  
  return 0;
}
