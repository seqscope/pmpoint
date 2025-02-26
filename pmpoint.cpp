#include "pmpoint.h"
#include "qgenlib/commands.h"
#include "qgenlib/qgen_utils.h"

int32_t cmd_summarize_pmtiles(int32_t argc, char **argv);
int32_t cmd_export_pmtiles(int32_t argc, char **argv);
int32_t cmd_count_tiles(int32_t argc, char **argv);

int32_t main(int32_t argc, char **argv)
{
  commandHelp.copyright_str = "Copyright (c) 2022-2025 by Hyun Min Kang";

  commandList cl;

  BEGIN_LONG_COMMANDS(longCommandlines)
  LONG_COMMAND_GROUP("Tools for PMTiles", NULL)
  LONG_COMMAND("summarize", &cmd_summarize_pmtiles, "Summarize a PMTIles file")
  LONG_COMMAND("export", &cmd_export_pmtiles, "Extract points from a PMTIles file")
  LONG_COMMAND("count-tiles", &cmd_count_tiles, "Count number of points in each tiles from a PMTiles file")
  END_LONG_COMMANDS();

  cl.Add(new longCommands("Available Commands", longCommandlines));

  if (argc < 2)
  {
    printf("[pmpoint] -- Utilities for handling PMTiles Points\n\n");
    fprintf(stderr, " Copyright (c) 2022-2024 by Hyun Min Kang\n");
    fprintf(stderr, " Licensed under the Apache License v2.0 http://www.apache.org/licenses/\n\n");
    fprintf(stderr, "To run a specific command      : %s [command] [options]\n", argv[0]);
    fprintf(stderr, "For detailed instructions, run : %s --help\n", argv[0]);
    cl.Status();
    return 1;
  }
  else
  {
    if (strcmp(argv[1], "--help") == 0)
    {
      cl.HelpMessage();
    }
    else
    {
      return cl.Read(argc, argv);
    }
  }
  return 0;
}
