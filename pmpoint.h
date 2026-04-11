#ifndef __PMPOINT_H
#define __PMPOINT_H

#include "ext/nlohmann/json.hpp"
#include <getopt.h>
#include <cstdio>
#include <vector>
#include <algorithm>
#include <utility>
#include <string>
#include <set>
#include <ctime>
#include <cmath>
#include <cassert>

int32_t cmd_build_mlt_point_pmtiles(int32_t argc, char** argv);

#include "qgenlib/params.h"
#include "qgenlib/qgen_error.h"
#include "qgenlib/hts_utils.h"

#include "version.h"

#endif // __PMPOINT_H
