cmake_minimum_required(VERSION 2.8...3.23)
project(pmpoint)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
#set(CMAKE_CXX_EXTENSIONS OFF)
if(CMAKE_COMPILER_IS_GNUCXX OR CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14")
elseif(MSVC)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /std:c++14")
endif()
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -pthread -std=c++14")

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/bin")

find_path(HTS_INCLUDE_DIRS htslib/bgzf.h ../htslib/bgzf.h HINTS /usr/lib/x86_64-linux-gnu/ /usr/lib/ /usr/lib64/ submodules/htslib ../htslib)
if(NOT HTS_INCLUDE_DIRS )
    message(FATAL_ERROR "libhts HTS_INCLUDE_DIRS not found")
endif()

message(STATUS ${HTS_INCLUDE_DIRS})

find_library(HTS_LIBRARIES libhts.a HINTS /usr/lib/x86_64-linux-gnu/ /usr/lib/ /usr/lib64/ submodules/htslib/ ../htslib/)
if(NOT HTS_LIBRARIES)
    message(FATAL_ERROR "libhts HTS_LIBRARIES not found")
endif()

message(STATUS ${HTS_LIBRARIES})

find_path(QGEN_INCLUDE_DIRS qgenlib/qgen_error.h  HINTS /usr/lib/x86_64-linux-gnu/ /usr/lib/ /usr/lib64/ submodules/qgenlib ../qgenlib)
if(NOT QGEN_INCLUDE_DIRS )
    message(FATAL_ERROR "libqgen QGEN_INCLUDE_DIRS not found")
endif()

message(STATUS ${QGEN_INCLUDE_DIRS})

find_library(QGEN_LIBRARIES libqgen.a  HINTS /usr/lib/x86_64-linux-gnu/ /usr/lib/ /usr/lib64/ submodules/qgenlib/lib ../qgenlib/lib/)
if(NOT QGEN_LIBRARIES)
    message(FATAL_ERROR "libqgen QGEN_LIBRARIES not found")
endif()

message(STATUS ${QGEN_LIBRARIES})

include_directories(${HTS_INCLUDE_DIRS} ${QGEN_INCLUDE_DIRS} ext)

set(SOURCE_FILES
    cmd_export_pmtiles.cpp
    cmd_summarize_pmtiles.cpp
    polygon.h
    pmt_utils.h
    pmt_utils.cpp
    mvt_pts.h
    mvt_pts.cpp
    pmt_pts.h
    pmt_pts.cpp
    pmpoint.cpp
    pmpoint.h
    ext/mapbox/geometry_io.hpp
    ext/mapbox/feature.hpp
    ext/mapbox/vector_tile.hpp
    ext/mapbox/variant_io.hpp
    ext/mapbox/recursive_wrapper.hpp
    ext/mapbox/vector_tile/version.hpp
    ext/mapbox/vector_tile/vector_tile_config.hpp
    ext/mapbox/variant.hpp
    ext/mapbox/variant_visitor.hpp
    ext/mapbox/optional.hpp
    ext/mapbox/geometry/line_string.hpp
    ext/mapbox/geometry/for_each_point.hpp
    ext/mapbox/geometry/envelope.hpp
    ext/mapbox/geometry/multi_polygon.hpp
    ext/mapbox/geometry/empty.hpp
    ext/mapbox/geometry/point.hpp
    ext/mapbox/geometry/polygon.hpp
    ext/mapbox/geometry/multi_point.hpp
    ext/mapbox/geometry/multi_line_string.hpp
    ext/mapbox/geometry/point_arithmetic.hpp
    ext/mapbox/geometry/geometry.hpp
    ext/mapbox/geometry/box.hpp
    ext/mapbox/geometry.hpp
    ext/protozero/pbf_reader.hpp
    ext/protozero/pbf_message.hpp
    ext/protozero/varint.hpp
    ext/protozero/pbf_writer.hpp
    ext/protozero/iterators.hpp
    ext/protozero/config.hpp
    ext/protozero/byteswap.hpp
    ext/protozero/version.hpp
    ext/protozero/pbf_builder.hpp
    ext/protozero/types.hpp
    ext/protozero/exception.hpp    
    ext/nlohmann/json.hpp
    ext/PMTiles/pmtiles.hpp
    )

SET(APP_EXE pmpoint)

#find_package(OpenMP REQUIRED)
ADD_EXECUTABLE(${APP_EXE} ${SOURCE_FILES})

find_library(ZLIB z HINTS /usr/lib/x86_64-linux-gnu/ /usr/lib/ /usr/lib64/)
if(NOT ZLIB)
    message(FATAL_ERROR "libz library not found")
endif()

find_library(CURLLIB curl HINTS /usr/lib/x86_64-linux-gnu/ /usr/lib/ /usr/lib64/)
if(NOT CURLLIB)
    message(FATAL_ERROR "libcurl library not found")
endif()

find_package (BZip2)
if (NOT BZIP2_FOUND)
    message(FATAL_ERROR "libbz2 library not found")
else()
    include_directories(${BZIP2_INCLUDE_DIRS})
    set(BZIP2 ${BZIP2_LIBRARIES})
endif()

find_library(LZMA lzma HINTS /usr/lib/x86_64-linux-gnu/ /usr/lib/ /usr/lib64/)
if(NOT LZMA)
    message(WARNING "liblzma library not found, if you specified --disable-lzma when compiling libhts, please ignore this warning")
    set(LZMA "")
else()
    set(LZMA ${LZMA_LIBRARIES})
endif()

find_library(DEFLATELIB deflate libdeflate.a HINTS /usr/lib/x86_64-linux-gnu/ /usr/lib/ /usr/lib64/ ../libdeflate/)
if(NOT DEFLATELIB)
    message(FATAL_ERROR "libdeflate library not found")
endif()

find_library(CRYPTOLIB crypto HINTS /usr/lib/x86_64-linux-gnu/ /usr/lib/ /usr/lib64/)
if(NOT CRYPTOLIB)
    message(FATAL_ERROR "libcrypto library not found")
endif()

find_package(Threads REQUIRED)

TARGET_LINK_LIBRARIES(${APP_EXE} PRIVATE ${QGEN_LIBRARIES} ${HTS_LIBRARIES} ${ZLIB} ${LZMA} ${BZIP2} ${CURLLIB} ${DEFLATELIB} ${CRYPTOLIB} Threads::Threads)
target_compile_options(${APP_EXE} PRIVATE)

install(TARGETS ${APP_EXE} RUNTIME DESTINATION bin)
