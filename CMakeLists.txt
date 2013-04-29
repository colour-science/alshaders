cmake_minimum_required(VERSION 2.8)

project(alShaders)

set(CMAKE_BUILD_TYPE RELEASE)
set(CMAKE_VERBOSE_MAKEFILE TRUE)
set(CMAKE_SKIP_RPATH TRUE)
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -fvisibility=hidden")
set(CMAKE_CXX_FLAGS_DEBUG "-g -fvisibility=hidden")

# pre-set the arnold and mtoa versions and root paths
set(ARNOLD_VERSION 4.0.10.2)
set(ARNOLD_ROOT /Users/anders/vfx/arnold/${ARNOLD_VERSION})

set(MTOA_VERSION 0.22.1)
set(MTOA_ROOT /Users/anders/vfx/mtoa/${MTOA_VERSION})

#set boost root path
set(BOOST_ROOT /Users/anders/vfx/boost/1.43.0)
set(BOOST_INCLUDE_PATH ${BOOST_ROOT}/include)

# check if we have a local cmake include file and include it if we do
# this is useful for overriding our arnold and mtoa locations
if(EXISTS ${CMAKE_SOURCE_DIR}/local.cmake)
    message(INFO " --- Using local.cmake")
    include(${CMAKE_SOURCE_DIR}/local.cmake)
else()
    message(INFO " --- No local.cmake found")
endif()

include_directories(${ARNOLD_ROOT}/include)
link_directories(${ARNOLD_ROOT}/bin)

set(MTOA_SHADERS ${MTOA_ROOT}/shaders)
set(MTOA_SCRIPTS ${MTOA_ROOT}/scripts)
set(MTOA_UI ${MTOA_ROOT}/scripts/mtoa/ui/ae)

include_directories(common)

#include boost for alPhotometric
include_directories(${BOOST_INCLUDE_PATH})

set(SUBDIRECTORIES 
        alSurface
        alRemap
        alCombine
        alNoise
        alLayer
        alPattern
        alCurvature
        alBlackbody
        alPhotometric
        alInputVector
        alHair
        common
)

foreach(SUBDIR ${SUBDIRECTORIES})
    add_subdirectory(${SUBDIR})
endforeach()

