cmake_minimum_required(VERSION 2.6)

project(alShaders)

set(ALS_MAJOR_VERSION 0)
set(ALS_MINOR_VERSION 4)
set(ALS_PATCH_VERSION 0b10)
set(ALS_VERSION "${ALS_MAJOR_VERSION}.${ALS_MINOR_VERSION}.${ALS_PATCH_VERSION}")

set(CMAKE_BUILD_TYPE RELEASE)
set(CMAKE_VERBOSE_MAKEFILE FALSE)
set(CMAKE_SKIP_RPATH TRUE)
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -fvisibility=hidden")
set(CMAKE_CXX_FLAGS_DEBUG "-g -fvisibility=hidden")

# Fix to stop crash on OSX < 10.9
set(CMAKE_OSX_DEPLOYMENT_TARGET 10.8)

# pre-set the arnold and mtoa versions and root paths
set(ARNOLD_VERSION 4.1.3.2)
set(ARNOLD_ROOT /Users/anders/vfx/arnold/${ARNOLD_VERSION})

set(MTOA_VERSION 1.0.0)
set(MTOA_ROOT /Users/anders/vfx/mtoa/${MTOA_VERSION}_maya2014)

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
# fix for windows paths
# FIXME: ugly double-if because OR wants to spit a warning
if (WIN32)
link_directories(${ARNOLD_ROOT}/lib)
endif()
if (WIN64)
link_directories(${ARNOLD_ROOT}/lib)
endif()


set(MTOA_SHADERS ${MTOA_ROOT}/shaders)
set(MTOA_SCRIPTS ${MTOA_ROOT}/scripts)
set(MTOA_ICONS ${MTOA_ROOT}/icons)
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
        alInputVector
        alHair
        alColorSpace
        common
)

foreach(SUBDIR ${SUBDIRECTORIES})
    add_subdirectory(${SUBDIR})
endforeach()

configure_file(
    "${PROJECT_SOURCE_DIR}/package.in.py"
    "${PROJECT_SOURCE_DIR}/package.py"
)

configure_file(
    "${PROJECT_SOURCE_DIR}/test.in.py"
    "${PROJECT_SOURCE_DIR}/test.py"
)

