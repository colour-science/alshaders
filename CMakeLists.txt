cmake_minimum_required(VERSION 2.6)

project(alShaders)

set(ALS_MAJOR_VERSION 0)
set(ALS_MINOR_VERSION 4)
set(ALS_PATCH_VERSION 0b5)
set(ALS_VERSION "${ALS_MAJOR_VERSION}.${ALS_MINOR_VERSION}.${ALS_PATCH_VERSION}")

set(CMAKE_BUILD_TYPE RELEASE)
set(CMAKE_VERBOSE_MAKEFILE FALSE)
set(CMAKE_SKIP_RPATH TRUE)
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -fvisibility=hidden")
set(CMAKE_CXX_FLAGS_DEBUG "-g -fvisibility=hidden")

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


if(DEFINED MTOA_ROOT)
message(INFO " --- Will compile for MtoA")
set(MTOA_SHADERS ${MTOA_ROOT}/shaders)
set(MTOA_SCRIPTS ${MTOA_ROOT}/scripts)
set(MTOA_ICONS ${MTOA_ROOT}/icons)
set(MTOA_UI ${MTOA_ROOT}/scripts/mtoa/ui/ae)
else()
message(INFO " --- Will not compile for MtoA. To compile for Mtoa set the MTOA_ROOT variable in local.cmake.")
endif()

if(DEFINED HTOA_ROOT)
message(INFO " --- Will compile for HtoA")
set(HTOA_SHADERS ${HTOA_ROOT}/shaders)
else()
message(INFO " --- Will not compile for HtoA. To compile for HtoA set the HTOA_ROOT variable in local.cmake.")
endif()

include_directories(common)
include_directories(/usr/local/Cellar/ilmbase/1.0.3/include/)

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
	alTriplanar
        alSwitch
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

