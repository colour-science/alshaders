cmake_minimum_required(VERSION 2.8)

project(alShaders)

set(CMAKE_BUILD_TYPE RELEASE)
set(CMAKE_CXX_FLAGS_RELEASE "-O3")
set(CMAKE_CXX_FLAGS_DEBUG "-ggdb")

set(ARNOLD_VERSION 4.0.10.2)
set(ARNOLD_ROOT /usr/people/anders-l/tools/arnold/${ARNOLD_VERSION})
include_directories(${ARNOLD_ROOT}/include)
link_directories(${ARNOLD_ROOT}/bin)

set(MTOA_VERSION 0.20.0)
set(MTOA_ROOT /usr/people/anders-l/tools/maya/2012/mtoa/${MTOA_VERSION})
set(MTOA_SHADERS ${MTOA_ROOT}/shaders)
set(MTOA_SCRIPTS ${MTOA_ROOT}/scripts)

include_directories(common)

set(SUBDIRECTORIES 
		alSurface
		alRemap
		alCombine
		alNoise
		alLayer
		alPattern
		alCurvature
		common
)

foreach(SUBDIR ${SUBDIRECTORIES})
	add_subdirectory(${SUBDIR})
endforeach()

