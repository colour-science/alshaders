cmake_minimum_required(VERSION 2.8)

project(alShaders)

set(ARNOLD_VERSION 4.0.10.2)
set(ARNOLD_ROOT /Users/anders/vfx/arnold/${ARNOLD_VERSION})
include_directories(${ARNOLD_ROOT}/include)
link_directories(${ARNOLD_ROOT}/bin)

set(MTOA_VERSION 0.20.0)
set(MTOA_ROOT /Users/anders/vfx/mtoa/${MTOA_VERSION})
set(MTOA_SHADERS ${MTOA_ROOT}/shaders)
set(MTOA_SCRIPTS ${MTOA_ROOT}/scripts)

include_directories(common)

add_subdirectory(alSurface)



