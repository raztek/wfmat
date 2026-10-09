cmake_minimum_required(VERSION 3.25)
project(wfmat VERSION 0.1.0 LANGUAGES CXX
        DESCRIPTION "Wavefront medial axis transform of regions bounded by lines and circular arcs")

option(WFMAT_BUILD_TESTS "Build the test suite" ${PROJECT_IS_TOP_LEVEL})
option(WFMAT_BUILD_TOOLS "Build wfmat-cli" ${PROJECT_IS_TOP_LEVEL})
option(WFMAT_WERROR "Treat compiler warnings as errors" OFF)
option(WFMAT_FETCH_DEPS "Always fetch dependencies instead of looking for installed ones" OFF)
option(WFMAT_BUILD_BENCHMARKS "Build the benchmarks (needs an installed Google Benchmark)" OFF)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
include(Dependencies)

# [LIB-03] Strict IEEE semantics: no contraction into FMA, never fast-math.
add_library(wfmat_options INTERFACE)
target_compile_features(wfmat_options INTERFACE cxx_std_20)
if(MSVC)
    target_compile_options(wfmat_options INTERFACE /fp:precise /W4 /permissive- /utf-8 $<$<BOOL:${WFMAT_WERROR}>:/WX>)
else()
    target_compile_options(wfmat_options INTERFACE -ffp-contract=off -fno-fast-math -Wall -Wextra -Wpedantic
                                                   $<$<BOOL:${WFMAT_WERROR}>:-Werror>)
endif()

add_library(wfmat
    src/geom.cpp
    src/cluster.cpp
    src/engine.cpp
    src/kernel.cpp
    src/mat.cpp
    src/region.cpp
    src/prepare.cpp
    src/validate.cpp)
add_library(wfmat::wfmat ALIAS wfmat)
target_include_directories(wfmat PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>)
target_link_libraries(wfmat PUBLIC tl::expected PRIVATE wfmat_boost wfmat_options)
target_compile_features(wfmat PUBLIC cxx_std_20)

add_library(wfmat_io
    src/io/json.cpp
    src/io/svg.cpp)
add_library(wfmat::io ALIAS wfmat_io)
target_link_libraries(wfmat_io PUBLIC wfmat PRIVATE nlohmann_json::nlohmann_json wfmat_options)

if(WFMAT_BUILD_TOOLS)
    add_executable(wfmat-cli tools/wfmat-cli.cpp)
    target_link_libraries(wfmat-cli PRIVATE wfmat_io wfmat_options)
endif()

if(WFMAT_BUILD_TESTS)
    enable_testing()
    add_subdirectory(tests)
endif()

# [VER-07] Benchmarks.
if(WFMAT_BUILD_BENCHMARKS)
    find_package(benchmark REQUIRED)
    add_executable(wfmat-bench bench/bench_mat.cpp)
    target_link_libraries(wfmat-bench PRIVATE wfmat benchmark::benchmark wfmat_options)
endif()
