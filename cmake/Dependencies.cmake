# [LIB-01], [LIB-03] Dependencies: installed packages first (vcpkg manifest or system), FetchContent otherwise.
include(FetchContent)

macro(wfmat_find_or_fetch)
    if(NOT WFMAT_FETCH_DEPS)
        find_package(${ARGN} QUIET)
    endif()
endmacro()

# Boost (header-only parts: Geometry, Multiprecision, Math).
add_library(wfmat_boost INTERFACE)
wfmat_find_or_fetch(Boost 1.81 CONFIG)
if(TARGET Boost::headers)
    target_link_libraries(wfmat_boost INTERFACE Boost::headers)
else()
    wfmat_find_or_fetch(boost_geometry 1.81 CONFIG)
    wfmat_find_or_fetch(boost_multiprecision 1.81 CONFIG)
    wfmat_find_or_fetch(boost_math 1.81 CONFIG)
    if(NOT (TARGET Boost::geometry AND TARGET Boost::multiprecision AND TARGET Boost::math))
        message(STATUS "wfmat: fetching Boost")
        set(BOOST_INCLUDE_LIBRARIES geometry multiprecision math)
        set(BOOST_ENABLE_CMAKE ON)
        FetchContent_Declare(Boost
            URL https://github.com/boostorg/boost/releases/download/boost-1.84.0/boost-1.84.0.tar.xz
            URL_HASH SHA256=2e64e5d79a738d0fa6fb546c6e5c2bd28f88d268a2a080546f74e5ff98f29d0e
            DOWNLOAD_EXTRACT_TIMESTAMP ON
            EXCLUDE_FROM_ALL)
        FetchContent_MakeAvailable(Boost)
    endif()
    target_link_libraries(wfmat_boost INTERFACE Boost::geometry Boost::multiprecision Boost::math)
endif()

# tl::expected (CC0): an installed package, else the vendored single header in third_party/.
wfmat_find_or_fetch(tl-expected 1.1 CONFIG)
if(NOT TARGET tl::expected)
    add_library(wfmat_tl_expected INTERFACE)
    target_include_directories(wfmat_tl_expected SYSTEM INTERFACE
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/third_party/tl-expected/include>)
    add_library(tl::expected ALIAS wfmat_tl_expected)
endif()

# nlohmann/json (MIT).
wfmat_find_or_fetch(nlohmann_json 3.11 CONFIG)
if(NOT TARGET nlohmann_json::nlohmann_json)
    message(STATUS "wfmat: fetching nlohmann/json")
    FetchContent_Declare(nlohmann_json
        URL https://github.com/nlohmann/json/releases/download/v3.11.3/json.tar.xz
        URL_HASH SHA256=d6c65aca6b1ed68e7a182f4757257b107ae403032760ed6ef121c9d55e81757d
        DOWNLOAD_EXTRACT_TIMESTAMP ON
        EXCLUDE_FROM_ALL)
    FetchContent_MakeAvailable(nlohmann_json)
endif()

# Catch2 v3 (BSL-1.0), tests only.
if(WFMAT_BUILD_TESTS)
    wfmat_find_or_fetch(Catch2 3.4 CONFIG)
    if(NOT TARGET Catch2::Catch2WithMain)
        message(STATUS "wfmat: fetching Catch2")
        FetchContent_Declare(Catch2
            GIT_REPOSITORY https://github.com/catchorg/Catch2.git
            GIT_TAG v3.5.4
            GIT_SHALLOW ON
            EXCLUDE_FROM_ALL)
        FetchContent_MakeAvailable(Catch2)
        list(APPEND CMAKE_MODULE_PATH "${catch2_SOURCE_DIR}/extras")
    endif()
endif()
