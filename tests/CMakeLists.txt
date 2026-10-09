include(Catch)

add_executable(wfmat_tests
    test_geom.cpp
    test_prepare.cpp
    test_json.cpp
    test_svg.cpp
    test_validate.cpp
    test_canonical.cpp
    test_kernel.cpp
    test_engine.cpp
    test_arcs.cpp
    test_degenerate.cpp)
target_link_libraries(wfmat_tests PRIVATE wfmat wfmat_io Catch2::Catch2WithMain nlohmann_json::nlohmann_json wfmat_boost_polygon wfmat_options)
target_compile_definitions(wfmat_tests PRIVATE
    WFMAT_TEST_DATA="${CMAKE_CURRENT_SOURCE_DIR}/data"
    WFMAT_TEST_OUTPUT="${CMAKE_CURRENT_BINARY_DIR}/output"
    $<$<CXX_COMPILER_ID:MSVC>:_CRT_SECURE_NO_WARNINGS>)  # std::getenv for the WFMAT_* test knobs
catch_discover_tests(wfmat_tests)
