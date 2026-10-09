include(Catch)

add_executable(wfmat_tests
    test_geom.cpp
    test_prepare.cpp
    test_json.cpp
    test_svg.cpp
    test_validate.cpp
    test_canonical.cpp)
target_link_libraries(wfmat_tests PRIVATE wfmat wfmat_io Catch2::Catch2WithMain wfmat_options)
target_compile_definitions(wfmat_tests PRIVATE
    WFMAT_TEST_DATA="${CMAKE_CURRENT_SOURCE_DIR}/data"
    WFMAT_TEST_OUTPUT="${CMAKE_CURRENT_BINARY_DIR}/output")
catch_discover_tests(wfmat_tests)
