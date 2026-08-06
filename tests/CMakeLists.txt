# ###############################################################################
# Tests
# ###############################################################################

set(test_sources
  box_cut_equivalence.cpp
  orient_sign_convention.cpp
  run_test.cpp
  utils/csv_reader.cpp
  utils/csv_reader.h
)


add_executable(envelope_unit_tests ${test_sources})
  # set_property(TARGET envelope_unit_tests PROPERTY COMPILE_WARNING_AS_ERROR ON)

################################################################################
# Required Libraries
################################################################################

target_link_libraries(envelope_unit_tests PUBLIC FastEnvelope::FastEnvelope)

include(catch2)
target_link_libraries(envelope_unit_tests PUBLIC Catch2::Catch2WithMain)


target_compile_definitions(
  envelope_unit_tests
  PUBLIC ENVELOPE_TEST_DATA_DIR="${CMAKE_CURRENT_SOURCE_DIR}/data/"
)

################################################################################
# Compiler options
################################################################################

target_compile_definitions(envelope_unit_tests PUBLIC CATCH_CONFIG_ENABLE_BENCHMARKING)

################################################################################
# Register tests
################################################################################

FetchContent_GetProperties(catch2)
list(APPEND CMAKE_MODULE_PATH ${catch2_SOURCE_DIR}/extras)
include(Catch)

# Register tests
set(PARSE_CATCH_TESTS_ADD_TO_CONFIGURE_DEPENDS ON)
catch_discover_tests(envelope_unit_tests ADD_TAGS_AS_LABELS EXTRA_ARGS --rng-seed 1)
