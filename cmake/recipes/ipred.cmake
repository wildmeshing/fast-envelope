if(TARGET indirectPredicates)
    return()
endif()

message(STATUS "Third-party: creating target 'indirectPredicates'")

# The exact/indirect predicates and the number types used to come from
# teseoch/Indirect_Predicates, which vendored numerics.h and built it as a library.
# Take them from MarcoAttene's upstream instead -- NFG (numbers) plus
# Indirect_Predicates (predicates) -- because VolumeRemesher moved to the same pair
# and the two cannot coexist in one binary.
#
# Both forks put `expansionObject` in the GLOBAL namespace, and they disagree on
# whether its members are static: teseoch's are non-static, MarcoAttene's are static.
# The Itanium ABI mangles the two identically, so a binary linking both keeps one
# definition and every call through the other convention has its arguments shifted by
# one register -- an immediate segfault. Sharing one copy is the only way for a
# downstream project to use both libraries.
#
# SOURCE_SUBDIR points at a directory with no CMakeLists.txt so MakeAvailable populates
# the sources without add_subdirectory()ing them: both are header-only here, and
# Indirect_Predicates' own CMakeLists builds a test executable we do not want. The
# declarations match VolumeRemesher's, so whichever project declares them first wins
# and the build ends up with exactly one copy.
include(FetchContent)
FetchContent_Declare(nfg
    GIT_REPOSITORY https://github.com/MarcoAttene/nfg.git
    GIT_TAG f1df171345e91eeccea3c8b1e3f703a9fe27fca3
    SOURCE_SUBDIR do-not-configure
)
FetchContent_Declare(indirect_predicates
    GIT_REPOSITORY https://github.com/MarcoAttene/Indirect_Predicates.git
    GIT_TAG 8d354b9affed68332d99c2a361e240543e0dba45
    SOURCE_SUBDIR do-not-configure
)
FetchContent_MakeAvailable(nfg indirect_predicates)

add_library(indirectPredicates INTERFACE)
target_include_directories(indirectPredicates INTERFACE
    ${nfg_SOURCE_DIR}/include
    ${indirect_predicates_SOURCE_DIR}/include
)
target_compile_features(indirectPredicates INTERFACE cxx_std_20)

# On ARM the exact predicates rely on SIMDe to emulate the x86 AVX2/FMA intrinsics
# NFG's numerics.h includes directly under __ARM_NEON.
if(CMAKE_SYSTEM_PROCESSOR MATCHES "^(arm64|aarch64|ARM64)$")
    FetchContent_Declare(simde
        GIT_REPOSITORY https://github.com/simd-everywhere/simde.git
        GIT_TAG v0.8.2
    )
    FetchContent_MakeAvailable(simde)
    target_include_directories(indirectPredicates INTERFACE ${simde_SOURCE_DIR}/simde)
endif()
