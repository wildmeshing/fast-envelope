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
# Indirect_Predicates' own CMakeLists builds a test executable we do not want.
#
# These declarations must stay IDENTICAL to VolumeRemesher's (its root CMakeLists.txt).
# FetchContent keeps the first declaration it sees and silently ignores later ones, so in a
# downstream that uses both -- wildmeshing-toolkit does -- whichever is declared first
# decides the version for everyone. Identical declarations make that harmless; differing
# ones make the include order load-bearing and mean one of the two projects is built and
# tested against predicates it does not get in production.
#
# Currently the wildmeshing fork: upstream 797c07c plus one fix, while
# MarcoAttene/Indirect_Predicates#15 is open. LNC/BPT/TBC reported a homogeneous
# denominator whose sign the interval filter could not decide as usable (getIntervalLambda
# returned an unconditional true, and the constructor never normalized), where SSI/LPI/TPI
# do both. The indirect predicates return sgn(det) of a determinant scaled by that
# denominator, so the answer came back as the true orientation times an unknown sign.
include(FetchContent)
FetchContent_Declare(nfg
    GIT_REPOSITORY https://github.com/MarcoAttene/nfg.git
    GIT_TAG f1df171345e91eeccea3c8b1e3f703a9fe27fca3
    SOURCE_SUBDIR do-not-configure
)
FetchContent_Declare(indirect_predicates
    GIT_REPOSITORY https://github.com/wildmeshing/Indirect_Predicates.git
    GIT_TAG 013c7c1c3315c8dd56cbf32dce7bc3041e778372
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
