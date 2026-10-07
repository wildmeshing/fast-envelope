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
# Currently the pins of VolumeRemesher since its PR #28, both upstream's latest. NFG 9b7635a
# includes MarcoAttene/NFG#4, which counts _M_X64 as SSE2: before it, MSVC files built
# without /arch:AVX2 (this library's, in a downstream that strips the flag) compiled the
# scalar interval_number layout while VolumeRemesher's compiled the SIMD one, and the two
# store the same bounds in opposite order. The one fix the wildmeshing fork of
# Indirect_Predicates carried (LNC/BPT/TBC reporting an undecided denominator sign as usable)
# is MarcoAttene/Indirect_Predicates#15, merged. Both libraries define everything inside
# namespaces -- NFG and IPs -- which is what the using-directives in
# indirectPredicates/ip_filtered.h are for.
include(FetchContent)
# NFG comes from wildmeshing/NFG#1 (05f99ea) for now, on top of 9b7635a: it adds
# bignatural::trimMemoryPool(), which lets a caller return what the thread-local bignatural
# pool grew into during an exact computation. VolumeRemesher pins the same commit (its #30);
# back to MarcoAttene/NFG once the change is merged there.
FetchContent_Declare(nfg
    GIT_REPOSITORY https://github.com/wildmeshing/NFG.git
    GIT_TAG 05f99ea5584d8ba180ea35f878b0c02f2f888fd6
    SOURCE_SUBDIR do-not-configure
)
FetchContent_Declare(indirect_predicates
    GIT_REPOSITORY https://github.com/MarcoAttene/Indirect_Predicates.git
    GIT_TAG a5cb2c4e6a2e8ae137f23ffcee91d113e8f3a1c0
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
