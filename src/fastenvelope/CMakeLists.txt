set(SOURCES
        AABB.h
        AABB.cpp
        FastEnvelope.cpp
        FastEnvelope.h
        FastEnvelope2D.cpp
        FastEnvelope2D.h
        Morton.h
        Morton.cpp
        Types.hpp
        common_algorithms.cpp
        common_algorithms.h

        indirectPredicates/ip_filtered_ex.h
        indirectPredicates/ip_filtered_ex.cpp
        indirectPredicates/ip_filtered.h
)

if(FAST_ENVELOPE_WITH_GMP)
        SET(SOURCES ${SOURCES} Rational.hpp)
endif()

prepend_current_path(SOURCES)
fast_envelope_set_source_group(${SOURCES})
target_sources(FastEnvelope PRIVATE ${SOURCES})
