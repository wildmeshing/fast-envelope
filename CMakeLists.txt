# Detects whether this is a top-level project
get_directory_property(HAS_PARENT PARENT_DIRECTORY)

if(HAS_PARENT)
    set(FAST_ENVELOPE_TOPLEVEL_PROJECT OFF)
else()
    set(FAST_ENVELOPE_TOPLEVEL_PROJECT ON)
endif()

# Check required CMake version
set(REQUIRED_CMAKE_VERSION "3.25.0")
if(FAST_ENVELOPE_TOPLEVEL_PROJECT)
    cmake_minimum_required(VERSION ${REQUIRED_CMAKE_VERSION})
    if(CMAKE_VERSION VERSION_GREATER_EQUAL "4.0.0")
        message(WARNING "CMake version is ${CMAKE_VERSION}, setting stuff for dependencies.")
        SET(CMAKE_POLICY_VERSION_MINIMUM ${REQUIRED_CMAKE_VERSION})
    endif()
else()
    # Don't use cmake_minimum_required here to avoid implicitly overriding parent policies
    if(${CMAKE_VERSION} VERSION_LESS ${REQUIRED_CMAKE_VERSION})
        message(FATAL_ERROR "CMake required version to build FastEnvelope is ${REQUIRED_CMAKE_VERSION}")
    endif()
endif()

# Enable ccache if available
find_program(CCACHE_PROGRAM ccache)
if(CCACHE_PROGRAM)
    option(WMTK_WITH_CCACHE "Enable ccache when building WMTK" ${FAST_ENVELOPE_TOPLEVEL_PROJECT})
else()
    option(WMTK_WITH_CCACHE "Enable ccache when building WMTK" OFF)
endif()
if(WMTK_WITH_CCACHE AND CCACHE_PROGRAM)
    set(ccacheEnv
        CCACHE_BASEDIR=${CMAKE_BINARY_DIR}
        CCACHE_SLOPPINESS=clang_index_store,include_file_ctime,include_file_mtime,locale,pch_defines,time_macros
    )
    foreach(lang IN ITEMS C CXX)
        set(CMAKE_${lang}_COMPILER_LAUNCHER
            ${CMAKE_COMMAND} -E env ${ccacheEnv} ${CCACHE_PROGRAM}
        )
    endforeach()
endif()

################################################################################
# CMake Policies
################################################################################

cmake_policy(SET CMP0054 NEW) # Only interpret if() arguments as variables or keywords when unquoted.
cmake_policy(SET CMP0076 NEW) # target_sources() command converts relative paths to absolute.
if(CMAKE_VERSION VERSION_GREATER_EQUAL "3.24")
  cmake_policy(SET CMP0135 NEW) # Set the timestamps of all extracted contents to the time of the extraction.
endif()

# ###############################################################################
project(FastEnvelope)

# ###############################################################################

if(INPUT_THIRD_PARTY_DIR)
	set(THIRD_PARTY_DIR ${CMAKE_CURRENT_SOURCE_DIR}/${INPUT_THIRD_PARTY_DIR}/)
else()
	set(THIRD_PARTY_DIR ${CMAKE_CURRENT_SOURCE_DIR}/3rdparty/)
endif()

list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake/recipes)
list(APPEND CMAKE_MODULE_PATH ${THIRD_PARTY_DIR}/Catch2/contrib)

# Color output
include(UseColors)

# Prepend function
include(PrependCurrentPath)

# FastEnvelope utils
include(FastEnvelopeUtils)

# Extra warnings
include(Warnings)

include(fast_envelope_copy_dll)

# Sort projects inside the solution
set_property(GLOBAL PROPERTY USE_FOLDERS ON)

# Generate position independent code by default
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

################################################################################


################################################################################

# FastEnvelope options
option(FAST_ENVELOPE_WITH_UNIT_TESTS   			"unit test project" 			ON)
option(FAST_ENVELOPE_ENABLE_TBB        			"Enable TBB"                    ON)
# option(FAST_ENVELOPE_USE_FLOAT          		"Use floats instead of double"  OFF)

option(FAST_ENVELOPE_WITH_GEOGRAM_PREDICATES        "Use whole Geogram for predicates"        ON)
option(FAST_ENVELOPE_WITH_GEOGRAM_PSM_PREDICATES    "Use Geogram predicates only"             OFF)
option(FAST_ENVELOPE_WITH_GMP						"Use gmp to have rational predicates"	  OFF)
option(FAST_ENVELOPE_WITH_TIMER   					"use timer in the library for debug" 	  OFF)

# Sanitizer options
option(SANITIZE_ADDRESS   "Sanitize Address"       OFF)
option(SANITIZE_MEMORY    "Sanitize Memory"        OFF)
option(SANITIZE_THREAD    "Sanitize Thread"        OFF)
option(SANITIZE_UNDEFINED "Sanitize Undefined"     OFF)

# Options for libigl modules
option(LIBIGL_USE_STATIC_LIBRARY          "Use libigl as static library" OFF)
option(LIBIGL_WITH_ANTTWEAKBAR            "Use AntTweakBar"              OFF)
option(LIBIGL_WITH_CGAL                   "Use CGAL"                     OFF)
option(LIBIGL_WITH_COMISO                 "Use CoMiso"                   OFF)
option(LIBIGL_WITH_CORK                   "Use Cork"                     OFF)
option(LIBIGL_WITH_EMBREE                 "Use Embree"                   OFF)
option(LIBIGL_WITH_LIM                    "Use LIM"                      OFF)
option(LIBIGL_WITH_MATLAB                 "Use Matlab"                   OFF)
option(LIBIGL_WITH_MOSEK                  "Use MOSEK"                    OFF)
option(LIBIGL_WITH_OPENGL                 "Use OpenGL"                   OFF)
option(LIBIGL_WITH_OPENGL_GLFW            "Use GLFW"                     OFF)
option(LIBIGL_WITH_OPENGL_GLFW_IMGUI      "Use ImGui"                    OFF)
option(LIBIGL_WITH_PNG                    "Use PNG"                      OFF)
option(LIBIGL_WITH_PYTHON                 "Use Python"                   OFF)
option(LIBIGL_WITH_TETGEN                 "Use Tetgen"                   OFF)
option(LIBIGL_WITH_TRIANGLE               "Use Triangle"                 OFF)
option(LIBIGL_WITH_VIEWER                 "Use OpenGL viewer"            OFF)
option(LIBIGL_WITH_XML                    "Use XML"                      OFF)

if(NOT ${FAST_ENVELOPE_WITH_GEOGRAM_PREDICATES} AND NOT ${FAST_ENVELOPE_WITH_GEOGRAM_PSM_PREDICATES})
	set(LIBIGL_WITH_PREDICATES ON)
else()
	set(LIBIGL_WITH_PREDICATES OFF)
endif()



# Sanitizers
if(FAST_ENVELOPE_WITH_SANITIZERS)
	list(APPEND CMAKE_MODULE_PATH ${THIRD_PARTY_DIR}/sanitizers-cmake/cmake)
endif()

# Setup dependencies
include(FastEnvelopeDependencies)

################################################################################
# FastEnvelope library
################################################################################

# add_library() can only be called without any source since CMake 3.11 ...
add_library(${PROJECT_NAME} src/common_algorithms.cpp)
# C++ standard
target_compile_features(${PROJECT_NAME} PUBLIC cxx_std_17)

# Public include directory for FastEnvelope
target_include_directories(${PROJECT_NAME} PUBLIC ${PROJECT_BINARY_DIR}/include)
if(MSVC)
	# Fix Windows issue with ambiguity of std::byte
	# target_compile_definitions(${PROJECT_NAME} PUBLIC _HAS_STD_BYTE=0)
	target_compile_definitions(${PROJECT_NAME} PUBLIC WIN32_LEAN_AND_MEAN)
endif()
# Extra warnings
target_link_libraries(${PROJECT_NAME} PRIVATE warnings::all)


# libigl
if(FAST_ENVELOPE_TOPLEVEL_PROJECT)
	include(libigl)
	target_link_libraries(${PROJECT_NAME} PUBLIC igl::core)
	if(LIBIGL_WITH_PREDICATES)
		target_link_libraries(${PROJECT_NAME} PUBLIC igl::predicates)
	endif()
else()
	if(LIBIGL_WITH_PREDICATES)
		include(libigl)
		target_link_libraries(${PROJECT_NAME} PUBLIC igl::core)
		target_link_libraries(${PROJECT_NAME} PUBLIC igl::predicates)
	else()
		# Assume eigen exists in the other lib
		target_link_libraries(${PROJECT_NAME} PUBLIC Eigen3::Eigen)
	endif()
endif()

if(FAST_ENVELOPE_WITH_SANITIZERS)
	add_sanitizers(${PROJECT_NAME})
endif()


################################################################################
# Required libraries
################################################################################


# TBB
if(FAST_ENVELOPE_ENABLE_TBB)
	include(onetbb)
endif()


# if(${FAST_ENVELOPE_USE_FLOAT})
# 	target_compile_definitions(${PROJECT_NAME} PUBLIC -DFAST_ENVELOPE_USE_FLOAT)
# endif()

include(cli11)
include(spdlog)

#spdlog
target_link_libraries(${PROJECT_NAME} PUBLIC spdlog::spdlog)
target_link_libraries(${PROJECT_NAME} PUBLIC CLI11::CLI11)

if(FAST_ENVELOPE_WITH_GEOGRAM_PSM_PREDICATES)
	# fast_envelope_download_geogram_predicates()

	add_library(GeoPredicates geogram_predicates/Predicates_psm.cpp)
	target_include_directories(GeoPredicates PUBLIC geogram_predicates/)
	target_compile_features(GeoPredicates PRIVATE ${CXX11_FEATURES})
	target_compile_definitions(GeoPredicates PRIVATE GEO_STATIC_LIBS)

	target_link_libraries(${PROJECT_NAME} PUBLIC GeoPredicates)
	target_compile_definitions(${PROJECT_NAME} PUBLIC ENVELOPE_WITH_GEO_PSM ENVELOPE_WITH_GEO)

elseif(FAST_ENVELOPE_WITH_GEOGRAM_PREDICATES)
	if(NOT TARGET geogram AND NOT TARGET geogram::geogram)
		# fast_envelope_download_geogram()
		include(geogram)
	endif()
	if(TARGET geogram)
		target_link_libraries(${PROJECT_NAME} PUBLIC geogram)
	elseif(TARGET geogram::geogram)
		target_link_libraries(${PROJECT_NAME} PUBLIC geogram::geogram)
	else()
		MESSAGE(FATAL_ERROR "Unable to get geogram or geogram::geogram targets")
	endif()
	target_compile_definitions(${PROJECT_NAME} PUBLIC ENVELOPE_WITH_GEO)
endif()



# Predicates
include("ipred")
target_link_libraries(${PROJECT_NAME} PRIVATE indirectPredicates)
if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
	message(WARNING "On CLang, there is no compiler flag which is required for ensuring the correctness of the algorithm.")
endif()
set(INDIRECTPREDICATES_HEADERS
	indirectPredicates/ip_filtered.h
	indirectPredicates/ip_filtered_ex.h
)
indirect_predicates_copy_headers(${INDIRECTPREDICATES_HEADERS})

# ip_filtered_ex.cpp holds the projection/normal helpers. It used to be pulled in by
# ip_filtered.cpp, which has been deleted, so compile it directly.
target_sources(${PROJECT_NAME} PRIVATE indirectPredicates/ip_filtered_ex.cpp)

if(FAST_ENVELOPE_WITH_TIMER)
	target_compile_definitions(${PROJECT_NAME} PRIVATE FAST_ENVELOPE_USE_TIMER)
endif()





################################################################################
# Subdirectories
################################################################################

# Sources
add_subdirectory(src)

if(FAST_ENVELOPE_WITH_GMP)
	find_package(GMP REQUIRED)

	target_link_libraries(${PROJECT_NAME} PUBLIC gmp::gmp)
	target_compile_definitions(${PROJECT_NAME} PUBLIC ENVELOPE_WITH_GMP)

endif()


################################################################################
# FastEnvelope binary
################################################################################

# Main executable
if(FAST_ENVELOPE_TOPLEVEL_PROJECT)
	# add an application library
	#add_subdirectory(app/benchmark)
	set(APP_SOURCES
		app/benchmark/MeshIO.hpp
		app/benchmark/MeshIO.cpp
		app/benchmark/mesh_AABB.h
		app/benchmark/mesh_AABB.cpp
		app/benchmark/AABBWrapper.cpp
		app/benchmark/AABBWrapper.h
		app/benchmark/csv_reader.h
		app/benchmark/csv_reader.cpp
		app/benchmark/getRSS.hpp
		app/benchmark/sampling.h
		app/benchmark/sampling.cpp
	)
	fast_envelope_copy_headers(${APP_SOURCES})
	add_library(${PROJECT_NAME}_app ${APP_SOURCES})
	target_include_directories(${PROJECT_NAME}_app PUBLIC ${PROJECT_BINARY_DIR}/include)
	target_include_directories(${PROJECT_NAME}_app PUBLIC app)

	# link libraries to the app library
	if(NOT FAST_ENVELOPE_WITH_GEOGRAM_PREDICATES)
		# Geogram
		# fast_envelope_download_geogram()
		include(geogram)
		target_link_libraries(${PROJECT_NAME}_app PUBLIC geogram)
	endif()
	target_link_libraries(${PROJECT_NAME}_app
		PUBLIC
			${PROJECT_NAME}

	)

	# executable
	add_executable(${PROJECT_NAME}_bin app/main.cpp)
	target_link_libraries(${PROJECT_NAME}_bin PUBLIC ${PROJECT_NAME}_app)

	target_compile_features(${PROJECT_NAME}_bin PUBLIC ${CXX14_FEATURES})
	if(FAST_ENVELOPE_WITH_SANITIZERS)
		add_sanitizers(${PROJECT_NAME}_bin)
	endif()


	if(FAST_ENVELOPE_ENABLE_TBB)
		target_link_libraries(${PROJECT_NAME}_bin PUBLIC TBB::tbb)
	endif()
	if(NOT (${CMAKE_VERSION} VERSION_LESS "3.6.0"))
		set_property(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY VS_STARTUP_PROJECT ${PROJECT_NAME}_bin)
	endif()
	fast_envelope_copy_dll(${PROJECT_NAME}_bin)

	# unit tests
	if(FAST_ENVELOPE_WITH_UNIT_TESTS)

		add_executable(${PROJECT_NAME}_test app/test.cpp)
		target_compile_definitions(${PROJECT_NAME}_test PUBLIC
    ENVELOPE_TEST_DATA_DIR="${CMAKE_CURRENT_SOURCE_DIR}/app/test_data/")
		target_link_libraries(${PROJECT_NAME}_test PUBLIC ${PROJECT_NAME}_app)

		target_compile_features(${PROJECT_NAME}_test PUBLIC ${CXX14_FEATURES})
		if(FAST_ENVELOPE_WITH_SANITIZERS)
			add_sanitizers(${PROJECT_NAME}_test)
		endif()


		if(FAST_ENVELOPE_ENABLE_TBB)
			target_link_libraries(${PROJECT_NAME}_test PUBLIC TBB::tbb)
		endif()
		if(NOT (${CMAKE_VERSION} VERSION_LESS "3.6.0"))
			set_property(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY VS_STARTUP_PROJECT ${PROJECT_NAME}_test)
		endif()

		fast_envelope_copy_dll(${PROJECT_NAME}_test)

		# Self-contained characterization test guarding the AABB box-cut hoist.
		# No library deps; runnable under ctest (exits non-zero on any mismatch).
		add_executable(${PROJECT_NAME}_test_boxcut tests/box_cut_equivalence.cpp)
		target_compile_features(${PROJECT_NAME}_test_boxcut PUBLIC ${CXX14_FEATURES})
		enable_testing()
		add_test(NAME box_cut_equivalence COMMAND ${PROJECT_NAME}_test_boxcut)
	endif()
	###################################################################
	#add_executable(${PROJECT_NAME}_bin app/main.cpp)
	#target_sources(FastEnvelope_bin PRIVATE ${SOURCES})
	#
	#target_compile_features(${PROJECT_NAME}_bin PUBLIC ${CXX14_FEATURES})
	#
	#
	#if(NOT FAST_ENVELOPE_WITH_GEOGRAM_PREDICATES)
	#	# Geogram
	#	fast_envelope_download_geogram()
	#	include(geogram)
	#
	#endif()
	#target_link_libraries(${PROJECT_NAME}_bin PUBLIC geogram)
	#target_link_libraries(${PROJECT_NAME}_bin
	#	PUBLIC
	#		${PROJECT_NAME}
	#		warnings::all
	#)
	#
	#if(FAST_ENVELOPE_WITH_SANITIZERS)
	#	add_sanitizers(${PROJECT_NAME}_bin)
	#endif()
	#
	#
	#if(FAST_ENVELOPE_ENABLE_TBB)
	#	target_link_libraries(${PROJECT_NAME}_bin PUBLIC TBB::tbb)
	#endif()
	#
	#if(NOT (${CMAKE_VERSION} VERSION_LESS "3.6.0"))
	#	set_property(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY VS_STARTUP_PROJECT ${PROJECT_NAME}_bin)
	#endif()

	##########################################################################################################



endif()
