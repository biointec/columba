# CMake minimum version and project name
cmake_minimum_required(VERSION 3.14)
project(columba)

# Enforce compiler restrictions
if(WIN32 AND CMAKE_COMPILER_IS_MINGW)
    message(FATAL_ERROR "MinGW-w64 is not supported. Please use a different compiler.")
endif()



# ---- Configuration Options ----

# Option to enable run-length compression (default: OFF)
option(RUN_LENGTH_COMPRESSION "Enable run-length compression" OFF)

# Option to enable developer mode (default: OFF)
option(DEVELOPER_MODE "Enable developer mode" OFF)
if(DEVELOPER_MODE)
    add_definitions(-DDEVELOPER_MODE)
    message(STATUS "Developer mode enabled.")
else()
    message(STATUS "Developer mode disabled.")
endif()

# Compiler flags for warnings and C++ standard
set(COMMON_FLAGS "-Wall -pedantic -mpopcnt -std=gnu++11")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${COMMON_FLAGS}")

# Set default 32-bit or 64-bit type based on RUN_LENGTH_COMPRESSION
if(RUN_LENGTH_COMPRESSION)
    set(THIRTY_TWO_DEFAULT OFF)
    add_definitions(-DRUN_LENGTH_COMPRESSION)
    message(STATUS "Run-length compression enabled: using the bidirectional move structure (b-move).")
else()
    set(THIRTY_TWO_DEFAULT ON)
    message(STATUS "Run-length compression disabled: using the bidirectional FM-index.")
endif()

option(THIRTY_TWO "Use 32-bit types" ${THIRTY_TWO_DEFAULT})
if(THIRTY_TWO)
    add_definitions(-DTHIRTY_TWO)
    message(STATUS "Using 32-bit types.")
else()
    message(STATUS "Using 64-bit types.")
endif()

# Always define ALPHABET as 5
add_definitions(-DALPHABET=5)

# ---- Restrict Windows to non-RLC mode ----
if(WIN32)
    if(RUN_LENGTH_COMPRESSION)
        message(FATAL_ERROR "Run-length compression is not supported on Windows. Please disable RUN_LENGTH_COMPRESSION.")
    endif()
endif()

# ---- Source Files ----

set(SOURCE_FILES
    src/indexhelpers.cpp
    src/bitparallelmatrix.cpp
    src/searchstrategy.cpp
    src/nucleotide.cpp
    src/search.cpp
    src/logger.cpp
    src/parameters/alignparameters.cpp
    src/indexinterface.cpp
)

# Add compression-specific or FM-index-specific source files
if(RUN_LENGTH_COMPRESSION)
    list(APPEND SOURCE_FILES src/bmove/bmove.cpp src/bmove/moverepr.cpp)
else()
    list(APPEND SOURCE_FILES src/fmindex/fmindex.cpp)
endif()

# Create the auxiliary library
add_library(AuxiliaryClasses ${SOURCE_FILES})

# ---- Update Git Submodules ----

# Update the Big-BWT submodule

if(RUN_LENGTH_COMPRESSION)
    # Find Git
    find_package(Git QUIET)

    set(BIG_BWT_USABLE TRUE)

    # Find Python 3.8 or greater
    find_package(Python 3.8 REQUIRED COMPONENTS Interpreter)

    if(NOT Python_FOUND)
        message(WARNING "Python 3.8 or greater is required for prefix-free parsing, but not found.")
        set(BIG_BWT_USABLE FALSE)
    else()
        message(STATUS "Found Python: ${Python_EXECUTABLE}")

        # Check if psutil is installed
        execute_process(
            COMMAND ${Python_EXECUTABLE} -c "import psutil"
            RESULT_VARIABLE PSUTIL_INSTALLED
            ERROR_QUIET
            OUTPUT_QUIET
        )

        if(NOT PSUTIL_INSTALLED EQUAL "0")
            message(WARNING "Python package 'psutil' is required for prefix-free parsing, but not installed. Please install it using pip (pip install psutil) if you want to use prefix-free parsing.")
            set(BIG_BWT_USABLE FALSE)
        else()
            message(STATUS "Python 3.8 or greater and psutil are installed.")
            set(BIG_BWT_USABLE TRUE)
        endif()
    endif()


    # Proceed with submodule update if usable
    if(GIT_FOUND AND EXISTS "${PROJECT_SOURCE_DIR}/.git")
        if(BIG_BWT_USABLE)
            message(STATUS "Updating Big-BWT submodule")
            execute_process(COMMAND ${GIT_EXECUTABLE} submodule update --init external/Big-BWT
                            WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                            RESULT_VARIABLE GIT_SUBMOD_RESULT)
            if(NOT GIT_SUBMOD_RESULT EQUAL "0")
                message(FATAL_ERROR "git submodule update --init failed with ${GIT_SUBMOD_RESULT}, please checkout submodules")
            endif()
        else()
            message(WARNING "Prefix-free parsing is not usable due to missing dependencies (Python 3.8 or psutil).")
        endif()
    endif()
endif()

# ---- Executables and Linking ----

add_executable(columba src/parallel.cpp src/fastq.cpp src/seqfile.cpp src/util.cpp)
target_link_libraries(columba AuxiliaryClasses)

if(RUN_LENGTH_COMPRESSION)
    add_executable(columba_build src/logger.cpp src/buildindex.cpp src/bmove/moverepr.cpp src/parameters/buildparameters.cpp)
    if(BIG_BWT_USABLE)
        add_definitions(-DBIG_BWT_USABLE)
        # Build the Big-BWT submodule if compression is enabled
        add_custom_command(
            OUTPUT ${CMAKE_CURRENT_SOURCE_DIR}/external/Big-BWT/bwtcheck
            COMMAND ${CMAKE_COMMAND} -E chdir ${CMAKE_CURRENT_SOURCE_DIR}/external/Big-BWT make
            COMMENT "Building Big-BWT executable"
        )
        add_custom_target(bigbwt_target ALL DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/external/Big-BWT/bwtcheck)
        add_dependencies(columba_build bigbwt_target)
    endif()
        # Copy the script to the build directory after columba_build is built
        add_custom_command(
            TARGET columba_build POST_BUILD
            COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_SOURCE_DIR}/src/bmove/columba_build_pfp.sh ${CMAKE_BINARY_DIR}/columba_build_pfp.sh
            COMMENT "Copying columba_build_pfp.sh to the build directory"
        )
    
else()
    add_executable(columba_build src/buildindex.cpp src/logger.cpp src/parameters/buildparameters.cpp )
endif()

# ---- External Dependencies ----

# FetchContent module for including external projects
include(FetchContent)

# libsais library (added to columba_build)

# Set specific build options for libsais
message(STATUS "Fetching libsais library")
set(LIBSAIS_USE_OPENMP OFF CACHE BOOL "Disable OpenMP support in libsais")
set(LIBSAIS_BUILD_SHARED_LIB OFF CACHE BOOL "Build libsais as a static library")

FetchContent_Declare(
  libsais
  GIT_REPOSITORY https://github.com/IlyaGrebnov/libsais.git
  GIT_TAG v2.8.3
)
FetchContent_MakeAvailable(libsais)
target_link_libraries(columba_build libsais)

# parallel-hashmap library
message(STATUS "Fetching parallel-hashmap library")
FetchContent_Declare(
    parallel-hashmap
    GIT_REPOSITORY https://github.com/greg7mdp/parallel-hashmap.git
    GIT_TAG e5b892baed478513adcb6425773cae1eda033057
)
FetchContent_MakeAvailable(parallel-hashmap)
include_directories(${parallel-hashmap_SOURCE_DIR})

# fmt library
message(STATUS "Fetching fmt library")
FetchContent_Declare(
  fmtlib
  GIT_REPOSITORY https://github.com/fmtlib/fmt.git
  GIT_TAG 10.2.1
)
FetchContent_MakeAvailable(fmtlib)
target_link_libraries(AuxiliaryClasses fmt::fmt)
target_link_libraries(columba fmt::fmt)

# FetchContent module for including external projects
include(FetchContent)

# Fetch fmtlib from GitHub
FetchContent_Declare(
  fmtlib
  GIT_REPOSITORY https://github.com/fmtlib/fmt.git
  GIT_TAG 10.2.1 # Adjust the version as needed
)

FetchContent_MakeAvailable(fmtlib)

target_link_libraries(AuxiliaryClasses fmt::fmt)
target_link_libraries(columba  fmt::fmt)
if(RUN_LENGTH_COMPRESSION)
	target_link_libraries(columba_build fmt::fmt)
endif()

# ---- Build Configuration ----

# Set default build type to Release if not specified
set(default_build_type "Release")
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    message(STATUS "Setting build type to '${default_build_type}' as none was specified.")
    set(CMAKE_BUILD_TYPE "${default_build_type}" CACHE STRING "Choose the type of build." FORCE)
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()

# ----Installed Dependencies and Linking ----

# Check for zlib
find_package(ZLIB)
if(ZLIB_FOUND)
    add_definitions(-DHAVE_ZLIB)
    include_directories(${ZLIB_INCLUDE_DIR})
    target_link_libraries(columba ${ZLIB_LIBRARY})
endif()

# Link Threads package
find_package(Threads REQUIRED)
target_link_libraries(columba Threads::Threads)

# SDLS library (if RUN_LENGTH_COMPRESSION is enabled)
if(RUN_LENGTH_COMPRESSION)
  # Set the path to your custom FindSDSL.cmake module
    set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/")
    find_package(SDSL REQUIRED)

    if(SDSL_FOUND)
        message(STATUS "Found SDSL library at ${SDSL_LIBRARY}")
        message(STATUS "SDSL include directory: ${SDSL_INCLUDE_DIR}")

        # Include directories
        include_directories(${SDSL_INCLUDE_DIR})
        # Link the SDSL library
        target_link_libraries(columba ${SDSL_LIBRARY})
        target_link_libraries(columba_build ${SDSL_LIBRARY})
    else()
        message(FATAL_ERROR "\nFATAL ERROR: SDSL library not found. Please install it and try again. If "
            "you installed SDSL in a non-standard location, you can point cmake to the "
            "installation location as follows: \ncmake -DSDSL_INCLUDE_DIR=<path-to-sdsl>/include .")
    endif()

endif()

# ---- Installation ----

# Install targets to the bin directory
install(TARGETS columba columba_build  DESTINATION bin)


