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

# Option to enable move structures for the phi operations (default: OFF)
option(PHI_MOVE "Enable phi move structures" OFF)
if(PHI_MOVE AND NOT RUN_LENGTH_COMPRESSION)
    message(WARNING "PHI_MOVE is not supported with Vanilla (non-RLC) mode and will be ignored.")
endif()

# Option to enable developer mode (default: OFF)
option(DEVELOPER_MODE "Enable developer mode" OFF)
if(DEVELOPER_MODE)
    add_definitions(-DDEVELOPER_MODE)
    message(STATUS "Developer mode enabled.")
else()
    message(STATUS "Developer mode disabled.")
endif()

# Compiler flags for warnings and C++ standard
set(COMMON_FLAGS "-Wall -pedantic -DCMAKE_POLICY_VERSION_MINIMUM=3.5 -std=gnu++11")
include(CheckCXXCompilerFlag)
check_cxx_compiler_flag("-mpopcnt" HAS_POPCNT)
if(HAS_POPCNT)
    set(COMMON_FLAGS "${COMMON_FLAGS} -mpopcnt")
endif()
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${COMMON_FLAGS}")

# Set default 32-bit or 64-bit type based on RUN_LENGTH_COMPRESSION
if(RUN_LENGTH_COMPRESSION)
    set(THIRTY_TWO_DEFAULT OFF)
    add_definitions(-DRUN_LENGTH_COMPRESSION)
    message(STATUS "Run-length compression enabled: using the bidirectional move structure (b-move).")

    if(PHI_MOVE)
        add_definitions(-DPHI_MOVE)
        message(STATUS "Phi move structures enabled.")
    else()
        message(STATUS "Phi move structures disabled.")
    endif()
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

find_package(Git QUIET)
message(STATUS "Checking and updating submodules...")
execute_process(
    COMMAND ${GIT_EXECUTABLE} submodule update --init --recursive
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    RESULT_VARIABLE GIT_SUBMOD_RESULT
)

if(NOT GIT_SUBMOD_RESULT EQUAL "0")
    message(FATAL_ERROR "Failed to initialize git submodules. Please run 'git submodule update --init --recursive' manually.")
endif()


# Update the Big-BWT submodule

if(RUN_LENGTH_COMPRESSION)
    # Find Git
    find_package(Git QUIET)

    set(BIG_BWT_USABLE TRUE)


    # Find Python 3.8 or greater
    find_package(Python3 3.8 COMPONENTS Interpreter)

    if (Python3_FOUND)

	message(STATUS "Python executable detected by CMake: ${Python3_EXECUTABLE}")
 
	# Check Python version
    	execute_process(
        	COMMAND ${Python3_EXECUTABLE} -c "import sys; print(f'{sys.version_info[0]} {sys.version_info[1]}')"
        	OUTPUT_VARIABLE PYTHON_VERSION
        	OUTPUT_STRIP_TRAILING_WHITESPACE
    	)

    	# Split the major and minor version
    	string(REPLACE " " ";" PYTHON_VERSION_LIST "${PYTHON_VERSION}")
    	list(GET PYTHON_VERSION_LIST 0 PYTHON_MAJOR_VERSION)
    	list(GET PYTHON_VERSION_LIST 1 PYTHON_MINOR_VERSION)

    	# Check if the Python version is >= 3.8
    	if (PYTHON_MAJOR_VERSION GREATER 3 OR (PYTHON_MAJOR_VERSION EQUAL 3 AND PYTHON_MINOR_VERSION GREATER_EQUAL 8))
        	message(STATUS "Found Python: ${Python3_EXECUTABLE} (Version ${PYTHON_MAJOR_VERSION}.${PYTHON_MINOR_VERSION})")

        	# Check if psutil is installed
        	execute_process(
            		COMMAND ${Python3_EXECUTABLE} -c "import psutil"
            		RESULT_VARIABLE PSUTIL_INSTALLED
            		ERROR_QUIET
            		OUTPUT_QUIET
        	)

        	if (NOT PSUTIL_INSTALLED EQUAL "0")
            		message(WARNING "Python package 'psutil' is required for prefix-free parsing, but not installed. Please install it using pip (pip install psutil) if you want to use prefix-free parsing.")
            		set(BIG_BWT_USABLE FALSE)
        	else()
            		message(STATUS "Python 3.8 or greater and psutil are installed.")
            		set(BIG_BWT_USABLE TRUE)
        	endif()

    	else()
        	message(WARNING "Found Python version ${PYTHON_MAJOR_VERSION}.${PYTHON_MINOR_VERSION}, but Python 3.8 or greater is required for prefix-free parsing.")
        	set(BIG_BWT_USABLE FALSE)
    	endif()
    else()
	    message(WARNING "Python not found, but is required for prefix-free parsing")
	    set(BIG_BWT_USABLE FALSE)
    endif()


    # Proceed with submodule update if usable
    if(GIT_FOUND AND EXISTS "${PROJECT_SOURCE_DIR}/.git")
        if(BIG_BWT_USABLE)
            message(STATUS "Updating Big-BWT submodule")
            execute_process(COMMAND ${GIT_EXECUTABLE} submodule update --init -q external/Big-BWT
                            WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                            RESULT_VARIABLE GIT_SUBMOD_RESULT)
            if(NOT GIT_SUBMOD_RESULT EQUAL "0")
                message(FATAL_ERROR "git submodule update --init failed with ${GIT_SUBMOD_RESULT}, please checkout submodules")
            endif()
        else()
            message(WARNING "Prefix-free parsing is not usable due to missing dependencies (Python 3.8 and psutil), if you do not want to use prefix-free parsing you can ignore this message.")
        endif()
    endif()
endif()

# ---- Executables and Linking ----

add_executable(columba src/parallel.cpp src/fastq.cpp src/seqfile.cpp src/util.cpp)
target_link_libraries(columba AuxiliaryClasses)

if(RUN_LENGTH_COMPRESSION)
    add_executable(columba_build src/logger.cpp src/buildindex.cpp src/bmove/moverepr.cpp src/parameters/buildparameters.cpp src/seqfile.cpp src/util.cpp)
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
    add_executable(columba_build src/buildindex.cpp src/logger.cpp src/parameters/buildparameters.cpp src/seqfile.cpp )
endif()

# ---- External Dependencies ----

# overwrite the message function to suppress messages if MESSAGE_QUIET is set
function(message)
    if (NOT MESSAGE_QUIET)
        _message(${ARGN})
    endif()
endfunction()
set(MESSAGE_QUIET OFF)

set(CMAKE_SUPPRESS_DEVELOPER_WARNINGS 1 CACHE INTERNAL "Suppress developer warnings")
set(CMAKE_WARN_DEPRECATED OFF CACHE BOOL "" FORCE)

# FetchContent module for including external projects
include(FetchContent)

# libsais library (added to columba_build)

if(THIRTY_TWO)
    # Set specific build options for libsais
    message(STATUS "Fetching libsais library")
    set(LIBSAIS_USE_OPENMP OFF CACHE BOOL "Disable OpenMP support in libsais")
    set(LIBSAIS_BUILD_SHARED_LIB OFF CACHE BOOL "Build libsais as a static library")

    FetchContent_Declare(
    libsais
    GIT_REPOSITORY https://github.com/IlyaGrebnov/libsais.git
    GIT_TAG v2.8.5
    )
    FetchContent_MakeAvailable(libsais)
    target_link_libraries(columba_build libsais)
else()
    find_package(Git QUIET)
    # update git submodule libdivsufsort
    message(STATUS "Updating libdivsufsort submodule")
    execute_process(COMMAND ${GIT_EXECUTABLE} submodule update --init  external/libdivsufsort
                    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                    RESULT_VARIABLE GIT_SUBMOD_RESULT)
    if(NOT GIT_SUBMOD_RESULT EQUAL "0")
        message(FATAL_ERROR "git submodule update --init failed with ${GIT_SUBMOD_RESULT}, please checkout submodules")
    endif()


    set(BUILD_DIVSUFSORT64 ON CACHE BOOL "Build divsufsort64")
    set(BUILD_EXAMPLES OFF CACHE BOOL "Build examples")
   

    
    set(MESSAGE_QUIET ON)
    add_subdirectory(external/libdivsufsort)
    set(MESSAGE_QUIET OFF)

     # If needed, add dependencies to ensure libdivsufsort is built first
    add_dependencies(columba_build divsufsort64)


    # Link against the libdivsufsort library
    target_link_libraries(columba_build  divsufsort64)


    # Specify include directories for the headers that are generated during build
    # Set the include path to the appropriate build directory where the headers are located
    target_include_directories(columba_build PRIVATE
        ${CMAKE_CURRENT_BINARY_DIR}/external/libdivsufsort/include
    )

   

endif()

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


# Enable IWYU (include-what-you-use) if desired
option(ENABLE_IWYU "Enable IWYU" OFF)
if(ENABLE_IWYU)
    find_program(IWYU_EXECUTABLE NAMES include-what-you-use iwyu)
    if(IWYU_EXECUTABLE)
        message(STATUS "Found IWYU: ${IWYU_EXECUTABLE}")
        set_property(TARGET AuxiliaryClasses PROPERTY CXX_INCLUDE_WHAT_YOU_USE ${IWYU_EXECUTABLE})
        set_property(TARGET columba PROPERTY CXX_INCLUDE_WHAT_YOU_USE ${IWYU_EXECUTABLE})
        set_property(TARGET columba_build PROPERTY CXX_INCLUDE_WHAT_YOU_USE ${IWYU_EXECUTABLE})
    else()
        message(FATAL_ERROR "IWYU not found. Please install IWYU or disable ENABLE_IWYU.")
    endif()
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
    target_link_libraries(columba_build ${ZLIB_LIBRARY})
endif()

# Link Threads package
find_package(Threads REQUIRED)
target_link_libraries(columba Threads::Threads)

# SDLS library (if RUN_LENGTH_COMPRESSION is enabled)
if(RUN_LENGTH_COMPRESSION)
  # Set the path to your custom FindSDSL.cmake module
    set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/")
    find_package(SDSL QUIET)

    if(SDSL_FOUND)
        message(STATUS "Found SDSL library at ${SDSL_LIBRARY}")
        message(STATUS "SDSL include directory: ${SDSL_INCLUDE_DIR}")

        include_directories(${SDSL_INCLUDE_DIR})
        target_link_libraries(columba ${SDSL_LIBRARY})
        target_link_libraries(columba_build ${SDSL_LIBRARY})
      

    else()
        message(STATUS "SDSL-lite not found, fetching and installing it locally...")

        include(ExternalProject)

        set(SDSL_INSTALL_DIR "${CMAKE_BINARY_DIR}/external/sdsl_install")

        ExternalProject_Add(sdsl-lite_ext
            GIT_REPOSITORY https://github.com/simongog/sdsl-lite.git
            GIT_TAG master
            CONFIGURE_COMMAND ""
            BUILD_COMMAND ""
            INSTALL_COMMAND bash <SOURCE_DIR>/install.sh ${SDSL_INSTALL_DIR}
            UPDATE_DISCONNECTED 1
            LOG_DOWNLOAD ON
            LOG_UPDATE ON
            LOG_CONFIGURE ON
            LOG_BUILD ON
            LOG_INSTALL ON
        )

        # Post-install cleanup of gtest that sdsl-lite installs
        add_custom_command(
            TARGET sdsl-lite_ext
            POST_BUILD
            COMMAND ${CMAKE_COMMAND} -E echo "Removing unnecessary gtest from SDSL-lite install..."
            COMMAND ${CMAKE_COMMAND} -E remove_directory "${SDSL_INSTALL_DIR}/include/gtest"
            COMMAND ${CMAKE_COMMAND} -E remove "${SDSL_INSTALL_DIR}/lib/libgtest.a"
            COMMAND ${CMAKE_COMMAND} -E remove "${SDSL_INSTALL_DIR}/lib/libgtest_main.a"
        )



        # Point to where it’ll be installed
        set(SDSL_INCLUDE_DIR "${SDSL_INSTALL_DIR}/include")
        set(SDSL_LIBRARY "${SDSL_INSTALL_DIR}/lib/libsdsl.a")

        include_directories(${SDSL_INCLUDE_DIR})
        target_link_libraries(columba ${SDSL_LIBRARY})
        target_link_libraries(columba_build ${SDSL_LIBRARY})

        # Ensure your targets wait until SDSL is built
        add_dependencies(AuxiliaryClasses sdsl-lite_ext)
        add_dependencies(columba sdsl-lite_ext)
        add_dependencies(columba_build sdsl-lite_ext)

    endif()
endif()

# ---- Installation ----

# Install targets to the bin directory
install(TARGETS columba columba_build DESTINATION bin)

