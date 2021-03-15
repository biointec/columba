cmake_minimum_required(VERSION 3.0)

project(fmindex)

add_executable(columba src/main.cpp src/fmindex.cpp  src/searchstrategy.cpp src/nucleotide.cpp)
add_executable(columba_build src/buildIndex.cpp )

SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -pedantic -mpopcnt -std=gnu++11")

#add_compile_definitions(ALPHABET=5)
add_definitions(-DALPHABET=5)

# set the module path
set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")

set(default_build_type "Release")
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "Setting build type to '${default_build_type}' as none was specified.")
  set(CMAKE_BUILD_TYPE "${default_build_type}" CACHE
      STRING "Choose the type of build." FORCE)
  # Set the possible values of build type for cmake-gui
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
    "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()

# set include path for Google's sparse hash table
find_package(SparseHash)
if (SPARSEHASH_FOUND)
    include_directories(${SPARSEHASH_INCLUDE_DIR})
else (SPARSEHASH_FOUND)
    message(FATAL_ERROR "\nFATAL ERROR: The required Google SparseHash package"
            " could not be found on this system.  Please refer to the Detox "
            "manual for the Google Sparsehash installation instructions.  If "
            "you installed Google Sparsehash in a non-standard location "
            "(e.g. somewhere in your homedir), you can point cmake to the "
            "installation location as follows: \ncmake "
            "-DSPARSEHASH_INCLUDE_DIR=<path-to-sparsehash>/include .")
endif(SPARSEHASH_FOUND)

# check if zlib is present
find_package(ZLIB)
if (ZLIB_FOUND)
    add_definitions(-DHAVE_ZLIB)
    include_directories(${ZLIB_INCLUDE_DIR})
endif(ZLIB_FOUND)

target_link_libraries(columba pthread)

if (ZLIB_FOUND)
        target_link_libraries(columba ${ZLIB_LIBRARY})
endif (ZLIB_FOUND)

install(TARGETS columba DESTINATION bin)
install(TARGETS columba_build DESTINATION bin)



