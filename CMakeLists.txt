cmake_minimum_required(VERSION 3.0)

project(fmindex)

add_executable(columba src/main.cpp src/fmindex.cpp src/bidirecfmindex.cpp src/customtypedefs.cpp src/nucleotide.cpp)
add_executable(columba_build src/buildIndex.cpp  src/customtypedefs.cpp)

SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -pedantic -mpopcnt ")

#add_compile_definitions(ALPHABET=5)
add_definitions(-DALPHABET=5)


set(default_build_type "Release")
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "Setting build type to '${default_build_type}' as none was specified.")
  set(CMAKE_BUILD_TYPE "${default_build_type}" CACHE
      STRING "Choose the type of build." FORCE)
  # Set the possible values of build type for cmake-gui
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
    "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()

install(TARGETS columba DESTINATION bin)
install(TARGETS columba_build DESTINATION bin)

