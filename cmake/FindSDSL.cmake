# Try to find SDSL
# Once done, this will define
#
# SDSL_FOUND - system has SDSL
# SDSL_INCLUDE_DIR - the SDSL include directories
# SDSL_LIBRARY - the SDSL libraries

if(SDSL_INCLUDE_DIR AND SDSL_LIBRARY)
    set(SDSL_FIND_QUIETLY TRUE)
endif(SDSL_INCLUDE_DIR AND SDSL_LIBRARY)

find_path(SDSL_INCLUDE_DIR
    NAMES sdsl
    HINTS
        ${CMAKE_BINARY_DIR}/external/sdsl_install
        ${CMAKE_BINARY_DIR}/external/sdsl_install/include
        ENV CPATH
        ENV HOME
        ENV C_INCLUDE_PATH
        ENV CPLUS_INCLUDE_PATH
        /usr/include
        /usr/local/include
    PATH_SUFFIXES include
)

find_library(SDSL_LIBRARY
    NAMES sdsl libsdsl.a
    HINTS
        ${CMAKE_BINARY_DIR}/external/sdsl_install
        ${CMAKE_BINARY_DIR}/external/sdsl_install/lib
        ENV LD_LIBRARY_PATH
        ENV HOME
        ENV LIBRARY_PATH
        /usr/lib
        /usr/local/lib
    PATH_SUFFIXES lib
)



# handle the QUIETLY and REQUIRED arguments and set SDSL_FOUND to TRUE if
# all listed variables are TRUE
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(SDSL DEFAULT_MSG SDSL_INCLUDE_DIR SDSL_LIBRARY)

mark_as_advanced(SDSL_INCLUDE_DIR SDSL_LIBRARY)
