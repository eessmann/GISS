# FindHYPRE.cmake
#
# Finds the HYPRE library
#
# This will define the following variables
#
#    HYPRE_FOUND
#
# and the following imported targets
#
#     HYPRE::HYPRE
#
# Author: Simon Praetorius <simon.praetorius@tu-dresden.de>


mark_as_advanced(HYPRE_FOUND HYPRE_VERSION HYPRE_LIBRARY HYPRE_HEADER HYPRE_INCLUDE_DIR HYPRE_LIB_DIR)

find_library(HYPRE_LIBRARY HYPRE
        HINTS
        ${HYPRE_ROOT}
        ${PETSC_LIBRARY_DIRS}
        ENV HYPRE_ROOT
        ENV EBROOTHYPRE
        PATH_SUFFIXES lib/ lib64/
        DOC "Library file for HYPRE"
        )

find_file(HYPRE_HEADER HYPRE.h
        HINTS
        ${HYPRE_ROOT}
        ${PETSC_INCLUDE_DIRS}
        ENV HYPRE_ROOT
        ENV EBROOTHYPRE
        PATH_SUFFIXES include/
        DOC "Headerfile hypre.h for HYPRE"
        )
if (HYPRE_HEADER)
    get_filename_component(HYPRE_INCLUDE_DIR ${HYPRE_HEADER} PATH)

    # extract HYPRE version number
    file(STRINGS ${HYPRE_INCLUDE_DIR}/HYPRE_config.h HYPRE_VERSION_LINE REGEX "#define HYPRE_RELEASE_VERSION")
    string(REGEX MATCH "([0-9]+[.][0-9]+[.][0-9]+)" HYPRE_VERSION_REGEX ${HYPRE_VERSION_LINE})
    set(HYPRE_VERSION ${CMAKE_MATCH_1})
endif (HYPRE_HEADER)

find_package(MPI QUIET COMPONENTS CXX)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(HYPRE
        REQUIRED_VARS HYPRE_LIBRARY HYPRE_HEADER MPI_FOUND
        VERSION_VAR HYPRE_VERSION
        )

if(HYPRE_FOUND AND MPI_FOUND AND NOT TARGET HYPRE::HYPRE)
    add_library(HYPRE::HYPRE INTERFACE IMPORTED GLOBAL)

    get_filename_component(HYPRE_LIB_DIR ${HYPRE_LIBRARY} PATH)
    file(GLOB HYPRE_LIBS ${HYPRE_LIB_DIR}/${CMAKE_SHARED_LIBRARY_PREFIX}HYPRE*${CMAKE_SHARED_LIBRARY_SUFFIX})
    list(APPEND HYPRE_LIBS ${HYPRE_LIBRARY})

    find_package(LAPACK QUIET)
    if (LAPACK_FOUND)
        list(APPEND HYPRE_LIBS ${LAPACK_LIBRARIES})
        set_target_properties(HYPRE::HYPRE PROPERTIES
                INTERFACE_LINK_OPTIONS "${LAPACK_LINKER_FLAGS}"
                )
    endif (LAPACK_FOUND)

    set_target_properties(HYPRE::HYPRE PROPERTIES
            INTERFACE_INCLUDE_DIRECTORIES "${HYPRE_INCLUDE_DIR}"
            INTERFACE_LINK_LIBRARIES "${HYPRE_LIBS}"
            )
    target_link_libraries(HYPRE::HYPRE INTERFACE MPI::MPI_CXX)
endif()
