cmake_minimum_required(VERSION 3.15...3.17)

# Temporary Pretty printing module
include(CMakePrintHelpers)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/modules")

# Based on Gerris Version numbers for compatibility
set(GFS_MAJOR_VERSION 2)
set(GFS_MINOR_VERSION 0)
set(GFS_MICRO_VERSION 1)
set(GFS_INTERFACE_AGE 0)
set(GFS_BINARY_AGE 0)
set(GFS_VERSION
        ${GFS_MAJOR_VERSION}.${GFS_MINOR_VERSION}.${GFS_MICRO_VERSION})

project(GISS
        VERSION ${GFS_VERSION}
        LANGUAGES C CXX)
include(cmake/StandardProjectSettings.cmake)

# Link this 'library' to set the c standard
add_library(project_options INTERFACE)
target_compile_features(project_options INTERFACE c_std_11 cxx_std_17)

if(CMAKE_CXX_COMPILER_ID MATCHES ".*Clang")
  option(ENABLE_BUILD_WITH_TIME_TRACE "Enable -ftime-trace to generate time tracing .json files on clang" OFF)
  if (ENABLE_BUILD_WITH_TIME_TRACE)
    add_compile_definitions(project_options INTERFACE -ftime-trace)
  endif()
endif()
 
# Link this 'library' to use the warnings specified in CompilerWarnings.cmake
add_library(project_warnings INTERFACE)

# enable cache system
include(cmake/Cache.cmake)

# standard compiler warnings
include(cmake/CompilerWarnings.cmake)
set_project_warnings(project_warnings)

# sanitizer options if supported by compiler
include(cmake/Sanitizers.cmake)
enable_sanitizers(project_options)

# enable doxygen
include(cmake/Doxygen.cmake)
enable_doxygen()

# allow for static analysis options
include(cmake/StaticAnalyzers.cmake)

option(BUILD_SHARED_LIBS "Enable compilation of shared libraries" OFF)
option(ENABLE_TESTING "Enable Test Builds" OFF)
option(ENABLE_DOCS "Enable Documentation Builds" OFF)
option(GISS_USE_MPI "Enable compilation with MPI" OFF)
option(GISS_USE_BULLET "Enable Bullet Physic Library" OFF)



# Setting up Conan dependencies based on our needs
# before loading Conan

include(cmake/Conan.cmake)
run_conan()

# Finding OpenMPI
if(GISS_USE_MPI)
  find_package(MPI REQUIRED)
  message(STATUS "Run: ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${MPIEXEC_MAX_NUMPROCS} ${MPIEXEC_PREFLAGS} EXECUTABLE ${MPIEXEC_POSTFLAGS} ARGS")
endif()

# Finding glib
find_package(glib REQUIRED)
if(glib_FOUND)
  message(STATUS "GLIB library found: ${glib_LIBRARIES_TARGETS}")
endif()

if(ENABLE_TESTING)
  enable_testing()
  message(
    "Building Tests and running Test-cases. This might take some time"
  )
  add_subdirectory(test)
endif()

if(ENABLE_DOCS)
  message(
    "Building Documentation"
  )
  add_subdirectory(doc)
endif()

# Core Gerris Solver
add_subdirectory(src)

# Modules including GISS
add_subdirectory(modules)

