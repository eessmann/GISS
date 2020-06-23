cmake_minimum_required(VERSION 3.15...3.17)

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

# Temporary Pretty printing module
include(CMakePrintHelpers)
include(GNUInstallDirs)

include(cmake/StandardProjectSettings.cmake)

# Link this 'library' to set the c standard
add_library(project_options INTERFACE)
target_compile_features(project_options INTERFACE c_std_11)
set(CMAKE_C_STANDARD 11)

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
option(GISS_USE_MPI "Enable compilation with MPI" ON)
option(GISS_USE_MODULES "Enable compilation of GFS Modules" ON)

# GFS Module Options
if(GISS_USE_MODULES)
  option(GISS_USE_BULLET "Enable Bullet Physic Library" OFF)
  option(GISS_USE_BUBBLES "Enable Bubble Module" ON)
  option(GISS_USE_MAPS "Enable Map Module" ON)
  option(GISS_USE_PART "Enable Particle Module" ON)
  option(GISS_USE_FFTW "Enable FFTW Module" ON)
  option(GISS_USE_TERRA "Enable Terrain Module" ON)
  option(GISS_USE_TIDE "Enable Tide Module" ON)
  option(GISS_USE_OKADA "Enable Okada Module" ON)
  option(GISS_USE_DF3 "Enable DF3 Module" ON)
  option(GISS_USE_ELECTRO "Enable Electrohydro Module" ON)
  option(GISS_USE_SKEW "Enable Skewsymmetric Module" ON)
  option(GISS_USE_TOPICS "Enable Topics Module" ON)
  option(GISS_USE_LAYERED "Enable Layered Module" ON)
  option(GISS_USE_CULVERT "Enable Culvert Module" ON)
  option(GISS_USE_STOKES "Enable Stokes Module" ON)
  option(GISS_USE_HYPRE "Enable HYPRE Module" OFF)
endif()

# Enabling Fortran for modules
if(GISS_USE_STOKES)
  message(WARNING "Enabling Fortran support, IPO will be disabled")
  set(ENABLE_IPO OFF)
  set(ENABLE_FORTRAN ON)
  enable_language(Fortran)
endif()


# Set up some extra Conan dependencies based on our needs
# before loading Conan
set(CONAN_EXTRA_REQUIRES "")
set(CONAN_EXTRA_OPTIONS "")

# Finding glib
set(CONAN_EXTRA_REQUIRES ${CONAN_EXTRA_REQUIRES}
        glib/2.64.0@bincrafters/stable)

set(CONAN_EXTRA_OPTIONS ${CONAN_EXTRA_OPTIONS}
        glib:shared=True)

# Finding gts
set(CONAN_EXTRA_REQUIRES ${CONAN_EXTRA_REQUIRES}
        gts/0.7.6@eessmann/conan-gts)

set(CONAN_EXTRA_OPTIONS ${CONAN_EXTRA_OPTIONS}
        gts:shared=False)

# Finding gsl
set(CONAN_EXTRA_REQUIRES ${CONAN_EXTRA_REQUIRES}
        gsl/2.6@eessmann/conan-gsl)

# FFTW Library
if(GISS_USE_FFTW)
  set(CONAN_EXTRA_REQUIRES ${CONAN_EXTRA_REQUIRES}
          fftw/3.3.8@bincrafters/stable)
  set(CONAN_EXTRA_OPTIONS ${CONAN_EXTRA_OPTIONS}
          fftw:openmp=True
          fftw:threads=True
          fftw:combinedthreads=True)
endif()

# PROJ Library
if(GISS_USE_MAPS)
  set(CONAN_EXTRA_REQUIRES ${CONAN_EXTRA_REQUIRES}
          proj/7.0.1)
endif()

# BULLET Library
if(GISS_USE_BULLET)
  set(CONAN_EXTRA_REQUIRES ${CONAN_EXTRA_REQUIRES}
          bullet3/2.87@bincrafters/stable)
endif()

# NETCDF Library
if(GISS_USE_TIDE)
  set(CONAN_EXTRA_REQUIRES ${CONAN_EXTRA_REQUIRES}
          netcdf-c/4.7.4@eessmann/conan-netcdf-c)
endif()

    include(cmake/Conan.cmake)
run_conan()

# Setting up OpenMPI
if(GISS_USE_MPI)
  find_package(MPI REQUIRED)
  message(STATUS "Run: ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${MPIEXEC_MAX_NUMPROCS} ${MPIEXEC_PREFLAGS} EXECUTABLE ${MPIEXEC_POSTFLAGS} ARGS")
endif()

# HYPRE library
if(GISS_USE_HYPRE)
  if(MPI_FOUND)
    find_package(HYPRE REQUIRED)
  else()
    set(GISS_USE_MPI ON)
    message(WARNING "HYPRE Requires MPI")
  endif()
endif()


# Core Gerris Solver
add_subdirectory(src)

# Gerris Modules
if(GISS_USE_MODULES)
  add_subdirectory(modules)
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

