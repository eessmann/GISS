include(CheckIncludeFiles)
include(CheckSymbolExists)

# Finding PkgConfig
find_package(PkgConfig)

# testing system include files
CHECK_INCLUDE_FILES(dlfcn.h HAVE_DLFCN_H LANGUAGE C)
CHECK_INCLUDE_FILES(fenv.h HAVE_FENV_H LANGUAGE C)
CHECK_INCLUDE_FILES(getopt.h HAVE_GETOPT_H LANGUAGE C)
CHECK_INCLUDE_FILES(inttypes.h HAVE_INTTYPES_H LANGUAGE C)
CHECK_INCLUDE_FILES(stdint.h HAVE_STDINT_H LANGUAGE C)
CHECK_INCLUDE_FILES(memory.h HAVE_MEMORY_H LANGUAGE C)
CHECK_INCLUDE_FILES(stdlib.h HAVE_STDLIB_H LANGUAGE C)
CHECK_INCLUDE_FILES(strings.h HAVE_STRINGS_H LANGUAGE C)
CHECK_INCLUDE_FILES(string.h HAVE_STRING_H LANGUAGE C)
CHECK_INCLUDE_FILES(sys/stat.h HAVE_SYS_STAT_H LANGUAGE C)
CHECK_INCLUDE_FILES(sys/types.h HAVE_SYS_TYPES_H LANGUAGE C)
CHECK_INCLUDE_FILES(unistd.h HAVE_UNISTD_H LANGUAGE C)

# testing for system functions
check_symbol_exists(getopt_long "getopt.h" HAVE_GETOPT_LONG)
check_symbol_exists(funopen "stdio.h" HAVE_FUNOPEN)
check_symbol_exists(open_memstream "stdio.h" HAVE_OPEN_MEMSTREAM)

# setting external program
set(HAVE_PKG_CONFIG ${PKG_CONFIG_FOUND})
set(HAVE_MPI ${MPI_FOUND})

# Generating Config header
configure_file(config.h.in ${CMAKE_CURRENT_BINARY_DIR}/config/config.h)

install(FILES ${CMAKE_CURRENT_BINARY_DIR}/config/config.h
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
