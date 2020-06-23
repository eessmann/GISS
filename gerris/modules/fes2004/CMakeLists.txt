# FES 2004 Tide Module
set(fes_sources
        fes2004_alloc.c
        fes2004_error.c
        fes2004_extraction.c
        fes2004_init.c
        fes2004_io.c
        fes2004_prediction.c
        fes2004_kernel.c)

set(fes_headers
        fes.h
        fes2004_lib.h)


add_library(fes)
target_sources(fes
        PRIVATE ${fes_sources}
        PUBLIC ${fes_header})
target_link_libraries(fes
        INTERFACE project_options project_warnings
        PUBLIC CONAN_PKG::gsl CONAN_PKG::netcdf-c m)