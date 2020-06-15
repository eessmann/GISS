# Setting up GFS library
# Gerris headers
set(gfs_headers
        ftt.h
        fluid.h
        variable.h
        output.h
        solid.h
        poisson.h
        advection.h
        boundary.h
        mpi_boundary.h
        timestep.h
        domain.h
        init.h
        refine.h
        event.h
        simulation.h
        graphic.h
        adaptive.h
        source.h
        tension.h
        vof.h
        utils.h
        ocean.h
        wave.h
        levelset.h
        isocube.h
        cartesian.h
        surface.h
        unstructured.h
        map.h
        river.h
        moving.h
        balance.h
        metric.h
        particle.h
        myc.h
        myc2d.h
        spatial.h
        function.h
        version.h)

set(gfs_sources
        ftt.c
        fluid.c
        variable.c
        output.c
        solid.c
        poisson.c
        advection.c
        boundary.c
        mpi_boundary.c
        timestep.c
        domain.c
        init.c
        refine.c
        event.c
        simulation.c
        graphic.c
        adaptive.c
        source.c
        tension.c
        vof.c
        utils.c
        ocean.c
        wave.c
        levelset.c
        cartesian.c
        surface.c
        unstructured.c
        map.c
        river.c
        moving.c
        balance.c
        metric.c
        particle.c
        open_memstream.c
        ftt_internal.c
        moving2.c
        verify.h)

find_package(GTS REQUIRED)

cmake_print_variables(GTS_LIBRARIES CMAKE_GTS_CXX_FLAGS GTS_INCLUDE_DIR GTS_LINK_DIRECTORIES GTS_EXE_LINKER_FLAGS GTS_DEFINITIONS)

add_library(gfs)
target_sources(gfs
        PRIVATE ${gfs_sources} config.h
        PUBLIC ${gfs_headers})

target_link_libraries(gfs
        PUBLIC glib::glib)

if(GISS_USE_MPI)
    target_link_libraries(gfs
            PUBLIC MPI::MPI_C)
endif()