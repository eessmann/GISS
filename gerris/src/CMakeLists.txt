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
        ftt_internal.h
        moving2.h
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


group 1 = "ftt.c", "ftt.h", "moving.c", "advection.h", "moving.h", "utils.h", \
"timestep.h", "poisson.h", "domain.h", "variable.h", "fluid.c", \
"fluid.h", "surface.c", "boundary.h", "ftt_internal.h", "surface.h", \
"ftt_internal.c", "moving2.h", "moving2.c"

group 2 = "river.h", "init.c", "map.h", "metric.c", "metric.h", "source.h", \
"refine.h", "simulation.c", "ocean.h", "ocean.c", "timestep.c", \
"tension.h", "poisson.c", "advection.c", "map.c", "particle.h", \
"particle.c"

group 3 = "domain.c", "mpi_boundary.h", "config.h", "solid.c", "version.h", \
"balance.c", "balance.h", "unstructured.h", "boundary.c", \
"open_memstream.c", "graphic.c", "isocube.h", "unstructured.c", \
"mpi_boundary.c", "verify.h"

group 4 = "adaptive.h", "river.c", "wave.h", "init.h", "simulation.h", \
"refine.c", "output.c", "solid.h", "event.h", "event.c", \
"adaptive.c", "graphic.h", "source.c", "wave.c", "output.h"

group 5 = "vof.h", "variable.c", "vof.c", "levelset.h", "tension.c", "myc.h", \
"levelset.c", "myc2d.h"

group 6 = "cartesian.c", "cartesian.h", "utils.c"

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