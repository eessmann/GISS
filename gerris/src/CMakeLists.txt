# Setting up GFS library
# Gerris headers

set(gfs_headers
        adaptive.h
        advection.h
        balance.h
        boundary.h
        cartesian.h
        domain.h
        event.h
        fluid.h
        ftt.h
        ftt_internal.h
        function.h
        graphic.h
        init.h
        isocube.h
        levelset.h
        map.h
        metric.h
        moving.h
        moving2.h
        mpi_boundary.h
        myc.h
        myc2d.h
        ocean.h
        output.h
        particle.h
        poisson.h
        refine.h
        river.h
        simulation.h
        solid.h
        source.h
        spatial.h
        surface.h
        tension.h
        timestep.h
        unstructured.h
        utils.h
        variable.h
        version.h
        vof.h
        wave.h)

set(gfs_sources
        adaptive.c
        advection.c
        balance.c
        boundary.c
        cartesian.c
        domain.c
        event.c
        fluid.c
        ftt.c
        ftt_internal.c
        graphic.c
        init.c
        levelset.c
        map.c
        metric.c
        moving.c
        moving2.c
        mpi_boundary.c
        ocean.c
        open_memstream.c
        output.c
        particle.c
        poisson.c
        refine.c
        river.c
        simulation.c
        solid.c
        source.c
        surface.c
        tension.c
        timestep.c
        unstructured.c
        utils.c
        variable.c
        verify.h
        vof.c
        wave.c)

add_library(gfs_base)
target_sources(gfs_base
        PRIVATE ${gfs_sources}
        PUBLIC ${gfs_headers})

target_link_libraries(gfs_base
        PUBLIC glib::glib)