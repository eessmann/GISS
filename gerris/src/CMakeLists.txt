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
        #verify.h
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
        #open_memstream.c
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
        vof.c
        wave.c)

set(pkgdatadir ${CMAKE_INSTALL_DATADIR}/gerris)
set(moduledir ${CMAKE_INSTALL_LIBDIR}/gerris)

# Gerris 3D library
add_library(gfs3D)
target_sources(gfs3D
        PUBLIC config.h ${gfs_headers} ${gfs_sources})

target_link_libraries(gfs3D
        INTERFACE project_options #project_warnings
        PUBLIC glib::glib GTS::GTS m)

if(MPI_FOUND)
    target_link_libraries(gfs3D PUBLIC MPI::MPI_CXX)
endif()

target_compile_definitions(gfs3D
        PUBLIC GFS_DATA_DIR="${pkgdatadir}" GFS_MODULES_DIR="${moduledir}")


# Gerris 2D library
add_library(gfs2D)
target_sources(gfs2D
        PUBLIC config.h ${gfs_headers} ${gfs_sources})

target_link_libraries(gfs2D
        INTERFACE project_options #project_warnings
        PUBLIC glib::glib GTS::GTS m)

if(MPI_FOUND)
    target_link_libraries(gfs2D PUBLIC MPI::MPI_CXX)
endif()

target_compile_definitions(gfs2D
        PUBLIC GFS_DATA_DIR="${pkgdatadir}" GFS_MODULES_DIR="${moduledir}" FTT_2D=1)

# Gerris 3D executable
add_executable(gerris3D gerris.c)
target_link_libraries(gerris3D
        INTERFACE project_options #project_warnings
        PRIVATE gfs3D)

# Gerris 2D executable
add_executable(gerris2D gerris.c)
target_link_libraries(gerris2D
        INTERFACE project_options #project_warnings
        PRIVATE gfs2D)

