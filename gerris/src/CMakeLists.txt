# Setting up GFS library
# Gerris headers
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

set(pkgdatadir ${CMAKE_INSTALL_DATADIR}/gerris PARENT_SCOPE)
set(moduledir ${CMAKE_INSTALL_LIBDIR}/gerris PARENT_SCOPE)

# Gerris 3D library
add_library(gfs3D ${gfs_sources})

target_link_libraries(gfs3D
        INTERFACE project_options project_warnings
        PUBLIC CONAN_PKG::glib CONAN_PKG::gts CONAN_PKG::gsl m)

target_include_directories(gfs3D
        PUBLIC ${CMAKE_SOURCE_DIR}/include/gerris ${CMAKE_BINARY_DIR}/config)

if(MPI_FOUND)
    target_link_libraries(gfs3D PUBLIC MPI::MPI_CXX)
endif()

target_compile_definitions(gfs3D
        PUBLIC GFS_DATA_DIR="${CMAKE_INSTALL_FULL_DATADIR}/gerris"
               GFS_MODULES_DIR="${CMAKE_INSTALL_FULL_LIBDIR}/gerris")
set_target_properties(gfs3D PROPERTIES PUBLIC_HEADER ${CMAKE_SOURCE_DIR}/include/gfs.h)

# Gerris 2D library
add_library(gfs2D ${gfs_sources})

target_link_libraries(gfs2D
        INTERFACE project_options project_warnings
        PUBLIC CONAN_PKG::glib CONAN_PKG::gts CONAN_PKG::gsl m)

target_include_directories(gfs2D
        PUBLIC ${CMAKE_SOURCE_DIR}/include/gerris ${CMAKE_BINARY_DIR}/config)

if(MPI_FOUND)
    target_link_libraries(gfs2D PUBLIC MPI::MPI_CXX)
endif()

target_compile_definitions(gfs2D
        PUBLIC  GFS_DATA_DIR="${CMAKE_INSTALL_FULL_DATADIR}/gerris"
                GFS_MODULES_DIR="${CMAKE_INSTALL_FULL_LIBDIR}/gerris"
                FTT_2D=1)
set_target_properties(gfs2D PROPERTIES PUBLIC_HEADER ${CMAKE_SOURCE_DIR}/include/gfs.h)

# Gerris 3D executable
add_executable(gerris3D gerris.c)
target_link_libraries(gerris3D
        INTERFACE project_options project_warnings
        PRIVATE gfs3D)
target_include_directories(gerris3D PRIVATE ${CMAKE_BINARY_DIR}/config)

if(MPI_FOUND)
    target_link_libraries(gerris3D PUBLIC MPI::MPI_CXX)
endif()

# Gerris 2D executable
add_executable(gerris2D gerris.c)
target_link_libraries(gerris2D
        INTERFACE project_options project_warnings
        PRIVATE gfs2D)

target_include_directories(gerris2D PRIVATE ${CMAKE_BINARY_DIR}/config)

if(MPI_FOUND)
    target_link_libraries(gerris2D PUBLIC MPI::MPI_CXX)
endif()

# Generating PC files
configure_file(gerris2D.pc.in ${CMAKE_CURRENT_BINARY_DIR}/gerris2D.pc)
configure_file(gerris3D.pc.in ${CMAKE_CURRENT_BINARY_DIR}/gerris3D.pc)

# Installing Gerris Solver and Library
install(TARGETS gerris2D gerris3D gfs2D gfs3D
        EXPORT GISS
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
        PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
        )

install(DIRECTORY ${CMAKE_SOURCE_DIR}/include/gerris
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/gerris2D.pc ${CMAKE_CURRENT_BINARY_DIR}/gerris3D.pc
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/pkgconfig)
install(PROGRAMS m4.awk
        DESTINATION ${CMAKE_INSTALL_DATADIR}/gerris)
