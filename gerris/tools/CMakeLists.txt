# gfs2oogl 3D and 2D
add_executable(gfs2oogl2D gfs2oogl.c)
target_link_libraries(gfs2oogl2D
        PRIVATE gfs2D)

add_executable(gfs2oogl3D gfs2oogl.c)
target_link_libraries(gfs2oogl3D
        PRIVATE gfs3D)
target_include_directories(gfs2oogl2D PRIVATE ${CMAKE_BINARY_DIR}/config)
target_include_directories(gfs2oogl3D PRIVATE ${CMAKE_BINARY_DIR}/config)


install(TARGETS gfs2oogl2D gfs2oogl3D DESTINATION ${CMAKE_INSTALL_BINDIR})

# streamanime
add_executable(streamanime streamanime.c)
target_link_libraries(streamanime
        PRIVATE gfs3D)
target_include_directories(streamanime PRIVATE ${CMAKE_BINARY_DIR}/config)

install(TARGETS streamanime DESTINATION ${CMAKE_INSTALL_BINDIR})

# ppmcombine
add_executable(ppmcombine ppmcombine.c)
target_link_libraries(ppmcombine
        PRIVATE gfs2D)
target_include_directories(ppmcombine PRIVATE ${CMAKE_BINARY_DIR}/config)

install(TARGETS ppmcombine DESTINATION ${CMAKE_INSTALL_BINDIR})

# gfscompare
add_executable(gfscompare2D gfscompare.c)
target_link_libraries(gfscompare2D
        PRIVATE gfs2D)
target_include_directories(gfscompare2D PRIVATE ${CMAKE_BINARY_DIR}/config)

add_executable(gfscompare3D gfscompare.c)
target_link_libraries(gfscompare3D
        PRIVATE gfs3D)
target_include_directories(gfscompare3D PRIVATE ${CMAKE_BINARY_DIR}/config)

install(TARGETS gfscompare2D gfscompare3D DESTINATION ${CMAKE_INSTALL_BINDIR})

# gfscombine
add_executable(gfscombine2D gfscombine.c)
target_link_libraries(gfscombine2D
        PRIVATE gfs2D)
target_include_directories(gfscombine2D PRIVATE ${CMAKE_BINARY_DIR}/config)

add_executable(gfscombine3D gfscombine.c)
target_link_libraries(gfscombine3D
        PRIVATE gfs3D)
target_include_directories(gfscombine3D PRIVATE ${CMAKE_BINARY_DIR}/config)

install(TARGETS gfscombine2D gfscombine3D DESTINATION ${CMAKE_INSTALL_BINDIR})

# gfsjoin
add_executable(gfsjoin2D gfsjoin2.c)
target_link_libraries(gfsjoin2D
        PRIVATE gfs2D)
target_include_directories(gfsjoin2D PRIVATE ${CMAKE_BINARY_DIR}/config)

add_executable(gfsjoin3D gfsjoin2.c)
target_link_libraries(gfsjoin3D
        PRIVATE gfs3D)
target_include_directories(gfsjoin3D PRIVATE ${CMAKE_BINARY_DIR}/config)

install(TARGETS gfsjoin2D gfsjoin3D DESTINATION ${CMAKE_INSTALL_BINDIR})

# shapes
add_executable(shapes shapes.c)
target_link_libraries(shapes
        PRIVATE CONAN_PKG::gts m)
target_include_directories(shapes PRIVATE ${CMAKE_BINARY_DIR}/config)

install(TARGETS shapes DESTINATION ${CMAKE_INSTALL_BINDIR})

# Tool Scripts
set(tool_scripts
    darcs2dist
    bat2gts
    ppm2mpeg
    ppm2theora
    ppm2video
    gfs2gfs
    gfsjoin
    gfsplot)

install(PROGRAMS ${tool_scripts}
        DESTINATION ${CMAKE_INSTALL_BINDIR})

# Generateing build_function
set(SH_INCLUDE_DIRS ${CONAN_INCLUDE_DIRS} ${MPI_C_INCLUDE_DIRS} ${HYPRE_INCLUDE_DIR})
string(REPLACE ";" " -I" SH_INCLUDE_STR "${SH_INCLUDE_DIRS}")

set(SH_LIBS ${CONAN_LIB_DIRS} ${MPI_C_LIBRARIES} ${HYPRE_LIBRARY})
string(REPLACE ";" " -L" SH_LIBS_STR "${SH_LIBS}")

configure_file(build_function.in ${CMAKE_CURRENT_BINARY_DIR}/build_function @ONLY)
install(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/build_function
        DESTINATION ${CMAKE_INSTALL_DATADIR}/gerris)