# Setting Fortran Linking Macros
if(ENABLE_FORTRAN)
    include(FortranCInterface)
    FortranCInterface_HEADER(${CMAKE_CURRENT_SOURCE_DIR}/FC.h MACRO_NAMESPACE "FC_")
endif()

# ktd Modules
add_subdirectory(kdt)

# rsurface2kdt Executable
add_subdirectory(RStarTree)

add_executable(rsurface2kdt)
target_sources(rsurface2kdt
        PRIVATE rsurface2kdt.c rsurface.c rsurface.h)
target_link_libraries(rsurface2kdt
        PRIVATE cSmRST kdt m)

# particulate_common
add_library(particulate_common3D)
target_sources(particulate_common3D
        PRIVATE particulatecommon.c
        PUBLIC  particulatecommon.h)

target_link_libraries(particulate_common3D
        INTERFACE project_options project_warnings
        PUBLIC gfs3D)

add_library(particulate_common2D)
target_sources(particulate_common2D
        PRIVATE particulatecommon.c
        PUBLIC  particulatecommon.h)

target_link_libraries(particulate_common2D
        PRIVATE project_options project_warnings
        PUBLIC gfs2D)

# particulates
if(GISS_USE_PART)
    # 3D particulates
    add_library(particulates3D SHARED)
    target_sources(particulates3D PRIVATE particulates.c)
    target_link_libraries(particulates3D
            INTERFACE project_options project_warnings
            PUBLIC particulate_common3D)

    # 2D particulates
    add_library(particulates2D SHARED)
    target_sources(particulates2D PRIVATE particulates.c)
    target_link_libraries(particulates2D
            INTERFACE project_options project_warnings
            PUBLIC particulate_common2D)
endif()


# bubbles
if(GISS_USE_BUBBLES)
    # 3D Bubbles
    add_library(bubbles3D SHARED)
    target_sources(bubbles3D PRIVATE bubbles.c)
    target_link_libraries(bubbles3D
            INTERFACE project_options project_warnings
            PUBLIC particulate_common3D)

    # 2D Bubbles
    add_library(bubbles2D SHARED)
    target_sources(bubbles2D PRIVATE bubbles.c)
    target_link_libraries(bubbles2D
            INTERFACE project_options project_warnings
            PUBLIC particulate_common2D)
endif()

# maps
if(GISS_USE_MAPS)
    # 3D Map
    add_library(map3D SHARED)
    target_sources(map3D PRIVATE map.c)
    target_link_libraries(map3D
            INTERFACE project_options project_warnings
            PUBLIC gfs3D CONAN_PKG::proj)

    # 2D Map
    add_library(map2D SHARED)
    target_sources(map2D PRIVATE map.c)
    target_link_libraries(map2D
            INTERFACE project_options project_warnings
            PUBLIC gfs2D CONAN_PKG::proj)
endif()

# fftw
if(GISS_USE_FFTW)
    # 3D FFT
    add_library(fft3D SHARED)
    target_sources(fft3D PRIVATE fft.c)
    target_link_libraries(fft3D
            INTERFACE project_options project_warnings
            PUBLIC gfs3D CONAN_PKG::fftw)

    # 2D Map
    add_library(fft2D SHARED)
    target_sources(fft2D PRIVATE fft.c)
    target_link_libraries(fft2D
            INTERFACE project_options project_warnings
            PUBLIC gfs2D CONAN_PKG::fftw)
endif()

# terrain
if(GISS_USE_TERRA)
    # 3D Terrain
    add_library(terrain3D SHARED)
    target_sources(terrain3D PRIVATE terrain.c)
    target_link_libraries(terrain3D
            INTERFACE project_options project_warnings
            PUBLIC gfs3D kdt CONAN_PKG::gsl)

    # 2D Terrain
    add_library(terrain2D SHARED)
    target_sources(terrain2D PRIVATE terrain.c)
    target_link_libraries(terrain2D
            INTERFACE project_options project_warnings
            PUBLIC gfs2D kdt CONAN_PKG::gsl)
endif()

# tide
if(GISS_USE_TIDE)
    # FES2004 Tide Model
    add_subdirectory(fes2004)

    # 3D Tide
    add_library(tide3D SHARED)
    target_sources(tide3D PRIVATE tide.c)
    target_link_libraries(tide3D
            INTERFACE project_options project_warnings
            PUBLIC gfs3D fes)

    # 2D Tide
    add_library(tide2D SHARED)
    target_sources(tide2D PRIVATE tide.c)
    target_link_libraries(tide2D
            INTERFACE project_options project_warnings
            PUBLIC gfs2D fes)
endif()

# okada
if(GISS_USE_OKADA)
    # 3D okada
    add_library(okada3D SHARED)
    target_sources(okada3D PRIVATE okada.c)
    target_link_libraries(okada3D
            INTERFACE project_options project_warnings
            PUBLIC gfs3D)

    # 2D okada
    add_library(okada2D SHARED)
    target_sources(okada2D PRIVATE okada.c)
    target_link_libraries(okada2D
            INTERFACE project_options project_warnings
            PUBLIC gfs2D)
endif()

# df3
if(GISS_USE_DF3)
    # Only 3D df3
    add_library(df33D SHARED)
    target_sources(df33D PRIVATE df3.c)
    target_link_libraries(df33D
            INTERFACE project_options project_warnings
            PUBLIC gfs3D)
endif()

# electrohydro
if(GISS_USE_ELECTRO)
    # 3D electrohydro
    add_library(electrohydro3D SHARED)
    target_sources(electrohydro3D PRIVATE electrohydro.c)
    target_link_libraries(electrohydro3D
            INTERFACE project_options project_warnings
            PUBLIC gfs3D)

    # 2D electrohydro
    add_library(electrohydro2D SHARED)
    target_sources(electrohydro2D PRIVATE electrohydro.c)
    target_link_libraries(electrohydro2D
            INTERFACE project_options project_warnings
            PUBLIC gfs2D)
endif()

# skewsymmetric
if(GISS_USE_SKEW)
    # 3D skewsymmetric
    add_library(skewsymmetric3D SHARED)
    target_sources(skewsymmetric3D PRIVATE skewsymmetric.c)
    target_link_libraries(skewsymmetric3D
            INTERFACE project_options project_warnings
            PUBLIC gfs3D)

    # 2D skewsymmetric
    add_library(skewsymmetric2D SHARED)
    target_sources(skewsymmetric2D PRIVATE skewsymmetric.c)
    target_link_libraries(skewsymmetric2D
            INTERFACE project_options project_warnings
            PUBLIC gfs2D)
endif()

# topics
if(GISS_USE_TOPICS)
    # 3D topics
    add_library(topics3D SHARED)
    target_sources(topics3D PRIVATE topics.c)
    target_link_libraries(topics3D
            INTERFACE project_options project_warnings
            PUBLIC gfs3D)

    # 2D topics
    add_library(topics2D SHARED)
    target_sources(topics2D PRIVATE topics.c)
    target_link_libraries(topics2D
            INTERFACE project_options project_warnings
            PUBLIC gfs2D)
endif()

# layered
if(GISS_USE_LAYERED)
    # 3D layered
    add_library(layered3D SHARED)
    target_sources(layered3D PRIVATE layered.c)
    target_link_libraries(layered3D
            INTERFACE project_options project_warnings
            PUBLIC gfs3D)

    # 2D layered
    add_library(layered2D SHARED)
    target_sources(layered2D PRIVATE layered.c)
    target_link_libraries(layered2D
            INTERFACE project_options project_warnings
            PUBLIC gfs2D)
endif()

# culvert
if(GISS_USE_CULVERT)
    # Culvert Model
    add_subdirectory(culvert)

    # 3D culvert
    add_library(culvert3D SHARED)
    target_sources(culvert3D PRIVATE culvert.c)
    target_link_libraries(culvert3D
            INTERFACE project_options project_warnings
            PUBLIC gfs3D boyd87 m)

    # 2D culvert
    add_library(culvert2D SHARED)
    target_sources(culvert2D PRIVATE culvert.c)
    target_link_libraries(culvert2D
            INTERFACE project_options project_warnings
            PUBLIC gfs2D boyd87 m)
endif()

# stokes
if(GISS_USE_STOKES)

    # Fortran Library
    add_library(CW263)
    target_sources(CW263 PUBLIC CW263.f)

    # 3D stokes
    add_library(stokes3D SHARED)
    target_sources(stokes3D PRIVATE stokes.c FC.h)
    set_property(TARGET stokes3D PROPERTY LINKER_LANGUAGE C)
    target_link_libraries(stokes3D
            PRIVATE CW263
            INTERFACE project_options project_warnings
            PUBLIC gfs3D)

    # 2D Tide
    add_library(stokes2D SHARED)
    target_sources(stokes2D PRIVATE stokes.c FC.h)
    set_property(TARGET stokes2D PROPERTY LINKER_LANGUAGE C)
    target_link_libraries(stokes2D
            PRIVATE CW263
            INTERFACE project_options project_warnings
            PUBLIC gfs2D)
endif()