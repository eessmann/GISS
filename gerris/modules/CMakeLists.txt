
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