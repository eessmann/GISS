option(GISS_USE_BUBBLES "Enable Gerris Bubble solver" ON)

# particulate_common
add_library(particulate_common3D)
target_sources(particulate_common3D
        PRIVATE particulatecommon.c
        PUBLIC  particulatecommon.h)

target_link_libraries(particulate_common3D
        PUBLIC gfs3D)

add_library(particulate_common2D)
target_sources(particulate_common2D
        PRIVATE particulatecommon.c
        PUBLIC  particulatecommon.h)

target_link_libraries(particulate_common2D
        PUBLIC gfs2D)

# particulates
add_library(particulates3D)
target_sources(particulates3D PRIVATE particulates.c)
target_link_libraries(particulates3D PUBLIC particulate_common3D)

add_library(particulates2D)
target_sources(particulates2D PRIVATE particulates.c)
target_link_libraries(particulates2D PUBLIC particulate_common2D)


# bubbles
if(GISS_USE_BUBBLES )
    # 3D Bubbles
    add_library(bubbles3D)
    target_sources(bubbles3D PRIVATE bubbles.c)
    target_link_libraries(bubbles3D PUBLIC particulate_common3D)

    # 2D Bubbles
    add_library(bubbles2D)
    target_sources(bubbles2D PRIVATE bubbles.c)
    target_link_libraries(bubbles2D PUBLIC particulate_common2D)
endif()