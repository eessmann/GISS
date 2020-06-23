# kdt Library
add_library(kdt)
target_sources(kdt
        PRIVATE kdt.c
        PUBLIC kdt.h)
target_link_libraries(kdt
        INTERFACE project_options project_warnings)
target_compile_definitions(kdt
        PUBLIC FILE_OFFSET_BITS=64)

# xyz2kdt Executable
add_executable(xyz2kdt)
target_sources(xyz2kdt PRIVATE xyz2kdt.c)
target_link_libraries(xyz2kdt
        PRIVATE kdt m
        INTERFACE project_options project_warnings)

# kdt2kdt Executable
add_executable(kdt2kdt)
target_sources(kdt2kdt PRIVATE kdt2kdt.c)
target_link_libraries(kdt2kdt
        PRIVATE kdt m
        INTERFACE project_options project_warnings)

# kdtquery Executable
add_executable(kdtquery)
target_sources(kdtquery PRIVATE kdtquery.c)
target_link_libraries(kdtquery
        PRIVATE kdt m
        INTERFACE project_options project_warnings)