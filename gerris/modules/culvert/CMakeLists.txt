# Culvert Module Libraries
add_library(boyd87)
target_sources(boyd87
        PRIVATE boyd87.c
        PUBLIC boyd87.h)
target_link_libraries(boyd87
        PUBLIC m
        INTERFACE project_options project_warnings)

# Culvert Executables
add_executable(check_outlet check_outlet.c)
target_link_libraries(check_outlet
        PRIVATE boyd87 m)

add_executable(check_inlet check_inlet.c)
target_link_libraries(check_inlet
        PRIVATE boyd87 m)

add_executable(discharge_table discharge_table.c)
target_link_libraries(discharge_table
        PRIVATE boyd87 m)
