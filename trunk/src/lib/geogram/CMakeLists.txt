aux_source_directories(SOURCES "Source Files"               .)
aux_source_directories(SOURCES "Source Files\\numerics"     numerics)
aux_source_directories(SOURCES "Source Files\\points"       points)
aux_source_directories(SOURCES "Source Files\\delaunay"     delaunay)

include_directories(${CMAKE_BINARY_DIR}/src/lib)

add_library(geogram ${SOURCES})

