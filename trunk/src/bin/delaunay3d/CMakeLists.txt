aux_source_directories(SOURCES "" .)
vor_add_executable(delaunay3d ${SOURCES})
target_link_libraries(delaunay3d geogram)
