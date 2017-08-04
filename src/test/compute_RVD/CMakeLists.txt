aux_source_directories(SOURCES "" .)
vor_add_executable(compute_RVD ${SOURCES})
target_link_libraries(compute_RVD geogram)


