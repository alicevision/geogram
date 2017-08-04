aux_source_directories(SOURCES "" .)
vor_add_executable(hexdom_pipeline ${SOURCES})
target_link_libraries(hexdom_pipeline exploragram geogram)


