aux_source_directories(SOURCES "" .)
vor_add_executable(test_opennl ${SOURCES})
target_link_libraries(test_opennl geogram)


