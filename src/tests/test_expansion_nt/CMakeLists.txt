aux_source_directories(SOURCES "" .)
vor_add_executable(test_expansion_nt ${SOURCES})
target_link_libraries(test_expansion_nt geogram)
