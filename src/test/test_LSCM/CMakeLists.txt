aux_source_directories(SOURCES "" .)
vor_add_executable(test_LSCM ${SOURCES})
target_link_libraries(test_LSCM geogram)
