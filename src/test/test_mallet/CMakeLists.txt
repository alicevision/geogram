aux_source_directories(SOURCES "" .)
vor_add_executable(test_mallet ${SOURCES})
target_link_libraries(test_mallet geogram)
