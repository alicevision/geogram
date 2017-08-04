aux_source_directories(SOURCES "" .)
vor_add_executable(test_MH ${SOURCES})
target_link_libraries(test_MH geogram)
