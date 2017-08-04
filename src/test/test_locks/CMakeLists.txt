aux_source_directories(SOURCES "" .)
vor_add_executable(test_locks ${SOURCES})
target_link_libraries(test_locks geogram)


