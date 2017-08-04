aux_source_directories(SOURCES "" .)
vor_add_executable(test_HLBFGS ${SOURCES})
target_link_libraries(test_HLBFGS geogram)


