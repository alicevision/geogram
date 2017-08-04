aux_source_directories(SOURCES "" .)
vor_add_executable(test_RVC ${SOURCES})
target_link_libraries(test_RVC geogram)


