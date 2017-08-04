aux_source_directories(SOURCES "" .)
vor_add_executable(vorpastat ${SOURCES})
target_link_libraries(vorpastat geogram)

