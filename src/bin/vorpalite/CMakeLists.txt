aux_source_directories(SOURCES "" .)
vor_add_executable(vorpalite ${SOURCES})
target_link_libraries(vorpalite geogram)

install_runtime_targets(vorpalite)

