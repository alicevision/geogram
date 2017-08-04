aux_source_directories(SOURCES "" .)
vor_add_executable(geodump ${SOURCES})
target_link_libraries(geodump geogram)

install_runtime_targets(geodump)

