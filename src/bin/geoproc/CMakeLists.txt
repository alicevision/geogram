include(${CMAKE_SOURCE_DIR}/cmake/opengl.cmake)

aux_source_directories(SOURCES "" .)
add_executable(geoproc ${SOURCES})
target_link_libraries(geoproc geogram_gfx geogram ${GLFW_LIBRARIES})
install_runtime_targets(geoproc)
