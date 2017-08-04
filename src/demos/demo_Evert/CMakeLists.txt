include(${CMAKE_SOURCE_DIR}/cmake/opengl.cmake)

aux_source_directories(SOURCES "" .)
add_executable(demo_evert ${SOURCES})
target_link_libraries(demo_evert geogram_gfx geogram ${GLFW_LIBRARIES})
