include(${CMAKE_SOURCE_DIR}/cmake/opengl.cmake)

aux_source_directories(SOURCES "" .)
add_executable(demo_GLUP ${SOURCES})
target_link_libraries(demo_GLUP geogram_gfx geogram ${GLFW_LIBRARIES})
