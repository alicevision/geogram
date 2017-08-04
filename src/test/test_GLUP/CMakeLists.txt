include(${PROJECT_SOURCE_DIR}/cmake/opengl.cmake)

aux_source_directories(SOURCES "" .)
add_executable(test_GLUP ${SOURCES})

target_link_libraries(test_GLUP geogram geogram_gfx ${GLFW_LIBRARIES})
