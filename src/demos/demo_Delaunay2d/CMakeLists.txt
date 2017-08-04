include(${CMAKE_SOURCE_DIR}/cmake/opengl.cmake)

aux_source_directories(SOURCES "" .)
add_executable(demo_Delaunay2d ${SOURCES})
target_link_libraries(demo_Delaunay2d geogram geogram_gfx ${GLFW_LIBRARIES})
