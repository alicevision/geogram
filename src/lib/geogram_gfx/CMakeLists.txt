include(${CMAKE_SOURCE_DIR}/cmake/opengl.cmake)

add_subdirectory(third_party)

aux_source_directories(SOURCES "Source Files"              .)
aux_source_directories(SOURCES "Source Files\\glut_viewer" glut_viewer)
aux_source_directories(SOURCES "Source Files\\glup_viewer" glup_viewer)
aux_source_directories(SOURCES "Source Files\\api"         api)
aux_source_directories(SOURCES "Source Files\\basic"       basic)
aux_source_directories(SOURCES "Source Files\\mesh"        mesh)
aux_source_directories(SOURCES "Source Files\\GLUP"        GLUP)

add_library(geogram_gfx ${SOURCES} $<TARGET_OBJECTS:geogram_gfx_third_party>)

if(VORPALINE_BUILD_DYNAMIC)
    # Must provide all symbols at link time
    target_link_libraries(geogram_gfx ${OPENGL_LIBRARIES})
    target_link_libraries(geogram_gfx geogram)
    target_link_libraries(geogram_gfx glfw)        
endif()

# Install the library
install_devkit_targets(geogram_gfx)

# There are no public header files to install so far
