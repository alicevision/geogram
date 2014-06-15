aux_source_directories(SOURCES "Source Files"            .)
aux_source_directories(SOURCES "Source Files\\numerics"     numerics)
aux_source_directories(SOURCES "Source Files\\points"       points)

include_directories(${CMAKE_BINARY_DIR}/src/lib)

add_library(geogram ${SOURCES} $<TARGET_OBJECTS:geogram_third_party>)
if(UNIX AND VORPALINE_BUILD_DYNAMIC)
    target_link_libraries(geogram pthread)
endif()

