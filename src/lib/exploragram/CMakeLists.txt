
aux_source_directories(SOURCES "Source Files"                    .)
aux_source_directories(SOURCES "Source Files\\hexdom"            hexdom)
aux_source_directories(SOURCES "Source Files\\optimal_transport" optimal_transport)

include_directories(${PROJECT_BINARY_DIR}/src/lib)

add_library(exploragram ${SOURCES})

set_target_properties(exploragram PROPERTIES
                      VERSION ${VORPALINE_VERSION}
                      SOVERSION ${VORPALINE_VERSION_MAJOR})

if(UNIX AND VORPALINE_BUILD_DYNAMIC)
    target_link_libraries(exploragram pthread dl)
endif()

if(WIN32)
    target_link_libraries(exploragram psapi)
endif()

target_link_libraries(exploragram geogram)
