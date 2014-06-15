aux_source_directories(SOURCES "" .)
include_directories(.)
vor_add_executable(mcc ${SOURCES})
target_link_libraries(mcc geogram)
