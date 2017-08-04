# Tetgen (by Hang Si) is currently included by default.
# Triangle (by Jonathan Shewchuk) is currently included by default.

set(GEOGRAM_WITH_TETGEN TRUE)
set(GEOGRAM_WITH_TRIANGLE TRUE)

add_subdirectory(third_party)

aux_source_directories(SOURCES "Source Files"            .)
aux_source_directories(SOURCES "Source Files\\basic"     basic)
aux_source_directories(SOURCES "Source Files\\numerics"  numerics)
aux_source_directories(SOURCES "Source Files\\mesh"      mesh)
aux_source_directories(SOURCES "Source Files\\delaunay"  delaunay)
aux_source_directories(SOURCES "Source Files\\voronoi"   voronoi)
aux_source_directories(SOURCES "Source Files\\points"    points)
aux_source_directories(SOURCES "Source Files\\api"       api)
aux_source_directories(SOURCES "Source Files\\NL"        NL)

if(GEOGRAM_WITH_TETGEN)
   add_definitions(-DGEOGRAM_WITH_TETGEN)
endif()

if(GEOGRAM_WITH_TRIANGLE)
   add_definitions(-DGEOGRAM_WITH_TRIANGLE)
endif()

# Parallel Delaunay now works on both Windows and Unix !
add_definitions(-DGEOGRAM_WITH_PDEL)

# Copy the default version info to the build directory
configure_file(basic/version.h.in version.h @ONLY)
list(APPEND SOURCES version.h)
set_source_files_properties(version.h PROPERTIES GENERATED true)

include_directories(${CMAKE_BINARY_DIR}/src/lib)

add_library(geogram ${SOURCES} $<TARGET_OBJECTS:geogram_third_party>)

if(UNIX AND VORPALINE_BUILD_DYNAMIC)
    target_link_libraries(geogram pthread dl)
endif()

if(WIN32)
    target_link_libraries(geogram psapi)
endif()

# Install the library
install_devkit_targets(geogram)

# Install include files for the standard devkit
install(
    DIRECTORY api
    DESTINATION include/geogram
    COMPONENT devkit
    FILES_MATCHING PATTERN *.h
)

# Install include files for the full devkit
install(
    DIRECTORY .
    DESTINATION include/geogram
    COMPONENT devkit-full
    FILES_MATCHING PATTERN *.h
    # Exclude all files related to licensing
    REGEX /license/ EXCLUDE
)

