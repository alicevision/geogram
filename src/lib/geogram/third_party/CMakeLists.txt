vor_reset_warning_level()

# Since we are in a CMake "object library" (i.e.,
# a set of object files compiled with same options),
# the automatically generated geogram_gfx_EXPORTS
# pre-processor definition of the parent library
# is not automatically inherited (this is why
# I set it manually)
if(BUILD_SHARED_LIBS)
   add_definitions(-Dgeogram_EXPORTS)
endif()

aux_source_directories(SOURCES "Source Files"           .)
aux_source_directories(SOURCES "Source Files\\HLBFGS"   HLBFGS)
aux_source_directories(SOURCES "Source Files\\LM6"      LM6)
aux_source_directories(SOURCES "Source Files\\rply"     rply)
aux_source_directories(SOURCES "Source Files\\shewchuk" shewchuk)
aux_source_directories(SOURCES "Source Files\\zlib"     zlib)
aux_source_directories(SOURCES "Source Files\\gzstream" gzstream)

if(GEOGRAM_WITH_TETGEN)
    aux_source_directories(SOURCES "Source Files\\tetgen" tetgen)
endif()

if(GEOGRAM_WITH_TRIANGLE)
    aux_source_directories(SOURCES "Source Files\\triangle" triangle)
endif()

add_library(geogram_third_party OBJECT ${SOURCES})

