
vor_reset_warning_level()

if(WIN32)
#  TODO: there are probably some MSVC warnings to 
#  deactivate here...
else()
   add_definitions(-Wno-unused -Wno-sign-compare -Wno-strict-aliasing)
endif()

if(GEOGRAM_WITH_GRAPHICS)
   if(NOT GEOGRAM_WITH_EMSCRIPTEN)
     add_subdirectory(glfw)
   endif()  
endif()

