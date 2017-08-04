
add_subdirectory(vorpastat)
add_subdirectory(vorpacomp)
add_subdirectory(geodump)
add_subdirectory(geocod)

if(GEOGRAM_WITH_GRAPHICS)
  add_subdirectory(vorpaview)
  add_subdirectory(geobox)
endif()

if(GEOGRAM_WITH_VORPALINE) 
  add_subdirectory(vorpaline)
else()
  add_subdirectory(vorpalite)
endif()

if(GEOGRAM_WITH_FPG)
  add_subdirectory(fpg)
endif()
