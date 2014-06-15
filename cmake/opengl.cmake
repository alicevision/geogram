#
# Searches for OpenGL package
# Sets variable OPENGL_LIBRARIES to the list of GL and GLU libraries
#

find_package(OpenGL)
if(OPENGL_FOUND)
    include_directories(${OPENGL_INCLUDE_DIR})
else()
    if(WIN32)
	list(APPEND OPENGL_LIBRARIES opengl32)
    else()
	list(APPEND OPENGL_LIBRARIES GL)
    endif()
endif()

if(OPENGL_GLU_FOUND)
    include_directories(${GLU_INCLUDE_DIR})
else()
    if(WIN32)
        list(APPEND OPENGL_LIBRARIES glu32)
    else()
	list(APPEND OPENGL_LIBRARIES GLU)
    endif()
endif()

if(WIN32)
    # freeglut uses some timer functions defined in winmm
    list(APPEND OPENGL_LIBRARIES winmm)

    # Vorpaline's drag and drop support in freeglut requires shell32   
    list(APPEND OPENGL_LIBRARIES shell32)   
else()
    find_package(X11)
    list(APPEND OPENGL_LIBRARIES X11)
endif()

