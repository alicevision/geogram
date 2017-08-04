##############################################################################
# Geogram/Vorpaline root CMakeList
##############################################################################

# CMake 2.8.11 is required for 2 reasons:
# - it is the first version that fully supports the specification of Visual
# Studio toolsets (v110_xp).
# - it is the version that supports the command string(TIMESTAMP ...)
cmake_minimum_required(VERSION 2.8.11)

# Note: geogram.cmake defines GEOGRAM_WITH_VORPALINE
# that we could have used instead,
# but geogram.cmake needs to be included after the project() 
# command, since project() resets CFLAGS and CXXFLAGS.

if(IS_DIRECTORY ${CMAKE_SOURCE_DIR}/src/lib/vorpalib)
    project(Vorpaline)
else()
    project(Geogram)
endif()

# Optional modules
# (can be overriden in CMakeOptions.txt)

if(NOT DEFINED GEOGRAM_WITH_GRAPHICS) 
    set(GEOGRAM_WITH_GRAPHICS TRUE)
endif()    

include(cmake/geogram.cmake)

set(VORPALINE_VERSION_MAJOR 1)
set(VORPALINE_VERSION_MINOR 0)
set(VORPALINE_VERSION_PATCH 0)
set(VORPALINE_VERSION ${VORPALINE_VERSION_MAJOR}.${VORPALINE_VERSION_MINOR}.${VORPALINE_VERSION_PATCH})


# Determine the current Build-OS (Build-platform without the compiler info)
string(REGEX REPLACE "-[^-]+$" "" VORPALINE_OS ${VORPALINE_PLATFORM})

# Determine the current build date
string(TIMESTAMP VORPALINE_BUILD_DATE "%Y-%m-%d %H:%M:%S")
string(TIMESTAMP YEAR "%Y")

# Determine the current build number
# This is set by Jenkins in environment variable BUILD_NUMBER
set(VORPALINE_BUILD_NUMBER $ENV{BUILD_NUMBER})


##############################################################################
# Get SVN revision info

if(GEOGRAM_WITH_VORPALINE)
   find_package(Subversion QUIET)
   if(NOT SUBVERSION_FOUND)
       message(WARNING "Subversion executable not found - cannot determine current revision")
   else()
       Subversion_WC_INFO(${PROJECT_SOURCE_DIR} Vorpaline)
       message(STATUS "Vorpaline revision is ${Vorpaline_WC_REVISION}")
       set(VORPALINE_SVN_REVISION ${Vorpaline_WC_REVISION})
   endif()
endif()   

##############################################################################
# Geogram/Vorpaline sources

add_subdirectory(src/lib/geogram)
add_subdirectory(src/bin/vorpastat)
add_subdirectory(src/bin/vorpacomp)
add_subdirectory(src/bin/geodump)

if(GEOGRAM_WITH_FPG)
   add_subdirectory(src/bin/fpg)
endif()

if(GEOGRAM_WITH_VORPALINE)
   add_subdirectory(src/lib/vorpalib)
   add_subdirectory(src/bin/vorpaline)
else()
   add_subdirectory(src/bin/vorpalite)
endif()

if(GEOGRAM_WITH_GRAPHICS)
   add_subdirectory(src/lib/geogram_gfx)
   add_subdirectory(src/bin/vorpaview)
   add_subdirectory(src/bin/geobox)   
endif()

add_subdirectory(src/lib/third_party)
add_subdirectory(src/test)

if(GEOGRAM_WITH_GRAPHICS)
   add_subdirectory(src/demos)
endif()

add_subdirectory(doc)
add_subdirectory(tests)


##############################################################################
# Cleanup from previous builds

file(REMOVE ${CMAKE_BINARY_DIR}/doc/LICENSE.txt)


##############################################################################
# Vorpaline installation

# Configure CPack

set(CPACK_PACKAGE_NAME ${CMAKE_PROJECT_NAME})
set(CPACK_SYSTEM_NAME ${VORPALINE_OS})
set(CPACK_PACKAGE_VENDOR "INRIA - ALICE")
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "A flexible mesh generator")

set(CPACK_PACKAGE_VERSION_MAJOR ${VORPALINE_VERSION_MAJOR})
set(CPACK_PACKAGE_VERSION_MINOR ${VORPALINE_VERSION_MINOR})
set(CPACK_PACKAGE_VERSION_PATCH ${VORPALINE_VERSION_PATCH})
set(CPACK_PACKAGE_VERSION ${VORPALINE_VERSION})
set(CPACK_PACKAGE_FILE_NAME ${CPACK_PACKAGE_NAME}-${CPACK_PACKAGE_VERSION})
set(CPACK_COMPONENT_INCLUDE_TOPLEVEL_DIRECTORY true)

if(WIN32)
    set(CPACK_GENERATOR ZIP)
else()
    set(CPACK_GENERATOR TGZ)
endif()

# Enable component-based packaging for archive generators (TGZ, ZIP)
set(CPACK_ARCHIVE_COMPONENT_INSTALL ON)

set(CPACK_COMPONENTS_ALL runtime devkit devkit-full doc-devkit doc-devkit-full)
set(CPACK_COMPONENTS_GROUPING "IGNORE")

set(CPACK_COMPONENT_RUNTIME_DISPLAY_NAME "Vorpaline Application")
set(CPACK_COMPONENT_RUNTIME_GROUP "Runtime")

set(CPACK_COMPONENT_DEVKIT_DISPLAY_NAME "Vorpaline Developer Kit")
set(CPACK_COMPONENT_DEVKIT_GROUP "Development")

set(CPACK_COMPONENT_DEVKIT-FULL_DISPLAY_NAME "Vorpaline Full Developer Kit")
set(CPACK_COMPONENT_DEVKIT-FULL_GROUP "Development")

set(CPACK_COMPONENT_DOC-DEVKIT_DISPLAY_NAME "Vorpaline API Developer Kit Documentation")
set(CPACK_COMPONENT_DOC-DEVKIT_GROUP "Documentation")

set(CPACK_COMPONENT_DOC-DEVKIT-FULL_DISPLAY_NAME "Vorpaline Full Developer Kit Documentation")
set(CPACK_COMPONENT_DOC-DEVKIT-FULL_GROUP "Documentation")

set(CPACK_COMPONENT_DOC-DEVKIT-INTERNAL_DISPLAY_NAME "Vorpaline Internal Developer Kit Documentation")
set(CPACK_COMPONENT_DOC-DEVKIT-INTERNAL_GROUP "Documentation")

# Copy the helper script to build individual packages to the binary directory
configure_file(
    tools/make_package.pl.in
    make_package.pl
    @ONLY
)

# This must always be last!
include(CPack)

