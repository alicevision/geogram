##############################################################################
# Geogram root CMakeList
##############################################################################

# Make sure that GEOGRAM_PLATFORM is defined

if(NOT DEFINED GEOGRAM_PLATFORM)
    message(FATAL_ERROR
        " CMake variable GEOGRAM_PLATFORM is not defined.\n"
        " Please run configure.{sh,bat} to setup the build tree."
    )
endif()

# CMake 2.8.11 is required for 2 reasons:
# - it is the first version that fully supports the specification of Visual
# Studio toolsets (v110_xp).
# - it is the version that supports the command string(TIMESTAMP ...)

cmake_minimum_required(VERSION 2.8.11)

project(Geogram)

set(GEOGRAM_VERSION_MAJOR 0)
set(GEOGRAM_VERSION_MINOR 9)
set(GEOGRAM_VERSION_PATCH 0)
set(GEOGRAM_VERSION ${GEOGRAM_VERSION_MAJOR}.${GEOGRAM_VERSION_MINOR}.${GEOGRAM_VERSION_PATCH})

# Determine the current Build-OS (Build-platform without the compiler info)
string(REGEX REPLACE "-[^-]+$" "" GEOGRAM_OS ${GEOGRAM_PLATFORM})

# Determine the current build date
string(TIMESTAMP GEOGRAM_BUILD_DATE "%Y-%m-%d %H:%M:%S")
string(TIMESTAMP YEAR "%Y")

# Determine the current build number
# This is set by Jenkins in environment variable BUILD_NUMBER
set(GEOGRAM_BUILD_NUMBER $ENV{BUILD_NUMBER})

##############################################################################
# Get SVN revision info

find_package(Subversion QUIET)
if(NOT SUBVERSION_FOUND)
    message(WARNING "Subversion executable not found - cannot determine current revision")
else()
    Subversion_WC_INFO(${PROJECT_SOURCE_DIR} Geogram)
    message(STATUS "Geogram revision is ${Geogram_WC_REVISION}")
    set(GEOGRAM_SVN_REVISION ${Geogram_WC_REVISION})
endif()

##############################################################################
# Build configuration

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)

if(EXISTS ${CMAKE_SOURCE_DIR}/CMakeOptions.txt)
   message(STATUS "Using local options file: ${CMAKE_SOURCE_DIR}/CMakeOptions.txt")
   include(${CMAKE_SOURCE_DIR}/CMakeOptions.txt)
endif()

include(${CMAKE_SOURCE_DIR}/cmake/utilities.cmake)
include(${CMAKE_SOURCE_DIR}/cmake/platforms/${GEOGRAM_PLATFORM}/config.cmake)

# This test is there to keep CMake happy about unused variable CMAKE_BUILD_TYPE
if(CMAKE_BUILD_TYPE STREQUAL "")
endif()

# Static versus dynamic builds
if(GEOGRAM_BUILD_DYNAMIC)
    set(BUILD_SHARED_LIBS TRUE)
    # Object files in OBJECT libraries are compiled in static mode, event if
    # BUILD_SHARED_LIBS is true! We must set CMAKE_POSITION_INDEPENDENT_CODE
    # to force compilation in dynamic mode.
    set(CMAKE_POSITION_INDEPENDENT_CODE TRUE)
else()
    set(BUILD_SHARED_LIBS FALSE)
endif()

##############################################################################
# Geogram sources

include_directories(${CMAKE_SOURCE_DIR}/src/lib)

add_subdirectory(src/lib/geogram)
add_subdirectory(src/bin/fpg)
add_subdirectory(src/bin/delaunay3d)

#add_subdirectory(doc)
#add_subdirectory(tests)

##############################################################################
# Cleanup from previous builds

file(REMOVE ${CMAKE_BINARY_DIR}/doc/LICENSE.txt)


##############################################################################
# Geogram installation

# Configure CPack

set(CPACK_SYSTEM_NAME ${GEOGRAM_OS})
SET(CPACK_PACKAGE_VENDOR "INRIA - ALICE")

set(CPACK_PACKAGE_VERSION_MAJOR ${GEOGRAM_VERSION_MAJOR})
set(CPACK_PACKAGE_VERSION_MINOR ${GEOGRAM_VERSION_MINOR})
set(CPACK_PACKAGE_VERSION_PATCH ${GEOGRAM_VERSION_PATCH})
set(CPACK_PACKAGE_VERSION ${GEOGRAM_VERSION})

if(WIN32)
    set(CPACK_GENERATOR ZIP)
else()
    set(CPACK_GENERATOR TGZ)
endif()

# Enable composant-based packaging for archive generators (TGZ, ZIP)
set(CPACK_ARCHIVE_COMPONENT_INSTALL ON)

set(CPACK_COMPONENTS_ALL runtime devkit devkit-full)
set(CPACK_COMPONENTS_GROUPING "IGNORE")

set(CPACK_COMPONENT_RUNTIME_DISPLAY_NAME "Geogram Application")
set(CPACK_COMPONENT_RUNTIME_GROUP "Runtime")

set(CPACK_COMPONENT_DEVKIT_DISPLAY_NAME "Geogram Developer Kit")
set(CPACK_COMPONENT_DEVKIT_GROUP "Development")

set(CPACK_COMPONENT_DEVKIT-FULL_DISPLAY_NAME "Geogram Full Developer Kit")
set(CPACK_COMPONENT_DEVKIT-FULL_GROUP "Development")

# This must always be last!
include(CPack)

