# Install script for directory: /home/levy/Programming/Vorpaline/trunk/src/lib/geogram_gfx

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/usr/local")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "Debug")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Install shared libraries without execute permission?
if(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)
  set(CMAKE_INSTALL_SO_NO_EXE "1")
endif()

if(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "devkit")
  foreach(file
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libgeogram_gfx.so.1.3.4"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libgeogram_gfx.so.1"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libgeogram_gfx.so"
      )
    if(EXISTS "${file}" AND
       NOT IS_SYMLINK "${file}")
      file(RPATH_CHECK
           FILE "${file}"
           RPATH "/usr/local/lib")
    endif()
  endforeach()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE SHARED_LIBRARY FILES
    "/home/levy/Programming/Vorpaline/trunk/lib/libgeogram_gfx.so.1.3.4"
    "/home/levy/Programming/Vorpaline/trunk/lib/libgeogram_gfx.so.1"
    "/home/levy/Programming/Vorpaline/trunk/lib/libgeogram_gfx.so"
    )
  foreach(file
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libgeogram_gfx.so.1.3.4"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libgeogram_gfx.so.1"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libgeogram_gfx.so"
      )
    if(EXISTS "${file}" AND
       NOT IS_SYMLINK "${file}")
      file(RPATH_CHANGE
           FILE "${file}"
           OLD_RPATH "/home/levy/Programming/Vorpaline/trunk/lib:"
           NEW_RPATH "/usr/local/lib")
      if(CMAKE_INSTALL_DO_STRIP)
        execute_process(COMMAND "/usr/bin/strip" "${file}")
      endif()
    endif()
  endforeach()
endif()

if(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "devkit-full")
  foreach(file
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libgeogram_gfx.so.1.3.4"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libgeogram_gfx.so.1"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libgeogram_gfx.so"
      )
    if(EXISTS "${file}" AND
       NOT IS_SYMLINK "${file}")
      file(RPATH_CHECK
           FILE "${file}"
           RPATH "/usr/local/lib")
    endif()
  endforeach()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE SHARED_LIBRARY FILES
    "/home/levy/Programming/Vorpaline/trunk/lib/libgeogram_gfx.so.1.3.4"
    "/home/levy/Programming/Vorpaline/trunk/lib/libgeogram_gfx.so.1"
    "/home/levy/Programming/Vorpaline/trunk/lib/libgeogram_gfx.so"
    )
  foreach(file
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libgeogram_gfx.so.1.3.4"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libgeogram_gfx.so.1"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libgeogram_gfx.so"
      )
    if(EXISTS "${file}" AND
       NOT IS_SYMLINK "${file}")
      file(RPATH_CHANGE
           FILE "${file}"
           OLD_RPATH "/home/levy/Programming/Vorpaline/trunk/lib:"
           NEW_RPATH "/usr/local/lib")
      if(CMAKE_INSTALL_DO_STRIP)
        execute_process(COMMAND "/usr/bin/strip" "${file}")
      endif()
    endif()
  endforeach()
endif()

if(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "devkit")
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/geogram1/geogram_gfx" TYPE DIRECTORY FILES "/home/levy/Programming/Vorpaline/trunk/src/lib/geogram_gfx/." FILES_MATCHING REGEX "/[^/]*\\.h$")
endif()

if(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "devkit-full")
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/geogram1/geogram_gfx" TYPE DIRECTORY FILES "/home/levy/Programming/Vorpaline/trunk/src/lib/geogram_gfx/." FILES_MATCHING REGEX "/[^/]*\\.h$")
endif()

if(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/pkgconfig" TYPE FILE FILES "/home/levy/Programming/Vorpaline/trunk/geogram_gfx.pc")
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for each subdirectory.
  include("/home/levy/Programming/Vorpaline/trunk/src/lib/geogram_gfx/third_party/cmake_install.cmake")

endif()

