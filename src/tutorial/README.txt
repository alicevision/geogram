Before compiling the tutorials, Geogram needs to be installed, either:

- to one of its defaults location:
   "/usr/local/"   (Linux, Mac OS/X)
   "Program Files" (Windows)

   This is the default if you build Geogram and use "make install" as root
(Linux, Mac OS/X) or build the "install" target from Visual C++ (Windows)
Windows note: to execute the "install" target, you will need to run
Visual C++ as administrator (else you are not allowed to write in
"Program Files").

- or it can be installed to a custom location, but then
   the environment variable GEOGRAM_INSTALL_PREFIX needs
   to be set.

To compile a tutorial:

Linux, Mac/OSX:
===============

  You got two alternatives, either using CMake (and the FindGeogram.cmake
file provided in each example), or using pkg-config (that uses the .pc
files installed in the lib/pkgconfig subdirectory of GEOGRAM_INSTALL_ROOT,
i.e. /usr/local/lib/pkgconfig by default)


 Using cmake
 -----------
  cd tutorials/<tutorial_name>
  cmake .
  make

 Using pkg-config
 ----------------
  cd tutorials/<tutorial_name>
  ./make_it.sh
  

Windows:
========
  start cmake_gui
  set sources directory as tutorials/<tutorial_name>
  push the "configure" button
  push the "generate" button
  open <tutorial_name>.dsw in Visual C++
  select "Release" (default is "Debug")
  build
  
  
