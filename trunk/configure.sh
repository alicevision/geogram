#!/bin/sh

# This file for Linux users, 
# launches CMake and creates configuration for
# Release and Debug modes.


echo
echo ============= Checking for CMake ============
echo

if (cmake --version); then
   echo "Found CMake"
else
   echo "Error: CMake not found, please install it (see http://www.cmake.org/)"
   exit 1
fi

# Parse command line arguments

cmake_options=
while [ -n "$1" ]; do
    case "$1" in
        --debug)
            cmake_options="$cmake_options --trace"
            shift
            ;;

        --with-*=*)
            cmake_option=`echo "$1" | sed 's/--with-\([^=]*\)=\(.*\)$/-DGEOGRAM_WITH_\U\1\E:STRING="\2"/'`
            cmake_options="$cmake_options $cmake_option"
            shift
            ;;

        --with-*)
            cmake_option=`echo "$1" | sed 's/--with-\(.*\)$/-DGEOGRAM_WITH_\U\1:BOOL=TRUE/'`
            cmake_options="$cmake_options $cmake_option"
            shift
            ;;

        --help)
            cat <<END
NAME
    configure.sh

SYNOPSIS
    Prepares the build environment for Geogram.
    
    - For Unix builds, the script creates 2 build trees for Debug and Release
    build in a 'build' sub directory under the project root.

    - For Windows builds, the script creates a single build tree that supports
    all cmake build types (Debug, Release, RelWithDebInfo, MinSizeRel)
    build in a 'build' sub directory under the project root.

USAGE
    configure.sh [options] build-platform

OPTIONS

    --help
        Prints this page.

    --with-gcov
        Builds the project for coverage analysis with gcov    

    --with-gprof
        Builds the project for performance analysis with gprof

    --with-asan
        Builds the project with Google's AddressSanitizer (dynamic builds only)
        See: http://code.google.com/p/address-sanitizer/

    --with-tsan
        Builds the project with Google's ThreadSanitizer (dynamic builds only)
        See: https://code.google.com/p/thread-sanitizer/

    --with-ddt=ddt-root-dir
        Builds the project for memory analysis with Allinea's DDT installed in
        the specified directory: ddt-root-dir

PLATFORM
    Build platform supported by Geogram. See cmake/platforms for supported
    platforms.
END
            exit
            ;;
            
        -*)
            echo "Error: unrecognized option: $1"
            return 1
            ;;
        *)
            break;
            ;;
    esac
done

# Check the current OS

os="$1"
if [ -z "$os" ]; then
    os=`uname -a`
    case "$os" in
        Linux*x86_64*)
            os=Linux64-gcc
            ;;
        Linux*amd64*)
            os=Linux64-gcc
            ;;
        Linux*i586*|Linux*i686*)
            os=Linux32-gcc
            ;;
        *)
            echo "Error: OS not supported: $os"
            exit 1
            ;;
    esac
fi

#  Import plaform specific environment

. cmake/platforms/$os/setvars.sh || exit 1

# Generate the Makefiles

for config in Release Debug; do
   platform=$os-$config
   echo
   echo ============= Creating makefiles for $platform ============
   echo

   build_dir=build/$platform
   mkdir -p $build_dir
   (cd $build_dir; cmake -DCMAKE_BUILD_TYPE:STRING=$config -DGEOGRAM_PLATFORM:STRING=$os $cmake_options ../../)
done

echo
echo ============== Geogram build configured ==================
echo

cat << EOF
To build geogram:
  - go to build/$os-Release or build/$os-Debug
  - run 'make' or 'cmake --build .'

Note: local configuration can be specified in CMakeOptions.txt
(see CMakeOptions.txt.sample for an example)
You'll need to re-run configure.sh if you create or modify CMakeOptions.txt

EOF
