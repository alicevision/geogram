# MakeMake.sh: generates a Makefile for compiling with MinGW
#
# Create geogram/build/MinGW directory
# cd geogram/build/MinGW
# ../../tools/MakeMake.sh > Makefile
# make

DATE=$(date)
SRCDIR=../../src/lib

CMAKELIST=$SRCDIR/../../CMakeLists.txt
VERSION_MAJOR=$(cat $CMAKELIST | grep 'set(VORPALINE_VERSION_MAJOR' | sed -e 's|[^0-9]||g')
VERSION_MINOR=$(cat $CMAKELIST | grep 'set(VORPALINE_VERSION_MINOR' | sed -e 's|[^0-9]||g')
VERSION_PATCH=$(cat $CMAKELIST | grep 'set(VORPALINE_VERSION_PATCH' | sed -e 's|[^0-9]||g')
VERSION=$VERSION_MAJOR.$VERSION_MINOR.$VERSION_PATCH

mkdir -p geogram
cat > geogram/version.h << EOF
#ifndef GEOGRAM_BASIC_VERSION
#define GEOGRAM_BASIC_VERSION

#define VORPALINE_VERSION_MAJOR "$VERSION_MAJOR"
#define VORPALINE_VERSION_MINOR "$VERSION_MINOR"
#define VORPALINE_VERSION_PATCH "$VERSION_PATCH"
#define VORPALINE_VERSION "$VERSION"
#define VORPALINE_BUILD_NUMBER ""
#define VORPALINE_BUILD_DATE "$DATE"
#define VORPALINE_SVN_REVISION "????"

#endif
EOF

cat << EOF
# This Makefile was automatically generated using
# geogram's MakeMake utility.

SRCDIR=$SRCDIR
CC=x86_64-w64-mingw32-gcc-win32
CXX=x86_64-w64-mingw32-g++-win32
AR=x86_64-w64-mingw32-ar
RANLIB=x86_64-w64-mingw32-ranlib

# Uncomment for DLL build (untested), see also 
# end of this file.
# EXPORTS=-Dgeogram_EXPORTS -Dexploragram_EXPORTS
EXPORTS=

COPT=-I\$(SRCDIR) -I. \$(EXPORTS)
CXXOPT=-I\$(SRCDIR) -I. -Wno-deprecated -std=c++11 \$(EXPORTS) 
LDXX=\$(CXX)
LDXXOPT=-static 

all: vorpalite_static.exe
EOF

CSRC=$(find $SRCDIR/geogram $SRCDIR/exploragram  -name "*.c" -print)
CPPSRC=$(find $SRCDIR/geogram $SRCDIR/exploragram -name "*.cpp" -print)

SRC2OBJ="sed -e 's|$SRCDIR/||g' -e 's|/|_|g' -e 's|\.cpp|\.o|g' -e 's|\.c|\.o|g'"

OBJ=$(echo $CSRC $CPPSRC | eval $SRC2OBJ)

echo
echo OBJ=$OBJ
echo

for i in $CSRC 
do
   obj=$(echo $i | eval $SRC2OBJ)
   echo $obj: $i
cat << EOF
	\$(CC) -c \$(COPT) $i -o $obj

EOF
done

for i in $CPPSRC 
do
   obj=$(echo $i | eval $SRC2OBJ)
   echo $obj: $i
cat << EOF
	\$(CXX) -c \$(CXXOPT) $i -o $obj

EOF
done

#============================================================================================

cat <<EOF

libexplogeogram_static.a: \$(OBJ)
	\$(AR) cq libexplogeogram_static.a \$(OBJ)
	\$(RANLIB) libexplogeogram_static.a

vorpalite.o: \$(SRCDIR)/../bin/vorpalite/main.cpp
	\$(CXX) -c \$(CXXOPT)  \$(SRCDIR)/../bin/vorpalite/main.cpp -o vorpalite.o

vorpalite_static.exe: vorpalite.o libexplogeogram_static.a
	\$(LDXX) \$(LDXXOPT) -o vorpalite_static.exe vorpalite.o libexplogeogram_static.a

# DLL build (untested)
#
#libexplogeogram_dll.a: \$(OBJ)
#	\$(CXX) -shared -o libexplogeogram_dll.dll \$(OBJ) -Wl,--out-implib,libexplogeogram_dll.a
#
#vorpalite_dynamic.exe: vorpalite.o libexplogeogram_dll.a
#	\$(CXX) \$(CXXOPT) -o vorpalite_dynamic.exe vorpalite.o libexplogeogram_dll.a

EOF
