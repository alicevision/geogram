#!/bin/sh
#
#      \\V (O |R |P /A |L |I |N |E 
# (C) Bruno Levy, INRIA - ALICE, 2012,2013
#
#   Confidential - proprietary software
#

export VORPALINE_BUILD_CONFIG="Debug"
export VORPALINE_SOURCE_DIR="/home/levy/Programming/Vorpaline/trunk"
export VORPALINE_BUILD_DIR="/home/levy/Programming/Vorpaline/trunk"
export VORPALINE_BIN_DIR="/home/levy/Programming/Vorpaline/trunk/bin"
export VORPALINE_LIB_DIR="/home/levy/Programming/Vorpaline/trunk/lib"

ulimit -c unlimited

args=
while [ -n "$1" ]; do
    case "$1" in
        --with-*=*)
            var=`echo "$1" | sed 's/--with-\([^=]*\)=\(.*\)$/VORPALINE_WITH_\U\1\E=\2/'`
            export "$var"
            shift
            ;;
        --with-*)
            var=`echo "$1" | sed 's/--with-\(.*\)$/VORPALINE_WITH_\U\1=1/'`
            export "$var"
            shift
            ;;
        *)
            args="$args $1"
            shift;
            ;;
    esac
done

