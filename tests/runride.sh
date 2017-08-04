#!/bin/sh
# Wrapper for the RobotFramework editor

ulimit -c unlimited

scriptdir=`dirname "$0"`
. "$scriptdir/testenv.sh" || exit 1
exec ride.py $args

