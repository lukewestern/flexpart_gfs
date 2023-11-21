#!/bin/bash
# By MB
# run FLEXPART, but provide a bit of information
echo "Welcome, running FLEXPART "
echo "Using defaults (/pathnames)"
cat /pathnames
echo "Mount volumes to change inputs"
echo "Git: $COMMIT"
echo "EXECUTING FLEXPART"
if [ $# -eq 1 ]; then
    echo "trying to execute: /src/$1"
    if [ -e /src/"$1" ]; then
        echo "Executing: /src/$1"
        /src/"$1" /pathnames
    else
        echo "Falling back to default, executing: /src/FLEXPART_ETA"
        /src/FLEXPART_ETA /pathnames
    fi
else
    echo "Executing: /src/FLEXPART_ETA"
    /src/FLEXPART_ETA /pathnames
fi
echo "FINISHED"
