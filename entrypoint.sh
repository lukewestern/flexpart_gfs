#!/bin/bash
# By MB
# run FLEXPART, but provide a bit of information

echo "Welcome, running FLEXPART "
echo "Using defaults:"
cat /pathnames
echo "Mount volumes to change inputs"
echo "Git: $COMMIT"
echo "EXECUTING FLEXPART"
/src/FLEXPART /pathnames
echo "FINISHED"