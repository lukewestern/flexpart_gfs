#!/bin/bash
# By MB
# run default tests
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m'
MM='FLEXPART AUTOMATIC OPTIONS TEST'

pause() {
    read -n1 -rsp $'Press any key to continue or Ctrl+C to exit...\n'
}

warning() {
    printf "%-68s[$YELLOW%10s$NC]\n" "$@" "SKIPPED"
    return 0
}

report() {
    if [ $? -eq 0 ]; then
        printf "%-68s[$GREEN%10s$NC]\n" "$@" "OK"
        return 0
    else
        printf "%-68s[$RED%10s$NC]\n" "$@" "FAILED"
        return 1
    fi
}

#
# initial conditions
#
test -f ../src/FLEXPART
report "[$MM] executable: ../src/FLEXPART" || exit 1
ln -s ../src/FLEXPART .
test -d ./default_options
report "[$MM] default options: ./default_options" || exit 1
cp -f ./default_options ./current
mkdir ./output/
#
# Different options tests
#
STATUS=0
sed "/LOUTRESTART=/c\ LOUTRESTART=   -1," ./default_options/COMMAND > ./current/COMMAND
./FLEXPART pathnames
STATUS=$($STATUS + $?)
exit $STATUS