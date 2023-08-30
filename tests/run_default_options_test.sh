#!/bin/bash
# By MB
# run default tests
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m'
MM='FLEXPART AUTOMATIC OPTIONS TEST'

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
warning "[$MM] $PWD"
#
# Change to directory of this script
#
cd $(dirname $0)
#
# Check for Flexpart executable build before
#
test -f ../src/FLEXPART
report "[$MM] executable: ../src/FLEXPART" || exit 1
ln -s ../src/FLEXPART .
test -d ./default_options
report "[$MM] default options: ./default_options" || exit 1
cp -rf ./default_options ./current
mkdir -p ./output/
#
# Different options tests
#
STATUS=0
TESTSRUN=0
#
# run Options test with restart output disabled
#
sed "/LOUTRESTART=/c\ LOUTRESTART=   -1," ./default_options/COMMAND > ./current/COMMAND
./FLEXPART pathnames
report "[$MM] TEST $TESTRUN (LOUTRESTART)"
STATUS=$((STATUS + $?))
TESTSRUN=$((TESTSRUN + 1))
# clean up
rm -rf ./current ./output/*
#
# 
#BACKWARD WET DEPOSITION
cp -rf ./default_options ./current
sed "/LDIRECT=/c\ LDIRECT=   -1," ./default_options/COMMAND > ./current/COMMAND
# change release
# 
# IND_RECEPTOR=          1, ! Unit to be used at the receptor; [0]no receptor [1]mass 2]mass mixing ratio 3]wet depo. 4]dry depo.
# ! Release start time in UTC HHMISS: HH hours, MI=minutes, SS=seconds
sed -i "s/ITIME1.*/ITIME1  =   030000,/" ./current/RELEASES
sed -i "s/ITIME2.*/ITIME2  =   030000,/" ./current/RELEASES
sed -i "s/IND_RECEPTOR.*/IND_RECEPTOR=  3," ./current/COMMAND
./FLEXPART pathnames
report "[$MM] TEST $TESTRUN (IND_RECEPTOR=3)"
STATUS=$((STATUS + $?))
TESTSRUN=$((TESTSRUN + 1))
# clean up
rm -rf ./current ./output/*
#
# 
#BACKWARD DRY DEPOSITION
cp -rf ./default_options ./current
sed "/LDIRECT=/c\ LDIRECT=   -1," ./default_options/COMMAND > ./current/COMMAND
sed -i "s/ITIME1.*/ITIME1  =   030000,/" ./current/RELEASES
sed -i "s/ITIME2.*/ITIME2  =   030000,/" ./current/RELEASES
sed -i "/IND_RECEPTOR.*/IND_RECEPTOR=  4,"./current/COMMAND
./FLEXPART pathnames
report "[$MM] TEST $TESTRUN (IND_RECEPTOR=4)"
STATUS=$((STATUS + $?))
TESTSRUN=$((TESTSRUN + 1))
# clean up
rm -rf ./current ./output/*
#
#
#OUTRESTART
cp -rf ./default_options ./current
sed -i "/LOUTRESTART.*/LOUTRESTART=  3600,"./current/COMMAND
./FLEXPART pathnames
report "[$MM] TEST $TESTRUN (LOUTRESTART=3600)"
STATUS=$((STATUS + $?))
TESTSRUN=$((TESTSRUN + 1))
# clean up
rm -rf ./current ./output/*
#
#
#CTL option
cp -rf ./default_options ./current
sed -i "/CTL.*/CTL=  10.0000,"./current/COMMAND
./FLEXPART pathnames
report "[$MM] TEST $TESTRUN (CTL=10)"
STATUS=$((STATUS + $?))
TESTSRUN=$((TESTSRUN + 1))
# clean up
rm -rf ./current ./output/*
#
#
#NETCDF output
cp -rf ./default_options ./current
sed -i "/IOUT.*/IOUT=  9,"./current/COMMAND
./FLEXPART pathnames
report "[$MM] TEST $TESTRUN (IOUT=9)"
STATUS=$((STATUS + $?))
TESTSRUN=$((TESTSRUN + 1))
# clean up
rm -rf ./current ./output/*
#
#
#NETCDF particle output
cp -rf ./default_options ./current
sed -i "/IPOUT.*/IPOUT=  1,"./current/COMMAND
./FLEXPART pathnames
report "[$MM] TEST $TESTRUN (IPOUT=1)"
STATUS=$((STATUS + $?))
TESTSRUN=$((TESTSRUN + 1))
# clean up
rm -rf ./current ./output/*
#
#
#NETCDF particle output at the end
cp -rf ./default_options ./current
sed -i "/IPOUT.*/IPOUT=  2,"./current/COMMAND
./FLEXPART pathnames
report "[$MM] TEST $TESTRUN (IPOUT=2)"
STATUS=$((STATUS + $?))
TESTSRUN=$((TESTSRUN + 1))
# clean up
rm -rf ./current ./output/*
#
#
#LSUBGRID
cp -rf ./default_options ./current
sed -i "/LSUBGRID.*/LSUBGRID=  1,"./current/COMMAND
./FLEXPART pathnames
report "[$MM] TEST $TESTRUN (LSUBGRID=1)"
STATUS=$((STATUS + $?))
TESTSRUN=$((TESTSRUN + 1))
# clean up
rm -rf ./current ./output/*
#
#
#AGESPECTRA
cp -rf ./default_options ./current
sed -i "/LAGESPECTRA.*/LAGESPECTRA=  1,"./current/COMMAND
./FLEXPART pathnames
report "[$MM] TEST $TESTRUN (LAGESPECTRA=1)"
STATUS=$((STATUS + $?))
TESTSRUN=$((TESTSRUN + 1))
# clean up
rm -rf ./current ./output/*
#
#
#IPIN
./FLEXPART pathnames
mv output/restart_20090101020000 output/restart.bin
cp -rf ./default_options ./current
sed -i "/IPIN.*/IPIN=  1,"./current/COMMAND
./FLEXPART pathnames
report "[$MM] TEST $TESTRUN (IPIN=1)"
STATUS=$((STATUS + $?))
TESTSRUN=$((TESTSRUN + 1))
# clean up
rm -rf ./current ./output/*
#
#
#IOUTPUTFOREACHRELEASE
cp -rf ./default_options ./current
sed -i "/IOUT.*/IOUT=  9,"./current/COMMAND
sed -i "/IOUTPUTFOREACHRELEASE.*/IOUTPUTFOREACHRELEASE=  1,"./current/COMMAND
./FLEXPART pathnames
report "[$MM] TEST $TESTRUN (IOUTPUTFOREACHRELEASE=1)"
STATUS=$((STATUS + $?))
TESTSRUN=$((TESTSRUN + 1))
#
#
#IFLUX
cp -rf ./default_options ./current
sed -i "/IFLUX.*/IFLUX=  1,"./current/COMMAND
./FLEXPART pathnames
report "[$MM] TEST $TESTRUN (IFLUX=1)"
STATUS=$((STATUS + $?))
TESTSRUN=$((TESTSRUN + 1))
# clean up
rm -rf ./current ./output/*
#
#
#MDOMAINFILL
cp -rf ./default_options ./current
sed -i "/MDOMAINFILL.*/MDOMAINFILL=  1,"./current/COMMAND
./FLEXPART pathnames
report "[$MM] TEST $TESTRUN (MDOMAINFILL=1)"
STATUS=$((STATUS + $?))
TESTSRUN=$((TESTSRUN + 1))
# clean up
rm -rf ./current ./output/*
#
#
#MQUASILAG
cp -rf ./default_options ./current
sed -i "/MQUASILAG.*/MQUASILAG=  1,"./current/COMMAND
./FLEXPART pathnames
report "[$MM] TEST $TESTRUN (MQUASILAG=1)"
STATUS=$((STATUS + $?))
TESTSRUN=$((TESTSRUN + 1))
# clean up
rm -rf ./current ./output/*
#
#
#NESTED_OUTPUT
cp -rf ./default_options ./current
sed -i "/NESTED_OUTPUT.*/NESTED_OUTPUT=  1,"./current/COMMAND
./FLEXPART pathnames
report "[$MM] TEST $TESTRUN (NESTED_OUTPUT=1)"
STATUS=$((STATUS + $?))
TESTSRUN=$((TESTSRUN + 1))
# clean up
rm -rf ./current ./output/*
#
#
#NESTED_INPUT
#
#
#LINIT_COND
cp -rf ./default_options ./current
sed -i "/LINIT_COND.*/LINIT_COND=  1,"./current/COMMAND
./FLEXPART pathnames
report "[$MM] TEST $TESTRUN (LINIT_COND=1)"
STATUS=$((STATUS + $?))
TESTSRUN=$((TESTSRUN + 1))
# clean up
rm -rf ./current ./output/*
#
#
#CBLFLAG
cp -rf ./default_options ./current
sed -i "/CTL.*/CTL=  10.0000,"./current/COMMAND
sed -i "/IFINE.*/IFINE=  10,"./current/COMMAND
sed -i "/CBLFLAG.*/CBLFLAG=  1,"./current/COMMAND
./FLEXPART pathnames
report "[$MM] TEST $TESTRUN (CBLFLAG=1)"
STATUS=$((STATUS + $?))
TESTSRUN=$((TESTSRUN + 1))
# clean up
rm -rf ./current ./output/*
# Test particale output (ipout=1)
# between openmp=1 and openmp=?? the locations should be identical
# no convection no turbulence
#

#
# Add comaprison with MASTER branch outputs
# add as a volume
#
#

#
# FINAL
#
echo "[$MM] Tests failed: $STATUS / $TESTSRUN"
#
# Return collective error status
#
exit $STATUS
