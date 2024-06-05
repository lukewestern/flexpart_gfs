
module purge
module load eccodes
export CPATH=/opt/tools/eccodes/2.28.0-gcc9.4/include:/opt/tools/netcdf-fortran/4.6.0-gcc9.4/include
export LIBRARY_PATH=$LIBRARY_PATH:/opt/tools/eccodes/2.28.0-gcc9.4/lib

rm -f *.o *.mod

# eta ncf - new default FLEXPART_ETA
FC=/opt/tools/gcc/9.4.0/bin/gfortran  make -f makefile_gfortran
#executable: FLEXPART_ETA

# meter bin - classic
#FC=/opt/tools/gcc/9.4.0/bin/gfortran  make -f makefile_gfortran eta=no ncf=no
# executable: FLEXPART_BIN
 
# meter ncf:
#FC=/opt/tools/gcc/9.4.0/bin/gfortran make -f makefile_gfortran eta=no
# executable: FLEXPART
 
# eta/bin
#FC=/opt/tools/gcc/9.4.0/bin/gfortran  make -f makefile_gfortran ncf=no
# executable: FLEXPART_ETA_BIN
