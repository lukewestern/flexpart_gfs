SHELL = /bin/bash
################################################################################
#  DESCRIPTION
#    Makefile for FLEXPART. Standard (serial) and parallel (MPI) version
#
#    Dependencies are resolved in this makefile, so parallel make is
#    possible ("make -j")
#
#    At NILU we have installed gcc-4.9.1 and libraries under user /homevip/flexpart
#    ("ROOT_DIR")
#    To use gfortran version 4.9, add "gcc=4.9" to the make command, e.g.
#       'make -j ecmwf gcc=4.9',
#    also set environment variable LD_LIBRARY_PATH to point to compiler libraries
#
#    Makefile was modified to produce unified executable for both ECMWF and GFS meteo data formats
#    gributils were included to detect format of meteo data
#
#    Cpp directives USE_MPIINPLACE were added to three source files. The effect of these directives 
#    are to enable the MPI_IN_PLACE option only if compiled with a -DUSE_MPIINPLACE directive. 
#    Otherwise, a safer option (which requires the allocation of another array) is used by default. 
#    In makefile added the -x f95-cpp-input flag for compiling of cpp directives.
# 
#  USAGE
#    Compile serial FLEXPART 
#      make [-j] serial
#
#  NETCDF OUTPUT
#    To remove support for output in netCDF format, append `ncf=no` to the
#    `make` command
# 
################################################################################

## PROGRAMS
# Unified executable names
# The same executable is used for both ECMWF and GFS metdata

# Serial processing executable
FLEXPART-SERIAL = FLEXPART

export OMP_NESTED=TRUE

#export TAU_MAKEFILE=/jetfs/home/lbakels/Codes/tau/x86_64/lib/Makefile.tau-pdt-openmp-opari
#export TAU_OPTIONS='-optCompInst -optRevert'
#export PROFILEDIR=/jetfs/home/lbakels/Profile2/
ifeq ($(shell hostname), jet01.img.univie.ac.at)
        export INTEL_LICENSE_FILE = /opt/intel/serverlicenses/l_5WHWBJ6D.lic
endif

        #ROOT_DIR = /homevip/flexpart/
        ROOT_DIR = . #MD

        F90 = /opt/intel/compilers_and_libraries_2020/linux/bin/intel64/ifort#/jetfs/home/lbakels/Codes/tau/x86_64/bin/tau_f90.sh#/opt/intel/compilers_and_libraries_2020/linux/bin/intel64/ifort

        #MD for jet

        INCPATH1 = /jetfs/spack/opt/spack/linux-rhel8-skylake_avx512/intel-20.0.2/eccodes-2.18.0-6tadpgreot7jf4yoaiqmqueiihhdcsxk/include
        LIBPATH1 = /jetfs/spack/opt/spack/linux-rhel8-skylake_avx512/intel-20.0.2/eccodes-2.18.0-6tadpgreot7jf4yoaiqmqueiihhdcsxk/lib
        INCPATH2 = /jetfs/spack/opt/spack/linux-rhel8-skylake_avx512/intel-20.0.2/netcdf-fortran-4.5.3-irdm5gqccsig2om7jqu376h6tlebh4bc/include
        LIBPATH2 = /jetfs/spack/opt/spack/linux-rhel8-skylake_avx512/intel-20.0.2/netcdf-fortran-4.5.3-irdm5gqccsig2om7jqu376h6tlebh4bc/lib
        LIBPATH3  = 

### Enable netCDF output?
ifeq ($(ncf), no)
	NCOPT = -UUSE_NCF
else
	NCOPT = -DUSE_NCF -lnetcdff
endif

# path to gributils used to detect meteodata format
VPATH = gributils/

## OPTIMIZATION LEVEL
O_LEV = 2 # [0,1,2,3,g,s,fast]
O_LEV_DBG = g # [0,g]
FUSER = -qopenmp -xSKYLAKE-AVX512#-traceback -check bounds #-xSKYLAKE-AVX512
## LIBRARIES
LIBS = -leccodes -leccodes_f90 -lm -ljasper $(NCOPT) #MD

FFLAGS   = -I$(INCPATH1) -I$(INCPATH2) -O$(O_LEV) -g -m64 -fpp -convert little_endian -assume byterecl -mcmodel=large -O$(O_LEV) $(NCOPT) $(FUSER)  #-convert little_endian -fpp -WB -check all

DBGFLAGS = -I$(INCPATH1) -I$(INCPATH2) -O$(O_LEV_DBG) -g3 -ggdb3 -cpp -m64 -mcmodel=large -fconvert=little-endian -frecord-marker=4 -fmessage-length=0 -flto=jobserver -O$(O_LEV_DBG) $(NCOPT) -fbacktrace   -Wall  -fdump-core $(FUSER)  #  -ffpe-trap=invalid,overflow,denormal,underflow,zero  -Warray-bounds -fcheck=all

LDFLAGS  = $(FFLAGS) -L$(LIBPATH1) -Wl,-rpath,$(LIBPATH1) $(LIBS) -L$(LIBPATH2) -Wl,-rpath,$(LIBPATH2) #-L$(LIBPATH3) #MD add LIBPATH2 and LIBPATH3 as rpath (for dynamic libraries)
LDDEBUG  = $(DBGFLAGS) -L$(LIBPATH1) -Wl,-rpath,$(LIBPATH1) $(LIBS) -L$(LIBPATH2) -Wl,-rpath,$(LIBPATH2) -Wl, -rpaht,$(LIBPATH3) #MD add LIBPATH2 and LIBPATH3 as rpath (for dynamic libraries)

MODOBJS = \
par_mod.o    		class_gribfile_mod.o \
particle_mod.o       interpol_mod.o \
coordinates_ecmwf.o  conv_mod.o \
windfields_mod.o     txt_output_mod.o \
turbulence_mod.o     cmapf_mod.o \
unc_mod.o            oh_mod.o \
xmass_mod.o          flux_mod.o \
point_mod.o          outg_mod.o \
mean_mod.o	 	random_mod.o \
com_mod.o		advance_mod.o \
readoptions.o  		binary_output_mod.o \
drydepo_mod.o 		wetdepo_mod.o

OBJECTS_NCF = netcdf_output_mod.o

OBJECTS = \
releaseparticles.o 	partoutput.o \
conccalc.o \
timemanager.o FLEXPART.o	\
readpartpositions.o \
domainfill.o \
initialize_particle.o		\
assignland.o \
caldate.o               \
coordtrafo.o            psih.o \
psim.o			\
pbl_profile.o \
juldate.o               \
calcpar.o               \
detectformat.o  \
qvsat.o                 \
plumetraj.o \
calcpv.o \
initial_cond_calc.o     \
get_settling.o	\
cbl.o 

ifeq ($(ncf), no)
	OBJECTS	:= $(OBJECTS)
else
	OBJECTS := $(OBJECTS) $(OBJECTS_NCF)
endif

%.o: %.mod

# serial executable
serial: $(FLEXPART-SERIAL)
serial: FC := $(F90)

$(FLEXPART-SERIAL): $(MODOBJS) $(OBJECTS)
	+$(FC) -o $@ $(MODOBJS) $(OBJECTS) $(LDFLAGS)

%.o: %.f90
	+$(FC) -c $(FFLAGS) $<

clean:
	\rm -f *.o *.mod

cleanall:
	\rm -f *.o *.mod $(FLEXPART-SERIAL)


.SUFFIXES = $(SUFFIXES) .f90

## DEPENDENCIES
advance_mod.o: cmapf_mod.o com_mod.o turbulence_mod.o interpol_mod.o par_mod.o \
	point_mod.o random_mod.o coordinates_ecmwf.o particle_mod.o
assignland.o: com_mod.o par_mod.o windfields_mod.o
binary_output_mod.o: com_mod.o mean_mod.o outg_mod.o par_mod.o point_mod.o unc_mod.o \
	windfields_mod.o
calcpar.o: com_mod.o par_mod.o drydepo_mod.o windfields_mod.o
calcpv.o: com_mod.o par_mod.o windfields_mod.o
caldate.o: par_mod.o
cbl.o: com_mod.o par_mod.o random_mod.o
cmapf_mod.o: par_mod.o
com_mod.o: par_mod.o
conccalc.o: com_mod.o outg_mod.o par_mod.o unc_mod.o coordinates_ecmwf.o particle_mod.o
conv_mod.o: par_mod.o com_mod.o flux_mod.o par_mod.o particle_mod.o class_gribfile_mod.o \
	random_mod.o interpol_mod.o coordinates_ecmwf.o windfields_mod.o
coordinates_ecmwf.o: par_mod.o com_mod.o interpol_mod.o particle_mod.o windfields_mod.o
coordtrafo.o: com_mod.o par_mod.o point_mod.o
detectformat.o: com_mod.o par_mod.o
domainfill.o: com_mod.o par_mod.o point_mod.o random_mod.o \
       interpol_mod.o coordinates_ecmwf.o particle_mod.o windfields_mod.o
drydepo_mod.o: com_mod.o par_mod.o unc_mod.o interpol_mod.o point_mod.o windfields_mod.o
FLEXPART.o: com_mod.o conv_mod.o netcdf_output_mod.o par_mod.o point_mod.o \
	random_mod.o readoptions.o binary_output_mod.o txt_output_mod.o windfields_mod.o
flux_mod.o: com_mod.o outg_mod.o par_mod.o coordinates_ecmwf.o particle_mod.o
get_settling.o: com_mod.o par_mod.o windfields_mod.o
initial_cond_calc.o: com_mod.o outg_mod.o par_mod.o unc_mod.o particle_mod.o \
	coordinates_ecmwf.o
initialize_particle.o: com_mod.o turbulence_mod.o interpol_mod.o par_mod.o random_mod.o \
	interpol_mod.o coordinates_ecmwf.o particle_mod.o windfields_mod.o
interpol_mod.o: par_mod.o com_mod.o turbulence_mod.o particle_mod.o windfields_mod.o
juldate.o: par_mod.o
mean_mod.o: par_mod.o
netcdf_output_mod.o: com_mod.o mean_mod.o outg_mod.o par_mod.o point_mod.o \
	unc_mod.o particle_mod.o windfields_mod.o
oh_mod.o: com_mod.o par_mod.o particle_mod.o windfields_mod.o
outg_mod.o: par_mod.o com_mod.o oh_mod.o par_mod.o unc_mod.o windfields_mod.o
particle_mod.o: com_mod.o par_mod.o
partoutput.o: com_mod.o par_mod.o interpol_mod.o coordinates_ecmwf.o netcdf_output_mod.o \
	particle_mod.o
pbl_profile.o: par_mod.o
plumetraj.o: com_mod.o mean_mod.o par_mod.o point_mod.o particle_mod.o coordinates_ecmwf.o \
	windfields_mod.o
psih.o: par_mod.o
psim.o: par_mod.o
readoptions.o: com_mod.o par_mod.o drydepo_mod.o
readpartpositions.o: com_mod.o par_mod.o random_mod.o particle_mod.o \
	particle_mod.o coordinates_ecmwf.o netcdf_output_mod.o
releaseparticles.o: com_mod.o par_mod.o point_mod.o random_mod.o xmass_mod.o \
	interpol_mod.o coordinates_ecmwf.o windfields_mod.o
timemanager.o: com_mod.o flux_mod.o netcdf_output_mod.o oh_mod.o outg_mod.o advance_mod.o\
	par_mod.o point_mod.o unc_mod.o xmass_mod.o particle_mod.o coordinates_ecmwf.o \
	binary_output_mod.o drydepo_mod.o wetdepo_mod.o conv_mod.o windfields_mod.o
turbulence_mod.o: par_mod.o com_mod.o random_mod.o particle_mod.o
txt_output_mod.o: com_mod.o par_mod.o point_mod.o outg_mod.o
unc_mod.o: par_mod.o
wetdepo_mod.o: com_mod.o par_mod.o point_mod.o particle_mod.o interpol_mod.o \
	coordinates_ecmwf.o unc_mod.o windfields_mod.o
windfields_mod.o: par_mod.o com_mod.o cmapf_mod.o point_mod.o class_gribfile_mod.o