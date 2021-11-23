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
par_mod.o    		com_mod.o \
particle_mod.o       interpol_mod.o \
coordinates_ecmwf.o  conv_mod.o \
windfields_mod.o \
turbulence_mod.o     cmapf_mod.o \
unc_mod.o            oh_mod.o \
xmass_mod.o          flux_mod.o \
point_mod.o          outg_mod.o \
mean_mod.o	 	random_mod.o \
class_gribfile_mod.o	advance_mod.o \
readoptions.o  		binary_output_mod.o	


## Serial versions (MPI version with same functionality and name '_mpi.f90' exists)
OBJECTS_SERIAL = \
	releaseparticles.o 	partoutput.o \
	conccalc.o \
	init_domainfill.o \
	timemanager.o FLEXPART.o	\
	readpartpositions.o \
	boundcond_domainfill.o	\
	redist.o		

OBJECTS_NCF = netcdf_output_mod.o

OBJECTS = \
initialize.o		\
writeprecip.o \
assignland.o\
part0.o 		gethourlyOH.o\
caldate.o               partdep.o \
coordtrafo.o            psih.o \
raerod.o 		drydepokernel.o \
erf.o                    \
ew.o			 \
get_vdep_prob.o   \
get_wetscav.o \
psim.o			outgrid_init.o 	\
outgrid_init_nest.o   	\
photo_O1D.o 		  \
getrb.o                 obukhov.o \
getrc.o                 \
getvdep.o \
richardson.o \
scalev.o                \
pbl_profile.o \
juldate.o               \
wetdepokernel.o \
calcpar.o               wetdepo.o \
windalign.o		\
calcpar_nests.o \
detectformat.o  \
calcfluxes.o            fluxoutput.o \
qvsat.o                 skplin.o \
sort2.o                 distance.o \
centerofmass.o          plumetraj.o \
openouttraj.o           calcpv.o \
distance2.o \
clustering.o 		openreceptors.o \
wetdepokernel_nest.o \
drydepokernel_nest.o    zenithangle.o \
ohreaction.o            getvdep_nests.o \
initial_cond_calc.o     initial_cond_output.o initial_cond_output_inversion.o \
dynamic_viscosity.o     get_settling.o	\
initialize_cbl_vel.o	re_initialize_particle.o \
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

$(FLEXPART-SERIAL): $(MODOBJS) $(OBJECTS) $(OBJECTS_SERIAL)
	+$(FC) -o $@ $(MODOBJS) $(OBJECTS) $(OBJECTS_SERIAL) $(LDFLAGS)

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
assignland.o: com_mod.o par_mod.o
binary_output_mod.o: com_mod.o mean_mod.o outg_mod.o par_mod.o point_mod.o unc_mod.o
boundcond_domainfill.o: com_mod.o par_mod.o point_mod.o random_mod.o \
	coordinates_ecmwf.o interpol_mod.o particle_mod.o
calcfluxes.o: com_mod.o flux_mod.o outg_mod.o par_mod.o coordinates_ecmwf.o particle_mod.o
calcpar.o: com_mod.o par_mod.o
calcpar_nests.o: com_mod.o par_mod.o
calcpv.o: com_mod.o par_mod.o
caldate.o: par_mod.o
cbl.o: com_mod.o par_mod.o
centerofmass.o: par_mod.o
clustering.o: par_mod.o particle_mod.o
cmapf_mod.o: par_mod.o
com_mod.o: par_mod.o
conccalc.o: com_mod.o outg_mod.o par_mod.o unc_mod.o coordinates_ecmwf.o particle_mod.o
conv_mod.o: par_mod.o com_mod.o flux_mod.o par_mod.o particle_mod.o class_gribfile_mod.o
coordinates_ecmwf.o: par_mod.o com_mod.o interpol_mod.o particle_mod.o
coordtrafo.o: com_mod.o par_mod.o point_mod.o
detectformat.o: com_mod.o par_mod.o
distance.o: par_mod.o
distance2.o: par_mod.o
drydepokernel.o: com_mod.o par_mod.o unc_mod.o
drydepokernel_nest.o: com_mod.o par_mod.o unc_mod.o
erf.o: par_mod.o
FLEXPART.o: com_mod.o conv_mod.o netcdf_output_mod.o par_mod.o point_mod.o \
	random_mod.o readoptions.o binary_output_mod.o
fluxoutput.o: com_mod.o flux_mod.o outg_mod.o par_mod.o
get_settling.o: com_mod.o par_mod.o
get_vdep_prob.o: com_mod.o interpol_mod.o par_mod.o point_mod.o
get_wetscav.o: com_mod.o par_mod.o point_mod.o interpol_mod.o particle_mod.o \
	coordinates_ecmwf.o
gethourlyOH.o: com_mod.o oh_mod.o par_mod.o
getrb.o: par_mod.o
getrc.o: com_mod.o par_mod.o
getvdep.o: com_mod.o par_mod.o
getvdep_nests.o: com_mod.o par_mod.o
grib2check.o: com_mod.o par_mod.o
init_domainfill.o: com_mod.o par_mod.o point_mod.o random_mod.o \
       interpol_mod.o coordinates_ecmwf.o particle_mod.o
initial_cond_calc.o: com_mod.o outg_mod.o par_mod.o unc_mod.o particle_mod.o \
	coordinates_ecmwf.o
initial_cond_output.o: com_mod.o outg_mod.o par_mod.o point_mod.o unc_mod.o
initial_cond_output_inversion.o: com_mod.o outg_mod.o par_mod.o point_mod.o \
	unc_mod.o
initialize.o: com_mod.o turbulence_mod.o interpol_mod.o par_mod.o random_mod.o \
	interpol_mod.o coordinates_ecmwf.o particle_mod.o
initialize_cbl_vel.o: com_mod.o par_mod.o random_mod.o
interpol_mod.o: par_mod.o com_mod.o turbulence_mod.o particle_mod.o
juldate.o: par_mod.o
mean_mod.o: par_mod.o
netcdf_output_mod.o: com_mod.o mean_mod.o outg_mod.o par_mod.o point_mod.o \
	unc_mod.o particle_mod.o
obukhov.o: par_mod.o
ohreaction.o: com_mod.o oh_mod.o par_mod.o particle_mod.o 
openouttraj.o: com_mod.o par_mod.o point_mod.o
openreceptors.o: com_mod.o par_mod.o
outg_mod.o: par_mod.o
outgrid_init.o: com_mod.o flux_mod.o oh_mod.o outg_mod.o par_mod.o unc_mod.o
outgrid_init_nest.o: com_mod.o outg_mod.o par_mod.o unc_mod.o
part0.o: par_mod.o
partdep.o: com_mod.o par_mod.o
particle_mod.o: com_mod.o par_mod.o
partoutput.o: com_mod.o par_mod.o interpol_mod.o coordinates_ecmwf.o netcdf_output_mod.o \
	particle_mod.o
pbl_profile.o: par_mod.o
plumetraj.o: com_mod.o mean_mod.o par_mod.o point_mod.o particle_mod.o coordinates_ecmwf.o
psih.o: par_mod.o
psim.o: par_mod.o
raerod.o: par_mod.o
re_initialize_particle.o: com_mod.o par_mod.o
readoptions.o: com_mod.o par_mod.o
readpartpositions.o: com_mod.o par_mod.o random_mod.o particle_mod.o \
	particle_mod.o coordinates_ecmwf.o
redist.o: com_mod.o conv_mod.o par_mod.o random_mod.o interpol_mod.o coordinates_ecmwf.o \
	particle_mod.o
releaseparticles.o: com_mod.o par_mod.o point_mod.o random_mod.o xmass_mod.o \
	interpol_mod.o coordinates_ecmwf.o
richardson.o: par_mod.o
scalev.o: par_mod.o
shift_field.o: par_mod.o
shift_field_0.o: par_mod.o
timemanager.o: com_mod.o flux_mod.o netcdf_output_mod.o oh_mod.o outg_mod.o advance_mod.o\
	par_mod.o point_mod.o unc_mod.o xmass_mod.o particle_mod.o coordinates_ecmwf.o \
	binary_output_mod.o
turbulence_mod.o: par_mod.o com_mod.o random_mod.o particle_mod.o
unc_mod.o: par_mod.o
wetdepo.o: com_mod.o par_mod.o point_mod.o particle_mod.o
wetdepokernel.o: com_mod.o par_mod.o unc_mod.o
wetdepokernel_nest.o: com_mod.o par_mod.o unc_mod.o
windfields_mod.o: par_mod.o com_mod.o conv_mod.o cmapf_mod.o
writeprecip.o: com_mod.o par_mod.o point_mod.o
zenithangle.o: par_mod.o