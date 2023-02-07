# Installation

## Download FLEXPART
There are two options to download _FLEXPART_:

  - **tar ball**

    You can download a tar ball with the latest release from the [_FLEXPART_ page](https://www.flexpart.eu/wiki/FpRoadmap) and then untar the file:

        $ tar -xvf <flexpart_vX.X.tar>

  - **git repository**

    If you have git installed, you can clone the latest version from our git repository:

        $ git clone --single-branch --branch master https://gitlab.phaidra.org/flexpart/flexpart

## Compiler
_FLEXPART_ is written in Fortran 95(or are there newer elements?). The following compilers can be used to compile _FLEXPART_:

  - GNU Fortran compiler (`gfortran`)
  - Intel Fortran compiler (`ifort`)
  - others?

For running _FLEXPART_ in parallel mode, a compiler supporting [OpenMP](https://www.openmp.org/) is required.

## Libraries
_FLEXPART_ uses the following libraries:

  - [ecCodes](https://confluence.ecmwf.int/display/ECC)
  - [netCDF](https://docs.unidata.ucar.edu/netcdf-fortran/current/)

These libraries are usually available as packages in most Linux distributions and MacOS package managers. For example:

  - In Debian/Ubuntu: `sudo apt install libeccodes-dev libnetcdff` 
  - In Fedora: `no idea`
  - In MacOS + Homebrew: `brew install eccodes netcdf`

## Parameters
Before compiling FLEXPART, you might want to change parameters defined in par_mod.f90

- *wind_coord_type*: for ECMWF meteorological data, you can set this to ETA to use the native eta coordinate system. Otherwise, set this to METER.
- *mesoscale_turbulence*: by default the mesocale turbulence is switched off, but can be switched on setting this variable to .true.
- *max_partoutput_filesize*: maximum output of each partoutput NetCDF-4 file in Mb before a new one is created.
- *max_numthreads_grid*: when using many openmp threads and gridded output (IOUT>0 in COMMAND option file), this variable sets a maximum on how many threads are used for doing the reductions on the grid. A high number can result in a significant increase in RAM usage.

## Compiling
_FLEXPART_ is compiled with [make](https://www.gnu.org/software/make/), which uses the makefile in the `src` subdirectory. In the makefile, make sure that the library and include paths point to the directories where `ecCodes` and `netCDF` are installed. Starting from the root directory, you can then compile _FLEXPART_ with the following steps:

    $ cd src/
    $ make -j -f <prefered_makefile>

