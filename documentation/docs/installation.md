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
_FLEXPART_ 11 is written in Fortran 2018. The following compilers can be used to compile _FLEXPART_:

  - GNU Fortran compiler version 8+ (`gfortran`)
  - Intel Fortran compiler (`ifort`)

For running _FLEXPART_ in parallel mode, a compiler supporting [OpenMP](https://www.openmp.org/) is required.

## Libraries
_FLEXPART_ uses the following libraries:

  - [ecCodes](https://confluence.ecmwf.int/display/ECC)
  - [netCDF](https://docs.unidata.ucar.edu/netcdf-fortran/current/)

These libraries are usually available as packages in most Linux distributions and MacOS package managers. For example:

  - In Debian/Ubuntu: `sudo apt install libeccodes-dev libnetcdff` 
  - In MacOS + Homebrew: `brew install eccodes netcdf`

## Parameters
Before compiling FLEXPART, you might want to change parameters defined in par_mod.f90

- `wind_coord_type`: for ECMWF meteorological data, you can set this to ETA to use the native eta coordinate system. Otherwise, set this to METER.
- `mesoscale_turbulence`: by default the mesocale turbulence is switched off, but can be switched on setting this variable to .true.
- `max_partoutput_filesize`: maximum output of each partoutput NetCDF-4 file in Mb before a new one is created.
- `max_numthreads_grid`: when using many openmp threads and gridded output (IOUT>0 in COMMAND option file), this variable sets a maximum on how many threads are used for doing the reductions on the grid. A high number can result in a significant increase in RAM usage.

## <a name="compiling"></a>Compiling FLEXPART
_FLEXPART_ is compiled with [make](https://www.gnu.org/software/make/), which uses the makefile in the `src` subdirectory. Starting from the root directory, you can then compile _FLEXPART_ with the following steps:

    $ cd src/
    $ make -j -f <prefered_makefile>

This will create the executable `FLEXPART`

### Compiling FLEXPART with eta coordinates (ECMWF only)

    $ make -j -f <prefered_makefile> eta=yes

This will create the executable `FLEXPART_ETA`

### Compiling FLEXPART without NetCDF libraries

    $ make -j -f <prefered_makefile> eta=<yes/no> ncf=no

This will create the executable `FLEXPART_BIN`, or `FLEXPART_ETA_BIN` for `eta=yes`.
### <a name="paths"></a>Paths
the makefile provided uses a predefined `CPATH` and `LIBRARY_PATH`, the default on most servers. If these paths are not available, you will need to edit the makefile, adding the correct paths manually.

### <a name="optimisation"></a>Optimisation
The default and recommended optimisation flags in the makefile are `-O3 -march=native -mtune=native -fopenmp`. If not requiring a parallelised version, the `-fopenmp` can be switched off without consequence.

### <a name="dependency"></a>Dependencies
At the bottom of the default makefile, one can find the list of modules an how they depend on each other in ascending order.