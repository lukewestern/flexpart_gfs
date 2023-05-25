![](https://www.flexpart.eu/chrome/site/flexpart_banner.png)
# Welcome to Flexpart - The Lagrangian particle dispersion model

This is the main development site @ University of Vienna.

Other references:

- [FLEXPART.eu](https://flexpart.eu)
- [FLEXPART@NILU](https://git.nilu.no/flexpart/flexpart)

## What is this repository for?

* This repository contains versions of the Lagrangian model FLEXPART
* Development versions
* Issues on the FLEXPART model
* Feature requests for future versions

## Getting started with Flexpart

The model is written in Fortran. It needs to be compiled for the architecture that runs it.

### 1.Configuration

  The Makefiles e.g. `src/makefile_gfortran` can use environmental variables:
   - `CPATH` for include directories
   - `LIBRARY_PATH` for libraries

   This is commonly used on HPC systems. e.g. JET and VSC. These paths will be added via `rpath` so statically linked.
   Otherwise edit the makefile with the paths to libraries and include files.

**Required Dependencies:**

 * [ecCodes](https://confluence.ecmwf.int/display/ECC) from ECMWF. 
 * [NetCDF](https://www.unidata.ucar.edu/software/netcdf/) (optional) from UCAR
 * Fortran Compiler with [Fortran 2008](https://fortranwiki.org/fortran/show/Fortran+2008+status) standard, e.g. GCC compiler `8+` or INTEL `19+` or INTEL-ONEAPI 
 * make utils

### 2.Compilation

Clone the git repository or download one of the [releases](https://gitlab.phaidra.org/flexpart/flexpart/-/releases)

```bash
# clone the repository to your directory or download one of our releases
git clone https://gitlab.phaidra.org/flexpart/flexpart.git
# change to the SRC directory
cd flexpart/src
# Remember to configure your libraries in the makefile or environmental variables
# use for example the GCC makefile
make -f makefile_gfortran
# this will create the FLEXPART executable
file ./FLEXPART
# Check its dependencies:
ldd ./FLEXPART
```
Now you are almost ready to run.

### 3.Deployment instructions 

   FLEXPART is a standalone executable  
   The necessary ECMWF wind fields can be obtained testing flex_ecmwf
   The AVAILABLE file works with the default ERA 5 retrieved winds
   In the winds are available in flex_ecmwf/work it should suffice to execute 
   `./src/FLEXPART` in the main directory  

## Contribution guidelines

* The version contributed should compile on a reference version of the system and compiler. 
   - `FLEXPART 10.4` used as reference gfortran 5.4 on Ubuntu 16.04
   - `FLEXPART 11` uses as reference gfortran 8.5.0 on AlmaLinux 8
* Code contribution including new features and bug fixes should be complemented with appropriate tests
   An essential test consists of a set of input files and directories that allow FLEXPART to run.
   A test can be accompanied by output files for verification
* Code review
* report issues via mail to [support](mailto:flexpart-support.img-wien@univie.ac.at)
* become an active developer and request a user account.
