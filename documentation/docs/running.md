# Running

<span style="color:red">
**WARNING:**
</span>
For many systems it is required to print the following in the command line or to set it in your submit script before executing FLEXPART compiled with OpenMP:
~~~
ulimit -s unlimited
~~~
Not setting this can result in a segmentation fault close to the start of your simulation.

Additionaly, one should set the following for optimal memory usage:
~~~
export OMP_PLACES=cores
export OMP_PROC_BIND=true
~~~
Not setting these can result in unnecessary slow performance of your application.

## Command line usage

FLEXPART accepts two command line options:

- `pathnames`, setting all appropriate paths, as explained in [Configuration](configuration.md#config).
- `-v <verbosity>`, currently not operational

## Exit code

The introduction of `error stop` in Fortran 2008 now garantees FLEXPART to only exit with code `0` for successful runs. Any other exit code indicates a failed run.

## Input data

To run FLEXPART, there are three important (sets) of files that need to be specified.
These are:

- the [**option files**](configuration.md#options), defining the set-up of the run,
- the [**pathnames file**](configuration.md#pathnames), defining the paths of where input and output are located, 
- the [**AVAILABLE file**](configuration.md#available), listing all available meteorological input

A full description can be found in [Configuration](configuration.md#config).

## OpenMP

Where most of FLEXPART's computational time is spent is very dependent on the specific problem to be solved and the set-up of FLEXPART. For example, when many particles are released from a single release point, initially most time is spent on particle trajectory computations. However, when a global high-resolution domain for the meteorological input data is used, significant time is spent on the convection computations on the grid. On the other hand, when few particles are used, computations on the gridded meteorological input data (e.g., coordinate transformations) are taking a large share. For this reason, we implemented OpenMP parallelisation throughout FLEXPART and tried to avoid bottlenecks at least for the most common set-ups.

We parallelised all particle based computations, apart from their initial release in the \texttt{releaseparticles} subroutine. On top of that, we parallelised the reading and computations on the meteorological fields, including the convection, wet and dry deposition, and the vertical coordinate transformation of the fields. Lastly, we parallelised the computations needed for the output, both for the gridded output and the particle dump.

One drawback of OpenMP parallelisation is that it is more difficult for users to make changes than in serial code, since they also are likely to have to update OpenMP regions. To minimise errors, we therefore strongly recommend users to make changes in the form of subroutines and functions and avoid the use of global variables.

## HPC systems

### SLURM example script

```
#!/bin/bash
#SBATCH --job-name=example
#SBATCH --output=example.log
#SBATCH --nodes=1 --ntasks-per-node=10 --ntasks-per-core=2 
#SBATCH --mem=30GB 
#SBATCH --time=20:00:00

export OMP_NUM_THREADS=10
export OMP_PLACES=cores
export OMP_PROC_BIND=true
ulimit -s unlimited

./FLEXPART_ETA pathnames
```