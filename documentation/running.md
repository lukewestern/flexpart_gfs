# Running
To run FLEXPART, there are three important (sets) of files that need to be specified.
These are:
- the **option files**, defining the set-up of the run,
- the **pathnames file**, defining the paths of where input and output are located, 
- and the **AVAILABLE file**, listing all available meteorological input

In addition to the regular input files listed above, a simulation can also be started using a NetCDF file listing all particles to be released. This option can be switched on by specifying IPIN=3 in the COMMAND option file. More information about how to use this option can be found in Silvia Bucci's paper (link)

## Option files
These files define the simulation settings. At the start of a simulation, a copy of each file will be written to the output directory defined in the **pathnames file**.
All option files should be presented as namelists (i.e. &OPTIONFILE).
### COMMAND
Sets the behaviour of the run (time range, backward or forward, output frequency, etc.). A table of all options is listed below
| Variable name | Description | Value **default** |
| ----------- | ----------- | ----------- |
| LDIRECT | Simulation direction in time | **1 (forward)** or -1 (backward) |
| IBDATE | Start date of the simulation | YYYYMMDD: YYYY=year, MM=month, DD=day |
| IBTIME | Start time of the simulation | HHMISS: HH=hours, MI=minutes, SS=seconds. UTC zone. |
| IEDATE | End date of the simulation | YYYYMMDD: YYYY=year, MM=month, DD=day |
| IETIME | End time of the simulation | HHMISS: HH=hours, MI=minutes, SS=seconds. UTC zone. |
| LOUTSTEP | Interval of model output. Average concentrations are calculated every LOUTSTEP (seconds) | **10800** |
| LOUTAVER | Concentration averaging interval, instantaneous for value of zero (seconds) | **10800** |
| LOUTSAMPLE | Numerical sampling rate of output, higher statistical accuracy with shorter intervals (seconds) | **900** |
| LOUTRESTART | Time interval when a restart file is written (seconds) | **999999999** |
| LSYNCTIME | All processes are synchronized to this time interval; all values above should be dividable by this number (seconds) | **900** |
| CTL | Factor by which particle transport time step in the ABL must be smaller than the Lagrangian timescale t l ; resulting time steps can be shorter than LSYNCTIME; LSYNCTIME is used if CTL < 0 | **-5.0** |
| IFINE | Additional reduction factor for time step used for vertical transport only considered if CTL > 1 | **4** |
| IOUT | Switch determining the gridded output type | 0 (no gridded output), **1 (forward: mass concentration; backwards: residence time)**, 2 (volume mixing ratio), 3 (1 and 2 combined), 4 (plume trajectories), 5 (1 and 4 combined), Add 8 for NetCDF output |
| IPOUT | Switch for particle position output | **0 (no particle output)**, 1 (particle output every LOUTSTEP), 2 (particle output at the end of the simulation) |
| LSUBGRID | Increase in ABL heights due to subgrid-scale orographic variations | **0 (off)**, 1 (on) |
| LCONVECTION | Switch for convection parameterization | 0 (off), **1 (on)** |
| LAGESPECTRA | Switch for calculation of age spectra (needs file AGECLASSES option file) | 0 (off), **1 (on)** |
| IPIN | Particle information input | **0 (using RELEASES option file)**, 1 (using restart.bin file), 2 (using previous partoutput file), 3 (self made initial conditions), 4 (restart.bin and self made initial conditions) |
| IOUTPUTFOREACHRELEASE | Switch for separate output fields for each location in the RELEASE file | 0 (no), **1 (yes)** |
| IFLUX | Output of mass fluxes through output grid box boundaries (northward, southward, eastward, westward, upward and downward) | 0 (off), **1 (on)** |
| MDOMAINFILL | Switch for domain-filling calculations: particles are initialized to reproduce air density or stratospheric ozone density; for limited-area simulations, particles are generated at the domain boundaries | **0 (no)**, 1 (like air density), 2 (stratospheric ozone tracer) |
| IND_SOURCE | Unit to be used at the source; see Seibert and Frank (2004); Eckhardt et al. (2017) | **1 (mass)**, 2 (mass mixing ratio) |
| IND_RECEPTOR | Unit to be used at the receptor; see Seibert and Frank (2004); Eckhardt et al. (2017) | 0 (no receptor), **1 (mass)**, 2 (mass mixing ratio), 3 (backward only: wet deposition),  4 (backward only: dry depostion) |
| MQUASILAG | Quasi-Lagrangian mode to track individual numbered particles | **0 (off)**, 1 (on) |
| NESTED_OUTPUT | Switch to produce output also for a nested domain | **0 (no)**, 1 (yes) |
| LINIT_COND | Switch to produce output sensitivity to initial conditions given in concentration or mixing ratio units (in backwards mode only) | **0 (no receptor)**, 1 (mass), 2 (mass mixing ratio) |
| SURF_ONLY | Output of SRR for fluxes only for the lowest model layer, most useful for backward runs when LINIT_COND set to 1 or 2 | **0 (no)**, 1 (yes) |
| CBLFLAG | Skewed rather than Gaussian turbulence in the convective ABL; when turned on, very short time steps should be used (see CTL and IFINE) | **0 (no)**, 1 (yes) |
### RELEASES
### SPECIES
### OUTGRID
### OUTGRID_NEST
### AGECLASSES
### RECEPTORS
### PARTOPTIONS

## Pathnames file

## AVAILABLE file

## Restarting a simulation
In case your simulation crashes or if you simply want to extend your simulation period, it is possible to run using the restart option (COMMAND option file: IPIN=1). You will need to decide if you will need this option before starting your initial simulation: LOUTRESTART in the COMMAND option file needs to be set to an appropriate time interval. For example, you can choose to set LOUTRESTART = 172800 s to get a new restart file ever 2 days. The restart files are written in binary and their name specifies the time within your simulation period they are written.

To run from one of these files, simply rename the desired restart_XXX.bin file to restart.bin, set IPIN=1 and you can restart your run from there.

WARNING: If you chose to use gridded data output (IOUT>0), then new data will be written to this file. If it is not desirable to overwrite a gridded data output file from a previous run, copy this file to another directory.
