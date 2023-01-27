# Running
To run FLEXPART, there are three important (sets) of files that need to be specified.
These are:
- the **option files**, defining the set-up of the run,
- the **pathnames file**, defining the paths of where input and output are located, 
- and the **AVAILABLE file**, listing all available meteorological input

## Option files
These files define the simulation settings. At the start of a simulation, a copy of each file will be written to the output directory defined in the **pathnames file**.
All option files should be presented as namelists (i.e. &OPTIONFILE).
### COMMAND
Sets the behaviour of the run (time range, backward or forward, output frequency, etc.).
| Variable name | Description | Value **default** |
| ----------- | ----------- | ----------- |
| LDIRECT | Simulation direction in time | **1 (forward)** or -1 (backward) |
| IBDATE | Start date of the simulation | YYYYMMDD: YYYY=year, MM=month, DD=day |
| IBTIME | Start time of the simulation | HHMISS: HH=hours, MI=minutes, SS=seconds. UTC zone. |
| IEDATE | End date of the simulation | YYYYMMDD: YYYY=year, MM=month, DD=day |
| IETIME | End time of the simulation | HHMISS: HH=hours, MI=minutes, SS=seconds. UTC zone. |
| LOUTSTEP | Interval of model output. Average concentrations are calculated every LOUTSTEP | **10800** s |
| LOUTAVER | Concentration averaging interval, instantaneous for value of zero | **10800** s |
| LOUTSAMPLE | Numerical sampling rate of output, higher statistical accuracy with shorter intervals | **900** s |
| LOUTRESTART | Time interval when a restart file is written | **999999999** |
| LSYNCTIME | All processes are synchronized to this time interval; all values above should be dividable by this number | **900** s |
| CTL |     | **-5.0** |
| IFINE |     | 4 |
| IOUT |     | 1 |
| IPOUT |     | 0 |
| LSUBGRID |     | 0 |
| LCONVECTION |     | 1 |
| LAGESPECTRA |     | 0 |
| IPIN |     | 0 |
| IOUTPUTFOREACHRELEASE |     | 1 |
| IFLUX |     | 0 |
| MDOMAINFILL |     | 0 |
| IND_SOURCE |     | 1 |
| IND_RECEPTOR |     | 1 |
| MQUASILAG |     | 0 |
| NESTED_OUTPUT |     | 0 |
| LINIT_COND |     | 0 |
| SURF_ONLY |     | 0 |
| CBLFLAG |     | 0 |
### RELEASES
### SPECIES
### OUTGRID
### OUTGRID_NEST
### AGECLASSES
### RECEPTORS
### PARTOPTIONS

## Pathnames file

## AVAILABLE file
