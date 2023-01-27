# Running
To run FLEXPART, there are three important (sets) of files that need to be specified.
These are:
- the **option files**, defining the set-up of the run,
- the **pathnames file**, defining the paths of where input and output are located, 
- and the **AVAILABLE file**, listing all available meteorological input

## Option files
These files define the simulation settings. At the start of a simulation, a copy of each file will be written to the output directory defined in the **pathnames file**.
All option files should be presented as namelists (i.e. &OPTIONFILE)
### COMMAND
Sets the behaviour of the run (time range, backward or forward, output frequency, etc.)
### RELEASES
### SPECIES
### OUTGRID
### OUTGRID_NEST
### AGECLASSES
### RECEPTORS
### PARTOPTIONS

## Pathnames file

## AVAILABLE file
