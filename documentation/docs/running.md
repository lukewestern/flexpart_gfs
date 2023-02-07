# Running FLEXPART
To run FLEXPART, there are three important (sets) of files that need to be specified.
These are:

- the [**option files**](running.md#options), defining the set-up of the run,
- the [**pathnames file**](running.md#pathnames), defining the paths of where input and output are located, 
- the [**AVAILABLE file**](running.md#available), listing all available meteorological input,

Of course, there is also the **par_mod.f90** file, which needs to be specified before compiling (see section compiling), but the parameters in this file are expected to not have to be changed between simulations.

In addition to the regular input files listed above, a simulation can also be started using a NetCDF file listing all particles to be released. This option can be switched on by specifying IPIN=3 in the COMMAND option file. More information about how to use this option can be found in Silvia Bucci's paper (link).

When wanting to restart a previous simulation, see [restarting a simulation](running.md#restart).

## <a name="options"></a>Option files
These files define the simulation settings. At the start of a simulation, a copy of each file will be written to the output directory defined in the [**pathnames file**](running.md#pathnames).
All option files should be presented as namelists (i.e. &OPTIONFILE). A template of these files can be found in the options/ directory within the repository.

Inside the `options/` directory a template of all option files can be found:

- [COMMAND](running.md#command)
- [RELEASES](running.md#releases)
- [SPECIES](running.md#species)
- [OUTGRID](running.md#outgrid)
- [OUTGRID_NESTED](running.md#outgrid_nested)
- [AGECLASSES](running.md#ageclasses)
- [RECEPTORS](running.md#receptors)
- [PARTOPTIONS](running.md#partoptions)

### <a name="command"></a>COMMAND
Sets the behaviour of the run (time range, backward or forward, output frequency, etc.). A table of all options is listed below.

- **Time variables**: Flexpart can be run in forward or backward mode. In forward mode, particles are being traced forward in time, while in backward more, the origin of particles are being traced, going backward in time. This can be set by the [LDIRECT](running.md#ldirect) variable. The start and end of the simulation are set by [IBDATE](running.md#IBDATE):[IBTIME](running.md#IBTIME) and [IEDATE](running.md#IEDATE):[IETIME](running.md#IETIME). [IEDATE](running.md#IEDATE):[IETIME](running.md#IETIME) is always at a later time than [IBDATE](running.md#IBDATE):[IBTIME](running.md#IBTIME), also for backwards simulations. Output variables can be written at specified times: [LOUTSTEP](running.md#LOUTSTEP), and restart files will be written at every [LOUTRESTART](running.md#LOUTRESTART) interval.
- **Numerical variables**: [LSYNCTIME](running.md#LSYNCTIME) and LOUTSAMPLE set the integration interval, smaller generally giving better results, although below a certain number, not much will be gained. With the [CTL](running.md#CTL) and [IFINE](running.md#IFINE) setting, you can make integration steps even smaller for the turbulence computations.
- **Output variables**: The output is written at every [LOUTSTEP](running.md#LOUTSTEP) interval. Both gridded data ([IOUT](running.md#IOUT)>0) and particle based data ([IPOUT](running.md#IPOUT)=1) can be written to NetCDF files (binary option for gridded data). Nested output can be set by the [NESTED_OUTPUT](running.md#NESTED_OUTPUT) switched. Note that for gridded output, the [OUTGRID](running.md#OUTGRID) for ([IOUT](running.md#IOUT)>0) and [OUTGRID_NESTED](running.md#OUTGRID_NESTED) (for [NESTED_OUTPUT](running.md#NESTED_OUTPUT)=1) option files should be specified. Other output variables can be set in the par_mod.f90 file. Namely, the size of the NetCDF files that contain the particle based data (max_partoutput_filesize). [IND_RECEPTOR](running.md#IND_RECEPTOR) can be set to get concentrations or mixing ratios at specified receptor points set in the RECEPTORS options file. For backward simulations, [IND_RECEPTOR](running.md#IND_RECEPTOR) can be used to get wet or dry deposition gridded data. [SURF_ONLY](running.md#SURF_ONLY) and [LINIT_COND](running.md#LINIT_COND) are only working for binary output. 
- **Input variables**: IPIN can be set to chose the input type: either initial conditions from particles come from the [RELEASES](running.md#releases) file ([IPIN](running.md#IPIN)=0), from restart files of a previous run ([IPIN](running.md#IPIN)=1),
from a particle netCDF file written in a previous run (only works when the correct fields in [PARTOPTIONS](running.md#PARTOPTIONS) are chosen) ([IPIN](running.md#IPIN)=2), or from user-defined initial particle conditions ([IPIN](running.md#IPIN)=3). [MDOMAINFILL](running.md#MDOMAINFILL) can be set to distribute particles according to the air density or stratospheric ozone density profiles. This option overwrites the vertical levels set in the [RELEASES](running.md#releases) option file.

| Variable name | Description | Value **default** |
| ----------- | ----------- | ----------- |
| <a name="ldirect"></a>LDIRECT | Simulation direction in time | **1 (forward)** or -1 (backward) |
| <a name="IBDATE"></a>IBDATE | Start date of the simulation | YYYYMMDD: YYYY=year, MM=month, DD=day |
| <a name="IBTIME"></a>IBTIME | Start time of the simulation | HHMISS: HH=hours, MI=minutes, SS=seconds. UTC zone. |
| <a name="IEDATE"></a>IEDATE | End date of the simulation | YYYYMMDD: YYYY=year, MM=month, DD=day |
| <a name="IETIME"></a>IETIME | End time of the simulation | HHMISS: HH=hours, MI=minutes, SS=seconds. UTC zone. |
| <a name="LOUTSTEP"></a>LOUTSTEP | Interval of model output. Average concentrations are calculated every LOUTSTEP (seconds) | **10800** |
| <a name="LOUTAVER"></a>LOUTAVER | Concentration averaging interval, instantaneous for value of zero (seconds) | **10800** |
| <a name="LOUTSAMPLE"></a>LOUTSAMPLE | Numerical sampling rate of output, higher statistical accuracy with shorter intervals (seconds) | **900** |
| <a name="LOUTRESTART"></a>LOUTRESTART | Time interval when a restart file is written (seconds) | **999999999** |
| <a name="LSYNCTIME"></a>LSYNCTIME | All processes are synchronized to this time interval; all values above should be dividable by this number (seconds) | **900** |
| <a name="CTL"></a>CTL | Factor by which particle transport time step in the ABL must be smaller than the Lagrangian timescale t l ; resulting time steps can be shorter than LSYNCTIME; LSYNCTIME is used if CTL < 0 | **-5.0** |
| <a name="IFINE"></a>IFINE | Additional reduction factor for time step used for vertical transport only considered if CTL > 1 | **4** |
| <a name="IOUT"></a>IOUT | Switch determining the gridded output type | 0 (no gridded output), **1 (forward: mass concentration; backwards: residence time)**, 2 (volume mixing ratio), 3 (1 and 2 combined), 4 (plume trajectories), 5 (1 and 4 combined), Add 8 for NetCDF output |
| <a name="IPOUT"></a>IPOUT | Switch for particle position output | **0 (no particle output)**, 1 (particle output every LOUTSTEP), 2 (particle output at the end of the simulation) |
| <a name="LSUBGRID"></a>LSUBGRID | Increase in ABL heights due to subgrid-scale orographic variations | **0 (off)**, 1 (on) |
| <a name="LCONVECTION"></a>LCONVECTION | Switch for convection parameterization | 0 (off), **1 (on)** |
| <a name="LAGESPECTRA"></a>LAGESPECTRA | Switch for calculation of age spectra (needs file [AGECLASSES](running.md#ageclasses) option file) | 0 (off), **1 (on)** |
| <a name="IPIN"></a>IPIN | Particle information input. Starting from [RELEASES](running.md#releases) option file, form restart.bin, or user-defined particle input data (see Silvia Bucci's stuff) | **0 (using RELEASES option file)**, 1 (using restart.bin file), 2 (using previous partoutput file), 3 (self made initial conditions), 4 (restart.bin and self made initial conditions) |
| <a name="IOUTPUTFOREACHRELEASE"></a>IOUTPUTFOREACHRELEASE | Switch for separate output fields for each location in the [RELEASES](running.md#releases) file | 0 (no), **1 (yes)** |
| <a name="IFLUX"></a>IFLUX | Output of mass fluxes through output grid box boundaries (northward, southward, eastward, westward, upward and downward) | 0 (off), **1 (on)** |
| <a name="MDOMAINFILL"></a>MDOMAINFILL | Switch for domain-filling calculations: particles are initialized to reproduce air density or stratospheric ozone density; for limited-area simulations, particles are generated at the domain boundaries | **0 (no)**, 1 (like air density), 2 (stratospheric ozone tracer) |
| <a name="IND_SOURCE"></a>IND_SOURCE | Unit to be used at the source; see Seibert and Frank (2004); Eckhardt et al. (2017) | **1 (mass)**, 2 (mass mixing ratio) |
| <a name="IND_RECEPTOR"></a>IND_RECEPTOR | Unit to be used at the receptor; see Seibert and Frank (2004); Eckhardt et al. (2017) | 0 (no receptor), **1 (mass)**, 2 (mass mixing ratio), 3 (backward only: wet deposition),  4 (backward only: dry depostion) |
| <a name="MQUASILAG"></a>MQUASILAG | Quasi-Lagrangian mode to track individual numbered particles | **0 (off)**, 1 (on) |
| <a name="NESTED_OUTPUT"></a>NESTED_OUTPUT | Switch to produce output also for a nested domain | **0 (no)**, 1 (yes) |
| <a name="LINIT_COND"></a>LINIT_COND | Switch to produce output sensitivity to initial conditions given in concentration or mixing ratio units (in backwards mode only) | **0 (no)**, 1 (mass), 2 (mass mixing ratio) |
| <a name="SURF_ONLY"></a>SURF_ONLY | Output of SRR for fluxes only for the lowest model layer, most useful for backward runs when LINIT_COND set to 1 or 2 | **0 (no)**, 1 (yes) |
| <a name="CBLFLAG"></a>CBLFLAG | Skewed rather than Gaussian turbulence in the convective ABL; when turned on, very short time steps should be used (see CTL and IFINE) | **0 (no)**, 1 (yes) |

### <a name="releases"></a>RELEASES
This file contains the information about the particles initial conditions: how many, where and when they will be released, their mass and what species they are (defined in the SPECIES files).

### <a name="species"></a>SPECIES
The subdirectory options/SPECIES/ needs to contain one or more files named SPECIES_nnn. For each species nnn listed in the header section of the RELEASES file, such a SPECIES_nnn file must exist. The parameters in the SPECIES_nnn file, contained in the namelist &SPECIES_PARAMS, set the species name and define the physicochemical properties of the species; they are described in Table 10. These are important for simulating radioactive or chemical decay, wet deposition (scavenging) for gases and aerosols, dry deposition for gases and aerosols, particle settling, and chemical reaction with the OH radical. Some parameters are only necessary for gas tracers and some are only necessary for aerosol tracers; thus, a namelist does not need to contain all parameters for both gases and particles. Optionally, since FLEXPART version 6.0, information about temporal emission variations can be added at the end of the file.

The following specifies the parameters associated with each physicochemical process simulated.

- Radioactive or chemical decay: set with pdecay; off if pdecay<0.
- Wet deposition for gases: set with pweta_gas, pwetb_gas (for below-cloud) and phenry (for in-cloud). Switch off for both in- and below-cloud if either pweta_gas or pwetb_gas is negative.
- Wet deposition for aerosols: set with pccn_aero, pin_aero for in-cloud scavenging and pcrain_aero, pcsnow_aero and pdquer for below-cloud scavenging.
- Dry deposition for aerosols: set with pdensity, pdquer, pndia, and psigma; off if pdensity < 0. 
- Dry deposition for gases: set with phenry, pf0 and preldiff; off if preldiff < 0. Alternatively, a constant dry deposition velocity pdryvel can be given. 
- Settling of particles: set with pdensity and pdquer.
- Shape of particles: set with PSHAPE, PASPECTRATIO, PLA, PIA, PSA, and PORIENT
- OH reaction: chemical reaction with the OH radical can be turned on by giving parameter pohcconst (cm 3 molecule −1 s −1 ), pohdconst (K) and pohnconst (no unit) positive values; defined by Eq. (13) in Pisso et al. (2019).
- Emission variation: emission variation during the hours (local time) of the day and during the days of the week can be specified. Factors should be 1.0 on average to obtain unbiased emissions overall. The area source factors (useful, e.g., for traffic emissions) are applied to emis sions with a lower release height below 0.5 m above ground level (a.g.l.) and the point source factors (useful, e.g., for power plant emissions) to emissions with a lower release height than 0.5 m a.g.l. Default values are 1.0.

| Variable name | Description |
| ----------- | ----------- |
|PSPECIES | Tracer name |
|PDECAY | Species half life |
|PWETA_GAS | Below-cloud scavenging (gases) - A (weta_gas) |
|PWETB_GAS | Below-cloud scavenging (gases) - B (wetb_gas) |
|PCRAIN_AERO | Below-cloud scavenging (particles) - Crain (crain_aero) |
|PCSNOW_AERO | Below-cloud scavenging (particles) - Csnow (csnow_aero) |
|PCCN_AERO | In-cloud scavenging (particles) - CCNeff (ccn_aero) |
|PIN_AERO | In-cloud scavenging (particles) - INeff (in_aero) |
|PDENSITY | Dry deposition (particles) - rho |
|PDQUER | Dry deposition (particles) - dquer (equivalent diameter for shape) |
|PDSIGMA | Dry deposition (particles) - dsig |
|PNDIA | Dry deposition (particles) - ndia |
|PDRYVEL | Alternative: dry deposition velocity |
|PRELDIFF | Dry deposition (gases) - D |
|PHENRY | Dry deposition (gases) - Henrys const. |
|PF0 | Dry deposition (gases) - f0 (reactivity) |
|PWEIGHTMOLAR | molweight |
|POHCCONST | OH Reaction rate - C [cm^3/molecule/sec] |
|POHDCONST | OH Reaction rate - D [K] |
|POHNCONST | OH Reaction rate - C [cm^3/molecule/sec] |
|PSHAPE | 0 for sphere, 1 any shape (defined by axes PLA,PIA,PSA), 2-cylinder, 3-cube, 4-tetrahedron, 5-octahedron, 6-ellipsoid |
|PASPECTRATIO | Aspect ratio of cylinders: works for PSHAPE=2 only |
|PLA | Longest axis in micrometer (Bagheri & Bonadonna 2016): only for PSHAPE=1 |
|PIA | Intermediate axis in micrometer: only for PSHAPE=1 |
|PSA | Smallest axis in micrometer: only for PSHAPE=1 |
|PORIENT | 0 for horizontal, 1 for random orientation of particles, 2 for an average between random and horizontal |

### <a name="outgrid"></a>OUTGRID

### <a name="outgrid_nest"></a>OUTGRID_NEST

### <a name="ageclasses"></a>AGECLASSES

### <a name="receptors"></a>RECEPTORS

### <a name="partoptions"></a>PARTOPTIONS
This option file is only necessary when requiring particle properties to be written out (IPOUT=1 in the COMMAND option file). In this file, the user can set what particle properties and interpolated fields they want to be written to files. At the moment, the available fields that can be written to file are:

- particle positions (longitude, latitude and height), 
- potential vorticity, 
- specific humidity, 
- density, temperature, 
- pressure, 
- particle mass, 
- separate cumulative wet and dry deposition masses, 
- settling velocity, 
- 3D velocities, 
- the height of the PBL, tropopause and topography.

Each property can also be printed out as an average instead of an instantaneous value. For example, if one makes internal time steps of 600 seconds each,
and writes properties to files every hour, the outputted value will be the average of the 6 previous values of the particle of the past hour. Note that this comes with an additional computational cost.

## <a name="pathnames"></a>Pathnames file

## <a name="available"></a>AVAILABLE file

## <a name="restart"></a>Restarting a simulation
In case your simulation crashes or if you simply want to extend your simulation period, it is possible to run using the restart option (COMMAND option file: IPIN=1). You will need to decide if you will need this option before starting your initial simulation: LOUTRESTART in the COMMAND option file needs to be set to an appropriate time interval. For example, you can choose to set LOUTRESTART = 172800 s to get a new restart file ever 2 days. The restart files are written in binary and their name specifies the time within your simulation period they are written.

To run from one of these files, simply rename the desired restart_XXX.bin file to restart.bin, set IPIN=1 and you can restart your run from there.

WARNING: If you chose to use gridded data output (IOUT>0), then new data will be written to this file. If it is not desirable to overwrite a gridded data output file from a previous run, copy this file to another directory.
