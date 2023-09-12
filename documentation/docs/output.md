# Output files

An overview of all possible output files is provided in the table below. Notice that not all these files are written out in every model run; the user settings control which files are produced. At the beginning of a run, FLEXPART records descriptive metadata to either the gridded NetCDF output file ([IOUT](running.md#IOUT)>8 or [LNETCDFOUT](running.md#LNETCDFOUT)=1), the particle output (IPOUT](running.md#iout)=1), or to a dedicated binary file header and a plain text file **header_txt** (with the exception of the orography data and release information). Corresponding files header_nest are produced if nested output is selected. 

When [IOUT](running.md#iout) is set to a non-zero value, FLEXPART produces gridded output. When requiring binary files ([IOUT](running.md#iout)<8), separate files are created for every output time step ([LOUTSTEP](running.md#loutstep)), species and domain (mother and, if requested, nest). The naming convention for these files is grid_type_date_nnn. When requiring NetCDF output ([IOUT](running.md#iout)>8), this information is all contained in one file for the mother grid, and one for a possible nested grid ([NESTED_OUTPUT](running.md#nested_output)). When [IPOUT](running.md#ipout) is switched on, this information is also contained in the **partoutput_nnn.nc** files.

| Name | Format | Switches | Description of contents |
| ---- | ------ | -------- | ----------------------- |
| ** *Header* ** | | | |
| **`header`** | binary | [IOUT](running.md#iout)<8 | run metadata + ancillary data, included in **grid_conc_date.nc** and **partoutput_date.nc** when [IOUT](running.md#iout)>8 and [IPOUT](running.md#ipout)=1, respectively. |
| **`header_nest`** | binary | [NESTED_OUTPUT](running.md#nested_output)=1 [IOUT](running.md#iout)<8 | run metadata + ancillary data, included in **grid_conc_date_nest.nc** and **partoutput_date_nest.nc** when [IOUT](running.md#iout)>8 and [IPOUT](running.md#ipout)=1, respectively. |
| **`header_txt`** | text | [IOUT](running.md#iout)<8 | human-readable run metadata (from [COMMAND](running.md#command)) |
| **`header_txt_releases`** | text | [IOUT](running.md#iout)<8 | human-readable run metadata (from [RELEASES](running.md#releases) or **part_ic.nc**) |
| **`dates`** | text | [IOUT](running.md#iout)<8 | time series: dates of output files |
| ** *Gridded data* ** | | | |
| **`grid_conc_date_nnn`** | binary | [LDIRECT](running.md#ldirect)=1 [IOUT](running.md#iout)=1,3,5 | 3D tracer mass density and 2D deposition |
| **`grid_pptv_date_nnn`** | binary | [LDIRECT](running.md#ldirect)=1 [IOUT](running.md#iout)=2,3 | 3D tracer volume mixing ratio and 2D deposition |
| **`grid_time_date_nnn`** | binary | [LDIRECT](running.md#ldirect)=-1 [IOUT](running.md#iout)=1 | 3D sensitivity of atmospheric receptor to emissions |
| **`grid_drydep_date_nnn`** | binary | [LDIRECT](running.md#ldirect)=-1 [IOUT](running.md#iout)=1 [IND_RECEPTOR](running.md#ind_receptor)=3 | 3D sensitivity of dry deposition receptor to emissions |
| **`grid_wetdep_date_nnn`** | binary | [LDIRECT](running.md#ldirect)=-1 [IOUT](running.md#iout)=1 [IND_RECEPTOR](running.md#ind_receptor)=4 | 3D sensitivity of wet deposition receptor to emissions |
| **`grid_conc_date.nc`** | NetCDF | [LDIRECT](running.md#ldirect)=1 [IOUT](running.md#iout)=9,10,11,13 | 3D tracer and 2D wet and dry deposition |
| **`grid_time_date.nc`** | NetCDF | [LDIRECT](running.md#ldirect)=-1 [IOUT](running.md#iout)=9 | 3D sensitivity of atmospheric receptor to emissions |
| **`grid_drydep_date.nc`** | NetCDF | [LDIRECT](running.md#ldirect)=-1 [IOUT](running.md#iout)=9 [IND_RECEPTOR](running.md#ind_receptor)=3 | 3D sensitivity of dry deposition receptor to emissions |
| **`grid_wetdep_date.nc`** | NetCDF | [LDIRECT](running.md#ldirect)=-1 [IOUT](running.md#iout)=9 [IND_RECEPTOR](running.md#ind_receptor)=4 | 3D sensitivity of wet deposition receptor to emissions |
| **`grid_initial_nnn`** | binary | [LDIRECT](running.md#ldirect)=-1 [LINIT_COND](running.md#linit_cond)>0 | 3D sensitivity of receptor concentrations and deposition to initial conditions |
| ** *Nested gridded data* ** | | | |
| **`grid_conc_nest_date_nnn`** | binary | [NESTED_OUTPUT](running.md#nested_output)=1 [LDIRECT](running.md#ldirect)=1 [IOUT](running.md#iout)=1,3,5 | 3D tracer mass density and 2D deposition |
| **`grid_pptv_nest_date_nnn`** | binary | [NESTED_OUTPUT](running.md#nested_output)=1 [LDIRECT](running.md#ldirect)=1 [IOUT](running.md#iout)=2,3 | 3D tracer volume mixing ratio and 2D deposition |
| **`grid_time_nest_date_nnn`** | binary | [NESTED_OUTPUT](running.md#nested_output)=1 [LDIRECT](running.md#ldirect)=-1 [IOUT](running.md#iout)=1 | 3D sensitivity of atmospheric receptor to emissions |
| **`grid_drydep_nest_date_nnn`** | binary | [NESTED_OUTPUT](running.md#nested_output)=1 [LDIRECT](running.md#ldirect)=-1 [IOUT](running.md#iout)=1 [IND_RECEPTOR](running.md#ind_receptor)=3 | 3D sensitivity of dry deposition receptor to emissions |
| **`grid_wetdep_nest_date_nnn`** | binary | [NESTED_OUTPUT](running.md#nested_output)=1 [LDIRECT](running.md#ldirect)=-1 [IOUT](running.md#iout)=1 [IND_RECEPTOR](running.md#ind_receptor)=4 | 3D sensitivity of wet deposition receptor to emissions |
| **`grid_conc_nest_date.nc`** | NetCDF | [NESTED_OUTPUT](running.md#nested_output)=1 [LDIRECT](running.md#ldirect)=1 [IOUT](running.md#iout)=9,10,11,13 | 3D tracer and 2D wet and dry deposition |
| **`grid_time_nest_date.nc`** | NetCDF | [NESTED_OUTPUT](running.md#nested_output)=1 [LDIRECT](running.md#ldirect)=-1 [IOUT](running.md#iout)=9 | 3D sensitivity of atmospheric receptor to emissions |
| **`grid_drydep_nest_date.nc`** | NetCDF | [NESTED_OUTPUT](running.md#nested_output)=1 [LDIRECT](running.md#ldirect)=-1 [IOUT](running.md#iout)=9 [IND_RECEPTOR](running.md#ind_receptor)=3 | 3D sensitivity of dry deposition receptor to emissions |
| **`grid_wetdep_nest_date.nc`** | NetCDF | [NESTED_OUTPUT](running.md#nested_output)=1 [LDIRECT](running.md#ldirect)=-1 [IOUT](running.md#iout)=9 [IND_RECEPTOR](running.md#ind_receptor)=4 | 3D sensitivity of wet deposition receptor to emissions |
| ** *Particle data* ** | | | |
| **`partoutput_date.nc`** | NetCDF | [IPOUT](running.md#ipout)=1,2,3 | Data at particle level. Output fields set in [PARTOPTIONS](running.md#partoptions) |
| **`partinit_date.nc`** | NetCDF | [IPOUT](running.md#ipout)=1,2,3 [IPIN](running.md#ipin)=1 | Initial particle data at t=0. Output fields set in [PARTOPTIONS](running.md#partoptions) |
| **`restart_date`** | binary | [LOUTRESTART](running.md#loutrestart)>=0 | Binary restart files are written to file at each [LOUTRESTART](running.md#loutrestart) interval, see [Restarting a simulation](running.md#restart) |