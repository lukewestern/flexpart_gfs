! SPDX-FileCopyrightText: FLEXPART 1998-2019, see flexpart_license.txt
! SPDX-License-Identifier: GPL-3.0-or-later

program flexpart

  !*****************************************************************************
  !                                                                            *
  !     This is the Lagrangian Particle Dispersion Model FLEXPART.             *
  !     The main program manages the reading of model run specifications, etc. *
  !     All actual computing is done within subroutine timemanager.            *
  !                                                                            *
  !     Author: A. Stohl                                                       *
  !                                                                            *
  !     18 May 1996                                                            *
  !                                                                            *
  !*****************************************************************************
  ! Changes:                                                                   *
  !   Unified ECMWF and GFS builds                                             *
  !   Marian Harustak, 12.5.2017                                               *
  !     - Added detection of metdata format using gributils routines           *
  !     - Distinguished calls to ecmwf/gfs gridcheck versions based on         *
  !       detected metdata format                                              *
  !     - Passed metdata format down to timemanager                            *
  !*****************************************************************************
  !                                                                            *
  ! Variables:                                                                 *
  !                                                                            *
  ! Constants:                                                                 *
  !                                                                            *
  !*****************************************************************************
  use omp_lib, only: OMP_GET_MAX_THREADS
  use par_mod
  use com_mod
  use random_mod, only: gasdev1
  use timemanager_mod
  use output_mod

  implicit none

  integer ::              &
    i,                    & ! loop variable for random numbers
    idummy=-320             ! dummy value used by the random routine
  character(len=256) ::   &
    inline_options          ! pathfile, flexversion, arg2

  ! Keeping track of the total running time of FLEXPART, printed out at the end.
  !*****************************************************************************
  CALL SYSTEM_CLOCK(count_clock, count_rate, count_max)
  s_total = (count_clock - count_clock0)/real(count_rate)

  ! Generate a large number of random numbers
  !******************************************
  do i=1,maxrand-1,2
    call gasdev1(idummy,rannumb(i),rannumb(i+1))
  end do
  call gasdev1(idummy,rannumb(maxrand),rannumb(maxrand-1))


  ! FLEXPART version string
  flexversion_major = '10' ! Major version number, also used for species file names
  flexversion='Version '//trim(flexversion_major)//'.4 (2019-11-12)'
  verbosity=0

  ! Read the pathnames where input/output files are stored
  !*******************************************************

  inline_options='none'
  select case (iargc())
  case (2)
    call getarg(1,arg1)
    pathfile=arg1
    call getarg(2,arg2)
    inline_options=arg2
  case (1)
    call getarg(1,arg1)
    pathfile=arg1
    if (arg1(1:1).eq.'-') then
      write(pathfile,'(a11)') './pathnames'
      inline_options=arg1 
    endif
  case (0)
    write(pathfile,'(a11)') './pathnames'
  end select
  
  ! Print the GPL License statement
  !*******************************************************
  print*,'Welcome to FLEXPART ', trim(flexversion)
  print*,'FLEXPART is free software released under the GNU General Public License.'
  write(*,*) 'FLEXPART is running with ', trim(wind_coord_type), 'coordinates.'
  ! Reading the number of threads available and print them for user
  !****************************************************************
#ifdef _OPENMP
    numthreads = OMP_GET_MAX_THREADS()
    numthreads = min(40,numthreads)
#else
    numthreads = 1
#endif

  if (numthreads.gt.1) then
    write(*,*)
    write(*,*) "*********** WARNING  *********************************"
    write(*,*) "* FLEXPART running in parallel mode                  *"
    write(*,*) "* Number of uncertainty classes in                   *"
    write(*,*) "* set to number of threads:", numthreads, ".         *"
    write(*,*) "******************************************************"
    write(*,*)
  endif

  ! Reading user specified options, allocating fields and checking bounds
  !**********************************************************************
  call read_options_and_initialise_flexpart

  ! Inform whether output kernel is used or not
  !*********************************************
  if (lroot) then
    if (.not.lusekerneloutput) then
      write(*,*) "Concentrations are calculated without using kernel"
    else
      write(*,*) "Concentrations are calculated using kernel"
    end if
  end if

  if (turboff) write(*,*) 'Turbulence switched off'

  ! For continuation of previous run or from user defined initial 
  ! conditions, read in particle positions
  !*************************************************************************
  call initialise_warm_start

  ! Calculate particle trajectories
  !********************************
  call timemanager

  CALL SYSTEM_CLOCK(count_clock, count_rate, count_max)
  s_total = (count_clock - count_clock0)/real(count_rate) - s_total
  
  write(*,*) 'Read wind fields: ', s_readwind, ' seconds'
  write(*,*) 'Write particle average files: ', s_writepartav, ' seconds'
  write(*,*) 'Write particle files: ', s_writepart, ' seconds'
  write(*,*) 'Total running time: ', s_total, ' seconds'
  write(*,*) 'CONGRATULATIONS: YOU HAVE SUCCESSFULLY COMPLETED A FLE&
       &XPART MODEL RUN!'

end program flexpart


subroutine read_options_and_initialise_flexpart

  use point_mod
  use par_mod
  use com_mod
  use conv_mod
  use class_gribfile
  use readoptions
  use windfields_mod
  use plume_mod
  use initialise_mod
  use drydepo_mod
  use getfields_mod
  use interpol_mod, only: interpol_allocate
  use outg_mod
  use binary_output_mod

  implicit none

  integer ::              &
    i,                    & ! loop variable for number of points
    inest                   ! loop variable for nested gridcells

  ! Read pathnames from file in working director that specify I/O directories
  !**************************************************************************
  call readpaths

  ! Read the user specifications for the current model run
  !*******************************************************
  call readcommand

  ! Read the age classes to be used
  !********************************
  call readageclasses

  ! Read, which wind fields are available within the modelling period
  !******************************************************************
  call readavailable

  if (ipout.ne.0) call readpartoptions

  ! Detect metdata format
  !**********************
  call detectformat

  if (metdata_format.eq.GRIBFILE_CENTRE_ECMWF) then
    print *,'ECMWF metdata detected'
  elseif (metdata_format.eq.GRIBFILE_CENTRE_NCEP) then
    print *,'NCEP metdata detected'
  else
    print *,'Unknown metdata format'
    stop
  endif

  ! Allocate memory for windfields
  !*******************************
  call windfields_allocate

  ! Read the model grid specifications,
  ! both for the mother domain and eventual nests
  !**********************************************
  if (metdata_format.eq.GRIBFILE_CENTRE_ECMWF) then
    call gridcheck_ecmwf
  else 
    call gridcheck_gfs    
  endif

  ! Set the upper level for where the convection will be working
  !*************************************************************
  call set_upperlevel_convect

  if (numbnests.ge.1) then
  ! If nested wind fields are used, allocate arrays
  !************************************************
    call windfields_nest_allocate
    call gridcheck_nests
  endif

  ! Read the output grid specifications if requested by user
  !*********************************************************
  if (iout.ne.0) then
    call readoutgrid

    if (nested_output.eq.1) then
      call readoutgrid_nest
    endif
  endif

  ! Read the receptor points for which extra concentrations are to be calculated
  !*****************************************************************************
  call readreceptors

  ! Read the physico-chemical species property table
  !*************************************************
  !SEC: now only needed SPECIES are read in readreleases.f
  !call readspecies

  ! Read the landuse inventory
  !***************************
  call readlanduse ! CHECK ETA

  ! Read the coordinates of the release locations
  !**********************************************
  call readreleases ! CHECK ETA
  ! Convert the release point coordinates from geografical to grid coordinates
  !***************************************************************************
  call coordtrafo(nxmin1,nymin1) ! CHECK ETA

  ! Read and compute surface resistances to dry deposition of gases
  !****************************************************************
  call readdepo ! CHECK ETA

  ! Allocate dry deposition fields if necessary
  !*********************************************
  call drydepo_allocate
  call convection_allocate
  call getfields_allocate
  call interpol_allocate

  ! Assign fractional cover of landuse classes to each ECMWF grid point
  !********************************************************************
  call assignland ! CHECK ETA

  ! Calculate volume, surface area, etc., of all output grid cells
  ! Allocate fluxes and OHfield if necessary
  !***************************************************************
  if (iout.ne.0) then
    call outgrid_init ! CHECK ETA
    if (nested_output.eq.1) call outgrid_init_nest ! CHECK ETA
  endif

  ! Read the OH field
  !******************
  if (OHREA.eqv..TRUE.) then
    call readOHfield ! CHECK ETA
  endif

  call openreceptors ! CHECK ETA
  if ((iout.eq.4).or.(iout.eq.5)) call openouttraj ! CHECK ETA

  ! Releases can only start and end at discrete times (multiples of lsynctime)
  !***************************************************************************
  do i=1,numpoint
    ireleasestart(i)=nint(real(ireleasestart(i))/real(lsynctime))*lsynctime
    ireleaseend(i)=nint(real(ireleaseend(i))/real(lsynctime))*lsynctime
  end do

  ! Initialize cloud-base mass fluxes for the convection scheme
  !************************************************************

  cbaseflux(0:nxmin1,0:nymin1)=0.
  do inest=1,numbnests
    cbasefluxn(0:nxn(inest)-1,0:nyn(inest)-1,inest)=0.
  end do

end subroutine read_options_and_initialise_flexpart


subroutine initialise_warm_start
  use com_mod
  use initialise_mod
  use netcdf_output_mod

  implicit none

  itime_init=0
  if ((ipin.eq.1).or.(ipin.eq.4)) then ! Restarting from restart.bin file
    call readrestart
  else if (ipin.eq.2) then ! Restarting from netcdf partoutput file
#ifdef USE_NCF
    call readpartpositions
#else
    stop 'Compile with netCDF if you want to use the ipin=2 option.'
#endif
  else if (ipin.eq.3) then ! User defined particle properties
  ! Reading initial conditions from netcdf file
#ifdef USE_NCF
    call readinitconditions_netcdf
#else
    stop 'Compile with netCDF if you want to use the ipin=3 option.'
#endif
  else
    numpart=0
    numparticlecount=0
  endif

end subroutine initialise_warm_start