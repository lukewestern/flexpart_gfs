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
  use point_mod
  use par_mod
  use com_mod
  use conv_mod
  use random_mod, only: gasdev1
  use class_gribfile
  use readoptions
  use windfields_mod
  use timemanager_mod

#ifdef USE_NCF
  use netcdf_output_mod, only: writeheader_netcdf
#endif
  use binary_output_mod
  use txt_output_mod

  implicit none

  integer ::              &
    i,j,                  & ! loop variables
    ix,jy,                & ! grid indices
    inest,                & ! loop variable for nested gridcells
    iopt,                 & ! temporarily storing inline options
    detectformat,         & ! integer function in detectformat.f90 (should this be here?)
    idummy=-320,          & ! dummy value used by the random routine
    metdata_format=GRIBFILE_CENTRE_UNKNOWN  ! storing the input data type (ECMWF/NCEP)
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

  ! Set the coordinate system. At the moment only ECMWF is possible. This bit needs
  ! to be a parameter that can be set at compile time. Throughout the code there
  ! will be select cases statements or ifdefs
  !*****************************************************************
  wind_coord_type='ETA'
  !wind_coord_type='METER'
  
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
  
  ! Read pathnames from file in working director that specify I/O directories
  !**************************************************************************
  call readpaths


  ! Read the user specifications for the current model run
  !*******************************************************
  call readcommand

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
 
  ! Initialize arrays in com_mod
  !*****************************
  ! call com_mod_allocate_part(maxpart)

  ! Read the age classes to be used
  !********************************
  call readageclasses

  ! Read, which wind fields are available within the modelling period
  !******************************************************************
  call readavailable

  ! Detect metdata format. GFS not supported, but can be added if converted
  ! to ECMWF eta coordinates or the appropriate parts in advance and the
  ! interpolation subroutines need to be changed.
  !**********************
  metdata_format = detectformat()

  if (metdata_format.eq.GRIBFILE_CENTRE_ECMWF) then
    print *,'ECMWF metdata detected'
  ! elseif (metdata_format.eq.GRIBFILE_CENTRE_NCEP) then
  !   print *,'NCEP metdata detected'
  else
    print *,'Unknown metdata format'
    stop
  endif

  ! If nested wind fields are used, allocate arrays
  !************************************************
  call com_mod_allocate_nests

  ! Read the model grid specifications,
  ! both for the mother domain and eventual nests
  !**********************************************
  call gridcheck_ecmwf
  call set_upperlevel_convect
  call gridcheck_nests

  ! Read the output grid specifications if requested by user
  !*********************************************************
  if (grid_output.eq.1) then
    call readoutgrid

    if (nested_output.eq.1) then
      call readoutgrid_nest
    endif
  endif

  ! Read the receptor points for which extra concentrations are to be calculated
  !*****************************************************************************
  call readreceptors ! CHECK ETA

  ! Read the physico-chemical species property table
  !*************************************************
  !SEC: now only needed SPECIES are read in readreleases.f
  !call readspecies

  ! Read the landuse inventory
  !***************************
  call readlanduse ! CHECK ETA

  ! Assign fractional cover of landuse classes to each ECMWF grid point
  !********************************************************************
  call assignland ! CHECK ETA

  ! Read the coordinates of the release locations
  !**********************************************
  call readreleases ! CHECK ETA

  ! Read and compute surface resistances to dry deposition of gases
  !****************************************************************
  call readdepo ! CHECK ETA

  ! Convert the release point coordinates from geografical to grid coordinates
  !***************************************************************************
  call coordtrafo ! CHECK ETA

  ! For continuation of previous run, read in particle positions
  !*************************************************************
  if (ipin.eq.1) then
    !stop "Convert readpartpositions to netcdf and then use particle module"
    call readpartpositions
  else
    numpart=0
    numparticlecount=0
  endif

  ! Calculate volume, surface area, etc., of all output grid cells
  ! Allocate fluxes and OHfield if necessary
  !***************************************************************
  if (grid_output.eq.1) then
    call outgrid_init ! CHECK ETA
    if (nested_output.eq.1) call outgrid_init_nest ! CHECK ETA
  endif

  ! Read the OH field
  !******************
  if (OHREA.eqv..TRUE.) then
    call readOHfield ! CHECK ETA
  endif

  ! Write basic information on the simulation to a file "header"
  ! and open files that are to be kept open throughout the simulation
  !******************************************************************
  if (grid_output.eq.1) then
#ifdef USE_NCF
    if (lnetcdfout.eq.1) then 
      call writeheader_netcdf(lnest=.false.)
    else 
      call writeheader
    end if

    if (nested_output.eq.1) then
      if (lnetcdfout.eq.1) then
        call writeheader_netcdf(lnest=.true.)
      else
        call writeheader_nest
      endif
    endif
#endif
  endif

  call writeheader ! CHECK ETA
  ! FLEXPART 9.2 ticket ?? write header in ASCII format 
  call writeheader_txt
  !if (nested_output.eq.1) call writeheader_nest
  if (nested_output.eq.1.and.surf_only.ne.1) call writeheader_nest
  if (nested_output.eq.1.and.surf_only.eq.1) call writeheader_nest_surf
  if (nested_output.ne.1.and.surf_only.eq.1) call writeheader_surf

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

  do jy=0,nymin1
    do ix=0,nxmin1
      cbaseflux(ix,jy)=0.
    end do
  end do
  do inest=1,numbnests
    do jy=0,nyn(inest)-1
      do ix=0,nxn(inest)-1
        cbasefluxn(ix,jy,inest)=0.
      end do
    end do
  end do

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
  ! Calculate particle trajectories
  !********************************
  call timemanager(metdata_format)

  CALL SYSTEM_CLOCK(count_clock, count_rate, count_max)
  s_total = (count_clock - count_clock0)/real(count_rate) - s_total
  
  write(*,*) 'Read wind fields: ', s_readwind, ' seconds'
  write(*,*) 'Write particle average files: ', s_writepartav, ' seconds'
  write(*,*) 'Write particle files: ', s_writepart, ' seconds'
  write(*,*) 'Total running time: ', s_total, ' seconds'
  write(*,*) 'CONGRATULATIONS: YOU HAVE SUCCESSFULLY COMPLETED A FLE&
       &XPART MODEL RUN!'

end program flexpart
