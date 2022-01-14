! SPDX-FileCopyrightText: FLEXPART 1998-2019, see flexpart_license.txt
! SPDX-License-Identifier: GPL-3.0-or-later

  !*****************************************************************************
  !                                                                            *
  !  This module handles all gridded netcdf output for concentration or        *
  !  residence time and wet and dry deposition output.                         *
  !                                                                            *
  !  - writeheader_netcdf generates file including all information previously  *
  !    stored in separate header files                                         *
  !  - concoutput_netcdf write concentration output and wet and dry deposition *
  !                                                                            *
  !     Author: D. Brunner                                                     *
  !                                                                            *
  !     12 April 2013                                                          *
  !                                                                            *
  ! HSO: 21 Oct 2014
  !  - added option to not writeout releases information by changing 
  !    switch write_releases
  !  - additional updates for FLEXPART 9.x
  ! 
  ! ESO 2016
  !  - Deposition fields can be calculated in double precision, see variable
  !    'dep_prec' in par_mod
  !  - Hardcoded options 'write_vol' and 'write_area' for grid cell
  !    volume and area
  !*****************************************************************************


module netcdf_output_mod

  use netcdf

  use point_mod, only: ireleasestart,ireleaseend,kindz,dx,xlon0,dy,ylat0,&
                       xpoint1,ypoint1,xpoint2,ypoint2,zpoint1,zpoint2,npart,xmass
  use outg_mod,  only: outheight,oroout,densityoutgrid,factor3d,volume,&
                       wetgrid,wetgridsigma,drygrid,drygridsigma,grid,gridsigma,&
                       area,arean,volumen,orooutn
  use par_mod,   only: dep_prec, sp, dp, maxspec, maxreceptor, nclassunc,&
                       unitoutrecept,unitoutreceptppt,unittmp
  use com_mod,   only: path,length,ldirect,ibdate,ibtime,iedate,ietime, &
                       loutstep,loutaver,loutsample,outlon0,outlat0,&
                       numxgrid,numygrid,dxout,dyout,numzgrid, &
                       outlon0n,outlat0n,dxoutn,dyoutn,numxgridn,numygridn, &
                       nspec,maxpointspec_act,species,numpoint,&
                       compoint,method,lsubgrid,lconvection,&
                       ind_source,ind_receptor,nageclass,lage,&
                       drydep,wetdep,decay,weta_gas,wetb_gas, numbnests, &
                       ccn_aero,in_aero, mintime, & ! wetc_in,wetd_in, &
                       reldiff,henry,f0,density,dquer,dsigma,dryvel,&
                       weightmolar,ohcconst,ohdconst,vsetaver,&

                       memind,xreceptor,yreceptor,numreceptor,creceptor,iout, &
                       itsplit, lsynctime, ctl, ifine, lagespectra, ipin, &
                       ioutputforeachrelease, iflux, mdomainfill, mquasilag, & 
                       nested_output, ipout, surf_only, linit_cond, &
                       flexversion,mpi_mode,DRYBKDEP,WETBKDEP,numpart,numpoint
  use windfields_mod, only: oro,rho,nxmax,height,nxmin1,nymin1,nz, &
                       ! for concoutput_netcdf and concoutput_nest_netcdf 
                       oron,rhon,xresoln,yresoln,xrn,xln,yrn,yln,nxn,nyn

  use mean_mod

  implicit none

  private

  public :: writeheader_netcdf,concoutput_surf_nest_netcdf,concoutput_netcdf,&
       &concoutput_nest_netcdf,concoutput_surf_netcdf,writeheader_partoutput,partoutput_netcdf,&
       open_partoutput_file,close_partoutput_file,readpartpositions_netcdf,create_particles_initialoutput,&
       write_particles_initialoutput

  !  include 'netcdf.inc'

  ! parameter for data compression (1-9, 9 = most aggressive)
  integer, parameter :: deflate_level = 5
  logical, parameter :: min_size = .false.   ! if set true, redundant fields (topography) are not written to minimize file size
  character(len=255), parameter :: institution = 'NILU'

  integer            :: tpointer=0,tpointer_part=0,ppointer_part=0,partinitpointer=0
  character(len=255) :: ncfname, ncfnamen, ncfname_part, ncfname_partinit!(maxpoint)

  ! netcdf dimension and variable IDs for main and nested output grid
  integer, dimension(maxspec) :: specID,specIDppt, wdspecID,ddspecID
  integer, dimension(maxspec) :: specIDn,specIDnppt, wdspecIDn,ddspecIDn
  integer                     :: timeID, timeIDn, timeIDpart
  integer, dimension(6)       :: dimids, dimidsn
  integer, dimension(5)       :: depdimids, depdimidsn

  !IDs for partoutput
  integer             :: partID
  integer             :: itramemID,topoID,pvID,qvID,rhoID
  integer             :: hmixID,trID,ttID,lonID,latID,levID,massID(maxspec)
  integer             :: partIDi,tIDi,lonIDi,latIDi,levIDi ! For initial particle outputs

  real :: eps
  !  private:: writemetadata, output_units, nf90_err

  ! switch output of release point information on/off
  logical, parameter :: write_releases = .true.

  ! switch output of grid cell volume and area on/off
  logical, parameter :: write_vol = .false.
  logical, parameter :: write_area = .false.

  ! coordinate transformation from internal to world coord
  real :: xp1,yp1,xp2,yp2

contains

!****************************************************************
! determine output units (see table 1 in Stohl et al., ACP 2005
!****************************************************************
subroutine output_units(units)
  implicit none
  character(len=15), intent(out) :: units
  if (ldirect.eq.1) then
     ! forward simulation
     if (ind_source.eq.1) then
        if (ind_receptor.eq.1) then
           units = 'ng m-3'   ! hes the kg in Tab1 is only indicating the units of the relase not the output
        else
           units = 'ng kg-1'
        endif
     else
        if (ind_receptor.eq.1) then
           units = 'ng m-3'
        else
           units = 'ng kg-1'
        endif
     endif
  else
     ! backward simulation
     if (ind_source.eq.1) then
        if (ind_receptor.eq.1) then
           units = 's'
        else 
           units = 's m3 kg-1'
        endif
     else
        if (ind_receptor.eq.1) then
           units = 's kg m-3'
        else
           units = 's'
        endif
     endif
  endif
end subroutine output_units


!****************************************************************
! write metadata to netCDF file 
!****************************************************************
subroutine writemetadata(ncid,lnest)
  
  implicit none 

  integer, intent(in) :: ncid
  logical, intent(in) :: lnest
  integer             :: status
  character           :: time*10,date*8,adate*8,atime*6
  character(5)        :: zone
  character(255)      :: login_name, host_name

  ! gather system information 
  call date_and_time(date,time,zone)
  call getlog(login_name)
  call hostnm(host_name)
  
  ! hes CF convention requires these attributes
  call nf90_err(nf90_put_att(ncid, nf90_global, 'Conventions', 'CF-1.6'))
  call nf90_err(nf90_put_att(ncid, nf90_global, 'title', 'FLEXPART model output'))
  call nf90_err(nf90_put_att(ncid, nf90_global, 'institution', trim(institution)))
  call nf90_err(nf90_put_att(ncid, nf90_global, 'source', trim(flexversion)//' model output'))
  call nf90_err(nf90_put_att(ncid, nf90_global, 'history', date(1:4)//'-'//date(5:6)// &
       '-'//date(7:8)//' '//time(1:2)//':'//time(3:4)//' '//zone//'  created by '//  &
       trim(login_name)//' on '//trim(host_name)))
  call nf90_err(nf90_put_att(ncid, nf90_global, 'references', &
       'Stohl et al., Atmos. Chem. Phys., 2005, doi:10.5194/acp-5-2461-200'))

  ! attributes describing model run
  !************************************************************************************

  if (lnest) then
     call nf90_err(nf90_put_att(ncid, nf90_global, 'outlon0', outlon0n))
     call nf90_err(nf90_put_att(ncid, nf90_global, 'outlat0', outlat0n))
     call nf90_err(nf90_put_att(ncid, nf90_global, 'dxout', dxoutn))
     call nf90_err(nf90_put_att(ncid, nf90_global, 'dyout', dyoutn))
  else
     call nf90_err(nf90_put_att(ncid, nf90_global, 'outlon0', outlon0))
     call nf90_err(nf90_put_att(ncid, nf90_global, 'outlat0', outlat0))
     call nf90_err(nf90_put_att(ncid, nf90_global, 'dxout', dxout))
     call nf90_err(nf90_put_att(ncid, nf90_global, 'dyout', dyout))
  endif
  !   vertical levels stored in grid structure

  ! COMMAND file settings
  call nf90_err(nf90_put_att(ncid, nf90_global, 'ldirect', ldirect))
  write(adate,'(i8.8)') ibdate
  write(atime,'(i6.6)') ibtime
  call nf90_err(nf90_put_att(ncid, nf90_global, 'ibdate', adate))
  call nf90_err(nf90_put_att(ncid, nf90_global, 'ibtime', atime))
  write(adate,'(i8.8)') iedate
  write(atime,'(i6.6)') ietime
  call nf90_err(nf90_put_att(ncid, nf90_global, 'iedate', adate))
  call nf90_err(nf90_put_att(ncid, nf90_global, 'ietime', atime))
  call nf90_err(nf90_put_att(ncid, nf90_global, 'loutstep', loutstep))
  call nf90_err(nf90_put_att(ncid, nf90_global, 'loutaver', loutaver))
  call nf90_err(nf90_put_att(ncid, nf90_global, 'loutsample', loutsample))
  call nf90_err(nf90_put_att(ncid, nf90_global, 'itsplit', itsplit))
  call nf90_err(nf90_put_att(ncid, nf90_global, 'lsynctime', lsynctime))
  call nf90_err(nf90_put_att(ncid, nf90_global, 'ctl', ctl))
  call nf90_err(nf90_put_att(ncid, nf90_global, 'ifine', ifine))
  call nf90_err(nf90_put_att(ncid, nf90_global, 'iout', iout))
  call nf90_err(nf90_put_att(ncid, nf90_global, 'ipout', ipout))
  call nf90_err(nf90_put_att(ncid, nf90_global, 'lsubgrid', lsubgrid))
  call nf90_err(nf90_put_att(ncid, nf90_global, 'lconvection', lconvection))
  call nf90_err(nf90_put_att(ncid, nf90_global, 'lagespectra', lagespectra))
  call nf90_err(nf90_put_att(ncid, nf90_global, 'ipin', ipin))
  call nf90_err(nf90_put_att(ncid, nf90_global, 'ioutputforeachrelease', ioutputforeachrelease))
  call nf90_err(nf90_put_att(ncid, nf90_global, 'iflux', iflux))
  call nf90_err(nf90_put_att(ncid, nf90_global, 'mdomainfill', mdomainfill))
  call nf90_err(nf90_put_att(ncid, nf90_global, 'ind_source', ind_source))
  call nf90_err(nf90_put_att(ncid, nf90_global, 'ind_receptor', ind_receptor))
  call nf90_err(nf90_put_att(ncid, nf90_global, 'mquasilag', mquasilag))
  call nf90_err(nf90_put_att(ncid, nf90_global, 'nested_output', nested_output))
  call nf90_err(nf90_put_att(ncid, nf90_global, 'surf_only', surf_only))
  call nf90_err(nf90_put_att(ncid, nf90_global, 'linit_cond', linit_cond))
end subroutine writemetadata


!****************************************************************
! netcdf error message handling
!****************************************************************
subroutine nf90_err(status)
  implicit none
  integer, intent (in) :: status
   if(status /= nf90_noerr) then
      print *, trim(nf90_strerror(status))
      stop 'Stopped'
    end if
end subroutine nf90_err


!****************************************************************
! Create netcdf file and write header/metadata information
! lnest = .false. : Create main output file
! lnest = .true.  : Create nested output file
!****************************************************************
subroutine writeheader_netcdf(lnest)

  implicit none

  logical, intent(in) :: lnest

  integer :: ncid, sID, wdsID, ddsID
  integer :: timeDimID, latDimID, lonDimID, levDimID
  integer :: nspecDimID, npointDimID, nageclassDimID, ncharDimID, pointspecDimID
  integer :: tID, lonID, latID, levID, poleID, lageID, oroID
  integer :: volID, areaID
  integer :: rellng1ID, rellng2ID, rellat1ID, rellat2ID, relzz1ID, relzz2ID
  integer :: relcomID, relkindzID, relstartID, relendID, relpartID, relxmassID
  integer :: nnx, nny 
  integer, dimension(6)       :: dIDs
  integer, dimension(5)       :: depdIDs
  character(len=255)          :: fname
  character(len=15)           :: units
  character(len=20)           :: fprefix
  character(len=3)            :: anspec
  CHARACTER                   :: adate*8,atime*6,timeunit*32
  !REAL, DIMENSION(1000)       :: coord
  real, allocatable, dimension(:) :: coord

  integer                     :: cache_size
  integer, dimension(6)       :: chunksizes
  integer, dimension(5)       :: dep_chunksizes

  integer                     :: i,ix,jy
  integer                     :: test_unit


  ! Check if output directory exists (the netcdf library will
  ! otherwise give an error which can look confusing). 
  ! *********************************************************************
  open(unit=unittmp,file=trim(path(2)(1:length(2)))//'test_dir.txt',status='replace',&
       &err=100)
  close (unittmp, status='delete')
  goto 101
100 write(*,FMT='(80("#"))') 
  write(*,*) 'ERROR: output directory ', trim(path(2)(1:length(2))), ' does not exist&
       & (or failed to write there).' 
  write(*,*) 'EXITING' 
  write(*,FMT='(80("#"))')
  stop
101 continue

  !************************
  ! Create netcdf file
  !************************

  if (ldirect.eq.1) then
     write(adate,'(i8.8)') ibdate
     write(atime,'(i6.6)') ibtime
     fprefix = 'grid_conc_'
  else
     write(adate,'(i8.8)') iedate
     write(atime,'(i6.6)') ietime
     fprefix = 'grid_time_'
  endif
  if (DRYBKDEP) fprefix='grid_drydep_'
  if (WETBKDEP) fprefix='grid_wetdep_'

  if (lnest) then
     fname = path(2)(1:length(2))//trim(fprefix)//adate//atime//'_nest.nc'
     ncfnamen = fname
     nnx = numxgridn
     nny = numygridn
  else
     fname = path(2)(1:length(2))//trim(fprefix)//adate//atime//'.nc'
     ncfname = fname
     nnx = numxgrid
     nny = numygrid
  endif

  cache_size = 16 * nnx * nny * numzgrid

  ! setting cache size in bytes. It is set to 4 times the largest data block that is written
  !   size_type x nx x ny x nz
  ! create file
  call nf90_err(nf90_create(trim(fname), cmode = nf90_hdf5, ncid = ncid, &
    cache_size = cache_size))  

  ! create dimensions:
  !*************************
  ! time
  call nf90_err(nf90_def_dim(ncid, 'time', nf90_unlimited, timeDimID))
  timeunit = 'seconds since '//adate(1:4)//'-'//adate(5:6)// &
     '-'//adate(7:8)//' '//atime(1:2)//':'//atime(3:4)

  ! lon
  call nf90_err(nf90_def_dim(ncid, 'longitude', nnx, lonDimID))
  ! lat
  call nf90_err(nf90_def_dim(ncid, 'latitude', nny, latDimID))
  ! level
  call nf90_err(nf90_def_dim(ncid, 'height', numzgrid, levDimID))
  ! number of species
  call nf90_err(nf90_def_dim(ncid, 'numspec', nspec, nspecDimID))
  ! number of release points
  call nf90_err(nf90_def_dim(ncid, 'pointspec', maxpointspec_act, pointspecDimID))
  ! number of age classes
  call nf90_err(nf90_def_dim(ncid, 'nageclass', nageclass, nageclassDimID))
  ! dimension for release point characters
  call nf90_err(nf90_def_dim(ncid, 'nchar', 45, ncharDimID))
  ! number of actual release points
  call nf90_err(nf90_def_dim(ncid, 'numpoint', numpoint, npointDimID))


  ! create variables
  !*************************

  ! time
  call nf90_err(nf90_def_var(ncid, 'time', nf90_int, (/ timeDimID /), tID))
  call nf90_err(nf90_put_att(ncid, tID, 'units', timeunit))
  call nf90_err(nf90_put_att(ncid, tID, 'calendar', 'proleptic_gregorian'))
  if (lnest) then
     timeIDn = tID
  else
     timeID = tID
  endif

  ! lon
  call nf90_err(nf90_def_var(ncid, 'longitude', nf90_float, (/ lonDimID /), lonID))
  call nf90_err(nf90_put_att(ncid, lonID, 'long_name', 'longitude in degree east'))
  call nf90_err(nf90_put_att(ncid, lonID, 'axis', 'Lon'))
  call nf90_err(nf90_put_att(ncid, lonID, 'units', 'degrees_east'))
  call nf90_err(nf90_put_att(ncid, lonID, 'standard_name', 'grid_longitude'))
  call nf90_err(nf90_put_att(ncid, lonID, 'description', 'grid cell centers'))

  ! lat
  call nf90_err(nf90_def_var(ncid, 'latitude', nf90_float, (/ latDimID /), latID))
  call nf90_err(nf90_put_att(ncid, latID, 'long_name', 'latitude in degree north'))
  call nf90_err(nf90_put_att(ncid, latID, 'axis', 'Lat'))
  call nf90_err(nf90_put_att(ncid, latID, 'units', 'degrees_north'))
  call nf90_err(nf90_put_att(ncid, latID, 'standard_name', 'grid_latitude'))
  call nf90_err(nf90_put_att(ncid, latID, 'description', 'grid cell centers'))


  ! height
  call nf90_err(nf90_def_var(ncid, 'height', nf90_float, (/ levDimID /), levID))
  ! call nf90_err(nf90_put_att(ncid, levID, 'axis', 'Z'))
  call nf90_err(nf90_put_att(ncid, levID, 'units', 'meters'))
  call nf90_err(nf90_put_att(ncid, levID, 'positive', 'up'))
  call nf90_err(nf90_put_att(ncid, levID, 'standard_name', 'height'))
  call nf90_err(nf90_put_att(ncid, levID, 'long_name', 'height above ground'))

  ! volume
  if (write_vol) call nf90_err(nf90_def_var(ncid, 'volume', nf90_float, &
       &(/ lonDimID, latDimID, levDimID /), volID))
  ! area 
  if (write_area) call nf90_err(nf90_def_var(ncid, 'area', nf90_float, &
       &(/ lonDimID, latDimID /), areaID))


  if (write_releases.eqv..true.) then
    ! release comment
    call nf90_err(nf90_def_var(ncid, 'RELCOM', nf90_char, (/ ncharDimID,npointDimID /), &
         relcomID))
    call nf90_err(nf90_put_att(ncid, relcomID, 'long_name', 'release point name'))
    ! release longitude 1
    call nf90_err(nf90_def_var(ncid, 'RELLNG1', nf90_float, (/ npointDimID /), rellng1ID))
    call nf90_err(nf90_put_att(ncid, rellng1ID, 'units', 'degrees_east'))
    call nf90_err(nf90_put_att(ncid, rellng1ID, 'long_name', &
         'release longitude lower left corner'))
    ! release longitude 2
    call nf90_err(nf90_def_var(ncid, 'RELLNG2', nf90_float, (/ npointDimID /), rellng2ID))
    call nf90_err(nf90_put_att(ncid, rellng2ID, 'units', 'degrees_east'))
    call nf90_err(nf90_put_att(ncid, rellng2ID, 'long_name', &
         'release longitude upper right corner'))
    ! release latitude 1
    call nf90_err(nf90_def_var(ncid, 'RELLAT1', nf90_float, (/ npointDimID /), rellat1ID))
    call nf90_err(nf90_put_att(ncid, rellat1ID, 'units', 'degrees_north'))
    call nf90_err(nf90_put_att(ncid, rellat1ID, 'long_name', &
         'release latitude lower left corner'))
    ! release latitude 2
    call nf90_err(nf90_def_var(ncid, 'RELLAT2', nf90_float, (/ npointDimID /), rellat2ID))
    call nf90_err(nf90_put_att(ncid, rellat2ID, 'units', 'degrees_north'))
    call nf90_err(nf90_put_att(ncid, rellat2ID, 'long_name', &
         'release latitude upper right corner'))

    ! hes: if rotated_ll it would be convenient also to write the the release points in rotated_coordinates

    ! release height bottom
    call nf90_err(nf90_def_var(ncid, 'RELZZ1', nf90_float, (/ npointDimID /), relzz1ID))
    call nf90_err(nf90_put_att(ncid, relzz1ID, 'units', 'meters'))
    call nf90_err(nf90_put_att(ncid, relzz1ID, 'long_name', 'release height bottom'))
    ! release height top
    call nf90_err(nf90_def_var(ncid, 'RELZZ2', nf90_float, (/ npointDimID /), relzz2ID))
    call nf90_err(nf90_put_att(ncid, relzz2ID, 'units', 'meters'))
    call nf90_err(nf90_put_att(ncid, relzz2ID, 'long_name', 'release height top'))
    ! release kind
    call nf90_err(nf90_def_var(ncid, 'RELKINDZ', nf90_int, (/ npointDimID /), relkindzID))
    call nf90_err(nf90_put_att(ncid, relkindzID, 'long_name', 'release kind'))
    ! release start
    call nf90_err(nf90_def_var(ncid, 'RELSTART', nf90_int, (/ npointDimID /), relstartID))
    call nf90_err(nf90_put_att(ncid, relstartID, 'units', 'seconds'))
    call nf90_err(nf90_put_att(ncid, relstartID, 'long_name', &
         'release start relative to simulation start'))
    ! release end
    call nf90_err(nf90_def_var(ncid, 'RELEND', nf90_int, (/ npointDimID /), relendID))
    call nf90_err(nf90_put_att(ncid, relendID, 'units', 'seconds'))
    call nf90_err(nf90_put_att(ncid, relendID, 'long_name', &
         'release end relative to simulation start'))
    ! release particles
    call nf90_err(nf90_def_var(ncid, 'RELPART', nf90_int, (/ npointDimID /), relpartID))
    call nf90_err(nf90_put_att(ncid, relpartID, 'long_name', 'number of release particles'))
    ! release particle masses
    call nf90_err(nf90_def_var(ncid, 'RELXMASS', nf90_float, (/ npointDimID, nspecDimID /), &
         relxmassID))
    call nf90_err(nf90_put_att(ncid, relxmassID, 'long_name', 'total release particle mass'))
  end if
 
  ! age classes
  call nf90_err(nf90_def_var(ncid, 'LAGE', nf90_int, (/ nageclassDimID /), lageID))
  call nf90_err(nf90_put_att(ncid, lageID, 'units', 'seconds'))
  call nf90_err(nf90_put_att(ncid, lageID, 'long_name', 'age class'))

  ! output orography
  if (.not. min_size) then
    call nf90_err(nf90_def_var(ncid, 'ORO', nf90_int, (/ lonDimID, latDimID /), oroID,  &
      deflate_level=deflate_level, chunksizes= (/ nnx, nny /)))
    call nf90_err(nf90_put_att(ncid, oroID, 'standard_name', 'surface altitude'))
    call nf90_err(nf90_put_att(ncid, oroID, 'long_name', 'outgrid surface altitude'))
    call nf90_err(nf90_put_att(ncid, oroID, 'units', 'm'))
  end if

  ! concentration output, wet and dry deposition variables (one per species)
  call output_units(units)

  dIDs = (/ londimid, latdimid, levdimid, timedimid, pointspecdimid, nageclassdimid /)
  depdIDs = (/ londimid, latdimid, timedimid, pointspecdimid, nageclassdimid /)
  if (lnest) then
     dimidsn    = dIDs
     depdimidsn = depdIDs
  else
     dimids    = dIDs
     depdimids = depdIDs
  endif

  ! set chunksizes according to largest written portion of data in an individual call to 
  ! nf90_put_var
  chunksizes = (/ nnx, nny, numzgrid, 1, 1, 1 /)
  dep_chunksizes = (/ nnx, nny, 1, 1, 1 /)

  do i = 1,nspec
     write(anspec,'(i3.3)') i

     ! concentration output
     if ((iout.eq.1).or.(iout.eq.3).or.(iout.eq.5)) then
        call nf90_err(nf90_def_var(ncid,'spec'//anspec//'_mr', nf90_float, dIDs, sID , &
             deflate_level = deflate_level,  &
             chunksizes = chunksizes ))
        call nf90_err(nf90_put_att(ncid, sID, 'units', units))
        call nf90_err(nf90_put_att(ncid, sID, 'long_name', species(i)))
        call nf90_err(nf90_put_att(ncid, sID, 'decay', decay(i)))
        call nf90_err(nf90_put_att(ncid, sID, 'weightmolar', weightmolar(i)))
  !        call nf90_err(nf90_put_att(ncid, sID, 'ohreact', ohreact(i)))
        call nf90_err(nf90_put_att(ncid, sID, 'ohcconst', ohcconst(i)))
        call nf90_err(nf90_put_att(ncid, sID, 'ohdconst', ohdconst(i)))
        call nf90_err(nf90_put_att(ncid, sID, 'vsetaver', vsetaver(i)))

        if (lnest) then
           specIDn(i) = sID
        else
           specID(i) = sID
        endif
     endif

     ! mixing ratio output
     if ((iout.eq.2).or.(iout.eq.3)) then
        call nf90_err(nf90_def_var(ncid,'spec'//anspec//'_pptv', nf90_float, dIDs, sID , &
             deflate_level = deflate_level,  &
             chunksizes = chunksizes ))
        call nf90_err(nf90_put_att(ncid, sID, 'units', 'pptv'))
        call nf90_err(nf90_put_att(ncid, sID, 'long_name', species(i)))
        call nf90_err(nf90_put_att(ncid, sID, 'decay', decay(i)))
        call nf90_err(nf90_put_att(ncid, sID, 'weightmolar', weightmolar(i)))
  !        call nf90_err(nf90_put_att(ncid, sID, 'ohreact', ohreact(i)))
        call nf90_err(nf90_put_att(ncid, sID, 'ohcconst', ohcconst(i)))
        call nf90_err(nf90_put_att(ncid, sID, 'ohdconst', ohdconst(i)))
        call nf90_err(nf90_put_att(ncid, sID, 'vsetaver', vsetaver(i)))

        if (lnest) then
           specIDnppt(i) = sID
        else
           specIDppt(i) = sID
        endif
     endif

     ! wet and dry deposition fields for forward runs
     if (wetdep) then
        call nf90_err(nf90_def_var(ncid,'WD_spec'//anspec, nf90_float, depdIDs, &
             wdsID, deflate_level = deflate_level, &
             chunksizes = dep_chunksizes))
        call nf90_err(nf90_put_att(ncid, wdsID, 'units', '1e-12 kg m-2'))
        call nf90_err(nf90_put_att(ncid, wdsID, 'weta_gas', weta_gas(i)))
        call nf90_err(nf90_put_att(ncid, wdsID, 'wetb_gas', wetb_gas(i)))
        call nf90_err(nf90_put_att(ncid, wdsID, 'ccn_aero', ccn_aero(i)))
        call nf90_err(nf90_put_att(ncid, wdsID, 'in_aero', in_aero(i)))
        ! call nf90_err(nf90_put_att(ncid, wdsID, 'wetc_in', wetc_in(i)))
        ! call nf90_err(nf90_put_att(ncid, wdsID, 'wetd_in', wetd_in(i)))
        call nf90_err(nf90_put_att(ncid, wdsID, 'dquer', dquer(i)))
        call nf90_err(nf90_put_att(ncid, wdsID, 'henry', henry(i)))
        if (lnest) then
           wdspecIDn(i) = wdsID
        else
           wdspecID(i) = wdsID
        endif
     endif
     if (drydep) then
        call nf90_err(nf90_def_var(ncid,'DD_spec'//anspec, nf90_float, depdIDs, &
             ddsID, deflate_level = deflate_level, &
             chunksizes = dep_chunksizes))
        call nf90_err(nf90_put_att(ncid, ddsID, 'units', '1e-12 kg m-2'))
        call nf90_err(nf90_put_att(ncid, ddsID, 'dryvel', dryvel(i)))
        call nf90_err(nf90_put_att(ncid, ddsID, 'reldiff', reldiff(i)))
        call nf90_err(nf90_put_att(ncid, ddsID, 'henry', henry(i)))
        call nf90_err(nf90_put_att(ncid, ddsID, 'f0', f0(i)))
        call nf90_err(nf90_put_att(ncid, ddsID, 'dquer', dquer(i)))
        call nf90_err(nf90_put_att(ncid, ddsID, 'density', density(i)))
        call nf90_err(nf90_put_att(ncid, ddsID, 'dsigma', dsigma(i)))
        if (lnest) then
           ddspecIDn(i) = ddsID
        else
           ddspecID(i) = ddsID
        endif
     endif
  end do


  ! global (metadata) attributes
  !*******************************
  call writemetadata(ncid,lnest)


  ! moves the file from define to data mode
  call nf90_err(nf90_enddef(ncid))

  !  ! hes: inquire var definition
  !  do i = 1,nspec
  !     write(anspec,'(i3.3)') i
  !
  !     ! concentration output
  !     if (iout.eq.1.or.iout.eq.3.or.iout.eq.5) then
  !        if (lnest) then
  !           sID = specIDn(i)
  !        else
  !           sID = specID(i)
  !        endif
  !        call nf90_err(nf90_inquire_variable(ncid, sID, chunksizes=inq_chunksizes))
  !        write(*,*) "Chunksizes for var "//anspec//": ", inq_chunksizes
  !     endif
  !  end do

  
  ! fill with data
  !******************************
  ! longitudes (grid cell centers)
  if (lnest) then
    if (.not.allocated(coord)) allocate(coord(numxgridn))
     do i = 1,numxgridn
        coord(i) = outlon0n + (i-0.5)*dxoutn
     enddo
     call nf90_err(nf90_put_var(ncid, lonID, coord(1:numxgridn)))
     deallocate(coord)
  else
    if (.not.allocated(coord)) allocate(coord(numxgrid))
     do i = 1,numxgrid
        coord(i) = outlon0 + (i-0.5)*dxout
     enddo
     call nf90_err(nf90_put_var(ncid, lonID, coord(1:numxgrid)))
     deallocate(coord)
  endif
  ! latitudes (grid cell centers)
  if (lnest) then
    if (.not.allocated(coord)) allocate(coord(numygridn))
     do i = 1,numygridn
        coord(i) = outlat0n + (i-0.5)*dyoutn
     enddo
     call nf90_err(nf90_put_var(ncid, latID, coord(1:numygridn)))
     deallocate(coord)
  else
    if (.not.allocated(coord)) allocate(coord(numygrid))
     do i = 1,numygrid
        coord(i) = outlat0 + (i-0.5)*dyout
     enddo
     call nf90_err(nf90_put_var(ncid, latID, coord(1:numygrid)))
     deallocate(coord)
  endif
  ! levels
  call nf90_err(nf90_put_var(ncid, levID, outheight(1:numzgrid)))

  ! volume
  if (write_vol) then
    if (lnest) then
      call nf90_err(nf90_put_var(ncid, volID, volumen(:,:,:)))
    else
      call nf90_err(nf90_put_var(ncid, volID, volume(:,:,:)))
    end if
  end if

  ! area
  if (write_area) then
    if (lnest) then
      call nf90_err(nf90_put_var(ncid, areaID, arean(:,:)))
    else
      call nf90_err(nf90_put_var(ncid, areaID, area(:,:)))
    end if
  end if

  if (write_releases.eqv..true.) then
    ! release point information
    do i = 1,numpoint
       call nf90_err(nf90_put_var(ncid, relstartID, ireleasestart(i), (/i/)))
       call nf90_err(nf90_put_var(ncid, relendID, ireleaseend(i), (/i/)))
       call nf90_err(nf90_put_var(ncid, relkindzID, kindz(i), (/i/)))
       xp1=xpoint1(i)*dx+xlon0
       yp1=ypoint1(i)*dy+ylat0
       xp2=xpoint2(i)*dx+xlon0
       yp2=ypoint2(i)*dy+ylat0
       call nf90_err(nf90_put_var(ncid, rellng1ID, xp1, (/i/)))
       call nf90_err(nf90_put_var(ncid, rellng2ID, xp2, (/i/)))
       call nf90_err(nf90_put_var(ncid, rellat1ID, yp1, (/i/)))
       call nf90_err(nf90_put_var(ncid, rellat2ID, yp2, (/i/)))
       call nf90_err(nf90_put_var(ncid, relzz1ID, zpoint1(i), (/i/)))
       call nf90_err(nf90_put_var(ncid, relzz2ID, zpoint2(i), (/i/)))
       call nf90_err(nf90_put_var(ncid, relpartID, npart(i), (/i/)))
       if (i .le. 1000) then
         call nf90_err(nf90_put_var(ncid, relcomID, compoint(i), (/1,i/), (/45,1/)))
       else
         call nf90_err(nf90_put_var(ncid, relcomID, 'NA', (/1,i/), (/45,1/)))
       endif 
       call nf90_err(nf90_put_var(ncid, relxmassID, xmass(i,1:nspec), (/i,1/), (/1,nspec/)))
    end do
  end if

  ! age classes
  call nf90_err(nf90_put_var(ncid, lageID, lage(1:nageclass)))

  ! orography 
  if (.not. min_size) then
    if (lnest) then
      call nf90_err(nf90_put_var(ncid, oroID, orooutn(0:(nnx-1), 0:(nny-1))))
    else
      call nf90_err(nf90_put_var(ncid, oroID, oroout(0:(nnx-1), 0:(nny-1))))
    endif
  end if

  call nf90_err(nf90_close(ncid))

  return
end subroutine writeheader_netcdf

subroutine concoutput_netcdf(itime,outnum,gridtotalunc,wetgridtotalunc,drygridtotalunc)
     
  !                          i     i          o             o
  !       o
  !*****************************************************************************
  !                                                                            *
  !     Output of the concentration grid and the receptor concentrations.      *
  !                                                                            *
  !     Author: A. Stohl                                                       *
  !                                                                            *
  !     24 May 1995                                                            *
  !                                                                            *
  !     13 April 1999, Major update: if output size is smaller, dump output in *
  !                    sparse matrix format; additional output of uncertainty  *
  !                                                                            *
  !     05 April 2000, Major update: output of age classes; output for backward*
  !                    runs is time spent in grid cell times total mass of     *
  !                    species.                                                *
  !                                                                            *
  !     17 February 2002, Appropriate dimensions for backward and forward runs *
  !                       are now specified in module par_mod                  *
  !                                                                            *
  !     June 2006, write grid in sparse matrix with a single write command     *
  !                in order to save disk space                                 *
  !                                                                            *
  !     2008 new sparse matrix format                                          *
  !                                                                            *
  !     February 2010, Dominik Brunner, Empa                                   *
  !                    Adapted for COSMO                                       *
  !                    Remark: calculation of density could be improved.       *
  !                    Currently, it is calculated for the lower left corner   *
  !                    of each output grid cell rather than for its center.    *
  !                    Furthermore, the average density could be calculated    *
  !                    from the difference in pressure at the top and bottom   *
  !                    of each cell rather than by interpolation.              *
  !                                                                            *
  !     April 2013, Dominik Brunner, Empa                                      *
  !                    Adapted for netcdf output                               *
  !                                                                            *
  !*****************************************************************************
  !                                                                            *
  ! Variables:                                                                 *
  ! outnum          number of samples                                          *
  ! ncells          number of cells with non-zero concentrations               *
  ! sparse          .true. if in sparse matrix format, else .false.            *
  ! tot_mu          1 for forward, initial mass mixing ration for backw. runs  *
  !                                                                            *
  !*****************************************************************************

  use unc_mod, only: gridunc,drygridunc,wetgridunc,drygridunc0,wetgridunc0

  implicit none

  integer, intent(in) :: itime
  real, intent(in)    :: outnum
  real(dep_prec),intent(out):: wetgridtotalunc,drygridtotalunc
  real, intent(out)   :: gridtotalunc
  real                :: densityoutrecept(maxreceptor)
  integer             :: ncid,kp,ks,kz,ix,jy,iix,jjy,kzz,kzzm1,ngrid
  integer             :: nage,i,l,jj
  real                :: tot_mu(maxspec,maxpointspec_act)
  real                :: halfheight,dz,dz1,dz2
  real                :: xl,yl,xlrot,ylrot,zagnd,zagndprev
  real(dep_prec)      :: auxgrid(nclassunc)
  real(dep_prec)      :: gridtotal,gridsigmatotal
  real(dep_prec)      :: wetgridtotal,wetgridsigmatotal
  real(dep_prec)      :: drygridtotal,drygridsigmatotal
  ! real(sp)            :: gridtotal,gridsigmatotal
  ! real(sp)            :: wetgridtotal,wetgridsigmatotal
  ! real(sp)            :: drygridtotal,drygridsigmatotal

  real, parameter     :: weightair=28.97

  eps=nxmax/3.e5

  ! open output file
  call nf90_err(nf90_open(trim(ncfname), nf90_write, ncid))

  ! write time
  tpointer = tpointer + 1
  call nf90_err(nf90_put_var( ncid, timeID, itime, (/ tpointer /)))
  
  ! For forward simulations, output fields have dimension MAXSPEC,
  ! for backward simulations, output fields have dimension MAXPOINT.
  ! Thus, make loops either about nspec, or about numpoint
  !*****************************************************************

  if (ldirect.eq.1) then
    do ks=1,nspec
      do kp=1,maxpointspec_act
        tot_mu(ks,kp)=1.0
      end do
    end do
  else
    do ks=1,nspec
      do kp=1,maxpointspec_act
        tot_mu(ks,kp)=xmass(kp,ks)
      end do
    end do
  endif


  !*******************************************************************
  ! Compute air density:
  ! brd134: we now take into account whether we are in the mother or in
  !    a nested domain (before only from mother domain)
  ! Determine center altitude of output layer, and interpolate density
  ! data to that altitude
  !*******************************************************************

  do kz=1,numzgrid
    if (kz.eq.1) then
      halfheight=outheight(1)/2.
    else
      halfheight=(outheight(kz)+outheight(kz-1))/2.
    endif
    do kzz=2,nz
      if ((height(kzz-1).lt.halfheight).and. &
           (height(kzz).gt.halfheight)) exit
    end do
    kzz=max(min(kzz,nz),2)
    dz1=halfheight-height(kzz-1)
    dz2=height(kzz)-halfheight
    dz=dz1+dz2

    do jy=0,numygrid-1
      do ix=0,numxgrid-1
        xl=outlon0+real(ix)*dxout
        yl=outlat0+real(jy)*dyout
        ! grid index in mother domain
        xl=(xl-xlon0)/dx
        yl=(yl-ylat0)/dx

        ngrid=0
        do jj=numbnests,1,-1
          if ( xl.gt.xln(jj)+eps .and. xl.lt.xrn(jj)-eps .and. &
                 yl.gt.yln(jj)+eps .and. yl.lt.yrn(jj)-eps ) then
            ngrid=jj
            exit 
          end if
        end do

        if (ngrid.eq.0) then
          iix=max(min(nint(xl),nxmin1),0) ! if output grid cell is outside mother domain
          jjy=max(min(nint(yl),nymin1),0)

          densityoutgrid(ix,jy,kz)=(rho(iix,jjy,kzz,memind(2))*dz1+ &
             rho(iix,jjy,kzz-1,memind(2))*dz2)/dz
        else
          xl=(xl-xln(ngrid))*xresoln(ngrid)
          yl=(yl-yln(ngrid))*yresoln(ngrid)
          iix=max(min(nint(xl),nxn(ngrid)-1),0)
          jjy=max(min(nint(yl),nyn(ngrid)-1),0)

          densityoutgrid(ix,jy,kz)=(rhon(iix,jjy,kzz,memind(2), ngrid)*dz1+ &
             rhon(iix,jjy,kzz-1,memind(2), ngrid)*dz2)/dz
        endif
      end do
    end do
  end do

  ! brd134: for receptor points no option for nests yet to specify density
  !    and also altitude zreceptor not considered yet (needs revision)
  if (numreceptor.gt.0) then 
    do i=1,numreceptor
      xl=xreceptor(i)
      yl=yreceptor(i)
      iix=max(min(nint(xl),nxmin1),0)
      jjy=max(min(nint(yl),nymin1),0)
      densityoutrecept(i)=rho(iix,jjy,1,memind(2))
    end do
  endif

  ! Output is different for forward and backward simulations
  if (ldirect.eq.1) then
     do kz=1,numzgrid
        do jy=0,numygrid-1
           do ix=0,numxgrid-1
              factor3d(ix,jy,kz)=1.e12/volume(ix,jy,kz)/outnum
           end do
        end do
     end do
  else
     do kz=1,numzgrid
        do jy=0,numygrid-1
           do ix=0,numxgrid-1
              factor3d(ix,jy,kz)=real(abs(loutaver))/outnum
           end do
        end do
     end do
  endif

  !*********************************************************************
  ! Determine the standard deviation of the mean concentration or mixing
  ! ratio (uncertainty of the output) and the dry and wet deposition
  !*********************************************************************

  gridtotal=0.
  gridsigmatotal=0.
  gridtotalunc=0.
  wetgridtotal=0._dep_prec
  wetgridsigmatotal=0._dep_prec
  wetgridtotalunc=0._dep_prec
  drygridtotal=0._dep_prec
  drygridsigmatotal=0._dep_prec
  drygridtotalunc=0._dep_prec

  do ks=1,nspec

    do kp=1,maxpointspec_act
      do nage=1,nageclass

        do jy=0,numygrid-1
          do ix=0,numxgrid-1

            ! WET DEPOSITION
            if ((wetdep).and.(ldirect.gt.0)) then
              if (mpi_mode.gt.0) then
                do l=1,nclassunc
                  auxgrid(l)=wetgridunc0(ix,jy,ks,kp,l,nage)
                end do
              else
                do l=1,nclassunc
                  auxgrid(l)=wetgridunc(ix,jy,ks,kp,l,nage)
                end do
              end if
              call mean(auxgrid,wetgrid(ix,jy), &
                   wetgridsigma(ix,jy),nclassunc)
              ! Multiply by number of classes to get total concentration
              wetgrid(ix,jy)=wetgrid(ix,jy)*real(nclassunc,kind=sp)
              wetgridtotal=wetgridtotal+wetgrid(ix,jy)
              ! Calculate standard deviation of the mean
              wetgridsigma(ix,jy)= &
                   wetgridsigma(ix,jy)* &
                   sqrt(real(nclassunc,kind=dep_prec))
              wetgridsigmatotal=wetgridsigmatotal+ &
                   wetgridsigma(ix,jy)
            endif

            ! DRY DEPOSITION
            if ((drydep).and.(ldirect.gt.0)) then
              if (mpi_mode.gt.0) then
                do l=1,nclassunc
                  auxgrid(l)=drygridunc0(ix,jy,ks,kp,l,nage)
                end do
              else
                do l=1,nclassunc
                  auxgrid(l)=drygridunc(ix,jy,ks,kp,l,nage)
                end do
              end if
              call mean(auxgrid,drygrid(ix,jy), &
                   drygridsigma(ix,jy),nclassunc)
              ! Multiply by number of classes to get total concentration
              drygrid(ix,jy)=drygrid(ix,jy)*real(nclassunc,kind=sp)
              drygridtotal=drygridtotal+drygrid(ix,jy)
              ! Calculate standard deviation of the mean
              drygridsigma(ix,jy)= &
                   drygridsigma(ix,jy)* &
                   sqrt(real(nclassunc, kind=dep_prec))
              drygridsigmatotal=drygridsigmatotal+ &
                   drygridsigma(ix,jy)
            endif

            ! CONCENTRATION OR MIXING RATIO
            do kz=1,numzgrid
              do l=1,nclassunc
                auxgrid(l)=gridunc(ix,jy,kz,ks,kp,l,nage)
              end do
              call mean(auxgrid,grid(ix,jy,kz), &
                   gridsigma(ix,jy,kz),nclassunc)
              ! Multiply by number of classes to get total concentration
              grid(ix,jy,kz)= &
                   grid(ix,jy,kz)*real(nclassunc)
              gridtotal=gridtotal+grid(ix,jy,kz)
              ! Calculate standard deviation of the mean
              gridsigma(ix,jy,kz)= &
                   gridsigma(ix,jy,kz)* &
                   sqrt(real(nclassunc))
              gridsigmatotal=gridsigmatotal+ &
                   gridsigma(ix,jy,kz)
            end do
          end do
        end do

  !       print*,gridtotal,maxpointspec_act

        !*******************************************************************
        ! Generate output: may be in concentration (ng/m3) or in mixing
        ! ratio (ppt) or both
        ! Output the position and the values alternated multiplied by
        ! 1 or -1, first line is number of values, number of positions
        ! For backward simulations, the unit is seconds, stored in grid_time
        !*******************************************************************

        ! Concentration output
        !*********************
        if ((iout.eq.1).or.(iout.eq.3).or.(iout.eq.5)) then

          ! Wet deposition
          if ((ldirect.eq.1).and.(WETDEP)) then
            call nf90_err(nf90_put_var(ncid,wdspecID(ks),1.e12*&
                 wetgrid(0:numxgrid-1,0:numygrid-1)/area(0:numxgrid-1,0:numygrid-1),&
                 (/ 1,1,tpointer,kp,nage /), (/ numxgrid,numygrid,1,1,1 /)))
          end if

          ! Dry deposition
          if ((ldirect.eq.1).and.(DRYDEP)) then
            call nf90_err(nf90_put_var(ncid,ddspecID(ks),1.e12*&
                 drygrid(0:numxgrid-1,0:numygrid-1)/area(0:numxgrid-1,0:numygrid-1),&
                 (/ 1,1,tpointer,kp,nage /), (/ numxgrid,numygrid,1,1,1 /)))
          endif

          ! Concentrations
          call nf90_err(nf90_put_var(ncid,specID(ks),grid(0:numxgrid-1,0:numygrid-1,&
             1:numzgrid)*factor3d(0:numxgrid-1,0:numygrid-1,1:numzgrid)/tot_mu(ks,kp),&
               (/ 1,1,1,tpointer,kp,nage /), (/ numxgrid,numygrid,numzgrid,1,1,1 /) ))
 
        endif !  concentration output

        ! Mixing ratio output
        !********************

        if ((iout.eq.2).or.(iout.eq.3)) then      ! mixing ratio

          ! Wet deposition
          if ((ldirect.eq.1).and.(WETDEP)) then
            call nf90_err(nf90_put_var(ncid,wdspecID(ks),1.e12*&
                 wetgrid(0:numxgrid-1,0:numygrid-1)/area(0:numxgrid-1,0:numygrid-1),&
                 (/ 1,1,tpointer,kp,nage /), (/ numxgrid,numygrid,1,1,1 /)))

          endif

          ! Dry deposition
          if ((ldirect.eq.1).and.(DRYDEP)) then
            call nf90_err(nf90_put_var(ncid,ddspecID(ks),1.e12*&
                 drygrid(0:numxgrid-1,0:numygrid-1)/area(0:numxgrid-1,0:numygrid-1),&
                 (/ 1,1,tpointer,kp,nage /), (/ numxgrid,numygrid,1,1,1 /)))
          endif

          ! Mixing ratios
          call nf90_err(nf90_put_var(ncid,specIDppt(ks),weightair/weightmolar(ks)*&
               grid(0:numxgrid-1,0:numygrid-1,1:numzgrid)*&
               factor3d(0:numxgrid-1,0:numygrid-1,1:numzgrid)/&
               densityoutgrid(0:numxgrid-1,0:numygrid-1,1:numzgrid),&
               (/ 1,1,1,tpointer,kp,nage /), (/ numxgrid,numygrid,numzgrid,1,1,1 /)))

        endif ! output for ppt

      end do
    end do

  end do

  ! Close netCDF file
  !**************************
  call nf90_err(nf90_close(ncid))


  if (gridtotal.gt.0.) gridtotalunc=gridsigmatotal/gridtotal
  if (wetgridtotal.gt.0.) wetgridtotalunc=wetgridsigmatotal/ &
       wetgridtotal
  if (drygridtotal.gt.0.) drygridtotalunc=real(drygridsigmatotal/ &
       drygridtotal, kind=dep_prec)

  ! Dump of receptor concentrations

  if (numreceptor.gt.0 .and. (iout.eq.2 .or. iout.eq.3)  ) then
    write(unitoutreceptppt) itime
    do ks=1,nspec
      write(unitoutreceptppt) (1.e12*creceptor(i,ks)/outnum* &
           weightair/weightmolar(ks)/densityoutrecept(i),i=1,numreceptor)
    end do
  endif

  ! Dump of receptor concentrations

  if (numreceptor.gt.0) then
    write(unitoutrecept) itime
    do ks=1,nspec
      write(unitoutrecept) (1.e12*creceptor(i,ks)/outnum,i=1,numreceptor)
    end do
  endif


  ! Reinitialization of grid
  !*************************
  if (numreceptor.gt.0) creceptor(1:numreceptor,1:nspec) = 0.
  gridunc(:,:,:,1:nspec,:,:,1:nageclass) = 0.  
end subroutine concoutput_netcdf

subroutine concoutput_surf_netcdf(itime,outnum,gridtotalunc,wetgridtotalunc,drygridtotalunc)

  use unc_mod, only: gridunc,drygridunc,wetgridunc,drygridunc0,wetgridunc0

  implicit none

  integer, intent(in) :: itime
  real, intent(in)    :: outnum
  real(sp), intent(out)   :: gridtotalunc
  real(dep_prec), intent(out)   :: wetgridtotalunc,drygridtotalunc

  print*,'Netcdf output for surface only not yet implemented'
end subroutine concoutput_surf_netcdf

subroutine concoutput_nest_netcdf(itime,outnum)
  !                               i     i 
  !*****************************************************************************
  !                                                                            *
  !     Output of the concentration grid and the receptor concentrations.      *
  !                                                                            *
  !     Author: A. Stohl                                                       *
  !                                                                            *
  !     24 May 1995                                                            *
  !                                                                            *
  !     13 April 1999, Major update: if output size is smaller, dump output in *
  !                    sparse matrix format; additional output of uncertainty  *
  !                                                                            *
  !     05 April 2000, Major update: output of age classes; output for backward*
  !                    runs is time spent in grid cell times total mass of     *
  !                    species.                                                *
  !                                                                            *
  !     17 February 2002, Appropriate dimensions for backward and forward runs *
  !                    are now specified in module par_mod                     *
  !                                                                            *
  !     June 2006, write grid in sparse matrix with a single write command     *
  !                    in order to save disk space                             *
  !                                                                            *
  !     2008 new sparse matrix format                                          *
  !                                                                            *
  !     19 February 2010, Dominik Brunner, Empa: Adapted for COSMO             *
  !                                                                            *
  !     April 2013, Dominik Brunner, Empa                                      *
  !                    Adapted for netcdf output                               *
  !                                                                            *
  !*****************************************************************************
  !                                                                            *
  ! Variables:                                                                 *
  ! itime           current simulation time                                    *
  ! outnum          number of samples                                          *
  !                                                                            *
  !*****************************************************************************

  use unc_mod, only: griduncn,drygriduncn,wetgriduncn,drygriduncn0,wetgriduncn0
 
  implicit none

  integer, intent(in) :: itime
  real, intent(in)    :: outnum
  real                :: densityoutrecept(maxreceptor)
  integer             :: ncid,kp,ks,kz,ix,jy,iix,jjy,kzz,kzzm1,ngrid
  integer             :: nage,i,l, jj
  real                :: tot_mu(maxspec,maxpointspec_act)
  real                :: halfheight,dz,dz1,dz2
  real                :: xl,yl,xlrot,ylrot,zagnd,zagndprev
  real(dep_prec)      :: auxgrid(nclassunc)
  real                :: gridtotal
  real, parameter     :: weightair=28.97

  eps=nxmax/3.e5

  ! open output file
  call nf90_err(nf90_open(trim(ncfnamen), nf90_write, ncid))

  ! write time (do not increase time counter here, done in main output domain)
  call nf90_err(nf90_put_var( ncid, timeID, itime, (/ tpointer /)))
  
  ! For forward simulations, output fields have dimension MAXSPEC,
  ! for backward simulations, output fields have dimension MAXPOINT.
  ! Thus, make loops either about nspec, or about numpoint
  !*****************************************************************

  if (ldirect.eq.1) then
    do ks=1,nspec
      do kp=1,maxpointspec_act
        tot_mu(ks,kp)=1.0
      end do
    end do
  else
    do ks=1,nspec
      do kp=1,maxpointspec_act
        tot_mu(ks,kp)=xmass(kp,ks)
      end do
    end do
  endif


  !*******************************************************************
  ! Compute air density:
  ! brd134: we now take into account whether we are in the mother or in
  !    a nested domain (before only from mother domain)
  ! Determine center altitude of output layer, and interpolate density
  ! data to that altitude
  !*******************************************************************

  do kz=1,numzgrid
    if (kz.eq.1) then
      halfheight=outheight(1)/2.
    else
      halfheight=(outheight(kz)+outheight(kz-1))/2.
    endif
    do kzz=2,nz
      if ((height(kzz-1).lt.halfheight).and. &
           (height(kzz).gt.halfheight)) exit
    end do
    kzz=max(min(kzz,nz),2)
    dz1=halfheight-height(kzz-1)
    dz2=height(kzz)-halfheight
    dz=dz1+dz2

    do jy=0,numygridn-1
      do ix=0,numxgridn-1
        xl=outlon0n+real(ix)*dxoutn
        yl=outlat0n+real(jy)*dyoutn
        xl=(xl-xlon0)/dx
        yl=(yl-ylat0)/dy

        ngrid=0
        do jj=numbnests,1,-1
          if ( xl.gt.xln(jj)+eps .and. xl.lt.xrn(jj)-eps .and. &
                 yl.gt.yln(jj)+eps .and. yl.lt.yrn(jj)-eps ) then
            ngrid=jj
            exit 
          end if
        end do

        if (ngrid.eq.0) then
          iix=max(min(nint(xl),nxmin1),0)
          jjy=max(min(nint(yl),nymin1),0)

          densityoutgrid(ix,jy,kz)=(rho(iix,jjy,kzz,memind(2))*dz1+ &
             rho(iix,jjy,kzz-1,memind(2))*dz2)/dz
        else
          xl=(xl-xln(ngrid))*xresoln(ngrid)
          yl=(yl-yln(ngrid))*yresoln(ngrid)
          iix=max(min(nint(xl),nxn(ngrid)-1),0)
          jjy=max(min(nint(yl),nyn(ngrid)-1),0)
          densityoutgrid(ix,jy,kz)=(rhon(iix,jjy,kzz,memind(2), ngrid)*dz1+ &
             rhon(iix,jjy,kzz-1,memind(2), ngrid)*dz2)/dz
        endif

      end do
    end do
  end do

  do i=1,numreceptor
    xl=xreceptor(i)
    yl=yreceptor(i)
    iix=max(min(nint(xl),nxmin1),0)
    jjy=max(min(nint(yl),nymin1),0)
    densityoutrecept(i)=rho(iix,jjy,1,memind(2))
  end do

  ! Output is different for forward and backward simulations
  if (ldirect.eq.1) then
     do kz=1,numzgrid
        do jy=0,numygridn-1
           do ix=0,numxgridn-1
              factor3d(ix,jy,kz)=1.e12/volumen(ix,jy,kz)/outnum
           end do
        end do
     end do
  else
     do kz=1,numzgrid
        do jy=0,numygridn-1
           do ix=0,numxgridn-1
              factor3d(ix,jy,kz)=real(abs(loutaver))/outnum
           end do
        end do
     end do
  endif

  !*********************************************************************
  ! Determine the standard deviation of the mean concentration or mixing
  ! ratio (uncertainty of the output) and the dry and wet deposition
  !*********************************************************************

  gridtotal=0.

  do ks=1,nspec

    do kp=1,maxpointspec_act
      do nage=1,nageclass

        do jy=0,numygridn-1
          do ix=0,numxgridn-1
            ! WET DEPOSITION
            if ((WETDEP).and.(ldirect.gt.0)) then
              if (mpi_mode.gt.0) then
                do l=1,nclassunc
                  auxgrid(l)=wetgriduncn0(ix,jy,ks,kp,l,nage)
                end do
              else
                do l=1,nclassunc
                  auxgrid(l)=wetgriduncn(ix,jy,ks,kp,l,nage)
                end do
              end if
              call mean(auxgrid,wetgrid(ix,jy), &
                   wetgridsigma(ix,jy),nclassunc)
              ! Multiply by number of classes to get total concentration
              wetgrid(ix,jy)=wetgrid(ix,jy)*real(nclassunc)
              ! Calculate standard deviation of the mean
              wetgridsigma(ix,jy)= &
                   wetgridsigma(ix,jy)* &
                   sqrt(real(nclassunc,kind=dep_prec))
            endif

            ! DRY DEPOSITION
            if ((DRYDEP).and.(ldirect.gt.0)) then
              if (mpi_mode.gt.0) then
                do l=1,nclassunc
                  auxgrid(l)=drygriduncn0(ix,jy,ks,kp,l,nage)
                end do
              else
                do l=1,nclassunc
                  auxgrid(l)=drygriduncn(ix,jy,ks,kp,l,nage)
                end do
              end if
              call mean(auxgrid,drygrid(ix,jy), &
                   drygridsigma(ix,jy),nclassunc)
              ! Multiply by number of classes to get total concentration
              drygrid(ix,jy)=drygrid(ix,jy)*real(nclassunc)
              ! Calculate standard deviation of the mean
              drygridsigma(ix,jy)= &
                   drygridsigma(ix,jy)* &
                   sqrt(real(nclassunc,kind=dep_prec))
            endif

            ! CONCENTRATION OR MIXING RATIO
            do kz=1,numzgrid
              do l=1,nclassunc
                auxgrid(l)=griduncn(ix,jy,kz,ks,kp,l,nage)
              end do
              call mean(auxgrid,grid(ix,jy,kz), &
                   gridsigma(ix,jy,kz),nclassunc)
              ! Multiply by number of classes to get total concentration
              grid(ix,jy,kz)= &
                   grid(ix,jy,kz)*real(nclassunc)
              gridtotal=gridtotal+grid(ix,jy,kz)
              ! Calculate standard deviation of the mean
              gridsigma(ix,jy,kz)= &
                   gridsigma(ix,jy,kz)* &
                   sqrt(real(nclassunc))
            end do
          end do
        end do

  !       print*,gridtotal,maxpointspec_act

        !*******************************************************************
        ! Generate output: may be in concentration (ng/m3) or in mixing
        ! ratio (ppt) or both
        ! Output the position and the values alternated multiplied by
        ! 1 or -1, first line is number of values, number of positions
        ! For backward simulations, the unit is seconds, stored in grid_time
        !*******************************************************************

        ! Concentration output
        !*********************
        if ((iout.eq.1).or.(iout.eq.3).or.(iout.eq.5)) then

          ! Wet deposition
          if ((ldirect.eq.1).and.(WETDEP)) then
             call nf90_err(nf90_put_var(ncid,wdspecIDn(ks),1.e12*&
                  wetgrid(0:numxgridn-1,0:numygridn-1)/arean(0:numxgridn-1,0:numygridn-1),&
                  (/ 1,1,tpointer,kp,nage /), (/ numxgridn,numygridn,1,1,1 /)))
          endif

          ! Dry deposition
          if ((ldirect.eq.1).and.(DRYDEP)) then
             call nf90_err(nf90_put_var(ncid,ddspecIDn(ks),1.e12*&
                  drygrid(0:numxgridn-1,0:numygridn-1)/arean(0:numxgridn-1,0:numygridn-1),&
                  (/ 1,1,tpointer,kp,nage /), (/ numxgridn,numygridn,1,1,1 /)))
          endif

          ! Concentrations
          call nf90_err(nf90_put_var(ncid,specIDn(ks),grid(0:numxgridn-1,0:numygridn-1,&
             1:numzgrid)*factor3d(0:numxgridn-1,0:numygridn-1,1:numzgrid)/tot_mu(ks,kp),&
               (/ 1,1,1,tpointer,kp,nage /), (/ numxgridn,numygridn,numzgrid,1,1,1 /)))
 
        endif !  concentration output

        ! Mixing ratio output
        !********************

        if ((iout.eq.2).or.(iout.eq.3)) then      ! mixing ratio

          ! Wet deposition
          if ((ldirect.eq.1).and.(WETDEP)) then
             call nf90_err(nf90_put_var(ncid,wdspecIDn(ks),1.e12*&
                  wetgrid(0:numxgridn-1,0:numygridn-1)/arean(0:numxgridn-1,0:numygridn-1),&
                  (/ 1,1,tpointer,kp,nage /), (/ numxgridn,numygridn,1,1,1 /)))
          endif

          ! Dry deposition
          if ((ldirect.eq.1).and.(DRYDEP)) then
             call nf90_err(nf90_put_var(ncid,ddspecIDn(ks),1.e12*&
                  drygrid(0:numxgridn-1,0:numygridn-1)/arean(0:numxgridn-1,0:numygridn-1),&
                  (/ 1,1,tpointer,kp,nage /), (/ numxgridn,numygridn,1,1,1 /)))
          endif

          ! Mixing ratios
          call nf90_err(nf90_put_var(ncid,specIDnppt(ks),weightair/weightmolar(ks)*&
               grid(0:numxgridn-1,0:numygridn-1,1:numzgrid)*&
               factor3d(0:numxgridn-1,0:numygridn-1,1:numzgrid)/&
               densityoutgrid(0:numxgridn-1,0:numygridn-1,1:numzgrid),&
               (/ 1,1,1,tpointer,kp,nage /), (/ numxgridn,numygridn,numzgrid,1,1,1 /)))

        endif ! output for ppt

      end do
    end do

  end do

  ! Close netCDF file
  !**************************
  call nf90_err(nf90_close(ncid))

  ! Reinitialization of grid
  !*************************

  if (numreceptor.gt.0) creceptor(1:numreceptor,1:nspec) = 0.
  griduncn(:,:,:,1:nspec,:,:,1:nageclass) = 0.  
end subroutine concoutput_nest_netcdf

subroutine concoutput_surf_nest_netcdf(itime,outnum)

  implicit none

  integer, intent(in) :: itime
  real, intent(in)    :: outnum

  print*,'Netcdf output for surface only not yet implemented'
end subroutine concoutput_surf_nest_netcdf

subroutine create_particles_initialoutput(itime,idate,itime_start,idate_start)

  implicit none

  integer, intent(in) :: itime,idate,itime_start,idate_start
  ! integer, intent(in) :: irelease
  integer             :: cache_size,ncid,j,totpart
  integer             :: partDimID
  character(len=11)   :: fprefix
  character(len=3)    :: anspec,arelease
  character           :: adate*8,atime*6,adate_start*8,atime_start*6,timeunit*32
  character(len=255)  :: fname_partoutput
  real                :: fillval

  write(adate,'(i8.8)') idate
  write(atime,'(i6.6)') itime
  write(adate_start,'(i8.8)') idate_start
  write(atime_start,'(i6.6)') itime_start
  ! write(arelease, '(i3.3)') irelease
  fprefix = 'partinit_'!rel'//arelease//'_'

  fname_partoutput = path(2)(1:length(2))//trim(fprefix)//adate//atime//'.nc'
  !ncfname_part(irelease) = fname_partoutput
  ncfname_partinit = fname_partoutput

  call nf90_err(nf90_create(trim(fname_partoutput), cmode = nf90_hdf5, ncid = ncid))!, &
    ! cache_size = cache_size))

  ! create dimensions:
  !*************************
  
  ! particle
  partinitpointer=0
  call nf90_err(nf90_def_dim(ncid, 'particle', nf90_unlimited, partDimID))

  ! create variables
  !*************************

  ! particles
  call nf90_err(nf90_def_var(ncid, 'particle', nf90_int, (/ partDimID/), partIDi))
  call nf90_err(nf90_put_att(ncid, partIDi, 'long_name', 'particle index'))

  fillval = -1.
  ! time
  timeunit = 'seconds since '//adate_start(1:4)//'-'//adate_start(5:6)// &
     '-'//adate_start(7:8)//' '//atime_start(1:2)//':'//atime_start(3:4)

  call nf90_err(nf90_def_var(ncid, 'time', nf90_int, (/ partDimID /), tIDi))
  call nf90_err(nf90_def_var_deflate(ncid,tIDi,shuffle=0,deflate=1,deflate_level=1))
  call nf90_err(nf90_put_att(ncid, tIDi, 'long_name', 'time of release'))
  call nf90_err(nf90_put_att(ncid, tIDi, '_FillValue', -1))
  call nf90_err(nf90_put_att(ncid, tIDi, 'axis', 't'))
  call nf90_err(nf90_put_att(ncid, tIDi, 'units', timeunit))
  call nf90_err(nf90_put_att(ncid, tIDi, 'calendar', 'proleptic_gregorian'))
  call nf90_err(nf90_put_att(ncid, tIDi, 'standard_name', 'time'))
  call nf90_err(nf90_put_att(ncid, tIDi, 'description', 'time of release'))


  ! lon
  call nf90_err(nf90_def_var(ncid, 'longitude', nf90_float, (/ partDimID /), lonIDi))
  call nf90_err(nf90_def_var_deflate(ncid,lonIDi,shuffle=0,deflate=1,deflate_level=1))
  call nf90_err(nf90_put_att(ncid, lonIDi, 'long_name', 'longitude in degree east'))
  call nf90_err(nf90_put_att(ncid, lonIDi, '_FillValue', fillval))
  call nf90_err(nf90_put_att(ncid, lonIDi, 'axis', 'Lon'))
  call nf90_err(nf90_put_att(ncid, lonIDi, 'units', 'degrees_east'))
  call nf90_err(nf90_put_att(ncid, lonIDi, 'standard_name', 'longitude'))
  call nf90_err(nf90_put_att(ncid, lonIDi, 'description', 'longitude of particles'))

  ! lat
  call nf90_err(nf90_def_var(ncid, 'latitude', nf90_float, (/ partDimID /), latIDi))
  call nf90_err(nf90_def_var_deflate(ncid,latIDi,shuffle=0,deflate=1,deflate_level=1))
  call nf90_err(nf90_put_att(ncid, latIDi, 'long_name', 'latitude in degree north'))
  call nf90_err(nf90_put_att(ncid, latIDi, 'axis', 'Lat'))
  call nf90_err(nf90_put_att(ncid, latIDi, '_FillValue', fillval))
  call nf90_err(nf90_put_att(ncid, latIDi, 'units', 'degrees_north'))
  call nf90_err(nf90_put_att(ncid, latIDi, 'standard_name', 'latitude'))
  call nf90_err(nf90_put_att(ncid, latIDi, 'description', 'latitude of particles'))

  ! height
  call nf90_err(nf90_def_var(ncid, 'height', nf90_float, (/ partDimID /), levIDi))
  call nf90_err(nf90_def_var_deflate(ncid,levIDi,shuffle=0,deflate=1,deflate_level=1))
  call nf90_err(nf90_put_att(ncid, levIDi, 'units', 'meters'))
  call nf90_err(nf90_put_att(ncid, levIDi, '_FillValue', fillval))
  call nf90_err(nf90_put_att(ncid, levIDi, 'positive', 'up'))
  call nf90_err(nf90_put_att(ncid, levIDi, 'standard_name', 'height'))
  call nf90_err(nf90_put_att(ncid, levIDi, 'long_name', 'height above ground'))

  ! moves the file from define to data mode
  call nf90_err(nf90_enddef(ncid))

  call nf90_err(nf90_close(ncid))
end subroutine create_particles_initialoutput

subroutine write_particles_initialoutput(itime,istart,iend)
   use particle_mod

   implicit none

   integer, intent(in) ::  &
      itime,               & ! time of particle release
      istart,              & ! index of first newly released particle
      iend                   ! index of last newly released partile
   integer, allocatable    :: partindices(:),releasetimes(:)
   integer :: newpart,ncid,j

   newpart = iend-istart
   call nf90_err(nf90_open(trim(ncfname_partinit), nf90_write, ncid))

   allocate ( partindices(newpart) )
  
   do j=1,newpart
      partindices(j)=j+partinitpointer
   end do 
   call nf90_err(nf90_put_var(ncid,partIDi,partindices,(/ partinitpointer+1 /),(/ newpart /)))
   deallocate (partindices)
   
   allocate ( releasetimes(newpart) )
   releasetimes=itime
   call nf90_err(nf90_put_var(ncid,tIDi,releasetimes,(/ partinitpointer+1 /),(/ newpart /)))
   deallocate (releasetimes)
   call nf90_err(nf90_put_var(ncid,lonIDi,xlon0+part(partinitpointer+1:iend)%xlon*dx, (/ partinitpointer+1 /),(/ newpart /)))
   call nf90_err(nf90_put_var(ncid,latIDi,ylat0+part(partinitpointer+1:iend)%ylat*dy, (/ partinitpointer+1 /),(/ newpart /)))
   call nf90_err(nf90_put_var(ncid,levIDi,part(partinitpointer+1:iend)%z, (/ partinitpointer+1 /),(/ newpart /)))

   call nf90_err(nf90_close(ncid))

   partinitpointer = partinitpointer+newpart
end subroutine write_particles_initialoutput

subroutine writeheader_partoutput(itime,idate,itime_start,idate_start)!,irelease)

  implicit none

  integer, intent(in) :: itime,idate,itime_start,idate_start
  ! integer, intent(in) :: irelease
  integer             :: cache_size,ncid,j,totpart
  integer             :: timeDimID,partDimID,tID,memDimID
  character(len=11)   :: fprefix
  character(len=3)    :: anspec,arelease
  character           :: adate*8,atime*6,adate_start*8,atime_start*6,timeunit*32
  character(len=255)  :: fname_partoutput
  real                :: fillval

  open(unit=unittmp,file=trim(path(2)(1:length(2)))//'test_dir.txt',status='replace',&
       &err=110)
  close (unittmp, status='delete')

  write(adate,'(i8.8)') idate
  write(atime,'(i6.6)') itime
  write(adate_start,'(i8.8)') idate_start
  write(atime_start,'(i6.6)') itime_start
  ! write(arelease, '(i3.3)') irelease
  fprefix = 'partoutput_'!rel'//arelease//'_'

  fname_partoutput = path(2)(1:length(2))//trim(fprefix)//adate//atime//'.nc'
  !ncfname_part(irelease) = fname_partoutput
  ncfname_part = fname_partoutput

  totpart=0
  do j=1,numpoint
    totpart = totpart+npart(j)
  end do
  !totpart = maxpart!max(numpart,totpart)
  !cache_size = 4 * 1 * (12+nspec)

  write(*,*) nspec,numpart,totpart

  call nf90_err(nf90_create(trim(fname_partoutput), cmode = nf90_hdf5, ncid = ncid))!, &
    ! cache_size = cache_size))

  ! create dimensions:
  !*************************
  ! time
  call nf90_err(nf90_def_dim(ncid, 'time', nf90_unlimited, timeDimID))
  timeunit = 'seconds since '//adate_start(1:4)//'-'//adate_start(5:6)// &
     '-'//adate_start(7:8)//' '//atime_start(1:2)//':'//atime_start(3:4)

  ! particle
  call nf90_err(nf90_def_dim(ncid, 'particle', nf90_unlimited, partDimID)) !totpart needs to be the actual number of particles


  ! create variables
  !*************************

  ! time
  tpointer_part=0
  call nf90_err(nf90_def_var(ncid, 'time', nf90_int, (/ timeDimID /), tID))
  call nf90_err(nf90_put_att(ncid, tID, 'units', timeunit))
  call nf90_err(nf90_put_att(ncid, tID, 'calendar', 'proleptic_gregorian'))

  timeIDpart=tID
  ! particles
  call nf90_err(nf90_def_var(ncid, 'particle', nf90_int, (/ partDimID/), partID))
  call nf90_err(nf90_put_att(ncid, partID, 'long_name', 'particle index'))

  fillval = -1.
  ! lon
  call nf90_err(nf90_def_var(ncid, 'longitude', nf90_float, (/ timeDimID,partDimID /), lonID))
  call nf90_err(nf90_def_var_chunking(ncid,lonID,NF90_CHUNKED,chunksizes=(/ 1,totpart /)))
  call nf90_err(nf90_def_var_deflate(ncid,lonID,shuffle=0,deflate=1,deflate_level=1))
  call nf90_err(nf90_put_att(ncid, lonID, 'long_name', 'longitude in degree east'))
  call nf90_err(nf90_put_att(ncid, lonID, '_FillValue', fillval))
  call nf90_err(nf90_put_att(ncid, lonID, 'axis', 'Lon'))
  call nf90_err(nf90_put_att(ncid, lonID, 'units', 'degrees_east'))
  call nf90_err(nf90_put_att(ncid, lonID, 'standard_name', 'longitude'))
  call nf90_err(nf90_put_att(ncid, lonID, 'description', 'longitude of particles'))

  ! lat
  call nf90_err(nf90_def_var(ncid, 'latitude', nf90_float, (/ timeDimID,partDimID /), latID))
  call nf90_err(nf90_def_var_chunking(ncid,latID,NF90_CHUNKED,chunksizes=(/ 1,totpart /)))
  call nf90_err(nf90_def_var_deflate(ncid,latID,shuffle=0,deflate=1,deflate_level=1))
  call nf90_err(nf90_put_att(ncid, latID, 'long_name', 'latitude in degree north'))
  call nf90_err(nf90_put_att(ncid, latID, 'axis', 'Lat'))
  call nf90_err(nf90_put_att(ncid, latID, '_FillValue', fillval))
  call nf90_err(nf90_put_att(ncid, latID, 'units', 'degrees_north'))
  call nf90_err(nf90_put_att(ncid, latID, 'standard_name', 'latitude'))
  call nf90_err(nf90_put_att(ncid, latID, 'description', 'latitude of particles'))

  ! height
  call nf90_err(nf90_def_var(ncid, 'height', nf90_float, (/ timeDimID,partDimID /), levID))
  call nf90_err(nf90_def_var_chunking(ncid,levID,NF90_CHUNKED,chunksizes=(/ 1,totpart /)))
  call nf90_err(nf90_def_var_deflate(ncid,levID,shuffle=0,deflate=1,deflate_level=1))
  call nf90_err(nf90_put_att(ncid, levID, 'units', 'meters'))
  call nf90_err(nf90_put_att(ncid, levID, '_FillValue', fillval))
  call nf90_err(nf90_put_att(ncid, levID, 'positive', 'up'))
  call nf90_err(nf90_put_att(ncid, levID, 'standard_name', 'height'))
  call nf90_err(nf90_put_att(ncid, levID, 'long_name', 'height above ground'))

  ! itramem
  ! call nf90_err(nf90_def_var(ncid, 'itramem', nf90_float, (/ timeDimID,partDimID /), itramemID))

  ! topo
  call nf90_err(nf90_def_var(ncid, 'topo', nf90_float, (/ timeDimID,partDimID /), topoID))
  call nf90_err(nf90_def_var_chunking(ncid,topoID,NF90_CHUNKED,chunksizes=(/ 1,totpart /)))
  call nf90_err(nf90_def_var_deflate(ncid,topoID,shuffle=0,deflate=1,deflate_level=1))
  call nf90_err(nf90_put_att(ncid, topoID, 'units', 'meters'))
  call nf90_err(nf90_put_att(ncid, topoID, '_FillValue', fillval))
  call nf90_err(nf90_put_att(ncid, topoID, 'standard_name', 'topography'))
  call nf90_err(nf90_put_att(ncid, topoID, 'long_name', 'topography above sealevel'))

  ! Potential vorticity
  call nf90_err(nf90_def_var(ncid, 'pv', nf90_float, (/ timeDimID,partDimID /), pvID))
  call nf90_err(nf90_def_var_chunking(ncid,pvID,NF90_CHUNKED,chunksizes=(/ 1,totpart /)))
  call nf90_err(nf90_def_var_deflate(ncid,pvID,shuffle=0,deflate=1,deflate_level=1))
  call nf90_err(nf90_put_att(ncid, pvID, 'units', '?'))
  call nf90_err(nf90_put_att(ncid, pvID, '_FillValue', fillval))
  call nf90_err(nf90_put_att(ncid, pvID, 'standard_name', 'potential_vorticity'))
  call nf90_err(nf90_put_att(ncid, pvID, 'long_name', 'potential vorticity'))
  
  ! Specific humidity
  call nf90_err(nf90_def_var(ncid, 'qv', nf90_float, (/ timeDimID,partDimID /), qvID))
  call nf90_err(nf90_def_var_chunking(ncid,qvID,NF90_CHUNKED,chunksizes=(/ 1,totpart /)))
  call nf90_err(nf90_def_var_deflate(ncid,qvID,shuffle=0,deflate=1,deflate_level=1))
  call nf90_err(nf90_put_att(ncid, qvID, 'units', '?'))
  call nf90_err(nf90_put_att(ncid, qvID, '_FillValue', fillval))
  call nf90_err(nf90_put_att(ncid, qvID, 'standard_name', 'specific_humidity'))
  call nf90_err(nf90_put_att(ncid, qvID, 'long_name', 'specific humidity'))

  ! Density
  call nf90_err(nf90_def_var(ncid, 'rho', nf90_float, (/ timeDimID,partDimID /), rhoID))
  call nf90_err(nf90_def_var_chunking(ncid,rhoID,NF90_CHUNKED,chunksizes=(/ 1,totpart /)))
  call nf90_err(nf90_def_var_deflate(ncid,rhoID,shuffle=0,deflate=1,deflate_level=1))
  call nf90_err(nf90_put_att(ncid, rhoID, 'units', 'kg/m3'))
  call nf90_err(nf90_put_att(ncid, rhoID, '_FillValue', fillval))
  call nf90_err(nf90_put_att(ncid, rhoID, 'positive', 'up'))
  call nf90_err(nf90_put_att(ncid, rhoID, 'standard_name', 'density'))
  call nf90_err(nf90_put_att(ncid, rhoID, 'long_name', 'density'))

  ! Mixing layer height
  call nf90_err(nf90_def_var(ncid, 'hmix', nf90_float, (/ timeDimID,partDimID /), hmixID))
  call nf90_err(nf90_def_var_chunking(ncid,hmixID,NF90_CHUNKED,chunksizes=(/ 1,totpart /)))
  call nf90_err(nf90_def_var_deflate(ncid,hmixID,shuffle=0,deflate=1,deflate_level=1))
  call nf90_err(nf90_put_att(ncid, hmixID, 'units', 'meters'))
  call nf90_err(nf90_put_att(ncid, hmixID, '_FillValue', fillval))
  call nf90_err(nf90_put_att(ncid, hmixID, 'positive', 'up'))
  call nf90_err(nf90_put_att(ncid, hmixID, 'standard_name', 'hmix'))
  call nf90_err(nf90_put_att(ncid, hmixID, 'long_name', 'height above ground of mixing layer'))

  ! tr
  call nf90_err(nf90_def_var(ncid, 'tr', nf90_float, (/ timeDimID,partDimID /), trID))
  call nf90_err(nf90_def_var_chunking(ncid,trID,NF90_CHUNKED,chunksizes=(/ 1,totpart /)))
  call nf90_err(nf90_def_var_deflate(ncid,trID,shuffle=0,deflate=1,deflate_level=1))
  call nf90_err(nf90_put_att(ncid, trID, 'units', 'meters'))
  call nf90_err(nf90_put_att(ncid, trID, '_FillValue', fillval))
  call nf90_err(nf90_put_att(ncid, trID, 'positive', 'up'))
  call nf90_err(nf90_put_att(ncid, trID, 'standard_name', 'htropo'))
  call nf90_err(nf90_put_att(ncid, trID, 'long_name', 'height above ground of tropopause'))

  ! Temperature
  call nf90_err(nf90_def_var(ncid, 'tt', nf90_float, (/ timeDimID,partDimID /), ttID))
  call nf90_err(nf90_def_var_chunking(ncid,ttID,NF90_CHUNKED,chunksizes=(/ 1,totpart /)))
  call nf90_err(nf90_def_var_deflate(ncid,ttID,shuffle=0,deflate=1,deflate_level=1))
  call nf90_err(nf90_put_att(ncid, ttID, 'units', 'K'))
  call nf90_err(nf90_put_att(ncid, ttID, '_FillValue', fillval))
  call nf90_err(nf90_put_att(ncid, ttID, 'positive', 'up'))
  call nf90_err(nf90_put_att(ncid, ttID, 'standard_name', 'temperature'))
  call nf90_err(nf90_put_att(ncid, ttID, 'long_name', 'temperature'))

  do j=1,nspec
    ! Masses
    write(anspec, '(i3.3)') j
    call nf90_err(nf90_def_var(ncid, 'mass'//anspec, nf90_float, (/ timeDimID,partDimID /), massID(j)))
    call nf90_err(nf90_def_var_chunking(ncid,massID(j),NF90_CHUNKED,chunksizes=(/ 1,totpart /)))
    call nf90_err(nf90_def_var_deflate(ncid,massID(j),shuffle=0,deflate=1,deflate_level=1))
    call nf90_err(nf90_put_att(ncid, massID(j), 'units', 'kg'))
    call nf90_err(nf90_put_att(ncid, massID(j), '_FillValue', fillval))
    call nf90_err(nf90_put_att(ncid, massID(j), 'positive', 'up'))
    call nf90_err(nf90_put_att(ncid, massID(j), 'standard_name', 'mass'//anspec))
    call nf90_err(nf90_put_att(ncid, massID(j), 'long_name', 'mass for nspec'//anspec))    
  end do

  ! global (metadata) attributes
  !*******************************
  call writemetadata(ncid,lnest=.false.)

  ! moves the file from define to data mode
  call nf90_err(nf90_enddef(ncid))

  call nf90_err(nf90_close(ncid))

  return
110 write(*,FMT='(80("#"))') 
  write(*,*) 'ERROR: output directory ', trim(path(2)(1:length(2))), ' does not exist&
       & (or failed to write there).' 
  write(*,*) 'EXITING' 
  write(*,FMT='(80("#"))')
  stop
end subroutine writeheader_partoutput

subroutine open_partoutput_file(ncid)!,irelease)
  
  implicit none 

  integer, intent(inout)         :: ncid
  !integer, intent(in)            :: irelease

  call nf90_err(nf90_open(trim(ncfname_part), nf90_write, ncid))
end subroutine open_partoutput_file

subroutine close_partoutput_file(ncid)
  
  implicit none 

  integer                        :: ncid

  call nf90_err(nf90_close(ncid))
end subroutine close_partoutput_file

subroutine partoutput_netcdf(itime,field,fieldname,imass,ncid)
  
  implicit none

  integer, intent(in)            :: itime,imass
  real, intent(in)               :: field(numpart)
  character(2), intent(in)       :: fieldname  ! input field to interpolate over
  integer, allocatable           :: partindices(:)
  integer                        :: ncid,newpart,j
  ! ! open output file
  ! call nf90_err(nf90_open(trim(ncfname_part), nf90_write, ncid))
  select case(fieldname)
    case('TI')
      ! write time
      tpointer_part = tpointer_part + 1
      call nf90_err(nf90_put_var(ncid, timeIDpart, itime, (/ tpointer_part /)))
    case('PA')
        newpart = numpart - ppointer_part
        
        if (tpointer_part.eq.1) then 
          allocate ( partindices(numpart) )
          do j=1,numpart 
            partindices(j)=j+ppointer_part
          end do 

          call nf90_err(nf90_put_var(ncid, partID,partindices, (/ 1 /),(/ numpart /)))

          deallocate (partindices)

          ppointer_part = numpart 

        else if (newpart.ge.0) then

          allocate ( partindices(newpart) )
          do j=1,newpart
            partindices(j)=j+ppointer_part
          end do

          call nf90_err(nf90_put_var(ncid, partID,partindices, (/ ppointer_part+1 /),(/ newpart /)))
          
          deallocate (partindices)

          ppointer_part = numpart
        endif 
    case('LO')
      call nf90_err(nf90_put_var(ncid,lonID,field, (/ tpointer_part,1 /),(/ 1,numpart /)))
    case('LA')
      call nf90_err(nf90_put_var(ncid,latID,field, (/ tpointer_part,1 /),(/ 1,numpart /)))
    case('ZZ')
      call nf90_err(nf90_put_var(ncid,levID,field, (/ tpointer_part,1 /),(/ 1,numpart /)))
    case('IT')
      call nf90_err(nf90_put_var(ncid,itramemID,field, (/ tpointer_part,1 /),(/ 1,numpart /)))
    case('TO')
      call nf90_err(nf90_put_var(ncid,topoID,field, (/ tpointer_part,1 /),(/ 1,numpart /)))
    case('PV')
      call nf90_err(nf90_put_var(ncid,pvID,field, (/ tpointer_part,1 /),(/ 1,numpart /)))
    case('QV')
      call nf90_err(nf90_put_var(ncid,qvID,field, (/ tpointer_part,1 /),(/ 1,numpart /)))
    case('RH')
      call nf90_err(nf90_put_var(ncid,rhoID,field, (/ tpointer_part,1 /),(/ 1,numpart /)))
    case('HM')
      call nf90_err(nf90_put_var(ncid,hmixID,field, (/ tpointer_part,1 /),(/ 1,numpart /)))
    case('TR')
      call nf90_err(nf90_put_var(ncid,trID,field, (/ tpointer_part,1 /),(/ 1,numpart /)))
    case('TT')
      call nf90_err(nf90_put_var(ncid,ttID,field, (/ tpointer_part,1 /),(/ 1,numpart /)))
    case('MA')
      call nf90_err(nf90_put_var(ncid,massID(imass),field, (/ tpointer_part,1 /),(/ 1,numpart /)))
  end select

  ! call nf90_err(nf90_close(ncid))
end subroutine partoutput_netcdf

subroutine readpartpositions_netcdf(ibtime,ibdate)
  use random_mod
  use particle_mod
  use date_mod

  implicit none 

  integer, intent(in) :: ibtime,ibdate
  integer             :: ncidend,tIDend,pIDend,tempIDend
  integer             :: tlen,plen,tend,i
  integer             :: idate_start,itime_start
  character           :: adate*8,atime*6,timeunit*32,adate_start*8,atime_start*6
  real(kind=dp)       :: julin,julcommand,julpartin

  integer :: idummy = -8

  write(adate,'(i8.8)') ibdate
  write(atime,'(i6.6)') ibtime
  
  if (mquasilag.ne.0) then 
    write(*,*) 'Combination of ipin, netcdf partoutput, and mquasilag!=0 does not work yet'
    stop 
  endif

  ! Open partoutput_end.nc file
  call nf90_err(nf90_open(trim('partoutput_end.nc'), mode=NF90_NOWRITE,ncid=ncidend))

  ! Take the positions of the particles at the last timestep in the file
  ! It needs to be the same as given in the COMMAND file, this is arbitrary
  ! and should be removed in the future for easier use

  ! First get the time dimension
  call nf90_err(nf90_inq_dimid(ncid=ncidend,name='time',dimid=tIDend))
  call nf90_err(nf90_inquire_dimension(ncid=ncidend,dimid=tIDend,len=tlen))

  ! Check if the time corresponds to the one given in the COMMAND file
  call nf90_err(nf90_inq_varid(ncid=ncidend,name='time',varid=tIDend))
  call nf90_err(nf90_get_att(ncid=ncidend,varid=tIDend,name='units',values=timeunit))
  call nf90_err(nf90_get_var(ncid=ncidend,varid=tIDend,values=tend,start=(/ tlen /)))!,count=(/ 1 /)))
  adate_start(1:4) = timeunit(15:18)
  adate_start(5:6) = timeunit(20:21)
  adate_start(7:8) = timeunit(23:24)
  atime_start = '000000'
  atime_start(1:2) = timeunit(26:27)
  atime_start(3:4) = timeunit(29:30)
  read(adate_start,*) idate_start
  read(atime_start,*) itime_start
  julin = juldate(idate_start,itime_start)+real(tend,kind=dp)/86400._dp
  julcommand = juldate(ibdate,ibtime)
  if (abs(julin-julcommand).gt.1.e-5) then 
    write(*,*) 'ERROR: The given starting time and date do not correspond to'
    write(*,*) 'the last timestep of partoutput_end.nc:'
    write(*,*) julin,julcommand,tend
    stop 
  endif

  ! Then the particle dimension
  call nf90_err(nf90_inq_dimid(ncid=ncidend,name='particle',dimid=pIDend))
  call nf90_err(nf90_inquire_dimension(ncid=ncidend,dimid=pIDend,len=plen))

  ! Now spawn the correct number of particles
  write(*,*) 'Npart:',plen
  call spawn_particles(0,plen)

  ! And give them the correct positions
  ! Longitude
  call nf90_err(nf90_inq_varid(ncid=ncidend,name='longitude',varid=tempIDend))
  call nf90_err(nf90_get_var(ncid=ncidend,varid=tempIDend,values=part(:)%xlon, & 
    start=(/ tlen, 1 /),count=(/ 1, plen /)))
  part(:)%xlon=(part(:)%xlon-xlon0)/dx
  ! Latitude
  call nf90_err(nf90_inq_varid(ncid=ncidend,name='latitude',varid=tempIDend))
  call nf90_err(nf90_get_var(ncid=ncidend,varid=tempIDend,values=part(:)%ylat, & 
    start=(/ tlen, 1 /),count=(/ 1, plen /)))
  part(:)%ylat=(part(:)%ylat-ylat0)/dx
  ! Height
  call nf90_err(nf90_inq_varid(ncid=ncidend,name='height',varid=tempIDend))
  call nf90_err(nf90_get_var(ncid=ncidend,varid=tempIDend,values=part(:)%z, & 
    start=(/ tlen, 1 /),count=(/ 1, plen /)))

  do i=1,plen
    part(i)%nclass=min(int(ran1(idummy)*real(nclassunc))+1, &
         nclassunc)    
    part(i)%idt=mintime
    part(i)%npoint=i
  end do

  call nf90_err(nf90_close(ncidend))
end subroutine readpartpositions_netcdf

end module netcdf_output_mod
