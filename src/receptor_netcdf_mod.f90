! SPDX-FileCopyrightText: FLEXPART 1998-2019, see flexpart_license.txt
! SPDX-License-Identifier: GPL-3.0-or-later

module receptor_netcdf_mod

  !*****************************************************************************
  !                                                                            *
  !    This module contains variables and subroutines for reading input        *
  !    for satellite receptors and for saving the output of all receptor       *
  !    types to netcdf file                                                    *
  !                                                                            *
  !    Author, Rona Thompson                                                   *
  !                                                                            *
  !*****************************************************************************

  use netcdf
  use par_mod
  use com_mod
  use point_mod
  use date_mod
  use netcdf_output_mod, only: nf90_err
  use windfields_mod,    only: prs, height, nzmax, nz

  implicit none

  ! general receptors
  integer :: rpointer
  integer :: nc_id
  integer :: recdim_id, nchardim_id
  integer :: timevar_id, recvar_id, reclonvar_id, reclatvar_id, recaltvar_id, recnamevar_id
  integer, dimension(:), allocatable :: concvar_id, uncvar_id
  integer :: nnvar_id, xkvar_id 
  ! satellite receptors
  integer :: spointer
  character(len=8), dimension(:), allocatable :: sat_name 
  integer :: ncsat_id
  integer :: satrecdim_id, sataltdim_id, satlayerdim_id
  integer :: sattimevar_id, satrecvar_id, satlonvar_id, satlatvar_id, sataltvar_id, satnamevar_id
  integer, dimension(:), allocatable :: satvar_id, satuncvar_id
  integer :: satnnvar_id, satxkvar_id

  contains

  subroutine receptorout_init

  !*****************************************************************************
  !                                                                            *
  !    Intitialize netcdf files for receptor concentrations                    *
  !                                                                            *
  !    Author: R. Thompson, Sep-2023                                           *
  !                                                                            *
  !                                                                            *
  !*****************************************************************************

    implicit none

    character(len=256) :: fn, timeunit
    character(len=8)   :: adate
    character(len=6)   :: atime
    integer            :: ks, n, ks_start
    character(8)       :: date
    character(10)      :: time
    character(5)       :: zone

    if (numreceptor.eq.0) then
      return
    endif

    if (llcmoutput) then
      ks_start=2
    else
      ks_start=1
    endif

    ! time at simulation start
    if (ldirect.eq.1) then
      write(adate,'(i8.8)') ibdate
      write(atime,'(i6.6)') ibtime
    else
      write(adate,'(i8.8)') iedate
      write(atime,'(i6.6)') ietime
    end if

    if ((iout.eq.1).or.(iout.eq.3).or.(iout.eq.5)) then
      ! concentration output file
      write(fn,'(A)') path(2)(1:length(2))//'receptor_conc.nc'
    else if (iout.eq.2) then
      ! mixing ratio output file
      ! note: for iout=3 (both conc and mixing ratio) only output conc at receptors
      write(fn,'(A)') path(2)(1:length(2))//'receptor_pptv.nc'
    endif

    ! initialization
    if (.not.allocated(concvar_id)) then
      allocate(concvar_id(nspec))
      allocate(uncvar_id(nspec))
    endif

    ! pointer for receptor dimension
    rpointer=1

    ! create file
    call nf90_err( nf90_create(trim(fn), nf90_clobber, nc_id) )

    ! define dimensions
    !******************

    call nf90_err( nf90_def_dim(nc_id, "rec", numreceptor, recdim_id) )
    call nf90_err( nf90_def_dim(nc_id, 'nchar', 16, nchardim_id) )

    ! define variables
    !*****************

    ! time
    timeunit = 'seconds since '//adate(1:4)//'-'//adate(5:6)// &
                  '-'//adate(7:8)//' '//atime(1:2)//':'//atime(3:4)
    call nf90_err( nf90_def_var(nc_id, 'time', nf90_int, (/ recdim_id /), timevar_id) )
    call nf90_err( nf90_put_att(nc_id, timevar_id, 'units', trim(timeunit)) )
    call nf90_err( nf90_put_att(nc_id, timevar_id, 'calendar', 'proleptic_gregorian') )

    ! receptors 
    call nf90_err( nf90_def_var(nc_id, "rec", nf90_float, (/ recdim_id /), recvar_id) )
    call nf90_err( nf90_put_att(nc_id, recvar_id, "longname", "receptors") )
    call nf90_err( nf90_put_att(nc_id, recvar_id, "units", "index") )

    ! receptor names 
    call nf90_err( nf90_def_var(nc_id, "receptorname", nf90_char, (/ nchardim_id, recdim_id /), &
                    recnamevar_id) )
    call nf90_err( nf90_put_att(nc_id, recnamevar_id, "longname", "receptor name") )

    ! receptor longitude 
    call nf90_err( nf90_def_var(nc_id, "lon", nf90_float, (/ recdim_id /), reclonvar_id) )
    call nf90_err( nf90_put_att(nc_id, reclonvar_id, "longname", "receptor longitude") )
    call nf90_err( nf90_put_att(nc_id, reclonvar_id, "units", "degrees_east") )

    ! receptor latitude
    call nf90_err( nf90_def_var(nc_id, "lat", nf90_float, (/ recdim_id /), reclatvar_id) )
    call nf90_err( nf90_put_att(nc_id, reclatvar_id, "longname", "receptor latitude") )
    call nf90_err( nf90_put_att(nc_id, reclatvar_id, "units", "degrees_north") )

    ! receptor altitude
    call nf90_err( nf90_def_var(nc_id, "lev", nf90_float, (/ recdim_id /), recaltvar_id) )
    call nf90_err( nf90_put_att(nc_id, recaltvar_id, "longname", "receptor altitude") )
    call nf90_err( nf90_put_att(nc_id, recaltvar_id, "units", "meters") )

    ! species specific variables
    do ks=ks_start,nspec
      write(*,*) 'species = ',trim(species(ks))
      ! concentration/mixing ratio variables
      call nf90_err( nf90_def_var(nc_id, trim(species(ks)), nf90_float, (/ recdim_id /), concvar_id(ks)) )
      call nf90_err( nf90_def_var(nc_id, trim(species(ks))//"_uncert", nf90_float, (/ recdim_id /), &
                       uncvar_id(ks)) )
      write(*,*) 'concvar_id(ks) = ',concvar_id(ks)
      if ((iout.eq.1).or.(iout.eq.3).or.(iout.eq.5)) then
        call nf90_err( nf90_put_att(nc_id, concvar_id(ks), "units", "ng/m3") )
        call nf90_err( nf90_put_att(nc_id, concvar_id(ks), "longname", "mean receptor concentration") )
        call nf90_err( nf90_put_att(nc_id, uncvar_id(ks), "units", "ng/m3") )
        call nf90_err( nf90_put_att(nc_id, uncvar_id(ks), "longname", "uncertainty receptor concentration") )
      else if ((iout.eq.2)) then
        call nf90_err( nf90_put_att(nc_id, concvar_id(ks), "units", "pptv") )
        call nf90_err( nf90_put_att(nc_id, concvar_id(ks), "longname", "mean receptor VMR") )
        call nf90_err( nf90_put_att(nc_id, uncvar_id(ks), "units", "pptv") )
        call nf90_err( nf90_put_att(nc_id, uncvar_id(ks), "longname", "uncertainty receptor VMR") )
      endif
    end do

    ! not species specific variables
    !  number of particles in receptor output
    call nf90_err( nf90_def_var(nc_id,"npart", nf90_float, (/ recdim_id /), nnvar_id) )
    call nf90_err( nf90_put_att(nc_id, nnvar_id, "units", "counts") )
    call nf90_err( nf90_put_att(nc_id, nnvar_id, "longname","number of particles at receptor") )
    !  average kernel weight at receptor
    call nf90_err( nf90_def_var(nc_id,"kernel", nf90_float, (/ recdim_id /), xkvar_id) )
    call nf90_err( nf90_put_att(nc_id, xkvar_id, "units", "") )
    call nf90_err( nf90_put_att(nc_id, xkvar_id, "longname", "average kernel weight at receptor") )

    ! write global attributes
    !************************

    call date_and_time(date,time,zone)
    call nf90_err( nf90_put_att(nc_id, nf90_global, 'Conventions', 'CF-1.6') )
    call nf90_err( nf90_put_att(nc_id, nf90_global, 'title', 'FLEXPART receptor output') )
    call nf90_err( nf90_put_att(nc_id, nf90_global, 'source', trim(flexversion)//' model output') )
    call nf90_err( nf90_put_att(nc_id, nf90_global, 'history', date(1:4)//'-'//date(5:6)//&
                     '-'//date(7:8)//' '//time(1:2)//':'//time(3:4)//' '//zone ) )
    call nf90_err( nf90_put_att(nc_id, nf90_global, 'references', &
          'Stohl et al., Atmos. Chem. Phys., 2005, doi:10.5194/acp-5-2461-200;'//&
          'Henne et al., in Lagrangian Modeling of the Atmosphere, 2012, doi:10.1029/2012GM001247') )

    ! end definition
    call nf90_err( nf90_enddef(nc_id) )

    ! write dimension variables 
    !**************************

    ! receptor index
    call nf90_err( nf90_put_var(nc_id, recvar_id, (/(n,n=1,numreceptor)/)) )

    ! close file
    call nf90_err( nf90_close(nc_id) )

  end subroutine receptorout_init


  subroutine satelliteout_init

  !*****************************************************************************
  !                                                                            *
  !    Intitialize netcdf files for satellite concentrations                   *
  !                                                                            *
  !    Author: R. Thompson, Oct-2023                                           *
  !                                                                            *
  !                                                                            *
  !*****************************************************************************

    implicit none

    character(len=256) :: fn, timeunit
    character(len=8)   :: adate
    character(len=6)   :: atime
    integer            :: ks, n, ks_start
    character(8)       :: date
    character(10)      :: time
    character(5)       :: zone

    if (numsatreceptor.eq.0) then
      return
    endif

    if (llcmoutput) then
      ks_start=2
    else
      ks_start=1
    endif

    ! time at simulation start
    if (ldirect.eq.1) then
      write(adate,'(i8.8)') ibdate
      write(atime,'(i6.6)') ibtime
    else
      write(adate,'(i8.8)') iedate
      write(atime,'(i6.6)') ietime
    end if

    ! mixing ratio output for satellites 
    write(fn,'(A)') path(2)(1:length(2))//'satellite_pptv.nc'

   ! initialization
    if (.not.allocated(satvar_id)) then
      allocate(satvar_id(nspec))
      allocate(satuncvar_id(nspec))
    endif

    ! pointer for satellite receptor dimension
    spointer=1

    ! create file
    call nf90_err( nf90_create(trim(fn), nf90_clobber, ncsat_id) )

    ! define dimensions
    !******************

    call nf90_err( nf90_def_dim(ncsat_id, "rec", numsatreceptor, satrecdim_id) )
    call nf90_err( nf90_def_dim(ncsat_id, "layer", (nlayermax-1), satlayerdim_id) )
    call nf90_err( nf90_def_dim(ncsat_id, "level", nlayermax, sataltdim_id) )
    call nf90_err( nf90_def_dim(ncsat_id, 'nchar', 16, nchardim_id) )

    ! define variables
    !*****************

    ! Note: did not include variables for kernel parameters as currently use same
    ! at all sites and did not include windspeed, ageclasses, and release points

    ! time
    timeunit = 'seconds since '//adate(1:4)//'-'//adate(5:6)// &
                  '-'//adate(7:8)//' '//atime(1:2)//':'//atime(3:4)
    call nf90_err( nf90_def_var(ncsat_id, 'time', nf90_int, (/ satrecdim_id /), sattimevar_id) )
    call nf90_err( nf90_put_att(ncsat_id, sattimevar_id, 'units', trim(timeunit)) )
    call nf90_err( nf90_put_att(ncsat_id, sattimevar_id, 'calendar', 'proleptic_gregorian') )

    ! receptor names 
    call nf90_err( nf90_def_var(ncsat_id, "receptorname", nf90_char, (/ nchardim_id, satrecdim_id /), &
                    satnamevar_id) )
    call nf90_err( nf90_put_att(ncsat_id, satnamevar_id, "longname", "receptor name") )

    ! receptors 
    call nf90_err( nf90_def_var(ncsat_id, "rec", nf90_float, (/ satrecdim_id /), satrecvar_id) )
    call nf90_err( nf90_put_att(ncsat_id, satrecvar_id, "longname", "receptors") )
    call nf90_err( nf90_put_att(ncsat_id, satrecvar_id, "units", "index") )

    ! receptor longitude 
    call nf90_err( nf90_def_var(ncsat_id, "lon", nf90_float, (/ satrecdim_id /), satlonvar_id) )
    call nf90_err( nf90_put_att(ncsat_id, satlonvar_id, "longname", "receptor longitude") )
    call nf90_err( nf90_put_att(ncsat_id, satlonvar_id, "units", "degrees_east") )

    ! receptor latitude
    call nf90_err( nf90_def_var(ncsat_id, "lat", nf90_float, (/ satrecdim_id /), satlatvar_id) )
    call nf90_err( nf90_put_att(ncsat_id, satlatvar_id, "longname", "receptor latitude") )
    call nf90_err( nf90_put_att(ncsat_id, satlatvar_id, "units", "degrees_north") )

    ! receptor altitude
    call nf90_err( nf90_def_var(ncsat_id, "alt", nf90_float, (/ satrecdim_id, sataltdim_id /), sataltvar_id) )
    call nf90_err( nf90_put_att(ncsat_id, sataltvar_id, "longname", "receptor altitude of levels") )
    call nf90_err( nf90_put_att(ncsat_id, sataltvar_id, "units", "meters") )

    ! species specific variables
    do ks=ks_start,nspec
      ! mixing ratio output for each layer of retrieval
      call nf90_err( nf90_def_var(ncsat_id, trim(species(ks)), nf90_float, (/ satrecdim_id, satlayerdim_id /), &
                       satvar_id(ks)) )
      call nf90_err( nf90_put_att(ncsat_id, satvar_id(ks), "units", "pptv") )
      call nf90_err( nf90_put_att(ncsat_id, satvar_id(ks), "longname", "mean VMR") )
      ! uncertainty output for each layer of retrieval
      call nf90_err( nf90_def_var(ncsat_id, trim(species(ks))//"_uncert", nf90_float, (/ satrecdim_id, satlayerdim_id /), &
                       satuncvar_id(ks)) )
      call nf90_err( nf90_put_att(ncsat_id, satuncvar_id(ks), "units", "pptv") )
      call nf90_err( nf90_put_att(ncsat_id, satuncvar_id(ks), "longname", "uncertainty VMR") )
    end do

    ! not species specific variables
    ! number of particles in receptor output
    call nf90_err( nf90_def_var(ncsat_id,"npart", nf90_float, (/ satrecdim_id, satlayerdim_id /), satnnvar_id) )
    call nf90_err( nf90_put_att(ncsat_id, satnnvar_id, "units", "counts") )
    call nf90_err( nf90_put_att(ncsat_id, satnnvar_id, "longname","number of particles at receptor") )
    ! average kernel weight at receptor
    call nf90_err( nf90_def_var(ncsat_id,"kernel", nf90_float, (/ satrecdim_id, satlayerdim_id /), satxkvar_id) )
    call nf90_err( nf90_put_att(ncsat_id, satxkvar_id, "units", "") )
    call nf90_err( nf90_put_att(ncsat_id, satxkvar_id, "longname", "average kernel weight at receptor") )

    ! write global attributes
    !************************

    call date_and_time(date,time,zone)
    call nf90_err( nf90_put_att(ncsat_id, nf90_global, 'Conventions', 'CF-1.6') )
    call nf90_err( nf90_put_att(ncsat_id, nf90_global, 'title', 'FLEXPART receptor output') )
    call nf90_err( nf90_put_att(ncsat_id, nf90_global, 'source', trim(flexversion)//' model output') )
    call nf90_err( nf90_put_att(ncsat_id, nf90_global, 'history', date(1:4)//'-'//date(5:6)//&
                     '-'//date(7:8)//' '//time(1:2)//':'//time(3:4)//' '//zone ) )
    call nf90_err( nf90_put_att(ncsat_id, nf90_global, 'references', &
          'Stohl et al., Atmos. Chem. Phys., 2005, doi:10.5194/acp-5-2461-200;'//&
          'Henne et al., in Lagrangian Modeling of the Atmosphere, 2012, doi:10.1029/2012GM001247') )

    ! end definition
    call nf90_err( nf90_enddef(ncsat_id) )

    ! write dimension variables 
    !**************************

    ! receptor index
    call nf90_err( nf90_put_var(ncsat_id, satrecvar_id, (/(n,n=1,numsatreceptor)/)) )

    ! close file
    call nf90_err( nf90_close(ncsat_id) )

  end subroutine satelliteout_init


  subroutine receptor_output_netcdf()

  !*****************************************************************************
  !                                                                            *
  !     This routine writes receptor concentrations to netcdf files            *
  !                                                                            *
  !     Author: R. Thompson                                                    *
  !     January 2024                                                           *
  !                                                                            *
  !*****************************************************************************

    implicit none

    integer :: ks, ks_start
    character(len=256) :: fn

    ! Open output files
    !******************

    if ((iout.eq.1).or.(iout.eq.3).or.(iout.eq.5)) then
      ! concentration output file
      write(fn,'(A)') path(2)(1:length(2))//'receptor_conc.nc'
    else if (iout.eq.2) then
      ! mixing ratio output file
      ! note for iout=3 (both conc and mixing ratio) only output conc at receptors
      write(fn,'(A)') path(2)(1:length(2))//'receptor_pptv.nc'
    endif
    call nf90_err( nf90_open(trim(fn), nf90_write, nc_id) )

    ! Get variable ids
    !*****************

    if (llcmoutput) then
      ks_start=2
    else
      ks_start=1
    endif

    print*, 'receptor_output_netcdf: nc_id  = ',nc_id

    do ks = ks_start,nspec
      print*, 'receptor_output_netcdf: species, concvar_id = ',trim(species(ks)), concvar_id(ks)
      call nf90_err( nf90_inq_varid(nc_id, trim(species(ks)), concvar_id(ks)) )
      call nf90_err( nf90_inq_varid(nc_id, trim(species(ks))//"_uncert", uncvar_id(ks)) )
      print*, 'receptor_output_netcdf: after read concvar_id = ', concvar_id(ks)
    end do
    call nf90_err( nf90_inq_varid(nc_id, "npart", nnvar_id) )
    call nf90_err( nf90_inq_varid(nc_id, "kernel", xkvar_id) )

  end subroutine receptor_output_netcdf


  subroutine satellite_output_netcdf()

  !*****************************************************************************
  !                                                                            *
  !     This routine writes receptor concentrations to netcdf files            *
  !                                                                            *
  !     Author: R. Thompson                                                    *
  !     January 2024                                                           *
  !                                                                            *
  !*****************************************************************************

    implicit none

    integer :: ks, ks_start
    character(len=256) :: fn

    ! Open output files
    !******************

    write(fn,'(A)') path(2)(1:length(2))//'satellite_pptv.nc'
    call nf90_err( nf90_open(trim(fn), nf90_write, ncsat_id) )

    ! Get variable ids
    !*****************

    if (llcmoutput) then
      ks_start=2
    else
      ks_start=1
    endif

    ! satellite receptors
    do ks = ks_start,nspec
      call nf90_err( nf90_inq_varid(ncsat_id, trim(species(ks)), satvar_id(ks)) )
      call nf90_err( nf90_inq_varid(ncsat_id, trim(species(ks))//"_uncert", satuncvar_id(ks)) )
    end do
    call nf90_err( nf90_inq_varid(ncsat_id, "npart", satnnvar_id) )
    call nf90_err( nf90_inq_varid(ncsat_id, "kernel", satxkvar_id) )

  end subroutine satellite_output_netcdf


  subroutine write_receptor_netcdf(crec,cunc,nnrec,xkrec,nrec)

  !*****************************************************************************
  !                                                                            *
  !     This routine writes receptor concentrations to netcdf files            *
  !                                                                            *
  !     Author: R. Thompson                                                    *
  !     January 2024                                                           *
  !                                                                            *
  !*****************************************************************************

    implicit none
  
    integer :: ks, ks_start, nrec
    real, dimension(nspec,nlayermax) :: crec, cunc
    real, dimension(nlayermax) :: nnrec, xkrec

    if (llcmoutput) then
      ks_start=2
    else
      ks_start=1
    endif

    ! species specific
    do ks=ks_start,nspec
      call nf90_err( nf90_put_var(nc_id, concvar_id(ks), crec(ks,1), (/rpointer/) ) )
      call nf90_err( nf90_put_var(nc_id, uncvar_id(ks), cunc(ks,1), (/rpointer/) ) )
    end do

    ! not species specific output
    ! receptor name
    call nf90_err( nf90_put_var(nc_id, recnamevar_id, receptorname(nrec), (/1,rpointer/), (/16,1/) ) )
    ! receptor time
    call nf90_err( nf90_put_var(nc_id, timevar_id, treceptor(nrec), (/rpointer/) ) )
    ! receptor longitude
    call nf90_err( nf90_put_var(nc_id, reclonvar_id, xreceptor(nrec)*dx+xlon0, (/rpointer/) ) )
    ! receptor latitude
    call nf90_err( nf90_put_var(nc_id, reclatvar_id, yreceptor(nrec)*dy+ylat0, (/rpointer/) ) )
    ! receptor latitude
    call nf90_err( nf90_put_var(nc_id, recaltvar_id, zreceptor(nrec), (/rpointer/) ) )
    ! average number of particles all timesteps for each receptor
    call nf90_err( nf90_put_var(nc_id, nnvar_id, nnrec(1), (/rpointer/) ) )
    ! average kernel all timesteps
    call nf90_err( nf90_put_var(nc_id, xkvar_id, xkrec(1), (/rpointer/) ) )

  end subroutine write_receptor_netcdf

  
  subroutine write_satellite_netcdf(crec,cunc,nnrec,xkrec,nrec)

  !*****************************************************************************
  !                                                                            *
  !     This routine writes receptor concentrations to netcdf files            *
  !                                                                            *
  !     Author: R. Thompson                                                    *
  !     January 2024                                                           *
  !                                                                            *
  !*****************************************************************************

    implicit none

    integer :: ks, ks_start, nrec
    real, dimension(nspec,nlayermax) :: crec, cunc
    real, dimension(nlayermax) :: nnrec, xkrec

    if (llcmoutput) then
      ks_start=2
    else
      ks_start=1
    endif

    ! species specific output
    do ks = ks_start,nspec
      call nf90_err( nf90_put_var(ncsat_id, satvar_id(ks), crec(ks,1:nlayermax-1), (/spointer,1/), (/1,nlayermax-1/) ) )
      call nf90_err( nf90_put_var(ncsat_id, satuncvar_id(ks), cunc(ks,1:nlayermax-1), (/spointer,1/), (/1,nlayermax-1/) ) )
    end do

    ! non-species specific
    ! receptor name
    call nf90_err( nf90_put_var(ncsat_id, satnamevar_id, satellitename(nrec), (/1,spointer/), (/16,1/) ) )
    ! receptor time 
    call nf90_err( nf90_put_var(ncsat_id, sattimevar_id, tsatellite(nrec), (/spointer/) ) )
    ! receptor longitude
    call nf90_err( nf90_put_var(ncsat_id, satlonvar_id, xsatellite(nrec)*dx+xlon0, (/spointer/) ) )
    ! receptor latitude
    call nf90_err( nf90_put_var(ncsat_id, satlatvar_id, ysatellite(nrec)*dy+ylat0, (/spointer/) ) )
    ! receptor altitude
    call nf90_err( nf90_put_var(ncsat_id, sataltvar_id, zsatellite(:,nrec), (/spointer,1/), (/1,(nlayermax)/) ) )
    ! average number of particles all timesteps for each receptor
    call nf90_err( nf90_put_var(ncsat_id, satnnvar_id, nnrec(1:nlayermax-1), (/spointer,1/), (/1,nlayermax-1/) ) )
    ! average kernel all timesteps
    call nf90_err( nf90_put_var(ncsat_id, satxkvar_id, xkrec(1:nlayermax-1), (/spointer,1/), (/1,nlayermax-1/) ) )

  end subroutine write_satellite_netcdf

  
  !*****************************************************************************

  subroutine close_receptor_netcdf

    call nf90_err( nf90_close(nc_id) )

  end subroutine close_receptor_netcdf

  !*****************************************************************************

  subroutine close_satellite_netcdf

    call nf90_err( nf90_close(ncsat_id) )

  end subroutine close_satellite_netcdf

  !*****************************************************************************


  subroutine readreceptors_satellite

  !*****************************************************************************
  !                                                                            *
  !     This routine reads the satellite retrieval information for which       *
  !     mixing ratios should be modelled                                       *
  !                                                                            *
  !     Author: R. Thompson                                                    *
  !     October 2023                                                           *
  !                                                                            *
  !*****************************************************************************

    implicit none

    character(len=256) :: file_name, ppath, pfile, strtmp1, strtmp2
    character(len=256), dimension(:), allocatable :: sat_file, sat_path
    character(len=8) :: anretr, psatname
    character(len=4) :: ayear
    character(len=2) :: amonth, aday
    integer :: jjjjmmdd, hhmmss, yyyy, mm, dd
    integer :: nc_id, dimid, varid
    real    :: xm, ym
    real(kind=dp) :: jul, jd
    real, allocatable, dimension(:,:) :: xpts, ypts, zpt1, zpt2
    integer, allocatable, dimension(:) :: npt, sdate, stime
    integer :: readerror, writeerror, stat
    integer :: j, nn, nr, nretr, nlayer, tmp_numsat
    integer,parameter :: unitreceptorout=2
    logical :: lexist
    real, allocatable, dimension(:) :: tmp_xsat, tmp_ysat, tmp_satarea
    real, allocatable, dimension(:,:) :: tmp_zsat
    integer, allocatable, dimension(:) :: tmp_tsat
    character(len=16), allocatable, dimension(:) :: tmp_satname

    ! declare namelist
    namelist /satellites/ psatname, ppath, pfile

    numsatreceptor=0

    ! For backward runs, do not allow receptor output
    !************************************************

    if (ldirect.lt.0) then
      return
    endif

    ! Open the SATELLITES file and read path and file info 
    !*****************************************************  

    open(unitreceptor,file=path(1)(1:length(1))//'SATELLITES',form='formatted',status='old',iostat=readerror)

    if (readerror.ne.0) then
      write(*,*) 'FLEXPART WARNING readreceptors_satellite: no satellite file found'
      return
    endif

    ! try namelist input
    read(unitreceptor,satellites,iostat=readerror)
    close(unitreceptor)

    if (readerror.ne.0) then
      write(*,*) ' #### FLEXPART ERROR in readreceptors_satellite: #### '
      write(*,*) ' #### error in namelist input #### '
      error stop
    endif ! only namelist input


    ! prepare namelist output if requested
    if (nmlout) then
      open(unitreceptorout,file=path(2)(1:length(2))//'SATELLITES.namelist',&
           &access='append',status='replace',iostat=writeerror)
      if (writeerror.ne.0) then
        write(*,*) ' #### FLEXPART ERROR readreceptors_satellite: cannot open file #### '
        write(*,*) ' #### '//trim(path(2)(1:length(2)))//'SATELLITES.namelist #### '
        error stop
      endif
    endif

    ! Get number of satellites
    !*************************

    open (unitreceptor,file=trim(path(1))//'SATELLITES',status='old',iostat=readerror)
    j=0
    do while (readerror.eq.0)
      read(unitreceptor,satellites,iostat=readerror)
      if (readerror.ne.0) exit
      j=j+1
    end do
    numsatellite=j
    write(*,*) 'Number of satellites: ',numsatellite
    close(unitreceptor)

    ! Allocate arrays
    !****************

    allocate(sat_name(numsatellite),sat_path(numsatellite),sat_file(numsatellite))
    allocate(nnsatlayer(numsatellite))

    ! Read satellite info
    !********************

    open(unitreceptor,file=path(1)(1:length(1))//'SATELLITES',form='formatted',status='old',iostat=readerror)
    j=0
    do while (readerror.eq.0)
      read(unitreceptor,satellites,iostat=readerror)
      if (readerror.ne.0) exit
      j=j+1
      write(*,*) 'readreceptors_satellite: psatname, ppath = ',trim(psatname),', ',trim(ppath)
      sat_name(j)=trim(psatname)
      sat_path(j)=trim(ppath)
      sat_file(j)=trim(pfile)
      ! namelist output
      if (nmlout) then
        write(unitreceptorout,nml=satellites)
      endif
    end do
    close(unitreceptor)
    close(unitreceptorout)

    ! Get number of satellite receptors
    !**********************************

    ! Loop over satellites and days
    jd=bdate
    tmp_numsat=0
    do while (jd.le.edate)
      do j=1,numsatellite
        ! get filename for current day
        call caldate(jd, jjjjmmdd, hhmmss)
        file_name=sat_file(j)
        call getfilename(jjjjmmdd, file_name)
        write(*,*) 'readreceptors_satellite: file_name = ',file_name
        inquire(file=trim(sat_path(j))//'/'//trim(file_name),exist=lexist)
        if (.not.lexist) then
          write(*,*) 'readreceptors_satellite: no retrievals file for '//ayear//amonth//aday
          cycle
        endif
        ! open file
        call nf90_err( nf90_open(trim(sat_path(j))//'/'//trim(file_name), nf90_nowrite, nc_id) )
        ! read dimensions
        call nf90_err( nf90_inq_dimid(nc_id, 'retrieval', dimid) )
        call nf90_err( nf90_inquire_dimension(nc_id, dimid, len=nretr) )
        call nf90_err( nf90_inq_dimid(nc_id, 'nlayer', dimid) )
        call nf90_err( nf90_inquire_dimension(nc_id, dimid, len=nlayer) )
        tmp_numsat=tmp_numsat+nretr
        nlayermax=max(nlayermax,nlayer)
        nnsatlayer(j)=nlayer
        call nf90_err(nf90_close(nc_id))
      end do
      jd=jd+1d0
    end do
    nlayermax=nlayermax+1 ! for levels
    print*, 'readreceptors_satellite: tmp_numsat = ',tmp_numsat

    ! Allocate temporary arrays
    !**************************

    allocate(tmp_xsat(tmp_numsat),tmp_ysat(tmp_numsat),&
             tmp_tsat(tmp_numsat),tmp_satarea(tmp_numsat),&
             tmp_zsat(nlayermax,tmp_numsat),tmp_satname(tmp_numsat))

    ! Read satellite retrievals info
    !*******************************

    ! Loop over days in simulation
    numsatreceptor=0
    jd=bdate
    do while (jd.le.edate)

      do j=1,numsatellite

        ! get filename for current day
        ! assumes same format as output from prep_satellite
        call caldate(jd, jjjjmmdd, hhmmss)
        print*, 'jjjjmmdd = ',jjjjmmdd
        file_name=sat_file(j)
        call getfilename(jjjjmmdd, file_name)
        write(*,*) 'readreceptors_satellite: file_name = ',file_name
        inquire(file=trim(sat_path(j))//'/'//trim(file_name),exist=lexist)
        if (.not.lexist) then
          write(*,*) 'readreceptors_satellite: no retrievals file for '//ayear//amonth//aday
          cycle
        endif

        ! open file
        call nf90_err( nf90_open(trim(sat_path(j))//'/'//trim(file_name), nf90_nowrite, nc_id) )

        ! read dimensions
        call nf90_err( nf90_inq_dimid(nc_id, 'retrieval', dimid) )
        call nf90_err( nf90_inquire_dimension(nc_id, dimid, len=nretr) )
        call nf90_err( nf90_inq_dimid(nc_id, 'nlayer', dimid) )
        call nf90_err( nf90_inquire_dimension(nc_id, dimid, len=nlayer) )

        ! allocate temporary variables
        allocate(sdate(nretr),stat=stat)
        if (stat.ne.0) write(*,*)'ERROR: could not allocate sdate'
        allocate(stime(nretr),stat=stat)
        if (stat.ne.0) write(*,*)'ERROR: could not allocate stime'
        allocate(xpts(nretr,4),stat=stat)
        if (stat.ne.0) write(*,*)'ERROR: could not allocate xpts'
        allocate(ypts(nretr,4),stat=stat)
        if (stat.ne.0) write(*,*)'ERROR: could not allocate ypts'
        allocate(zpt1(nretr,nlayer),stat=stat)
        if (stat.ne.0) write(*,*)'ERROR: could not allocate zpt1'
        allocate(zpt2(nretr,nlayer),stat=stat)
        if (stat.ne.0) write(*,*)'ERROR: could not allocate zpt2'

        ! read coordinate variables      
        call nf90_err( nf90_inq_varid(nc_id,'idate',varid) )
        call nf90_err( nf90_get_var(nc_id,varid,sdate) )
        call nf90_err( nf90_inq_varid(nc_id,'itime',varid) )
        call nf90_err( nf90_get_var(nc_id,varid,stime) )
        call nf90_err( nf90_inq_varid(nc_id,'xpoints',varid) )
        call nf90_err( nf90_get_var(nc_id,varid,xpts) )
        call nf90_err( nf90_inq_varid(nc_id,'ypoints',varid) )
        call nf90_err( nf90_get_var(nc_id,varid,ypts) )
        call nf90_err( nf90_inq_varid(nc_id,'zpoint1',varid) )
        call nf90_err( nf90_get_var(nc_id,varid,zpt1) )
        call nf90_err( nf90_inq_varid(nc_id,'zpoint2',varid) )
        call nf90_err( nf90_get_var(nc_id,varid,zpt2) )
        call nf90_err( nf90_close(nc_id) )

        ! write to coordinates receptor variables
        do nr=1,nretr
          ! time of retrieval from bdate in seconds
          jul=juldate(sdate(nr),stime(nr))
          if ((jul.lt.bdate).or.(jul.gt.edate)) cycle
          numsatreceptor=numsatreceptor+1
          write(anretr,'(I8.8)') numsatreceptor
          tmp_satname(numsatreceptor)=trim(sat_name(j))//'_'//anretr
          tmp_tsat(numsatreceptor)=int((jul-bdate)*24.*3600.) ! time in sec
          ! transform to grid coordinates
          tmp_xsat(numsatreceptor)=(sum(xpts(nr,:))/4.-xlon0)/dx
          tmp_ysat(numsatreceptor)=(sum(ypts(nr,:))/4.-ylat0)/dy
          ! vertical coordinates layer boundaries in Pa
          tmp_zsat(1:nlayer,numsatreceptor)=zpt1(nr,:)
          tmp_zsat(nlayer+1,numsatreceptor)=zpt2(nr,nlayer)
          ! area for mixing ratio calc (used if ind_samp = -1)
          xm=r_earth*cos((sum(ypts(nr,:))/4.)*pi/180.)*dx/180.*pi
          ym=r_earth*dy/180.*pi
          tmp_satarea(numsatreceptor)=xm*ym
          !! test
          if (nr.le.10) then
            write(*,*) 'readreceptors_satellite:'
            write(*,*) 'numsatreceptor,nr = ',numsatreceptor, nr
            write(*,*) 'satellitename = ',tmp_satname(numsatreceptor)
            write(*,*) 'tsatellite = ',tmp_tsat(numsatreceptor)
            write(*,*) 'xsatellite, ysatellite = ',tmp_xsat(numsatreceptor), tmp_ysat(numsatreceptor)
            write(*,*) 'xmidpt, ymidpt = ',sum(xpts(nr,:))/4., sum(ypts(nr,:))/4.
            write(*,*) 'zsatellite = ',tmp_zsat(:,numsatreceptor)
          endif
          !!
        end do ! nretr

        deallocate(sdate, stime, xpts, ypts, zpt1, zpt2)

      end do ! numsatellite

      jd=jd+1d0

    end do ! jd

    write(*,*) 'readreceptors_satellite: numsatreceptor = ',numsatreceptor

    ! Reallocate satellite arrays to actual size
    !*******************************************

    allocate(xsatellite(numsatreceptor),ysatellite(numsatreceptor),&
           tsatellite(numsatreceptor),satellitearea(numsatreceptor),&
           zsatellite(nlayermax,numsatreceptor),satellitename(numsatreceptor))

    xsatellite=tmp_xsat(1:numsatreceptor)
    ysatellite=tmp_ysat(1:numsatreceptor)
    tsatellite=tmp_tsat(1:numsatreceptor)
    satellitearea=tmp_satarea(1:numsatreceptor)
    zsatellite=tmp_zsat(:,1:numsatreceptor)
    satellitename=tmp_satname(1:numsatreceptor)

    deallocate(tmp_xsat,tmp_ysat,tmp_tsat,tmp_satarea,tmp_zsat,tmp_satname)
    deallocate(sat_path,sat_file)

  end subroutine readreceptors_satellite


  subroutine getfilename(jjjjmmdd,file_name)

  !*****************************************************************************
  !                                                                            *
  !    Get actual filename based on dates from generic name                    *
  !                                                                            *
  !    Author: R. Thompson, Sep-2023                                           *
  !                                                                            *
  !                                                                            *
  !*****************************************************************************

    implicit none

    character(len=256) :: file_name, strtmp1, strtmp2
    character(len=4) :: ayear
    character(len=2) :: amonth, aday
    integer :: jjjjmmdd, yyyy, mm, dd, nn

    yyyy=jjjjmmdd/10000
    mm=(jjjjmmdd-yyyy*10000)/100
    dd=jjjjmmdd-yyyy*10000-mm*100
    write(ayear,'(I4)') yyyy
    write(amonth,'(I2.2)') mm
    write(aday,'(I2.2)') dd
    nn=index(file_name,'YYYY',back=.false.)
    if (nn.ne.0) then
      strtmp1=file_name(1:nn-1)
      nn=index(file_name,'YYYY',back=.true.)
      strtmp2=file_name(nn+4:len_trim(file_name))
      file_name=trim(strtmp1)//ayear//trim(strtmp2)
    endif
    nn=index(file_name,'MM',back=.false.)
    if (nn.ne.0) then
      strtmp1=file_name(1:nn-1)
      nn=index(file_name,'MM',back=.true.)
      strtmp2=file_name(nn+2:len_trim(file_name))
      file_name=trim(strtmp1)//amonth//trim(strtmp2)
    endif
    nn=index(file_name,'DD',back=.false.)
    if (nn.ne.0) then
      strtmp1=file_name(1:nn-1)
      nn=index(file_name,'DD',back=.true.)
      strtmp2=file_name(nn+2:len_trim(file_name))
      file_name=trim(strtmp1)//aday//trim(strtmp2)
    endif

  end subroutine getfilename


  subroutine verttransform_satellite

  !*****************************************************************************
  !                                                                            *
  !     This routine transforms the vertical coordinate of the satellite       *
  !     receptors from pressure to height above ground                         *
  !                                                                            *
  !     Author: R. Thompson                                                    *
  !     October 2023                                                           *
  !                                                                            *
  !*****************************************************************************

    implicit none

    integer :: nr, nn, nl, kz, nchar
    integer :: numsatlayer
    integer :: ix, jy, ixp, jyp
    real    :: p1, p2, p3, p4, ddx, ddy, rddx, rddy, dz1, dz2
    real    :: press, pressold, zpres

    print*, 'height = ',height

!$OMP PARALLEL &
!$OMP PRIVATE(nr,nn,nchar,numsatlayer,nl,zpres,ix,jy,ddx,ddy,ixp,jyp,rddx,rddy,&
!$OMP p1,p2,p3,p4,kz,press,pressold,dz1,dz2)

!$OMP DO
    do nr=1,numsatreceptor

      !! test
      if (nr.lt.20) write(*,*) 'zsatellite before transform = ',zsatellite(:,nr)
      !!

      ! get actual number vertical layers for this retrieval
      do nn=1,numsatellite
        nchar=len_trim(sat_name(nn))
        if (satellitename(nr)(1:nchar).eq.trim(sat_name(nn))) then
          numsatlayer=nnsatlayer(nn)
          exit
        endif
      end do

      do nl=1,(numsatlayer+1)

        ! pressure of level in Pa 
        zpres=zsatellite(nl,nr)
        ix=int(xsatellite(nr))
        jy=int(ysatellite(nr))
        ddx=xsatellite(nr)-real(ix)
        ddy=ysatellite(nr)-real(jy)
        ixp=ix+1
        jyp=jy+1
        rddx=1.-ddx
        rddy=1.-ddy
        p1=rddx*rddy
        p2=ddx*rddy
        p3=rddx*ddy
        p4=ddx*ddy

        ! use pressure from second field regardless of time stamp
        ! think this is adequate otherwise need to transform just
        ! for retrievals in recoutaverage interval before calling receptorcalc
        do kz=1,nz
          press=p1*prs(ix,jy,kz,2) &
                 +p2*prs(ixp,jy,kz,2) &
                 +p3*prs(ix,jyp,kz,2) &
                 +p4*prs(ixp,jyp,kz,2)
          if (kz.eq.1) pressold=press
          if (press.lt.zpres) then
            if (kz.eq.1) then
              zsatellite(nl,nr)=height(1)/2.
            else
              dz1=pressold-zpres
              dz2=zpres-press
              zsatellite(nl,nr)=(height(kz-1)*dz2+height(kz)*dz1)/(dz1+dz2)
            endif
            exit  ! found height end loop over nz
          endif
          pressold=press
        end do ! nz
        if ((kz.gt.nz).and.(zpres.le.press)) then
          ! ztra1 press less than top of windfield press
          zsatellite(nl,nr)=0.5*(height(nz-1)+height(nz))
        endif

      end do ! numsatlayer

      !! test
      if (nr.lt.20)  write(*,*) 'zsatellite after transform = ',zsatellite(:,nr)
      !!

    end do ! numsatreceptor
!$OMP END DO
!$OMP END PARALLEL

  end subroutine verttransform_satellite

end module receptor_netcdf_mod

