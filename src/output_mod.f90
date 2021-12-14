! SPDX-FileCopyrightText: FLEXPART 1998-2019, see flexpart_license.txt
! SPDX-License-Identifier: GPL-3.0-or-later

module output_mod
  
  use com_mod
  use par_mod
  use date_mod
#ifdef USE_NCF  
  use netcdf_output_mod
#endif
  use binary_output_mod
  use txt_output_mod

  implicit none

contains

subroutine initialise_output(itime,filesize)
  implicit none
  
  integer, intent(in) :: itime
  real, intent(inout) :: filesize
#ifdef USE_NCF
  real(kind=dp) ::          &
    jul
  integer ::                &
    jjjjmmdd,ihmmss,i
#endif

  ! Writing header information to either binary or NetCDF format
  if (itime.eq.0) then
    if (grid_output.eq.1) then
#ifdef USE_NCF
      if (lnetcdfout.eq.1) then 
        call writeheader_netcdf(lnest=.false.)
      else 
        call writeheader_binary
      end if

      if (nested_output.eq.1) then
        if (lnetcdfout.eq.1) then
          call writeheader_netcdf(lnest=.true.)
        else
          call writeheader_binary_nest
        endif
      endif
#endif
    endif

    call writeheader_binary ! CHECK ETA
    ! FLEXPART 9.2 ticket ?? write header in ASCII format 
    call writeheader_txt
    !if (nested_output.eq.1) call writeheader_nest
    if (nested_output.eq.1.and.surf_only.ne.1) call writeheader_binary_nest
    if (nested_output.eq.1.and.surf_only.eq.1) call writeheader_binary_nest_surf
    if (nested_output.ne.1.and.surf_only.eq.1) call writeheader_binary_surf

    ! NetCDF only: Create file for storing initial particle positions.
#ifdef USE_NCF
    if (mdomainfill.eq.0) then
      if (ldirect.eq.1) then
        call create_particles_initialoutput(ibtime,ibdate,ibtime,ibdate)
      else
        call create_particles_initialoutput(ietime,iedate,ietime,iedate)
      endif
    endif
    ! Create header files for files that store the particle output
    if (ipout.ge.1) then
      if (ldirect.eq.1) then
        call writeheader_partoutput(ibtime,ibdate,ibtime,ibdate)
      else 
        call writeheader_partoutput(ietime,iedate,ietime,iedate)
      endif
    endif
#endif
  
  ! In case the particle output file is becoming larger than the maximum set
  ! in par_mod, create a new one while keeping track of the filesize.
  else if ((mod(itime,ipoutfac*loutstep).eq.0).and.(ipout.ge.1)) then
#ifdef USE_NCF
    if (filesize.ge.max_partoutput_filesize) then 
      jul=bdate+real(itime,kind=dp)/86400._dp
      call caldate(jul,jjjjmmdd,ihmmss)
      if (ldirect.eq.1) then 
        call writeheader_partoutput(ihmmss,jjjjmmdd,ibtime,ibdate)
      else 
        call writeheader_partoutput(ihmmss,jjjjmmdd,ietime,iedate)
      endif
      filesize = 0.
    endif
    do i=1,numpoint
      filesize = filesize + npart(i)*13.*4./1000000.
    end do
#endif
  endif
end subroutine initialise_output

subroutine finalise_output(itime)
  ! Complete the calculation of initial conditions for particles not yet terminated
  
  implicit none 

  integer, intent(in) :: itime
  integer :: j

  do j=1,numpart
    if (linit_cond.ge.1) call initial_cond_calc(itime,j)
  end do

  if (ipout.eq.2) call output_particles(itime)!,active_per_rel)     ! dump particle positions

  if (linit_cond.ge.1) then
    if(linversionout.eq.1) then
      call initial_cond_output_inversion(itime)   ! dump initial cond. field
    else
      call initial_cond_output(itime)   ! dump initial cond. fielf
    endif
  endif
end subroutine finalise_output

subroutine output_particles(itime)
  !                        i
  !*****************************************************************************
  !                                                                            *
  !     Dump all particle positions                                            *
  !                                                                            *
  !     Author: A. Stohl                                                       *
  !                                                                            *
  !     12 March 1999                                                          *
  !                                                                            *
  !*****************************************************************************
  !                                                                            *
  ! Variables:                                                                 *
  !                                                                            *
  !*****************************************************************************

  use interpol_mod
  use coordinates_ecmwf
  use particle_mod
#ifdef USE_NCF
  use netcdf
  use netcdf_output_mod, only: partoutput_netcdf,open_partoutput_file,close_partoutput_file
  use omp_lib, only: OMP_GET_THREAD_NUM
#endif

  implicit none

  real(kind=dp) :: jul
  integer :: itime,i,j,jjjjmmdd,ihmmss
  real :: tr(2),hm(2)
  character :: adate*8,atime*6

  real :: xlon(numpart),ylat(numpart),ztemp1,ztemp2
  real :: tti(numpart),rhoi(numpart),pvi(numpart),qvi(numpart)
  real :: topo(numpart),hmixi(numpart),tri(numpart),ztemp(numpart)
  real :: masstemp(numpart,nspec)

#ifdef USE_NCF
  integer  :: ncid, mythread, thread_divide(12),mass_divide(nspec)
#endif

  ! Some variables needed for temporal interpolation
  !*************************************************
  call find_time_variables(itime)

!$OMP PARALLEL PRIVATE(i,tr,hm)
!$OMP DO
  do i=1,numpart
  ! Take only valid particles
  !**************************
    xlon(i)=-1.
    ylat(i)=-1.
    tti(i)=-1.
    rhoi(i)=-1.
    pvi(i)=-1.
    qvi(i)=-1.
    topo(i)=-1.
    hmixi(i)=-1.
    tri(i)=-1.
    ztemp(i)=-1.
    do j=1,nspec
      masstemp(i,j)=-1.
    end do
    if (part(i)%alive) then
      xlon(i)=xlon0+part(i)%xlon*dx
      ylat(i)=ylat0+part(i)%ylat*dy

  !*****************************************************************************
  ! Interpolate several variables (PV, specific humidity, etc.) to particle position
  !*****************************************************************************
      call determine_grid_coordinates(real(part(i)%xlon),real(part(i)%ylat))
      call find_grid_distances(real(part(i)%xlon),real(part(i)%ylat))
  ! Topography
  !***********
      call bilinear_horizontal_interpolation_2dim(oro,topo(i))

      ! First set dz1out from interpol_mod to -1 so it only is calculated once per particle
      !************************************************************************************
      dz1out=-1
      ! Potential vorticity
      call interpol_partoutput_value('PV',pvi(i),i)
      ! Specific humidity
      call interpol_partoutput_value('QV',qvi(i),i)
      ! Temperature
      call interpol_partoutput_value('TT',tti(i),i)
      ! Density
      call interpol_partoutput_value('RH',rhoi(i),i)
      ! Reset dz1out
      !*************
      dz1out=-1

  ! Tropopause and PBL height
  !**************************
  ! Tropopause
      call bilinear_horizontal_interpolation(tropopause,tr,1,1)
      call temporal_interpolation(tr(1),tr(2),tri(i))
  ! PBL height
      call bilinear_horizontal_interpolation(hmix,hm,1,1)
      call temporal_interpolation(hm(1),hm(2),hmixi(i))


  ! Convert eta z coordinate to meters if necessary
  !************************************************
      call update_zeta_to_z(itime, i)
      ztemp(i)=part(i)%z

  ! Assign the masses
  !******************
      do j=1,nspec
        masstemp(i,j)=part(i)%mass(j)
      end do
    endif 
  end do

!$OMP END DO
!$OMP END PARALLEL
  if (numpart.gt.0) then
    write(*,*) 'topo: ', topo(1), 'z:', part(1)%zeta,part(1)%z
    write(*,*) 'xtra,xeta: ', part(1)%xlon
    write(*,*) 'ytra,yeta: ', part(1)%ylat
    write(*,*) pvi(1),qvi(1),tti(1),rhoi(1),part(1)%alive,&
      count%alive,count%spawned,count%terminated
  endif

  ! Determine current calendar date, needed for the file name
  !**********************************************************

  jul=bdate+real(itime,kind=dp)/86400._dp
  call caldate(jul,jjjjmmdd,ihmmss)
  write(adate,'(i8.8)') jjjjmmdd
  write(atime,'(i6.6)') ihmmss

  if (lnetcdfout.eq.1) then
  ! open output file
    call open_partoutput_file(ncid)

    ! Dividing the openmp threads for writing
    j=0
    do i=1,10
      if (j.eq.numthreads) j = 0
      thread_divide(i) = j
      j = j + 1
    end do
    do i=1,nspec
      if (j.eq.numthreads) j = 0
      mass_divide(i) = j
      j = j + 1
    end do

    ! First allocate the time and particle dimention within the netcdf file
    call partoutput_netcdf(itime,xlon,'TI',j,ncid)
    call partoutput_netcdf(itime,xlon,'PA',j,ncid)

    ! Fill the fields in parallel
    if (numpart.gt.0) then
!$OMP PARALLEL PRIVATE(j,mythread)
#ifdef USE_NCF
      mythread = omp_get_thread_num()
      if (mythread.eq.thread_divide(1)) call partoutput_netcdf(itime,xlon,'LO',j,ncid)
      if (mythread.eq.thread_divide(2)) call partoutput_netcdf(itime,ylat,'LA',j,ncid)
      if (mythread.eq.thread_divide(3)) call partoutput_netcdf(itime,ztemp,'ZZ',j,ncid)
      !if (mythread.eq.thread_divide(12)) call partoutput_netcdf_int(itime,itramem(1:numpart),'IT',j,ncid)
      if (mythread.eq.thread_divide(4)) call partoutput_netcdf(itime,topo,'TO',j,ncid)
      if (mythread.eq.thread_divide(5)) call partoutput_netcdf(itime,pvi,'PV',j,ncid)
      if (mythread.eq.thread_divide(6)) call partoutput_netcdf(itime,qvi,'QV',j,ncid)
      if (mythread.eq.thread_divide(7)) call partoutput_netcdf(itime,rhoi,'RH',j,ncid)
      if (mythread.eq.thread_divide(8)) call partoutput_netcdf(itime,hmixi,'HM',j,ncid)
      if (mythread.eq.thread_divide(9)) call partoutput_netcdf(itime,tri,'TR',j,ncid)
      if (mythread.eq.thread_divide(10)) call partoutput_netcdf(itime,tti,'TT',j,ncid)
      do j=1,nspec
        if (mythread.eq.mass_divide(j)) call partoutput_netcdf(itime,masstemp(:,j),'MA',j,ncid)
      end do
#endif
!$OMP END PARALLEL
    endif
    call close_partoutput_file(ncid)
  else
    ! Open output file and write the output
    !**************************************

    if (ipout.eq.1.or.ipout.eq.3) then
      open(unitpartout,file=path(2)(1:length(2))//'partposit_'//adate// &
           atime,form='unformatted')
    else
      open(unitpartout,file=path(2)(1:length(2))//'partposit_end', &
           form='unformatted')
    endif

    ! Write current time to file
    !***************************

    write(unitpartout) itime
    do i=1,numpart
    ! Take only valid particles
    !**************************

      if (part(i)%alive) then
    ! Write the output
    !*****************      
        write(unitpartout) part(i)%npoint,xlon(i),ylat(i),part(i)%z, &
             part(i)%tstart,topo(i),pvi(i),qvi(i),rhoi(i),hmixi(i),tri(i),tti(i), &
             (part(i)%mass(j),j=1,nspec)
      endif
    end do


    write(unitpartout) -99999,-9999.9,-9999.9,-9999.9,-99999, &
         -9999.9,-9999.9,-9999.9,-9999.9,-9999.9,-9999.9,-9999.9, &
         (-9999.9,j=1,nspec)


    close(unitpartout)
  endif
end subroutine output_particles

subroutine output_concentrations(itime,loutstart,loutend,loutnext,outnum)
  use unc_mod
  use outg_mod
  use par_mod
  use com_mod
#ifdef USE_NCF
  use netcdf_output_mod, only: concoutput_netcdf,concoutput_nest_netcdf,&
       &concoutput_surf_netcdf,concoutput_surf_nest_netcdf
#endif
  use binary_output_mod 

  implicit none

  integer,intent(in) ::     &
    itime                     ! time index
  integer,intent(inout) ::  &
    loutstart,loutend,      & ! concentration calculation starting and ending time
    loutnext
  real,intent(inout) ::     &
    outnum                    ! concentration calculation sample number
  real(sp) ::               &
    gridtotalunc              ! concentration calculation related
  real(dep_prec) ::         &
    wetgridtotalunc,        & ! concentration calculation related
    drygridtotalunc           ! concentration calculation related
  real ::                   &
    weight                    ! concentration calculation sample weight

  ! Is the time within the computation interval, if not, return
  !************************************************************
  if ((ldirect*itime.lt.ldirect*loutstart).or.(ldirect*itime.gt.ldirect*loutend)) then
    return
  endif

  ! If we are exactly at the start or end of the concentration averaging interval,
  ! give only half the weight to this sample
  !*****************************************************************************
  if (mod(itime-loutstart,loutsample).eq.0) then
    if ((itime.eq.loutstart).or.(itime.eq.loutend)) then
      weight=0.5
    else
      weight=1.0
    endif
    outnum=outnum+weight
    call conccalc(itime,weight)
  endif

  ! If it is not time yet to write outputs, return
  !***********************************************
  if ((itime.ne.loutend).or.(outnum.le.0)) then
    return
  endif

  ! Output and reinitialization of grid
  ! If necessary, first sample of new grid is also taken
  !*****************************************************
  if ((iout.le.3.).or.(iout.eq.5)) then
    if (surf_only.ne.1) then 
#ifdef USE_NCF
      call concoutput_netcdf(itime,outnum,gridtotalunc,wetgridtotalunc,drygridtotalunc)
#else
      call concoutput(itime,outnum,gridtotalunc,wetgridtotalunc,drygridtotalunc)
#endif
    else
#ifdef USE_NCF
      call concoutput_surf_netcdf(itime,outnum,gridtotalunc,wetgridtotalunc,drygridtotalunc)
#else
      if (linversionout.eq.1) then
        call concoutput_inversion(itime,outnum,gridtotalunc,wetgridtotalunc,drygridtotalunc)
      else
        call concoutput_surf(itime,outnum,gridtotalunc,wetgridtotalunc,drygridtotalunc)
      endif
#endif
    endif

    if (nested_output .eq. 1) then
#ifdef USE_NCF
      if (surf_only.ne.1) then
        call concoutput_nest_netcdf(itime,outnum)
      else 
        call concoutput_surf_nest_netcdf(itime,outnum)
      endif
#else
      if (surf_only.ne.1) then
        call concoutput_nest(itime,outnum)
      else 
        if(linversionout.eq.1) then
          call concoutput_inversion_nest(itime,outnum)
        else 
          call concoutput_surf_nest(itime,outnum)
        endif
      endif
#endif
    endif
    outnum=0.
  endif

  write(*,45) itime,numpart,gridtotalunc,wetgridtotalunc,drygridtotalunc

45      format(i13,' Seconds simulated: ',i13, ' Particles:    Uncertainty: ',3f7.3)

  loutnext=loutnext+loutstep
  loutstart=loutnext-loutaver/2
  loutend=loutnext+loutaver/2
  if (itime.eq.loutstart) then
    weight=0.5
    outnum=outnum+weight
    call conccalc(itime,weight)
  endif
end subroutine output_concentrations

end module output_mod