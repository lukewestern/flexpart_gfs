! SPDX-FileCopyrightText: FLEXPART 1998-2019, see flexpart_license.txt
! SPDX-License-Identifier: GPL-3.0-or-later

subroutine partoutput(itime)!,active_per_rel)
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

  use par_mod
  use com_mod
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
  !integer :: ix,jy,ixp,jyp,indexh,m,il,ind,indz,indzp
  !real :: xlon,ylat,ztemp
  !real :: topo,hmixi,qvi,tri
  !real :: tti,rhoi,pvi
  real :: tr(2),hm(2)
  character :: adate*8,atime*6

  real :: xlon(numpart),ylat(numpart)
  real :: tti(numpart),rhoi(numpart),pvi(numpart),qvi(numpart)
  real :: topo(numpart),hmixi(numpart),tri(numpart)
  !logical  :: active_per_rel(maxpoint)

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
      call update_zcoord(itime, i)
    endif 
  end do

!$OMP END DO
!$OMP END PARALLEL

  write(*,*) 'topo: ', topo(1), 'z:', part(1)%zeta,part(1)%z
  write(*,*) 'xtra,xeta: ', part(1)%xlon
  write(*,*) 'ytra,yeta: ', part(1)%ylat
  write(*,*) pvi(1),qvi(1),tti(1),rhoi(1),part(1)%alive,&
    count%alive,count%spawned,count%terminated


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
!$OMP PARALLEL PRIVATE(j,mythread)
#ifdef USE_NCF
    mythread = omp_get_thread_num()
    if (mythread.eq.thread_divide(1)) call partoutput_netcdf(itime,xlon,'LO',j,ncid)
    if (mythread.eq.thread_divide(2)) call partoutput_netcdf(itime,ylat,'LA',j,ncid)
    if (mythread.eq.thread_divide(3)) call partoutput_netcdf(itime,part(1:numpart)%z,'ZZ',j,ncid)
    !if (mythread.eq.thread_divide(12)) call partoutput_netcdf_int(itime,itramem(1:numpart),'IT',j,ncid)
    if (mythread.eq.thread_divide(4)) call partoutput_netcdf(itime,topo,'TO',j,ncid)
    if (mythread.eq.thread_divide(5)) call partoutput_netcdf(itime,pvi,'PV',j,ncid)
    if (mythread.eq.thread_divide(6)) call partoutput_netcdf(itime,qvi,'QV',j,ncid)
    if (mythread.eq.thread_divide(7)) call partoutput_netcdf(itime,rhoi,'RH',j,ncid)
    if (mythread.eq.thread_divide(8)) call partoutput_netcdf(itime,hmixi,'HM',j,ncid)
    if (mythread.eq.thread_divide(9)) call partoutput_netcdf(itime,tri,'TR',j,ncid)
    if (mythread.eq.thread_divide(10)) call partoutput_netcdf(itime,tti,'TT',j,ncid)
    do j=1,nspec
      if (mythread.eq.mass_divide(j)) call partoutput_netcdf(itime,part(1:numpart)%mass(j),'MA',j,ncid)
    end do
#endif
!$OMP END PARALLEL
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

end subroutine partoutput
