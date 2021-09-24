! SPDX-FileCopyrightText: FLEXPART 1998-2019, see flexpart_license.txt
! SPDX-License-Identifier: GPL-3.0-or-later


subroutine partpos_average(itime,j)


!**********************************************************************
! This subroutine averages particle quantities, to be used for particle
! dump (in partoutput.f90). Averaging is done over output interval.
!**********************************************************************

  use par_mod
  use com_mod
  use interpol_mod
  use coordinates_ecmwf

  implicit none

  integer :: itime,j
  real :: xlon,ylat,x,y,z
  real :: topo,hm(2),hmixi,pvi,qvi
  real :: tti,rhoi,ttemp
  real :: uui,vvi
  real :: tr(2),tri!,energy



 ! Some variables needed for temporal interpolation
  !*************************************************
  call find_time_variables(itime)

  xlon=xlon0+part(j)%xlon*dx
  ylat=ylat0+part(j)%ylat*dy

  !*****************************************************************************
  ! Interpolate several variables (PV, specific humidity, etc.) to particle position
  !*****************************************************************************

  call determine_grid_coordinates(real(part(j)%xlon),real(part(j)%ylat))
  call find_grid_distances(real(part(j)%xlon),real(part(j)%ylat))

  ! Topography
  !***********
  call bilinear_horizontal_interpolation_2dim(oro,topo)

  ! Potential vorticity, specific humidity, temperature, and density
  !*****************************************************************
  ! First set dz1out from interpol_mod to -1 so it only is calculated once per particle
  !************************************************************************************
  dz1out=-1
  ! Potential vorticity
  call interpol_partoutput_value('PV',pvi,j)
  ! Specific humidity
  call interpol_partoutput_value('QV',qvi,j)
  ! Temperature
  call interpol_partoutput_value('TT',tti,j)
  ! U wind
  call interpol_partoutput_value('UU',uui,j)
  ! V wind
  call interpol_partoutput_value('VV',vvi,j)
  ! Density
  call interpol_partoutput_value('RH',rhoi,j)
  ! Reset dz1out
  !*************
  dz1out=-1
  ! Tropopause and PBL height
  !**************************
  ! Tropopause
  call bilinear_horizontal_interpolation(tropopause,tr,1,1)
  call temporal_interpolation(tr(1),tr(2),tri)
  ! PBL height
  call bilinear_horizontal_interpolation(hmix,hm,1,1)
  call temporal_interpolation(hm(1),hm(2),hmixi)

  ! Convert eta z coordinate to meters if necessary. Can be moved to output only
  !************************************************
  if (wind_coord_type.eq.'ETA') then
    call update_zcoord(itime,j)
  endif

  ! energy=tti*cpa+(ztemp1+topo)*9.81+qvi*2501000.+(uui**2+vvi**2)/2.

  ! Add new values to sum and increase counter by one
  !**************************************************

  npart_av(j)=npart_av(j)+1

  ! Calculate Cartesian 3D coordinates suitable for averaging
  !**********************************************************

  xlon=xlon*pi180
  ylat=ylat*pi180
  x = cos(ylat)*sin(xlon)
  y = -1.*cos(ylat)*cos(xlon)
  z = sin(ylat)


  if (j.eq.1) then
    write(*,*) 'topo: ', topo, 'z:', part(j)%zeta,part(j)%z
    write(*,*) 'xtra,xeta: ', part(j)%xlon
    write(*,*) 'ytra,yeta: ', part(j)%ylat
    write(*,*) pvi,qvi,tti,uui,vvi,rhoi
  endif

  part_av_cartx(j)=part_av_cartx(j)+x
  part_av_carty(j)=part_av_carty(j)+y
  part_av_cartz(j)=part_av_cartz(j)+z
  part_av_z(j)=part(j)%z!part_av_z(j)+ztemp1
  part_av_topo(j)=part_av_topo(j)+topo
  part_av_pv(j)=part_av_pv(j)+pvi
  part_av_qv(j)=part_av_qv(j)+qvi
  part_av_tt(j)=part_av_tt(j)+tti
  part_av_uu(j)=part_av_uu(j)+uui
  part_av_vv(j)=part_av_vv(j)+vvi
  part_av_rho(j)=part_av_rho(j)+rhoi
  part_av_tro(j)=part_av_tro(j)+tri
  part_av_hmix(j)=part_av_hmix(j)+hmixi
  ! part_av_energy(j)=part_av_energy(j)+energy

return
end subroutine partpos_average 
