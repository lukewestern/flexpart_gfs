

module coordinates_ecmwf

  use par_mod
  use com_mod
  use windfields_mod

contains

subroutine update_zeta_to_z(itime, ipart)
  use particle_mod
  implicit none 

  integer, intent(in) ::          &
    itime,                        & ! time index
    ipart                           ! particle index

  if (.not. wind_coord_type.eq.'ETA') return

  if (part(ipart)%etaupdate) return

  call zeta_to_z(itime,part(ipart)%xlon,part(ipart)%ylat,part(ipart)%zeta,part(ipart)%z)
  part(ipart)%etaupdate = .true.
  part(ipart)%meterupdate = .true.
end subroutine update_zeta_to_z

subroutine update_z_to_zeta(itime, ipart)
  use particle_mod
  implicit none 

  integer, intent(in) ::          &
    itime,                        & ! time index
    ipart                           ! particle index

  if (.not. part(ipart)%alive) return
  if (.not. wind_coord_type.eq.'ETA') return

  if (part(ipart)%meterupdate) return

  call z_to_zeta(itime,part(ipart)%xlon,part(ipart)%ylat,part(ipart)%z,part(ipart)%zeta)
  part(ipart)%etaupdate = .true.
  part(ipart)%meterupdate = .true.
end subroutine update_z_to_zeta

subroutine z_to_zeta(itime,xt,yt,zold,zteta)
  !                        i    i   o  o  o
  !        o       o       o    i  i  i   o
  !*****************************************************************************
  ! Converting z from eta coordinates to meters                                *
  !                                                                            *
  !*****************************************************************************
  !                                                                            *
  ! Variables:                                                                 *
  ! itime [s]          current temporal position                               *
  ! xteta,yteta,zteta                   spatial position of trajectory         *
  !                                                                            *
  ! etauvheight defined in windfields: half model heights for ETA coordinates  *
  ! Constants:                                                                 *
  !                                                                            *
  !*****************************************************************************
  use interpol_mod

  implicit none
  integer, intent(in) ::          &
    itime                           ! time index
  integer ::                      &
    i,m,indexh                  ! loop indices
  real(kind=dp), intent(in) ::    &
    xt,yt                           ! particle position
  real(kind=dp), intent(in) ::    &
    zold                            ! particle verticle position in eta coordinates
  real(kind=dp), intent(inout) :: &
    zteta                           ! converted output z in meters
  real ::                         &
    frac,                         & ! fraction between z levels
    ztemp1,ztemp2,                & ! z positions of the two encompassing levels
    ttemp_old,ttemp1(2),ttemp_new,& ! storing virtual temperature
    psint1(2),psint                 ! pressure of encompassing levels

  call determine_grid_coordinates(real(xt),real(yt))
  call find_grid_distances(real(xt),real(yt))
  call find_time_variables(itime)

  ! Integration method as used in the original verttransform_ecmwf.f90
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  ztemp1 = 0.
  do i=2,nz-1

    call bilinear_horizontal_interpolation(etauvheight,ttemp1,i,nzmax)
    call temporal_interpolation(ttemp1(1),ttemp1(2),ztemp2)

    if (ztemp2.gt.real(zold)) then
      frac = (real(zold)-ztemp1)/(ztemp2-ztemp1)
      exit
    else if (i.eq.nz-1) then
      frac = 1.
      exit
    endif
    ttemp_old=ttemp_new
    ztemp1=ztemp2
  end do

  zteta=real(uvheight(i-1)*(1.-frac)+uvheight(i)*frac,kind=dp)
end subroutine z_to_zeta

subroutine zeta_to_z(itime,xt,yt,zteta,ztout)
  !                        i    i   o  o  o
  !        o       o       o    i  i  i   o
  !*****************************************************************************
  ! Converting z from eta coordinates to meters                                *
  !                                                                            *
  !*****************************************************************************
  !                                                                            *
  ! Variables:                                                                 *
  ! itime [s]          current temporal position                               *
  ! xt,yt,zteta                   spatial position of trajectory               *
  !                                                                            *
  !                                                                            *
  !                                                                            *
  !*****************************************************************************
  use interpol_mod

  implicit none
  integer, intent(in) ::          &
    itime                           ! time index
  integer ::                      &
    i,j,k,m,indexh                  ! loop indices
  real(kind=dp), intent(in) ::    &
    xt,yt                           ! particle position
  real(kind=dp), intent(in) ::    &
    zteta                           ! particle verticle position in eta coordinates
  real(kind=dp), intent(inout) :: &
    ztout                           ! converted output z in meters
  real(kind=dp) ::                &
    frac                            ! fraction between z levels
  real ::                         &
    ztemp1,ztemp2,                & ! z positions of the two encompassing levels
    ttemp_old,ttemp1(2),ttemp_new,& ! storing virtual temperature
    psint1(2),psint                 ! pressure of encompassing levels
 

  ! Convert eta z coordinate to meters
  !***********************************
  call determine_grid_coordinates(real(xt),real(yt))
  call find_grid_distances(real(xt),real(yt))
  call find_time_variables(itime)

  k=nz-1
  frac=1.
  do k=2,nz-1
    if (zteta.ge.real(uvheight(k),kind=dp)) then
      frac=(zteta-real(uvheight(k-1),kind=dp))/(real(uvheight(k)-uvheight(k-1),kind=dp))
      exit
    endif
  end do

  call bilinear_horizontal_interpolation(etauvheight,ttemp1,k-1,nzmax)
  call temporal_interpolation(ttemp1(1),ttemp1(2),ztemp1)

  call bilinear_horizontal_interpolation(etauvheight,ttemp1,k,nzmax)
  call temporal_interpolation(ttemp1(1),ttemp1(2),ztemp2)

  ztout = real(ztemp1,kind=dp)*(1.-frac)+real(ztemp2,kind=dp)*frac
end subroutine zeta_to_z

end module coordinates_ecmwf