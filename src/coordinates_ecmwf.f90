

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
  if (.not. part(ipart)%alive) return  
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

  if (.not. wind_coord_type.eq.'ETA') return
  if (.not. part(ipart)%alive) return
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
    i,m,k,n                         ! loop indices
  real(kind=dp), intent(in) ::    &
    xt,yt                           ! particle position
  real(kind=dp), intent(in) ::    &
    zold                            ! particle verticle position in eta coordinates
  real(kind=dp), intent(inout) :: &
    zteta                           ! converted output z in meters
  real ::                         &
    frac,                         & ! fraction between z levels
    ztemp1,ztemp2,                & ! z positions of the two encompassing levels
    ttemp1(2),                    & ! storing virtual temperature
    psint1(2),psint                 ! pressure of encompassing levels
  real ::                         &
    prx,pr1,pr2     ! pressure of encompassing levels

  if (.not. logarithmic_interpolation) then
    call z_to_zeta_lin(itime,xt,yt,zold,zteta)
    return
  endif

  call determine_grid_coordinates(real(xt),real(yt))
  call find_grid_distances(real(xt),real(yt))
  call find_time_variables(itime)

  ! Integration method as used in the original verttransform_ecmwf.f90
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! First estimate the level it is at, to reduce computation time
  n=nz-3
  do i=2,nz-1
    if ((etauvheight(ix,jy,i,memind(1)).gt.real(zold)) .or. &
        (etauvheight(ixp,jy,i,memind(1)).gt.real(zold)) .or. & 
        (etauvheight(ix,jyp,i,memind(1)).gt.real(zold)) .or. &
        (etauvheight(ixp,jyp,i,memind(1)).gt.real(zold))) then
      n=i-2
      exit
    endif
  end do
  n=max(n,2)

  ztemp1 = 0.
  do i=n,nz-1
    k=i
    do m=1,2
      call horizontal_interpolation(etauvheight,ttemp1(m),i,memind(m),nzmax)
    end do
    call temporal_interpolation(ttemp1(1),ttemp1(2),ztemp2)

    if (ztemp2.gt.real(zold)) then
      !frac = (real(zold)-ztemp1)/(ztemp2-ztemp1)
      exit
    else if (i.eq.nz-1) then
      frac = 1.
      exit
    endif
    ztemp1=ztemp2
  end do

  if (k.lt.nz-1) then 
    do m=1,2
      call horizontal_interpolation(ps,psint1(m),1,memind(m),1)
    end do
    call temporal_interpolation(psint1(1),psint1(2),psint)  
    pr1=akz(k-1) + bkz(k-1)*psint
    pr2=akz(k) + bkz(k)*psint

    prx=pr1/exp(log(pr2/pr1)/(ztemp2-ztemp1)*ztemp1) * exp(log(pr2/pr1)/(ztemp2-ztemp1)*real(zold))
    frac=(prx-pr1)/(pr2 - pr1)
  endif

  zteta=real(uvheight(k-1)*(1.-frac)+uvheight(k)*frac,kind=dp)
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
    ttemp1(2),                    & ! storing virtual temperature
    psint1(2),psint,prx,pr1,pr2     ! pressure of encompassing levels
 
  if (.not. logarithmic_interpolation) then
    call zeta_to_z_lin(itime,xt,yt,zteta,ztout)
    return
  endif

  ! Convert eta z coordinate to meters
  !***********************************
  call determine_grid_coordinates(real(xt),real(yt))
  call find_grid_distances(real(xt),real(yt))
  call find_time_variables(itime)

  k=nz-1
  frac=1.
  do i=2,nz-1
    k=i
    if (zteta.ge.real(uvheight(k),kind=dp)) then
      frac=(zteta-real(uvheight(k-1),kind=dp))/(real(uvheight(k)-uvheight(k-1),kind=dp))
      exit
    endif
  end do

  do m=1,2
    call horizontal_interpolation(ps,psint1(m),1,memind(m),1)
  end do
  call temporal_interpolation(psint1(1),psint1(2),psint)  
  pr1=akz(k-1) + bkz(k-1)*psint
  pr2=akz(k) + bkz(k)*psint
  prx=pr1*(1.-frac) + pr2*frac
  
  do m=1,2
    call horizontal_interpolation(etauvheight,ttemp1(m),k-1,memind(m),nzmax)
  end do
  call temporal_interpolation(ttemp1(1),ttemp1(2),ztemp1)

  do m=1,2
    call horizontal_interpolation(etauvheight,ttemp1(m),k,memind(m),nzmax)
  end do
  call temporal_interpolation(ttemp1(1),ttemp1(2),ztemp2)
  
  if ((pr2.eq.0).or.(pr1.eq.0)) then
    ztout = real(ztemp1,kind=dp)*(1.-frac)+real(ztemp2,kind=dp)*frac
    return
  endif
  
  ztout = ztemp1 + (ztemp2-ztemp1)/log(pr2/pr1)*log(prx/pr1)
end subroutine zeta_to_z

subroutine z_to_zeta_lin(itime,xt,yt,zold,zteta)
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
    i,m,k,n                         ! loop indices
  real(kind=dp), intent(in) ::    &
    xt,yt                           ! particle position
  real(kind=dp), intent(in) ::    &
    zold                            ! particle verticle position in eta coordinates
  real(kind=dp), intent(inout) :: &
    zteta                           ! converted output z in meters
  real ::                         &
    frac,                         & ! fraction between z levels
    ztemp1,ztemp2,                & ! z positions of the two encompassing levels
    ttemp1(2),                    & ! storing virtual temperature
    psint1(2),psint                 ! pressure of encompassing levels
  real ::                         &
    prx,pr1,pr2     ! pressure of encompassing levels

  call determine_grid_coordinates(real(xt),real(yt))
  call find_grid_distances(real(xt),real(yt))
  call find_time_variables(itime)

  ! Integration method as used in the original verttransform_ecmwf.f90
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! First estimate the level it is at, to reduce computation time
  n=nz-3
  do i=2,nz-1
    if ((etauvheight(ix,jy,i,memind(1)).gt.real(zold)) .or. &
        (etauvheight(ixp,jy,i,memind(1)).gt.real(zold)) .or. &
        (etauvheight(ix,jyp,i,memind(1)).gt.real(zold)) .or. &
        (etauvheight(ixp,jyp,i,memind(1)).gt.real(zold))) then
      n=i-2
      exit
    endif
  end do
  n=max(n,2)

  ztemp1 = 0.
  do i=n,nz-1
    k=i
    do m=1,2
      call horizontal_interpolation(etauvheight,ttemp1(m),i,memind(m),nzmax)
    end do
    call temporal_interpolation(ttemp1(1),ttemp1(2),ztemp2)

    if (ztemp2.gt.real(zold)) then
      frac = (real(zold)-ztemp1)/(ztemp2-ztemp1)
      exit
    else if (i.eq.nz-1) then
      frac = 1.
      exit
    endif
    ztemp1=ztemp2
  end do

  zteta=real(uvheight(k-1)*(1.-frac)+uvheight(k)*frac,kind=dp)
end subroutine z_to_zeta_lin

subroutine zeta_to_z_lin(itime,xt,yt,zteta,ztout)
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
    ttemp1(2),                    & ! storing virtual temperature
    psint1(2),psint                 ! pressure of encompassing levels
 

  ! Convert eta z coordinate to meters
  !***********************************
  call determine_grid_coordinates(real(xt),real(yt))
  call find_grid_distances(real(xt),real(yt))
  call find_time_variables(itime)

  k=nz-1
  frac=1.
  do i=2,nz-1
    k=i
    if (zteta.ge.real(uvheight(k),kind=dp)) then
      frac=(zteta-real(uvheight(k-1),kind=dp))/(real(uvheight(k)-uvheight(k-1),kind=dp))
      exit
    endif
  end do
  
  do m=1,2
    call horizontal_interpolation(etauvheight,ttemp1(m),k-1,memind(m),nzmax)
  end do
  call temporal_interpolation(ttemp1(1),ttemp1(2),ztemp1)

  do m=1,2
    call horizontal_interpolation(etauvheight,ttemp1(m),k,memind(m),nzmax)
  end do
  call temporal_interpolation(ttemp1(1),ttemp1(2),ztemp2)
  
  ztout = real(ztemp1,kind=dp)*(1.-frac)+real(ztemp2,kind=dp)*frac
end subroutine zeta_to_z_lin

end module coordinates_ecmwf