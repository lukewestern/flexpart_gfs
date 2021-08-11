! SPDX-FileCopyrightText: FLEXPART 1998-2019, see flexpart_license.txt
! SPDX-License-Identifier: GPL-3.0-or-later

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
  !                                                                            *
  ! Constants:                                                                 *
  !                                                                            *
  !*****************************************************************************

  use par_mod
  use com_mod
  use interpol_mod

  implicit none
  integer, intent(in) ::          &
    itime                           ! time index
  integer ::                      &
    i,m,indexh                  ! loop indices
  real(kind=dp), intent(in) ::    &
    xt,yt                           ! particle position
  real, intent(in) ::             &
    zold                            ! particle verticle position in eta coordinates
  real, intent(inout) ::          &
    zteta                           ! converted output z in meters
  real ::                         &
    ttemp_old,ttemp1(2),ttemp_new,& ! storing virtual temperature
    ew,                           & ! why does this function need to be declared here?
    ztemp1,ztemp2,                & ! z positions of the two encompassing levels
    frac,                         & ! fraction between z levels
    psint1(2),psint                 ! pressure of encompassing levels

  ix=int(xt)
  jy=int(yt)
  ixp=ix+1
  jyp=jy+1

  ddx=xt-real(ix)
  ddy=yt-real(jy)
  rddx=1.-ddx
  rddy=1.-ddy
  p1=rddx*rddy
  p2=ddx*rddy
  p3=rddx*ddy
  p4=ddx*ddy

  dt1=real(itime-memtime(1))
  dt2=real(memtime(2)-itime)
  dtt=1./(dt1+dt2)


  do m=1,2
    indexh=memind(m)
    psint1(m)=p1*ps(ix ,jy ,1,indexh) &
          +p2*ps(ixp,jy ,1,indexh) &
          +p3*ps(ix ,jyp,1,indexh) &
          +p4*ps(ixp,jyp,1,indexh)
    ttemp1(m)=p1*tt2(ix ,jy ,1,indexh)*(1.+0.378*ew(td2(ix,jy,1,indexh))/ps(ix,jy,1,indexh)) &
          +p2*tt2(ixp,jy ,1,indexh)*(1.+0.378*ew(td2(ixp,jy,1,indexh))/ps(ixp,jy,1,indexh)) &
          +p3*tt2(ix ,jyp,1,indexh)*(1.+0.378*ew(td2(ix,jyp,1,indexh))/ps(ix,jyp,1,indexh)) &
          +p4*tt2(ixp,jyp,1,indexh)*(1.+0.378*ew(td2(ixp,jyp,1,indexh))/ps(ixp,jyp,1,indexh))
  end do
  psint=(psint1(1)*dt2+psint1(2)*dt1)*dtt

  ! Integration method as used in the original verttransform_ecmwf.f90
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  ttemp_old=(ttemp1(1)*dt2+ttemp1(2)*dt1)*dtt
  ztemp1 = 0.
  do i=2,nz-1
    do m=1,2
      indexh=memind(m)
      ttemp1(m)=p1*tteta(ix ,jy ,i,indexh)*(1.+0.608*qveta(ix,jy,i,indexh)) &
            +p2*tteta(ixp,jy ,i,indexh)*(1.+0.608*qveta(ixp,jy,i,indexh)) &
            +p3*tteta(ix ,jyp,i,indexh)*(1.+0.608*qveta(ix,jyp,i,indexh)) &
            +p4*tteta(ixp,jyp,i,indexh)*(1.+0.608*qveta(ixp,jyp,i,indexh))
    end do
    ttemp_new=(ttemp1(1)*dt2+ttemp1(2)*dt1)*dtt

    if (abs(ttemp_new-ttemp_old).gt.0.2) then
      ztemp2=ztemp1+r_air/ga*log((akz(i-1)+bkz(i-1)*psint)/(akz(i)+bkz(i)*psint))* &
        (ttemp_new-ttemp_old)/log(ttemp_new/ttemp_old)
    else
      ztemp2=ztemp1+r_air/ga*log((akz(i-1)+bkz(i-1)*psint)/(akz(i)+bkz(i)*psint))*ttemp_new
    endif

    if (ztemp2.gt.zold) then
      frac = (zold-ztemp1)/(ztemp2-ztemp1)
      exit
    else if (i.eq.nz-1) then
      frac = 1.
      exit
    endif
    ttemp_old=ttemp_new
    ztemp1=ztemp2
  end do

  zteta=uvheight(i-1)*(1.-frac)+uvheight(i)*frac

end subroutine z_to_zeta