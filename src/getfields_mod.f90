! SPDX-FileCopyrightText: FLEXPART 1998-2019, see flexpart_license.txt
! SPDX-License-Identifier: GPL-3.0-or-later
module getfields_mod
  
  use par_mod 
  use com_mod
  use windfields_mod

  implicit none

  private :: obukhov,richardson,scalev,calcpar,calcpar_nests,calcpv,calcpv_nests

  public :: getfields
contains

subroutine getfields(itime,nstop,metdata_format)
  !                       i     o
  !*****************************************************************************
  !                                                                            *
  !  This subroutine manages the 3 data fields to be kept in memory.           *
  !  During the first time step of petterssen it has to be fulfilled that the  *
  !  first data field must have |wftime|<itime, i.e. the absolute value of     *
  !  wftime must be smaller than the absolute value of the current time in [s].*
  !  The other 2 fields are the next in time after the first one.              *
  !  Pointers (memind) are used, because otherwise one would have to resort the*
  !  wind fields, which costs a lot of computing time. Here only the pointers  *
  !  are resorted.                                                             *
  !                                                                            *
  !     Author: A. Stohl                                                       *
  !                                                                            *
  !     29 April 1994                                                          *
  !                                                                            *
  !*****************************************************************************
  !  Changes, Bernd C. Krueger, Feb. 2001:                                     *
  !   Variables tth,qvh,tthn,qvhn (on eta coordinates) in common block.        *
  !   Function of nstop extended.                                              *
  !                                                                            *
  !   Unified ECMWF and GFS builds                                             *
  !   Marian Harustak, 12.5.2017                                               *
  !     - Added passing of metdata_format as it was needed by called routines  *
  !*****************************************************************************
  !                                                                            *
  ! Variables:                                                                 *
  ! lwindinterval [s]    time difference between the two wind fields read in   *
  ! indj                 indicates the number of the wind field to be read in  *
  ! indmin               remembers the number of wind fields already treated   *
  ! memind(2)            pointer, on which place the wind fields are stored    *
  ! memtime(2) [s]       times of the wind fields, which are kept in memory    *
  ! itime [s]            current time since start date of trajectory calcu-    *
  !                      lation                                                *
  ! nstop                > 0, if trajectory has to be terminated               *
  ! nx,ny,nuvz,nwz       field dimensions in x,y and z direction               *
  ! uu(0:nxmax,0:nymax,nuvzmax,2)   wind components in x-direction [m/s]       *
  ! vv(0:nxmax,0:nymax,nuvzmax,2)   wind components in y-direction [m/s]       *
  ! ww(0:nxmax,0:nymax,nwzmax,2)    wind components in z-direction [deltaeta/s]*
  ! tt(0:nxmax,0:nymax,nuvzmax,2)   temperature [K]                            *
  ! ps(0:nxmax,0:nymax,2)           surface pressure [Pa]                      *
  ! metdata_format     format of metdata (ecmwf/gfs)                           *
  !                                                                            *
  ! Constants:                                                                 *
  ! idiffmax             maximum allowable time difference between 2 wind      *
  !                      fields                                                *
  !                                                                            *
  !*****************************************************************************

  use class_gribfile

  implicit none

  integer :: indj,itime,nstop,memaux
  integer :: metdata_format

  real :: uuh(0:nxmax-1,0:nymax-1,nuvzmax)
  real :: vvh(0:nxmax-1,0:nymax-1,nuvzmax)
  real :: pvh(0:nxmax-1,0:nymax-1,nuvzmax)
  real :: wwh(0:nxmax-1,0:nymax-1,nwzmax)
  real :: uuhn(0:nxmaxn-1,0:nymaxn-1,nuvzmax,maxnests)
  real :: vvhn(0:nxmaxn-1,0:nymaxn-1,nuvzmax,maxnests)
  real :: pvhn(0:nxmaxn-1,0:nymaxn-1,nuvzmax,maxnests)
  real :: wwhn(0:nxmaxn-1,0:nymaxn-1,nwzmax,maxnests)
  ! RLT added partial pressure water vapor 
  real :: pwater(0:nxmax-1,0:nymax-1,nzmax,numwfmem)
  integer :: kz, ix
  character(len=100) :: rowfmt

  integer :: indmin = 1


  ! Check, if wind fields are available for the current time step
  !**************************************************************

  nstop=0
  if ((ldirect*wftime(1).gt.ldirect*itime).or. &
       (ldirect*wftime(numbwf).lt.ldirect*itime)) then
    write(*,*) 'FLEXPART WARNING: NO WIND FIELDS ARE AVAILABLE.'
    write(*,*) 'A TRAJECTORY HAS TO BE TERMINATED.'
    nstop=4
    return
  endif


  if ((ldirect*memtime(1).le.ldirect*itime).and. &
       (ldirect*memtime(2).gt.ldirect*itime)) then

  ! The right wind fields are already in memory -> don't do anything
  !*****************************************************************

    continue

  else if ((ldirect*memtime(2).le.ldirect*itime).and. &
       (memtime(2).ne.999999999)) then


  ! Current time is after 2nd wind field
  ! -> Resort wind field pointers, so that current time is between 1st and 2nd
  !***************************************************************************

    memaux=memind(1)
    memind(1)=memind(2)
    memind(2)=memaux
    memtime(1)=memtime(2)


  ! Read a new wind field and store it on place memind(2)
  !******************************************************

    do indj=indmin,numbwf-1
      if (ldirect*wftime(indj+1).gt.ldirect*itime) then
        if (metdata_format.eq.GRIBFILE_CENTRE_ECMWF) then
          call SYSTEM_CLOCK(count_clock, count_rate, count_max)
          s_temp = (count_clock - count_clock0)/real(count_rate)
          call readwind_ecmwf(indj+1,memind(2),uuh,vvh,wwh)
          call SYSTEM_CLOCK(count_clock, count_rate, count_max)
          s_readwind = s_readwind + ((count_clock - count_clock0)/real(count_rate)-s_temp)
        else
          call readwind_gfs(indj+1,memind(2),uuh,vvh,wwh)
        end if
        call readwind_nests(indj+1,memind(2),uuhn,vvhn,wwhn)
        call calcpar(memind(2),uuh,vvh,pvh,metdata_format)
        call calcpar_nests(memind(2),uuhn,vvhn,pvhn,metdata_format)
        if (metdata_format.eq.GRIBFILE_CENTRE_ECMWF) then
          call verttransform_ecmwf(memind(2),uuh,vvh,wwh,pvh)
        else
          call verttransform_gfs(memind(2),uuh,vvh,wwh,pvh)
        end if
        call verttransform_nests(memind(2),uuhn,vvhn,wwhn,pvhn)
        memtime(2)=wftime(indj+1)
        nstop = 1
        exit
      endif
    end do
    indmin=indj

    if (WETBKDEP) then
      call writeprecip(itime,memind(1))
    endif

  else

  ! No wind fields, which can be used, are currently in memory
  ! -> read both wind fields
  !***********************************************************

    do indj=indmin,numbwf-1
      if ((ldirect*wftime(indj).le.ldirect*itime).and. &
           (ldirect*wftime(indj+1).gt.ldirect*itime)) then
        memind(1)=1
        if (metdata_format.eq.GRIBFILE_CENTRE_ECMWF) then
          call SYSTEM_CLOCK(count_clock, count_rate, count_max)
          s_temp = (count_clock - count_clock0)/real(count_rate)
          call readwind_ecmwf(indj,memind(1),uuh,vvh,wwh)
          call SYSTEM_CLOCK(count_clock, count_rate, count_max)
          s_readwind = s_readwind + ((count_clock - count_clock0)/real(count_rate)-s_temp)
        else
          call readwind_gfs(indj,memind(1),uuh,vvh,wwh)
        end if
        call readwind_nests(indj,memind(1),uuhn,vvhn,wwhn)
        call calcpar(memind(1),uuh,vvh,pvh,metdata_format)
        call calcpar_nests(memind(1),uuhn,vvhn,pvhn,metdata_format)
        if (metdata_format.eq.GRIBFILE_CENTRE_ECMWF) then
          call verttransform_ecmwf(memind(1),uuh,vvh,wwh,pvh)
        else
          call verttransform_gfs(memind(1),uuh,vvh,wwh,pvh)
        end if
        call verttransform_nests(memind(1),uuhn,vvhn,wwhn,pvhn)
        memtime(1)=wftime(indj)
        memind(2)=2
        if (metdata_format.eq.GRIBFILE_CENTRE_ECMWF) then
          call SYSTEM_CLOCK(count_clock, count_rate, count_max)
          s_temp = (count_clock - count_clock0)/real(count_rate)
          call readwind_ecmwf(indj+1,memind(2),uuh,vvh,wwh)
          call SYSTEM_CLOCK(count_clock, count_rate, count_max)
          s_readwind = s_readwind + ((count_clock - count_clock0)/real(count_rate)-s_temp)
        else
          call readwind_gfs(indj+1,memind(2),uuh,vvh,wwh)
        end if
        call readwind_nests(indj+1,memind(2),uuhn,vvhn,wwhn)
        call calcpar(memind(2),uuh,vvh,pvh,metdata_format)
        call calcpar_nests(memind(2),uuhn,vvhn,pvhn,metdata_format)
        if (metdata_format.eq.GRIBFILE_CENTRE_ECMWF) then
          call verttransform_ecmwf(memind(2),uuh,vvh,wwh,pvh)
        else
          call verttransform_gfs(memind(2),uuh,vvh,wwh,pvh)
        end if
        call verttransform_nests(memind(2),uuhn,vvhn,wwhn,pvhn)
        memtime(2)=wftime(indj+1)
        nstop = 1
        exit
      endif
    end do
    indmin=indj

    if (WETBKDEP) then
      call writeprecip(itime,memind(1))
    endif

  end if

  ! RLT calculate dry air density
  pwater=qv*prs/((r_air/r_water)*(1.-qv)+qv)
  rho_dry=(prs-pwater)/(r_air*tt)

  lwindinterv=abs(memtime(2)-memtime(1))

  if (lwindinterv.gt.idiffmax) nstop=3
end subroutine getfields

subroutine calcpv(n,uuh,vvh,pvh)
  !               i  i   i   o
  !*****************************************************************************
  !                                                                            *
  !  Calculation of potential vorticity on 3-d grid.                           *
  !                                                                            *
  !     Author: P. James                                                       *
  !     3 February 2000                                                        *
  !                                                                            *
  !     Adaptation to FLEXPART, A. Stohl, 1 May 2000                           *
  !                                                                            *
  !*****************************************************************************
  !                                                                            *
  ! Variables:                                                                 *
  ! n                  temporal index for meteorological fields (1 to 2)       *
  !                                                                            *
  ! Constants:                                                                 *
  !                                                                            *
  !*****************************************************************************

  implicit none

  integer :: n,ix,jy,i,j,k,kl,ii,jj,klvrp,klvrm,klpt,kup,kdn,kch
  integer :: jyvp,jyvm,ixvp,ixvm,jumpx,jumpy,jux,juy,ivrm,ivrp,ivr
  integer :: nlck
  real :: vx(2),uy(2),phi,tanphi,cosphi,dvdx,dudy,f
  real :: theta,thetap,thetam,dthetadp,dt1,dt2,dt
  real :: pvavr,ppml(0:nxmax-1,0:nymax-1,nuvzmax),ppmk(0:nxmax-1,0:nymax-1,nuvzmax)
  real :: thup,thdn
  real,parameter :: eps=1.e-5, p0=101325
  real :: uuh(0:nxmax-1,0:nymax-1,nuvzmax)
  real :: vvh(0:nxmax-1,0:nymax-1,nuvzmax)
  real :: pvh(0:nxmax-1,0:nymax-1,nuvzmax)

  ! Set number of levels to check for adjacent theta
  nlck=nuvz/3
  !
  ! Loop over entire grid
  !**********************
  do kl=1,nuvz
    do jy=0,nymin1
      do ix=0,nxmin1
         ppml(ix,jy,kl)=akz(kl)+bkz(kl)*ps(ix,jy,1,n)
      enddo
    enddo
  enddo

!  ppmk(:,:,1:nuvz)=(100000./ppml(:,:,1:nuvz))**kappa
  ppmk(0:nxmin1,0:nymin1,1:nuvz)=(100000./ppml(0:nxmin1,0:nymin1,1:nuvz))**kappa
!$OMP PARALLEL PRIVATE(jy,ix,kl,phi,f,tanphi,cosphi,jyvp,jyvm,jumpy,juy, &
!$OMP ixvp,ixvm,jumpx,ivrp,ivrm,jux,theta,klvrp,klvrm,klpt,thetap,thetam,dthetadp, &
!$OMP ii,i,ivr,kdn,kch,kup,thdn,thup,dt1,dt2,dt,vx,k,dvdx, &
!$OMP jj,j,uy,dudy)
!$OMP DO
  do jy=0,nymin1
    if (sglobal.and.jy.eq.0) goto 10
    if (nglobal.and.jy.eq.nymin1) goto 10
    phi = (ylat0 + jy * dy) * pi / 180.
    f = 0.00014585 * sin(phi)
    tanphi = tan(phi)
    cosphi = cos(phi)
  ! Provide a virtual jy+1 and jy-1 in case we are on domain edge (Lat)
    jyvp=jy+1
    jyvm=jy-1
    if (jy.eq.0) jyvm=0
    if (jy.eq.nymin1) jyvp=nymin1
  ! Define absolute gap length
    jumpy=2
    if (jy.eq.0.or.jy.eq.nymin1) jumpy=1
    if (sglobal.and.jy.eq.1) then
      jyvm=1
      jumpy=1
    end if
    if (nglobal.and.jy.eq.ny-2) then
      jyvp=ny-2
      jumpy=1
    end if
    juy=jumpy
  !
    do ix=0,nxmin1
  ! Provide a virtual ix+1 and ix-1 in case we are on domain edge (Long)
      ixvp=ix+1
      ixvm=ix-1
      jumpx=2
      if (xglobal) then
         ivrp=ixvp
         ivrm=ixvm
         if (ixvm.lt.0) ivrm=ixvm+nxmin1
         if (ixvp.ge.nx) ivrp=ixvp-nx+1
      else
        if (ix.eq.0) ixvm=0
        if (ix.eq.nxmin1) ixvp=nxmin1
        ivrp=ixvp
        ivrm=ixvm
  ! Define absolute gap length
        if (ix.eq.0.or.ix.eq.nxmin1) jumpx=1
      end if
      jux=jumpx
  !
  ! Loop over the vertical
  !***********************

      do kl=1,nuvz
        theta=tth(ix,jy,kl,n)*ppmk(ix,jy,kl)
        klvrp=kl+1
        klvrm=kl-1
        klpt=kl
  ! If top or bottom level, dthetadp is evaluated between the current
  ! level and the level inside, otherwise between level+1 and level-1
  !
        if (klvrp.gt.nuvz) klvrp=nuvz
        if (klvrm.lt.1) klvrm=1
        thetap=tth(ix,jy,klvrp,n)*ppmk(ix,jy,klvrp)
        thetam=tth(ix,jy,klvrm,n)*ppmk(ix,jy,klvrm)
        dthetadp=(thetap-thetam)/(ppml(ix,jy,klvrp)-ppml(ix,jy,klvrm))

  ! Compute vertical position at pot. temperature surface on subgrid
  ! and the wind at that position
  !*****************************************************************
  ! a) in x direction
        ii=0
        do i=ixvm,ixvp,jumpx
          ivr=i
          if (xglobal) then
             if (i.lt.0) ivr=ivr+nxmin1
             if (i.ge.nx) ivr=ivr-nx+1
          end if
          ii=ii+1
  ! Search adjacent levels for current theta value
  ! Spiral out from current level for efficiency
          kup=klpt-1
          kdn=klpt
          kch=0
40        continue
  ! Upward branch
          kup=kup+1
          if (kch.ge.nlck) goto 21     ! No more levels to check,
  !                                       ! and no values found
          if (kup.ge.nuvz) goto 41
          kch=kch+1
          k=kup
          thdn=tth(ivr,jy,k,n)*ppmk(ivr,jy,k)
          thup=tth(ivr,jy,k+1,n)*ppmk(ivr,jy,k+1)


          if (((thdn.ge.theta).and.(thup.le.theta)).or. &
          ((thdn.le.theta).and.(thup.ge.theta))) then
            dt1=abs(theta-thdn)
            dt2=abs(theta-thup)
            dt=dt1+dt2
            if (dt.lt.eps) then   ! Avoid division by zero error
              dt1=0.5             ! G.W., 10.4.1996
              dt2=0.5
              dt=1.0
            endif
            vx(ii)=(vvh(ivr,jy,k)*dt2+vvh(ivr,jy,k+1)*dt1)/dt
            goto 20
          endif
41        continue
  ! Downward branch
          kdn=kdn-1
          if (kdn.lt.1) goto 40
          kch=kch+1
          k=kdn
          thdn=tth(ivr,jy,k,n)*ppmk(ivr,jy,k)
          thup=tth(ivr,jy,k+1,n)*ppmk(ivr,jy,k+1)

          if (((thdn.ge.theta).and.(thup.le.theta)).or. &
          ((thdn.le.theta).and.(thup.ge.theta))) then
            dt1=abs(theta-thdn)
            dt2=abs(theta-thup)
            dt=dt1+dt2
            if (dt.lt.eps) then   ! Avoid division by zero error
              dt1=0.5             ! G.W., 10.4.1996
              dt2=0.5
              dt=1.0
            endif
            vx(ii)=(vvh(ivr,jy,k)*dt2+vvh(ivr,jy,k+1)*dt1)/dt
            goto 20
          endif
          goto 40
  ! This section used when no values were found
21        continue
  ! Must use vv at current level and long. jux becomes smaller by 1
          vx(ii)=vvh(ix,jy,kl)
          jux=jux-1
  ! Otherwise OK
20        continue
        end do
        if (jux.gt.0) then
          dvdx=(vx(2)-vx(1))/real(jux)/(dx*pi/180.)
        else
          dvdx=vvh(ivrp,jy,kl)-vvh(ivrm,jy,kl)
          dvdx=dvdx/real(jumpx)/(dx*pi/180.)
  ! Only happens if no equivalent theta value
  ! can be found on either side, hence must use values
  ! from either side, same pressure level.
        end if

  ! b) in y direction

        jj=0
        do j=jyvm,jyvp,jumpy
          jj=jj+1
  ! Search adjacent levels for current theta value
  ! Spiral out from current level for efficiency
          kup=klpt-1
          kdn=klpt
          kch=0
70        continue
  ! Upward branch
          kup=kup+1
          if (kch.ge.nlck) goto 51     ! No more levels to check,
  !                                     ! and no values found
          if (kup.ge.nuvz) goto 71
          kch=kch+1
          k=kup
          thdn=tth(ix,j,k,n)*ppmk(ix,j,k)
          thup=tth(ix,j,k+1,n)*ppmk(ix,j,k+1)
          if (((thdn.ge.theta).and.(thup.le.theta)).or. &
          ((thdn.le.theta).and.(thup.ge.theta))) then
            dt1=abs(theta-thdn)
            dt2=abs(theta-thup)
            dt=dt1+dt2
            if (dt.lt.eps) then   ! Avoid division by zero error
              dt1=0.5             ! G.W., 10.4.1996
              dt2=0.5
              dt=1.0
            endif
            uy(jj)=(uuh(ix,j,k)*dt2+uuh(ix,j,k+1)*dt1)/dt
            goto 50
          endif
71        continue
  ! Downward branch
          kdn=kdn-1
          if (kdn.lt.1) goto 70
          kch=kch+1
          k=kdn
          thdn=tth(ix,j,k,n)*ppmk(ix,j,k)
          thup=tth(ix,j,k+1,n)*ppmk(ix,j,k+1)
          if (((thdn.ge.theta).and.(thup.le.theta)).or. &
          ((thdn.le.theta).and.(thup.ge.theta))) then
            dt1=abs(theta-thdn)
            dt2=abs(theta-thup)
            dt=dt1+dt2
            if (dt.lt.eps) then   ! Avoid division by zero error
              dt1=0.5             ! G.W., 10.4.1996
              dt2=0.5
              dt=1.0
            endif
            uy(jj)=(uuh(ix,j,k)*dt2+uuh(ix,j,k+1)*dt1)/dt
            goto 50
          endif
          goto 70
  ! This section used when no values were found
51        continue
  ! Must use uu at current level and lat. juy becomes smaller by 1
          uy(jj)=uuh(ix,jy,kl)
          juy=juy-1
  ! Otherwise OK
50        continue
        end do
      if (juy.gt.0) then
      dudy=(uy(2)-uy(1))/real(juy)/(dy*pi/180.)
      else
      dudy=uuh(ix,jyvp,kl)-uuh(ix,jyvm,kl)
      dudy=dudy/real(jumpy)/(dy*pi/180.)
      end if
  !
      pvh(ix,jy,kl)=dthetadp*(f+(dvdx/cosphi-dudy &
           +uuh(ix,jy,kl)*tanphi)/r_earth)*(-1.e6)*9.81


  !
  ! Resest jux and juy
      jux=jumpx
      juy=jumpy
      end do
    end do
10  continue
  end do
!$OMP END DO 
!$OMP END PARALLEL
  !
  ! Fill in missing PV values on poles, if present
  ! Use mean PV of surrounding latitude ring
  !
  if (sglobal) then
     do kl=1,nuvz
        pvavr=0.
        do ix=0,nxmin1
           pvavr=pvavr+pvh(ix,1,kl)
        end do
        pvavr=pvavr/real(nx)
        jy=0
        do ix=0,nxmin1
           pvh(ix,jy,kl)=pvavr
        end do
     end do
  end if
  if (nglobal) then
     do kl=1,nuvz
        pvavr=0.
        do ix=0,nxmin1
           pvavr=pvavr+pvh(ix,ny-2,kl)
        end do
        pvavr=pvavr/real(nx)
        jy=nymin1
        do ix=0,nxmin1
           pvh(ix,jy,kl)=pvavr
        end do
     end do
  end if

end subroutine calcpv

subroutine calcpv_nests(l,n,uuhn,vvhn,pvhn)
  !                     i i  i    i    o
  !*****************************************************************************
  !                                                                            *
  !  Calculation of potential vorticity on 3-d nested grid                     *
  !                                                                            *
  !     Author: P. James                                                       *
  !     22 February 2000                                                       *
  !                                                                            *
  !*****************************************************************************
  !                                                                            *
  ! Variables:                                                                 *
  ! n                  temporal index for meteorological fields (1 to 2)       *
  ! l                  index of current nest                                   *
  !                                                                            *
  ! Constants:                                                                 *
  !                                                                            *
  !*****************************************************************************
  
  implicit none

  integer :: n,ix,jy,i,j,k,kl,ii,jj,klvrp,klvrm,klpt,kup,kdn,kch
  integer :: jyvp,jyvm,ixvp,ixvm,jumpx,jumpy,jux,juy,ivrm,ivrp,ivr
  integer :: nlck,l
  real :: vx(2),uy(2),phi,tanphi,cosphi,dvdx,dudy,f
  real :: theta,thetap,thetam,dthetadp,dt1,dt2,dt
  real :: ppml(0:nxmaxn-1,0:nymaxn-1,nuvzmax),ppmk(0:nxmaxn-1,0:nymaxn-1,nuvzmax)
  real :: thup,thdn
  real,parameter :: eps=1.e-5,p0=101325
  real :: uuhn(0:nxmaxn-1,0:nymaxn-1,nuvzmax,maxnests)
  real :: vvhn(0:nxmaxn-1,0:nymaxn-1,nuvzmax,maxnests)
  real :: pvhn(0:nxmaxn-1,0:nymaxn-1,nuvzmax,maxnests)

  ! Set number of levels to check for adjacent theta
  nlck=nuvz/3
  !
  ! Loop over entire grid
  !**********************
  do kl=1,nuvz
    do jy=0,nyn(l)-1
      do ix=0,nxn(l)-1
         ppml(ix,jy,kl)=akz(kl)+bkz(kl)*psn(ix,jy,1,n,l)
      enddo
    enddo
  enddo
  !  ppmk=(100000./ppml)**kappa
  ppmk(0:nxn(l)-1,0:nyn(l)-1,1:nuvz)=(100000./ppml(0:nxn(l)-1,0:nyn(l)-1,1:nuvz))**kappa

  do jy=0,nyn(l)-1
    phi = (ylat0n(l) + jy * dyn(l)) * pi / 180.
    f = 0.00014585 * sin(phi)
    tanphi = tan(phi)
    cosphi = cos(phi)
  ! Provide a virtual jy+1 and jy-1 in case we are on domain edge (Lat)
      jyvp=jy+1
      jyvm=jy-1
      if (jy.eq.0) jyvm=0
      if (jy.eq.nyn(l)-1) jyvp=nyn(l)-1
  ! Define absolute gap length
      jumpy=2
      if (jy.eq.0.or.jy.eq.nyn(l)-1) jumpy=1
      juy=jumpy
  !
    do ix=0,nxn(l)-1
  ! Provide a virtual ix+1 and ix-1 in case we are on domain edge (Long)
      ixvp=ix+1
      ixvm=ix-1
      jumpx=2
      if (ix.eq.0) ixvm=0
      if (ix.eq.nxn(l)-1) ixvp=nxn(l)-1
      ivrp=ixvp
      ivrm=ixvm
  ! Define absolute gap length
      if (ix.eq.0.or.ix.eq.nxn(l)-1) jumpx=1
      jux=jumpx
  !
  ! Loop over the vertical
  !***********************

      do kl=1,nuvz
        theta=tthn(ix,jy,kl,n,l)*ppmk(ix,jy,kl)
        klvrp=kl+1
        klvrm=kl-1
        klpt=kl
  ! If top or bottom level, dthetadp is evaluated between the current
  ! level and the level inside, otherwise between level+1 and level-1
  !
        if (klvrp.gt.nuvz) klvrp=nuvz
        if (klvrm.lt.1) klvrm=1
        thetap=tthn(ix,jy,klvrp,n,l)*ppmk(ix,jy,klvrp)
        thetam=tthn(ix,jy,klvrm,n,l)*ppmk(ix,jy,klvrm)
        dthetadp=(thetap-thetam)/(ppml(ix,jy,klvrp)-ppml(ix,jy,klvrm))

  ! Compute vertical position at pot. temperature surface on subgrid
  ! and the wind at that position
  !*****************************************************************
  ! a) in x direction
        ii=0
        do i=ixvm,ixvp,jumpx
          ivr=i
          ii=ii+1
  ! Search adjacent levels for current theta value
  ! Spiral out from current level for efficiency
          kup=klpt-1
          kdn=klpt
          kch=0
40        continue
  ! Upward branch
          kup=kup+1
          if (kch.ge.nlck) goto 21     ! No more levels to check,
  !                                       ! and no values found
          if (kup.ge.nuvz) goto 41
          kch=kch+1
          k=kup
          thdn=tthn(ivr,jy,k,n,l)*ppmk(ivr,jy,k)
          thup=tthn(ivr,jy,k+1,n,l)*ppmk(ivr,jy,k+1)

      if (((thdn.ge.theta).and.(thup.le.theta)).or. &
           ((thdn.le.theta).and.(thup.ge.theta))) then
              dt1=abs(theta-thdn)
              dt2=abs(theta-thup)
              dt=dt1+dt2
              if (dt.lt.eps) then   ! Avoid division by zero error
                dt1=0.5             ! G.W., 10.4.1996
                dt2=0.5
                dt=1.0
              endif
    vx(ii)=(vvhn(ivr,jy,k,l)*dt2+vvhn(ivr,jy,k+1,l)*dt1)/dt
              goto 20
            endif
41        continue
  ! Downward branch
          kdn=kdn-1
          if (kdn.lt.1) goto 40
          kch=kch+1
          k=kdn
          thdn=tthn(ivr,jy,k,n,l)*ppmk(ivr,jy,k)
          thup=tthn(ivr,jy,k+1,n,l)*ppmk(ivr,jy,k+1)
      if (((thdn.ge.theta).and.(thup.le.theta)).or. &
           ((thdn.le.theta).and.(thup.ge.theta))) then
              dt1=abs(theta-thdn)
              dt2=abs(theta-thup)
              dt=dt1+dt2
              if (dt.lt.eps) then   ! Avoid division by zero error
                dt1=0.5             ! G.W., 10.4.1996
                dt2=0.5
                dt=1.0
              endif
    vx(ii)=(vvhn(ivr,jy,k,l)*dt2+vvhn(ivr,jy,k+1,l)*dt1)/dt
              goto 20
            endif
            goto 40
  ! This section used when no values were found
21      continue
  ! Must use vv at current level and long. jux becomes smaller by 1
        vx(ii)=vvhn(ix,jy,kl,l)
        jux=jux-1
  ! Otherwise OK
20        continue
        end do
      if (jux.gt.0) then
      dvdx=(vx(2)-vx(1))/real(jux)/(dxn(l)*pi/180.)
      else
      dvdx=vvhn(ivrp,jy,kl,l)-vvhn(ivrm,jy,kl,l)
      dvdx=dvdx/real(jumpx)/(dxn(l)*pi/180.)
  ! Only happens if no equivalent theta value
  ! can be found on either side, hence must use values
  ! from either side, same pressure level.
      end if

  ! b) in y direction

        jj=0
        do j=jyvm,jyvp,jumpy
          jj=jj+1
  ! Search adjacent levels for current theta value
  ! Spiral out from current level for efficiency
          kup=klpt-1
          kdn=klpt
          kch=0
70        continue
  ! Upward branch
          kup=kup+1
          if (kch.ge.nlck) goto 51     ! No more levels to check,
  !                                     ! and no values found
          if (kup.ge.nuvz) goto 71
          kch=kch+1
          k=kup
          thdn=tthn(ix,j,k,n,l)*ppmk(ix,j,k)
          thup=tthn(ix,j,k+1,n,l)*ppmk(ix,j,k+1)
      if (((thdn.ge.theta).and.(thup.le.theta)).or. &
           ((thdn.le.theta).and.(thup.ge.theta))) then
              dt1=abs(theta-thdn)
              dt2=abs(theta-thup)
              dt=dt1+dt2
              if (dt.lt.eps) then   ! Avoid division by zero error
                dt1=0.5             ! G.W., 10.4.1996
                dt2=0.5
                dt=1.0
              endif
        uy(jj)=(uuhn(ix,j,k,l)*dt2+uuhn(ix,j,k+1,l)*dt1)/dt
              goto 50
            endif
71        continue
  ! Downward branch
          kdn=kdn-1
          if (kdn.lt.1) goto 70
          kch=kch+1
          k=kdn
          thdn=tthn(ix,j,k,n,l)*ppmk(ix,j,k)
          thup=tthn(ix,j,k+1,n,l)*ppmk(ix,j,k+1)
      if (((thdn.ge.theta).and.(thup.le.theta)).or. &
           ((thdn.le.theta).and.(thup.ge.theta))) then
              dt1=abs(theta-thdn)
              dt2=abs(theta-thup)
              dt=dt1+dt2
              if (dt.lt.eps) then   ! Avoid division by zero error
                dt1=0.5             ! G.W., 10.4.1996
                dt2=0.5
                dt=1.0
              endif
        uy(jj)=(uuhn(ix,j,k,l)*dt2+uuhn(ix,j,k+1,l)*dt1)/dt
              goto 50
            endif
            goto 70
  ! This section used when no values were found
51      continue
  ! Must use uu at current level and lat. juy becomes smaller by 1
        uy(jj)=uuhn(ix,jy,kl,l)
        juy=juy-1
  ! Otherwise OK
50        continue
        end do
      if (juy.gt.0) then
      dudy=(uy(2)-uy(1))/real(juy)/(dyn(l)*pi/180.)
      else
      dudy=uuhn(ix,jyvp,kl,l)-uuhn(ix,jyvm,kl,l)
      dudy=dudy/real(jumpy)/(dyn(l)*pi/180.)
      end if
  !
      pvhn(ix,jy,kl,l)=dthetadp*(f+(dvdx/cosphi-dudy &
           +uuhn(ix,jy,kl,l)*tanphi)/r_earth)*(-1.e6)*9.81

  !
  ! Resest jux and juy
      jux=jumpx
      juy=jumpy
      end do
    end do
  end do
  !
end subroutine calcpv_nests

subroutine calcpar(n,uuh,vvh,pvh,metdata_format)
  !                   i  i   i   o
  !*****************************************************************************
  !                                                                            *
  !     Computation of several boundary layer parameters needed for the        *
  !     dispersion calculation and calculation of dry deposition velocities.   *
  !     All parameters are calculated over the entire grid.                    *
  !                                                                            *
  !     Author: A. Stohl                                                       *
  !                                                                            *
  !     21 May 1995                                                            *
  !                                                                            *
  ! ------------------------------------------------------------------         *
  !     Petra Seibert, Feb 2000:                                               *
  !     convection scheme:                                                     *
  !     new variables in call to richardson                                    *
  !                                                                            *
  !   CHANGE 17/11/2005 Caroline Forster NCEP GFS version                      *
  !                                                                            *
  !   Changes, Bernd C. Krueger, Feb. 2001:                                    *
  !    Variables tth and qvh (on eta coordinates) in common block              *
  !                                                                            *
  !   Unified ECMWF and GFS builds                                             *
  !   Marian Harustak, 12.5.2017                                               *
  !     - Merged calcpar and calcpar_gfs into one routine using if-then        *
  !       for meteo-type dependent code                                        *
  !*****************************************************************************

  !*****************************************************************************
  !                                                                            *
  ! Variables:                                                                 *
  ! n                  temporal index for meteorological fields (1 to 3)       *
  ! uuh                                                                        *
  ! vvh                                                                        *
  ! pvh                                                                        *
  ! metdata_format     format of metdata (ecmwf/gfs)                           * 
  !                                                                            *
  ! Constants:                                                                 *
  !                                                                            *
  !                                                                            *
  ! Functions:                                                                 *
  ! scalev             computation of ustar                                    *
  ! obukhov            computatio of Obukhov length                            *
  !                                                                            *
  !*****************************************************************************

  use par_mod
  use com_mod
  use class_gribfile
  use drydepo_mod
  use windfields_mod
  use qvsat_mod

  implicit none

  integer :: metdata_format
  integer :: n,ix,jy,i,kz,lz,kzmin,llev,loop_start
  real :: ttlev(nuvzmax),qvlev(nuvzmax),ol,hmixplus
  real :: ulev(nuvzmax),vlev(nuvzmax),rh,vd(maxspec),subsceff,ylat
  real :: altmin,tvold,pold,zold,pint,tv,zlev(nuvzmax),hmixdummy,akzdummy
  real :: uuh(0:nxmax-1,0:nymax-1,nuvzmax)
  real :: vvh(0:nxmax-1,0:nymax-1,nuvzmax)
  real :: pvh(0:nxmax-1,0:nymax-1,nuvzmax)
  real,parameter :: const=r_air/ga

  !write(*,*) 'in calcpar writting snowheight'
  !***********************************
  ! for test: write out snow depths

  ! open(4,file='slandusetest',form='formatted')
  ! do 5 ix=0,nxmin1
  !5       write (4,*) (sd(ix,jy,1,n),jy=0,nymin1)
  !  close(4)


  ! Loop over entire grid
  !**********************

  ! openmp change. Issue with vdep (with O2+openmp), this needs to be fixed before this can be switched on.
  ! !$OMP PARALLEL PRIVATE(jy,ix,ulev,vlev,ttlev,qvlev,llev,ylat,ol,i,hmixplus, &
  ! !$OMP subsceff,vd,kz,lz,zlev,rh,kzmin,pold,zold,tvold,pint,tv,loop_start)

  ! !$OMP DO
  do jy=0,nymin1

  ! Set minimum height for tropopause
  !**********************************

    ylat=ylat0+real(jy)*dy
    if ((ylat.ge.-20.).and.(ylat.le.20.)) then
      altmin = 5000.
    else
      if ((ylat.gt.20.).and.(ylat.lt.40.)) then
        altmin=2500.+(40.-ylat)*125.
      else if ((ylat.gt.-40.).and.(ylat.lt.-20.)) then
        altmin=2500.+(40.+ylat)*125.
      else
        altmin=2500.
      endif
    endif

    do ix=0,nxmin1

  ! 1) Calculation of friction velocity
  !************************************

      ustar(ix,jy,1,n)=scalev(ps(ix,jy,1,n),tt2(ix,jy,1,n), &
           td2(ix,jy,1,n),surfstr(ix,jy,1,n))
      if (ustar(ix,jy,1,n).le.1.e-8) ustar(ix,jy,1,n)=1.e-8

  ! 2) Calculation of inverse Obukhov length scale
  !***********************************************

      if (metdata_format.eq.GRIBFILE_CENTRE_NCEP) then
        ! NCEP version: find first level above ground
        llev = 0
        do i=1,nuvz
          if (ps(ix,jy,1,n).lt.akz(i)) llev=i
        end do
        llev = llev+1
        if (llev.gt.nuvz) llev = nuvz-1
        ! NCEP version

        ! calculate inverse Obukhov length scale with tth(llev)
        ol=obukhov(ps(ix,jy,1,n),tt2(ix,jy,1,n),td2(ix,jy,1,n), &
             tth(ix,jy,llev,n),ustar(ix,jy,1,n),sshf(ix,jy,1,n),akm,bkm,akz(llev),metdata_format)
      else
        llev=0
        ol=obukhov(ps(ix,jy,1,n),tt2(ix,jy,1,n),td2(ix,jy,1,n), &
            tth(ix,jy,2,n),ustar(ix,jy,1,n),sshf(ix,jy,1,n),akm,bkm,akzdummy,metdata_format)
      end if

      if (ol.ne.0.) then
        oli(ix,jy,1,n)=1./ol
      else
        oli(ix,jy,1,n)=99999.
      endif


  ! 3) Calculation of convective velocity scale and mixing height
  !**************************************************************

      do i=1,nuvz
        ulev(i)=uuh(ix,jy,i)
        vlev(i)=vvh(ix,jy,i)
        ttlev(i)=tth(ix,jy,i,n)
        qvlev(i)=qvh(ix,jy,i,n)
      end do

      if (metdata_format.eq.GRIBFILE_CENTRE_NCEP) then
        ! NCEP version hmix has been read in in readwind.f, is therefore not calculated here
      call richardson(ps(ix,jy,1,n),ustar(ix,jy,1,n),ttlev,qvlev, &
           ulev,vlev,nuvz,akz,bkz,sshf(ix,jy,1,n),tt2(ix,jy,1,n), &
             td2(ix,jy,1,n),hmixdummy,wstar(ix,jy,1,n),hmixplus,metdata_format)
      else
        call richardson(ps(ix,jy,1,n),ustar(ix,jy,1,n),ttlev,qvlev, &
             ulev,vlev,nuvz,akz,bkz,sshf(ix,jy,1,n),tt2(ix,jy,1,n), &
             td2(ix,jy,1,n),hmix(ix,jy,1,n),wstar(ix,jy,1,n),hmixplus,metdata_format)
      end if

      if(lsubgrid.eq.1) then
        subsceff=min(excessoro(ix,jy),hmixplus)
      else
        subsceff=0.0
      endif
  !
  ! CALCULATE HMIX EXCESS ACCORDING TO SUBGRIDSCALE VARIABILITY AND STABILITY
  !
      hmix(ix,jy,1,n)=hmix(ix,jy,1,n)+subsceff
      hmix(ix,jy,1,n)=max(hmixmin,hmix(ix,jy,1,n)) ! set minimum PBL height
      hmix(ix,jy,1,n)=min(hmixmax,hmix(ix,jy,1,n)) ! set maximum PBL height

  ! 4) Calculation of dry deposition velocities
  !********************************************

      if (DRYDEP) then
  ! Sabine Eckhardt, Dec 06: use new index for z0 for water depending on
  ! windspeed
        z0(7)=0.016*ustar(ix,jy,1,n)*ustar(ix,jy,1,n)/ga

  ! Calculate relative humidity at surface
  !***************************************
        rh=ew(td2(ix,jy,1,n),ps(ix,jy,1,n))/ew(tt2(ix,jy,1,n),ps(ix,jy,1,n))

        call getvdep(n,ix,jy,ustar(ix,jy,1,n), &
             tt2(ix,jy,1,n),ps(ix,jy,1,n),1./oli(ix,jy,1,n), &
             ssr(ix,jy,1,n),rh,lsprec(ix,jy,1,n)+convprec(ix,jy,1,n), &
             sd(ix,jy,1,n),vd)

        do i=1,nspec
          vdep(ix,jy,i,n)=vd(i)
        end do

      endif

  !******************************************************
  ! Calculate height of thermal tropopause (Hoinka, 1997)
  !******************************************************

  ! 1) Calculate altitudes of model levels
  !***************************************

      tvold=tt2(ix,jy,1,n)*(1.+0.378*ew(td2(ix,jy,1,n),ps(ix,jy,1,n))/ &
           ps(ix,jy,1,n))
      pold=ps(ix,jy,1,n)
      zold=0.
      if (metdata_format.eq.GRIBFILE_CENTRE_ECMWF) then
        loop_start=2
      else
        loop_start=llev
      end if
      do kz=loop_start,nuvz
        pint=akz(kz)+bkz(kz)*ps(ix,jy,1,n)  ! pressure on model layers
        tv=tth(ix,jy,kz,n)*(1.+0.608*qvh(ix,jy,kz,n))

        if (abs(tv-tvold).gt.0.2) then
          zlev(kz)=zold+const*log(pold/pint)*(tv-tvold)/log(tv/tvold)
        else
          zlev(kz)=zold+const*log(pold/pint)*tv
        endif
        tvold=tv
        pold=pint
        zold=zlev(kz)
      end do

  ! 2) Define a minimum level kzmin, from which upward the tropopause is
  !    searched for. This is to avoid inversions in the lower troposphere
  !    to be identified as the tropopause
  !************************************************************************

      if (metdata_format.eq.GRIBFILE_CENTRE_ECMWF) then
        !LB, The CTM version has 2 (as bugfix), so I changed it 2 from 1 to try out
        loop_start=2
      else
        loop_start=llev
      end if

      do kz=loop_start,nuvz
        if (zlev(kz).ge.altmin) then
          kzmin=kz
          exit
        endif
      end do

  ! 3) Search for first stable layer above minimum height that fulfills the
  !    thermal tropopause criterion
  !************************************************************************

      do kz=kzmin,nuvz
        do lz=kz+1,nuvz
          if ((zlev(lz)-zlev(kz)).gt.2000.) then
            if (((tth(ix,jy,kz,n)-tth(ix,jy,lz,n))/ &
                 (zlev(lz)-zlev(kz))).lt.0.002) then
              tropopause(ix,jy,1,n)=zlev(kz)
              goto 51
            endif
            goto 50
          endif
        end do
50      continue
      end do
51    continue


    end do
  end do
  ! !$OMP END DO
  ! !$OMP END PARALLEL
  ! openmp change end

  ! Calculation of potential vorticity on 3-d grid
  !***********************************************

  call calcpv(n,uuh,vvh,pvh)
end subroutine calcpar

subroutine calcpar_nests(n,uuhn,vvhn,pvhn,metdata_format)
  !                         i  i    i    o
  !*****************************************************************************
  !                                                                            *
  !     Computation of several boundary layer parameters needed for the        *
  !     dispersion calculation and calculation of dry deposition velocities.   *
  !     All parameters are calculated over the entire grid.                    *
  !     This routine is similar to calcpar, but is used for the nested grids.  *
  !                                                                            *
  !     Author: A. Stohl                                                       *
  !                                                                            *
  !     8 February 1999                                                        *
  !                                                                            *
  ! ------------------------------------------------------------------         *
  !     Petra Seibert, Feb 2000:                                               *
  !     convection scheme:                                                     *
  !     new variables in call to richardson                                    *
  !                                                                            *
  !*****************************************************************************
  !  Changes, Bernd C. Krueger, Feb. 2001:                                     *
  !   Variables tth and qvh (on eta coordinates) in common block               *
  !                                                                            *
  !   Unified ECMWF and GFS builds                                             *
  !   Marian Harustak, 12.5.2017                                               *
  !     - Added passing of metdata_format as it was needed by called routines  *
  !*****************************************************************************
  !                                                                            *
  ! Variables:                                                                 *
  ! n                  temporal index for meteorological fields (1 to 3)       *
  ! metdata_format     format of metdata (ecmwf/gfs)                           *
  !                                                                            *
  ! Constants:                                                                 *
  !                                                                            *
  !                                                                            *
  ! Functions:                                                                 *
  ! scalev             computation of ustar                                    *
  ! obukhov            computatio of Obukhov length                            *
  !                                                                            *
  !*****************************************************************************

  use par_mod
  use com_mod
  use drydepo_mod
  use qvsat_mod

  implicit none

  integer :: metdata_format
  integer :: n,ix,jy,i,l,kz,lz,kzmin
  real :: ttlev(nuvzmax),qvlev(nuvzmax),ol,hmixplus,dummyakzllev
  real :: ulev(nuvzmax),vlev(nuvzmax),rh,vd(maxspec),subsceff,ylat
  real :: altmin,tvold,pold,zold,pint,tv,zlev(nuvzmax)
  real :: uuhn(0:nxmaxn-1,0:nymaxn-1,nuvzmax,maxnests)
  real :: vvhn(0:nxmaxn-1,0:nymaxn-1,nuvzmax,maxnests)
  real :: pvhn(0:nxmaxn-1,0:nymaxn-1,nuvzmax,maxnests)
  real,parameter :: const=r_air/ga


  ! Loop over all nests
  !********************

  do l=1,numbnests

  ! Loop over entire grid
  !**********************

  do jy=0,nyn(l)-1

  ! Set minimum height for tropopause
  !**********************************

    ylat=ylat0n(l)+real(jy)*dyn(l)
    if ((ylat.ge.-20.).and.(ylat.le.20.)) then
      altmin = 5000.
    else
      if ((ylat.gt.20.).and.(ylat.lt.40.)) then
        altmin=2500.+(40.-ylat)*125.
      else if ((ylat.gt.-40.).and.(ylat.lt.-20.)) then
        altmin=2500.+(40.+ylat)*125.
      else
        altmin=2500.
      endif
    endif

    do ix=0,nxn(l)-1

  ! 1) Calculation of friction velocity
  !************************************

      ustarn(ix,jy,1,n,l)=scalev(psn(ix,jy,1,n,l),tt2n(ix,jy,1,n,l), &
           td2n(ix,jy,1,n,l),surfstrn(ix,jy,1,n,l))
      if (ustarn(ix,jy,1,n,l).le.1.e-8) ustarn(ix,jy,1,n,l)=1.e-8

  ! 2) Calculation of inverse Obukhov length scale
  !***********************************************

      ol=obukhov(psn(ix,jy,1,n,l),tt2n(ix,jy,1,n,l), &
           td2n(ix,jy,1,n,l),tthn(ix,jy,2,n,l),ustarn(ix,jy,1,n,l), &
           sshfn(ix,jy,1,n,l),akm,bkm,dummyakzllev,metdata_format)
      if (ol.ne.0.) then
        olin(ix,jy,1,n,l)=1./ol
      else
        olin(ix,jy,1,n,l)=99999.
      endif


  ! 3) Calculation of convective velocity scale and mixing height
  !**************************************************************

      do i=1,nuvz
        ulev(i)=uuhn(ix,jy,i,l)
        vlev(i)=vvhn(ix,jy,i,l)
        ttlev(i)=tthn(ix,jy,i,n,l)
        qvlev(i)=qvhn(ix,jy,i,n,l)
      end do

      call richardson(psn(ix,jy,1,n,l),ustarn(ix,jy,1,n,l),ttlev, &
           qvlev,ulev,vlev,nuvz,akz,bkz,sshfn(ix,jy,1,n,l), &
           tt2n(ix,jy,1,n,l),td2n(ix,jy,1,n,l),hmixn(ix,jy,1,n,l), &
           wstarn(ix,jy,1,n,l),hmixplus,metdata_format)

      if(lsubgrid.eq.1) then
        subsceff=min(excessoron(ix,jy,l),hmixplus)
      else
        subsceff=0.0
      endif
  !
  ! CALCULATE HMIX EXCESS ACCORDING TO SUBGRIDSCALE VARIABILITY AND STABILITY
  !
      hmixn(ix,jy,1,n,l)=hmixn(ix,jy,1,n,l)+subsceff
      hmixn(ix,jy,1,n,l)=max(hmixmin,hmixn(ix,jy,1,n,l)) ! minim PBL height
      hmixn(ix,jy,1,n,l)=min(hmixmax,hmixn(ix,jy,1,n,l)) ! maxim PBL height


  ! 4) Calculation of dry deposition velocities
  !********************************************

      if (DRYDEP) then
        ! z0(4)=0.016*ustarn(ix,jy,1,n,l)*ustarn(ix,jy,1,n,l)/ga
        ! z0(9)=0.016*ustarn(ix,jy,1,n,l)*ustarn(ix,jy,1,n,l)/ga
        z0(7)=0.016*ustarn(ix,jy,1,n,l)*ustarn(ix,jy,1,n,l)/ga

  ! Calculate relative humidity at surface
  !***************************************
        rh=ew(td2n(ix,jy,1,n,l),psn(ix,jy,1,n,l))/ew(tt2n(ix,jy,1,n,l),psn(ix,jy,1,n,l))

        call getvdep_nests(n,ix,jy,ustarn(ix,jy,1,n,l), &
             tt2n(ix,jy,1,n,l),psn(ix,jy,1,n,l),1./olin(ix,jy,1,n,l), &
             ssrn(ix,jy,1,n,l),rh,lsprecn(ix,jy,1,n,l)+ &
             convprecn(ix,jy,1,n,l),sdn(ix,jy,1,n,l),vd,l)

        do i=1,nspec
          vdepn(ix,jy,i,n,l)=vd(i)
        end do

      endif

  !******************************************************
  ! Calculate height of thermal tropopause (Hoinka, 1997)
  !******************************************************

  ! 1) Calculate altitudes of ECMWF model levels
  !*********************************************

      tvold=tt2n(ix,jy,1,n,l)*(1.+0.378*ew(td2n(ix,jy,1,n,l),psn(ix,jy,1,n,l))/ &
           psn(ix,jy,1,n,l))
      pold=psn(ix,jy,1,n,l)
      zold=0.
      do kz=2,nuvz
        pint=akz(kz)+bkz(kz)*psn(ix,jy,1,n,l)  ! pressure on model layers
        tv=tthn(ix,jy,kz,n,l)*(1.+0.608*qvhn(ix,jy,kz,n,l))

        if (abs(tv-tvold).gt.0.2) then
         zlev(kz)=zold+const*log(pold/pint)*(tv-tvold)/log(tv/tvold)
        else
          zlev(kz)=zold+const*log(pold/pint)*tv
        endif
        tvold=tv
        pold=pint
        zold=zlev(kz)
      end do

  ! 2) Define a minimum level kzmin, from which upward the tropopause is
  !    searched for. This is to avoid inversions in the lower troposphere
  !    to be identified as the tropopause
  !************************************************************************

      do kz=1,nuvz
        if (zlev(kz).ge.altmin) then
          kzmin=kz
          exit
        endif
      end do

  ! 3) Search for first stable layer above minimum height that fulfills the
  !    thermal tropopause criterion
  !************************************************************************

      do kz=kzmin,nuvz
        do lz=kz+1,nuvz
          if ((zlev(lz)-zlev(kz)).gt.2000.) then
            if (((tthn(ix,jy,kz,n,l)-tthn(ix,jy,lz,n,l))/ &
                 (zlev(lz)-zlev(kz))).lt.0.002) then
              tropopausen(ix,jy,1,n,l)=zlev(kz)
              goto 51
            endif
            goto 50
          endif
        end do
50      continue
      end do
51    continue


    end do
  end do

  ! Calculation of potential vorticity on 3-d grid
  !***********************************************

  call calcpv_nests(l,n,uuhn,vvhn,pvhn)

  end do
end subroutine calcpar_nests

real function obukhov(ps,tsurf,tdsurf,tlev,ustar,hf,akm,bkm,plev,metdata_format)

  !********************************************************************
  !                                                                   *
  !                       Author: G. WOTAWA                           *
  !                       Date:   1994-06-27                          *
  !                                                                   *
  !     This program calculates Obukhov scale height from surface     *
  !     meteorological data and sensible heat flux.                   *
  !                                                                   *
  !********************************************************************
  !                                                                   *
  !  Update: A. Stohl, 2000-09-25, avoid division by zero by          *
  !  setting ustar to minimum value                                   *
  !  CHANGE: 17/11/2005 Caroline Forster NCEP GFS version             *
  !                                                                   *
  !   Unified ECMWF and GFS builds                                    *
  !   Marian Harustak, 12.5.2017                                      *
  !     - Merged obukhov and obukhov_gfs into one routine using       *
  !       if-then for meteo-type dependent code                       *
  !                                                                   *
  !********************************************************************
  !                                                                   *
  !     INPUT:                                                        *
  !                                                                   *
  !     ps      surface pressure [Pa]                                 *
  !     tsurf   surface temperature [K]                               *
  !     tdsurf  surface dew point [K]                                 *
  !     tlev    temperature first model level [K]                     *
  !     ustar   scale velocity [m/s]                                  *
  !     hf      surface sensible heat flux [W/m2]                     *
  !     akm     ECMWF vertical discretization parameter               *
  !     bkm     ECMWF vertical discretization parameter               *
  !     plev                                                          *
  !     metdata_format format of metdata (ecmwf/gfs)                  *
  !                                                                   *
  !********************************************************************

  use par_mod
  use class_gribfile
  use qvsat_mod

  implicit none

  integer :: metdata_format
  real :: akm(nwzmax),bkm(nwzmax)
  real :: ps,tsurf,tdsurf,tlev,ustar,hf,e,tv,rhoa,plev
  real :: ak1,bk1,theta,thetastar


  e=ew(tdsurf,ps)                           ! vapor pressure
  tv=tsurf*(1.+0.378*e/ps)               ! virtual temperature
  rhoa=ps/(r_air*tv)                      ! air density
  if (metdata_format.eq.GRIBFILE_CENTRE_ECMWF) then
  ak1=(akm(1)+akm(2))/2.
  bk1=(bkm(1)+bkm(2))/2.
  plev=ak1+bk1*ps                        ! Pressure level 1
  end if
  theta=tlev*(100000./plev)**(r_air/cpa) ! potential temperature
  if (ustar.le.0.) ustar=1.e-8
  thetastar=hf/(rhoa*cpa*ustar)           ! scale temperature
  if(abs(thetastar).gt.1.e-10) then
     obukhov=theta*ustar**2/(karman*ga*thetastar)
  else
     obukhov=9999                        ! zero heat flux
  endif
  if (obukhov.gt. 9999.) obukhov= 9999.
  if (obukhov.lt.-9999.) obukhov=-9999.
end function obukhov

subroutine richardson(psurf,ust,ttlev,qvlev,ulev,vlev,nuvz, &
       akz,bkz,hf,tt2,td2,h,wst,hmixplus,metdata_format)
  !                        i    i    i     i    i    i    i
  ! i   i  i   i   i  o  o     o
  !****************************************************************************
  !                                                                           *
  !     Calculation of mixing height based on the critical Richardson number. *
  !     Calculation of convective time scale.                                 *
  !     For unstable conditions, one iteration is performed. An excess        *
  !     temperature (dependent on hf and wst) is calculated, added to the     *
  !     temperature at the lowest model level. Then the procedure is repeated.*
  !                                                                           *
  !     Author: A. Stohl                                                      *
  !                                                                           *
  !     22 August 1996                                                        *
  !                                                                           *
  !     Literature:                                                           *
  !     Vogelezang DHP and Holtslag AAM (1996): Evaluation and model impacts  *
  !     of alternative boundary-layer height formulations. Boundary-Layer     *
  !     Meteor. 81, 245-269.                                                  *
  !                                                                           *
  !****************************************************************************
  !                                                                           *
  !     Update: 1999-02-01 by G. Wotawa                                       *
  !                                                                           *
  !     Two meter level (temperature, humidity) is taken as reference level   *
  !     instead of first model level.                                         *
  !     New input variables tt2, td2 introduced.                              *
  !                                                                           *
  !     CHANGE: 17/11/2005 Caroline Forster NCEP GFS version                  * 
  !                                                                           *
  !     Unified ECMWF and GFS builds                                          *
  !     Marian Harustak, 12.5.2017                                            *
  !       - Merged richardson and richardson_gfs into one routine using       *
  !         if-then for meteo-type dependent code                             *
  !                                                                           *
  !****************************************************************************
  !                                                                           *
  ! Variables:                                                                *
  ! h                          mixing height [m]                              *
  ! hf                         sensible heat flux                             *
  ! psurf                      surface pressure at point (xt,yt) [Pa]         *
  ! tv                         virtual temperature                            *
  ! wst                        convective velocity scale                      *
  ! metdata_format             format of metdata (ecmwf/gfs)                  *
  !                                                                           *
  ! Constants:                                                                *
  ! ric                        critical Richardson number                     *
  !                                                                           *
  !****************************************************************************

  use par_mod
  use class_gribfile
  use qvsat_mod

  implicit none

  integer :: metdata_format
  integer :: i,k,nuvz,iter,llev,loop_start
  real :: tv,tvold,zref,z,zold,pint,pold,theta,thetaref,ri
  real :: akz(nuvz),bkz(nuvz),ulev(nuvz),vlev(nuvz),hf,wst,tt2,td2
  real :: psurf,ust,ttlev(nuvz),qvlev(nuvz),h,excess
  real :: thetaold,zl,ul,vl,thetal,ril,hmixplus,wspeed,bvfsq,bvf
  real :: rh,rhold,rhl,theta1,theta2,zl1,zl2,thetam
  real,parameter    :: const=r_air/ga, ric=0.25, b=100., bs=8.5
  integer,parameter :: itmax=3

  excess=0.0
  iter=0

  if (metdata_format.eq.GRIBFILE_CENTRE_NCEP) then
    ! NCEP version: find first model level above ground
    !**************************************************

     llev = 0
     do i=1,nuvz
       if (psurf.lt.akz(i)) llev=i
     end do
     llev = llev+1
    ! sec llev should not be 1!
     if (llev.eq.1) llev = 2
     if (llev.gt.nuvz) llev = nuvz-1
    ! NCEP version
  end if


  ! Compute virtual temperature and virtual potential temperature at
  ! reference level (2 m)
  !*****************************************************************

30   iter=iter+1

  pold=psurf
  tvold=tt2*(1.+0.378*ew(td2,psurf)/psurf)
  zold=2.0
  zref=zold
  rhold=ew(td2,psurf)/ew(tt2,psurf)

  thetaref=tvold*(100000./pold)**(r_air/cpa)+excess
  thetaold=thetaref


  ! Integrate z up to one level above zt
  !*************************************
  if (metdata_format.eq.GRIBFILE_CENTRE_ECMWF) then
    loop_start=2
  else
    loop_start=llev
  end if
  do k=loop_start,nuvz
    pint=akz(k)+bkz(k)*psurf  ! pressure on model layers
    tv=ttlev(k)*(1.+0.608*qvlev(k))

    if (abs(tv-tvold).gt.0.2) then
      z=zold+const*log(pold/pint)*(tv-tvold)/log(tv/tvold)
    else
      z=zold+const*log(pold/pint)*tv
    endif

    theta=tv*(100000./pint)**(r_air/cpa)
  ! Petra
    rh = qvlev(k) / f_qvsat( pint, ttlev(k) )


  ! Calculate Richardson number at each level
  !****************************************

    ri=ga/thetaref*(theta-thetaref)*(z-zref)/ &
         max(((ulev(k)-ulev(2))**2+(vlev(k)-vlev(2))**2+b*ust**2),0.1)

  !  addition of second condition: MH should not be placed in an
  !  unstable layer (PS / Feb 2000)
    if (ri.gt.ric .and. thetaold.lt.theta) goto 20

    tvold=tv
    pold=pint
    rhold=rh
    thetaold=theta
    zold=z
  end do
  k=k-1 ! ESO: make sure k <= nuvz (ticket #139)
20   continue

  ! Determine Richardson number between the critical levels
  !********************************************************

  zl1=zold
  theta1=thetaold
  do i=1,20
    zl=zold+real(i)/20.*(z-zold)
    ul=ulev(k-1)+real(i)/20.*(ulev(k)-ulev(k-1))
    vl=vlev(k-1)+real(i)/20.*(vlev(k)-vlev(k-1))
    thetal=thetaold+real(i)/20.*(theta-thetaold)
    rhl=rhold+real(i)/20.*(rh-rhold)
    ril=ga/thetaref*(thetal-thetaref)*(zl-zref)/ &
         max(((ul-ulev(2))**2+(vl-vlev(2))**2+b*ust**2),0.1)
    zl2=zl
    theta2=thetal
    if (ril.gt.ric) exit
    zl1=zl
    theta1=thetal
  end do

  h=zl
  thetam=0.5*(theta1+theta2)
  wspeed=sqrt(ul**2+vl**2)                    ! Wind speed at z=hmix
  bvfsq=(ga/thetam)*(theta2-theta1)/(zl2-zl1) ! Brunt-Vaisala frequency
                                              ! at z=hmix

  ! Under stable conditions, limit the maximum effect of the subgrid-scale topography
  ! by the maximum lifting possible from the available kinetic energy
  !*****************************************************************************

  if(bvfsq.le.0.) then
    hmixplus=9999.
  else
    bvf=sqrt(bvfsq)
    hmixplus=wspeed/bvf*convke
  endif


  ! Calculate convective velocity scale
  !************************************

  if (hf.lt.0.) then
    wst=(-h*ga/thetaref*hf/cpa)**0.333
    excess=-bs*hf/cpa/wst
    if (iter.lt.itmax) goto 30
  else
    wst=0.
  endif

end subroutine richardson

real function scalev(ps,t,td,stress)

  !********************************************************************
  !                                                                   *
  !                       Author: G. WOTAWA                           *
  !                       Date:   1994-06-27                          *
  !                       Update: 1996-05-21 A. Stohl                 *
  !                                                                   *
  !********************************************************************
  !                                                                   *
  !     This Programm calculates scale velocity ustar from surface    *
  !     stress and air density.                                       *
  !                                                                   *
  !********************************************************************
  !                                                                   *
  !     INPUT:                                                        *
  !                                                                   *
  !     ps      surface pressure [Pa]                                 *
  !     t       surface temperature [K]                               *
  !     td      surface dew point [K]                                 *
  !     stress  surface stress [N/m2]                                 *
  !                                                                   *
  !********************************************************************
  
  use par_mod
  use qvsat_mod

  implicit none

  real :: ps,t,td,e,tv,rhoa,stress

  e=ew(td,ps)                       ! vapor pressure
  tv=t*(1.+0.378*e/ps)           ! virtual temperature
  rhoa=ps/(r_air*tv)              ! air density
  scalev=sqrt(abs(stress)/rhoa)

end function scalev

end module getfields_mod