! SPDX-FileCopyrightText: FLEXPART 1998-2019, see flexpart_license.txt
! SPDX-License-Identifier: GPL-3.0-or-later

module verttransform_mod
  use par_mod
  use com_mod
  use windfields_mod

  implicit none

contains

subroutine verttransform_ecmwf(n,uuh,vvh,wwh,pvh)
  !                              i  i   i   i   i
  !*****************************************************************************
  !                                                                            *
  !     This subroutine transforms temperature, dew point temperature and      *
  !     wind components from eta to meter coordinates.                         *
  !     The vertical wind component is transformed from Pa/s to m/s using      *
  !     the conversion factor pinmconv.                                        *
  !     In addition, this routine calculates vertical density gradients        *
  !     needed for the parameterization of the turbulent velocities.           *
  !                                                                            *
  !     Author: A. Stohl, G. Wotawa                                            *
  !                                                                            *
  !     12 August 1996                                                         *
  !     Update: 16 January 1998                                                *
  !                                                                            *
  !     Major update: 17 February 1999                                         *
  !     by G. Wotawa                                                           *
  !                                                                            *
  !     - Vertical levels for u, v and w are put together                      *
  !     - Slope correction for vertical velocity: Modification of calculation  *
  !       procedure                                                            *
  !                                                                            *
  !*****************************************************************************
  !  Changes, Bernd C. Krueger, Feb. 2001:
  !   Variables tth and qvh (on eta coordinates) from common block
  !
  ! Sabine Eckhardt, March 2007
  ! added the variable cloud for use with scavenging - descr. in com_mod
  !
  ! Unified ECMWF and GFS builds
  ! Marian Harustak, 12.5.2017 
  !     - Renamed from verttransform to verttransform_ecmwf
  !*****************************************************************************
  ! Date: 2017-05-30 modification of a bug in ew. Don Morton (CTBTO project)   *
  !*****************************************************************************
  !                                                                            *
  ! Variables:                                                                 *
  ! nx,ny,nz                        field dimensions in x,y and z direction    *
  ! clouds(0:nxmax,0:nymax,0:nzmax,numwfmem) cloud field for wet deposition    *
  ! uu(0:nxmax,0:nymax,nzmax,numwfmem)     wind components in x-direction [m/s]*
  ! vv(0:nxmax,0:nymax,nzmax,numwfmem)     wind components in y-direction [m/s]*
  ! ww(0:nxmax,0:nymax,nzmax,numwfmem)     wind components in z-direction      *
  !                                          [deltaeta/s]                      *
  ! tt(0:nxmax,0:nymax,nzmax,numwfmem)     temperature [K]                     *
  ! pv(0:nxmax,0:nymax,nzmax,numwfmem)     potential voriticity (pvu)          *
  ! ps(0:nxmax,0:nymax,numwfmem)           surface pressure [Pa]               *
  !                                                                            *
  !*****************************************************************************

  use par_mod
  use com_mod
  use cmapf_mod, only: cc2gll
  use qvsat_mod

  implicit none

  integer, intent(in) :: n
  real,intent(in),dimension(0:nxmax-1,0:nymax-1,nuvzmax) :: uuh,vvh,pvh
  real,intent(in),dimension(0:nxmax-1,0:nymax-1,nwzmax) :: wwh

  real,dimension(0:nxmax-1,0:nymax-1,nuvzmax) :: rhoh
  real,dimension(0:nxmax-1,0:nymax-1,nzmax) :: pinmconv
  ! RLT added pressure
  real,dimension(0:nxmax-1,0:nymax-1,nuvzmax) :: prsh

  logical :: init = .true.

  !*************************************************************************
  ! If verttransform is called the first time, initialize heights of the   *
  ! z levels in meter. The heights are the heights of model levels, where  *
  ! u,v,T and qv are given, and of the interfaces, where w is given. So,   *
  ! the vertical resolution in the z system is doubled. As reference point,*
  ! the lower left corner of the grid is used.                             *
  ! Unlike in the eta system, no difference between heights for u,v and    *
  ! heights for w exists.                                                  *
  !*************************************************************************


  !eso measure CPU time
  !  call mpif_mtime('verttransform',0)

  if (init) then

  ! Search for a point with high surface pressure (i.e. not above significant topography)
  ! Then, use this point to construct a reference z profile, to be used at all times
  !*****************************************************************************
    call initialise_ecmwf_verttransform(n)

  ! Do not repeat initialization of the Cartesian z grid
  !*****************************************************

    init=.false.
  endif


  ! Compute heights of eta levels and their respective pressure and density fields
  !*******************************************************************************
  call verttransform_ecmwf_heights(n,rhoh,pinmconv,prsh)

  ! Transform the wind fields to the internal coordinate system and save the native ETA 
  ! fields when case wind_coord_type==ETA
  !*************************************************************
  call verttransform_ecmwf_transform_windfields(n,uuh,vvh,wwh,pvh,rhoh,prsh,pinmconv)

  ! If north or south pole is in the domain, calculate wind velocities in polar
  ! stereographic coordinates
  !*******************************************************************
  call verttransform_ecmwf_stereographic(n)

  ! Create cloud fields
  !*********************
  call verttransform_ecmwf_clouds(n,readclouds,sumclouds,nxmin1,nymin1,clouds(:,:,:,n), &
    cloudsh(:,:,n),clw(:,:,:,n),ctwc(:,:,n),clwc(:,:,:,n),ciwc(:,:,:,n),lsprec(:,:,1,n), &
    convprec(:,:,1,n),rho(:,:,:,n),tt(:,:,:,n),qv(:,:,:,n))
end subroutine verttransform_ecmwf

subroutine initialise_ecmwf_verttransform(n)
  use par_mod
  use com_mod
  use qvsat_mod
  use initialise_mod
  use output_mod

  implicit none

  integer, intent(in) :: n
  real,dimension(0:nxmax-1,0:nymax-1) ::  tvold,pold,pint,tv
  integer :: ix,jy,kz,ixm,jym
  real,parameter :: const=r_air/ga

  if (ipin.eq.1) then
    call read_heightlevels(height,nmixz)
    return
  endif

  loop1: do jy=0,nymin1
    do ix=0,nxmin1
      if (ps(ix,jy,1,n).gt.100000.) then
        ixm=ix
        jym=jy
        exit loop1
      endif
    end do
  end do loop1

  tvold(ixm,jym)=tt2(ixm,jym,1,n)*(1.+0.378*ew(td2(ixm,jym,1,n),ps(ixm,jym,1,n))/ &
       ps(ixm,jym,1,n))
  pold(ixm,jym)=ps(ixm,jym,1,n)
  height(1)=0.

  do kz=2,nuvz
    pint=akz(kz)+bkz(kz)*ps(ixm,jym,1,n)
    tv=tth(ixm,jym,kz,n)*(1.+0.608*qvh(ixm,jym,kz,n))

    if (abs(tv(ixm,jym)-tvold(ixm,jym)).gt.0.2) then
      height(kz)= &
           height(kz-1)+const*log(pold(ixm,jym)/pint(ixm,jym))* &
           (tv(ixm,jym)-tvold(ixm,jym))/log(tv(ixm,jym)/tvold(ixm,jym))
    else
      height(kz)=height(kz-1)+ &
           const*log(pold(ixm,jym)/pint(ixm,jym))*tv(ixm,jym)
    endif

    tvold(ixm,jym)=tv(ixm,jym)
    pold(ixm,jym)=pint(ixm,jym)
  end do

  ! Determine highest levels that can be within PBL
  !************************************************

  do kz=1,nz
    if (height(kz).gt.hmixmax) then
      nmixz=kz
      exit
    endif
  end do

  call output_heightlevels(height,nmixz)
end subroutine initialise_ecmwf_verttransform

subroutine verttransform_ecmwf_heights(n,rhoh,pinmconv,prsh)
  use par_mod
  use com_mod
  use qvsat_mod

  implicit none

  integer, intent(in) :: n
  real,intent(inout),dimension(0:nxmax-1,0:nymax-1,nuvzmax) :: rhoh
  real,intent(inout),dimension(0:nxmax-1,0:nymax-1,nzmax) :: pinmconv
  ! RLT added pressure
  real,intent(inout),dimension(0:nxmax-1,0:nymax-1,nuvzmax) :: prsh
  real,dimension(0:nxmax-1,0:nymax-1) :: tvold,pold,pint,tv

  integer :: ix,jy,kz,iz
  real,parameter :: const=r_air/ga


  ! Loop over the whole grid
  !*************************

  do jy=0,nymin1
    do ix=0,nxmin1
      tvold(ix,jy)=tt2(ix,jy,1,n)*(1.+0.378*ew(td2(ix,jy,1,n),ps(ix,jy,1,n))/ &
           ps(ix,jy,1,n))
    enddo
  enddo
  pold=ps(:,:,1,n)
  etauvheight(:,:,1,n)=0.
  etawheight(:,:,1,n)=0.
  rhoh(:,:,1)=pold/(r_air*tvold)
  ! RLT add pressure
  prsh(:,:,1)=ps(:,:,1,n)

  ! Compute heights in meter of the half-model levels over the whole grid
  !*****************************************************************
  do kz=2,nuvz
    pint=akz(kz)+bkz(kz)*ps(:,:,1,n)
    ! RLT add pressure
    prsh(:,:,kz)=pint
    tv=tth(:,:,kz,n)*(1.+0.608*qvh(:,:,kz,n))
    rhoh(:,:,kz)=pint(:,:)/(r_air*tv)

    where (abs(tv-tvold).gt.0.2) 
      etauvheight(:,:,kz,n)=etauvheight(:,:,kz-1,n)+const*log(pold/pint)* &
           (tv-tvold)/log(tv/tvold)
    elsewhere
      etauvheight(:,:,kz,n)=etauvheight(:,:,kz-1,n)+const*log(pold/pint)*tv
    endwhere

    tvold=tv
    pold=pint
  end do

  ! Compute heights in meter of the model levels over the whole grid
  !*****************************************************************
  do kz=2,nwz-1
    etawheight(:,:,kz,n)=(etauvheight(:,:,kz+1,n)+etauvheight(:,:,kz,n))/2.
  end do
  etawheight(:,:,nwz,n)=etawheight(:,:,nwz-1,n)+ &
       etauvheight(:,:,nuvz,n)-etauvheight(:,:,nuvz-1,n)

  ! pinmconv=(h2-h1)/(p2-p1)

  pinmconv(:,:,1)=(etauvheight(:,:,2,n))/ &
       ((aknew(2)+bknew(2)*ps(:,:,1,n))- &
       (aknew(1)+bknew(1)*ps(:,:,1,n)))
  do kz=2,nz-1
    pinmconv(:,:,kz)=(etauvheight(:,:,kz+1,n)-etauvheight(:,:,kz-1,n))/ &
         ((aknew(kz+1)+bknew(kz+1)*ps(:,:,1,n))- &
         (aknew(kz-1)+bknew(kz-1)*ps(:,:,1,n)))
  end do
  pinmconv(:,:,nz)=(etauvheight(:,:,nz,n)-etauvheight(:,:,nz-1,n))/ &
       ((aknew(nz)+bknew(nz)*ps(:,:,1,n))- &
       (aknew(nz-1)+bknew(nz-1)*ps(:,:,1,n)))
end subroutine verttransform_ecmwf_heights

subroutine verttransform_ecmwf_transform_windfields(n,uuh,vvh,wwh,pvh,rhoh,prsh,pinmconv)
  use par_mod
  use com_mod
  use qvsat_mod

  implicit none

  integer,intent(in) :: n
  real,intent(in),dimension(0:nxmax-1,0:nymax-1,nuvzmax) :: uuh,vvh,pvh
  real,intent(in),dimension(0:nxmax-1,0:nymax-1,nwzmax) :: wwh
  real,intent(in),dimension(0:nxmax-1,0:nymax-1,nuvzmax) :: rhoh
  real,intent(in),dimension(0:nxmax-1,0:nymax-1,nzmax) :: pinmconv
  ! RLT added pressure
  real,intent(in),dimension(0:nxmax-1,0:nymax-1,nuvzmax) :: prsh

  real,dimension(0:nxmax-1,0:nymax-1) ::  dpdeta

  real,dimension(0:nymax-1) :: cosf

  integer,dimension(0:nxmax-1,0:nymax-1) :: idx

  integer :: ix,jy,kz,iz,kmin,ixp,jyp,ix1,jy1
  real :: dz1,dz2,dz
  real :: xlon,ylat,xlonr,dzdx,dzdy
  real :: dzdx1,dzdx2,dzdy1,dzdy2

  ! Levels, where u,v,t and q are given
  !************************************
  uu(:,:,1,n)=uuh(:,:,1)
  vv(:,:,1,n)=vvh(:,:,1)
  tt(:,:,1,n)=tth(:,:,1,n)
  qv(:,:,1,n)=qvh(:,:,1,n)
  !hg adding the cloud water 
  if (readclouds) then
    clwc(:,:,1,n)=clwch(:,:,1,n)
    if (.not.sumclouds) ciwc(:,:,1,n)=ciwch(:,:,1,n)
  end if
  !hg 
  pv(:,:,1,n)=pvh(:,:,1)
  rho(:,:,1,n)=rhoh(:,:,1)
  ! RLT add pressure
  prs(:,:,1,n)=prsh(:,:,1)

  uu(:,:,nz,n)=uuh(:,:,nuvz)
  vv(:,:,nz,n)=vvh(:,:,nuvz)
  tt(:,:,nz,n)=tth(:,:,nuvz,n)
  qv(:,:,nz,n)=qvh(:,:,nuvz,n)
  !hg adding the cloud water
  if (readclouds) then
    clwc(:,:,nz,n)=clwch(:,:,nuvz,n)
    if (.not.sumclouds) ciwc(:,:,nz,n)=ciwch(:,:,nuvz,n)
  end if
  !hg
  pv(:,:,nz,n)=pvh(:,:,nuvz)
  rho(:,:,nz,n)=rhoh(:,:,nuvz)
  ! RLT
  prs(:,:,nz,n)=prsh(:,:,nuvz)

  kmin=2
  idx=kmin
  do iz=2,nz-1
    do jy=0,nymin1
      do ix=0,nxmin1
        if(height(iz).gt.etauvheight(ix,jy,nuvz,n)) then
          uu(ix,jy,iz,n)=uu(ix,jy,nz,n)
          vv(ix,jy,iz,n)=vv(ix,jy,nz,n)
          tt(ix,jy,iz,n)=tt(ix,jy,nz,n)
          qv(ix,jy,iz,n)=qv(ix,jy,nz,n)
  !hg adding the cloud water
          if (readclouds) then
            clwc(ix,jy,iz,n)=clwc(ix,jy,nz,n)
            if (.not.sumclouds) ciwc(ix,jy,iz,n)=ciwc(ix,jy,nz,n)
          end if
  !hg
          pv(ix,jy,iz,n)=pv(ix,jy,nz,n)
          rho(ix,jy,iz,n)=rho(ix,jy,nz,n)
  ! RLT
          prs(ix,jy,iz,n)=prs(ix,jy,nz,n)
        else
          innuvz: do kz=idx(ix,jy),nuvz
            if (idx(ix,jy) .le. kz .and. (height(iz).gt.etauvheight(ix,jy,kz-1,n)).and. &
                 (height(iz).le.etauvheight(ix,jy,kz,n))) then
              idx(ix,jy)=kz
              exit innuvz
            endif
          enddo innuvz
        endif
      enddo
    enddo
    do jy=0,nymin1
      do ix=0,nxmin1
        if(height(iz).le.etauvheight(ix,jy,nuvz,n)) then
          kz=idx(ix,jy)
          dz1=height(iz)-etauvheight(ix,jy,kz-1,n)
          dz2=etauvheight(ix,jy,kz,n)-height(iz)
          dz=dz1+dz2
          uu(ix,jy,iz,n)=(uuh(ix,jy,kz-1)*dz2+uuh(ix,jy,kz)*dz1)/dz
          vv(ix,jy,iz,n)=(vvh(ix,jy,kz-1)*dz2+vvh(ix,jy,kz)*dz1)/dz
          tt(ix,jy,iz,n)=(tth(ix,jy,kz-1,n)*dz2 &
               +tth(ix,jy,kz,n)*dz1)/dz
          qv(ix,jy,iz,n)=(qvh(ix,jy,kz-1,n)*dz2 &
               +qvh(ix,jy,kz,n)*dz1)/dz
  !hg adding the cloud water
          if (readclouds) then
            clwc(ix,jy,iz,n)=(clwch(ix,jy,kz-1,n)*dz2+clwch(ix,jy,kz,n)*dz1)/dz
            if (.not.sumclouds) &
                 &ciwc(ix,jy,iz,n)=(ciwch(ix,jy,kz-1,n)*dz2+ciwch(ix,jy,kz,n)*dz1)/dz
          end if
  !hg
          pv(ix,jy,iz,n)=(pvh(ix,jy,kz-1)*dz2+pvh(ix,jy,kz)*dz1)/dz
          rho(ix,jy,iz,n)=(rhoh(ix,jy,kz-1)*dz2+rhoh(ix,jy,kz)*dz1)/dz
  ! RLT add pressure
          prs(ix,jy,iz,n)=(prsh(ix,jy,kz-1)*dz2+prsh(ix,jy,kz)*dz1)/dz
        endif
      enddo
    enddo
  enddo


  ! Levels, where w is given
  !*************************

  ww(:,:,1,n)=wwh(:,:,1)*pinmconv(:,:,1)
  ww(:,:,nz,n)=wwh(:,:,nwz)*pinmconv(:,:,nz)
  kmin=2
  idx=kmin
  do iz=2,nz
    do jy=0,nymin1
      do ix=0,nxmin1
        inn:         do kz=idx(ix,jy),nwz
          if(idx(ix,jy) .le. kz .and. height(iz).gt.etawheight(ix,jy,kz-1,n).and. &
               height(iz).le.etawheight(ix,jy,kz,n)) then
            idx(ix,jy)=kz
            exit inn
          endif
        enddo inn
      enddo
    enddo
    do jy=0,nymin1
      do ix=0,nxmin1
        kz=idx(ix,jy)
        dz1=height(iz)-etawheight(ix,jy,kz-1,n)
        dz2=etawheight(ix,jy,kz,n)-height(iz)
        dz=dz1+dz2
        ww(ix,jy,iz,n)=(wwh(ix,jy,kz-1)*pinmconv(ix,jy,kz-1)*dz2 &
             +wwh(ix,jy,kz)*pinmconv(ix,jy,kz)*dz1)/dz
      enddo
    enddo
  enddo

  ! Compute density gradients at intermediate levels
  !*************************************************

  drhodz(:,:,1,n)=(rho(:,:,2,n)-rho(:,:,1,n))/ &
       (height(2)-height(1))
  do kz=2,nz-1
    drhodz(:,:,kz,n)=(rho(:,:,kz+1,n)-rho(:,:,kz-1,n))/ &
         (height(kz+1)-height(kz-1))
  end do
  drhodz(:,:,nz,n)=drhodz(:,:,nz-1,n)


  !****************************************************************
  ! Compute slope of eta levels in windward direction and resulting
  ! vertical wind correction
  !****************************************************************

  do jy=1,ny-2
    cosf(jy)=1./cos((real(jy)*dy+ylat0)*pi180)
  enddo

  kmin=2
  idx=kmin
  do iz=2,nz-1
    do jy=1,ny-2
      do ix=1,nx-2

        inneta: do kz=idx(ix,jy),nz
          if (idx(ix,jy) .le. kz .and. (height(iz).gt.etauvheight(ix,jy,kz-1,n)).and. &
               (height(iz).le.etauvheight(ix,jy,kz,n))) then
            idx(ix,jy)=kz
            exit inneta
          endif
        enddo inneta
      enddo
    enddo

    do jy=1,ny-2
      do ix=1,nx-2
        kz=idx(ix,jy)
        dz1=height(iz)-etauvheight(ix,jy,kz-1,n)
        dz2=etauvheight(ix,jy,kz,n)-height(iz)
        dz=dz1+dz2
        ix1=ix-1
        jy1=jy-1
        ixp=ix+1
        jyp=jy+1

        dzdx1=(etauvheight(ixp,jy,kz-1,n)-etauvheight(ix1,jy,kz-1,n))/2.
        dzdx2=(etauvheight(ixp,jy,kz,n)-etauvheight(ix1,jy,kz,n))/2.
        dzdx=(dzdx1*dz2+dzdx2*dz1)/dz

        dzdy1=(etauvheight(ix,jyp,kz-1,n)-etauvheight(ix,jy1,kz-1,n))/2.
        dzdy2=(etauvheight(ix,jyp,kz,n)-etauvheight(ix,jy1,kz,n))/2.
        dzdy=(dzdy1*dz2+dzdy2*dz1)/dz

        ww(ix,jy,iz,n)=ww(ix,jy,iz,n)+(dzdx*uu(ix,jy,iz,n)*dxconst*cosf(jy)+dzdy*vv(ix,jy,iz,n)*dyconst)

      end do

    end do
  end do

  ! Keep original fields if wind_coord_type==ETA
  if (wind_coord_type.eq.'ETA') then
    uueta(:,:,:,n) = uuh(:,:,:)
    vveta(:,:,:,n) = vvh(:,:,:)
    tteta(:,:,:,n) = tth(:,:,:,n)
    qveta(:,:,:,n) = qvh(:,:,:,n)
    pveta(:,:,:,n) = pvh(:,:,:)
    rhoeta(:,:,:,n) = rhoh(:,:,:)
    prseta(:,:,:,n) = prsh(:,:,:)
    drhodzeta(:,:,1,n)=(rhoeta(:,:,2,n)-rhoeta(:,:,1,n))/ &
         (height(2)-height(1))
    do kz=2,nz-1
      drhodzeta(:,:,kz,n)=(rhoeta(:,:,kz+1,n)-rhoeta(:,:,kz-1,n))/ &
           (height(kz+1)-height(kz-1)) ! Note that this is still in SI units and not in eta
    end do
    drhodzeta(:,:,nz,n)=drhodzeta(:,:,nz-1,n)
    tvirtual(:,:,:,n)=tteta(:,:,:,n)* &
      ((qveta(:,:,:,n)+0.622)/(0.622*qveta(:,:,:,n)+0.622)) ! eq A11 from Mid-latitude atmospheric dynamics by Jonathan E. Martin
    !tvirtual(:,:,:,n)=tteta(:,:,:,n)*(1.+0.608*qveta(:,:,:,n))
    do jy=0,ny-1
      do ix=0,nx-1
        tvirtual(ix,jy,1,n)=tt2(ix,jy,1,n)*(1.+0.378*ew(td2(ix,jy,1,n),ps(ix,jy,1,n))/ps(ix,jy,1,n))
      end do 
    end do


    ! Convert w from Pa/s to eta/s, following FLEXTRA
    !************************************************
    do kz=1,nuvz-1
      if (kz.eq.1) then
        dpdeta=(akm(kz+1)-akm(kz)+(bkm(kz+1)-bkm(kz))*ps(:,:,1,n))/ &
          (wheight(kz+1)-wheight(kz))
      else if (kz.eq.nuvz-1) then
        dpdeta=(akm(kz)-akm(kz-1)+(bkm(kz)-bkm(kz-1))*ps(:,:,1,n))/ &
          (wheight(kz)-wheight(kz-1))
      else
        dpdeta=(akm(kz+1)-akm(kz-1)+(bkm(kz+1)-bkm(kz-1))*ps(:,:,1,n))/ &
          (wheight(kz+1)-wheight(kz-1))
      endif
      wweta(:,:,kz,n)=wwh(:,:,kz)/dpdeta
    end do
  endif
end subroutine verttransform_ecmwf_transform_windfields

subroutine verttransform_ecmwf_stereographic(n)
  use par_mod
  use com_mod
  use cmapf_mod, only: cc2gll
  use qvsat_mod

  implicit none

  integer, intent(in) :: n

  integer :: ix,jy,iz
  real :: xlon,ylat,xlonr
  real :: uuaux,vvaux,uupolaux,vvpolaux,ddpol,ffpol,wdummy

  if (nglobal) then
    do iz=1,nz
      do jy=int(switchnorthg)-2,nymin1
        ylat=ylat0+real(jy)*dy
        do ix=0,nxmin1
          xlon=xlon0+real(ix)*dx
          call cc2gll(northpolemap,ylat,xlon,uu(ix,jy,iz,n), &
               vv(ix,jy,iz,n),uupol(ix,jy,iz,n), &
               vvpol(ix,jy,iz,n))
          if (wind_coord_type.eq.'ETA') then
            call cc2gll(northpolemap,ylat,xlon,uueta(ix,jy,iz,n), &
                 vveta(ix,jy,iz,n),uupoleta(ix,jy,iz,n), &
                 vvpoleta(ix,jy,iz,n))
          endif
        end do
      end do
    end do


    do iz=1,nz

  ! CALCULATE FFPOL, DDPOL FOR CENTRAL GRID POINT
  !
  !   AMSnauffer Nov 18 2004 Added check for case vv=0
  !
      xlon=xlon0+real(nx/2-1)*dx
      xlonr=xlon*pi/180.
      ffpol=sqrt(uu(nx/2-1,nymin1,iz,n)**2+ &
           vv(nx/2-1,nymin1,iz,n)**2)
      if (vv(nx/2-1,nymin1,iz,n).lt.0.) then
        ddpol=atan(uu(nx/2-1,nymin1,iz,n)/ &
             vv(nx/2-1,nymin1,iz,n))-xlonr
      else if (vv(nx/2-1,nymin1,iz,n).gt.0.) then
        ddpol=pi+atan(uu(nx/2-1,nymin1,iz,n)/ &
             vv(nx/2-1,nymin1,iz,n))-xlonr
      else
        ddpol=pi/2-xlonr
      endif
      if(ddpol.lt.0.) ddpol=2.0*pi+ddpol
      if(ddpol.gt.2.0*pi) ddpol=ddpol-2.0*pi

  ! CALCULATE U,V FOR 180 DEG, TRANSFORM TO POLAR STEREOGRAPHIC GRID
      xlon=180.0
      xlonr=xlon*pi/180.
      ylat=90.0
      uuaux=-ffpol*sin(xlonr+ddpol)
      vvaux=-ffpol*cos(xlonr+ddpol)
      call cc2gll(northpolemap,ylat,xlon,uuaux,vvaux,uupolaux, &
           vvpolaux)

      jy=nymin1
      do ix=0,nxmin1
        uupol(ix,jy,iz,n)=uupolaux
        vvpol(ix,jy,iz,n)=vvpolaux
      end do
    end do

    if (wind_coord_type.eq.'ETA') then    
      do iz=1,nz

        xlon=xlon0+real(nx/2-1)*dx
        xlonr=xlon*pi/180.
        ffpol=sqrt(uueta(nx/2-1,nymin1,iz,n)**2+ &
             vveta(nx/2-1,nymin1,iz,n)**2)
        if (vveta(nx/2-1,nymin1,iz,n).lt.0.) then
          ddpol=atan(uueta(nx/2-1,nymin1,iz,n)/ &
               vveta(nx/2-1,nymin1,iz,n))-xlonr
        else if (vveta(nx/2-1,nymin1,iz,n).gt.0.) then
          ddpol=pi+atan(uueta(nx/2-1,nymin1,iz,n)/ &
               vveta(nx/2-1,nymin1,iz,n))-xlonr
        else
          ddpol=pi/2-xlonr
        endif
        if(ddpol.lt.0.) ddpol=2.0*pi+ddpol
        if(ddpol.gt.2.0*pi) ddpol=ddpol-2.0*pi

  ! CALCULATE U,V FOR 180 DEG, TRANSFORM TO POLAR STEREOGRAPHIC GRID
        xlon=180.0
        xlonr=xlon*pi/180.
        ylat=90.0
        uuaux=-ffpol*sin(xlonr+ddpol)
        vvaux=-ffpol*cos(xlonr+ddpol)
        call cc2gll(northpolemap,ylat,xlon,uuaux,vvaux,uupolaux, &
             vvpolaux)

        jy=nymin1
        do ix=0,nxmin1
          uupoleta(ix,jy,iz,n)=uupolaux
          vvpoleta(ix,jy,iz,n)=vvpolaux
        end do
      end do
    endif


  ! Fix: Set W at pole to the zonally averaged W of the next equator-
  ! ward parallel of latitude

    do iz=1,nz
      wdummy=0.
      jy=ny-2
      do ix=0,nxmin1
        wdummy=wdummy+ww(ix,jy,iz,n)
      end do
      wdummy=wdummy/real(nx)
      jy=nymin1
      do ix=0,nxmin1
        ww(ix,jy,iz,n)=wdummy
      end do
    end do

    if (wind_coord_type.eq.'ETA') then
      do iz=1,nz
        wdummy=0.
        jy=ny-2
        do ix=0,nxmin1
          wdummy=wdummy+wweta(ix,jy,iz,n)
        end do
        wdummy=wdummy/real(nx)
        jy=nymin1
        do ix=0,nxmin1
          wweta(ix,jy,iz,n)=wdummy
        end do
      end do
    endif

  endif


  ! If south pole is in the domain, calculate wind velocities in polar
  ! stereographic coordinates
  !*******************************************************************

  if (sglobal) then
    do iz=1,nz
      do jy=0,int(switchsouthg)+3
        ylat=ylat0+real(jy)*dy
        do ix=0,nxmin1
          xlon=xlon0+real(ix)*dx
          call cc2gll(southpolemap,ylat,xlon,uu(ix,jy,iz,n), &
               vv(ix,jy,iz,n),uupol(ix,jy,iz,n), &
               vvpol(ix,jy,iz,n))
          if (wind_coord_type.eq.'ETA') then
            call cc2gll(southpolemap,ylat,xlon,uueta(ix,jy,iz,n), &
                 vveta(ix,jy,iz,n),uupoleta(ix,jy,iz,n), &
                 vvpoleta(ix,jy,iz,n))
          endif
        end do
      end do
    end do

    do iz=1,nz

  ! CALCULATE FFPOL, DDPOL FOR CENTRAL GRID POINT
  !
  !   AMSnauffer Nov 18 2004 Added check for case vv=0
  !
      xlon=xlon0+real(nx/2-1)*dx
      xlonr=xlon*pi/180.
      ffpol=sqrt(uu(nx/2-1,0,iz,n)**2+ &
           vv(nx/2-1,0,iz,n)**2)
      if (vv(nx/2-1,0,iz,n).lt.0.) then
        ddpol=atan(uu(nx/2-1,0,iz,n)/ &
             vv(nx/2-1,0,iz,n))+xlonr
      else if (vv(nx/2-1,0,iz,n).gt.0.) then
        ddpol=pi+atan(uu(nx/2-1,0,iz,n)/ &
             vv(nx/2-1,0,iz,n))+xlonr
      else
        ddpol=pi/2-xlonr
      endif
      if(ddpol.lt.0.) ddpol=2.0*pi+ddpol
      if(ddpol.gt.2.0*pi) ddpol=ddpol-2.0*pi

  ! CALCULATE U,V FOR 180 DEG, TRANSFORM TO POLAR STEREOGRAPHIC GRID
      xlon=180.0
      xlonr=xlon*pi/180.
      ylat=-90.0
      uuaux=+ffpol*sin(xlonr-ddpol)
      vvaux=-ffpol*cos(xlonr-ddpol)
      call cc2gll(northpolemap,ylat,xlon,uuaux,vvaux,uupolaux, &
           vvpolaux)

      jy=0
      do ix=0,nxmin1
        uupol(ix,jy,iz,n)=uupolaux
        vvpol(ix,jy,iz,n)=vvpolaux
      end do
    end do

    if (wind_coord_type.eq.'ETA') then
      do iz=1,nz
  ! CALCULATE FFPOL, DDPOL FOR CENTRAL GRID POINT
  !
  !   AMSnauffer Nov 18 2004 Added check for case vv=0
  !
        xlon=xlon0+real(nx/2-1)*dx
        xlonr=xlon*pi/180.
        ffpol=sqrt(uueta(nx/2-1,0,iz,n)**2+ &
             vveta(nx/2-1,0,iz,n)**2)
        if (vveta(nx/2-1,0,iz,n).lt.0.) then
          ddpol=atan(uueta(nx/2-1,0,iz,n)/ &
               vveta(nx/2-1,0,iz,n))+xlonr
        else if (vveta(nx/2-1,0,iz,n).gt.0.) then
          ddpol=pi+atan(uueta(nx/2-1,0,iz,n)/ &
               vveta(nx/2-1,0,iz,n))+xlonr
        else
          ddpol=pi/2-xlonr
        endif
        if(ddpol.lt.0.) ddpol=2.0*pi+ddpol
        if(ddpol.gt.2.0*pi) ddpol=ddpol-2.0*pi

  ! CALCULATE U,V FOR 180 DEG, TRANSFORM TO POLAR STEREOGRAPHIC GRID
        xlon=180.0
        xlonr=xlon*pi/180.
        ylat=-90.0
        uuaux=+ffpol*sin(xlonr-ddpol)
        vvaux=-ffpol*cos(xlonr-ddpol)
        call cc2gll(northpolemap,ylat,xlon,uuaux,vvaux,uupolaux, &
             vvpolaux)

        jy=0
        do ix=0,nxmin1
          uupoleta(ix,jy,iz,n)=uupolaux
          vvpoleta(ix,jy,iz,n)=vvpolaux
        end do
      end do
    endif

  ! Fix: Set W at pole to the zonally averaged W of the next equator-
  ! ward parallel of latitude

    do iz=1,nz
      wdummy=0.
      jy=1
      do ix=0,nxmin1
        wdummy=wdummy+ww(ix,jy,iz,n)
      end do
      wdummy=wdummy/real(nx)
      jy=0
      do ix=0,nxmin1
        ww(ix,jy,iz,n)=wdummy
      end do
    end do

    if (wind_coord_type.eq.'ETA') then
      do iz=1,nz
        wdummy=0.
        jy=1
        do ix=0,nxmin1
          wdummy=wdummy+wweta(ix,jy,iz,n)
        end do
        wdummy=wdummy/real(nx)
        jy=0
        do ix=0,nxmin1
          wweta(ix,jy,iz,n)=wdummy
        end do
      end do
    endif
  endif
end subroutine verttransform_ecmwf_stereographic

subroutine verttransform_ecmwf_clouds(n,lreadclouds,lsumclouds,nxlim,nylim,clouds_tmp,cloudsh_tmp,&
  clw_tmp,ctwc_tmp,clwc_tmp,ciwc_tmp,lsprec_tmp,convprec_tmp,rho_tmp,tt_tmp,qv_tmp)
  use par_mod
  use com_mod
  use cmapf_mod, only: cc2gll
  use qvsat_mod

  implicit none

  logical,intent(in) :: lreadclouds,lsumclouds
  integer, intent(in) :: nxlim,nylim
  integer, intent(in) :: n
  integer(kind=1),intent(inout) :: clouds_tmp(0:nxlim,0:nylim,nzmax)
  integer,intent(inout) :: cloudsh_tmp(0:nxlim,0:nylim)
  real,intent(inout) :: clw_tmp(0:nxlim,0:nylim,nzmax)
  real,intent(inout) :: ctwc_tmp(0:nxlim,0:nylim)
  real,intent(inout) :: clwc_tmp(0:nxlim,0:nylim,nzmax)
  real,intent(in) :: ciwc_tmp(0:nxlim,0:nylim,nzmax)
  real,intent(in) :: lsprec_tmp(0:nxlim,0:nylim),convprec_tmp(0:nxlim,0:nylim)
  real,intent(in),dimension(0:nxlim,0:nylim,nzmax) :: rho_tmp,tt_tmp,qv_tmp

  integer,dimension(0:nxmax-1,0:nymax-1) :: rain_cloud_above

  integer :: ix,jy,kz,kz_inv
  real :: pressure,rh,lsp,convp,cloudh_min,prec

  !***********************************************************************************  
  if (lreadclouds) then !HG METHOD
  ! The method is loops all grids vertically and constructs the 3D matrix for clouds
  ! Cloud top and cloud bottom gid cells are assigned as well as the total column 
  ! cloud water. For precipitating grids, the type and whether it is in or below 
  ! cloud scavenging are assigned with numbers 2-5 (following the old metod).
  ! Distinction is done for lsp and convp though they are treated the same in regards
  ! to scavenging. Also clouds that are not precipitating are defined which may be 
  ! to include future cloud processing by non-precipitating-clouds. 
  !***********************************************************************************
    !write(*,*) 'Global ECMWF fields: using cloud water'
    clw_tmp(0:nxlim,0:nylim,:)=0.0
  !    icloud_stats(:,:,:,n)=0.0
    ctwc_tmp(:,:)=0.0
    clouds_tmp(0:nxlim,0:nylim,:)=0
  ! If water/ice are read separately into clwc and ciwc, store sum in clwc
    if (.not.lsumclouds) then 
      clwc_tmp(0:nxlim,0:nylim,:) = clwc_tmp(0:nxlim,0:nylim,:) + ciwc_tmp(:,:,:)
    end if
    do jy=0,nylim
      do ix=0,nxlim
        lsp=lsprec_tmp(ix,jy)
        convp=convprec_tmp(ix,jy)
        prec=lsp+convp
  !        tot_cloud_h=0
  ! Find clouds in the vertical
        do kz=1, nz-1 !go from top to bottom
          if (clwc_tmp(ix,jy,kz).gt.0) then      
  ! assuming rho is in kg/m3 and hz in m gives: kg/kg * kg/m3 *m3/kg /m = m2/m3 
            clw_tmp(ix,jy,kz)=(clwc_tmp(ix,jy,kz)*rho_tmp(ix,jy,kz))*(height(kz+1)-height(kz))
  !            tot_cloud_h=tot_cloud_h+(height(kz+1)-height(kz)) 
            
  !            icloud_stats(ix,jy,4,n)= icloud_stats(ix,jy,4,n)+clw(ix,jy,kz,n)          ! Column cloud water [m3/m3]
            ctwc_tmp(ix,jy) = ctwc_tmp(ix,jy)+clw_tmp(ix,jy,kz)
  !            icloud_stats(ix,jy,3,n)= min(height(kz+1),height(kz))                     ! Cloud BOT height stats      [m]
            cloudh_min=min(height(kz+1),height(kz))
          endif
        end do

  ! If Precipitation. Define removal type in the vertical
        if ((lsp.gt.0.01).or.(convp.gt.0.01)) then ! cloud and precipitation

          do kz=nz,2,-1 !go Bottom up!
            if (clw_tmp(ix,jy,kz).gt. 0) then ! is in cloud
              cloudsh_tmp(ix,jy)=cloudsh_tmp(ix,jy)+height(kz)-height(kz-1) 
              clouds_tmp(ix,jy,kz)=1                               ! is a cloud
              if (lsp.ge.convp) then
                clouds_tmp(ix,jy,kz)=3                            ! lsp in-cloud
              else
                clouds_tmp(ix,jy,kz)=2                             ! convp in-cloud
              endif                                              ! convective or large scale
            elseif((clw_tmp(ix,jy,kz).le.0) .and. (cloudh_min.ge.height(kz))) then ! is below cloud
              if (lsp.ge.convp) then
                clouds_tmp(ix,jy,kz)=5                             ! lsp dominated washout
              else
                clouds_tmp(ix,jy,kz)=4                             ! convp dominated washout
              endif                                              ! convective or large scale 
            endif

            if (height(kz).ge. 19000) then                        ! set a max height for removal
              clouds_tmp(ix,jy,kz)=0
            endif !clw>0
          end do !nz
        endif ! precipitation
      end do
    end do

  ! eso: copy the relevant data to clw4 to reduce amount of communicated data for MPI
  !    ctwc(:,:,n) = icloud_stats(:,:,4,n)

  !**************************************************************************
  else       ! use old definitions
  !**************************************************************************
  !   create a cloud and rainout/washout field, clouds occur where rh>80%
  !   total cloudheight is stored at level 0
    !write(*,*) 'Global fields: using cloud water from Parameterization'
    do jy=0,nylim
      do ix=0,nxlim
  ! OLD METHOD
        rain_cloud_above(ix,jy)=0
        lsp=lsprec_tmp(ix,jy)
        convp=convprec_tmp(ix,jy)
        cloudsh_tmp(ix,jy)=0
        do kz_inv=1,nz-1
          kz=nz-kz_inv+1
          pressure=rho_tmp(ix,jy,kz)*r_air*tt_tmp(ix,jy,kz)
          rh=qv_tmp(ix,jy,kz)/f_qvsat(pressure,tt_tmp(ix,jy,kz))
          clouds_tmp(ix,jy,kz)=0
          if (rh.gt.0.8) then ! in cloud
            if ((lsp.gt.0.01).or.(convp.gt.0.01)) then ! cloud and precipitation
              rain_cloud_above(ix,jy)=1
              cloudsh_tmp(ix,jy)=cloudsh_tmp(ix,jy)+ &
                   height(kz)-height(kz-1)
              if (lsp.ge.convp) then
                clouds_tmp(ix,jy,kz)=3 ! lsp dominated rainout
              else
                clouds_tmp(ix,jy,kz)=2 ! convp dominated rainout
              endif
            else ! no precipitation
              clouds_tmp(ix,jy,kz)=1 ! cloud
            endif
          else ! no cloud
            if (rain_cloud_above(ix,jy).eq.1) then ! scavenging
              if (lsp.ge.convp) then
                clouds_tmp(ix,jy,kz)=5 ! lsp dominated washout
              else
                clouds_tmp(ix,jy,kz)=4 ! convp dominated washout
              endif
            endif
          endif
        end do
  !END OLD METHOD
      end do
    end do
  endif !readclouds
end subroutine verttransform_ecmwf_clouds

subroutine verttransform_gfs(n,uuh,vvh,wwh,pvh)
  !                      i  i   i   i   i
  !*****************************************************************************
  !                                                                            *
  !     This subroutine transforms temperature, dew point temperature and      *
  !     wind components from eta to meter coordinates.                         *
  !     The vertical wind component is transformed from Pa/s to m/s using      *
  !     the conversion factor pinmconv.                                        *
  !     In addition, this routine calculates vertical density gradients        *
  !     needed for the parameterization of the turbulent velocities.           *
  !                                                                            *
  !     Author: A. Stohl, G. Wotawa                                            *
  !                                                                            *
  !     12 August 1996                                                         *
  !     Update: 16 January 1998                                                *
  !                                                                            *
  !     Major update: 17 February 1999                                         *
  !     by G. Wotawa                                                           *
  !     CHANGE 17/11/2005 Caroline Forster, NCEP GFS version                   *
  !                                                                            *
  !   - Vertical levels for u, v and w are put together                        *
  !   - Slope correction for vertical velocity: Modification of calculation    *
  !     procedure                                                              *
  !                                                                            *
  !*****************************************************************************
  !  Changes, Bernd C. Krueger, Feb. 2001:
  !   Variables tth and qvh (on eta coordinates) from common block
  !
  !   Unified ECMWF and GFS builds                                      
  !   Marian Harustak, 12.5.2017                                        
  !     - Renamed routine from verttransform to verttransform_gfs
  !
  !*****************************************************************************
  !                                                                            *
  ! Variables:                                                                 *
  ! nx,ny,nz                        field dimensions in x,y and z direction    *
  ! uu(0:nxmax,0:nymax,nzmax,2)     wind components in x-direction [m/s]       *
  ! vv(0:nxmax,0:nymax,nzmax,2)     wind components in y-direction [m/s]       *
  ! ww(0:nxmax,0:nymax,nzmax,2)     wind components in z-direction [deltaeta/s]*
  ! tt(0:nxmax,0:nymax,nzmax,2)     temperature [K]                            *
  ! pv(0:nxmax,0:nymax,nzmax,2)     potential voriticity (pvu)                 *
  ! ps(0:nxmax,0:nymax,2)           surface pressure [Pa]                      *
  ! clouds(0:nxmax,0:nymax,0:nzmax,2) cloud field for wet deposition           *
  !                                                                            *
  !*****************************************************************************

  use par_mod
  use com_mod
  use cmapf_mod
  use qvsat_mod

  implicit none

  integer :: ix,jy,kz,iz,n,kmin,kl,klp,ix1,jy1,ixp,jyp,ixm,jym
  integer :: rain_cloud_above,kz_inv
  real :: pressure
  real :: rh,lsp,cloudh_min,convp,prec
  real :: rhoh(nuvzmax),pinmconv(nzmax)
  real :: pint,tv,tvold,pold,dz1,dz2,dz,ui,vi
  real :: xlon,ylat,xlonr,dzdx,dzdy
  real :: dzdx1,dzdx2,dzdy1,dzdy2,cosf
  real :: uuaux,vvaux,uupolaux,vvpolaux,ddpol,ffpol,wdummy
  real :: uuh(0:nxmax-1,0:nymax-1,nuvzmax)
  real :: vvh(0:nxmax-1,0:nymax-1,nuvzmax)
  real :: pvh(0:nxmax-1,0:nymax-1,nuvzmax)
  real :: wwh(0:nxmax-1,0:nymax-1,nwzmax)
  real :: wzlev(nwzmax),uvwzlev(0:nxmax-1,0:nymax-1,nzmax)
  real,parameter :: const=r_air/ga

  ! NCEP version
  integer :: llev, i

  logical :: init = .true.


  !*************************************************************************
  ! If verttransform is called the first time, initialize heights of the   *
  ! z levels in meter. The heights are the heights of model levels, where  *
  ! u,v,T and qv are given.                                                *
  !*************************************************************************

  if (init) then

  ! Search for a point with high surface pressure (i.e. not above significant topography)
  ! Then, use this point to construct a reference z profile, to be used at all times
  !*****************************************************************************

    y_loop: do jy=0,nymin1
      x_loop: do ix=0,nxmin1
        if (ps(ix,jy,1,n).gt.100000.) then
          ixm=ix
          jym=jy
          exit y_loop
        endif
      end do x_loop
    end do y_loop

    tvold=tt2(ixm,jym,1,n)*(1.+0.378*ew(td2(ixm,jym,1,n),ps(ixm,jym,1,n))/ &
      ps(ixm,jym,1,n))
    pold=ps(ixm,jym,1,n)
    height(1)=0.

    do kz=2,nuvz
      pint=akz(kz)+bkz(kz)*ps(ixm,jym,1,n)
      tv=tth(ixm,jym,kz,n)*(1.+0.608*qvh(ixm,jym,kz,n))

      if (abs(tv-tvold).gt.0.2) then
        height(kz)=height(kz-1)+const*log(pold/pint)* &
        (tv-tvold)/log(tv/tvold)
      else
        height(kz)=height(kz-1)+const*log(pold/pint)*tv
      endif

      tvold=tv
      pold=pint
    end do


  ! Determine highest levels that can be within PBL
  !************************************************

    do kz=1,nz
      if (height(kz).gt.hmixmax) then
        nmixz=kz
        exit
      endif
    end do
  
  ! Do not repeat initialization of the Cartesian z grid
  !*****************************************************

    init=.false.

  endif


  ! Loop over the whole grid
  !*************************

  do jy=0,nymin1
    do ix=0,nxmin1

  ! NCEP version: find first level above ground
      llev = 0
      do i=1,nuvz
        if (ps(ix,jy,1,n).lt.akz(i)) llev=i
      end do
       llev = llev+1
       if (llev.gt.nuvz-2) llev = nuvz-2
  !     if (llev.eq.nuvz-2) write(*,*) 'verttransform
  !    +WARNING: LLEV eq NUZV-2'
  ! NCEP version


  ! compute height of pressure levels above ground
  !***********************************************

      tvold=tth(ix,jy,llev,n)*(1.+0.608*qvh(ix,jy,llev,n))
      pold=akz(llev)
      wzlev(llev)=0.
      uvwzlev(ix,jy,llev)=0.
      rhoh(llev)=pold/(r_air*tvold)

      do kz=llev+1,nuvz
        pint=akz(kz)+bkz(kz)*ps(ix,jy,1,n)
        tv=tth(ix,jy,kz,n)*(1.+0.608*qvh(ix,jy,kz,n))
        rhoh(kz)=pint/(r_air*tv)

        if (abs(tv-tvold).gt.0.2) then
          uvwzlev(ix,jy,kz)=uvwzlev(ix,jy,kz-1)+const*log(pold/pint)* &
          (tv-tvold)/log(tv/tvold)
        else
          uvwzlev(ix,jy,kz)=uvwzlev(ix,jy,kz-1)+const*log(pold/pint)*tv
        endif
        wzlev(kz)=uvwzlev(ix,jy,kz)

        tvold=tv
        pold=pint
      end do

  ! pinmconv=(h2-h1)/(p2-p1)

      pinmconv(llev)=(uvwzlev(ix,jy,llev+1)-uvwzlev(ix,jy,llev))/ &
           ((aknew(llev+1)+bknew(llev+1)*ps(ix,jy,1,n))- &
           (aknew(llev)+bknew(llev)*ps(ix,jy,1,n)))
      do kz=llev+1,nz-1
        pinmconv(kz)=(uvwzlev(ix,jy,kz+1)-uvwzlev(ix,jy,kz-1))/ &
             ((aknew(kz+1)+bknew(kz+1)*ps(ix,jy,1,n))- &
             (aknew(kz-1)+bknew(kz-1)*ps(ix,jy,1,n)))
      end do
      pinmconv(nz)=(uvwzlev(ix,jy,nz)-uvwzlev(ix,jy,nz-1))/ &
           ((aknew(nz)+bknew(nz)*ps(ix,jy,1,n))- &
           (aknew(nz-1)+bknew(nz-1)*ps(ix,jy,1,n)))


  ! Levels, where u,v,t and q are given
  !************************************

      uu(ix,jy,1,n)=uuh(ix,jy,llev)
      vv(ix,jy,1,n)=vvh(ix,jy,llev)
      tt(ix,jy,1,n)=tth(ix,jy,llev,n)
      qv(ix,jy,1,n)=qvh(ix,jy,llev,n)
  ! IP & SEC, 201812 add clouds
      if (readclouds) then
         clwc(ix,jy,1,n)=clwch(ix,jy,llev,n)
      endif 
      pv(ix,jy,1,n)=pvh(ix,jy,llev)
      rho(ix,jy,1,n)=rhoh(llev)
      pplev(ix,jy,1,n)=akz(llev)
      uu(ix,jy,nz,n)=uuh(ix,jy,nuvz)
      vv(ix,jy,nz,n)=vvh(ix,jy,nuvz)
      tt(ix,jy,nz,n)=tth(ix,jy,nuvz,n)
      qv(ix,jy,nz,n)=qvh(ix,jy,nuvz,n)
  ! IP & SEC, 201812 add clouds
      if (readclouds) then
         clwc(ix,jy,nz,n)=clwch(ix,jy,nuvz,n)
      endif
      pv(ix,jy,nz,n)=pvh(ix,jy,nuvz)
      rho(ix,jy,nz,n)=rhoh(nuvz)
      pplev(ix,jy,nz,n)=akz(nuvz)
      kmin=llev+1
      do iz=2,nz-1
        do kz=kmin,nuvz
          if(height(iz).gt.uvwzlev(ix,jy,nuvz)) then
            uu(ix,jy,iz,n)=uu(ix,jy,nz,n)
            vv(ix,jy,iz,n)=vv(ix,jy,nz,n)
            tt(ix,jy,iz,n)=tt(ix,jy,nz,n)
            qv(ix,jy,iz,n)=qv(ix,jy,nz,n)
  ! IP & SEC, 201812 add clouds
            if (readclouds) then
               clwc(ix,jy,iz,n)=clwc(ix,jy,nz,n)
            endif
            pv(ix,jy,iz,n)=pv(ix,jy,nz,n)
            rho(ix,jy,iz,n)=rho(ix,jy,nz,n)
            pplev(ix,jy,iz,n)=pplev(ix,jy,nz,n)
            exit
          endif
          if ((height(iz).gt.uvwzlev(ix,jy,kz-1)).and. &
          (height(iz).le.uvwzlev(ix,jy,kz))) then
            dz1=height(iz)-uvwzlev(ix,jy,kz-1)
            dz2=uvwzlev(ix,jy,kz)-height(iz)
            dz=dz1+dz2
            uu(ix,jy,iz,n)=(uuh(ix,jy,kz-1)*dz2+uuh(ix,jy,kz)*dz1)/dz
            vv(ix,jy,iz,n)=(vvh(ix,jy,kz-1)*dz2+vvh(ix,jy,kz)*dz1)/dz
            tt(ix,jy,iz,n)=(tth(ix,jy,kz-1,n)*dz2 &
            +tth(ix,jy,kz,n)*dz1)/dz
            qv(ix,jy,iz,n)=(qvh(ix,jy,kz-1,n)*dz2 &
            +qvh(ix,jy,kz,n)*dz1)/dz
  ! IP & SEC, 201812 add clouds
            if (readclouds) then
               clwc(ix,jy,iz,n)=(clwch(ix,jy,kz-1,n)*dz2 &
               +clwch(ix,jy,kz,n)*dz1)/dz
            endif
            pv(ix,jy,iz,n)=(pvh(ix,jy,kz-1)*dz2+pvh(ix,jy,kz)*dz1)/dz
            rho(ix,jy,iz,n)=(rhoh(kz-1)*dz2+rhoh(kz)*dz1)/dz
            pplev(ix,jy,iz,n)=(akz(kz-1)*dz2+akz(kz)*dz1)/dz
          endif
        end do
      end do


  ! Levels, where w is given
  !*************************

      ww(ix,jy,1,n)=wwh(ix,jy,llev)*pinmconv(llev)
      ww(ix,jy,nz,n)=wwh(ix,jy,nwz)*pinmconv(nz)
      kmin=llev+1
      do iz=2,nz
        do kz=kmin,nwz
          if ((height(iz).gt.wzlev(kz-1)).and. &
          (height(iz).le.wzlev(kz))) then
            dz1=height(iz)-wzlev(kz-1)
            dz2=wzlev(kz)-height(iz)
            dz=dz1+dz2
            ww(ix,jy,iz,n)=(wwh(ix,jy,kz-1)*pinmconv(kz-1)*dz2 &
            +wwh(ix,jy,kz)*pinmconv(kz)*dz1)/dz
          endif
        end do
      end do


  ! Compute density gradients at intermediate levels
  !*************************************************

      drhodz(ix,jy,1,n)=(rho(ix,jy,2,n)-rho(ix,jy,1,n))/ &
           (height(2)-height(1))
      do kz=2,nz-1
        drhodz(ix,jy,kz,n)=(rho(ix,jy,kz+1,n)-rho(ix,jy,kz-1,n))/ &
        (height(kz+1)-height(kz-1))
      end do
      drhodz(ix,jy,nz,n)=drhodz(ix,jy,nz-1,n)

    end do
  end do


  !****************************************************************
  ! Compute slope of eta levels in windward direction and resulting
  ! vertical wind correction
  !****************************************************************

  do jy=1,ny-2
    cosf=cos((real(jy)*dy+ylat0)*pi180)
    do ix=1,nx-2

  ! NCEP version: find first level above ground
      llev = 0
      do i=1,nuvz
       if (ps(ix,jy,1,n).lt.akz(i)) llev=i
      end do
       llev = llev+1
       if (llev.gt.nuvz-2) llev = nuvz-2
  !     if (llev.eq.nuvz-2) write(*,*) 'verttransform
  !    +WARNING: LLEV eq NUZV-2'
  ! NCEP version

      kmin=llev+1
      do iz=2,nz-1

        ui=uu(ix,jy,iz,n)*dxconst/cosf
        vi=vv(ix,jy,iz,n)*dyconst

        do kz=kmin,nz
          if ((height(iz).gt.uvwzlev(ix,jy,kz-1)).and. &
          (height(iz).le.uvwzlev(ix,jy,kz))) then
            dz1=height(iz)-uvwzlev(ix,jy,kz-1)
            dz2=uvwzlev(ix,jy,kz)-height(iz)
            dz=dz1+dz2
            kl=kz-1
            klp=kz
            exit
          endif
        end do

        ix1=ix-1
        jy1=jy-1
        ixp=ix+1
        jyp=jy+1

        dzdx1=(uvwzlev(ixp,jy,kl)-uvwzlev(ix1,jy,kl))/2.
        dzdx2=(uvwzlev(ixp,jy,klp)-uvwzlev(ix1,jy,klp))/2.
        dzdx=(dzdx1*dz2+dzdx2*dz1)/dz

        dzdy1=(uvwzlev(ix,jyp,kl)-uvwzlev(ix,jy1,kl))/2.
        dzdy2=(uvwzlev(ix,jyp,klp)-uvwzlev(ix,jy1,klp))/2.
        dzdy=(dzdy1*dz2+dzdy2*dz1)/dz

        ww(ix,jy,iz,n)=ww(ix,jy,iz,n)+(dzdx*ui+dzdy*vi)

      end do

    end do
  end do


  ! If north pole is in the domain, calculate wind velocities in polar
  ! stereographic coordinates
  !*******************************************************************

  if (nglobal) then
    do jy=int(switchnorthg)-2,nymin1
      ylat=ylat0+real(jy)*dy
      do ix=0,nxmin1
        xlon=xlon0+real(ix)*dx
        do iz=1,nz
          call cc2gll(northpolemap,ylat,xlon,uu(ix,jy,iz,n), &
               vv(ix,jy,iz,n),uupol(ix,jy,iz,n), &
               vvpol(ix,jy,iz,n))
        end do
      end do
    end do


    do iz=1,nz

  ! CALCULATE FFPOL, DDPOL FOR CENTRAL GRID POINT
      xlon=xlon0+real(nx/2-1)*dx
      xlonr=xlon*pi/180.
      ffpol=sqrt(uu(nx/2-1,nymin1,iz,n)**2+vv(nx/2-1,nymin1,iz,n)**2)
      if (vv(nx/2-1,nymin1,iz,n).lt.0.) then
        ddpol=atan(uu(nx/2-1,nymin1,iz,n)/vv(nx/2-1,nymin1,iz,n))-xlonr
      elseif (vv(nx/2-1,nymin1,iz,n).gt.0.) then
        ddpol=pi+atan(uu(nx/2-1,nymin1,iz,n)/ &
        vv(nx/2-1,nymin1,iz,n))-xlonr
      else
        ddpol=pi/2-xlonr
      endif
      if(ddpol.lt.0.) ddpol=2.0*pi+ddpol
      if(ddpol.gt.2.0*pi) ddpol=ddpol-2.0*pi

  ! CALCULATE U,V FOR 180 DEG, TRANSFORM TO POLAR STEREOGRAPHIC GRID
      xlon=180.0
      xlonr=xlon*pi/180.
      ylat=90.0
      uuaux=-ffpol*sin(xlonr+ddpol)
      vvaux=-ffpol*cos(xlonr+ddpol)
      call cc2gll(northpolemap,ylat,xlon,uuaux,vvaux,uupolaux,vvpolaux)
      jy=nymin1
      do ix=0,nxmin1
        uupol(ix,jy,iz,n)=uupolaux
        vvpol(ix,jy,iz,n)=vvpolaux
      end do
    end do


  ! Fix: Set W at pole to the zonally averaged W of the next equator-
  ! ward parallel of latitude

    do iz=1,nz
      wdummy=0.
      jy=ny-2
      do ix=0,nxmin1
        wdummy=wdummy+ww(ix,jy,iz,n)
      end do
      wdummy=wdummy/real(nx)
      jy=nymin1
      do ix=0,nxmin1
        ww(ix,jy,iz,n)=wdummy
      end do
    end do

  endif


  ! If south pole is in the domain, calculate wind velocities in polar
  ! stereographic coordinates
  !*******************************************************************

  if (sglobal) then
    do jy=0,int(switchsouthg)+3
      ylat=ylat0+real(jy)*dy
      do ix=0,nxmin1
        xlon=xlon0+real(ix)*dx
        do iz=1,nz
          call cc2gll(southpolemap,ylat,xlon,uu(ix,jy,iz,n), &
          vv(ix,jy,iz,n),uupol(ix,jy,iz,n),vvpol(ix,jy,iz,n))
        end do
      end do
    end do

    do iz=1,nz

  ! CALCULATE FFPOL, DDPOL FOR CENTRAL GRID POINT
      xlon=xlon0+real(nx/2-1)*dx
      xlonr=xlon*pi/180.
      ffpol=sqrt(uu(nx/2-1,0,iz,n)**2+vv(nx/2-1,0,iz,n)**2)
      if(vv(nx/2-1,0,iz,n).lt.0.) then
        ddpol=atan(uu(nx/2-1,0,iz,n)/vv(nx/2-1,0,iz,n))+xlonr
      elseif (vv(nx/2-1,0,iz,n).gt.0.) then
        ddpol=pi+atan(uu(nx/2-1,0,iz,n)/vv(nx/2-1,0,iz,n))-xlonr
      else
        ddpol=pi/2-xlonr
      endif
      if(ddpol.lt.0.) ddpol=2.0*pi+ddpol
      if(ddpol.gt.2.0*pi) ddpol=ddpol-2.0*pi

  ! CALCULATE U,V FOR 180 DEG, TRANSFORM TO POLAR STEREOGRAPHIC GRID
      xlon=180.0
      xlonr=xlon*pi/180.
      ylat=-90.0
      uuaux=+ffpol*sin(xlonr-ddpol)
      vvaux=-ffpol*cos(xlonr-ddpol)
      call cc2gll(northpolemap,ylat,xlon,uuaux,vvaux,uupolaux,vvpolaux)

      jy=0
      do ix=0,nxmin1
        uupol(ix,jy,iz,n)=uupolaux
        vvpol(ix,jy,iz,n)=vvpolaux
      end do
    end do


  ! Fix: Set W at pole to the zonally averaged W of the next equator-
  ! ward parallel of latitude

    do iz=1,nz
      wdummy=0.
      jy=1
      do ix=0,nxmin1
        wdummy=wdummy+ww(ix,jy,iz,n)
      end do
      wdummy=wdummy/real(nx)
      jy=0
      do ix=0,nxmin1
        ww(ix,jy,iz,n)=wdummy
      end do
    end do
  endif



  !***********************************************************************************
  ! IP & SEC, 201812 GFS clouds read
  if (readclouds) then
  ! The method is loops all grids vertically and constructs the 3D matrix for clouds
  ! Cloud top and cloud bottom gid cells are assigned as well as the total column
  ! cloud water. For precipitating grids, the type and whether it is in or below
  ! cloud scavenging are assigned with numbers 2-5 (following the old metod).
  ! Distinction is done for lsp and convp though they are treated the same in regards
  ! to scavenging. Also clouds that are not precipitating are defined which may be
  ! to include future cloud processing by non-precipitating-clouds.
  !***********************************************************************************
    write(*,*) 'Global NCEP fields: using cloud water'
    clw(:,:,:,n)=0.0
    ctwc(:,:,n)=0.0
    clouds(:,:,:,n)=0
  ! If water/ice are read separately into clwc and ciwc, store sum in clwc
    do jy=0,nymin1
      do ix=0,nxmin1
        lsp=lsprec(ix,jy,1,n)
        convp=convprec(ix,jy,1,n)
        prec=lsp+convp
  ! Find clouds in the vertical
        do kz=1, nz-1 !go from top to bottom
          if (clwc(ix,jy,kz,n).gt.0) then
  ! assuming rho is in kg/m3 and hz in m gives: kg/kg * kg/m3 *m3/kg /m = m2/m3
            clw(ix,jy,kz,n)=(clwc(ix,jy,kz,n)*rho(ix,jy,kz,n))*(height(kz+1)-height(kz))
            ctwc(ix,jy,n) = ctwc(ix,jy,n)+clw(ix,jy,kz,n)
            cloudh_min=min(height(kz+1),height(kz))
          endif
        end do

  ! If Precipitation. Define removal type in the vertical
        if ((lsp.gt.0.01).or.(convp.gt.0.01)) then ! cloud and precipitation

          do kz=nz,2,-1 !go Bottom up!
            if (clw(ix,jy,kz,n).gt. 0) then ! is in cloud
              cloudsh(ix,jy,n)=cloudsh(ix,jy,n)+height(kz)-height(kz-1)
              clouds(ix,jy,kz,n)=1                               ! is a cloud
              if (lsp.ge.convp) then
                clouds(ix,jy,kz,n)=3                             ! lsp in-cloud
              else
                clouds(ix,jy,kz,n)=2                             ! convp in-cloud
              endif                                              ! convective or large scale
            elseif((clw(ix,jy,kz,n).le.0) .and. (cloudh_min.ge.height(kz))) then ! is below cloud
              if (lsp.ge.convp) then
                clouds(ix,jy,kz,n)=5                             ! lsp dominated washout
              else
                clouds(ix,jy,kz,n)=4                             ! convp dominated washout
              endif                                              ! convective or large scale
            endif

            if (height(kz).ge. 19000) then                        ! set a max height for removal
              clouds(ix,jy,kz,n)=0
            endif !clw>0
          end do !nz
        endif ! precipitation
      end do
    end do
  else
  write(*,*) 'Global NCEP fields: using cloud water from Parameterization'
  !   write (*,*) 'initializing clouds, n:',n,nymin1,nxmin1,nz
  !   create a cloud and rainout/washout field, clouds occur where rh>80%
  !   total cloudheight is stored at level 0
  do jy=0,nymin1
    do ix=0,nxmin1
      rain_cloud_above=0
      lsp=lsprec(ix,jy,1,n)
      convp=convprec(ix,jy,1,n)
      cloudsh(ix,jy,n)=0
      do kz_inv=1,nz-1
         kz=nz-kz_inv+1
         pressure=rho(ix,jy,kz,n)*r_air*tt(ix,jy,kz,n)
         rh=qv(ix,jy,kz,n)/f_qvsat(pressure,tt(ix,jy,kz,n))
         clouds(ix,jy,kz,n)=0
         if (rh.gt.0.8) then ! in cloud
           if ((lsp.gt.0.01).or.(convp.gt.0.01)) then ! cloud and precipitation
              rain_cloud_above=1
              cloudsh(ix,jy,n)=cloudsh(ix,jy,n)+height(kz)-height(kz-1)
              if (lsp.ge.convp) then
                 clouds(ix,jy,kz,n)=3 ! lsp dominated rainout
              else
                 clouds(ix,jy,kz,n)=2 ! convp dominated rainout
              endif
           else ! no precipitation
             clouds(ix,jy,kz,n)=1 ! cloud
           endif
         else ! no cloud
           if (rain_cloud_above.eq.1) then ! scavenging
             if (lsp.ge.convp) then
               clouds(ix,jy,kz,n)=5 ! lsp dominated washout
             else
               clouds(ix,jy,kz,n)=4 ! convp dominated washout
             endif
           endif
         endif
      end do
    end do
  end do
  endif  ! IP & SEC 201812, GFS clouds read
end subroutine verttransform_gfs

subroutine verttransform_nests(n,uuhn,vvhn,wwhn,pvhn)
  !                            i   i    i    i   i
  !*****************************************************************************
  !                                                                            *
  !     This subroutine transforms temperature, dew point temperature and      *
  !     wind components from eta to meter coordinates.                         *
  !     The vertical wind component is transformed from Pa/s to m/s using      *
  !     the conversion factor pinmconv.                                        *
  !     In addition, this routine calculates vertical density gradients        *
  !     needed for the parameterization of the turbulent velocities.           *
  !     It is similar to verttransform, but makes the transformations for      *
  !     the nested grids.                                                      *
  !                                                                            *
  !     Author: A. Stohl, G. Wotawa                                            *
  !                                                                            *
  !     12 August 1996                                                         *
  !     Update: 16 January 1998                                                *
  !                                                                            *
  !     Major update: 17 February 1999                                         *
  !     by G. Wotawa                                                           *
  !                                                                            *
  !     - Vertical levels for u, v and w are put together                      *
  !     - Slope correction for vertical velocity: Modification of calculation  *
  !       procedure                                                            *
  !                                                                            *
  !*****************************************************************************
  !  Changes, Bernd C. Krueger, Feb. 2001:       (marked "C-cv")
  !   Variables tthn and qvhn (on eta coordinates) from common block
  !*****************************************************************************
  ! Sabine Eckhardt, March 2007
  ! add the variable cloud for use with scavenging - descr. in com_mod
  !*****************************************************************************
  ! ESO, 2016
  ! -note that divide-by-zero occurs when nxmaxn,nymaxn etc. are larger than 
  !  the actual field dimensions
  !*****************************************************************************
  ! Date: 2017-05-30 modification of a bug in ew. Don Morton (CTBTO project)   *
  !*****************************************************************************
  !                                                                            *
  ! Variables:                                                                 *
  ! nxn,nyn,nuvz,nwz                field dimensions in x,y and z direction    *
  ! uun                             wind components in x-direction [m/s]       *
  ! vvn                             wind components in y-direction [m/s]       *
  ! wwn                             wind components in z-direction [deltaeta/s]*
  ! ttn                             temperature [K]                            *
  ! pvn                             potential vorticity (pvu)                  *
  ! psn                             surface pressure [Pa]                      *
  !                                                                            *
  !*****************************************************************************

  use par_mod
  use com_mod
  use qvsat_mod

  implicit none

  real,intent(in),dimension(0:nxmaxn-1,0:nymaxn-1,nuvzmax,maxnests) :: uuhn,vvhn,pvhn
  real,intent(in),dimension(0:nxmaxn-1,0:nymaxn-1,nwzmax,maxnests) :: wwhn

  real,dimension(0:nxmaxn-1,0:nymaxn-1,nuvzmax) :: rhohn,uvzlev,wzlev
  real,dimension(0:nxmaxn-1,0:nymaxn-1,nzmax) :: pinmconv
  real,dimension(0:nxmaxn-1,0:nymaxn-1) :: tvold,pold,pint,tv,dpdeta
  real,dimension(0:nymaxn-1) :: cosf

  integer,dimension(0:nxmaxn-1,0:nymaxn-1) :: rain_cloud_above, idx

  integer :: ix,jy,kz,iz,n,l,kmin,kl,klp,ix1,jy1,ixp,jyp,kz_inv
  real :: pressure,rh,lsp,convp,cloudh_min,prec

  real :: dz1,dz2,dz
  real :: dzdx,dzdy
  real :: dzdx1,dzdx2,dzdy1,dzdy2
  real,parameter :: const=r_air/ga
  real :: tot_cloud_h
  integer :: nxm1, nym1

  !  real,parameter :: precmin = 0.002 ! minimum prec in mm/h for cloud diagnostics

  ! Loop over all nests
  !********************

  do l=1,numbnests
    nxm1=nxn(l)-1 
    nym1=nyn(l)-1 

  ! Loop over the whole grid
  !*************************

    do jy=0,nyn(l)-1
      do ix=0,nxn(l)-1
        tvold(ix,jy)=tt2n(ix,jy,1,n,l)*(1.+0.378*ew(td2n(ix,jy,1,n,l),psn(ix,jy,1,n,l))/ &
             psn(ix,jy,1,n,l))
      end do
    end do

    pold(0:nxm1,0:nym1)=psn(0:nxm1,0:nym1,1,n,l)
    uvzlev(0:nxm1,0:nym1,1)=0.
    wzlev(0:nxm1,0:nym1,1)=0.
    rhohn(0:nxm1,0:nym1,1)=pold(0:nxm1,0:nym1)/(r_air*tvold(0:nxm1,0:nym1))

  ! Compute heights of eta levels
  !******************************

    do kz=2,nuvz
      pint(0:nxm1,0:nym1)=akz(kz)+bkz(kz)*psn(0:nxm1,0:nym1,1,n,l)
      tv(0:nxm1,0:nym1)=tthn(0:nxm1,0:nym1,kz,n,l)*(1.+0.608*qvhn(0:nxm1,0:nym1,kz,n,l))
      rhohn(0:nxm1,0:nym1,kz)=pint(0:nxm1,0:nym1)/(r_air*tv(0:nxm1,0:nym1))

      where (abs(tv(0:nxm1,0:nym1)-tvold(0:nxm1,0:nym1)).gt.0.2) 
        uvzlev(0:nxm1,0:nym1,kz)=uvzlev(0:nxm1,0:nym1,kz-1)+const*&
             &log(pold(0:nxm1,0:nym1)/pint(0:nxm1,0:nym1))* &
             (tv(0:nxm1,0:nym1)-tvold(0:nxm1,0:nym1))/&
             &log(tv(0:nxm1,0:nym1)/tvold(0:nxm1,0:nym1))
      elsewhere
        uvzlev(0:nxm1,0:nym1,kz)=uvzlev(0:nxm1,0:nym1,kz-1)+const*&
             &log(pold(0:nxm1,0:nym1)/pint(0:nxm1,0:nym1))*tv(0:nxm1,0:nym1)
      endwhere

      tvold(0:nxm1,0:nym1)=tv(0:nxm1,0:nym1)
      pold(0:nxm1,0:nym1)=pint(0:nxm1,0:nym1)

    end do

    do kz=2,nwz-1
      wzlev(0:nxm1,0:nym1,kz)=(uvzlev(0:nxm1,0:nym1,kz+1)+uvzlev(0:nxm1,0:nym1,kz))/2.
    end do
    wzlev(0:nxm1,0:nym1,nwz)=wzlev(0:nxm1,0:nym1,nwz-1)+ &
         uvzlev(0:nxm1,0:nym1,nuvz)-uvzlev(0:nxm1,0:nym1,nuvz-1)


    pinmconv(0:nxm1,0:nym1,1)=(uvzlev(0:nxm1,0:nym1,2))/ &
         ((aknew(2)+bknew(2)*psn(0:nxm1,0:nym1,1,n,l))- &
         (aknew(1)+bknew(1)*psn(0:nxm1,0:nym1,1,n,l)))
    do kz=2,nz-1
      pinmconv(0:nxm1,0:nym1,kz)=(uvzlev(0:nxm1,0:nym1,kz+1)-uvzlev(0:nxm1,0:nym1,kz-1))/ &
           ((aknew(kz+1)+bknew(kz+1)*psn(0:nxm1,0:nym1,1,n,l))- &
           (aknew(kz-1)+bknew(kz-1)*psn(0:nxm1,0:nym1,1,n,l)))
    end do
    pinmconv(0:nxm1,0:nym1,nz)=(uvzlev(0:nxm1,0:nym1,nz)-uvzlev(0:nxm1,0:nym1,nz-1))/ &
         ((aknew(nz)+bknew(nz)*psn(0:nxm1,0:nym1,1,n,l))- &
         (aknew(nz-1)+bknew(nz-1)*psn(0:nxm1,0:nym1,1,n,l)))

  ! Levels, where u,v,t and q are given
  !************************************

    uun(0:nxm1,0:nym1,1,n,l)=uuhn(0:nxm1,0:nym1,1,l)
    vvn(0:nxm1,0:nym1,1,n,l)=vvhn(0:nxm1,0:nym1,1,l)
    ttn(0:nxm1,0:nym1,1,n,l)=tthn(0:nxm1,0:nym1,1,n,l)
    qvn(0:nxm1,0:nym1,1,n,l)=qvhn(0:nxm1,0:nym1,1,n,l)
    if (readclouds_nest(l)) then
      clwcn(0:nxm1,0:nym1,1,n,l)=clwchn(0:nxm1,0:nym1,1,n,l)
      if (.not.sumclouds_nest(l)) ciwcn(0:nxm1,0:nym1,1,n,l)=ciwchn(0:nxm1,0:nym1,1,n,l)
    end if
    pvn(0:nxm1,0:nym1,1,n,l)=pvhn(0:nxm1,0:nym1,1,l)
    rhon(0:nxm1,0:nym1,1,n,l)=rhohn(0:nxm1,0:nym1,1)

    uun(0:nxm1,0:nym1,nz,n,l)=uuhn(0:nxm1,0:nym1,nuvz,l)
    vvn(0:nxm1,0:nym1,nz,n,l)=vvhn(0:nxm1,0:nym1,nuvz,l)
    ttn(0:nxm1,0:nym1,nz,n,l)=tthn(0:nxm1,0:nym1,nuvz,n,l)
    qvn(0:nxm1,0:nym1,nz,n,l)=qvhn(0:nxm1,0:nym1,nuvz,n,l)
    if (readclouds_nest(l)) then
      clwcn(0:nxm1,0:nym1,nz,n,l)=clwchn(0:nxm1,0:nym1,nuvz,n,l)
      if (.not.sumclouds_nest(l)) ciwcn(0:nxm1,0:nym1,nz,n,l)=ciwchn(0:nxm1,0:nym1,nuvz,n,l)
    end if
    pvn(0:nxm1,0:nym1,nz,n,l)=pvhn(0:nxm1,0:nym1,nuvz,l)
    rhon(0:nxm1,0:nym1,nz,n,l)=rhohn(0:nxm1,0:nym1,nuvz)


    kmin=2
    idx=kmin
    do iz=2,nz-1
      do jy=0,nyn(l)-1
        do ix=0,nxn(l)-1
          if(height(iz).gt.uvzlev(ix,jy,nuvz)) then
            uun(ix,jy,iz,n,l)=uun(ix,jy,nz,n,l)
            vvn(ix,jy,iz,n,l)=vvn(ix,jy,nz,n,l)
            ttn(ix,jy,iz,n,l)=ttn(ix,jy,nz,n,l)
            qvn(ix,jy,iz,n,l)=qvn(ix,jy,nz,n,l)
  !hg adding the cloud water
            if (readclouds_nest(l)) then
              clwcn(ix,jy,iz,n,l)=clwcn(ix,jy,nz,n,l)
              if (.not.sumclouds_nest(l)) ciwcn(ix,jy,iz,n,l)=ciwcn(ix,jy,nz,n,l)
            end if
  !hg
            pvn(ix,jy,iz,n,l)=pvn(ix,jy,nz,n,l)
            rhon(ix,jy,iz,n,l)=rhon(ix,jy,nz,n,l)
          else
            innuvz: do kz=idx(ix,jy),nuvz
              if (idx(ix,jy) .le. kz .and. (height(iz).gt.uvzlev(ix,jy,kz-1)).and. &
                   (height(iz).le.uvzlev(ix,jy,kz))) then
                idx(ix,jy)=kz
                exit innuvz
              endif
            enddo innuvz
          endif
        enddo
      enddo
      do jy=0,nyn(l)-1
        do ix=0,nxn(l)-1
          if(height(iz).le.uvzlev(ix,jy,nuvz)) then
            kz=idx(ix,jy)
            dz1=height(iz)-uvzlev(ix,jy,kz-1)
            dz2=uvzlev(ix,jy,kz)-height(iz)
            dz=dz1+dz2
            uun(ix,jy,iz,n,l)=(uuhn(ix,jy,kz-1,l)*dz2+uuhn(ix,jy,kz,l)*dz1)/dz
            vvn(ix,jy,iz,n,l)=(vvhn(ix,jy,kz-1,l)*dz2+vvhn(ix,jy,kz,l)*dz1)/dz
            ttn(ix,jy,iz,n,l)=(tthn(ix,jy,kz-1,n,l)*dz2 &
                 +tthn(ix,jy,kz,n,l)*dz1)/dz
            qvn(ix,jy,iz,n,l)=(qvhn(ix,jy,kz-1,n,l)*dz2 &
                 +qvhn(ix,jy,kz,n,l)*dz1)/dz
  !hg adding the cloud water
            if (readclouds_nest(l)) then
              clwcn(ix,jy,iz,n,l)=(clwchn(ix,jy,kz-1,n,l)*dz2+clwchn(ix,jy,kz,n,l)*dz1)/dz
              if (.not.sumclouds_nest(l)) &
                   &ciwcn(ix,jy,iz,n,l)=(ciwchn(ix,jy,kz-1,n,l)*dz2+ciwchn(ix,jy,kz,n,l)*dz1)/dz
            end if
  !hg
            pvn(ix,jy,iz,n,l)=(pvhn(ix,jy,kz-1,l)*dz2+pvhn(ix,jy,kz,l)*dz1)/dz
            rhon(ix,jy,iz,n,l)=(rhohn(ix,jy,kz-1)*dz2+rhohn(ix,jy,kz)*dz1)/dz
          endif
        enddo
      enddo
    enddo

  ! Levels, where w is given
  !*************************

    wwn(0:nxm1,0:nym1,1,n,l)=wwhn(0:nxm1,0:nym1,1,l)*pinmconv(0:nxm1,0:nym1,1)
    wwn(0:nxm1,0:nym1,nz,n,l)=wwhn(0:nxm1,0:nym1,nwz,l)*pinmconv(0:nxm1,0:nym1,nz)
    kmin=2
    idx=kmin
    do iz=2,nz
      do jy=0,nym1
        do ix=0,nxm1
          inn:         do kz=idx(ix,jy),nwz
            if(idx(ix,jy) .le. kz .and. height(iz).gt.wzlev(ix,jy,kz-1).and. &
                 height(iz).le.wzlev(ix,jy,kz)) then
              idx(ix,jy)=kz
              exit inn
            endif
          enddo inn
        enddo
      enddo
      do jy=0,nym1
        do ix=0,nxm1
          kz=idx(ix,jy)
          dz1=height(iz)-wzlev(ix,jy,kz-1)
          dz2=wzlev(ix,jy,kz)-height(iz)
          dz=dz1+dz2
          wwn(ix,jy,iz,n,l)=(wwhn(ix,jy,kz-1,l)*pinmconv(ix,jy,kz-1)*dz2 &
               +wwhn(ix,jy,kz,l)*pinmconv(ix,jy,kz)*dz1)/dz
        enddo
      enddo
    enddo

  ! Compute density gradients at intermediate levels
  !*************************************************

    drhodzn(0:nxm1,0:nym1,1,n,l)=(rhon(0:nxm1,0:nym1,2,n,l)-rhon(0:nxm1,0:nym1,1,n,l))/ &
         (height(2)-height(1))
    do kz=2,nz-1
      drhodzn(0:nxm1,0:nym1,kz,n,l)=(rhon(0:nxm1,0:nym1,kz+1,n,l)-rhon(0:nxm1,0:nym1,kz-1,n,l))/ &
           (height(kz+1)-height(kz-1))
    end do
    drhodzn(0:nxm1,0:nym1,nz,n,l)=drhodzn(0:nxm1,0:nym1,nz-1,n,l)


  !****************************************************************
  ! Compute slope of eta levels in windward direction and resulting
  ! vertical wind correction
  !****************************************************************

    do jy=1,nyn(l)-2
      cosf(jy)=1./cos((real(jy)*dyn(l)+ylat0n(l))*pi180)
    end do

    kmin=2
    idx=kmin
    do iz=2,nz-1
      do jy=1,nyn(l)-2
        do ix=1,nxn(l)-2

          inneta: do kz=idx(ix,jy),nz
            if (idx(ix,jy) .le. kz .and. (height(iz).gt.uvzlev(ix,jy,kz-1)).and. &
                 (height(iz).le.uvzlev(ix,jy,kz))) then
              idx(ix,jy)=kz
              exit inneta
            endif
          enddo inneta
        enddo
      enddo

      do jy=1,nyn(l)-2
        do ix=1,nxn(l)-2
          kz=idx(ix,jy)
          dz1=height(iz)-uvzlev(ix,jy,kz-1)
          dz2=uvzlev(ix,jy,kz)-height(iz)
          dz=dz1+dz2
          ix1=ix-1
          jy1=jy-1
          ixp=ix+1
          jyp=jy+1

          dzdx1=(uvzlev(ixp,jy,kz-1)-uvzlev(ix1,jy,kz-1))/2.
          dzdx2=(uvzlev(ixp,jy,kz)-uvzlev(ix1,jy,kz))/2.
          dzdx=(dzdx1*dz2+dzdx2*dz1)/dz

          dzdy1=(uvzlev(ix,jyp,kz-1)-uvzlev(ix,jy1,kz-1))/2.
          dzdy2=(uvzlev(ix,jyp,kz)-uvzlev(ix,jy1,kz))/2.
          dzdy=(dzdy1*dz2+dzdy2*dz1)/dz

          wwn(ix,jy,iz,n,l)=wwn(ix,jy,iz,n,l)+(dzdx*uun(ix,jy,iz,n,l)*dxconst*xresoln(l)*cosf(jy)+&
               &dzdy*vvn(ix,jy,iz,n,l)*dyconst*yresoln(l))

        end do

      end do
    end do

  ! Keep original fields if wind_coord_type==ETA
    if (wind_coord_type.eq.'ETA') then
      uuetan(0:nxm1,0:nym1,1:nz-1,n,l) = uuhn(0:nxm1,0:nym1,1:nz-1,l)
      vvetan(0:nxm1,0:nym1,1:nz-1,n,l) = vvhn(0:nxm1,0:nym1,1:nz-1,l)
      ttetan(0:nxm1,0:nym1,1:nz-1,n,l) = tthn(0:nxm1,0:nym1,1:nz-1,n,l)
      qvetan(0:nxm1,0:nym1,1:nz-1,n,l) = qvhn(0:nxm1,0:nym1,1:nz-1,n,l)
      pvetan(0:nxm1,0:nym1,1:nz-1,n,l) = pvhn(0:nxm1,0:nym1,1:nz-1,l)
      rhoetan(0:nxm1,0:nym1,1:nz-1,n,l) = rhohn(0:nxm1,0:nym1,1:nz-1)

      uuetan(0:nxm1,0:nym1,nz,n,l)=uuhn(0:nxm1,0:nym1,nuvz,l)
      vvetan(0:nxm1,0:nym1,nz,n,l)=vvhn(0:nxm1,0:nym1,nuvz,l)
      ttetan(0:nxm1,0:nym1,nz,n,l)=tthn(0:nxm1,0:nym1,nuvz,n,l)
      qvetan(0:nxm1,0:nym1,nz,n,l)=qvhn(0:nxm1,0:nym1,nuvz,n,l)
      pvetan(0:nxm1,0:nym1,nz,n,l)=pvhn(0:nxm1,0:nym1,nuvz,l)
      rhoetan(0:nxm1,0:nym1,nz,n,l)=rhohn(0:nxm1,0:nym1,nuvz)
      

      drhodzetan(0:nxm1,0:nym1,1,n,l)=(rhoetan(0:nxm1,0:nym1,2,n,l)-rhoetan(0:nxm1,0:nym1,1,n,l))/ &
           (height(2)-height(1))
      do kz=2,nz-1
        drhodzetan(0:nxm1,0:nym1,kz,n,l)=(rhoetan(0:nxm1,0:nym1,kz+1,n,l)-rhoetan(0:nxm1,0:nym1,kz-1,n,l))/ &
             (height(kz+1)-height(kz-1)) ! Note that this is still in SI units and not in eta
      end do
      drhodzetan(0:nxm1,0:nym1,nz,n,l)=drhodzetan(0:nxm1,0:nym1,nz-1,n,l)

      ! Convert w from Pa/s to eta/s, following FLEXTRA
      !************************************************
      do kz=1,nuvz-1
        if (kz.eq.1) then
          dpdeta=(akm(kz+1)-akm(kz)+(bkm(kz+1)-bkm(kz))*ps(:,:,1,n))/ &
            (wheight(kz+1)-wheight(kz))
        else if (kz.eq.nuvz-1) then
          dpdeta=(akm(kz)-akm(kz-1)+(bkm(kz)-bkm(kz-1))*ps(:,:,1,n))/ &
            (wheight(kz)-wheight(kz-1))
        else
          dpdeta=(akm(kz+1)-akm(kz-1)+(bkm(kz+1)-bkm(kz-1))*ps(:,:,1,n))/ &
            (wheight(kz+1)-wheight(kz-1))
        endif
        wwetan(:,:,kz,n,l)=wwhn(:,:,kz,l)/dpdeta
      end do      
    endif
 

    ! Create cloud fields
    !*********************
    call verttransform_ecmwf_clouds(n,readclouds_nest(l),sumclouds_nest(l),nxm1,nym1,        &
      cloudsn(0:nxm1,0:nym1,:,n,l),cloudshn(0:nxm1,0:nym1,n,l),clwn(0:nxm1,0:nym1,:,n,l),  &
      ctwcn(0:nxm1,0:nym1,n,l),clwcn(0:nxm1,0:nym1,:,n,l),ciwcn(0:nxm1,0:nym1,:,n,l),        &
      lsprecn(0:nxm1,0:nym1,1,n,l),convprecn(0:nxm1,0:nym1,1,n,l),rhon(0:nxm1,0:nym1,:,n,l), &
      ttn(0:nxm1,0:nym1,:,n,l),qvn(0:nxm1,0:nym1,:,n,l))
      
  end do ! end loop over nests
end subroutine verttransform_nests

end module verttransform_mod