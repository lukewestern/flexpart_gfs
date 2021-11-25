! SPDX-FileCopyrightText: FLEXPART 1998-2019, see flexpart_license.txt
! SPDX-License-Identifier: GPL-3.0-or-later

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

  implicit none

  integer :: metdata_format
  integer :: n,ix,jy,i,kz,lz,kzmin,llev,loop_start
  real :: ttlev(nuvzmax),qvlev(nuvzmax),obukhov,scalev,ol,hmixplus
  real :: ulev(nuvzmax),vlev(nuvzmax),ew,rh,vd(maxspec),subsceff,ylat
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

  implicit none

  integer :: metdata_format
  integer :: n,ix,jy,i,l,kz,lz,kzmin
  real :: ttlev(nuvzmax),qvlev(nuvzmax),obukhov,scalev,ol,hmixplus,dummyakzllev
  real :: ulev(nuvzmax),vlev(nuvzmax),ew,rh,vd(maxspec),subsceff,ylat
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
        rh=ew(td2n(ix,jy,1,n,l))/ew(tt2n(ix,jy,1,n,l))

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

      tvold=tt2n(ix,jy,1,n,l)*(1.+0.378*ew(td2n(ix,jy,1,n,l))/ &
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

  implicit none

  integer :: metdata_format
  real :: akm(nwzmax),bkm(nwzmax)
  real :: ps,tsurf,tdsurf,tlev,ustar,hf,e,ew,tv,rhoa,plev
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

  implicit none

  integer :: metdata_format
  integer :: i,k,nuvz,iter,llev,loop_start
  real :: tv,tvold,zref,z,zold,pint,pold,theta,thetaref,ri
  real :: akz(nuvz),bkz(nuvz),ulev(nuvz),vlev(nuvz),hf,wst,tt2,td2,ew
  real :: psurf,ust,ttlev(nuvz),qvlev(nuvz),h,excess
  real :: thetaold,zl,ul,vl,thetal,ril,hmixplus,wspeed,bvfsq,bvf
  real :: f_qvsat,rh,rhold,rhl,theta1,theta2,zl1,zl2,thetam
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

  implicit none

  real :: ps,t,td,e,ew,tv,rhoa,stress

  e=ew(td,ps)                       ! vapor pressure
  tv=t*(1.+0.378*e/ps)           ! virtual temperature
  rhoa=ps/(r_air*tv)              ! air density
  scalev=sqrt(abs(stress)/rhoa)

end function scalev