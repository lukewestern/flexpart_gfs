! SPDX-FileCopyrightText: FLEXPART 1998-2019, see flexpart_license.txt
! SPDX-License-Identifier: GPL-3.0-or-later

subroutine advance(itime,nrelpoint,ldt,up,vp,wp, &
       usigold,vsigold,wsigold,nstop,xt,yt,zt,zteta,prob,icbt,pp)
  !                     i    i  i/oi/oi/o
  !  i/o     i/o     i/o     o  i/oi/oi/o i/o  i/o
  !*****************************************************************************
  !                                                                            *
  !  Calculation of turbulent particle trajectories utilizing a                *
  !  zero-acceleration scheme, which is corrected by a numerically more        *
  !  accurate Petterssen scheme whenever possible.                             *
  !                                                                            *
  !  Particle positions are read in, incremented, and returned to the calling  *
  !  program.                                                                  *
  !                                                                            *
  !  In different regions of the atmosphere (PBL vs. free troposphere),        *
  !  different parameters are needed for advection, parameterizing turbulent   *
  !  velocities, etc. For efficiency, different interpolation routines have    *
  !  been written for these different cases, with the disadvantage that there  *
  !  exist several routines doing almost the same. They all share the          *
  !  included file 'interpol_mod'. The following                               *
  !  interpolation routines are used:                                          *
  !                                                                            *
  !  interpol_all(_nests)     interpolates everything (called inside the PBL)  *
  !  interpol_misslev(_nests) if a particle moves vertically in the PBL,       *
  !                           additional parameters are interpolated if it     *
  !                           crosses a model level                            *
  !  interpol_wind(_nests)    interpolates the wind and determines the         *
  !                           standard deviation of the wind (called outside   *
  !                           PBL) also interpolates potential vorticity       *
  !  interpol_wind_short(_nests) only interpolates the wind (needed for the    *
  !                           Petterssen scheme)                               *
  !  interpol_vdep(_nests)    interpolates deposition velocities               *
  !                                                                            *
  !                                                                            *
  !     Author: A. Stohl                                                       *
  !                                                                            *
  !     16 December 1997                                                       *
  !                                                                            *
  !  Changes:                                                                  *
  !                                                                            *
  !  8 April 2000: Deep convection parameterization                            *
  !                                                                            *
  !  May 2002: Petterssen scheme introduced                                    *
  !                                                                            *
  !*****************************************************************************
  !                                                                            *
  ! Variables:                                                                 *
  ! icbt               1 if particle not transferred to forbidden state,       *
  !                    else -1                                                 *
  ! dawsave            accumulated displacement in along-wind direction        *
  ! dcwsave            accumulated displacement in cross-wind direction        *
  ! dxsave             accumulated displacement in longitude                   *
  ! dysave             accumulated displacement in latitude                    *
  ! h [m]              Mixing height                                           *
  ! lwindinterv [s]    time interval between two wind fields                   *
  ! itime [s]          time at which this subroutine is entered                *
  ! itimec [s]         actual time, which is incremented in this subroutine    *
  ! href [m]           height for which dry deposition velocity is calculated  *
  ! ladvance [s]       Total integration time period                           *
  ! ldirect            1 forward, -1 backward                                  *
  ! ldt [s]            Time step for the next integration                      *
  ! lsynctime [s]      Synchronisation interval of FLEXPART                    *
  ! ngrid              index which grid is to be used                          *
  ! nrand              index for a variable to be picked from rannumb          *
  ! nstop              if > 1 particle has left domain and must be stopped     *
  ! prob               probability of absorption due to dry deposition         *
  ! rannumb(maxrand)   normally distributed random variables                   *
  ! rhoa               air density                                             *
  ! rhograd            vertical gradient of the air density                    *
  ! up,vp,wp           random velocities due to turbulence (along wind, cross  *
  !                    wind, vertical wind                                     *
  ! usig,vsig,wsig     mesoscale wind fluctuations                             *
  ! usigold,vsigold,wsigold  like usig, etc., but for the last time step       *
  ! vdepo              Deposition velocities for all species                   *
  ! xt,yt,zt           Particle position                                       *
  !                                                                            *
  !*****************************************************************************

  use point_mod
  use par_mod
  use com_mod
  use interpol_mod
  use hanna_mod
  use cmapf_mod
  use random_mod, only: ran3

  ! openmp change
  use omp_lib, only: OMP_GET_THREAD_NUM
  ! openmp change end

  implicit none
  real, parameter ::              &
    eps=nxmax/3.e5,               &
    eps2=1.e-9,                   &
    eps3=tiny(1.0),               &
    eps_eta=1.e-4
  integer, intent(in) ::          &
    itime,                        & ! time index
    nrelpoint,                    & ! particle index
    pp                              ! temporary, will be removed
  integer, intent(inout) ::       &
    nstop,                        & ! flag to stop particle if it leaves the domain
    ldt                             ! next timestep
  integer(kind=2), intent(inout) :: &
    icbt                            ! flag for forbidden state particle
  integer ::                      &
    itimec,                       &
    i,j,k,                        & ! loop variables
    nrand,                        & ! random number used for turbulence
    loop,                         & ! loop variable for time in the PBL
    memindnext,                   & ! seems useless
    mind,                         & ! windfield index
    ngr,                          & ! temporary new grid index of moved particle
    nix,njy,                      & ! nexted grid indices
    ks,nsp,                       & ! loop variables for vertical levels
    flagrein,                     & ! flag used in CBL scheme
    thread,                       & ! number of openmp threads (probably can be removed)
    idummy = -7                     ! used in random number routines
  real, intent(inout) ::          &
    zt,                           & ! z particle position in meters, want to keep this local to advance in future
    zteta,                        & ! z particle position in eta coordinates
    up,vp,wp,                     & ! random turbulence velocities
    usigold,vsigold,wsigold         ! old mesoscale wind fluctuations
  real(kind=dp), intent(inout) :: &
    xt, yt                          ! particle positions on grid
  real ::                         &
    xts,yts,                      & ! local 'real' copy of the particle position
    xtn,ytn,                      & ! nested particle position
    weight,                       & ! transition above the tropopause
    dz,dz1,dz2,                   & ! values used for interpolating between vertical levels
    xlon,ylat,xpol,ypol,          & ! temporarily storing new particle positions
    gridsize,cosfact,             & ! used to compute new positions of particles
    ru,rv,rw,                     & ! used for computing turbulence
    dt,                           & ! real(ldt)
    ux,vy,                        & ! random turbulent velocities above PBL
    tropop,                       & ! height of troposphere
    prob(maxspec),                & ! dry deposition ground absorption probability
    dxsave,dysave,                & ! accumulated displacement in long and lat
    dawsave,dcwsave,              & ! accumulated displacement in wind directions
    uold,vold,wold,woldeta,       & ! 
    r,rs,                         & ! mesoscale related
    vdepo(maxspec),               & ! deposition velocities for all species
    h1(2),                        & ! mixing height
    rhoa,                         & ! air density, used in CBL
    rhograd,                      & ! vertical gradient of the air density, used in CBL
    delz,                         & ! change in vertical position due to turbulence
    dtf,rhoaux,dtftlw,ath,bth,    & ! CBL related
    ptot_lhh,Q_lhh,phi_lhh,       & ! CBL related
    old_wp_buf,dcas,dcas1,        & ! CBL related
    del_test,                     & ! CBL related
    uxscale,wpscale,              & ! factor used in calculating turbulent perturbations above PBL
    ztemp,                        & ! temporarily storing z position
    settling = 0.                   ! settling velo

  ! openmp change
  save idummy
!$OMP THREADPRIVATE(idummy)  
!$    if (idummy.eq.-7) then
!$      thread = OMP_GET_THREAD_NUM()
!$      idummy = idummy - thread
!$    endif
  ! openmp change end 

  nstop=0
  do i=1,nmixz
    indzindicator(i)=.true.
  end do

  indzeta=nz-1
  indzpeta=nz
  do i=2,nz
    if (wheight(i).lt.zteta) then
      indzeta=i-1
      indzpeta=i
      exit
    endif
  end do

  induv=nz-1
  indpuv=nz
  do i=2,nz
    if (uvheight(i).lt.zteta) then
      induv=i-1
      indpuv=i
      exit
    endif
  end do
  
  if (DRYDEP) then    ! reset probability for deposition
    do ks=1,nspec
      depoindicator(ks)=.true.
      prob(ks)=0.
    end do
  endif

  dxsave=0.           ! reset position displacements
  dysave=0.           ! due to mean wind
  dawsave=0.          ! and turbulent wind
  dcwsave=0.

  itimec=itime

  nrand=int(ran3(idummy)*real(maxrand-1))+1

  ! Determine whether lat/long grid or polarstereographic projection
  ! is to be used
  ! Furthermore, determine which nesting level to be used
  !*****************************************************************

  if (nglobal.and.(yt.gt.switchnorthg)) then
    ngrid=-1
  else if (sglobal.and.(yt.lt.switchsouthg)) then
    ngrid=-2
  else
    ngrid=0
    do j=numbnests,1,-1
      if ((xt.gt.xln(j)+eps).and.(xt.lt.xrn(j)-eps).and. &
           (yt.gt.yln(j)+eps).and.(yt.lt.yrn(j)-eps)) then
        ngrid=j
        exit
      endif
    end do
  endif

  !***************************
  ! Interpolate necessary data
  !***************************

  if (abs(itime-memtime(1)).lt.abs(itime-memtime(2))) then
    memindnext=1
  else
    memindnext=2
  endif

  ! Determine nested grid coordinates
  !**********************************

  if (ngrid.gt.0) then
    xtn=(xt-xln(ngrid))*xresoln(ngrid)
    ytn=(yt-yln(ngrid))*yresoln(ngrid)
    ix=int(xtn)
    jy=int(ytn)
    nix=nint(xtn)
    njy=nint(ytn)
  else
    ix=int(xt)
    jy=int(yt)
    nix=nint(xt)
    njy=nint(yt)
  endif
  ixp=ix+1
  jyp=jy+1

 ! Determine the lower left corner and its distance to the current position
  !*************************************************************************

  ddx=xt-real(ix)
  ddy=yt-real(jy)
  rddx=1.-ddx
  rddy=1.-ddy
  p1=rddx*rddy
  p2=ddx*rddy
  p3=rddx*ddy
  p4=ddx*ddy

 ! Calculate variables for time interpolation
  !*******************************************

  dt1=real(itime-memtime(1))
  dt2=real(memtime(2)-itime)
  dtt=1./(dt1+dt2)

! eso: Temporary fix for particle exactly at north pole
  if (jyp >= nymax) then
    ! write(*,*) 'WARNING: advance.f90 jyp >= nymax. xt,yt:',xt,yt
    jyp=jyp-1
  end if

  if (jyp >= nymax) then
    jyp=jyp-1
  end if
  ! Compute maximum mixing height around particle position
  !*******************************************************
  
  ! Convert z(eta) to z(m) for the turbulence scheme, w(m/s) 
  ! is computed in verttransform_ecmwf.f90
  call zeta_to_z(itime,xt,yt,zteta,zt)

  h=0.
  if (ngrid.le.0) then
    do k=1,2
       mind=memind(k) ! eso: compatibility with 3-field version
       if (interpolhmix) then
             h1(k)=p1*hmix(ix ,jy ,1,mind) &
                 + p2*hmix(ixp,jy ,1,mind) &
                 + p3*hmix(ix ,jyp,1,mind) &
                 + p4*hmix(ixp,jyp,1,mind)
        else
          do j=jy,jyp
            do i=ix,ixp
               if (hmix(i,j,1,mind).gt.h) h=hmix(i,j,1,mind)
            end do
          end do
        endif
    end do
    tropop=tropopause(nix,njy,1,1)
  else
    do k=1,2
      mind=memind(k)
      do j=jy,jyp
        do i=ix,ixp
          if (hmixn(i,j,1,mind,ngrid).gt.h) h=hmixn(i,j,1,mind,ngrid)
        end do
      end do
    end do
    tropop=tropopausen(nix,njy,1,1,ngrid)
  endif

  if (interpolhmix) h=(h1(1)*dt2+h1(2)*dt1)*dtt 
  zeta=zt/h

  !*************************************************************
  ! If particle is in the PBL, interpolate once and then make a
  ! time loop until end of interval is reached
  !*************************************************************
  ! In the PBL we use meters instead of eta coordinates for the vertical transport
  if (zeta.le.1.) then

  ! BEGIN TIME LOOP
  !================
    loop=0
    pbl_loop : do
      loop=loop+1
      if (method.eq.1) then
        ldt=min(ldt,abs(lsynctime-itimec+itime))
        itimec=itimec+ldt*ldirect
      else
        ldt=abs(lsynctime)
        itimec=itime+lsynctime
      endif
      dt=real(ldt)

      zeta=zt/h


      if (loop.eq.1) then
        if (ngrid.le.0) then
          xts=real(xt)
          yts=real(yt)
          call interpol_all(itime,xts,yts,zt,zteta)
        else
          call interpol_all_nests(itime,xtn,ytn,zt)
        endif

      else


  ! Determine the level below the current position for u,v,rho
  !***********************************************************

        do i=2,nz
          if (height(i).gt.zt) then
            indz=i-1
            indzp=i
            exit
          endif
        end do

        indzeta=nz-1
        indzpeta=nz
        do i=2,nz
          if (wheight(i).lt.zteta) then
            indzeta=i-1
            indzpeta=i
            exit
          endif
        end do

        induv=nz-1
        indpuv=nz
        do i=2,nz
          if (uvheight(i).lt.zteta) then
            induv=i-1
            indpuv=i
            exit
          endif
        end do

  ! If one of the levels necessary is not yet available,
  ! calculate it
  !*****************************************************

        do i=indzeta,indzpeta
          if (indzindicator(i)) then
            if (ngrid.le.0) then
              call interpol_misslev(i)
            else
              call interpol_misslev_nests(i)
            endif
          endif
        end do
      endif


  ! Vertical interpolation of u,v,w,rho and drhodz
  !***********************************************

  ! Vertical distance to the level below and above current position
  ! both in terms of (u,v) and (w) fields
  !****************************************************************
      ! Keep w(m) for computing the turbulence
      dz=1./(height(indzp)-height(indz))
      dz1=(zt-height(indz))*dz
      dz2=(height(indzp)-zt)*dz
      w=dz1*wprof(indzp)+dz2*wprof(indz)

      dz=1./(uvheight(indpuv)-uvheight(induv))
      dz1=(zteta-uvheight(induv))*dz
      dz2=(uvheight(indpuv)-zteta)*dz
      u=dz1*uprof(indpuv)+dz2*uprof(induv)
      v=dz1*vprof(indpuv)+dz2*vprof(induv)
      rhoa=dz1*rhoprof(indpuv)+dz2*rhoprof(induv)
      rhograd=dz1*rhogradprof(indpuv)+dz2*rhogradprof(induv)

      ! Compute w(eta) and use where possible
      dz=1./(wheight(indzpeta)-wheight(indzeta))
      dz1=(zteta-wheight(indzeta))*dz
      dz2=(wheight(indzpeta)-zteta)*dz
      weta=dz1*wprofeta(indzpeta)+ dz2*wprofeta(indzeta)

  ! Compute the turbulent disturbances
  ! Determine the sigmas and the timescales
  !****************************************

      if (turbswitch) then
        call hanna(zt)
      else
        call hanna1(zt)
      endif

  !*****************************************
  ! Determine the new diffusivity velocities
  !*****************************************

  ! Horizontal components
  !**********************

      if (nrand+1.gt.maxrand) nrand=1
      if (dt/tlu.lt..5) then
        up=(1.-dt/tlu)*up+rannumb(nrand)*sigu*sqrt(2.*dt/tlu)
      else
        ru=exp(-dt/tlu)
        up=ru*up+rannumb(nrand)*sigu*sqrt(1.-ru**2)
      endif
      if (dt/tlv.lt..5) then
        vp=(1.-dt/tlv)*vp+rannumb(nrand+1)*sigv*sqrt(2.*dt/tlv)
      else
        rv=exp(-dt/tlv)
        vp=rv*vp+rannumb(nrand+1)*sigv*sqrt(1.-rv**2)
      endif
      nrand=nrand+2


      if (nrand+ifine.gt.maxrand) nrand=1
      rhoaux=rhograd/rhoa
      dtf=dt*fine

      dtftlw=dtf/tlw

  ! Loop over ifine short time steps for vertical component
  !********************************************************

      do i=1,ifine

  ! Determine the drift velocity and density correction velocity
  !*************************************************************

        if (turbswitch) then
          if (dtftlw.lt..5) then
  !*************************************************************
  !************** CBL options added by mc see routine cblf90 ***
            ! LB needs to be checked if this works with openmp
            if (cblflag.eq.1) then  !modified by mc
              if (-h/ol.gt.5) then  !modified by mc
              !if (ol.lt.0.) then   !modified by mc  
              !if (ol.gt.0.) then   !modified by mc : for test
                  !print  *,zt,wp,ath,bth,tlw,dtf,'prima'
                  flagrein=0
                  nrand=nrand+1
                  old_wp_buf=wp
                  call cbl(wp,zt,ust,wst,h,rhoa,rhograd,sigw,dsigwdz,tlw,ptot_lhh,Q_lhh,phi_lhh,ath,bth,ol,flagrein) !inside the routine for inverse time
                  wp=(wp+ath*dtf+bth*rannumb(nrand)*sqrt(dtf))*real(icbt) 
                  ! wp=(wp+ath*dtf+bth*gasdev2(mydum)*sqrt(dtf))*real(icbt) 
                  delz=wp*dtf
                  if (flagrein.eq.1) then
                      call re_initialize_particle(zt,ust,wst,h,sigw,old_wp_buf,nrand,ol)
                      wp=old_wp_buf
                      delz=wp*dtf
                      nan_count=nan_count+1
                  end if
                  !print  *,zt,wp,ath,bth,tlw,dtf,rannumb(nrand+i),icbt
                  !pause                  
              else 
                  nrand=nrand+1
                  old_wp_buf=wp
                  ath=-wp/tlw+sigw*dsigwdz+wp*wp/sigw*dsigwdz+sigw*sigw/rhoa*rhograd  !1-note for inverse time should be -wp/tlw*ldirect+... calculated for wp=-wp
                                                                                      !2-but since ldirect =-1 for inverse time and this must be calculated for (-wp) and
                                                                                      !3-the gaussian pdf is symmetric (i.e. pdf(w)=pdf(-w) ldirect can be discarded
                  bth=sigw*rannumb(nrand)*sqrt(2.*dtftlw)
                  wp=(wp+ath*dtf+bth)*real(icbt)  
                  delz=wp*dtf
                  del_test=(1.-wp)/wp !catch infinity value
                  if (isnan(wp).or.isnan(del_test)) then 
                      nrand=nrand+1                      
                      wp=sigw*rannumb(nrand)
                      delz=wp*dtf
                      nan_count2=nan_count2+1
                      !print *,'NaN coutner equal to:', nan_count,'reduce ifine if this number became a non-negligible fraction of the particle number'
                  end if  
              end if
  !******************** END CBL option *******************************            
  !*******************************************************************            
            else
                 wp=((1.-dtftlw)*wp+rannumb(nrand+i)*sqrt(2.*dtftlw) &
                 +dtf*(dsigwdz+rhoaux*sigw))*real(icbt) 
                 delz=wp*sigw*dtf
            end if
          else
            rw=exp(-dtftlw)
            wp=(rw*wp+rannumb(nrand+i)*sqrt(1.-rw**2) &
                 +tlw*(1.-rw)*(dsigwdz+rhoaux*sigw))*real(icbt)
            delz=wp*sigw*dtf
          endif
          
        else
          rw=exp(-dtftlw)
          wp=(rw*wp+rannumb(nrand+i)*sqrt(1.-rw**2)*sigw &
               +tlw*(1.-rw)*(dsigw2dz+rhoaux*sigw**2))*real(icbt)
          delz=wp*dtf
        endif

        if (turboff) then
!sec switch off turbulence
          up=0.0
          vp=0.0
          wp=0.0
          delz=0.
        endif

  !****************************************************
  ! Compute turbulent vertical displacement of particle
  !****************************************************

        if (abs(delz).gt.h) delz=mod(delz,h)

  ! Determine if particle transfers to a "forbidden state" below the ground
  ! or above the mixing height
  !************************************************************************

        if (delz.lt.-zt) then         ! reflection at ground
          icbt=-1
          zt=-zt-delz
        else if (delz.gt.(h-zt)) then ! reflection at h
          icbt=-1
          zt=-zt-delz+2.*h
        else                         ! no reflection
          icbt=1
          zt=zt+delz
        endif

        if (i.ne.ifine) then
          zeta=zt/h
          call hanna_short(zt)
        endif

      end do
      if (cblflag.ne.1) nrand=nrand+i

  ! Determine time step for next integration
  !*****************************************

      if (turbswitch) then
        ldt=int(min(tlw,h/max(2.*abs(wp*sigw),1.e-5), &
             0.5/abs(dsigwdz))*ctl)
      else
        ldt=int(min(tlw,h/max(2.*abs(wp),1.e-5))*ctl)
      endif
      ldt=max(ldt,mintime)


  ! If particle represents only a single species, add gravitational settling
  ! velocity. The settling velocity is zero for gases, or if particle
  ! represents more than one species
  !*************************************************************************

      if (mdomainfill.eq.0) then
        if (lsettling) then
          do nsp=1,nspec
            if (xmass(nrelpoint,nsp).gt.eps3) exit
          end do
          if (nsp.gt.nspec) then
            nsp=nspec
          end if
          if (density(nsp).gt.0.) then
            call get_settling(itime,real(xt),real(yt),zt,nsp,settling)  !bugfix
            w=w+settling
          end if
        end if
      endif

  ! Horizontal displacements during time step dt are small real values compared
  ! to the position; adding the two, would result in large numerical errors.
  ! Thus, displacements are accumulated during lsynctime and are added to the
  ! position at the end
  !****************************************************************************

      dxsave=dxsave+u*dt
      dysave=dysave+v*dt
      dawsave=dawsave+up*dt
      dcwsave=dcwsave+vp*dt
      ! How can I change the w to w(eta) efficiently?
      zt=zt+w*dt*real(ldirect)

      ! HSO/AL: Particle managed to go over highest level -> interpolation error in goto 700
      !          alias interpol_wind (division by zero)
      if (zt.ge.height(nz)) zt=height(nz)-100.*eps

      if (zt.gt.h) then
        call z_to_zeta(itime,xt,yt,zt,zteta)
        if (itimec.eq.itime+lsynctime) goto 99
        goto 700    ! complete the current interval above PBL
      endif
      
  ! Determine probability of deposition
  !************************************

      if ((DRYDEP).and.(zt.lt.2.*href)) then
        do ks=1,nspec
          if (DRYDEPSPEC(ks)) then
            if (depoindicator(ks)) then
              if (ngrid.le.0) then
                call interpol_vdep(ks,vdepo(ks))
              else
                call interpol_vdep_nests(ks,vdepo(ks))
              endif
            endif
  ! correction by Petra Seibert, 10 April 2001
  !   this formulation means that prob(n) = 1 - f(0)*...*f(n)
  !   where f(n) is the exponential term
               prob(ks)=1.+(prob(ks)-1.)* &
                    exp(-vdepo(ks)*abs(dt)/(2.*href))
          endif
        end do
      endif

      if (zt.lt.0.) zt=min(h-eps2,-1.*zt)    ! if particle below ground -> reflection


      if (itimec.eq.(itime+lsynctime)) then
        usig=0.5*(usigprof(indpuv)+usigprof(induv))
        vsig=0.5*(vsigprof(indpuv)+vsigprof(induv))
        wsig=0.5*(wsigprof(indzp)+wsigprof(indz))
        wsigeta=0.5*(wsigprofeta(indzpeta)+wsigprofeta(indzeta))
        ! Converting the z position that changed through turbulence motions to eta coords
        call z_to_zeta(itime,xt,yt,zt,zteta)
        goto 99  ! finished
      endif
    end do pbl_loop

  ! END TIME LOOP
  !==============
  endif

  !**********************************************************
  ! For all particles that are outside the PBL, make a single
  ! time step. Only horizontal turbulent disturbances are
  ! calculated. Vertical disturbances are reset.
  !**********************************************************


  ! Interpolate the wind
  !*********************

700   continue
  if (ngrid.le.0) then
    xts=real(xt)
    yts=real(yt)
    call interpol_wind(itime,xts,yts,zt,zteta,pp)
  else
    call interpol_wind_nests(itime,xtn,ytn,zt)
  endif

  ! Compute everything for above the PBL

  ! Assume constant, uncorrelated, turbulent perturbations
  ! In the stratosphere, use a small vertical diffusivity d_strat,
  ! in the troposphere, use a larger horizontal diffusivity d_trop.
  ! Turbulent velocity scales are determined based on sqrt(d_trop/dt)
  !******************************************************************

  ldt=abs(lsynctime-itimec+itime)
  dt=real(ldt)

  if (zt.lt.tropop) then  ! in the troposphere
    uxscale=sqrt(2.*d_trop/dt)
    if (nrand+1.gt.maxrand) nrand=1
    ux=rannumb(nrand)*uxscale
    vy=rannumb(nrand+1)*uxscale
    nrand=nrand+2
    wp=0.
  else if (zt.lt.tropop+1000.) then     ! just above the tropopause: make transition
    weight=(zt-tropop)/1000.
    uxscale=sqrt(2.*d_trop/dt*(1.-weight))
    if (nrand+2.gt.maxrand) nrand=1
    ux=rannumb(nrand)*uxscale
    vy=rannumb(nrand+1)*uxscale
    wpscale=sqrt(2.*d_strat/dt*weight)
    wp=rannumb(nrand+2)*wpscale+d_strat/1000.
    nrand=nrand+3
  else                 ! in the stratosphere
    if (nrand.gt.maxrand) nrand=1
    ux=0.
    vy=0.
    wpscale=sqrt(2.*d_strat/dt)
    wp=rannumb(nrand)*wpscale
    nrand=nrand+1
  endif

  if (turboff) then
!sec switch off turbulence
    ux=0.0
    vy=0.0
    wp=0.0
  endif

  ! If particle represents only a single species, add gravitational settling
  ! velocity. The settling velocity is zero for gases
  !*************************************************************************
  ! Does not work in eta coordinates yet
  if (mdomainfill.eq.0) then
    if (lsettling) then
      do nsp=1,nspec
        if (xmass(nrelpoint,nsp).gt.eps3) exit
      end do
      if (nsp.gt.nspec) then
        nsp=nspec
      end if
      ! LB needs to be checked if this works with openmp and change to eta coords
      if (density(nsp).gt.0.) then
        call get_settling(itime,real(xt),real(yt),zt,nsp,settling)  !bugfix
        w=w+settling
        zt=zt+settling*dt*real(ldirect)
      end if
    endif
  end if


  ! Calculate position at time step itime+lsynctime
  !************************************************
  dxsave=dxsave+(u+ux)*dt
  dysave=dysave+(v+vy)*dt
  zt=zt+(wp)*dt*real(ldirect)
  if (zt.lt.0.) zt=min(h-eps2,-1.*zt)    ! if particle below ground -> reflection
  call z_to_zeta(itime,xt,yt,zt,zteta)

  zteta=zteta+(weta)*dt*real(ldirect)

  if (zteta.ge.1.) zteta=1.-(zteta-1.)
  if (zteta.eq.1.) zteta=zteta-eps_eta
  ! if (zteta.ge.uvheight(2)) zteta=uvheight(2) -(zteta - uvheight(2))


99   continue



  !****************************************************************
  ! Add mesoscale random disturbances
  ! This is done only once for the whole lsynctime interval to save
  ! computation time
  !****************************************************************


  ! Mesoscale wind velocity fluctuations are obtained by scaling
  ! with the standard deviation of the grid-scale winds surrounding
  ! the particle location, multiplied by a factor turbmesoscale.
  ! The autocorrelation time constant is taken as half the
  ! time interval between wind fields
  !****************************************************************

  r=exp(-2.*real(abs(lsynctime))/real(lwindinterv))
  rs=sqrt(1.-r**2)
  if (nrand+2.gt.maxrand) nrand=1
  usigold=r*usigold+rs*rannumb(nrand)*usig*turbmesoscale
  vsigold=r*vsigold+rs*rannumb(nrand+1)*vsig*turbmesoscale
  wsigold=r*wsigold+rs*rannumb(nrand+2)*wsigeta*turbmesoscale

  dxsave=dxsave+usigold*real(lsynctime)
  dysave=dysave+vsigold*real(lsynctime)
  zteta=zteta+wsigold*real(lsynctime)

  if (zteta.ge.1.) zteta=1.-(zteta-1.)
  if (zteta.eq.1.) zteta=zteta-eps_eta
  ! if (zteta.ge.uvheight(2)) zteta=uvheight(2) -(zteta - uvheight(2))

  !*************************************************************
  ! Transform along and cross wind components to xy coordinates,
  ! add them to u and v, transform u,v to grid units/second
  ! and calculate new position
  !*************************************************************

  call windalign(dxsave,dysave,dawsave,dcwsave,ux,vy)
  dxsave=dxsave+ux   ! comment by mc: comment this line to stop the particles horizontally for test reasons 
  dysave=dysave+vy
  if (ngrid.ge.0) then
    cosfact=dxconst/cos((yt*dy+ylat0)*pi180)
    xt=xt+real(dxsave*cosfact*real(ldirect),kind=dp)
    yt=yt+real(dysave*dyconst*real(ldirect),kind=dp)
  else if (ngrid.eq.-1) then      ! around north pole
    xlon=xlon0+real(xt)*dx                                !comment by mc: compute old particle position
    ylat=ylat0+real(yt)*dy
    call cll2xy(northpolemap,ylat,xlon,xpol,ypol)   !convert old particle position in polar stereographic
    gridsize=1000.*cgszll(northpolemap,ylat,xlon)   !calculate size in m of grid element in polar stereographic coordinate
    dxsave=dxsave/gridsize                          !increment from meter to grdi unit
    dysave=dysave/gridsize
    xpol=xpol+dxsave*real(ldirect)                  !position in grid unit polar stereographic
    ypol=ypol+dysave*real(ldirect)
    call cxy2ll(northpolemap,xpol,ypol,ylat,xlon)  !convert to lat long coordinate
    xt=real((xlon-xlon0)/dx,kind=dp)                             !convert to grid units in lat long coordinate, comment by mc
    yt=real((ylat-ylat0)/dy,kind=dp)
  else if (ngrid.eq.-2) then    ! around south pole
    xlon=xlon0+real(xt)*dx
    ylat=ylat0+real(yt)*dy
    call cll2xy(southpolemap,ylat,xlon,xpol,ypol)
    gridsize=1000.*cgszll(southpolemap,ylat,xlon)
    dxsave=dxsave/gridsize
    dysave=dysave/gridsize
    xpol=xpol+dxsave*real(ldirect)
    ypol=ypol+dysave*real(ldirect)
    call cxy2ll(southpolemap,xpol,ypol,ylat,xlon)
    xt=real((xlon-xlon0)/dx,kind=dp)
    yt=real((ylat-ylat0)/dy,kind=dp)
  endif

  ! If global data are available, use cyclic boundary condition
  !************************************************************
  if (xglobal) then
    if (xt.ge.real(nxmin1,kind=dp)) xt=xt-real(nxmin1,kind=dp)
    if (xt.lt.0.) xt=xt+real(nxmin1,kind=dp)
    if (xt.le.real(eps,kind=dp)) xt=real(eps,kind=dp)
    if (abs(xt-real(nxmin1)).le.eps) xt=real(nxmin1-eps,kind=dp)
  endif

  ! HSO/AL: Prevent particles from disappearing at the pole
  !******************************************************************

  if ( yt.lt.0. ) then
    xt=mod(xt+180.,360.)
    yt=-yt
  else if ( yt.gt.real(nymin1,kind=dp) ) then
    xt=mod(xt+180.,360.)
    yt=2.*real(nymin1,kind=dp)-yt
  endif

  ! Check position: If trajectory outside model domain, terminate it
  !*****************************************************************

  if ((xt.lt.0.).or.(xt.ge.real(nxmin1)).or.(yt.lt.0.).or. &
       (yt.gt.real(nymin1))) then
    nstop=3
    return
  endif

  ! If particle above highest model level, set it back into the domain
  !*******************************************************************

  if (zteta.le.uvheight(nz)) zteta=uvheight(nz)+eps_eta


  !************************************************************************
  ! Now we could finish, as this was done in FLEXPART versions up to 4.0.
  ! However, truncation errors of the advection can be significantly
  ! reduced by doing one iteration of the Petterssen scheme, if this is
  ! possible.
  ! Note that this is applied only to the grid-scale winds, not to
  ! the turbulent winds.
  !************************************************************************

  ! The Petterssen scheme can only applied with long time steps (only then u
  ! is the "old" wind as required by the scheme); otherwise do nothing
  !*************************************************************************

  if (ldt.ne.abs(lsynctime)) return

  ! The Petterssen scheme can only be applied if the ending time of the time step
  ! (itime+ldt*ldirect) is still between the two wind fields held in memory;
  ! otherwise do nothing
  !******************************************************************************

  if (abs(itime+ldt*ldirect).gt.abs(memtime(2))) return

  ! Apply it also only if starting and ending point of current time step are on
  ! the same grid; otherwise do nothing
  !*****************************************************************************
  if (nglobal.and.(yt.gt.switchnorthg)) then
    ngr=-1
  else if (sglobal.and.(yt.lt.switchsouthg)) then
    ngr=-2
  else
    ngr=0
    do j=numbnests,1,-1
      if ((xt.gt.xln(j)+eps).and.(xt.lt.xrn(j)-eps).and. &
           (yt.gt.yln(j)+eps).and.(yt.lt.yrn(j)-eps)) then
        ngr=j
        exit
      endif
    end do
  endif

  if (ngr.ne.ngrid) return

  ! Determine nested grid coordinates
  !**********************************

  if (ngrid.gt.0) then
    xtn=(xt-xln(ngrid))*xresoln(ngrid)
    ytn=(yt-yln(ngrid))*yresoln(ngrid)
    ix=int(xtn)
    jy=int(ytn)
  else
    ix=int(xt)
    jy=int(yt)
  endif
  ixp=ix+1
  jyp=jy+1

  ! Memorize the old wind
  !**********************

  uold=u
  vold=v
  !wold=w

  woldeta=weta
  ! Interpolate wind at new position and time
  !******************************************

  if (ngrid.le.0) then
    xts=real(xt)
    yts=real(yt)
    if (ix.gt.720) then
      write(*,*) 'ERROR! ', pp,xt,yt
    endif
    call interpol_wind_short(itime+ldt*ldirect,xts,yts,zt,zteta)
  else
    call interpol_wind_short_nests(itime+ldt*ldirect,xtn,ytn,zt)
  endif

  if (mdomainfill.eq.0) then
    if (lsettling) then
      do nsp=1,nspec
        if (xmass(nrelpoint,nsp).gt.eps3) exit
      end do
      if (nsp.gt.nspec) then
        nsp=nspec
      end if
      if (density(nsp).gt.0.) then
        call zeta_to_z(itime,xt,yt,zteta,zt)
        call get_settling(itime+ldt,real(xt),real(yt),zt,nsp,settling) !bugfix
        call z_to_zeta(itime,xt,yt,zt+settling*real(ldt*ldirect),ztemp)
        weta=weta+(ztemp-zteta)/real(ldt*ldirect)
      end if
    endif
  end if


  ! Determine the difference vector between new and old wind
  ! (use half of it to correct position according to Petterssen)
  !*************************************************************

  u=(u-uold)/2.
  v=(v-vold)/2.
  !w=(w-wold)/2.

  weta=(weta-woldeta)/2.

  ! Finally, correct the old position
  !**********************************
  zteta=zteta+weta*real(ldt*ldirect)

  if (zteta.ge.1.) zteta=1.-(zteta-1.)
  if (zteta.eq.1.) zteta=zteta-eps_eta
  ! if (zteta.ge.uvheight(2)) zteta=uvheight(2) -(zteta - uvheight(2))

  if (ngrid.ge.0) then
    cosfact=dxconst/cos((real(yt)*dy+ylat0)*pi180)
    xt=xt+real(u*cosfact*real(ldt*ldirect),kind=dp)
    yt=yt+real(v*dyconst*real(ldt*ldirect),kind=dp)
  else if (ngrid.eq.-1) then      ! around north pole
    xlon=xlon0+real(xt)*dx
    ylat=ylat0+real(yt)*dy
    call cll2xy(northpolemap,ylat,xlon,xpol,ypol)
    gridsize=1000.*cgszll(northpolemap,ylat,xlon)
    u=u/gridsize
    v=v/gridsize
    xpol=xpol+u*real(ldt*ldirect)
    ypol=ypol+v*real(ldt*ldirect)
    call cxy2ll(northpolemap,xpol,ypol,ylat,xlon)
    xt=real((xlon-xlon0)/dx,kind=dp)
    yt=real((ylat-ylat0)/dy,kind=dp)
  else if (ngrid.eq.-2) then    ! around south pole
    xlon=xlon0+real(xt)*dx
    ylat=ylat0+real(yt)*dy
    call cll2xy(southpolemap,ylat,xlon,xpol,ypol)
    gridsize=1000.*cgszll(southpolemap,ylat,xlon)
    u=u/gridsize
    v=v/gridsize
    xpol=xpol+u*real(ldt*ldirect)
    ypol=ypol+v*real(ldt*ldirect)
    call cxy2ll(southpolemap,xpol,ypol,ylat,xlon)
    xt=real((xlon-xlon0)/dx,kind=dp)
    yt=real((ylat-ylat0)/dy,kind=dp)
  endif

  ! If global data are available, use cyclic boundary condition
  !************************************************************

  if (xglobal) then
    if (xt.ge.real(nxmin1,kind=dp)) xt=xt-real(nxmin1,kind=dp)
    if (xt.lt.0.) xt=xt+real(nxmin1,kind=dp)
    if (xt.le.eps) xt=real(eps,kind=dp)
    if (abs(xt-real(nxmin1,kind=dp)).le.eps) xt=real(nxmin1-eps,kind=dp)
  endif

  ! HSO/AL: Prevent particles from disappearing at the pole
  !******************************************************************
  if ( yt.lt.0. ) then
    xt=mod(xt+180.,360.)
    yt=-yt
  else if ( yt.gt.real(nymin1,kind=dp) ) then
    xt=mod(xt+180.,360.)
    yt=2.*real(nymin1,kind=dp)-yt
  endif

  ! Check position: If trajectory outside model domain, terminate it
  !*****************************************************************
  if ((xt.lt.0.).or.(xt.ge.real(nxmin1,kind=dp)).or.(yt.lt.0.).or. &
       (yt.gt.real(nymin1,kind=dp))) then
    nstop=3
    return
  endif

  ! If particle above highest model level, set it back into the domain
  !*******************************************************************

  if (zteta.le.uvheight(nz)) zteta=uvheight(nz)+eps_eta


end subroutine advance

