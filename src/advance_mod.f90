! SPDX-FileCopyrightText: FLEXPART 1998-2019, see flexpart_license.txt
! SPDX-License-Identifier: GPL-3.0-or-later
module advance_mod
  use point_mod
  use par_mod
  use com_mod
  use interpol_mod
  use cmapf_mod
  use random_mod, only: ran3
  use coordinates_ecmwf
  use particle_mod
  use turbulence_mod
  use settling_mod

  implicit none 
    real, parameter ::              &
      eps2=1.e-9,                   &
      eps3=tiny(1.0),               &
      eps_eta=1.e-4
    real ::                         &
      eps
  private :: advance_abovePBL,advance_PBL,advance_PettersonCorrection,&
    advance_updateXY,advance_adjusttopheight
contains
  
subroutine advance(itime,ipart)
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
  ! xt,yt,zt           Particle position                                       *
  !                                                                            *
  !*****************************************************************************

  ! openmp change
  use omp_lib, only: OMP_GET_THREAD_NUM
  ! openmp change end

  implicit none
  integer, intent(in) ::          &
    itime,                        & ! time index
    ipart                           ! particle index
  integer ::                      &
    itimec,                       &
    i,j,k,                        & ! loop variables
    nrand,                        & ! random number used for turbulence
    memindnext,                   & ! seems useless
    ngr,                          & ! temporary new grid index of moved particle
    nsp,                          & ! loop variables for number of species
    thread,                       & ! number of openmp threads (probably can be removed)
    idummy = -7                     ! used in random number routines
  real ::                         &
    ux,vy,                        & ! random turbulent velocities above PBL
    tropop,                       & ! height of troposphere
    dxsave,dysave,                & ! accumulated displacement in long and lat
    dawsave,dcwsave                 ! accumulated displacement in wind directions
  logical ::                      &
    abovePBL                        ! flag that will be set to 'true' if computation needs to be completed above PBL
  !type(particle) :: part
  ! openmp change
  save idummy
!$OMP THREADPRIVATE(idummy)  
!$    if (idummy.eq.-7) then
!$      thread = OMP_GET_THREAD_NUM()
!$      idummy = idummy - thread
!$    endif
  ! openmp change end 

  eps=nxmax/3.e5

  part(ipart)%nstop=.false.
  do i=1,nmixz
    indzindicator(i)=.true.
  end do
  
  if (DRYDEP) then    ! reset probability for deposition
    depoindicator=.true.
    part(ipart)%prob=0.
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
  ! call find_ngrid(part(ipart)%xlon,part(ipart)%ylat)
  if (nglobal.and.(part(ipart)%ylat.gt.switchnorthg)) then
    ngrid=-1
  else if (sglobal.and.(part(ipart)%ylat.lt.switchsouthg)) then
    ngrid=-2
  else
    ngrid=0
    do j=numbnests,1,-1
      if ((part(ipart)%xlon.gt.xln(j)+eps).and.(part(ipart)%xlon.lt.xrn(j)-eps).and. &
           (part(ipart)%ylat.gt.yln(j)+eps).and.(part(ipart)%ylat.lt.yrn(j)-eps)) then
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
  ! Determine the lower left corner and its distance to the current position
  ! Calculate variables for time interpolation
  !*******************************************
  call initialise_interpol_mod(itime,real(part(ipart)%xlon),real(part(ipart)%ylat),&
    real(part(ipart)%z),real(part(ipart)%zeta))

  ! Compute maximum mixing height around particle position
  !*******************************************************
  
  ! Convert z(eta) to z(m) for the turbulence scheme, w(m/s) 
  ! is computed in verttransform_ecmwf.f90

  call update_zeta_to_z(itime,ipart)

  ! Compute the height of the troposphere and the PBL at the x-y location of the particle
  call interpol_htropo_hmix(tropop,h)
  zeta=real(part(ipart)%z)/h

  !*************************************************************
  ! If particle is in the PBL, interpolate once and then make a
  ! time loop until end of interval is reached
  !*************************************************************
  ! In the PBL we use meters instead of eta coordinates for the vertical transport
  abovePBL=.true.
  if (zeta.le.1.) then
    abovePBL=.false.
    call advance_PBL(itime,itimec,&
      dxsave,dysave,dawsave,dcwsave,abovePBL,nrand,ipart)
  endif 

  !**********************************************************
  ! For all particles that are outside the PBL, make a single
  ! time step. Only horizontal turbulent disturbances are
  ! calculated. Vertical disturbances are reset.
  !**********************************************************

  ! Interpolate the wind
  !*********************
  if (abovePBL) then 
    call advance_abovePBL(itime,itimec,&
      dxsave,dysave,ux,vy,tropop,nrand,ipart)
  endif ! Above PBL computation

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
  if (.not. turboff) then ! mesoscale turbulence is found to give issues, so turned off
    ! call turbulence_mesoscale(nrand,dxsave,dysave,ipart,usig,vsig,wsig,wsigeta,eps_eta)

    !*************************************************************
    ! Transform along and cross wind components to xy coordinates,
    ! add them to u and v, transform u,v to grid units/second
    ! and calculate new position
    !*************************************************************

    call windalign(dxsave,dysave,dawsave,dcwsave,ux,vy)
    dxsave=dxsave+ux   ! comment by mc: comment this line to stop the particles horizontally for test reasons 
    dysave=dysave+vy
  endif

  call advance_updateXY(dxsave,dysave,ipart)

  ! If particle above highest model level, set it back into the domain
  !*******************************************************************
  call advance_adjusttopheight(ipart)
  
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

  if (part(ipart)%idt.ne.abs(lsynctime)) return

  ! The Petterssen scheme can only be applied if the ending time of the time step
  ! (itime+ldt*ldirect) is still between the two wind fields held in memory;
  ! otherwise do nothing
  !******************************************************************************

  if (abs(itime+part(ipart)%idt*ldirect).gt.abs(memtime(2))) return

  ! Apply it also only if starting and ending point of current time step are on
  ! the same grid; otherwise do nothing
  !*****************************************************************************
  ! ngr = ngrid 
  ! call find_ngrid(part(ipart)%xlon,part(ipart)%ylat)

  if (nglobal.and.(real(part(ipart)%ylat).gt.switchnorthg)) then
    ngr=-1
  else if (sglobal.and.(real(part(ipart)%ylat).lt.switchsouthg)) then
    ngr=-2
  else
    ngr=0
    do j=numbnests,1,-1
      if ((real(part(ipart)%xlon).gt.xln(j)+eps).and.(real(part(ipart)%xlon).lt.xrn(j)-eps).and. &
           (real(part(ipart)%ylat).gt.yln(j)+eps).and.(real(part(ipart)%ylat).lt.yrn(j)-eps)) then
        ngr=j
        exit
      endif
    end do
  endif

  if (ngr.ne.ngrid) return

  call advance_PettersonCorrection(itime,ipart)
end subroutine advance

subroutine advance_abovePBL(itime,itimec,dxsave,dysave,&
  ux,vy,tropop,nrand,ipart)

  implicit none
  integer, intent(in) ::          &
    itime,                        & ! time index
    ipart                           ! particle index
  integer, intent(inout) ::       &
    itimec,                       & ! next timestep
    nrand                           ! random number used for turbulence
  real, intent(in) ::             &
    tropop                          ! height of troposphere
  real, intent(inout) ::          &
    ux,vy,                        & ! random turbulent velocities above PBL
    dxsave,dysave                   ! accumulated displacement in long and lat
  real ::                         &
    dt,                           & ! real(ldt)
    xts,yts,zts,ztseta,           & ! local 'real' copy of the particle position
    wp,                           & ! random turbulence velocities
    settling = 0.                   ! settling velocity
  integer ::                      &
    nsp                             ! loop variables for number of species

  zts=real(part(ipart)%z)
  ztseta=real(part(ipart)%zeta)
  xts=real(part(ipart)%xlon)
  yts=real(part(ipart)%ylat)
  call interpol_wind(itime,xts,yts,zts,ztseta,ipart)

  ! Compute everything for above the PBL

  ! Assume constant, uncorrelated, turbulent perturbations
  ! In the stratosphere, use a small vertical diffusivity d_strat,
  ! in the troposphere, use a larger horizontal diffusivity d_trop.
  ! Turbulent velocity scales are determined based on sqrt(d_trop/dt)
  !******************************************************************

  part(ipart)%idt=abs(lsynctime-itimec+itime)
  dt=real(part(ipart)%idt)

  if (.not.turboff) then
    call turbulence_stratosphere(dt,nrand,ux,vy,wp,tropop,zts)
  else
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
        if (xmass(part(ipart)%npoint,nsp).gt.eps3) exit
      end do
      if (nsp.gt.nspec) then
        nsp=nspec
      end if
      ! LB needs to be checked if this works with openmp and change to eta coords
      if (density(nsp).gt.0.) then
        call get_settling(itime,xts,yts,zts,nsp,settling)  !bugfix
        w=w+settling
        if (wind_coord_type.eq.'ETA') then 
          call update_zeta_to_z(itime,ipart)
          call update_z(ipart,settling*dt*real(ldirect))
        endif
      end if
    endif
  end if

  ! Calculate position at time step itime+lsynctime
  !************************************************
  dxsave=dxsave+(u+ux)*dt
  dysave=dysave+(v+vy)*dt
 
  select case (wind_coord_type)
    case ('ETA')
      if ((.not.turboff).or.(lsettling)) then
        call update_zeta_to_z(itime,ipart)
        call update_z(ipart,wp*dt*real(ldirect))
        if (part(ipart)%z.lt.0.) call set_z(ipart,min(h-eps2,-1.*part(ipart)%z))  ! if particle below ground -> reflection
        call update_z_to_zeta(itime,ipart)
      endif
      call update_zeta(ipart,weta*dt*real(ldirect))
      if (part(ipart)%zeta.ge.1.) call set_zeta(ipart,1.-(part(ipart)%zeta-1.))
      if (part(ipart)%zeta.eq.1.) call update_zeta(ipart,-eps_eta)
    case ('METER')
      call update_z(ipart,(w+wp)*dt*real(ldirect))
      if (part(ipart)%z.lt.0.) call set_z(ipart,min(h-eps2,-1.*part(ipart)%z))
    case default
      call update_z(ipart,(w+wp)*dt*real(ldirect))
      if (part(ipart)%z.lt.0.) call set_z(ipart,min(h-eps2,-1.*part(ipart)%z))
  end select
end subroutine advance_abovePBL

subroutine advance_PBL(itime,itimec,&
  dxsave,dysave,dawsave,dcwsave,abovePBL,nrand,ipart)
  use drydepo_mod, only: drydepo_probability

  implicit none

  logical, intent(inout) ::       &
    abovePBL                        ! flag that will be set to 'true' if computation needs to be completed above PBL
  integer, intent(in) ::          &
    itime,                        & ! time index
    ipart                           ! particle index
  real, intent(inout) ::          &
    dxsave,dysave,                & ! accumulated displacement in long and lat
    dawsave,dcwsave                 ! accumulated displacement in wind directions
  integer, intent(inout) ::       &
    itimec,                       & ! next timestep
    nrand                           ! random number used for turbulence
  real ::                         &
    dt,                           & ! real(ldt)
    xts,yts,zts,ztseta,           & ! local 'real' copy of the particle position
    rhoa,                         & ! air density, used in CBL
    rhograd,                      & ! vertical gradient of the air density, used in CBL
    settling = 0.                   ! settling velocity
  integer ::                      &
    loop,                         & ! loop variable for time in the PBL
    nsp                             ! loop variable for species

  eps=nxmax/3.e5

  ! BEGIN TIME LOOP
  !================
  loop=0
  pbl_loop : do
    loop=loop+1
    if (method.eq.1) then
      part(ipart)%idt=min(part(ipart)%idt,abs(lsynctime-itimec+itime))
      itimec=itimec+part(ipart)%idt*ldirect
    else
      part(ipart)%idt=abs(lsynctime)
      itimec=itime+lsynctime
    endif
    dt=real(part(ipart)%idt)
    call update_zeta_to_z(itime,ipart)
    xts=real(part(ipart)%xlon)
    yts=real(part(ipart)%ylat)
    zts=real(part(ipart)%z)
    ztseta=real(part(ipart)%zeta)

    zeta=real(part(ipart)%z)/h
    if (loop.eq.1) then
      if (ngrid.le.0) then
        xts=real(part(ipart)%xlon)
        yts=real(part(ipart)%ylat)
        call interpol_all(itime,xts,yts,zts,ztseta)
      else
        call interpol_all_nests(itime,xtn,ytn,zts,ztseta)
      endif

    else
      ! Determine the level below the current position for u,v,rho
      !***********************************************************
      call find_z_level(zts,ztseta) ! Not sure if zteta levels are necessary here

      ! If one of the levels necessary is not yet available,
      ! calculate it
      !*****************************************************
      call interpol_misslev()
    endif


  ! Vertical interpolation of u,v,w,rho and drhodz
  !***********************************************

  ! Vertical distance to the level below and above current position
  ! both in terms of (u,v) and (w) fields
  !****************************************************************
    call interpol_mixinglayer(zts,ztseta,rhoa,rhograd)

  ! Compute the turbulent disturbances
  ! Determine the sigmas and the timescales 
  !****************************************
    if (.not.turboff) then
      call turbulence_boundarylayer(ipart,nrand,dt,zts,rhoa,rhograd) ! Note: zts and nrand get updated
    ! Determine time step for next integration
    !*****************************************
      if (turbswitch) then
        part(ipart)%idt=int(min(tlw,h/max(2.*abs(part(ipart)%turbvel%w*sigw),1.e-5), &
             0.5/abs(dsigwdz))*ctl)
      else
        part(ipart)%idt=int(min(tlw,h/max(2.*abs(part(ipart)%turbvel%w),1.e-5))*ctl)
      endif
    else
      part(ipart)%turbvel%u=0.0
      part(ipart)%turbvel%v=0.0
      part(ipart)%turbvel%w=0.0
    endif

    part(ipart)%idt=max(part(ipart)%idt,mintime)


  ! If particle represents only a single species, add gravitational settling
  ! velocity. The settling velocity is zero for gases, or if particle
  ! represents more than one species
  !*************************************************************************

    if (mdomainfill.eq.0) then
      if (lsettling) then
        do nsp=1,nspec
          if (xmass(part(ipart)%npoint,nsp).gt.eps3) exit
        end do
        if (nsp.gt.nspec) then
          nsp=nspec
        end if
        if (density(nsp).gt.0.) then
          call get_settling(itime,xts,yts,zts,nsp,settling)  !bugfix
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
    dawsave=dawsave+part(ipart)%turbvel%u*dt
    dcwsave=dcwsave+part(ipart)%turbvel%v*dt
    ! How can I change the w to w(eta) efficiently?
    select case (wind_coord_type)
      case ('ETA')
        call update_z(ipart,w*dt*real(ldirect))
        zts=real(part(ipart)%z)
        ! HSO/AL: Particle managed to go over highest level -> interpolation error in goto 700
        !          alias interpol_wind (division by zero)
        if (zts.ge.height(nz)) call set_z(ipart,height(nz)-100.*eps)
        zts=real(part(ipart)%z)
      case ('METER')
        call update_z(ipart,w*dt*real(ldirect))
        call advance_adjusttopheight(ipart)
    end select

    
    if (zts.gt.h) then
      call update_z_to_zeta(itime,ipart)
      if (itimec.ne.itime+lsynctime) abovePBL=.true. ! complete the current interval above PBL
      return 
    endif
    
  ! Determine probability of deposition
  !************************************
    call drydepo_probability(part(ipart)%prob,dt,zts)

    if (zts.lt.0.) call set_z(ipart,min(h-eps2,-1.*part(ipart)%z))    ! if particle below ground -> reflection


    if (itimec.eq.(itime+lsynctime)) then
      call interpol_average()
      ! Converting the z position that changed through turbulence motions to eta coords
      call update_z_to_zeta(itime,ipart)
      return  ! finished
    endif
  end do pbl_loop
  call update_z_to_zeta(itime,ipart)
end subroutine advance_PBL

subroutine advance_PettersonCorrection(itime,ipart)

  implicit none 

  integer, intent(in) ::          &
    itime,                        & ! time index
    ipart                           ! particle index
  integer ::                      &
    nsp                             ! loop variables for number of species
  real ::                         &
    xts,yts,zts,ztseta,           & ! local 'real' copy of the particle position
    uold,vold,wold,woldeta,       & !
    settling = 0.                   ! settling velocity
  real(kind=dp) ::                &
    ztemp                           ! temporarily storing z position


  xts=real(part(ipart)%xlon)
  yts=real(part(ipart)%ylat)
  zts=real(part(ipart)%z)
  ztseta=real(part(ipart)%zeta)


  ! Determine nested grid coordinates
  !**********************************
  call determine_grid_coordinates(xts,yts)

  ! Memorize the old wind
  !**********************

  uold=u
  vold=v

  select case (wind_coord_type)
    case ('ETA')
      woldeta=weta
    case ('METER')
      wold=w
    case default
      wold=w
  end select

  ! Interpolate wind at new position and time
  !******************************************

  call interpol_wind_short(itime+part(ipart)%idt*ldirect,xts,yts,zts,ztseta)

  if (mdomainfill.eq.0) then
    if (lsettling) then
      do nsp=1,nspec
        if (xmass(part(ipart)%npoint,nsp).gt.eps3) exit
      end do
      if (nsp.gt.nspec) then
        nsp=nspec
      end if
      if (density(nsp).gt.0.) then
        select case (wind_coord_type)

          case ('ETA')
            call update_zeta_to_z(itime+part(ipart)%idt,ipart)
            call get_settling(itime+part(ipart)%idt,xts,yts,real(part(ipart)%z),nsp,settling) !bugfix
            call z_to_zeta(itime+part(ipart)%idt,part(ipart)%xlon,part(ipart)%ylat,&
              part(ipart)%z+real(settling*real(part(ipart)%idt*ldirect),kind=dp),ztemp)
            weta=weta+(real(ztemp)-ztseta)/real(part(ipart)%idt*ldirect)

          case ('METER')
            call get_settling(itime+part(ipart)%idt,xts,yts,zts,nsp,settling)
            w=w+settling

          case default 
            call get_settling(itime+part(ipart)%idt,xts,yts,zts,nsp,settling)
            w=w+settling
        end select            
      end if
    endif
  end if

  ! Determine the difference vector between new and old wind
  ! (use half of it to correct position according to Petterssen)
  !*************************************************************

  u=(u-uold)/2.
  v=(v-vold)/2.

  select case (wind_coord_type)
    case ('ETA')
      weta=(weta-woldeta)/2.
      call update_zeta(ipart,weta*real(part(ipart)%idt*ldirect))
      if (part(ipart)%zeta.ge.1.) call set_zeta(ipart,1.-(part(ipart)%zeta-1.))
      if (part(ipart)%zeta.eq.1.) call update_zeta(ipart,-eps_eta)

    case ('METER')
      w=(w-wold)/2.
      call update_z(ipart,w*real(part(ipart)%idt*ldirect))
      if (part(ipart)%z.lt.0.) call set_z(ipart,min(h-eps2,-1.*part(ipart)%z))    ! if particle below ground -> reflection

    case default 
      w=(w-wold)/2.
      call update_z(ipart,w*real(part(ipart)%idt*ldirect))
      if (part(ipart)%z.lt.0.) call set_z(ipart,min(h-eps2,-1.*part(ipart)%z)) 
  end select  

  ! Finally, correct the old position
  !**********************************
  call advance_updateXY(u*part(ipart)%idt,v*part(ipart)%idt,ipart)

  ! If particle above highest model level, set it back into the domain
  !*******************************************************************
  call advance_adjusttopheight(ipart)
end subroutine advance_PettersonCorrection

subroutine advance_updateXY(xchange,ychange,ipart)

  implicit none

  integer, intent(in) ::          &
    ipart                           ! particle number
  real, intent(in) ::             &
    xchange,ychange                 ! change in position
  real ::                         &
    xlon,ylat,xpol,ypol,          & ! temporarily storing new particle positions
    gridsize,cosfact                ! used to compute new positions of particles

  eps=nxmax/3.e5

  if (ngrid.ge.0) then
    cosfact=dxconst/cos((real(part(ipart)%ylat)*dy+ylat0)*pi180)
    call update_xlon(ipart,real(xchange*cosfact*real(ldirect),kind=dp))
    call update_ylat(ipart,real(ychange*dyconst*real(ldirect),kind=dp))
  else if (ngrid.eq.-1) then      ! around north pole
    xlon=xlon0+real(part(ipart)%xlon)*dx            !comment by mc: compute old particle position
    ylat=ylat0+real(part(ipart)%ylat)*dy
    call cll2xy(northpolemap,ylat,xlon,xpol,ypol)   !convert old particle position in polar stereographic
    gridsize=1000.*cgszll(northpolemap,ylat,xlon)   !calculate size in m of grid element in polar stereographic coordinate
    xpol=xpol+xchange/gridsize*real(ldirect)        !position in grid unit polar stereographic
    ypol=ypol+ychange/gridsize*real(ldirect)
    call cxy2ll(northpolemap,xpol,ypol,ylat,xlon)   !convert to lat long coordinate
    call set_xlon(ipart,real((xlon-xlon0)/dx,kind=dp))!convert to grid units in lat long coordinate, comment by mc
    call set_ylat(ipart,real((ylat-ylat0)/dy,kind=dp))
  else if (ngrid.eq.-2) then    ! around south pole
    xlon=xlon0+real(part(ipart)%xlon)*dx
    ylat=ylat0+real(part(ipart)%ylat)*dy
    call cll2xy(southpolemap,ylat,xlon,xpol,ypol)
    gridsize=1000.*cgszll(southpolemap,ylat,xlon)
    xpol=xpol+xchange/gridsize*real(ldirect)
    ypol=ypol+ychange/gridsize*real(ldirect)
    call cxy2ll(southpolemap,xpol,ypol,ylat,xlon)
    call set_xlon(ipart,real((xlon-xlon0)/dx,kind=dp))
    call set_ylat(ipart,real((ylat-ylat0)/dy,kind=dp))
  endif

  ! If global data are available, use cyclic boundary condition
  !************************************************************
  if (xglobal) then
    if (part(ipart)%xlon.ge.real(nxmin1,kind=dp)) call update_xlon(ipart,-real(nxmin1,kind=dp))
    if (part(ipart)%xlon.lt.0.) call update_xlon(ipart,real(nxmin1,kind=dp))
    if (part(ipart)%xlon.le.real(eps,kind=dp)) call set_xlon(ipart,real(eps,kind=dp))
    if (abs(part(ipart)%xlon-real(nxmin1,kind=dp)).le.eps) call set_xlon(ipart,real(nxmin1-eps,kind=dp))
  endif

  ! HSO/AL: Prevent particles from disappearing at the pole
  !******************************************************************
  if ( part(ipart)%ylat.lt.0. ) then
    call set_xlon(ipart,mod(part(ipart)%xlon+180.,360.))
    call set_ylat(ipart,-part(ipart)%ylat)
  else if ( part(ipart)%ylat.gt.real(nymin1,kind=dp) ) then
    call set_xlon(ipart,mod(part(ipart)%xlon+180.,360.))
    call set_ylat(ipart,2.*real(nymin1,kind=dp)-part(ipart)%ylat)
  endif

  ! Check position: If trajectory outside model domain, terminate it
  !*****************************************************************
  if ((part(ipart)%xlon.lt.0.).or.(part(ipart)%xlon.ge.real(nxmin1,kind=dp)).or.(part(ipart)%ylat.lt.0.).or. &
       (part(ipart)%ylat.gt.real(nymin1,kind=dp))) then
    part(ipart)%nstop=.true.
    return
  endif
end subroutine advance_updateXY

subroutine advance_adjusttopheight(ipart)

  implicit none 

  integer, intent(in) ::          &
    ipart                           ! particle index

  eps=nxmax/3.e5

  select case (wind_coord_type)
    case ('ETA')
      if (part(ipart)%zeta.le.real(uvheight(nz),kind=dp)) then
        call set_zeta(ipart,uvheight(nz)+eps_eta)
      endif
    case ('METER')
      if (part(ipart)%z.ge.real(height(nz),kind=dp)) call set_z(ipart,height(nz)-100.*eps)
    case default
      if (part(ipart)%z.ge.real(height(nz),kind=dp)) call set_z(ipart,height(nz)-100.*eps)
  end select
end subroutine advance_adjusttopheight

end module advance_mod
