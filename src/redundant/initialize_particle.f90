! SPDX-FileCopyrightText: FLEXPART 1998-2019, see flexpart_license.txt
! SPDX-License-Identifier: GPL-3.0-or-later

subroutine initialize_particle(itime,ipart)
  !                        i    i   o  o  o
  !        o       o       o    i  i  i   o
  !*****************************************************************************
  !                                                                            *
  !  Calculation of trajectories utilizing a zero-acceleration scheme. The time*
  !  step is determined by the Courant-Friedrichs-Lewy (CFL) criterion. This   *
  !  means that the time step must be so small that the displacement within    *
  !  this time step is smaller than 1 grid distance. Additionally, a temporal  *
  !  CFL criterion is introduced: the time step must be smaller than the time  *
  !  interval of the wind fields used for interpolation.                       *
  !  For random walk simulations, these are the only time step criteria.       *
  !  For the other options, the time step is also limited by the Lagrangian    *
  !  time scale.                                                               *
  !                                                                            *
  !     Author: A. Stohl                                                       *
  !                                                                            *
  !     16 December 1997                                                       *
  !                                                                            *
  !  Literature:                                                               *
  !                                                                            *
  !*****************************************************************************
  !                                                                            *
  ! Variables:                                                                 *
  ! h [m]              Mixing height                                           *
  ! lwindinterv [s]    time interval between two wind fields                   *
  ! itime [s]          current temporal position                               *
  ! ldt [s]            Suggested time step for next integration                *
  ! ladvance [s]       Total integration time period                           *
  ! rannumb(maxrand)   normally distributed random variables                   *
  ! usig,vsig,wsig     uncertainties of wind velocities due to interpolation   *
  ! xt,yt,zt           Next time step's spatial position of trajectory         *
  !                                                                            *
  !                                                                            *
  ! Constants:                                                                 *
  ! cfl                factor, by which the time step has to be smaller than   *
  !                    the spatial CFL-criterion                               *
  ! cflt               factor, by which the time step has to be smaller than   *
  !                    the temporal CFL-criterion                              *
  !                                                                            *
  !*****************************************************************************

  use par_mod
  use com_mod
  use windfields_mod
  use interpol_mod
  use turbulence_mod
  use random_mod, only: ran3
  use interpol_mod
  use coordinates_ecmwf
  use particle_mod

  use omp_lib

  implicit none

  integer,intent(in) ::  &
    itime,               &
    ipart
  integer :: i,j,k,m,indexh
  integer :: nrand
  real :: dz,dz1,dz2,wp
  real :: ttemp,dummy1,dummy2
  real :: xt,yt,zt,zteta
  integer :: thread
  save idummy

  integer :: idummy = -7

!$OMP THREADPRIVATE(idummy)
!$    if (idummy.eq.-7) then
!$      thread = OMP_GET_THREAD_NUM()
!$      idummy = idummy - thread
!$    endif 

  part(ipart)%icbt=1           ! initialize particle to no "reflection"

  nrand=int(ran3(idummy)*real(maxrand-1))+1

  xt = real(part(ipart)%xlon)
  yt = real(part(ipart)%ylat)
  zt = real(part(ipart)%z)
  zteta = real(part(ipart)%zeta)

  !******************************
  ! 2. Interpolate necessary data
  !******************************

  ! Compute maximum mixing height around particle position
  !*******************************************************
  call determine_grid_coordinates(xt,yt)
  
  h=max(hmix(ix ,jy,1,memind(1)), &
       hmix(ixp,jy ,1,memind(1)), &
       hmix(ix ,jyp,1,memind(1)), &
       hmix(ixp,jyp,1,memind(1)), &
       hmix(ix ,jy ,1,memind(2)), &
       hmix(ixp,jy ,1,memind(2)), &
       hmix(ix ,jyp,1,memind(2)), &
       hmix(ixp,jyp,1,memind(2)))

  zeta=zt/h


  !*************************************************************
  ! If particle is in the PBL, interpolate once and then make a
  ! time loop until end of interval is reached
  !*************************************************************

  if (zeta.le.1.) then

    call interpol_all(itime,xt,yt,zt,zteta)

  ! Vertical interpolation of u,v,w,rho and drhodz
  !***********************************************

  ! Vertical distance to the level below and above current position
  ! both in terms of (u,v) and (w) fields
  !****************************************************************
    call interpol_mixinglayer(zt,zteta,dummy1,dummy2)

  ! Compute the turbulent disturbances

  ! Determine the sigmas and the timescales
  !****************************************

    if (turbswitch) then
      call hanna(zt)
    else
      call hanna1(zt)
    endif


  ! Determine the new diffusivity velocities
  !*****************************************

    if (nrand+2.gt.maxrand) nrand=1
    part(ipart)%turbvel%u=rannumb(nrand)*sigu
    part(ipart)%turbvel%v=rannumb(nrand+1)*sigv
    part(ipart)%turbvel%w=rannumb(nrand+2)
    if (.not.turbswitch) then     ! modified by mc
      part(ipart)%turbvel%w=part(ipart)%turbvel%w*sigw
    else if (cblflag.eq.1) then   ! modified by mc
      if(-h/ol.gt.5) then
  !if (ol.lt.0.) then
  !if (ol.gt.0.) then !by mc : only for test correct is lt.0
        call initialize_cbl_vel(idummy,zt,ust,wst,h,sigw,part(ipart)%turbvel%w,ol)
      else
        part(ipart)%turbvel%w=part(ipart)%turbvel%w*sigw
      end if
    end if


  ! Determine time step for next integration
  !*****************************************

    if (turbswitch) then
      part(ipart)%idt=int(min(tlw,h/max(2.*abs(part(ipart)%turbvel%w*sigw),1.e-5), &
           0.5/abs(dsigwdz),600.)*ctl)
    else
      part(ipart)%idt=int(min(tlw,h/max(2.*abs(part(ipart)%turbvel%w),1.e-5),600.)*ctl)
    endif
    part(ipart)%idt=max(part(ipart)%idt,mintime)

    call interpol_average()
    ! usig=(usigprof(indzp)+usigprof(indz))/2.
    ! vsig=(vsigprof(indzp)+vsigprof(indz))/2.
    ! wsig=(wsigprof(indzp)+wsigprof(indz))/2.

    ! wsigeta=(wsigprofeta(indzpeta)+wsigprofeta(indzeta))/2.

  else



  !**********************************************************
  ! For all particles that are outside the PBL, make a single
  ! time step. Only horizontal turbulent disturbances are
  ! calculated. Vertical disturbances are reset.
  !**********************************************************


  ! Interpolate the wind
  !*********************

    call interpol_wind(itime,xt,yt,zt,zteta,10)


  ! Compute everything for above the PBL

  ! Assume constant turbulent perturbations
  !****************************************

    part(ipart)%idt=abs(lsynctime)

    if (nrand+1.gt.maxrand) nrand=1
    part(ipart)%turbvel%u=rannumb(nrand)*0.3
    part(ipart)%turbvel%v=rannumb(nrand+1)*0.3
    nrand=nrand+2
    part(ipart)%turbvel%w=0.
    sigw=0.

  endif

  !****************************************************************
  ! Add mesoscale random disturbances
  ! This is done only once for the whole lsynctime interval to save
  ! computation time
  !****************************************************************


  ! It is assumed that the average interpolation error is 1/2 sigma
  ! of the surrounding points, autocorrelation time constant is
  ! 1/2 of time interval between wind fields
  !****************************************************************

  if (nrand+2.gt.maxrand) nrand=1
  part(ipart)%mesovel%u=rannumb(nrand)*usig
  part(ipart)%mesovel%v=rannumb(nrand+1)*vsig
  select case (wind_coord_type)
    case ('ETA')
      part(ipart)%mesovel%w=rannumb(nrand+2)*wsigeta
    case ('METER')
      part(ipart)%mesovel%w=rannumb(nrand+2)*wsig
    case default
      part(ipart)%mesovel%w=rannumb(nrand+2)*wsig
  end select  

end subroutine initialize_particle
