! SPDX-FileCopyrightText: FLEXPART 1998-2019, see flexpart_license.txt
! SPDX-License-Identifier: GPL-3.0-or-later

subroutine timemanager(metdata_format)

  !*****************************************************************************
  !                                                                            *
  ! Handles the computation of trajectories, i.e. determines which             *
  ! trajectories have to be computed at what time.                             *
  ! Manages dry+wet deposition routines, radioactive decay and the computation *
  ! of concentrations.                                                         *
  !                                                                            *
  !     Author: A. Stohl                                                       *
  !                                                                            *
  !     20 May 1996                                                            *
  !                                                                            *
  !*****************************************************************************
  !  Changes, Bernd C. Krueger, Feb. 2001:                                     *
  !        Call of convmix when new windfield is read                          *
  !------------------------------------                                        *
  !  Changes Petra Seibert, Sept 2002                                          *
  !     fix wet scavenging problem                                             *
  !     Code may not be correct for decay of deposition!                       *
  !  Changes Petra Seibert, Nov 2002                                           *
  !     call convection BEFORE new fields are read in BWD mode                 *
  !  Changes Caroline Forster, Feb 2005                                        *
  !   new interface between flexpart and convection scheme                     *
  !   Emanuel's latest subroutine convect43c.f is used                         *
  !  Changes Stefan Henne, Harald Sodemann, 2013-2014                          *
  !   added netcdf output code                                                 *
  !  Changes Espen Sollum 2014                                                 *
  !   For compatibility with MPI version,                                      *
  !   variables uap,ucp,uzp,us,vs,ws,cbt now in module com_mod                 *
  !  Unified ECMWF and GFS builds                                              *
  !   Marian Harustak, 12.5.2017                                               *
  !   - Added passing of metdata_format as it was needed by called routines    *
  !*****************************************************************************
  !                                                                            *
  ! Variables:                                                                 *
  ! DEP                .true. if either wet or dry deposition is switched on   *
  ! decay(maxspec) [1/s] decay constant for radioactive decay                  *
  ! DRYDEP             .true. if dry deposition is switched on                 *
  ! ideltas [s]        modelling period                                        *
  ! itime [s]          actual temporal position of calculation                 *
  ! ldeltat [s]        time since computation of radioact. decay of depositions*
  ! loutaver [s]       averaging period for concentration calculations         *
  ! loutend [s]        end of averaging for concentration calculations         *
  ! loutnext [s]       next time at which output fields shall be centered      *
  ! loutsample [s]     sampling interval for averaging of concentrations       *
  ! loutstart [s]      start of averaging for concentration calculations       *
  ! loutstep [s]       time interval for which concentrations shall be         *
  !                    calculated                                              *
  ! npoint             index, which starting point the trajectory has          *
  !                    starting positions of trajectories                      *
  ! nstop              serves as indicator for fate of particles               *
  !                    in the particle loop                                    *
  ! nstop1             serves as indicator for wind fields (see getfields)     *
  ! outnum             number of samples for each concentration calculation    *
  ! outnum             number of samples for each concentration calculation    *
  ! prob               probability of absorption at ground due to dry          *
  !                    deposition                                              *
  ! WETDEP             .true. if wet deposition is switched on                 *
  ! weight             weight for each concentration sample (1/2 or 1)         *
  ! metdata_format     format of metdata (ecmwf/gfs)                           *
  !                                                                            *
  !*****************************************************************************
  ! openmp change
  use omp_lib
  ! openmp change end
  use unc_mod
  use point_mod
  use xmass_mod
  use flux_mod
  use outg_mod
  use oh_mod
  use par_mod
  use com_mod
#ifdef USE_NCF
  use netcdf_output_mod, only: concoutput_netcdf,concoutput_nest_netcdf,&
       &concoutput_surf_netcdf,concoutput_surf_nest_netcdf,writeheader_partoutput
#endif
  use coordinates_ecmwf
  use particle_mod
  use conv_mod
  use windfields_mod
  use advance_mod, only: advance

  implicit none
  real, parameter ::        &
    e_inv = 1.0/exp(1.0)  
  integer, intent(in) ::    &
    metdata_format            ! Data type of the windfields
  integer ::                &
    j,i,                    & ! loop variable
    ks,                     & ! loop variable species
    kp,                     & ! loop variable for maxpointspec_act
    l,                      & ! loop variable over nclassunc
    n,                      & ! loop variable over particles
    itime=0,                & ! time index
    nstop1,                 & ! windfield existence flag
    loutnext,               & ! following timestep
    loutstart,loutend,      & ! concentration calculation starting and ending time
    ix,jy,                  & ! gridcell indices
    ldeltat,                & ! radioactive decay time
    itage,nage,             & ! related to age classes
    idummy,                 & ! used for the random routines
    i_nan=0,ii_nan,total_nan_intl=0, &  !added by mc to check instability in CBL scheme 
    thread                    ! openmp change (not sure if necessary)
  ! logical ::                &
  !   active_per_rel(maxpoint)  ! are there particles active in each release
#ifdef USE_NCF
  real ::                   &
    filesize!(maxpoint)        ! Keeping track of the size of the particledump output, so it can be splitted
  real(kind=dp) ::          &
    jul
  integer ::                &
    jjjjmmdd,ihmmss
#endif
  real ::                   &
    outnum,                 & ! concentration calculation sample number
    weight,                 & ! concentration calculation sample weight
    prob_rec(maxspec),      & ! dry deposition related
    decfact,                & ! radioactive decay factor
    wetscav,                & ! wet scavenging
    xmassfract,             & ! dry deposition related
    grfraction(3)             ! wet deposition related
  real(sp) ::               &
    gridtotalunc              ! concentration calculation related
  real(dep_prec) ::         &
    drydeposit(maxspec),    & ! dry deposition related
    wetgridtotalunc,        & ! concentration calculation related
    drygridtotalunc           ! concentration calculation related

  integer :: npart_alive=0

  ! First output for time 0
  !************************

  loutnext=loutstep/2
  outnum=0.
  loutstart=loutnext-loutaver/2
  loutend=loutnext+loutaver/2

  !**********************************************************************
  ! Loop over the whole modelling period in time steps of mintime seconds
  !**********************************************************************

  write(*,46) float(itime)/3600,itime,numpart

#ifdef USE_NCF
  filesize=0.
  ! active_per_rel=.false.
#endif

  do itime=0,ideltas,lsynctime

  ! Computation of wet deposition, OH reaction and mass transfer
  ! between two species every lsynctime seconds
  ! maybe wet depo frequency can be relaxed later but better be on safe side
  ! wetdepo must be called BEFORE new fields are read in but should not
  ! be called in the very beginning before any fields are loaded, or
  ! before particles are in the system
  ! Code may not be correct for decay of deposition
  ! changed by Petra Seibert 9/02
  !********************************************************************

    if (WETDEP .and. itime .ne. 0 .and. numpart .gt. 0) then
      call wetdepo(itime,lsynctime,loutnext)
    endif

    if (OHREA .and. itime .ne. 0 .and. numpart .gt. 0) &
      call ohreaction(itime,lsynctime,loutnext)

  ! compute convection for backward runs
  !*************************************

    if ((ldirect.eq.-1).and.(lconvection.eq.1).and.(itime.lt.0)) then    
      call convmix(itime,metdata_format)
    endif

  ! Get necessary wind fields if not available
  !*******************************************
    call getfields(itime,nstop1,metdata_format)
    if (nstop1.gt.1) stop 'NO METEO FIELDS AVAILABLE'

  ! In case of ETA coordinates being read from file, convert the z positions
  !*************************************************************************
    if ((ipin.eq.1).and.(itime.eq.0).and.(wind_coord_type.eq.'ETA')) then 
      if (numpart.le.0) stop 'Something is going wrong reading the old particle file!'
!$OMP PARALLEL PRIVATE(i)
!$OMP DO
      do i=1,numpart
        call update_z_to_zeta(itime, i)
      end do
!$OMP END DO
!$OMP END PARALLEL
    endif

  ! Get hourly OH fields if not available 
  !****************************************************
    if (OHREA) then
      call gethourlyOH(itime)
    endif
        
  ! Release particles
  !******************
    if (mdomainfill.ge.1) then
      if (itime.eq.0) then   
        call init_domainfill
      else 
        call boundcond_domainfill(itime,loutend)
      endif
    else
      call releaseparticles(itime)
    endif

#ifdef USE_NCF
    if (ipout.ge.1) then
      if (itime.eq.0) then
        if (ldirect.eq.1) then
          call writeheader_partoutput(ibtime,ibdate,ibtime,ibdate)
        else 
          call writeheader_partoutput(ietime,iedate,ietime,iedate)
        endif
      else if (mod(itime,ipoutfac*loutstep).eq.0) then
        if (filesize.ge.max_partoutput_filesize) then 
          jul=bdate+real(itime,kind=dp)/86400._dp
          call caldate(jul,jjjjmmdd,ihmmss)
          if (ldirect.eq.1) then 
            call writeheader_partoutput(ihmmss,jjjjmmdd,ibtime,ibdate)
          else 
            call writeheader_partoutput(ihmmss,jjjjmmdd,ietime,iedate)
          endif
          filesize = 0.
        endif
        do i=1,numpoint
          filesize = filesize + npart(i)*13.*4./1000000.
        end do
      endif
    endif
#endif

  ! Compute convective mixing for forward runs
  ! for backward runs it is done before next windfield is read in
  !**************************************************************
    if ((ldirect.eq.1).and.(lconvection.eq.1)) then
      call convmix(itime,metdata_format)
    endif

  ! If middle of averaging period of output fields is reached, accumulated
  ! deposited mass radioactively decays
  !***********************************************************************
  ! This should go in a subroutine
    if (DEP.and.(itime.eq.loutnext).and.(ldirect.gt.0)) then
      do ks=1,nspec
      do kp=1,maxpointspec_act
        if (decay(ks).gt.0.) then
          do nage=1,nageclass
            do l=1,nclassunc
  ! Mother output grid
              do jy=0,numygrid-1
                do ix=0,numxgrid-1
                  wetgridunc(ix,jy,ks,kp,l,nage)= &
                       wetgridunc(ix,jy,ks,kp,l,nage)* &
                       exp(-1.*outstep*decay(ks))
                  drygridunc(ix,jy,ks,kp,l,nage)= &
                       drygridunc(ix,jy,ks,kp,l,nage)* &
                       exp(-1.*outstep*decay(ks))
                end do
              end do
  ! Nested output grid
              if (nested_output.eq.1) then
                do jy=0,numygridn-1
                  do ix=0,numxgridn-1
                    wetgriduncn(ix,jy,ks,kp,l,nage)= &
                         wetgriduncn(ix,jy,ks,kp,l,nage)* &
                         exp(-1.*outstep*decay(ks))
                    drygriduncn(ix,jy,ks,kp,l,nage)= &
                         drygriduncn(ix,jy,ks,kp,l,nage)* &
                         exp(-1.*outstep*decay(ks))
                  end do
                end do
              endif
            end do
          end do
        endif
      end do
      end do
    endif

  ! Check whether concentrations are to be calculated
  !**************************************************
  ! Put all of the concentration stuff in a subroutine
    if ((ldirect*itime.ge.ldirect*loutstart).and. &
         (ldirect*itime.le.ldirect*loutend)) then ! add to grid
      if (mod(itime-loutstart,loutsample).eq.0) then

  ! If we are exactly at the start or end of the concentration averaging interval,
  ! give only half the weight to this sample
  !*****************************************************************************

        if ((itime.eq.loutstart).or.(itime.eq.loutend)) then
          weight=0.5
        else
          weight=1.0
        endif
        outnum=outnum+weight
        call conccalc(itime,weight)
      endif

  ! Output and reinitialization of grid
  ! If necessary, first sample of new grid is also taken
  !*****************************************************

      if ((itime.eq.loutend).and.(outnum.gt.0.)) then
if (grid_output.eq.1) then
        if ((iout.le.3.).or.(iout.eq.5)) then
          if (surf_only.ne.1) then 
            if (lnetcdfout.eq.1) then 
#ifdef USE_NCF
              call concoutput_netcdf(itime,outnum,gridtotalunc,wetgridtotalunc,drygridtotalunc)
#endif
            else 
              call concoutput(itime,outnum,gridtotalunc,wetgridtotalunc,drygridtotalunc)
            endif
          else  
            if (verbosity.eq.1) then
              print*,'call concoutput_surf '
              call system_clock(count_clock)
              write(*,*) 'system clock',count_clock - count_clock0   
            endif
            if (lnetcdfout.eq.1) then
#ifdef USE_NCF
              call concoutput_surf_netcdf(itime,outnum,gridtotalunc,wetgridtotalunc,drygridtotalunc)
#endif
            else
              if (linversionout.eq.1) then
                call concoutput_inversion(itime,outnum,gridtotalunc,wetgridtotalunc,drygridtotalunc)
                if (verbosity.eq.1) then
                  print*,'called concoutput_inversion'
                  call system_clock(count_clock)
                  write(*,*) 'system clock',count_clock - count_clock0 
                endif
              else
                call concoutput_surf(itime,outnum,gridtotalunc,wetgridtotalunc,drygridtotalunc)
              endif
              if (verbosity.eq.1) then
                print*,'called concoutput_surf '
                call system_clock(count_clock)
                write(*,*) 'system clock',count_clock - count_clock0   
              endif
            endif
          endif

          if (nested_output .eq. 1) then
            if (lnetcdfout.eq.0) then
              if (surf_only.ne.1) then
                call concoutput_nest(itime,outnum)
              else 
                if(linversionout.eq.1) then
                  call concoutput_inversion_nest(itime,outnum)
                else 
                call concoutput_surf_nest(itime,outnum)
              endif
              endif
            else
#ifdef USE_NCF
              if (surf_only.ne.1) then
                call concoutput_nest_netcdf(itime,outnum)
              else 
                call concoutput_surf_nest_netcdf(itime,outnum)
              endif
#endif
            endif
          endif
          outnum=0.
        endif
endif
        if ((iout.eq.4).or.(iout.eq.5)) call plumetraj(itime)
        if (iflux.eq.1) call fluxoutput(itime)
        write(*,45) itime,numpart,gridtotalunc,wetgridtotalunc,drygridtotalunc
 

45      format(i13,' Seconds simulated: ',i13, ' Particles:    Uncertainty: ',3f7.3)
46      format(' Simulated ',f7.1,' hours (',i13,' s), ',i13, ' particles')
        if (ipout.ge.1) then
          if (mod(itime,ipoutfac*loutstep).eq.0) then
            call SYSTEM_CLOCK(count_clock, count_rate, count_max)
            s_temp = (count_clock - count_clock0)/real(count_rate)
            call partoutput(itime)!,active_per_rel) ! dump particle positions
            call SYSTEM_CLOCK(count_clock, count_rate, count_max)
            s_writepart = s_writepart + ((count_clock - count_clock0)/real(count_rate)-s_temp)
          endif
        endif
        loutnext=loutnext+loutstep
        loutstart=loutnext-loutaver/2
        loutend=loutnext+loutaver/2
        if (itime.eq.loutstart) then
          weight=0.5
          outnum=outnum+weight
          call conccalc(itime,weight)
        endif

      endif
    endif


    if (itime.eq.ideltas) exit         ! almost finished

  ! Compute interval since radioactive decay of deposited mass was computed
  !************************************************************************

    if (itime.lt.loutnext) then
      ldeltat=itime-(loutnext-loutstep)
    else                                  ! first half of next interval
      ldeltat=itime-loutnext
    endif


  ! Loop over all particles
  !************************
  ! Various variables for testing reason of CBL scheme, by mc
    well_mixed_vector=0. !erase vector to test well mixed condition: modified by mc
    well_mixed_norm=0.   !erase normalization to test well mixed condition: modified by mc
    avg_ol=0.
    avg_wst=0.
    avg_h=0.
    avg_air_dens=0.  !erase vector to obtain air density at particle positions: modified by mc
  !-----------------------------------------------------------------------------

  ! openmp change
  ! LB, openmp following CTM version, need to be very careful due to big differences
  ! between the openmp loop in this and the CTM version
!$OMP PARALLEL PRIVATE(prob_rec,ks,thread,j)

#if (defined _OPENMP)
        thread = OMP_GET_THREAD_NUM()
#endif

!$OMP DO
    do j=1,numpart

  ! If integration step is due, do it
  !**********************************
      if (.not. part(j)%alive) cycle

  ! Initialize newly released particle
  !***********************************
      if ((part(j)%tstart.eq.itime).or.(itime.eq.0)) then
        call update_zeta_to_z(itime, j)
        call initialize(itime,part(j)%idt, &
            part(j)%turbvel%u,part(j)%turbvel%v,part(j)%turbvel%w, &
            part(j)%mesovel%u,part(j)%mesovel%v,part(j)%mesovel%w, &
            part(j)%xlon,part(j)%ylat,part(j)%z, &
            part(j)%zeta,part(j)%icbt)
      endif

  ! Memorize particle positions
  !****************************
      part(j)%xlon_prev=part(j)%xlon
      part(j)%ylat_prev=part(j)%ylat
      part(j)%z_prev=part(j)%z

  ! RECEPTOR: dry/wet depovel
  !****************************
  ! Before the particle is moved 
  ! the calculation of the scavenged mass shall only be done once after release
  ! xscav_frac1 was initialised with a negative value

      if  (DRYBKDEP) then
        do ks=1,nspec
          if  ((xscav_frac1(j,ks).lt.0)) then
            call update_zeta_to_z(itime,j)
            call get_vdep_prob(itime,part(j)%xlon,part(j)%ylat,part(j)%z,prob_rec)
            if (DRYDEPSPEC(ks)) then        ! dry deposition
              xscav_frac1(j,ks)=prob_rec(ks)
            else
              part(j)%mass(ks)=0.
              xscav_frac1(j,ks)=0.
            endif
          endif
        enddo
      endif

  ! Integrate Lagevin equation for lsynctime seconds
  !*************************************************

      call advance(itime,j)

    end do 

!$OMP END DO
!$OMP END PARALLEL

    do j=1,numpart
  ! If integration step is due, do it
  !**********************************
      if (.not. part(j)%alive) cycle

  ! Determine age class of the particle
      itage=abs(itime-part(j)%tstart)
      do nage=1,nageclass
        if (itage.lt.lage(nage)) exit
      end do
  ! Calculate the gross fluxes across layer interfaces
  !***************************************************

      if (iflux.eq.1) call calcfluxes(itime,nage,j,part(j)%xlon_prev, &
        part(j)%ylat_prev,part(j)%z_prev) !OMP reduction necessary for flux array


  ! Determine, when next time step is due
  ! If trajectory is terminated, mark it
  !**************************************
      if (part(j)%nstop) then
        if (linit_cond.ge.1) call initial_cond_calc(itime,j) !OMP reduction necessary for init_cond
        call terminate_particle(j)
      else

  ! Dry deposition and radioactive decay for each species
  ! Also check maximum (of all species) of initial mass remaining on the particle;
  ! if it is below a threshold value, terminate particle
  !*****************************************************************************

        xmassfract=0.
        do ks=1,nspec
          if (decay(ks).gt.0.) then             ! radioactive decay
            decfact=exp(-real(abs(lsynctime))*decay(ks))
          else
            decfact=1.
          endif

          if (DRYDEPSPEC(ks)) then        ! dry deposition
            drydeposit(ks)=part(j)%mass(ks)*part(j)%prob(ks)*decfact
            part(j)%mass(ks)=part(j)%mass(ks)*(1.-part(j)%prob(ks))*decfact
            if (decay(ks).gt.0.) then   ! correct for decay (see wetdepo)
              drydeposit(ks)=drydeposit(ks)* &
                   exp(real(abs(ldeltat))*decay(ks))
            endif
          else                           ! no dry deposition
            part(j)%mass(ks)=part(j)%mass(ks)*decfact
          endif

  ! Skip check on mass fraction when npoint represents particle number
          if (mdomainfill.eq.0.and.mquasilag.eq.0) then
            if (xmass(part(j)%npoint,ks).gt.0.) then
                 xmassfract=max(xmassfract,real(npart(part(j)%npoint))* &
                 part(j)%mass(ks)/xmass(part(j)%npoint,ks))
            endif

          else
            xmassfract=1.0
          end if
        end do

        if (xmassfract.lt.minmass) then   ! terminate all particles carrying less mass
          call terminate_particle(j)
        endif

!        Sabine Eckhardt, June 2008
!        don't create depofield for backward runs
        if (DRYDEP.AND.(ldirect.eq.1)) then !OMP reduction necessary for drygridunc

          if (ioutputforeachrelease.eq.1) then
              kp=part(j)%npoint
          else
              kp=1
          endif

          call drydepokernel(part(j)%nclass,drydeposit,real(part(j)%xlon), &
               real(part(j)%ylat),nage,kp)
          if (nested_output.eq.1) call drydepokernel_nest( &
               part(j)%nclass,drydeposit,real(part(j)%xlon),real(part(j)%ylat), &
               nage,kp)
        endif

  ! Terminate trajectories that are older than maximum allowed age
  !***************************************************************

        if ((part(j)%alive).and.(abs(itime-part(j)%tstart).ge.lage(nageclass))) then
          if (linit_cond.ge.1) call initial_cond_calc(itime+lsynctime,j)
          call terminate_particle(j)
        endif
      endif

    end do !loop over particles

  ! openmp change end

  ! Counter of "unstable" particle velocity during a time scale of
  ! maximumtl=20 minutes (defined in com_mod)
  !***************************************************************
    
    total_nan_intl=0
    i_nan=i_nan+1 ! added by mc to count nan during a time of maxtl (i.e. maximum tl fixed here to 20 minutes, see com_mod)
    sum_nan_count(i_nan)=nan_count
    if (i_nan > maxtl/lsynctime) i_nan=1 !lsynctime must be <= maxtl
    do ii_nan=1, (maxtl/lsynctime) 
      total_nan_intl=total_nan_intl+sum_nan_count(ii_nan)
    end do
  ! Output to keep track of the numerical instabilities in CBL simulation and if
  ! they are compromising the final result (or not)
    if (cblflag.eq.1) print *,j,itime,'nan_synctime',nan_count,'nan_tl',total_nan_intl  
          
  end do


  ! Complete the calculation of initial conditions for particles not yet terminated
  !*****************************************************************************

  do j=1,numpart
    if (linit_cond.ge.1) call initial_cond_calc(itime,j)
  end do

  if (ipout.eq.2) call partoutput(itime)!,active_per_rel)     ! dump particle positions

  if (linit_cond.ge.1) then
    if(linversionout.eq.1) then
      call initial_cond_output_inversion(itime)   ! dump initial cond. field
    else
      call initial_cond_output(itime)   ! dump initial cond. fielf
    endif
  endif

  ! De-allocate memory and end
  !***************************
  call deallocate_all_particles()
  if (iflux.eq.1) then
      deallocate(flux)
  endif
  if (OHREA) then
      deallocate(OH_field,OH_hourly,lonOH,latOH,altOH)
  endif
  if (grid_output.eq.1) then
  if (ldirect.gt.0) then
  deallocate(drygridunc,wetgridunc)
  endif
  deallocate(gridunc)
  endif
  deallocate(xpoint1,xpoint2,ypoint1,ypoint2,zpoint1,zpoint2,xmass)
  deallocate(ireleasestart,ireleaseend,npart,kindz)
  deallocate(xmasssave)
  if (nested_output.eq.1) then
     deallocate(orooutn, arean, volumen)
     if (ldirect.gt.0) then
     deallocate(griduncn,drygriduncn,wetgriduncn)
     endif
  endif
  if (grid_output.eq.1) then
  deallocate(outheight,outheighthalf)
  deallocate(oroout, area, volume)
  endif

end subroutine timemanager

