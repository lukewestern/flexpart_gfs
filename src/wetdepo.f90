! SPDX-FileCopyrightText: FLEXPART 1998-2019, see flexpart_license.txt
! SPDX-License-Identifier: GPL-3.0-or-later

subroutine wetdepo(itime,ltsample,loutnext)
!                  i      i        i
!*****************************************************************************
!                                                                            *
! Calculation of wet deposition using the concept of scavenging coefficients.*
! For lack of detailed information, washout and rainout are jointly treated. *
! It is assumed that precipitation does not occur uniformly within the whole *
! grid cell, but that only a fraction of the grid cell experiences rainfall. *
! This fraction is parameterized from total cloud cover and rates of large   *
! scale and convective precipitation.                                        *
!                                                                            *
!    Author: A. Stohl                                                        *
!                                                                            *
!    1 December 1996                                                         *
!                                                                            *
! Correction by Petra Seibert, Sept 2002:                                    *
! use centred precipitation data for integration                             *
! Code may not be correct for decay of deposition!                           *
!                                                                            *
!*****************************************************************************
!                                                                            *
! Variables:                                                                 *
! ix,jy              indices of output grid cell for each particle           *
! itime [s]          actual simulation time [s]                              *
! jpart              particle index                                          *
! ldeltat [s]        interval since radioactive decay was computed           *
! loutnext [s]       time for which gridded deposition is next output        *
! loutstep [s]       interval at which gridded deposition is output          *
! ltsample [s]       interval over which mass is deposited                   *
! wetdeposit         mass that is wet deposited                              *
! wetgrid            accumulated deposited mass on output grid               *
! wetscav            scavenging coefficient                                  *
!                                                                            *
! Constants:                                                                 *
!                                                                            *
!*****************************************************************************

  use point_mod
  use par_mod
  use com_mod
  use unc_mod, only:wetgridunc
  use particle_mod

  implicit none

  integer :: jpart,itime,ltsample,loutnext,ldeltat
  integer :: itage,nage
  integer :: ks, kp
  integer(selected_int_kind(16)), dimension(nspec) :: blc_count, inc_count
  real :: grfraction(3),wetscav
  real :: wetdeposit(maxspec),restmass
  real,parameter :: smallnum = tiny(0.0) ! smallest number that can be handled

! Compute interval since radioactive decay of deposited mass was computed
!************************************************************************

  if (itime.le.loutnext) then
    ldeltat=itime-(loutnext-loutstep)
  else                                  ! first half of next interval
    ldeltat=itime-loutnext
  endif

! Loop over all particles
!************************

  blc_count(:)=0
  inc_count(:)=0

! OMP doesn't work yet, a reduction is necessary for the kernel function
  do jpart=1,numpart

    ! Check if memory has been deallocated
    if (.not. is_particle_allocated(jpart)) cycle

    ! Check if particle is still allive
    if (.not. part(jpart)%alive) cycle

! Determine age class of the particle - nage is used for the kernel
!******************************************************************
     itage=abs(itime-part(jpart)%tstart)
     do nage=1,nageclass
       if (itage.lt.lage(nage)) exit
     end do

    do ks=1,nspec      ! loop over species

      if (WETDEPSPEC(ks).eqv..false.) cycle 

!**************************************************
! CALCULATE DEPOSITION 
!**************************************************
!       wetscav=0.
       
!        write(*,*) ks,dquer(ks), crain_aero(ks),csnow_aero(ks)
!       if (((dquer(ks).le.0.).and.(weta_gas(ks).gt.0..or.wetb_gas(ks).gt.0.)) &
!          .or. &
!          ((dquer(ks).gt.0.).and.(crain_aero(ks).gt.0..or.csnow_aero(ks).gt.0.).or. &
!            (ccn_aero(ks).gt0) .or. (in_aero(ks).gt.0) .or. (henry(ks).gt.0)))  then

      call get_wetscav(itime,ltsample,loutnext,jpart,ks,grfraction,inc_count,blc_count,wetscav) ! OMP carefully check

      ! WETBKDEP moved here from timemanager.f90
      ! xscav_frac1...factor used to weight the wet depostion along the back-trajectory
      ! (based on wetdepo at receptor)
      ! jpart...particle index
      ! ks...species index
      ! grfraction(1)...fraction of grid, for which precipitation occurs; set in get_wetscav
      ! zpoint2(npoint(jpart))...forced to be 20000m for wetdepo backward
      ! zpoint1(npoint(jpart))...forced to be 0m for wetdepo backward
      ! zpoint1,zpoint2...height range, over which release takes place; set in RELEASE file
      ! as z1 and z2 => but for wetdepo backward always set to 0 and 20000m
      if (WETBKDEP) then
        if ((xscav_frac1(jpart,ks).lt.-0.1)) then   ! particle marked as starting particle
          if (wetscav.gt.0.) then
             xscav_frac1(jpart,ks)=wetscav*(zpoint2(part(jpart)%npoint)-&
             zpoint1(part(jpart)%npoint))*grfraction(1)
             ! apl3 => print out particle number, wetscav, xscav_frac1
!              ! apl_8
!              if (ks.eq.1) write(*,924) 'wetdepo:',itime,jpart,wetscav,&
!                    zpoint2(npoint(jpart)),zpoint1(npoint(jpart)),&
!                    grfraction(1),xscav_frac1(jpart,ks)
! 924 format(a,1x,2i9,5es13.5)
          else
            part(jpart)%mass(ks)=0.
            xscav_frac1(jpart,ks)=0.
          endif
        endif
      endif

      if (wetscav.gt.0.) then
        wetdeposit(ks)=part(jpart)%mass(ks)* &
             (1.-exp(-wetscav*abs(ltsample)))*grfraction(1)  ! wet deposition
      else ! if no scavenging
        wetdeposit(ks)=0.
      endif
 
      restmass = part(jpart)%mass(ks)-wetdeposit(ks)
      if (ioutputforeachrelease.eq.1) then
        kp=part(jpart)%npoint
      else
        kp=1
      endif
      if (restmass .gt. smallnum) then
        part(jpart)%mass(ks)=restmass
!   depostatistic
!   wetdepo_sum(ks,kp)=wetdepo_sum(ks,kp)+wetdeposit(ks)
!   depostatistic
      else
        part(jpart)%mass(ks)=0.
      endif
!   Correct deposited mass to the last time step when radioactive decay of
!   gridded deposited mass was calculated
      if (decay(ks).gt.0.) then
        wetdeposit(ks)=wetdeposit(ks)*exp(abs(ldeltat)*decay(ks))
      endif

!    endif ! no deposition
    end do ! loop over species

! Sabine Eckhardt, June 2008 create deposition runs only for forward runs
! Add the wet deposition to accumulated amount on output grid and nested output grid
!*****************************************************************************

    if (ldirect.eq.1) then !OMP reduction necessary for wetgridunc
      call wetdepokernel(part(jpart)%nclass,wetdeposit,real(part(jpart)%xlon), &
           real(part(jpart)%ylat),nage,kp)
      if (nested_output.eq.1) call wetdepokernel_nest(part(jpart)%nclass, &
           wetdeposit,real(part(jpart)%xlon),real(part(jpart)%ylat),nage,kp)
    endif

  end do ! all particles

! count the total number of below-cloud and in-cloud occurences:
  tot_blc_count(1:nspec)=tot_blc_count(1:nspec)+blc_count(1:nspec)
  tot_inc_count(1:nspec)=tot_inc_count(1:nspec)+inc_count(1:nspec)

end subroutine wetdepo
