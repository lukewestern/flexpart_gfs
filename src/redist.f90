! SPDX-FileCopyrightText: FLEXPART 1998-2019, see flexpart_license.txt
! SPDX-License-Identifier: GPL-3.0-or-later

subroutine redist (itime,ipart,ktop,ipconv)

  !**************************************************************************
  ! Do the redistribution of particles due to convection
  ! This subroutine is called for each particle which is assigned
  ! a new vertical position randomly, based on the convective redistribution
  ! matrix
  !**************************************************************************

  ! Petra Seibert, Feb 2001, Apr 2001, May 2001, Jan 2002, Nov 2002 and
  ! Andreas Frank, Nov 2002

  ! Caroline Forster:  November 2004 - February 2005

  use par_mod
  use com_mod
  use conv_mod
  use random_mod
  use omp_lib
  use interpol_mod
  use coordinates_ecmwf

  implicit none

  real,parameter :: const=r_air/ga
  integer :: ipart, ktop,ipconv,itime
  integer :: k, kz, levnew, levold
  real :: uvzlev(nuvzmax)
  real :: wsub(nuvzmax)
  real :: totlevmass, wsubpart
  real :: temp_levold,temp_levold1
  real :: sub_levold,sub_levold1
  real :: pint, pold, rn, tv, tvold, dlevfrac
  real :: ew,ztold,ffraction
  real :: tv1, tv2, dlogp, dz, dz1, dz2
  save :: iseed, uvzlev
  integer :: iseed = -88

!$OMP THREADPRIVATE(iseed,uvzlev)
!$  if (iseed.eq.-88) then
!$    iseed = iseed - OMP_GET_THREAD_NUM()
!$  endif

  ! ipart   ... number of particle to be treated

  ipconv=1

  !  determine vertical grid position of particle in the eta system
  !****************************************************************
  select case (wind_coord_type)

    case ('ETA')
      ztold = ztra1eta(abs(ipart))
      ! find old particle grid position
      levold = nconvtop
      do kz = 2, nconvtop
        if (uvheight(kz) .le. ztold ) then
          levold = kz-1
          exit
        endif
      end do

    case ('METER')

      ! determine height of the eta half-levels (uvzlev)
      ! do that only once for each grid column
      ! i.e. when ktop.eq.1
      !**************************************************************

      if (ktop .le. 1) then

        tvold=tt2conv*(1.+0.378*ew(td2conv)/psconv)
        pold=psconv
        uvzlev(1)=0.

        pint = phconv(2)
      !  determine next virtual temperatures
        tv1 = tconv(1)*(1.+0.608*qconv(1))
        tv2 = tconv(2)*(1.+0.608*qconv(2))
      !  interpolate virtual temperature to half-level
        tv = tv1 + (tv2-tv1)*(pconv(1)-phconv(2))/(pconv(1)-pconv(2))
        tv = tv1 + (tv2-tv1)*(pconv(1)-phconv(2))/(pconv(1)-pconv(2))
        if (abs(tv-tvold).gt.0.2) then
          uvzlev(2) = uvzlev(1) + &
               const*log(pold/pint)* &
               (tv-tvold)/log(tv/tvold)
        else
          uvzlev(2) = uvzlev(1)+ &
               const*log(pold/pint)*tv
        endif
        tvold=tv
        tv1=tv2
        pold=pint

      ! integrate profile (calculation of height agl of eta layers) as required
        do kz = 3, nconvtop+1
      !    note that variables defined in calcmatrix.f (pconv,tconv,qconv)
      !    start at the first real ECMWF model level whereas kz and
      !    thus uvzlev(kz) starts at the surface. uvzlev is defined at the
      !    half-levels (between the tconv, qconv etc. values !)
      !    Thus, uvzlev(kz) is the lower boundary of the tconv(kz) cell.
          pint = phconv(kz)
      !    determine next virtual temperatures
          tv2 = tconv(kz)*(1.+0.608*qconv(kz))
      !    interpolate virtual temperature to half-level
          tv = tv1 + (tv2-tv1)*(pconv(kz-1)-phconv(kz))/ &
               (pconv(kz-1)-pconv(kz))
          tv = tv1 + (tv2-tv1)*(pconv(kz-1)-phconv(kz))/ &
               (pconv(kz-1)-pconv(kz))
          if (abs(tv-tvold).gt.0.2) then
            uvzlev(kz) = uvzlev(kz-1) + &
                 const*log(pold/pint)* &
                 (tv-tvold)/log(tv/tvold)
          else
            uvzlev(kz) = uvzlev(kz-1)+ &
                 const*log(pold/pint)*tv
          endif
          tvold=tv
          tv1=tv2
          pold=pint


        end do

        ktop = 2

      endif
      
      ztold = ztra1(abs(ipart))
      ! find old particle grid position
      levold = nconvtop
      do kz = 2, nconvtop
        if (uvzlev(kz) .ge. ztold ) then
          levold = kz-1
          exit
        endif
      end do

    case default
      write(*,*) 'The wind_coord_type is not defined in redist.f90'
      stop

  end select

  ! If the particle is above the potentially convective domain, it will be skipped
  if (levold.ne.nconvtop) then

    ! now redistribute particles
    !****************************

    !  Choose a random number and find corresponding level of destination
    !  Random numbers to be evenly distributed in [0,1]

    rn = ran3(iseed)

    ! initialize levnew

    levnew = levold

    ffraction = 0.
    totlevmass=dpr(levold)/ga
    loop1: do k = 1,nconvtop
    ! for backward runs use the transposed matrix
     if (ldirect.eq.1) then
       ffraction=ffraction+fmassfrac(levold,k) &
            /totlevmass
     else
       ffraction=ffraction+fmassfrac(k,levold) &
            /totlevmass
     endif
     if (rn.le.ffraction) then
       levnew=k
    ! avoid division by zero or a too small number
    ! if division by zero or a too small number happens the
    ! particle is assigned to the center of the grid cell
       if (ffraction.gt.1.e-20) then
        if (ldirect.eq.1) then
          dlevfrac = (ffraction-rn) / fmassfrac(levold,k) * totlevmass
        else
          dlevfrac = (ffraction-rn) / fmassfrac(k,levold) * totlevmass
        endif
       else
         dlevfrac = 0.5
       endif
       exit loop1
     endif
    end do loop1

    ! now assign new position to particle
    select case (wind_coord_type)

      case ('ETA')
        if (levnew.le.nconvtop) then
          if (levnew.eq.levold) then
            ztra1eta(abs(ipart)) = ztold
          else
            dlogp = (1.-dlevfrac)* &
               (uvheight(levnew+1)-uvheight(levnew))
            ztra1eta(abs(ipart)) = uvheight(levnew)+dlogp
            if (ztra1eta(abs(ipart)).ge.1.) ztra1eta(abs(ipart))=1.-(ztra1eta(abs(ipart))-1.)
            if (ztra1eta(abs(ipart)).eq.1.) ztra1eta(abs(ipart))=ztra1eta(abs(ipart))-1.e-5
            if (ipconv.gt.0) ipconv=-1
          endif
        endif

      case ('METER')
        if (levnew.le.nconvtop) then
          if (levnew.eq.levold) then
            ztra1(abs(ipart)) = ztold
          else
            dlogp = (1.-dlevfrac)* &
                 (log(phconv(levnew+1))-log(phconv(levnew)))
            pint = log(phconv(levnew))+dlogp
            dz1 = pint - log(phconv(levnew))
            dz2 = log(phconv(levnew+1)) - pint
            dz = dz1 + dz2
            ztra1(abs(ipart)) = (uvzlev(levnew)*dz2+uvzlev(levnew+1)*dz1)/dz
            if (ztra1(abs(ipart)).lt.0.) &
                ztra1(abs(ipart))=-1.*ztra1(abs(ipart))
            if (ipconv.gt.0) ipconv=-1
          endif
        endif

      case default
        write(*,*) 'The wind_coord_type is not defined in redist.f90'
        stop

    end select

    ! displace particle according to compensating subsidence
    ! this is done to those particles, that were not redistributed
    ! by the matrix
    !**************************************************************

    if (levnew.le.nconvtop.and.levnew.eq.levold) then

      ! determine compensating vertical velocity at the levels
      ! above and below the particel position
      ! increase compensating subsidence by the fraction that
      ! is displaced by convection to this level

      if (levold.gt.1) then
       temp_levold = tconv(levold-1) + &
            (tconv(levold)-tconv(levold-1)) &
            *(pconv(levold-1)-phconv(levold))/ &
            (pconv(levold-1)-pconv(levold))
       sub_levold = sub(levold)/(1.-ga*sub(levold)/dpr(levold))
       wsub(levold)=-1.*sub_levold*r_air*temp_levold/(phconv(levold))
      else
       wsub(levold)=0.
      endif

       temp_levold1 = tconv(levold) + &
            (tconv(levold+1)-tconv(levold)) &
            *(pconv(levold)-phconv(levold+1))/ &
            (pconv(levold)-pconv(levold+1))
       sub_levold1 = sub(levold+1)/(1.-ga*sub(levold+1)/dpr(levold+1))
       wsub(levold+1)=-1.*sub_levold1*r_air*temp_levold1/ &
            (phconv(levold+1))

      ! interpolate wsub to the vertical particle position
      select case (wind_coord_type)
        case ('ETA')
          ztold = ztra1eta(abs(ipart))
          dz1 = ztold - uvheight(levold)
          dz2 = uvheight(levold+1) - ztold
          dz = dz1 + dz2

          ! Convert z(eta) to z(m) in order to add subsidence
          call zeta_to_z(itime,xtra1(abs(ipart)),ytra1(abs(ipart)), &
            ztra1eta(abs(ipart)),ztra1(abs(ipart)))

          ztold=ztra1(abs(ipart))
          wsubpart = (dz2*wsub(levold)+dz1*wsub(levold+1))/dz
          
          ztra1(abs(ipart)) = ztold+wsubpart*real(lsynctime)

          if (ztra1(abs(ipart)).lt.0.) then
             ztra1(abs(ipart))=-1.*ztra1(abs(ipart))
          endif

          ! Convert new z(m) back to z(eta)
          call z_to_zeta(itime,xtra1(abs(ipart)),ytra1(abs(ipart)), &
            ztra1(abs(ipart)),ztra1eta(abs(ipart)))

        case ('METER')
          ztold = ztra1(abs(ipart))
          dz1 = ztold - uvzlev(levold)
          dz2 = uvzlev(levold+1) - ztold
          dz = dz1 + dz2

          wsubpart = (dz2*wsub(levold)+dz1*wsub(levold+1))/dz
          ztra1(abs(ipart)) = ztold+wsubpart*real(lsynctime)
          if (ztra1(abs(ipart)).lt.0.) then
             ztra1(abs(ipart))=-1.*ztra1(abs(ipart))
          endif

        case default
          write(*,*) 'The wind_coord_type is not defined in redist.f90'
          stop
      end select
    endif      !(levnew.le.nconvtop.and.levnew.eq.levold)
  endif
  ! Maximum altitude .5 meter below uppermost model level
  !*******************************************************

  select case (wind_coord_type)
    case ('ETA')
      if (ztra1(abs(ipart)) .gt. height(nz)-0.5) &
           ztra1(abs(ipart)) = height(nz)-0.5

      if (ztra1eta(abs(ipart)) .lt. uvheight(nz)) &
           ztra1eta(abs(ipart)) = uvheight(nz)+1.e-4
      if (ztra1eta(abs(ipart)).ge.1.) ztra1eta(abs(ipart))=1.-(ztra1eta(abs(ipart))-1.)
      if (ztra1eta(abs(ipart)).eq.1.) ztra1eta(abs(ipart))=ztra1eta(abs(ipart))-1.e-4
    case ('METER')
      if (ztra1(abs(ipart)) .gt. height(nz)-0.5) &
           ztra1(abs(ipart)) = height(nz)-0.5
    case default
      write(*,*) 'The wind_coord_type is not defined in redist.f90'
      stop
  end select

end subroutine redist
