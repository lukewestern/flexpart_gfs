! SPDX-FileCopyrightText: FLEXPART 1998-2019, see flexpart_license.txt
! SPDX-License-Identifier: GPL-3.0-or-later

module oh_mod

  !includes OH concentration field as well as the height information
  !for this field

  implicit none

  integer :: nxOH,nyOH,nzOH
  real, allocatable, dimension(:) :: lonOH,latOH,altOH
  real, allocatable, dimension(:,:,:,:) :: OH_hourly
  real, allocatable, dimension (:,:,:,:) :: OH_field
  real, dimension(2) :: memOHtime
  real, dimension(360,180,12) :: jrate_average
  real, dimension(360) :: lonjr
  real, dimension(180) :: latjr

contains

subroutine ohreaction(itime,ltsample,loutnext)
  !                     i      i        i
  !*****************************************************************************
  !                                                                            *
  !                                                                            *
  !    Author: R.L. Thompson                                                   *
  !                                                                            *
  !    Nov 2014                                                                *
  !                                                                            *
  !                                                                            *
  !*****************************************************************************
  ! Variables:                                                                 *
  ! ix,jy                indices of output grid cell for each particle         *
  ! itime [s]            actual simulation time [s]                            *
  ! jpart                particle index                                        *
  ! ldeltat [s]          interval since radioactive decay was computed         *
  ! loutnext [s]         time for which gridded deposition is next output      *
  ! loutstep [s]         interval at which gridded deposition is output        *
  ! oh_average [molecule/cm^3] OH Concentration                                *
  ! ltsample [s]         interval over which mass is deposited                 *
  !                                                                            *
  !*****************************************************************************
  use par_mod
  use com_mod
  use particle_mod

  implicit none

  integer :: jpart,itime,ltsample,loutnext,ldeltat,j,k,ix,jy!,ijx,jjy
  integer :: ngrid,interp_time,m,n,h,indz,i!,ia,il
  integer :: jjjjmmdd,hhmmss,OHx,OHy,OHz
  real, dimension(nzOH) :: altOHtop
  real :: xlon,ylat 
  real :: xtn,ytn
  real :: restmass,ohreacted,oh_average
  real :: ohrate,temp 
  real, parameter :: smallnum = tiny(0.0) ! smallest number that can be handled
  real(kind=dp) :: jul

  ! Compute interval since radioactive decay of deposited mass was computed
  !************************************************************************

  if (itime.le.loutnext) then
    ldeltat=itime-(loutnext-loutstep)
  else                                  ! first half of next interval
    ldeltat=itime-loutnext
  endif

  jul=bdate+real(itime,kind=dp)/86400.
  call caldate(jul,jjjjmmdd,hhmmss)
  m=(jjjjmmdd-(jjjjmmdd/10000)*10000)/100
  h=hhmmss/10000

  ! Loop over particles
  !*****************************************
  ! !$OMP PARALLEL PRIVATE(jpart,xtn,ytn,j,k,ix,jy,interp_time, &
  ! !$OMP n,indz,i,xlon,ylat,OHx,OHy,OHz,oh_average,temp,ohrate, &
  ! !$OMP restmass,ohreacted,altOHtop,ngrid)

  
  do jpart=1,numpart

    ! Determine which nesting level to be used
    ngrid=0
    do j=numbnests,1,-1
      if ((part(jpart)%xlon.gt.xln(j)).and.(part(jpart)%xlon.lt.xrn(j)).and. &
           (part(jpart)%ylat.gt.yln(j)).and.(part(jpart)%ylat.lt.yrn(j))) then
        ngrid=j
        exit
      endif
    end do

    ! Determine nested grid coordinates
    if (ngrid.gt.0) then
      xtn=(part(jpart)%xlon-xln(ngrid))*xresoln(ngrid)
      ytn=(part(jpart)%ylat-yln(ngrid))*yresoln(ngrid)
      ix=int(xtn)
      jy=int(ytn)
    else
      ix=int(part(jpart)%xlon)
      jy=int(part(jpart)%ylat)
    endif

    interp_time=nint(itime-0.5*ltsample)
    n=2
    if(abs(memtime(1)-interp_time).lt.abs(memtime(2)-interp_time)) n=1

    do i=2,nz
      if (height(i).gt.part(jpart)%z) then
        indz=i-1
        exit
      endif
    end do

    ! Get OH from nearest grid-cell and specific month 
    !*************************************************

    ! world coordinates
    xlon=part(jpart)%xlon*dx+xlon0
    if (xlon.gt.180) then
       xlon=xlon-360
    endif
    ylat=part(jpart)%ylat*dy+ylat0

    ! get position in the OH field
    OHx=minloc(abs(lonOH-xlon),dim=1,mask=abs(lonOH-xlon).eq.minval(abs(lonOH-xlon)))
    OHy=minloc(abs(latOH-ylat),dim=1,mask=abs(latOH-ylat).eq.minval(abs(latOH-ylat)))

    ! get the level of the OH field for the particle
    ! z is the z-coord of the trajectory above model orography in metres
    ! altOH is the height of the centre of the level in the OH field above orography
    do i=2,nzOH
      altOHtop(i-1)=altOH(i)+0.5*(altOH(i)-altOH(i-1))
    end do
    altOHtop(nzOH)=altOH(nzOH)+0.5*(altOH(nzOH)-altOH(nzOH-1))
    OHz=minloc(abs(altOHtop-part(jpart)%z),dim=1,mask=abs(altOHtop-part(jpart)%z)&
            &.eq.minval(abs(altOHtop-part(jpart)%z)))

    ! Interpolate between hourly OH fields to current time
    !*****************************************************

    oh_average=OH_hourly(OHx,OHy,OHz,1)+&
               &(OH_hourly(OHx,OHy,OHz,2)-OH_hourly(OHx,OHy,OHz,1))*&
               &(itime-memOHtime(1))/(memOHtime(2)-memOHtime(1))

    if (oh_average.gt.smallnum) then

      ! Computation of the OH reaction
      !**********************************************************

      temp=tt(ix,jy,indz,n)

      do k=1,nspec                                
        if (ohcconst(k).gt.0.) then
          ohrate=ohcconst(k)*temp**ohnconst(k)*exp(-ohdconst(k)/temp)*oh_average
          ! new particle mass
          restmass = part(jpart)%mass(k)*exp(-1*ohrate*abs(ltsample))
          if (restmass .gt. smallnum) then
            part(jpart)%mass(k)=restmass
          else
            part(jpart)%mass(k)=0.
          endif
          ohreacted=part(jpart)%mass(k)*(1-exp(-1*ohrate*abs(ltsample)))
          if (jpart.eq.535) write(*,*) 'ohreaction', part(jpart)%mass(k),k
        else
          ohreacted=0.
        endif
      end do
    endif  ! oh_average.gt.smallnum 

  end do  !continue loop over all particles

  ! !$OMP END DO
  ! !$OMP END PARALLEL
end subroutine ohreaction

subroutine gethourlyOH(itime)
  !                     i     
  !*****************************************************************************
  !                                                                            *
  !                                                                            *
  !    Author: R.L. Thompson                                                   *
  !                                                                            *
  !    Nov 2014                                                                *
  !                                                                            *
  !                                                                            *
  !*****************************************************************************
  ! Variables:                                                                 *
  !                                                                            *
  !*****************************************************************************
  use par_mod
  use com_mod

  implicit none
  
  integer :: itime
  integer :: ix,jy,kz,m1,m2
  integer :: ijx,jjy
  integer :: jjjjmmdd,hhmmss
  real :: sza,jrate,photo_O1D,zenithangle
  real(kind=dp) :: jul1,jul2


  ! Check hourly OH field is available for the current time step
  !**************************************************************

  if ((ldirect*memOHtime(1).le.ldirect*itime).and. &
       (ldirect*memOHtime(2).gt.ldirect*itime)) then

  ! The right OH fields are already in memory -> don't do anything
  !****************************************************************

    continue

  else if ((ldirect*memOHtime(2).le.ldirect*itime).and. &
       (memOHtime(2).ne.0.)) then

    ! Current time is after 2nd OH field
    !************************************

    memOHtime(1)=memOHtime(2)
    memOHtime(2)=memOHtime(1)+ldirect*3600.
    OH_hourly(:,:,:,1)=OH_hourly(:,:,:,2)

    ! Compute new hourly value of OH
    !**********************************************************

    jul2=bdate+memOHtime(2)/86400._dp  ! date for next hour
    call caldate(jul2,jjjjmmdd,hhmmss)
    m2=(jjjjmmdd-(jjjjmmdd/10000)*10000)/100


    do kz=1,nzOH
      do jy=1,nyOH
        do ix=1,nxOH
          ijx=minloc(abs(lonjr-lonOH(ix)),dim=1,mask=abs(lonjr-lonOH(ix)).eq.minval(abs(lonjr-lonOH(ix))))
          jjy=minloc(abs(latjr-latOH(jy)),dim=1,mask=abs(latjr-latOH(jy)).eq.minval(abs(latjr-latOH(jy))))
          ! calculate solar zenith angle in degrees (sza) 
          sza=zenithangle(latOH(jy),lonOH(ix),jul2)
          ! calculate J(O1D) (jrate)
          jrate=photo_O1D(sza)
          ! apply hourly correction to OH
          if(jrate_average(ijx,jjy,m2).gt.0.) then
            OH_hourly(ix,jy,kz,2)=OH_field(ix,jy,kz,m2)*jrate/jrate_average(ijx,jjy,m2)
          else
            OH_hourly(ix,jy,kz,2)=0.
          endif
          !! for testing !!
          ! if(jy.eq.36.and.ix.eq.36.and.kz.eq.1) then
          !   write(999,fmt='(F6.3)') jrate/jrate_average(ijx,jjy,m2)
          ! endif
          ! if(jy.eq.11.and.ix.eq.36.and.kz.eq.1) then
          !   write(998,fmt='(F6.3)') jrate/jrate_average(ijx,jjy,m2)
          ! endif
        end do
      end do
    end do

  else

    ! No OH fields in memory -> compute both hourly OH fields
    !**********************************************************

    jul1=bdate  ! begin date of simulation (julian)
    call caldate(jul1,jjjjmmdd,hhmmss)
    m1=(jjjjmmdd-(jjjjmmdd/10000)*10000)/100
    memOHtime(1)=0.

    jul2=bdate+ldirect*real(1./24.,kind=dp)  ! date for next hour
    call caldate(jul2,jjjjmmdd,hhmmss)
    m2=(jjjjmmdd-(jjjjmmdd/10000)*10000)/100
    memOHtime(2)=ldirect*3600.


    do kz=1,nzOH
      do jy=1,nyOH
        do ix=1,nxOH
          ijx=minloc(abs(lonjr-lonOH(ix)),dim=1,mask=abs(lonjr-lonOH(ix)).eq.minval(abs(lonjr-lonOH(ix))))
          jjy=minloc(abs(latjr-latOH(jy)),dim=1,mask=abs(latjr-latOH(jy)).eq.minval(abs(latjr-latOH(jy))))
          ! calculate solar zenith angle in degrees (sza), beginning 
          sza=zenithangle(latOH(jy),lonOH(ix),jul1)
          ! calculate J(O1D) (jrate), beginning
          jrate=photo_O1D(sza)
          ! apply hourly correction to OH
          if(jrate_average(ijx,jjy,m1).gt.0.) then
            OH_hourly(ix,jy,kz,1)=OH_field(ix,jy,kz,m1)*jrate/jrate_average(ijx,jjy,m1)
          else
            OH_hourly(ix,jy,kz,1)=0.
          endif
          ! calculate solar zenith angle in degrees (sza), after 1-hour 
          sza=zenithangle(latOH(jy),lonOH(ix),jul2)
          ! calculate J(O1D) (jrate), after 1-hour
          jrate=photo_O1D(sza)
          ! apply hourly correction to OH
          if(jrate_average(ijx,jjy,m2).gt.0.) then
            OH_hourly(ix,jy,kz,2)=OH_field(ix,jy,kz,m2)*jrate/jrate_average(ijx,jjy,m2)
          else
            OH_hourly(ix,jy,kz,2)=0.
          endif
        end do
      end do
    end do

  endif
end subroutine gethourlyOH

end module oh_mod
