! SPDX-FileCopyrightText: FLEXPART 1998-2019, see flexpart_license.txt
! SPDX-License-Identifier: GPL-3.0-or-later

subroutine partoutput(itime)
  !                        i
  !*****************************************************************************
  !                                                                            *
  !     Dump all particle positions                                            *
  !                                                                            *
  !     Author: A. Stohl                                                       *
  !                                                                            *
  !     12 March 1999                                                          *
  !                                                                            *
  !*****************************************************************************
  !                                                                            *
  ! Variables:                                                                 *
  !                                                                            *
  !*****************************************************************************

  use par_mod
  use com_mod
#ifdef USE_NCF
  use netcdf
  use netcdf_output_mod, only: partoutput_netcdf,open_partoutput_file,close_partoutput_file
  use omp_lib, only: OMP_GET_THREAD_NUM
#endif

  implicit none

  real(kind=dp) :: jul
  integer :: itime,i,j,jjjjmmdd,ihmmss
  integer :: ix,jy,ixp,jyp,indexh,m,il,ind,indz,indzp
  real :: dt1,dt2,dtt,ddx,ddy,rddx,rddy,p1,p2,p3,p4,dz1,dz2,dz
  real :: hm(2),pv1(2),pvprof(2),qv1(2),qvprof(2)
  real :: tt1(2),ttprof(2),rho1(2),rhoprof(2)
  real :: tr(2)
  character :: adate*8,atime*6

  real :: xlon(numpart),ylat(numpart)
  real :: tti(numpart),rhoi(numpart),pvi(numpart),qvi(numpart)
  real :: topo(numpart),hmixi(numpart),tri(numpart)

#ifdef USE_NCF
  integer  :: ncid
#endif


  ! Some variables needed for temporal interpolation
  !*************************************************

  dt1=real(itime-memtime(1))
  dt2=real(memtime(2)-itime)
  dtt=1./(dt1+dt2)


  do i=1,numpart
  ! Take only valid particles
  !**************************
    xlon(i)=-1.
    ylat(i)=-1.
    tti(i)=-1.
    rhoi(i)=-1.
    pvi(i)=-1.
    qvi(i)=-1.
    topo(i)=-1.
    hmixi(i)=-1.
    tri(i)=-1.
    if (itra1(i).eq.itime) then
      xlon(i)=xlon0+xtra1(i)*dx
      ylat(i)=ylat0+ytra1(i)*dy

  !*****************************************************************************
  ! Interpolate several variables (PV, specific humidity, etc.) to particle position
  !*****************************************************************************

      ix=xtra1(i)
      jy=ytra1(i)
      ixp=ix+1
      jyp=jy+1
      ddx=xtra1(i)-real(ix)
      ddy=ytra1(i)-real(jy)
      rddx=1.-ddx
      rddy=1.-ddy
      p1=rddx*rddy
      p2=ddx*rddy
      p3=rddx*ddy
      p4=ddx*ddy

! eso: Temporary fix for particle exactly at north pole
      if (jyp >= nymax) then
      !  write(*,*) 'WARNING: conccalc.f90 jyp >= nymax'
        jyp=jyp-1
      end if

  ! Topography
  !***********

      topo(i)=p1*oro(ix ,jy) &
           + p2*oro(ixp,jy) &
           + p3*oro(ix ,jyp) &
           + p4*oro(ixp,jyp)


  ! Potential vorticity, specific humidity, temperature, and density
  !*****************************************************************

      do il=2,nz
        if (height(il).gt.ztra1(i)) then
          indz=il-1
          indzp=il
          exit
        endif
      end do

      dz1=ztra1(i)-height(indz)
      dz2=height(indzp)-ztra1(i)
      dz=1./(dz1+dz2)


      do ind=indz,indzp
        do m=1,2
          indexh=memind(m)

  ! Potential vorticity
          pv1(m)=p1*pv(ix ,jy ,ind,indexh) &
               +p2*pv(ixp,jy ,ind,indexh) &
               +p3*pv(ix ,jyp,ind,indexh) &
               +p4*pv(ixp,jyp,ind,indexh)
  ! Specific humidity
          qv1(m)=p1*qv(ix ,jy ,ind,indexh) &
               +p2*qv(ixp,jy ,ind,indexh) &
               +p3*qv(ix ,jyp,ind,indexh) &
               +p4*qv(ixp,jyp,ind,indexh)
  ! Temperature
          tt1(m)=p1*tt(ix ,jy ,ind,indexh) &
               +p2*tt(ixp,jy ,ind,indexh) &
               +p3*tt(ix ,jyp,ind,indexh) &
               +p4*tt(ixp,jyp,ind,indexh)
  ! Density
          rho1(m)=p1*rho(ix ,jy ,ind,indexh) &
               +p2*rho(ixp,jy ,ind,indexh) &
               +p3*rho(ix ,jyp,ind,indexh) &
               +p4*rho(ixp,jyp,ind,indexh)
        end do
        pvprof(ind-indz+1)=(pv1(1)*dt2+pv1(2)*dt1)*dtt
        qvprof(ind-indz+1)=(qv1(1)*dt2+qv1(2)*dt1)*dtt
        ttprof(ind-indz+1)=(tt1(1)*dt2+tt1(2)*dt1)*dtt
        rhoprof(ind-indz+1)=(rho1(1)*dt2+rho1(2)*dt1)*dtt
      end do
      pvi(i)=(dz1*pvprof(2)+dz2*pvprof(1))*dz
      qvi(i)=(dz1*qvprof(2)+dz2*qvprof(1))*dz
      tti(i)=(dz1*ttprof(2)+dz2*ttprof(1))*dz
      rhoi(i)=(dz1*rhoprof(2)+dz2*rhoprof(1))*dz

  ! Tropopause and PBL height
  !**************************

      do m=1,2
        indexh=memind(m)

  ! Tropopause
        tr(m)=p1*tropopause(ix ,jy ,1,indexh) &
             + p2*tropopause(ixp,jy ,1,indexh) &
             + p3*tropopause(ix ,jyp,1,indexh) &
             + p4*tropopause(ixp,jyp,1,indexh)

  ! PBL height
        hm(m)=p1*hmix(ix ,jy ,1,indexh) &
             + p2*hmix(ixp,jy ,1,indexh) &
             + p3*hmix(ix ,jyp,1,indexh) &
             + p4*hmix(ixp,jyp,1,indexh)
      end do

      hmixi(i)=(hm(1)*dt2+hm(2)*dt1)*dtt
      tri(i)=(tr(1)*dt2+tr(2)*dt1)*dtt

    endif 
  end do
  ! Determine current calendar date, needed for the file name
  !**********************************************************

  jul=bdate+real(itime,kind=dp)/86400._dp
  call caldate(jul,jjjjmmdd,ihmmss)
  write(adate,'(i8.8)') jjjjmmdd
  write(atime,'(i6.6)') ihmmss

  if (lnetcdfout.eq.1) then
#ifdef USE_NCF
  ! open output file
    call open_partoutput_file(ncid)
    ! First allocate the time and particle dimention within the netcdf file
    call partoutput_netcdf(itime,xlon,'TI',j,ncid)
    call partoutput_netcdf(itime,xlon,'PA',j,ncid)
    call partoutput_netcdf(itime,xlon,'LO',j,ncid)
    call partoutput_netcdf(itime,ylat,'LA',j,ncid)
    call partoutput_netcdf(itime,ztra1,'ZZ',j,ncid)
    call partoutput_netcdf(itime,topo,'TO',j,ncid)
    call partoutput_netcdf(itime,pvi,'PV',j,ncid)
    call partoutput_netcdf(itime,qvi,'QV',j,ncid)
    call partoutput_netcdf(itime,rhoi,'RH',j,ncid)
    call partoutput_netcdf(itime,hmixi,'HM',j,ncid)
    call partoutput_netcdf(itime,tri,'TR',j,ncid)
    call partoutput_netcdf(itime,tti,'TT',j,ncid)
    do j=1,nspec
      call partoutput_netcdf(itime,xmass1(:,j),'MA',j,ncid)
    end do
    call close_partoutput_file(ncid)
#endif
  else
    ! Open output file and write the output
    !**************************************

    if (ipout.eq.1.or.ipout.eq.3) then
      open(unitpartout,file=path(2)(1:length(2))//'partposit_'//adate// &
           atime,form='unformatted')
    else
      open(unitpartout,file=path(2)(1:length(2))//'partposit_end', &
           form='unformatted')
    endif

    ! Write current time to file
    !***************************

    write(unitpartout) itime
    do i=1,numpart
    ! Take only valid particles
    !**************************

      if (itra1(i).eq.itime) then
    ! Write the output
    !*****************      
        write(unitpartout) npoint(i),xlon(i),ylat(i),ztra1(i), &
             itramem(i),topo(i),pvi(i),qvi(i),rhoi(i),hmixi(i),tri(i),tti(i), &
             (xmass1(i,j),j=1,nspec)
      endif
    end do


    write(unitpartout) -99999,-9999.9,-9999.9,-9999.9,-99999, &
         -9999.9,-9999.9,-9999.9,-9999.9,-9999.9,-9999.9,-9999.9, &
         (-9999.9,j=1,nspec)


    close(unitpartout)
  endif

end subroutine partoutput
