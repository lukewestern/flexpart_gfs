! SPDX-FileCopyrightText: FLEXPART 1998-2019, see flexpart_license.txt
! SPDX-License-Identifier: GPL-3.0-or-later

module output_mod
  
  use com_mod
  use par_mod
  use date_mod
#ifdef USE_NCF  
  use netcdf_output_mod
#endif
  use binary_output_mod
  use txt_output_mod

  implicit none

contains

subroutine initialise_output(itime,filesize)
  implicit none
  
  integer, intent(in) :: itime
  real, intent(inout) :: filesize
#ifdef USE_NCF
  real(kind=dp) ::          &
    jul
  integer ::                &
    jjjjmmdd,ihmmss,i
#endif

  ! Writing header information to either binary or NetCDF format
  if (itime.eq.0) then
    if (iout.ne.0) then
#ifdef USE_NCF
      if (lnetcdfout.eq.1) then 
        call writeheader_netcdf(lnest=.false.)
      else 
        call writeheader_binary
      end if

      if (nested_output.eq.1) then
        if (lnetcdfout.eq.1) then
          call writeheader_netcdf(lnest=.true.)
        else
          call writeheader_binary_nest
        endif
      endif
#endif
    endif

    call writeheader_binary ! CHECK ETA
    ! FLEXPART 9.2 ticket ?? write header in ASCII format 
    call writeheader_txt
    !if (nested_output.eq.1) call writeheader_nest
    if (nested_output.eq.1.and.surf_only.ne.1) call writeheader_binary_nest
    if (nested_output.eq.1.and.surf_only.eq.1) call writeheader_binary_nest_surf
    if (nested_output.ne.1.and.surf_only.eq.1) call writeheader_binary_surf

    ! NetCDF only: Create file for storing initial particle positions.
#ifdef USE_NCF
    if (mdomainfill.eq.0) then
      if (ldirect.eq.1) then
        call create_particles_initialoutput(ibtime,ibdate,ibtime,ibdate)
      else
        call create_particles_initialoutput(ietime,iedate,ietime,iedate)
      endif
    endif
    ! Create header files for files that store the particle output
    if (ipout.ge.1) then
      if (ldirect.eq.1) then
        call writeheader_partoutput(ibtime,ibdate,ibtime,ibdate)
      else 
        call writeheader_partoutput(ietime,iedate,ietime,iedate)
      endif
    endif
#endif
  
  ! In case the particle output file is becoming larger than the maximum set
  ! in par_mod, create a new one while keeping track of the filesize.
  else if ((mod(itime,ipoutfac*loutstep).eq.0).and.(ipout.ge.1)) then
#ifdef USE_NCF
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
#endif
  endif
end subroutine initialise_output

subroutine finalise_output(itime)
  ! Complete the calculation of initial conditions for particles not yet terminated
  
  implicit none 

  integer, intent(in) :: itime
  integer :: j

  do j=1,numpart
    if (linit_cond.ge.1) call initial_cond_calc(itime,j)
  end do

  if (ipout.eq.2) call output_particles(itime)!,active_per_rel)     ! dump particle positions

  if (linit_cond.ge.1) then
    if(linversionout.eq.1) then
      call initial_cond_output_inversion(itime)   ! dump initial cond. field
    else
      call initial_cond_output(itime)   ! dump initial cond. fielf
    endif
  endif
end subroutine finalise_output

subroutine output_particles(itime)
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
  use particle_mod
  use windfields_mod
  implicit none

  real(kind=dp) :: jul
  integer :: itime,i,j,jjjjmmdd,ihmmss
  integer :: ix,jy,ixp,jyp,indexh,m,il,ind,indz,indzp
  real :: dt1,dt2,dtt,ddx,ddy,rddx,rddy,p1,p2,p3,p4,dz1,dz2,dz
  real :: hm(2),pv1(2),pvprof(2),qv1(2),qvprof(2),pr1(2),prprof(2)
  real :: tt1(2),ttprof(2),rho1(2),rhoprof(2)
  real :: tr(2)
  character :: adate*8,atime*6

  real :: xlon(numpart),ylat(numpart),ztemp(numpart)
  real :: tti(numpart),rhoi(numpart),pvi(numpart),qvi(numpart),pri(numpart)
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
    pri(i)=-1.
    ztemp(i)=-1.
    !if (part(i)%itra1.eq.itime) then
      xlon(i)=xlon0+part(i)%xlon*dx
      ylat(i)=ylat0+part(i)%ylat*dy

  !*****************************************************************************
  ! Interpolate several variables (PV, specific humidity, etc.) to particle position
  !*****************************************************************************

      ix=part(i)%xlon
      jy=part(i)%ylat
      ixp=ix+1
      jyp=jy+1
      ddx=part(i)%xlon-real(ix)
      ddy=part(i)%ylat-real(jy)
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
        if (height(il).gt.part(i)%z) then
          indz=il-1
          indzp=il
          exit
        endif
      end do

      dz1=part(i)%z-height(indz)
      dz2=height(indzp)-part(i)%z
      dz=1./(dz1+dz2)
      ztemp(i)=part(i)%z

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
  ! Pressure
          pr1(m)=p1*prs(ix ,jy ,ind,indexh) &
               +p2*prs(ixp,jy ,ind,indexh) &
               +p3*prs(ix ,jyp,ind,indexh) &
               +p4*prs(ixp,jyp,ind,indexh)
        end do
        pvprof(ind-indz+1)=(pv1(1)*dt2+pv1(2)*dt1)*dtt
        qvprof(ind-indz+1)=(qv1(1)*dt2+qv1(2)*dt1)*dtt
        ttprof(ind-indz+1)=(tt1(1)*dt2+tt1(2)*dt1)*dtt
        rhoprof(ind-indz+1)=(rho1(1)*dt2+rho1(2)*dt1)*dtt
        prprof(ind-indz+1)=(pr1(1)*dt2+pr1(2)*dt1)*dtt
      end do
      pvi(i)=(dz1*pvprof(2)+dz2*pvprof(1))*dz
      qvi(i)=(dz1*qvprof(2)+dz2*qvprof(1))*dz
      tti(i)=(dz1*ttprof(2)+dz2*ttprof(1))*dz
      rhoi(i)=(dz1*rhoprof(2)+dz2*rhoprof(1))*dz
      pri(i)=(dz1*prprof(2)+dz2*prprof(1))*dz

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

    !endif 
  end do
  ! Determine current calendar date, needed for the file name
  !**********************************************************
  if (numpart.gt.0) then
    write(*,*) 'topo: ', topo(1), 'z:', part(1)%z
    write(*,*) 'xtra: ', xlon(1)
    write(*,*) 'ytra: ', ylat(1)
    !write(*,*) 'mass: ', xmass1(1,1)
    write(*,*) pvi(1),qvi(1),tti(1),rhoi(1)
  endif
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
    call partoutput_netcdf(itime,ztemp,'ZZ',j,ncid)
    call partoutput_netcdf(itime,topo,'TO',j,ncid)
    call partoutput_netcdf(itime,pvi,'PV',j,ncid)
    call partoutput_netcdf(itime,qvi,'QV',j,ncid)
    call partoutput_netcdf(itime,rhoi,'RH',j,ncid)
    call partoutput_netcdf(itime,hmixi,'HM',j,ncid)
    call partoutput_netcdf(itime,tri,'TR',j,ncid)
    call partoutput_netcdf(itime,tti,'TT',j,ncid)
    call partoutput_netcdf(itime,pri,'PR',j,ncid)
    do j=1,nspec
     ! call partoutput_netcdf(itime,xmass1(:,j),'MA',j,ncid)
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

      !if (itra1(i).eq.itime) then
    ! Write the output
    !*****************      
        !write(unitpartout) npoint(i),xlon(i),ylat(i),ztra1(i), &
        !     itramem(i),topo(i),pvi(i),qvi(i),rhoi(i),hmixi(i),tri(i),tti(i), &
        !     (xmass1(i,j),j=1,nspec)
      !endif
    end do


    write(unitpartout) -99999,-9999.9,-9999.9,-9999.9,-99999, &
         -9999.9,-9999.9,-9999.9,-9999.9,-9999.9,-9999.9,-9999.9, &
         (-9999.9,j=1,nspec)


    close(unitpartout)
  endif
end subroutine output_particles

subroutine output_concentrations(itime,loutstart,loutend,loutnext,outnum)
  use unc_mod
  use outg_mod
  use par_mod
  use com_mod
#ifdef USE_NCF
  use netcdf_output_mod, only: concoutput_netcdf,concoutput_nest_netcdf,&
       &concoutput_surf_netcdf,concoutput_surf_nest_netcdf
#endif
  use binary_output_mod 

  implicit none

  integer,intent(in) ::     &
    itime                     ! time index
  integer,intent(inout) ::  &
    loutstart,loutend,      & ! concentration calculation starting and ending time
    loutnext
  real,intent(inout) ::     &
    outnum                    ! concentration calculation sample number
  real(sp) ::               &
    gridtotalunc              ! concentration calculation related
  real(dep_prec) ::         &
    wetgridtotalunc,        & ! concentration calculation related
    drygridtotalunc           ! concentration calculation related
  real ::                   &
    weight                    ! concentration calculation sample weight

  ! Is the time within the computation interval, if not, return
  !************************************************************
  if ((ldirect*itime.lt.ldirect*loutstart).or.(ldirect*itime.gt.ldirect*loutend)) then
    return
  endif

  ! If we are exactly at the start or end of the concentration averaging interval,
  ! give only half the weight to this sample
  !*****************************************************************************
  if (mod(itime-loutstart,loutsample).eq.0) then
    if ((itime.eq.loutstart).or.(itime.eq.loutend)) then
      weight=0.5
    else
      weight=1.0
    endif
    outnum=outnum+weight
    call conccalc(itime,weight)
  endif

  ! If it is not time yet to write outputs, return
  !***********************************************
  if ((itime.ne.loutend).or.(outnum.le.0)) then
    return
  endif

  ! Output and reinitialization of grid
  ! If necessary, first sample of new grid is also taken
  !*****************************************************
  if ((iout.le.3.).or.(iout.eq.5)) then
    if (surf_only.ne.1) then 
#ifdef USE_NCF
      call concoutput_netcdf(itime,outnum,gridtotalunc,wetgridtotalunc,drygridtotalunc)
#else
      call concoutput(itime,outnum,gridtotalunc,wetgridtotalunc,drygridtotalunc)
#endif
    else
#ifdef USE_NCF
      call concoutput_surf_netcdf(itime,outnum,gridtotalunc,wetgridtotalunc,drygridtotalunc)
#else
      if (linversionout.eq.1) then
        call concoutput_inversion(itime,outnum,gridtotalunc,wetgridtotalunc,drygridtotalunc)
      else
        call concoutput_surf(itime,outnum,gridtotalunc,wetgridtotalunc,drygridtotalunc)
      endif
#endif
    endif

    if (nested_output .eq. 1) then
#ifdef USE_NCF
      if (surf_only.ne.1) then
        call concoutput_nest_netcdf(itime,outnum)
      else 
        call concoutput_surf_nest_netcdf(itime,outnum)
      endif
#else
      if (surf_only.ne.1) then
        call concoutput_nest(itime,outnum)
      else 
        if(linversionout.eq.1) then
          call concoutput_inversion_nest(itime,outnum)
        else 
          call concoutput_surf_nest(itime,outnum)
        endif
      endif
#endif
    endif
    outnum=0.
  endif

  write(*,45) itime,numpart,gridtotalunc,wetgridtotalunc,drygridtotalunc

45      format(i13,' Seconds simulated: ',i13, ' Particles:    Uncertainty: ',3f7.3)

  loutnext=loutnext+loutstep
  loutstart=loutnext-loutaver/2
  loutend=loutnext+loutaver/2
  if (itime.eq.loutstart) then
    weight=0.5
    outnum=outnum+weight
    call conccalc(itime,weight)
  endif
end subroutine output_concentrations

subroutine conccalc(itime,weight)
  !                      i     i
  !*****************************************************************************
  !                                                                            *
  !     Calculation of the concentrations on a regular grid using volume       *
  !     sampling                                                               *
  !                                                                            *
  !     Author: A. Stohl                                                       *
  !                                                                            *
  !     24 May 1996                                                            *
  !                                                                            *
  !     April 2000: Update to calculate age spectra                            *
  !                 Bug fix to avoid negative conc. at the domain boundaries,  *
  !                 as suggested by Petra Seibert                              *
  !                                                                            *
  !     2 July 2002: re-order if-statements in order to optimize CPU time      *
  !                                                                            *
  !                                                                            *
  !*****************************************************************************
  !                                                                            *
  ! Variables:                                                                 *
  ! nspeciesdim     = nspec for forward runs, 1 for backward runs              *
  !                                                                            *
  !*****************************************************************************

  use unc_mod
  use outg_mod
  use par_mod
  use com_mod
  use omp_lib, only: OMP_GET_THREAD_NUM
  use interpol_mod, only: interpol_density,ix,jy,ixp,jyp,ddx,ddy
  use coordinates_ecmwf
  use particle_mod

  implicit none

  integer,intent(in) :: itime
  real,intent(in) :: weight
  integer :: itage,i,kz,ks,n,nage
  integer :: il,ind,indz,indzp,nrelpointer
  real :: hx,hy,hz,h,xd,yd,zd,xkern,r2,c(maxspec)
  real :: rhoi
  real :: xl,yl,wx,wy,w
  real,parameter :: factor=.596831, hxmax=6.0, hymax=4.0, hzmax=150.
  !  integer xscav_count

  ! For forward simulations, make a loop over the number of species;
  ! for backward simulations, make an additional loop over the
  ! releasepoints
  !***************************************************************************
  !  xscav_count=0
  do i=1,numpart
    if (.not.part(i)%alive) cycle

  ! Determine age class of the particle
    itage=abs(itime-part(i)%tstart)
    do nage=1,nageclass
      if (itage.lt.lage(nage)) exit
    end do

  !  if (xscav_frac1(i,1).lt.0) xscav_count=xscav_count+1
           
  ! For special runs, interpolate the air density to the particle position
  !************************************************************************
  !***********************************************************************
  !AF IND_SOURCE switches between different units for concentrations at the source
  !Af    NOTE that in backward simulations the release of particles takes place
  !Af    at the receptor and the sampling at the source.
  !Af          1="mass"
  !Af          2="mass mixing ratio"
  !Af IND_RECEPTOR switches between different units for concentrations at the receptor
  !Af          1="mass"
  !Af          2="mass mixing ratio"

  !Af switches for the conccalcfile:
  !AF IND_SAMP =  0 : xmass * 1
  !Af IND_SAMP = -1 : xmass / rho

  !Af ind_samp is defined in readcommand.f

    if ( ind_samp .eq. -1 ) then
      call update_zeta_to_z(itime,i)
      call interpol_density(i,rhoi)
    elseif (ind_samp.eq.0) then 
      rhoi = 1.
    endif

  !****************************************************************************
  ! 1. Evaluate grid concentrations using a uniform kernel of bandwidths dx, dy
  !****************************************************************************


  ! For backward simulations, look from which release point the particle comes from
  ! For domain-filling trajectory option, npoint contains a consecutive particle
  ! number, not the release point information. Therefore, nrelpointer is set to 1
  ! for the domain-filling option.
  !*****************************************************************************

    if ((ioutputforeachrelease.eq.0).or.(mdomainfill.eq.1)) then
       nrelpointer=1
    else
       nrelpointer=part(i)%npoint
    endif

    do kz=1,numzgrid                ! determine height of cell
      if (outheight(kz).gt.part(i)%z) exit
    end do

    if (kz.le.numzgrid) then           ! inside output domain


  !********************************
  ! Do everything for mother domain
  !********************************

      xl=(part(i)%xlon*dx+xoutshift)/dxout
      yl=(part(i)%ylat*dy+youtshift)/dyout
      ix=int(xl)
      if (xl.lt.0.) ix=ix-1
      jy=int(yl)
      if (yl.lt.0.) jy=jy-1



  ! For particles aged less than 3 hours, attribute particle mass to grid cell
  ! it resides in rather than use the kernel, in order to avoid its smoothing effect.
  ! For older particles, use the uniform kernel.
  ! If a particle is close to the domain boundary, do not use the kernel either.
  !*****************************************************************************

      if ((.not.lusekerneloutput).or.(itage.lt.10800).or. &
           (xl.lt.0.5).or.(yl.lt.0.5).or. &
           (xl.gt.real(numxgrid-1)-0.5).or. &
           (yl.gt.real(numygrid-1)-0.5)) then             ! no kernel, direct attribution to grid cell
        if ((ix.ge.0).and.(jy.ge.0).and.(ix.le.numxgrid-1).and. &
             (jy.le.numygrid-1)) then
          if (DRYBKDEP.or.WETBKDEP) then
            do ks=1,nspec
              gridunc(ix,jy,kz,ks,nrelpointer,part(i)%nclass,nage)= &
                   gridunc(ix,jy,kz,ks,nrelpointer,part(i)%nclass,nage)+ &
                   part(i)%mass(ks)/rhoi*weight*max(xscav_frac1(i,ks),0.0)
            end do
          else
            if (lparticlecountoutput) then
              do ks=1,nspec
                gridunc(ix,jy,kz,ks,nrelpointer,part(i)%nclass,nage)= &
                     gridunc(ix,jy,kz,ks,nrelpointer,part(i)%nclass,nage)+1
              end do
            else
              do ks=1,nspec
                gridunc(ix,jy,kz,ks,nrelpointer,part(i)%nclass,nage)= &
                     gridunc(ix,jy,kz,ks,nrelpointer,part(i)%nclass,nage)+ &
                     part(i)%mass(ks)/rhoi*weight
              end do
            end if
          endif
        endif

      else                                 ! attribution via uniform kernel 

        ddx=xl-real(ix)                   ! distance to left cell border
        ddy=yl-real(jy)                   ! distance to lower cell border
        if (ddx.gt.0.5) then
          ixp=ix+1
          wx=1.5-ddx
        else
          ixp=ix-1
          wx=0.5+ddx
        endif

        if (ddy.gt.0.5) then
          jyp=jy+1
          wy=1.5-ddy
        else
          jyp=jy-1
          wy=0.5+ddy
        endif


  ! Determine mass fractions for four grid points
  !**********************************************

        if ((ix.ge.0).and.(ix.le.numxgrid-1)) then
          if ((jy.ge.0).and.(jy.le.numygrid-1)) then
            w=wx*wy
            if (DRYBKDEP.or.WETBKDEP) then
               do ks=1,nspec
                 gridunc(ix,jy,kz,ks,nrelpointer,part(i)%nclass,nage)= &
                   gridunc(ix,jy,kz,ks,nrelpointer,part(i)%nclass,nage)+ &
                   part(i)%mass(ks)/rhoi*w*weight*max(xscav_frac1(i,ks),0.0)
               end do
            else
               do ks=1,nspec
                 gridunc(ix,jy,kz,ks,nrelpointer,part(i)%nclass,nage)= &
                   gridunc(ix,jy,kz,ks,nrelpointer,part(i)%nclass,nage)+ &
                   part(i)%mass(ks)/rhoi*weight*w
               end do
            endif
          endif

          if ((jyp.ge.0).and.(jyp.le.numygrid-1)) then
            w=wx*(1.-wy)
            if (DRYBKDEP.or.WETBKDEP) then
              do ks=1,nspec
                 gridunc(ix,jyp,kz,ks,nrelpointer,part(i)%nclass,nage)= &
                   gridunc(ix,jyp,kz,ks,nrelpointer,part(i)%nclass,nage)+ &
                   part(i)%mass(ks)/rhoi*weight*w*max(xscav_frac1(i,ks),0.0)
               end do
             else
              do ks=1,nspec
                 gridunc(ix,jyp,kz,ks,nrelpointer,part(i)%nclass,nage)= &
                   gridunc(ix,jyp,kz,ks,nrelpointer,part(i)%nclass,nage)+ &
                   part(i)%mass(ks)/rhoi*weight*w
               end do
             endif
          endif
        endif !ix ge 0


        if ((ixp.ge.0).and.(ixp.le.numxgrid-1)) then
          if ((jyp.ge.0).and.(jyp.le.numygrid-1)) then
            w=(1.-wx)*(1.-wy)
            if (DRYBKDEP.or.WETBKDEP) then
               do ks=1,nspec
                 gridunc(ixp,jyp,kz,ks,nrelpointer,part(i)%nclass,nage)= &
                   gridunc(ixp,jyp,kz,ks,nrelpointer,part(i)%nclass,nage)+ &
                   part(i)%mass(ks)/rhoi*w*weight*max(xscav_frac1(i,ks),0.0)
               end do
            else
               do ks=1,nspec
                 gridunc(ixp,jyp,kz,ks,nrelpointer,part(i)%nclass,nage)= &
                   gridunc(ixp,jyp,kz,ks,nrelpointer,part(i)%nclass,nage)+ &
                   part(i)%mass(ks)/rhoi*weight*w
               end do
            endif
          endif

          if ((jy.ge.0).and.(jy.le.numygrid-1)) then
            w=(1.-wx)*wy
            if (DRYBKDEP.or.WETBKDEP) then
               do ks=1,nspec
                 gridunc(ixp,jy,kz,ks,nrelpointer,part(i)%nclass,nage)= &
                   gridunc(ixp,jy,kz,ks,nrelpointer,part(i)%nclass,nage)+ &
                   part(i)%mass(ks)/rhoi*weight*w*max(xscav_frac1(i,ks),0.0)
               end do
            else
               do ks=1,nspec
                 gridunc(ixp,jy,kz,ks,nrelpointer,part(i)%nclass,nage)= &
                   gridunc(ixp,jy,kz,ks,nrelpointer,part(i)%nclass,nage)+ &
                   part(i)%mass(ks)/rhoi*weight*w
               end do
            endif
          endif
        endif !ixp ge 0
     endif

  !************************************
  ! Do everything for the nested domain
  !************************************

      if (nested_output.eq.1) then
        xl=(part(i)%xlon*dx+xoutshiftn)/dxoutn
        yl=(part(i)%ylat*dy+youtshiftn)/dyoutn
        ix=int(xl)
        if (xl.lt.0.) ix=ix-1
        jy=int(yl)
        if (yl.lt.0.) jy=jy-1


  ! For particles aged less than 3 hours, attribute particle mass to grid cell
  ! it resides in rather than use the kernel, in order to avoid its smoothing effect.
  ! For older particles, use the uniform kernel.
  ! If a particle is close to the domain boundary, do not use the kernel either.
  !*****************************************************************************

        if ((itage.lt.10800).or.(xl.lt.0.5).or.(yl.lt.0.5).or. &
             (xl.gt.real(numxgridn-1)-0.5).or. &
             (yl.gt.real(numygridn-1)-0.5).or.((.not.lusekerneloutput))) then
  ! no kernel, direct attribution to grid cell
          if ((ix.ge.0).and.(jy.ge.0).and.(ix.le.numxgridn-1).and. &
               (jy.le.numygridn-1)) then
            if (DRYBKDEP.or.WETBKDEP) then
               do ks=1,nspec
                 griduncn(ix,jy,kz,ks,nrelpointer,part(i)%nclass,nage)= &
                   griduncn(ix,jy,kz,ks,nrelpointer,part(i)%nclass,nage)+ &
                   part(i)%mass(ks)/rhoi*weight*max(xscav_frac1(i,ks),0.0)
               end do
            else
              if (lparticlecountoutput) then
                do ks=1,nspec
                  griduncn(ix,jy,kz,ks,nrelpointer,part(i)%nclass,nage)= &
                       griduncn(ix,jy,kz,ks,nrelpointer,part(i)%nclass,nage)+1
                end do
              else
                do ks=1,nspec
                  griduncn(ix,jy,kz,ks,nrelpointer,part(i)%nclass,nage)= &
                       griduncn(ix,jy,kz,ks,nrelpointer,part(i)%nclass,nage)+ &
                       part(i)%mass(ks)/rhoi*weight
                end do
              endif
            endif
          endif
          
        else                                 ! attribution via uniform kernel

          ddx=xl-real(ix)                   ! distance to left cell border
          ddy=yl-real(jy)                   ! distance to lower cell border
          if (ddx.gt.0.5) then
            ixp=ix+1
            wx=1.5-ddx
          else
            ixp=ix-1
            wx=0.5+ddx
          endif

          if (ddy.gt.0.5) then
            jyp=jy+1
            wy=1.5-ddy
          else
            jyp=jy-1
            wy=0.5+ddy
          endif


  ! Determine mass fractions for four grid points
  !**********************************************

          if ((ix.ge.0).and.(ix.le.numxgridn-1)) then
            if ((jy.ge.0).and.(jy.le.numygridn-1)) then
              w=wx*wy
              if (DRYBKDEP.or.WETBKDEP) then
                 do ks=1,nspec
                   griduncn(ix,jy,kz,ks,nrelpointer,part(i)%nclass,nage)= &
                     griduncn(ix,jy,kz,ks,nrelpointer,part(i)%nclass,nage)+ &
                     part(i)%mass(ks)/rhoi*weight*w*max(xscav_frac1(i,ks),0.0)
                 end do
              else
                do ks=1,nspec
                   griduncn(ix,jy,kz,ks,nrelpointer,part(i)%nclass,nage)= &
                     griduncn(ix,jy,kz,ks,nrelpointer,part(i)%nclass,nage)+ &
                     part(i)%mass(ks)/rhoi*weight*w
                 end do
              endif
            endif

            if ((jyp.ge.0).and.(jyp.le.numygridn-1)) then
              w=wx*(1.-wy)
              if (DRYBKDEP.or.WETBKDEP) then
                 do ks=1,nspec
                   griduncn(ix,jyp,kz,ks,nrelpointer,part(i)%nclass,nage)= &
                     griduncn(ix,jyp,kz,ks,nrelpointer,part(i)%nclass,nage)+ &
                     part(i)%mass(ks)/rhoi*weight*w*max(xscav_frac1(i,ks),0.0)
                 end do
              else
                 do ks=1,nspec
                   griduncn(ix,jyp,kz,ks,nrelpointer,part(i)%nclass,nage)= &
                     griduncn(ix,jyp,kz,ks,nrelpointer,part(i)%nclass,nage)+ &
                     part(i)%mass(ks)/rhoi*weight*w
                 end do
              endif
            endif
          endif


          if ((ixp.ge.0).and.(ixp.le.numxgridn-1)) then
            if ((jyp.ge.0).and.(jyp.le.numygridn-1)) then
              w=(1.-wx)*(1.-wy)
              if (DRYBKDEP.or.WETBKDEP) then
                 do ks=1,nspec
                   griduncn(ixp,jyp,kz,ks,nrelpointer,part(i)%nclass,nage)= &
                     griduncn(ixp,jyp,kz,ks,nrelpointer,part(i)%nclass,nage)+ &
                     part(i)%mass(ks)/rhoi*weight*w*max(xscav_frac1(i,ks),0.0)
                 end do
              else
                 do ks=1,nspec
                   griduncn(ixp,jyp,kz,ks,nrelpointer,part(i)%nclass,nage)= &
                     griduncn(ixp,jyp,kz,ks,nrelpointer,part(i)%nclass,nage)+ &
                     part(i)%mass(ks)/rhoi*weight*w
                 end do
              endif
            endif

            if ((jy.ge.0).and.(jy.le.numygridn-1)) then
              w=(1.-wx)*wy
              if (DRYBKDEP.or.WETBKDEP) then
                 do ks=1,nspec
                   griduncn(ixp,jy,kz,ks,nrelpointer,part(i)%nclass,nage)= &
                     griduncn(ixp,jy,kz,ks,nrelpointer,part(i)%nclass,nage)+ &
                     part(i)%mass(ks)/rhoi*weight*w*max(xscav_frac1(i,ks),0.0)
                 end do
              else
                 do ks=1,nspec
                    griduncn(ixp,jy,kz,ks,nrelpointer,part(i)%nclass,nage)= &
                     griduncn(ixp,jy,kz,ks,nrelpointer,part(i)%nclass,nage)+ &
                     part(i)%mass(ks)/rhoi*weight*w
                 end do
              endif
            endif
          endif
        endif
      endif
    endif
  end do
  !  write(*,*) 'xscav count:',xscav_count

  !***********************************************************************
  ! 2. Evaluate concentrations at receptor points, using the kernel method
  !***********************************************************************

  do n=1,numreceptor


  ! Reset concentrations
  !*********************

    do ks=1,nspec
      c(ks)=0.
    end do


  ! Estimate concentration at receptor
  !***********************************

    do i=1,numpart

      if (.not. part(i)%alive) cycle
      itage=abs(itime-part(i)%tstart)

      hz=min(50.+0.3*sqrt(real(itage)),hzmax)
      zd=part(i)%z/hz
      if (zd.gt.1.) cycle          ! save computing time, leave loop

      hx=min((0.29+2.222e-3*sqrt(real(itage)))*dx+ &
           real(itage)*1.2e-5,hxmax)                     ! 80 km/day
      xd=(part(i)%xlon-xreceptor(n))/hx
      if (xd*xd.gt.1.) cycle       ! save computing time, leave loop

      hy=min((0.18+1.389e-3*sqrt(real(itage)))*dy+ &
           real(itage)*7.5e-6,hymax)                     ! 80 km/day
      yd=(part(i)%ylat-yreceptor(n))/hy
      if (yd*yd.gt.1.) cycle       ! save computing time, leave loop
      h=hx*hy*hz

      r2=xd*xd+yd*yd+zd*zd
      if (r2.lt.1.) then
        xkern=factor*(1.-r2)
        do ks=1,nspec
          c(ks)=c(ks)+part(i)%mass(ks)*xkern/h
        end do
      endif
    end do

    do ks=1,nspec
      creceptor(n,ks)=creceptor(n,ks)+2.*weight*c(ks)/receptorarea(n)
    end do
  end do
end subroutine conccalc

end module output_mod
