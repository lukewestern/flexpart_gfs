! SPDX-FileCopyrightText: FLEXPART 1998-2019, see flexpart_license.txt
! SPDX-License-Identifier: GPL-3.0-or-later

module receptor_mod

  !*****************************************************************************
  !                                                                            *
  !    This module contains variables and subroutines for  calculating         *
  !    receptor concentrations and writing these to file                       *
  !                                                                            *
  !*****************************************************************************

  use par_mod
  use com_mod
  use point_mod
  use particle_mod
  use date_mod
  use windfields_mod,      only: rho, prs, height, nzmax, nz !, qv
#if USE_NCF
  use receptor_netcdf_mod, only: nc_id, concvar_id, uncvar_id, &
                                recnamevar_id, timevar_id, &
                                reclonvar_id, reclatvar_id, recaltvar_id, &
                                nnvar_id, xkvar_id, rpointer, &
                                ncsat_id, satvar_id, satuncvar_id, &
                                satnamevar_id, sattimevar_id, &
                                satlonvar_id, satlatvar_id, sataltvar_id, &
                                satnnvar_id, satxkvar_id, spointer, sat_name, &
                                receptor_output_netcdf, write_receptor_netcdf, &
                                satellite_output_netcdf, write_satellite_netcdf, &
                                close_receptor_netcdf, close_satellite_netcdf
#endif
  use binary_output_mod,  only: receptorout_init_binary, write_receptor_binary, &
                                satelliteout_init_binary, write_satellite_binary

  implicit none

  contains

 
  subroutine alloc_receptor

  !*****************************************************************************
  !                                                                            *
  !    This routine allocates variables for receptor concentrations            *
  !                                                                            *
  !*****************************************************************************

    implicit none

    if (numreceptor.gt.0) then
      allocate( creceptor(nspec,maxrecsample) )
      allocate( crecuncert(nspec,maxrecsample) )
      allocate( nnreceptor(maxrecsample) )
      allocate( xkreceptor(maxrecsample) )
      creceptor(:,:)=0.
      crecuncert(:,:)=0.
      nnreceptor(:)=0.
      xkreceptor(:)=0.
    endif

    if (numsatreceptor.gt.0) then
      allocate( csatellite(nspec,nlayermax,maxrecsample) )
      allocate( csatuncert(nspec,nlayermax,maxrecsample) )
      allocate( nnsatellite(nlayermax,maxrecsample) )
      allocate( xksatellite(nlayermax,maxrecsample) )
      csatellite(:,:,:)=0.
      csatuncert(:,:,:)=0.
      nnsatellite(:,:)=0.
      xksatellite(:,:)=0.
    endif

  end subroutine alloc_receptor


  subroutine receptoroutput(itime,lrecoutstart,lrecoutend,lrecoutnext,recoutnum,recoutnumsat)

  !*****************************************************************************
  !                                                                            *
  !    Output concentrations at receptors                                      *
  !                                                                            *
  !    Author: R. Thompson, Sep-2023                                           *
  !                                                                            *
  !                                                                            *
  !*****************************************************************************

    implicit none

    integer   :: itime, lrecoutstart, lrecoutend, lrecoutnext
    real, dimension(maxrecsample) :: recoutnum
    real, dimension(nlayermax,maxrecsample) :: recoutnumsat
    real      :: weight
    real(kind=dp) :: jul
    character(len=256) :: fn, fnsat
    character :: adate*8, atime*6
    integer   :: jjjjmmdd, ihmmss
    integer   :: numsatlayer, nchar, ks_start
    integer   :: n, nn, ix, jy, ixp, jyp, indz, indzp, il, ind, ks, k, kz
    real      :: ddx, ddy, rddx, rddy, p1, p2, p3, p4, dz1, dz2, dz
    real      :: rhoi, zmid
    real, dimension(2) :: rho_p
    real, dimension(:), allocatable   :: densityoutrecept, nnrec, xkrec
    real, dimension(:,:), allocatable :: crec, cunc
    real, parameter :: weightair=28.97
    integer, parameter :: unitoutrecdates=109
    logical :: lexist

    if (llcmoutput) then
      ks_start=2
    else
      ks_start=1
    endif

    if (mod(itime-lrecoutstart,lrecoutsample).eq.0) then
      if ((itime.eq.lrecoutstart).or.(itime.eq.lrecoutend)) then
        weight=0.5
      else
        weight=1.0
      endif
      write(*,*) 'calling receptorcalc at itime = ',itime
      call receptorcalc(itime,weight,lrecoutstart,lrecoutend,recoutnum,recoutnumsat)
    endif

    if ((itime.eq.lrecoutend).and.(any(recoutnum.gt.0.).or.any(recoutnumsat.gt.0.))) then

      ! output receptor concentrations
      !*******************************

      if ((iout.le.3.).or.(iout.eq.5)) then
        write(*,*) 'calling receptoroutput at itime = ',itime

        ! Determine current date for output in dates_receptor file
        !**********************************************************

        jul=bdate+dble(float(itime))/86400.
        call caldate(jul,jjjjmmdd,ihmmss)
        write(adate,'(i8.8)') jjjjmmdd
        write(atime,'(i6.6)') ihmmss

        inquire(file=path(2)(1:length(2))//'dates_receptors',exist=lexist)
        if (.not.lexist) then
          ! initialize dates output file
          open(unitoutrecdates,file=path(2)(1:length(2))//'dates_receptors',STATUS='REPLACE')
          close(unitoutrecdates)
        endif
        open(unitoutrecdates,file=path(2)(1:length(2))//'dates_receptors', &
              ACCESS='APPEND', STATUS='OLD')
        write(unitoutrecdates,'(a)') adate//atime
        close(unitoutrecdates)

        ! For netcdf output open files and get variable info
        !***************************************************

        if (lnetcdfout.eq.1) then
#ifdef USE_NCF
          if (numreceptor.gt.0) then
            print*, 'before: nc_id, concvar_id = ',nc_id, concvar_id
            call receptor_output_netcdf()
            print*, 'after: nc_id, concvar_id = ',nc_id, concvar_id
          endif
          if (numsatreceptor.gt.0) then
            print*, 'before: ncsat_id, satvar_id = ',ncsat_id, satvar_id
            call satellite_output_netcdf()
            print*, 'after: ncsat_id, satvar_id = ',ncsat_id, satvar_id
          endif
#endif
        endif

        ! Initialize variables
        !*********************

        allocate(densityoutrecept(nlayermax),nnrec(nlayermax),xkrec(nlayermax))
        allocate(crec(nspec,nlayermax),cunc(nspec,nlayermax))
        nnrec(:)=0.
        xkrec(:)=0.
        crec(:,:)=0.
        cunc(:,:)=0.

        ! Loop over general receptors
        !**************************** 

!$OMP PARALLEL &
!$OMP PRIVATE(n,k,ix,jy,ixp,jyp,ddx,ddy,rddx,rddy,p1,p2,p3,p4,indz,indzp,il,dz1,dz2,dz, &
!$OMP         ind,rho_p,rhoi,densityoutrecept,ks,crec,cunc,nnrec,xkrec) &
!$OMP SHARED(rpointer) 

!$OMP DO
        do n=1,numreceptor

          if ((treceptor(n).lt.lrecoutstart).or. &
               (treceptor(n).ge.lrecoutend)) cycle  ! skip if not in current sampling time interval

          ! find indice to creceptor, xkreceptor, and nnreceptor
          do k=1,maxrecsample
            if (cpointer(k).eq.n) exit
          end do

          !! test
          write(*,*) 'receptoroutput: lrecoutstart, lrecoutend, treceptor(n) = ',lrecoutstart, lrecoutend, treceptor(n)
          write(*,*) 'receptoroutput: n, cpointer(k), rpointer = ',n, cpointer(k), rpointer

          if (.not.llcmoutput) then

            ! Compute air density
            !*********************

            ix=int(xreceptor(n))
            jy=int(yreceptor(n))
            ixp=ix+1
            jyp=jy+1
            ddx=xreceptor(n)-float(ix)
            ddy=yreceptor(n)-float(jy)
            rddx=1.-ddx
            rddy=1.-ddy
            p1=rddx*rddy
            p2=ddx*rddy
            p3=rddx*ddy
            p4=ddx*ddy

            indz=nzmax-1
            indzp=nzmax
            do il=2,nzmax
              if (height(il).gt.zreceptor(n)) then
                indz=il-1
                indzp=il
                exit
              endif
            end do

            dz1=zreceptor(n)-height(indz)
            dz2=height(indzp)-zreceptor(n)
            dz=1./(dz1+dz2)

            ! Take density from 2nd wind field in memory
            !********************************************

            do ind=indz,indzp
              ! assume moist air density
              rho_p(ind-indz+1)=p1*rho(ix ,jy ,ind,2) + &
                                p2*rho(ixp,jy ,ind,2) + &
                                p3*rho(ix ,jyp,ind,2) + &
                                p4*rho(ixp,jyp,ind,2)
              ! dry air density             
!              rho_p(ind-indz+1)= &
!                 p1*rho(ix ,jy ,ind,2)*(1. - qv(ix ,jy ,ind,2)) + &
!                 p2*rho(ixp,jy ,ind,2)*(1. - qv(ixp,jy ,ind,2)) + &
!                 p3*rho(ix ,jyp,ind,2)*(1. - qv(ix ,jyp,ind,2)) + &
!                 p4*rho(ixp,jyp,ind,2)*(1. - qv(ixp,jyp,ind,2))
            end do
            rhoi=(dz1*rho_p(2)+dz2*rho_p(1))*dz
            densityoutrecept(1)=rhoi

          else 

            ! no normalization by air density
            densityoutrecept(1)=1.

          endif ! llcmoutput

          ! Write receptor output
          !**********************

          if (recoutnum(k).gt.0.) then
            do ks = ks_start,nspec
              ! write mass concentration
              if ((iout.eq.1).or.(iout.eq.3).or.(iout.eq.5)) then
                ! concentration (ng/m3)
                crec(ks,1)=creceptor(ks,k)*1.E12/recoutnum(k)
                cunc(ks,1)=crecuncert(ks,k)*1.E12/recoutnum(k)
              else if (iout.eq.2) then
                ! mixing ratio (ppt)
                ! note: for iout=3 (both conc and mixing ratio) only output conc at receptors
                crec(ks,1)=creceptor(ks,k)*1.E12*weightair/weightmolar(ks)/ &
                              densityoutrecept(1)/recoutnum(k)
                cunc(ks,1)=crecuncert(ks,k)*1.E12*weightair/weightmolar(ks)/ &
                              densityoutrecept(1)/recoutnum(k)
              endif
            end do
            nnrec(1)=nnreceptor(k)/recoutnum(k)
            xkrec(1)=xkreceptor(k)/recoutnum(k)
          else
            ! no particles found in kernel for this receptor
            crec(:,1)=-999.
            cunc(:,1)=-999.
            nnrec(1)=-999.
            xkrec(1)=-999.
          endif

!$OMP CRITICAL
          ! write receptor output this time interval 
          if (lnetcdfout.eq.1) then
#if USE_NCF
            call write_receptor_netcdf(crec,cunc,nnrec,xkrec,n)
#endif
          else
            call write_receptor_binary(crec,cunc,nnrec,xkrec,n)
          endif

          ! advance output time index
          rpointer=rpointer+1
!$OMP END CRITICAL

        end do ! numreceptor
!$OMP END DO
!$OMP END PARALLEL

        !! test
        write(*,*) 'receptoroutput: rpointer = ',rpointer

        ! Loop over satellites
        !*********************

!$OMP PARALLEL &
!$OMP PRIVATE(n,k,ix,jy,ixp,jyp,ddx,ddy,rddx,rddy,p1,p2,p3,p4,kz,numsatlayer,zmid,indz,indzp, &
!$OMP         il,dz1,dz2,dz,ind,rho_p,rhoi,densityoutrecept,ks,crec,cunc,nnrec,xkrec) &
!$OMP SHARED(spointer) 

!$OMP DO
        do n=1,numsatreceptor

          if ((tsatellite(n).lt.lrecoutstart).or. &
               (tsatellite(n).ge.lrecoutend)) cycle  ! skip if not in current sampling time interval

          ! find indice to csatellite, xksatellite, and nnsatellite
          do k=1,maxrecsample
            if (csatpointer(k).eq.n) exit
          end do
          !! test
          write(*,*) 'receptoroutput: n, csatpointer(k) = ',n, csatpointer(k)

          ! get actual number vertical layers for this retrieval
          do nn=1,numsatellite
            nchar=len_trim(sat_name(nn))
            if (satellitename(n)(1:nchar).eq.trim(sat_name(nn))) then
              numsatlayer=nnsatlayer(nn)
              exit
            endif
          end do

          if (.not.llcmoutput) then

            ! Compute air density
            !*********************

            ix=int(xsatellite(n))
            jy=int(ysatellite(n))
            ixp=ix+1
            jyp=jy+1
            ddx=xsatellite(n)-float(ix)
            ddy=ysatellite(n)-float(jy)
            rddx=1.-ddx
            rddy=1.-ddy
            p1=rddx*rddy
            p2=ddx*rddy
            p3=rddx*ddy
            p4=ddx*ddy

            do kz=1,numsatlayer
  
              zmid=0.5*(zsatellite(kz,n)+zsatellite(kz+1,n))
              indz=nzmax-1
              indzp=nzmax
              do il=2,nzmax
                if (height(il).gt.zmid) then
                  indz=il-1
                  indzp=il
                  exit
                endif
              end do

              dz1=zmid-height(indz)
              dz2=height(indzp)-zmid
              dz=1./(dz1+dz2)

              ! Take density from 2nd wind field in memory
              !********************************************

              do ind=indz,indzp
                ! assume moist air density
                rho_p(ind-indz+1)=p1*rho(ix ,jy ,ind,2) + &
                                  p2*rho(ixp,jy ,ind,2) + &
                                  p3*rho(ix ,jyp,ind,2) + &
                                  p4*rho(ixp,jyp,ind,2)
                ! dry air density             
!                rho_p(ind-indz+1)= &
!                   p1*rho(ix ,jy ,ind,2)*(1. - qv(ix ,jy ,ind,2)) + &
!                   p2*rho(ixp,jy ,ind,2)*(1. - qv(ixp,jy ,ind,2)) + &
!                   p3*rho(ix ,jyp,ind,2)*(1. - qv(ix ,jyp,ind,2)) + &
!                   p4*rho(ixp,jyp,ind,2)*(1. - qv(ixp,jyp,ind,2))
              end do
              rhoi=(dz1*rho_p(2)+dz2*rho_p(1))*dz
              densityoutrecept(kz)=rhoi
 
            end do ! numsatlayer  

          else

            ! no normalization by density
            densityoutrecept(:)=1.

          endif ! llcmoutput

          ! Write receptor output
          !**********************

          do kz=1,numsatlayer
            if (recoutnumsat(kz,k).gt.0.) then
              do ks = ks_start,nspec
                ! write mass concentration
                if ((iout.eq.1).or.(iout.eq.3).or.(iout.eq.5)) then
                  ! concentration (ng/m3)
                  crec(ks,kz)=csatellite(ks,kz,k)*1.E12/recoutnumsat(kz,k)
                  cunc(ks,kz)=csatuncert(ks,kz,k)*1.E12/recoutnumsat(kz,k)
                else if (iout.eq.2) then
                  ! mixing ratio (ppt)
                  ! note: for iout=3 (both conc and mixing ratio) only output conc at receptors
                  crec(ks,kz)=csatellite(ks,kz,k)*1.E12*weightair/weightmolar(ks)/ &
                                densityoutrecept(kz)/recoutnumsat(kz,k)
                  cunc(ks,kz)=csatuncert(ks,kz,k)*1.E12*weightair/weightmolar(ks)/ &
                                densityoutrecept(kz)/recoutnumsat(kz,k)
                endif
              end do
              nnrec(kz)=nnsatellite(kz,k)/recoutnumsat(kz,k)
              xkrec(kz)=xksatellite(kz,k)/recoutnumsat(kz,k)
            else
              crec(:,kz)=-999.
              cunc(:,kz)=-999.
              nnrec(kz)=-999.
              xkrec(kz)=-999.
            endif
          end do

!$OMP CRITICAL
          ! write satellite output this time interval
          if (lnetcdfout.eq.1) then
#if USE_NCF
            call write_satellite_netcdf(crec,cunc,nnrec,xkrec,n)
#endif  
          else
            call write_satellite_binary(crec,cunc,nnrec,xkrec,n)
          endif

          ! advance output time index
          spointer=spointer+1
!$OMP END CRITICAL

        end do ! numsatreceptor
!$OMP END DO
!$OMP END PARALLEL

        ! close files
        if (lnetcdfout.eq.1) then
#ifdef USE_NCF
          if (numreceptor.gt.0) then
            call close_receptor_netcdf
          endif
          if (numsatreceptor.gt.0) then
            call close_satellite_netcdf
          endif
#endif
        endif

        ! Reinitialization
        !*****************

        deallocate(densityoutrecept,nnrec,xkrec)
        deallocate(crec,cunc)

        creceptor(:,:)=0.
        crecuncert(:,:)=0.
        nnreceptor(:)=0.
        xkreceptor(:)=0.
        if (numsatreceptor.gt.0) then
          csatellite(:,:,:)=0.
          csatuncert(:,:,:)=0.
          nnsatellite(:,:)=0.
          xksatellite(:,:)=0.
        endif

        ! End of receptor output
        !***********************

        recoutnum(:)=0.
        recoutnumsat(:,:)=0.

      endif ! output receptor conc

      ! Update output timesteps
      !************************
      lrecoutnext=lrecoutnext+lrecoutstep
      lrecoutstart=lrecoutnext-lrecoutaver/2
      lrecoutend=lrecoutnext+lrecoutaver/2

      if (itime.eq.lrecoutstart) then
        weight=0.5
        call receptorcalc(itime,weight,lrecoutstart,lrecoutend,recoutnum,recoutnumsat)
      endif

    endif ! calculate receptors

  end subroutine receptoroutput


  subroutine receptorcalc(itime,weight,lrecoutstart,lrecoutend,recoutnum,recoutnumsat)

  !*****************************************************************************
  !                                                                            *
  !     Calculation of the concentrations at receptor points using the         *
  !     kernel method                                                          *
  !                                                                            *
  !     Author: A. Stohl                                                       *
  !                                                                            *
  !     Modifications:                                                         * 
  !       Sep-2023, R. Thompson: added option for domain filling mode          *
  !                                                                            *
  !*****************************************************************************

    implicit none

    integer :: itime, itage, lrecoutstart, lrecoutend
    real, dimension(maxrecsample) :: recoutnum
    real, dimension(nlayermax,maxrecsample) :: recoutnumsat
    real    :: weight
    real    :: xd, yd, zd, hx, hy, hz, h, r2, xkern
    real, dimension(nspec) :: conc 
    real, dimension(:,:), allocatable :: unc 
    integer :: n, j, k, ks, kz, i, nn, ks_start, nchar, numsatlayer
    real    :: hxmax, hymax, hzmax, hxsat, hysat, rec_ff, xksum, eta, zmid
    real, parameter :: factor=.596831  ! 15/(8*pi)
    real, parameter :: zref=2000.      ! normalizing height for calculating eta
    integer, parameter :: jmax=4000    ! max number of particles in kernel

    ! initialization
    !***************

    if (llcmoutput) then
      ks_start=2
    else
      ks_start=1
    endif
    allocate(unc(nspec,jmax))

    ! hxmax and hymax in degrees, hzmax in metres
    !********************************************

    if (mdomainfill.ne.1) then
      hxmax=6.0
      hymax=4.0
      hzmax=150.
    else
      hxmax=2.0
      hymax=1.5
      hzmax=300.
      hxsat=1.75
      hysat=1.25
    endif

    ! convert h-values to grid coordinates
    !*************************************

    hxmax=hxmax/dx
    hymax=hymax/dy
    hxsat=hxsat/dx
    hysat=hysat/dy

    write(*,*) 'hxmax, hymax = ',hxmax, hymax
    write(*,*) 'hxsat, hysat = ',hxsat, hysat

    ! Loop over receptors
    !********************

    ! pointer for creceptor, xkreceptor and nnreceptor
    cpointer(:)=0
    k=0

    do n=1,numreceptor

      if ((treceptor(n).lt.lrecoutstart).or. &
           (treceptor(n).gt.lrecoutend)) cycle  ! skip if not in current sampling time interval

      !! test
      write(*,*) 'receptorcalc: lrecoutstart, lrecoutend, treceptor(n) = ',lrecoutstart, lrecoutend, treceptor(n)

      ! update pointer
      k=k+1
      if (k.gt.maxrecsample) then
        write(*,*) 'FLEXPART ERROR in receptorcalc: maxrecsample too small'
        error stop
      endif
      cpointer(k)=n

      !! test
      write(*,*) 'receptorcalc: n, receptorname(n) = ',n, receptorname(n)

      ! Reset concentrations for new receptor
      conc(:)=0.
      unc(:,:)=0.
      xksum=0.
      j=0

      if (zreceptor(n).lt.hzmax) then
        rec_ff=0.5 + 0.9375*(zreceptor(n)/hzmax) - &
                 0.625*(zreceptor(n)/hzmax)**3 + &
                 0.1875*(zreceptor(n)/hzmax)**5
        rec_ff=1./rec_ff
      else
        rec_ff=1.
      endif
      h=hxmax*hymax*hzmax

      !! test
      write(*,*) 'receptorcalc: rec_ff = ',rec_ff

      ! Estimate concentration at receptor
      !***********************************

!$OMP PARALLEL &
!$OMP PRIVATE(i,itage,hz,zd,hx,xd,hy,yd,h,r2,xkern,ks) &
!$OMP SHARED(j,unc,conc) &
!$OMP REDUCTION(+:xksum,xkreceptor,nnreceptor)

!$OMP DO
      do i=1,count%alive

        if (mdomainfill.ne.1) then

          ! not domain filling run so consider age of particle
          itage=abs(itime-part(i)%tstart)

          hz=min(50.+0.3*sqrt(real(itage)),hzmax)
          zd=part(i)%z/hz
          if (zd.gt.1.) cycle           ! save computing time

          hx=min((0.29+2.222e-3*sqrt(real(itage)))*dx+ &
               real(itage)*1.2e-5,hxmax)                     ! 80 km/day
          xd=(part(i)%xlon-xreceptor(n))/hx
          if (xd*xd.gt.1.) cycle        ! save computing time

          hy=min((0.18+1.389e-3*sqrt(real(itage)))*dy+ &
               real(itage)*7.5e-6,hymax)                     ! 80 km/day
          yd=(part(i)%ylat-yreceptor(n))/hy
          if (yd*yd.gt.1.) cycle        ! save computing time
          h=hx*hy*hz

          r2=xd*xd+yd*yd+zd*zd
          if (r2.ge.1.) cycle           ! save computing time
          xkern=2.*factor*(1.-r2)       ! parabolic kernel

        else

          ! domain filling run 
          xd=(part(i)%xlon-xreceptor(n))/hxmax
          if (xd*xd.gt.1) cycle         ! save computing time

          yd=(part(i)%ylat-yreceptor(n))/hymax
          if (yd*yd.gt.1.) cycle        ! save computing time

          zd=(part(i)%z-zreceptor(n))/hzmax
          if (zd*zd.gt.1.) cycle        ! save computing time

          r2=xd*xd+yd*yd+zd*zd
          if (r2.ge.1.) cycle           ! save computing time
          xkern=rec_ff*factor*(1.-r2)   ! parabolic kernel

        endif ! mdomainfill

        ! counter of particles used in conc calculation
        ! note: use of atomic may not be efficient, alternative split unc, j over threads
        ! and combine all threads after end of parallel section, need to save j from each thread
!$OMP ATOMIC
        j=j+1
        if (j.gt.jmax) then
          write(*,*) 'FLEXPART ERROR in receptorcalc: size of jmax too small'
          stop
        endif

        do ks=ks_start,nspec
          if (llcmoutput) then
            ! special case LCM output use mass ratio species to airtracer
            ! species 1 is always airtracer
            conc(ks)=conc(ks) + mass(i,ks)/mass(i,1) * &
                      weight * xkern
            unc(ks,j)=mass(i,ks)/mass(i,1)
          else
            ! normal case
            conc(ks)=conc(ks) + mass(i,ks) * &
                       weight * xkern/h/receptorarea(n)
            unc(ks,j)=mass(i,ks)/h/receptorarea(n)
          endif
        end do
        nnreceptor(k)=nnreceptor(k) + 1.
        xkreceptor(k)=xkreceptor(k) + xkern
        xksum=xksum + xkern

      end do ! count%alive
!$OMP END DO
!$OMP END PARALLEL

      do ks=ks_start,nspec
        if (conc(ks).gt.0.) then
          ! only do if conc could be calculated for this receptor at this time
          creceptor(ks,k)=creceptor(ks,k) + conc(ks)/xksum
          crecuncert(ks,k)=crecuncert(ks,k) + &
                               sqrt(sum((unc(ks,1:j)-conc(ks)/xksum/weight)**2)/real(j))
        endif
      end do
      if (any(conc(:).gt.0.)) then
        recoutnum(k)=recoutnum(k)+weight
      endif

      !! test
      write(*,*) 'receptorcalc: j, conc, xksum = ',j, conc(2), xksum
      write(*,*) 'receptorcalc: n, k, cpointer(k) = ',n, k, cpointer(k)
      write(*,*) 'receptorcalc: nnreceptor(k) = ',nnreceptor(k)
      write(*,*) 'receptorcalc: xkreceptor(k) = ',xkreceptor(k)
      write(*,*) 'receptorcalc: creceptor(2,k) = ',creceptor(2,k)
      write(*,*) 'receptorcalc: crecuncert(2,k) = ',crecuncert(2,k)

    end do ! numreceptor

    ! Loop over satellites
    !*********************

    ! pointer for csatellite, xksatellite and nnsatellite
    csatpointer(:)=0
    k=0

    do n=1,numsatreceptor

      if ((tsatellite(n).lt.lrecoutstart).or. &
           (tsatellite(n).gt.lrecoutend)) cycle  ! skip if not in current sampling time interval
   
       !! test
      write(*,*) 'receptorcalc: lrecoutstart, lrecoutend, tsatellite(n) = ',lrecoutstart, lrecoutend, tsatellite(n)

      ! update pointer
      k=k+1
      if (k.gt.maxrecsample) then
        write(*,*) 'FLEXPART ERROR in receptorcalc: maxrecsample too small'
        stop
      endif
      csatpointer(k)=n

      ! get actual number vertical layers for this retrieval
      do nn=1,numsatellite
        nchar=len_trim(sat_name(nn))
        if (satellitename(n)(1:nchar).eq.trim(sat_name(nn))) then
          numsatlayer=nnsatlayer(nn)
          exit
        endif
      end do

      !! test
      write(*,*) 'receptorcalc: n, satellitename(n) = ',n, satellitename(n)

      do kz=1,numsatlayer

        ! Reset concentrations for new receptor
        conc(:)=0.
        unc(:,:)=0.
        xksum=0.
        j=0

        ! height of layer
        hz=zsatellite(kz+1,n)-zsatellite(kz,n)

        ! rare cases when 2 satellite pressure layers fall in same meteo layer
        ! set minimum height between layers 
        hz=max(200.,hz)

        ! midpoint of layer
        zmid=zsatellite(kz,n)+0.5*hz

        ! factor by which to expand horizontal threshhold
        ! as expect particle density to decrease with altitude
!        eta=max(1.,sqrt(zsatellite(kz,n)/zref))
        eta=max(1.,sqrt(0.5*zmid/zref))

        h=hxsat*hysat*(0.5*hz)

        !! test
        write(*,*) 'zmid, hz, eta =',zmid, hz, eta

!$OMP PARALLEL &
!$OMP PRIVATE(i,zd,xd,yd,r2,xkern,ks) &
!$OMP SHARED(j,unc,conc) &
!$OMP REDUCTION(+:xksum,xksatellite,nnsatellite)

!$OMP DO
        do i=1,count%alive

          ! sample satellite retrievals for domain-filling and 
          ! non-domain-filling modes in the same way

          zd=(part(i)%z-zmid)/(0.5*hz)
          if (zd*zd.gt.1.) cycle        ! save computing time

          xd=(part(i)%xlon-xsatellite(n))/(eta*hxsat)
          if (xd*xd.gt.1) cycle         ! save computing time

          yd=(part(i)%ylat-ysatellite(n))/(eta*hysat)
          if (yd*yd.gt.1.) cycle        ! save computing time

          r2=xd*xd+yd*yd+zd*zd
          if (r2.ge.1.) cycle           ! save computing time

          xkern=factor*(1.-r2)          ! parabolic kernel

          ! counter of particles used in conc calculation
!$OMP ATOMIC
          j=j+1
          if (j.gt.jmax) then
            write(*,*) 'FLEXPART ERROR in receptorcalc: size of jmax too small'
            stop
          endif

          do ks=ks_start,nspec
            if (llcmoutput) then
              ! special case LCM output use mass ratio species to airtracer
              ! species 1 is always airtracer
              conc(ks)=conc(ks) + mass(i,ks)/mass(i,1) * &
                        weight * xkern
              unc(ks,j)=mass(i,ks)/mass(i,1)
            else
              ! normal case
              conc(ks)=conc(ks) + mass(i,ks) * &
                         weight * xkern/h/satellitearea(n)
              unc(ks,j)=mass(i,ks)/h/satellitearea(n)
            endif
          end do
          nnsatellite(kz,k)=nnsatellite(kz,k) + 1.
          xksatellite(kz,k)=xksatellite(kz,k) + xkern
          xksum=xksum + xkern

        end do ! count%alive
!$OMP END DO
!$OMP END PARALLEL

        do ks=ks_start,nspec
          if (conc(ks).gt.0.) then
            ! only do if conc could be calculated for this receptor and time
            csatellite(ks,kz,k)=csatellite(ks,kz,k) + conc(ks)/xksum
            csatuncert(ks,kz,k)=csatuncert(ks,kz,k) + &
                                 sqrt(sum((unc(ks,1:j)-conc(ks)/xksum/weight)**2)/real(j))
          endif
        end do
        if (any(conc(:).gt.0.)) then
          recoutnumsat(kz,k)=recoutnumsat(kz,k)+weight
        endif

        !! test
        write(*,*) 'receptorcalc: for satellite: j, conc, xksum = ',j, conc(2), xksum
        write(*,*) 'receptorcalc: n, k, csatpointer(k) = ',n, k, csatpointer(k)
        write(*,*) 'receptorcalc: nnsatellite(:,k) = ',nnsatellite(:,k)
        write(*,*) 'receptorcalc: xksatellite(:,k) = ',xksatellite(:,k)

      end do ! numsatlayer

    end do ! numsatreceptor

  end subroutine receptorcalc


end module receptor_mod

