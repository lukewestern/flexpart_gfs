! SPDX-FileCopyrightText: FLEXPART 1998-2019, see flexpart_license.txt
! SPDX-License-Identifier: GPL-3.0-or-later

  !*****************************************************************************
  !                                                                            *
  !    This module contains variables and subroutines for  calculating         *
  !    chemical loss of species                                                *
  !                                                                            *
  !*****************************************************************************

module chemistry_mod

  use netcdf
  use par_mod
  use com_mod
  use date_mod
  use totals_mod, only: chem_loss
  use netcdf_output_mod, only: nf90_err

  ! reagent field variables 

  implicit none

  integer                                 :: nxCL,nyCL,nzCL
  integer                                 :: nxjr, nyjr
  real, allocatable, dimension(:)         :: lonCL,latCL,altCL
  real                                    :: dxCL,dyCL,dxjr,dyjr
  real, allocatable, dimension(:,:,:,:,:) :: CL_field           ! chemical loss fields at input resolution
  real, allocatable, dimension(:,:,:,:)   :: clfield_cur        ! chemical loss fields current hour
  integer, dimension(2)                   :: memCLtime          ! time of fields in memory (sec)  
  integer                                 :: curCLhour          ! current hour since start of simulation
  real(kind=dp), dimension(2)             :: CL_time            ! julian date of fields in memory
  real, allocatable, dimension(:,:,:)     :: jrate_average      ! monthly average actinic flux 
  real, allocatable, dimension(:)         :: lonjr,latjr

  private :: zenithangle, photo_O1D

  contains

  subroutine readreagents()

  !*****************************************************************************
  !                                                                            *
  !   This routine reads names and input paths of chemical reagents            *
  !                                                                            *
  !   Author: R. Thompson, Sep 2023                                            *
  !                                                                            *
  !*****************************************************************************
  !                                                                            *
  ! Variables:                                                                 *
  ! preagent              Reagent name                                         *
  ! preag_path            Path to reagent fields                               *
  ! phourly               Logical for hourly interpolate (0 = no, 1 = yes)     *
  !                                                                            *
  !*****************************************************************************

    implicit none

    character(len=256) :: preag_path
    character(len=10)  :: preagent
    integer :: phourly
    integer,parameter :: unitreagents=2, unitreagentsout=3
    integer :: readerror
    integer :: j

    ! declare namelist
    namelist /reagent_params/ &
         preagent, preag_path, phourly

    preagent="" ! read failure indicator value
    preag_path=""
    phourly=0
    reag_hourly(:)=0

    open(unitreagents,file=path(1)(1:length(1))//'REAGENTS',status='old',form='formatted',iostat=readerror)
    if (readerror.ne.0) then
      ! no REAGENTS file
      nreagent=0
      go to 999
    endif

    ! prepare namelist output if requested
    if (nmlout.and.lroot) then
      open(unitreagentsout,file=path(2)(1:length(2))//'REAGENTS.namelist',access='append',status='replace',iostat=readerror)
      if (readerror.ne.0) then
        write(*,*) '#### FLEXPART MODEL ERROR CANNOT CREATE FILE:   ####'
        write(*,*) '#### ',path(2)(1:length(2))//'REAGENTS.namelist ####'
      endif
    endif

    ! try namelist input
    read(unitreagents,reagent_params,iostat=readerror)
    close(unitreagents)
    if (readerror.ne.0) then
      ! not namelist format
      nreagent=0
      go to 999
    endif

    ! namelist input
    open(unitreagents,file=path(1)(1:length(1))//'REAGENTS',status='old',form='formatted',iostat=readerror)
    j=0
    do while (readerror.eq.0)
      j=j+1
      if (j.gt.maxreagent) then
        write(*,*) ' #### FLEXPART MODEL ERROR! TOO MANY REAGENTS #### '
        write(*,*) ' #### MAXIMUM NUMBER IS ',maxreagent,'        #### '
        write(*,*) ' #### PLEASE MAKE CHANGES IN FILE REAGENTS    #### '
        stop
      endif
      read(unitreagents,reagent_params,iostat=readerror)
      reagents(j)=preagent
      reag_path(j)=preag_path
      reag_hourly(j)=phourly
      ! namelist output
      if (nmlout.and.lroot) then
        write(unitreagentsout,nml=reagent_params)
      endif
    end do
    nreagent=j-1
    if (lroot) then
      write(*,*) 'Number of reagents: ',nreagent
      write(*,*) 'Reagent names: ',reagents(1:nreagent)
    endif
    close(unitreagents)
    close(unitreagentsout)

999 continue ! no reagents file

  end subroutine readreagents

  subroutine getchemfield(itime)

  !*****************************************************************************
  !                                                                            *
  !    This routine reads the chemical reagent fields into memory              *
  !                                                                            *
  !    Author: Rona Thompson, Sep 2023                                         *
  !                                                                            *
  !*****************************************************************************
  !                                                                            *
  ! Variables:                                                                 *
  !                                                                            *
  ! itime            time since start of simulation in sec                     *
  !                                                                            *
  !*****************************************************************************

    implicit none

    integer          :: itime
    integer          :: jjjjmmdd, hhmmss, mm, eomday
    integer          :: nr, memid
    character(len=2) :: amonth
    character(len=256):: CL_name, jr_name
    logical          :: lexist

    do nr=1,nreagent

      if(lroot) write(*,*) 'Getting fields for reagent: ',reagents(nr)
      print*, 'ldirect*memCLtime(1) = ',ldirect*memCLtime(1)
      print*, 'ldirect*memCLtime(2) = ',ldirect*memCLtime(2)

      ! Check fields are available for the current time step
      !*****************************************************

      if ((ldirect*memCLtime(1).le.ldirect*itime).and. &
           (ldirect*memCLtime(2).gt.ldirect*itime)) then

        ! The rightfields are already in memory -> don't do anything
        exit

      else if ((ldirect*memCLtime(2).le.ldirect*itime).and. &
           (memCLtime(2).ne.0.)) then

        ! Current time is after 2nd chem field
        !*************************************

        write(*,*) 'Updating CL field for: ',trim(reagents(nr))

        CL_field(:,:,:,nr,1)=CL_field(:,:,:,nr,2)
        memCLtime(1)=memCLtime(2)  ! time in sec
        CL_time(1)=CL_time(2)      ! julian date
        memid=2

        ! Read new chem field and store in 2nd position
        !**********************************************

        ! julian date of next chem field assuming monthly fields
        call caldate(CL_time(1), jjjjmmdd, hhmmss)
        eomday=calceomday(jjjjmmdd/100)
        memCLtime(2)=memCLtime(1)+ldirect*eomday*24*3600   ! time in sec
        CL_time(2)=CL_time(1)+real(ldirect*eomday,kind=dp) ! julian date 
        !! test
        write(*,*) 'memCLtime = ',memCLtime(1), memCLtime(2)
        write(*,*) 'CL_time = ',CL_time(1), CL_time(2)
        !!
        call caldate(CL_time(2), jjjjmmdd,hhmmss)
        mm=int((jjjjmmdd-(jjjjmmdd/10000)*10000)/100)
        write(amonth,'(I2.2)') mm
        write(CL_name,'(A)') trim(reag_path(nr))//trim(reagents(nr))//'_'//amonth//'.nc'
        inquire(file=CL_name,exist=lexist)
        if (lexist) then
          write(*,*) 'Reading CL field: ',CL_name
          call readchemfield(CL_name, memid, nr)
        else
          write(*,*) '#### FLEXPART ERROR                ####'
          write(*,*) '#### CHEMISTRY FIELD NOT FOUND     ####'
          write(*,*) '#### '//trim(CL_name)//' ####'
          error stop
        endif

        if (reag_hourly(nr).gt.0) then

          ! Read average jrates and store in 2nd position
          !**********************************************

          write(jr_name,'(A)') trim(reag_path(nr))//'jrate_average.nc'
          inquire(file=jr_name,exist=lexist)
          if (lexist) then
            write(*,*) 'Reading jrate field: ',jr_name
            call readjrate(jr_name, memid, mm)
          else
            write(*,*) '#### FLEXPART ERROR                ####'
            write(*,*) '#### JRATE_AVERAGE NOT FOUND       ####'
            write(*,*) '#### '//trim(jr_name)//' ####'
            error stop
          endif

        endif

      else

        ! No chem field in memory that can be used
        !******************************************  

        write(*,*) 'Reading two new CL fields for: ',trim(reagents(nr))

        ! read both fields into memory
        do memid=1,2
          if (memid.eq.1) then
            CL_time(memid)=bdate+real(ldirect*itime,kind=dp)/86400._dp
          else
            CL_time(memid)=CL_time(memid-1)+real(ldirect*eomday,kind=dp)
          endif
          memCLtime(memid)=int((CL_time(memid)-bdate)*86400._dp)*ldirect
          call caldate(CL_time(memid), jjjjmmdd, hhmmss)
          mm=int((jjjjmmdd-(jjjjmmdd/10000)*10000)/100)
          eomday=calceomday(jjjjmmdd/100)
          write(amonth,'(I2.2)') mm
          write(CL_name,'(A)') trim(reag_path(nr))//trim(reagents(nr))//'_'//amonth//'.nc'
          inquire(file=CL_name,exist=lexist)
          if (lexist) then
            write(*,*) 'Reading CL field: ',CL_name
            call readchemfield(CL_name, memid, nr)
          else
            write(*,*) '#### FLEXPART ERROR                ####'
            write(*,*) '#### CHEMISTRY FIELD NOT FOUND     ####'
            write(*,*) '#### '//trim(CL_name)//' ####'
            error stop
          endif
          if (reag_hourly(nr).gt.0) then
            ! Read average jrates into memory 
            write(jr_name,'(A)') trim(reag_path(nr))//'jrate_average.nc'
            inquire(file=jr_name,exist=lexist)
            if (lexist) then
              write(*,*) 'Reading jrate field: ',jr_name
              call readjrate(jr_name, memid, mm)
            else
              write(*,*) '#### FLEXPART ERROR                ####'
              write(*,*) '#### JRATE_AVERAGE NOT FOUND       ####'
              write(*,*) '#### '//trim(jr_name)//' ####'
              error stop
            endif
          endif ! reag_hourly
        end do ! memid

      endif ! update hourly fields

    end do ! nreagent

  end subroutine getchemfield

  subroutine readchemfield(CL_name,memid,nr)

  !*****************************************************************************
  !                                                                            *
  !    Reads the chemical reagent fields into memory                           *
  !                                                                            *
  !    Author: Rona Thompson, Sep 2023                                         *
  !                                                                            *
  !*****************************************************************************
  !                                                                            *
  ! Variables:                                                                 *
  !                                                                            *
  ! CL_name           name of chemical reagent file                            *
  ! memid             time index to chemical field variable                    *
  ! nr                reagent index to chemical field variable                 *
  !                                                                            * 
  !*****************************************************************************

    implicit none

    character(len=256) :: CL_name
    integer            :: memid,nr
    integer            :: ncid,dimid,varid

    ! Read netcdf file
    !******************

    ! open file
    call nf90_err( nf90_open(trim(CL_name),nf90_NOWRITE,ncid) )

    ! longitude
    call nf90_err( nf90_inq_dimid(ncid,'lon',dimid) )
    call nf90_err( nf90_inquire_dimension(ncid,dimid,len=nxCL) )
    if (.not.allocated(lonCL)) allocate(lonCL(nxCL))
    call nf90_err( nf90_inq_varid(ncid,'lon',varid) )
    call nf90_err( nf90_get_var(ncid,varid,lonCL) )
    dxCL=abs(lonCL(2)-lonCL(1))

    ! latitude
    call nf90_err( nf90_inq_dimid(ncid,'lat',dimid) )
    call nf90_err( nf90_inquire_dimension(ncid,dimid,len=nyCL) )
    if (.not.allocated(latCL)) allocate(latCL(nyCL))
    call nf90_err( nf90_inq_varid(ncid,'lat',varid) )
    call nf90_err( nf90_get_var(ncid,varid,latCL) )
    dyCL=abs(latCL(2)-latCL(1))

    ! altitude
    call nf90_err( nf90_inq_dimid(ncid,'lev',dimid) )
    call nf90_err( nf90_inquire_dimension(ncid,dimid,len=nzCL) )
    if (.not.allocated(altCL)) allocate(altCL(nzCL))
    call nf90_err( nf90_inq_varid(ncid,'lev',varid) )
    call nf90_err( nf90_get_var(ncid,varid,altCL) )

    ! chemical field
    call nf90_err( nf90_inq_varid(ncid,trim(reagents(nr)),varid) )
    if (.not.allocated(CL_field)) then
      allocate(CL_field(nxCL,nyCL,nzCL,nreagent,2))
      allocate(clfield_cur(nxCL,nyCL,nzCL,nreagent))
    endif
    call nf90_err( nf90_get_var(ncid,varid,CL_field(:,:,:,nr,memid)) )

    ! close file
    call nf90_err( nf90_close(ncid) )

    return

  end subroutine readchemfield

  subroutine getchemhourly(itime)

  !*****************************************************************************
  !                                                                            *
  !    This routine interpolates the chemistry fields to current hour and      *
  !    if required using information on solar zenith angle                     *
  !                                                                            *
  !    Author: Rona Thompson, Mar 2024                                         *
  !                                                                            *
  !*****************************************************************************
  !                                                                            *
  ! Variables:                                                                 *
  ! itime [s]            actual simulation time [s]                            *
  !                                                                            *
  !*****************************************************************************

    implicit none

    integer          :: itime, curhour, interp_time
    integer          :: nr, kz, ix, jy
    real             :: dt1, dt2, dtt, sza, jrate, jrate_cur
    integer          :: jjjjmmdd, hhmmss
    integer          :: jrx, jry
    real(kind=dp)    :: jul
    !! testing
    character(len=4) :: atime
    character(len=20) :: frmt
    real, dimension(nxjr,nyjr) :: jscalar
    
    jscalar(:,:)=1.

    ! current hour of simulation
    curhour=itime/3600
    print*, 'getchemhourly: curhour, curCLhour = ',curhour, curCLhour

    jul=bdate+real(itime,kind=dp)/86400._dp
    call caldate(jul,jjjjmmdd,hhmmss)

    if ((ldirect*curCLhour.eq.ldirect*curhour).and.(ldirect*itime.gt.0)) then
      ! chemical field is for current hour -> don't do anything
      return
    else
      ! interpolate to middle of hour
      curCLhour=curhour
      interp_time=curhour*3600+1800
      dt1=float(interp_time-memCLtime(1))
      dt2=float(memCLtime(2)-interp_time)
      dtt=1./(dt1+dt2)
      !! testing
      print*, 'getchemhourly: dt1, dt2, dtt = ',dt1,dt2,dtt 
      do nr=1,nreagent
        write(*,*) 'Interpolating fields for reagent: ',reagents(nr)
        clfield_cur(:,:,:,nr)=(dt2*CL_field(:,:,:,nr,1) + dt1*CL_field(:,:,:,nr,2))*dtt
        if (reag_hourly(nr).gt.0) then
          ! use actinic flux (jrate) for interpolation
!$OMP PARALLEL &
!$OMP PRIVATE(kz,jy,jry,ix,jrx,sza,jrate,jrate_cur) 
!$OMP DO
          do kz=1,nzCL
            do jy=1,nyCL  
              ! jrate_average dimensions given as grid midpoints
              jry=int((latCL(jy)-(latjr(1)-0.5*dyjr))/dyjr)+1
              !! testing
              if (kz.eq.1.and.jy.lt.10) print*, 'latCL, latjr = ',latCL(jy), latjr(jry)
              do ix=1,nxCL
                ! jrate_average dimensions given as grid midpoints
                jrx=int((lonCL(ix)-(lonjr(1)-0.5*dxjr))/dxjr)+1
                !! testing
                if (kz.eq.1.and.jy.lt.10.and.ix.lt.10) print*, 'lonCL, lonjr = ',lonCL(ix), lonjr(jrx)
                ! solar zenith angle
                sza=zenithangle(latjr(jry),lonjr(jrx),jjjjmmdd,hhmmss)
                ! calculate J(O1D) (jrate)
                jrate=photo_O1D(sza)
                jrate_cur=(dt2*jrate_average(jrx,jry,1) + &
                           dt1*jrate_average(jrx,jry,2))*dtt
                ! apply hourly correction to chem field
                if(jrate_cur.gt.0.) then
                  clfield_cur(ix,jy,kz,nr)=clfield_cur(ix,jy,kz,nr)*jrate/jrate_cur
                  !! testing
                  if (kz.eq.1) jscalar(ix,jy)=jrate/jrate_cur
                endif
              end do
            end do
          end do
!$OMP END DO
!$OMP END PARALLEL
        endif ! reag_hourly
      end do ! nreagent
    endif ! curhour

    !! testing
    write(frmt,fmt='(A,I4,A)') '(',nxCL,'(E14.6))'
    write(atime,fmt='(I4.4)') curhour
    open(600,file=path(2)(1:length(2))//'clfield_'//atime//'.txt',action='write',status='replace')
    do kz=1,nzCL
      do jy=1,nyCL
        write(600,fmt=frmt) clfield_cur(:,jy,kz,1)
      end do
    end do
    close(600)
    write(frmt,fmt='(A,I4,A)') '(',nxjr,'(E14.6))'
    open(600,file=path(2)(1:length(2))//'jscalar_'//atime//'.txt',action='write',status='replace')
    do jy=1,nyjr
      write(600,fmt=frmt) jscalar(:,jy)
    end do
    close(600)
    if (itime.eq.0) then
      write(frmt,fmt='(A,I4,A)') '(',nxjr,'(E14.6))'
      open(600,file=path(2)(1:length(2))//'javerage_'//atime//'.txt',action='write',status='replace')
      do jy=1,nyjr
        write(600,fmt=frmt) jrate_average(:,jy,1)
      end do
    endif
    close(600)

  end subroutine getchemhourly

  subroutine chemreaction(itime)

  !*****************************************************************************
  !                                                                            *
  !    This routine computes the chemical loss of each species and updates     *
  !    the particle mass                                                       *
  !                                                                            *
  !    Author: Rona Thompson, Sep 2023                                         *
  !    updated for v11, Jan 2024                                               *
  !                                                                            *
  !*****************************************************************************
  !                                                                            *
  ! Variables:                                                                 *
  ! itime [s]            actual simulation time [s]                            *
  !                                                                            *
  !*****************************************************************************

    use particle_mod,   only: count, part, mass
    use point_mod,      only: xlon0, ylat0, dx, dy
    use windfields_mod, only: xresoln, yresoln, xln, yln, xrn, yrn, height, tt, nz
    use omp_lib

    implicit none

    integer               :: jpart,itime
    integer               :: ii,i,j,ks,ix,jy,kz,jrx,jry,nr
    integer               :: ngrid,interp_time,n,indz
    integer               :: jjjjmmdd,hhmmss,mm,hh,m1,m2
    integer               :: clx,cly,clz,clzm,ithread
    real, dimension(nzCL) :: altCLtop
    real, dimension(2)    :: cl_tmp
    real                  :: xlon,ylat
    real                  :: xtn,ytn
    real                  :: dt1,dt2,dtt,dtt1,dtt2,dttt,dz1,dz2,dzz
    real                  :: restmass,clreacted,cl_cur
    real                  :: jrate,jrate_cur,sza
    real                  :: clrate,temp
    real, parameter       :: smallnum = tiny(0.0) ! smallest number that can be handled
    real(kind=dp)         :: jul
    real                  :: lonjrx,latjry
#ifdef _OPENMP
    real(kind=dp), allocatable, dimension(:,:,:) :: chem_loss_tmp
#endif

    ! use middle of synchronisation time step
    interp_time=nint(itime+0.5*lsynctime)
    dtt1=float(interp_time-memtime(1))
    dtt2=float(memtime(2)-interp_time)
    dttt=1./(dtt1+dtt2)

    ! initialization
    chem_loss(:,:)=0d0

#ifdef _OPENMP
    allocate( chem_loss_tmp(nreagent,nspec,numthreads) )
    chem_loss_tmp(:,:,:) = 0d0
#endif 
    ! Loop over particles
    !*****************************************

!$OMP PARALLEL &
!$OMP PRIVATE(ii,jpart,ngrid,j,xtn,ytn,ix,jy, &
!$OMP   xlon,ylat,clx,cly,clz,clzm,kz,altCLtop,dz1,dz2,dzz,nr,i, &
!$OMP   cl_cur,indz,temp,ks,clrate,restmass,clreacted,ithread)
! !$OMP REDUCTION(+:chem_loss) 

#ifdef _OPENMP
    ithread = OMP_GET_THREAD_NUM()+1 ! Starts with 1
#else
    ithread = 1
#endif

!$OMP DO
    do ii=1,count%alive

      jpart=count%ialive(ii)
      ! Determine which nesting level to be used
      ngrid=0
      do j=numbnests,1,-1 ! check if need +/- dxn below
        if ((part(jpart)%xlon.gt.xln(j)).and.(part(jpart)%xlon.lt.xrn(j)).and. &
             (part(jpart)%ylat.gt.yln(j)).and.(part(jpart)%ylat.lt.yrn(j))) then
          ngrid=j
          exit
        endif
      end do

      ! Determine nested grid coordinates
      if (ngrid.gt.0) then
        xtn=(real(part(jpart)%xlon)-xln(ngrid))*xresoln(ngrid)
        ytn=(real(part(jpart)%ylat)-yln(ngrid))*yresoln(ngrid)
        ix=int(xtn)
        jy=int(ytn)
      else
        ix=int(part(jpart)%xlon)
        jy=int(part(jpart)%ylat)
      endif

      ! Get CL from nearest grid-cell for current hour
      !***********************************************

      ! world coordinates
      xlon=real(part(jpart)%xlon)*dx+xlon0
      if (xlon.gt.180.) then
        xlon=xlon-360.
      endif
      ylat=real(part(jpart)%ylat)*dy+ylat0
      if (ylat.gt.90.) then
        ylat=90.
      endif
      ! get position in the chem field
      ! assumes chem field dimensions given as grid midpoints
      clx=min(nxCL, int((xlon-(lonCL(1)-0.5*dxCL))/dxCL)+1)
      cly=min(nyCL, int((ylat-(latCL(1)-0.5*dyCL))/dyCL)+1)

      ! get the level of the chem field for the particle
      ! z is the z-coord of the trajectory above model orography in metres
      ! altCL is the height of the centre of the level in the chem field above orography
      do kz=2,nzCL
        altCLtop(kz-1)=altCL(kz-1)+0.5*(altCL(kz)-altCL(kz-1))
      end do
      altCLtop(nzCL)=altCL(nzCL)+0.5*(altCL(nzCL)-altCL(nzCL-1))
      clzm=nzCL-1
      do clz=1,nzCL
        if (real(part(jpart)%z).lt.altCLtop(clz)) then
          clzm=clz-1
          exit
        endif
      end do
      clz=min(nzCL,clz)
      if (clzm.eq.0 ) then
        dz1=1.
        dz2=1.
        clzm=clz
      else
        dz1=real(part(jpart)%z)-altCL(clzm)
        dz2=altCL(clz)-real(part(jpart)%z)
      endif
      if (dz1.eq.(-1.*dz2)) then
        dz1=1.
        dz2=1.
      endif
      dzz=1./(dz1+dz2)

      do nr=1,nreagent

        ! chem reagent for particle time and location 
        cl_cur=(dz2*clfield_cur(clx,cly,clzm,nr) + &
                dz1*clfield_cur(clx,cly,clz,nr))*dzz

        ! Compute chemical loss
        !**********************

        if (cl_cur.gt.smallnum) then
          indz=nz
          do kz=2,nz
            if (height(kz).gt.part(jpart)%z) then
              indz=kz-1
              exit
            endif
          end do
          temp=(tt(ix,jy,indz,1)*dtt2 &
                 + tt(ix,jy,indz,2)*dtt1)*dttt
          do ks=1,nspec
            if (reaccconst(nr,ks).gt.0.) then
              ! k = CT^Nexp(-D/temp)[reagent]
              clrate=reaccconst(nr,ks)*(temp**reacnconst(nr,ks))*exp(-1.*reacdconst(nr,ks)/temp)*cl_cur
              ! new particle mass
              restmass=mass(jpart,ks)*exp(-1.*clrate*abs(lsynctime))
              if (restmass.gt.smallnum) then
                clreacted=mass(jpart,ks)-restmass
                mass(jpart,ks)=restmass
              else
                clreacted=mass(jpart,ks)
                mass(jpart,ks)=0.
              endif
#ifdef _OPENMP
              chem_loss_tmp(nr,ks,ithread)=chem_loss_tmp(nr,ks,ithread)+real(clreacted,kind=dp)
#else
              chem_loss(nr,ks)=chem_loss(nr,ks)+real(clreacted,kind=dp)
#endif
            endif
          end do  ! nspec
        endif   

      end do  ! nreagent

    end do  ! loop over all particles
!$OMP END DO

!$OMP END PARALLEL
  
#ifdef _OPENMP
  do ithread=1,numthreads
    chem_loss(:,:) = chem_loss(:,:)+chem_loss_tmp(:,:,ithread)
  end do 
  deallocate( chem_loss_tmp )
#endif

  end subroutine chemreaction


  subroutine readjrate(jr_name,memid,mm)

  !*****************************************************************************
  !                                                                            *
  !    Reads the average actinic flux rates into memory                        *
  !                                                                            *
  !    Author: Rona Thompson, Sep 2023                                         *
  !                                                                            *
  !*****************************************************************************
  !                                                                            *
  ! Variables:                                                                 *
  !                                                                            *
  ! jr_name           name of the file                                         *
  ! memid             time index to chemical field variable                    *
  ! mm                month to read                                            *
  !                                                                            * 
  !*****************************************************************************

    implicit none

    character(len=256) :: jr_name
    integer            :: memid,mm
    integer            :: ncid,dimid,varid

    ! Read netcdf file
    !******************

    ! open file
    call nf90_err( nf90_open(trim(jr_name),nf90_NOWRITE,ncid) )

    ! longitude
    call nf90_err( nf90_inq_dimid(ncid,'longitude',dimid) )
    call nf90_err( nf90_inquire_dimension(ncid,dimid,len=nxjr) )
    if (.not.allocated(lonjr)) allocate(lonjr(nxjr))
    call nf90_err( nf90_inq_varid(ncid,'longitude',varid) )
    call nf90_err( nf90_get_var(ncid,varid,lonjr) )
    dxjr=abs(lonjr(2)-lonjr(1))

    ! latitude
    call nf90_err( nf90_inq_dimid(ncid,'latitude',dimid) )
    call nf90_err( nf90_inquire_dimension(ncid,dimid,len=nyjr) )
    if (.not.allocated(latjr)) allocate(latjr(nyjr))
    call nf90_err( nf90_inq_varid(ncid,'latitude',varid) )
    call nf90_err( nf90_get_var(ncid,varid,latjr) )
    dyjr=abs(latjr(2)-latjr(1))

    ! jrate field
    call nf90_err( nf90_inq_varid(ncid,'jrate',varid) )
    if (.not.allocated(jrate_average)) then
      allocate(jrate_average(nxjr,nyjr,2))
    endif
    call nf90_err( nf90_get_var(ncid,varid,jrate_average(:,:,memid),start=(/1,1,mm/),count=(/nxjr,nyjr,1/)) )

    ! close file
    call nf90_err( nf90_close(ncid) )

    return

  end subroutine readjrate

  function photo_O1D(sza) result(jrate)

  !*****************************************************************************
  !                                                                            *
  !    Calculates J(O1D) photolysis rate based on solar zenith angle           *
  !                                                                            *
  !    Author: A. Stohl, Nov 2014                                              *
  !                                                                            *
  !*****************************************************************************
  !                                                                            *
  !    INPUT:                                                                  *
  !    sza        solar zenith angle (degrees)                                 *
  !                                                                            *
  !    OUTPUT:                                                                 *
  !    photo_O1D  J(O1D) photoylsis rate                                       *
  !                                                                            *
  !*****************************************************************************

    implicit none

    real, intent(in)  :: sza
    real :: jrate
    integer :: iz,ik
    real :: z1,z2,zg,f1,f2,dummy
    real :: photo_NO2
    integer, parameter :: nzenith=11
    real, parameter :: pi=3.1415927
    real, dimension(nzenith) :: zangle,fact_photo

    ! zangle: zenith angles for which fact_photo is tabulated
    ! fact_photo: conversion of photolysis rate of NO2 to photolysis 
    !     rate of O3 into O1D as a function of solar zenith angle

    zangle=(/0.,10.,20.,30.,40.,50.,60.,70.,78.,86.,90.0001/)
    fact_photo=(/0.4616E-02,0.4478E-02,0.4131E-02,0.3583E-02,0.2867E-02,&
      &0.2081E-02,0.1235E-02,0.5392E-03,0.2200E-03,0.1302E-03,0.0902E-03/)

    if (sza.lt.90.) then
      do iz=1,nzenith-1
        if(sza.ge.zangle(iz)) ik=iz
      end do
      z1=1./cos(zangle(ik)*pi/180.)
      z2=1./cos(zangle(ik+1)*pi/180.)
      zg=1./cos(sza*pi/180.)
      dummy=(zg-z1)/(z2-z1)
      f1=alog(fact_photo(ik))
      f2=alog(fact_photo(ik+1))
      photo_NO2=1.45e-2*exp(-0.4/cos(sza*pi/180.))
      jrate=photo_NO2*exp(f1+(f2-f1)*dummy)
    else
      jrate=0.
    endif

  end function photo_O1D


  function zenithangle(ylat,xlon,jjjjmmdd,hhmmss) result(sza)

  !*****************************************************************************
  !                                                                            *
  !    This function returns the sinus of solar elevation as a function        *
  !    of geographic longitde, latitude and GMT-Time.                          *
  !                                                                            *
  !    Author: G. Wotawa, Nov 1993                                             *
  !                                                                            *
  !*****************************************************************************
  !                                                                            *
  !     INPUT:                                                                 *
  !                                                                            *
  !            ylat          geographical latitude  [DEG]                      *
  !            xlon          geographical longitude [DEG]                      *
  !            jjjj          Year                                              *
  !            mm            Month                                             *
  !            dd            Day                                               *
  !            hh            Hour                                              *
  !            minute        Minute                                            *
  !                                                                            *
  !*****************************************************************************

    implicit none

    integer,       intent(in)  :: jjjjmmdd,hhmmss
    real,          intent(in)  :: ylat,xlon
    real    :: sza
    integer :: jjjj,mm,id,iu,minute 
    integer :: ndaynum
    real :: sinsol,solelev
    real :: rnum,rylat,ttime,dekl,rdekl,eq
    real,parameter :: pi=3.1415927

    jjjj=jjjjmmdd/10000
    mm=jjjjmmdd/100-jjjj*100
    id=jjjjmmdd-jjjj*10000-mm*100
    iu=hhmmss/10000
    minute=hhmmss/100-100*iu

    ndaynum=31*(mm-1)+id
    if(mm.gt.2) ndaynum=ndaynum-int(0.4*mm+2.3)
    if((mm.gt.2).and.(jjjj/4*4.eq.jjjj)) ndaynum=ndaynum+1

    rnum=2.*pi*ndaynum/365.
    rylat=pi*ylat/180.
    ttime=real(iu)+real(minute)/60.

    dekl=0.396+3.631*sin(rnum)+0.038*sin(2.*rnum)+0.077*sin(3.*rnum)- &
         22.97*cos(rnum)-0.389*cos(2.*rnum)-0.158*cos(3.*rnum)
    rdekl=pi*dekl/180.

    eq=(0.003-7.343*sin(rnum)-9.47*sin(2.*rnum)- &
         0.329*sin(3.*rnum)-0.196*sin(4.*rnum)+ &
         0.552*cos(rnum)-3.020*cos(2.*rnum)- &
         0.076*cos(3.*rnum)-0.125*cos(4.*rnum))/60.
    sinsol=sin(rylat)*sin(rdekl)+cos(rylat)*cos(rdekl)* &
         cos((ttime-12.+xlon/15.+eq)*pi/12.)
    ! Calculate the maximum solar elevation on that day
    !sinsol=sin(rylat)*sin(rdekl)+cos(rylat)*cos(rdekl)*
    !    &       cos((eq)*pi/12.)
    solelev=asin(sinsol)*180./pi
    sza=90.-solelev

  end function zenithangle

  function find_minloc(n,vect,val) result(ind)

  !*****************************************************************************
  !                                                                            *
  !    This function returns index of a vector corresponding to the closest    *
  !    value to variable val                                                   *
  !                                                                            *
  !    Author: R. Thompson, Mar 2024                                           *
  !                                                                            *
  !*****************************************************************************

    implicit none

    integer,            intent(in) :: n
    real, dimension(n), intent(in) :: vect
    real,               intent(in) :: val
    integer                        :: i, ind
    real                           :: diff, diff_min

    diff_min=huge(1.0)

    ind=0
    do i=1,n
      diff=abs(vect(i)-val)
      if ( diff.lt.diff_min ) then
        diff_min=diff
        ind=i
      endif
    end do    

  end function find_minloc

end module chemistry_mod

