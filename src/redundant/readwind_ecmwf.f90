! SPDX-FileCopyrightText: FLEXPART 1998-2019, see flexpart_license.txt
! SPDX-License-Identifier: GPL-3.0-or-later

subroutine readwind_ecmwf(indj,n,uuh,vvh,wwh)

!**********************************************************************
!                                                                     *
!             TRAJECTORY MODEL SUBROUTINE READWIND                    *
!                                                                     *
!**********************************************************************
!                                                                     *
!             AUTHOR:      G. WOTAWA                                  *
!             DATE:        1997-08-05                                 *
!             LAST UPDATE: 2000-10-17, Andreas Stohl                  *
!             CHANGE: 11/01/2008, Harald Sodemann, GRIB1/2 input with *
!                                 ECMWF grib_api                      *
!             CHANGE: 03/12/2008, Harald Sodemann, update to f90 with *
!                                 ECMWF grib_api                      *
!                                                                     *
!**********************************************************************
!  Changes, Bernd C. Krueger, Feb. 2001:
!   Variables tth and qvh (on eta coordinates) in common block
!
!   Unified ECMWF and GFS builds                                      
!   Marian Harustak, 12.5.2017                                        
!     - Renamed from readwind to readwind_ecmwf                     
!
!**********************************************************************
!                                                                     *
! DESCRIPTION:                                                        *
!                                                                     *
! READING OF ECMWF METEOROLOGICAL FIELDS FROM INPUT DATA FILES. THE   *
! INPUT DATA FILES ARE EXPECTED TO BE AVAILABLE IN GRIB CODE          *
!                                                                     *
! INPUT:                                                              *
! indj               indicates number of the wind field to be read in *
! n                  temporal index for meteorological fields (1 to 3)*
!                                                                     *
! IMPORTANT VARIABLES FROM COMMON BLOCK:                              *
!                                                                     *
! wfname             File name of data to be read in                  *
! nx,ny,nuvz,nwz     expected field dimensions                        *
! nlev_ec            number of vertical levels ecmwf model            *
! uu,vv,ww           wind fields                                      *
! tt,qv              temperature and specific humidity                *
! ps                 surface pressure                                 *
!                                                                     *
!**********************************************************************

  use grib_api
  use par_mod
  use com_mod

  implicit none

!  include 'grib_api.h'

!HSO  parameters for grib_api
  integer :: ifile
  integer :: iret
  integer, dimension(:), allocatable   :: igrib
  integer :: nfield, ii, arsize
  integer :: gribVer,parCat,parNum,typSurf,valSurf,discipl,parId
  integer :: gotGrid
!HSO  end

  real(kind=4) :: uuh(0:nxmax-1,0:nymax-1,nuvzmax)
  real(kind=4) :: vvh(0:nxmax-1,0:nymax-1,nuvzmax)
  real(kind=4) :: wwh(0:nxmax-1,0:nymax-1,nwzmax)
  integer :: indj,i,j,k,n,levdiff2,iumax,iwmax!,ifield
  integer :: kz

! VARIABLES AND ARRAYS NEEDED FOR GRIB DECODING

! dimension of isec2 at least (22+n), where n is the number of parallels or
! meridians in a quasi-regular (reduced) Gaussian or lat/long grid

! dimension of zsec2 at least (10+nn), where nn is the number of vertical
! coordinate parameters

  integer :: isec1(56),isec2(22+nxmax+nymax)
  real(kind=4), allocatable, dimension(:) :: zsec4
!  real(kind=4) :: zsec4(jpunp)
  real(kind=4) :: xaux,yaux,xaux0,yaux0
  real(kind=8) :: xauxin,yauxin
  real,parameter :: eps=1.e-4
  real(kind=4) :: nsss(0:nxmax-1,0:nymax-1),ewss(0:nxmax-1,0:nymax-1)
  real :: plev1,pmean,tv,fu,hlev1,ff10m,fflev1,conversion_factor
  integer :: stat

  logical :: hflswitch,strswitch!,readcloud

!HSO  grib api error messages
  character(len=24) :: gribErrorMsg = 'Error reading grib file'
  character(len=20) :: gribFunction = 'readwind'

  hflswitch=.false.
  strswitch=.false.
!ZHG test the grib fields that have lcwc without using them
!  readcloud=.false.

  levdiff2=nlev_ec-nwz+1
  iumax=0
  iwmax=0

!
! OPENING OF DATA FILE (GRIB CODE)
!
5 call grib_open_file(ifile,path(3)(1:length(3)) &
       //trim(wfname(indj)),'r',iret)
  if (iret.ne.GRIB_SUCCESS) then
    goto 888   ! ERROR DETECTED
  endif

  call grib_count_in_file(ifile,nfield)

  ! allocate memory for grib handles
  allocate(igrib(nfield), stat=stat)
  if (stat.ne.0) stop "Could not allocate igrib"
  ! initialise
  igrib(:) = -1

  do ii = 1,nfield
    call grib_new_from_file(ifile, igrib(ii), iret)
  end do

  call grib_close_file(ifile)

!turn on support for multi fields messages */
!call grib_multi_support_on

  gotGrid=0

!$OMP PARALLEL DEFAULT(none) &
!$OMP SHARED (nfield, igrib, gribFunction, nxfield, ny, nlev_ec, dx, xlon0, ylat0, &
!$OMP   n, tth, uuh, vvh, iumax, qvh, ps, wwh, iwmax, sd, msl, tcc, u10, v10, tt2, &
!$OMP   td2, lsprec, convprec, sshf, hflswitch, ssr, ewss, nsss, strswitch, oro,   &
!$OMP   excessoro, lsm, nymin1,ciwch,clwch,readclouds,sumclouds) & 
!$OMP PRIVATE(ii, gribVer, iret, isec1, discipl, parCat, parNum, parId,typSurf, valSurf, &
!$OMP   zsec4, isec2, gribErrorMsg, xauxin, yauxin, xaux, yaux, xaux0,  &
!$OMP   yaux0, k, arsize, stat, conversion_factor)  &
!$OMP REDUCTION(+:gotGrid)
!
! GET NEXT FIELDS
!
  ! allocate memory for reading from grib
  allocate(zsec4(nxfield*ny), stat=stat)
  if (stat.ne.0) stop "Could not allocate zsec4"

!$OMP DO SCHEDULE(static)

  fieldloop : do ii=1,nfield

!first see if we read GRIB1 or GRIB2
  call grib_get_int(igrib(ii),'editionNumber',gribVer,iret)
  call grib_check(iret,gribFunction,gribErrorMsg)

  if (gribVer.eq.1) then ! GRIB Edition 1

    !print*,'GRiB Edition 1'
    !read the grib2 identifiers
    call grib_get_int(igrib(ii),'indicatorOfParameter',isec1(6),iret)
    call grib_check(iret,gribFunction,gribErrorMsg)
    call grib_get_int(igrib(ii),'level',isec1(8),iret)
    call grib_check(iret,gribFunction,gribErrorMsg)

!change code for etadot to code for omega
    if (isec1(6).eq.77) then
      isec1(6)=135
    endif

    conversion_factor=1.

  else

!print*,'GRiB Edition 2'
!read the grib2 identifiers
    call grib_get_int(igrib(ii),'discipline',discipl,iret)
    call grib_check(iret,gribFunction,gribErrorMsg)
    call grib_get_int(igrib(ii),'parameterCategory',parCat,iret)
    call grib_check(iret,gribFunction,gribErrorMsg)
    call grib_get_int(igrib(ii),'parameterNumber',parNum,iret)
    call grib_check(iret,gribFunction,gribErrorMsg)
    call grib_get_int(igrib(ii),'typeOfFirstFixedSurface',typSurf,iret)
    call grib_check(iret,gribFunction,gribErrorMsg)
    call grib_get_int(igrib(ii),'level',valSurf,iret)
    call grib_check(iret,gribFunction,gribErrorMsg)
    call grib_get_int(igrib(ii),'paramId',parId,iret)
    call grib_check(iret,gribFunction,gribErrorMsg)

!print*,discipl,parCat,parNum,typSurf,valSurf

!convert to grib1 identifiers
    isec1(6)=-1
    isec1(7)=-1
    isec1(8)=-1
    isec1(8)=valSurf     ! level
    conversion_factor=1.
    if ((parCat.eq.0).and.(parNum.eq.0).and.(typSurf.eq.105)) then ! T
      isec1(6)=130         ! indicatorOfParameter
    elseif ((parCat.eq.2).and.(parNum.eq.2).and.(typSurf.eq.105)) then ! U
      isec1(6)=131         ! indicatorOfParameter
    elseif ((parCat.eq.2).and.(parNum.eq.3).and.(typSurf.eq.105)) then ! V
      isec1(6)=132         ! indicatorOfParameter
    elseif ((parCat.eq.1).and.(parNum.eq.0).and.(typSurf.eq.105)) then ! Q
      isec1(6)=133         ! indicatorOfParameter
! ESO Cloud water is in a) fields CLWC and CIWC, *or* b) field QC 
    elseif ((parCat.eq.1).and.(parNum.eq.83).and.(typSurf.eq.105)) then ! clwc
      isec1(6)=246         ! indicatorOfParameter
    elseif ((parCat.eq.1).and.(parNum.eq.84).and.(typSurf.eq.105)) then ! ciwc
      isec1(6)=247         ! indicatorOfParameter
! ESO qc(=clwc+ciwc):
    elseif ((parCat.eq.201).and.(parNum.eq.31).and.(typSurf.eq.105)) then ! qc
      isec1(6)=201031         ! indicatorOfParameter
    elseif ((parCat.eq.3).and.(parNum.eq.0).and.(typSurf.eq.1)) then !SP
      isec1(6)=134         ! indicatorOfParameter
    elseif ((parCat.eq.2).and.(parNum.eq.32)) then ! W, actually eta dot
      isec1(6)=135         ! indicatorOfParameter
    elseif ((parCat.eq.128).and.(parNum.eq.77)) then ! W, actually eta dot
      isec1(6)=135         ! indicatorOfParameter
    elseif ((parCat.eq.3).and.(parNum.eq.0).and.(typSurf.eq.101)) then !SLP
      isec1(6)=151         ! indicatorOfParameter
    elseif ((parCat.eq.2).and.(parNum.eq.2).and.(typSurf.eq.103)) then ! 10U
      isec1(6)=165         ! indicatorOfParameter
    elseif ((parCat.eq.2).and.(parNum.eq.3).and.(typSurf.eq.103)) then ! 10V
      isec1(6)=166         ! indicatorOfParameter
    elseif ((parCat.eq.0).and.(parNum.eq.0).and.(typSurf.eq.103)) then ! 2T
      isec1(6)=167         ! indicatorOfParameter
    elseif ((parCat.eq.0).and.(parNum.eq.6).and.(typSurf.eq.103)) then ! 2D
      isec1(6)=168         ! indicatorOfParameter
    elseif ((parCat.eq.1).and.(parNum.eq.11).and.(typSurf.eq.1)) then ! SD
      isec1(6)=141         ! indicatorOfParameter
      conversion_factor=1000.
    elseif ((parCat.eq.6).and.(parNum.eq.1) .or. parId .eq. 164) then ! CC
      isec1(6)=164         ! indicatorOfParameter
    elseif ((parCat.eq.1).and.(parNum.eq.9) .or. parId .eq. 142) then ! LSP
      isec1(6)=142         ! indicatorOfParameter
    elseif ((parCat.eq.1).and.(parNum.eq.10)) then ! CP
      isec1(6)=143         ! indicatorOfParameter
      conversion_factor=1000.
    elseif ((parCat.eq.0).and.(parNum.eq.11).and.(typSurf.eq.1)) then ! SHF
      isec1(6)=146         ! indicatorOfParameter
    elseif ((parCat.eq.4).and.(parNum.eq.9).and.(typSurf.eq.1)) then ! SR
      isec1(6)=176         ! indicatorOfParameter
!    elseif ((parCat.eq.2).and.(parNum.eq.17) .or. parId .eq. 180) then ! EWSS --wrong
    elseif ((parCat.eq.2).and.(parNum.eq.38) .or. parId .eq. 180) then ! EWSS --correct
      isec1(6)=180         ! indicatorOfParameter
!    elseif ((parCat.eq.2).and.(parNum.eq.18) .or. parId .eq. 181) then ! NSSS --wrong
    elseif ((parCat.eq.2).and.(parNum.eq.37) .or. parId .eq. 181) then ! NSSS --correct
      isec1(6)=181         ! indicatorOfParameter
    elseif ((parCat.eq.3).and.(parNum.eq.4)) then ! ORO
      isec1(6)=129         ! indicatorOfParameter
    elseif ((parCat.eq.3).and.(parNum.eq.7) .or. parId .eq. 160) then ! SDO
      isec1(6)=160         ! indicatorOfParameter
    elseif ((discipl.eq.2).and.(parCat.eq.0).and.(parNum.eq.0).and. &
         (typSurf.eq.1)) then ! LSM
      isec1(6)=172         ! indicatorOfParameter
    elseif (parNum.eq.152) then 
      isec1(6)=152         ! avoid warning for lnsp    
    else
      print*,'***WARNING: undefined GRiB2 message found!',discipl, &
           parCat,parNum,typSurf
    endif
    if(parId .ne. isec1(6) .and. parId .ne. 77) then
      write(*,*) 'parId',parId, 'isec1(6)',isec1(6)
!    stop
    endif

  endif

!HSO  get the size and data of the values array
  if (isec1(6).ne.-1) then
    call grib_get_real4_array(igrib(ii),'values',zsec4,iret)
    call grib_check(iret,gribFunction,gribErrorMsg)
  endif

!HSO  get the required fields from section 2 in a gribex compatible manner
  if (ii.eq.1) then
    call grib_get_int(igrib(ii),'numberOfPointsAlongAParallel',isec2(2),iret)
    call grib_check(iret,gribFunction,gribErrorMsg)
    call grib_get_int(igrib(ii),'numberOfPointsAlongAMeridian',isec2(3),iret)
    call grib_check(iret,gribFunction,gribErrorMsg)
    call grib_get_int(igrib(ii),'numberOfVerticalCoordinateValues',isec2(12))
    call grib_check(iret,gribFunction,gribErrorMsg)
! CHECK GRID SPECIFICATIONS
    if(isec2(2).ne.nxfield) stop 'READWIND: NX NOT CONSISTENT'
    if(isec2(3).ne.ny) stop 'READWIND: NY NOT CONSISTENT'
    if(isec2(12)/2-1.ne.nlev_ec) &
         stop 'READWIND: VERTICAL DISCRETIZATION NOT CONSISTENT'
  endif ! ifield

!$OMP CRITICAL
!HSO  get the second part of the grid dimensions only from GRiB1 messages
  if (isec1(6) .eq. 167 .and. (gotGrid.eq.0)) then
    call grib_get_real8(igrib(ii),'longitudeOfFirstGridPointInDegrees', &
         xauxin,iret)
    call grib_check(iret,gribFunction,gribErrorMsg)
    call grib_get_real8(igrib(ii),'latitudeOfLastGridPointInDegrees', &
         yauxin,iret)
    call grib_check(iret,gribFunction,gribErrorMsg)
    if (xauxin.gt.180.) xauxin=xauxin-360.0
    if (xauxin.lt.-180.) xauxin=xauxin+360.0

    xaux=xauxin+real(nxshift)*dx
    yaux=yauxin
    if (xaux.gt.180.) xaux=xaux-360.0
    if(abs(xaux-xlon0).gt.eps) &
         stop 'READWIND: LOWER LEFT LONGITUDE NOT CONSISTENT'
    if(abs(yaux-ylat0).gt.eps) &
         stop 'READWIND: LOWER LEFT LATITUDE NOT CONSISTENT'
    gotGrid=1
  endif ! gotGrid
!$OMP END CRITICAL

  k=isec1(8)
  select case(isec1(6)) 
!! TEMPERATURE
    case(130) 
      do j=0,nymin1
        do i=0,nxfield-1
          tth(i,j,nlev_ec-k+2,n) = zsec4(nxfield*(ny-j-1)+i+1)
        end do 
      end do 
!! U VELOCITY
    case(131) 
      do j=0,nymin1
        do i=0,nxfield-1
          uuh(i,j,nlev_ec-k+2) = zsec4(nxfield*(ny-j-1)+i+1)
        end do 
      end do 
!$OMP CRITICAL
      iumax=max(iumax,nlev_ec-k+1)
!$OMP END CRITICAL
!! V VELOCITY
    case(132)
      do j=0,nymin1
        do i=0,nxfield-1
          vvh(i,j,nlev_ec-k+2) = zsec4(nxfield*(ny-j-1)+i+1)
        end do 
      end do
!! SPEC. HUMIDITY
    case(133)
      do j=0,nymin1
        do i=0,nxfield-1
          qvh(i,j,nlev_ec-k+2,n) = zsec4(nxfield*(ny-j-1)+i+1)
          if (qvh(i,j,nlev_ec-k+2,n) .lt. 0.) &
               qvh(i,j,nlev_ec-k+2,n) = 0.
    !        this is necessary because the gridded data may contain
    !        spurious negative values
        end do 
      end do
!! SURF. PRESS.
    case(134) 
      do j=0,nymin1
        do i=0,nxfield-1
          ps(i,j,1,n) = zsec4(nxfield*(ny-j-1)+i+1)
        end do 
      end do
!! W VELOCITY
    case(135)
      do j=0,nymin1
        do i=0,nxfield-1
          wwh(i,j,nlev_ec-k+1) = zsec4(nxfield*(ny-j-1)+i+1)
        end do 
      end do 
!$OMP CRITICAL
      iwmax=max(iwmax,nlev_ec-k+1)
!$OMP END CRITICAL
!! SNOW DEPTH
    case(141) 
      do j=0,nymin1
        do i=0,nxfield-1
          sd(i,j,1,n)= zsec4(nxfield*(ny-j-1)+i+1)/conversion_factor
        end do 
      end do
!! SEA LEVEL PRESS.
    case(151)
      do j=0,nymin1
        do i=0,nxfield-1
          msl(i,j,1,n) = zsec4(nxfield*(ny-j-1)+i+1)
        end do 
      end do
!! CLOUD COVER
    case(164)
      do j=0,nymin1
        do i=0,nxfield-1
          tcc(i,j,1,n) =  zsec4(nxfield*(ny-j-1)+i+1)
        end do 
      end do 
!! 10 M U VELOCITY
    case(165) 
      do j=0,nymin1
        do i=0,nxfield-1
          u10(i,j,1,n)= zsec4(nxfield*(ny-j-1)+i+1)
        end do 
      end do 
!! 10 M V VELOCITY
    case(166) 
      do j=0,nymin1
        do i=0,nxfield-1
          v10(i,j,1,n) = zsec4(nxfield*(ny-j-1)+i+1)
        end do 
      end do 
!! 2 M TEMPERATURE
    case(167)
      do j=0,nymin1
        do i=0,nxfield-1
          tt2(i,j,1,n) = zsec4(nxfield*(ny-j-1)+i+1)
        end do 
      end do 
!! 2 M DEW POINT
    case(168)
      do j=0,nymin1
        do i=0,nxfield-1
          td2(i,j,1,n) = zsec4(nxfield*(ny-j-1)+i+1)
        end do 
      end do
!! LARGE SCALE PREC.
    case(142)
      do j=0,nymin1
        do i=0,nxfield-1
          lsprec(i,j,1,n)=zsec4(nxfield*(ny-j-1)+i+1)
          if (lsprec(i,j,1,n).lt.0.) lsprec(i,j,1,n)=0.
        end do 
      end do
!! CONVECTIVE PREC.
    case(143)
      do j=0,nymin1
        do i=0,nxfield-1
          convprec(i,j,1,n)=zsec4(nxfield*(ny-j-1)+i+1)/conversion_factor
          if (convprec(i,j,1,n).lt.0.) convprec(i,j,1,n)=0.
        end do 
      end do
!! SENS. HEAT FLUX
    case(146)
      do j=0,nymin1
        do i=0,nxfield-1
          sshf(i,j,1,n) = zsec4(nxfield*(ny-j-1)+i+1)
!$OMP CRITICAL
          if(zsec4(nxfield*(ny-j-1)+i+1).ne.0.) &  
            hflswitch=.true.    ! Heat flux available
!$OMP END CRITICAL
        end do 
      end do
!! SOLAR RADIATION
    case(176)
      do j=0,nymin1
        do i=0,nxfield-1
          ssr(i,j,1,n)=zsec4(nxfield*(ny-j-1)+i+1)
          if (ssr(i,j,1,n).lt.0.) ssr(i,j,1,n)=0.
        end do 
      end do
!! EW SURFACE STRESS
    case(180)
      do j=0,nymin1
        do i=0,nxfield-1
          ewss(i,j) = zsec4(nxfield*(ny-j-1)+i+1)
!$OMP CRITICAL
          if (zsec4(nxfield*(ny-j-1)+i+1).ne.0.) strswitch=.true.    ! stress available
!$OMP END CRITICAL
        end do 
      end do 
!! NS SURFACE STRESS
    case(181)
     do j=0,nymin1
       do i=0,nxfield-1
         nsss(i,j) = zsec4(nxfield*(ny-j-1)+i+1)
!$OMP CRITICAL
         if (zsec4(nxfield*(ny-j-1)+i+1).ne.0.) strswitch=.true.    ! stress available
!$OMP END CRITICAL
       end do 
     end do 
!! ECMWF OROGRAPHY
    case(129)
      do j=0,nymin1
        do i=0,nxfield-1
          oro(i,j) = zsec4(nxfield*(ny-j-1)+i+1)/ga
        end do 
      end do
!! STANDARD DEVIATION OF OROGRAPHY
    case(160)
      do j=0,nymin1
        do i=0,nxfield-1
          excessoro(i,j) = zsec4(nxfield*(ny-j-1)+i+1)
        end do 
      end do 
!! ECMWF LAND SEA MASK
    case(172)
      do j=0,nymin1
        do i=0,nxfield-1
          lsm(i,j) = zsec4(nxfield*(ny-j-1)+i+1)
        end do 
      end do
!! CLWC  Cloud liquid water content [kg/kg]
    case(246)    
      do j=0,nymin1
        do i=0,nxfield-1
          clwch(i,j,nlev_ec-k+2,n)=zsec4(nxfield*(ny-j-1)+i+1)
        end do
      end do
!$OMP CRITICAL
      readclouds=.true.
      sumclouds=.false.
!$OMP END CRITICAL
!! CIWC  Cloud ice water content
    case(247)
      do j=0,nymin1
        do i=0,nxfield-1
          ciwch(i,j,nlev_ec-k+2,n)=zsec4(nxfield*(ny-j-1)+i+1)
        end do
      end do
!ZHG end
!ESO read qc (=clwc+ciwc)
!! QC  Cloud liquid water content [kg/kg]
    case(201031)
      do j=0,nymin1
        do i=0,nxfield-1      
          clwch(i,j,nlev_ec-k+2,n)=zsec4(nxfield*(ny-j-1)+i+1)
        end do
      end do
!$OMP CRITICAL
      readclouds=.true.
      sumclouds=.false.
!$OMP END CRITICAL

  end select
  ! do j=0,nymin1
  !   do i=0,nxfield-1
  !     k=isec1(8)
      ! if(isec1(6).eq.130) tth(i,j,nlev_ec-k+2,n)= &!! TEMPERATURE
      !      zsec4(nxfield*(ny-j-1)+i+1)
      ! if(isec1(6).eq.131) uuh(i,j,nlev_ec-k+2)= &!! U VELOCITY
      !      zsec4(nxfield*(ny-j-1)+i+1)
      ! if(isec1(6).eq.132) vvh(i,j,nlev_ec-k+2)= &!! V VELOCITY
      !      zsec4(nxfield*(ny-j-1)+i+1)
!       if(isec1(6).eq.133) then                      !! SPEC. HUMIDITY
!         qvh(i,j,nlev_ec-k+2,n)=zsec4(nxfield*(ny-j-1)+i+1)
!         if (qvh(i,j,nlev_ec-k+2,n) .lt. 0.) &
!              qvh(i,j,nlev_ec-k+2,n) = 0.
! !        this is necessary because the gridded data may contain
! !        spurious negative values
!       endif
      ! if(isec1(6).eq.134) ps(i,j,1,n)= &!! SURF. PRESS.
      !      zsec4(nxfield*(ny-j-1)+i+1)

      ! if(isec1(6).eq.135) wwh(i,j,nlev_ec-k+1)= &!! W VELOCITY
      !      zsec4(nxfield*(ny-j-1)+i+1)
      ! if(isec1(6).eq.141) sd(i,j,1,n)= &!! SNOW DEPTH
      !      zsec4(nxfield*(ny-j-1)+i+1)/conversion_factor
      ! if(isec1(6).eq.151) msl(i,j,1,n)= &!! SEA LEVEL PRESS.
      !      zsec4(nxfield*(ny-j-1)+i+1)
      ! if(isec1(6).eq.164) tcc(i,j,1,n)= &!! CLOUD COVER
      !      zsec4(nxfield*(ny-j-1)+i+1)
      ! if(isec1(6).eq.165) u10(i,j,1,n)= &!! 10 M U VELOCITY
      !      zsec4(nxfield*(ny-j-1)+i+1)
      ! if(isec1(6).eq.166) v10(i,j,1,n)= &!! 10 M V VELOCITY
      !      zsec4(nxfield*(ny-j-1)+i+1)
      ! if(isec1(6).eq.167) tt2(i,j,1,n)= &!! 2 M TEMPERATURE
      !      zsec4(nxfield*(ny-j-1)+i+1)
      ! if(isec1(6).eq.168) td2(i,j,1,n)= &!! 2 M DEW POINT
      !      zsec4(nxfield*(ny-j-1)+i+1)
      ! if(isec1(6).eq.142) then                      !! LARGE SCALE PREC.
      !   lsprec(i,j,1,n)=zsec4(nxfield*(ny-j-1)+i+1)
      !   if (lsprec(i,j,1,n).lt.0.) lsprec(i,j,1,n)=0.
      ! endif
      ! if(isec1(6).eq.143) then                      !! CONVECTIVE PREC.
      !   convprec(i,j,1,n)=zsec4(nxfield*(ny-j-1)+i+1)/conversion_factor
      !   if (convprec(i,j,1,n).lt.0.) convprec(i,j,1,n)=0.
      ! endif
      ! if(isec1(6).eq.146) sshf(i,j,1,n)= &!! SENS. HEAT FLUX
      !      zsec4(nxfield*(ny-j-1)+i+1)
      ! if((isec1(6).eq.146).and.(zsec4(nxfield*(ny-j-1)+i+1).ne.0.)) &
      !      hflswitch=.true.    ! Heat flux available
      ! if(isec1(6).eq.176) then                      !! SOLAR RADIATION
      !   ssr(i,j,1,n)=zsec4(nxfield*(ny-j-1)+i+1)
      !   if (ssr(i,j,1,n).lt.0.) ssr(i,j,1,n)=0.
      ! endif
!       if(isec1(6).eq.180) ewss(i,j)= &!! EW SURFACE STRESS
!            zsec4(nxfield*(ny-j-1)+i+1)
!       if(isec1(6).eq.181) nsss(i,j)= &!! NS SURFACE STRESS
!            zsec4(nxfield*(ny-j-1)+i+1)
!       if(((isec1(6).eq.180).or.(isec1(6).eq.181)).and. &
!            (zsec4(nxfield*(ny-j-1)+i+1).ne.0.)) strswitch=.true.    ! stress available
! !sec        strswitch=.true.
      ! if(isec1(6).eq.129) oro(i,j)= &!! ECMWF OROGRAPHY
      !      zsec4(nxfield*(ny-j-1)+i+1)/ga
      ! if(isec1(6).eq.160) excessoro(i,j)= &!! STANDARD DEVIATION OF OROGRAPHY
      !      zsec4(nxfield*(ny-j-1)+i+1)
      ! if(isec1(6).eq.172) lsm(i,j)= &!! ECMWF LAND SEA MASK
      !      zsec4(nxfield*(ny-j-1)+i+1)
      ! if(isec1(6).eq.131) iumax=max(iumax,nlev_ec-k+1)
      ! if(isec1(6).eq.135) iwmax=max(iwmax,nlev_ec-k+1)
!ZHG READING CLOUD FIELDS ASWELL
! ESO TODO: add check for if one of clwc/ciwc missing (error),
! also if all 3 cw fields present, use qc and disregard the others
!       if(isec1(6).eq.246) then  !! CLWC  Cloud liquid water content [kg/kg]
!         clwch(i,j,nlev_ec-k+2,n)=zsec4(nxfield*(ny-j-1)+i+1)
!         readclouds=.true.
!         sumclouds=.false.
!       endif
!       if(isec1(6).eq.247) then  !! CIWC  Cloud ice water content
!         ciwch(i,j,nlev_ec-k+2,n)=zsec4(nxfield*(ny-j-1)+i+1)
!       endif
! !ZHG end
! !ESO read qc (=clwc+ciwc)
!       if(isec1(6).eq.201031) then  !! QC  Cloud liquid water content [kg/kg]
!         clwch(i,j,nlev_ec-k+2,n)=zsec4(nxfield*(ny-j-1)+i+1)
!         readclouds=.true.
!         sumclouds=.true.
!       endif

!     end do
!   end do

  call grib_release(igrib(ii))

  end do fieldloop
!$OMP END DO
  deallocate(zsec4)
!$OMP END PARALLEL

  deallocate(igrib)
!
! CLOSING OF INPUT DATA FILE
!

! 50 call grib_close_file(ifile)

!error message if no fields found with correct first longitude in it
  if (gotGrid.eq.0) then
    print*,'***ERROR: input file needs to contain GRiB1 formatted'// &
         'messages'
    stop
  endif

  if(levdiff2.eq.0) then
    iwmax=nlev_ec+1
    do i=0,nxmin1
      do j=0,nymin1
        wwh(i,j,nlev_ec+1)=0.
      end do
    end do
  endif

! For global fields, assign the leftmost data column also to the rightmost
! data column; if required, shift whole grid by nxshift grid points
!*************************************************************************

  if (xglobal) then
    call shift_field_0(ewss,nxfield,ny)
    call shift_field_0(nsss,nxfield,ny)
    call shift_field_0(oro,nxfield,ny)
    call shift_field_0(excessoro,nxfield,ny)
    call shift_field_0(lsm,nxfield,ny)
    call shift_field(ps,nxfield,ny,1,1,2,n)
    call shift_field(sd,nxfield,ny,1,1,2,n)
    call shift_field(msl,nxfield,ny,1,1,2,n)
    call shift_field(tcc,nxfield,ny,1,1,2,n)
    call shift_field(u10,nxfield,ny,1,1,2,n)
    call shift_field(v10,nxfield,ny,1,1,2,n)
    call shift_field(tt2,nxfield,ny,1,1,2,n)
    call shift_field(td2,nxfield,ny,1,1,2,n)
    call shift_field(lsprec,nxfield,ny,1,1,2,n)
    call shift_field(convprec,nxfield,ny,1,1,2,n)
    call shift_field(sshf,nxfield,ny,1,1,2,n)
    call shift_field(ssr,nxfield,ny,1,1,2,n)
    call shift_field(tth,nxfield,ny,nuvzmax,nuvz,2,n)
    call shift_field(qvh,nxfield,ny,nuvzmax,nuvz,2,n)
    call shift_field(uuh,nxfield,ny,nuvzmax,nuvz,1,1)
    call shift_field(vvh,nxfield,ny,nuvzmax,nuvz,1,1)
    call shift_field(wwh,nxfield,ny,nwzmax,nwz,1,1)
!ZHG
    call shift_field(clwch,nxfield,ny,nuvzmax,nuvz,2,n)
    if (.not.sumclouds) call shift_field(ciwch,nxfield,ny,nuvzmax,nuvz,2,n)
!ZHG end

  endif

  do i=0,nxmin1
    do j=0,nymin1
      surfstr(i,j,1,n)=sqrt(ewss(i,j)**2+nsss(i,j)**2)
    end do
  end do

  if ((.not.hflswitch).or.(.not.strswitch)) then
    write(*,*) 'WARNING: No flux data contained in GRIB file ', &
         wfname(indj)

! CALCULATE USTAR AND SSHF USING THE PROFILE METHOD
! As ECMWF has increased the model resolution, such that now the first model
! level is at about 10 m (where 10-m wind is given), use the 2nd ECMWF level
! (3rd model level in FLEXPART) for the profile method
!***************************************************************************

    do i=0,nxmin1
      do j=0,nymin1
        plev1=akz(3)+bkz(3)*ps(i,j,1,n)
        pmean=0.5*(ps(i,j,1,n)+plev1)
        tv=tth(i,j,3,n)*(1.+0.61*qvh(i,j,3,n))
        fu=-r_air*tv/ga/pmean
        hlev1=fu*(plev1-ps(i,j,1,n))   ! HEIGTH OF FIRST MODEL LAYER
        ff10m= sqrt(u10(i,j,1,n)**2+v10(i,j,1,n)**2)
        fflev1=sqrt(uuh(i,j,3)**2+vvh(i,j,3)**2)
        call pbl_profile(ps(i,j,1,n),td2(i,j,1,n),hlev1, &
             tt2(i,j,1,n),tth(i,j,3,n),ff10m,fflev1, &
             surfstr(i,j,1,n),sshf(i,j,1,n))
        if(sshf(i,j,1,n).gt.200.) sshf(i,j,1,n)=200.
        if(sshf(i,j,1,n).lt.-400.) sshf(i,j,1,n)=-400.
      end do
    end do
  endif


! Assign 10 m wind to model level at eta=1.0 to have one additional model
! level at the ground
! Specific humidity is taken the same as at one level above
! Temperature is taken as 2 m temperature
!**************************************************************************

  do i=0,nxmin1
    do j=0,nymin1
      uuh(i,j,1)=u10(i,j,1,n)
      vvh(i,j,1)=v10(i,j,1,n)
      qvh(i,j,1,n)=qvh(i,j,2,n)
      tth(i,j,1,n)=tt2(i,j,1,n)
    end do
  end do

  if(iumax.ne.nuvz-1) stop 'READWIND: NUVZ NOT CONSISTENT'
  if(iwmax.ne.nwz)    stop 'READWIND: NWZ NOT CONSISTENT'

  return

888 write(*,*) ' #### FLEXPART MODEL ERROR! WINDFIELD         #### '
  write(*,*) ' #### ',wfname(indj),'                    #### '
  write(*,*) ' #### IS NOT GRIB FORMAT !!!                  #### '
  stop 'Execution terminated'

end subroutine readwind_ecmwf

