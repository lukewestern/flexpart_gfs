! SPDX-FileCopyrightText: FLEXPART 1998-2019, see flexpart_license.txt
! SPDX-License-Identifier: GPL-3.0-or-later

module windfields_mod
  use par_mod
  use com_mod

	implicit none

  ! Fixed fields, unchangeable with time
  !*************************************

  real :: oro(0:nxmax-1,0:nymax-1)   ! orography of the ECMWF model
  real :: excessoro(0:nxmax-1,0:nymax-1)   ! excess orography mother domain
  real :: lsm(0:nxmax-1,0:nymax-1)   ! land sea mask of the ECMWF model
  

  ! 3d fields
  !**********
  real :: uueta(0:nxmax-1,0:nymax-1,nzmax,numwfmem)
  real :: vveta(0:nxmax-1,0:nymax-1,nzmax,numwfmem)
  real :: uupoleta(0:nxmax-1,0:nymax-1,nzmax,numwfmem)
  real :: vvpoleta(0:nxmax-1,0:nymax-1,nzmax,numwfmem)
  real :: wweta(0:nxmax-1,0:nymax-1,nzmax,numwfmem)
  real :: tteta(0:nxmax-1,0:nymax-1,nzmax,numwfmem)
  real :: qveta(0:nxmax-1,0:nymax-1,nzmax,numwfmem)
  real :: pveta(0:nxmax-1,0:nymax-1,nzmax,numwfmem)
  real :: rhoeta(0:nxmax-1,0:nymax-1,nzmax,numwfmem)
  real :: drhodzeta(0:nxmax-1,0:nymax-1,nzmax,numwfmem)
  real :: tvirtual(0:nxmax-1,0:nymax-1,nzmax,numwfmem)

  ! uu,vv,ww [m/2]       wind components in x,y and z direction
  real :: uu(0:nxmax-1,0:nymax-1,nzmax,numwfmem)
  real :: vv(0:nxmax-1,0:nymax-1,nzmax,numwfmem)
  ! uupol,vvpol [m/s]    wind components in polar stereographic projection
  real :: uupol(0:nxmax-1,0:nymax-1,nzmax,numwfmem)
  real :: vvpol(0:nxmax-1,0:nymax-1,nzmax,numwfmem)
  real :: ww(0:nxmax-1,0:nymax-1,nzmax,numwfmem)
  real :: tt(0:nxmax-1,0:nymax-1,nzmax,numwfmem) ! temperature data [K]
  real :: qv(0:nxmax-1,0:nymax-1,nzmax,numwfmem) ! specific humidity data

  ! ZHG adding cloud water 
  real :: clwc(0:nxmax-1,0:nymax-1,nzmax,numwfmem)=0.0 !liquid   [kg/kg]
  real :: ciwc(0:nxmax-1,0:nymax-1,nzmax,numwfmem)=0.0 !ice      [kg/kg]
  real :: clw(0:nxmax-1,0:nymax-1,nzmax,numwfmem)=0.0  !combined [m3/m3]
  ! RLT add pressure and dry air density
  real :: prs(0:nxmax-1,0:nymax-1,nzmax,numwfmem) ! air pressure
  real :: rho_dry(0:nxmax-1,0:nymax-1,nzmax,numwfmem)
  real :: pv(0:nxmax-1,0:nymax-1,nzmax,numwfmem) ! potential vorticity
  real :: rho(0:nxmax-1,0:nymax-1,nzmax,numwfmem) ! air density [kg/m3]
  real :: drhodz(0:nxmax-1,0:nymax-1,nzmax,numwfmem) ! vertical air density gradient [kg/m2]
  ! tth,qvh              tth,qvh on original eta levels
  real :: tth(0:nxmax-1,0:nymax-1,nuvzmax,numwfmem)
  real :: qvh(0:nxmax-1,0:nymax-1,nuvzmax,numwfmem)
  real :: clwch(0:nxmax-1,0:nymax-1,nuvzmax,numwfmem)=0.0
  real :: ciwch(0:nxmax-1,0:nymax-1,nuvzmax,numwfmem)=0.0

  real :: pplev(0:nxmax-1,0:nymax-1,nuvzmax,numwfmem)
  real :: ctwc(0:nxmax-1,0:nymax-1,numwfmem) ! ESO: =icloud_stats(:,:,4,:) total cloud water content
  
  ! clouds:   no cloud, no precipitation   0
  !      cloud, no precipitation      1
  !      rainout  conv/lsp dominated  2/3
  !      washout  conv/lsp dominated  4/5
  ! PS 2013
  !c icloudbot (m)        cloud bottom height
  !c icloudthck (m)       cloud thickness     

  !scavenging NIK, PS
  integer(kind=1) :: clouds(0:nxmax-1,0:nymax-1,nzmax,numwfmem)
  integer :: cloudsh(0:nxmax-1,0:nymax-1,numwfmem)

  ! 2d fields
  !**********
  real :: ps(0:nxmax-1,0:nymax-1,1,numwfmem)   ! surface pressure
  real :: sd(0:nxmax-1,0:nymax-1,1,numwfmem)   ! snow depth
  real :: msl(0:nxmax-1,0:nymax-1,1,numwfmem)   ! mean sea level pressure
  real :: tcc(0:nxmax-1,0:nymax-1,1,numwfmem)   ! total cloud cover
  real :: u10(0:nxmax-1,0:nymax-1,1,numwfmem)   ! 10 meter u
  real :: v10(0:nxmax-1,0:nymax-1,1,numwfmem)   ! 10 meter v
  real :: tt2(0:nxmax-1,0:nymax-1,1,numwfmem)   ! 2 meter temperature
  real :: td2(0:nxmax-1,0:nymax-1,1,numwfmem)   ! 2 meter dew point
  real :: lsprec(0:nxmax-1,0:nymax-1,1,numwfmem)   ! large scale total precipitation [mm/h]
  real :: convprec(0:nxmax-1,0:nymax-1,1,numwfmem)   ! convective precipitation [mm/h]
  real :: sshf(0:nxmax-1,0:nymax-1,1,numwfmem)   ! surface sensible heat flux
  real :: ssr(0:nxmax-1,0:nymax-1,1,numwfmem)   ! surface solar radiation
  real :: surfstr(0:nxmax-1,0:nymax-1,1,numwfmem)   ! surface stress
  real :: ustar(0:nxmax-1,0:nymax-1,1,numwfmem)   ! friction velocity [m/s]
  real :: wstar(0:nxmax-1,0:nymax-1,1,numwfmem)   ! convective velocity scale [m/s]
  real :: hmix(0:nxmax-1,0:nymax-1,1,numwfmem)   ! mixing height [m]
  real :: tropopause(0:nxmax-1,0:nymax-1,1,numwfmem)   ! altitude of thermal tropopause [m]
  real :: oli(0:nxmax-1,0:nymax-1,1,numwfmem)   ! inverse Obukhov length (1/L) [m]

  ! Save uv and w field heights for ETA coordinate system
  real,dimension(0:nxmax-1,0:nymax-1,nuvzmax,numwfmem) :: etauvheight,etawheight
contains

subroutine gridcheck_ecmwf

  !**********************************************************************
  !                                                                     *
  !             FLEXPART MODEL SUBROUTINE GRIDCHECK                     *
  !                                                                     *
  !**********************************************************************
  !                                                                     *
  !             AUTHOR:      G. WOTAWA                                  *
  !             DATE:        1997-08-06                                 *
  !             LAST UPDATE: 1997-10-10                                 *
  !                                                                     *
  !             Update:      1999-02-08, global fields allowed, A. Stohl*
  !             CHANGE: 11/01/2008, Harald Sodemann, GRIB1/2 input with *
  !                                 ECMWF grib_api                      *
  !             CHANGE: 03/12/2008, Harald Sodemann, update to f90 with *
  !                                 ECMWF grib_api                      *
  !                                                                     *
  !   Unified ECMWF and GFS builds                                      *
  !   Marian Harustak, 12.5.2017                                        *
  !     - Renamed from gridcheck to gridcheck_ecmwf                     *
  !                                                                     *
  !**********************************************************************
  !                                                                     *
  ! DESCRIPTION:                                                        *
  !                                                                     *
  ! THIS SUBROUTINE DETERMINES THE GRID SPECIFICATIONS (LOWER LEFT      *
  ! LONGITUDE, LOWER LEFT LATITUDE, NUMBER OF GRID POINTS, GRID DIST-   *
  ! ANCE AND VERTICAL DISCRETIZATION OF THE ECMWF MODEL) FROM THE       *
  ! GRIB HEADER OF THE FIRST INPUT FILE. THE CONSISTANCY (NO CHANGES    *
  ! WITHIN ONE FLEXPART RUN) IS CHECKED IN THE ROUTINE "READWIND" AT    *
  ! ANY CALL.                                                           *
  !                                                                     *
  ! XLON0                geographical longitude of lower left gridpoint *
  ! YLAT0                geographical latitude of lower left gridpoint  *
  ! NX                   number of grid points x-direction              *
  ! NY                   number of grid points y-direction              *
  ! DX                   grid distance x-direction                      *
  ! DY                   grid distance y-direction                      *
  ! NUVZ                 number of grid points for horizontal wind      *
  !                      components in z direction                      *
  ! NWZ                  number of grid points for vertical wind        *
  !                      component in z direction                       *
  ! sizesouth, sizenorth give the map scale (i.e. number of virtual grid*
  !                      points of the polar stereographic grid):       *
  !                      used to check the CFL criterion                *
  ! UVHEIGHT(1)-         heights of gridpoints where u and v are        *
  ! UVHEIGHT(NUVZ)       given                                          *
  ! WHEIGHT(1)-          heights of gridpoints where w is given         *
  ! WHEIGHT(NWZ)                                                        *
  !                                                                     *
  !**********************************************************************

  use grib_api
  use cmapf_mod, only: stlmbr,stcm2p

  implicit none

  !HSO  parameters for grib_api
  integer :: ifile
  integer :: iret
  integer :: igrib
  integer :: gotGrid
  real(kind=4) :: xaux1,xaux2,yaux1,yaux2
  real(kind=8) :: xaux1in,xaux2in,yaux1in,yaux2in
  integer :: gribVer,parCat,parNum,typSurf,valSurf,discipl,parId
  !HSO  end
  integer :: ix,jy,i,ifn,ifield,j,k,iumax,iwmax,numskip
  real :: sizesouth,sizenorth,xauxa,conversion_factor

  ! VARIABLES AND ARRAYS NEEDED FOR GRIB DECODING

  ! dimension of isec2 at least (22+n), where n is the number of parallels or
  ! meridians in a quasi-regular (reduced) Gaussian or lat/long grid

  ! dimension of zsec2 at least (10+nn), where nn is the number of vertical
  ! coordinate parameters

  integer :: isec1(56),isec2(22+nxmax+nymax)
  real(kind=4) :: zsec2(60+2*nuvzmax),zsec4(jpunp)
  character(len=1) :: opt

  !HSO  grib api error messages
  character(len=24) :: gribErrorMsg = 'Error reading grib file'
  character(len=20) :: gribFunction = 'gridcheck'


  iumax=0
  iwmax=0

  if(ideltas.gt.0) then
    ifn=1
  else
    ifn=numbwf
  endif
  !
  ! OPENING OF DATA FILE (GRIB CODE)
  !
5 call grib_open_file(ifile,path(3)(1:length(3)) &
       //trim(wfname(ifn)),'r',iret)
  if (iret.ne.GRIB_SUCCESS) then
    goto 999   ! ERROR DETECTED
  endif
  !turn on support for multi fields messages
  !call grib_multi_support_on

  gotGrid=0
  ifield=0
10 ifield=ifield+1

  !
  ! GET NEXT FIELDS
  !
  call grib_new_from_file(ifile,igrib,iret)
  if (iret.eq.GRIB_END_OF_FILE )  then
    goto 30    ! EOF DETECTED
  elseif (iret.ne.GRIB_SUCCESS) then
    goto 999   ! ERROR DETECTED
  endif

  !first see if we read GRIB1 or GRIB2
  call grib_get_int(igrib,'editionNumber',gribVer,iret)
  call grib_check(iret,gribFunction,gribErrorMsg)

  if (gribVer.eq.1) then ! GRIB Edition 1

    !print*,'GRiB Edition 1'
    !read the grib2 identifiers
    call grib_get_int(igrib,'indicatorOfParameter',isec1(6),iret)
    call grib_check(iret,gribFunction,gribErrorMsg)
    call grib_get_int(igrib,'level',isec1(8),iret)
    call grib_check(iret,gribFunction,gribErrorMsg)

    !change code for etadot to code for omega
    if (isec1(6).eq.77) then
      isec1(6)=135
    endif

    !print*,isec1(6),isec1(8)

  else

    !print*,'GRiB Edition 2'
    !read the grib2 identifiers
    call grib_get_int(igrib,'discipline',discipl,iret)
    call grib_check(iret,gribFunction,gribErrorMsg)
    call grib_get_int(igrib,'parameterCategory',parCat,iret)
    call grib_check(iret,gribFunction,gribErrorMsg)
    call grib_get_int(igrib,'parameterNumber',parNum,iret)
    call grib_check(iret,gribFunction,gribErrorMsg)
    call grib_get_int(igrib,'typeOfFirstFixedSurface',typSurf,iret)
    call grib_check(iret,gribFunction,gribErrorMsg)
    call grib_get_int(igrib,'level',valSurf,iret)
    call grib_check(iret,gribFunction,gribErrorMsg)
    call grib_get_int(igrib,'paramId',parId,iret)
    call grib_check(iret,gribFunction,gribErrorMsg)

    !print*,discipl,parCat,parNum,typSurf,valSurf

    !convert to grib1 identifiers
    isec1(6)=-1
    isec1(7)=-1
    isec1(8)=-1
    isec1(8)=valSurf     ! level
    if ((parCat.eq.0).and.(parNum.eq.0).and.(typSurf.eq.105)) then ! T
      isec1(6)=130         ! indicatorOfParameter
    elseif ((parCat.eq.2).and.(parNum.eq.2).and.(typSurf.eq.105)) then ! U
      isec1(6)=131         ! indicatorOfParameter
    elseif ((parCat.eq.2).and.(parNum.eq.3).and.(typSurf.eq.105)) then ! V
      isec1(6)=132         ! indicatorOfParameter
    elseif ((parCat.eq.1).and.(parNum.eq.0).and.(typSurf.eq.105)) then ! Q
      isec1(6)=133         ! indicatorOfParameter
      !ZHG FOR CLOUDS FROM GRIB
    elseif ((parCat.eq.1).and.(parNum.eq.83).and.(typSurf.eq.105)) then ! clwc
      isec1(6)=246         ! indicatorOfParameter
    elseif ((parCat.eq.1).and.(parNum.eq.84).and.(typSurf.eq.105)) then ! ciwc
      isec1(6)=247         ! indicatorOfParameter
      !ZHG end
      ! ESO qc(=clwc+ciwc)
    elseif ((parCat.eq.201).and.(parNum.eq.31).and.(typSurf.eq.105)) then ! qc
      isec1(6)=201031      ! indicatorOfParameter
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
    elseif ((parCat.eq.6).and.(parNum.eq.1) .or. parId .eq. 164) then ! CC
      isec1(6)=164         ! indicatorOfParameter
    elseif ((parCat.eq.1).and.(parNum.eq.9) .or. parId .eq. 142) then ! LSP
      isec1(6)=142         ! indicatorOfParameter
    elseif ((parCat.eq.1).and.(parNum.eq.10)) then ! CP
      isec1(6)=143         ! indicatorOfParameter
    elseif ((parCat.eq.0).and.(parNum.eq.11).and.(typSurf.eq.1)) then ! SHF
      isec1(6)=146         ! indicatorOfParameter
    elseif ((parCat.eq.4).and.(parNum.eq.9).and.(typSurf.eq.1)) then ! SR
      isec1(6)=176         ! indicatorOfParameter
    elseif ((parCat.eq.2).and.(parNum.eq.17) .or. parId .eq. 180) then ! EWSS
      isec1(6)=180         ! indicatorOfParameter
    elseif ((parCat.eq.2).and.(parNum.eq.18) .or. parId .eq. 181) then ! NSSS
      isec1(6)=181         ! indicatorOfParameter
    elseif ((parCat.eq.3).and.(parNum.eq.4)) then ! ORO
      isec1(6)=129         ! indicatorOfParameter
    elseif ((parCat.eq.3).and.(parNum.eq.7) .or. parId .eq. 160) then ! SDO
      isec1(6)=160         ! indicatorOfParameter
    elseif ((discipl.eq.2).and.(parCat.eq.0).and.(parNum.eq.0).and. &
         (typSurf.eq.1)) then ! LSM
      isec1(6)=172         ! indicatorOfParameter
    else
      print*,'***ERROR: undefined GRiB2 message found!',discipl, &
           parCat,parNum,typSurf
    endif
    if(parId .ne. isec1(6) .and. parId .ne. 77) then
      write(*,*) 'parId',parId, 'isec1(6)',isec1(6)
      !    stop
    endif

  endif

  CALL grib_get_int(igrib,'numberOfPointsAlongAParallel', &
       isec2(2),iret)
  ! nx=isec2(2)
  ! WRITE(*,*) nx,nxmax
  IF (isec2(2).GT.nxmax) THEN
    WRITE(*,*) 'FLEXPART error: Too many grid points in x direction.'
    WRITE(*,*) 'Reduce resolution of wind fields.'
    WRITE(*,*) 'Or change parameter settings in file ecmwf_mod.'
    WRITE(*,*) isec2(2),nxmax
  !    STOP
  ENDIF

  !get the size and data of the values array
  if (isec1(6).ne.-1) then
    call grib_get_real4_array(igrib,'values',zsec4,iret)
    call grib_check(iret,gribFunction,gribErrorMsg)
  endif

  if (ifield.eq.1) then

    !HSO  get the required fields from section 2 in a gribex compatible manner
    call grib_get_int(igrib,'numberOfPointsAlongAParallel', &
         isec2(2),iret)
    call grib_check(iret,gribFunction,gribErrorMsg)
    call grib_get_int(igrib,'numberOfPointsAlongAMeridian', &
         isec2(3),iret)
    call grib_check(iret,gribFunction,gribErrorMsg)
    call grib_get_real8(igrib,'longitudeOfFirstGridPointInDegrees', &
         xaux1in,iret)
    call grib_check(iret,gribFunction,gribErrorMsg)
    call grib_get_int(igrib,'numberOfVerticalCoordinateValues', &
         isec2(12),iret)
    call grib_check(iret,gribFunction,gribErrorMsg)

    !  get the size and data of the vertical coordinate array
    call grib_get_real4_array(igrib,'pv',zsec2,iret)
    call grib_check(iret,gribFunction,gribErrorMsg)

    nxfield=isec2(2)
    ny=isec2(3)
    nlev_ec=isec2(12)/2-1
  endif

  !HSO  get the second part of the grid dimensions only from GRiB1 messages
  if (isec1(6) .eq. 167 .and. (gotGrid.eq.0)) then
    call grib_get_real8(igrib,'longitudeOfLastGridPointInDegrees', &
         xaux2in,iret)
    call grib_check(iret,gribFunction,gribErrorMsg)
    call grib_get_real8(igrib,'latitudeOfLastGridPointInDegrees', &
         yaux1in,iret)
    call grib_check(iret,gribFunction,gribErrorMsg)
    call grib_get_real8(igrib,'latitudeOfFirstGridPointInDegrees', &
         yaux2in,iret)
    call grib_check(iret,gribFunction,gribErrorMsg)
    xaux1=xaux1in
    xaux2=xaux2in
    yaux1=yaux1in
    yaux2=yaux2in
    if (xaux1.gt.180.) xaux1=xaux1-360.0
    if (xaux2.gt.180.) xaux2=xaux2-360.0
    if (xaux1.lt.-180.) xaux1=xaux1+360.0
    if (xaux2.lt.-180.) xaux2=xaux2+360.0
    if (xaux2.lt.xaux1) xaux2=xaux2+360.0
    xlon0=xaux1
    ylat0=yaux1
    dx=(xaux2-xaux1)/real(nxfield-1)
    dy=(yaux2-yaux1)/real(ny-1)
    dxconst=180./(dx*r_earth*pi)
    dyconst=180./(dy*r_earth*pi)
    gotGrid=1
    ! Check whether fields are global
    ! If they contain the poles, specify polar stereographic map
    ! projections using the stlmbr- and stcm2p-calls
    !***********************************************************

    xauxa=abs(xaux2+dx-360.-xaux1)
    if (xauxa.lt.0.001) then
      nx=nxfield+1                 ! field is cyclic
      xglobal=.true.
      if (abs(nxshift).ge.nx) &
           stop 'nxshift in file par_mod is too large'
      xlon0=xlon0+real(nxshift)*dx
    else
      nx=nxfield
      xglobal=.false.
      if (nxshift.ne.0) &
           stop 'nxshift (par_mod) must be zero for non-global domain'
    endif
    nxmin1=nx-1
    nymin1=ny-1
    if (xlon0.gt.180.) xlon0=xlon0-360.
    xauxa=abs(yaux1+90.)
    if (xglobal.and.xauxa.lt.0.001) then
      sglobal=.true.               ! field contains south pole
      ! Enhance the map scale by factor 3 (*2=6) compared to north-south
      ! map scale
      sizesouth=6.*(switchsouth+90.)/dy
      call stlmbr(southpolemap,-90.,0.)
      call stcm2p(southpolemap,0.,0.,switchsouth,0.,sizesouth, &
           sizesouth,switchsouth,180.)
      switchsouthg=(switchsouth-ylat0)/dy
    else
      sglobal=.false.
      switchsouthg=999999.
    endif
    xauxa=abs(yaux2-90.)
    if (xglobal.and.xauxa.lt.0.001) then
      nglobal=.true.               ! field contains north pole
      ! Enhance the map scale by factor 3 (*2=6) compared to north-south
      ! map scale
      sizenorth=6.*(90.-switchnorth)/dy
      call stlmbr(northpolemap,90.,0.)
      call stcm2p(northpolemap,0.,0.,switchnorth,0.,sizenorth, &
           sizenorth,switchnorth,180.)
      switchnorthg=(switchnorth-ylat0)/dy
    else
      nglobal=.false.
      switchnorthg=999999.
    endif
    if (nxshift.lt.0) &
         stop 'nxshift (par_mod) must not be negative'
    if (nxshift.ge.nxfield) stop 'nxshift (par_mod) too large'
  endif ! gotGrid

  if (nx.gt.nxmax) then
    write(*,*) 'FLEXPART error: Too many grid points in x direction.'
    write(*,*) 'Reduce resolution of wind fields.'
    write(*,*) 'Or change parameter settings in file par_mod.'
    write(*,*) nx,nxmax
    stop
  endif

  if (ny.gt.nymax) then
    write(*,*) 'FLEXPART error: Too many grid points in y direction.'
    write(*,*) 'Reduce resolution of wind fields.'
    write(*,*) 'Or change parameter settings in file par_mod.'
    write(*,*) ny,nymax
    stop
  endif

  k=isec1(8)
  if(isec1(6).eq.131) iumax=max(iumax,nlev_ec-k+1)
  if(isec1(6).eq.135) iwmax=max(iwmax,nlev_ec-k+1)

  if(isec1(6).eq.129) then
    do jy=0,ny-1
      do ix=0,nxfield-1
        oro(ix,jy)=zsec4(nxfield*(ny-jy-1)+ix+1)/ga
      end do
    end do
  endif
  if(isec1(6).eq.172) then
    do jy=0,ny-1
      do ix=0,nxfield-1
        lsm(ix,jy)=zsec4(nxfield*(ny-jy-1)+ix+1)
      end do
    end do
  endif
  if(isec1(6).eq.160) then
    do jy=0,ny-1
      do ix=0,nxfield-1
        excessoro(ix,jy)=zsec4(nxfield*(ny-jy-1)+ix+1)
      end do
    end do
  endif

  call grib_release(igrib)
  goto 10                      !! READ NEXT LEVEL OR PARAMETER
  !
  ! CLOSING OF INPUT DATA FILE
  !

30 call grib_close_file(ifile)

  !error message if no fields found with correct first longitude in it
  if (gotGrid.eq.0) then
    print*,'***ERROR: input file needs to contain GRiB1 formatted'// &
         'messages'
    stop
  endif

  nuvz=iumax
  nwz =iwmax
  if(nuvz.eq.nlev_ec) nwz=nlev_ec+1

  if (nuvz+1.gt.nuvzmax) then
    write(*,*) 'FLEXPART error: Too many u,v grid points in z '// &
         'direction.'
    write(*,*) 'Reduce resolution of wind fields.'
    write(*,*) 'Or change parameter settings in file par_mod.'
    write(*,*) nuvz+1,nuvzmax
    stop
  endif

  if (nwz.gt.nwzmax) then
    write(*,*) 'FLEXPART error: Too many w grid points in z '// &
         'direction.'
    write(*,*) 'Reduce resolution of wind fields.'
    write(*,*) 'Or change parameter settings in file par_mod.'
    write(*,*) nwz,nwzmax
    stop
  endif

  ! If desired, shift all grids by nxshift grid cells
  !**************************************************

  if (xglobal) then
    call shift_field_0(oro,nxfield,ny)
    call shift_field_0(lsm,nxfield,ny)
    call shift_field_0(excessoro,nxfield,ny)
  endif

  ! Output of grid info
  !********************

  if (lroot) then
    write(*,'(a,2i7)') ' Vertical levels in ECMWF data: ', &
         nuvz+1,nwz
    write(*,*)
    write(*,'(a)') ' Mother domain:'
    write(*,'(a,f10.5,a,f10.5,a,f10.5)') '  Longitude range: ', &
         xlon0,' to ',xlon0+(nx-1)*dx,'   Grid distance: ',dx
    write(*,'(a,f10.5,a,f10.5,a,f10.5)') '  Latitude range : ', &
         ylat0,' to ',ylat0+(ny-1)*dy,'   Grid distance: ',dy
    write(*,*)
  end if

  ! CALCULATE VERTICAL DISCRETIZATION OF ECMWF MODEL
  ! PARAMETER akm,bkm DESCRIBE THE HYBRID "ETA" COORDINATE SYSTEM

  numskip=nlev_ec-nuvz  ! number of ecmwf model layers not used
  ! by trajectory model
  !do 8940 i=1,244
  !   write (*,*) 'zsec2:',i,ifield,zsec2(i),numskip
  !940  continue
  !   stop
  ! SEC SEC SEC
  ! for unknown reason zsec 1 to 10 is filled in this version
  ! compared to the old one
  ! SEC SEC SE
  do i=1,nwz
    j=numskip+i
    k=nlev_ec+1+numskip+i
    akm(nwz-i+1)=zsec2(j)
    !   write (*,*) 'ifield:',ifield,k,j,zsec2(10+j)
    bkm(nwz-i+1)=zsec2(k)
    wheight(nwz-i+1)=akm(nwz-i+1)/101325.+bkm(nwz-i+1) ! From FLEXTRA
  end do

  !
  ! CALCULATION OF AKZ, BKZ
  ! AKZ,BKZ: model discretization parameters at the center of each model
  !     layer
  !
  ! Assign the 10 m winds to an artificial model level with akz=0 and bkz=1.0,
  ! i.e. ground level
  !*****************************************************************************

  akz(1)=0.
  bkz(1)=1.0
  uvheight(1)=1.
  do i=1,nuvz
    uvheight(i+1)=0.5*(wheight(i+1)+wheight(i)) ! From FLEXTRA
    akz(i+1)=0.5*(akm(i+1)+akm(i))
    bkz(i+1)=0.5*(bkm(i+1)+bkm(i))
  end do
  ! exuvheight=wheight
  nuvz=nuvz+1

  ! NOTE: In FLEXPART versions up to 4.0, the number of model levels was doubled
  ! upon the transformation to z levels. In order to save computer memory, this is
  ! not done anymore in the standard version. However, this option can still be
  ! switched on by replacing the following lines with those below, that are
  ! currently commented out. For this, similar changes are necessary in
  ! verttransform.f and verttranform_nests.f
  !*****************************************************************************

  nz=nuvz
  if (nz.gt.nzmax) stop 'nzmax too small'
  do i=1,nuvz
    aknew(i)=akz(i)
    bknew(i)=bkz(i)
  end do

  ! Switch on following lines to use doubled vertical resolution
  !*************************************************************
  !nz=nuvz+nwz-1
  !if (nz.gt.nzmax) stop 'nzmax too small'
  !do 100 i=1,nwz
  !  aknew(2*(i-1)+1)=akm(i)
  !00     bknew(2*(i-1)+1)=bkm(i)
  !do 110 i=2,nuvz
  !  aknew(2*(i-1))=akz(i)
  !10     bknew(2*(i-1))=bkz(i)
  ! End doubled vertical resolution
  return

999 write(*,*)
  write(*,*) ' ###########################################'// &
       '###### '
  write(*,*) '       TRAJECTORY MODEL SUBROUTINE GRIDCHECK:'
  write(*,*) ' CAN NOT OPEN INPUT DATA FILE '//wfname(ifn)
  write(*,*) ' ###########################################'// &
       '###### '
  write(*,*)
  write(*,'(a)') '!!! PLEASE INSERT A NEW CD-ROM AND   !!!'
  write(*,'(a)') '!!! PRESS ANY KEY TO CONTINUE...     !!!'
  write(*,'(a)') '!!! ...OR TERMINATE FLEXPART PRESSING!!!'
  write(*,'(a)') '!!! THE "X" KEY...                   !!!'
  write(*,*)
  read(*,'(a)') opt
  if(opt.eq.'X') then
    stop
  else
    goto 5
  endif
end subroutine gridcheck_ecmwf

subroutine gridcheck_gfs

  !**********************************************************************
  !                                                                     *
  !             FLEXPART MODEL SUBROUTINE GRIDCHECK                     *
  !                                                                     *
  !**********************************************************************
  !                                                                     *
  !             AUTHOR:      G. WOTAWA                                  *
  !             DATE:        1997-08-06                                 *
  !             LAST UPDATE: 1997-10-10                                 *
  !                                                                     *
  !             Update:      1999-02-08, global fields allowed, A. Stohl*
  !             CHANGE: 17/11/2005, Caroline Forster, GFS data          *
  !             CHANGE: 11/01/2008, Harald Sodemann, GRIB1/2 input with *
  !                                 ECMWF grib_api                      *
  !             CHANGE: 03/12/2008, Harald Sodemann, update to f90 with *
  !                                 ECMWF grib_api                      *
  !                                                                     *
  !   Unified ECMWF and GFS builds                                      *
  !   Marian Harustak, 12.5.2017                                        *
  !     - Renamed routine from gridcheck to gridcheck_gfs               *
  !                                                                     *
  !**********************************************************************
  !                                                                     *
  ! DESCRIPTION:                                                        *
  !                                                                     *
  ! THIS SUBROUTINE DETERMINES THE GRID SPECIFICATIONS (LOWER LEFT      *
  ! LONGITUDE, LOWER LEFT LATITUDE, NUMBER OF GRID POINTS, GRID DIST-   *
  ! ANCE AND VERTICAL DISCRETIZATION OF THE ECMWF MODEL) FROM THE       *
  ! GRIB HEADER OF THE FIRST INPUT FILE. THE CONSISTANCY (NO CHANGES    *
  ! WITHIN ONE FLEXPART RUN) IS CHECKED IN THE ROUTINE "READWIND" AT    *
  ! ANY CALL.                                                           *
  !                                                                     *
  ! XLON0                geographical longitude of lower left gridpoint *
  ! YLAT0                geographical latitude of lower left gridpoint  *
  ! NX                   number of grid points x-direction              *
  ! NY                   number of grid points y-direction              *
  ! DX                   grid distance x-direction                      *
  ! DY                   grid distance y-direction                      *
  ! NUVZ                 number of grid points for horizontal wind      *
  !                      components in z direction                      *
  ! NWZ                  number of grid points for vertical wind        *
  !                      component in z direction                       *
  ! sizesouth, sizenorth give the map scale (i.e. number of virtual grid*
  !                      points of the polar stereographic grid):       *
  !                      used to check the CFL criterion                *
  ! UVHEIGHT(1)-         heights of gridpoints where u and v are        *
  ! UVHEIGHT(NUVZ)       given                                          *
  ! WHEIGHT(1)-          heights of gridpoints where w is given         *
  ! WHEIGHT(NWZ)                                                        *
  !                                                                     *
  !**********************************************************************

  use grib_api
  use cmapf_mod, only: stlmbr,stcm2p

  implicit none

  !HSO  parameters for grib_api
  integer :: ifile
  integer :: iret
  integer :: igrib
  real(kind=4) :: xaux1,xaux2,yaux1,yaux2
  real(kind=8) :: xaux1in,xaux2in,yaux1in,yaux2in
  integer :: gribVer,parCat,parNum,typSurf,valSurf,discipl
  !HSO  end
  integer :: ix,jy,i,ifn,ifield,j,k,iumax,iwmax,numskip
  real :: sizesouth,sizenorth,xauxa,pint
  real :: akm_usort(nwzmax)
  real,parameter :: eps=0.0001

  ! NCEP GFS
  real :: pres(nwzmax), help

  integer :: i179,i180,i181

  ! VARIABLES AND ARRAYS NEEDED FOR GRIB DECODING

  integer :: isec1(8),isec2(3)
  real(kind=4) :: zsec4(jpunp)
  character(len=1) :: opt

  !HSO  grib api error messages
  character(len=24) :: gribErrorMsg = 'Error reading grib file'
  character(len=20) :: gribFunction = 'gridcheckwind_gfs'
  !
  if (numbnests.ge.1) then
  write(*,*) ' ###########################################'
  write(*,*) ' FLEXPART ERROR SUBROUTINE GRIDCHECK:'
  write(*,*) ' NO NESTED WINDFIELDAS ALLOWED FOR GFS!      '
  write(*,*) ' ###########################################'
  stop
  endif

  iumax=0
  iwmax=0

  if(ideltas.gt.0) then
    ifn=1
  else
    ifn=numbwf
  endif
  !
  ! OPENING OF DATA FILE (GRIB CODE)
  !
5   call grib_open_file(ifile,path(3)(1:length(3)) &
         //trim(wfname(ifn)),'r',iret)
  if (iret.ne.GRIB_SUCCESS) then
    goto 999   ! ERROR DETECTED
  endif
  !turn on support for multi fields messages
  call grib_multi_support_on

  ifield=0
10   ifield=ifield+1
  !
  ! GET NEXT FIELDS
  !
  call grib_new_from_file(ifile,igrib,iret)
  if (iret.eq.GRIB_END_OF_FILE )  then
    goto 30    ! EOF DETECTED
  elseif (iret.ne.GRIB_SUCCESS) then
    goto 999   ! ERROR DETECTED
  endif

  !first see if we read GRIB1 or GRIB2
  call grib_get_int(igrib,'editionNumber',gribVer,iret)
  call grib_check(iret,gribFunction,gribErrorMsg)

  if (gribVer.eq.1) then ! GRIB Edition 1

  !read the grib1 identifiers
  call grib_get_int(igrib,'indicatorOfParameter',isec1(6),iret)
  call grib_check(iret,gribFunction,gribErrorMsg)
  call grib_get_int(igrib,'indicatorOfTypeOfLevel',isec1(7),iret)
  call grib_check(iret,gribFunction,gribErrorMsg)
  call grib_get_int(igrib,'level',isec1(8),iret)
  call grib_check(iret,gribFunction,gribErrorMsg)

  !get the size and data of the values array
  call grib_get_real4_array(igrib,'values',zsec4,iret)
  call grib_check(iret,gribFunction,gribErrorMsg)

  else ! GRIB Edition 2

  !read the grib2 identifiers
  call grib_get_int(igrib,'discipline',discipl,iret)
  call grib_check(iret,gribFunction,gribErrorMsg)
  call grib_get_int(igrib,'parameterCategory',parCat,iret)
  call grib_check(iret,gribFunction,gribErrorMsg)
  call grib_get_int(igrib,'parameterNumber',parNum,iret)
  call grib_check(iret,gribFunction,gribErrorMsg)
  call grib_get_int(igrib,'typeOfFirstFixedSurface',typSurf,iret)
  call grib_check(iret,gribFunction,gribErrorMsg)
  call grib_get_int(igrib,'scaledValueOfFirstFixedSurface', &
       valSurf,iret)
  call grib_check(iret,gribFunction,gribErrorMsg)

  !convert to grib1 identifiers
  isec1(6)=-1
  isec1(7)=-1
  isec1(8)=-1
  if ((parCat.eq.2).and.(parNum.eq.2).and.(typSurf.eq.100)) then ! U
    isec1(6)=33          ! indicatorOfParameter
    isec1(7)=100         ! indicatorOfTypeOfLevel
    isec1(8)=valSurf/100 ! level, convert to hPa
  elseif ((parCat.eq.3).and.(parNum.eq.5).and.(typSurf.eq.1)) then ! TOPO
    isec1(6)=7           ! indicatorOfParameter
    isec1(7)=1           ! indicatorOfTypeOfLevel
    isec1(8)=0
  elseif ((parCat.eq.0).and.(parNum.eq.0).and.(typSurf.eq.1) &
       .and.(discipl.eq.2)) then ! LSM
    isec1(6)=81          ! indicatorOfParameter
    isec1(7)=1           ! indicatorOfTypeOfLevel
    isec1(8)=0
  endif

  if (isec1(6).ne.-1) then
  !  get the size and data of the values array
    call grib_get_real4_array(igrib,'values',zsec4,iret)
    call grib_check(iret,gribFunction,gribErrorMsg)
  endif

  endif ! gribVer

  if(ifield.eq.1) then

  !get the required fields from section 2
  !store compatible to gribex input
  call grib_get_int(igrib,'numberOfPointsAlongAParallel', &
       isec2(2),iret)
  call grib_check(iret,gribFunction,gribErrorMsg)
  call grib_get_int(igrib,'numberOfPointsAlongAMeridian', &
       isec2(3),iret)
  call grib_check(iret,gribFunction,gribErrorMsg)
  call grib_get_real8(igrib,'longitudeOfFirstGridPointInDegrees', &
       xaux1in,iret)
  call grib_check(iret,gribFunction,gribErrorMsg)
  call grib_get_real8(igrib,'longitudeOfLastGridPointInDegrees', &
       xaux2in,iret)
  call grib_check(iret,gribFunction,gribErrorMsg)
  call grib_get_real8(igrib,'latitudeOfLastGridPointInDegrees', &
       yaux1in,iret)
  call grib_check(iret,gribFunction,gribErrorMsg)
  call grib_get_real8(igrib,'latitudeOfFirstGridPointInDegrees', &
       yaux2in,iret)
  call grib_check(iret,gribFunction,gribErrorMsg)

  ! Fix for flexpart.eu ticket #48
  if (xaux2in.lt.0) xaux2in = 359.0

  xaux1=xaux1in
  xaux2=xaux2in
  yaux1=yaux1in
  yaux2=yaux2in

  nxfield=isec2(2)
  ny=isec2(3)
  if((abs(xaux1).lt.eps).and.(xaux2.ge.359)) then ! NCEP DATA FROM 0 TO
    xaux1=-179.0                             ! 359 DEG EAST ->
    xaux2=-179.0+360.-360./real(nxfield)    ! TRANSFORMED TO -179
  endif                                      ! TO 180 DEG EAST
  if (xaux1.gt.180) xaux1=xaux1-360.0
  if (xaux2.gt.180) xaux2=xaux2-360.0
  if (xaux1.lt.-180) xaux1=xaux1+360.0
  if (xaux2.lt.-180) xaux2=xaux2+360.0
  if (xaux2.lt.xaux1) xaux2=xaux2+360.
  xlon0=xaux1
  ylat0=yaux1
  dx=(xaux2-xaux1)/real(nxfield-1)
  dy=(yaux2-yaux1)/real(ny-1)
  dxconst=180./(dx*r_earth*pi)
  dyconst=180./(dy*r_earth*pi)
  !HSO end edits


  ! Check whether fields are global
  ! If they contain the poles, specify polar stereographic map
  ! projections using the stlmbr- and stcm2p-calls
  !***********************************************************

    xauxa=abs(xaux2+dx-360.-xaux1)
    if (xauxa.lt.0.001) then
      nx=nxfield+1                 ! field is cyclic
      xglobal=.true.
      if (abs(nxshift).ge.nx) &
           stop 'nxshift in file par_mod is too large'
      xlon0=xlon0+real(nxshift)*dx
    else
      nx=nxfield
      xglobal=.false.
      if (nxshift.ne.0) &
           stop 'nxshift (par_mod) must be zero for non-global domain'
    endif
    nxmin1=nx-1
    nymin1=ny-1
    if (xlon0.gt.180.) xlon0=xlon0-360.
    xauxa=abs(yaux1+90.)
    if (xglobal.and.xauxa.lt.0.001) then
      sglobal=.true.               ! field contains south pole
  ! Enhance the map scale by factor 3 (*2=6) compared to north-south
  ! map scale
      sizesouth=6.*(switchsouth+90.)/dy
      call stlmbr(southpolemap,-90.,0.)
      call stcm2p(southpolemap,0.,0.,switchsouth,0.,sizesouth, &
           sizesouth,switchsouth,180.)
      switchsouthg=(switchsouth-ylat0)/dy
    else
      sglobal=.false.
      switchsouthg=999999.
    endif
    xauxa=abs(yaux2-90.)
    if (xglobal.and.xauxa.lt.0.001) then
      nglobal=.true.               ! field contains north pole
  ! Enhance the map scale by factor 3 (*2=6) compared to north-south
  ! map scale
      sizenorth=6.*(90.-switchnorth)/dy
      call stlmbr(northpolemap,90.,0.)
      call stcm2p(northpolemap,0.,0.,switchnorth,0.,sizenorth, &
           sizenorth,switchnorth,180.)
      switchnorthg=(switchnorth-ylat0)/dy
    else
      nglobal=.false.
      switchnorthg=999999.
    endif
  endif ! ifield.eq.1

  if (nxshift.lt.0) stop 'nxshift (par_mod) must not be negative'
  if (nxshift.ge.nxfield) stop 'nxshift (par_mod) too large'

  ! NCEP ISOBARIC LEVELS
  !*********************

  if((isec1(6).eq.33).and.(isec1(7).eq.100)) then ! check for U wind
    iumax=iumax+1
    pres(iumax)=real(isec1(8))*100.0
  endif


  i179=nint(179./dx)
  if (dx.lt.0.7) then
    i180=nint(180./dx)+1    ! 0.5 deg data
  else
    i180=nint(179./dx)+1    ! 1 deg data
  endif
  i181=i180+1


  ! NCEP TERRAIN
  !*************

  if((isec1(6).eq.007).and.(isec1(7).eq.001)) then
    do jy=0,ny-1
      do ix=0,nxfield-1
        help=zsec4(nxfield*(ny-jy-1)+ix+1)
        if(ix.le.i180) then
          oro(i179+ix,jy)=help
          excessoro(i179+ix,jy)=0.0 ! ISOBARIC SURFACES: SUBGRID TERRAIN DISREGARDED
        else
          oro(ix-i181,jy)=help
          excessoro(ix-i181,jy)=0.0 ! ISOBARIC SURFACES: SUBGRID TERRAIN DISREGARDED
        endif
      end do
    end do
  endif

  ! NCEP LAND SEA MASK
  !*******************

  if((isec1(6).eq.081).and.(isec1(7).eq.001)) then
    do jy=0,ny-1
      do ix=0,nxfield-1
        help=zsec4(nxfield*(ny-jy-1)+ix+1)
        if(ix.le.i180) then
          lsm(i179+ix,jy)=help
        else
          lsm(ix-i181,jy)=help
        endif
      end do
    end do
  endif

  call grib_release(igrib)

  goto 10                      !! READ NEXT LEVEL OR PARAMETER
  !
  ! CLOSING OF INPUT DATA FILE
  !

  ! HSO
30   continue
  call grib_close_file(ifile)
  ! HSO end edits

  nuvz=iumax
  nwz =iumax
  nlev_ec=iumax

  if (nx.gt.nxmax) then
   write(*,*) 'FLEXPART error: Too many grid points in x direction.'
    write(*,*) 'Reduce resolution of wind fields.'
    write(*,*) 'Or change parameter settings in file par_mod.'
    write(*,*) nx,nxmax
    stop
  endif

  if (ny.gt.nymax) then
   write(*,*) 'FLEXPART error: Too many grid points in y direction.'
    write(*,*) 'Reduce resolution of wind fields.'
    write(*,*) 'Or change parameter settings in file par_mod.'
    write(*,*) ny,nymax
    stop
  endif

  if (nuvz.gt.nuvzmax) then
    write(*,*) 'FLEXPART error: Too many u,v grid points in z '// &
         'direction.'
    write(*,*) 'Reduce resolution of wind fields.'
    write(*,*) 'Or change parameter settings in file par_mod.'
    write(*,*) nuvz,nuvzmax
    stop
  endif

  if (nwz.gt.nwzmax) then
    write(*,*) 'FLEXPART error: Too many w grid points in z '// &
         'direction.'
    write(*,*) 'Reduce resolution of wind fields.'
    write(*,*) 'Or change parameter settings in file par_mod.'
    write(*,*) nwz,nwzmax
    stop
  endif

  ! If desired, shift all grids by nxshift grid cells
  !**************************************************

  if (xglobal) then
    call shift_field_0(oro,nxfield,ny)
    call shift_field_0(lsm,nxfield,ny)
    call shift_field_0(excessoro,nxfield,ny)
  endif

  ! Output of grid info
  !********************

  if (lroot) then
    write(*,*)
    write(*,*)
    write(*,'(a,2i7)') 'Vertical levels in NCEP data: ', &
         nuvz,nwz
    write(*,*)
    write(*,'(a)') 'Mother domain:'
    write(*,'(a,f10.2,a1,f10.2,a,f10.2)') '  Longitude range: ', &
         xlon0,' to ',xlon0+(nx-1)*dx,'   Grid distance: ',dx
    write(*,'(a,f10.2,a1,f10.2,a,f10.2)') '  Latitude range : ', &
         ylat0,' to ',ylat0+(ny-1)*dy,'   Grid distance: ',dy
    write(*,*)
  end if

  ! CALCULATE VERTICAL DISCRETIZATION OF ECMWF MODEL
  ! PARAMETER akm,bkm DESCRIBE THE HYBRID "ETA" COORDINATE SYSTEM

  numskip=nlev_ec-nuvz  ! number of ecmwf model layers not used
                        ! by trajectory model
  do i=1,nwz
    j=numskip+i
    k=nlev_ec+1+numskip+i
    akm_usort(nwz-i+1)=pres(nwz-i+1)
    bkm(nwz-i+1)=0.0
  end do

  !******************************
  ! change Sabine Eckhardt: akm should always be in descending order ... readwind adapted!
  !******************************
      do i=1,nwz
         if (akm_usort(1).gt.akm_usort(2)) then
            akm(i)=akm_usort(i)
         else
            akm(i)=akm_usort(nwz-i+1)
         endif
      end do

  !
  ! CALCULATION OF AKZ, BKZ
  ! AKZ,BKZ: model discretization parameters at the center of each model
  !     layer
  !
  ! Assign the 10 m winds to an artificial model level with akz=0 and bkz=1.0,
  ! i.e. ground level
  !*****************************************************************************

  do i=1,nuvz
     akz(i)=akm(i)
     bkz(i)=bkm(i)
  end do

  ! NOTE: In FLEXPART versions up to 4.0, the number of model levels was doubled
  ! upon the transformation to z levels. In order to save computer memory, this is
  ! not done anymore in the standard version. However, this option can still be
  ! switched on by replacing the following lines with those below, that are
  ! currently commented out. For this, similar changes are necessary in
  ! verttransform.f and verttranform_nests.f
  !*****************************************************************************

  nz=nuvz
  if (nz.gt.nzmax) stop 'nzmax too small'
  do i=1,nuvz
    aknew(i)=akz(i)
    bknew(i)=bkz(i)
  end do

  ! Switch on following lines to use doubled vertical resolution
  !*************************************************************
  !nz=nuvz+nwz-1
  !if (nz.gt.nzmax) stop 'nzmax too small'
  !do 100 i=1,nwz
  !  aknew(2*(i-1)+1)=akm(i)
  !00     bknew(2*(i-1)+1)=bkm(i)
  !do 110 i=2,nuvz
  !  aknew(2*(i-1))=akz(i)
  !10     bknew(2*(i-1))=bkz(i)
  ! End doubled vertical resolution
  return

999   write(*,*)
  write(*,*) ' ###########################################'// &
       '###### '
  write(*,*) '       TRAJECTORY MODEL SUBROUTINE GRIDCHECK:'
  write(*,*) ' CAN NOT OPEN INPUT DATA FILE '//wfname(ifn)
  write(*,*) ' ###########################################'// &
       '###### '
  write(*,*)
  write(*,'(a)') '!!! PLEASE INSERT A NEW CD-ROM AND   !!!'
  write(*,'(a)') '!!! PRESS ANY KEY TO CONTINUE...     !!!'
  write(*,'(a)') '!!! ...OR TERMINATE FLEXPART PRESSING!!!'
  write(*,'(a)') '!!! THE "X" KEY...                   !!!'
  write(*,*)
  read(*,'(a)') opt
  if(opt.eq.'X') then
    stop
  else
    goto 5
  endif

end subroutine gridcheck_gfs

subroutine gridcheck_nests

  !*****************************************************************************
  !                                                                            *
  !     This routine checks the grid specification for the nested model        *
  !     domains. It is similar to subroutine gridcheck, which checks the       *
  !     mother domain.                                                         *
  !                                                                            *
  !     Authors: A. Stohl, G. Wotawa                                           *
  !                                                                            *
  !     8 February 1999                                                        *
  !                                                                            *
  !*****************************************************************************
  !  CHANGE: 11/01/2008, Harald Sodemann, GRIB1/2 input with ECMWF grib_api    *
  !  CHANGE: 03/12/2008, Harald Sodemann, change to f90 grib_api               *
  !*****************************************************************************

  use grib_api

  implicit none

  !HSO  parameters for grib_api
  integer :: ifile
  integer :: iret
  integer :: igrib
  integer :: gribVer,parCat,parNum,typSurf,valSurf,discipl
  integer :: parID !added by mc for making it consistent with new gridcheck.f90
  integer :: gotGrib
  !HSO  end
  integer :: i,j,k,l,ifn,ifield,iumax,iwmax,numskip,nlev_ecn
  integer :: nuvzn,nwzn
  real :: akmn(nwzmax),bkmn(nwzmax),akzn(nuvzmax),bkzn(nuvzmax)
  real(kind=4) :: xaux1,xaux2,yaux1,yaux2
  real(kind=8) :: xaux1in,xaux2in,yaux1in,yaux2in
  real :: conversion_factor !added by mc to make it consistent with new gridchek.f90

  ! VARIABLES AND ARRAYS NEEDED FOR GRIB DECODING

  ! dimension of isec2 at least (22+n), where n is the number of parallels or
  ! meridians in a quasi-regular (reduced) Gaussian or lat/long grid

  ! dimension of zsec2 at least (10+nn), where nn is the number of vertical
  ! coordinate parameters

  integer :: isec1(56),isec2(22+nxmaxn+nymaxn)
  real(kind=4) :: zsec2(60+2*nuvzmax),zsec4(jpunp)

  !HSO  grib api error messages
  character(len=24) :: gribErrorMsg = 'Error reading grib file'
  character(len=20) :: gribFunction = 'gridcheck_nests'

  xresoln(0)=1.       ! resolution enhancement for mother grid
  yresoln(0)=1.       ! resolution enhancement for mother grid

  ! Loop about all nesting levels
  !******************************

  do l=1,numbnests

    iumax=0
    iwmax=0

    if(ideltas.gt.0) then
      ifn=1
    else
      ifn=numbwf
    endif
  !
  ! OPENING OF DATA FILE (GRIB CODE)
  !
  ifile=0
  igrib=0
  iret=0

5   call grib_open_file(ifile,path(numpath+2*(l-1)+1) &
         (1:length(numpath+2*(l-1)+1))//trim(wfnamen(l,ifn)),'r',iret)
  if (iret.ne.GRIB_SUCCESS) then
    goto 999   ! ERROR DETECTED
  endif
  !turn on support for multi fields messages
  !call grib_multi_support_on

  gotGrib=0
  ifield=0
10   ifield=ifield+1

  !
  ! GET NEXT FIELDS
  !
  call grib_new_from_file(ifile,igrib,iret)
  if (iret.eq.GRIB_END_OF_FILE)  then
    goto 30    ! EOF DETECTED
  elseif (iret.ne.GRIB_SUCCESS) then
    goto 999   ! ERROR DETECTED
  endif

  !first see if we read GRIB1 or GRIB2
  call grib_get_int(igrib,'editionNumber',gribVer,iret)
  call grib_check(iret,gribFunction,gribErrorMsg)

  if (gribVer.eq.1) then ! GRIB Edition 1

  !print*,'GRiB Edition 1'
  !read the grib2 identifiers
  call grib_get_int(igrib,'indicatorOfParameter',isec1(6),iret)
  call grib_check(iret,gribFunction,gribErrorMsg)
  call grib_get_int(igrib,'level',isec1(8),iret)
  call grib_check(iret,gribFunction,gribErrorMsg)

  !change code for etadot to code for omega
  if (isec1(6).eq.77) then
    isec1(6)=135
  endif

  !print*,isec1(6),isec1(8)

  else

  !print*,'GRiB Edition 2'
  !read the grib2 identifiers
  call grib_get_int(igrib,'discipline',discipl,iret)
  call grib_check(iret,gribFunction,gribErrorMsg)
  call grib_get_int(igrib,'parameterCategory',parCat,iret)
  call grib_check(iret,gribFunction,gribErrorMsg)
  call grib_get_int(igrib,'parameterNumber',parNum,iret)
  call grib_check(iret,gribFunction,gribErrorMsg)
  call grib_get_int(igrib,'typeOfFirstFixedSurface',typSurf,iret)
  call grib_check(iret,gribFunction,gribErrorMsg)
  call grib_get_int(igrib,'level',valSurf,iret)
  call grib_check(iret,gribFunction,gribErrorMsg)
  call grib_get_int(igrib,'paramId',parId,iret) !added by mc to make it consisitent with new grid_check.f90
  call grib_check(iret,gribFunction,gribErrorMsg) !added by mc to make it consisitent with new  grid_check.f90

  !print*,discipl,parCat,parNum,typSurf,valSurf

  !convert to grib1 identifiers
  isec1(6)=-1
  isec1(7)=-1
  isec1(8)=-1
  isec1(8)=valSurf     ! level
  if ((parCat.eq.0).and.(parNum.eq.0).and.(typSurf.eq.105)) then ! T
    isec1(6)=130         ! indicatorOfParameter
  elseif ((parCat.eq.2).and.(parNum.eq.2).and.(typSurf.eq.105)) then ! U
    isec1(6)=131         ! indicatorOfParameter
  elseif ((parCat.eq.2).and.(parNum.eq.3).and.(typSurf.eq.105)) then ! V
    isec1(6)=132         ! indicatorOfParameter
  elseif ((parCat.eq.1).and.(parNum.eq.0).and.(typSurf.eq.105)) then ! Q
    isec1(6)=133         ! indicatorOfParameter
  elseif ((parCat.eq.1).and.(parNum.eq.83).and.(typSurf.eq.105)) then ! clwc
    isec1(6)=246         ! indicatorOfParameter
  elseif ((parCat.eq.1).and.(parNum.eq.84).and.(typSurf.eq.105)) then ! ciwc
    isec1(6)=247         ! indicatorOfParameter
  !ZHG end
  ! ESO qc(=clwc+ciwc)
  elseif ((parCat.eq.201).and.(parNum.eq.31).and.(typSurf.eq.105)) then ! qc
    isec1(6)=201031      ! indicatorOfParameter
  elseif ((parCat.eq.3).and.(parNum.eq.0).and.(typSurf.eq.1)) then !SP
    isec1(6)=134         ! indicatorOfParameter
  elseif ((parCat.eq.2).and.(parNum.eq.32)) then ! W, actually eta dot
    isec1(6)=135         ! indicatorOfParameter
  elseif ((parCat.eq.128).and.(parNum.eq.77)) then ! W, actually eta dot !added bymc to make it consistent with new gridcheck.f90
    isec1(6)=135         ! indicatorOfParameter    !
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
  elseif ((parCat.eq.6).and.(parNum.eq.1) .or. parId .eq. 164) then ! CC !added by mc to make it consistent with new gridchek.f90
    isec1(6)=164         ! indicatorOfParameter
 elseif ((parCat.eq.1).and.(parNum.eq.9) .or. parId .eq. 142) then ! LSP !added by mc to make it consistent with new gridchek.f90
    isec1(6)=142         ! indicatorOfParameter
  elseif ((parCat.eq.1).and.(parNum.eq.10)) then ! CP
    isec1(6)=143         ! indicatorOfParameter
  elseif ((parCat.eq.0).and.(parNum.eq.11).and.(typSurf.eq.1)) then ! SHF
    isec1(6)=146         ! indicatorOfParameter
  elseif ((parCat.eq.4).and.(parNum.eq.9).and.(typSurf.eq.1)) then ! SR
    isec1(6)=176         ! indicatorOfParameter
  elseif ((parCat.eq.2).and.(parNum.eq.17) .or. parId .eq. 180) then ! EWSS !added by mc to make it consistent with new gridchek.f90
    isec1(6)=180         ! indicatorOfParameter
  elseif ((parCat.eq.2).and.(parNum.eq.18) .or. parId .eq. 181) then ! NSSS !added by mc to make it consistent with new gridchek.f90
    isec1(6)=181         ! indicatorOfParameter
  elseif ((parCat.eq.3).and.(parNum.eq.4)) then ! ORO
    isec1(6)=129         ! indicatorOfParameter
  elseif ((parCat.eq.3).and.(parNum.eq.7) .or. parId .eq. 160) then ! SDO !added by mc to make it consistent with new gridchek.f90
    isec1(6)=160         ! indicatorOfParameter
  elseif ((discipl.eq.2).and.(parCat.eq.0).and.(parNum.eq.0).and. &
       (typSurf.eq.1)) then ! LSM
    isec1(6)=172         ! indicatorOfParameter
  else
    print*,'***ERROR: undefined GRiB2 message found!',discipl, &
         parCat,parNum,typSurf
  endif
  if(parId .ne. isec1(6) .and. parId .ne. 77) then !added by mc to make it consistent with new gridchek.f90
    write(*,*) 'parId',parId, 'isec1(6)',isec1(6)
  !    stop
  endif

  endif

  !get the size and data of the values array
  if (isec1(6).ne.-1) then
    call grib_get_real4_array(igrib,'values',zsec4,iret)
    call grib_check(iret,gribFunction,gribErrorMsg)
  endif

  !HSO  get the required fields from section 2 in a gribex compatible manner
  if (ifield.eq.1) then
    call grib_get_int(igrib,'numberOfPointsAlongAParallel', &
         isec2(2),iret)
    call grib_check(iret,gribFunction,gribErrorMsg)
    call grib_get_int(igrib,'numberOfPointsAlongAMeridian', &
         isec2(3),iret)
    call grib_check(iret,gribFunction,gribErrorMsg)
    call grib_get_int(igrib,'numberOfVerticalCoordinateValues', &
         isec2(12),iret)
    call grib_check(iret,gribFunction,gribErrorMsg)
  !HSO    get the size and data of the vertical coordinate array
    call grib_get_real4_array(igrib,'pv',zsec2,iret)
    call grib_check(iret,gribFunction,gribErrorMsg)

    nxn(l)=isec2(2)
    nyn(l)=isec2(3)
    nlev_ecn=isec2(12)/2-1
  endif ! ifield

  if (nxn(l).gt.nxmaxn) then
  write(*,*) 'FLEXPART error: Too many grid points in x direction.'
  write(*,*) 'Reduce resolution of wind fields (file GRIDSPEC)'
  write(*,*) 'for nesting level ',l
  write(*,*) 'Or change parameter settings in file par_mod.'
  write(*,*) nxn(l),nxmaxn
  stop
  endif

  if (nyn(l).gt.nymaxn) then
  write(*,*) 'FLEXPART error: Too many grid points in y direction.'
  write(*,*) 'Reduce resolution of wind fields (file GRIDSPEC)'
  write(*,*) 'for nesting level ',l
  write(*,*) 'Or change parameter settings in file par_mod.'
  write(*,*) nyn(l),nymaxn
  stop
  endif

  !HSO  get the second part of the grid dimensions only from GRiB1 messages
  if (isec1(6) .eq. 167 .and. (gotGrib.eq.0)) then !added by mc to make it consistent with new gridchek.f90 note that gotGrid must be changed in gotGrib!!
    call grib_get_real8(igrib,'longitudeOfFirstGridPointInDegrees', & !comment by mc: note that this was in the (if (ifield.eq.1) ..end above in gridchek.f90 see line 257
         xaux1in,iret)
    call grib_check(iret,gribFunction,gribErrorMsg)
    call grib_get_real8(igrib,'longitudeOfLastGridPointInDegrees', &
         xaux2in,iret)
    call grib_check(iret,gribFunction,gribErrorMsg)
    call grib_get_real8(igrib,'latitudeOfLastGridPointInDegrees', &
         yaux1in,iret)
    call grib_check(iret,gribFunction,gribErrorMsg)
    call grib_get_real8(igrib,'latitudeOfFirstGridPointInDegrees', &
         yaux2in,iret)
    call grib_check(iret,gribFunction,gribErrorMsg)
    xaux1=xaux1in
    xaux2=xaux2in
    yaux1=yaux1in
    yaux2=yaux2in
    if(xaux1.gt.180.) xaux1=xaux1-360.0
    if(xaux2.gt.180.) xaux2=xaux2-360.0
    if(xaux1.lt.-180.) xaux1=xaux1+360.0
    if(xaux2.lt.-180.) xaux2=xaux2+360.0
    if (xaux2.lt.xaux1) xaux2=xaux2+360.0
    xlon0n(l)=xaux1
    ylat0n(l)=yaux1
    dxn(l)=(xaux2-xaux1)/real(nxn(l)-1)
    dyn(l)=(yaux2-yaux1)/real(nyn(l)-1)
    gotGrib=1 !commetn by mc note tahthere gotGRIB is used instead of gotGrid!!!
  endif ! ifield.eq.1

  k=isec1(8)
  if(isec1(6).eq.131) iumax=max(iumax,nlev_ec-k+1)
  if(isec1(6).eq.135) iwmax=max(iwmax,nlev_ec-k+1)

  if(isec1(6).eq.129) then
    do j=0,nyn(l)-1
      do i=0,nxn(l)-1
        oron(i,j,l)=zsec4(nxn(l)*(nyn(l)-j-1)+i+1)/ga
      end do
    end do
  endif
  if(isec1(6).eq.172) then
    do j=0,nyn(l)-1
      do i=0,nxn(l)-1
        lsmn(i,j,l)=zsec4(nxn(l)*(nyn(l)-j-1)+i+1)/ga
      end do
    end do
  endif
  if(isec1(6).eq.160) then
    do j=0,nyn(l)-1
      do i=0,nxn(l)-1
        excessoron(i,j,l)=zsec4(nxn(l)*(nyn(l)-j-1)+i+1)/ga
      end do
    end do
  endif

  call grib_release(igrib)
  goto 10                 !! READ NEXT LEVEL OR PARAMETER
  !
  ! CLOSING OF INPUT DATA FILE
  !

30   call grib_close_file(ifile)

  !error message if no fields found with correct first longitude in it
  if (gotGrib.eq.0) then
    print*,'***ERROR: input file needs to contain GRiB1 formatted'// &
         'messages'
    stop
  endif

  nuvzn=iumax
  nwzn=iwmax
  if(nuvzn.eq.nlev_ec) nwzn=nlev_ecn+1

  if ((nuvzn.gt.nuvzmax).or.(nwzn.gt.nwzmax)) then
  write(*,*) 'FLEXPART error: Nested wind fields have too many'// &
       'vertical levels.'
  write(*,*) 'Problem was encountered for nesting level ',l
  stop
  endif


  ! Output of grid info
  !********************

  write(*,'(a,i2,a)') ' Nested domain ',l,':'
  write(*,'(a,f10.5,a,f10.5,a,f10.5)') '  Longitude range: ', &
       xlon0n(l),' to ',xlon0n(l)+(nxn(l)-1)*dxn(l), &
       '   Grid distance: ',dxn(l)
  write(*,'(a,f10.5,a,f10.5,a,f10.5)') '  Latitude range : ', &
       ylat0n(l),' to ',ylat0n(l)+(nyn(l)-1)*dyn(l), &
       '   Grid distance: ',dyn(l)
  write(*,*)

  ! Determine, how much the resolutions in the nests are enhanced as
  ! compared to the mother grid
  !*****************************************************************

  xresoln(l)=dx/dxn(l)
  yresoln(l)=dy/dyn(l)

  ! Determine the mother grid coordinates of the corner points of the
  ! nested grids
  ! Convert first to geographical coordinates, then to grid coordinates
  !********************************************************************

  xaux1=xlon0n(l)
  xaux2=xlon0n(l)+real(nxn(l)-1)*dxn(l)
  yaux1=ylat0n(l)
  yaux2=ylat0n(l)+real(nyn(l)-1)*dyn(l)

  xln(l)=(xaux1-xlon0)/dx
  xrn(l)=(xaux2-xlon0)/dx
  yln(l)=(yaux1-ylat0)/dy
  yrn(l)=(yaux2-ylat0)/dy


  if ((xln(l).lt.0.).or.(yln(l).lt.0.).or. &
       (xrn(l).gt.real(nxmin1)).or.(yrn(l).gt.real(nymin1))) then
    write(*,*) 'Nested domain does not fit into mother domain'
    write(*,*) 'For global mother domain fields, you can shift'
    write(*,*) 'shift the mother domain into x-direction'
    write(*,*) 'by setting nxshift (file par_mod) to a'
    write(*,*) 'positive value. Execution is terminated.'
    stop
  endif


  ! CALCULATE VERTICAL DISCRETIZATION OF ECMWF MODEL
  ! PARAMETER akm,bkm DESCRIBE THE HYBRID "ETA" COORDINATE SYSTEM

  numskip=nlev_ecn-nuvzn ! number of ecmwf model layers not used by FLEXPART
  do i=1,nwzn
    j=numskip+i
    k=nlev_ecn+1+numskip+i
    akmn(nwzn-i+1)=zsec2(j)
    bkmn(nwzn-i+1)=zsec2(k)
  end do

  !
  ! CALCULATION OF AKZ, BKZ
  ! AKZ,BKZ: model discretization parameters at the center of each model
  !     layer
  !
  ! Assign the 10 m winds to an artificial model level with akz=0 and bkz=1.0,
  ! i.e. ground level
  !*****************************************************************************

    akzn(1)=0.
    bkzn(1)=1.0
    do i=1,nuvzn
      akzn(i+1)=0.5*(akmn(i+1)+akmn(i))
      bkzn(i+1)=0.5*(bkmn(i+1)+bkmn(i))
    end do
    nuvzn=nuvzn+1

  ! Check, whether the heights of the model levels of the nested
  ! wind fields are consistent with those of the mother domain.
  ! If not, terminate model run.
  !*************************************************************

    do i=1,nuvz
      if ((akzn(i).ne.akz(i)).or.(bkzn(i).ne.bkz(i))) then
  write(*,*) 'FLEXPART error: The wind fields of nesting level',l
  write(*,*) 'are not consistent with the mother domain:'
  write(*,*) 'Differences in vertical levels detected.'
        stop
      endif
    end do

    do i=1,nwz
      if ((akmn(i).ne.akm(i)).or.(bkmn(i).ne.bkm(i))) then
  write(*,*) 'FLEXPART error: The wind fields of nesting level',l
  write(*,*) 'are not consistent with the mother domain:'
  write(*,*) 'Differences in vertical levels detected.'
        stop
      endif
    end do

  end do

  return

999   write(*,*)
  write(*,*) ' ###########################################'// &
       '###### '
  write(*,*) '                FLEXPART SUBROUTINE GRIDCHECK:'
  write(*,*) ' CAN NOT OPEN INPUT DATA FILE '//wfnamen(l,ifn)
  write(*,*) ' FOR NESTING LEVEL ',k
  write(*,*) ' ###########################################'// &
       '###### '
  stop

end subroutine gridcheck_nests

subroutine pbl_profile(ps,td2m,zml1,t2m,tml1,u10m,uml1,stress,hf)

  !********************************************************************
  !                                                                   *
  !                    G. WOTAWA, 1995-07-07                          *
  !                                                                   *
  !********************************************************************
  !                                                                   *
  ! DESCRIPTION: CALCULATION OF FRICTION VELOCITY AND SURFACE SENS-   *
  !              IBLE HEAT FLUX USING THE PROFILE METHOD (BERKOVICZ   *
  !              AND PRAHM, 1982)                                     *
  !                                                                   *
  ! Output now is surface stress instead of ustar                     *
  !                                                                   *
  !                                                                   *
  !********************************************************************
  !                                                                   *
  ! INPUT:                                                            *
  !                                                                   *
  !                                                                   *
  ! ps      surface pressure(Pa)                                      *
  ! td2m    two metre dew point(K)                                    *
  ! zml1    heigth of first model level (m)                           *
  ! t2m     two metre temperature (K)                                 *
  ! tml1    temperature first model level (K)                         *
  ! u10m    ten metre wind speed (ms-1)                               *
  ! uml1    wind speed first model level (ms-1)                       *
  !                                                                   *
  !********************************************************************
  !                                                                   *
  ! OUTPUT:                                                           *
  !                                                                   *
  ! stress  surface stress (i.e., friction velocity (ms-1) squared    *
  !                         multiplied with air density)              *
  ! hf      surface sensible heat flux (Wm-2)                         *
  !                                                                   *
  !********************************************************************
  ! ustar   friction velocity (ms-1)                                  *
  ! maxiter maximum number of iterations                              *
  !********************************************************************

  use par_mod
  use qvsat_mod
  use stability_correction

  implicit none

  integer :: iter
  real :: ps,td2m,rhoa,zml1,t2m,tml1,u10m,uml1,ustar,hf
  real :: al,alold,aldiff,tmean,crit
  real :: deltau,deltat,thetastar,e,tv,stress
  integer,parameter :: maxiter=10
  real,parameter    :: r1=0.74

  e=ew(td2m,ps)               ! vapor pressure
  tv=t2m*(1.+0.378*e/ps)   ! virtual temperature
  rhoa=ps/(r_air*tv)       ! air density

  deltau=uml1-u10m         !! Wind Speed difference between
                           !! Model level 1 and 10 m

  if(deltau.le.0.001) then    !! Monin-Obukhov Theory not
    al=9999.               !! applicable --> Set dummy values
    ustar=0.01
    stress=ustar*ustar*rhoa
    hf=0.0
    return
  endif
  deltat=tml1-t2m+0.0098*(zml1-2.)  !! Potential temperature difference
                                    !! between model level 1 and 10 m

  if(abs(deltat).le.0.03) then    !! Neutral conditions
    hf=0.0
    al=9999.
    ustar=(vonkarman*deltau)/ &
         (log(zml1/10.)-psim(zml1,al)+psim(10.,al))
    stress=ustar*ustar*rhoa
    return
  endif

  tmean=0.5*(t2m+tml1)
  crit=(0.0219*tmean*(zml1-2.0)*deltau**2)/ &
       (deltat*(zml1-10.0)**2)
  if((deltat.gt.0).and.(crit.le.1.)) then
                                    !! Successive approximation will
    al=50.                          !! not converge
    ustar=(vonkarman*deltau)/ &
         (log(zml1/10.)-psim(zml1,al)+psim(10.,al))
    thetastar=(vonkarman*deltat/r1)/ &
         (log(zml1/2.)-psih(zml1,al)+psih(2.,al))
    hf=rhoa*cpa*ustar*thetastar
    stress=ustar*ustar*rhoa
    return
  endif

  al=9999.                 ! Start iteration assuming neutral conditions
  do iter=1,maxiter
    alold=al
    ustar=(vonkarman*deltau)/ &
         (log(zml1/10.)-psim(zml1,al)+psim(10.,al))
    thetastar=(vonkarman*deltat/r1)/ &
         (log(zml1/2.)-psih(zml1,al)+psih(2.,al))
    al=(tmean*ustar**2)/(ga*vonkarman*thetastar)
    aldiff=abs((al-alold)/alold)
    if(aldiff.lt.0.01) exit  !! Successive approximation successful
  end do
  hf=rhoa*cpa*ustar*thetastar
  if(al.gt.9999.) al=9999.
  if(al.lt.-9999.) al=-9999.

  stress=ustar*ustar*rhoa
end subroutine pbl_profile

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

  implicit none

  !  include 'grib_api.h'

  !HSO  parameters for grib_api
  integer :: ifile
  integer :: iret
  integer, dimension(:), allocatable   :: igrib
  integer :: nfield, ii, arsize
  integer :: gribVer,parCat,parNum,typSurf,valSurf,discipl,parId
  integer :: gotGrid
  ! HSO  end

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
  call grib_open_file(ifile,path(3)(1:length(3)) &
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

subroutine readwind_gfs(indj,n,uuh,vvh,wwh)

  !***********************************************************************
  !*                                                                     *
  !*             TRAJECTORY MODEL SUBROUTINE READWIND                    *
  !*                                                                     *
  !***********************************************************************
  !*                                                                     *
  !*             AUTHOR:      G. WOTAWA                                  *
  !*             DATE:        1997-08-05                                 *
  !*             LAST UPDATE: 2000-10-17, Andreas Stohl                  *
  !*             CHANGE: 01/02/2001, Bernd C. Krueger, Variables tth and *
  !*                     qvh (on eta coordinates) in common block        *
  !*             CHANGE: 16/11/2005, Caroline Forster, GFS data          *
  !*             CHANGE: 11/01/2008, Harald Sodemann, Input of GRIB1/2   *
  !*                     data with the ECMWF grib_api library            *
  !*             CHANGE: 03/12/2008, Harald Sodemann, update to f90 with *
  !*                                 ECMWF grib_api                      *
  !                                                                      *
  !   Unified ECMWF and GFS builds                                       *
  !   Marian Harustak, 12.5.2017                                         *
  !     - Renamed routine from readwind to readwind_gfs                  *
  !*                                                                     *
  !***********************************************************************
  !*                                                                     *
  !* DESCRIPTION:                                                        *
  !*                                                                     *
  !* READING OF ECMWF METEOROLOGICAL FIELDS FROM INPUT DATA FILES. THE   *
  !* INPUT DATA FILES ARE EXPECTED TO BE AVAILABLE IN GRIB CODE          *
  !*                                                                     *
  !* INPUT:                                                              *
  !* indj               indicates number of the wind field to be read in *
  !* n                  temporal index for meteorological fields (1 to 3)*
  !*                                                                     *
  !* IMPORTANT VARIABLES FROM COMMON BLOCK:                              *
  !*                                                                     *
  !* wfname             File name of data to be read in                  *
  !* nx,ny,nuvz,nwz     expected field dimensions                        *
  !* nlev_ec            number of vertical levels ecmwf model            *
  !* uu,vv,ww           wind fields                                      *
  !* tt,qv              temperature and specific humidity                *
  !* ps                 surface pressure                                 *
  !*                                                                     *
  !***********************************************************************

  use grib_api
  use qvsat_mod

  implicit none

  !HSO  new parameters for grib_api
  integer :: ifile
  integer :: iret
  integer :: igrib
  integer :: gribVer,parCat,parNum,typSurf,valSurf,discipl
  !HSO end edits
  real :: uuh(0:nxmax-1,0:nymax-1,nuvzmax)
  real :: vvh(0:nxmax-1,0:nymax-1,nuvzmax)
  real :: wwh(0:nxmax-1,0:nymax-1,nwzmax)
  integer :: ii,indj,i,j,k,n,levdiff2,ifield,iumax,iwmax

  ! NCEP
  integer :: numpt,numpu,numpv,numpw,numprh,numpclwch
  real :: help, temp
  real :: elev
  real :: ulev1(0:nxmax-1,0:nymax-1),vlev1(0:nxmax-1,0:nymax-1)
  real :: tlev1(0:nxmax-1,0:nymax-1)
  real :: qvh2(0:nxmax-1,0:nymax-1)

  integer :: i179,i180,i181

  ! VARIABLES AND ARRAYS NEEDED FOR GRIB DECODING
  !HSO kept isec1, isec2 and zsec4 for consistency with gribex GRIB input

  integer :: isec1(8),isec2(3)
  real(kind=4) :: zsec4(jpunp)
  real(kind=4) :: xaux,yaux,xaux0,yaux0
  real(kind=8) :: xauxin,yauxin
  real,parameter :: eps=1.e-4
  real(kind=4) :: ewss(0:nxmax-1,0:nymax-1),nsss(0:nxmax-1,0:nymax-1)
  real :: plev1,hlev1,ff10m,fflev1

  logical :: hflswitch,strswitch

  !HSO  for grib api error messages
  character(len=24) :: gribErrorMsg = 'Error reading grib file'
  character(len=20) :: gribFunction = 'readwind_gfs'
  character(len=20) :: shortname


  hflswitch=.false.
  strswitch=.false.
  levdiff2=nlev_ec-nwz+1
  iumax=0
  iwmax=0


  ! OPENING OF DATA FILE (GRIB CODE)

  !HSO
  call grib_open_file(ifile,path(3)(1:length(3)) &
         //trim(wfname(indj)),'r',iret)
  if (iret.ne.GRIB_SUCCESS) then
    goto 888   ! ERROR DETECTED
  endif
  !turn on support for multi fields messages
  call grib_multi_support_on

  numpt=0
  numpu=0
  numpv=0
  numpw=0
  numprh=0
  numpclwch=0
  ifield=0
10   ifield=ifield+1
  !
  ! GET NEXT FIELDS
  !
  call grib_new_from_file(ifile,igrib,iret)
  if (iret.eq.GRIB_END_OF_FILE)  then
    goto 50    ! EOF DETECTED
  elseif (iret.ne.GRIB_SUCCESS) then
    goto 888   ! ERROR DETECTED
  endif

  !first see if we read GRIB1 or GRIB2
  call grib_get_int(igrib,'editionNumber',gribVer,iret)
  !  call grib_check(iret,gribFunction,gribErrorMsg)

  if (gribVer.eq.1) then ! GRIB Edition 1

  !read the grib1 identifiers
  call grib_get_int(igrib,'indicatorOfParameter',isec1(6),iret)
  !  call grib_check(iret,gribFunction,gribErrorMsg)
  call grib_get_int(igrib,'indicatorOfTypeOfLevel',isec1(7),iret)
  !  call grib_check(iret,gribFunction,gribErrorMsg)
  call grib_get_int(igrib,'level',isec1(8),iret)
  !  call grib_check(iret,gribFunction,gribErrorMsg)

  else ! GRIB Edition 2

  !read the grib2 identifiers
  call grib_get_string(igrib,'shortName',shortname,iret)

  call grib_get_int(igrib,'discipline',discipl,iret)
  !  call grib_check(iret,gribFunction,gribErrorMsg)
  call grib_get_int(igrib,'parameterCategory',parCat,iret)
  !  call grib_check(iret,gribFunction,gribErrorMsg)
  call grib_get_int(igrib,'parameterNumber',parNum,iret)
  !  call grib_check(iret,gribFunction,gribErrorMsg)
  call grib_get_int(igrib,'typeOfFirstFixedSurface',typSurf,iret)
  !  call grib_check(iret,gribFunction,gribErrorMsg)
  call grib_get_int(igrib,'scaledValueOfFirstFixedSurface', &
       valSurf,iret)
  !  call grib_check(iret,gribFunction,gribErrorMsg)
  
  !  write(*,*) 'Field: ',ifield,parCat,parNum,typSurf,shortname
  !convert to grib1 identifiers
  isec1(6)=-1
  isec1(7)=-1
  isec1(8)=-1
  if ((parCat.eq.0).and.(parNum.eq.0).and.(typSurf.eq.100)) then ! T
    isec1(6)=11          ! indicatorOfParameter
    isec1(7)=100         ! indicatorOfTypeOfLevel
    isec1(8)=valSurf/100 ! level, convert to hPa
  elseif ((parCat.eq.2).and.(parNum.eq.2).and.(typSurf.eq.100)) then ! U
    isec1(6)=33          ! indicatorOfParameter
    isec1(7)=100         ! indicatorOfTypeOfLevel
    isec1(8)=valSurf/100 ! level, convert to hPa
  elseif ((parCat.eq.2).and.(parNum.eq.3).and.(typSurf.eq.100)) then ! V
    isec1(6)=34          ! indicatorOfParameter
    isec1(7)=100         ! indicatorOfTypeOfLevel
    isec1(8)=valSurf/100 ! level, convert to hPa
  elseif ((parCat.eq.2).and.(parNum.eq.8).and.(typSurf.eq.100)) then ! W
    isec1(6)=39          ! indicatorOfParameter
    isec1(7)=100         ! indicatorOfTypeOfLevel
    isec1(8)=valSurf/100 ! level, convert to hPa
  elseif ((parCat.eq.1).and.(parNum.eq.1).and.(typSurf.eq.100)) then ! RH
    isec1(6)=52          ! indicatorOfParameter
    isec1(7)=100         ! indicatorOfTypeOfLevel
    isec1(8)=valSurf/100 ! level, convert to hPa
  elseif ((parCat.eq.1).and.(parNum.eq.1).and.(typSurf.eq.103)) then ! RH2
    isec1(6)=52          ! indicatorOfParameter
    isec1(7)=105         ! indicatorOfTypeOfLevel
    isec1(8)=2
  elseif ((parCat.eq.0).and.(parNum.eq.0).and.(typSurf.eq.103)) then ! T2
    isec1(6)=11          ! indicatorOfParameter
    isec1(7)=105         ! indicatorOfTypeOfLevel
    isec1(8)=2
  elseif ((parCat.eq.2).and.(parNum.eq.2).and.(typSurf.eq.103)) then ! U10
    isec1(6)=33          ! indicatorOfParameter
    isec1(7)=105         ! indicatorOfTypeOfLevel
    isec1(8)=10
  elseif ((parCat.eq.2).and.(parNum.eq.3).and.(typSurf.eq.103)) then ! V10
    isec1(6)=34          ! indicatorOfParameter
    isec1(7)=105         ! indicatorOfTypeOfLevel
    isec1(8)=10
  elseif ((parCat.eq.1).and.(parNum.eq.22).and.(typSurf.eq.100)) then ! CLWMR Cloud Mixing Ratio [kg/kg]:
    isec1(6)=153         ! indicatorOfParameter
    isec1(7)=100         ! indicatorOfTypeOfLevel
    isec1(8)=valSurf/100 ! level, convert to hPa
  elseif ((parCat.eq.3).and.(parNum.eq.1).and.(typSurf.eq.101)) then ! SLP
    isec1(6)=2           ! indicatorOfParameter
    isec1(7)=102         ! indicatorOfTypeOfLevel
    isec1(8)=0
  elseif ((parCat.eq.3).and.(parNum.eq.0).and.(typSurf.eq.1)) then ! SP
    isec1(6)=1           ! indicatorOfParameter
    isec1(7)=1           ! indicatorOfTypeOfLevel
    isec1(8)=0
  elseif ((parCat.eq.1).and.(parNum.eq.13).and.(typSurf.eq.1)) then ! SNOW
    isec1(6)=66          ! indicatorOfParameter
    isec1(7)=1           ! indicatorOfTypeOfLevel
    isec1(8)=0
  elseif ((parCat.eq.0).and.(parNum.eq.0).and.(typSurf.eq.104)) then ! T sigma 0
    isec1(6)=11          ! indicatorOfParameter
    isec1(7)=107         ! indicatorOfTypeOfLevel
    isec1(8)=0.995       ! lowest sigma level
  elseif ((parCat.eq.2).and.(parNum.eq.2).and.(typSurf.eq.104)) then ! U sigma 0
    isec1(6)=33          ! indicatorOfParameter
    isec1(7)=107         ! indicatorOfTypeOfLevel
    isec1(8)=0.995       ! lowest sigma level
  elseif ((parCat.eq.2).and.(parNum.eq.3).and.(typSurf.eq.104)) then ! V sigma 0
    isec1(6)=34          ! indicatorOfParameter
    isec1(7)=107         ! indicatorOfTypeOfLevel
    isec1(8)=0.995       ! lowest sigma level
  elseif ((parCat.eq.3).and.(parNum.eq.5).and.(typSurf.eq.1)) then ! TOPO
    isec1(6)=7           ! indicatorOfParameter
    isec1(7)=1           ! indicatorOfTypeOfLevel
    isec1(8)=0
  elseif ((parCat.eq.0).and.(parNum.eq.0).and.(typSurf.eq.1) &
       .and.(discipl.eq.2)) then ! LSM
    isec1(6)=81          ! indicatorOfParameter
    isec1(7)=1           ! indicatorOfTypeOfLevel
    isec1(8)=0
  elseif ((parCat.eq.3).and.(parNum.eq.196).and.(typSurf.eq.1)) then ! BLH
    isec1(6)=221         ! indicatorOfParameter
    isec1(7)=1           ! indicatorOfTypeOfLevel
    isec1(8)=0
  elseif ((parCat.eq.1).and.(parNum.eq.7).and.(typSurf.eq.1)) then ! LSP/TP
    isec1(6)=62          ! indicatorOfParameter
    isec1(7)=1           ! indicatorOfTypeOfLevel
    isec1(8)=0
  elseif ((parCat.eq.1).and.(parNum.eq.196).and.(typSurf.eq.1)) then ! CP
    isec1(6)=63          ! indicatorOfParameter
    isec1(7)=1           ! indicatorOfTypeOfLevel
    isec1(8)=0
  endif

  endif ! gribVer

  if (isec1(6).ne.-1) then
  !  get the size and data of the values array
    call grib_get_real4_array(igrib,'values',zsec4,iret)
  !    call grib_check(iret,gribFunction,gribErrorMsg)
  endif

  if(ifield.eq.1) then

  !get the required fields from section 2
  !store compatible to gribex input
  call grib_get_int(igrib,'numberOfPointsAlongAParallel', &
       isec2(2),iret)
  !  call grib_check(iret,gribFunction,gribErrorMsg)
  call grib_get_int(igrib,'numberOfPointsAlongAMeridian', &
       isec2(3),iret)
  !  call grib_check(iret,gribFunction,gribErrorMsg)
  call grib_get_real8(igrib,'longitudeOfFirstGridPointInDegrees', &
       xauxin,iret)
  !  call grib_check(iret,gribFunction,gribErrorMsg)
  call grib_get_real8(igrib,'latitudeOfLastGridPointInDegrees', &
       yauxin,iret)
  !  call grib_check(iret,gribFunction,gribErrorMsg)
  xaux=xauxin+real(nxshift)*dx
  yaux=yauxin

  ! CHECK GRID SPECIFICATIONS

    if(isec2(2).ne.nxfield) stop 'READWIND: NX NOT CONSISTENT'
    if(isec2(3).ne.ny) stop 'READWIND: NY NOT CONSISTENT'
    if(xaux.eq.0.) xaux=-179.0     ! NCEP DATA
    xaux0=xlon0
    yaux0=ylat0
    if(xaux.lt.0.) xaux=xaux+360.
    if(yaux.lt.0.) yaux=yaux+360.
    if(xaux0.lt.0.) xaux0=xaux0+360.
    if(yaux0.lt.0.) yaux0=yaux0+360.
    if(abs(xaux-xaux0).gt.eps) &
         stop 'READWIND: LOWER LEFT LONGITUDE NOT CONSISTENT'
    if(abs(yaux-yaux0).gt.eps) &
         stop 'READWIND: LOWER LEFT LATITUDE NOT CONSISTENT'
  endif
  !HSO end of edits

  i179=nint(179./dx)
  if (dx.lt.0.7) then
    i180=nint(180./dx)+1    ! 0.5 deg data
  else
    i180=nint(179./dx)+1    ! 1 deg data
  endif
  i181=i180+1

  if (isec1(6).ne.-1) then

  do j=0,nymin1
    do i=0,nxfield-1
      if((isec1(6).eq.011).and.(isec1(7).eq.100)) then
  ! TEMPERATURE
         if((i.eq.0).and.(j.eq.0)) then
            do ii=1,nuvz
              if ((isec1(8)*100.0).eq.akz(ii)) numpt=ii
            end do
        endif
        help=zsec4(nxfield*(ny-j-1)+i+1)
        if(i.le.i180) then
          tth(i179+i,j,numpt,n)=help
        else
          tth(i-i181,j,numpt,n)=help
        endif
      endif
      if((isec1(6).eq.033).and.(isec1(7).eq.100)) then
  ! U VELOCITY
         if((i.eq.0).and.(j.eq.0)) then
            do ii=1,nuvz
              if ((isec1(8)*100.0).eq.akz(ii)) numpu=ii
            end do
        endif
        help=zsec4(nxfield*(ny-j-1)+i+1)
        if(i.le.i180) then
          uuh(i179+i,j,numpu)=help
        else
          uuh(i-i181,j,numpu)=help
        endif
      endif
      if((isec1(6).eq.034).and.(isec1(7).eq.100)) then
  ! V VELOCITY
         if((i.eq.0).and.(j.eq.0)) then
            do ii=1,nuvz
              if ((isec1(8)*100.0).eq.akz(ii)) numpv=ii
            end do
        endif
        help=zsec4(nxfield*(ny-j-1)+i+1)
        if(i.le.i180) then
          vvh(i179+i,j,numpv)=help
        else
          vvh(i-i181,j,numpv)=help
        endif
      endif
      if((isec1(6).eq.052).and.(isec1(7).eq.100)) then
  ! RELATIVE HUMIDITY -> CONVERT TO SPECIFIC HUMIDITY LATER
         if((i.eq.0).and.(j.eq.0)) then
            do ii=1,nuvz
              if ((isec1(8)*100.0).eq.akz(ii)) numprh=ii
            end do
        endif
        help=zsec4(nxfield*(ny-j-1)+i+1)
        if(i.le.i180) then
          qvh(i179+i,j,numprh,n)=help
        else
          qvh(i-i181,j,numprh,n)=help
        endif
      endif
      if((isec1(6).eq.001).and.(isec1(7).eq.001)) then
  ! SURFACE PRESSURE
        help=zsec4(nxfield*(ny-j-1)+i+1)
        if(i.le.i180) then
          ps(i179+i,j,1,n)=help
        else
          ps(i-i181,j,1,n)=help
        endif
      endif
      if((isec1(6).eq.039).and.(isec1(7).eq.100)) then
  ! W VELOCITY
         if((i.eq.0).and.(j.eq.0)) then
            do ii=1,nuvz
              if ((isec1(8)*100.0).eq.akz(ii)) numpw=ii
            end do
        endif
        help=zsec4(nxfield*(ny-j-1)+i+1)
        if(i.le.i180) then
          wwh(i179+i,j,numpw)=help
        else
          wwh(i-i181,j,numpw)=help
        endif
      endif
      if((isec1(6).eq.066).and.(isec1(7).eq.001)) then
  ! SNOW DEPTH
        help=zsec4(nxfield*(ny-j-1)+i+1)
        if(i.le.i180) then
          sd(i179+i,j,1,n)=help
        else
          sd(i-i181,j,1,n)=help
        endif
      endif
      if((isec1(6).eq.002).and.(isec1(7).eq.102)) then
  ! MEAN SEA LEVEL PRESSURE
        help=zsec4(nxfield*(ny-j-1)+i+1)
        if(i.le.i180) then
          msl(i179+i,j,1,n)=help
        else
          msl(i-i181,j,1,n)=help
        endif
      endif
      if((isec1(6).eq.071).and.(isec1(7).eq.244)) then
  ! TOTAL CLOUD COVER
        help=zsec4(nxfield*(ny-j-1)+i+1)
        if(i.le.i180) then
          tcc(i179+i,j,1,n)=help
        else
          tcc(i-i181,j,1,n)=help
        endif
      endif
      if((isec1(6).eq.033).and.(isec1(7).eq.105).and. &
           (isec1(8).eq.10)) then
  ! 10 M U VELOCITY
        help=zsec4(nxfield*(ny-j-1)+i+1)
        if(i.le.i180) then
          u10(i179+i,j,1,n)=help
        else
          u10(i-i181,j,1,n)=help
        endif
      endif
      if((isec1(6).eq.034).and.(isec1(7).eq.105).and. &
           (isec1(8).eq.10)) then
  ! 10 M V VELOCITY
        help=zsec4(nxfield*(ny-j-1)+i+1)
        if(i.le.i180) then
          v10(i179+i,j,1,n)=help
        else
          v10(i-i181,j,1,n)=help
        endif
      endif
      if((isec1(6).eq.011).and.(isec1(7).eq.105).and. &
           (isec1(8).eq.02)) then
  ! 2 M TEMPERATURE
        help=zsec4(nxfield*(ny-j-1)+i+1)
        if(i.le.i180) then
          tt2(i179+i,j,1,n)=help
        else
          tt2(i-i181,j,1,n)=help
        endif
      endif
      if((isec1(6).eq.017).and.(isec1(7).eq.105).and. &
           (isec1(8).eq.02)) then
  ! 2 M DEW POINT TEMPERATURE
        help=zsec4(nxfield*(ny-j-1)+i+1)
        if(i.le.i180) then
          td2(i179+i,j,1,n)=help
        else
          td2(i-i181,j,1,n)=help
        endif
      endif
      if((isec1(6).eq.062).and.(isec1(7).eq.001)) then
  ! LARGE SCALE PREC.
        help=zsec4(nxfield*(ny-j-1)+i+1)
        if(i.le.i180) then
          lsprec(i179+i,j,1,n)=help
        else
          lsprec(i-i181,j,1,n)=help
        endif
      endif
      if((isec1(6).eq.063).and.(isec1(7).eq.001)) then
  ! CONVECTIVE PREC.
        help=zsec4(nxfield*(ny-j-1)+i+1)
        if(i.le.i180) then
          convprec(i179+i,j,1,n)=help
        else
          convprec(i-i181,j,1,n)=help
        endif
      endif
      if((isec1(6).eq.007).and.(isec1(7).eq.001)) then
  ! TOPOGRAPHY
        help=zsec4(nxfield*(ny-j-1)+i+1)
        if(i.le.i180) then
          oro(i179+i,j)=help
          excessoro(i179+i,j)=0.0 ! ISOBARIC SURFACES: SUBGRID TERRAIN DISREGARDED
        else
          oro(i-i181,j)=help
          excessoro(i-i181,j)=0.0 ! ISOBARIC SURFACES: SUBGRID TERRAIN DISREGARDED
        endif
      endif
      if((isec1(6).eq.081).and.(isec1(7).eq.001)) then
  ! LAND SEA MASK
        help=zsec4(nxfield*(ny-j-1)+i+1)
        if(i.le.i180) then
          lsm(i179+i,j)=help
        else
          lsm(i-i181,j)=help
        endif
      endif
      if((isec1(6).eq.221).and.(isec1(7).eq.001)) then
  ! MIXING HEIGHT
        help=zsec4(nxfield*(ny-j-1)+i+1)
        if(i.le.i180) then
          hmix(i179+i,j,1,n)=help
        else
          hmix(i-i181,j,1,n)=help
        endif
      endif
      if((isec1(6).eq.052).and.(isec1(7).eq.105).and. &
           (isec1(8).eq.02)) then
  ! 2 M RELATIVE HUMIDITY
        help=zsec4(nxfield*(ny-j-1)+i+1)
        if(i.le.i180) then
          qvh2(i179+i,j)=help
        else
          qvh2(i-i181,j)=help
        endif
      endif
      if((isec1(6).eq.011).and.(isec1(7).eq.107)) then
  ! TEMPERATURE LOWEST SIGMA LEVEL
        help=zsec4(nxfield*(ny-j-1)+i+1)
        if(i.le.i180) then
          tlev1(i179+i,j)=help
        else
          tlev1(i-i181,j)=help
        endif
      endif
      if((isec1(6).eq.033).and.(isec1(7).eq.107)) then
  ! U VELOCITY LOWEST SIGMA LEVEL
        help=zsec4(nxfield*(ny-j-1)+i+1)
        if(i.le.i180) then
          ulev1(i179+i,j)=help
        else
          ulev1(i-i181,j)=help
        endif
      endif
      if((isec1(6).eq.034).and.(isec1(7).eq.107)) then
  ! V VELOCITY LOWEST SIGMA LEVEL
        help=zsec4(nxfield*(ny-j-1)+i+1)
        if(i.le.i180) then
          vlev1(i179+i,j)=help
        else
          vlev1(i-i181,j)=help
        endif
      endif
  ! SEC & IP 12/2018 read GFS clouds
      if((isec1(6).eq.153).and.(isec1(7).eq.100)) then  !! CLWCR  Cloud liquid water content [kg/kg] 
         if((i.eq.0).and.(j.eq.0)) then
            do ii=1,nuvz
              if ((isec1(8)*100.0).eq.akz(ii)) numpclwch=ii
            end do
        endif
        help=zsec4(nxfield*(ny-j-1)+i+1)
        if(i.le.i180) then
          clwch(i179+i,j,numpclwch,n)=help
        else
          clwch(i-i181,j,numpclwch,n)=help
        endif
        readclouds=.true.
        sumclouds=.true.
  !        readclouds=.false.
  !       sumclouds=.false.
      endif


    end do
  end do

  endif

  if((isec1(6).eq.33).and.(isec1(7).eq.100)) then
  ! NCEP ISOBARIC LEVELS
    iumax=iumax+1
  endif

  call grib_release(igrib)
  goto 10                      !! READ NEXT LEVEL OR PARAMETER
  !
  ! CLOSING OF INPUT DATA FILE
  !

  !HSO close grib file
50   continue
  call grib_close_file(ifile)

  ! SENS. HEAT FLUX
  sshf(:,:,1,n)=0.0     ! not available from gfs.tccz.pgrbfxx files
  hflswitch=.false.    ! Heat flux not available
  ! SOLAR RADIATIVE FLUXES
  ssr(:,:,1,n)=0.0      ! not available from gfs.tccz.pgrbfxx files
  ! EW SURFACE STRESS
  ewss=0.0         ! not available from gfs.tccz.pgrbfxx files
  ! NS SURFACE STRESS
  nsss=0.0         ! not available from gfs.tccz.pgrbfxx files
  strswitch=.false.    ! stress not available

  ! CONVERT TP TO LSP (GRIB2 only)
  if (gribVer.eq.2) then
    do j=0,nymin1
    do i=0,nxfield-1
     if(i.le.i180) then
     if (convprec(i179+i,j,1,n).lt.lsprec(i179+i,j,1,n)) then ! neg precip would occur
         lsprec(i179+i,j,1,n)= &
              lsprec(i179+i,j,1,n)-convprec(i179+i,j,1,n)
     else
         lsprec(i179+i,j,1,n)=0
     endif
     else
     if (convprec(i-i181,j,1,n).lt.lsprec(i-i181,j,1,n)) then
          lsprec(i-i181,j,1,n)= &
               lsprec(i-i181,j,1,n)-convprec(i-i181,j,1,n)
     else
          lsprec(i-i181,j,1,n)=0
     endif
     endif
    enddo
    enddo
  endif
  !HSO end edits


  ! TRANSFORM RH TO SPECIFIC HUMIDITY

  do j=0,ny-1
    do i=0,nxfield-1
      do k=1,nuvz
        help=qvh(i,j,k,n)
        temp=tth(i,j,k,n)
        plev1=akm(k)+bkm(k)*ps(i,j,1,n)
        elev=ew(temp,plev1)*help/100.0
        qvh(i,j,k,n)=xmwml*(elev/(plev1-((1.0-xmwml)*elev)))
      end do
    end do
  end do

  ! CALCULATE 2 M DEW POINT FROM 2 M RELATIVE HUMIDITY
  ! USING BOLTON'S (1980) FORMULA
  ! BECAUSE td2 IS NOT AVAILABLE FROM NCEP GFS DATA

  do j=0,ny-1
    do i=0,nxfield-1
        help=qvh2(i,j)
        temp=tt2(i,j,1,n)
        plev1=akm(k)+bkm(k)*ps(i,j,1,n)
        elev=ew(temp,plev1)/100.*help/100.   !vapour pressure in hPa
        td2(i,j,1,n)=243.5/(17.67/log(elev/6.112)-1)+273.
        if (help.le.0.) td2(i,j,1,n)=tt2(i,j,1,n)
    end do
  end do

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
    call shift_field_0(ulev1,nxfield,ny)
    call shift_field_0(vlev1,nxfield,ny)
    call shift_field_0(tlev1,nxfield,ny)
    call shift_field_0(qvh2,nxfield,ny)
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
    call shift_field(hmix,nxfield,ny,1,1,2,n)
    call shift_field(tth,nxfield,ny,nuvzmax,nuvz,2,n)
    call shift_field(qvh,nxfield,ny,nuvzmax,nuvz,2,n)
    call shift_field(uuh,nxfield,ny,nuvzmax,nuvz,1,1)
    call shift_field(vvh,nxfield,ny,nuvzmax,nuvz,1,1)
    call shift_field(wwh,nxfield,ny,nwzmax,nwz,1,1)
  ! IP & SEC adding GFS Clouds 20181205
    call shift_field(clwch,nxfield,ny,nuvzmax,nuvz,2,n)
  endif

  do i=0,nxmin1
    do j=0,nymin1
  ! Convert precip. from mm/s -> mm/hour
      convprec(i,j,1,n)=convprec(i,j,1,n)*3600.
      lsprec(i,j,1,n)=lsprec(i,j,1,n)*3600.
      surfstr(i,j,1,n)=sqrt(ewss(i,j)**2+nsss(i,j)**2)
    end do
  end do

  if ((.not.hflswitch).or.(.not.strswitch)) then
  !  write(*,*) 'WARNING: No flux data contained in GRIB file ',
  !    +  wfname(indj)

  ! CALCULATE USTAR AND SSHF USING THE PROFILE METHOD
  !***************************************************************************

    do i=0,nxmin1
      do j=0,nymin1
        hlev1=30.0                     ! HEIGHT OF FIRST MODEL SIGMA LAYER
        ff10m= sqrt(u10(i,j,1,n)**2+v10(i,j,1,n)**2)
        fflev1=sqrt(ulev1(i,j)**2+vlev1(i,j)**2)
        call pbl_profile(ps(i,j,1,n),td2(i,j,1,n),hlev1, &
             tt2(i,j,1,n),tlev1(i,j),ff10m,fflev1, &
             surfstr(i,j,1,n),sshf(i,j,1,n))
        if(sshf(i,j,1,n).gt.200.) sshf(i,j,1,n)=200.
        if(sshf(i,j,1,n).lt.-400.) sshf(i,j,1,n)=-400.
      end do
    end do
  endif

  if(iumax.ne.nuvz) stop 'READWIND: NUVZ NOT CONSISTENT'
  if(iumax.ne.nwz)    stop 'READWIND: NWZ NOT CONSISTENT'

  return
888   write(*,*) ' #### FLEXPART MODEL ERROR! WINDFIELD         #### '
  write(*,*) ' #### ',wfname(indj),'                    #### '
  write(*,*) ' #### IS NOT GRIB FORMAT !!!                  #### '
  stop 'Execution terminated'
999   write(*,*) ' #### FLEXPART MODEL ERROR! WINDFIELD         #### '
  write(*,*) ' #### ',wfname(indj),'                    #### '
  write(*,*) ' #### CANNOT BE OPENED !!!                    #### '
  stop 'Execution terminated'

end subroutine readwind_gfs

subroutine readwind_nests(indj,n,uuhn,vvhn,wwhn)
  !                           i   i  o    o    o
  !*****************************************************************************
  !                                                                            *
  !     This routine reads the wind fields for the nested model domains.       *
  !     It is similar to subroutine readwind, which reads the mother domain.   *
  !                                                                            *
  !     Authors: A. Stohl, G. Wotawa                                           *
  !                                                                            *
  !     8 February 1999                                                        *
  !                                                                            *
  !     Last update: 17 October 2000, A. Stohl                                 *
  !                                                                            *
  !*****************************************************************************
  !  Changes, Bernd C. Krueger, Feb. 2001:                                     *
  !        Variables tthn and qvhn (on eta coordinates) in common block        *
  !  CHANGE: 11/01/2008, Harald Sodemann, GRIB1/2 input with ECMWF grib_api    *
  !  CHANGE: 03/12/2008, Harald Sodemann, update to f90 with ECMWF grib_api    *
  !*****************************************************************************

  use grib_api

  implicit none

  !HSO  parameters for grib_api
  integer :: ifile
  integer :: iret
  integer :: igrib
  integer :: gribVer,parCat,parNum,typSurf,valSurf,discipl
  integer :: parId !!added by mc for making it consistent with new readwind.f90
  integer :: gotGrid
  !HSO  end

  real :: uuhn(0:nxmaxn-1,0:nymaxn-1,nuvzmax,maxnests)
  real :: vvhn(0:nxmaxn-1,0:nymaxn-1,nuvzmax,maxnests)
  real :: wwhn(0:nxmaxn-1,0:nymaxn-1,nwzmax,maxnests)
  integer :: indj,i,j,k,n,levdiff2,ifield,iumax,iwmax,l

  ! VARIABLES AND ARRAYS NEEDED FOR GRIB DECODING

  ! dimension of isec2 at least (22+n), where n is the number of parallels or
  ! meridians in a quasi-regular (reduced) Gaussian or lat/long grid

  ! dimension of zsec2 at least (10+nn), where nn is the number of vertical
  ! coordinate parameters

  integer :: isec1(56),isec2(22+nxmaxn+nymaxn)
  real(kind=4) :: zsec4(jpunp)
  real(kind=4) :: xaux,yaux
  real(kind=8) :: xauxin,yauxin
  real,parameter :: eps=1.e-4
  real :: ewss(0:nxmaxn-1,0:nymaxn-1),nsss(0:nxmaxn-1,0:nymaxn-1)
  real :: plev1,pmean,tv,fu,hlev1,ff10m,fflev1
  real :: conversion_factor !added by mc to make it consistent with new gridchek.f90

  logical :: hflswitch,strswitch

  !HSO  grib api error messages
  character(len=24) :: gribErrorMsg = 'Error reading grib file'
  character(len=20) :: gribFunction = 'readwind_nests'

  do l=1,numbnests
    hflswitch=.false.
    strswitch=.false.
    levdiff2=nlev_ec-nwz+1
    iumax=0
    iwmax=0

    ifile=0
    igrib=0
    iret=0

  !
  ! OPENING OF DATA FILE (GRIB CODE)
  !

5   call grib_open_file(ifile,path(numpath+2*(l-1)+1) &
         (1:length(numpath+2*(l-1)+1))//trim(wfnamen(l,indj)),'r')
  if (iret.ne.GRIB_SUCCESS) then
    goto 888   ! ERROR DETECTED
  endif
  !turn on support for multi fields messages */
  !call grib_multi_support_on

    gotGrid=0
    ifield=0
10   ifield=ifield+1
  !
  ! GET NEXT FIELDS
  !
  call grib_new_from_file(ifile,igrib,iret)
  if (iret.eq.GRIB_END_OF_FILE)  then
    goto 50    ! EOF DETECTED
  elseif (iret.ne.GRIB_SUCCESS) then
    goto 888   ! ERROR DETECTED
  endif

  !first see if we read GRIB1 or GRIB2
  call grib_get_int(igrib,'editionNumber',gribVer,iret)
  call grib_check(iret,gribFunction,gribErrorMsg)

  if (gribVer.eq.1) then ! GRIB Edition 1

  !print*,'GRiB Edition 1'
  !read the grib2 identifiers
  call grib_get_int(igrib,'indicatorOfParameter',isec1(6),iret)
  call grib_check(iret,gribFunction,gribErrorMsg)
  call grib_get_int(igrib,'level',isec1(8),iret)
  call grib_check(iret,gribFunction,gribErrorMsg)

  !change code for etadot to code for omega
  if (isec1(6).eq.77) then
    isec1(6)=135
  endif

  conversion_factor=1.


  else

  !print*,'GRiB Edition 2'
  !read the grib2 identifiers
  call grib_get_int(igrib,'discipline',discipl,iret)
  call grib_check(iret,gribFunction,gribErrorMsg)
  call grib_get_int(igrib,'parameterCategory',parCat,iret)
  call grib_check(iret,gribFunction,gribErrorMsg)
  call grib_get_int(igrib,'parameterNumber',parNum,iret)
  call grib_check(iret,gribFunction,gribErrorMsg)
  call grib_get_int(igrib,'typeOfFirstFixedSurface',typSurf,iret)
  call grib_check(iret,gribFunction,gribErrorMsg)
  call grib_get_int(igrib,'level',valSurf,iret)
  call grib_check(iret,gribFunction,gribErrorMsg)
  call grib_get_int(igrib,'paramId',parId,iret) !added by mc to make it consisitent with new readwind.f90
  call grib_check(iret,gribFunction,gribErrorMsg) !added by mc to make it consisitent with new readwind.f90

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
  elseif ((parCat.eq.2).and.(parNum.eq.32)) then ! W, actually eta dot !
    isec1(6)=135         ! indicatorOfParameter
  elseif ((parCat.eq.128).and.(parNum.eq.77)) then ! W, actually eta dot !added by mc to make it consisitent with new readwind.f90
    isec1(6)=135         ! indicatorOfParameter    !added by mc to make it consisitent with new readwind.f90
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
    conversion_factor=1000. !added by mc to make it consisitent with new readwind.f90
  elseif ((parCat.eq.6).and.(parNum.eq.1) .or. parId .eq. 164) then ! CC !added by mc to make it consisitent with new readwind.f90
    isec1(6)=164         ! indicatorOfParameter
  elseif ((parCat.eq.1).and.(parNum.eq.9) .or. parId .eq. 142) then ! LSP !added by mc to make it consisitent with new readwind.f90
    isec1(6)=142         ! indicatorOfParameter
  elseif ((parCat.eq.1).and.(parNum.eq.10)) then ! CP
    isec1(6)=143         ! indicatorOfParameter
    conversion_factor=1000. !added by mc to make it consisitent with new readwind.f90
  elseif ((parCat.eq.0).and.(parNum.eq.11).and.(typSurf.eq.1)) then ! SHF
    isec1(6)=146         ! indicatorOfParameter
  elseif ((parCat.eq.4).and.(parNum.eq.9).and.(typSurf.eq.1)) then ! SR
    isec1(6)=176         ! indicatorOfParameter
  elseif ((parCat.eq.2).and.(parNum.eq.38) .or. parId .eq. 180) then ! EWSS !added by mc to make it consisitent with new readwind.f90
    isec1(6)=180         ! indicatorOfParameter
  elseif ((parCat.eq.2).and.(parNum.eq.37) .or. parId .eq. 181) then ! NSSS !added by mc to make it consisitent with new readwind.f90
    isec1(6)=181         ! indicatorOfParameter
  elseif ((parCat.eq.3).and.(parNum.eq.4)) then ! ORO
    isec1(6)=129         ! indicatorOfParameter
   elseif ((parCat.eq.3).and.(parNum.eq.7) .or. parId .eq. 160) then ! SDO !added by mc to make it consisitent with new readwind.f90
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
  if(parId .ne. isec1(6) .and. parId .ne. 77) then !added by mc to make it consisitent with new readwind.f90
    write(*,*) 'parId',parId, 'isec1(6)',isec1(6)  !
  !    stop
  endif

  endif

  !HSO  get the size and data of the values array
  if (isec1(6).ne.-1) then
    call grib_get_real4_array(igrib,'values',zsec4,iret)
    call grib_check(iret,gribFunction,gribErrorMsg)
  endif

  !HSO  get the required fields from section 2 in a gribex compatible manner
  if(ifield.eq.1) then
  call grib_get_int(igrib,'numberOfPointsAlongAParallel', &
       isec2(2),iret)
  call grib_check(iret,gribFunction,gribErrorMsg)
  call grib_get_int(igrib,'numberOfPointsAlongAMeridian', &
       isec2(3),iret)
  call grib_check(iret,gribFunction,gribErrorMsg)
  call grib_get_int(igrib,'numberOfVerticalCoordinateValues', &
       isec2(12))
  call grib_check(iret,gribFunction,gribErrorMsg)
  ! CHECK GRID SPECIFICATIONS
  if(isec2(2).ne.nxn(l)) stop &
  'READWIND: NX NOT CONSISTENT FOR A NESTING LEVEL'
  if(isec2(3).ne.nyn(l)) stop &
  'READWIND: NY NOT CONSISTENT FOR A NESTING LEVEL'
  if(isec2(12)/2-1.ne.nlev_ec) stop 'READWIND: VERTICAL DISCRET&
       &IZATION NOT CONSISTENT FOR A NESTING LEVEL'
  endif ! ifield

  !HSO  get the second part of the grid dimensions only from GRiB1 messages
 if (isec1(6) .eq. 167 .and. (gotGrid.eq.0)) then ! !added by mc to make it consisitent with new readwind.f90
    call grib_get_real8(igrib,'longitudeOfFirstGridPointInDegrees', &
         xauxin,iret)
    call grib_check(iret,gribFunction,gribErrorMsg)
    call grib_get_real8(igrib,'latitudeOfLastGridPointInDegrees', &
         yauxin,iret)
    call grib_check(iret,gribFunction,gribErrorMsg)
    if (xauxin.gt.180.) xauxin=xauxin-360.0
    if (xauxin.lt.-180.) xauxin=xauxin+360.0

    xaux=xauxin
    yaux=yauxin
    if (abs(xaux-xlon0n(l)).gt.eps) &
    stop 'READWIND: LOWER LEFT LONGITUDE NOT CONSISTENT FOR A NESTING LEVEL'
    if (abs(yaux-ylat0n(l)).gt.eps) &
    stop 'READWIND: LOWER LEFT LATITUDE NOT CONSISTENT FOR A NESTING LEVEL'
    gotGrid=1
  endif

    do j=0,nyn(l)-1
      do i=0,nxn(l)-1
        k=isec1(8)
        if(isec1(6).eq.130) tthn(i,j,nlev_ec-k+2,n,l)= &!! TEMPERATURE
             zsec4(nxn(l)*(nyn(l)-j-1)+i+1)
        if(isec1(6).eq.131) uuhn(i,j,nlev_ec-k+2,l)= &!! U VELOCITY
             zsec4(nxn(l)*(nyn(l)-j-1)+i+1)
        if(isec1(6).eq.132) vvhn(i,j,nlev_ec-k+2,l)= &!! V VELOCITY
             zsec4(nxn(l)*(nyn(l)-j-1)+i+1)
        if(isec1(6).eq.133) then                         !! SPEC. HUMIDITY
          qvhn(i,j,nlev_ec-k+2,n,l)=zsec4(nxn(l)*(nyn(l)-j-1)+i+1)
          if (qvhn(i,j,nlev_ec-k+2,n,l) .lt. 0.) &
               qvhn(i,j,nlev_ec-k+2,n,l) = 0.
  !          this is necessary because the gridded data may contain
  !          spurious negative values
        endif
        if(isec1(6).eq.134) psn(i,j,1,n,l)= &!! SURF. PRESS.
             zsec4(nxn(l)*(nyn(l)-j-1)+i+1)
        if(isec1(6).eq.135) wwhn(i,j,nlev_ec-k+1,l)= &!! W VELOCITY
             zsec4(nxn(l)*(nyn(l)-j-1)+i+1)
        if(isec1(6).eq.141) sdn(i,j,1,n,l)= &!! SNOW DEPTH
             zsec4(nxn(l)*(nyn(l)-j-1)+i+1)/conversion_factor !added by mc to make it consisitent with new readwind.f90!
        if(isec1(6).eq.151) msln(i,j,1,n,l)= &!! SEA LEVEL PRESS.
             zsec4(nxn(l)*(nyn(l)-j-1)+i+1)
        if(isec1(6).eq.164) tccn(i,j,1,n,l)= &!! CLOUD COVER
             zsec4(nxn(l)*(nyn(l)-j-1)+i+1)
        if(isec1(6).eq.165) u10n(i,j,1,n,l)= &!! 10 M U VELOCITY
             zsec4(nxn(l)*(nyn(l)-j-1)+i+1)
        if(isec1(6).eq.166) v10n(i,j,1,n,l)= &!! 10 M V VELOCITY
             zsec4(nxn(l)*(nyn(l)-j-1)+i+1)
        if(isec1(6).eq.167) tt2n(i,j,1,n,l)= &!! 2 M TEMPERATURE
             zsec4(nxn(l)*(nyn(l)-j-1)+i+1)
        if(isec1(6).eq.168) td2n(i,j,1,n,l)= &!! 2 M DEW POINT
             zsec4(nxn(l)*(nyn(l)-j-1)+i+1)
        if(isec1(6).eq.142) then                         !! LARGE SCALE PREC.
          lsprecn(i,j,1,n,l)=zsec4(nxn(l)*(nyn(l)-j-1)+i+1)
          if (lsprecn(i,j,1,n,l).lt.0.) lsprecn(i,j,1,n,l)=0.
        endif
        if(isec1(6).eq.143) then                         !! CONVECTIVE PREC.
          convprecn(i,j,1,n,l)=zsec4(nxn(l)*(nyn(l)-j-1)+i+1)/conversion_factor !added by mc to make it consisitent with new readwind.f90
          if (convprecn(i,j,1,n,l).lt.0.) convprecn(i,j,1,n,l)=0.
        endif
        if(isec1(6).eq.146) sshfn(i,j,1,n,l)= &!! SENS. HEAT FLUX
             zsec4(nxn(l)*(nyn(l)-j-1)+i+1)
        if((isec1(6).eq.146).and. &
             (zsec4(nxn(l)*(nyn(l)-j-1)+i+1).ne.0.)) hflswitch=.true.    ! Heat flux available
        if(isec1(6).eq.176) then                         !! SOLAR RADIATION
          ssrn(i,j,1,n,l)=zsec4(nxn(l)*(nyn(l)-j-1)+i+1)
          if (ssrn(i,j,1,n,l).lt.0.) ssrn(i,j,1,n,l)=0.
        endif
        if(isec1(6).eq.180) ewss(i,j)= &!! EW SURFACE STRESS
             zsec4(nxn(l)*(nyn(l)-j-1)+i+1)
        if(isec1(6).eq.181) nsss(i,j)= &!! NS SURFACE STRESS
             zsec4(nxn(l)*(nyn(l)-j-1)+i+1)
        if(((isec1(6).eq.180).or.(isec1(6).eq.181)).and. &
             (zsec4(nxn(l)*(nyn(l)-j-1)+i+1).ne.0.)) strswitch=.true.    ! stress available
        if(isec1(6).eq.129) oron(i,j,l)= &!! ECMWF OROGRAPHY
             zsec4(nxn(l)*(nyn(l)-j-1)+i+1)/ga
        if(isec1(6).eq.160) excessoron(i,j,l)= &!! STANDARD DEVIATION OF OROGRAPHY
             zsec4(nxn(l)*(nyn(l)-j-1)+i+1)
        if(isec1(6).eq.172) lsmn(i,j,l)= &!! ECMWF LAND SEA MASK
             zsec4(nxn(l)*(nyn(l)-j-1)+i+1)
        if(isec1(6).eq.131) iumax=max(iumax,nlev_ec-k+1)
        if(isec1(6).eq.135) iwmax=max(iwmax,nlev_ec-k+1)

  ! ESO TODO:
  ! -add check for if one of clwc/ciwc missing (error),
  !    also if all 3 cw fields present, use qc and disregard the others
        if(isec1(6).eq.246) then  !! CLWC  Cloud liquid water content [kg/kg]
          clwchn(i,j,nlev_ec-k+2,n,l)=zsec4(nxn(l)*(nyn(l)-j-1)+i+1)
          readclouds_nest(l)=.true.
          sumclouds_nest(l)=.false.
        endif
        if(isec1(6).eq.247) then  !! CIWC  Cloud ice water content
          ciwchn(i,j,nlev_ec-k+2,n,l)=zsec4(nxn(l)*(nyn(l)-j-1)+i+1)
        endif
  !ZHG end
  !ESO read qc (=clwc+ciwc)
        if(isec1(6).eq.201031) then  !! QC  Cloud liquid water content [kg/kg]
          clwchn(i,j,nlev_ec-k+2,n,l)=zsec4(nxn(l)*(nyn(l)-j-1)+i+1)
          readclouds_nest(l)=.true.
          sumclouds_nest(l)=.true.
        endif


      end do
    end do

  call grib_release(igrib)
  goto 10                      !! READ NEXT LEVEL OR PARAMETER
  !
  ! CLOSING OF INPUT DATA FILE
  !
50   call grib_close_file(ifile)

  !error message if no fields found with correct first longitude in it
  if (gotGrid.eq.0) then
    print*,'***ERROR: input file needs to contain GRiB1 formatted'// &
         'messages'
    stop
  endif

  if(levdiff2.eq.0) then
    iwmax=nlev_ec+1
    do i=0,nxn(l)-1
      do j=0,nyn(l)-1
        wwhn(i,j,nlev_ec+1,l)=0.
      end do
    end do
  endif

  do i=0,nxn(l)-1
    do j=0,nyn(l)-1
      surfstrn(i,j,1,n,l)=sqrt(ewss(i,j)**2+nsss(i,j)**2)
    end do
  end do

  if ((.not.hflswitch).or.(.not.strswitch)) then
    write(*,*) 'WARNING: No flux data contained in GRIB file ', &
         wfnamen(l,indj)

  ! CALCULATE USTAR AND SSHF USING THE PROFILE METHOD
  ! As ECMWF has increased the model resolution, such that now the first model
  ! level is at about 10 m (where 10-m wind is given), use the 2nd ECMWF level
  ! (3rd model level in FLEXPART) for the profile method
  !***************************************************************************

    do i=0,nxn(l)-1
      do j=0,nyn(l)-1
        plev1=akz(3)+bkz(3)*psn(i,j,1,n,l)
        pmean=0.5*(psn(i,j,1,n,l)+plev1)
        tv=tthn(i,j,3,n,l)*(1.+0.61*qvhn(i,j,3,n,l))
        fu=-r_air*tv/ga/pmean
        hlev1=fu*(plev1-psn(i,j,1,n,l))   ! HEIGTH OF FIRST MODEL LAYER
        ff10m= sqrt(u10n(i,j,1,n,l)**2+v10n(i,j,1,n,l)**2)
        fflev1=sqrt(uuhn(i,j,3,l)**2+vvhn(i,j,3,l)**2)
        call pbl_profile(psn(i,j,1,n,l),td2n(i,j,1,n,l),hlev1, &
             tt2n(i,j,1,n,l),tthn(i,j,3,n,l),ff10m,fflev1, &
             surfstrn(i,j,1,n,l),sshfn(i,j,1,n,l))
        if(sshfn(i,j,1,n,l).gt.200.) sshfn(i,j,1,n,l)=200.
        if(sshfn(i,j,1,n,l).lt.-400.) sshfn(i,j,1,n,l)=-400.
      end do
    end do
  endif


  ! Assign 10 m wind to model level at eta=1.0 to have one additional model
  ! level at the ground
  ! Specific humidity is taken the same as at one level above
  ! Temperature is taken as 2 m temperature
  !**************************************************************************

    do i=0,nxn(l)-1
      do j=0,nyn(l)-1
        uuhn(i,j,1,l)=u10n(i,j,1,n,l)
        vvhn(i,j,1,l)=v10n(i,j,1,n,l)
        qvhn(i,j,1,n,l)=qvhn(i,j,2,n,l)
        tthn(i,j,1,n,l)=tt2n(i,j,1,n,l)
      end do
    end do

    if(iumax.ne.nuvz-1) stop &
         'READWIND: NUVZ NOT CONSISTENT FOR A NESTING LEVEL'
    if(iwmax.ne.nwz) stop &
         'READWIND: NWZ NOT CONSISTENT FOR A NESTING LEVEL'

  end do

  return
888   write(*,*) ' #### FLEXPART MODEL ERROR! WINDFIELD         #### '
  write(*,*) ' #### ',wfnamen(l,indj),' FOR NESTING LEVEL  #### '
  write(*,*) ' #### ',l,' IS NOT GRIB FORMAT !!!           #### '
  stop 'Execution terminated'


999   write(*,*) ' #### FLEXPART MODEL ERROR! WINDFIELD         #### '
  write(*,*) ' #### ',wfnamen(l,indj),'                    #### '
  write(*,*) ' #### CANNOT BE OPENED FOR NESTING LEVEL ',l,'####'

end subroutine readwind_nests

subroutine shift_field_0(field,nxf,nyf)
  !                          i/o   i   i
  !*****************************************************************************
  !                                                                            *
  !  This subroutine shifts global fields by nxshift grid cells, in order to   *
  !  facilitate all sorts of nested wind fields, or output grids, which,       *
  !  without shifting, would overlap with the domain "boundary".               *
  !                                                                            *
  !    Author: A. Stohl                                                        *
  !                                                                            *
  !    3 July 2002                                                             *
  !                                                                            *
  !*****************************************************************************
  !                                                                            *
  ! Variables:                                                                 *
  !                                                                            *
  ! Constants:                                                                 *
  !                                                                            *
  !*****************************************************************************

  implicit none

  integer :: nxf,nyf,ix,jy,ixs
  real :: field(0:nxmax-1,0:nymax-1),xshiftaux(0:nxmax-1)

  ! Loop over y and z
  !******************

  do jy=0,nyf-1

  ! Shift the data
  !***************

    if (nxshift.ne.0) then
      do ix=0,nxf-1
        if (ix.ge.nxshift) then
          ixs=ix-nxshift
        else
          ixs=nxf-nxshift+ix
        endif
        xshiftaux(ixs)=field(ix,jy)
      end do
      do ix=0,nxf-1
        field(ix,jy)=xshiftaux(ix)
      end do
    endif

  ! Repeat the westernmost grid cells at the easternmost domain "boundary"
  !***********************************************************************

    field(nxf,jy)=field(0,jy)
  end do

  return
end subroutine shift_field_0

subroutine shift_field(field,nxf,nyf,nzfmax,nzf,nmax,n)
  !                        i/o   i   i    i     i   i   i
  !*****************************************************************************
  !                                                                            *
  !  This subroutine shifts global fields by nxshift grid cells, in order to   *
  !  facilitate all sorts of nested wind fields, or output grids, which,       *
  !  without shifting, would overlap with the domain "boundary".               *
  !                                                                            *
  !    Author: A. Stohl                                                        *
  !                                                                            *
  !    3 July 2002                                                             *
  !                                                                            *
  !*****************************************************************************
  !                                                                            *
  ! Variables:                                                                 *
  !                                                                            *
  ! Constants:                                                                 *
  !                                                                            *
  !*****************************************************************************

  implicit none

  integer :: nxf,nyf,nzf,n,ix,jy,kz,ixs,nzfmax,nmax
  real :: field(0:nxmax-1,0:nymax-1,nzfmax,nmax),xshiftaux(0:nxmax-1)

  ! Loop over y and z
  !******************

  do kz=1,nzf
    do jy=0,nyf-1

  ! Shift the data
  !***************

      if (nxshift.ne.0) then
        do ix=0,nxf-1
          if (ix.ge.nxshift) then
            ixs=ix-nxshift
          else
            ixs=nxf-nxshift+ix
          endif
          xshiftaux(ixs)=field(ix,jy,kz,n)
        end do
        do ix=0,nxf-1
          field(ix,jy,kz,n)=xshiftaux(ix)
        end do
      endif

  ! Repeat the westernmost grid cells at the easternmost domain "boundary"
  !***********************************************************************

      field(nxf,jy,kz,n)=field(0,jy,kz,n)
    end do
  end do
end subroutine shift_field

subroutine verttransform_ecmwf(n,uuh,vvh,wwh,pvh)
  !                              i  i   i   i   i
  !*****************************************************************************
  !                                                                            *
  !     This subroutine transforms temperature, dew point temperature and      *
  !     wind components from eta to meter coordinates.                         *
  !     The vertical wind component is transformed from Pa/s to m/s using      *
  !     the conversion factor pinmconv.                                        *
  !     In addition, this routine calculates vertical density gradients        *
  !     needed for the parameterization of the turbulent velocities.           *
  !                                                                            *
  !     Author: A. Stohl, G. Wotawa                                            *
  !                                                                            *
  !     12 August 1996                                                         *
  !     Update: 16 January 1998                                                *
  !                                                                            *
  !     Major update: 17 February 1999                                         *
  !     by G. Wotawa                                                           *
  !                                                                            *
  !     - Vertical levels for u, v and w are put together                      *
  !     - Slope correction for vertical velocity: Modification of calculation  *
  !       procedure                                                            *
  !                                                                            *
  !*****************************************************************************
  !  Changes, Bernd C. Krueger, Feb. 2001:
  !   Variables tth and qvh (on eta coordinates) from common block
  !
  ! Sabine Eckhardt, March 2007
  ! added the variable cloud for use with scavenging - descr. in com_mod
  !
  ! Unified ECMWF and GFS builds
  ! Marian Harustak, 12.5.2017 
  !     - Renamed from verttransform to verttransform_ecmwf
  !*****************************************************************************
  ! Date: 2017-05-30 modification of a bug in ew. Don Morton (CTBTO project)   *
  !*****************************************************************************
  !                                                                            *
  ! Variables:                                                                 *
  ! nx,ny,nz                        field dimensions in x,y and z direction    *
  ! clouds(0:nxmax,0:nymax,0:nzmax,numwfmem) cloud field for wet deposition    *
  ! uu(0:nxmax,0:nymax,nzmax,numwfmem)     wind components in x-direction [m/s]*
  ! vv(0:nxmax,0:nymax,nzmax,numwfmem)     wind components in y-direction [m/s]*
  ! ww(0:nxmax,0:nymax,nzmax,numwfmem)     wind components in z-direction      *
  !                                          [deltaeta/s]                      *
  ! tt(0:nxmax,0:nymax,nzmax,numwfmem)     temperature [K]                     *
  ! pv(0:nxmax,0:nymax,nzmax,numwfmem)     potential voriticity (pvu)          *
  ! ps(0:nxmax,0:nymax,numwfmem)           surface pressure [Pa]               *
  !                                                                            *
  !*****************************************************************************

  use par_mod
  use com_mod
  use cmapf_mod, only: cc2gll
  use qvsat_mod

  implicit none

  real,intent(in),dimension(0:nxmax-1,0:nymax-1,nuvzmax) :: uuh,vvh,pvh
  real,intent(in),dimension(0:nxmax-1,0:nymax-1,nwzmax) :: wwh

  real,dimension(0:nxmax-1,0:nymax-1,nuvzmax) :: rhoh
  real,dimension(0:nxmax-1,0:nymax-1,nzmax) :: pinmconv
  ! RLT added pressure
  real,dimension(0:nxmax-1,0:nymax-1,nuvzmax) :: prsh
  real,dimension(0:nxmax-1,0:nymax-1) ::  tvold,pold,pint,tv,dpdeta
  real,dimension(0:nymax-1) :: cosf

  integer,dimension(0:nxmax-1,0:nymax-1) :: rain_cloud_above,idx

  integer :: ix,jy,kz,iz,n,kmin,ix1,jy1,ixp,jyp,ixm,jym,kz_inv
  real :: pressure,rh,lsp,convp,cloudh_min,prec
  real :: dz1,dz2,dz
  real :: xlon,ylat,xlonr,dzdx,dzdy
  real :: dzdx1,dzdx2,dzdy1,dzdy2
  real :: uuaux,vvaux,uupolaux,vvpolaux,ddpol,ffpol,wdummy
  real,parameter :: const=r_air/ga
  real,parameter :: precmin = 0.002 ! minimum prec in mm/h for cloud diagnostics

  logical :: init = .true.

  !*************************************************************************
  ! If verttransform is called the first time, initialize heights of the   *
  ! z levels in meter. The heights are the heights of model levels, where  *
  ! u,v,T and qv are given, and of the interfaces, where w is given. So,   *
  ! the vertical resolution in the z system is doubled. As reference point,*
  ! the lower left corner of the grid is used.                             *
  ! Unlike in the eta system, no difference between heights for u,v and    *
  ! heights for w exists.                                                  *
  !*************************************************************************


  !eso measure CPU time
  !  call mpif_mtime('verttransform',0)

  if (init) then

  ! Search for a point with high surface pressure (i.e. not above significant topography)
  ! Then, use this point to construct a reference z profile, to be used at all times
  !*****************************************************************************

    loop1: do jy=0,nymin1
      do ix=0,nxmin1
        if (ps(ix,jy,1,n).gt.100000.) then
          ixm=ix
          jym=jy
          exit loop1
        endif
      end do
    end do loop1

    tvold(ixm,jym)=tt2(ixm,jym,1,n)*(1.+0.378*ew(td2(ixm,jym,1,n),ps(ixm,jym,1,n))/ &
         ps(ixm,jym,1,n))
    pold(ixm,jym)=ps(ixm,jym,1,n)
    height(1)=0.

    do kz=2,nuvz
      pint=akz(kz)+bkz(kz)*ps(ixm,jym,1,n)
      tv=tth(ixm,jym,kz,n)*(1.+0.608*qvh(ixm,jym,kz,n))

      if (abs(tv(ixm,jym)-tvold(ixm,jym)).gt.0.2) then
        height(kz)= &
             height(kz-1)+const*log(pold(ixm,jym)/pint(ixm,jym))* &
             (tv(ixm,jym)-tvold(ixm,jym))/log(tv(ixm,jym)/tvold(ixm,jym))
      else
        height(kz)=height(kz-1)+ &
             const*log(pold(ixm,jym)/pint(ixm,jym))*tv(ixm,jym)
      endif

      tvold(ixm,jym)=tv(ixm,jym)
      pold(ixm,jym)=pint(ixm,jym)
    end do

  ! Determine highest levels that can be within PBL
  !************************************************

    do kz=1,nz
      if (height(kz).gt.hmixmax) then
        nmixz=kz
        exit
      endif
    end do

  ! Do not repeat initialization of the Cartesian z grid
  !*****************************************************

    init=.false.
  endif


  ! Loop over the whole grid
  !*************************

  do jy=0,nymin1
    do ix=0,nxmin1
      tvold(ix,jy)=tt2(ix,jy,1,n)*(1.+0.378*ew(td2(ix,jy,1,n),ps(ix,jy,1,n))/ &
           ps(ix,jy,1,n))
    enddo
  enddo
  pold=ps(:,:,1,n)
  etauvheight(:,:,1,n)=0.
  etawheight(:,:,1,n)=0.
  rhoh(:,:,1)=pold/(r_air*tvold)
  ! RLT add pressure
  prsh(:,:,1)=ps(:,:,1,n)


  ! Compute heights of eta levels
  !******************************

  do kz=2,nuvz
    pint=akz(kz)+bkz(kz)*ps(:,:,1,n)
    ! RLT add pressure
    prsh(:,:,kz)=pint
    tv=tth(:,:,kz,n)*(1.+0.608*qvh(:,:,kz,n))
    rhoh(:,:,kz)=pint(:,:)/(r_air*tv)

    where (abs(tv-tvold).gt.0.2) 
      etauvheight(:,:,kz,n)=etauvheight(:,:,kz-1,n)+const*log(pold/pint)* &
           (tv-tvold)/log(tv/tvold)
    elsewhere
      etauvheight(:,:,kz,n)=etauvheight(:,:,kz-1,n)+const*log(pold/pint)*tv
    endwhere

    tvold=tv
    pold=pint
  end do


  do kz=2,nwz-1
    etawheight(:,:,kz,n)=(etauvheight(:,:,kz+1,n)+etauvheight(:,:,kz,n))/2.
  end do
  etawheight(:,:,nwz,n)=etawheight(:,:,nwz-1,n)+ &
       etauvheight(:,:,nuvz,n)-etauvheight(:,:,nuvz-1,n)

  ! pinmconv=(h2-h1)/(p2-p1)

  pinmconv(:,:,1)=(etauvheight(:,:,2,n))/ &
       ((aknew(2)+bknew(2)*ps(:,:,1,n))- &
       (aknew(1)+bknew(1)*ps(:,:,1,n)))
  do kz=2,nz-1
    pinmconv(:,:,kz)=(etauvheight(:,:,kz+1,n)-etauvheight(:,:,kz-1,n))/ &
         ((aknew(kz+1)+bknew(kz+1)*ps(:,:,1,n))- &
         (aknew(kz-1)+bknew(kz-1)*ps(:,:,1,n)))
  end do
  pinmconv(:,:,nz)=(etauvheight(:,:,nz,n)-etauvheight(:,:,nz-1,n))/ &
       ((aknew(nz)+bknew(nz)*ps(:,:,1,n))- &
       (aknew(nz-1)+bknew(nz-1)*ps(:,:,1,n)))

  ! Levels, where u,v,t and q are given
  !************************************


  uu(:,:,1,n)=uuh(:,:,1)
  vv(:,:,1,n)=vvh(:,:,1)
  tt(:,:,1,n)=tth(:,:,1,n)
  qv(:,:,1,n)=qvh(:,:,1,n)
  !hg adding the cloud water 
  if (readclouds) then
    clwc(:,:,1,n)=clwch(:,:,1,n)
    if (.not.sumclouds) ciwc(:,:,1,n)=ciwch(:,:,1,n)
  end if
  !hg 
  pv(:,:,1,n)=pvh(:,:,1)
  rho(:,:,1,n)=rhoh(:,:,1)
  ! RLT add pressure
  prs(:,:,1,n)=prsh(:,:,1)

  uu(:,:,nz,n)=uuh(:,:,nuvz)
  vv(:,:,nz,n)=vvh(:,:,nuvz)
  tt(:,:,nz,n)=tth(:,:,nuvz,n)
  qv(:,:,nz,n)=qvh(:,:,nuvz,n)
  !hg adding the cloud water
  if (readclouds) then
    clwc(:,:,nz,n)=clwch(:,:,nuvz,n)
    if (.not.sumclouds) ciwc(:,:,nz,n)=ciwch(:,:,nuvz,n)
  end if
  !hg
  pv(:,:,nz,n)=pvh(:,:,nuvz)
  rho(:,:,nz,n)=rhoh(:,:,nuvz)
  ! RLT
  prs(:,:,nz,n)=prsh(:,:,nuvz)

  kmin=2
  idx=kmin
  do iz=2,nz-1
    do jy=0,nymin1
      do ix=0,nxmin1
        if(height(iz).gt.etauvheight(ix,jy,nuvz,n)) then
          uu(ix,jy,iz,n)=uu(ix,jy,nz,n)
          vv(ix,jy,iz,n)=vv(ix,jy,nz,n)
          tt(ix,jy,iz,n)=tt(ix,jy,nz,n)
          qv(ix,jy,iz,n)=qv(ix,jy,nz,n)
  !hg adding the cloud water
          if (readclouds) then
            clwc(ix,jy,iz,n)=clwc(ix,jy,nz,n)
            if (.not.sumclouds) ciwc(ix,jy,iz,n)=ciwc(ix,jy,nz,n)
          end if
  !hg
          pv(ix,jy,iz,n)=pv(ix,jy,nz,n)
          rho(ix,jy,iz,n)=rho(ix,jy,nz,n)
  ! RLT
          prs(ix,jy,iz,n)=prs(ix,jy,nz,n)
        else
          innuvz: do kz=idx(ix,jy),nuvz
            if (idx(ix,jy) .le. kz .and. (height(iz).gt.etauvheight(ix,jy,kz-1,n)).and. &
                 (height(iz).le.etauvheight(ix,jy,kz,n))) then
              idx(ix,jy)=kz
              exit innuvz
            endif
          enddo innuvz
        endif
      enddo
    enddo
    do jy=0,nymin1
      do ix=0,nxmin1
        if(height(iz).le.etauvheight(ix,jy,nuvz,n)) then
          kz=idx(ix,jy)
          dz1=height(iz)-etauvheight(ix,jy,kz-1,n)
          dz2=etauvheight(ix,jy,kz,n)-height(iz)
          dz=dz1+dz2
          uu(ix,jy,iz,n)=(uuh(ix,jy,kz-1)*dz2+uuh(ix,jy,kz)*dz1)/dz
          vv(ix,jy,iz,n)=(vvh(ix,jy,kz-1)*dz2+vvh(ix,jy,kz)*dz1)/dz
          tt(ix,jy,iz,n)=(tth(ix,jy,kz-1,n)*dz2 &
               +tth(ix,jy,kz,n)*dz1)/dz
          qv(ix,jy,iz,n)=(qvh(ix,jy,kz-1,n)*dz2 &
               +qvh(ix,jy,kz,n)*dz1)/dz
  !hg adding the cloud water
          if (readclouds) then
            clwc(ix,jy,iz,n)=(clwch(ix,jy,kz-1,n)*dz2+clwch(ix,jy,kz,n)*dz1)/dz
            if (.not.sumclouds) &
                 &ciwc(ix,jy,iz,n)=(ciwch(ix,jy,kz-1,n)*dz2+ciwch(ix,jy,kz,n)*dz1)/dz
          end if
  !hg
          pv(ix,jy,iz,n)=(pvh(ix,jy,kz-1)*dz2+pvh(ix,jy,kz)*dz1)/dz
          rho(ix,jy,iz,n)=(rhoh(ix,jy,kz-1)*dz2+rhoh(ix,jy,kz)*dz1)/dz
  ! RLT add pressure
          prs(ix,jy,iz,n)=(prsh(ix,jy,kz-1)*dz2+prsh(ix,jy,kz)*dz1)/dz
        endif
      enddo
    enddo
  enddo


  ! Levels, where w is given
  !*************************

  ww(:,:,1,n)=wwh(:,:,1)*pinmconv(:,:,1)
  ww(:,:,nz,n)=wwh(:,:,nwz)*pinmconv(:,:,nz)
  kmin=2
  idx=kmin
  do iz=2,nz
    do jy=0,nymin1
      do ix=0,nxmin1
        inn:         do kz=idx(ix,jy),nwz
          if(idx(ix,jy) .le. kz .and. height(iz).gt.etawheight(ix,jy,kz-1,n).and. &
               height(iz).le.etawheight(ix,jy,kz,n)) then
            idx(ix,jy)=kz
            exit inn
          endif
        enddo inn
      enddo
    enddo
    do jy=0,nymin1
      do ix=0,nxmin1
        kz=idx(ix,jy)
        dz1=height(iz)-etawheight(ix,jy,kz-1,n)
        dz2=etawheight(ix,jy,kz,n)-height(iz)
        dz=dz1+dz2
        ww(ix,jy,iz,n)=(wwh(ix,jy,kz-1)*pinmconv(ix,jy,kz-1)*dz2 &
             +wwh(ix,jy,kz)*pinmconv(ix,jy,kz)*dz1)/dz
      enddo
    enddo
  enddo

  ! Compute density gradients at intermediate levels
  !*************************************************

  drhodz(:,:,1,n)=(rho(:,:,2,n)-rho(:,:,1,n))/ &
       (height(2)-height(1))
  do kz=2,nz-1
    drhodz(:,:,kz,n)=(rho(:,:,kz+1,n)-rho(:,:,kz-1,n))/ &
         (height(kz+1)-height(kz-1))
  end do
  drhodz(:,:,nz,n)=drhodz(:,:,nz-1,n)

  !    end do
  !  end do


  !****************************************************************
  ! Compute slope of eta levels in windward direction and resulting
  ! vertical wind correction
  !****************************************************************

  do jy=1,ny-2
    cosf(jy)=1./cos((real(jy)*dy+ylat0)*pi180)
  enddo

  kmin=2
  idx=kmin
  do iz=2,nz-1
    do jy=1,ny-2
      do ix=1,nx-2

        inneta: do kz=idx(ix,jy),nz
          if (idx(ix,jy) .le. kz .and. (height(iz).gt.etauvheight(ix,jy,kz-1,n)).and. &
               (height(iz).le.etauvheight(ix,jy,kz,n))) then
            idx(ix,jy)=kz
            exit inneta
          endif
        enddo inneta
      enddo
    enddo

    do jy=1,ny-2
      do ix=1,nx-2
        kz=idx(ix,jy)
        dz1=height(iz)-etauvheight(ix,jy,kz-1,n)
        dz2=etauvheight(ix,jy,kz,n)-height(iz)
        dz=dz1+dz2
        ix1=ix-1
        jy1=jy-1
        ixp=ix+1
        jyp=jy+1

        dzdx1=(etauvheight(ixp,jy,kz-1,n)-etauvheight(ix1,jy,kz-1,n))/2.
        dzdx2=(etauvheight(ixp,jy,kz,n)-etauvheight(ix1,jy,kz,n))/2.
        dzdx=(dzdx1*dz2+dzdx2*dz1)/dz

        dzdy1=(etauvheight(ix,jyp,kz-1,n)-etauvheight(ix,jy1,kz-1,n))/2.
        dzdy2=(etauvheight(ix,jyp,kz,n)-etauvheight(ix,jy1,kz,n))/2.
        dzdy=(dzdy1*dz2+dzdy2*dz1)/dz

        ww(ix,jy,iz,n)=ww(ix,jy,iz,n)+(dzdx*uu(ix,jy,iz,n)*dxconst*cosf(jy)+dzdy*vv(ix,jy,iz,n)*dyconst)

      end do

    end do
  end do

  ! Keep original fields if wind_coord_type==ETA
  if (wind_coord_type.eq.'ETA') then
    uueta(:,:,:,n) = uuh(:,:,:)
    vveta(:,:,:,n) = vvh(:,:,:)
    tteta(:,:,:,n) = tth(:,:,:,n)
    qveta(:,:,:,n) = qvh(:,:,:,n)
    pveta(:,:,:,n) = pvh(:,:,:)
    rhoeta(:,:,:,n) = rhoh(:,:,:)
    drhodzeta(:,:,1,n)=(rhoeta(:,:,2,n)-rhoeta(:,:,1,n))/ &
         (height(2)-height(1))
    do kz=2,nz-1
      drhodzeta(:,:,kz,n)=(rhoeta(:,:,kz+1,n)-rhoeta(:,:,kz-1,n))/ &
           (height(kz+1)-height(kz-1)) ! Note that this is still in SI units and not in eta
    end do
    drhodzeta(:,:,nz,n)=drhodzeta(:,:,nz-1,n)
    tvirtual(:,:,:,n)=tteta(:,:,:,n)* &
      ((qveta(:,:,:,n)+0.622)/(0.622*qveta(:,:,:,n)+0.622)) ! eq A11 from Mid-latitude atmospheric dynamics by Jonathan E. Martin
    !tvirtual(:,:,:,n)=tteta(:,:,:,n)*(1.+0.608*qveta(:,:,:,n))
    do jy=0,ny-1
      do ix=0,nx-1
        tvirtual(ix,jy,1,n)=tt2(ix,jy,1,n)*(1.+0.378*ew(td2(ix,jy,1,n),ps(ix,jy,1,n))/ps(ix,jy,1,n))
      end do 
    end do


    ! Convert w from Pa/s to eta/s, following FLEXTRA
    !************************************************
    do kz=1,nuvz-1
      if (kz.eq.1) then
        dpdeta=(akm(kz+1)-akm(kz)+(bkm(kz+1)-bkm(kz))*ps(:,:,1,n))/ &
          (wheight(kz+1)-wheight(kz))
      else if (kz.eq.nuvz-1) then
        dpdeta=(akm(kz)-akm(kz-1)+(bkm(kz)-bkm(kz-1))*ps(:,:,1,n))/ &
          (wheight(kz)-wheight(kz-1))
      else
        dpdeta=(akm(kz+1)-akm(kz-1)+(bkm(kz+1)-bkm(kz-1))*ps(:,:,1,n))/ &
          (wheight(kz+1)-wheight(kz-1))
      endif
      wweta(:,:,kz,n)=wwh(:,:,kz)/dpdeta
    end do
  endif

  ! If north pole is in the domain, calculate wind velocities in polar
  ! stereographic coordinates
  !*******************************************************************

  if (nglobal) then
    do iz=1,nz
      do jy=int(switchnorthg)-2,nymin1
        ylat=ylat0+real(jy)*dy
        do ix=0,nxmin1
          xlon=xlon0+real(ix)*dx
          call cc2gll(northpolemap,ylat,xlon,uu(ix,jy,iz,n), &
               vv(ix,jy,iz,n),uupol(ix,jy,iz,n), &
               vvpol(ix,jy,iz,n))
          if (wind_coord_type.eq.'ETA') then
            call cc2gll(northpolemap,ylat,xlon,uueta(ix,jy,iz,n), &
                 vveta(ix,jy,iz,n),uupoleta(ix,jy,iz,n), &
                 vvpoleta(ix,jy,iz,n))
          endif
        end do
      end do
    end do


    do iz=1,nz

  ! CALCULATE FFPOL, DDPOL FOR CENTRAL GRID POINT
  !
  !   AMSnauffer Nov 18 2004 Added check for case vv=0
  !
      xlon=xlon0+real(nx/2-1)*dx
      xlonr=xlon*pi/180.
      ffpol=sqrt(uu(nx/2-1,nymin1,iz,n)**2+ &
           vv(nx/2-1,nymin1,iz,n)**2)
      if (vv(nx/2-1,nymin1,iz,n).lt.0.) then
        ddpol=atan(uu(nx/2-1,nymin1,iz,n)/ &
             vv(nx/2-1,nymin1,iz,n))-xlonr
      else if (vv(nx/2-1,nymin1,iz,n).gt.0.) then
        ddpol=pi+atan(uu(nx/2-1,nymin1,iz,n)/ &
             vv(nx/2-1,nymin1,iz,n))-xlonr
      else
        ddpol=pi/2-xlonr
      endif
      if(ddpol.lt.0.) ddpol=2.0*pi+ddpol
      if(ddpol.gt.2.0*pi) ddpol=ddpol-2.0*pi

  ! CALCULATE U,V FOR 180 DEG, TRANSFORM TO POLAR STEREOGRAPHIC GRID
      xlon=180.0
      xlonr=xlon*pi/180.
      ylat=90.0
      uuaux=-ffpol*sin(xlonr+ddpol)
      vvaux=-ffpol*cos(xlonr+ddpol)
      call cc2gll(northpolemap,ylat,xlon,uuaux,vvaux,uupolaux, &
           vvpolaux)

      jy=nymin1
      do ix=0,nxmin1
        uupol(ix,jy,iz,n)=uupolaux
        vvpol(ix,jy,iz,n)=vvpolaux
      end do
    end do

    if (wind_coord_type.eq.'ETA') then    
      do iz=1,nz

        xlon=xlon0+real(nx/2-1)*dx
        xlonr=xlon*pi/180.
        ffpol=sqrt(uueta(nx/2-1,nymin1,iz,n)**2+ &
             vveta(nx/2-1,nymin1,iz,n)**2)
        if (vveta(nx/2-1,nymin1,iz,n).lt.0.) then
          ddpol=atan(uueta(nx/2-1,nymin1,iz,n)/ &
               vveta(nx/2-1,nymin1,iz,n))-xlonr
        else if (vveta(nx/2-1,nymin1,iz,n).gt.0.) then
          ddpol=pi+atan(uueta(nx/2-1,nymin1,iz,n)/ &
               vveta(nx/2-1,nymin1,iz,n))-xlonr
        else
          ddpol=pi/2-xlonr
        endif
        if(ddpol.lt.0.) ddpol=2.0*pi+ddpol
        if(ddpol.gt.2.0*pi) ddpol=ddpol-2.0*pi

  ! CALCULATE U,V FOR 180 DEG, TRANSFORM TO POLAR STEREOGRAPHIC GRID
        xlon=180.0
        xlonr=xlon*pi/180.
        ylat=90.0
        uuaux=-ffpol*sin(xlonr+ddpol)
        vvaux=-ffpol*cos(xlonr+ddpol)
        call cc2gll(northpolemap,ylat,xlon,uuaux,vvaux,uupolaux, &
             vvpolaux)

        jy=nymin1
        do ix=0,nxmin1
          uupoleta(ix,jy,iz,n)=uupolaux
          vvpoleta(ix,jy,iz,n)=vvpolaux
        end do
      end do
    endif


  ! Fix: Set W at pole to the zonally averaged W of the next equator-
  ! ward parallel of latitude

    do iz=1,nz
      wdummy=0.
      jy=ny-2
      do ix=0,nxmin1
        wdummy=wdummy+ww(ix,jy,iz,n)
      end do
      wdummy=wdummy/real(nx)
      jy=nymin1
      do ix=0,nxmin1
        ww(ix,jy,iz,n)=wdummy
      end do
    end do

    if (wind_coord_type.eq.'ETA') then
      do iz=1,nz
        wdummy=0.
        jy=ny-2
        do ix=0,nxmin1
          wdummy=wdummy+wweta(ix,jy,iz,n)
        end do
        wdummy=wdummy/real(nx)
        jy=nymin1
        do ix=0,nxmin1
          wweta(ix,jy,iz,n)=wdummy
        end do
      end do
    endif

  endif


  ! If south pole is in the domain, calculate wind velocities in polar
  ! stereographic coordinates
  !*******************************************************************

  if (sglobal) then
    do iz=1,nz
      do jy=0,int(switchsouthg)+3
        ylat=ylat0+real(jy)*dy
        do ix=0,nxmin1
          xlon=xlon0+real(ix)*dx
          call cc2gll(southpolemap,ylat,xlon,uu(ix,jy,iz,n), &
               vv(ix,jy,iz,n),uupol(ix,jy,iz,n), &
               vvpol(ix,jy,iz,n))
          if (wind_coord_type.eq.'ETA') then
            call cc2gll(southpolemap,ylat,xlon,uueta(ix,jy,iz,n), &
                 vveta(ix,jy,iz,n),uupoleta(ix,jy,iz,n), &
                 vvpoleta(ix,jy,iz,n))
          endif
        end do
      end do
    end do

    do iz=1,nz

  ! CALCULATE FFPOL, DDPOL FOR CENTRAL GRID POINT
  !
  !   AMSnauffer Nov 18 2004 Added check for case vv=0
  !
      xlon=xlon0+real(nx/2-1)*dx
      xlonr=xlon*pi/180.
      ffpol=sqrt(uu(nx/2-1,0,iz,n)**2+ &
           vv(nx/2-1,0,iz,n)**2)
      if (vv(nx/2-1,0,iz,n).lt.0.) then
        ddpol=atan(uu(nx/2-1,0,iz,n)/ &
             vv(nx/2-1,0,iz,n))+xlonr
      else if (vv(nx/2-1,0,iz,n).gt.0.) then
        ddpol=pi+atan(uu(nx/2-1,0,iz,n)/ &
             vv(nx/2-1,0,iz,n))+xlonr
      else
        ddpol=pi/2-xlonr
      endif
      if(ddpol.lt.0.) ddpol=2.0*pi+ddpol
      if(ddpol.gt.2.0*pi) ddpol=ddpol-2.0*pi

  ! CALCULATE U,V FOR 180 DEG, TRANSFORM TO POLAR STEREOGRAPHIC GRID
      xlon=180.0
      xlonr=xlon*pi/180.
      ylat=-90.0
      uuaux=+ffpol*sin(xlonr-ddpol)
      vvaux=-ffpol*cos(xlonr-ddpol)
      call cc2gll(northpolemap,ylat,xlon,uuaux,vvaux,uupolaux, &
           vvpolaux)

      jy=0
      do ix=0,nxmin1
        uupol(ix,jy,iz,n)=uupolaux
        vvpol(ix,jy,iz,n)=vvpolaux
      end do
    end do

    if (wind_coord_type.eq.'ETA') then
      do iz=1,nz
  ! CALCULATE FFPOL, DDPOL FOR CENTRAL GRID POINT
  !
  !   AMSnauffer Nov 18 2004 Added check for case vv=0
  !
        xlon=xlon0+real(nx/2-1)*dx
        xlonr=xlon*pi/180.
        ffpol=sqrt(uueta(nx/2-1,0,iz,n)**2+ &
             vveta(nx/2-1,0,iz,n)**2)
        if (vveta(nx/2-1,0,iz,n).lt.0.) then
          ddpol=atan(uueta(nx/2-1,0,iz,n)/ &
               vveta(nx/2-1,0,iz,n))+xlonr
        else if (vveta(nx/2-1,0,iz,n).gt.0.) then
          ddpol=pi+atan(uueta(nx/2-1,0,iz,n)/ &
               vveta(nx/2-1,0,iz,n))+xlonr
        else
          ddpol=pi/2-xlonr
        endif
        if(ddpol.lt.0.) ddpol=2.0*pi+ddpol
        if(ddpol.gt.2.0*pi) ddpol=ddpol-2.0*pi

  ! CALCULATE U,V FOR 180 DEG, TRANSFORM TO POLAR STEREOGRAPHIC GRID
        xlon=180.0
        xlonr=xlon*pi/180.
        ylat=-90.0
        uuaux=+ffpol*sin(xlonr-ddpol)
        vvaux=-ffpol*cos(xlonr-ddpol)
        call cc2gll(northpolemap,ylat,xlon,uuaux,vvaux,uupolaux, &
             vvpolaux)

        jy=0
        do ix=0,nxmin1
          uupoleta(ix,jy,iz,n)=uupolaux
          vvpoleta(ix,jy,iz,n)=vvpolaux
        end do
      end do
    endif

  ! Fix: Set W at pole to the zonally averaged W of the next equator-
  ! ward parallel of latitude

    do iz=1,nz
      wdummy=0.
      jy=1
      do ix=0,nxmin1
        wdummy=wdummy+ww(ix,jy,iz,n)
      end do
      wdummy=wdummy/real(nx)
      jy=0
      do ix=0,nxmin1
        ww(ix,jy,iz,n)=wdummy
      end do
    end do

    if (wind_coord_type.eq.'ETA') then
      do iz=1,nz
        wdummy=0.
        jy=1
        do ix=0,nxmin1
          wdummy=wdummy+wweta(ix,jy,iz,n)
        end do
        wdummy=wdummy/real(nx)
        jy=0
        do ix=0,nxmin1
          wweta(ix,jy,iz,n)=wdummy
        end do
      end do
    endif
  endif


  !***********************************************************************************  
  if (readclouds) then !HG METHOD
  ! The method is loops all grids vertically and constructs the 3D matrix for clouds
  ! Cloud top and cloud bottom gid cells are assigned as well as the total column 
  ! cloud water. For precipitating grids, the type and whether it is in or below 
  ! cloud scavenging are assigned with numbers 2-5 (following the old metod).
  ! Distinction is done for lsp and convp though they are treated the same in regards
  ! to scavenging. Also clouds that are not precipitating are defined which may be 
  ! to include future cloud processing by non-precipitating-clouds. 
  !***********************************************************************************
    write(*,*) 'Global ECMWF fields: using cloud water'
    clw(:,:,:,n)=0.0
  !    icloud_stats(:,:,:,n)=0.0
    ctwc(:,:,n)=0.0
    clouds(:,:,:,n)=0
  ! If water/ice are read separately into clwc and ciwc, store sum in clwc
    if (.not.sumclouds) then 
      clwc(:,:,:,n) = clwc(:,:,:,n) + ciwc(:,:,:,n)
    end if
    do jy=0,nymin1
      do ix=0,nxmin1
        lsp=lsprec(ix,jy,1,n)
        convp=convprec(ix,jy,1,n)
        prec=lsp+convp
  !        tot_cloud_h=0
  ! Find clouds in the vertical
        do kz=1, nz-1 !go from top to bottom
          if (clwc(ix,jy,kz,n).gt.0) then      
  ! assuming rho is in kg/m3 and hz in m gives: kg/kg * kg/m3 *m3/kg /m = m2/m3 
            clw(ix,jy,kz,n)=(clwc(ix,jy,kz,n)*rho(ix,jy,kz,n))*(height(kz+1)-height(kz))
  !            tot_cloud_h=tot_cloud_h+(height(kz+1)-height(kz)) 
            
  !            icloud_stats(ix,jy,4,n)= icloud_stats(ix,jy,4,n)+clw(ix,jy,kz,n)          ! Column cloud water [m3/m3]
            ctwc(ix,jy,n) = ctwc(ix,jy,n)+clw(ix,jy,kz,n)
  !            icloud_stats(ix,jy,3,n)= min(height(kz+1),height(kz))                     ! Cloud BOT height stats      [m]
            cloudh_min=min(height(kz+1),height(kz))
          endif
        end do

  ! If Precipitation. Define removal type in the vertical
        if ((lsp.gt.0.01).or.(convp.gt.0.01)) then ! cloud and precipitation

          do kz=nz,2,-1 !go Bottom up!
            if (clw(ix,jy,kz,n).gt. 0) then ! is in cloud
              cloudsh(ix,jy,n)=cloudsh(ix,jy,n)+height(kz)-height(kz-1) 
              clouds(ix,jy,kz,n)=1                               ! is a cloud
              if (lsp.ge.convp) then
                clouds(ix,jy,kz,n)=3                            ! lsp in-cloud
              else
                clouds(ix,jy,kz,n)=2                             ! convp in-cloud
              endif                                              ! convective or large scale
            elseif((clw(ix,jy,kz,n).le.0) .and. (cloudh_min.ge.height(kz))) then ! is below cloud
              if (lsp.ge.convp) then
                clouds(ix,jy,kz,n)=5                             ! lsp dominated washout
              else
                clouds(ix,jy,kz,n)=4                             ! convp dominated washout
              endif                                              ! convective or large scale 
            endif

            if (height(kz).ge. 19000) then                        ! set a max height for removal
              clouds(ix,jy,kz,n)=0
            endif !clw>0
          end do !nz
        endif ! precipitation
      end do
    end do

  ! eso: copy the relevant data to clw4 to reduce amount of communicated data for MPI
  !    ctwc(:,:,n) = icloud_stats(:,:,4,n)

  !**************************************************************************
  else       ! use old definitions
  !**************************************************************************
  !   create a cloud and rainout/washout field, clouds occur where rh>80%
  !   total cloudheight is stored at level 0
    write(*,*) 'Global fields: using cloud water from Parameterization'
    do jy=0,nymin1
      do ix=0,nxmin1
  ! OLD METHOD
        rain_cloud_above(ix,jy)=0
        lsp=lsprec(ix,jy,1,n)
        convp=convprec(ix,jy,1,n)
        cloudsh(ix,jy,n)=0
        do kz_inv=1,nz-1
          kz=nz-kz_inv+1
          pressure=rho(ix,jy,kz,n)*r_air*tt(ix,jy,kz,n)
          rh=qv(ix,jy,kz,n)/f_qvsat(pressure,tt(ix,jy,kz,n))
          clouds(ix,jy,kz,n)=0
          if (rh.gt.0.8) then ! in cloud
            if ((lsp.gt.0.01).or.(convp.gt.0.01)) then ! cloud and precipitation
              rain_cloud_above(ix,jy)=1
              cloudsh(ix,jy,n)=cloudsh(ix,jy,n)+ &
                   height(kz)-height(kz-1)
              if (lsp.ge.convp) then
                clouds(ix,jy,kz,n)=3 ! lsp dominated rainout
              else
                clouds(ix,jy,kz,n)=2 ! convp dominated rainout
              endif
            else ! no precipitation
              clouds(ix,jy,kz,n)=1 ! cloud
            endif
          else ! no cloud
            if (rain_cloud_above(ix,jy).eq.1) then ! scavenging
              if (lsp.ge.convp) then
                clouds(ix,jy,kz,n)=5 ! lsp dominated washout
              else
                clouds(ix,jy,kz,n)=4 ! convp dominated washout
              endif
            endif
          endif
        end do
  !END OLD METHOD
      end do
    end do
  endif !readclouds
end subroutine verttransform_ecmwf

subroutine verttransform_gfs(n,uuh,vvh,wwh,pvh)
  !                      i  i   i   i   i
  !*****************************************************************************
  !                                                                            *
  !     This subroutine transforms temperature, dew point temperature and      *
  !     wind components from eta to meter coordinates.                         *
  !     The vertical wind component is transformed from Pa/s to m/s using      *
  !     the conversion factor pinmconv.                                        *
  !     In addition, this routine calculates vertical density gradients        *
  !     needed for the parameterization of the turbulent velocities.           *
  !                                                                            *
  !     Author: A. Stohl, G. Wotawa                                            *
  !                                                                            *
  !     12 August 1996                                                         *
  !     Update: 16 January 1998                                                *
  !                                                                            *
  !     Major update: 17 February 1999                                         *
  !     by G. Wotawa                                                           *
  !     CHANGE 17/11/2005 Caroline Forster, NCEP GFS version                   *
  !                                                                            *
  !   - Vertical levels for u, v and w are put together                        *
  !   - Slope correction for vertical velocity: Modification of calculation    *
  !     procedure                                                              *
  !                                                                            *
  !*****************************************************************************
  !  Changes, Bernd C. Krueger, Feb. 2001:
  !   Variables tth and qvh (on eta coordinates) from common block
  !
  !   Unified ECMWF and GFS builds                                      
  !   Marian Harustak, 12.5.2017                                        
  !     - Renamed routine from verttransform to verttransform_gfs
  !
  !*****************************************************************************
  !                                                                            *
  ! Variables:                                                                 *
  ! nx,ny,nz                        field dimensions in x,y and z direction    *
  ! uu(0:nxmax,0:nymax,nzmax,2)     wind components in x-direction [m/s]       *
  ! vv(0:nxmax,0:nymax,nzmax,2)     wind components in y-direction [m/s]       *
  ! ww(0:nxmax,0:nymax,nzmax,2)     wind components in z-direction [deltaeta/s]*
  ! tt(0:nxmax,0:nymax,nzmax,2)     temperature [K]                            *
  ! pv(0:nxmax,0:nymax,nzmax,2)     potential voriticity (pvu)                 *
  ! ps(0:nxmax,0:nymax,2)           surface pressure [Pa]                      *
  ! clouds(0:nxmax,0:nymax,0:nzmax,2) cloud field for wet deposition           *
  !                                                                            *
  !*****************************************************************************

  use par_mod
  use com_mod
  use cmapf_mod
  use qvsat_mod

  implicit none

  integer :: ix,jy,kz,iz,n,kmin,kl,klp,ix1,jy1,ixp,jyp,ixm,jym
  integer :: rain_cloud_above,kz_inv
  real :: pressure
  real :: rh,lsp,cloudh_min,convp,prec
  real :: rhoh(nuvzmax),pinmconv(nzmax)
  real :: pint,tv,tvold,pold,dz1,dz2,dz,ui,vi
  real :: xlon,ylat,xlonr,dzdx,dzdy
  real :: dzdx1,dzdx2,dzdy1,dzdy2,cosf
  real :: uuaux,vvaux,uupolaux,vvpolaux,ddpol,ffpol,wdummy
  real :: uuh(0:nxmax-1,0:nymax-1,nuvzmax)
  real :: vvh(0:nxmax-1,0:nymax-1,nuvzmax)
  real :: pvh(0:nxmax-1,0:nymax-1,nuvzmax)
  real :: wwh(0:nxmax-1,0:nymax-1,nwzmax)
  real :: wzlev(nwzmax),uvwzlev(0:nxmax-1,0:nymax-1,nzmax)
  real,parameter :: const=r_air/ga

  ! NCEP version
  integer :: llev, i

  logical :: init = .true.


  !*************************************************************************
  ! If verttransform is called the first time, initialize heights of the   *
  ! z levels in meter. The heights are the heights of model levels, where  *
  ! u,v,T and qv are given.                                                *
  !*************************************************************************

  if (init) then

  ! Search for a point with high surface pressure (i.e. not above significant topography)
  ! Then, use this point to construct a reference z profile, to be used at all times
  !*****************************************************************************

    do jy=0,nymin1
      do ix=0,nxmin1
        if (ps(ix,jy,1,n).gt.100000.) then
          ixm=ix
          jym=jy
          goto 3
        endif
      end do
    end do
3   continue


    tvold=tt2(ixm,jym,1,n)*(1.+0.378*ew(td2(ixm,jym,1,n),ps(ixm,jym,1,n))/ &
      ps(ixm,jym,1,n))
    pold=ps(ixm,jym,1,n)
    height(1)=0.

    do kz=2,nuvz
      pint=akz(kz)+bkz(kz)*ps(ixm,jym,1,n)
      tv=tth(ixm,jym,kz,n)*(1.+0.608*qvh(ixm,jym,kz,n))

      if (abs(tv-tvold).gt.0.2) then
        height(kz)=height(kz-1)+const*log(pold/pint)* &
        (tv-tvold)/log(tv/tvold)
      else
        height(kz)=height(kz-1)+const*log(pold/pint)*tv
      endif

      tvold=tv
      pold=pint
    end do


  ! Determine highest levels that can be within PBL
  !************************************************

    do kz=1,nz
      if (height(kz).gt.hmixmax) then
        nmixz=kz
        goto 9
      endif
    end do
9   continue
  
  ! Do not repeat initialization of the Cartesian z grid
  !*****************************************************

    init=.false.

  endif


  ! Loop over the whole grid
  !*************************

  do jy=0,nymin1
    do ix=0,nxmin1

  ! NCEP version: find first level above ground
      llev = 0
      do i=1,nuvz
        if (ps(ix,jy,1,n).lt.akz(i)) llev=i
      end do
       llev = llev+1
       if (llev.gt.nuvz-2) llev = nuvz-2
  !     if (llev.eq.nuvz-2) write(*,*) 'verttransform
  !    +WARNING: LLEV eq NUZV-2'
  ! NCEP version


  ! compute height of pressure levels above ground
  !***********************************************

      tvold=tth(ix,jy,llev,n)*(1.+0.608*qvh(ix,jy,llev,n))
      pold=akz(llev)
      wzlev(llev)=0.
      uvwzlev(ix,jy,llev)=0.
      rhoh(llev)=pold/(r_air*tvold)

      do kz=llev+1,nuvz
        pint=akz(kz)+bkz(kz)*ps(ix,jy,1,n)
        tv=tth(ix,jy,kz,n)*(1.+0.608*qvh(ix,jy,kz,n))
        rhoh(kz)=pint/(r_air*tv)

        if (abs(tv-tvold).gt.0.2) then
          uvwzlev(ix,jy,kz)=uvwzlev(ix,jy,kz-1)+const*log(pold/pint)* &
          (tv-tvold)/log(tv/tvold)
        else
          uvwzlev(ix,jy,kz)=uvwzlev(ix,jy,kz-1)+const*log(pold/pint)*tv
        endif
        wzlev(kz)=uvwzlev(ix,jy,kz)

        tvold=tv
        pold=pint
      end do

  ! pinmconv=(h2-h1)/(p2-p1)

      pinmconv(llev)=(uvwzlev(ix,jy,llev+1)-uvwzlev(ix,jy,llev))/ &
           ((aknew(llev+1)+bknew(llev+1)*ps(ix,jy,1,n))- &
           (aknew(llev)+bknew(llev)*ps(ix,jy,1,n)))
      do kz=llev+1,nz-1
        pinmconv(kz)=(uvwzlev(ix,jy,kz+1)-uvwzlev(ix,jy,kz-1))/ &
             ((aknew(kz+1)+bknew(kz+1)*ps(ix,jy,1,n))- &
             (aknew(kz-1)+bknew(kz-1)*ps(ix,jy,1,n)))
      end do
      pinmconv(nz)=(uvwzlev(ix,jy,nz)-uvwzlev(ix,jy,nz-1))/ &
           ((aknew(nz)+bknew(nz)*ps(ix,jy,1,n))- &
           (aknew(nz-1)+bknew(nz-1)*ps(ix,jy,1,n)))


  ! Levels, where u,v,t and q are given
  !************************************

      uu(ix,jy,1,n)=uuh(ix,jy,llev)
      vv(ix,jy,1,n)=vvh(ix,jy,llev)
      tt(ix,jy,1,n)=tth(ix,jy,llev,n)
      qv(ix,jy,1,n)=qvh(ix,jy,llev,n)
  ! IP & SEC, 201812 add clouds
      if (readclouds) then
         clwc(ix,jy,1,n)=clwch(ix,jy,llev,n)
      endif 
      pv(ix,jy,1,n)=pvh(ix,jy,llev)
      rho(ix,jy,1,n)=rhoh(llev)
      pplev(ix,jy,1,n)=akz(llev)
      uu(ix,jy,nz,n)=uuh(ix,jy,nuvz)
      vv(ix,jy,nz,n)=vvh(ix,jy,nuvz)
      tt(ix,jy,nz,n)=tth(ix,jy,nuvz,n)
      qv(ix,jy,nz,n)=qvh(ix,jy,nuvz,n)
  ! IP & SEC, 201812 add clouds
      if (readclouds) then
         clwc(ix,jy,nz,n)=clwch(ix,jy,nuvz,n)
      endif
      pv(ix,jy,nz,n)=pvh(ix,jy,nuvz)
      rho(ix,jy,nz,n)=rhoh(nuvz)
      pplev(ix,jy,nz,n)=akz(nuvz)
      kmin=llev+1
      do iz=2,nz-1
        do kz=kmin,nuvz
          if(height(iz).gt.uvwzlev(ix,jy,nuvz)) then
            uu(ix,jy,iz,n)=uu(ix,jy,nz,n)
            vv(ix,jy,iz,n)=vv(ix,jy,nz,n)
            tt(ix,jy,iz,n)=tt(ix,jy,nz,n)
            qv(ix,jy,iz,n)=qv(ix,jy,nz,n)
  ! IP & SEC, 201812 add clouds
            if (readclouds) then
               clwc(ix,jy,iz,n)=clwc(ix,jy,nz,n)
            endif
            pv(ix,jy,iz,n)=pv(ix,jy,nz,n)
            rho(ix,jy,iz,n)=rho(ix,jy,nz,n)
            pplev(ix,jy,iz,n)=pplev(ix,jy,nz,n)
            goto 30
          endif
          if ((height(iz).gt.uvwzlev(ix,jy,kz-1)).and. &
          (height(iz).le.uvwzlev(ix,jy,kz))) then
            dz1=height(iz)-uvwzlev(ix,jy,kz-1)
            dz2=uvwzlev(ix,jy,kz)-height(iz)
            dz=dz1+dz2
            uu(ix,jy,iz,n)=(uuh(ix,jy,kz-1)*dz2+uuh(ix,jy,kz)*dz1)/dz
            vv(ix,jy,iz,n)=(vvh(ix,jy,kz-1)*dz2+vvh(ix,jy,kz)*dz1)/dz
            tt(ix,jy,iz,n)=(tth(ix,jy,kz-1,n)*dz2 &
            +tth(ix,jy,kz,n)*dz1)/dz
            qv(ix,jy,iz,n)=(qvh(ix,jy,kz-1,n)*dz2 &
            +qvh(ix,jy,kz,n)*dz1)/dz
  ! IP & SEC, 201812 add clouds
            if (readclouds) then
               clwc(ix,jy,iz,n)=(clwch(ix,jy,kz-1,n)*dz2 &
               +clwch(ix,jy,kz,n)*dz1)/dz
            endif
            pv(ix,jy,iz,n)=(pvh(ix,jy,kz-1)*dz2+pvh(ix,jy,kz)*dz1)/dz
            rho(ix,jy,iz,n)=(rhoh(kz-1)*dz2+rhoh(kz)*dz1)/dz
            pplev(ix,jy,iz,n)=(akz(kz-1)*dz2+akz(kz)*dz1)/dz
          endif
        end do
30      continue
      end do


  ! Levels, where w is given
  !*************************

      ww(ix,jy,1,n)=wwh(ix,jy,llev)*pinmconv(llev)
      ww(ix,jy,nz,n)=wwh(ix,jy,nwz)*pinmconv(nz)
      kmin=llev+1
      do iz=2,nz
        do kz=kmin,nwz
          if ((height(iz).gt.wzlev(kz-1)).and. &
          (height(iz).le.wzlev(kz))) then
            dz1=height(iz)-wzlev(kz-1)
            dz2=wzlev(kz)-height(iz)
            dz=dz1+dz2
            ww(ix,jy,iz,n)=(wwh(ix,jy,kz-1)*pinmconv(kz-1)*dz2 &
            +wwh(ix,jy,kz)*pinmconv(kz)*dz1)/dz
          endif
        end do
      end do


  ! Compute density gradients at intermediate levels
  !*************************************************

      drhodz(ix,jy,1,n)=(rho(ix,jy,2,n)-rho(ix,jy,1,n))/ &
           (height(2)-height(1))
      do kz=2,nz-1
        drhodz(ix,jy,kz,n)=(rho(ix,jy,kz+1,n)-rho(ix,jy,kz-1,n))/ &
        (height(kz+1)-height(kz-1))
      end do
      drhodz(ix,jy,nz,n)=drhodz(ix,jy,nz-1,n)

    end do
  end do


  !****************************************************************
  ! Compute slope of eta levels in windward direction and resulting
  ! vertical wind correction
  !****************************************************************

  do jy=1,ny-2
    cosf=cos((real(jy)*dy+ylat0)*pi180)
    do ix=1,nx-2

  ! NCEP version: find first level above ground
      llev = 0
      do i=1,nuvz
       if (ps(ix,jy,1,n).lt.akz(i)) llev=i
      end do
       llev = llev+1
       if (llev.gt.nuvz-2) llev = nuvz-2
  !     if (llev.eq.nuvz-2) write(*,*) 'verttransform
  !    +WARNING: LLEV eq NUZV-2'
  ! NCEP version

      kmin=llev+1
      do iz=2,nz-1

        ui=uu(ix,jy,iz,n)*dxconst/cosf
        vi=vv(ix,jy,iz,n)*dyconst

        do kz=kmin,nz
          if ((height(iz).gt.uvwzlev(ix,jy,kz-1)).and. &
          (height(iz).le.uvwzlev(ix,jy,kz))) then
            dz1=height(iz)-uvwzlev(ix,jy,kz-1)
            dz2=uvwzlev(ix,jy,kz)-height(iz)
            dz=dz1+dz2
            kl=kz-1
            klp=kz
            goto 47
          endif
        end do

47      ix1=ix-1
        jy1=jy-1
        ixp=ix+1
        jyp=jy+1

        dzdx1=(uvwzlev(ixp,jy,kl)-uvwzlev(ix1,jy,kl))/2.
        dzdx2=(uvwzlev(ixp,jy,klp)-uvwzlev(ix1,jy,klp))/2.
        dzdx=(dzdx1*dz2+dzdx2*dz1)/dz

        dzdy1=(uvwzlev(ix,jyp,kl)-uvwzlev(ix,jy1,kl))/2.
        dzdy2=(uvwzlev(ix,jyp,klp)-uvwzlev(ix,jy1,klp))/2.
        dzdy=(dzdy1*dz2+dzdy2*dz1)/dz

        ww(ix,jy,iz,n)=ww(ix,jy,iz,n)+(dzdx*ui+dzdy*vi)

      end do

    end do
  end do


  ! If north pole is in the domain, calculate wind velocities in polar
  ! stereographic coordinates
  !*******************************************************************

  if (nglobal) then
    do jy=int(switchnorthg)-2,nymin1
      ylat=ylat0+real(jy)*dy
      do ix=0,nxmin1
        xlon=xlon0+real(ix)*dx
        do iz=1,nz
          call cc2gll(northpolemap,ylat,xlon,uu(ix,jy,iz,n), &
               vv(ix,jy,iz,n),uupol(ix,jy,iz,n), &
               vvpol(ix,jy,iz,n))
        end do
      end do
    end do


    do iz=1,nz

  ! CALCULATE FFPOL, DDPOL FOR CENTRAL GRID POINT
      xlon=xlon0+real(nx/2-1)*dx
      xlonr=xlon*pi/180.
      ffpol=sqrt(uu(nx/2-1,nymin1,iz,n)**2+vv(nx/2-1,nymin1,iz,n)**2)
      if (vv(nx/2-1,nymin1,iz,n).lt.0.) then
        ddpol=atan(uu(nx/2-1,nymin1,iz,n)/vv(nx/2-1,nymin1,iz,n))-xlonr
      elseif (vv(nx/2-1,nymin1,iz,n).gt.0.) then
        ddpol=pi+atan(uu(nx/2-1,nymin1,iz,n)/ &
        vv(nx/2-1,nymin1,iz,n))-xlonr
      else
        ddpol=pi/2-xlonr
      endif
      if(ddpol.lt.0.) ddpol=2.0*pi+ddpol
      if(ddpol.gt.2.0*pi) ddpol=ddpol-2.0*pi

  ! CALCULATE U,V FOR 180 DEG, TRANSFORM TO POLAR STEREOGRAPHIC GRID
      xlon=180.0
      xlonr=xlon*pi/180.
      ylat=90.0
      uuaux=-ffpol*sin(xlonr+ddpol)
      vvaux=-ffpol*cos(xlonr+ddpol)
      call cc2gll(northpolemap,ylat,xlon,uuaux,vvaux,uupolaux,vvpolaux)
      jy=nymin1
      do ix=0,nxmin1
        uupol(ix,jy,iz,n)=uupolaux
        vvpol(ix,jy,iz,n)=vvpolaux
      end do
    end do


  ! Fix: Set W at pole to the zonally averaged W of the next equator-
  ! ward parallel of latitude

    do iz=1,nz
      wdummy=0.
      jy=ny-2
      do ix=0,nxmin1
        wdummy=wdummy+ww(ix,jy,iz,n)
      end do
      wdummy=wdummy/real(nx)
      jy=nymin1
      do ix=0,nxmin1
        ww(ix,jy,iz,n)=wdummy
      end do
    end do

  endif


  ! If south pole is in the domain, calculate wind velocities in polar
  ! stereographic coordinates
  !*******************************************************************

  if (sglobal) then
    do jy=0,int(switchsouthg)+3
      ylat=ylat0+real(jy)*dy
      do ix=0,nxmin1
        xlon=xlon0+real(ix)*dx
        do iz=1,nz
          call cc2gll(southpolemap,ylat,xlon,uu(ix,jy,iz,n), &
          vv(ix,jy,iz,n),uupol(ix,jy,iz,n),vvpol(ix,jy,iz,n))
        end do
      end do
    end do

    do iz=1,nz

  ! CALCULATE FFPOL, DDPOL FOR CENTRAL GRID POINT
      xlon=xlon0+real(nx/2-1)*dx
      xlonr=xlon*pi/180.
      ffpol=sqrt(uu(nx/2-1,0,iz,n)**2+vv(nx/2-1,0,iz,n)**2)
      if(vv(nx/2-1,0,iz,n).lt.0.) then
        ddpol=atan(uu(nx/2-1,0,iz,n)/vv(nx/2-1,0,iz,n))+xlonr
      elseif (vv(nx/2-1,0,iz,n).gt.0.) then
        ddpol=pi+atan(uu(nx/2-1,0,iz,n)/vv(nx/2-1,0,iz,n))-xlonr
      else
        ddpol=pi/2-xlonr
      endif
      if(ddpol.lt.0.) ddpol=2.0*pi+ddpol
      if(ddpol.gt.2.0*pi) ddpol=ddpol-2.0*pi

  ! CALCULATE U,V FOR 180 DEG, TRANSFORM TO POLAR STEREOGRAPHIC GRID
      xlon=180.0
      xlonr=xlon*pi/180.
      ylat=-90.0
      uuaux=+ffpol*sin(xlonr-ddpol)
      vvaux=-ffpol*cos(xlonr-ddpol)
      call cc2gll(northpolemap,ylat,xlon,uuaux,vvaux,uupolaux,vvpolaux)

      jy=0
      do ix=0,nxmin1
        uupol(ix,jy,iz,n)=uupolaux
        vvpol(ix,jy,iz,n)=vvpolaux
      end do
    end do


  ! Fix: Set W at pole to the zonally averaged W of the next equator-
  ! ward parallel of latitude

    do iz=1,nz
      wdummy=0.
      jy=1
      do ix=0,nxmin1
        wdummy=wdummy+ww(ix,jy,iz,n)
      end do
      wdummy=wdummy/real(nx)
      jy=0
      do ix=0,nxmin1
        ww(ix,jy,iz,n)=wdummy
      end do
    end do
  endif



  !***********************************************************************************
  ! IP & SEC, 201812 GFS clouds read
  if (readclouds) then
  ! The method is loops all grids vertically and constructs the 3D matrix for clouds
  ! Cloud top and cloud bottom gid cells are assigned as well as the total column
  ! cloud water. For precipitating grids, the type and whether it is in or below
  ! cloud scavenging are assigned with numbers 2-5 (following the old metod).
  ! Distinction is done for lsp and convp though they are treated the same in regards
  ! to scavenging. Also clouds that are not precipitating are defined which may be
  ! to include future cloud processing by non-precipitating-clouds.
  !***********************************************************************************
    write(*,*) 'Global NCEP fields: using cloud water'
    clw(:,:,:,n)=0.0
    ctwc(:,:,n)=0.0
    clouds(:,:,:,n)=0
  ! If water/ice are read separately into clwc and ciwc, store sum in clwc
    do jy=0,nymin1
      do ix=0,nxmin1
        lsp=lsprec(ix,jy,1,n)
        convp=convprec(ix,jy,1,n)
        prec=lsp+convp
  ! Find clouds in the vertical
        do kz=1, nz-1 !go from top to bottom
          if (clwc(ix,jy,kz,n).gt.0) then
  ! assuming rho is in kg/m3 and hz in m gives: kg/kg * kg/m3 *m3/kg /m = m2/m3
            clw(ix,jy,kz,n)=(clwc(ix,jy,kz,n)*rho(ix,jy,kz,n))*(height(kz+1)-height(kz))
            ctwc(ix,jy,n) = ctwc(ix,jy,n)+clw(ix,jy,kz,n)
            cloudh_min=min(height(kz+1),height(kz))
          endif
        end do

  ! If Precipitation. Define removal type in the vertical
        if ((lsp.gt.0.01).or.(convp.gt.0.01)) then ! cloud and precipitation

          do kz=nz,2,-1 !go Bottom up!
            if (clw(ix,jy,kz,n).gt. 0) then ! is in cloud
              cloudsh(ix,jy,n)=cloudsh(ix,jy,n)+height(kz)-height(kz-1)
              clouds(ix,jy,kz,n)=1                               ! is a cloud
              if (lsp.ge.convp) then
                clouds(ix,jy,kz,n)=3                             ! lsp in-cloud
              else
                clouds(ix,jy,kz,n)=2                             ! convp in-cloud
              endif                                              ! convective or large scale
            elseif((clw(ix,jy,kz,n).le.0) .and. (cloudh_min.ge.height(kz))) then ! is below cloud
              if (lsp.ge.convp) then
                clouds(ix,jy,kz,n)=5                             ! lsp dominated washout
              else
                clouds(ix,jy,kz,n)=4                             ! convp dominated washout
              endif                                              ! convective or large scale
            endif

            if (height(kz).ge. 19000) then                        ! set a max height for removal
              clouds(ix,jy,kz,n)=0
            endif !clw>0
          end do !nz
        endif ! precipitation
      end do
    end do
  else
  write(*,*) 'Global NCEP fields: using cloud water from Parameterization'
  !   write (*,*) 'initializing clouds, n:',n,nymin1,nxmin1,nz
  !   create a cloud and rainout/washout field, clouds occur where rh>80%
  !   total cloudheight is stored at level 0
  do jy=0,nymin1
    do ix=0,nxmin1
      rain_cloud_above=0
      lsp=lsprec(ix,jy,1,n)
      convp=convprec(ix,jy,1,n)
      cloudsh(ix,jy,n)=0
      do kz_inv=1,nz-1
         kz=nz-kz_inv+1
         pressure=rho(ix,jy,kz,n)*r_air*tt(ix,jy,kz,n)
         rh=qv(ix,jy,kz,n)/f_qvsat(pressure,tt(ix,jy,kz,n))
         clouds(ix,jy,kz,n)=0
         if (rh.gt.0.8) then ! in cloud
           if ((lsp.gt.0.01).or.(convp.gt.0.01)) then ! cloud and precipitation
              rain_cloud_above=1
              cloudsh(ix,jy,n)=cloudsh(ix,jy,n)+height(kz)-height(kz-1)
              if (lsp.ge.convp) then
                 clouds(ix,jy,kz,n)=3 ! lsp dominated rainout
              else
                 clouds(ix,jy,kz,n)=2 ! convp dominated rainout
              endif
           else ! no precipitation
             clouds(ix,jy,kz,n)=1 ! cloud
           endif
         else ! no cloud
           if (rain_cloud_above.eq.1) then ! scavenging
             if (lsp.ge.convp) then
               clouds(ix,jy,kz,n)=5 ! lsp dominated washout
             else
               clouds(ix,jy,kz,n)=4 ! convp dominated washout
             endif
           endif
         endif
      end do
    end do
  end do
  endif  ! IP & SEC 201812, GFS clouds read
end subroutine verttransform_gfs

subroutine verttransform_nests(n,uuhn,vvhn,wwhn,pvhn)
  !                            i   i    i    i   i
  !*****************************************************************************
  !                                                                            *
  !     This subroutine transforms temperature, dew point temperature and      *
  !     wind components from eta to meter coordinates.                         *
  !     The vertical wind component is transformed from Pa/s to m/s using      *
  !     the conversion factor pinmconv.                                        *
  !     In addition, this routine calculates vertical density gradients        *
  !     needed for the parameterization of the turbulent velocities.           *
  !     It is similar to verttransform, but makes the transformations for      *
  !     the nested grids.                                                      *
  !                                                                            *
  !     Author: A. Stohl, G. Wotawa                                            *
  !                                                                            *
  !     12 August 1996                                                         *
  !     Update: 16 January 1998                                                *
  !                                                                            *
  !     Major update: 17 February 1999                                         *
  !     by G. Wotawa                                                           *
  !                                                                            *
  !     - Vertical levels for u, v and w are put together                      *
  !     - Slope correction for vertical velocity: Modification of calculation  *
  !       procedure                                                            *
  !                                                                            *
  !*****************************************************************************
  !  Changes, Bernd C. Krueger, Feb. 2001:       (marked "C-cv")
  !   Variables tthn and qvhn (on eta coordinates) from common block
  !*****************************************************************************
  ! Sabine Eckhardt, March 2007
  ! add the variable cloud for use with scavenging - descr. in com_mod
  !*****************************************************************************
  ! ESO, 2016
  ! -note that divide-by-zero occurs when nxmaxn,nymaxn etc. are larger than 
  !  the actual field dimensions
  !*****************************************************************************
  ! Date: 2017-05-30 modification of a bug in ew. Don Morton (CTBTO project)   *
  !*****************************************************************************
  !                                                                            *
  ! Variables:                                                                 *
  ! nxn,nyn,nuvz,nwz                field dimensions in x,y and z direction    *
  ! uun                             wind components in x-direction [m/s]       *
  ! vvn                             wind components in y-direction [m/s]       *
  ! wwn                             wind components in z-direction [deltaeta/s]*
  ! ttn                             temperature [K]                            *
  ! pvn                             potential vorticity (pvu)                  *
  ! psn                             surface pressure [Pa]                      *
  !                                                                            *
  !*****************************************************************************

  use par_mod
  use com_mod
  use qvsat_mod

  implicit none

  real,intent(in),dimension(0:nxmaxn-1,0:nymaxn-1,nuvzmax,maxnests) :: uuhn,vvhn,pvhn
  real,intent(in),dimension(0:nxmaxn-1,0:nymaxn-1,nwzmax,maxnests) :: wwhn

  real,dimension(0:nxmaxn-1,0:nymaxn-1,nuvzmax) :: rhohn,uvzlev,wzlev
  real,dimension(0:nxmaxn-1,0:nymaxn-1,nzmax) :: pinmconv
  real,dimension(0:nxmaxn-1,0:nymaxn-1) :: tvold,pold,pint,tv
  real,dimension(0:nymaxn-1) :: cosf

  integer,dimension(0:nxmaxn-1,0:nymaxn-1) :: rain_cloud_above, idx

  integer :: ix,jy,kz,iz,n,l,kmin,kl,klp,ix1,jy1,ixp,jyp,kz_inv
  real :: pressure,rh,lsp,convp,cloudh_min,prec

  real :: dz1,dz2,dz
  real :: dzdx,dzdy
  real :: dzdx1,dzdx2,dzdy1,dzdy2
  real,parameter :: const=r_air/ga
  real :: tot_cloud_h
  integer :: nxm1, nym1

  !  real,parameter :: precmin = 0.002 ! minimum prec in mm/h for cloud diagnostics

  ! Loop over all nests
  !********************

  do l=1,numbnests
    nxm1=nxn(l)-1 
    nym1=nyn(l)-1 

  ! Loop over the whole grid
  !*************************

    do jy=0,nyn(l)-1
      do ix=0,nxn(l)-1
        tvold(ix,jy)=tt2n(ix,jy,1,n,l)*(1.+0.378*ew(td2n(ix,jy,1,n,l),psn(ix,jy,1,n,l))/ &
             psn(ix,jy,1,n,l))
      end do
    end do

    pold(0:nxm1,0:nym1)=psn(0:nxm1,0:nym1,1,n,l)
    uvzlev(0:nxm1,0:nym1,1)=0.
    wzlev(0:nxm1,0:nym1,1)=0.
    rhohn(0:nxm1,0:nym1,1)=pold(0:nxm1,0:nym1)/(r_air*tvold(0:nxm1,0:nym1))

  ! Compute heights of eta levels
  !******************************

    do kz=2,nuvz
      pint(0:nxm1,0:nym1)=akz(kz)+bkz(kz)*psn(0:nxm1,0:nym1,1,n,l)
      tv(0:nxm1,0:nym1)=tthn(0:nxm1,0:nym1,kz,n,l)*(1.+0.608*qvhn(0:nxm1,0:nym1,kz,n,l))
      rhohn(0:nxm1,0:nym1,kz)=pint(0:nxm1,0:nym1)/(r_air*tv(0:nxm1,0:nym1))

      where (abs(tv(0:nxm1,0:nym1)-tvold(0:nxm1,0:nym1)).gt.0.2) 
        uvzlev(0:nxm1,0:nym1,kz)=uvzlev(0:nxm1,0:nym1,kz-1)+const*&
             &log(pold(0:nxm1,0:nym1)/pint(0:nxm1,0:nym1))* &
             (tv(0:nxm1,0:nym1)-tvold(0:nxm1,0:nym1))/&
             &log(tv(0:nxm1,0:nym1)/tvold(0:nxm1,0:nym1))
      elsewhere
        uvzlev(0:nxm1,0:nym1,kz)=uvzlev(0:nxm1,0:nym1,kz-1)+const*&
             &log(pold(0:nxm1,0:nym1)/pint(0:nxm1,0:nym1))*tv(0:nxm1,0:nym1)
      endwhere

      tvold(0:nxm1,0:nym1)=tv(0:nxm1,0:nym1)
      pold(0:nxm1,0:nym1)=pint(0:nxm1,0:nym1)

    end do

    do kz=2,nwz-1
      wzlev(0:nxm1,0:nym1,kz)=(uvzlev(0:nxm1,0:nym1,kz+1)+uvzlev(0:nxm1,0:nym1,kz))/2.
    end do
    wzlev(0:nxm1,0:nym1,nwz)=wzlev(0:nxm1,0:nym1,nwz-1)+ &
         uvzlev(0:nxm1,0:nym1,nuvz)-uvzlev(0:nxm1,0:nym1,nuvz-1)


    pinmconv(0:nxm1,0:nym1,1)=(uvzlev(0:nxm1,0:nym1,2))/ &
         ((aknew(2)+bknew(2)*psn(0:nxm1,0:nym1,1,n,l))- &
         (aknew(1)+bknew(1)*psn(0:nxm1,0:nym1,1,n,l)))
    do kz=2,nz-1
      pinmconv(0:nxm1,0:nym1,kz)=(uvzlev(0:nxm1,0:nym1,kz+1)-uvzlev(0:nxm1,0:nym1,kz-1))/ &
           ((aknew(kz+1)+bknew(kz+1)*psn(0:nxm1,0:nym1,1,n,l))- &
           (aknew(kz-1)+bknew(kz-1)*psn(0:nxm1,0:nym1,1,n,l)))
    end do
    pinmconv(0:nxm1,0:nym1,nz)=(uvzlev(0:nxm1,0:nym1,nz)-uvzlev(0:nxm1,0:nym1,nz-1))/ &
         ((aknew(nz)+bknew(nz)*psn(0:nxm1,0:nym1,1,n,l))- &
         (aknew(nz-1)+bknew(nz-1)*psn(0:nxm1,0:nym1,1,n,l)))

  ! Levels, where u,v,t and q are given
  !************************************

    uun(0:nxm1,0:nym1,1,n,l)=uuhn(0:nxm1,0:nym1,1,l)
    vvn(0:nxm1,0:nym1,1,n,l)=vvhn(0:nxm1,0:nym1,1,l)
    ttn(0:nxm1,0:nym1,1,n,l)=tthn(0:nxm1,0:nym1,1,n,l)
    qvn(0:nxm1,0:nym1,1,n,l)=qvhn(0:nxm1,0:nym1,1,n,l)
    if (readclouds_nest(l)) then
      clwcn(0:nxm1,0:nym1,1,n,l)=clwchn(0:nxm1,0:nym1,1,n,l)
      if (.not.sumclouds_nest(l)) ciwcn(0:nxm1,0:nym1,1,n,l)=ciwchn(0:nxm1,0:nym1,1,n,l)
    end if
    pvn(0:nxm1,0:nym1,1,n,l)=pvhn(0:nxm1,0:nym1,1,l)
    rhon(0:nxm1,0:nym1,1,n,l)=rhohn(0:nxm1,0:nym1,1)

    uun(0:nxm1,0:nym1,nz,n,l)=uuhn(0:nxm1,0:nym1,nuvz,l)
    vvn(0:nxm1,0:nym1,nz,n,l)=vvhn(0:nxm1,0:nym1,nuvz,l)
    ttn(0:nxm1,0:nym1,nz,n,l)=tthn(0:nxm1,0:nym1,nuvz,n,l)
    qvn(0:nxm1,0:nym1,nz,n,l)=qvhn(0:nxm1,0:nym1,nuvz,n,l)
    if (readclouds_nest(l)) then
      clwcn(0:nxm1,0:nym1,nz,n,l)=clwchn(0:nxm1,0:nym1,nuvz,n,l)
      if (.not.sumclouds_nest(l)) ciwcn(0:nxm1,0:nym1,nz,n,l)=ciwchn(0:nxm1,0:nym1,nuvz,n,l)
    end if
    pvn(0:nxm1,0:nym1,nz,n,l)=pvhn(0:nxm1,0:nym1,nuvz,l)
    rhon(0:nxm1,0:nym1,nz,n,l)=rhohn(0:nxm1,0:nym1,nuvz)


    kmin=2
    idx=kmin
    do iz=2,nz-1
      do jy=0,nyn(l)-1
        do ix=0,nxn(l)-1
          if(height(iz).gt.uvzlev(ix,jy,nuvz)) then
            uun(ix,jy,iz,n,l)=uun(ix,jy,nz,n,l)
            vvn(ix,jy,iz,n,l)=vvn(ix,jy,nz,n,l)
            ttn(ix,jy,iz,n,l)=ttn(ix,jy,nz,n,l)
            qvn(ix,jy,iz,n,l)=qvn(ix,jy,nz,n,l)
  !hg adding the cloud water
            if (readclouds_nest(l)) then
              clwcn(ix,jy,iz,n,l)=clwcn(ix,jy,nz,n,l)
              if (.not.sumclouds_nest(l)) ciwcn(ix,jy,iz,n,l)=ciwcn(ix,jy,nz,n,l)
            end if
  !hg
            pvn(ix,jy,iz,n,l)=pvn(ix,jy,nz,n,l)
            rhon(ix,jy,iz,n,l)=rhon(ix,jy,nz,n,l)
          else
            innuvz: do kz=idx(ix,jy),nuvz
              if (idx(ix,jy) .le. kz .and. (height(iz).gt.uvzlev(ix,jy,kz-1)).and. &
                   (height(iz).le.uvzlev(ix,jy,kz))) then
                idx(ix,jy)=kz
                exit innuvz
              endif
            enddo innuvz
          endif
        enddo
      enddo
      do jy=0,nyn(l)-1
        do ix=0,nxn(l)-1
          if(height(iz).le.uvzlev(ix,jy,nuvz)) then
            kz=idx(ix,jy)
            dz1=height(iz)-uvzlev(ix,jy,kz-1)
            dz2=uvzlev(ix,jy,kz)-height(iz)
            dz=dz1+dz2
            uun(ix,jy,iz,n,l)=(uuhn(ix,jy,kz-1,l)*dz2+uuhn(ix,jy,kz,l)*dz1)/dz
            vvn(ix,jy,iz,n,l)=(vvhn(ix,jy,kz-1,l)*dz2+vvhn(ix,jy,kz,l)*dz1)/dz
            ttn(ix,jy,iz,n,l)=(tthn(ix,jy,kz-1,n,l)*dz2 &
                 +tthn(ix,jy,kz,n,l)*dz1)/dz
            qvn(ix,jy,iz,n,l)=(qvhn(ix,jy,kz-1,n,l)*dz2 &
                 +qvhn(ix,jy,kz,n,l)*dz1)/dz
  !hg adding the cloud water
            if (readclouds_nest(l)) then
              clwcn(ix,jy,iz,n,l)=(clwchn(ix,jy,kz-1,n,l)*dz2+clwchn(ix,jy,kz,n,l)*dz1)/dz
              if (.not.sumclouds_nest(l)) &
                   &ciwcn(ix,jy,iz,n,l)=(ciwchn(ix,jy,kz-1,n,l)*dz2+ciwchn(ix,jy,kz,n,l)*dz1)/dz
            end if
  !hg
            pvn(ix,jy,iz,n,l)=(pvhn(ix,jy,kz-1,l)*dz2+pvhn(ix,jy,kz,l)*dz1)/dz
            rhon(ix,jy,iz,n,l)=(rhohn(ix,jy,kz-1)*dz2+rhohn(ix,jy,kz)*dz1)/dz
          endif
        enddo
      enddo
    enddo

  ! Levels, where w is given
  !*************************

    wwn(0:nxm1,0:nym1,1,n,l)=wwhn(0:nxm1,0:nym1,1,l)*pinmconv(0:nxm1,0:nym1,1)
    wwn(0:nxm1,0:nym1,nz,n,l)=wwhn(0:nxm1,0:nym1,nwz,l)*pinmconv(0:nxm1,0:nym1,nz)
    kmin=2
    idx=kmin
    do iz=2,nz
      do jy=0,nym1
        do ix=0,nxm1
          inn:         do kz=idx(ix,jy),nwz
            if(idx(ix,jy) .le. kz .and. height(iz).gt.wzlev(ix,jy,kz-1).and. &
                 height(iz).le.wzlev(ix,jy,kz)) then
              idx(ix,jy)=kz
              exit inn
            endif
          enddo inn
        enddo
      enddo
      do jy=0,nym1
        do ix=0,nxm1
          kz=idx(ix,jy)
          dz1=height(iz)-wzlev(ix,jy,kz-1)
          dz2=wzlev(ix,jy,kz)-height(iz)
          dz=dz1+dz2
          wwn(ix,jy,iz,n,l)=(wwhn(ix,jy,kz-1,l)*pinmconv(ix,jy,kz-1)*dz2 &
               +wwhn(ix,jy,kz,l)*pinmconv(ix,jy,kz)*dz1)/dz
        enddo
      enddo
    enddo

  ! Compute density gradients at intermediate levels
  !*************************************************

    drhodzn(0:nxm1,0:nym1,1,n,l)=(rhon(0:nxm1,0:nym1,2,n,l)-rhon(0:nxm1,0:nym1,1,n,l))/ &
         (height(2)-height(1))
    do kz=2,nz-1
      drhodzn(0:nxm1,0:nym1,kz,n,l)=(rhon(0:nxm1,0:nym1,kz+1,n,l)-rhon(0:nxm1,0:nym1,kz-1,n,l))/ &
           (height(kz+1)-height(kz-1))
    end do
    drhodzn(0:nxm1,0:nym1,nz,n,l)=drhodzn(0:nxm1,0:nym1,nz-1,n,l)


  !****************************************************************
  ! Compute slope of eta levels in windward direction and resulting
  ! vertical wind correction
  !****************************************************************

    do jy=1,nyn(l)-2
      cosf(jy)=1./cos((real(jy)*dyn(l)+ylat0n(l))*pi180)
    end do

    kmin=2
    idx=kmin
    do iz=2,nz-1
      do jy=1,nyn(l)-2
        do ix=1,nxn(l)-2

          inneta: do kz=idx(ix,jy),nz
            if (idx(ix,jy) .le. kz .and. (height(iz).gt.uvzlev(ix,jy,kz-1)).and. &
                 (height(iz).le.uvzlev(ix,jy,kz))) then
              idx(ix,jy)=kz
              exit inneta
            endif
          enddo inneta
        enddo
      enddo

      do jy=1,nyn(l)-2
        do ix=1,nxn(l)-2
          kz=idx(ix,jy)
          dz1=height(iz)-uvzlev(ix,jy,kz-1)
          dz2=uvzlev(ix,jy,kz)-height(iz)
          dz=dz1+dz2
          ix1=ix-1
          jy1=jy-1
          ixp=ix+1
          jyp=jy+1

          dzdx1=(uvzlev(ixp,jy,kz-1)-uvzlev(ix1,jy,kz-1))/2.
          dzdx2=(uvzlev(ixp,jy,kz)-uvzlev(ix1,jy,kz))/2.
          dzdx=(dzdx1*dz2+dzdx2*dz1)/dz

          dzdy1=(uvzlev(ix,jyp,kz-1)-uvzlev(ix,jy1,kz-1))/2.
          dzdy2=(uvzlev(ix,jyp,kz)-uvzlev(ix,jy1,kz))/2.
          dzdy=(dzdy1*dz2+dzdy2*dz1)/dz

          wwn(ix,jy,iz,n,l)=wwn(ix,jy,iz,n,l)+(dzdx*uun(ix,jy,iz,n,l)*dxconst*xresoln(l)*cosf(jy)+&
               &dzdy*vvn(ix,jy,iz,n,l)*dyconst*yresoln(l))

        end do

      end do
    end do


  !***********************************************************************************  
    if (readclouds_nest(l)) then !HG METHOD
  ! The method is loops all grids vertically and constructs the 3D matrix for clouds
  ! Cloud top and cloud bottom gid cells are assigned as well as the total column 
  ! cloud water. For precipitating grids, the type and whether it is in or below 
  ! cloud scavenging are assigned with numbers 2-5 (following the old metod).
  ! Distinction is done for lsp and convp though they are treated the same in regards
  ! to scavenging. Also clouds that are not precipitating are defined which may be 
  ! to include future cloud processing by non-precipitating-clouds. 
  !***********************************************************************************
      write(*,*) 'Nested ECMWF fields: using cloud water'
      clwn(0:nxm1,0:nym1,:,n,l)=0.0
  !    icloud_stats(0:nxm1,0:nym1,:,n)=0.0
      ctwcn(0:nxm1,0:nym1,n,l)=0.0
      cloudsn(0:nxm1,0:nym1,:,n,l)=0
  ! If water/ice are read separately into clwc and ciwc, store sum in clwcn
      if (.not.sumclouds_nest(l)) then 
        clwcn(0:nxm1,0:nym1,:,n,l) = clwcn(0:nxm1,0:nym1,:,n,l) + ciwcn(0:nxm1,0:nym1,:,n,l)
      end if
      do jy=0,nym1
        do ix=0,nxm1
          lsp=lsprecn(ix,jy,1,n,l)
          convp=convprecn(ix,jy,1,n,l)
          prec=lsp+convp
          tot_cloud_h=0
  ! Find clouds in the vertical
          do kz=1, nz-1 !go from top to bottom
            if (clwcn(ix,jy,kz,n,l).gt.0) then      
  ! assuming rho is in kg/m3 and hz in m gives: kg/kg * kg/m3 *m3/kg /m = m2/m3 
              clwn(ix,jy,kz,n,l)=(clwcn(ix,jy,kz,n,l)*rhon(ix,jy,kz,n,l))*(height(kz+1)-height(kz))
              tot_cloud_h=tot_cloud_h+(height(kz+1)-height(kz)) 
              ctwcn(ix,jy,n,l) = ctwcn(ix,jy,n,l)+clwn(ix,jy,kz,n,l)
  !            icloud_stats(ix,jy,4,n)= icloud_stats(ix,jy,4,n)+clw(ix,jy,kz,n)          ! Column cloud water [m3/m3]
  !           icloud_stats(ix,jy,3,n)= min(height(kz+1),height(kz))                     ! Cloud BOT height stats      [m]
              cloudh_min=min(height(kz+1),height(kz))
            endif
          end do

  ! If Precipitation. Define removal type in the vertical
          if ((lsp.gt.0.01).or.(convp.gt.0.01)) then ! cloud and precipitation

            do kz=nz,1,-1 !go Bottom up!
              if (clwn(ix,jy,kz,n,l).gt.0.0) then ! is in cloud
                cloudshn(ix,jy,n,l)=cloudshn(ix,jy,n,l)+height(kz)-height(kz-1) 
                cloudsn(ix,jy,kz,n,l)=1                               ! is a cloud
                if (lsp.ge.convp) then
                  cloudsn(ix,jy,kz,n,l)=3                            ! lsp in-cloud
                else
                  cloudsn(ix,jy,kz,n,l)=2                             ! convp in-cloud
                endif                                              ! convective or large scale
              else if((clwn(ix,jy,kz,n,l).le.0.0) .and. (cloudh_min.ge.height(kz))) then ! is below cloud
                if (lsp.ge.convp) then
                  cloudsn(ix,jy,kz,n,l)=5                             ! lsp dominated washout
                else
                  cloudsn(ix,jy,kz,n,l)=4                             ! convp dominated washout
                endif                                              ! convective or large scale 
              endif

              if (height(kz).ge. 19000) then                        ! set a max height for removal
                cloudsn(ix,jy,kz,n,l)=0
              endif !clw>0
            end do !nz
          endif ! precipitation
        end do
      end do
  !********************************************************************
    else ! old method:
  !********************************************************************
      write(*,*) 'Nested fields: using cloud water from Parameterization'
      do jy=0,nyn(l)-1
        do ix=0,nxn(l)-1
          rain_cloud_above(ix,jy)=0
          lsp=lsprecn(ix,jy,1,n,l)
          convp=convprecn(ix,jy,1,n,l)
          cloudshn(ix,jy,n,l)=0
          do kz_inv=1,nz-1
            kz=nz-kz_inv+1
            pressure=rhon(ix,jy,kz,n,l)*r_air*ttn(ix,jy,kz,n,l)
            rh=qvn(ix,jy,kz,n,l)/f_qvsat(pressure,ttn(ix,jy,kz,n,l))
            cloudsn(ix,jy,kz,n,l)=0
            if (rh.gt.0.8) then ! in cloud
              if ((lsp.gt.0.01).or.(convp.gt.0.01)) then
                rain_cloud_above(ix,jy)=1
                cloudshn(ix,jy,n,l)=cloudshn(ix,jy,n,l)+ &
                     height(kz)-height(kz-1)
                if (lsp.ge.convp) then
                  cloudsn(ix,jy,kz,n,l)=3 ! lsp dominated rainout
                else
                  cloudsn(ix,jy,kz,n,l)=2 ! convp dominated rainout
                endif
              else ! no precipitation
                cloudsn(ix,jy,kz,n,l)=1 ! cloud
              endif
            else ! no cloud
              if (rain_cloud_above(ix,jy).eq.1) then ! scavenging
                if (lsp.ge.convp) then
                  cloudsn(ix,jy,kz,n,l)=5 ! lsp dominated washout
                else
                  cloudsn(ix,jy,kz,n,l)=4 ! convp dominated washout
                endif
              endif
            endif
          end do
        end do
      end do
    end if

  end do ! end loop over nests
end subroutine verttransform_nests

subroutine writeprecip(itime,imem)

  !*****************************************************************************
  !                                                                            *
  !  This routine produces a file containing total precipitation for each      *
  !  releases point.                                                           *
  !                                                                            *
  !     Author: S. Eckhardt                                                    * 
  !     7 Mai 2017                                                             *
  !*****************************************************************************

  use point_mod
  use par_mod
  use com_mod
  use date_mod

  implicit none

  integer :: jjjjmmdd,ihmmss,itime,i
  real(kind=dp) :: jul
  character :: adate*8,atime*6

  integer :: ix,jy,imem
  real :: xp1,yp1

  
  if (itime.eq.0) then
      open(unitprecip,file=path(2)(1:length(2))//'wetscav_precip.txt', &
       form='formatted',err=998)
  else
      open(unitprecip,file=path(2)(1:length(2))//'wetscav_precip.txt', &
       ACCESS='APPEND',form='formatted',err=998)
  endif

  jul=bdate+real(itime,kind=dp)/86400._dp
  call caldate(jul,jjjjmmdd,ihmmss)
  write(adate,'(i8.8)') jjjjmmdd
  write(atime,'(i6.6)') ihmmss

  do i=1,numpoint
    xp1=xpoint1(i)*dx+xlon0 !lat, long (real) coord
    yp1=ypoint1(i)*dy+ylat0 !lat, long (real) coord
    ix=int((xpoint1(i)+xpoint2(i))/2.)
    jy=int((ypoint1(i)+ypoint2(i))/2.)
    write(unitprecip,*)  jjjjmmdd, ihmmss, & 
           xp1,yp1,lsprec(ix,jy,1,imem),convprec(ix,jy,1,imem) !time is the same as in the ECMWF windfield
! units mm/h, valid for the time given in the windfield
  end do

  close(unitprecip)

  return


998   write(*,*) ' #### FLEXPART MODEL ERROR!   THE FILE         #### '
  write(*,*) ' #### '//path(2)(1:length(2))//'header_txt'//' #### '
  write(*,*) ' #### CANNOT BE OPENED. IF A FILE WITH THIS    #### '
  write(*,*) ' #### NAME ALREADY EXISTS, DELETE IT AND START #### '
  write(*,*) ' #### THE PROGRAM AGAIN.                       #### '
  stop
end subroutine writeprecip

end module windfields_mod