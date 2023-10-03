! SPDX-FileCopyrightText: FLEXPART 1998-2019, see flexpart_license.txt
! SPDX-License-Identifier: GPL-3.0-or-later

  !*****************************************************************************
  !                                                                            *
  ! This module stores all meteorological input data and contains routines     *
  ! reading and allocating this data                                           *
  !                                                                            *
  ! L. Bakels 2023: Dynamical allocation of all windfields                     *
  !                                                                            *
  !*****************************************************************************

module windfields_mod
  use par_mod
  use com_mod
  use point_mod
  use pbl_profile_mod

  implicit none

 !******************************************************************************
 ! Variables associated with the ECMWF meteorological input data ("wind fields")
 !******************************************************************************

  integer ::            &
    numbwf                ! actual number of wind fields

  integer,allocatable,dimension(:) ::            &
    wftime         ! times relative to beginning time of wind fields [s]

  character(len=255),allocatable,dimension(:) :: &
    wfname        ! file names of wind fields

  ! Nested equivalents
  !*******************
  character(len=255),allocatable,dimension(:,:) :: &
    wfnamen                                          ! nested wind field names

  !Windfield parameters
  !********************
  integer :: nxmax,nymax,nuvzmax,nwzmax,nzmax !Size of windfield

  ! Fixed fields, unchangeable with time
  !*************************************
  real, allocatable,dimension(:,:) :: &
    oro,                              & ! orography of the ECMWF model
    excessoro,                        & ! excess orography mother domain
    lsm                                 ! land sea mask of the ECMWF model

  ! Nested fields, unchangeable with time
  !**************************************
  real, allocatable,dimension(:,:,:) :: &
    oron,                               & ! orography of the ECMWF model
    excessoron,                         & ! excess orography mother domain
    lsmn                                  ! land sea mask of the ECMWF model

  ! 3d fields necessary for eta coordinates option
  !************************************************
  real, allocatable,dimension(:,:,:,:) :: &
    uueta,vveta,                          & ! wind components on half model levels in x and y direction [m/s]
    wweta,                                & ! wind component on model levels in z direction [eta/s]
    uupoleta,vvpoleta,                    & ! wind components on half model levels in polar stereographic projection [m/s]
    tteta,                                & ! temperature data on half model levels [K]
    pveta,                                & ! potential vorticity on half model levels
    rhoeta,                               & ! air density on half model levels [kg/m3]
    prseta,                               & ! air pressure on half model levels
    drhodzeta,                            & ! vertical air density gradient on half model levels [kg/m2]
    !tvirtual,                             & ! Virtual temperature on half model levels
    etauvheight,etawheight                  ! Saved half model and model heights for ETA coordinate system [m]

  ! 3d fields
  !**********
  real, allocatable,dimension(:,:,:,:) :: &
    uu,vv,ww,                             & ! wind components in x,y and z direction [m/s]
    uupol,vvpol,                          & ! wind components in polar stereographic projection [m/s]
    tt,tth,                               & ! temperature data on internal and half model levels [K]
    qv,qvh,                               & ! specific humidity data on internal and half model levels
    pv,                                   & ! potential vorticity
    rho,                                  & ! air density [kg/m3]
    drhodz,                               & ! vertical air density gradient [kg/m2]
    pplev,                                & ! Pressure on half model levels
    prs,                                  & ! air pressure RLT
    rho_dry                                 ! dry air density RLT Only printed out in binary mode???

  ! Cloud properties
  ! clouds:   no cloud, no precipitation   0
  !      cloud, no precipitation      1
  !      rainout  conv/lsp dominated  2/3
  !      washout  conv/lsp dominated  4/5
  ! PS 2013
  !*****************************************
  real, allocatable,dimension(:,:,:,:) :: &
    clwc,                                 & ! liquid   [kg/kg] ZHG
    ciwc,                                 & ! ice      [kg/kg] ZHG
    clw,                                  & ! combined [m3/m3] ZHG
    clwch,                                & ! original eta level liquid [kg/kg] ZHG
    ciwch                                   ! original eta level ice [kg/kg] ZHG
  real, allocatable,dimension(:,:,:) ::   &
    ctwc                                    ! ESO: =icloud_stats(:,:,4,:) total cloud water content
  integer(kind=1),allocatable,dimension(:,:,:,:) :: &
    clouds                                  ! scavenging NIK, PS
  integer,allocatable,dimension(:,:,:) :: &
    cloudsh                                 ! scavenging NIK, PS

  ! 3d nested fields
  !*****************
  real,allocatable,dimension(:,:,:,:,:) :: &
    uun, vvn, wwn,                         & ! wind components in x,y and z direction [m/s]
    ttn, tthn,                             & ! temperature data on internal and half model levels [K]
    qvn, qvhn,                             & ! specific humidity data on internal and half model levels
    pvn,                                   & ! potential vorticity
    rhon,                                  & ! air density [kg/m3]
    prsn,                                   & ! air pressure RLT
    drhodzn                                  ! vertical air density gradient [kg/m2] 

  ! ETA equivalents
  real,allocatable,dimension(:,:,:,:,:) :: &
    uuetan,vvetan,                         & ! wind components on half model levels in x and y direction [m/s]
    wwetan,                                & ! wind component on model levels in z direction [eta/s]
    ttetan,                                & ! temperature data on half model levels [K]
    pvetan,                                & ! potential vorticity on half model levels
    rhoetan,                               & ! air density on half model levels [kg/m3]
    prsetan,                               & ! air pressure on half model levels
    drhodzetan,                            & ! vertical air density gradient on half model levels [kg/m2]
    !tvirtualn,                             & ! Virtual temperature on half model levels
    etauvheightn,etawheightn                 ! Saved half model and model heights for ETA coordinate system [m]


  ! Nested cloud properties
  real,allocatable,dimension(:,:,:,:,:) :: &
    clwcn,                                 & ! liquid   [kg/kg] ZHG
    ciwcn,                                 & ! ice      [kg/kg] ZHG
    clwn,                                  & ! combined [m3/m3] ZHG
    clwchn,                                & ! original eta level liquid [kg/kg] ZHG
    ciwchn                                   ! original eta level ice [kg/kg] ZHG
  real,allocatable,dimension(:,:,:,:) ::   &
    ctwcn                                    ! ESO: =icloud_stats(:,:,4,:) total cloud water content
  integer(kind=1),allocatable,dimension(:,:,:,:,:) :: &
    cloudsn                                  ! scavenging NIK, PS
  integer,allocatable,dimension(:,:,:,:) :: &
    cloudshn                                 ! scavenging NIK, PS

  ! 2d fields
  !**********
  real, allocatable,dimension(:,:,:,:) :: &
    ps,                                 & ! surface pressure
    sd,                                 & ! snow depth
    msl,                                & ! mean sea level pressure
    tcc,                                & ! total cloud cover
    u10,                                & ! 10 meter u
    v10,                                & ! 10 meter v
    tt2,                                & ! 2 meter temperature
    td2,                                & ! 2 meter dew point
    lsprec,                             & ! large scale total precipitation [mm/h]
    convprec,                           & ! convective precipitation [mm/h]
    sshf,                               & ! surface sensible heat flux
    ssr,                                & ! surface solar radiation
    sfcstress,                            & ! surface stress
    ustar,                              & ! friction velocity [m/s]
    wstar,                              & ! convective velocity scale [m/s]
    hmix,                               & ! mixing height [m]
    tropopause,                         & ! altitude of thermal tropopause [m]
    oli                                   ! inverse Obukhov length (1/L) [m]

  ! 2d nested fields
  !*******************
  real, allocatable,dimension(:,:,:,:,:) :: &
    psn,                                    & ! surface pressure
    sdn,                                    & ! snow depth
    msln,                                   & ! mean sea level pressure
    tccn,                                   & ! total cloud cover
    u10n,                                   & ! 10 meter u
    v10n,                                   & ! 10 meter v
    tt2n,                                   & ! 2 meter temperature
    td2n,                                   & ! 2 meter dew point
    lsprecn,                                & ! large scale total precipitation [mm/h]
    convprecn,                              & ! convective precipitation [mm/h]
    sshfn,                                  & ! surface sensible heat flux
    ssrn,                                   & ! surface solar radiation
    sfcstressn,                               & ! surface stress
    ustarn,                                 & ! friction velocity [m/s]
    wstarn,                                 & ! convective velocity scale [m/s]
    hmixn,                                  & ! mixing height [m]
    tropopausen,                            & ! altitude of thermal tropopause [m]
    olin,                                   & ! inverse Obukhov length (1/L) [m]
    vdepn                                     !

  integer :: metdata_format  ! storing the input data type (ECMWF/NCEP)

  !****************************************************************************
  ! Variables defining actual size and geographical location of the wind fields
  !****************************************************************************

  integer ::   &
    nx,ny,nz,  & ! actual dimensions of wind fields in x,y and z direction, respectively
    nxmin1,    & ! nx-1
    nymin1,    & ! ny-1
    nxfield,   & ! same as nx for limited area fields, but for global fields nx=nxfield+1
    nuvz,nwz,  & ! vertical dimension of original ECMWF data (u,v components/ w components(staggered grid))
    nmixz,     & ! number of levels up to maximum PBL height (3500 m)
    nlev_ec      ! number of levels ECMWF model
  real ::      &
    dxconst,   & ! auxiliary variables for utransform
    dyconst      ! auxiliary variables for vtransform
  integer :: nconvlevmax !maximum number of levels for convection
  integer :: na ! parameter used in Emanuel's convect subroutine

  !*************************************************
  ! Variables used for vertical model discretization
  !*************************************************
  real,allocatable,dimension(:) :: &
    height,                        & ! heights of all levels [m]
    wheight,                       & ! model level heights [m]
    uvheight,                      & ! half-model level heights [m]
    akm,bkm,                       & ! coefficients which regulate vertical discretization of ecmwf model levels
    akz,bkz,                       & ! model discretization coefficients at the centre of the layers
    aknew,bknew                      ! model discretization coefficients at the interpolated levels

  !*********************************************************************
  ! Variables characterizing size and location of the nested wind fields
  !*********************************************************************

  integer,allocatable,dimension(:) :: &
    nxn,nyn                             ! actual dimensions of nested wind fields in x and y direction
  real,allocatable,dimension(:) ::    &
    dxn,dyn,                          & ! grid distances in x,y direction for the nested grids
    xlon0n,                           & ! geographical longitude of lower left grid point of nested wind fields
    ylat0n                              ! geographical latitude of lower left grid point of nested wind fields

  !*************************************************
  ! Certain auxiliary variables needed for the nests
  !*************************************************

  real,allocatable,dimension(:) :: &
    xresoln,yresoln,               & ! Factors by which the resolutions in the nests
                                     ! are enhanced compared to mother grid
    xln,yln,xrn,yrn                  ! Corner points of nested grids in grid coordinates of mother grid

contains

subroutine detectformat

  !*****************************************************************************
  !                                                                            *
  !   This routine reads the 1st file with windfields to determine             *
  !   the format.                                                              *
  !                                                                            *
  !     Authors: M. Harustak                                                   *
  !                                                                            *
  !     6 May 2015                                                             *
  !                                                                            *
  !   Unified ECMWF and GFS builds                                             *
  !   Marian Harustak, 12.5.2017                                               *
  !     - Added routine to FP10 Flexpart distribution                          *
  !*****************************************************************************
  !                                                                            *
  ! Variables:                                                                 *
  ! fname                file name of file to check                            *
  !                                                                            *
  !*****************************************************************************

  use par_mod
  use com_mod
  use class_gribfile_mod


  implicit none

  character(len=255) :: filename

  ! If no file is available
  if ( numbwf.le.0 ) then
    print*,'No wind file available'
    metdata_format = GRIBFILE_CENTRE_UNKNOWN
    return
  endif

  ! construct filename
  filename = path(3)(1:length(3)) // trim(wfname(1))
 
  ! get format
  metdata_format = gribfile_centre(TRIM(filename))
end subroutine detectformat

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
  integer :: gotGrid,stat
  real(kind=4) :: xaux1,xaux2,yaux1,yaux2
  real(kind=8) :: xaux1in,xaux2in,yaux1in,yaux2in
  integer :: gribVer,parCat,parNum,typSurf,valSurf,discipl,parId
  !HSO  end
  integer :: ix,jy,i,ifn,ifield,j,k,iumax,iwmax,numskip
  real :: sizesouth,sizenorth,xauxa

  ! VARIABLES AND ARRAYS NEEDED FOR GRIB DECODING

  ! dimension of isec2 at least (22+n), where n is the number of parallels or
  ! meridians in a quasi-regular (reduced) Gaussian or lat/long grid

  ! dimension of zsec2 at least (10+nn), where nn is the number of vertical
  ! coordinate parameters

  integer :: isec1(56), isec2(12)
  real(kind=4),allocatable,dimension(:) :: zsec2,zsec4
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
  do while(.true.)
    ifield=ifield+1
    !
    ! GET NEXT FIELDS
    !
    call grib_new_from_file(ifile,igrib,iret)
    if (iret.eq.GRIB_END_OF_FILE ) then
      exit   ! EOF DETECTED
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

    call grib_get_int(igrib,'numberOfPointsAlongAParallel', &
         isec2(2),iret)
    ! ! nx=isec2(2)
    ! ! WRITE(*,*) nx,nxmax
    ! if (isec2(2).gt.nxmax) then
    !   WRITE(*,*) 'FLEXPART error: Too many grid points in x direction.'
    !   WRITE(*,*) 'Reduce resolution of wind fields.'
    !   WRITE(*,*) 'Or change parameter settings in file ecmwf_mod.'
    !   WRITE(*,*) isec2(2),nxmax
    ! !    STOP
    ! endif

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

      nxfield=isec2(2)
      ny=isec2(3)
      nlev_ec=isec2(12)/2-1

      allocate(zsec2(60+2*nlev_ec), stat=stat)
      if (stat.ne.0) error stop "Could not allocate zsec2"
      allocate(zsec4(16*nxfield*ny), stat=stat)
      if (stat.ne.0) error stop "Could not allocate zsec4"
      !  get the size and data of the vertical coordinate array
      call grib_get_real4_array(igrib,'pv',zsec2,iret)
      call grib_check(iret,gribFunction,gribErrorMsg)
    endif

    !get the size and data of the values array
    if (isec1(6).ne.-1) then
      call grib_get_real4_array(igrib,'values',zsec4,iret)
      call grib_check(iret,gribFunction,gribErrorMsg)
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
      xaux1=real(xaux1in)
      xaux2=real(xaux2in)
      yaux1=real(yaux1in)
      yaux2=real(yaux2in)
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
          error stop 'nxshift in file par_mod is too large'
        xlon0=xlon0+real(nxshift)*dx
      else
        nx=nxfield
        xglobal=.false.
        if (nxshift.ne.0) &
          error stop 'nxshift (par_mod) must be zero for non-global domain'
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
        error stop 'nxshift (par_mod) must not be negative'
      if (nxshift.ge.nxfield) error stop 'nxshift (par_mod) too large'
    endif ! gotGrid

    ! if (nx.gt.nxmax) then
    !   write(*,*) 'FLEXPART error: Too many grid points in x direction.'
    !   write(*,*) 'Reduce resolution of wind fields.'
    !   write(*,*) 'Or change parameter settings in file par_mod.'
    !   write(*,*) nx,nxmax
    !   error stop
    ! endif

    ! if (ny.gt.nymax) then
    !   write(*,*) 'FLEXPART error: Too many grid points in y direction.'
    !   write(*,*) 'Reduce resolution of wind fields.'
    !   write(*,*) 'Or change parameter settings in file par_mod.'
    !   write(*,*) ny,nymax
    !   error stop
    ! endif

    k=isec1(8)
    if(isec1(6).eq.131) iumax=max(iumax,nlev_ec-k+1)
    if(isec1(6).eq.135) iwmax=max(iwmax,nlev_ec-k+1)

    if (isec1(6) .eq. 167) then
      ! Asking grid values and allocate memory to read windfields
      nxmax=nxfield
      if (xglobal) then
        nxmax=nxfield+1
      endif
      nymax=ny
      nwzmax=iwmax+1
      nuvzmax=iumax+1
      nzmax=nuvzmax
      nconvlevmax=iumax
      na=nuvzmax
      ! Temporary nxmax and nymax
      call alloc_fixedfields
      write(*,*) 'grid dim:',nxmax,nymax,nwzmax,nuvzmax,nconvlevmax,na
    endif

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
  end do                      !! READ NEXT LEVEL OR PARAMETER
  !
  ! CLOSING OF INPUT DATA FILE
  !
  call grib_close_file(ifile)

  ! Allocate memory for windfields
  !*******************************
  call alloc_windfields

  !error message if no fields found with correct first longitude in it
  if (gotGrid.eq.0) then
    print*,'***ERROR: input file needs to contain GRiB1 formatted'// &
         'messages'
    error stop
  endif

  nuvz=iumax
  nwz =iwmax
  if(nuvz.eq.nlev_ec) nwz=nlev_ec+1

  ! if (nuvz+1.gt.nuvzmax) then
  !   write(*,*) 'FLEXPART error: Too many u,v grid points in z '// &
  !        'direction.'
  !   write(*,*) 'Reduce resolution of wind fields.'
  !   write(*,*) 'Or change parameter settings in file par_mod.'
  !   write(*,*) nuvz+1,nuvzmax
  !   stop
  ! endif

  ! if (nwz.gt.nwzmax) then
  !   write(*,*) 'FLEXPART error: Too many w grid points in z '// &
  !        'direction.'
  !   write(*,*) 'Reduce resolution of wind fields.'
  !   write(*,*) 'Or change parameter settings in file par_mod.'
  !   write(*,*) nwz,nwzmax
  !   stop
  ! endif

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
  akm=0
  bkm=0
  akz=0
  bkz=0
  do i=1,nwz ! LB: should start counting fom 0 to get the top level?
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
  ! verttransform.f and verttranform_nest.f
  !*****************************************************************************

  nz=nuvz
  if (nz.gt.nzmax) error stop 'nzmax too small'
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
    error stop
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
  integer :: igrib,stat
  real(kind=4) :: xaux1,xaux2,yaux1,yaux2
  real(kind=8) :: xaux1in,xaux2in,yaux1in,yaux2in
  integer :: gribVer,parCat,parNum,typSurf,valSurf,discipl
  !HSO  end
  integer :: ix,jy,i,ifn,ifield,j,k,iumax,iwmax,numskip
  real :: sizesouth,sizenorth,xauxa
  real,allocatable,dimension(:) :: akm_usort,pres,tmppres
  real,parameter :: eps=0.0001

  ! NCEP GFS
  real :: help

  integer :: i179,i180,i181

  ! VARIABLES AND ARRAYS NEEDED FOR GRIB DECODING

  integer :: isec1(8),isec2(3)
  real(kind=4),allocatable,dimension(:) :: zsec4
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
  error stop
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
  do
    ifield=ifield+1
    !
    ! GET NEXT FIELDS
    !
    call grib_new_from_file(ifile,igrib,iret)
    if (iret.eq.GRIB_END_OF_FILE )  then
      exit    ! EOF DETECTED
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

      xaux1=real(xaux1in)
      xaux2=real(xaux2in)
      yaux1=real(yaux1in)
      yaux2=real(yaux2in)

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
          error stop 'nxshift in file par_mod is too large'
        xlon0=xlon0+real(nxshift)*dx
      else
        nx=nxfield
        xglobal=.false.
        if (nxshift.ne.0) &
          error stop 'nxshift (par_mod) must be zero for non-global domain'
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
      ! Set nxmax and nymax and allocate the fields for oro lsm and excessoro
      nxmax=nx
      nymax=ny
      call alloc_fixedfields
      if (gribVer.eq.1) then
        allocate( zsec4(16*nxmax*nymax),stat=stat) ! GRIB Edition 1
        if (stat.ne.0) error stop "Could not allocate zsec2"
      endif
    endif ! ifield.eq.1

    if (nxshift.lt.0) error stop 'nxshift (par_mod) must not be negative'
    if (nxshift.ge.nxfield) error stop 'nxshift (par_mod) too large'

    if (gribVer.eq.1) then ! GRIB Edition 1
      if (isec1(6).ne.-1) then
      !  get the size and data of the values array
        call grib_get_real4_array(igrib,'values',zsec4,iret)
        call grib_check(iret,gribFunction,gribErrorMsg)
      endif
    endif

    ! NCEP ISOBARIC LEVELS
    !*********************

    if((isec1(6).eq.33).and.(isec1(7).eq.100)) then ! check for U wind
      iumax=iumax+1
      allocate( tmppres(iumax), stat=stat)
      if (stat.ne.0) error stop "Could not allocate tmppres"
      if (iumax.gt.1) tmppres(1:iumax)=pres
      pres(iumax)=real(isec1(8))*100.0
      call move_alloc(tmppres,pres)
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
    ! IP 8/5/23 allocate fields missing for GFS reading 
    call alloc_fixedfields
    ! IP 8/5/23
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

  end do                      !! READ NEXT LEVEL OR PARAMETER
  !
  ! CLOSING OF INPUT DATA FILE
  !

  ! HSO
  call grib_close_file(ifile)
  ! HSO end edits

  nuvz=iumax
  nwz =iumax
  nlev_ec=iumax

  ! Allocate memory for windfields
  !*******************************
  nwzmax=nwz
  nuvzmax=nuvz
  nzmax=nuvz
  call alloc_windfields

  if (nx.gt.nxmax) then
    write(*,*) 'FLEXPART error: Too many grid points in x direction.'
    write(*,*) 'Reduce resolution of wind fields.'
    write(*,*) 'Or change parameter settings in file par_mod.'
    write(*,*) nx,nxmax
    error stop
  endif

  if (ny.gt.nymax) then
    write(*,*) 'FLEXPART error: Too many grid points in y direction.'
    write(*,*) 'Reduce resolution of wind fields.'
    write(*,*) 'Or change parameter settings in file par_mod.'
    write(*,*) ny,nymax
    error stop
  endif

  if (nuvz.gt.nuvzmax) then
    write(*,*) 'FLEXPART error: Too many u,v grid points in z '// &
         'direction.'
    write(*,*) 'Reduce resolution of wind fields.'
    write(*,*) 'Or change parameter settings in file par_mod.'
    write(*,*) nuvz,nuvzmax
    error stop
  endif

  if (nwz.gt.nwzmax) then
    write(*,*) 'FLEXPART error: Too many w grid points in z '// &
         'direction.'
    write(*,*) 'Reduce resolution of wind fields.'
    write(*,*) 'Or change parameter settings in file par_mod.'
    write(*,*) nwz,nwzmax
    error stop
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

  allocate( akm_usort(nwzmax), stat=stat)
  if (stat.ne.0) error stop "Could not allocate akm_usort"
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
  ! verttransform.f and verttranform_nest.f
  !*****************************************************************************

  nz=nuvz
  if (nz.gt.nzmax) error stop 'nzmax too small'
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
  deallocate(  akm_usort,pres,zsec4 )
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
    error stop
  else
    goto 5
  endif

end subroutine gridcheck_gfs

subroutine gridcheck_nest

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
  integer :: igrib,stat
  integer :: gribVer,parCat,parNum,typSurf,valSurf,discipl
  integer :: parID !added by mc for making it consistent with new gridcheck.f90
  integer :: gotGrib
  !HSO  end
  integer :: i,j,k,l,ifn,ifield,iumax,iwmax,numskip,nlev_ecn
  integer :: nuvzn,nwzn
  real :: akmn(nwzmax),bkmn(nwzmax),akzn(nuvzmax),bkzn(nuvzmax)
  real(kind=4) :: xaux1,xaux2,yaux1,yaux2
  real(kind=8) :: xaux1in,xaux2in,yaux1in,yaux2in

  ! VARIABLES AND ARRAYS NEEDED FOR GRIB DECODING

  ! dimension of isec2 at least (22+n), where n is the number of parallels or
  ! meridians in a quasi-regular (reduced) Gaussian or lat/long grid

  ! dimension of zsec2 at least (10+nn), where nn is the number of vertical
  ! coordinate parameters

  integer :: isec1(56),isec2(12) !(22+nxmaxn+nymaxn)
  real(kind=4),allocatable,dimension(:) :: zsec2,zsec4

  !HSO  grib api error messages
  character(len=24) :: gribErrorMsg = 'Error reading grib file'
  character(len=20) :: gribFunction = 'gridcheck_nest'

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

    call grib_open_file(ifile,path(numpath+2*(l-1)+1) &
           (1:length(numpath+2*(l-1)+1))//trim(wfnamen(l,ifn)),'r',iret)
    if (iret.ne.GRIB_SUCCESS) then
      goto 999   ! ERROR DETECTED
    endif
    !turn on support for multi fields messages
    !call grib_multi_support_on

    gotGrib=0
    ifield=0
    do
      ifield=ifield+1

      !
      ! GET NEXT FIELDS
      !
      call grib_new_from_file(ifile,igrib,iret)
      if (iret.eq.GRIB_END_OF_FILE)  then
        exit    ! EOF DETECTED
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

        nxn(l)=isec2(2)
        nyn(l)=isec2(3)
        nlev_ecn=isec2(12)/2-1

        if (nxn(l).gt.nxmaxn) nxmaxn=nxn(l)
        if (nyn(l).gt.nymaxn) nymaxn=nyn(l)
        allocate( zsec2(60+2*nlev_ecn), stat=stat)
        if (stat.ne.0) error stop "Could not allocate zsec2"
        allocate( zsec4(16*nxmaxn*nymaxn), stat=stat)
        if (stat.ne.0) error stop "Could not allocate zsec4"

        !HSO    get the size and data of the vertical coordinate array
        call grib_get_real4_array(igrib,'pv',zsec2,iret)
        call grib_check(iret,gribFunction,gribErrorMsg)

        call alloc_fixedfields_nest
        write(*,*) 'Dimensions nest:',nxmaxn,nymaxn,nlev_ecn
      endif ! ifield

      !get the size and data of the values array
      if (isec1(6).ne.-1) then
        call grib_get_real4_array(igrib,'values',zsec4,iret)
        call grib_check(iret,gribFunction,gribErrorMsg)
      endif

      if (nxn(l).gt.nxmaxn) then
        write(*,*) 'FLEXPART error: Too many grid points in x direction.'
        write(*,*) 'Reduce resolution of wind fields (file GRIDSPEC)'
        write(*,*) 'for nesting level ',l
        write(*,*) 'Or change parameter settings in file par_mod.'
        write(*,*) nxn(l),nxmaxn
        error stop
      endif

      if (nyn(l).gt.nymaxn) then
        write(*,*) 'FLEXPART error: Too many grid points in y direction.'
        write(*,*) 'Reduce resolution of wind fields (file GRIDSPEC)'
        write(*,*) 'for nesting level ',l
        write(*,*) 'Or change parameter settings in file par_mod.'
        write(*,*) nyn(l),nymaxn
        error stop
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
        xaux1=real(xaux1in)
        xaux2=real(xaux2in)
        yaux1=real(yaux1in)
        yaux2=real(yaux2in)
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
    end do                 !! READ NEXT LEVEL OR PARAMETER
    !
    ! CLOSING OF INPUT DATA FILE
    !

    call grib_close_file(ifile)

    !error message if no fields found with correct first longitude in it
    if (gotGrib.eq.0) then
      print*,'***ERROR: input file needs to contain GRiB1 formatted'// &
           'messages'
      error stop
    endif

    nuvzn=iumax
    nwzn=iwmax
    if(nuvzn.eq.nlev_ec) nwzn=nlev_ecn+1

    if ((nuvzn.gt.nuvzmax).or.(nwzn.gt.nwzmax)) then
      write(*,*) 'FLEXPART error: Nested wind fields have too many '// &
           'vertical levels.'
      write(*,*) 'Problem was encountered for nesting level ',l
      error stop
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
      error stop
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
        error stop
      endif
    end do

    do i=1,nwz
      if ((akmn(i).ne.akm(i)).or.(bkmn(i).ne.bkm(i))) then
        write(*,*) 'FLEXPART error: The wind fields of nesting level',l
        write(*,*) 'are not consistent with the mother domain:'
        write(*,*) 'Differences in vertical levels detected.'
        error stop
      endif
    end do

    deallocate( zsec2,zsec4 )
  end do

  ! Allocate memory for windfields using nxmax, nymaxn, numbnest
  !*************************************************************
  if (numbnests.ge.1) then
    ! If nested wind fields are used, allocate arrays
    !************************************************
    call alloc_windfields_nest
  endif


  return

999   write(*,*)
  write(*,*) ' ###########################################'// &
       '###### '
  write(*,*) '                FLEXPART SUBROUTINE GRIDCHECK:'
  write(*,*) ' CAN NOT OPEN INPUT DATA FILE '//wfnamen(l,ifn)
  write(*,*) ' FOR NESTING LEVEL ',k
  write(*,*) ' ###########################################'// &
       '###### '
  error stop

end subroutine gridcheck_nest

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
  !  L. Bakels, 2021: OpenMP parallelisation (following CTM version)
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

  ! VARIABLES AND ARRAYS NEEDED FOR GRIB DECODING

  ! dimension of isec2 at least (22+n), where n is the number of parallels or
  ! meridians in a quasi-regular (reduced) Gaussian or lat/long grid

  ! dimension of zsec2 at least (10+nn), where nn is the number of vertical
  ! coordinate parameters

  integer :: isec1(56)
  integer,allocatable,dimension(:) :: isec2
  real(kind=4),allocatable,dimension(:) :: zsec4
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
  if (stat.ne.0) error stop "Could not allocate igrib"
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
!$OMP   excessoro, lsm, nymin1,ciwch,clwch,readclouds,sumclouds, nxshift,nxmax,nymax) & 
!$OMP PRIVATE(ii, gribVer, iret, isec1, discipl, parCat, parNum, parId,typSurf, valSurf, &
!$OMP   zsec4, isec2, gribErrorMsg, xauxin, yauxin, xaux, yaux, xaux0,  &
!$OMP   yaux0, k, arsize, stat, conversion_factor)  &
!$OMP REDUCTION(+:gotGrid)
  !
  ! GET NEXT FIELDS
  !
  allocate( isec2(22+nxmax+nymax),zsec4( 16*nxmax*nymax ),stat=stat)
  if (stat.ne.0) write(*,*)'ERROR: could not allocate isec2 or zsec4'

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
    if(isec2(2).ne.nxfield) error stop 'READWIND: NX NOT CONSISTENT'
    if(isec2(3).ne.ny) error stop 'READWIND: NY NOT CONSISTENT'
    if(isec2(12)/2-1.ne.nlev_ec) &
      error stop 'READWIND: VERTICAL DISCRETIZATION NOT CONSISTENT'
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

    xaux=real(xauxin)+real(nxshift)*dx
    yaux=real(yauxin)
    if (xaux.gt.180.) xaux=xaux-360.0
    if(abs(xaux-xlon0).gt.eps) &
      error stop 'READWIND: LOWER LEFT LONGITUDE NOT CONSISTENT'
    if(abs(yaux-ylat0).gt.eps) &
      error stop 'READWIND: LOWER LEFT LATITUDE NOT CONSISTENT'
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

  deallocate( zsec4,isec2 )
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
    error stop
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
      if ((ewss(i,j).eq.0.).and.(nsss(i,j).eq.0.)) then
        if ((i.ne.0).and.(j.ne.0).and.(i.ne.nxmin1).and.(j.ne.nymin1)) then
          ewss(i,j)=(ewss(i-1,j-1)+ewss(i+1,j+1)+ewss(i+1,j)+ewss(i-1,j)+ &
                     ewss(i,j+1)+ewss(i,j-1)+ewss(i-1,j+1)+ewss(i+1,j-1))/8.
          nsss(i,j)=(nsss(i-1,j-1)+nsss(i+1,j+1)+nsss(i+1,j)+nsss(i-1,j)+ &
                     nsss(i,j+1)+nsss(i,j-1)+nsss(i-1,j+1)+nsss(i+1,j-1))/8.
        else if ((i.eq.0).and.(j.eq.0)) then
          ewss(i,j)=(ewss(i+1,j+1)+ewss(i+1,j)+ewss(i,j+1))/3.
          nsss(i,j)=(nsss(i+1,j+1)+nsss(i+1,j)+nsss(i,j+1))/3.
        else if ((i.eq.nxmin1).and.(j.eq.nymin1)) then
          ewss(i,j)=(ewss(i-1,j-1)+ewss(i-1,j)+ewss(i,j-1))/3.
          nsss(i,j)=(nsss(i-1,j-1)+nsss(i-1,j)+nsss(i,j-1))/3.
        else if ((i.eq.0).and.(j.eq.nymin1)) then
          ewss(i,j)=(ewss(i+1,j-1)+ewss(i+1,j)+ewss(i,j-1))/3.
          nsss(i,j)=(nsss(i+1,j-1)+nsss(i+1,j)+nsss(i,j-1))/3.
        else if ((i.eq.nxmin1).and.(j.eq.0)) then
          ewss(i,j)=(ewss(i-1,j+1)+ewss(i-1,j)+ewss(i,j+1))/3.
          nsss(i,j)=(nsss(i-1,j+1)+nsss(i-1,j)+nsss(i,j+1))/3.
        else if (i.eq.0) then
          ewss(i,j)=(ewss(i+1,j+1)+ewss(i+1,j)+ewss(i,j+1)+ewss(i,j-1)+ewss(i+1,j-1))/5.
          nsss(i,j)=(nsss(i+1,j+1)+nsss(i+1,j)+nsss(i,j+1)+nsss(i,j-1)+nsss(i+1,j-1))/5.
        else if (i.eq.nxmin1) then
          ewss(i,j)=(ewss(i-1,j+1)+ewss(i-1,j)+ewss(i,j+1)+ewss(i,j-1)+ewss(i-1,j-1))/5.
          nsss(i,j)=(nsss(i-1,j+1)+nsss(i-1,j)+nsss(i,j+1)+nsss(i,j-1)+nsss(i-1,j-1))/5.
        else if (j.eq.0) then
          ewss(i,j)=(ewss(i+1,j+1)+ewss(i+1,j)+ewss(i-1,j)+ewss(i,j+1)+ewss(i-1,j+1))/5.
          nsss(i,j)=(nsss(i+1,j+1)+nsss(i+1,j)+nsss(i-1,j)+nsss(i,j+1)+nsss(i-1,j+1))/5.
        else if (j.eq.nymin1) then
          ewss(i,j)=(ewss(i+1,j-1)+ewss(i+1,j)+ewss(i-1,j)+ewss(i,j-1)+ewss(i-1,j-1))/5.
          nsss(i,j)=(nsss(i+1,j-1)+nsss(i+1,j)+nsss(i-1,j)+nsss(i,j-1)+nsss(i-1,j-1))/5.
        endif
      endif
      sfcstress(i,j,1,n)=sqrt(ewss(i,j)**2+nsss(i,j)**2)
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
             sfcstress(i,j,1,n),sshf(i,j,1,n))
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

  if(iumax.ne.nuvz-1) error stop 'READWIND: NUVZ NOT CONSISTENT'
  if(iwmax.ne.nwz)    error stop 'READWIND: NWZ NOT CONSISTENT'

  return

888 write(*,*) ' #### FLEXPART MODEL ERROR! WINDFIELD         #### '
  write(*,*) ' #### ',wfname(indj),'                    #### '
  write(*,*) ' #### IS NOT GRIB FORMAT !!!                  #### '
  error stop 'Execution terminated'

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
  real(kind=4) :: zsec4(16*nxmax*nymax)
  real(kind=4) :: xaux,yaux,xaux0,yaux0
  real(kind=8) :: xauxin,yauxin
  real,parameter :: eps=1.e-4
  real(kind=4) :: ewss(0:nxmax-1,0:nymax-1),nsss(0:nxmax-1,0:nymax-1)
  real :: plev1,hlev1,ff10m,fflev1

  logical :: hflswitch,strswitch

  !HSO  for grib api error messages
  character(len=24) :: gribErrorMsg = 'Error reading grib file'
  ! character(len=20) :: gribFunction = 'readwind_gfs'
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
  do
    ifield=ifield+1
    !
    ! GET NEXT FIELDS
    !
    call grib_new_from_file(ifile,igrib,iret)
    if (iret.eq.GRIB_END_OF_FILE)  then
      exit   ! EOF DETECTED
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
    xaux=real(xauxin)+real(nxshift)*dx
    yaux=real(yauxin)

    ! CHECK GRID SPECIFICATIONS

      if(isec2(2).ne.nxfield) error stop 'READWIND: NX NOT CONSISTENT'
      if(isec2(3).ne.ny) error stop 'READWIND: NY NOT CONSISTENT'
      if(xaux.eq.0.) xaux=-179.0     ! NCEP DATA
      xaux0=xlon0
      yaux0=ylat0
      if(xaux.lt.0.) xaux=xaux+360.
      if(yaux.lt.0.) yaux=yaux+360.
      if(xaux0.lt.0.) xaux0=xaux0+360.
      if(yaux0.lt.0.) yaux0=yaux0+360.
      if(abs(xaux-xaux0).gt.eps) &
           error stop 'READWIND: LOWER LEFT LONGITUDE NOT CONSISTENT'
      if(abs(yaux-yaux0).gt.eps) &
           error stop 'READWIND: LOWER LEFT LATITUDE NOT CONSISTENT'
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
  end do                      !! READ NEXT LEVEL OR PARAMETER
  !
  ! CLOSING OF INPUT DATA FILE
  !

  !HSO close grib file
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
      sfcstress(i,j,1,n)=sqrt(ewss(i,j)**2+nsss(i,j)**2)
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
             sfcstress(i,j,1,n),sshf(i,j,1,n))
        if(sshf(i,j,1,n).gt.200.) sshf(i,j,1,n)=200.
        if(sshf(i,j,1,n).lt.-400.) sshf(i,j,1,n)=-400.
      end do
    end do
  endif

  if(iumax.ne.nuvz) error stop 'READWIND: NUVZ NOT CONSISTENT'
  if(iumax.ne.nwz) error stop 'READWIND: NWZ NOT CONSISTENT'

  return
888   write(*,*) ' #### FLEXPART MODEL ERROR! WINDFIELD         #### '
  write(*,*) ' #### ',wfname(indj),'                    #### '
  write(*,*) ' #### IS NOT GRIB FORMAT !!!                  #### '
  error stop 'Execution terminated'

end subroutine readwind_gfs

subroutine readwind_nest(indj,n,uuhn,vvhn,wwhn)
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

  real :: uuhn(0:nxmaxn-1,0:nymaxn-1,nuvzmax,numbnests)
  real :: vvhn(0:nxmaxn-1,0:nymaxn-1,nuvzmax,numbnests)
  real :: wwhn(0:nxmaxn-1,0:nymaxn-1,nwzmax,numbnests)
  integer :: indj,i,j,k,n,levdiff2,ifield,iumax,iwmax,l

  ! VARIABLES AND ARRAYS NEEDED FOR GRIB DECODING

  ! dimension of isec2 at least (22+n), where n is the number of parallels or
  ! meridians in a quasi-regular (reduced) Gaussian or lat/long grid

  ! dimension of zsec2 at least (10+nn), where nn is the number of vertical
  ! coordinate parameters

  integer :: isec1(56),isec2(22+nxmaxn+nymaxn)
  real(kind=4) :: zsec4(16*nxmaxn*nymaxn)
  real(kind=4) :: xaux,yaux
  real(kind=8) :: xauxin,yauxin
  real,parameter :: eps=1.e-4
  real :: ewss(0:nxmaxn-1,0:nymaxn-1),nsss(0:nxmaxn-1,0:nymaxn-1)
  real :: plev1,pmean,tv,fu,hlev1,ff10m,fflev1
  real :: conversion_factor !added by mc to make it consistent with new gridchek.f90

  logical :: hflswitch,strswitch

  !HSO  grib api error messages
  character(len=24) :: gribErrorMsg = 'Error reading grib file'
  character(len=20) :: gribFunction = 'readwind_nest'

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

    call grib_open_file(ifile,path(numpath+2*(l-1)+1) &
         (1:length(numpath+2*(l-1)+1))//trim(wfnamen(l,indj)),'r')
    if (iret.ne.GRIB_SUCCESS) then
      goto 888   ! ERROR DETECTED
    endif
    !turn on support for multi fields messages */
    !call grib_multi_support_on

    gotGrid=0
    ifield=0
    do
      ifield=ifield+1
      !
      ! GET NEXT FIELDS
      !
      call grib_new_from_file(ifile,igrib,iret)
      if (iret.eq.GRIB_END_OF_FILE)  then
        exit    ! EOF DETECTED
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
      if(isec2(2).ne.nxn(l)) error stop &
      'READWIND: NX NOT CONSISTENT FOR A NESTING LEVEL'
      if(isec2(3).ne.nyn(l)) error stop &
      'READWIND: NY NOT CONSISTENT FOR A NESTING LEVEL'
      if(isec2(12)/2-1.ne.nlev_ec) error stop 'READWIND: VERTICAL DISCRET&
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

        xaux=real(xauxin)
        yaux=real(yauxin)
        if (abs(xaux-xlon0n(l)).gt.eps) &
        error stop 'READWIND: LOWER LEFT LONGITUDE NOT CONSISTENT FOR A NESTING LEVEL'
        if (abs(yaux-ylat0n(l)).gt.eps) &
        error stop 'READWIND: LOWER LEFT LATITUDE NOT CONSISTENT FOR A NESTING LEVEL'
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
    end do                      !! READ NEXT LEVEL OR PARAMETER
  !
  ! CLOSING OF INPUT DATA FILE
  !
    call grib_close_file(ifile)

    !error message if no fields found with correct first longitude in it
    if (gotGrid.eq.0) then
      print*,'***ERROR: input file needs to contain GRiB1 formatted'// &
           'messages'
      error stop
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
        sfcstressn(i,j,1,n,l)=sqrt(ewss(i,j)**2+nsss(i,j)**2)
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
               sfcstressn(i,j,1,n,l),sshfn(i,j,1,n,l))
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

    if(iumax.ne.nuvz-1) error stop &
         'READWIND: NUVZ NOT CONSISTENT FOR A NESTING LEVEL'
    if(iwmax.ne.nwz) error stop &
         'READWIND: NWZ NOT CONSISTENT FOR A NESTING LEVEL'

  end do

  return
888   write(*,*) ' #### FLEXPART MODEL ERROR! WINDFIELD         #### '
  write(*,*) ' #### ',wfnamen(l,indj),' FOR NESTING LEVEL  #### '
  write(*,*) ' #### ',l,' IS NOT GRIB FORMAT !!!           #### '
  error stop 'Execution terminated'


999   write(*,*) ' #### FLEXPART MODEL ERROR! WINDFIELD         #### '
  write(*,*) ' #### ',wfnamen(l,indj),'                    #### '
  write(*,*) ' #### CANNOT BE OPENED FOR NESTING LEVEL ',l,'####'

end subroutine readwind_nest

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

subroutine alloc_fixedfields
  implicit none
  integer :: stat

  allocate(oro(0:nxmax-1,0:nymax-1),stat=stat)
  if (stat.ne.0) error stop "Could not allocate oro"
  allocate(excessoro(0:nxmax-1,0:nymax-1),stat=stat)
  if (stat.ne.0) error stop "Could not allocate excessoro"
  allocate(lsm(0:nxmax-1,0:nymax-1),stat=stat)
  if (stat.ne.0) error stop "Could not allocate lsm"
  allocate(pv(0:nxmax-1,0:nymax-1,nzmax,numwfmem),stat=stat)
  if (stat.ne.0) error stop "Could not allocate pv"
end subroutine alloc_fixedfields

subroutine alloc_fixedfields_nest
  implicit none 
  integer :: stat

  allocate(oron(0:nxmaxn-1,0:nymaxn-1,numbnests),stat=stat)
  if (stat.ne.0) error stop "Could not allocate oron"
  allocate(excessoron(0:nxmaxn-1,0:nymaxn-1,numbnests),stat=stat)
  if (stat.ne.0) error stop "Could not allocate excessoron"
  allocate(lsmn(0:nxmaxn-1,0:nymaxn-1,numbnests),stat=stat)
  if (stat.ne.0) error stop "Could not allocate lsmn"
end subroutine alloc_fixedfields_nest

subroutine alloc_windfields
  implicit none
  integer :: stat
  ! Eta coordinates
  !****************
#ifdef ETA
  allocate(uueta(0:nxmax-1,0:nymax-1,nzmax,numwfmem),stat=stat)
  if (stat.ne.0) error stop "Could not allocate uueta"
  allocate(vveta(0:nxmax-1,0:nymax-1,nzmax,numwfmem),stat=stat)
  if (stat.ne.0) error stop "Could not allocate vveta"
  allocate(wweta(0:nxmax-1,0:nymax-1,nzmax,numwfmem),stat=stat)
  if (stat.ne.0) error stop "Could not allocate wweta"
  allocate(uupoleta(0:nxmax-1,0:nymax-1,nzmax,numwfmem),stat=stat)
  if (stat.ne.0) error stop "Could not allocate uupoleta"
  allocate(vvpoleta(0:nxmax-1,0:nymax-1,nzmax,numwfmem),stat=stat)
  if (stat.ne.0) error stop "Could not allocate vvpoleta"
  allocate(tteta(0:nxmax-1,0:nymax-1,nzmax,numwfmem),stat=stat)
  if (stat.ne.0) error stop "Could not allocate tteta"
  allocate(pveta(0:nxmax-1,0:nymax-1,nzmax,numwfmem),stat=stat)
  if (stat.ne.0) error stop "Could not allocate pveta"
  allocate(prseta(0:nxmax-1,0:nymax-1,nzmax,numwfmem),stat=stat)
  if (stat.ne.0) error stop "Could not allocate prseta"
  allocate(rhoeta(0:nxmax-1,0:nymax-1,nzmax,numwfmem),stat=stat)
  if (stat.ne.0) error stop "Could not allocate rhoeta"
  allocate(drhodzeta(0:nxmax-1,0:nymax-1,nzmax,numwfmem),stat=stat)
  if (stat.ne.0) error stop "Could not allocate drhodzeta"
#endif
  allocate(etauvheight(0:nxmax-1,0:nymax-1,nuvzmax,numwfmem),stat=stat)
  if (stat.ne.0) error stop "Could not allocate etauvheight"
  allocate(etawheight(0:nxmax-1,0:nymax-1,nuvzmax,numwfmem),stat=stat)
  if (stat.ne.0) error stop "Could not allocate etawheight"

  ! Intrinsic coordinates
  !**********************
  allocate(uu(0:nxmax-1,0:nymax-1,nzmax,numwfmem),stat=stat)
  if (stat.ne.0) error stop "Could not allocate uu"
  allocate(vv(0:nxmax-1,0:nymax-1,nzmax,numwfmem),stat=stat)
  if (stat.ne.0) error stop "Could not allocate vv"
  allocate(ww(0:nxmax-1,0:nymax-1,nzmax,numwfmem),stat=stat)
  if (stat.ne.0) error stop "Could not allocate ww"
  allocate(uupol(0:nxmax-1,0:nymax-1,nzmax,numwfmem),stat=stat)
  if (stat.ne.0) error stop "Could not allocate uupol"
  allocate(vvpol(0:nxmax-1,0:nymax-1,nzmax,numwfmem),stat=stat)
  if (stat.ne.0) error stop "Could not allocate vvpol"
  allocate(tt(0:nxmax-1,0:nymax-1,nzmax,numwfmem),stat=stat)
  if (stat.ne.0) error stop "Could not allocate tt"
  allocate(tth(0:nxmax-1,0:nymax-1,nuvzmax,numwfmem),stat=stat)
  if (stat.ne.0) error stop "Could not allocate tth"
  allocate(qv(0:nxmax-1,0:nymax-1,nzmax,numwfmem),stat=stat)
  if (stat.ne.0) error stop "Could not allocate qv"
  allocate(qvh(0:nxmax-1,0:nymax-1,nuvzmax,numwfmem),stat=stat)
  if (stat.ne.0) error stop "Could not allocate qvh"
  allocate(rho(0:nxmax-1,0:nymax-1,nzmax,numwfmem),stat=stat)
  if (stat.ne.0) error stop "Could not allocate rho"
  allocate(drhodz(0:nxmax-1,0:nymax-1,nzmax,numwfmem),stat=stat)
  if (stat.ne.0) error stop "Could not allocate drhodz"
  allocate(pplev(0:nxmax-1,0:nymax-1,nuvzmax,numwfmem),stat=stat)
  if (stat.ne.0) error stop "Could not allocate pplev"
  allocate(prs(0:nxmax-1,0:nymax-1,nzmax,numwfmem),stat=stat)
  if (stat.ne.0) error stop "Could not allocate prs"
  allocate(rho_dry(0:nxmax-1,0:nymax-1,nzmax,numwfmem),stat=stat)
  if (stat.ne.0) error stop "Could not allocate rho_dry"

  ! Cloud data
  !***********
  allocate(clwc(0:nxmax-1,0:nymax-1,nzmax,numwfmem),stat=stat)
  if (stat.ne.0) error stop "Could not allocate clwc"
  allocate(ciwc(0:nxmax-1,0:nymax-1,nzmax,numwfmem),stat=stat)
  if (stat.ne.0) error stop "Could not allocate ciwc"
  allocate(clw(0:nxmax-1,0:nymax-1,nzmax,numwfmem),stat=stat)
  if (stat.ne.0) error stop "Could not allocate clw"
  allocate(clwch(0:nxmax-1,0:nymax-1,nuvzmax,numwfmem),stat=stat)
  if (stat.ne.0) error stop "Could not allocate clwch"
  allocate(ciwch(0:nxmax-1,0:nymax-1,nuvzmax,numwfmem),stat=stat)
  if (stat.ne.0) error stop "Could not allocate ciwch"
  clwc=0.0
  ciwc=0.0
  clw=0.0
  clwch=0.0
  ciwch=0.0
  allocate(ctwc(0:nxmax-1,0:nymax-1,numwfmem),stat=stat)
  if (stat.ne.0) error stop "Could not allocate ctwc"
  allocate(cloudsh(0:nxmax-1,0:nymax-1,numwfmem),stat=stat)
  if (stat.ne.0) error stop "Could not allocate cloudsh"
  allocate(clouds(0:nxmax-1,0:nymax-1,nzmax,numwfmem),stat=stat)
  if (stat.ne.0) error stop "Could not allocate clouds"

  ! 2d fields
  !**********
  allocate(ps(0:nxmax-1,0:nymax-1,1,numwfmem),stat=stat)
  if (stat.ne.0) error stop "Could not allocate ps"
  allocate(sd(0:nxmax-1,0:nymax-1,1,numwfmem),stat=stat)
  if (stat.ne.0) error stop "Could not allocate sd"
  allocate(msl(0:nxmax-1,0:nymax-1,1,numwfmem),stat=stat)
  if (stat.ne.0) error stop "Could not allocate msl"
  allocate(tcc(0:nxmax-1,0:nymax-1,1,numwfmem),stat=stat)
  if (stat.ne.0) error stop "Could not allocate tcc"
  allocate(u10(0:nxmax-1,0:nymax-1,1,numwfmem),stat=stat)
  if (stat.ne.0) error stop "Could not allocate u10"
  allocate(v10(0:nxmax-1,0:nymax-1,1,numwfmem),stat=stat)
  if (stat.ne.0) error stop "Could not allocate v10"
  allocate(tt2(0:nxmax-1,0:nymax-1,1,numwfmem),stat=stat)
  if (stat.ne.0) error stop "Could not allocate tt2"
  allocate(td2(0:nxmax-1,0:nymax-1,1,numwfmem),stat=stat)
  if (stat.ne.0) error stop "Could not allocate td2"
  allocate(lsprec(0:nxmax-1,0:nymax-1,1,numwfmem),stat=stat)
  if (stat.ne.0) error stop "Could not allocate lsprec"
  allocate(convprec(0:nxmax-1,0:nymax-1,1,numwfmem),stat=stat)
  if (stat.ne.0) error stop "Could not allocate convprec"
  allocate(sshf(0:nxmax-1,0:nymax-1,1,numwfmem),stat=stat)
  if (stat.ne.0) error stop "Could not allocate sshf"
  allocate(ssr(0:nxmax-1,0:nymax-1,1,numwfmem),stat=stat)
  if (stat.ne.0) error stop "Could not allocate ssr"
  allocate(sfcstress(0:nxmax-1,0:nymax-1,1,numwfmem),stat=stat)
  if (stat.ne.0) error stop "Could not allocate sfcstress"
  allocate(ustar(0:nxmax-1,0:nymax-1,1,numwfmem),stat=stat)
  if (stat.ne.0) error stop "Could not allocate ustar"
  allocate(wstar(0:nxmax-1,0:nymax-1,1,numwfmem),stat=stat)
  if (stat.ne.0) error stop "Could not allocate wstar"
  allocate(hmix(0:nxmax-1,0:nymax-1,1,numwfmem),stat=stat)
  if (stat.ne.0) error stop "Could not allocate hmix"
  allocate(tropopause(0:nxmax-1,0:nymax-1,1,numwfmem),stat=stat)
  if (stat.ne.0) error stop "Could not allocate tropopause"
  allocate(oli(0:nxmax-1,0:nymax-1,1,numwfmem),stat=stat)
  if (stat.ne.0) error stop "Could not allocate oli"

  ! Vertical descritisation arrays
  !*******************************
  allocate(height(nzmax),wheight(nzmax),uvheight(nzmax),stat=stat)
  if (stat.ne.0) error stop "Could not allocate height arrays"
  allocate(akm(nwzmax),bkm(nwzmax),akz(nuvzmax),bkz(nuvzmax), &
    aknew(nzmax),bknew(nzmax),stat=stat)
  if (stat.ne.0) error stop "Could not allocate model level parameters"
end subroutine alloc_windfields

subroutine alloc_windfields_nest
  !*******************************************************************************    
  ! Dynamic allocation of arrays
  !
  ! For nested wind fields. 
  ! 
  !*******************************************************************************
  implicit none 
  integer :: stat

  allocate(uun(0:nxmaxn-1,0:nymaxn-1,nzmax,numwfmem,numbnests),stat=stat)
  if (stat.ne.0) error stop "Could not allocate uun"
  allocate(vvn(0:nxmaxn-1,0:nymaxn-1,nzmax,numwfmem,numbnests),stat=stat)
  if (stat.ne.0) error stop "Could not allocate vvn"
  allocate(wwn(0:nxmaxn-1,0:nymaxn-1,nzmax,numwfmem,numbnests),stat=stat)
  if (stat.ne.0) error stop "Could not allocate wwn"
  allocate(ttn(0:nxmaxn-1,0:nymaxn-1,nzmax,numwfmem,numbnests),stat=stat)
  if (stat.ne.0) error stop "Could not allocate ttn"
  allocate(qvn(0:nxmaxn-1,0:nymaxn-1,nzmax,numwfmem,numbnests),stat=stat)
  if (stat.ne.0) error stop "Could not allocate qvn"
  allocate(pvn(0:nxmaxn-1,0:nymaxn-1,nzmax,numwfmem,numbnests),stat=stat)
  if (stat.ne.0) error stop "Could not allocate pvn"
  allocate(clwcn(0:nxmaxn-1,0:nymaxn-1,nzmax,numwfmem,numbnests),stat=stat)
  if (stat.ne.0) error stop "Could not allocate clwcn"
  allocate(ciwcn(0:nxmaxn-1,0:nymaxn-1,nzmax,numwfmem,numbnests),stat=stat)
  if (stat.ne.0) error stop "Could not allocate ciwcn"
  allocate(clwn(0:nxmaxn-1,0:nymaxn-1,nzmax,numwfmem,numbnests),stat=stat)
  if (stat.ne.0) error stop "Could not allocate clwn"

  ! ETA equivalents
#ifdef ETA
  allocate(uuetan(0:nxmaxn-1,0:nymaxn-1,nzmax,numwfmem,numbnests),stat=stat)
  if (stat.ne.0) error stop "Could not allocate uuetan"
  allocate(vvetan(0:nxmaxn-1,0:nymaxn-1,nzmax,numwfmem,numbnests),stat=stat)
  if (stat.ne.0) error stop "Could not allocate vvetan"
  allocate(wwetan(0:nxmaxn-1,0:nymaxn-1,nzmax,numwfmem,numbnests),stat=stat)
  if (stat.ne.0) error stop "Could not allocate wwetan"
  allocate(ttetan(0:nxmaxn-1,0:nymaxn-1,nzmax,numwfmem,numbnests),stat=stat)
  if (stat.ne.0) error stop "Could not allocate ttetan"
  allocate(pvetan(0:nxmaxn-1,0:nymaxn-1,nzmax,numwfmem,numbnests),stat=stat)
  if (stat.ne.0) error stop "Could not allocate pvetan"
  allocate(prsetan(0:nxmaxn-1,0:nymaxn-1,nzmax,numwfmem,numbnests) ,stat=stat)
  if (stat.ne.0) error stop "Could not allocate prsetan"
  allocate(rhoetan(0:nxmaxn-1,0:nymaxn-1,nzmax,numwfmem,numbnests),stat=stat)
  if (stat.ne.0) error stop "Could not allocate rhoetan"
  allocate(drhodzetan(0:nxmaxn-1,0:nymaxn-1,nzmax,numwfmem,numbnests),stat=stat)
  if (stat.ne.0) error stop "Could not allocate drhodzetan"
#endif
  allocate(etauvheightn(0:nxmaxn-1,0:nymaxn-1,nuvzmax,numwfmem,numbnests),stat=stat)
  if (stat.ne.0) error stop "Could not allocate etauvheightn"
  allocate(etawheightn(0:nxmaxn-1,0:nymaxn-1,nuvzmax,numwfmem,numbnests),stat=stat)
  if (stat.ne.0) error stop "Could not allocate etawheightn"

  allocate(cloudsn(0:nxmaxn-1,0:nymaxn-1,nzmax,numwfmem,numbnests),stat=stat)
  if (stat.ne.0) error stop "Could not allocate cloudsn"
  allocate(cloudshn(0:nxmaxn-1,0:nymaxn-1,numwfmem,numbnests),stat=stat)
  if (stat.ne.0) error stop "Could not allocate cloudshn"
  allocate(prsn(0:nxmaxn-1,0:nymaxn-1,nzmax,numwfmem,numbnests),stat=stat)
  if (stat.ne.0) error stop "Could not allocate prsn"
  allocate(rhon(0:nxmaxn-1,0:nymaxn-1,nzmax,numwfmem,numbnests),stat=stat)
  if (stat.ne.0) error stop "Could not allocate rhon"
  allocate(drhodzn(0:nxmaxn-1,0:nymaxn-1,nzmax,numwfmem,numbnests),stat=stat)
  if (stat.ne.0) error stop "Could not allocate drhodzn"
  allocate(tthn(0:nxmaxn-1,0:nymaxn-1,nuvzmax,numwfmem,numbnests),stat=stat)
  if (stat.ne.0) error stop "Could not allocate tthn"
  allocate(qvhn(0:nxmaxn-1,0:nymaxn-1,nuvzmax,numwfmem,numbnests),stat=stat)
  if (stat.ne.0) error stop "Could not allocate qvhn"
  allocate(clwchn(0:nxmaxn-1,0:nymaxn-1,nuvzmax,numwfmem,numbnests),stat=stat)
  if (stat.ne.0) error stop "Could not allocate clwchn"
  allocate(ciwchn(0:nxmaxn-1,0:nymaxn-1,nuvzmax,numwfmem,numbnests),stat=stat)
  if (stat.ne.0) error stop "Could not allocate ciwchn"
  allocate(ctwcn(0:nxmaxn-1,0:nymaxn-1,numwfmem,numbnests),stat=stat)
  if (stat.ne.0) error stop "Could not allocate ctwcn"

  ! 2d fields
  !***********
  allocate(psn(0:nxmaxn-1,0:nymaxn-1,1,numwfmem,numbnests),stat=stat)
  if (stat.ne.0) error stop "Could not allocate psn"
  allocate(sdn(0:nxmaxn-1,0:nymaxn-1,1,numwfmem,numbnests),stat=stat)
  if (stat.ne.0) error stop "Could not allocate sdn"
  allocate(msln(0:nxmaxn-1,0:nymaxn-1,1,numwfmem,numbnests),stat=stat)
  if (stat.ne.0) error stop "Could not allocate msln"
  allocate(tccn(0:nxmaxn-1,0:nymaxn-1,1,numwfmem,numbnests),stat=stat)
  if (stat.ne.0) error stop "Could not allocate tccn"
  allocate(u10n(0:nxmaxn-1,0:nymaxn-1,1,numwfmem,numbnests),stat=stat)
  if (stat.ne.0) error stop "Could not allocate u10n"
  allocate(v10n(0:nxmaxn-1,0:nymaxn-1,1,numwfmem,numbnests),stat=stat)
  if (stat.ne.0) error stop "Could not allocate v10n"
  allocate(tt2n(0:nxmaxn-1,0:nymaxn-1,1,numwfmem,numbnests),stat=stat)
  if (stat.ne.0) error stop "Could not allocate tt2n"
  allocate(td2n(0:nxmaxn-1,0:nymaxn-1,1,numwfmem,numbnests),stat=stat)
  if (stat.ne.0) error stop "Could not allocate td2n"
  allocate(lsprecn(0:nxmaxn-1,0:nymaxn-1,1,numwfmem,numbnests),stat=stat)
  if (stat.ne.0) error stop "Could not allocate lsprecn"
  allocate(convprecn(0:nxmaxn-1,0:nymaxn-1,1,numwfmem,numbnests),stat=stat)
  if (stat.ne.0) error stop "Could not allocate convprecn"
  allocate(sshfn(0:nxmaxn-1,0:nymaxn-1,1,numwfmem,numbnests),stat=stat)
  if (stat.ne.0) error stop "Could not allocate sshfn"
  allocate(ssrn(0:nxmaxn-1,0:nymaxn-1,1,numwfmem,numbnests),stat=stat)
  if (stat.ne.0) error stop "Could not allocate ssrn"
  allocate(sfcstressn(0:nxmaxn-1,0:nymaxn-1,1,numwfmem,numbnests),stat=stat)
  if (stat.ne.0) error stop "Could not allocate sfcstressn"
  allocate(ustarn(0:nxmaxn-1,0:nymaxn-1,1,numwfmem,numbnests),stat=stat)
  if (stat.ne.0) error stop "Could not allocate ustarn"
  allocate(wstarn(0:nxmaxn-1,0:nymaxn-1,1,numwfmem,numbnests),stat=stat)
  if (stat.ne.0) error stop "Could not allocate wstarn"
  allocate(hmixn(0:nxmaxn-1,0:nymaxn-1,1,numwfmem,numbnests),stat=stat)
  if (stat.ne.0) error stop "Could not allocate hmixn"
  allocate(tropopausen(0:nxmaxn-1,0:nymaxn-1,1,numwfmem,numbnests),stat=stat)
  if (stat.ne.0) error stop "Could not allocate tropopausen"
  allocate(olin(0:nxmaxn-1,0:nymaxn-1,1,numwfmem,numbnests),stat=stat)
  if (stat.ne.0) error stop "Could not allocate olin"
  allocate(vdepn(0:nxmaxn-1,0:nymaxn-1,maxspec,numwfmem,numbnests),stat=stat)
  if (stat.ne.0) error stop "Could not allocate vdepn"

  ! Initialise
  !************  
  clwcn(:,:,:,:,:)=0.
  ciwcn(:,:,:,:,:)=0.
  clwchn(:,:,:,:,:)=0.
  ciwchn(:,:,:,:,:)=0.
end subroutine alloc_windfields_nest

subroutine alloc_nest_properties
  implicit none 
  integer :: stat

  allocate(nxn(numbnests),stat=stat)
  if (stat.ne.0) error stop "Could not allocate nxn"
  allocate(nyn(numbnests),stat=stat)
  if (stat.ne.0) error stop "Could not allocate nyn"
  allocate(dxn(numbnests),stat=stat)
  if (stat.ne.0) error stop "Could not allocate dxn"
  allocate(dyn(numbnests),stat=stat)
  if (stat.ne.0) error stop "Could not allocate dyn"
  allocate(xlon0n(numbnests),stat=stat)
  if (stat.ne.0) error stop "Could not allocate xlon0n"
  allocate(ylat0n(numbnests),stat=stat)
  if (stat.ne.0) error stop "Could not allocate ylat0n"

  allocate(xresoln(0:numbnests),stat=stat)
  if (stat.ne.0) error stop "Could not allocate xresoln"
  allocate(yresoln(0:numbnests),stat=stat)
  if (stat.ne.0) error stop "Could not allocate yresoln"
  allocate(xln(numbnests),stat=stat)
  if (stat.ne.0) error stop "Could not allocate xln"
  allocate(yln(numbnests),stat=stat)
  if (stat.ne.0) error stop "Could not allocate yln"
  allocate(xrn(numbnests),stat=stat)
  if (stat.ne.0) error stop "Could not allocate xrn"
  allocate(yrn(numbnests),stat=stat)
  if (stat.ne.0) error stop "Could not allocate yrn"
end subroutine alloc_nest_properties

subroutine dealloc_windfields_nest
  
  deallocate(wfnamen)

  deallocate(nxn,nyn,dxn,dyn,xlon0n,ylat0n)

  deallocate(oron,excessoron,lsmn)

  deallocate(uun,vvn,wwn,ttn,qvn,pvn,clwcn,ciwcn,clwn,cloudsn, &
    cloudshn,rhon,prsn,drhodzn,tthn,qvhn,clwchn,ciwchn,ctwcn)

#ifdef ETA
  deallocate(uuetan,vvetan,wwetan,ttetan,pvetan,prsetan,rhoetan, &
    drhodzetan,etauvheightn,etawheightn)
#endif

  deallocate(psn,sdn,msln,tccn,u10n,v10n,tt2n,td2n,lsprecn,convprecn, &
    sshfn,ssrn,sfcstressn,ustarn,wstarn,hmixn,tropopausen,olin,vdepn)

  deallocate(xresoln,yresoln,xln,yln,xrn,yrn)
end subroutine dealloc_windfields_nest

subroutine dealloc_windfields
  implicit none

  deallocate(wftime,wfname)
  deallocate(oro,excessoro,lsm)

#ifdef ETA
  deallocate(uueta,vveta,wweta,uupoleta,vvpoleta,tteta,pveta, &
    prseta,rhoeta,drhodzeta)
#endif

  deallocate(etauvheight,etawheight)
  
  deallocate(uu,vv,ww,uupol,vvpol,tt,tth,qv,qvh,pv,rho,drhodz,pplev,prs,rho_dry)

  deallocate(clwc,ciwc,clw,clwch,ciwch,ctwc,cloudsh,clouds)

  deallocate(ps,sd,msl,tcc,u10,v10,tt2,td2,lsprec,convprec,sshf,ssr,sfcstress, &
    ustar,wstar,hmix,tropopause,oli)

  deallocate(height,wheight,uvheight,akm,bkm,akz,bkz,aknew,bknew)
end subroutine dealloc_windfields

end module windfields_mod
