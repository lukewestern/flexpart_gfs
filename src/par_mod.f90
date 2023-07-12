! SPDX-FileCopyrightText: FLEXPART 1998-2019, see flexpart_license.txt
! SPDX-License-Identifier: GPL-3.0-or-later

!*******************************************************************************
!   Include file for calculation of particle trajectories (Program FLEXPART)   *
!        This file contains the parameter statements used in FLEXPART          *
!                                                                              *
!        Author: A. Stohl                                                      *
!                                                                              *
!        1997                                                                  *
!                                                                              *
!        Update 15 August 2013 IP                                              *
!                                                                              *
!        Anne Tipka, Petra Seibert, 2021-02: implement new interpolation       *
!           for precipitation according to #295 using 2 additional fields      *
!                                                                              *
!*******************************************************************************

module par_mod

  implicit none

  !****************************************************************
  ! Parameters defining KIND parameter for double/single precision
  !****************************************************************

  integer,parameter :: dp=selected_real_kind(P=15)
  integer,parameter :: sp=selected_real_kind(6)

  !****************************************************************
  ! dep_prec sets the precision for deposition calculations (sp or 
  ! dp). sp is default, dp can be used for increased precision.
  !****************************************************************

  integer,parameter :: dep_prec=dp

  !****************************************************************
  ! Set to F to disable use of kernel for concentrations/deposition
  !****************************************************************

  logical, parameter :: lusekerneloutput=.true.

  !*********************************************************************
  ! Set to T to change output units to number of particles per grid cell
  !*********************************************************************
  logical, parameter :: lparticlecountoutput=.false.

  !***********************************************************
  ! number of directories/files used for FLEXPART input/output
  !***********************************************************

  integer,parameter :: numpath=4

  ! numpath                 Number of different pathnames for input/output files


  !*****************************
  ! Physical and other constants
  !*****************************

  real,parameter :: pi=3.14159265, r_earth=6.371e6, r_air=287.05, ga=9.81
  real,parameter :: cpa=1004.6, kappa=0.286, pi180=pi/180., vonkarman=0.4
  ! additional constants RLT Aug-2017
  real,parameter :: rgas=8.31447 
  real,parameter :: r_water=461.495

  ! pi                      number "pi"
  ! pi180                   pi/180.
  ! r_earth                 radius of earth [m]
  ! r_air                   individual gas constant for dry air [J/kg/K]
  ! ga                      gravity acceleration of earth [m/s**2]
  ! cpa                     specific heat for dry air
  ! kappa                   exponent of formula for potential temperature
  ! vonkarman               von Karman constant
  ! rgas                    universal gas constant [J/mol/K]
  ! r_water                 specific gas constant for water vapor [J/kg/K]

  real,parameter :: karman=0.40, href=15., convke=2.0
  real,parameter :: hmixmin=100., hmixmax=4500.
  !real,parameter :: d_trop=50., d_strat=0.1
  real :: d_trop=50., d_strat=0.1, fturbmeso=0.16 ! turbulence factors can change for different runs
  real,parameter :: rho_water=1000. !ZHG 2015 [kg/m3]
  real,parameter :: ratio_incloud=6.2   !ZHG MAR2016
  real,parameter :: wet_a=1.e-5, wet_b=0.8 !AT

  ! karman            Karman's constant
  ! href [m]          Reference height for dry deposition
  ! konvke            Relative share of kinetic energy used for parcel lifting
  ! hmixmin,hmixmax   Minimum and maximum allowed PBL height
  ! fturbmeso         the factor by which standard deviations of winds at grid
  !                   points surrounding the particle positions are scaled to
  !                   yield the scales for the mesoscale wind velocity fluctuations
  ! d_trop [m2/s]     Turbulent diffusivity for horiz components in the troposphere
  ! d_strat [m2/s]    Turbulent diffusivity for vertical component in the stratosphere
  ! ratio_incloud     ZHG MAR2016
  ! wet_a, wet_b      for wetscav=wet_a*prec**wet_b if no cloud found, but precipitation occurs

  
  real,parameter :: xmwml=18.016/28.960
                  ! ratio of molar weights of water vapor and dry air


  !****************************************************
  ! Constants related to the stratospheric ozone tracer
  !****************************************************

  real,parameter :: ozonescale=60., pvcrit=2.

  ! ozonescale              ppbv O3 per PV unit
  ! pvcrit                  PV level of the tropopause



  !********************
  ! Some time constants
  !********************

  integer,parameter :: idiffnorm=10800, idiffmax=2*idiffnorm, minstep=1

  ! idiffnorm [s]           normal time interval between two wind fields
  ! idiffmax [s]            maximum time interval between two wind fields
  ! minstep [s]             minimum time step to be used within FLEXPART


  !*****************************************************************
  ! Parameters for polar stereographic projection close to the poles
  !*****************************************************************

  real,parameter :: switchnorth=75., switchsouth=-75.

  ! switchnorth    use polar stereographic grid north of switchnorth
  ! switchsouth    use polar stereographic grid south of switchsouth


  !*********************************************
  ! Maximum dimensions of the input mother grids
  !*********************************************
  
  ! ECMWF
! integer,parameter :: nxmax=361,nymax=181,nuvzmax=92,nwzmax=92,nzmax=92,nxshift=359 ! 1.0 deg 92 levels
!  integer,parameter :: nxmax=361,nymax=181,nuvzmax=138,nwzmax=138,nzmax=138,nxshift=0 ! 1.0 deg 138 levels
!   integer,parameter :: nxmax=361,nymax=181,nuvzmax=138,nwzmax=138,nzmax=138,nxshift=359 ! 1.0 deg 138 levels
  integer,parameter :: nxmax=721,nymax=361,nuvzmax=138,nwzmax=138,nzmax=138!,nxshift=359  ! 0.5 deg 138 levels
!  integer,parameter :: nxmax=181,nymax=91,nuvzmax=92,nwzmax=92,nzmax=92,nxshift=0  ! CERA 2.0 deg 92 levels

! GFS
!   integer,parameter :: nxmax=361,nymax=181,nuvzmax=138,nwzmax=138,nzmax=138
!   integer,parameter :: nxshift=0 ! shift not fixed for the executable 

  !*********************************
  ! Parmaters for GRIB file decoding
  !*********************************

  ! integer,parameter :: jpack=4*nxmax*nymax, jpunp=4*jpack
  integer,parameter :: jpack=4*361*181, jpunp=4*jpack
  ! jpack,jpunp             maximum dimensions needed for GRIB file decoding


  !*********************************************
  ! Maximum dimensions of the nested input grids
  !*********************************************

  integer,parameter :: maxnests=0,nxmaxn=0,nymaxn=0

  ! nxmax,nymax        maximum dimension of wind fields in x and y
  !                    direction, respectively
  ! nuvzmax,nwzmax     maximum dimension of (u,v) and (w) wind fields in z
  !                    direction (for fields on eta levels)
  ! nzmax              maximum dimension of wind fields in z direction
  !                    for the transformed Cartesian coordinates

  
  integer,parameter :: nconvlevmax = nuvzmax-1
  integer,parameter :: na = nconvlevmax+1

  ! ntracermax         maximum number of tracer species in convection
  ! nconvlevmax        maximum number of levels for convection
  ! na                 parameter used in Emanuel's convect subroutine


  !**************************************
  ! Maximum dimensions of the output grid
  !**************************************

  !integer,parameter :: maxageclass=1,maxzgrid=10,nclassunc=1
  integer,parameter :: maxageclass=1,nclassunc=1

  ! nclassunc               number of classes used to calculate the uncertainty
  !                         of the output
  ! maxageclass             maximum number of age classes used for output

  ! Sabine Eckhardt, June, 2008
  ! the dimensions of the OUTGRID are now set dynamically during runtime
  ! maxxgrid,maxygrid,maxzgrid    maximum dimensions in x,y,z direction
  ! maxxgridn,maxygridn           maximum dimension of the nested grid
  !integer maxxgrid,maxygrid,maxzgrid,maxxgridn,maxygridn
  !integer,parameter :: maxxgrid=361,maxygrid=181,maxxgridn=0,maxygridn=0)

  integer,parameter :: maxreceptor=20

  ! maxreceptor             maximum number of receptor points


  !**************************************************
  ! Maximum number of particles, species, and similar
  !**************************************************

  !integer,parameter :: maxpart=5000000
  integer,parameter :: maxspec=5

  real,parameter :: minmassfrac=0.0

  ! maxpart      Maximum number of particles
  ! maxspec      Maximum number of chemical species per release
  ! minmassfrac  Terminate particles carrying a lower fraction
  !                compared to their initial mass

  ! maxpoint is also set dynamically during runtime
  ! maxpoint     Maximum number of release locations

  ! ---------
  ! Sabine Eckhardt: change of landuse inventary numclass=13
  integer,parameter :: maxwf=1000000, maxtable=1000, numclass=13, maxndia=100
  integer,parameter :: numpf=1 ! number of precip fields original =1, new=3(AT and PS, #295)
  integer,parameter :: numwfmem=2 ! Serial version/MPI with 2 fields
  !integer,parameter :: numwfmem=3 ! MPI with 3 fields

  ! maxwf        maximum number of wind fields to be used for simulation
  ! maxtable     Maximum number of chemical species that can be tabulated 
  ! numclass     Number of landuse classes available to FLEXPART
  ! maxndia      Maximum number of diameter classes of particles
  ! numpf        Number of precipitation fields (1 standard, 3 #295)
  ! numwfmem     Number of windfields kept in memory. 2 for serial version, 
  !              2 or 3 for MPI version

  !**************************************************************************
  ! dimension of the OH field
  !**************************************************************************
  
  integer,parameter :: maxxOH=72, maxyOH=46, maxzOH=7
  
  !**************************************************************************
  ! aerosol below-cloud scavenging removal polynomial constants for rain & snow
  !**************************************************************************

  ! for bcscheme = 1
  ! rain (Laakso et al 2003)
  real, parameter :: bclr(6) = &
      (/274.35758, 332839.59273, 226656.57259, 58005.91340, 6588.38582, 0.244984/)
  ! snow (Kyro et al 2009)
  real, parameter :: bcls(6) = (/22.7, 0.0, 0.0, 1321.0, 381.0, 0.0/) 
  
  ! for bcscheme = 2 & 3
  ! AT (after Wang et al 2014)
  ! rain
  real, parameter :: bclr_a(4) = &
      (/-6.2609, 0.682, 0.8676, 0.1282/)
  real, parameter :: bclr_b(7) = &
      (/-14.707, 51.043, -97.306, 97.946, -53.923, 15.311, -1.751/)
  real, parameter :: bclr_c(2) = &
      (/0.723, 0.0303/)
  real, parameter :: bclr_e(7) = &
      (/-0.6492, 9.3483, -21.929, 25.317, -15.395, 4.7242, -0.5766/)
  ! snow 
  real, parameter :: bcls_a(7) = &
      (/-4.426, 1.394, -1.202, -3.2942, -1.9521, -0.4904, -0.0457/)
  real, parameter :: bcls_b(7) = &
      (/-4.3521, -0.7828, 12.768, -19.864, 13.618, -4.4350, 0.5551/)
  real, parameter :: bcls_c(7) = &
      (/0.5664, 0.0085, -0.1948, -0.6532, -0.5462, -0.1778, -0.0201/) 
  real, parameter :: bcls_e(7) = &
      (/0.5689, -0.0923, 0.0402, 1.4523, -2.078, 1.05, -0.1821/)


  !**************************************************************************
  ! Maximum number of particles to be released in a single atmospheric column
  ! for the domain-filling trajectories option
  !**************************************************************************

  integer,parameter :: maxcolumn=3000


  !*********************************
  ! Dimension of random number field
  !*********************************

  integer,parameter :: maxrand=6000000

  ! maxrand                 number of random numbers used
  

  !*****************************************************
  ! Number of clusters to be used for plume trajectories
  !*****************************************************

  integer,parameter :: ncluster=5

  !************************************
  ! Unit numbers for input/output files
  !************************************

  integer,parameter :: unitpath=1, unitcommand=1, unitageclasses=1, unitgrid=1
  integer,parameter :: unitavailab=1, unitreleases=88, unitpartout=93
  integer,parameter :: unitpartout_average=105, unitpartoptions=106
  integer,parameter :: unitrestart=106,unitheightlevels=107
  integer,parameter :: unitpartin=93, unitflux=98, unitouttraj=96
  integer,parameter :: unitvert=1, unitoro=1, unitpoin=1, unitreceptor=1
  integer,parameter :: unitreceptorout=2  
  integer,parameter :: unitoutgrid=97, unitoutgridppt=99, unitoutinfo=1
  integer,parameter :: unitspecies=1, unitoutrecept=91, unitoutreceptppt=92
  integer,parameter :: unitlsm=1, unitsfcdata=1, unitland=1, unitwesely=1
  integer,parameter :: unitOH=1
  integer,parameter :: unitdates=94, unitheader=90,unitheader_txt=100
  integer,parameter :: unitshortpart=95, unitprecip=101
  integer,parameter :: unitboundcond=89
  integer,parameter :: unittmp=101
! RLT
  integer,parameter :: unitoutfactor=102

!******************************************************
! integer code for missing values, used in wet scavenging (PS, 2012)
!******************************************************

  integer,parameter ::  icmv=-9999

!*******************************************************************************
! Maximum output of each partoutput NetCDF-4 file in Mb 
! before a new one is created
!*******************************************************************************

  integer,parameter :: max_partoutput_filesize=30000

  ! Set maximum number of threads for doing grid computations. 
  ! Recommended to set this to max 16
  ! High numbers create more overhead and a larger memory footprint
  !***********************************************************************
  integer,parameter :: max_numthreads_grid=1
  ! Set the coordinate system. At the moment only ECMWF is possible. This bit
  ! needs to be a parameter that can be set at compile time. 
  ! Throughout the code there will be SELECT CASE statements or IFDEFs
  !*******************************************************************
  
  !character(len=256),parameter :: wind_coord_type='ETA'
  character(len=256),parameter :: wind_coord_type='METER'

  ! This flag sets all vertical interpolation to logarithmic instead of linear
  !***************************************************************************
  logical,parameter :: log_interpol=.false.

  ! mesoscale turbulence is found to give issues, so turned off by default
  !***********************************************************************
  logical,parameter :: lmesoscale_turb=.false.

  ! Threshold equivalent diameter for interaction with surface sublayer 
  ! resistance (below 10 meters) in micrometer. Above this diameter there
  ! is no interaction
  !**********************************************************************
  real,parameter :: d_thresheqv=20
  
end module par_mod
