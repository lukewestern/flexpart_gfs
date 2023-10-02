! SPDX-FileCopyrightText: FLEXPART 1998-2019, see flexpart_license.txt
! SPDX-License-Identifier: GPL-3.0-or-later

!*****************************************************************************
!                                                                            *
!   L. Bakels 2022: This module contains all subroutines for                 *
!                   reading option files                                     *
!                                                                            *
!*****************************************************************************

module readoptions_mod
  use par_mod
  use com_mod
  use date_mod
  use point_mod
  use windfields_mod

  implicit none

contains

subroutine readageclasses

  !*****************************************************************************
  !                                                                            *
  !     This routine reads the age classes to be used for the current model    *
  !     run.                                                                   *
  !                                                                            *
  !     Author: A. Stohl                                                       *
  !     20 March 2000                                                          *
  !     HSO, 1 July 2014                                                       *
  !     Added optional namelist input                                          *
  !                                                                            *
  !*****************************************************************************
  !                                                                            *
  ! Variables:                                                                 *
  !                                                                            *
  ! Constants:                                                                 *
  !                                                                            *
  !*****************************************************************************

  implicit none

  integer :: i

  ! namelist help variables
  integer :: ios,stat
  ! namelist declaration
  namelist /nage/ &
    nageclass
  namelist /ageclass/ &
    lage

  nageclass=-1 ! preset to negative value to identify failed namelist input

  ! If age spectra calculation is switched off, set number of age classes
  ! to 1 and maximum age to a large number
  !**********************************************************************

  if (lagespectra.ne.1) then
    nageclass=1
    allocate( lage(nageclass),stat=stat)
    if (stat.ne.0) error stop "Could not allocate lage"
    lage(nageclass)=999999999
    return
  endif

  ! If age spectra claculation is switched on,
  ! open the AGECLASSSES file and read user options
  !************************************************

  open(unitageclasses,file=path(1)(1:length(1))//'AGECLASSES', &
    form='formatted',status='old',err=999)

  ! try to read in as a namelist

  read(unitageclasses,nage,iostat=ios)
  allocate( lage(nageclass),stat=stat)
  if (stat.ne.0) error stop "Could not allocate lage"
  read(unitageclasses,ageclass,iostat=ios)
  close(unitageclasses)

  if ((nageclass.lt.0).or.(ios.ne.0)) then
    open(unitageclasses,file=path(1)(1:length(1))//'AGECLASSES', &
      status='old',err=999)
    do i=1,13
      read(unitageclasses,*)
    end do
    read(unitageclasses,*) nageclass
    allocate( lage(nageclass),stat=stat)
    if (stat.ne.0) error stop "Could not allocate lage"
    read(unitageclasses,*) lage(1)
    do i=2,nageclass
      read(unitageclasses,*) lage(i)
    end do
    close(unitageclasses)
  endif

  ! write ageclasses file in namelist format to output directory if requested
  if (nmlout.and.lroot) then
    open(unitageclasses,file=path(2)(1:length(2))//'AGECLASSES.namelist', &
      err=1000)
    write(unitageclasses,nml=nage)
    write(unitageclasses,nml=ageclass)
    close(unitageclasses)
  endif

  if (lage(1).le.0) then
    write(*,*) ' #### FLEXPART MODEL ERROR! AGE OF FIRST      #### '
    write(*,*) ' #### CLASS MUST BE GREATER THAN ZERO. CHANGE #### '
    write(*,*) ' #### SETTINGS IN FILE AGECLASSES.            #### '
    error stop 'First age class must be larger than zero'
  endif

  do i=2,nageclass
    if (lage(i).le.lage(i-1)) then
      write(*,*) ' #### FLEXPART MODEL ERROR! AGE CLASSES     #### '
      write(*,*) ' #### MUST BE GIVEN IN TEMPORAL ORDER.      #### '
      write(*,*) ' #### CHANGE SETTINGS IN FILE AGECLASSES.   #### '
      error stop 'Age classes must be in temporal order'
    endif
  end do

  return

999   write(*,*) ' #### FLEXPART MODEL ERROR! FILE "AGECLASSES" #### '
  write(*,*) ' #### CANNOT BE OPENED IN THE DIRECTORY       #### '
  write(*,'(a)') path(1)(1:length(1))
  error stop 'AGECLASSES cannot be opened'

1000  write(*,*) ' #### FLEXPART MODEL ERROR! FILE "AGECLASSES" #### '
  write(*,*) ' #### CANNOT BE OPENED IN THE DIRECTORY       #### '
  write(*,'(a)') path(2)(1:length(2))
  error stop 'AGECLASSES cannot be opened'
end subroutine readageclasses

subroutine readavailable

  !*****************************************************************************
  !                                                                            *
  !   This routine reads the dates and times for which windfields are          *
  !   available.                                                               *
  !                                                                            *
  !     Authors: A. Stohl                                                      *
  !                                                                            *
  !     6 February 1994                                                        *
  !     8 February 1999, Use of nested fields, A. Stohl                        *
  !                                                                            *
  !*****************************************************************************
  !                                                                            *
  ! Variables:                                                                 *
  ! bdate                beginning date as Julian date                         *
  ! beg                  beginning date for windfields                         *
  ! endl                  ending date for windfields                            *
  ! fname                filename of wind field, help variable                 *
  ! ideltas [s]          duration of modelling period                          *
  ! idiff                time difference between 2 wind fields                 *
  ! idiffnorm            normal time difference between 2 wind fields          *
  ! idiffmax [s]         maximum allowable time between 2 wind fields          *
  ! jul                  julian date, help variable                            *
  ! numbwf               actual number of wind fields                          *
  ! wfname(maxwf)        file names of needed wind fields                      *
  ! wfspec(maxwf)        file specifications of wind fields (e.g., if on disc) *
  ! wftime(maxwf) [s]times of wind fields relative to beginning time           *
  ! wfname1,wfspec1,wftime1 = same as above, but only local (help variables)   *
  !                                                                            *
  ! Constants:                                                                 *
  ! maxwf                maximum number of wind fields                         *
  ! unitavailab          unit connected to file AVAILABLE                      *
  !                                                                            *
  !*****************************************************************************

  implicit none

  integer :: i,idiff,ldat,ltim,wftime1(maxwf),numbwfn(maxnests),k
  integer :: wftime1n(maxnests,maxwf),wftimen(maxnests,maxwf)
  logical :: lwarntd=.true.
  real(kind=dp) :: jul,beg,endl
  character(len=255) :: fname,spec,wfname1(maxwf),wfspec1(maxwf)
  character(len=255) :: wfname1n(maxnests,maxwf)
  character(len=40) :: wfspec1n(maxnests,maxwf)


  ! Windfields are only used, if they are within the modelling period.
  ! However, 1 additional day at the beginning and at the end is used for
  ! interpolation. -> Compute beginning and ending date for the windfields.
  !************************************************************************

  if (ideltas.gt.0) then         ! forward trajectories
    beg=bdate-1._dp
    endl=bdate+real(ideltas,kind=dp)/86400._dp+real(idiffmax,kind=dp)/ &
         86400._dp
  else                           ! backward trajectories
    beg=bdate+real(ideltas,kind=dp)/86400._dp-real(idiffmax,kind=dp)/ &
         86400._dp
    endl=bdate+1._dp
  endif

  ! Open the wind field availability file and read available wind fields
  ! within the modelling period.
  !*********************************************************************

  open(unitavailab,file=path(4)(1:length(4)),status='old', &
       err=999)

  do i=1,3
    read(unitavailab,*)
  end do

  numbwf=0
100   read(unitavailab,'(i8,1x,i6,2(6x,a255))',end=99) &
           ldat,ltim,fname,spec
    jul=juldate(ldat,ltim)
    if ((jul.ge.beg).and.(jul.le.endl)) then
      numbwf=numbwf+1
      if (numbwf.gt.maxwf) then      ! check exceedance of dimension
       write(*,*) 'Number of wind fields needed is too great.'
       write(*,*) 'Reduce modelling period (file "COMMAND") or'
       write(*,*) 'reduce number of wind fields (file "AVAILABLE").'
       error stop 'Number of wind fields is larger than maxwf (set in par_mod)'
      endif

      wfname1(numbwf)=fname(1:index(fname,' '))
      wfspec1(numbwf)=spec
      wftime1(numbwf)=nint((jul-bdate)*86400._dp)
    endif
    goto 100       ! next wind field

99   continue

  close(unitavailab)

  ! Open the wind field availability file and read available wind fields
  ! within the modelling period (nested grids)
  !*********************************************************************

  do k=1,numbnests
  !print*,length(numpath+2*(k-1)+1),length(numpath+2*(k-1)+2),length(4),length(3)
  !print*,path(numpath+2*(k-1)+2)(1:length(numpath+2*(k-1)+2))
    open(unitavailab,file=path(numpath+2*(k-1)+2) &
         (1:length(numpath+2*(k-1)+2)),status='old',err=998)

    do i=1,3
      read(unitavailab,*)
    end do

    numbwfn(k)=0
700   read(unitavailab,'(i8,1x,i6,2(6x,a255))',end=699) ldat, &
           ltim,fname,spec
      jul=juldate(ldat,ltim)
      if ((jul.ge.beg).and.(jul.le.endl)) then
        numbwfn(k)=numbwfn(k)+1
        if (numbwfn(k).gt.maxwf) then      ! check exceedance of dimension
       write(*,*) 'Number of nested wind fields is too great.'
       write(*,*) 'Reduce modelling period (file "COMMAND") or'
       write(*,*) 'reduce number of wind fields (file "AVAILABLE").'
          error stop 'Number of nested wind fields is larger than maxwf (set in par_mod)'
        endif

        wfname1n(k,numbwfn(k))=fname
        wfspec1n(k,numbwfn(k))=spec
        wftime1n(k,numbwfn(k))=nint((jul-bdate)*86400._dp)
      endif
      goto 700       ! next wind field

699   continue

    close(unitavailab)
  end do


  ! Check wind field times of file AVAILABLE (expected to be in temporal order)
  !****************************************************************************

  if (numbwf.eq.0) then
    write(*,*) ' #### FLEXPART MODEL ERROR! NO WIND FIELDS    #### '
    write(*,*) ' #### AVAILABLE FOR SELECTED TIME PERIOD.     #### '
    error stop 'No wind fields available for selected time period'
  endif

  do i=2,numbwf
    if (wftime1(i).le.wftime1(i-1)) then
      write(*,*) 'FLEXPART ERROR: FILE AVAILABLE IS CORRUPT.'
      write(*,*) 'THE WIND FIELDS ARE NOT IN TEMPORAL ORDER.'
      write(*,*) 'PLEASE CHECK FIELD ',wfname1(i)
      error stop 'AVAILABLE file error'
    endif
  end do

  ! Check wind field times of file AVAILABLE for the nested fields
  ! (expected to be in temporal order)
  !***************************************************************

  do k=1,numbnests
    if (numbwfn(k).eq.0) then
      write(*,*) '#### FLEXPART MODEL ERROR! NO WIND FIELDS  ####'
      write(*,*) '#### AVAILABLE FOR SELECTED TIME PERIOD.   ####'
      error stop 'No nested wind fields available for selected time period'
    endif

    do i=2,numbwfn(k)
      if (wftime1n(k,i).le.wftime1n(k,i-1)) then
      write(*,*) 'FLEXPART ERROR: FILE AVAILABLE IS CORRUPT. '
      write(*,*) 'THE NESTED WIND FIELDS ARE NOT IN TEMPORAL ORDER.'
      write(*,*) 'PLEASE CHECK FIELD ',wfname1n(k,i)
      write(*,*) 'AT NESTING LEVEL ',k
      error stop 'nested AVAILABLE file error'
      endif
    end do

  end do


  ! For backward trajectories, reverse the order of the windfields
  !***************************************************************

  if (ideltas.ge.0) then
    do i=1,numbwf
      wfname(i)=wfname1(i)
      wfspec(i)=wfspec1(i)
      wftime(i)=wftime1(i)
    end do
    do k=1,numbnests
      do i=1,numbwfn(k)
        wfnamen(k,i)=wfname1n(k,i)
        wfspecn(k,i)=wfspec1n(k,i)
        wftimen(k,i)=wftime1n(k,i)
      end do
    end do
  else
    do i=1,numbwf
      wfname(numbwf-i+1)=wfname1(i)
      wfspec(numbwf-i+1)=wfspec1(i)
      wftime(numbwf-i+1)=wftime1(i)
    end do
    do k=1,numbnests
      do i=1,numbwfn(k)
        wfnamen(k,numbwfn(k)-i+1)=wfname1n(k,i)
        wfspecn(k,numbwfn(k)-i+1)=wfspec1n(k,i)
        wftimen(k,numbwfn(k)-i+1)=wftime1n(k,i)
      end do
    end do
  endif

  ! Check the time difference between the wind fields. If it is big,
  ! write a warning message. If it is too big, terminate the trajectory.
  !*********************************************************************

  do i=2,numbwf
    idiff=abs(wftime(i)-wftime(i-1))
    if (idiff.gt.idiffmax.and.lroot) then
      write(*,*) 'FLEXPART WARNING: TIME DIFFERENCE BETWEEN TWO'
      write(*,*) 'WIND FIELDS IS TOO BIG FOR TRANSPORT CALCULATION.&
           &'
      write(*,*) 'THEREFORE, TRAJECTORIES HAVE TO BE SKIPPED.'
    else if (idiff.gt.idiffnorm.and.lroot.and.lwarntd) then
      write(*,*) 'FLEXPART WARNING: TIME DIFFERENCE BETWEEN TWO'
      write(*,*) 'WIND FIELDS IS BIG. THIS MAY CAUSE A DEGRADATION'
      write(*,*) 'OF SIMULATION QUALITY.'
      lwarntd=.false. ! only issue this warning once
    endif
  end do

  do k=1,numbnests
    if (numbwfn(k).ne.numbwf) then
      write(*,*) 'FLEXPART ERROR: THE AVAILABLE FILES FOR THE'
      write(*,*) 'NESTED WIND FIELDS ARE NOT CONSISTENT WITH'
      write(*,*) 'THE AVAILABLE FILE OF THE MOTHER DOMAIN.  '
      write(*,*) 'ERROR AT NEST LEVEL: ',k
      error stop 'Not the same number of nested wind field files as mother domain'
    endif
    do i=1,numbwf
      if (wftimen(k,i).ne.wftime(i)) then
        write(*,*) 'FLEXPART ERROR: THE AVAILABLE FILES FOR THE'
        write(*,*) 'NESTED WIND FIELDS ARE NOT CONSISTENT WITH'
        write(*,*) 'THE AVAILABLE FILE OF THE MOTHER DOMAIN.  '
        write(*,*) 'ERROR AT NEST LEVEL: ',k
        error stop 'Not the same time period of nested wind field files as mother domain'
      endif
    end do
  end do

  ! Reset the times of the wind fields that are kept in memory to no time
  !**********************************************************************

  do i=1,2
    memind(i)=i
    memtime(i)=999999999
  end do

  return

998   write(*,*) ' #### FLEXPART MODEL ERROR! AVAILABLE FILE   #### '
  write(*,'(a)') '     '//path(numpath+2*(k-1)+2) &
       (1:length(numpath+2*(k-1)+2))
  write(*,*) ' #### CANNOT BE OPENED             #### '
  error stop 'Nested AVAILABLE file cannot be opened'

999   write(*,*) ' #### FLEXPART MODEL ERROR! AVAILABLE FILE #### '
  write(*,'(a)') '     '//path(4)(1:length(4))
  write(*,*) ' #### CANNOT BE OPENED           #### '
  error stop 'AVAILABLE file cannot be opened'
end subroutine readavailable

subroutine readcommand

  !*****************************************************************************
  !                                                                            *
  !     This routine reads the user specifications for the current model run.  *
  !                                                                            *
  !     Author: A. Stohl                                                       *
  !                                                                            *
  !     18 May 1996                                                            *
  !     HSO, 1 July 2014                                                       *
  !     Added optional namelist input                                          *
  !                                                                            *
  !*****************************************************************************
  !                                                                            *
  ! Variables:                                                                 *
  ! bdate                beginning date as Julian date                         *
  ! ctl                  factor by which time step must be smaller than        *
  !                      Lagrangian time scale                                 *
  ! ibdate,ibtime        beginnning date and time (YYYYMMDD, HHMISS)           *
  ! ideltas [s]          modelling period                                      *
  ! iedate,ietime        ending date and time (YYYYMMDD, HHMISS)               *
  ! ifine                reduction factor for vertical wind time step          *
  ! outputforeachrel     for forward runs it is possible either to create      *
  !                      one outputfield or several for each releasepoint      *
  ! iflux                switch to turn on (1)/off (0) flux calculations       *
  ! iout                 1 for conc. (residence time for backward runs) output,*
  !                      2 for mixing ratio output, 3 both, 4 for plume        *
  !                      trajectory output, 5 = options 1 and 4                *
  ! ipin                 1 continue simulation with restart.bin file,          *
  !                      2 continue simulaion with dumped particle data, 0 no  *
  !                      3 use self-defined initial conditions in netcdf       *
  !                      4 initial run using option 3, restart from restart.bin*
  ! ipout                0 no particle dump, 1 every output time, 3 only at end*
  ! ipoutfac             increase particle dump interval by factor (default 1) *
  ! loutaver [s]         concentration output is an average over loutaver      *
  !                      seconds                                               *
  ! loutsample [s]       average is computed from samples taken every [s]      *
  !                      seconds                                               *
  ! loutstep [s]         time interval of concentration output                 *
  ! lsynctime [s]        synchronisation time interval for all particles       *
  ! lagespectra          switch to turn on (1)/off (0) calculation of age      *
  !                      spectra                                               *
  ! lconvection          value of either 0 and 1 indicating mixing by          *
  !                      convection                                            *
  !                      = 0 .. no convection                                  *
  !                      + 1 .. parameterisation of mixing by subgrid-scale    *
  !                              convection = on                               *
  ! lsubgrid             switch to turn on (1)/off (0) subgrid topography      *
  !                      parameterization                                      *
  ! method               method used to compute the particle pseudovelocities  *
  ! mdomainfill          1 use domain-filling option, 0 not, 2 use strat. O3   *
  !                                                                            *
  ! Constants:                                                                 *
  ! unitcommand          unit connected to file COMMAND                        *
  !                                                                            *
  !*****************************************************************************

  implicit none

  character(len=50) :: line
  integer :: ios
  integer :: lturbulence_meso

  namelist /command/ &
  ldirect, &
  ibdate,ibtime, &
  iedate,ietime, &
  loutstep, &
  loutaver, &
  loutsample, &
  loutrestart, &
  lsynctime, &
  ctl, &
  ifine, &
  iout, &
  ipout, &
  ipoutfac, &
  lsubgrid, &
  lconvection, &
  lturbulence, &
  lturbulence_meso, &
  lagespectra, &
  ipin, &
  ioutputforeachrelease, &
  iflux, &
  mdomainfill, &
  ind_source, &
  ind_receptor, &
  mquasilag, &
  nested_output, &
  linit_cond, &
  lnetcdfout, &
  sfc_only, &
  surf_only, &
  cblflag, &
  linversionout, &
  ohfields_path, &
  d_trop, &
  d_strat, &
  nxshift, &
  maxthreadgrid, &
  maxfilesize, &
  logvertinterp

  ! Presetting namelist command
  ldirect=0
  ibdate=20000101
  ibtime=0
  iedate=20000102
  ietime=0
  loutstep=10800
  loutaver=10800
  loutsample=900
  loutrestart=-1
  lsynctime=900
  ctl=-5.0
  ifine=4
  iout=3
  ipout=0
  ipoutfac=1
  lsubgrid=1
  lconvection=1
  lturbulence=1
  lturbulence_meso=0
  lagespectra=0
  ipin=0
  ioutputforeachrelease=1
  iflux=1
  mdomainfill=0
  ind_source=1
  ind_receptor=1
  mquasilag=0
  nested_output=0
  linit_cond=0
  lnetcdfout=1
  sfc_only=0
  surf_only=-1
  cblflag=0 ! if using old-style COMMAND file, set to 1 here to use mc cbl routine
  linversionout=0
  ohfields_path="../../flexin/"
  nxshift=-9999
  maxthreadgrid=1
  maxfilesize=10000
  logvertinterp=0

  !Af set release-switch
  WETBKDEP=.false.
  DRYBKDEP=.false.

  ! Open the command file and read user options
  ! Namelist input first: try to read as namelist file
  !**************************************************************************
  open(unitcommand,file=path(1)(1:length(1))//'COMMAND',status='old', &
    form='formatted',err=999)

  ! try namelist input (default)
  read(unitcommand,command,iostat=ios)
  if (ios.ne.0) then
        backspace(unitcommand)
        read(unitcommand,fmt='(A)') line
        if (lroot) write(*,*) &
            'Invalid line in COMMAND reads: '//trim(line)
  end if
  
  close(unitcommand)

  if (lturbulence_meso.ne.0) lmesoscale_turb=.true.
  if (logvertinterp.ne.0) log_interpol=.true.

  if (surf_only.ne.-1) then
    write(*,*) 'WARNING: SURF_ONLY in COMMAND will be deprecated and renamed SFC_ONLY'
    write(*,*) 'Continuing with SURF_ONLY...'
    sfc_only=surf_only
  endif
  ! distinguish namelist from fixed text input
  if ((ios.ne.0).or.(ldirect.eq.0)) then ! parse as text file format
    if (lroot) write(*,*) 'COMMAND either having unrecognised entries, &
      &or in old format, please update to namelist format.'
      error stop 'COMMAND has unrecognised entries or is not a namelist'
  endif ! input format

  ! write command file in namelist format to output directory if requested
  if (nmlout.and.lroot) then
    open(unitcommand,file=path(2)(1:length(2))//'COMMAND.namelist',err=1000)
    write(unitcommand,nml=command)
    close(unitcommand)
  endif

  ifine=max(ifine,1)

  ! Determine how Markov chain is formulated (for w or for w/sigw)
  !***************************************************************
  if (cblflag.eq.1) then ! added by mc to properly set parameters for CBL simulations 
    turbswitch=.true.
    if (lsynctime .gt. maxtl) lsynctime=maxtl  !maxtl defined in com_mod.f90
    if (ctl.lt.5) then
      print *,'WARNING: CBL flag active the ratio of TLu/dt has been set to 5'
      ctl=5.
    end if
    if (ifine*ctl.lt.50) then
      ifine=int(50./ctl)+1

      print *,'WARNING: CBL flag active the ratio of TLW/dt was < 50, &
        &ifine has been re-set to',ifine
  !pause
    endif
    print *,'WARNING: CBL flag active the ratio of TLW/dt is ',ctl*ifine
    print *,'WARNING: CBL flag active lsynctime is ',lsynctime
  else                    !added by mc
    if (ctl.ge.0.1) then
      turbswitch=.true.
    else
      turbswitch=.false.
      ifine=1
    endif
  endif                   !added by mc
  fine=1./real(ifine)
  ctl=1./ctl

  ! Set the switches required for the various options for input/output units
  !*************************************************************************
  !AF Set the switches IND_REL and IND_SAMP for the release and sampling
  !Af switches for the releasefile:
  !Af IND_REL =  1 : xmass * rho
  !Af IND_REL =  0 : xmass * 1

  !Af switches for the conccalcfile:
  !AF IND_SAMP =  0 : xmass * 1
  !Af IND_SAMP = -1 : xmass / rho

  !AF IND_SOURCE switches between different units for concentrations at the source
  !Af   NOTE that in backward simulations the release of computational particles
  !Af   takes place at the "receptor" and the sampling of particles at the "source".
  !Af          1 = mass units
  !Af          2 = mass mixing ratio units
  !Af IND_RECEPTOR switches between different units for concentrations at the receptor
  !            0 = no receptors
  !Af          1 = mass units
  !Af          2 = mass mixing ratio units
  !            3 = wet deposition in outputfield
  !            4 = dry deposition in outputfield

  if ( ldirect .eq. 1 ) then  ! FWD-Run
  !Af set release-switch
     if (ind_source .eq. 1 ) then !mass
        ind_rel = 0
     else ! mass mix
        ind_rel = 1
     endif
  !Af set sampling switch
     if (ind_receptor .le. 1) then !mass
        ind_samp = 0
     else ! mass mix
        ind_samp = -1
     endif
  elseif (ldirect .eq. -1 ) then !BWD-Run
  !Af set sampling switch
     if (ind_source .eq. 1 ) then !mass
        ind_samp = -1
     else ! mass mix
        ind_samp = 0
     endif
     select case (ind_receptor)
     case (1)  !  1 .. concentration at receptor
        ind_rel = 1
     case (2)  !  2 .. mixing ratio at receptor
        ind_rel = 0
     case (3)  ! 3 .. wet deposition in outputfield 
        ind_rel = 3
        if (lroot) then
          write(*,*) ' #### FLEXPART WET DEPOSITION BACKWARD MODE    #### '
          write(*,*) ' #### Releaseheight is forced to 0 - 20km      #### '
          write(*,*) ' #### Release is performed above ground lev    #### '
        end if
         WETBKDEP=.true.
         !allocate(xscav_frac1(maxpart,maxspec))
     case (4)  ! 4 .. dry deposition in outputfield
         ind_rel = 4
         if (lroot) then
           write(*,*) ' #### FLEXPART DRY DEPOSITION BACKWARD MODE    #### '
           write(*,*) ' #### Releaseheight is forced to 0 - 2*href    #### '
           write(*,*) ' #### Release is performed above ground lev    #### '
         end if
         DRYBKDEP=.true.
         !allocate(xscav_frac1(maxpart,maxspec))
     end select
  endif

  !*************************************************************
  ! Check whether valid options have been chosen in file COMMAND
  !*************************************************************

  ! Check options for initial condition output: Switch off for forward runs
  !************************************************************************

  if (ldirect.eq.1) linit_cond=0
  if ((linit_cond.lt.0).or.(linit_cond.gt.2)) then
    write(*,*) ' #### FLEXPART MODEL ERROR! INVALID OPTION    #### '
    write(*,*) ' #### FOR LINIT_COND IN FILE "COMMAND".       #### '
    error stop 'LINIT_COND in COMMAND is invalid'
  endif

  ! Check input dates
  !******************

  if (iedate.lt.ibdate) then
    write(*,*) ' #### FLEXPART MODEL ERROR! BEGINNING DATE    #### '
    write(*,*) ' #### IS LARGER THAN ENDING DATE. CHANGE      #### '
    write(*,*) ' #### EITHER POINT 2 OR POINT 3 IN FILE       #### '
    write(*,*) ' #### "COMMAND".                              #### '
    error stop 'COMMAND: beginning date is larger than ending date'
  else if (iedate.eq.ibdate) then
    if (ietime.lt.ibtime) then
    write(*,*) ' #### FLEXPART MODEL ERROR! BEGINNING TIME    #### '
    write(*,*) ' #### IS LARGER THAN ENDING TIME. CHANGE      #### '
    write(*,*) ' #### EITHER POINT 2 OR POINT 3 IN FILE       #### '
    write(*,*) ' #### "COMMAND".                              #### '
    error stop 'COMMAND: beginning time is larger than ending time'
    endif
  endif

! #ifndef USE_NCF
!   if ((loutrestart.ne.-1).or.(ipin.ne.0)) then
!     write(*,*) ' WARNING: restart option set with intervals'
!     write(*,*) ' LOUTRESTART', loutrestart
!     write(*,*) ' not possible when using binary gridded output'
!     write(*,*) ' ==> RESTART FUNCTION SWITCHED OFF!'
!   endif
!   if (ipin.ne.0) then 
!     write(*,*) ' ERROR: restart option not possible using binary'
!     write(*,*) ' output.'
!     write(*,*) ' Please only use IPIN>0 when compiling and running using'
!     write(*,*) ' netcdf output. '
!   endif
! #else
!   if ((surf_only.eq.1).or.(linversionout.eq.1)) then
!     write(*,*) ' ERROR: NetCDF output for surface only or for inversions'
!     write(*,*) ' is not yet implemented. Please compile without NetCDF.'
!     error stop 'Surface only option is not supported for NetCDF'
!   endif
! #endif

  ! Determine kind of dispersion method
  !************************************

  if (ctl.gt.0.) then
    method=1
    mintime=minstep
  else
    method=0
    mintime=lsynctime
  endif

  ! Check for netcdf output switch
  !*******************************
  if (iout.gt.8) then
    lnetcdfout=1
    iout = iout -8
  endif
  if (lnetcdfout.eq.1) then
#ifndef USE_NCF
    write(*,*) 'WARNING: netcdf output not activated during compile time &
      &but switched on in COMMAND or set to default value 1.'
    write(*,*) 'Please recompile with netcdf library (`make [...] ncf=yes`) &
      &when requiring NetCDF output.'
    write(*,*) 'LNETCDFOUT set to 0.'
    lnetcdfout = 0
#endif
  else
#ifdef USE_NCF
    write(*,*) 'WARNING: Executable compiled using NetCDF libraries, but &
      &BINARY output is requested. If this was unintended, please add 8 &
      &to IOUT or set LOUTNETCDF=1 in the COMMAND file.'
#endif
  endif

  if ((lnetcdfout.eq.1).and.((sfc_only.eq.1).or.(linversionout.eq.1))) then
    write(*,*) ' WARNING: NetCDF output for surface only or for inversions'
    write(*,*) ' is not yet implemented. LNETCDFOUT set to 0.'
    lnetcdfout=0
  endif

#ifndef USE_NCF
  if (ipout.ne.0) then
    write(*,*) 'ERROR: NETCDF missing! Please recompile with the netcdf'
    write(*,*) 'library if you want the particle dump or set IPOUT=0.'
    error stop 'FLEXPART not compiled with NetCDF'
  endif
#endif

  ! Check whether a valid option for gridded model output has been chosen
  !**********************************************************************

  if (iout.eq.0) then
    write(*,*) 'WARNING: IOUT set to zero, no gridded information will be &
      &written to file'
  else if ((iout.lt.0).or.(iout.gt.5)) then
    write(*,*) ' #### FLEXPART MODEL ERROR! FILE COMMAND:     #### '
    write(*,*) ' #### IOUT MUST BE 1, 2, 3, 4 OR 5 FOR        #### '
    write(*,*) ' #### STANDARD FLEXPART OUTPUT OR  9 - 13     #### '
    write(*,*) ' #### FOR NETCDF OUTPUT                       #### '
    error stop
  endif

  !AF check consistency between units and volume mixing ratio
  if ( ((iout.eq.2).or.(iout.eq.3)).and. &
       (ind_source.gt.1 .or.ind_receptor.gt.1) ) then
    write(*,*) ' #### FLEXPART MODEL ERROR! FILE COMMAND:     #### '
    write(*,*) ' #### VOLUME MIXING RATIO ONLY SUPPORTED      #### '
    write(*,*) ' #### FOR MASS UNITS (at the moment)          #### '
    error stop
  endif


  ! For quasilag output for each release is forbidden
  !*****************************************************************************

  if ((ioutputforeachrelease.eq.1).and.(mquasilag.eq.1)) then
      write(*,*) '#### FLEXPART MODEL ERROR! FILE COMMAND:     ####'
      write(*,*) '#### OUTPUTFOREACHRELEASE AND QUASILAGRANGIAN####'
      write(*,*) '#### MODE IS NOT POSSIBLE   !                ####'
      error stop
  endif


  ! For quasilag backward is forbidden
  !*****************************************************************************

  if ((ldirect.lt.0).and.(mquasilag.eq.1)) then
      write(*,*) '#### FLEXPART MODEL ERROR! FILE COMMAND:     ####'
      write(*,*) '#### FOR BACKWARD RUNS, QUASILAGRANGIAN MODE ####'
      write(*,*) '#### IS NOT POSSIBLE   !                     ####'
      error stop
  endif


  ! For backward runs one releasefield for all releases makes no sense,
  ! For quasilag and domainfill ioutputforechrelease is forbidden
  !*****************************************************************************

  if ((ldirect.lt.0).and.(ioutputforeachrelease.eq.0)) then
      write(*,*) '#### FLEXPART MODEL ERROR! FILE COMMAND:     ####'
      write(*,*) '#### FOR BACKWARD RUNS, IOUTPUTFOREACHRLEASE ####'
      write(*,*) '#### MUST BE SET TO ONE!                     ####'
      error stop
  endif


  ! For backward runs one releasefield for all releases makes no sense,
  ! and is "forbidden"
  !*****************************************************************************

  if ((mdomainfill.eq.1).and.(ioutputforeachrelease.eq.1)) then
      write(*,*) '#### FLEXPART MODEL ERROR! FILE COMMAND:     ####'
      write(*,*) '#### FOR DOMAIN FILLING RUNS OUTPUT FOR      ####'
      write(*,*) '#### EACH RELEASE IS FORBIDDEN !             ####'
      error stop
  endif

  ! Inversion output format only for backward runs
  !*****************************************************************************
  
  if ((linversionout.eq.1).and.(ldirect.eq.1)) then
      write(*,*) '#### FLEXPART MODEL ERROR! FILE COMMAND:     ####'
      write(*,*) '#### INVERSION OUTPUT FORMAT ONLY FOR        ####'
      write(*,*) '#### BACKWARD RUNS                           ####'
      error stop
  endif


  ! For domain-filling trajectories, a plume centroid trajectory makes no sense,
  ! For backward runs, only residence time output (iout=1) or plume trajectories (iout=4),
  ! or both (iout=5) makes sense; other output options are "forbidden"
  !*****************************************************************************

  if (ldirect.lt.0) then
    if ((iout.eq.2).or.(iout.eq.3)) then
      write(*,*) '#### FLEXPART MODEL ERROR! FILE COMMAND:     ####'
      write(*,*) '#### FOR BACKWARD RUNS, IOUT MUST BE 1,4,OR 5####'
      error stop
    endif
  endif


  ! For domain-filling trajectories, a plume centroid trajectory makes no sense,
  ! and is "forbidden"
  !*****************************************************************************

  if (mdomainfill.ge.1) then
    if ((iout.eq.4).or.(iout.eq.5)) then
      write(*,*) '#### FLEXPART MODEL ERROR! FILE COMMAND:     ####'
      write(*,*) '#### FOR DOMAIN-FILLING TRAJECTORY OPTION,   ####'
      write(*,*) '#### IOUT MUST NOT BE SET TO 4 OR 5.         ####'
      error stop
    endif
  endif

  ! Check whether a valid options for particle dump has been chosen
  !****************************************************************

  if ((ipout.ne.0).and.(ipout.ne.1).and.(ipout.ne.2).and.(ipout.ne.3)) then
    write(*,*) ' #### FLEXPART MODEL ERROR! FILE COMMAND:     #### '
    write(*,*) ' #### IPOUT MUST BE 0, 1, 2 OR 3!             #### '
    error stop
  endif

  ! Check whether input and output settings don't contradict
  !*********************************************************
  ! if (((iout.eq.4).or.(iout.eq.5)).and.((ipin.eq.3).or.(ipin.eq.4))) then
  !   write(*,*) ' #### FLEXPART MODEL ERROR! FILE COMMAND:     #### '
  !   write(*,*) ' #### IOUT CANNOT BE 4 or 5 (plume) WHEN      #### '
  !   write(*,*) ' #### READING FROM part_ic.nc (ipin=4/5)      #### '
  !   error stop
  ! endif    

  if(lsubgrid.ne.1.and.verbosity.eq.0) then
    write(*,*) '             ----------------               '
    write(*,*) ' INFORMATION: SUBGRIDSCALE TERRAIN EFFECT IS'
    write(*,*) ' NOT PARAMETERIZED DURING THIS SIMULATION.  '
    write(*,*) '             ----------------               '
  endif


  ! Check whether convection scheme is either turned on or off
  !***********************************************************

  if ((lconvection.ne.0).and.(lconvection.ne.1)) then
    write(*,*) ' #### FLEXPART MODEL ERROR! FILE COMMAND:     #### '
    write(*,*) ' #### LCONVECTION MUST BE SET TO EITHER 1 OR 0#### '
    error stop
  endif


  ! Check whether synchronisation interval is sufficiently short
  !*************************************************************

  if (lsynctime.gt.(idiffnorm/2)) then
    write(*,*) ' #### FLEXPART MODEL ERROR! SYNCHRONISATION   #### '
    write(*,*) ' #### TIME IS TOO LONG. MAKE IT SHORTER.      #### '
    write(*,*) ' #### MINIMUM HAS TO BE: ', idiffnorm/2
    error stop
  endif


  ! Check consistency of the intervals, sampling periods, etc., for model output
  !*****************************************************************************

  if (loutaver.eq.0) then
    write(*,*) ' #### FLEXPART MODEL ERROR! TIME AVERAGE OF   #### '
    write(*,*) ' #### CONCENTRATION FIELD OUTPUT MUST NOT BE  #### '
    write(*,*) ' #### ZERO.                                   #### '
    write(*,*) ' #### CHANGE INPUT IN FILE COMMAND.           #### '
    error stop
  endif

  if (loutaver.gt.loutstep) then
    write(*,*) ' #### FLEXPART MODEL ERROR! TIME AVERAGE OF   #### '
    write(*,*) ' #### CONCENTRATION FIELD OUTPUT MUST NOT BE  #### '
    write(*,*) ' #### GREATER THAN INTERVAL OF OUTPUT.        #### '
    write(*,*) ' #### CHANGE INPUT IN FILE COMMAND.           #### '
    error stop
  endif

  if (loutsample.gt.loutaver) then
    write(*,*) ' #### FLEXPART MODEL ERROR! SAMPLING TIME OF  #### '
    write(*,*) ' #### CONCENTRATION FIELD OUTPUT MUST NOT BE  #### '
    write(*,*) ' #### GREATER THAN TIME AVERAGE OF OUTPUT.    #### '
    write(*,*) ' #### CHANGE INPUT IN FILE COMMAND.           #### '
    error stop
  endif

  if (mod(loutaver,lsynctime).ne.0) then
    write(*,*) ' #### FLEXPART MODEL ERROR! AVERAGING TIME OF #### '
    write(*,*) ' #### CONCENTRATION FIELD MUST BE A MULTIPLE  #### '
    write(*,*) ' #### OF THE SYNCHRONISATION INTERVAL         #### '
    error stop
  endif

  if ((loutaver/lsynctime).lt.2) then
    write(*,*) ' #### FLEXPART MODEL ERROR! AVERAGING TIME OF #### '
    write(*,*) ' #### CONCENTRATION FIELD MUST BE AT LEAST    #### '
    write(*,*) ' #### TWICE THE SYNCHRONISATION INTERVAL      #### '
    error stop
  endif

  if (mod(loutstep,lsynctime).ne.0) then
    write(*,*) ' #### FLEXPART MODEL ERROR! INTERVAL BETWEEN  #### '
    write(*,*) ' #### CONCENTRATION FIELDS MUST BE A MULTIPLE #### '
    write(*,*) ' #### OF THE SYNCHRONISATION INTERVAL         #### '
    error stop
  endif

  if ((loutstep/lsynctime).lt.2) then
    write(*,*) ' #### FLEXPART MODEL ERROR! INTERVAL BETWEEN  #### '
    write(*,*) ' #### CONCENTRATION FIELDS MUST BE AT LEAST   #### '
    write(*,*) ' #### TWICE THE SYNCHRONISATION INTERVAL      #### '
    error stop
  endif

  if (mod(loutsample,lsynctime).ne.0) then
    write(*,*) ' #### FLEXPART MODEL ERROR! SAMPLING TIME OF  #### '
    write(*,*) ' #### CONCENTRATION FIELD MUST BE A MULTIPLE  #### '
    write(*,*) ' #### OF THE SYNCHRONISATION INTERVAL         #### '
    error stop
  endif

  if ((mquasilag.eq.1).and.(iout.ge.4)) then
    write(*,*) ' #### FLEXPART MODEL ERROR! CONFLICTING       #### '
    write(*,*) ' #### OPTIONS: IF MQUASILAG=1, PLUME          #### '
    write(*,*) ' #### TRAJECTORY OUTPUT IS IMPOSSIBLE.        #### '
    error stop
  endif

  ! Compute modeling time in seconds and beginning date in Julian date
  !*******************************************************************

  outstep=real(abs(loutstep))
  if (ldirect.eq.1) then
    bdate=juldate(ibdate,ibtime)
    edate=juldate(iedate,ietime)
    ideltas=nint((edate-bdate)*86400.)
  else if (ldirect.eq.-1) then
    loutaver=-1*loutaver
    loutstep=-1*loutstep
    loutsample=-1*loutsample
    lsynctime=-1*lsynctime
    bdate=juldate(iedate,ietime)
    edate=juldate(ibdate,ibtime)
    ideltas=nint((edate-bdate)*86400.)
  else
    write(*,*) ' #### FLEXPART MODEL ERROR! DIRECTION IN      #### '
    write(*,*) ' #### FILE "COMMAND" MUST BE EITHER -1 OR 1.  #### '
    error stop
  endif

  return

999   write(*,*) ' #### FLEXPART MODEL ERROR! FILE "COMMAND"    #### '
  write(*,*) ' #### CANNOT BE OPENED IN THE DIRECTORY       #### '
  write(*,'(a)') path(1)(1:length(1))
  error stop

1000   write(*,*) ' #### FLEXPART MODEL ERROR! FILE "COMMAND"    #### '
  write(*,*) ' #### CANNOT BE OPENED IN THE DIRECTORY       #### '
  write(*,'(a)') path(2)(1:length(2))
  error stop
end subroutine readcommand

subroutine readdepo

  !*****************************************************************************
  !                                                                            *
  !  Reads dry deposition parameters needed by the procedure of Wesely (1989). *
  !  Wesely (1989): Parameterization of surface resistances to gaseous         *
  !  dry deposition in regional-scale numerical models.                        *
  !  Atmos. Environ. 23, 1293-1304.                                            *
  !                                                                            *
  !                                                                            *
  !      AUTHOR: Andreas Stohl, 19 May 1995                                    *
  !                                                                            *
  !*****************************************************************************
  !                                                                            *
  ! Variables:                                                                 *
  !                                                                            *
  ! rcl(maxspec,5,9) [s/m] Lower canopy resistance                             *
  ! rgs(maxspec,5,9) [s/m] Ground resistance                                   *
  ! rlu(maxspec,5,9) [s/m] Leaf cuticular resistance                           *
  ! rm(maxspec) [s/m]      Mesophyll resistance, set in readreleases           *
  ! ri(maxspec) [s/m]      Stomatal resistance                                 *
  !                                                                            *
  ! Constants:                                                                 *
  !                                                                            *
  !*****************************************************************************

  implicit none

  ! FOR THIS SUBROUTINE, numclass=9 IS ASSUMED
  !*******************************************

  real :: rluh(5,numclass),rgssh(5,numclass),rgsoh(5,numclass)
  real :: rclsh(5,numclass),rcloh(5,numclass)
  integer :: i,j,ic
  logical :: ios, ios2
  character(12) :: file_sfcdepo


  ! Read deposition constants related with landuse and seasonal category
  !*********************************************************************
  file_sfcdepo='sfcdepo.t'
  inquire(file=path(1)(1:length(1))//trim(file_sfcdepo),exist=ios)
  if (.not. ios) then
    file_sfcdepo='surfdepo.t'
    inquire(file=path(1)(1:length(1))//trim(file_sfcdepo),exist=ios2)
    if (ios2) then
      write(*,*) 'WARNING: surfdepo.t should be renamed to sfcdepo.t'
      write(*,*) 'Continuing with surfdepo.t'
    endif
  endif

  open(unitwesely,file=path(1)(1:length(1))//trim(file_sfcdepo), &
       status='old',err=999)

  do i=1,16
    read(unitwesely,*)
  end do
  do i=1,5
    read(unitwesely,*)
    read(unitwesely,'(8x,13f8.0)') (ri(i,j),j=1,numclass)
    read(unitwesely,'(8x,13f8.0)') (rluh(i,j),j=1,numclass)
    read(unitwesely,'(8x,13f8.0)') (rac(i,j),j=1,numclass)
    read(unitwesely,'(8x,13f8.0)') (rgssh(i,j),j=1,numclass)
    read(unitwesely,'(8x,13f8.0)') (rgsoh(i,j),j=1,numclass)
    read(unitwesely,'(8x,13f8.0)') (rclsh(i,j),j=1,numclass)
    read(unitwesely,'(8x,13f8.0)') (rcloh(i,j),j=1,numclass)
  end do

  ! TEST
  ! do 31 i=1,5
  !   ri(i,13)=ri(i,5)
  !   rluh(i,13)=rluh(i,5)
  !   rac(i,13)=rac(i,5)
  !   rgssh(i,13)=rgssh(i,5)
  !   rgsoh(i,13)=rgsoh(i,5)
  !   rclsh(i,13)=rclsh(i,5)
  !   rcloh(i,13)=rcloh(i,5)
  !31             continue
  ! TEST
  ! Sabine Eckhardt, Dec 06, set resistances of 9999 to 'infinite' (1E25)
  do i=1,5
    do j=1,numclass
      if    (ri(i,j).eq.9999.)    ri(i,j)=1.E25
      if  (rluh(i,j).eq.9999.)  rluh(i,j)=1.E25
      if   (rac(i,j).eq.9999.)   rac(i,j)=1.E25
      if (rgssh(i,j).eq.9999.) rgssh(i,j)=1.E25
      if (rgsoh(i,j).eq.9999.) rgsoh(i,j)=1.E25
      if (rclsh(i,j).eq.9999.) rclsh(i,j)=1.E25
      if (rcloh(i,j).eq.9999.) rcloh(i,j)=1.E25
    end do
  end do



  do i=1,5
    do j=1,numclass
      ri(i,j)=max(ri(i,j),0.001)
      rluh(i,j)=max(rluh(i,j),0.001)
      rac(i,j)=max(rac(i,j),0.001)
      rgssh(i,j)=max(rgssh(i,j),0.001)
      rgsoh(i,j)=max(rgsoh(i,j),0.001)
      rclsh(i,j)=max(rclsh(i,j),0.001)
      rcloh(i,j)=max(rcloh(i,j),0.001)
    end do
  end do
  close(unitwesely)


  ! Compute additional parameters
  !******************************

  do ic=1,nspec
    if (reldiff(ic).gt.0.) then     ! gas is dry deposited
      do i=1,5
        do j=1,numclass
          rlu(ic,i,j)=rluh(i,j)/(1.e-5*henry(ic)+f0(ic))
          rgs(ic,i,j)=1./(henry(ic)/(10.e5*rgssh(i,j))+f0(ic)/ &
               rgsoh(i,j))
          rcl(ic,i,j)=1./(henry(ic)/(10.e5*rclsh(i,j))+f0(ic)/ &
               rcloh(i,j))
        end do
      end do
    endif
  end do


  return

999   write(*,*) '### FLEXPART ERROR! FILE              ###'
  write(*,*) '### sfcdepo.t DOES NOT EXIST.        ###'
  error stop
end subroutine readdepo

subroutine readOHfield

  !*****************************************************************************
  !                                                                            *
  ! Reads the OH field into memory                                             *
  !                                                                            *
  ! AUTHOR: R.L. Thompson, Nov 2014                                            *
  !                                                                            *
  ! UPDATES:                                                                   *
  !   03/2018 SEC: Converted original netCDF files to binary format            *
  !*****************************************************************************
  !                                                                            *
  ! Variables:                                                                 *
  !                                                                            *
  ! path(numpath)              contains the path names                         *
  ! lonOH(nxOH)                longitude of OH fields                          *
  ! latOH(nyOH)                latitude of OH fields                           *
  ! altOH(nzOH)                altitude of OH fields                           *
  ! etaOH(nzOH)                eta-levels of OH fields                         *
  ! OH_field(nxOH,nyOH,nzOH,m) OH concentration (molecules/cm3)                *
  !                                                                            *
  !                                                                            *
  !*****************************************************************************

  use ohr_mod

  implicit none

  integer :: i,j,k,l,ierr,stat
  real, dimension(:), allocatable :: etaOH

  !  real, parameter :: gasct=8.314   ! gas constant
  !  real, parameter :: mct=0.02894   ! kg mol-1
  !  real, parameter :: g=9.80665     ! m s-2
  !  real, parameter :: lrate=0.0065  ! K m-1
  real, parameter :: scalehgt=7000. ! scale height in metres


  open(unitOH,file=trim(ohfields_path) &
       //'OH_FIELDS/OH_variables.bin',status='old', &
       form='UNFORMATTED', iostat=ierr, convert='little_endian')

  if(ierr.ne.0) then
    write(*,*) 'Cannot read binary OH fields in ', &
      trim(ohfields_path)//'OH_FIELDS/OH_variables.bin'
    error stop
  endif

  read(unitOH) nxOH
  read(unitOH) nyOH
  read(unitOH) nzOH
  write(*,*) nxOH,nyOH,nzOH

  ! allocate variables
  allocate(lonOH(nxOH),stat=stat)
  if (stat.ne.0) error stop "Could not allocate lonOH"
  allocate(latOH(nyOH),stat=stat)
  if (stat.ne.0) error stop "Could not allocate latOH"
  allocate(etaOH(nzOH),stat=stat)
  if (stat.ne.0) error stop "Could not allocate etaOH"
  allocate(altOH(nzOH),stat=stat)
  if (stat.ne.0) error stop "Could not allocate altOH"
  allocate(OH_field(nxOH,nyOH,nzOH,12),stat=stat)
  if (stat.ne.0) error stop "Could not allocate OH_field"
  allocate(OH_hourly(nxOH,nyOH,nzOH,2),stat=stat)
  if (stat.ne.0) error stop "Could not allocate OH_hourly"

  read(unitOH) (lonjr(i),i=1,360)
  read(unitOH) (latjr(i),i=1,180)
  read(unitOH) (((jrate_average(i,j,k),i=1,360),j=1,180),k=1,12)
  read(unitOH) (lonOH(i),i=1,nxOH)
  read(unitOH) (latOH(i),i=1,nyOH)
  read(unitOH) (lonOH(i),i=1,nxOH)

  read(unitOH) (altOH(i),i=1,nzOH)
  read(unitOH) ((((OH_field(i,j,k,l),i=1,nxOH),j=1,nyOH),k=1,nzOH),l=1,12)
  read(unitOH) ((((OH_hourly(i,j,k,l),i=1,nxOH),j=1,nyOH),k=1,nzOH),l=1,2)
end subroutine readOHfield

subroutine readlanduse

  !*****************************************************************************
  !                                                                            *
  !      Reads the landuse inventory into memory and relates it to Leaf Area   *
  !      Index and roughness length.                                           *
  !                                                                            *
  !      AUTHOR: Andreas Stohl, 10 January 1994                                *
  !                                                                            *
  !*****************************************************************************
  !                                                                            *
  ! Variables:                                                                 *
  ! i                       loop indices                                       *
  ! landinvent(1200,600,13) area fractions of 13 landuse categories            *
  ! LENGTH(numpath)         length of the path names                           *
  ! PATH(numpath)           contains the path names                            *
  ! unitland                unit connected with landuse inventory              *
  !                                                                            *
  ! -----                                                                      *
  ! Sabine Eckhardt, Dec 06 - new landuse inventary                            *
  ! after                                                                      *
  ! Belward, A.S., Estes, J.E., and Kline, K.D., 1999,                         *
  ! The IGBP-DIS 1-Km Land-Cover Data Set DISCover:                            *
  ! A Project Overview: Photogrammetric Engineering and Remote Sensing,        *
  ! v. 65, no. 9, p. 1013-1020                                                 *
  !                                                                            *
  ! LANDUSE CATEGORIES:                                                        *
  !                                                                            *
  ! 1   Urban land                                                             *
  ! 2   Agricultural land                                                      *
  ! 3   Range land                                                             *
  ! 4   Deciduous forest                                                       *
  ! 5   Coniferous forest                                                      *
  ! 6   Mixed forest including wetland                                         *
  ! 7   water, both salt and fresh                                             *
  ! 8   barren land mostly desert                                              *
  ! 9   nonforested wetland                                                    *
  ! 10  mixed agricultural and range land                                      *
  ! 11  rocky open areas with low growing shrubs                               *
  ! 12  ice                                                                    *
  ! 13  rainforest                                                             *
  !                                                                            *
  !*****************************************************************************

  use drydepo_mod
  
  implicit none

  integer :: ix,jy,i,k,lu_cat,lu_perc
  integer(kind=1) :: ilr
  integer(kind=1) :: ilr_buffer(2160000)
  integer :: il,irecread
  real :: rlr, r2lr
  logical :: ios,ios2
  character(12) :: file_sfcdata


  ! Read landuse inventory
  !***********************
  ! The landuse information is saved in a compressed format and written
  ! out by records of the length of 1 BYTE. Each grid cell consists of 3
  ! Bytes, which include 3 landuse categories (val 1-13 and 16 percentage
  ! categories) So one half byte is used to store the Landusecat the other
  ! for the percentageclass in 6.25 steps (100/6.25=16)
  ! e.g.
  ! 4 3  percentage 4 = 4*6.25 => 25% landuse class 3
  ! 2 1  percentage 2 = 2*6.25 => 13% landuse class 1
  ! 1 12 percentage 1 = 1*6.26 => 6.25% landuse class 12

  open(unitland,file=path(1)(1:length(1))//'IGBP_int1.dat',status='old', &
    form='UNFORMATTED', err=998, convert='little_endian')
  read (unitland) (ilr_buffer(i),i=1,2160000)
  close(unitland)

  irecread=1
  do ix=1,1200
    do jy=1,600
  ! the 3 most abundant landuse categories in the inventory
  ! first half byte contains the landuse class
  ! second half byte contains the respective percentage
      do k=1,3
  ! 1 byte is read
        ilr=ilr_buffer(irecread)
  !      ilr=0
        irecread=irecread+1
  ! as only signed integer values exist an unsigned value is constructed
        if (ilr.lt.0) then
           il=ilr+256
        else
           il=ilr
        endif
  ! dividing by 16 has the effect to get rid of the right half of the byte
  ! so just the left half remains, this corresponds to a shift right of 4
  ! bits
        rlr=real(il)/16.
        lu_cat=int(rlr)
  ! the left half of the byte is substracted from the whole in order to
  ! get only the right half of the byte
        r2lr=rlr-int(rlr)
  ! shift left by 4
        lu_perc=int(r2lr*16.)
        landinvent(ix,jy,k)=lu_cat
        landinvent(ix,jy,k+3)=lu_perc
  ! if ((jy.lt.10).and.(ix.lt.10)) write(*,*) 'reading: ',ix,jy,lu_cat,lu_perc
      end do
    end do
  end do

  ! Read relation landuse,z0
  !*****************************
  file_sfcdata='sfcdata.t'
  inquire(file=path(1)(1:length(1))//trim(file_sfcdata),exist=ios)
  if (.not. ios) then
    file_sfcdata='surfdata.t'
    inquire(file=path(1)(1:length(1))//trim(file_sfcdata),exist=ios2)
    if (ios2) then
      write(*,*) 'WARNING: surfdata.t should be renamed to sfcdata.t'
      write(*,*) 'Continuing with surfdata.t'
    endif
  endif

  open(unitsfcdata,file=path(1)(1:length(1))//trim(file_sfcdata), &
            status='old',err=999)

  do i=1,4
    read(unitsfcdata,*)
  end do
  do i=1,numclass
    read(unitsfcdata,'(45x,f15.3)') z0(i)
  end do
  close(unitsfcdata)

  return

  ! Issue error messages
  !*********************
998 write(*,*) ' #### FLEXPART ERROR! FILE                     ####'
  write(*,*)   ' #### ', path(1)(1:length(1))//'IGBP_int1.dat'
  write(*,*)   " #### (LANDUSE INVENTORY) COULD NOT BE OPENED  ####"
  error stop

999 write(*,*) ' #### FLEXPART ERROR! FILE              ####'
  write(*,*)   ' #### ', path(1)(1:length(1))//'sfcdata.t'
  write(*,*)   ' #### DOES NOT EXIST.                   ####'
  error stop
end subroutine readlanduse

subroutine readoutgrid

  !*****************************************************************************
  !                                                                            *
  !     This routine reads the user specifications for the output grid.        *
  !                                                                            *
  !     Author: A. Stohl                                                       *
  !                                                                            *
  !     4 June 1996                                                            *
  !     HSO, 1 July 2014
  !     Added optional namelist input
  !                                                                            *
  !*****************************************************************************
  !                                                                            *
  ! Variables:                                                                 *
  ! dxout,dyout          grid distance                                         *
  ! numxgrid,numygrid,numzgrid    grid dimensions                              *
  ! outlon0,outlat0      lower left corner of grid                             *
  ! outheight(maxzgrid)  height levels of output grid [m]                      *
  !                                                                            *
  ! Constants:                                                                 *
  ! unitoutgrid          unit connected to file OUTGRID                        *
  !                                                                            *
  !*****************************************************************************

  use outgrid_mod

  implicit none

  integer :: i,j,stat
  real :: outhelp,xr,xr1,yr,yr1
  real,parameter :: eps=1.e-4

  ! namelist variables
  integer, parameter :: maxoutlev=500
  integer :: ios
  real,allocatable, dimension (:) :: outheights

  ! declare namelist
  namelist /outgrid/ &
    outlon0,outlat0, &
    numxgrid,numygrid, &
    dxout,dyout, &
    outheights

  ! allocate large array for reading input
  allocate(outheights(maxoutlev),stat=stat)
  if (stat.ne.0) write(*,*)'ERROR: could not allocate outheights'

  ! helps identifying failed namelist input
  dxout=-1.0
  outheights=-1.0

  ! Open the OUTGRID file and read output grid specifications
  !**********************************************************

  open(unitoutgrid,file=path(1)(1:length(1))//'OUTGRID',status='old', &
    form='formatted',err=999)

  ! try namelist input
  read(unitoutgrid,outgrid,iostat=ios)
  close(unitoutgrid)

  if ((dxout.le.0).or.(ios.ne.0)) then

    ios=1

    open(unitoutgrid,file=path(1)(1:length(1))//'OUTGRID',status='old',err=999)

    call skplin(5,unitoutgrid)

    ! 1.  Read horizontal grid specifications
    !****************************************

    call skplin(3,unitoutgrid)
    read(unitoutgrid,'(4x,f11.4)') outlon0
    call skplin(3,unitoutgrid)
    read(unitoutgrid,'(4x,f11.4)') outlat0
    call skplin(3,unitoutgrid)
    read(unitoutgrid,'(4x,i5)') numxgrid
    call skplin(3,unitoutgrid)
    read(unitoutgrid,'(4x,i5)') numygrid
    call skplin(3,unitoutgrid)
    read(unitoutgrid,'(4x,f12.5)') dxout
    call skplin(3,unitoutgrid)
    read(unitoutgrid,'(4x,f12.5)') dyout

  endif

  ! Check validity of output grid (shall be within model domain)
  !*************************************************************

  xr=outlon0+real(numxgrid)*dxout
  yr=outlat0+real(numygrid)*dyout
  xr1=xlon0+real(nxmin1)*dx
  yr1=ylat0+real(nymin1)*dy
  if ((outlon0+eps.lt.xlon0).or.(outlat0+eps.lt.ylat0) &
       .or.(xr.gt.xr1+eps).or.(yr.gt.yr1+eps)) then
    write(*,*) outlon0,outlat0
    write(*,*) xr1,yr1,xlon0,ylat0,xr,yr,dxout,dyout
    write(*,*) ' #### FLEXPART MODEL ERROR! PART OF OUTPUT    ####'
    write(*,*) ' #### GRID IS OUTSIDE MODEL DOMAIN. CHANGE    ####'
    write(*,*) ' #### FILE OUTGRID IN DIRECTORY               ####'
    write(*,'(a)') path(1)(1:length(1))
    error stop
  endif

  ! 2. Count Vertical levels of output grid
  !****************************************

  if (ios.ne.0) then
    j=0
100 j=j+1
    do i=1,3
      read(unitoutgrid,*,end=99)
    end do
    read(unitoutgrid,'(4x,f7.1)',end=99) outhelp
    if (outhelp.eq.0.) goto 99
    goto 100
99  numzgrid=j-1
  else
    do i=1,maxoutlev
      if (outheights(i).lt.0) exit
    end do
    numzgrid=i-1
  end if

  allocate(outheight(numzgrid),stat=stat)
  if (stat.ne.0) write(*,*)'ERROR: could not allocate outheight'
  allocate(outheighthalf(numzgrid),stat=stat)
  if (stat.ne.0) write(*,*)'ERROR: could not allocate outheighthalf'

  ! 2. Vertical levels of output grid
  !**********************************

  if (ios.ne.0) then

    rewind(unitoutgrid)
    call skplin(29,unitoutgrid)

    do j=1,numzgrid
      do i=1,3
        read(unitoutgrid,*)
      end do
      read(unitoutgrid,'(4x,f7.1)') outhelp
      outheight(j)=outhelp
      outheights(j)=outhelp
    end do
    close(unitoutgrid)

  else

    do j=1,numzgrid
      outheight(j)=outheights(j)
    end do

  endif

  ! write outgrid file in namelist format to output directory if requested
  if (nmlout.and.lroot) then
    ! reallocate outheights with actually required dimension for namelist writing
    deallocate(outheights)
    allocate(outheights(numzgrid),stat=stat)
    if (stat.ne.0) write(*,*)'ERROR: could not allocate outheights'

    do j=1,numzgrid
      outheights(j)=outheight(j)
    end do

    open(unitoutgrid,file=path(2)(1:length(2))//'OUTGRID.namelist',err=1000)
    write(unitoutgrid,nml=outgrid)
    close(unitoutgrid)
  endif

  ! Check whether vertical levels are specified in ascending order
  !***************************************************************

  do j=2,numzgrid
    if (outheight(j).le.outheight(j-1)) then
    write(*,*) ' #### FLEXPART MODEL ERROR! YOUR SPECIFICATION#### '
    write(*,*) ' #### OF OUTPUT LEVELS IS CORRUPT AT LEVEL    #### '
    write(*,*) ' #### ',j,'                              #### '
    write(*,*) ' #### PLEASE MAKE CHANGES IN FILE OUTGRID.    #### '
    endif
  end do

  ! Determine the half levels, i.e. middle levels of the output grid
  !*****************************************************************

  outheighthalf(1)=outheight(1)/2.
  do j=2,numzgrid
    outheighthalf(j)=(outheight(j-1)+outheight(j))/2.
  end do

  xoutshift=xlon0-outlon0
  youtshift=ylat0-outlat0

  allocate(oroout(0:numxgrid-1,0:numygrid-1),stat=stat)
  if (stat.ne.0) write(*,*)'ERROR: could not allocate oroout'
  allocate(area(0:numxgrid-1,0:numygrid-1),stat=stat)
  if (stat.ne.0) write(*,*)'ERROR: could not allocate area'
  allocate(volume(0:numxgrid-1,0:numygrid-1,numzgrid),stat=stat)
  if (stat.ne.0) write(*,*)'ERROR: could not allocate volume'
  allocate(areaeast(0:numxgrid-1,0:numygrid-1,numzgrid),stat=stat)
  if (stat.ne.0) write(*,*)'ERROR: could not allocate areaeast'
  allocate(areanorth(0:numxgrid-1,0:numygrid-1,numzgrid),stat=stat)
  if (stat.ne.0) write(*,*)'ERROR: could not allocate areanorth'
  return

999 write(*,*) ' #### FLEXPART MODEL ERROR! FILE "OUTGRID"    #### '
  write(*,*) ' #### CANNOT BE OPENED IN THE DIRECTORY       #### '
  write(*,'(a)') path(1)(1:length(1))
  error stop

1000 write(*,*) ' #### FLEXPART MODEL ERROR! FILE "OUTGRID"    #### '
  write(*,*) ' #### CANNOT BE OPENED IN THE DIRECTORY       #### '
  write(*,'(a)') path(2)(1:length(2))
  error stop
end subroutine readoutgrid

subroutine readoutgrid_nest

  !*****************************************************************************
  !                                                                            *
  !     This routine reads the user specifications for the output nest.        *
  !                                                                            *
  !     Author: A. Stohl                                                       *
  !                                                                            *
  !     4 June 1996                                                            *
  !                                                                            *
  !*****************************************************************************
  !                                                                            *
  ! Variables:                                                                 *
  ! dxoutn,dyoutn        grid distances of output nest                         *
  ! numxgridn,numygridn,numzgrid    nest dimensions                            *
  ! outlon0n,outlat0n    lower left corner of nest                             *
  ! outheight(maxzgrid)  height levels of output grid [m]                      *
  !                                                                            *
  ! Constants:                                                                 *
  ! unitoutgrid          unit connected to file OUTGRID                        *
  !                                                                            *
  !*****************************************************************************

  use outgrid_mod

  implicit none

  integer :: stat
  real :: xr,xr1,yr,yr1
  real,parameter :: eps=1.e-4

  integer :: ios

  ! declare namelist
  namelist /outgridn/ &
    outlon0n,outlat0n, &
    numxgridn,numygridn, &
    dxoutn,dyoutn

  ! helps identifying failed namelist input
  dxoutn=-1.0

  ! Open the OUTGRID file and read output grid specifications
  !**********************************************************

  open(unitoutgrid,file=path(1)(1:length(1))//'OUTGRID_NEST',form='formatted',status='old',err=999)

  ! try namelist input
  read(unitoutgrid,outgridn,iostat=ios)
  close(unitoutgrid)

  if ((dxoutn.le.0).or.(ios.ne.0)) then

    open(unitoutgrid,file=path(1)(1:length(1))//'OUTGRID_NEST',status='old',err=999)
    call skplin(5,unitoutgrid)

    ! 1.  Read horizontal grid specifications
    !****************************************

    call skplin(3,unitoutgrid)
    read(unitoutgrid,'(4x,f11.4)') outlon0n
    call skplin(3,unitoutgrid)
    read(unitoutgrid,'(4x,f11.4)') outlat0n
    call skplin(3,unitoutgrid)
    read(unitoutgrid,'(4x,i5)') numxgridn
    call skplin(3,unitoutgrid)
    read(unitoutgrid,'(4x,i5)') numygridn
    call skplin(3,unitoutgrid)
    read(unitoutgrid,'(4x,f12.5)') dxoutn
    call skplin(3,unitoutgrid)
    read(unitoutgrid,'(4x,f12.5)') dyoutn

    close(unitoutgrid)
  endif

  ! write outgrid_nest file in namelist format to output directory if requested
  if (nmlout.and.lroot) then
    open(unitoutgrid,file=path(2)(1:length(2))//'OUTGRID_NEST.namelist',err=1000)
    write(unitoutgrid,nml=outgridn)
    close(unitoutgrid)
  endif

  allocate(orooutn(0:numxgridn-1,0:numygridn-1),stat=stat)
  if (stat.ne.0) write(*,*)'ERROR: could not allocate orooutn'
  allocate(arean(0:numxgridn-1,0:numygridn-1),stat=stat)
  if (stat.ne.0) write(*,*)'ERROR: could not allocate arean'
  allocate(volumen(0:numxgridn-1,0:numygridn-1,numzgrid),stat=stat)
  if (stat.ne.0) write(*,*)'ERROR: could not allocate volumen'

  ! Check validity of output grid (shall be within model domain)
  !*************************************************************

  xr=outlon0n+real(numxgridn)*dxoutn
  yr=outlat0n+real(numygridn)*dyoutn
  xr1=xlon0+real(nxmin1)*dx
  yr1=ylat0+real(nymin1)*dy
  if ((outlon0n+eps.lt.xlon0).or.(outlat0n+eps.lt.ylat0) &
       .or.(xr.gt.xr1+eps).or.(yr.gt.yr1+eps)) then
    write(*,*) ' #### FLEXPART MODEL ERROR! PART OF OUTPUT    ####'
    write(*,*) ' #### NEST IS OUTSIDE MODEL DOMAIN. CHANGE    ####'
    write(*,*) ' #### FILE OUTGRID IN DIRECTORY               ####'
    write(*,'(a)') path(1)(1:length(1))
    error stop
  endif

  xoutshiftn=xlon0-outlon0n
  youtshiftn=ylat0-outlat0n
  return

999 write(*,*) ' #### FLEXPART MODEL ERROR! FILE "OUTGRID"    #### '
  write(*,*) ' #### CANNOT BE OPENED IN THE DIRECTORY       #### '
  write(*,'(a)') path(1)(1:length(1))
  error stop

1000 write(*,*) ' #### FLEXPART MODEL ERROR! FILE "OUTGRID"    #### '
  write(*,*) ' #### CANNOT BE OPENED IN THE DIRECTORY       #### '
  write(*,'(a)') path(2)(1:length(2))
  error stop
end subroutine readoutgrid_nest

subroutine readpaths

  !*****************************************************************************
  !                                                                            *
  !     Reads the pathnames, where input/output files are expected to be.      *
  !     The file pathnames must be available in the current working directory. *
  !                                                                            *
  !     Author: A. Stohl                                                       *
  !                                                                            *
  !     1 February 1994                                                        *
  !     last modified                                                          *
  !     HS, 7.9.2012                                                           *
  !     option to give pathnames file as command line option                   *
  !                                                                            *
  !*****************************************************************************
  !                                                                            *
  ! Variables:                                                                 *
  ! length(numpath)    lengths of the path names                               *
  ! path(numpath)      pathnames of input/output files                         *
  !                                                                            *
  ! Constants:                                                                 *
  ! numpath            number of pathnames to be read in                       *
  !                                                                            *
  !*****************************************************************************

  implicit none

  integer   :: i
  character(256) :: string_test 
  character(1) :: character_test 

  ! Read the pathname information stored in unitpath
  !*************************************************

  open(unitpath,file=trim(pathfile),status='old',err=999)

  do i=1,numpath
    read(unitpath,'(a)',err=998) path(i)
    length(i)=index(path(i),' ')-1

    
    string_test = path(i)
    character_test = string_test(length(i):length(i))
    !print*, 'character_test,  string_test ', character_test,  string_test 
      if ((character_test .NE. '/') .AND. (i .LT. 4))  then
         print*, 'WARNING: path not ending in /' 
         print*, path(i)
         path(i) = string_test(1:length(i)) // '/'
         length(i)=length(i)+1
         print*, 'fix: padded with /' 
         print*, path(i)
         print*, 'length(i) increased 1' 
      endif
  end do

  ! Check whether any nested subdomains are to be used
  !***************************************************

  do i=1,maxnests
  ! ESO 2016 Added 'end'/'err' in case user forgot '====' at end of file and
  ! maxnests > numbnests
    read(unitpath,'(a)', end=30, err=30) path(numpath+2*(i-1)+1)
    read(unitpath,'(a)', end=30, err=30) path(numpath+2*(i-1)+2)
    if (path(numpath+2*(i-1)+1)(1:5).eq.'=====') goto 30
    length(numpath+2*(i-1)+1)=index(path(numpath+2*(i-1)+1),' ')-1
    length(numpath+2*(i-1)+2)=index(path(numpath+2*(i-1)+2),' ')-1
  end do


  ! Determine number of available nested domains
  !*********************************************

30  numbnests=i-1

  close(unitpath)
  return

998   write(*,*) ' #### TRAJECTORY MODEL ERROR! ERROR WHILE     #### '
  write(*,*) ' #### READING FILE PATHNAMES.                 #### '
  error stop

999   write(*,*) ' #### TRAJECTORY MODEL ERROR! FILE "pathnames"#### '
  write(*,*) ' #### CANNOT BE OPENED IN THE CURRENT WORKING #### '
  write(*,*) ' #### DIRECTORY.                              #### '
  error stop
end subroutine readpaths

subroutine readreceptors

  !*****************************************************************************
  !                                                                            *
  !     This routine reads the user specifications for the receptor points.    *
  !                                                                            *
  !     Author: A. Stohl                                                       *
  !     1 August 1996                                                          *
  !                                                                            *
  !     HSO, 14 August 2013: Added optional namelist input
  !     PS, 2/2015: access= -> position=
  !     PS, 6/2015: variable names, simplify code
  !     PS, 3/2023: remove position=append, makes no sense for new file        *
  !                                                                            *
  !*****************************************************************************
  !                                                                            *
  ! Variables:                                                                 *
  ! receptorarea   area of dx*dy at location of receptor                       *
  ! receptorname   names of receptors                                          *
  ! xreceptor,yreceptor  coordinates of receptor points                        *
  !                                                                            *
  ! Constants:                                                                 *
  ! unitreceptor         unit connected to file RECEPTORS                      *
  !                                                                            *
  !*****************************************************************************

  implicit none

  integer :: j
  real :: xm,ym
  character(len=16) :: receptor

  integer :: ios,stat
  real :: lon,lat   ! for namelist input, lon/lat are used instead of x,y
  real,allocatable,dimension(:) :: tmpxrec,tmpyrec,tmprecarea
  character(len=16),allocatable,dimension(:) :: tmprecname
  ! declare namelist
  namelist /receptors/ receptor, lon, lat

  numreceptor=0 ! Initialise numreceptor
!CPS I comment this out - why should we not have receptor output in bwd runs? 
  ! For backward runs, do not allow receptor output. Thus, set number of
  ! receptors to zero
  !  if (ldirect.lt.0) then
  !    numreceptor=0
  !    return
  !  endif


  ! Open the RECEPTORS file and read output grid specifications
  !************************************************************

  open (unitreceptor,file=trim(path(1))//'RECEPTORS',form='formatted', &
    status='old',err=999)

  lon = -999.
  lat = -999.

  ! try namelist input
  read(unitreceptor,receptors,iostat=ios)
  close (unitreceptor)

  ! prepare namelist output if requested
  if (nmlout) open(unitreceptorout,file=trim(path(2))// &
    'RECEPTORS.namelist',err=1000)


  if (ios .ne. 0) then ! read as regular text file

    goto 991 ! wrong variable name in namelist

  elseif (ios .eq. 0) then ! read as namelist

    if (nint(lon) .eq. -999 .or. nint(lat) .eq. -999) goto 993

    ! PS: reopen file otherwise first receptor is skipped!
    open (unitreceptor,file=trim(path(1))//'RECEPTORS',status='old',err=999)

    ! Read the names and coordinates of the receptors
    !************************************************
    j=0
    do while (ios .eq. 0)
      j=j+1
      read(unitreceptor,receptors,iostat=ios)
      if (ios .eq. 0) then
        numreceptor = j
        allocate( tmprecname(j),tmpxrec(j),tmpyrec(j),tmprecarea(j),stat=stat)
        if (stat.ne.0) write(*,*)'ERROR: could not allocate tmp arrays in readrec'
        if (j.gt.1) then
          tmprecname(1:j-1)=receptorname
          tmpxrec(1:j-1)=xreceptor
          tmpyrec(1:j-1)=yreceptor
          tmprecarea(1:j-1)=receptorarea
        endif
        tmprecname(j)=receptor
        tmpxrec(j)=(lon-xlon0)/dx       ! transform to grid coordinates
        tmpyrec(j)=(lat-ylat0)/dy
        xm=r_earth*cos(lat*pi/180.)*dx/180.*pi
        ym=r_earth*dy/180.*pi
        tmprecarea(j)=xm*ym

        call move_alloc(tmprecname,receptorname)
        call move_alloc(tmpxrec,xreceptor)
        call move_alloc(tmpyrec,yreceptor)
        call move_alloc(tmprecarea,receptorarea)

      ! write receptors in namelist format to output directory if requested
        if (nmlout) write(unitreceptorout,nml=receptors)
!        if (nmlout) write(unitreceptorout,nml=nml_receptors)
      elseif (ios .gt. 0) then
        write(*,*) ' ### FLEXPART MODEL ERROR! Error in RECEPTORS namelist ###'
        stop 'Error in RECEPTORS namelist'
      ! else
      !   write (*,*) 'receptor read in nml format, ios<0', ios
      !   write (*,receptors)
      endif
    end do ! end nml receptors reading loop

  else ! ios<0 = EOF, read as conventional input file

    open (unitreceptor,file=trim(path(1))//'RECEPTORS',status='old',err=999)
    call skplin(5,unitreceptor)

    ! Read the names and coordinates of the receptors
    !************************************************

    j=1
100 continue
      read(unitreceptor,*,end=99)
      read(unitreceptor,*,end=99)
      read(unitreceptor,*,end=99)
      read(unitreceptor,'(4x,a16)',end=99) receptor
      call skplin(3,unitreceptor)
      read(unitreceptor,'(4x,f11.4)',end=99) lon
      call skplin(3,unitreceptor)
      read(unitreceptor,'(4x,f11.4)',end=99) lat
      if (lon.eq.0. .and. lat.eq.0. .and. &
         (receptor .eq. '                ')) then
        write(*,*) 'WARNING: looks like empty receptor at south pole;'// &
          ' will be skipped'
        j=j-1
        goto 100
      endif

      numreceptor = j
      allocate( tmprecname(j),tmpxrec(j),tmpyrec(j),tmprecarea(j),stat=stat)
      if (stat.ne.0) write(*,*)'ERROR: could not allocate tmp arrays in readrec'
      if (j.gt.1) then
        tmprecname(1:j-1)=receptorname
        tmpxrec(1:j-1)=xreceptor
        tmpyrec(1:j-1)=yreceptor
        tmprecarea(1:j-1)=receptorarea
      endif
      tmprecname(j)=receptor
      tmpxrec(j)=(lon-xlon0)/dx       ! transform to grid coordinates
      tmpyrec(j)=(lat-ylat0)/dy
      xm=r_earth*cos(lat*pi/180.)*dx/180.*pi
      ym=r_earth*dy/180.*pi
      tmprecarea(j)=xm*ym

      call move_alloc(tmprecname,receptorname)
      call move_alloc(tmpxrec,xreceptor)
      call move_alloc(tmpyrec,yreceptor)
      call move_alloc(tmprecarea,receptorarea)
      ! write receptors file in namelist format to output directory if requested
      if (nmlout) write(unitreceptorout,nml=receptors)
!      if (nmlout) write(unitreceptorout,nml=nml_receptors)
    goto 100

99  continue

  endif ! end no-nml / nml bloc

  close (unitreceptor)
  if (nmlout) close (unitreceptorout)

  write(*,*) 'Number of receptors: ',numreceptor
  return

991 continue
  write(*,*) '#### FLEXPART ERROR: wrong variable names present'
  write(*,*) '#### in namelist in file RECEPTORS'
  write(*,*) '#### note that in v11+ coordinate names are lon and lat'

  error stop

993 continue
  write(*,*) '#### FLEXPART ERROR: namelist in file RECEPTORS'
  write(*,*) '#### first receptor point did not contain lon and/or lat'
  write(*,*) '#### Check your namelist!'
  error stop

999 write(*,*) 'INFORMATION: input file RECEPTORS cannot be opened'
    write(*,*) 'in directory '//trim(path(1))
    write(*,*) 'Continuing without RECEPTOR'

  numreceptor=0
  return

1000 write(*,*) ' #### FLEXPART MODEL ERROR! File "RECEPTORS"      #### '
  write(*,*)    ' #### cannot be opened in the output directory    #### '
  write(*,'(a)') ' #### '//trim(path(2))
  write(*,*)    ' #### either write perm missing or old file exists ###'

  error stop

end subroutine readreceptors


subroutine readreleases

  !*****************************************************************************
  !                                                                            *
  !     This routine reads the release point specifications for the current    *
  !     model run. Several release points can be used at the same time.        *
  !                                                                            *
  !     Author: A. Stohl                                                       *
  !                                                                            *
  !     18 May 1996                                                            *
  !                                                                            *
  !     Update: 29 January 2001                                                *
  !     Release altitude can be either in magl or masl                         *
  !     HSO, 12 August 2013
  !     Added optional namelist input
  !                                                                            *
  !*****************************************************************************
  !                                                                            *
  ! Variables:                                                                 *
  ! decay               decay constant of species                              *
  ! dquer [um]          mean particle diameters                                *
  ! dsigma              e.g. dsigma=10 or dsigma=0.1 means that 68% of the mass*
  !                     are between 0.1*dquer and 10*dquer                     *
  ! ireleasestart, ireleaseend [s] starting time and ending time of each       *
  !                     release                                                *
  ! kindz               1: zpoint is in m agl, 2: zpoint is in m asl, 3: zpoint*
  !                     is in hPa                                              *
  ! npart               number of particles to be released                     *
  ! nspec               number of species to be released                       *
  ! density [kg/m3]     density of the particles                               *
  ! rm [s/m]            Mesophyll resistance                                   *
  ! species             name of species                                        *
  ! xmass               total mass of each species                             *
  ! xpoint1,ypoint1     geograf. coordinates of lower left corner of release   *
  !                     area                                                   *
  ! xpoint2,ypoint2     geograf. coordinates of upper right corner of release  *
  !                     area                                                   *
  ! weta_gas, wetb_gas  parameters for below-cloud scavenging (gas)            *
  ! crain_aero, csnow_aero  parameters for below-cloud scavenging (aerosol)    *
  ! ccn_aero, in_aero   parameters for in-cloud scavenging (aerosol)           *
  ! zpoint1,zpoint2     height range, over which release takes place           *
  ! num_min_discrete    if less, release cannot be randomized and happens at   *
  !                     time mid-point of release interval                     *
  ! lroot               true if serial version, or if MPI and root process     *
  !                                                                            *
  !*****************************************************************************

  use point_mod
  use xmass_mod
  use drydepo_mod

  implicit none

  integer :: numpartmax,i,j,id1,it1,id2,it2,stat,irel,ispc,nsettle
  integer,parameter :: num_min_discrete=100
  real :: releaserate,cun
  real(kind=dp) :: jul1,jul2,julm
  real,parameter :: eps2=1.e-9
  character(len=50) :: line

  ! help variables for namelist reading
  integer :: numpoints, parts, ios, nspec_init
  integer*2 :: zkind
  integer :: idate1, itime1, idate2, itime2
  real :: lon1,lon2,lat1,lat2,z1,z2
  character*40 :: comment
  integer,parameter :: unitreleasesout=2
  real,allocatable, dimension (:) :: mass
  integer,allocatable, dimension (:) :: specnum_rel,specnum_rel2
  real,allocatable,dimension(:) :: vsh,fracth,schmih

  ! declare namelists
  namelist /releases_ctrl/ &
       nspec, &
       specnum_rel

  namelist /release/ &
       idate1, itime1, &
       idate2, itime2, &
       lon1, lon2, &
       lat1, lat2, &
       z1, z2, &
       zkind, &
       mass, &
       parts, &
       comment

  numpoint=0
  nspec_init=50 ! necessary to allocate specnum_rel, would be cleaner to change RELEASES 
                ! in one extra namelist

  allocate(specnum_rel(nspec_init),stat=stat)
  if (stat.ne.0) write(*,*)'ERROR: could not allocate specnum_rel'

  ! presetting namelist releases_ctrl
  nspec = -1  ! use negative value to determine failed namelist input
  specnum_rel = 0

  !sec, read release to find how many releasepoints should be allocated
  open(unitreleases,file=path(1)(1:length(1))//'RELEASES',status='old', &
    form='formatted',err=999)

  ! check if namelist input provided
  read(unitreleases,releases_ctrl,iostat=ios)

  if (nspec.gt.nspec_init) then
    write(*,*) 'Requested number of species:',nspec
    error stop 'More than 50 species at a time is not possible when using &
     & RELEASES. You can do this with the part_ic.nc option (IPIN=3).'
  end if
  ! Allocate fields with maxspec
  maxspec=nspec
  call alloc_com()
  ! allocate with maxspec for first input loop
  allocate(mass(maxspec),stat=stat)
  if (stat.ne.0) write(*,*)'ERROR: could not allocate mass'

  if (ios.ne.0) then
        backspace(unitreleases)
        read(unitreleases,fmt='(A)') line
        if (lroot) write(*,*) &
             'Invalid line in RELEASES: '//trim(line)
            !cgz; Check if the number of species in RELEASES_CTRL is larger than the maximum number of species in par_mod
        if ((lroot) .and. nspec.gt.maxspec) goto 994
  end if
    
  ! prepare namelist output if requested
  if (nmlout.and.lroot) then
    open(unitreleasesout,file=path(2)(1:length(2))//'RELEASES.namelist', &
      access='append',status='replace',err=1000)
  endif

  if ((ios.ne.0).or.(nspec.lt.0)) then
    if (lroot) write(*,*) 'RELEASE either having unrecognised entries, &
      &or in old format, please update to namelist format.'
      error stop
  else
    if ((ipin.ne.3).and.(ipin.ne.4)) then ! Not necessary to read releases when using part_ic.nc
      ios=0
      do while (ios.eq.0) 
        idate1=-1
        read(unitreleases,release,iostat=ios)
        if ((idate1.lt.0).or.(ios.ne.0)) then
          ios=1
        else
          numpoint=numpoint+1
        endif
      end do
      ios=0
    else
      numpoint=1
    endif
  endif

  rewind(unitreleases)

  if (nspec.gt.maxspec) goto 994

  ! allocate arrays of matching size for number of species (namelist output)
  deallocate(mass)
  allocate(mass(nspec),stat=stat)
  if (stat.ne.0) write(*,*)'ERROR: could not allocate mass'
  allocate(specnum_rel2(nspec),stat=stat)
  if (stat.ne.0) write(*,*)'ERROR: could not allocate specnum_rel2'
  specnum_rel2=specnum_rel(1:nspec)
  deallocate(specnum_rel) 
  ! eso: BUG, crashes here for nspec=12 and maxspec=6,
  ! TODO: catch error and exit
  allocate(specnum_rel(nspec),stat=stat)
  if (stat.ne.0) write(*,*)'ERROR: could not allocate specnum_rel'
  specnum_rel=specnum_rel2
  deallocate(specnum_rel2)

  !allocate memory for numpoint releaspoints
  allocate(ireleasestart(numpoint),stat=stat)
  if (stat.ne.0) write(*,*)'ERROR: could not allocate ireleasestart'
  allocate(ireleaseend(numpoint),stat=stat)
  if (stat.ne.0) write(*,*)'ERROR: could not allocate ireleaseend'
  allocate(xpoint1(numpoint),stat=stat)
  if (stat.ne.0) write(*,*)'ERROR: could not allocate xpoint1'
  allocate(xpoint2(numpoint),stat=stat)
  if (stat.ne.0) write(*,*)'ERROR: could not allocate xpoint2'
  allocate(ypoint1(numpoint),stat=stat)
  if (stat.ne.0) write(*,*)'ERROR: could not allocate ypoint1'
  allocate(ypoint2(numpoint),stat=stat)
  if (stat.ne.0) write(*,*)'ERROR: could not allocate ypoint2'
  allocate(zpoint1(numpoint),stat=stat)
  if (stat.ne.0) write(*,*)'ERROR: could not allocate zpoint1'
  allocate(zpoint2(numpoint),stat=stat)
  if (stat.ne.0) write(*,*)'ERROR: could not allocate zpoint2'
  allocate(kindz(numpoint),stat=stat)
  if (stat.ne.0) write(*,*)'ERROR: could not allocate kindz'
  allocate(xmass(numpoint,maxspec),stat=stat)
  if (stat.ne.0) write(*,*)'ERROR: could not allocate xmass'
  allocate(rho_rel(numpoint),stat=stat)
  if (stat.ne.0) write(*,*)'ERROR: could not allocate rho_rel'
  allocate(npart(numpoint),stat=stat)
  if (stat.ne.0) write(*,*)'ERROR: could not allocate npart'
  allocate(xmasssave(numpoint),stat=stat)
  if (stat.ne.0) write(*,*)'ERROR: could not allocate xmasssave'

  if (lroot) write (*,*) 'Releasepoints : ', numpoint

  do i=1,numpoint
    xmasssave(i)=0.
  end do

  !now save the information
  DEP=.false.
  DRYDEP=.false.
  WETDEP=.false.
  OHREA=.false.
  do i=1,maxspec
    DRYDEPSPEC(i)=.false.
    WETDEPSPEC(i)=.false.
  end do

  ! namelist output
  if (nmlout.and.lroot) then
    write(unitreleasesout,nml=releases_ctrl)
  endif

  do i=1,nspec
    call readspecies(specnum_rel(i),i)
  end do

  ! Allocate fields that depend on ndia
  call alloc_com_ndia

  do i=1,nspec

  ! Allocate temporary memory necessary for the different diameter bins
  !********************************************************************
    allocate( vsh(ndia(i)),fracth(ndia(i)),schmih(ndia(i)),stat=stat)
    if (stat.ne.0) error stop "Could not allocate vsh,fracth,schmih"

  ! Molecular weight
  !*****************

    if (((iout.eq.2).or.(iout.eq.3)).and.(weightmolar(i).lt.0.)) then
      write(*,*) 'For mixing ratio output, valid molar weight'
      write(*,*) 'must be specified for all simulated species.'
      write(*,*) 'Check table SPECIES or choose concentration'
      write(*,*) 'output instead if molar weight is not known.'
      error stop
    endif

  ! Radioactive decay
  !******************

    decay(i)=0.693147/decay(i) !conversion half life to decay constant


  ! Dry deposition of gases
  !************************

    if (reldiff(i).gt.0.) rm(i)=1./(henry(i)/3000.+100.*f0(i))    ! mesophyll resistance

  ! Dry deposition of particles
  !****************************

    vsetaver(i)=0.
    cunningham(i)=0.
    dquer(i)=dquer(i)*1000000.         ! Conversion m to um
    if (density(i).gt.0.) then         ! Additional parameters
      call part0(dquer(i),dsigma(i),density(i),ndia(i),fracth,schmih,cun,vsh)
      do j=1,ndia(i)
        fract(i,j)=fracth(j)
        schmi(i,j)=schmih(j)
        vset(i,j)=vsh(j)
        cunningham(i)=cunningham(i)+cun*fract(i,j)
        vsetaver(i)=vsetaver(i)-vset(i,j)*fract(i,j)
      end do
      if (lroot) write(*,*) 'Average settling velocity: ',i,vsetaver(i)
    endif

  ! Dry deposition for constant deposition velocity
  !************************************************

    dryvel(i)=dryvel(i)*0.01         ! conversion to m/s

  ! Check if wet deposition or OH reaction shall be calculated
  !***********************************************************

  ! ESO 04.2016 check for below-cloud scavenging (gas or aerosol)
    if ((dquer(i).le.0..and.(weta_gas(i).gt.0. .or. wetb_gas(i).gt.0.)) .or. &
         &(dquer(i).gt.0. .and. (crain_aero(i) .gt. 0. .or. csnow_aero(i).gt.0.)))  then
      WETDEP=.true.
      WETDEPSPEC(i)=.true.
      if (lroot) then
        write (*,*) '  Below-cloud scavenging: ON'
  !  write (*,*) 'Below-cloud scavenging coefficients: ',weta(i),i
      end if
    else
      if (lroot) write (*,*) '  Below-cloud scavenging: OFF'
    endif

  ! NIK 31.01.2013 + 10.12.2013 + 15.02.2015
    if (dquer(i).gt.0..and.(ccn_aero(i).gt.0. .or. in_aero(i).gt.0.))  then
      WETDEP=.true.
      WETDEPSPEC(i)=.true.
      if (lroot) then
        write (*,*) '  In-cloud scavenging: ON'
  !        write (*,*) 'In-cloud scavenging coefficients: ',&
  !           &ccn_aero(i),in_aero(i),i !,wetc_in(i), wetd_in(i),i
      end if
    else
      if (lroot) write (*,*) '  In-cloud scavenging: OFF' 
    endif

    if (ohcconst(i).gt.0.) then
      OHREA=.true.
      if (lroot) write (*,*) '  OHreaction switched on: ',ohcconst(i),i
    endif

    if ((reldiff(i).gt.0.).or.(density(i).gt.0.).or.(dryvel(i).gt.0.)) then
      DRYDEP=.true.
      DRYDEPSPEC(i)=.true.
    endif

    deallocate(vsh,fracth,schmih)
  end do ! end loop over species

  if (WETDEP.or.DRYDEP) DEP=.true.

  ! Not necessary to read releases when using part_ic.nc
  !*****************************************************
  if ((ipin.eq.3).or.(ipin.eq.4)) then
    maxpointspec_act=1
    return 
  endif
  
  ! Read specifications for each release point
  !*******************************************
  numpoints=numpoint
  numpoint=0
  numpartmax=0
  releaserate=0.
101 numpoint=numpoint+1

  if (numpoint.gt.numpoints) goto 250
  zkind = 1
  mass = 0
  parts = 0
  comment = ' '
  read(unitreleases,release,iostat=ios)
  id1=idate1
  it1=itime1
  id2=idate2
  it2=itime2
  xpoint1(numpoint)=lon1
  xpoint2(numpoint)=lon2
  ypoint1(numpoint)=lat1
  ypoint2(numpoint)=lat2
  zpoint1(numpoint)=z1
  zpoint2(numpoint)=z2
  kindz(numpoint)=zkind
  do i=1,nspec
    xmass(numpoint,i)=mass(i)
  end do
  npart(numpoint)=parts
  compoint(min(1001,numpoint))=comment

! namelist output
  if (nmlout.and.lroot) then
    write(unitreleasesout,nml=release)
  endif

  ! If a release point contains no particles, stop and issue error message
  !***********************************************************************

  if (npart(numpoint).eq.0) then
    write(*,*) 'FLEXPART MODEL ERROR'
    write(*,*) 'RELEASES file is corrupt.'
    write(*,*) 'At least for one release point, there are zero'
    write(*,*) 'particles released. Make changes to RELEASES.'
    error stop
  endif

  ! If FLEXPART is run for backward deposition force zpoint
  !*********************************************************************
  if (WETBKDEP) then
    zpoint1(numpoint)=0.
    zpoint2(numpoint)=20000.
    kindz(numpoint)=1
  endif
  if (DRYBKDEP) then
    zpoint1(numpoint)=0.
    zpoint2(numpoint)=2.*href
    kindz(numpoint)=1
  endif


  ! Check whether x coordinates of release point are within model domain
  !*********************************************************************

  if (xpoint1(numpoint).lt.xlon0) &
       xpoint1(numpoint)=xpoint1(numpoint)+360.
  if (xpoint1(numpoint).gt.xlon0+(nxmin1)*dx) &
       xpoint1(numpoint)=xpoint1(numpoint)-360.
  if (xpoint2(numpoint).lt.xlon0) &
       xpoint2(numpoint)=xpoint2(numpoint)+360.
  if (xpoint2(numpoint).gt.xlon0+(nxmin1)*dx) &
       xpoint2(numpoint)=xpoint2(numpoint)-360.

  ! Determine relative beginning and ending times of particle release
  !******************************************************************

  jul1=juldate(id1,it1)
  jul2=juldate(id2,it2)
  julm=(jul1+jul2)/2.
  if (jul1.gt.jul2) then
    write(*,*) 'FLEXPART MODEL ERROR'
    write(*,*) 'Release stops before it begins.'
    write(*,*) 'Make changes to file RELEASES.'
    error stop
  endif
  if (mdomainfill.eq.0) then   ! no domain filling
    if (ldirect.eq.1) then
      if (((jul1.lt.bdate).or.(jul2.gt.edate)).and.(ipin.eq.0)) then
        write(*,*) 'FLEXPART MODEL ERROR'
        write(*,*) 'Release starts before simulation begins or ends'
        write(*,*) 'after simulation stops.'
        write(*,*) 'Make files COMMAND and RELEASES consistent.'
        error stop
      endif
      if (npart(numpoint).gt.num_min_discrete) then
        ireleasestart(numpoint)=int((jul1-bdate)*86400.)
        ireleaseend(numpoint)=int((jul2-bdate)*86400.)
      else
        ireleasestart(numpoint)=int((julm-bdate)*86400.)
        ireleaseend(numpoint)=int((julm-bdate)*86400.)
      endif
    else if (ldirect.eq.-1) then
      if (((jul1.lt.edate).or.(jul2.gt.bdate)).and.(ipin.eq.0)) then
        write(*,*) 'FLEXPART MODEL ERROR'
        write(*,*) 'Release starts before simulation begins or ends'
        write(*,*) 'after simulation stops.'
        write(*,*) 'Make files COMMAND and RELEASES consistent.'
        error stop
      endif
      if (npart(numpoint).gt.num_min_discrete) then
        ireleasestart(numpoint)=int((jul1-bdate)*86400.)
        ireleaseend(numpoint)=int((jul2-bdate)*86400.)
      else
        ireleasestart(numpoint)=int((julm-bdate)*86400.)
        ireleaseend(numpoint)=int((julm-bdate)*86400.)
      endif
    endif
  endif


  ! Determine the release rate (particles per second) and total number
  ! of particles released during the simulation
  !*******************************************************************

  if (ireleasestart(numpoint).ne.ireleaseend(numpoint)) then
    releaserate=releaserate+real(npart(numpoint))/ &
         real(ireleaseend(numpoint)-ireleasestart(numpoint))
  else
    releaserate=99999999.
  endif
  numpartmax=numpartmax+npart(numpoint)
  goto 101

250 close(unitreleases)

  if (nmlout.and.lroot) then
    close(unitreleasesout)
  endif

  !if (lroot) write (*,*) 'Particles allocated (maxpart)  : ',maxpart
  if (lroot) write (*,*) 'Particles released (numpartmax): ',numpartmax
  numpoint=numpoint-1

  if (ioutputforeachrelease.eq.1) then
    maxpointspec_act=numpoint
  else
    maxpointspec_act=1
  endif

  ! Disable settling if more than 1 species at any release point
  ! or if MQUASILAG and more than one species
  !*************************************************************

  if (mquasilag.ne.0) then
    if (nspec.gt.1) lsettling=.false.
  else
    do irel=1,numpoint
      nsettle=0
      do ispc=1,nspec
        if (xmass(irel,ispc).gt.eps2) nsettle=nsettle+1
      end do
      if (nsettle.gt.1) lsettling=.false.
    end do
  end if

  if (lroot) then
    if (.not.lsettling) then 
      write(*,*) 'WARNING: more than 1 species per release point, settling &
           &disabled'
    end if
  end if

  ! Check, whether the total number of particles may exceed totally allowed
  ! number of particles at some time during the simulation
  !************************************************************************

  ! if (releaserate.gt. &
  !      0.99*real(maxpart)/real(lage(nageclass))) then
  !   if (numpartmax.gt.maxpart.and.lroot) then
  !     write(*,*) '#####################################################'
  !     write(*,*) '#### FLEXPART MODEL SUBROUTINE READRELEASES:     ####'
  !     write(*,*) '####                                             ####'
  !     write(*,*) '####WARNING - TOTAL NUMBER OF PARTICLES SPECIFIED####'
  !     write(*,*) '#### IN FILE "RELEASES" MAY AT SOME POINT DURING ####'
  !     write(*,*) '#### THE SIMULATION EXCEED THE MAXIMUM ALLOWED   ####'
  !     write(*,*) '#### NUMBER (MAXPART).IF RELEASES DO NOT OVERLAP,####'
  !     write(*,*) '#### FLEXPART CAN POSSIBLY COMPLETE SUCCESSFULLY.####'
  !     write(*,*) '#### HOWEVER, FLEXPART MAY HAVE TO STOP          ####'
  !     write(*,*) '#### AT SOME TIME DURING THE SIMULATION. PLEASE  ####'
  !     write(*,*) '#### MAKE SURE THAT YOUR SETTINGS ARE CORRECT.   ####'
  !     write(*,*) '#####################################################'
  !     write(*,*) 'Maximum release rate may be: ',releaserate, &
  !          ' particles per second'
  !     write(*,*) 'Maximum allowed release rate is: ', &
  !          real(maxpart)/real(lage(nageclass)),' particles per second'
  !     write(*,*) &
  !          'Total number of particles released during the simulation is: ', &
  !          numpartmax
  !     write(*,*) 'Maximum allowed number of particles is: ',maxpart
  !   endif
  ! endif


  if (lroot) then
    write(*,FMT='(A,ES14.7)') ' Total mass released:', sum(xmass(1:numpoint,1:nspec))
  end if

  return

994 write(*,*) '#####################################################'
  write(*,*) '#### FLEXPART MODEL SUBROUTINE READRELEASES:     ####'
  write(*,*) '####                                             ####'
  write(*,*) '#### ERROR - MAXIMUM NUMBER OF EMITTED SPECIES IS####'
  write(*,*) '#### TOO LARGE. PLEASE REDUCE NUMBER OF SPECIES. ####'
  write(*,*) '#####################################################'
  error stop

999 write(*,*) '#####################################################'
  write(*,*) '   FLEXPART MODEL SUBROUTINE READRELEASES: '
  write(*,*)
  write(*,*) 'FATAL ERROR - FILE CONTAINING PARTICLE RELEASE POINTS'
  write(*,*) 'POINTS IS NOT AVAILABLE OR YOU ARE NOT'
  write(*,*) 'PERMITTED FOR ANY ACCESS'
  write(*,*) '#####################################################'
  error stop

1000 write(*,*) ' #### FLEXPART MODEL ERROR! FILE "RELEASES"    #### '
  write(*,*) ' #### CANNOT BE OPENED IN THE DIRECTORY       #### '
  write(*,'(a)') path(2)(1:length(2))
  error stop
end subroutine readreleases

subroutine readspecies(id_spec,pos_spec)

  !*****************************************************************************
  !                                                                            *
  !     This routine reads names and physical constants of chemical species/   *
  !     radionuclides given in the parameter pos_spec                          *
  !                                                                            *
  !   Author: A. Stohl                                                         *
  !                                                                            *
  !   11 July 1996                                                             *
  !                                                                            *
  !   Changes:                                                                 *
  !   N. Kristiansen, 31.01.2013: Including parameters for in-cloud scavenging *
  !                                                                            *
  !   HSO, 13 August 2013
  !   added optional namelist input
  !                                                                            *
  !*****************************************************************************
  !                                                                            *
  ! Variables:                                                                 *
  ! decaytime(maxtable)   half time for radiological decay                     *
  ! specname(maxtable)    names of chemical species, radionuclides             *
  ! weta_gas, wetb_gas    Parameters for below-cloud scavenging of gasses      *
  ! crain_aero,csnow_aero Parameters for below-cloud scavenging of aerosols    *
  ! ccn_aero,in_aero      Parameters for in-cloud scavenging of aerosols       *
  ! ohcconst              OH reaction rate constant C                          *
  ! ohdconst              OH reaction rate constant D                          *
  ! ohnconst              OH reaction rate constant n                          *
  ! id_spec               SPECIES number as referenced in RELEASE file         *
  ! id_pos                position where SPECIES data shall be stored          *
  ! ni                    Number of diameter classes of particles              *                                                              *
  ! Constants:                                                                 *
  !                                                                            *
  !*****************************************************************************

  implicit none

  integer :: i, pos_spec,j
  integer :: idow,ihour,id_spec
  character(len=3) :: aspecnumb

  character(len=16) :: pspecies
  character(len=50) :: line
  real :: pdecay, pweta_gas, pwetb_gas, preldiff, phenry, pf0, pdensity, pdquer
  real :: pdsigma, pdryvel, pweightmolar, pohcconst, pohdconst, pohnconst
  real :: pcrain_aero, pcsnow_aero, pccn_aero, pin_aero
  real :: parea_dow(7), parea_hour(24), ppoint_dow(7), ppoint_hour(24)
  integer :: pndia
  integer :: ios
  integer :: pshape,porient
  ! Daria Tatsii: species shape properties
  real ::pla,pia,psa,f,e,paspectratio

  ! declare namelist
  namelist /species_params/ &
       pspecies, pdecay, pweta_gas, pwetb_gas, &
       pcrain_aero, pcsnow_aero, pccn_aero, pin_aero, &
       preldiff, phenry, pf0, pdensity, pdquer, &
       pdsigma, pndia, pdryvel, pweightmolar, pohcconst, pohdconst, pohnconst, &
       parea_dow, parea_hour, ppoint_dow, ppoint_hour, &
       pshape, paspectratio, pla, pia, psa, porient

  pspecies="" ! read failure indicator value
  pdecay=-999.9
  pweta_gas=-9.9E-09
  pwetb_gas=0.0
  pcrain_aero=-9.9E-09
  pcsnow_aero=-9.9E-09
  pccn_aero=-9.9E-09
  pin_aero=-9.9E-09
  preldiff=-9.9
  phenry=0.0
  pf0=0.0
  pdensity=-9.9E09
  pdquer=0.0
  pdsigma=0.0
  pndia=1
  pdryvel=-9.99
  pohcconst=-9.99
  pohdconst=-9.9E-09
  pohnconst=2.0
  pweightmolar=-999.9
  parea_dow=-999.9
  parea_hour=-999.9
  ppoint_dow=-999.9
  ppoint_hour=-999.9
  pshape=0 ! 0 for sphere, 1 for other shapes
  paspectratio=-1.
  pla=-1. ! longest axis in micrometer
  pia=-1. ! Intermediate axis
  psa=-1. ! Smallest axis
  porient=0 ! 0 for horizontal, 1 for random


  do j=1,24           ! initialize everything to no variation
    parea_hour(j)=1.
    ppoint_hour(j)=1.
    area_hour(pos_spec,j)=1.
    point_hour(pos_spec,j)=1.
  end do
  do j=1,7
    parea_dow(j)=1.
    ppoint_dow(j)=1.
    area_dow(pos_spec,j)=1.
    point_dow(pos_spec,j)=1.
  end do

  ! Open the SPECIES file and read species names and properties
  !************************************************************
  specnum(pos_spec)=id_spec
  write(aspecnumb,'(i3.3)') specnum(pos_spec)
  open(unitspecies,file=path(1)(1:length(1))//'SPECIES/SPECIES_'//aspecnumb, &
    status='old',form='formatted',err=998)
  write(*,*) 'reading SPECIES',specnum(pos_spec)

  ASSSPEC=.FALSE.

  ! try namelist input
  read(unitspecies,species_params,iostat=ios)
  !CGZ add check on which line of species file problem occurs
  if (ios.ne.0) then
    backspace(unitspecies)
    read(unitspecies,fmt='(A)') line
    if (lroot) write(*,*) &
        'Invalid line in species: '//trim(line)
  end if
  close(unitspecies)

  if ((len(trim(pspecies)).eq.0).or.(ios.ne.0)) then ! no namelist found
    if (lroot) write(*,*) "SPECIES file not in NAMELIST format, attempting to &
         &read as fixed format"

    ios=1

    open(unitspecies,file=path(1)(1:length(1))//'SPECIES/SPECIES_'//aspecnumb, &
      status='old',err=998)

    do i=1,6
      read(unitspecies,*)
    end do

    read(unitspecies,'(a10)',end=22) species(pos_spec)
  !  write(*,*) species(pos_spec)
    read(unitspecies,'(f18.1)',end=22) decay(pos_spec)
  !  write(*,*) decay(pos_spec)
    read(unitspecies,'(e18.1)',end=22) weta_gas(pos_spec)
  !  write(*,*) weta_gas(pos_spec)
    read(unitspecies,'(f18.2)',end=22) wetb_gas(pos_spec)
  !  write(*,*) wetb_gas(pos_spec)
    read(unitspecies,'(e18.1)',end=22) crain_aero(pos_spec)
  !  write(*,*) crain_aero(pos_spec)
    read(unitspecies,'(f18.2)',end=22) csnow_aero(pos_spec)
  !  write(*,*) csnow_aero(pos_spec)
  !*** NIK 31.01.2013: including in-cloud scavening parameters
    read(unitspecies,'(e18.1)',end=22) ccn_aero(pos_spec)
  !  write(*,*) ccn_aero(pos_spec)
    read(unitspecies,'(f18.2)',end=22) in_aero(pos_spec)
  !  write(*,*) in_aero(pos_spec)
    read(unitspecies,'(f18.1)',end=22) reldiff(pos_spec)
  !  write(*,*) reldiff(pos_spec)
    read(unitspecies,'(e18.1)',end=22) henry(pos_spec)
  !  write(*,*) henry(pos_spec)
    read(unitspecies,'(f18.1)',end=22) f0(pos_spec)
  !  write(*,*) f0(pos_spec)
    read(unitspecies,'(e18.1)',end=22) density(pos_spec)
  !  write(*,*) density(pos_spec)
    read(unitspecies,'(e18.1)',end=22) dquer(pos_spec)
  !  write(*,*) 'dquer(pos_spec):', dquer(pos_spec)
    read(unitspecies,'(e18.1)',end=22) dsigma(pos_spec)
  !  write(*,*) dsigma(pos_spec)
    read(unitspecies,'(i16)',end=22) ndia(pos_spec)
    ! write(*,*) ndia(pos_spec)
    read(unitspecies,'(f18.2)',end=22) dryvel(pos_spec)
  !  write(*,*) dryvel(pos_spec)
    read(unitspecies,'(f18.2)',end=22) weightmolar(pos_spec)
  !  write(*,*) weightmolar(pos_spec)
    read(unitspecies,'(e18.2)',end=22) ohcconst(pos_spec)
  !  write(*,*) ohcconst(pos_spec)
    read(unitspecies,'(f8.2)',end=22) ohdconst(pos_spec)
  !  write(*,*) ohdconst(pos_spec)
    read(unitspecies,'(f8.2)',end=22) ohnconst(pos_spec)
  !  write(*,*) ohnconst(pos_spec)

  ! Read in daily and day-of-week variation of emissions, if available
  !*******************************************************************

    read(unitspecies,*,end=22)
    do j=1,24     ! 24 hours, starting with 0-1 local time
      read(unitspecies,*) ihour,area_hour(pos_spec,j),point_hour(pos_spec,j)
    end do
    read(unitspecies,*)
    do j=1,7      ! 7 days of the week, starting with Monday
      read(unitspecies,*) idow,area_dow(pos_spec,j),point_dow(pos_spec,j)
    end do

    pspecies=species(pos_spec)
    pdecay=decay(pos_spec)
    pweta_gas=weta_gas(pos_spec)
    pwetb_gas=wetb_gas(pos_spec)
    pcrain_aero=crain_aero(pos_spec)
    pcsnow_aero=csnow_aero(pos_spec)
    pccn_aero=ccn_aero(pos_spec)
    pin_aero=in_aero(pos_spec)
    preldiff=reldiff(pos_spec)
    phenry=henry(pos_spec)
    pf0=f0(pos_spec)
    pdensity=density(pos_spec)
    pdquer=dquer(pos_spec)
    pdsigma=dsigma(pos_spec)
    pndia=ndia(pos_spec)
    pdryvel=dryvel(pos_spec)
    pweightmolar=weightmolar(pos_spec)
    pohcconst=ohcconst(pos_spec)
    pohdconst=ohdconst(pos_spec)
    pohnconst=ohnconst(pos_spec)


    do j=1,24     ! 24 hours, starting with 0-1 local time
      parea_hour(j)=area_hour(pos_spec,j)
      ppoint_hour(j)=point_hour(pos_spec,j)
    end do
    do j=1,7      ! 7 days of the week, starting with Monday
      parea_dow(j)=area_dow(pos_spec,j)
      ppoint_dow(j)=point_dow(pos_spec,j)
    end do

  else ! namelist available

    species(pos_spec)=pspecies
    decay(pos_spec)=pdecay
    weta_gas(pos_spec)=pweta_gas
    wetb_gas(pos_spec)=pwetb_gas
    crain_aero(pos_spec)=pcrain_aero
    csnow_aero(pos_spec)=pcsnow_aero
    ccn_aero(pos_spec)=pccn_aero
    in_aero(pos_spec)=pin_aero
    reldiff(pos_spec)=preldiff
    henry(pos_spec)=phenry
    f0(pos_spec)=pf0
    density(pos_spec)=pdensity
    dquer(pos_spec)=pdquer
    dsigma(pos_spec)=pdsigma
    ndia(pos_spec)=pndia
    dryvel(pos_spec)=pdryvel
    weightmolar(pos_spec)=pweightmolar
    ohcconst(pos_spec)=pohcconst
    ohdconst(pos_spec)=pohdconst
    ohnconst(pos_spec)=pohnconst
    ishape(pos_spec)=pshape
    orient(pos_spec)=porient


    ! Daria Tatsii 2023: compute particle shape dimensions
    if (ishape(pos_spec).ge.1) then ! Compute shape according to given axes
      select case (ishape(pos_spec))
        case (1)
          write(*,*) "Particle shape USER-DEFINED for particle", id_spec
          if ((psa.le.0.0).or.(pia.le.0.0).or.(pla.le.0.0)) then
            write(*,*) "#### ERROR: Shape=1 (user-defined) is chosen, &
              &but no valid axes are provided."
            write(*,*) "#### SPECIES file requires SA, IA, and LA parameter &
              &greater than zero."
            error stop
          endif
          write(*,*) "SA,IA,LA:",psa,pia,pla
        case (2) ! Cylinders (fibers) !
          if (paspectratio.le.0.0) then
            write(*,*) "#### ERROR: Shape=2 cylinder is chosen, but no valid apect ratio is provided."
            write(*,*) "#### SPECIES file requires ASPECTRATIO parameter greater than zero."
            error stop
          endif
          psa=(((dquer(pos_spec)**3.0)*2.0)/ &
                  (3.0*paspectratio))**(1.0/3.0)
          pia=psa
          pla=psa*paspectratio
          write(*,*) "Particle shape CYLINDER for particle", id_spec
          write(*,*) "SA,IA,LA:",psa,pia,pla
        case (3) ! Cubes !
          write(*,*) "Particle shape CUBE for particle", id_spec
          psa=((dquer(pos_spec)**3)*pi/6.0)**(1.0/3.0)
          pia=(2.0**0.5)*psa
          pla=(3.0**0.5)*psa
          if ((psa.le.0.0).or.(pia.le.0.0).or.(pla.le.0.0)) then
            write(*,*) "#### ERROR: Shape=3 (user-defined) is chosen, but no valid axes are provided."
            write(*,*) "#### SPECIES file requires SA, IA, and LA parameter greater than zero."
            error stop
          endif
          write(*,*) "SA,IA,LA:",psa,pia,pla
        case (4) ! Tetrahedrons !
          write(*,*) "Particle shape TETRAHEDRON for particle", id_spec
          pla=((dquer(pos_spec)**3)*pi*2**(0.5))**(1.0/3.0)
          pia=pla*((3.0/4.0)**(0.5))
          psa=pla*((2.0/3.0)**(0.5))
          if ((psa.le.0.0).or.(pia.le.0.0).or.(pla.le.0.0)) then
            write(*,*) "#### ERROR: Shape=4 (user-defined) is chosen, but no valid axes are provided."
            write(*,*) "#### SPECIES file requires SA, IA, and LA parameter greater than zero."
            error stop
          endif
          write(*,*) "SA,IA,LA:",psa,pia,pla
        case (5) ! Octahedrons !
          write(*,*) "Particle shape OCTAHEDRON for particle", id_spec
          psa=dquer(pos_spec)*(pi/(2.0*2.0**(0.5)))**3
          pia=psa
          pla=psa*(2.0**(0.5))
          if ((psa.le.0.0).or.(pia.le.0.0).or.(pla.le.0.0)) then
            write(*,*) "#### ERROR: Shape=5 (user-defined) is chosen, but no valid axes are provided."
            write(*,*) "#### SPECIES file requires SA, IA, and LA parameter greater than zero."
            error stop
          endif
          write(*,*) "SA,IA,LA:",psa,pia,pla
        case (6) ! Ellipsoids !
          write(*,*) "Particle shape ELLIPSOID for particle", id_spec
          psa=dquer(pos_spec)/(2.0**(1.0/3.0)) 
          pia=psa
          pla=2*pia
          if ((psa.le.0.0).or.(pia.le.0.0).or.(pla.le.0.0)) then
            write(*,*) "#### ERROR: Shape=6 (user-defined) is chosen, but no valid axes are provided."
            write(*,*) "#### SPECIES file requires SA, IA, and LA parameter greater than zero."
            error stop
          endif
          write(*,*) "SA,IA,LA:",psa,pia,pla
      end select

      ! When using the shape option, dquer is the sphere equivalent diameter
      
      f=psa/pia
      e=pia/pla
      ! Drag coefficient scheme by Bagheri & Bonadonna, 2016 to define settling velocities of other shapes (by D.Tatsii)
      if ((ishape(pos_spec).eq.2).or.((ishape(pos_spec).eq.1).and. &
        (pia.eq.psa).and.(pla.ge.20.0*pia))) then

        Fn(pos_spec)=f*f*e  ! simplified equation, validated by experiments with fibers
        Fs(pos_spec)=f*e**(1.3)   ! simplified equation, validated by experiments with fibers
      else
        Fn(pos_spec)=f*f*e*((dquer(pos_spec))**3)/(psa*pia*pla) ! Newton's regime  
        Fs(pos_spec)=f*e**(1.3)*(dquer(pos_spec)**3/(psa*pia*pla)) ! Stokes' regime
      endif
      ! Pre-compute ks and kn values needed for horizontal and average orientation (B&B Figure 12 k_(s,max))
      ks1(pos_spec)=(Fs(pos_spec)**(1./3.) + Fs(pos_spec)**(-1./3.))/2.
      ks2(pos_spec)=0.5*((Fs(pos_spec)**0.05)+(Fs(pos_spec)**(-0.36)))
      kn2(pos_spec)=10.**(alpha2*(-log10(Fn(pos_spec)))**beta2)

    else ! Spheres
      write(*,*) "Particle shape SPHERE for particle", id_spec
    endif

    do j=1,24     ! 24 hours, starting with 0-1 local time
      area_hour(pos_spec,j)=parea_hour(j)
      point_hour(pos_spec,j)=ppoint_hour(j)
    end do
    do j=1,7      ! 7 days of the week, starting with Monday
      area_dow(pos_spec,j)=parea_dow(j)
      point_dow(pos_spec,j)=ppoint_dow(j)
    end do
  endif

  i=pos_spec

  !NIK 16.02.2015
  ! Check scavenging parameters given in SPECIES file

  if (lroot) then
  ! ZHG 2016.04.07 Start of changes
    write(*,*) ' '
    if (dquer(pos_spec) .gt.0)  write(*,'(a,i3,a,a,a)')       ' SPECIES: ', &
         id_spec,'  ', species(pos_spec),'  (AEROSOL) '
    if (dquer(pos_spec) .le.0)  write(*,'(a,i3,a,a,a)')       ' SPECIES: ', &
         id_spec,'  ', species(pos_spec),'  (GAS) '

  ! Particles
  !**********
    if (dquer(pos_spec).gt.0) then
      if (ccn_aero(pos_spec) .gt. 0) then
        write(*,'(a,f5.2)') '  Particle CCN  efficiency (CCNeff):', ccn_aero(pos_spec)
      else 
        write(*,'(a)')      '  Particle CCN  efficiency (CCNeff):   OFF'
      endif
      if (in_aero(pos_spec) .gt. 0) then
        write(*,'(a,f5.2)') '  Particle  IN  efficiency (INeff) :', in_aero(pos_spec)
      else
        write(*,'(a)')      '  Particle  IN  efficiency (INeff) :   OFF'
      endif
      if (crain_aero(pos_spec) .gt. 0) then
        write(*,'(a,f5.2)') '  Particle Rain efficiency (Crain) :', crain_aero(pos_spec)
      else
        write(*,'(a)')      '  Particle Rain efficiency (Crain) :   OFF'
      endif
      if (csnow_aero(pos_spec) .gt. 0) then
        write(*,'(a,f5.2)') '  Particle Snow efficiency (Csnow) :', csnow_aero(pos_spec)
      else
        write(*,'(a)')      '  Particle Snow efficiency (Csnow) :   OFF'
      end if
      if (density(pos_spec) .gt. 0) then
        write(*,'(a)') '  Dry deposition is turned         :   ON'
        if (reldiff(pos_spec).gt.0) then
          error stop 'density>0 (SPECIES is a particle) implies reldiff <=0  '
        endif
      else
        if (reldiff(pos_spec).le.0) then
          error stop 'density<=0 (SPECIES is a gas) implies reldiff >0  '
        endif      
        write(*,'(a)') '  Dry deposition is (density<0)    :   OFF'
      end if
      if (crain_aero(pos_spec).gt.10.0 .or. csnow_aero(pos_spec).gt.10.0 .or. &
           &ccn_aero(pos_spec).gt.1.0 .or. in_aero(pos_spec).gt.1.0) then
        write(*,*) '*******************************************'
        write(*,*) ' WARNING: Particle Scavenging parameter likely out of range '
        write(*,*) '       Likely   range for Crain    0.0-10'
        write(*,*) '       Likely   range for Csnow    0.0-10'
        write(*,*) '       Physical range for CCNeff   0.0-1'
        write(*,*) '       Physical range for INeff    0.0-1'
        write(*,*) '*******************************************'
      end if
    else
  ! Gas
  !****
      if (weta_gas(pos_spec) .gt. 0 .and. wetb_gas(pos_spec).gt.0) then
        write(*,*)          '  Wet removal for gases      is turned: ON'
        write(*,*)          '  Gas below-cloud scavenging parameter A  ', &
             &weta_gas(pos_spec)
        write(*,'(a,f5.2)') '  Gas below-cloud scavenging parameter B  ', &
             &wetb_gas(pos_spec)
      else
        write(*,*)          '  Wet removal for gases      is turned: OFF '
      end if
      if (reldiff(i).gt.0.) then
        write(*,*)          '  Dry deposition for gases   is turned: ON '
      else
        write(*,*)          '  Dry deposition for gases   is turned: OFF '
      end if
      if (weta_gas(pos_spec).gt.0.) then !if wet deposition is turned on
        if (weta_gas(pos_spec).gt.1E-04 .or. weta_gas(pos_spec).lt.1E-09 .or. &
             &wetb_gas(pos_spec).gt.0.8 .or. wetb_gas(pos_spec).lt.0.4) then
          write(*,*) '*******************************************'
          write(*,*) ' WARNING: Gas below-cloud scavengig is out of likely range'
          write(*,*) '          Likely range for A is 1E-04 to 1E-08'
          write(*,*) '          Likely range for B is 0.60  to 0.80 ' 
          write(*,*) '*******************************************'
        end if
      endif

      if (((weta_gas(pos_spec).gt.0).or.(wetb_gas(pos_spec).gt.0)).and.&
           &(henry(pos_spec).le.0)) then
        if (dquer(pos_spec).le.0) goto 996 ! no particle, no henry set
      endif
    end if
  end if

  if (ndia(pos_spec).gt.maxndia) then
    maxndia=ndia(pos_spec)
    ! write(*,*) 'NDIA in SPECIES file', pos_spec, 'set to', ndia(pos_spec), 'larger than maxndia', &
    !   maxndia, 'set in par_mod.f90'
  endif
  !  if (dsigma(i).eq.0.) dsigma(i)=1.0001   ! avoid floating exception
  if (dquer(i).gt.0 .and. dsigma(i).le.1.) then !dsigma(i)=1.0001   ! avoid floating exception
    !write(*,*) '#### FLEXPART MODEL ERROR!                      ####'
    write(*,*) '#### FLEXPART MODEL WARNING                     ####'
    write(*,*) '#### in SPECIES_',aspecnumb, '                             ####'
    write(*,*) '#### from v10.4 dsigma has to be larger than 1  ####'  
    write(*,*) '#### to adapt older SPECIES files,              ####' 
    write(*,*) '#### if dsigma was < 1                          ####' 
    write(*,*) '#### use the reciprocal of the old dsigma       ####'  
    if (.not.debug_mode) then 
      error stop
    else
       write(*,*) 'debug mode: continue'
    endif
  endif

  if ((reldiff(i).gt.0.).and.(density(i).gt.0.)) then
    write(*,*) '#### FLEXPART MODEL ERROR! FILE "SPECIES"    ####'
    write(*,*) '#### IS CORRUPT. SPECIES CANNOT BE BOTH      ####'
    write(*,*) '#### PARTICLE AND GAS.                       ####'
    write(*,*) '#### SPECIES NUMBER',aspecnumb
    error stop
  endif

22 close(unitspecies)

! namelist output if requested
  if (nmlout.and.lroot) then
    open(unitspecies,file=path(2)(1:length(2))//'SPECIES_'//aspecnumb//'.namelist',access='append',status='replace',err=1000)
    write(unitspecies,nml=species_params)
    close(unitspecies)
  endif

  return

996 write(*,*) '#####################################################'
  write(*,*) '#### FLEXPART MODEL ERROR!                      #### '
  write(*,*) '#### WET DEPOSITION SWITCHED ON, BUT NO HENRYS  #### '
  write(*,*) '#### CONSTANT IS SET                            ####'
  write(*,*) '#### PLEASE MODIFY SPECIES DESCR. FILE!        #### '
  write(*,*) '#####################################################'
  error stop

998 write(*,*) '#####################################################'
  write(*,*) '#### FLEXPART MODEL ERROR!                      #### '
  write(*,*) '#### THE SPECIES FILE FOR SPECIES ', id_spec
  write(*,*) '#### CANNOT BE FOUND: CREATE FILE'
  write(*,*) '#### ',path(1)(1:length(1)),'SPECIES/SPECIES_',aspecnumb
  write(*,*) '#####################################################'
  error stop

1000 write(*,*) ' #### FLEXPART MODEL ERROR! FILE "SPECIES_',aspecnumb,'.namelist'
  write(*,*) ' #### CANNOT BE OPENED IN THE DIRECTORY       #### '
  write(*,'(a)') path(2)(1:length(2))
  error stop
end subroutine readspecies

subroutine readpartoptions

  !*****************************************************************************
  !                                                                            *
  !     This routine reads the age classes to be used for the current model    *
  !     run.                                                                   *
  !                                                                            *
  !     Author: A. Stohl                                                       *
  !     20 March 2000                                                          *
  !     HSO, 1 July 2014                                                       *
  !     Added optional namelist input                                          *
  !                                                                            *
  !*****************************************************************************
  !                                                                            *
  ! Variables:                                                                 *
  !                                                                            *
  ! Constants:                                                                 *
  !                                                                            *
  !*****************************************************************************

  implicit none

  integer :: np,stat

  ! namelist help variables
  integer :: ios

  logical ::                      &
    longitude=.false.,            &
    longitude_average=.false.,    &
    latitude=.false.,             &
    latitude_average=.false.,     &
    height=.false.,               &
    height_average=.false.,       &
    pv=.false.,                   &
    pv_average=.false.,           &
    qv=.false.,                   &
    qv_average=.false.,           &
    density=.false.,              &
    density_average=.false.,      &
    temperature=.false.,          &
    temperature_average=.false.,  &
    pressure=.false.,             &
    pressure_average=.false.,     &
    mixingheight=.false.,         &
    mixingheight_average=.false., &
    tropopause=.false.,           &
    tropopause_average=.false.,   &
    topography=.false.,           &
    topography_average=.false.,   &
    mass=.false.,                 &
    mass_average=.false.,         &
    u=.false.,                    &
    u_average=.false.,            &
    v=.false.,                    &
    v_average=.false.,            &
    w=.false.,                    &
    w_average=.false.,            &
    vsettling=.false.,            &
    vsettling_average=.false.,    &
    wetdeposition=.false.,        &
    drydeposition=.false.

  ! namelist declaration
  namelist /partoptions/  &
    longitude,            &
    longitude_average,    &
    latitude,             &
    latitude_average,     &
    height,               &
    height_average,       &
    pv,                   &
    pv_average,           &
    qv,                   &
    qv_average,           &
    density,              &
    density_average,      &
    temperature,          &
    temperature_average,  &
    pressure,             &
    pressure_average,     &
    mixingheight,         &
    mixingheight_average, &
    tropopause,           &
    tropopause_average,   &
    topography,           &
    topography_average,   &
    mass,                 &
    mass_average,         &
    u,                    &
    u_average,            &
    v,                    &
    v_average,            &
    w,                    &
    w_average,            &
    vsettling,            &
    vsettling_average,    &
    wetdeposition,        &
    drydeposition

  ! If age spectra claculation is switched on,
  ! open the AGECLASSSES file and read user options
  !************************************************

  open(unitpartoptions,file=path(1)(1:length(1))//'PARTOPTIONS',form='formatted',status='old',err=9999)

  ! try to read in as a namelist
  read(unitpartoptions,partoptions,iostat=ios)
  close(unitpartoptions)

  if (ios.ne.0) then
    write(*,*) 'Namelist error in PARTOPTIONS file', trim(path(1)(1:length(1))//'PARTOPTIONS')
    error stop
  endif
  allocate( partopt(num_partopt),stat=stat)
  if (stat.ne.0) error stop "Could not allocate partopt"
  ! Save values in particle options derived type
  !*********************************************
  partopt(1)%long_name='longitude'
  partopt(1)%name='LO'
  partopt(1)%print=longitude

  partopt(2)%long_name='longitude_average'
  partopt(2)%name='lo'
  partopt(2)%print=longitude_average
  partopt(2)%average=.true.

  partopt(3)%long_name='latitude'
  partopt(3)%name='LA'
  partopt(3)%print=latitude

  partopt(4)%long_name='latitude_average'
  partopt(4)%name='la'
  partopt(4)%print=latitude_average
  partopt(4)%average=.true.

  partopt(5)%long_name='height'
  partopt(5)%name='ZZ'
  partopt(5)%print=height

  partopt(6)%long_name='height_average'
  partopt(6)%name='zz'
  partopt(6)%print=height_average
  partopt(6)%average=.true.

  partopt(7)%long_name='pv'
  partopt(7)%name='PV'
  partopt(7)%print=pv

  partopt(8)%long_name='pv_average'
  partopt(8)%name='pv'
  partopt(8)%print=pv_average
  partopt(8)%average=.true.

  partopt(9)%long_name='qv'
  partopt(9)%name='QV'
  partopt(9)%print=qv

  partopt(10)%long_name='qv_average'
  partopt(10)%name='qv'
  partopt(10)%print=qv_average
  partopt(10)%average=.true.

  partopt(11)%long_name='density'
  partopt(11)%name='RH'
  partopt(11)%print=density

  partopt(12)%long_name='density_average'
  partopt(12)%name='rh'
  partopt(12)%print=density_average
  partopt(12)%average=.true.

  partopt(13)%long_name='temperature'
  partopt(13)%name='TT'
  partopt(13)%print=temperature

  partopt(14)%long_name='temperature_average'
  partopt(14)%name='tt'
  partopt(14)%print=temperature_average
  partopt(14)%average=.true.

  partopt(15)%long_name='pressure'
  partopt(15)%name='PR'
  partopt(15)%print=pressure

  partopt(16)%long_name='pressure_average'
  partopt(16)%name='pr'
  partopt(16)%print=pressure_average
  partopt(16)%average=.true.

  partopt(17)%long_name='mixingheight'
  partopt(17)%name='HM'
  partopt(17)%print=mixingheight

  partopt(18)%long_name='mixingheight_average'
  partopt(18)%name='hm'
  partopt(18)%print=mixingheight_average
  partopt(18)%average=.true.

  partopt(19)%long_name='tropopause'
  partopt(19)%name='TR'
  partopt(19)%print=tropopause

  partopt(20)%long_name='tropopause_average'
  partopt(20)%name='tr'
  partopt(20)%print=tropopause_average
  partopt(20)%average=.true.

  partopt(21)%long_name='topography'
  partopt(21)%name='TO'
  partopt(21)%print=topography

  partopt(22)%long_name='topography_average'
  partopt(22)%name='to'
  partopt(22)%print=topography_average
  partopt(22)%average=.true.

  partopt(23)%long_name='mass'
  partopt(23)%name='MA'
  partopt(23)%print=mass

  partopt(24)%long_name='mass_average'
  partopt(24)%name='ma'
  partopt(24)%print=mass_average
  partopt(24)%average=.true.
  
  partopt(25)%long_name='u'
  partopt(25)%name='UU'
  partopt(25)%print=u

  partopt(26)%long_name='u_average'
  partopt(26)%name='uu'
  partopt(26)%print=u_average
  partopt(26)%average=.true.

  partopt(27)%long_name='v'
  partopt(27)%name='VV'
  partopt(27)%print=v

  partopt(28)%long_name='v_average'
  partopt(28)%name='vv'
  partopt(28)%print=v_average
  partopt(28)%average=.true.

  partopt(29)%long_name='w'
  partopt(29)%name='WW'
  partopt(29)%print=w

  partopt(30)%long_name='w_average'
  partopt(30)%name='ww'
  partopt(30)%print=w_average
  partopt(30)%average=.true.

  partopt(31)%long_name='vsettling'
  partopt(31)%name='VS'
  partopt(31)%print=vsettling

  partopt(32)%long_name='vsettling_average'
  partopt(32)%name='vs'
  partopt(32)%print=vsettling_average
  partopt(32)%average=.true.

  partopt(33)%long_name='wetdeposition'
  partopt(33)%name='WD'
  partopt(33)%print=wetdeposition

  partopt(34)%long_name='drydeposition'
  partopt(34)%name='DD'
  partopt(34)%print=drydeposition
  ! Numbers are assigned to the averaged fields for proper
  ! allocation and reading in particle_mod and output_mod
  !******************************************************
  n_average=0
  do np=1,num_partopt
    if (partopt(np)%average .and. partopt(np)%print) then 
      n_average=n_average+1
      partopt(np)%i_average = n_average
      if ((partopt(np)%name.eq.'MA') .or. (partopt(np)%name.eq.'ma')) then
        n_average=n_average + (maxspec-1)
      endif
    endif
  end do

  ! write partoptions file in namelist format to output directory if requested
  if (nmlout.and.lroot) then
    open(unitpartoptions,file=path(2)(1:length(2))//'PARTOPTIONS.namelist',err=10000)
    write(unitpartoptions,nml=partoptions)
    close(unitpartoptions)
  endif


  ! Restart files, when using in combination with averaged particle output, 
  ! need to be synchronised to prevent false averages in the first step of
  ! the new run
  !************************************************************************
  if ((ipout.ne.0).and.(n_average.gt.0).and.(loutrestart.ne.-1)) then
    if (mod(loutrestart,ipoutfac*loutstep).ne.0) then
      write(*,*) '### FLEXPART MODEL ERROR! FILE COMMAND:     ###'
      write(*,*) '### LOUTRESTART NEEDS TO BE DIVISABLE BY    ###'
      write(*,*) '### LOUTSTEP*IPOUTFAC.                      ###'
      error stop
    endif
  endif

  return

9999   write(*,*) ' #### FLEXPART MODEL ERROR! FILE "PARTOPTIONS" #### '
  write(*,*) ' #### CANNOT BE OPENED IN THE DIRECTORY       #### '
  write(*,'(a)') path(1)(1:length(1))
  error stop

10000  write(*,*) ' #### FLEXPART MODEL ERROR! FILE "PARTOPTIONS" #### '
  write(*,*) ' #### CANNOT BE OPENED IN THE DIRECTORY       #### '
  write(*,'(a)') path(2)(1:length(2))
  error stop
end subroutine readpartoptions

subroutine skplin(nlines,iunit)
  !                    i      i
  !*****************************************************************************
  !                                                                            *
  !   This routine reads nlines from unit iunit and discards them              *
  !                                                                            *
  !   Authors: Petra Seibert                                                   *
  !                                                                            *
  !   31 Dec 1998                                                              *
  !                                                                            *
  !*****************************************************************************
  !                                                                            *
  ! Variables:                                                                 *
  !                                                                            *
  ! iunit   unit number from which lines are to be skipped                     *
  ! nlines  number of lines to be skipped                                      *
  !                                                                            *
  !*****************************************************************************

  implicit none

  integer :: i,iunit, nlines

  do i=1,nlines
    read(iunit,*)
  end do

end subroutine skplin

end module readoptions_mod
