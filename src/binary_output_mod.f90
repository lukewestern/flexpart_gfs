! SPDX-FileCopyrightText: FLEXPART 1998-2019, see flexpart_license.txt
! SPDX-License-Identifier: GPL-3.0-or-later
module binary_output_mod
  use point_mod
  use outg_mod
  use par_mod
  use com_mod
  use date_mod
  use windfields_mod

  implicit none

contains

subroutine writeheader

  !*****************************************************************************
  !                                                                            *
  !  This routine produces a file header containing basic information on the   *
  !  settings of the FLEXPART run.                                             *
  !  The header file is essential and must be read in by any postprocessing    *
  !  program before reading in the output data.                                *
  !                                                                            *
  !     Author: A. Stohl                                                       *
  !                                                                            *
  !     7 August 2002                                                          *
  !                                                                            *
  !*****************************************************************************
  !                                                                            *
  !  Modified to remove TRIM around the output of flexversion so that          *
  !  it will be a constant length (defined in com_mod.f90) in output header    *
  !                                                                            *
  !     Don Morton, Boreal Scientific Computing                                *
  !     07 May 2017                                                            *
  !                                                                            *
  !*****************************************************************************
  !                                                                            *
  ! Variables:                                                                 *
  !                                                                            *
  ! xlon                   longitude                                           *
  ! xl                     model x coordinate                                  *
  ! ylat                   latitude                                            *
  ! yl                     model y coordinate                                  *
  !                                                                            *
  !*****************************************************************************

  use point_mod
  use outg_mod
  use par_mod
  use com_mod

  implicit none

  integer :: jjjjmmdd,ihmmss,i,ix,jy,j
  real :: xp1,yp1,xp2,yp2


  !************************
  ! Open header output file
  !************************

  open(unitheader,file=path(2)(1:length(2))//'header', &
       form='unformatted',err=998)


  ! Write the header information
  !*****************************

  if (ldirect.eq.1) then
    write(unitheader) ibdate,ibtime, flexversion
  else
    write(unitheader) iedate,ietime, flexversion
  endif

  ! Write info on output interval, averaging time, sampling time
  !*************************************************************

  write(unitheader) loutstep,loutaver,loutsample

  ! Write information on output grid setup
  !***************************************

  write(unitheader) outlon0,outlat0,numxgrid,numygrid, &
       dxout,dyout
  write(unitheader) numzgrid,(outheight(i),i=1,numzgrid)

  call caldate(bdate,jjjjmmdd,ihmmss)
  write(unitheader) jjjjmmdd,ihmmss

  ! Write number of species, and name for each species (+extra name for depositions)
  ! Indicate the dimension of the fields (i.e., 1 for deposition fields, numzgrid for
  ! concentration fields
  !*****************************************************************************

  write(unitheader) 3*nspec,maxpointspec_act
  do i=1,nspec
    write(unitheader) 1,'WD_'//species(i)(1:7)
    write(unitheader) 1,'DD_'//species(i)(1:7)
    write(unitheader) numzgrid,species(i)
  end do

  ! Write information on release points: total number, then for each point:
  ! start, end, coordinates, # of particles, name, mass
  !************************************************************************

  write(unitheader) numpoint
  do i=1,numpoint
    write(unitheader) ireleasestart(i),ireleaseend(i),kindz(i)
    xp1=xpoint1(i)*dx+xlon0
    yp1=ypoint1(i)*dy+ylat0
    xp2=xpoint2(i)*dx+xlon0
    yp2=ypoint2(i)*dy+ylat0
    write(unitheader) xp1,yp1,xp2,yp2,zpoint1(i),zpoint2(i)
    write(unitheader) npart(i),1
    if (numpoint.le.1000) then
      write(unitheader) compoint(i)
    else
      write(unitheader) compoint(1001)
    endif
    do j=1,nspec
      write(unitheader) xmass(i,j)
      write(unitheader) xmass(i,j)
      write(unitheader) xmass(i,j)
    end do
  end do

  ! Write information on some model switches
  !*****************************************

  write(unitheader) method,lsubgrid,lconvection, &
       ind_source,ind_receptor

  ! Write age class information
  !****************************

  write(unitheader) nageclass,(lage(i),i=1,nageclass)


  ! Write topography to output file
  !********************************

  do ix=0,numxgrid-1
    write(unitheader) (oroout(ix,jy),jy=0,numygrid-1)
  end do
  close(unitheader)

  return


998   write(*,*) ' #### FLEXPART MODEL ERROR!   THE FILE         #### '
  write(*,*) ' #### '//path(2)(1:length(2))//'header'//' #### '
  write(*,*) ' #### CANNOT BE OPENED. IF A FILE WITH THIS    #### '
  write(*,*) ' #### NAME ALREADY EXISTS, DELETE IT AND START #### '
  write(*,*) ' #### THE PROGRAM AGAIN.                       #### '
  stop

end subroutine writeheader

subroutine writeheader_nest

  !*****************************************************************************
  !                                                                            *
  !  This routine produces a file header containing basic information on the   *
  !  settings of the FLEXPART run.                                             *
  !  The header file is essential and must be read in by any postprocessing    *
  !  program before reading in the output data.                                *
  !                                                                            *
  !     Author: A. Stohl                                                       *
  !                                                                            *
  !     7 August 2002                                                          *
  !                                                                            *
  !*****************************************************************************
  !                                                                            *
  !  Modified to remove TRIM around the output of flexversion so that          *
  !  it will be a constant length (defined in com_mod.f90) in output header    *
  !                                                                            *
  !     Don Morton, Boreal Scientific Computing                                *
  !     07 May 2017                                                            *
  !                                                                            *
  !*****************************************************************************
  !                                                                            *
  ! Variables:                                                                 *
  !                                                                            *
  ! xlon                   longitude                                           *
  ! xl                     model x coordinate                                  *
  ! ylat                   latitude                                            *
  ! yl                     model y coordinate                                  *
  !                                                                            *
  !*****************************************************************************

  implicit none

  integer :: jjjjmmdd,ihmmss,i,ix,jy,j
  real :: xp1,yp1,xp2,yp2


  !************************
  ! Open header output file
  !************************

  open(unitheader,file=path(2)(1:length(2))//'header_nest', &
       form='unformatted',err=998)


  ! Write the header information
  !*****************************

  if (ldirect.eq.1) then
    write(unitheader) ibdate,ibtime, flexversion
  else
    write(unitheader) iedate,ietime, flexversion
  endif

  ! Write info on output interval, averaging time, sampling time
  !*************************************************************

  write(unitheader) loutstep,loutaver,loutsample

  ! Write information on output grid setup
  !***************************************

  write(unitheader) outlon0n,outlat0n,numxgridn,numygridn, &
       dxoutn,dyoutn
  write(unitheader) numzgrid,(outheight(i),i=1,numzgrid)

  call caldate(bdate,jjjjmmdd,ihmmss)
  write(unitheader) jjjjmmdd,ihmmss

  ! Write number of species, and name for each species (+extra name for depositions)
  ! Indicate the dimension of the fields (i.e., 1 for deposition fields, numzgrid for
  ! concentration fields
  !*****************************************************************************

  write(unitheader) 3*nspec,maxpointspec_act
  do i=1,nspec
    write(unitheader) 1,'WD_'//species(i)(1:7)
    write(unitheader) 1,'DD_'//species(i)(1:7)
    write(unitheader) numzgrid,species(i)
  end do

  ! Write information on release points: total number, then for each point:
  ! start, end, coordinates, # of particles, name, mass
  !************************************************************************

  write(unitheader) numpoint
  do i=1,numpoint
    write(unitheader) ireleasestart(i),ireleaseend(i),kindz(i)
    xp1=xpoint1(i)*dx+xlon0
    yp1=ypoint1(i)*dy+ylat0
    xp2=xpoint2(i)*dx+xlon0
    yp2=ypoint2(i)*dy+ylat0
    write(unitheader) xp1,yp1,xp2,yp2,zpoint1(i),zpoint2(i)
    write(unitheader) npart(i),1
    if (numpoint.le.1000) then
      write(unitheader) compoint(i)
    else
      write(unitheader) compoint(1001)
   endif
    do j=1,nspec
      write(unitheader) xmass(i,j)
      write(unitheader) xmass(i,j)
      write(unitheader) xmass(i,j)
    end do
  end do

  ! Write information on some model switches
  !*****************************************

  write(unitheader) method,lsubgrid,lconvection, &
       ind_source,ind_receptor

  ! Write age class information
  !****************************

  write(unitheader) nageclass,(lage(i),i=1,nageclass)


  ! Write topography to output file
  !********************************

  do ix=0,numxgridn-1
    write(unitheader) (orooutn(ix,jy),jy=0,numygridn-1)
  end do
  close(unitheader)

  return


998   write(*,*) ' #### FLEXPART MODEL ERROR!   THE FILE         #### '
  write(*,*) ' #### '//path(2)(1:length(2))//'header'//' #### '
  write(*,*) ' #### CANNOT BE OPENED. IF A FILE WITH THIS    #### '
  write(*,*) ' #### NAME ALREADY EXISTS, DELETE IT AND START #### '
  write(*,*) ' #### THE PROGRAM AGAIN.                       #### '
  stop

end subroutine writeheader_nest

subroutine writeheader_nest_surf

  !*****************************************************************************
  !                                                                            *
  !  This routine produces a file header containing basic information on the   *
  !  settings of the FLEXPART run.                                             *
  !  The header file is essential and must be read in by any postprocessing    *
  !  program before reading in the output data.                                *
  !                                                                            *
  !     Author: A. Stohl                                                       *
  !                                                                            *
  !     7 August 2002                                                          *
  !                                                                            *
  !*****************************************************************************
  !                                                                            *
  !  Modified to remove TRIM around the output of flexversion so that          *
  !  it will be a constant length (defined in com_mod.f90) in output header    *
  !                                                                            *
  !     Don Morton, Boreal Scientific Computing                                *
  !     07 May 2017                                                            *
  !                                                                            *
  !*****************************************************************************
  !                                                                            *
  ! Variables:                                                                 *
  !                                                                            *
  ! xlon                   longitude                                           *
  ! xl                     model x coordinate                                  *
  ! ylat                   latitude                                            *
  ! yl                     model y coordinate                                  *
  !                                                                            *
  !*****************************************************************************

  implicit none

  integer :: jjjjmmdd,ihmmss,i,ix,jy,j
  real :: xp1,yp1,xp2,yp2


  !************************
  ! Open header output file
  !************************

  open(unitheader,file=path(2)(1:length(2))//'header_nest_grid_time', &
       form='unformatted',err=998)


  ! Write the header information
  !*****************************

  if (ldirect.eq.1) then
    write(unitheader) ibdate,ibtime,flexversion
  else
    write(unitheader) iedate,ietime,flexversion
  endif

  ! Write info on output interval, averaging time, sampling time
  !*************************************************************

  write(unitheader) loutstep,loutaver,loutsample

  ! Write information on output grid setup
  !***************************************

  write(unitheader) outlon0n,outlat0n,numxgridn,numygridn, &
       dxoutn,dyoutn
  write(unitheader) 1,(outheight(1),i=1,1)

  call caldate(bdate,jjjjmmdd,ihmmss)
  write(unitheader) jjjjmmdd,ihmmss

  ! Write number of species, and name for each species (+extra name for depositions)
  ! Indicate the dimension of the fields (i.e., 1 for deposition fields, numzgrid for
  ! concentration fields
  !*****************************************************************************

  write(unitheader) 3*nspec,maxpointspec_act
  do i=1,nspec
    write(unitheader) 1,'WD_'//species(i)(1:7)
    write(unitheader) 1,'DD_'//species(i)(1:7)
    write(unitheader) 1,species(i)
  end do

  ! Write information on release points: total number, then for each point:
  ! start, end, coordinates, # of particles, name, mass
  !************************************************************************

  write(unitheader) numpoint
  do i=1,numpoint
    write(unitheader) ireleasestart(i),ireleaseend(i),kindz(i)
    xp1=xpoint1(i)*dx+xlon0
    yp1=ypoint1(i)*dy+ylat0
    xp2=xpoint2(i)*dx+xlon0
    yp2=ypoint2(i)*dy+ylat0
    write(unitheader) xp1,yp1,xp2,yp2,zpoint1(i),zpoint2(i)
    write(unitheader) npart(i),1
    if (numpoint.le.1000) then
      write(unitheader) compoint(i)
    else
      write(unitheader) compoint(1001)
   endif
    do j=1,nspec
      write(unitheader) xmass(i,j)
      write(unitheader) xmass(i,j)
      write(unitheader) xmass(i,j)
    end do
  end do

  ! Write information on some model switches
  !*****************************************

  write(unitheader) method,lsubgrid,lconvection, &
       ind_source,ind_receptor

  ! Write age class information
  !****************************

  write(unitheader) nageclass,(lage(i),i=1,nageclass)


  ! Write topography to output file
  !********************************

  do ix=0,numxgridn-1
    write(unitheader) (orooutn(ix,jy),jy=0,numygridn-1)
  end do
  close(unitheader)

  return


998   write(*,*) ' #### FLEXPART MODEL ERROR!   THE FILE         #### '
  write(*,*) ' #### '//path(2)(1:length(2))//'header'//' #### '
  write(*,*) ' #### CANNOT BE OPENED. IF A FILE WITH THIS    #### '
  write(*,*) ' #### NAME ALREADY EXISTS, DELETE IT AND START #### '
  write(*,*) ' #### THE PROGRAM AGAIN.                       #### '
  stop

end subroutine writeheader_nest_surf

subroutine writeheader_surf

  !*****************************************************************************
  !                                                                            *
  !  This routine produces a file header containing basic information on the   *
  !  settings of the FLEXPART run.                                             *
  !  The header file is essential and must be read in by any postprocessing    *
  !  program before reading in the output data.                                *
  !                                                                            *
  !     Author: A. Stohl                                                       *
  !                                                                            *
  !     7 August 2002                                                          *
  !                                                                            *
  !*****************************************************************************
  !                                                                            *
  !  Modified to remove TRIM around the output of flexversion so that          *
  !  it will be a constant length (defined in com_mod.f90) in output header    *
  !                                                                            *
  !     Don Morton, Boreal Scientific Computing                                *
  !     07 May 2017                                                            *
  !                                                                            *
  !*****************************************************************************
  !                                                                            *
  ! Variables:                                                                 *
  !                                                                            *
  ! xlon                   longitude                                           *
  ! xl                     model x coordinate                                  *
  ! ylat                   latitude                                            *
  ! yl                     model y coordinate                                  *
  !                                                                            *
  !*****************************************************************************

  implicit none

  integer :: jjjjmmdd,ihmmss,i,ix,jy,j
  real :: xp1,yp1,xp2,yp2


  !************************
  ! Open header output file
  !************************

  open(unitheader,file=path(2)(1:length(2))//'header_grid_time', &
       form='unformatted',err=998)


  ! Write the header information
  !*****************************

  if (ldirect.eq.1) then
    write(unitheader) ibdate,ibtime, flexversion
  else
    write(unitheader) iedate,ietime, flexversion
  endif

  ! Write info on output interval, averaging time, sampling time
  !*************************************************************

  write(unitheader) loutstep,loutaver,loutsample

  ! Write information on output grid setup
  !***************************************

  write(unitheader) outlon0,outlat0,numxgrid,numygrid, &
       dxout,dyout
  write(unitheader) 1,(outheight(1),i=1,1)

  call caldate(bdate,jjjjmmdd,ihmmss)
  write(unitheader) jjjjmmdd,ihmmss

  ! Write number of species, and name for each species (+extra name for depositions)
  ! Indicate the dimension of the fields (i.e., 1 for deposition fields, numzgrid for
  ! concentration fields
  !*****************************************************************************

  write(unitheader) 3*nspec,maxpointspec_act
  do i=1,nspec
    write(unitheader) 1,'WD_'//species(i)(1:7)
    write(unitheader) 1,'DD_'//species(i)(1:7)
    write(unitheader) 1,species(i)
  end do

  ! Write information on release points: total number, then for each point:
  ! start, end, coordinates, # of particles, name, mass
  !************************************************************************

  write(unitheader) numpoint
  do i=1,numpoint
    write(unitheader) ireleasestart(i),ireleaseend(i),kindz(i)
    xp1=xpoint1(i)*dx+xlon0
    yp1=ypoint1(i)*dy+ylat0
    xp2=xpoint2(i)*dx+xlon0
    yp2=ypoint2(i)*dy+ylat0
    write(unitheader) xp1,yp1,xp2,yp2,zpoint1(i),zpoint2(i)
    write(unitheader) npart(i),1
    if (numpoint.le.1000) then
      write(unitheader) compoint(i)
    else
      write(unitheader) compoint(1001)
    endif
    do j=1,nspec
      write(unitheader) xmass(i,j)
      write(unitheader) xmass(i,j)
      write(unitheader) xmass(i,j)
    end do
  end do

  ! Write information on some model switches
  !*****************************************

  write(unitheader) method,lsubgrid,lconvection, &
       ind_source,ind_receptor

  ! Write age class information
  !****************************

  write(unitheader) nageclass,(lage(i),i=1,nageclass)


  ! Write topography to output file
  !********************************

  do ix=0,numxgrid-1
    write(unitheader) (oroout(ix,jy),jy=0,numygrid-1)
  end do
  close(unitheader)

  return


998   write(*,*) ' #### FLEXPART MODEL ERROR!   THE FILE         #### '
  write(*,*) ' #### '//path(2)(1:length(2))//'header'//' #### '
  write(*,*) ' #### CANNOT BE OPENED. IF A FILE WITH THIS    #### '
  write(*,*) ' #### NAME ALREADY EXISTS, DELETE IT AND START #### '
  write(*,*) ' #### THE PROGRAM AGAIN.                       #### '
  stop

end subroutine writeheader_surf

subroutine openreceptors

  !*****************************************************************************
  !                                                                            *
  !  This routine opens the receptor output files and writes out the receptor  *
  !  names and the receptor locations. The receptor output files are not       *
  !  closed, but kept open throughout the simulation. Concentrations are       *
  !  continuously  dumped to these files.                                      *
  !                                                                            *
  !     Author: A. Stohl                                                       *
  !                                                                            *
  !     7 August 2002                                                          *
  !                                                                            *
  !*****************************************************************************
  !                                                                            *
  ! Variables:                                                                 *
  ! numreceptor            actual number of receptor points specified          *
  ! receptornames          names of the receptor points                        *
  ! xreceptor,yreceptor    coordinates of the receptor points                  *
  !                                                                            *
  !*****************************************************************************

  use par_mod
  use com_mod

  implicit none

  integer :: j

  ! Open output file for receptor points and write out a short header
  ! containing receptor names and locations
  !******************************************************************

  if (numreceptor.ge.1) then           ! do it only if receptors are specified

  ! Concentration output
  !*********************

    if ((iout.eq.1).or.(iout.eq.3).or.(iout.eq.5)) then
      open(unitoutrecept,file=path(2)(1:length(2))//'receptor_conc', &
           form='unformatted',err=997)
      write(unitoutrecept) (receptorname(j),j=1,numreceptor)
      write(unitoutrecept) (xreceptor(j)*dx+xlon0, &
           yreceptor(j)*dy+ylat0,j=1,numreceptor)
    endif

  ! Mixing ratio output
  !********************

    if ((iout.eq.2).or.(iout.eq.3)) then
      open(unitoutreceptppt,file=path(2)(1:length(2))//'receptor_pptv', &
           form='unformatted',err=998)
      write(unitoutreceptppt) (receptorname(j),j=1,numreceptor)
      write(unitoutreceptppt) (xreceptor(j)*dx+xlon0, &
           yreceptor(j)*dy+ylat0,j=1,numreceptor)
    endif
  endif

  return
997   write(*,*) ' #### FLEXPART MODEL ERROR! THE FILE           #### '
  write(*,*) ' ####              receptor_conc               #### '
  write(*,*) ' #### CANNOT BE OPENED.                        #### '
  stop

998   write(*,*) ' #### FLEXPART MODEL ERROR! THE FILE           #### '
  write(*,*) ' ####              receptor_pptv               #### '
  write(*,*) ' #### CANNOT BE OPENED.                        #### '
  stop
end subroutine openreceptors


subroutine concoutput(itime,outnum,gridtotalunc,wetgridtotalunc, &
     drygridtotalunc)
  !                        i     i          o             o
  !       o
  !*****************************************************************************
  !                                                                            *
  !     Output of the concentration grid and the receptor concentrations.      *
  !                                                                            *
  !     Author: A. Stohl                                                       *
  !                                                                            *
  !     24 May 1995                                                            *
  !                                                                            *
  !     13 April 1999, Major update: if output size is smaller, dump output    *
  !                    in sparse matrix format; additional output of           *
  !                    uncertainty                                             *
  !                                                                            *
  !     05 April 2000, Major update: output of age classes; output for backward*
  !                    runs is time spent in grid cell times total mass of     *
  !                    species.                                                *
  !                                                                            *
  !     17 February 2002, Appropriate dimensions for backward and forward runs *
  !                       are now specified in file par_mod                    *
  !                                                                            *
  !     June 2006, write grid in sparse matrix with a single write command     *
  !                in order to save disk space                                 *
  !                                                                            *
  !     2008 new sparse matrix format                                          *
  !                                                                            *
  !*****************************************************************************
  !                                                                            *
  ! Variables:                                                                 *
  ! outnum          number of samples                                          *
  ! ncells          number of cells with non-zero concentrations               *
  ! sparse          .true. if in sparse matrix format, else .false.            *
  ! tot_mu          1 for forward, initial mass mixing ration for backw. runs  *
  !                                                                            *
  !*****************************************************************************

  use unc_mod
  use mean_mod

  implicit none

  real(kind=dp) :: jul
  integer :: itime,i,ix,jy,kz,ks,kp,l,iix,jjy,kzz,nage,jjjjmmdd,ihmmss
  integer :: sp_count_i,sp_count_r
  real :: sp_fact
  real :: outnum,densityoutrecept(maxreceptor),xl,yl
  ! RLT
  real :: densitydryrecept(maxreceptor)
  real :: factor_dryrecept(maxreceptor)

  real(dep_prec) :: auxgrid(nclassunc)
  real(sp) :: gridtotal,gridsigmatotal,gridtotalunc
  real(dep_prec) :: wetgridtotal,wetgridsigmatotal,wetgridtotalunc
  real(dep_prec) :: drygridtotal,drygridsigmatotal,drygridtotalunc
  real :: halfheight,dz,dz1,dz2,tot_mu(maxspec,maxpointspec_act)
  real,parameter :: smallnum = tiny(0.0) ! smallest number that can be handled
  real,parameter :: weightair=28.97
  logical :: sp_zer
  logical,save :: init=.true.
  character :: adate*8,atime*6
  character(len=3) :: anspec
  integer :: mind
  character(LEN=8),save :: file_stat='REPLACE'
  logical :: ldates_file
  logical :: lexist
  integer :: ierr
  character(LEN=100) :: dates_char

  ! Determine current calendar date, needed for the file name
  !**********************************************************

  jul=bdate+real(itime,kind=dp)/86400._dp
  call caldate(jul,jjjjmmdd,ihmmss)
  write(adate,'(i8.8)') jjjjmmdd
  write(atime,'(i6.6)') ihmmss

  ! Overwrite existing dates file on first call, later append to it
  ! This fixes a bug where the dates file kept growing across multiple runs

  ! If 'dates' file exists in output directory, make a backup
  inquire(file=path(2)(1:length(2))//'dates', exist=ldates_file)
  if (ldates_file.and.init) then
    open(unit=unitdates, file=path(2)(1:length(2))//'dates',form='formatted', &
         &access='sequential', status='old', action='read', iostat=ierr)
    open(unit=unittmp, file=path(2)(1:length(2))//'dates.bak', access='sequential', &
         &status='replace', action='write', form='formatted', iostat=ierr)
    do while (.true.)
      read(unitdates, '(a)', iostat=ierr) dates_char
      if (ierr.ne.0) exit
      write(unit=unittmp, fmt='(a)', iostat=ierr, advance='yes') trim(dates_char)
    end do
    close(unit=unitdates)
    close(unit=unittmp)
  end if

  open(unitdates,file=path(2)(1:length(2))//'dates', ACCESS='APPEND', STATUS=file_stat)
  write(unitdates,'(a)') adate//atime
  close(unitdates)  

  ! Overwrite existing dates file on first call, later append to it
  ! This fixes a bug where the dates file kept growing across multiple runs
  IF (init) THEN
    file_stat='OLD'
    init=.false.
  END IF


  ! For forward simulations, output fields have dimension MAXSPEC,
  ! for backward simulations, output fields have dimension MAXPOINT.
  ! Thus, make loops either about nspec, or about numpoint
  !*****************************************************************


  if (ldirect.eq.1) then
    do ks=1,nspec
      do kp=1,maxpointspec_act
        tot_mu(ks,kp)=1
      end do
    end do
  else
    do ks=1,nspec
      do kp=1,maxpointspec_act
        tot_mu(ks,kp)=xmass(kp,ks)
      end do
    end do
  endif


  !*******************************************************************
  ! Compute air density: sufficiently accurate to take it
  ! from coarse grid at some time
  ! Determine center altitude of output layer, and interpolate density
  ! data to that altitude
  !*******************************************************************

  mind=memind(2)
  do kz=1,numzgrid
    if (kz.eq.1) then
      halfheight=outheight(1)/2.
    else
      halfheight=(outheight(kz)+outheight(kz-1))/2.
    endif
    do kzz=2,nz
      if ((height(kzz-1).lt.halfheight).and. &
           (height(kzz).gt.halfheight)) goto 46
    end do
46  kzz=max(min(kzz,nz),2)
    dz1=halfheight-height(kzz-1)
    dz2=height(kzz)-halfheight
    dz=dz1+dz2
    do jy=0,numygrid-1
      do ix=0,numxgrid-1
        xl=outlon0+real(ix)*dxout
        yl=outlat0+real(jy)*dyout
        xl=(xl-xlon0)/dx
        yl=(yl-ylat0)/dy !v9.1.1 
        iix=max(min(nint(xl),nxmin1),0)
        jjy=max(min(nint(yl),nymin1),0)
  ! densityoutgrid(ix,jy,kz)=(rho(iix,jjy,kzz,2)*dz1+ &
  !      rho(iix,jjy,kzz-1,2)*dz2)/dz
        densityoutgrid(ix,jy,kz)=(rho(iix,jjy,kzz,mind)*dz1+ &
             rho(iix,jjy,kzz-1,mind)*dz2)/dz
  ! RLT
        densitydrygrid(ix,jy,kz)=(rho_dry(iix,jjy,kzz,mind)*dz1+ &
             rho_dry(iix,jjy,kzz-1,mind)*dz2)/dz  
      end do
    end do
  end do

  do i=1,numreceptor
    xl=xreceptor(i)
    yl=yreceptor(i)
    iix=max(min(nint(xl),nxmin1),0)
    jjy=max(min(nint(yl),nymin1),0)
  !densityoutrecept(i)=rho(iix,jjy,1,2)
    densityoutrecept(i)=rho(iix,jjy,1,mind)
  ! RLT
    densitydryrecept(i)=rho_dry(iix,jjy,1,mind)
  end do

  ! RLT
  ! conversion factor for output relative to dry air
  factor_drygrid=densityoutgrid/densitydrygrid
  factor_dryrecept=densityoutrecept/densitydryrecept

  ! Output is different for forward and backward simulations
  do kz=1,numzgrid
    do jy=0,numygrid-1
      do ix=0,numxgrid-1
        if (ldirect.eq.1) then
          factor3d(ix,jy,kz)=1.e12/volume(ix,jy,kz)/outnum
        else
          factor3d(ix,jy,kz)=real(abs(loutaver))/outnum
        endif
      end do
    end do
  end do

  !*********************************************************************
  ! Determine the standard deviation of the mean concentration or mixing
  ! ratio (uncertainty of the output) and the dry and wet deposition
  !*********************************************************************

  gridtotal=0.
  gridsigmatotal=0.
  gridtotalunc=0.
  wetgridtotal=0.
  wetgridsigmatotal=0.
  wetgridtotalunc=0.
  drygridtotal=0.
  drygridsigmatotal=0.
  drygridtotalunc=0.

  do ks=1,nspec

    write(anspec,'(i3.3)') ks

    if (DRYBKDEP.or.WETBKDEP) then !scavdep output
      if (DRYBKDEP) & 
      open(unitoutgrid,file=path(2)(1:length(2))//'grid_drydep_'//adate// &
           atime//'_'//anspec,form='unformatted')
      if (WETBKDEP) & 
      open(unitoutgrid,file=path(2)(1:length(2))//'grid_wetdep_'//adate// &
           atime//'_'//anspec,form='unformatted')
      write(unitoutgrid) itime
    else
      if ((iout.eq.1).or.(iout.eq.3).or.(iout.eq.5)) then
        if (ldirect.eq.1) then
          open(unitoutgrid,file=path(2)(1:length(2))//'grid_conc_'//adate// &
             atime//'_'//anspec,form='unformatted')
        else
          open(unitoutgrid,file=path(2)(1:length(2))//'grid_time_'//adate// &
             atime//'_'//anspec,form='unformatted')
        endif
        write(unitoutgrid) itime
      endif
      if ((iout.eq.2).or.(iout.eq.3)) then      ! mixing ratio
        open(unitoutgridppt,file=path(2)(1:length(2))//'grid_pptv_'//adate// &
           atime//'_'//anspec,form='unformatted')
        write(unitoutgridppt) itime
      endif
    endif ! if deposition output

    do kp=1,maxpointspec_act
      do nage=1,nageclass

        do jy=0,numygrid-1
          do ix=0,numxgrid-1

  ! WET DEPOSITION
            if ((WETDEP).and.(ldirect.gt.0)) then
              do l=1,nclassunc
                auxgrid(l)=wetgridunc(ix,jy,ks,kp,l,nage)
              end do
              call mean(auxgrid,wetgrid(ix,jy), &
                   wetgridsigma(ix,jy),nclassunc)
  ! Multiply by number of classes to get total concentration
              wetgrid(ix,jy)=wetgrid(ix,jy) &
                   *nclassunc
              wetgridtotal=wetgridtotal+wetgrid(ix,jy)
  ! Calculate standard deviation of the mean
              wetgridsigma(ix,jy)= &
                   wetgridsigma(ix,jy)* &
                   sqrt(real(nclassunc))
              wetgridsigmatotal=wetgridsigmatotal+ &
                   wetgridsigma(ix,jy)
            endif

  ! DRY DEPOSITION
            if ((DRYDEP).and.(ldirect.gt.0)) then
              do l=1,nclassunc
                auxgrid(l)=drygridunc(ix,jy,ks,kp,l,nage)
              end do
              call mean(auxgrid,drygrid(ix,jy), &
                   drygridsigma(ix,jy),nclassunc)
  ! Multiply by number of classes to get total concentration
              drygrid(ix,jy)=drygrid(ix,jy)* &
                   nclassunc
              drygridtotal=drygridtotal+drygrid(ix,jy)
  ! Calculate standard deviation of the mean
              drygridsigma(ix,jy)= &
                   drygridsigma(ix,jy)* &
                   sqrt(real(nclassunc))
              drygridsigmatotal=drygridsigmatotal+ &
                   drygridsigma(ix,jy)
            endif

  ! CONCENTRATION OR MIXING RATIO
            do kz=1,numzgrid
              do l=1,nclassunc
                auxgrid(l)=gridunc(ix,jy,kz,ks,kp,l,nage)
              end do
              call mean(auxgrid,grid(ix,jy,kz), &
                   gridsigma(ix,jy,kz),nclassunc)
  ! Multiply by number of classes to get total concentration
              grid(ix,jy,kz)= &
                   grid(ix,jy,kz)*nclassunc
              gridtotal=gridtotal+grid(ix,jy,kz)
  ! Calculate standard deviation of the mean
              gridsigma(ix,jy,kz)= &
                   gridsigma(ix,jy,kz)* &
                   sqrt(real(nclassunc))
              gridsigmatotal=gridsigmatotal+ &
                   gridsigma(ix,jy,kz)
            end do
          end do
        end do

  !*******************************************************************
  ! Generate output: may be in concentration (ng/m3) or in mixing
  ! ratio (ppt) or both
  ! Output the position and the values alternated multiplied by
  ! 1 or -1, first line is number of values, number of positions
  ! For backward simulations, the unit is seconds, stored in grid_time
  !*******************************************************************

  ! Concentration output
  !*********************
        if ((iout.eq.1).or.(iout.eq.3).or.(iout.eq.5)) then

  ! Wet deposition
          sp_count_i=0
          sp_count_r=0
          sp_fact=-1.
          sp_zer=.true.
          if ((ldirect.eq.1).and.(WETDEP)) then
            do jy=0,numygrid-1
              do ix=0,numxgrid-1
  !oncentraion greater zero
                if (wetgrid(ix,jy).gt.smallnum) then
                  if (sp_zer.eqv..true.) then ! first non zero value
                    sp_count_i=sp_count_i+1
                    sparse_dump_i(sp_count_i)=ix+jy*numxgrid
                    sp_zer=.false.
                    sp_fact=sp_fact*(-1.)
                  endif
                  sp_count_r=sp_count_r+1
                  sparse_dump_r(sp_count_r)= &
                       sp_fact*1.e12*wetgrid(ix,jy)/area(ix,jy)
                else ! concentration is zero
                  sp_zer=.true.
                endif
              end do
            end do
          else
            sp_count_i=0
            sp_count_r=0
          endif
          write(unitoutgrid) sp_count_i
          write(unitoutgrid) (sparse_dump_i(i),i=1,sp_count_i)
          write(unitoutgrid) sp_count_r
          write(unitoutgrid) (sparse_dump_r(i),i=1,sp_count_r)

  ! Dry deposition
          sp_count_i=0
          sp_count_r=0
          sp_fact=-1.
          sp_zer=.true.
          if ((ldirect.eq.1).and.(DRYDEP)) then
            do jy=0,numygrid-1
              do ix=0,numxgrid-1
                if (drygrid(ix,jy).gt.smallnum) then
                  if (sp_zer.eqv..true.) then ! first non zero value
                    sp_count_i=sp_count_i+1
                    sparse_dump_i(sp_count_i)=ix+jy*numxgrid
                    sp_zer=.false.
                    sp_fact=sp_fact*(-1.)
                  endif
                  sp_count_r=sp_count_r+1
                  sparse_dump_r(sp_count_r)= &
                       sp_fact* &
                       1.e12*drygrid(ix,jy)/area(ix,jy)
                else ! concentration is zero
                  sp_zer=.true.
                endif
              end do
            end do
          else
            sp_count_i=0
            sp_count_r=0
          endif
          write(unitoutgrid) sp_count_i
          write(unitoutgrid) (sparse_dump_i(i),i=1,sp_count_i)
          write(unitoutgrid) sp_count_r
          write(unitoutgrid) (sparse_dump_r(i),i=1,sp_count_r)



  ! Concentrations
          sp_count_i=0
          sp_count_r=0
          sp_fact=-1.
          sp_zer=.true.
          do kz=1,numzgrid
            do jy=0,numygrid-1
              do ix=0,numxgrid-1
                if (grid(ix,jy,kz).gt.smallnum) then
                  if (sp_zer.eqv..true.) then ! first non zero value
                    sp_count_i=sp_count_i+1
                    sparse_dump_i(sp_count_i)= &
                         ix+jy*numxgrid+kz*numxgrid*numygrid
                    sp_zer=.false.
                    sp_fact=sp_fact*(-1.)
                  endif
                  sp_count_r=sp_count_r+1
                  if (lparticlecountoutput) then
                    sparse_dump_r(sp_count_r)= &
                         sp_fact* &
                         grid(ix,jy,kz)
                  else
                    sparse_dump_r(sp_count_r)= &
                         sp_fact* &
                         grid(ix,jy,kz)* &
                         factor3d(ix,jy,kz)/tot_mu(ks,kp)
                  end if

                else ! concentration is zero
                  sp_zer=.true.
                endif
              end do
            end do
          end do
          write(unitoutgrid) sp_count_i
          write(unitoutgrid) (sparse_dump_i(i),i=1,sp_count_i)
          write(unitoutgrid) sp_count_r
          write(unitoutgrid) (sparse_dump_r(i),i=1,sp_count_r)



        endif !  concentration output

  ! Mixing ratio output
  !********************

        if ((iout.eq.2).or.(iout.eq.3)) then      ! mixing ratio

  ! Wet deposition
          sp_count_i=0
          sp_count_r=0
          sp_fact=-1.
          sp_zer=.true.
          if ((ldirect.eq.1).and.(WETDEP)) then
            do jy=0,numygrid-1
              do ix=0,numxgrid-1
                if (wetgrid(ix,jy).gt.smallnum) then
                  if (sp_zer.eqv..true.) then ! first non zero value
                    sp_count_i=sp_count_i+1
                    sparse_dump_i(sp_count_i)= &
                         ix+jy*numxgrid
                    sp_zer=.false.
                    sp_fact=sp_fact*(-1.)
                  endif
                  sp_count_r=sp_count_r+1
                  sparse_dump_r(sp_count_r)= &
                       sp_fact* &
                       1.e12*wetgrid(ix,jy)/area(ix,jy)
                else ! concentration is zero
                  sp_zer=.true.
                endif
              end do
            end do
          else
            sp_count_i=0
            sp_count_r=0
          endif
          write(unitoutgridppt) sp_count_i
          write(unitoutgridppt) (sparse_dump_i(i),i=1,sp_count_i)
          write(unitoutgridppt) sp_count_r
          write(unitoutgridppt) (sparse_dump_r(i),i=1,sp_count_r)


  ! Dry deposition
          sp_count_i=0
          sp_count_r=0
          sp_fact=-1.
          sp_zer=.true.
          if ((ldirect.eq.1).and.(DRYDEP)) then
            do jy=0,numygrid-1
              do ix=0,numxgrid-1
                if (drygrid(ix,jy).gt.smallnum) then
                  if (sp_zer.eqv..true.) then ! first non zero value
                    sp_count_i=sp_count_i+1
                    sparse_dump_i(sp_count_i)= &
                         ix+jy*numxgrid
                    sp_zer=.false.
                    sp_fact=sp_fact*(-1)
                  endif
                  sp_count_r=sp_count_r+1
                  sparse_dump_r(sp_count_r)= &
                       sp_fact* &
                       1.e12*drygrid(ix,jy)/area(ix,jy)
                else ! concentration is zero
                  sp_zer=.true.
                endif
              end do
            end do
          else
            sp_count_i=0
            sp_count_r=0
          endif
          write(unitoutgridppt) sp_count_i
          write(unitoutgridppt) (sparse_dump_i(i),i=1,sp_count_i)
          write(unitoutgridppt) sp_count_r
          write(unitoutgridppt) (sparse_dump_r(i),i=1,sp_count_r)


  ! Mixing ratios
          sp_count_i=0
          sp_count_r=0
          sp_fact=-1.
          sp_zer=.true.
          do kz=1,numzgrid
            do jy=0,numygrid-1
              do ix=0,numxgrid-1
                if (grid(ix,jy,kz).gt.smallnum) then
                  if (sp_zer.eqv..true.) then ! first non zero value
                    sp_count_i=sp_count_i+1
                    sparse_dump_i(sp_count_i)= &
                         ix+jy*numxgrid+kz*numxgrid*numygrid
                    sp_zer=.false.
                    sp_fact=sp_fact*(-1.)
                  endif
                  sp_count_r=sp_count_r+1
                  sparse_dump_r(sp_count_r)= &
                       sp_fact* &
                       1.e12*grid(ix,jy,kz) &
                       /volume(ix,jy,kz)/outnum* &
                       weightair/weightmolar(ks)/densityoutgrid(ix,jy,kz)
                else ! concentration is zero
                  sp_zer=.true.
                endif
              end do
            end do
          end do
          write(unitoutgridppt) sp_count_i
          write(unitoutgridppt) (sparse_dump_i(i),i=1,sp_count_i)
          write(unitoutgridppt) sp_count_r
          write(unitoutgridppt) (sparse_dump_r(i),i=1,sp_count_r)

        endif ! output for ppt

      end do
    end do

    close(unitoutgridppt)
    close(unitoutgrid)

  end do

  ! RLT Aug 2017
  ! Write out conversion factor for dry air
  inquire(file=path(2)(1:length(2))//'factor_drygrid',exist=lexist)
  if (lexist) then
    ! open and append
    open(unitoutfactor,file=path(2)(1:length(2))//'factor_drygrid',form='unformatted',&
            status='old',action='write',access='append')
  else
    ! create new
    open(unitoutfactor,file=path(2)(1:length(2))//'factor_drygrid',form='unformatted',&
            status='new',action='write')
  endif
  sp_count_i=0
  sp_count_r=0
  sp_fact=-1.
  sp_zer=.true.
  do kz=1,numzgrid
    do jy=0,numygrid-1
      do ix=0,numxgrid-1
        if (factor_drygrid(ix,jy,kz).gt.(1.+smallnum).or.factor_drygrid(ix,jy,kz).lt.(1.-smallnum)) then
          if (sp_zer.eqv..true.) then ! first value not equal to one
            sp_count_i=sp_count_i+1
            sparse_dump_i(sp_count_i)= &
                  ix+jy*numxgrid+kz*numxgrid*numygrid
            sp_zer=.false.
            sp_fact=sp_fact*(-1.)
          endif
          sp_count_r=sp_count_r+1
          sparse_dump_r(sp_count_r)= &
               sp_fact*factor_drygrid(ix,jy,kz)
        else ! factor is one
          sp_zer=.true.
        endif
      end do
    end do
  end do
  write(unitoutfactor) sp_count_i
  write(unitoutfactor) (sparse_dump_i(i),i=1,sp_count_i)
  write(unitoutfactor) sp_count_r
  write(unitoutfactor) (sparse_dump_r(i),i=1,sp_count_r)
  close(unitoutfactor)


  if (gridtotal.gt.0.) gridtotalunc=gridsigmatotal/gridtotal
  if (wetgridtotal.gt.0.) wetgridtotalunc=wetgridsigmatotal/ &
       wetgridtotal
  if (drygridtotal.gt.0.) drygridtotalunc=drygridsigmatotal/ &
       drygridtotal

  ! Dump of receptor concentrations

  if (numreceptor.gt.0 .and. (iout.eq.2 .or. iout.eq.3)  ) then
    write(unitoutreceptppt) itime
    do ks=1,nspec
      write(unitoutreceptppt) (1.e12*creceptor(i,ks)/outnum* &
           weightair/weightmolar(ks)/densityoutrecept(i),i=1,numreceptor)
    end do
  endif

  ! Dump of receptor concentrations

  if (numreceptor.gt.0) then
    write(unitoutrecept) itime
    do ks=1,nspec
      write(unitoutrecept) (1.e12*creceptor(i,ks)/outnum, &
           i=1,numreceptor)
    end do
  endif

  ! RLT Aug 2017
  ! Write out conversion factor for dry air
  if (numreceptor.gt.0) then
    inquire(file=path(2)(1:length(2))//'factor_dryreceptor',exist=lexist)
     if (lexist) then
     ! open and append
      open(unitoutfactor,file=path(2)(1:length(2))//'factor_dryreceptor',form='unformatted',&
              status='old',action='write',access='append')
    else
      ! create new
      open(unitoutfactor,file=path(2)(1:length(2))//'factor_dryreceptor',form='unformatted',&
              status='new',action='write')
    endif
    write(unitoutfactor) itime
    write(unitoutfactor) (factor_dryrecept(i),i=1,numreceptor)
    close(unitoutfactor)
  endif

  creceptor(:,:)=0.
  gridunc(:,:,:,:,:,:,:)=0.
end subroutine concoutput

subroutine concoutput_nest(itime,outnum)
  !                        i     i
  !*****************************************************************************
  !                                                                            *
  !     Output of the concentration grid and the receptor concentrations.      *
  !                                                                            *
  !     Author: A. Stohl                                                       *
  !                                                                            *
  !     24 May 1995                                                            *
  !                                                                            *
  !     13 April 1999, Major update: if output size is smaller, dump output    *
  !                    in sparse matrix format; additional output of           *
  !                    uncertainty                                             *
  !                                                                            *
  !     05 April 2000, Major update: output of age classes; output for backward*
  !                    runs is time spent in grid cell times total mass of     *
  !                    species.                                                *
  !                                                                            *
  !     17 February 2002, Appropriate dimensions for backward and forward runs *
  !                       are now specified in file par_mod                    *
  !                                                                            *
  !     June 2006, write grid in sparse matrix with a single write command     *
  !                in order to save disk space                                 *
  !                                                                            *
  !     2008 new sparse matrix format                                          *
  !                                                                            *
  !*****************************************************************************
  !                                                                            *
  ! Variables:                                                                 *
  ! outnum          number of samples                                          *
  ! ncells          number of cells with non-zero concentrations               *
  ! sparse          .true. if in sparse matrix format, else .false.            *
  ! tot_mu          1 for forward, initial mass mixing ration for backw. runs  *
  !                                                                            *
  !*****************************************************************************
  use unc_mod
  use mean_mod

  implicit none

  real(kind=dp) :: jul
  integer :: itime,i,ix,jy,kz,ks,kp,l,iix,jjy,kzz,nage,jjjjmmdd,ihmmss
  integer :: sp_count_i,sp_count_r
  real :: sp_fact
  real :: outnum,densityoutrecept(maxreceptor),xl,yl
  ! RLT
  real :: densitydryrecept(maxreceptor)
  real :: factor_dryrecept(maxreceptor)

  !real densityoutgrid(0:numxgrid-1,0:numygrid-1,numzgrid),
  !    +grid(0:numxgrid-1,0:numygrid-1,numzgrid,maxspec,maxpointspec_act,
  !    +    maxageclass)
  !real wetgrid(0:numxgrid-1,0:numygrid-1,maxspec,maxpointspec_act,
  !    +       maxageclass)
  !real drygrid(0:numxgrid-1,0:numygrid-1,maxspec,
  !    +       maxpointspec_act,maxageclass)
  !real gridsigma(0:numxgrid-1,0:numygrid-1,numzgrid,maxspec,
  !    +       maxpointspec_act,maxageclass),
  !    +     drygridsigma(0:numxgrid-1,0:numygrid-1,maxspec,
  !    +     maxpointspec_act,maxageclass),
  !    +     wetgridsigma(0:numxgrid-1,0:numygrid-1,maxspec,
  !    +     maxpointspec_act,maxageclass)
  !real factor(0:numxgrid-1,0:numygrid-1,numzgrid)
  !real sparse_dump_r(numxgrid*numygrid*numzgrid)
  !integer sparse_dump_i(numxgrid*numygrid*numzgrid)

  !real sparse_dump_u(numxgrid*numygrid*numzgrid)
  real(dep_prec) :: auxgrid(nclassunc)
  real :: halfheight,dz,dz1,dz2,tot_mu(maxspec,maxpointspec_act)
  real,parameter :: smallnum = tiny(0.0) ! smallest number that can be handled
  real,parameter :: weightair=28.97
  logical :: sp_zer
  character :: adate*8,atime*6
  character(len=3) :: anspec
  logical :: lexist
  integer :: mind 


  ! Determine current calendar date, needed for the file name
  !**********************************************************

  jul=bdate+real(itime,kind=dp)/86400._dp
  call caldate(jul,jjjjmmdd,ihmmss)
  write(adate,'(i8.8)') jjjjmmdd
  write(atime,'(i6.6)') ihmmss


  ! For forward simulations, output fields have dimension MAXSPEC,
  ! for backward simulations, output fields have dimension MAXPOINT.
  ! Thus, make loops either about nspec, or about numpoint
  !*****************************************************************


    if (ldirect.eq.1) then
       do ks=1,nspec
         do kp=1,maxpointspec_act
           tot_mu(ks,kp)=1
         end do
       end do
   else
      do ks=1,nspec
             do kp=1,maxpointspec_act
               tot_mu(ks,kp)=xmass(kp,ks)
             end do
      end do
    endif


  !*******************************************************************
  ! Compute air density: sufficiently accurate to take it
  ! from coarse grid at some time
  ! Determine center altitude of output layer, and interpolate density
  ! data to that altitude
  !*******************************************************************

  mind=memind(2)
  do kz=1,numzgrid
    if (kz.eq.1) then
      halfheight=outheight(1)/2.
    else
      halfheight=(outheight(kz)+outheight(kz-1))/2.
    endif
    do kzz=2,nz
      if ((height(kzz-1).lt.halfheight).and. &
           (height(kzz).gt.halfheight)) goto 46
    end do
46   kzz=max(min(kzz,nz),2)
    dz1=halfheight-height(kzz-1)
    dz2=height(kzz)-halfheight
    dz=dz1+dz2
    do jy=0,numygridn-1
      do ix=0,numxgridn-1
        xl=outlon0n+real(ix)*dxoutn
        yl=outlat0n+real(jy)*dyoutn
        xl=(xl-xlon0)/dx
        yl=(yl-ylat0)/dy
        iix=max(min(nint(xl),nxmin1),0)
        jjy=max(min(nint(yl),nymin1),0)
        ! densityoutgrid(ix,jy,kz)=(rho(iix,jjy,kzz,2)*dz1+ &
        !      rho(iix,jjy,kzz-1,2)*dz2)/dz
        densityoutgrid(ix,jy,kz)=(rho(iix,jjy,kzz,mind)*dz1+ &
             rho(iix,jjy,kzz-1,mind)*dz2)/dz
  ! RLT
        densitydrygrid(ix,jy,kz)=(rho_dry(iix,jjy,kzz,mind)*dz1+ &
             rho_dry(iix,jjy,kzz-1,mind)*dz2)/dz
      end do
    end do
  end do

  do i=1,numreceptor
    xl=xreceptor(i)
    yl=yreceptor(i)
    iix=max(min(nint(xl),nxmin1),0)
    jjy=max(min(nint(yl),nymin1),0)
    !densityoutrecept(i)=rho(iix,jjy,1,2)
    densityoutrecept(i)=rho(iix,jjy,1,mind)
  ! RLT
    densitydryrecept(i)=rho_dry(iix,jjy,1,mind)
  end do

  ! RLT
  ! conversion factor for output relative to dry air
  factor_drygrid=densityoutgrid/densitydrygrid
  factor_dryrecept=densityoutrecept/densitydryrecept

  ! Output is different for forward and backward simulations
    do kz=1,numzgrid
      do jy=0,numygridn-1
        do ix=0,numxgridn-1
          if (ldirect.eq.1) then
            factor3d(ix,jy,kz)=1.e12/volumen(ix,jy,kz)/outnum
          else
            factor3d(ix,jy,kz)=real(abs(loutaver))/outnum
          endif
        end do
      end do
    end do

  !*********************************************************************
  ! Determine the standard deviation of the mean concentration or mixing
  ! ratio (uncertainty of the output) and the dry and wet deposition
  !*********************************************************************

  do ks=1,nspec

  write(anspec,'(i3.3)') ks
 
  if (DRYBKDEP.or.WETBKDEP) then !scavdep output
      if (DRYBKDEP) &
      open(unitoutgrid,file=path(2)(1:length(2))//'grid_drydep_nest_'//adate// &
           atime//'_'//anspec,form='unformatted')
      if (WETBKDEP) &
      open(unitoutgrid,file=path(2)(1:length(2))//'grid_wetdep_nest_'//adate// &
           atime//'_'//anspec,form='unformatted')
      write(unitoutgrid) itime
  else
     if ((iout.eq.1).or.(iout.eq.3).or.(iout.eq.5)) then
       if (ldirect.eq.1) then
         open(unitoutgrid,file=path(2)(1:length(2))//'grid_conc_nest_' &
           //adate// &
           atime//'_'//anspec,form='unformatted')
       else
         open(unitoutgrid,file=path(2)(1:length(2))//'grid_time_nest_' &
           //adate// &
           atime//'_'//anspec,form='unformatted')
       endif
       write(unitoutgrid) itime
     endif
  endif

  if ((iout.eq.2).or.(iout.eq.3)) then      ! mixing ratio
   open(unitoutgridppt,file=path(2)(1:length(2))//'grid_pptv_nest_' &
        //adate// &
        atime//'_'//anspec,form='unformatted')

    write(unitoutgridppt) itime
  endif

  do kp=1,maxpointspec_act
  do nage=1,nageclass

    do jy=0,numygridn-1
      do ix=0,numxgridn-1

  ! WET DEPOSITION
        if ((WETDEP).and.(ldirect.gt.0)) then
            do l=1,nclassunc
              auxgrid(l)=wetgriduncn(ix,jy,ks,kp,l,nage)
            end do
            call mean(auxgrid,wetgrid(ix,jy), &
                 wetgridsigma(ix,jy),nclassunc)
  ! Multiply by number of classes to get total concentration
            wetgrid(ix,jy)=wetgrid(ix,jy) &
                 *nclassunc
  ! Calculate standard deviation of the mean
            wetgridsigma(ix,jy)= &
                 wetgridsigma(ix,jy)* &
                 sqrt(real(nclassunc))
        endif

  ! DRY DEPOSITION
        if ((DRYDEP).and.(ldirect.gt.0)) then
            do l=1,nclassunc
              auxgrid(l)=drygriduncn(ix,jy,ks,kp,l,nage)
            end do
            call mean(auxgrid,drygrid(ix,jy), &
                 drygridsigma(ix,jy),nclassunc)
  ! Multiply by number of classes to get total concentration
            drygrid(ix,jy)=drygrid(ix,jy)* &
                 nclassunc
  ! Calculate standard deviation of the mean
            drygridsigma(ix,jy)= &
                 drygridsigma(ix,jy)* &
                 sqrt(real(nclassunc))
        endif

  ! CONCENTRATION OR MIXING RATIO
        do kz=1,numzgrid
            do l=1,nclassunc
              auxgrid(l)=griduncn(ix,jy,kz,ks,kp,l,nage)
            end do
            call mean(auxgrid,grid(ix,jy,kz), &
                 gridsigma(ix,jy,kz),nclassunc)
  ! Multiply by number of classes to get total concentration
            grid(ix,jy,kz)= &
                 grid(ix,jy,kz)*nclassunc
  ! Calculate standard deviation of the mean
            gridsigma(ix,jy,kz)= &
                 gridsigma(ix,jy,kz)* &
                 sqrt(real(nclassunc))
        end do
      end do
    end do


  !*******************************************************************
  ! Generate output: may be in concentration (ng/m3) or in mixing
  ! ratio (ppt) or both
  ! Output the position and the values alternated multiplied by
  ! 1 or -1, first line is number of values, number of positions
  ! For backward simulations, the unit is seconds, stored in grid_time
  !*******************************************************************

  ! Concentration output
  !*********************
  if ((iout.eq.1).or.(iout.eq.3).or.(iout.eq.5)) then

  ! Wet deposition
         sp_count_i=0
         sp_count_r=0
         sp_fact=-1.
         sp_zer=.true.
         if ((ldirect.eq.1).and.(WETDEP)) then
         do jy=0,numygridn-1
            do ix=0,numxgridn-1
  !oncentraion greater zero
              if (wetgrid(ix,jy).gt.smallnum) then
                 if (sp_zer.eqv..true.) then ! first non zero value
                    sp_count_i=sp_count_i+1
                    sparse_dump_i(sp_count_i)=ix+jy*numxgridn
                    sp_zer=.false.
                    sp_fact=sp_fact*(-1.)
                 endif
                 sp_count_r=sp_count_r+1
                 sparse_dump_r(sp_count_r)= &
                      sp_fact*1.e12*wetgrid(ix,jy)/arean(ix,jy)
  !                sparse_dump_u(sp_count_r)=
  !+                1.e12*wetgridsigma(ix,jy,ks,kp,nage)/area(ix,jy)
              else ! concentration is zero
                  sp_zer=.true.
              endif
            end do
         end do
         else
            sp_count_i=0
            sp_count_r=0
         endif
         write(unitoutgrid) sp_count_i
         write(unitoutgrid) (sparse_dump_i(i),i=1,sp_count_i)
         write(unitoutgrid) sp_count_r
         write(unitoutgrid) (sparse_dump_r(i),i=1,sp_count_r)
  !       write(unitoutgrid) sp_count_u
  !       write(unitoutgrid) (sparse_dump_u(i),i=1,sp_count_r)

  ! Dry deposition
         sp_count_i=0
         sp_count_r=0
         sp_fact=-1.
         sp_zer=.true.
         if ((ldirect.eq.1).and.(DRYDEP)) then
          do jy=0,numygridn-1
            do ix=0,numxgridn-1
              if (drygrid(ix,jy).gt.smallnum) then
                 if (sp_zer.eqv..true.) then ! first non zero value
                    sp_count_i=sp_count_i+1
                    sparse_dump_i(sp_count_i)=ix+jy*numxgridn
                    sp_zer=.false.
                    sp_fact=sp_fact*(-1.)
                 endif
                 sp_count_r=sp_count_r+1
                 sparse_dump_r(sp_count_r)= &
                      sp_fact* &
                      1.e12*drygrid(ix,jy)/arean(ix,jy)
  !                sparse_dump_u(sp_count_r)=
  !+                1.e12*drygridsigma(ix,jy,ks,kp,nage)/area(ix,jy)
              else ! concentration is zero
                  sp_zer=.true.
              endif
            end do
          end do
         else
            sp_count_i=0
            sp_count_r=0
         endif
         write(unitoutgrid) sp_count_i
         write(unitoutgrid) (sparse_dump_i(i),i=1,sp_count_i)
         write(unitoutgrid) sp_count_r
         write(unitoutgrid) (sparse_dump_r(i),i=1,sp_count_r)
  !       write(*,*) sp_count_u
  !       write(unitoutgrid) (sparse_dump_u(i),i=1,sp_count_r)



  ! Concentrations
         sp_count_i=0
         sp_count_r=0
         sp_fact=-1.
         sp_zer=.true.
          do kz=1,numzgrid
            do jy=0,numygridn-1
              do ix=0,numxgridn-1
                if (grid(ix,jy,kz).gt.smallnum) then
                  if (sp_zer.eqv..true.) then ! first non zero value
                    sp_count_i=sp_count_i+1
                    sparse_dump_i(sp_count_i)= &
                         ix+jy*numxgridn+kz*numxgridn*numygridn
                    sp_zer=.false.
                    sp_fact=sp_fact*(-1.)
                   endif
                   sp_count_r=sp_count_r+1
                   sparse_dump_r(sp_count_r)= &
                        sp_fact* &
                        grid(ix,jy,kz)* &
                        factor3d(ix,jy,kz)/tot_mu(ks,kp)
  !                 if ((factor(ix,jy,kz)/tot_mu(ks,kp)).eq.0)
  !    +              write (*,*) factor(ix,jy,kz),tot_mu(ks,kp),ks,kp
  !                sparse_dump_u(sp_count_r)=
  !+               ,gridsigma(ix,jy,kz,ks,kp,nage)*
  !+               factor(ix,jy,kz)/tot_mu(ks,kp)
              else ! concentration is zero
                  sp_zer=.true.
              endif
              end do
            end do
          end do
         write(unitoutgrid) sp_count_i
         write(unitoutgrid) (sparse_dump_i(i),i=1,sp_count_i)
         write(unitoutgrid) sp_count_r
         write(unitoutgrid) (sparse_dump_r(i),i=1,sp_count_r)
  !       write(unitoutgrid) sp_count_u
  !       write(unitoutgrid) (sparse_dump_u(i),i=1,sp_count_r)



    endif !  concentration output

  ! Mixing ratio output
  !********************

  if ((iout.eq.2).or.(iout.eq.3)) then      ! mixing ratio

  ! Wet deposition
         sp_count_i=0
         sp_count_r=0
         sp_fact=-1.
         sp_zer=.true.
         if ((ldirect.eq.1).and.(WETDEP)) then
          do jy=0,numygridn-1
            do ix=0,numxgridn-1
                if (wetgrid(ix,jy).gt.smallnum) then
                  if (sp_zer.eqv..true.) then ! first non zero value
                    sp_count_i=sp_count_i+1
                    sparse_dump_i(sp_count_i)= &
                         ix+jy*numxgridn
                    sp_zer=.false.
                    sp_fact=sp_fact*(-1.)
                 endif
                 sp_count_r=sp_count_r+1
                 sparse_dump_r(sp_count_r)= &
                      sp_fact* &
                      1.e12*wetgrid(ix,jy)/arean(ix,jy)
  !                sparse_dump_u(sp_count_r)=
  !    +            ,1.e12*wetgridsigma(ix,jy,ks,kp,nage)/area(ix,jy)
              else ! concentration is zero
                  sp_zer=.true.
              endif
            end do
          end do
         else
           sp_count_i=0
           sp_count_r=0
         endif
         write(unitoutgridppt) sp_count_i
         write(unitoutgridppt) (sparse_dump_i(i),i=1,sp_count_i)
         write(unitoutgridppt) sp_count_r
         write(unitoutgridppt) (sparse_dump_r(i),i=1,sp_count_r)
  !       write(unitoutgridppt) sp_count_u
  !       write(unitoutgridppt) (sparse_dump_u(i),i=1,sp_count_r)


  ! Dry deposition
         sp_count_i=0
         sp_count_r=0
         sp_fact=-1.
         sp_zer=.true.
         if ((ldirect.eq.1).and.(DRYDEP)) then
          do jy=0,numygridn-1
            do ix=0,numxgridn-1
                if (drygrid(ix,jy).gt.smallnum) then
                  if (sp_zer.eqv..true.) then ! first non zero value
                    sp_count_i=sp_count_i+1
                    sparse_dump_i(sp_count_i)= &
                         ix+jy*numxgridn
                    sp_zer=.false.
                    sp_fact=sp_fact*(-1)
                 endif
                 sp_count_r=sp_count_r+1
                 sparse_dump_r(sp_count_r)= &
                      sp_fact* &
                      1.e12*drygrid(ix,jy)/arean(ix,jy)
  !                sparse_dump_u(sp_count_r)=
  !    +            ,1.e12*drygridsigma(ix,jy,ks,kp,nage)/area(ix,jy)
              else ! concentration is zero
                  sp_zer=.true.
              endif
            end do
          end do
         else
           sp_count_i=0
           sp_count_r=0
         endif
         write(unitoutgridppt) sp_count_i
         write(unitoutgridppt) (sparse_dump_i(i),i=1,sp_count_i)
         write(unitoutgridppt) sp_count_r
         write(unitoutgridppt) (sparse_dump_r(i),i=1,sp_count_r)
  !       write(unitoutgridppt) sp_count_u
  !       write(unitoutgridppt) (sparse_dump_u(i),i=1,sp_count_r)


  ! Mixing ratios
         sp_count_i=0
         sp_count_r=0
         sp_fact=-1.
         sp_zer=.true.
          do kz=1,numzgrid
            do jy=0,numygridn-1
              do ix=0,numxgridn-1
                if (grid(ix,jy,kz).gt.smallnum) then
                  if (sp_zer.eqv..true.) then ! first non zero value
                    sp_count_i=sp_count_i+1
                    sparse_dump_i(sp_count_i)= &
                         ix+jy*numxgridn+kz*numxgridn*numygridn
                    sp_zer=.false.
                    sp_fact=sp_fact*(-1.)
                 endif
                 sp_count_r=sp_count_r+1
                 sparse_dump_r(sp_count_r)= &
                      sp_fact* &
                      1.e12*grid(ix,jy,kz) &
                      /volumen(ix,jy,kz)/outnum* &
                      weightair/weightmolar(ks)/densityoutgrid(ix,jy,kz)
  !                sparse_dump_u(sp_count_r)=
  !+              ,1.e12*gridsigma(ix,jy,kz,ks,kp,nage)/volume(ix,jy,kz)/
  !+              outnum*weightair/weightmolar(ks)/
  !+              densityoutgrid(ix,jy,kz)
              else ! concentration is zero
                  sp_zer=.true.
              endif
              end do
            end do
          end do
         write(unitoutgridppt) sp_count_i
         write(unitoutgridppt) (sparse_dump_i(i),i=1,sp_count_i)
         write(unitoutgridppt) sp_count_r
         write(unitoutgridppt) (sparse_dump_r(i),i=1,sp_count_r)
  !       write(unitoutgridppt) sp_count_u
  !       write(unitoutgridppt) (sparse_dump_u(i),i=1,sp_count_r)

      endif ! output for ppt

  end do
  end do

    close(unitoutgridppt)
    close(unitoutgrid)

  end do

  ! RLT Aug 2017
  ! Write out conversion factor for dry air
  inquire(file=path(2)(1:length(2))//'factor_drygrid_nest',exist=lexist)
  if (lexist) then
    ! open and append
    open(unitoutfactor,file=path(2)(1:length(2))//'factor_drygrid_nest',form='unformatted',&
            status='old',action='write',access='append')
  else
    ! create new
    open(unitoutfactor,file=path(2)(1:length(2))//'factor_drygrid_nest',form='unformatted',&
            status='new',action='write')
  endif
  sp_count_i=0
  sp_count_r=0
  sp_fact=-1.
  sp_zer=.true.
  do kz=1,numzgrid
    do jy=0,numygridn-1
      do ix=0,numxgridn-1
        if (factor_drygrid(ix,jy,kz).gt.(1.+smallnum).or.factor_drygrid(ix,jy,kz).lt.(1.-smallnum)) then
          if (sp_zer.eqv..true.) then ! first value not equal to one
            sp_count_i=sp_count_i+1
            sparse_dump_i(sp_count_i)= &
                  ix+jy*numxgridn+kz*numxgridn*numygridn
            sp_zer=.false.
            sp_fact=sp_fact*(-1.)
          endif
          sp_count_r=sp_count_r+1
          sparse_dump_r(sp_count_r)= &
               sp_fact*factor_drygrid(ix,jy,kz)
        else ! factor is one
          sp_zer=.true.
        endif
      end do
    end do
  end do
  write(unitoutfactor) sp_count_i
  write(unitoutfactor) (sparse_dump_i(i),i=1,sp_count_i)
  write(unitoutfactor) sp_count_r
  write(unitoutfactor) (sparse_dump_r(i),i=1,sp_count_r)
  close(unitoutfactor)

  creceptor(:,:)=0.
  griduncn(:,:,:,:,:,:,:)=0.
end subroutine concoutput_nest

subroutine concoutput_inversion(itime,outnum,gridtotalunc,wetgridtotalunc, &
     drygridtotalunc)
  !                        i     i          o             o
  !       o
  !*****************************************************************************
  !                                                                            *
  !     Output of the concentration grid and the receptor concentrations.      *
  !                                                                            *
  !     Author: A. Stohl                                                       *
  !                                                                            *
  !     24 May 1995                                                            *
  !                                                                            *
  !     13 April 1999, Major update: if output size is smaller, dump output    *
  !                    in sparse matrix format; additional output of           *
  !                    uncertainty                                             *
  !                                                                            *
  !     05 April 2000, Major update: output of age classes; output for backward*
  !                    runs is time spent in grid cell times total mass of     *
  !                    species.                                                *
  !                                                                            *
  !     17 February 2002, Appropriate dimensions for backward and forward runs *
  !                       are now specified in file par_mod                    *
  !                                                                            *
  !     June 2006, write grid in sparse matrix with a single write command     *
  !                in order to save disk space                                 *
  !                                                                            *
  !     2008 new sparse matrix format                                          *
  !
  !     January 2017,  Separate files by release but include all timesteps
  !                                                                            *
  !*****************************************************************************
  !                                                                            *
  ! Variables:                                                                 *
  ! outnum          number of samples                                          *
  ! ncells          number of cells with non-zero concentrations               *
  ! sparse          .true. if in sparse matrix format, else .false.            *
  ! tot_mu          1 for forward, initial mass mixing ration for backw. runs  *
  !                                                                            *
  !*****************************************************************************

  use unc_mod
  use mean_mod

  implicit none

  real(kind=dp) :: jul
  integer :: itime,i,ix,jy,kz,ks,kp,l,iix,jjy,kzz,nage,jjjjmmdd,ihmmss
  integer :: sp_count_i,sp_count_r
  real :: sp_fact
  real :: outnum,densityoutrecept(maxreceptor),xl,yl
  ! RLT
  real :: densitydryrecept(maxreceptor)
  real :: factor_dryrecept(maxreceptor)


  real(dep_prec) :: auxgrid(nclassunc)
  real(sp) :: gridtotal,gridsigmatotal,gridtotalunc
  real(dep_prec) :: wetgridtotal,wetgridsigmatotal,wetgridtotalunc
  real(dep_prec) :: drygridtotal,drygridsigmatotal,drygridtotalunc
  real :: halfheight,dz,dz1,dz2,tot_mu(maxspec,maxpointspec_act)
  real,parameter :: smallnum = tiny(0.0) ! smallest number that can be handled
  real,parameter :: weightair=28.97
  logical :: sp_zer
  character :: adate*8,atime*6
  character(len=3) :: anspec
  logical :: lexist
  character :: areldate*8,areltime*6
  logical,save :: lstart=.true.
  logical,save,allocatable,dimension(:) :: lstartrel
  integer :: ierr
  character(LEN=100) :: dates_char
  integer, parameter :: unitrelnames=654


  if(lstart) then
    allocate(lstartrel(maxpointspec_act))
    lstartrel(:)=.true.
  endif
  print*, 'lstartrel = ',lstartrel

  if (verbosity.eq.1) then
    print*,'inside concoutput_inversion '
    CALL SYSTEM_CLOCK(count_clock)
    WRITE(*,*) 'SYSTEM_CLOCK',count_clock - count_clock0   
  endif

  ! Determine current calendar date
  !**********************************************************

  jul=bdate+real(itime,kind=dp)/86400._dp
  call caldate(jul,jjjjmmdd,ihmmss)
  write(adate,'(i8.8)') jjjjmmdd
  write(atime,'(i6.6)') ihmmss


  ! Prepare output files for dates
  !**********************************************************

  ! Overwrite existing dates file on first call, later append to it
  ! If 'dates' file exists in output directory, copy to new file dates.old
  inquire(file=path(2)(1:length(2))//'dates', exist=lexist)
  if (lexist.and.lstart) then
    ! copy contents of existing dates file to dates.old
    print*, 'warning: replacing old dates file'
    open(unit=unitdates, file=path(2)(1:length(2))//'dates',form='formatted', &
         &access='sequential', status='old', action='read', iostat=ierr)
    open(unit=unittmp, file=path(2)(1:length(2))//'dates.old', access='sequential', &
         &status='replace', action='write', form='formatted', iostat=ierr)
    do while (.true.)
      read(unitdates, '(a)', iostat=ierr) dates_char
      if (ierr.ne.0) exit
      write(unit=unittmp, fmt='(a)', iostat=ierr, advance='yes') trim(dates_char)
    end do
    close(unit=unitdates)
    close(unit=unittmp)
    ! create new dates file
    open(unit=unitdates, file=path(2)(1:length(2))//'dates',form='formatted', &
         &access='sequential', status='replace', iostat=ierr)
    close(unit=unitdates)
  endif

  open(unitdates,file=path(2)(1:length(2))//'dates', ACCESS='APPEND')
  write(unitdates,'(a)') adate//atime
  close(unitdates)

  !CGZ: Make a filename with names of releases
  if (lstart) then
    open(unit=unitrelnames, file=path(2)(1:length(2))//'releases_out',form='formatted', &
         &access='sequential', status='replace', iostat=ierr)
    close(unitrelnames)
  endif

  print*, 'after creating dates files: lstart = ',lstart
  !  print*, 'outnum:',outnum
  !  print*, 'datetime:',adate//atime


  ! For forward simulations, output fields have dimension MAXSPEC,
  ! for backward simulations, output fields have dimension MAXPOINT.
  ! Thus, make loops either about nspec, or about numpoint
  !*****************************************************************


  if (ldirect.eq.1) then
    do ks=1,nspec
      do kp=1,maxpointspec_act
        tot_mu(ks,kp)=1
      end do
    end do
  else
    do ks=1,nspec
      do kp=1,maxpointspec_act
        tot_mu(ks,kp)=xmass(kp,ks)
      end do
    end do
  endif


  if (verbosity.eq.1) then
    print*,'concoutput_inversion 2'
    CALL SYSTEM_CLOCK(count_clock)
    WRITE(*,*) 'SYSTEM_CLOCK',count_clock - count_clock0   
  endif

  !*******************************************************************
  ! Compute air density: sufficiently accurate to take it
  ! from coarse grid at some time
  ! Determine center altitude of output layer, and interpolate density
  ! data to that altitude
  !*******************************************************************

  do kz=1,numzgrid
    if (kz.eq.1) then
      halfheight=outheight(1)/2.
    else
      halfheight=(outheight(kz)+outheight(kz-1))/2.
    endif
    do kzz=2,nz
      if ((height(kzz-1).lt.halfheight).and. &
           (height(kzz).gt.halfheight)) goto 46
    end do
46  kzz=max(min(kzz,nz),2)
    dz1=halfheight-height(kzz-1)
    dz2=height(kzz)-halfheight
    dz=dz1+dz2
    do jy=0,numygrid-1
      do ix=0,numxgrid-1
        xl=outlon0+real(ix)*dxout
        yl=outlat0+real(jy)*dyout
        xl=(xl-xlon0)/dx
        yl=(yl-ylat0)/dy
        iix=max(min(nint(xl),nxmin1),0)
        jjy=max(min(nint(yl),nymin1),0)
        densityoutgrid(ix,jy,kz)=(rho(iix,jjy,kzz,2)*dz1+ &
             rho(iix,jjy,kzz-1,2)*dz2)/dz
  ! RLT
        densitydrygrid(ix,jy,kz)=(rho_dry(iix,jjy,kzz,2)*dz1+ &
             rho_dry(iix,jjy,kzz-1,2)*dz2)/dz
      end do
    end do
  end do

  do i=1,numreceptor
    xl=xreceptor(i)
    yl=yreceptor(i)
    iix=max(min(nint(xl),nxmin1),0)
    jjy=max(min(nint(yl),nymin1),0)
    densityoutrecept(i)=rho(iix,jjy,1,2)
  ! RLT
    densitydryrecept(i)=rho_dry(iix,jjy,1,2)
  end do

  ! RLT
  ! conversion factor for output relative to dry air
  factor_drygrid=densityoutgrid/densitydrygrid
  factor_dryrecept=densityoutrecept/densitydryrecept

  ! Output is different for forward and backward simulations
  do kz=1,numzgrid
    do jy=0,numygrid-1
      do ix=0,numxgrid-1
        if (ldirect.eq.1) then
          factor3d(ix,jy,kz)=1.e12/volume(ix,jy,kz)/outnum
        else
          factor3d(ix,jy,kz)=real(abs(loutaver))/outnum
        endif
      end do
    end do
  end do

  !*********************************************************************
  ! Determine the standard deviation of the mean concentration or mixing
  ! ratio (uncertainty of the output) and the dry and wet deposition
  !*********************************************************************

  if (verbosity.eq.1) then
    print*,'concoutput_inversion 3 (sd)'
    CALL SYSTEM_CLOCK(count_clock)
    WRITE(*,*) 'SYSTEM_CLOCK',count_clock - count_clock0   
  endif
  gridtotal=0.
  gridsigmatotal=0.
  gridtotalunc=0.
  wetgridtotal=0.
  wetgridsigmatotal=0.
  wetgridtotalunc=0.
  drygridtotal=0.
  drygridsigmatotal=0.
  drygridtotalunc=0.

  do ks=1,nspec

    write(anspec,'(i3.3)') ks

    ! loop over releases
    do kp=1,maxpointspec_act

      print*, 'itime = ',itime
      !print*, 'lage(1) = ',lage(1)
      print*, 'ireleasestart(kp) = ',ireleasestart(kp)
      print*, 'ireleaseend(kp) = ',ireleaseend(kp)

      ! check itime is within release and backward trajectory length
      if (nageclass.eq.1) then
        if ((itime.gt.ireleaseend(kp)).or.(itime.lt.(ireleasestart(kp)-lage(1)))) then
          go to 10
        endif
      endif

      ! calculate date of release for filename
      jul=bdate+real(ireleasestart(kp),kind=dp)/86400._dp    ! this is the current day
      call caldate(jul,jjjjmmdd,ihmmss)     
      write(areldate,'(i8.8)') jjjjmmdd
      write(areltime,'(i6.6)') ihmmss
      print*, 'areldate/areltime = ',areldate//areltime

      ! calculate date of field
      jul=bdate+real(itime,kind=dp)/86400._dp
      call caldate(jul,jjjjmmdd,ihmmss)
      write(adate,'(i8.8)') jjjjmmdd
      write(atime,'(i6.6)') ihmmss

      if ((iout.eq.1).or.(iout.eq.3).or.(iout.eq.5)) then
        if (ldirect.eq.1) then
          ! concentrations
          inquire(file=path(2)(1:length(2))//'grid_conc_'//areldate// &
                  areltime//'_'//anspec,exist=lexist)
          if(lexist.and..not.lstartrel(kp)) then
            ! open and append to existing file
            open(unitoutgrid,file=path(2)(1:length(2))//'grid_conc_'//areldate// &
                 areltime//'_'//anspec,form='unformatted',status='old',action='write',access='append')
          else
            ! open new file
            open(unitoutgrid,file=path(2)(1:length(2))//'grid_conc_'//areldate// &
                 areltime//'_'//anspec,form='unformatted',status='replace',action='write')
          endif
        else
          ! residence times
          inquire(file=path(2)(1:length(2))//'grid_time_'//areldate// &
                  areltime//'_'//anspec,exist=lexist)
          if(lexist.and..not.lstartrel(kp)) then
            ! open and append to existing file
            open(unitoutgrid,file=path(2)(1:length(2))//'grid_time_'//areldate// &
                 areltime//'_'//anspec,form='unformatted',status='old',action='write',access='append')
          else
            ! open new file
            open(unitoutgrid,file=path(2)(1:length(2))//'grid_time_'//areldate// &
                 areltime//'_'//anspec,form='unformatted',status='replace',action='write')
            ! add part of the filename to a file (similar to dates) for easier post-processing
            open(unit=unitrelnames, file=path(2)(1:length(2))//'releases_out',form='formatted', &
                 & access='APPEND', iostat=ierr)
            write(unitrelnames,'(a)') areldate//areltime//'_'//anspec
            close(unitrelnames)
          endif
        endif
        write(unitoutgrid) jjjjmmdd
        write(unitoutgrid) ihmmss
      endif

      if ((iout.eq.2).or.(iout.eq.3)) then      
        ! mixing ratio
        inquire(file=path(2)(1:length(2))//'grid_pptv_'//areldate// &
                areltime//'_'//anspec,exist=lexist)
        if(lexist.and..not.lstartrel(kp)) then
          ! open and append to existing file
          open(unitoutgridppt,file=path(2)(1:length(2))//'grid_pptv_'//areldate// &
               areltime//'_'//anspec,form='unformatted',status='old',action='write',access='append')
        else
          ! open new file
          open(unitoutgridppt,file=path(2)(1:length(2))//'grid_pptv_'//areldate// &
               areltime//'_'//anspec,form='unformatted',status='replace',action='write')
        endif   
        write(unitoutgridppt) jjjjmmdd
        write(unitoutgridppt) ihmmss
      endif

      lstartrel(kp)=.false.

      do nage=1,nageclass

        do jy=0,numygrid-1
          do ix=0,numxgrid-1

  ! CONCENTRATION OR MIXING RATIO
            do kz=1,numzgrid
              do l=1,nclassunc
                auxgrid(l)=gridunc(ix,jy,kz,ks,kp,l,nage)
              end do
              call mean(auxgrid,grid(ix,jy,kz), &
                   gridsigma(ix,jy,kz),nclassunc)
  ! Multiply by number of classes to get total concentration
              grid(ix,jy,kz)= &
                   grid(ix,jy,kz)*nclassunc
              gridtotal=gridtotal+grid(ix,jy,kz)
  ! Calculate standard deviation of the mean
              gridsigma(ix,jy,kz)= &
                   gridsigma(ix,jy,kz)* &
                   sqrt(real(nclassunc))
              gridsigmatotal=gridsigmatotal+ &
                   gridsigma(ix,jy,kz)
            end do
          end do
        end do


  !*******************************************************************
  ! Generate output: may be in concentration (ng/m3) or in mixing
  ! ratio (ppt) or both
  ! Output the position and the values alternated multiplied by
  ! 1 or -1, first line is number of values, number of positions
  ! For backward simulations, the unit is seconds, stored in grid_time
  !*******************************************************************

        if (verbosity.eq.1) then
          print*,'concoutput_inversion 4 (output)'
          CALL SYSTEM_CLOCK(count_clock)
          WRITE(*,*) 'SYSTEM_CLOCK',count_clock - count_clock0   
        endif

  ! Concentration output
  !*********************

        if ((iout.eq.1).or.(iout.eq.3).or.(iout.eq.5)) then

          if (verbosity.eq.1) then
            print*,'concoutput_inversion (Wet deposition)'
            CALL SYSTEM_CLOCK(count_clock)
            WRITE(*,*) 'SYSTEM_CLOCK',count_clock - count_clock0   
          endif

          if (verbosity.eq.1) then
            print*,'concoutput_inversion (Concentrations)'
            CALL SYSTEM_CLOCK(count_clock)
            WRITE(*,*) 'SYSTEM_CLOCK',count_clock - count_clock0   
          endif

  ! Concentrations

  ! surf_only write only 1st layer 

          sp_count_i=0
          sp_count_r=0
          sp_fact=-1.
          sp_zer=.true.
          do kz=1,1
            do jy=0,numygrid-1
              do ix=0,numxgrid-1
                if (grid(ix,jy,kz).gt.smallnum) then
                  if (sp_zer.eqv..true.) then ! first non zero value
                    sp_count_i=sp_count_i+1
                    sparse_dump_i(sp_count_i)= &
                         ix+jy*numxgrid+kz*numxgrid*numygrid
                    sp_zer=.false.
                    sp_fact=sp_fact*(-1.)
                  endif
                  sp_count_r=sp_count_r+1
                  sparse_dump_r(sp_count_r)= &
                       sp_fact* &
                       grid(ix,jy,kz)* &
                       factor3d(ix,jy,kz)/tot_mu(ks,kp)
                  sparse_dump_u(sp_count_r)= &
                       gridsigma(ix,jy,kz)* &
                       factor3d(ix,jy,kz)/tot_mu(ks,kp)

                else ! concentration is zero
                  sp_zer=.true.
                endif
              end do
            end do
          end do
          write(unitoutgrid) sp_count_i
          write(unitoutgrid) (sparse_dump_i(i),i=1,sp_count_i)
          write(unitoutgrid) sp_count_r
          write(unitoutgrid) (sparse_dump_r(i),i=1,sp_count_r)

        endif !  concentration output

  ! Mixing ratio output
  !********************

        if ((iout.eq.2).or.(iout.eq.3)) then      ! mixing ratio

  ! Mixing ratios

  ! surf_only write only 1st layer 

          sp_count_i=0
          sp_count_r=0
          sp_fact=-1.
          sp_zer=.true.
          do kz=1,1
            do jy=0,numygrid-1
              do ix=0,numxgrid-1
                if (grid(ix,jy,kz).gt.smallnum) then
                  if (sp_zer.eqv..true.) then ! first non zero value
                    sp_count_i=sp_count_i+1
                    sparse_dump_i(sp_count_i)= &
                         ix+jy*numxgrid+kz*numxgrid*numygrid
                    sp_zer=.false.
                    sp_fact=sp_fact*(-1.)
                  endif
                  sp_count_r=sp_count_r+1
                  sparse_dump_r(sp_count_r)= &
                       sp_fact* &
                       1.e12*grid(ix,jy,kz) &
                       /volume(ix,jy,kz)/outnum* &
                       weightair/weightmolar(ks)/densityoutgrid(ix,jy,kz)
                  sparse_dump_u(sp_count_r)= &
                       1.e12*gridsigma(ix,jy,kz)/volume(ix,jy,kz)/ &
                       outnum*weightair/weightmolar(ks)/ &
                       densityoutgrid(ix,jy,kz)
                else ! concentration is zero
                  sp_zer=.true.
                endif
              end do
            end do
          end do
          write(unitoutgridppt) sp_count_i
          write(unitoutgridppt) (sparse_dump_i(i),i=1,sp_count_i)
          write(unitoutgridppt) sp_count_r
          write(unitoutgridppt) (sparse_dump_r(i),i=1,sp_count_r)

        endif ! output for ppt

      end do  ! nageclass

      close(unitoutgridppt)
      close(unitoutgrid)

      ! itime is outside range
10    continue

    end do  ! maxpointspec_act

  end do  ! nspec

  ! RLT Aug 2017
  ! Write out conversion factor for dry air
  inquire(file=path(2)(1:length(2))//'factor_drygrid',exist=lexist)
  if (lexist.and..not.lstart) then
    ! open and append
    open(unitoutfactor,file=path(2)(1:length(2))//'factor_drygrid',form='unformatted',&
            status='old',action='write',access='append')
  else
    ! create new
    open(unitoutfactor,file=path(2)(1:length(2))//'factor_drygrid',form='unformatted',&
            status='replace',action='write')
  endif
  sp_count_i=0
  sp_count_r=0
  sp_fact=-1.
  sp_zer=.true.
  do kz=1,1
    do jy=0,numygrid-1
      do ix=0,numxgrid-1
        if (factor_drygrid(ix,jy,kz).gt.(1.+smallnum).or.factor_drygrid(ix,jy,kz).lt.(1.-smallnum)) then
          if (sp_zer.eqv..true.) then ! first value not equal to one
            sp_count_i=sp_count_i+1
            sparse_dump_i(sp_count_i)= &
                  ix+jy*numxgrid+kz*numxgrid*numygrid
            sp_zer=.false.
            sp_fact=sp_fact*(-1.)
          endif
          sp_count_r=sp_count_r+1
          sparse_dump_r(sp_count_r)= &
               sp_fact*factor_drygrid(ix,jy,kz)
        else ! factor is one
          sp_zer=.true.
        endif
      end do
    end do
  end do
  write(unitoutfactor) sp_count_i
  write(unitoutfactor) (sparse_dump_i(i),i=1,sp_count_i)
  write(unitoutfactor) sp_count_r
  write(unitoutfactor) (sparse_dump_r(i),i=1,sp_count_r)
  close(unitoutfactor)


  if (gridtotal.gt.0.) gridtotalunc=gridsigmatotal/gridtotal

  ! Dump of receptor concentrations

  if (numreceptor.gt.0 .and. (iout.eq.2 .or. iout.eq.3)  ) then
    write(unitoutreceptppt) itime
    do ks=1,nspec
      write(unitoutreceptppt) (1.e12*creceptor(i,ks)/outnum* &
           weightair/weightmolar(ks)/densityoutrecept(i),i=1,numreceptor)
    end do
  endif

  ! Dump of receptor concentrations

  if (numreceptor.gt.0) then
    write(unitoutrecept) itime
    do ks=1,nspec
      write(unitoutrecept) (1.e12*creceptor(i,ks)/outnum, &
           i=1,numreceptor)
    end do
  endif

  ! RLT Aug 2017
  ! Write out conversion factor for dry air
  if (numreceptor.gt.0) then
    inquire(file=path(2)(1:length(2))//'factor_dryreceptor',exist=lexist)
     if (lexist.and..not.lstart) then
     ! open and append
      open(unitoutfactor,file=path(2)(1:length(2))//'factor_dryreceptor',form='unformatted',&
              status='old',action='write',access='append')
    else
      ! create new
      open(unitoutfactor,file=path(2)(1:length(2))//'factor_dryreceptor',form='unformatted',&
              status='replace',action='write')
    endif
    write(unitoutfactor) itime
    write(unitoutfactor) (factor_dryrecept(i),i=1,numreceptor)
    close(unitoutfactor)
  endif

  ! reset lstart
  if (lstart) then
    lstart=.false.
  endif
  print*, 'after writing output files: lstart = ',lstart


  ! Reinitialization of grid
  !*************************

  do ks=1,nspec
    do kp=1,maxpointspec_act
      do i=1,numreceptor
        creceptor(i,ks)=0.
      end do
      do jy=0,numygrid-1
        do ix=0,numxgrid-1
          do l=1,nclassunc
            do nage=1,nageclass
              do kz=1,numzgrid
                gridunc(ix,jy,kz,ks,kp,l,nage)=0.
              end do
            end do
          end do
        end do
      end do
    end do
  end do
end subroutine concoutput_inversion

subroutine concoutput_inversion_nest(itime,outnum)
  !                        i     i
  !*****************************************************************************
  !                                                                            *
  !     Output of the concentration grid and the receptor concentrations.      *
  !                                                                            *
  !     Author: A. Stohl                                                       *
  !                                                                            *
  !     24 May 1995                                                            *
  !                                                                            *
  !     13 April 1999, Major update: if output size is smaller, dump output    *
  !                    in sparse matrix format; additional output of           *
  !                    uncertainty                                             *
  !                                                                            *
  !     05 April 2000, Major update: output of age classes; output for backward*
  !                    runs is time spent in grid cell times total mass of     *
  !                    species.                                                *
  !                                                                            *
  !     17 February 2002, Appropriate dimensions for backward and forward runs *
  !                       are now specified in file par_mod                    *
  !                                                                            *
  !     June 2006, write grid in sparse matrix with a single write command     *
  !                in order to save disk space                                 *
  !                                                                            *
  !     2008 new sparse matrix format                                          *
  !
  !     January 2017,  Separate files by release but include all timesteps     *
  !                                                                            *
  !*****************************************************************************
  !                                                                            *
  ! Variables:                                                                 *
  ! outnum          number of samples                                          *
  ! ncells          number of cells with non-zero concentrations               *
  ! sparse          .true. if in sparse matrix format, else .false.            *
  ! tot_mu          1 for forward, initial mass mixing ration for backw. runs  *
  !                                                                            *
  !*****************************************************************************

  use unc_mod
  use mean_mod

  implicit none

  real(kind=dp) :: jul
  integer :: itime,i,ix,jy,kz,ks,kp,l,iix,jjy,kzz,nage,jjjjmmdd,ihmmss
  integer :: sp_count_i,sp_count_r
  real :: sp_fact
  real :: outnum,densityoutrecept(maxreceptor),xl,yl
  ! RLT
  real :: densitydryrecept(maxreceptor)
  real :: factor_dryrecept(maxreceptor)

  real(dep_prec) :: auxgrid(nclassunc)
  real :: halfheight,dz,dz1,dz2,tot_mu(maxspec,maxpointspec_act)
  real,parameter :: smallnum = tiny(0.0) ! smallest number that can be handled
  real,parameter :: weightair=28.97
  logical :: sp_zer
  logical,save :: lnstart=.true.
  logical,save,allocatable,dimension(:) :: lnstartrel
  character :: adate*8,atime*6
  character(len=3) :: anspec
  logical :: lexist
  character :: areldate*8,areltime*6

  if(lnstart) then
    allocate(lnstartrel(maxpointspec_act))
    lnstartrel(:)=.true.
  endif
  print*, 'lnstartrel = ',lnstartrel

  ! Determine current calendar date, needed for the file name
  !**********************************************************

  jul=bdate+real(itime,kind=dp)/86400._dp
  call caldate(jul,jjjjmmdd,ihmmss)
  write(adate,'(i8.8)') jjjjmmdd
  write(atime,'(i6.6)') ihmmss

  print*, 'outnum:',outnum
  print*, 'datetime:',adate//atime

  ! For forward simulations, output fields have dimension MAXSPEC,
  ! for backward simulations, output fields have dimension MAXPOINT.
  ! Thus, make loops either about nspec, or about numpoint
  !*****************************************************************


    if (ldirect.eq.1) then
       do ks=1,nspec
         do kp=1,maxpointspec_act
           tot_mu(ks,kp)=1
         end do
       end do
   else
      do ks=1,nspec
             do kp=1,maxpointspec_act
               tot_mu(ks,kp)=xmass(kp,ks)
             end do
      end do
    endif


  !*******************************************************************
  ! Compute air density: sufficiently accurate to take it
  ! from coarse grid at some time
  ! Determine center altitude of output layer, and interpolate density
  ! data to that altitude
  !*******************************************************************

  do kz=1,numzgrid
    if (kz.eq.1) then
      halfheight=outheight(1)/2.
    else
      halfheight=(outheight(kz)+outheight(kz-1))/2.
    endif
    do kzz=2,nz
      if ((height(kzz-1).lt.halfheight).and. &
           (height(kzz).gt.halfheight)) goto 46
    end do
46   kzz=max(min(kzz,nz),2)
    dz1=halfheight-height(kzz-1)
    dz2=height(kzz)-halfheight
    dz=dz1+dz2
    do jy=0,numygridn-1
      do ix=0,numxgridn-1
        xl=outlon0n+real(ix)*dxoutn
        yl=outlat0n+real(jy)*dyoutn
        xl=(xl-xlon0)/dx
        yl=(yl-ylat0)/dy
        iix=max(min(nint(xl),nxmin1),0)
        jjy=max(min(nint(yl),nymin1),0)
        densityoutgrid(ix,jy,kz)=(rho(iix,jjy,kzz,2)*dz1+ &
             rho(iix,jjy,kzz-1,2)*dz2)/dz
  ! RLT
        densitydrygrid(ix,jy,kz)=(rho_dry(iix,jjy,kzz,2)*dz1+ &
             rho_dry(iix,jjy,kzz-1,2)*dz2)/dz
      end do
    end do
  end do

  do i=1,numreceptor
    xl=xreceptor(i)
    yl=yreceptor(i)
    iix=max(min(nint(xl),nxmin1),0)
    jjy=max(min(nint(yl),nymin1),0)
    densityoutrecept(i)=rho(iix,jjy,1,2)
  ! RLT
    densitydryrecept(i)=rho_dry(iix,jjy,1,2)
  end do

  ! RLT
  ! conversion factor for output relative to dry air
  factor_drygrid=densityoutgrid/densitydrygrid
  factor_dryrecept=densityoutrecept/densitydryrecept

  ! Output is different for forward and backward simulations
    do kz=1,numzgrid
      do jy=0,numygridn-1
        do ix=0,numxgridn-1
          if (ldirect.eq.1) then
            factor3d(ix,jy,kz)=1.e12/volumen(ix,jy,kz)/outnum
          else
            factor3d(ix,jy,kz)=real(abs(loutaver))/outnum
          endif
        end do
      end do
    end do

  !*********************************************************************
  ! Determine the standard deviation of the mean concentration or mixing
  ! ratio (uncertainty of the output) and the dry and wet deposition
  !*********************************************************************

  do ks=1,nspec

  write(anspec,'(i3.3)') ks

    do kp=1,maxpointspec_act

      print*, 'itime = ',itime
      print*, 'lage(1) = ',lage(1)
      print*, 'ireleasestart(kp) = ',ireleasestart(kp)
      print*, 'ireleaseend(kp) = ',ireleaseend(kp)

      ! check itime is within release and backward trajectory length
      if (nageclass.eq.1) then
        if ((itime.gt.ireleaseend(kp)).or.(itime.lt.(ireleasestart(kp)-lage(1)))) then
          go to 10
        endif
      endif

      ! calculate date of release
      jul=bdate+real(ireleasestart(kp),kind=dp)/86400._dp    ! this is the current day
      call caldate(jul,jjjjmmdd,ihmmss)
      write(areldate,'(i8.8)') jjjjmmdd
      write(areltime,'(i6.6)') ihmmss
      print*, areldate//areltime

      ! calculate date of field
      jul=bdate+real(itime,kind=dp)/86400._dp
      call caldate(jul,jjjjmmdd,ihmmss)
      write(adate,'(i8.8)') jjjjmmdd
      write(atime,'(i6.6)') ihmmss
      print*, adate//atime

      if ((iout.eq.1).or.(iout.eq.3).or.(iout.eq.5)) then
        if (ldirect.eq.1) then
          ! concentrations
          inquire(file=path(2)(1:length(2))//'grid_conc_nest_'//areldate// &
                  areltime//'_'//anspec,exist=lexist)
          if(lexist.and..not.lnstartrel(kp)) then
            ! open and append to existing file
            open(unitoutgrid,file=path(2)(1:length(2))//'grid_conc_nest_'//areldate// &
                 areltime//'_'//anspec,form='unformatted',status='old',action='write',access='append')
          else
            ! open new file
            open(unitoutgrid,file=path(2)(1:length(2))//'grid_conc_nest_'//areldate// &
                 areltime//'_'//anspec,form='unformatted',status='replace',action='write')
          endif
        else
          ! residence times
          inquire(file=path(2)(1:length(2))//'grid_time_nest_'//areldate// &
                  areltime//'_'//anspec,exist=lexist)
          if(lexist.and..not.lnstartrel(kp)) then
            ! open and append to existing file
            open(unitoutgrid,file=path(2)(1:length(2))//'grid_time_nest_'//areldate// &
                 areltime//'_'//anspec,form='unformatted',status='old',action='write',access='append')
          else
            ! open new file
            open(unitoutgrid,file=path(2)(1:length(2))//'grid_time_nest_'//areldate// &
                 areltime//'_'//anspec,form='unformatted',status='replace',action='write')
          endif
        endif
        write(unitoutgrid) jjjjmmdd
        write(unitoutgrid) ihmmss
      endif

      if ((iout.eq.2).or.(iout.eq.3)) then
        ! mixing ratio
        inquire(file=path(2)(1:length(2))//'grid_pptv_nest_'//areldate// &
                areltime//'_'//anspec,exist=lexist)
        if(lexist.and..not.lnstartrel(kp)) then
          ! open and append to existing file
          open(unitoutgridppt,file=path(2)(1:length(2))//'grid_pptv_nest_'//areldate// &
               areltime//'_'//anspec,form='unformatted',status='old',action='write',access='append')
        else
          ! open new file
          open(unitoutgridppt,file=path(2)(1:length(2))//'grid_pptv_nest_'//areldate// &
               areltime//'_'//anspec,form='unformatted',status='replace',action='write')
        endif
        write(unitoutgridppt) jjjjmmdd
        write(unitoutgridppt) ihmmss
      endif

      lnstartrel(kp)=.false.

      do nage=1,nageclass

        do jy=0,numygridn-1
          do ix=0,numxgridn-1

  ! CONCENTRATION OR MIXING RATIO
            do kz=1,numzgrid
              do l=1,nclassunc
                auxgrid(l)=griduncn(ix,jy,kz,ks,kp,l,nage)
              end do
              call mean(auxgrid,grid(ix,jy,kz), &
                     gridsigma(ix,jy,kz),nclassunc)
  ! Multiply by number of classes to get total concentration
              grid(ix,jy,kz)= &
                   grid(ix,jy,kz)*nclassunc
  ! Calculate standard deviation of the mean
              gridsigma(ix,jy,kz)= &
                   gridsigma(ix,jy,kz)* &
                   sqrt(real(nclassunc))
            end do
          end do
        end do


  !*******************************************************************
  ! Generate output: may be in concentration (ng/m3) or in mixing
  ! ratio (ppt) or both
  ! Output the position and the values alternated multiplied by
  ! 1 or -1, first line is number of values, number of positions
  ! For backward simulations, the unit is seconds, stored in grid_time
  !*******************************************************************

  ! Concentration output
  !*********************

        if ((iout.eq.1).or.(iout.eq.3).or.(iout.eq.5)) then

  ! Concentrations

  ! surf_only write only 1st layer 

         sp_count_i=0
         sp_count_r=0
         sp_fact=-1.
         sp_zer=.true.
          do kz=1,1
            do jy=0,numygridn-1
              do ix=0,numxgridn-1
                if (grid(ix,jy,kz).gt.smallnum) then
                  if (sp_zer.eqv..true.) then ! first non zero value
                    sp_count_i=sp_count_i+1
                    sparse_dump_i(sp_count_i)= &
                         ix+jy*numxgridn+kz*numxgridn*numygridn
                    sp_zer=.false.
                    sp_fact=sp_fact*(-1.)
                   endif
                   sp_count_r=sp_count_r+1
                   sparse_dump_r(sp_count_r)= &
                        sp_fact* &
                        grid(ix,jy,kz)* &
                        factor3d(ix,jy,kz)/tot_mu(ks,kp)
  !                 if ((factor(ix,jy,kz)/tot_mu(ks,kp)).eq.0)
  !    +              write (*,*) factor(ix,jy,kz),tot_mu(ks,kp),ks,kp
                   sparse_dump_u(sp_count_r)= &
                        gridsigma(ix,jy,kz)* &
                        factor3d(ix,jy,kz)/tot_mu(ks,kp)
              else ! concentration is zero
                  sp_zer=.true.
              endif
              end do
            end do
          end do
         write(unitoutgrid) sp_count_i
         write(unitoutgrid) (sparse_dump_i(i),i=1,sp_count_i)
         write(unitoutgrid) sp_count_r
         write(unitoutgrid) (sparse_dump_r(i),i=1,sp_count_r)

      endif !  concentration output

  ! Mixing ratio output
  !********************

      if ((iout.eq.2).or.(iout.eq.3)) then      ! mixing ratio


  ! Mixing ratios

    ! surf_only write only 1st layer 

         sp_count_i=0
         sp_count_r=0
         sp_fact=-1.
         sp_zer=.true.
          do kz=1,1
            do jy=0,numygridn-1
              do ix=0,numxgridn-1
                if (grid(ix,jy,kz).gt.smallnum) then
                  if (sp_zer.eqv..true.) then ! first non zero value
                    sp_count_i=sp_count_i+1
                    sparse_dump_i(sp_count_i)= &
                         ix+jy*numxgridn+kz*numxgridn*numygridn
                    sp_zer=.false.
                    sp_fact=sp_fact*(-1.)
                 endif
                 sp_count_r=sp_count_r+1
                 sparse_dump_r(sp_count_r)= &
                      sp_fact* &
                      1.e12*grid(ix,jy,kz) &
                      /volumen(ix,jy,kz)/outnum* &
                      weightair/weightmolar(ks)/densityoutgrid(ix,jy,kz)
                 sparse_dump_u(sp_count_r)= &
                      1.e12*gridsigma(ix,jy,kz)/volumen(ix,jy,kz)/ &
                      outnum*weightair/weightmolar(ks)/ &
                      densityoutgrid(ix,jy,kz)
              else ! concentration is zero
                  sp_zer=.true.
              endif
              end do
            end do
          end do
          write(unitoutgridppt) sp_count_i
          write(unitoutgridppt) (sparse_dump_i(i),i=1,sp_count_i)
          write(unitoutgridppt) sp_count_r
          write(unitoutgridppt) (sparse_dump_r(i),i=1,sp_count_r)

        endif ! output for ppt
 
      end do ! nageclass

      close(unitoutgridppt)
      close(unitoutgrid)

      ! itime is outside range
10    continue

    end do ! maxpointspec_act

  end do ! nspec


  ! RLT Aug 2017
  ! Write out conversion factor for dry air
  inquire(file=path(2)(1:length(2))//'factor_drygrid_nest',exist=lexist)
  if (lexist.and..not.lnstart) then
    ! open and append
    open(unitoutfactor,file=path(2)(1:length(2))//'factor_drygrid_nest',form='unformatted',&
            status='old',action='write',access='append')
  else
    ! create new
    open(unitoutfactor,file=path(2)(1:length(2))//'factor_drygrid_nest',form='unformatted',&
            status='replace',action='write')
  endif
  sp_count_i=0
  sp_count_r=0
  sp_fact=-1.
  sp_zer=.true.
  do kz=1,1
    do jy=0,numygridn-1
      do ix=0,numxgridn-1
        if (factor_drygrid(ix,jy,kz).gt.(1.+smallnum).or.factor_drygrid(ix,jy,kz).lt.(1.-smallnum)) then
          if (sp_zer.eqv..true.) then ! first value not equal to one
            sp_count_i=sp_count_i+1
            sparse_dump_i(sp_count_i)= &
                  ix+jy*numxgridn+kz*numxgridn*numygridn
            sp_zer=.false.
            sp_fact=sp_fact*(-1.)
          endif
          sp_count_r=sp_count_r+1
          sparse_dump_r(sp_count_r)= &
               sp_fact*factor_drygrid(ix,jy,kz)
        else ! factor is one
          sp_zer=.true.
        endif
      end do
    end do
  end do
  write(unitoutfactor) sp_count_i
  write(unitoutfactor) (sparse_dump_i(i),i=1,sp_count_i)
  write(unitoutfactor) sp_count_r
  write(unitoutfactor) (sparse_dump_r(i),i=1,sp_count_r)
  close(unitoutfactor)

  ! reset lnstart
  if (lnstart) then
    lnstart=.false.
  endif

  ! Reinitialization of grid
  !*************************

  do ks=1,nspec
  do kp=1,maxpointspec_act
    do i=1,numreceptor
      creceptor(i,ks)=0.
    end do
    do jy=0,numygridn-1
      do ix=0,numxgridn-1
        do l=1,nclassunc
          do nage=1,nageclass
            do kz=1,numzgrid
              griduncn(ix,jy,kz,ks,kp,l,nage)=0.
            end do
          end do
        end do
      end do
    end do
  end do
  end do
end subroutine concoutput_inversion_nest

subroutine concoutput_surf(itime,outnum,gridtotalunc,wetgridtotalunc, &
     drygridtotalunc)
  !                        i     i          o             o
  !       o
  !*****************************************************************************
  !                                                                            *
  !     Output of the concentration grid and the receptor concentrations.      *
  !                                                                            *
  !     Author: A. Stohl                                                       *
  !                                                                            *
  !     24 May 1995                                                            *
  !                                                                            *
  !     13 April 1999, Major update: if output size is smaller, dump output    *
  !                    in sparse matrix format; additional output of           *
  !                    uncertainty                                             *
  !                                                                            *
  !     05 April 2000, Major update: output of age classes; output for backward*
  !                    runs is time spent in grid cell times total mass of     *
  !                    species.                                                *
  !                                                                            *
  !     17 February 2002, Appropriate dimensions for backward and forward runs *
  !                       are now specified in file par_mod                    *
  !                                                                            *
  !     June 2006, write grid in sparse matrix with a single write command     *
  !                in order to save disk space                                 *
  !                                                                            *
  !     2008 new sparse matrix format                                          *
  !                                                                            *
  !*****************************************************************************
  !                                                                            *
  ! Variables:                                                                 *
  ! outnum          number of samples                                          *
  ! ncells          number of cells with non-zero concentrations               *
  ! sparse          .true. if in sparse matrix format, else .false.            *
  ! tot_mu          1 for forward, initial mass mixing ration for backw. runs  *
  !                                                                            *
  !*****************************************************************************

  use unc_mod
  use mean_mod

  implicit none

  real(kind=dp) :: jul
  integer :: itime,i,ix,jy,kz,ks,kp,l,iix,jjy,kzz,nage,jjjjmmdd,ihmmss
  integer :: sp_count_i,sp_count_r
  real :: sp_fact
  real :: outnum,densityoutrecept(maxreceptor),xl,yl
  ! RLT
  real :: densitydryrecept(maxreceptor)
  real :: factor_dryrecept(maxreceptor)

  real(dep_prec) :: auxgrid(nclassunc)
  real(sp) :: gridtotal,gridsigmatotal,gridtotalunc
  real(dep_prec) :: wetgridtotal,wetgridsigmatotal,wetgridtotalunc
  real(dep_prec) :: drygridtotal,drygridsigmatotal,drygridtotalunc
  real :: halfheight,dz,dz1,dz2,tot_mu(maxspec,maxpointspec_act)
  real,parameter :: smallnum = tiny(0.0) ! smallest number that can be handled
  real,parameter :: weightair=28.97
  logical :: sp_zer
  character :: adate*8,atime*6
  character(len=3) :: anspec
  logical :: lexist


  if (verbosity.eq.1) then
    print*,'inside concoutput_surf '
    CALL SYSTEM_CLOCK(count_clock)
    WRITE(*,*) 'SYSTEM_CLOCK',count_clock - count_clock0   
  endif

  ! Determine current calendar date, needed for the file name
  !**********************************************************

  jul=bdate+real(itime,kind=dp)/86400._dp
  call caldate(jul,jjjjmmdd,ihmmss)
  write(adate,'(i8.8)') jjjjmmdd
  write(atime,'(i6.6)') ihmmss

  open(unitdates,file=path(2)(1:length(2))//'dates', ACCESS='APPEND')
  write(unitdates,'(a)') adate//atime
  close(unitdates)

  ! For forward simulations, output fields have dimension MAXSPEC,
  ! for backward simulations, output fields have dimension MAXPOINT.
  ! Thus, make loops either about nspec, or about numpoint
  !*****************************************************************


  if (ldirect.eq.1) then
    do ks=1,nspec
      do kp=1,maxpointspec_act
        tot_mu(ks,kp)=1
      end do
    end do
  else
    do ks=1,nspec
      do kp=1,maxpointspec_act
        tot_mu(ks,kp)=xmass(kp,ks)
      end do
    end do
  endif


  if (verbosity.eq.1) then
    print*,'concoutput_surf 2'
    CALL SYSTEM_CLOCK(count_clock)
    WRITE(*,*) 'SYSTEM_CLOCK',count_clock - count_clock0   
  endif

  !*******************************************************************
  ! Compute air density: sufficiently accurate to take it
  ! from coarse grid at some time
  ! Determine center altitude of output layer, and interpolate density
  ! data to that altitude
  !*******************************************************************

  do kz=1,numzgrid
    if (kz.eq.1) then
      halfheight=outheight(1)/2.
    else
      halfheight=(outheight(kz)+outheight(kz-1))/2.
    endif
    do kzz=2,nz
      if ((height(kzz-1).lt.halfheight).and. &
           (height(kzz).gt.halfheight)) goto 46
    end do
46  kzz=max(min(kzz,nz),2)
    dz1=halfheight-height(kzz-1)
    dz2=height(kzz)-halfheight
    dz=dz1+dz2
    do jy=0,numygrid-1
      do ix=0,numxgrid-1
        xl=outlon0+real(ix)*dxout
        yl=outlat0+real(jy)*dyout
        xl=(xl-xlon0)/dx
        yl=(yl-ylat0)/dy
        iix=max(min(nint(xl),nxmin1),0)
        jjy=max(min(nint(yl),nymin1),0)
        densityoutgrid(ix,jy,kz)=(rho(iix,jjy,kzz,2)*dz1+ &
             rho(iix,jjy,kzz-1,2)*dz2)/dz
  ! RLT
        densitydrygrid(ix,jy,kz)=(rho_dry(iix,jjy,kzz,2)*dz1+ &
             rho_dry(iix,jjy,kzz-1,2)*dz2)/dz
      end do
    end do
  end do

  do i=1,numreceptor
    xl=xreceptor(i)
    yl=yreceptor(i)
    iix=max(min(nint(xl),nxmin1),0)
    jjy=max(min(nint(yl),nymin1),0)
    densityoutrecept(i)=rho(iix,jjy,1,2)
  ! RLT
    densitydryrecept(i)=rho_dry(iix,jjy,1,2)
  end do

  ! RLT
  ! conversion factor for output relative to dry air
  factor_drygrid=densityoutgrid/densitydrygrid
  factor_dryrecept=densityoutrecept/densitydryrecept

  ! Output is different for forward and backward simulations
  do kz=1,numzgrid
    do jy=0,numygrid-1
      do ix=0,numxgrid-1
        if (ldirect.eq.1) then
          factor3d(ix,jy,kz)=1.e12/volume(ix,jy,kz)/outnum
        else
          factor3d(ix,jy,kz)=real(abs(loutaver))/outnum
        endif
      end do
    end do
  end do

  !*********************************************************************
  ! Determine the standard deviation of the mean concentration or mixing
  ! ratio (uncertainty of the output) and the dry and wet deposition
  !*********************************************************************

  if (verbosity.eq.1) then
    print*,'concoutput_surf 3 (sd)'
    CALL SYSTEM_CLOCK(count_clock)
    WRITE(*,*) 'SYSTEM_CLOCK',count_clock - count_clock0   
  endif
  gridtotal=0.
  gridsigmatotal=0.
  gridtotalunc=0.
  wetgridtotal=0.
  wetgridsigmatotal=0.
  wetgridtotalunc=0.
  drygridtotal=0.
  drygridsigmatotal=0.
  drygridtotalunc=0.

  do ks=1,nspec

    write(anspec,'(i3.3)') ks
    if ((iout.eq.1).or.(iout.eq.3).or.(iout.eq.5)) then
      if (ldirect.eq.1) then
        open(unitoutgrid,file=path(2)(1:length(2))//'grid_conc_'//adate// &
             atime//'_'//anspec,form='unformatted')
      else
        open(unitoutgrid,file=path(2)(1:length(2))//'grid_time_'//adate// &
             atime//'_'//anspec,form='unformatted')
      endif
      write(unitoutgrid) itime
    endif

    if ((iout.eq.2).or.(iout.eq.3)) then      ! mixing ratio
      open(unitoutgridppt,file=path(2)(1:length(2))//'grid_pptv_'//adate// &
           atime//'_'//anspec,form='unformatted')

      write(unitoutgridppt) itime
    endif

    do kp=1,maxpointspec_act
      do nage=1,nageclass

        do jy=0,numygrid-1
          do ix=0,numxgrid-1

  ! WET DEPOSITION
            if ((WETDEP).and.(ldirect.gt.0)) then
              do l=1,nclassunc
                auxgrid(l)=wetgridunc(ix,jy,ks,kp,l,nage)
              end do
              call mean(auxgrid,wetgrid(ix,jy), &
                   wetgridsigma(ix,jy),nclassunc)
  ! Multiply by number of classes to get total concentration
              wetgrid(ix,jy)=wetgrid(ix,jy) &
                   *nclassunc
              wetgridtotal=wetgridtotal+wetgrid(ix,jy)
  ! Calculate standard deviation of the mean
              wetgridsigma(ix,jy)= &
                   wetgridsigma(ix,jy)* &
                   sqrt(real(nclassunc))
              wetgridsigmatotal=wetgridsigmatotal+ &
                   wetgridsigma(ix,jy)
            endif

  ! DRY DEPOSITION
            if ((DRYDEP).and.(ldirect.gt.0)) then
              do l=1,nclassunc
                auxgrid(l)=drygridunc(ix,jy,ks,kp,l,nage)
              end do
              call mean(auxgrid,drygrid(ix,jy), &
                   drygridsigma(ix,jy),nclassunc)
  ! Multiply by number of classes to get total concentration
              drygrid(ix,jy)=drygrid(ix,jy)* &
                   nclassunc
              drygridtotal=drygridtotal+drygrid(ix,jy)
  ! Calculate standard deviation of the mean
              drygridsigma(ix,jy)= &
                   drygridsigma(ix,jy)* &
                   sqrt(real(nclassunc))
125           drygridsigmatotal=drygridsigmatotal+ &
                   drygridsigma(ix,jy)
            endif

  ! CONCENTRATION OR MIXING RATIO
            do kz=1,numzgrid
              do l=1,nclassunc
                auxgrid(l)=gridunc(ix,jy,kz,ks,kp,l,nage)
              end do
              call mean(auxgrid,grid(ix,jy,kz), &
                   gridsigma(ix,jy,kz),nclassunc)
  ! Multiply by number of classes to get total concentration
              grid(ix,jy,kz)= &
                   grid(ix,jy,kz)*nclassunc
              gridtotal=gridtotal+grid(ix,jy,kz)
  ! Calculate standard deviation of the mean
              gridsigma(ix,jy,kz)= &
                   gridsigma(ix,jy,kz)* &
                   sqrt(real(nclassunc))
              gridsigmatotal=gridsigmatotal+ &
                   gridsigma(ix,jy,kz)
            end do
          end do
        end do


  !*******************************************************************
  ! Generate output: may be in concentration (ng/m3) or in mixing
  ! ratio (ppt) or both
  ! Output the position and the values alternated multiplied by
  ! 1 or -1, first line is number of values, number of positions
  ! For backward simulations, the unit is seconds, stored in grid_time
  !*******************************************************************

        if (verbosity.eq.1) then
          print*,'concoutput_surf 4 (output)'
          CALL SYSTEM_CLOCK(count_clock)
          WRITE(*,*) 'SYSTEM_CLOCK',count_clock - count_clock0   
        endif

  ! Concentration output
  !*********************

        if ((iout.eq.1).or.(iout.eq.3).or.(iout.eq.5)) then

          if (verbosity.eq.1) then
            print*,'concoutput_surf (Wet deposition)'
            CALL SYSTEM_CLOCK(count_clock)
            WRITE(*,*) 'SYSTEM_CLOCK',count_clock - count_clock0   
          endif

  ! Wet deposition
          sp_count_i=0
          sp_count_r=0
          sp_fact=-1.
          sp_zer=.true.
          if ((ldirect.eq.1).and.(WETDEP)) then
            do jy=0,numygrid-1
              do ix=0,numxgrid-1
  ! concentraion greater zero
                if (wetgrid(ix,jy).gt.smallnum) then
                  if (sp_zer.eqv..true.) then ! first non zero value
                    sp_count_i=sp_count_i+1
                    sparse_dump_i(sp_count_i)=ix+jy*numxgrid
                    sp_zer=.false.
                    sp_fact=sp_fact*(-1.)
                  endif
                  sp_count_r=sp_count_r+1
                  sparse_dump_r(sp_count_r)= &
                       sp_fact*1.e12*wetgrid(ix,jy)/area(ix,jy)
                  sparse_dump_u(sp_count_r)= &
                       1.e12*wetgridsigma(ix,jy)/area(ix,jy)
                else ! concentration is zero
                  sp_zer=.true.
                endif
              end do
            end do
          else
            sp_count_i=0
            sp_count_r=0
          endif
          write(unitoutgrid) sp_count_i
          write(unitoutgrid) (sparse_dump_i(i),i=1,sp_count_i)
          write(unitoutgrid) sp_count_r
          write(unitoutgrid) (sparse_dump_r(i),i=1,sp_count_r)

          if (verbosity.eq.1) then
            print*,'concoutput_surf (Dry deposition)'
            CALL SYSTEM_CLOCK(count_clock)
            WRITE(*,*) 'SYSTEM_CLOCK',count_clock - count_clock0   
          endif
  ! Dry deposition
          sp_count_i=0
          sp_count_r=0
          sp_fact=-1.
          sp_zer=.true.
          if ((ldirect.eq.1).and.(DRYDEP)) then
            do jy=0,numygrid-1
              do ix=0,numxgrid-1
                if (drygrid(ix,jy).gt.smallnum) then
                  if (sp_zer.eqv..true.) then ! first non zero value
                    sp_count_i=sp_count_i+1
                    sparse_dump_i(sp_count_i)=ix+jy*numxgrid
                    sp_zer=.false.
                    sp_fact=sp_fact*(-1.)
                  endif
                  sp_count_r=sp_count_r+1
                  sparse_dump_r(sp_count_r)= &
                       sp_fact* &
                       1.e12*drygrid(ix,jy)/area(ix,jy)
                  sparse_dump_u(sp_count_r)= &
                       1.e12*drygridsigma(ix,jy)/area(ix,jy)
                else ! concentration is zero
                  sp_zer=.true.
                endif
              end do
            end do
          else
            sp_count_i=0
            sp_count_r=0
          endif
          write(unitoutgrid) sp_count_i
          write(unitoutgrid) (sparse_dump_i(i),i=1,sp_count_i)
          write(unitoutgrid) sp_count_r
          write(unitoutgrid) (sparse_dump_r(i),i=1,sp_count_r)

          if (verbosity.eq.1) then
            print*,'concoutput_surf (Concentrations)'
            CALL SYSTEM_CLOCK(count_clock)
            WRITE(*,*) 'SYSTEM_CLOCK',count_clock - count_clock0   
          endif

  ! Concentrations

  ! surf_only write only 1st layer 

          sp_count_i=0
          sp_count_r=0
          sp_fact=-1.
          sp_zer=.true.
          do kz=1,1
            do jy=0,numygrid-1
              do ix=0,numxgrid-1
                if (grid(ix,jy,kz).gt.smallnum) then
                  if (sp_zer.eqv..true.) then ! first non zero value
                    sp_count_i=sp_count_i+1
                    sparse_dump_i(sp_count_i)= &
                         ix+jy*numxgrid+kz*numxgrid*numygrid
                    sp_zer=.false.
                    sp_fact=sp_fact*(-1.)
                  endif
                  sp_count_r=sp_count_r+1
                  sparse_dump_r(sp_count_r)= &
                       sp_fact* &
                       grid(ix,jy,kz)* &
                       factor3d(ix,jy,kz)/tot_mu(ks,kp)
                  sparse_dump_u(sp_count_r)= &
                       gridsigma(ix,jy,kz)* &
                       factor3d(ix,jy,kz)/tot_mu(ks,kp)
                else ! concentration is zero
                  sp_zer=.true.
                endif
              end do
            end do
          end do
          write(unitoutgrid) sp_count_i
          write(unitoutgrid) (sparse_dump_i(i),i=1,sp_count_i)
          write(unitoutgrid) sp_count_r
          write(unitoutgrid) (sparse_dump_r(i),i=1,sp_count_r)
        endif !  concentration output

  ! Mixing ratio output
  !********************

        if ((iout.eq.2).or.(iout.eq.3)) then      ! mixing ratio

  ! Wet deposition
          sp_count_i=0
          sp_count_r=0
          sp_fact=-1.
          sp_zer=.true.
          if ((ldirect.eq.1).and.(WETDEP)) then
            do jy=0,numygrid-1
              do ix=0,numxgrid-1
                if (wetgrid(ix,jy).gt.smallnum) then
                  if (sp_zer.eqv..true.) then ! first non zero value
                    sp_count_i=sp_count_i+1
                    sparse_dump_i(sp_count_i)= &
                         ix+jy*numxgrid
                    sp_zer=.false.
                    sp_fact=sp_fact*(-1.)
                  endif
                  sp_count_r=sp_count_r+1
                  sparse_dump_r(sp_count_r)= &
                       sp_fact* &
                       1.e12*wetgrid(ix,jy)/area(ix,jy)
                  sparse_dump_u(sp_count_r)= &
                       1.e12*wetgridsigma(ix,jy)/area(ix,jy)
                else ! concentration is zero
                  sp_zer=.true.
                endif
              end do
            end do
          else
            sp_count_i=0
            sp_count_r=0
          endif
          write(unitoutgridppt) sp_count_i
          write(unitoutgridppt) (sparse_dump_i(i),i=1,sp_count_i)
          write(unitoutgridppt) sp_count_r
          write(unitoutgridppt) (sparse_dump_r(i),i=1,sp_count_r)

  ! Dry deposition
          sp_count_i=0
          sp_count_r=0
          sp_fact=-1.
          sp_zer=.true.
          if ((ldirect.eq.1).and.(DRYDEP)) then
            do jy=0,numygrid-1
              do ix=0,numxgrid-1
                if (drygrid(ix,jy).gt.smallnum) then
                  if (sp_zer.eqv..true.) then ! first non zero value
                    sp_count_i=sp_count_i+1
                    sparse_dump_i(sp_count_i)= &
                         ix+jy*numxgrid
                    sp_zer=.false.
                    sp_fact=sp_fact*(-1)
                  endif
                  sp_count_r=sp_count_r+1
                  sparse_dump_r(sp_count_r)= &
                       sp_fact* &
                       1.e12*drygrid(ix,jy)/area(ix,jy)
                  sparse_dump_u(sp_count_r)= &
                       1.e12*drygridsigma(ix,jy)/area(ix,jy)
                else ! concentration is zero
                  sp_zer=.true.
                endif
              end do
            end do
          else
            sp_count_i=0
            sp_count_r=0
          endif
          write(unitoutgridppt) sp_count_i
          write(unitoutgridppt) (sparse_dump_i(i),i=1,sp_count_i)
          write(unitoutgridppt) sp_count_r
          write(unitoutgridppt) (sparse_dump_r(i),i=1,sp_count_r)

  ! Mixing ratios

  ! surf_only write only 1st layer 

          sp_count_i=0
          sp_count_r=0
          sp_fact=-1.
          sp_zer=.true.
          do kz=1,1
            do jy=0,numygrid-1
              do ix=0,numxgrid-1
                if (grid(ix,jy,kz).gt.smallnum) then
                  if (sp_zer.eqv..true.) then ! first non zero value
                    sp_count_i=sp_count_i+1
                    sparse_dump_i(sp_count_i)= &
                         ix+jy*numxgrid+kz*numxgrid*numygrid
                    sp_zer=.false.
                    sp_fact=sp_fact*(-1.)
                  endif
                  sp_count_r=sp_count_r+1
                  sparse_dump_r(sp_count_r)= &
                       sp_fact* &
                       1.e12*grid(ix,jy,kz) &
                       /volume(ix,jy,kz)/outnum* &
                       weightair/weightmolar(ks)/densityoutgrid(ix,jy,kz)
                  sparse_dump_u(sp_count_r)= &
                       1.e12*gridsigma(ix,jy,kz)/volume(ix,jy,kz)/ &
                       outnum*weightair/weightmolar(ks)/ &
                       densityoutgrid(ix,jy,kz)
                else ! concentration is zero
                  sp_zer=.true.
                endif
              end do
            end do
          end do
          write(unitoutgridppt) sp_count_i
          write(unitoutgridppt) (sparse_dump_i(i),i=1,sp_count_i)
          write(unitoutgridppt) sp_count_r
          write(unitoutgridppt) (sparse_dump_r(i),i=1,sp_count_r)
        endif ! output for ppt

      end do
    end do

    close(unitoutgridppt)
    close(unitoutgrid)

  end do

  ! RLT Aug 2017
  ! Write out conversion factor for dry air
  inquire(file=path(2)(1:length(2))//'factor_drygrid',exist=lexist)
  if (lexist) then
    ! open and append
    open(unitoutfactor,file=path(2)(1:length(2))//'factor_drygrid',form='unformatted',&
            status='old',action='write',access='append')
  else
    ! create new
    open(unitoutfactor,file=path(2)(1:length(2))//'factor_drygrid',form='unformatted',&
            status='new',action='write')
  endif
  sp_count_i=0
  sp_count_r=0
  sp_fact=-1.
  sp_zer=.true.
  do kz=1,1
    do jy=0,numygrid-1
      do ix=0,numxgrid-1
        if (factor_drygrid(ix,jy,kz).gt.(1.+smallnum).or.factor_drygrid(ix,jy,kz).lt.(1.-smallnum)) then
          if (sp_zer.eqv..true.) then ! first value not equal to one
            sp_count_i=sp_count_i+1
            sparse_dump_i(sp_count_i)= &
                  ix+jy*numxgrid+kz*numxgrid*numygrid
            sp_zer=.false.
            sp_fact=sp_fact*(-1.)
          endif
          sp_count_r=sp_count_r+1
          sparse_dump_r(sp_count_r)= &
               sp_fact*factor_drygrid(ix,jy,kz)
        else ! factor is one
          sp_zer=.true.
        endif
      end do
    end do
  end do
  write(unitoutfactor) sp_count_i
  write(unitoutfactor) (sparse_dump_i(i),i=1,sp_count_i)
  write(unitoutfactor) sp_count_r
  write(unitoutfactor) (sparse_dump_r(i),i=1,sp_count_r)
  close(unitoutfactor)


  if (gridtotal.gt.0.) gridtotalunc=gridsigmatotal/gridtotal
  if (wetgridtotal.gt.0.) wetgridtotalunc=wetgridsigmatotal/ &
       wetgridtotal
  if (drygridtotal.gt.0.) drygridtotalunc=drygridsigmatotal/ &
       drygridtotal

  ! Dump of receptor concentrations

  if (numreceptor.gt.0 .and. (iout.eq.2 .or. iout.eq.3)  ) then
    write(unitoutreceptppt) itime
    do ks=1,nspec
      write(unitoutreceptppt) (1.e12*creceptor(i,ks)/outnum* &
           weightair/weightmolar(ks)/densityoutrecept(i),i=1,numreceptor)
    end do
  endif

  ! Dump of receptor concentrations

  if (numreceptor.gt.0) then
    write(unitoutrecept) itime
    do ks=1,nspec
      write(unitoutrecept) (1.e12*creceptor(i,ks)/outnum, &
           i=1,numreceptor)
    end do
  endif

  ! RLT Aug 2017
  ! Write out conversion factor for dry air
  if (numreceptor.gt.0) then
    inquire(file=path(2)(1:length(2))//'factor_dryreceptor',exist=lexist)
     if (lexist) then
     ! open and append
      open(unitoutfactor,file=path(2)(1:length(2))//'factor_dryreceptor',form='unformatted',&
              status='old',action='write',access='append')
    else
      ! create new
      open(unitoutfactor,file=path(2)(1:length(2))//'factor_dryreceptor',form='unformatted',&
              status='new',action='write')
    endif
    write(unitoutfactor) itime
    write(unitoutfactor) (factor_dryrecept(i),i=1,numreceptor)
    close(unitoutfactor)
  endif

  ! Reinitialization of grid
  !*************************

  do ks=1,nspec
    do kp=1,maxpointspec_act
      do i=1,numreceptor
        creceptor(i,ks)=0.
      end do
      do jy=0,numygrid-1
        do ix=0,numxgrid-1
          do l=1,nclassunc
            do nage=1,nageclass
              do kz=1,numzgrid
                gridunc(ix,jy,kz,ks,kp,l,nage)=0.
              end do
            end do
          end do
        end do
      end do
    end do
  end do
end subroutine concoutput_surf

subroutine concoutput_surf_nest(itime,outnum)
  !                        i     i
  !*****************************************************************************
  !                                                                            *
  !     Output of the concentration grid and the receptor concentrations.      *
  !                                                                            *
  !     Author: A. Stohl                                                       *
  !                                                                            *
  !     24 May 1995                                                            *
  !                                                                            *
  !     13 April 1999, Major update: if output size is smaller, dump output    *
  !                    in sparse matrix format; additional output of           *
  !                    uncertainty                                             *
  !                                                                            *
  !     05 April 2000, Major update: output of age classes; output for backward*
  !                    runs is time spent in grid cell times total mass of     *
  !                    species.                                                *
  !                                                                            *
  !     17 February 2002, Appropriate dimensions for backward and forward runs *
  !                       are now specified in file par_mod                    *
  !                                                                            *
  !     June 2006, write grid in sparse matrix with a single write command     *
  !                in order to save disk space                                 *
  !                                                                            *
  !     2008 new sparse matrix format                                          *
  !                                                                            *
  !*****************************************************************************
  !                                                                            *
  ! Variables:                                                                 *
  ! outnum          number of samples                                          *
  ! ncells          number of cells with non-zero concentrations               *
  ! sparse          .true. if in sparse matrix format, else .false.            *
  ! tot_mu          1 for forward, initial mass mixing ration for backw. runs  *
  !                                                                            *
  !*****************************************************************************

  use unc_mod
  use mean_mod

  implicit none

  real(kind=dp) :: jul
  integer :: itime,i,ix,jy,kz,ks,kp,l,iix,jjy,kzz,nage,jjjjmmdd,ihmmss
  integer :: sp_count_i,sp_count_r
  real :: sp_fact
  real :: outnum,densityoutrecept(maxreceptor),xl,yl
  ! RLT
  real :: densitydryrecept(maxreceptor)
  real :: factor_dryrecept(maxreceptor)

  real(dep_prec) :: auxgrid(nclassunc)
  real :: halfheight,dz,dz1,dz2,tot_mu(maxspec,maxpointspec_act)
  real,parameter :: smallnum = tiny(0.0) ! smallest number that can be handled
  real,parameter :: weightair=28.97
  logical :: sp_zer
  character :: adate*8,atime*6
  character(len=3) :: anspec
  logical :: lexist

  ! Determine current calendar date, needed for the file name
  !**********************************************************

  jul=bdate+real(itime,kind=dp)/86400._dp
  call caldate(jul,jjjjmmdd,ihmmss)
  write(adate,'(i8.8)') jjjjmmdd
  write(atime,'(i6.6)') ihmmss


  ! For forward simulations, output fields have dimension MAXSPEC,
  ! for backward simulations, output fields have dimension MAXPOINT.
  ! Thus, make loops either about nspec, or about numpoint
  !*****************************************************************


    if (ldirect.eq.1) then
       do ks=1,nspec
         do kp=1,maxpointspec_act
           tot_mu(ks,kp)=1
         end do
       end do
   else
      do ks=1,nspec
             do kp=1,maxpointspec_act
               tot_mu(ks,kp)=xmass(kp,ks)
             end do
      end do
    endif


  !*******************************************************************
  ! Compute air density: sufficiently accurate to take it
  ! from coarse grid at some time
  ! Determine center altitude of output layer, and interpolate density
  ! data to that altitude
  !*******************************************************************

  do kz=1,numzgrid
    if (kz.eq.1) then
      halfheight=outheight(1)/2.
    else
      halfheight=(outheight(kz)+outheight(kz-1))/2.
    endif
    do kzz=2,nz
      if ((height(kzz-1).lt.halfheight).and. &
           (height(kzz).gt.halfheight)) goto 46
    end do
46   kzz=max(min(kzz,nz),2)
    dz1=halfheight-height(kzz-1)
    dz2=height(kzz)-halfheight
    dz=dz1+dz2
    do jy=0,numygridn-1
      do ix=0,numxgridn-1
        xl=outlon0n+real(ix)*dxoutn
        yl=outlat0n+real(jy)*dyoutn
        xl=(xl-xlon0)/dx
        yl=(yl-ylat0)/dy
        iix=max(min(nint(xl),nxmin1),0)
        jjy=max(min(nint(yl),nymin1),0)
        densityoutgrid(ix,jy,kz)=(rho(iix,jjy,kzz,2)*dz1+ &
             rho(iix,jjy,kzz-1,2)*dz2)/dz
  ! RLT
        densitydrygrid(ix,jy,kz)=(rho_dry(iix,jjy,kzz,2)*dz1+ &
             rho_dry(iix,jjy,kzz-1,2)*dz2)/dz
      end do
    end do
  end do

    do i=1,numreceptor
      xl=xreceptor(i)
      yl=yreceptor(i)
      iix=max(min(nint(xl),nxmin1),0)
      jjy=max(min(nint(yl),nymin1),0)
      densityoutrecept(i)=rho(iix,jjy,1,2)
  ! RLT
    densitydryrecept(i)=rho_dry(iix,jjy,1,2)
    end do

  ! RLT
  ! conversion factor for output relative to dry air
  factor_drygrid=densityoutgrid/densitydrygrid
  factor_dryrecept=densityoutrecept/densitydryrecept

  ! Output is different for forward and backward simulations
    do kz=1,numzgrid
      do jy=0,numygridn-1
        do ix=0,numxgridn-1
          if (ldirect.eq.1) then
            factor3d(ix,jy,kz)=1.e12/volumen(ix,jy,kz)/outnum
          else
            factor3d(ix,jy,kz)=real(abs(loutaver))/outnum
          endif
        end do
      end do
    end do

  !*********************************************************************
  ! Determine the standard deviation of the mean concentration or mixing
  ! ratio (uncertainty of the output) and the dry and wet deposition
  !*********************************************************************

  do ks=1,nspec

  write(anspec,'(i3.3)') ks
  if ((iout.eq.1).or.(iout.eq.3).or.(iout.eq.5)) then
    if (ldirect.eq.1) then
      open(unitoutgrid,file=path(2)(1:length(2))//'grid_conc_nest_' &
           //adate// &
           atime//'_'//anspec,form='unformatted')
    else
      open(unitoutgrid,file=path(2)(1:length(2))//'grid_time_nest_' &
           //adate// &
           atime//'_'//anspec,form='unformatted')
    endif
     write(unitoutgrid) itime
   endif

  if ((iout.eq.2).or.(iout.eq.3)) then      ! mixing ratio
   open(unitoutgridppt,file=path(2)(1:length(2))//'grid_pptv_nest_' &
        //adate// &
        atime//'_'//anspec,form='unformatted')

    write(unitoutgridppt) itime
  endif

  do kp=1,maxpointspec_act
  do nage=1,nageclass

    do jy=0,numygridn-1
      do ix=0,numxgridn-1

  ! WET DEPOSITION
        if ((WETDEP).and.(ldirect.gt.0)) then
            do l=1,nclassunc
              auxgrid(l)=wetgriduncn(ix,jy,ks,kp,l,nage)
            end do
            call mean(auxgrid,wetgrid(ix,jy), &
                 wetgridsigma(ix,jy),nclassunc)
  ! Multiply by number of classes to get total concentration
            wetgrid(ix,jy)=wetgrid(ix,jy) &
                 *nclassunc
  ! Calculate standard deviation of the mean
            wetgridsigma(ix,jy)= &
                 wetgridsigma(ix,jy)* &
                 sqrt(real(nclassunc))
        endif

  ! DRY DEPOSITION
        if ((DRYDEP).and.(ldirect.gt.0)) then
            do l=1,nclassunc
              auxgrid(l)=drygriduncn(ix,jy,ks,kp,l,nage)
            end do
            call mean(auxgrid,drygrid(ix,jy), &
                 drygridsigma(ix,jy),nclassunc)
  ! Multiply by number of classes to get total concentration
            drygrid(ix,jy)=drygrid(ix,jy)* &
                 nclassunc
  ! Calculate standard deviation of the mean
            drygridsigma(ix,jy)= &
                 drygridsigma(ix,jy)* &
                 sqrt(real(nclassunc))
        endif

  ! CONCENTRATION OR MIXING RATIO
        do kz=1,numzgrid
            do l=1,nclassunc
              auxgrid(l)=griduncn(ix,jy,kz,ks,kp,l,nage)
            end do
            call mean(auxgrid,grid(ix,jy,kz), &
                 gridsigma(ix,jy,kz),nclassunc)
  ! Multiply by number of classes to get total concentration
            grid(ix,jy,kz)= &
                 grid(ix,jy,kz)*nclassunc
  ! Calculate standard deviation of the mean
            gridsigma(ix,jy,kz)= &
                 gridsigma(ix,jy,kz)* &
                 sqrt(real(nclassunc))
        end do
      end do
    end do


  !*******************************************************************
  ! Generate output: may be in concentration (ng/m3) or in mixing
  ! ratio (ppt) or both
  ! Output the position and the values alternated multiplied by
  ! 1 or -1, first line is number of values, number of positions
  ! For backward simulations, the unit is seconds, stored in grid_time
  !*******************************************************************

  ! Concentration output
  !*********************
  if ((iout.eq.1).or.(iout.eq.3).or.(iout.eq.5)) then

  ! Wet deposition
         sp_count_i=0
         sp_count_r=0
         sp_fact=-1.
         sp_zer=.true.
         if ((ldirect.eq.1).and.(WETDEP)) then
         do jy=0,numygridn-1
            do ix=0,numxgridn-1
  !oncentraion greater zero
              if (wetgrid(ix,jy).gt.smallnum) then
                 if (sp_zer.eqv..true.) then ! first non zero value
                    sp_count_i=sp_count_i+1
                    sparse_dump_i(sp_count_i)=ix+jy*numxgridn
                    sp_zer=.false.
                    sp_fact=sp_fact*(-1.)
                 endif
                 sp_count_r=sp_count_r+1
                 sparse_dump_r(sp_count_r)= &
                      sp_fact*1.e12*wetgrid(ix,jy)/arean(ix,jy)
                 sparse_dump_u(sp_count_r)= &
                      1.e12*wetgridsigma(ix,jy)/area(ix,jy)
              else ! concentration is zero
                  sp_zer=.true.
              endif
            end do
         end do
         else
            sp_count_i=0
            sp_count_r=0
         endif
         write(unitoutgrid) sp_count_i
         write(unitoutgrid) (sparse_dump_i(i),i=1,sp_count_i)
         write(unitoutgrid) sp_count_r
         write(unitoutgrid) (sparse_dump_r(i),i=1,sp_count_r)
  ! Dry deposition
         sp_count_i=0
         sp_count_r=0
         sp_fact=-1.
         sp_zer=.true.
         if ((ldirect.eq.1).and.(DRYDEP)) then
          do jy=0,numygridn-1
            do ix=0,numxgridn-1
              if (drygrid(ix,jy).gt.smallnum) then
                 if (sp_zer.eqv..true.) then ! first non zero value
                    sp_count_i=sp_count_i+1
                    sparse_dump_i(sp_count_i)=ix+jy*numxgridn
                    sp_zer=.false.
                    sp_fact=sp_fact*(-1.)
                 endif
                 sp_count_r=sp_count_r+1
                 sparse_dump_r(sp_count_r)= &
                      sp_fact* &
                      1.e12*drygrid(ix,jy)/arean(ix,jy)
                 sparse_dump_u(sp_count_r)= &
                      1.e12*drygridsigma(ix,jy)/area(ix,jy)
              else ! concentration is zero
                  sp_zer=.true.
              endif
            end do
          end do
         else
            sp_count_i=0
            sp_count_r=0
         endif
         write(unitoutgrid) sp_count_i
         write(unitoutgrid) (sparse_dump_i(i),i=1,sp_count_i)
         write(unitoutgrid) sp_count_r
         write(unitoutgrid) (sparse_dump_r(i),i=1,sp_count_r)

  ! Concentrations

  ! if surf_only write only 1st layer 

         if(surf_only.eq.1) then
         sp_count_i=0
         sp_count_r=0
         sp_fact=-1.
         sp_zer=.true.
          do kz=1,1
            do jy=0,numygridn-1
              do ix=0,numxgridn-1
                if (grid(ix,jy,kz).gt.smallnum) then
                  if (sp_zer.eqv..true.) then ! first non zero value
                    sp_count_i=sp_count_i+1
                    sparse_dump_i(sp_count_i)= &
                         ix+jy*numxgridn+kz*numxgridn*numygridn
                    sp_zer=.false.
                    sp_fact=sp_fact*(-1.)
                   endif
                   sp_count_r=sp_count_r+1
                   sparse_dump_r(sp_count_r)= &
                        sp_fact* &
                        grid(ix,jy,kz)* &
                        factor3d(ix,jy,kz)/tot_mu(ks,kp)
  !                 if ((factor(ix,jy,kz)/tot_mu(ks,kp)).eq.0)
  !    +              write (*,*) factor(ix,jy,kz),tot_mu(ks,kp),ks,kp
                   sparse_dump_u(sp_count_r)= &
                        gridsigma(ix,jy,kz)* &
                        factor3d(ix,jy,kz)/tot_mu(ks,kp)
              else ! concentration is zero
                  sp_zer=.true.
              endif
              end do
            end do
          end do
         write(unitoutgrid) sp_count_i
         write(unitoutgrid) (sparse_dump_i(i),i=1,sp_count_i)
         write(unitoutgrid) sp_count_r
         write(unitoutgrid) (sparse_dump_r(i),i=1,sp_count_r)
         else

  ! write full vertical resolution

         sp_count_i=0
         sp_count_r=0
         sp_fact=-1.
         sp_zer=.true.
          do kz=1,numzgrid
            do jy=0,numygridn-1
              do ix=0,numxgridn-1
                if (grid(ix,jy,kz).gt.smallnum) then
                  if (sp_zer.eqv..true.) then ! first non zero value
                    sp_count_i=sp_count_i+1
                    sparse_dump_i(sp_count_i)= &
                         ix+jy*numxgridn+kz*numxgridn*numygridn
                    sp_zer=.false.
                    sp_fact=sp_fact*(-1.)
                   endif
                   sp_count_r=sp_count_r+1
                   sparse_dump_r(sp_count_r)= &
                        sp_fact* &
                        grid(ix,jy,kz)* &
                        factor3d(ix,jy,kz)/tot_mu(ks,kp)
                   sparse_dump_u(sp_count_r)= &
                        gridsigma(ix,jy,kz)* &
                        factor3d(ix,jy,kz)/tot_mu(ks,kp)
              else ! concentration is zero
                  sp_zer=.true.
              endif
              end do
            end do
          end do
         write(unitoutgrid) sp_count_i
         write(unitoutgrid) (sparse_dump_i(i),i=1,sp_count_i)
         write(unitoutgrid) sp_count_r
         write(unitoutgrid) (sparse_dump_r(i),i=1,sp_count_r)
         endif ! surf_only


    endif !  concentration output

  ! Mixing ratio output
  !********************

  if ((iout.eq.2).or.(iout.eq.3)) then      ! mixing ratio

  ! Wet deposition
         sp_count_i=0
         sp_count_r=0
         sp_fact=-1.
         sp_zer=.true.
         if ((ldirect.eq.1).and.(WETDEP)) then
          do jy=0,numygridn-1
            do ix=0,numxgridn-1
                if (wetgrid(ix,jy).gt.smallnum) then
                  if (sp_zer.eqv..true.) then ! first non zero value
                    sp_count_i=sp_count_i+1
                    sparse_dump_i(sp_count_i)= &
                         ix+jy*numxgridn
                    sp_zer=.false.
                    sp_fact=sp_fact*(-1.)
                 endif
                 sp_count_r=sp_count_r+1
                 sparse_dump_r(sp_count_r)= &
                      sp_fact* &
                      1.e12*wetgrid(ix,jy)/arean(ix,jy)
                 sparse_dump_u(sp_count_r)= &
                      1.e12*wetgridsigma(ix,jy)/area(ix,jy)
              else ! concentration is zero
                  sp_zer=.true.
              endif
            end do
          end do
         else
           sp_count_i=0
           sp_count_r=0
         endif
         write(unitoutgridppt) sp_count_i
         write(unitoutgridppt) (sparse_dump_i(i),i=1,sp_count_i)
         write(unitoutgridppt) sp_count_r
         write(unitoutgridppt) (sparse_dump_r(i),i=1,sp_count_r)

  ! Dry deposition
         sp_count_i=0
         sp_count_r=0
         sp_fact=-1.
         sp_zer=.true.
         if ((ldirect.eq.1).and.(DRYDEP)) then
          do jy=0,numygridn-1
            do ix=0,numxgridn-1
                if (drygrid(ix,jy).gt.smallnum) then
                  if (sp_zer.eqv..true.) then ! first non zero value
                    sp_count_i=sp_count_i+1
                    sparse_dump_i(sp_count_i)= &
                         ix+jy*numxgridn
                    sp_zer=.false.
                    sp_fact=sp_fact*(-1)
                 endif
                 sp_count_r=sp_count_r+1
                 sparse_dump_r(sp_count_r)= &
                      sp_fact* &
                      1.e12*drygrid(ix,jy)/arean(ix,jy)
                 sparse_dump_u(sp_count_r)= &
                      1.e12*drygridsigma(ix,jy)/area(ix,jy)
              else ! concentration is zero
                  sp_zer=.true.
              endif
            end do
          end do
         else
           sp_count_i=0
           sp_count_r=0
         endif
         write(unitoutgridppt) sp_count_i
         write(unitoutgridppt) (sparse_dump_i(i),i=1,sp_count_i)
         write(unitoutgridppt) sp_count_r
         write(unitoutgridppt) (sparse_dump_r(i),i=1,sp_count_r)

  ! Mixing ratios

    ! if surf_only write only 1st layer 

         if(surf_only.eq.1) then
         sp_count_i=0
         sp_count_r=0
         sp_fact=-1.
         sp_zer=.true.
          do kz=1,1
            do jy=0,numygridn-1
              do ix=0,numxgridn-1
                if (grid(ix,jy,kz).gt.smallnum) then
                  if (sp_zer.eqv..true.) then ! first non zero value
                    sp_count_i=sp_count_i+1
                    sparse_dump_i(sp_count_i)= &
                         ix+jy*numxgridn+kz*numxgridn*numygridn
                    sp_zer=.false.
                    sp_fact=sp_fact*(-1.)
                 endif
                 sp_count_r=sp_count_r+1
                 sparse_dump_r(sp_count_r)= &
                      sp_fact* &
                      1.e12*grid(ix,jy,kz) &
                      /volumen(ix,jy,kz)/outnum* &
                      weightair/weightmolar(ks)/densityoutgrid(ix,jy,kz)
                 sparse_dump_u(sp_count_r)= &
                      1.e12*gridsigma(ix,jy,kz)/volumen(ix,jy,kz)/ &
                      outnum*weightair/weightmolar(ks)/ &
                      densityoutgrid(ix,jy,kz)
              else ! concentration is zero
                  sp_zer=.true.
              endif
              end do
            end do
          end do
         write(unitoutgridppt) sp_count_i
         write(unitoutgridppt) (sparse_dump_i(i),i=1,sp_count_i)
         write(unitoutgridppt) sp_count_r
         write(unitoutgridppt) (sparse_dump_r(i),i=1,sp_count_r)
         else

  ! write full vertical resolution

         sp_count_i=0
         sp_count_r=0
         sp_fact=-1.
         sp_zer=.true.
          do kz=1,numzgrid
            do jy=0,numygridn-1
              do ix=0,numxgridn-1
                if (grid(ix,jy,kz).gt.smallnum) then
                  if (sp_zer.eqv..true.) then ! first non zero value
                    sp_count_i=sp_count_i+1
                    sparse_dump_i(sp_count_i)= &
                         ix+jy*numxgridn+kz*numxgridn*numygridn
                    sp_zer=.false.
                    sp_fact=sp_fact*(-1.)
                 endif
                 sp_count_r=sp_count_r+1
                 sparse_dump_r(sp_count_r)= &
                      sp_fact* &
                      1.e12*grid(ix,jy,kz) &
                      /volumen(ix,jy,kz)/outnum* &
                      weightair/weightmolar(ks)/densityoutgrid(ix,jy,kz)
                 sparse_dump_u(sp_count_r)= &
                      1.e12*gridsigma(ix,jy,kz)/volumen(ix,jy,kz)/ &
                      outnum*weightair/weightmolar(ks)/ &
                      densityoutgrid(ix,jy,kz)
              else ! concentration is zero
                  sp_zer=.true.
              endif
              end do
            end do
          end do
         write(unitoutgridppt) sp_count_i
         write(unitoutgridppt) (sparse_dump_i(i),i=1,sp_count_i)
         write(unitoutgridppt) sp_count_r
         write(unitoutgridppt) (sparse_dump_r(i),i=1,sp_count_r)
         endif ! surf_only

      endif ! output for ppt

  end do
  end do

    close(unitoutgridppt)
    close(unitoutgrid)

  end do

  ! RLT Aug 2017
  ! Write out conversion factor for dry air
  inquire(file=path(2)(1:length(2))//'factor_drygrid_nest',exist=lexist)
  if (lexist) then
    ! open and append
    open(unitoutfactor,file=path(2)(1:length(2))//'factor_drygrid_nest',form='unformatted',&
            status='old',action='write',access='append')
  else
    ! create new
    open(unitoutfactor,file=path(2)(1:length(2))//'factor_drygrid_nest',form='unformatted',&
            status='new',action='write')
  endif
  sp_count_i=0
  sp_count_r=0
  sp_fact=-1.
  sp_zer=.true.
  do kz=1,1
    do jy=0,numygridn-1
      do ix=0,numxgridn-1
        if (factor_drygrid(ix,jy,kz).gt.(1.+smallnum).or.factor_drygrid(ix,jy,kz).lt.(1.-smallnum)) then
          if (sp_zer.eqv..true.) then ! first value not equal to one
            sp_count_i=sp_count_i+1
            sparse_dump_i(sp_count_i)= &
                  ix+jy*numxgridn+kz*numxgridn*numygridn
            sp_zer=.false.
            sp_fact=sp_fact*(-1.)
          endif
          sp_count_r=sp_count_r+1
          sparse_dump_r(sp_count_r)= &
               sp_fact*factor_drygrid(ix,jy,kz)
        else ! factor is one
          sp_zer=.true.
        endif
      end do
    end do
  end do
  write(unitoutfactor) sp_count_i
  write(unitoutfactor) (sparse_dump_i(i),i=1,sp_count_i)
  write(unitoutfactor) sp_count_r
  write(unitoutfactor) (sparse_dump_r(i),i=1,sp_count_r)
  close(unitoutfactor)



  ! Reinitialization of grid
  !*************************

  do ks=1,nspec
  do kp=1,maxpointspec_act
    do i=1,numreceptor
      creceptor(i,ks)=0.
    end do
    do jy=0,numygridn-1
      do ix=0,numxgridn-1
        do l=1,nclassunc
          do nage=1,nageclass
            do kz=1,numzgrid
              griduncn(ix,jy,kz,ks,kp,l,nage)=0.
            end do
          end do
        end do
      end do
    end do
  end do
  end do
end subroutine concoutput_surf_nest

subroutine initial_cond_output(itime)
  !                                 i
  !*****************************************************************************
  !                                                                            *
  !     Output of the initial condition sensitivity field.                     *
  !                                                                            *
  !     Author: A. Stohl                                                       *
  !                                                                            *
  !     24 May 1995                                                            *
  !                                                                            *
  !     13 April 1999, Major update: if output size is smaller, dump output    *
  !                    in sparse matrix format; additional output of           *
  !                    uncertainty                                             *
  !                                                                            *
  !     05 April 2000, Major update: output of age classes; output for backward*
  !                    runs is time spent in grid cell times total mass of     *
  !                    species.                                                *
  !                                                                            *
  !     17 February 2002, Appropriate dimensions for backward and forward runs *
  !                       are now specified in file par_mod                    *
  !                                                                            *
  !     June 2006, write grid in sparse matrix with a single write command     *
  !                in order to save disk space                                 *
  !                                                                            *
  !     2008 new sparse matrix format                                          *
  !                                                                            *
  !*****************************************************************************
  !                                                                            *
  ! Variables:                                                                 *
  ! ncells          number of cells with non-zero concentrations               *
  ! sparse          .true. if in sparse matrix format, else .false.            *
  !                                                                            *
  !*****************************************************************************

  use unc_mod

  implicit none

  integer :: itime,i,ix,jy,kz,ks,kp,sp_count_i,sp_count_r
  real :: sp_fact,fact_recept
  real,parameter :: smallnum = tiny(0.0) ! smallest number that can be handled
  logical :: sp_zer
  character(len=3) :: anspec


  !*********************************************************************
  ! Determine the standard deviation of the mean concentration or mixing
  ! ratio (uncertainty of the output) and the dry and wet deposition
  !*********************************************************************

  do ks=1,nspec

    write(anspec,'(i3.3)') ks
    open(97,file=path(2)(1:length(2))//'grid_initial'// &
         '_'//anspec,form='unformatted')
    write(97) itime

    do kp=1,maxpointspec_act

      if (ind_rel.eq.1) then
        fact_recept=rho_rel(kp)
      else
        fact_recept=1.
      endif

  !*******************************************************************
  ! Generate output: may be in concentration (ng/m3) or in mixing
  ! ratio (ppt) or both
  ! Output the position and the values alternated multiplied by
  ! 1 or -1, first line is number of values, number of positions
  ! For backward simulations, the unit is seconds, stored in grid_time
  !*******************************************************************

  ! Write out dummy "wet and dry deposition" fields, to keep same format
  ! as for concentration output
      sp_count_i=0
      sp_count_r=0
      write(97) sp_count_i
      write(97) (sparse_dump_i(i),i=1,sp_count_i)
      write(97) sp_count_r
      write(97) (sparse_dump_r(i),i=1,sp_count_r)
      write(97) sp_count_i
      write(97) (sparse_dump_i(i),i=1,sp_count_i)
      write(97) sp_count_r
      write(97) (sparse_dump_r(i),i=1,sp_count_r)


  ! Write out sensitivity to initial conditions
      sp_count_i=0
      sp_count_r=0
      sp_fact=-1.
      sp_zer=.true.
      do kz=1,numzgrid
        do jy=0,numygrid-1
          do ix=0,numxgrid-1
            if (init_cond(ix,jy,kz,ks,kp).gt.smallnum) then
              if (sp_zer.eqv..true.) then ! first non zero value
                sp_count_i=sp_count_i+1
                sparse_dump_i(sp_count_i)= &
                     ix+jy*numxgrid+kz*numxgrid*numygrid
                sp_zer=.false.
                sp_fact=sp_fact*(-1.)
              endif
              sp_count_r=sp_count_r+1
              sparse_dump_r(sp_count_r)=sp_fact* &
                   init_cond(ix,jy,kz,ks,kp)/xmass(kp,ks)*fact_recept
            else ! concentration is zero
              sp_zer=.true.
            endif
          end do
        end do
      end do
      write(97) sp_count_i
      write(97) (sparse_dump_i(i),i=1,sp_count_i)
      write(97) sp_count_r
      write(97) (sparse_dump_r(i),i=1,sp_count_r)


    end do

    close(97)

  end do
end subroutine initial_cond_output

subroutine initial_cond_output_inversion(itime)
  !                                 i
  !*****************************************************************************
  !                                                                            *
  !     Output of the initial condition sensitivity field.                     *
  !                                                                            *
  !     Author: A. Stohl                                                       *
  !                                                                            *
  !     24 May 1995                                                            *
  !                                                                            *
  !     13 April 1999, Major update: if output size is smaller, dump output    *
  !                    in sparse matrix format; additional output of           *
  !                    uncertainty                                             *
  !                                                                            *
  !     05 April 2000, Major update: output of age classes; output for backward*
  !                    runs is time spent in grid cell times total mass of     *
  !                    species.                                                *
  !                                                                            *
  !     17 February 2002, Appropriate dimensions for backward and forward runs *
  !                       are now specified in file par_mod                    *
  !                                                                            *
  !     June 2006, write grid in sparse matrix with a single write command     *
  !                in order to save disk space                                 *
  !                                                                            *
  !     2008 new sparse matrix format                                          *
  !                                                                            *
  !*****************************************************************************
  !                                                                            *
  ! Variables:                                                                 *
  ! ncells          number of cells with non-zero concentrations               *
  ! sparse          .true. if in sparse matrix format, else .false.            *
  !                                                                            *
  !*****************************************************************************

  use unc_mod

  implicit none

  integer :: itime,i,ix,jy,kz,ks,kp,sp_count_i,sp_count_r
  integer :: jjjjmmdd, ihmmss
  real(kind=dp) :: jul
  real :: sp_fact,fact_recept
  real,parameter :: smallnum = tiny(0.0) ! smallest number that can be handled
  logical :: sp_zer,lexist
  logical,save :: listart=.true.
  logical,save,allocatable,dimension(:) :: listartrel
  character :: adate*8,atime*6
  character :: areldate*8,areltime*6
  character(len=3) :: anspec

  if(listart) then
    allocate(listartrel(maxpointspec_act))
    listartrel(:)=.true.
  endif
  print*, 'listartrel = ',listartrel

  !*********************************************************************
  ! Determine the standard deviation of the mean concentration or mixing
  ! ratio (uncertainty of the output) and the dry and wet deposition
  !*********************************************************************

  do ks=1,nspec

    write(anspec,'(i3.3)') ks

    do kp=1,maxpointspec_act

      ! calculate date of release
      jul=bdate+real(ireleasestart(kp),kind=dp)/86400._dp    ! this is the current day
      call caldate(jul,jjjjmmdd,ihmmss)
      write(areldate,'(i8.8)') jjjjmmdd
      write(areltime,'(i6.6)') ihmmss
      print*, areldate//areltime

      ! calculate date of field
      jul=bdate+real(itime,kind=dp)/86400._dp
      call caldate(jul,jjjjmmdd,ihmmss)
      write(adate,'(i8.8)') jjjjmmdd
      write(atime,'(i6.6)') ihmmss
      print*, adate//atime

      inquire(file=path(2)(1:length(2))//'grid_initial_'//areldate// &
         areltime//'_'//anspec,exist=lexist)
      if(lexist.and..not.listartrel(kp)) then
        ! open and append to existing file
        open(97,file=path(2)(1:length(2))//'grid_initial_'//areldate// &
             areltime//'_'//anspec,form='unformatted',status='old',action='write',access='append') 
      else
        ! open new file
        open(97,file=path(2)(1:length(2))//'grid_initial_'//areldate// &
             areltime//'_'//anspec,form='unformatted',status='replace',action='write')
      endif
      write(97) jjjjmmdd
      write(97) ihmmss

      listartrel(kp)=.false.

      if (ind_rel.eq.1) then
        fact_recept=rho_rel(kp)
      else
        fact_recept=1.
      endif

  !*******************************************************************
  ! Generate output: may be in concentration (ng/m3) or in mixing
  ! ratio (ppt) or both
  ! Output the position and the values alternated multiplied by
  ! 1 or -1, first line is number of values, number of positions
  ! For backward simulations, the unit is seconds, stored in grid_time
  !*******************************************************************

  ! Write out dummy "wet and dry deposition" fields, to keep same format
  ! as for concentration output

  ! Write out sensitivity to initial conditions
      sp_count_i=0
      sp_count_r=0
      sp_fact=-1.
      sp_zer=.true.
      do kz=1,numzgrid
        do jy=0,numygrid-1
          do ix=0,numxgrid-1
            if (init_cond(ix,jy,kz,ks,kp).gt.smallnum) then
              if (sp_zer.eqv..true.) then ! first non zero value
                sp_count_i=sp_count_i+1
                sparse_dump_i(sp_count_i)= &
                     ix+jy*numxgrid+kz*numxgrid*numygrid
                sp_zer=.false.
                sp_fact=sp_fact*(-1.)
              endif
              sp_count_r=sp_count_r+1
              sparse_dump_r(sp_count_r)=sp_fact* &
                   init_cond(ix,jy,kz,ks,kp)/xmass(kp,ks)*fact_recept
            else ! concentration is zero
              sp_zer=.true.
            endif
          end do
        end do
      end do
      write(97) sp_count_i
      write(97) (sparse_dump_i(i),i=1,sp_count_i)
      write(97) sp_count_r
      write(97) (sparse_dump_r(i),i=1,sp_count_r)

      close(97)

    end do

  end do

  ! reset listart
  if (listart) then
    listart=.false.
  endif
end subroutine initial_cond_output_inversion

end module binary_output_mod