! SPDX-FileCopyrightText: FLEXPART 1998-2019, see flexpart_license.txt
! SPDX-License-Identifier: GPL-3.0-or-later

module outg_mod

  use par_mod, only: dep_prec, sp

  implicit none

  real,allocatable, dimension (:) :: outheight
  real,allocatable, dimension (:) :: outheighthalf
  real,allocatable, dimension (:,:) :: oroout
  real,allocatable, dimension (:,:) :: orooutn
  real,allocatable, dimension (:,:) :: area
  real,allocatable, dimension (:,:) :: arean
  real,allocatable, dimension (:,:,:) :: volume
  real,allocatable, dimension (:,:,:) :: volumen
  real,allocatable, dimension (:,:,:) :: areaeast
  real,allocatable, dimension (:,:,:) :: areanorth
  real,allocatable, dimension (:,:,:) :: densityoutgrid
  real,allocatable, dimension (:,:,:) :: densitydrygrid ! added RLT 
  real,allocatable, dimension (:,:,:) :: factor_drygrid ! added RLT 
  real,allocatable, dimension (:,:,:) :: factor3d
  real,allocatable, dimension (:,:,:) :: grid
  real(sp),allocatable, dimension (:,:) :: wetgrid
  real(sp),allocatable, dimension (:,:) :: drygrid
  real,allocatable, dimension (:,:,:) :: gridsigma
  real(dep_prec),allocatable, dimension (:,:) :: drygridsigma
  real(dep_prec),allocatable, dimension (:,:) :: wetgridsigma
  real,allocatable, dimension (:) :: sparse_dump_r
  real,allocatable, dimension (:) :: sparse_dump_u
  integer,allocatable, dimension (:) :: sparse_dump_i

  real,allocatable, dimension (:,:,:,:,:,:,:) :: flux
  real,allocatable, dimension (:,:,:,:,:,:,:,:) :: flux_omp

  real,allocatable, dimension (:,:,:,:,:) :: init_cond
  real,allocatable, dimension (:,:,:,:,:,:) :: init_cond_omp

  !1 fluxw west - east
  !2 fluxe east - west
  !3 fluxs south - north
  !4 fluxn north - south
  !5 fluxu upward
  !6 fluxd downward
  !real,allocatable, dimension (:,:,:) :: areanorth
  !real,allocatable, dimension (:,:,:) :: areaeast
contains

subroutine outgrid_init
  !
  !*****************************************************************************
  !                                                                            *
  !  This routine initializes the output grids                                 *
  !                                                                            *
  !     Author: A. Stohl                                                       *
  !                                                                            *
  !     7 August 2002                                                          *
  !                                                                            *
  !*****************************************************************************
  !                                                                            *
  ! Variables:                                                                 *
  !                                                                            *
  ! area               surface area of all output grid cells                   *
  ! areaeast           eastward facing wall area of all output grid cells      *
  ! areanorth          northward facing wall area of all output grid cells     *
  ! volume             volumes of all output grid cells                        *
  !                                                                            *
  !*****************************************************************************

  use oh_mod
  use unc_mod
  use par_mod
  use com_mod
  use windfields_mod

  implicit none

  integer :: ix,jy,kz,i,nage,l,iix,jjy,ixp,jyp,i1,j1,j,ngrid
  integer :: ks,kp,stat
  real :: ylat,gridarea,ylatp,ylatm,hzone,cosfactm,cosfactp
  real :: xlon,xl,yl,ddx,ddy,rddx,rddy,p1,p2,p3,p4,xtn,ytn,oroh
  real :: eps

  eps=nxmax/3.e5

  ! Compute surface area and volume of each grid cell: area, volume;
  ! and the areas of the northward and eastward facing walls: areaeast, areanorth
  !***********************************************************************
  do jy=0,numygrid-1
    ylat=outlat0+(real(jy)+0.5)*dyout
    ylatp=ylat+0.5*dyout
    ylatm=ylat-0.5*dyout
    if ((ylatm.lt.0).and.(ylatp.gt.0.)) then
      hzone=dyout*r_earth*pi180
    else

  ! Calculate area of grid cell with formula M=2*pi*R*h*dx/360,
  ! see Netz, Formeln der Mathematik, 5. Auflage (1983), p.90
  !************************************************************

      cosfactp=cos(ylatp*pi180)
      cosfactm=cos(ylatm*pi180)
      if (cosfactp.lt.cosfactm) then
        hzone=sqrt(1-cosfactp**2)- &
             sqrt(1-cosfactm**2)
        hzone=hzone*r_earth
      else
        hzone=sqrt(1-cosfactm**2)- &
             sqrt(1-cosfactp**2)
        hzone=hzone*r_earth
      endif
    endif

  ! Surface are of a grid cell at a latitude ylat
  !**********************************************

    gridarea=2.*pi*r_earth*hzone*dxout/360.

    do ix=0,numxgrid-1
      area(ix,jy)=gridarea

  ! Volume = area x box height
  !***************************

      volume(ix,jy,1)=area(ix,jy)*outheight(1)
      areaeast(ix,jy,1)=dyout*r_earth*pi180*outheight(1)
      areanorth(ix,jy,1)=cos(ylat*pi180)*dxout*r_earth*pi180* &
           outheight(1)
      do kz=2,numzgrid
        areaeast(ix,jy,kz)=dyout*r_earth*pi180* &
             (outheight(kz)-outheight(kz-1))
        areanorth(ix,jy,kz)=cos(ylat*pi180)*dxout*r_earth*pi180* &
             (outheight(kz)-outheight(kz-1))
        volume(ix,jy,kz)=area(ix,jy)*(outheight(kz)-outheight(kz-1))
      end do
    end do
  end do




  !******************************************************************
  ! Determine average height of model topography in output grid cells
  !******************************************************************

  ! Loop over all output grid cells
  !********************************

  do jjy=0,numygrid-1
    do iix=0,numxgrid-1
      oroh=0.

  ! Take 100 samples of the topography in every grid cell
  !******************************************************

      do j1=1,10
        ylat=outlat0+(real(jjy)+real(j1)/10.-0.05)*dyout
        yl=(ylat-ylat0)/dy
        do i1=1,10
          xlon=outlon0+(real(iix)+real(i1)/10.-0.05)*dxout
          xl=(xlon-xlon0)/dx

  ! Determine the nest we are in
  !*****************************

          ngrid=0
          do j=numbnests,1,-1
            if ((xl.gt.xln(j)+eps).and.(xl.lt.xrn(j)-eps).and. &
                 (yl.gt.yln(j)+eps).and.(yl.lt.yrn(j)-eps)) then
              ngrid=j
              exit
            endif
          end do

  ! Determine (nested) grid coordinates and auxiliary parameters used for interpolation
  !*****************************************************************************

          if (ngrid.gt.0) then
            xtn=(xl-xln(ngrid))*xresoln(ngrid)
            ytn=(yl-yln(ngrid))*yresoln(ngrid)
            ix=int(xtn)
            jy=int(ytn)
            ddy=ytn-real(jy)
            ddx=xtn-real(ix)
          else
            ix=int(xl)
            jy=int(yl)
            ddy=yl-real(jy)
            ddx=xl-real(ix)
          endif
          ixp=ix+1
          jyp=jy+1
          rddx=1.-ddx
          rddy=1.-ddy
          p1=rddx*rddy
          p2=ddx*rddy
          p3=rddx*ddy
          p4=ddx*ddy

          if (ngrid.gt.0) then
            oroh=oroh+p1*oron(ix ,jy ,ngrid) &
                 + p2*oron(ixp,jy ,ngrid) &
                 + p3*oron(ix ,jyp,ngrid) &
                 + p4*oron(ixp,jyp,ngrid)
          else
            oroh=oroh+p1*oro(ix ,jy) &
                 + p2*oro(ixp,jy) &
                 + p3*oro(ix ,jyp) &
                 + p4*oro(ixp,jyp)
          endif
        end do
      end do

  ! Divide by the number of samples taken
  !**************************************

      oroout(iix,jjy)=oroh/100.
    end do
  end do

  ! if necessary allocate flux fields
  if (iflux.eq.1) then
    allocate(flux(6,0:numxgrid-1,0:numygrid-1,numzgrid, &
         1:nspec,1:maxpointspec_act,1:nageclass),stat=stat)
    if (stat.ne.0) write(*,*)'ERROR: could not allocate flux array '
#ifdef _OPENMP
    allocate(flux_omp(6,0:numxgrid-1,0:numygrid-1,numzgrid, &
         1:nspec,1:maxpointspec_act,1:nageclass,numthreads))
    if (stat.ne.0) write(*,*)'ERROR: could not allocate flux_omp array '
#endif
  endif

  ! gridunc,griduncn        uncertainty of outputted concentrations
  allocate(gridunc(0:numxgrid-1,0:numygrid-1,numzgrid,maxspec, &
       maxpointspec_act,nclassunc,maxageclass),stat=stat)
  if (stat.ne.0) write(*,*)'ERROR: could not allocate gridunc'
#ifdef _OPENMP
  allocate(gridunc_omp(0:numxgrid-1,0:numygrid-1,numzgrid,maxspec, &
       maxpointspec_act,nclassunc,maxageclass,numthreads_grid),stat=stat)
  if (stat.ne.0) then
    write(*,*)'ERROR: could not allocate gridunc_omp',
    write(*,*)'increase the memory or reduce max_numthreads_grid in par_mod.f90.'
    stop
#endif
  if (ldirect.gt.0) then
    allocate(wetgridunc(0:numxgrid-1,0:numygrid-1,maxspec, &
         maxpointspec_act,nclassunc,maxageclass),stat=stat)
    if (stat.ne.0) write(*,*)'ERROR: could not allocate wetgridunc'
    allocate(drygridunc(0:numxgrid-1,0:numygrid-1,maxspec, &
         maxpointspec_act,nclassunc,maxageclass),stat=stat)
    if (stat.ne.0) write(*,*)'ERROR: could not allocate drygridunc'
#ifdef _OPENMP
    allocate(wetgridunc_omp(0:numxgrid-1,0:numygrid-1,maxspec, &
         maxpointspec_act,nclassunc,maxageclass,numthreads_grid),stat=stat)
    if (stat.ne.0) write(*,*)'ERROR: could not allocate wetgridunc_omp'
    allocate(drygridunc_omp(0:numxgrid-1,0:numygrid-1,maxspec, &
         maxpointspec_act,nclassunc,maxageclass,numthreads_grid),stat=stat)
    if (stat.ne.0) write(*,*)'ERROR: could not allocate drygridunc_omp'
#endif
  endif

#ifdef USE_MPIINPLACE
#else
! Extra field for totals at MPI root process
  if (lroot.and.mpi_mode.gt.0) then
! If MPI_IN_PLACE option is not used in mpi_mod.f90::mpif_tm_reduce_grid(),
! then an aux array is needed for parallel grid reduction
    allocate(gridunc0(0:numxgrid-1,0:numygrid-1,numzgrid,maxspec, &
         maxpointspec_act,nclassunc,maxageclass),stat=stat)
    if (stat.ne.0) write(*,*)'ERROR: could not allocate gridunc0'
  else if (.not.lroot.and.mpi_mode.gt.0) then
    allocate(gridunc0(1,1,1,1,1,1,1),stat=stat)
    if (stat.ne.0) write(*,*)'ERROR: could not allocate gridunc0'
  end if
#endif
  if (ldirect.gt.0) then
    if (lroot.and.mpi_mode.gt.0) then
      allocate(wetgridunc0(0:numxgrid-1,0:numygrid-1,maxspec, &
           maxpointspec_act,nclassunc,maxageclass),stat=stat)
      if (stat.ne.0) write(*,*)'ERROR: could not allocate wetgridunc0'
      allocate(drygridunc0(0:numxgrid-1,0:numygrid-1,maxspec, &
           maxpointspec_act,nclassunc,maxageclass),stat=stat)
      if (stat.ne.0) write(*,*)'ERROR: could not allocate drygridunc0'

  ! allocate a dummy to avoid compilator complaints
    else if (.not.lroot.and.mpi_mode.gt.0) then
      allocate(wetgridunc0(1,1,1,1,1,1),stat=stat)
      allocate(drygridunc0(1,1,1,1,1,1),stat=stat)
    end if
  end if

  !write (*,*) 'Dimensions for fields', numxgrid,numygrid, &
  !     maxspec,maxpointspec_act,nclassunc,maxageclass

  if (lroot) then
    write (*,*) 'Allocating fields for global output (x,y): ', numxgrid,numygrid
    write (*,*) 'Allocating fields for nested output (x,y): ', numxgridn,numygridn 
  end if

  ! allocate fields for concoutput with maximum dimension of outgrid
  ! and outgrid_nest

  allocate(gridsigma(0:max(numxgrid,numxgridn)-1, &
       0:max(numygrid,numygridn)-1,numzgrid),stat=stat)
    if (stat.ne.0) write(*,*)'ERROR: could not allocate gridunc'
  allocate(grid(0:max(numxgrid,numxgridn)-1, &
       0:max(numygrid,numygridn)-1,numzgrid),stat=stat)
    if (stat.ne.0) write(*,*)'ERROR: could not allocate gridunc'
  allocate(densityoutgrid(0:max(numxgrid,numxgridn)-1, &
       0:max(numygrid,numygridn)-1,numzgrid),stat=stat)
    if (stat.ne.0) write(*,*)'ERROR: could not allocate gridunc'
  ! RLT
  allocate(densitydrygrid(0:max(numxgrid,numxgridn)-1, &
       0:max(numygrid,numygridn)-1,numzgrid),stat=stat)
    if (stat.ne.0) write(*,*)'ERROR: could not allocate gridunc'
  allocate(factor_drygrid(0:max(numxgrid,numxgridn)-1, &
       0:max(numygrid,numygridn)-1,numzgrid),stat=stat)
    if (stat.ne.0) write(*,*)'ERROR: could not allocate gridunc'

  allocate(factor3d(0:max(numxgrid,numxgridn)-1, &
       0:max(numygrid,numygridn)-1,numzgrid),stat=stat)
    if (stat.ne.0) write(*,*)'ERROR: could not allocate gridunc'
  allocate(sparse_dump_r(max(numxgrid,numxgridn)* &
       max(numygrid,numygridn)*numzgrid),stat=stat)
    if (stat.ne.0) write(*,*)'ERROR: could not allocate gridunc'

   allocate(sparse_dump_u(max(numxgrid,numxgridn)* &
       max(numygrid,numygridn)*numzgrid),stat=stat)
        if (stat.ne.0) write(*,*)'ERROR: could not allocate gridunc'

  allocate(sparse_dump_i(max(numxgrid,numxgridn)* &
       max(numygrid,numygridn)*numzgrid),stat=stat)
    if (stat.ne.0) write(*,*)'ERROR: could not allocate gridunc'

  ! deposition fields are only allocated for forward runs
  if (ldirect.gt.0) then
     allocate(wetgridsigma(0:max(numxgrid,numxgridn)-1, &
          0:max(numygrid,numygridn)-1),stat=stat)
     if (stat.ne.0) write(*,*)'ERROR: could not allocate gridunc'
     allocate(drygridsigma(0:max(numxgrid,numxgridn)-1, &
          0:max(numygrid,numygridn)-1),stat=stat)
     if (stat.ne.0) write(*,*)'ERROR: could not allocate gridunc'
     allocate(wetgrid(0:max(numxgrid,numxgridn)-1, &
          0:max(numygrid,numygridn)-1),stat=stat)
     if (stat.ne.0) write(*,*)'ERROR: could not allocate gridunc'
     allocate(drygrid(0:max(numxgrid,numxgridn)-1, &
          0:max(numygrid,numygridn)-1),stat=stat)
     if (stat.ne.0) write(*,*)'ERROR: could not allocate gridunc'
  endif

  ! Initial condition field

  if (linit_cond.gt.0) then
    allocate(init_cond(0:numxgrid-1,0:numygrid-1,numzgrid,maxspec, &
         maxpointspec_act),stat=stat)
    if (stat.ne.0) write(*,*)'ERROR: could not allocate init_cond'
#ifdef _OPENMP
    allocate(init_cond_omp(0:numxgrid-1,0:numygrid-1,numzgrid,maxspec, &
         maxpointspec_act,numthreads),stat=stat)
    if (stat.ne.0) write(*,*)'ERROR: could not allocate init_cond_omp'
#endif
  endif

  !************************
  ! Initialize output grids
  !************************

  do ks=1,nspec
    do kp=1,maxpointspec_act
      if (numreceptor.gt.0) then
        do i=1,numreceptor
      ! Receptor points
          creceptor(i,ks)=0.
        end do
      endif
      do nage=1,nageclass
        do jy=0,numygrid-1
          do ix=0,numxgrid-1
            do l=1,nclassunc
    ! Deposition fields
              if (ldirect.gt.0) then
              wetgridunc(ix,jy,ks,kp,l,nage)=0.
              drygridunc(ix,jy,ks,kp,l,nage)=0.
#ifdef _OPENMP
              wetgridunc_omp(ix,jy,ks,kp,l,nage,:)=0.
              drygridunc_omp(ix,jy,ks,kp,l,nage,:)=0.
#endif
              endif
              do kz=1,numzgrid
                if (iflux.eq.1) then
    ! Flux fields
                   do i=1,5
                     flux(i,ix,jy,kz,ks,kp,nage)=0.
#ifdef _OPENMP
                     flux_omp(i,ix,jy,kz,ks,kp,nage,:)=0.
#endif
                   end do
                endif
    ! Initial condition field
                if ((l.eq.1).and.(nage.eq.1).and.(linit_cond.gt.0)) then
                  init_cond(ix,jy,kz,ks,kp)=0.
#ifdef _OPENMP
                  init_cond_omp(ix,jy,kz,ks,kp,:)=0.
#endif
                endif
    ! Concentration fields
                gridunc(ix,jy,kz,ks,kp,l,nage)=0.
#ifdef _OPENMP
                gridunc_omp(ix,jy,kz,ks,kp,l,nage,:)=0.
#endif
              end do
            end do
          end do
        end do
      end do
    end do
  end do
end subroutine outgrid_init

subroutine outgrid_init_nest
  !
  !*****************************************************************************
  !                                                                            *
  !  This routine calculates, for each grid cell of the output nest, the       *
  !  volume and the surface area.                                              *
  !                                                                            *
  !     Author: A. Stohl                                                       *
  !                                                                            *
  !    30 August 2004                                                          *
  !                                                                            *
  !*****************************************************************************
  !                                                                            *
  ! Variables:                                                                 *
  !                                                                            *
  ! arean              surface area of all output nest cells                   *
  ! volumen            volumes of all output nest cells                        *
  !                                                                            *
  !*****************************************************************************

  use unc_mod
  use par_mod
  use com_mod
  use windfields_mod

  implicit none

  integer :: ix,jy,kz,ks,kp,nage,l,iix,jjy,ixp,jyp,i1,j1,j,ngrid
  integer :: stat
  real :: ylat,gridarea,ylatp,ylatm,hzone,cosfactm,cosfactp
  real :: xlon,xl,yl,ddx,ddy,rddx,rddy,p1,p2,p3,p4,xtn,ytn,oroh
  real :: eps

  eps=nxmax/3.e5

  ! gridunc,griduncn        uncertainty of outputted concentrations
  allocate(griduncn(0:numxgridn-1,0:numygridn-1,numzgrid,maxspec, &
       maxpointspec_act,nclassunc,maxageclass),stat=stat)
  if (stat.ne.0) write(*,*)'ERROR:could not allocate nested gridunc'
#ifdef _OPENMP
  allocate(griduncn_omp(0:numxgridn-1,0:numygridn-1,numzgrid,maxspec, &
       maxpointspec_act,nclassunc,maxageclass,numthreads_grid),stat=stat)
  if (stat.ne.0) write(*,*)'ERROR:could not allocate nested gridunc_omp'
#endif

  if (ldirect.gt.0) then
    allocate(wetgriduncn(0:numxgridn-1,0:numygridn-1,maxspec, &
         maxpointspec_act,nclassunc,maxageclass),stat=stat)
    if (stat.ne.0) write(*,*)'ERROR:could not allocate nested gridunc'
    allocate(drygriduncn(0:numxgridn-1,0:numygridn-1,maxspec, &
         maxpointspec_act,nclassunc,maxageclass),stat=stat)
    if (stat.ne.0) write(*,*)'ERROR:could not allocate nested gridunc'
#ifdef _OPENMP
    allocate(wetgriduncn_omp(0:numxgridn-1,0:numygridn-1,maxspec, &
         maxpointspec_act,nclassunc,maxageclass,numthreads_grid),stat=stat)
    if (stat.ne.0) write(*,*)'ERROR:could not allocate nested wetgridunc_omp'
    allocate(drygriduncn_omp(0:numxgridn-1,0:numygridn-1,maxspec, &
         maxpointspec_act,nclassunc,maxageclass,numthreads_grid),stat=stat)
    if (stat.ne.0) write(*,*)'ERROR:could not allocate nested drygriduncn_omp'
#endif
  endif

#ifdef USE_MPIINPLACE
#else
  ! Extra field for totals at MPI root process
  if (lroot.and.mpi_mode.gt.0) then
  ! If MPI_IN_PLACE option is not used in mpi_mod.f90::mpif_tm_reduce_grid_nest(),
  ! then an aux array is needed for parallel grid reduction
    allocate(griduncn0(0:numxgridn-1,0:numygridn-1,numzgrid,maxspec, &
         maxpointspec_act,nclassunc,maxageclass),stat=stat)
    if (stat.ne.0) write(*,*)'ERROR:could not allocate nested gridunc'
  ! allocate a dummy to avoid compilator complaints
  else if (.not.lroot.and.mpi_mode.gt.0) then
    allocate(griduncn0(1,1,1,1,1,1,1),stat=stat)
  end if
#endif
  if (ldirect.gt.0) then
    if (lroot.and.mpi_mode.gt.0) then
      allocate(wetgriduncn0(0:numxgridn-1,0:numygridn-1,maxspec, &
           maxpointspec_act,nclassunc,maxageclass),stat=stat)
      if (stat.ne.0) write(*,*)'ERROR:could not allocate nested gridunc'
      allocate(drygriduncn0(0:numxgridn-1,0:numygridn-1,maxspec, &
           maxpointspec_act,nclassunc,maxageclass),stat=stat)
      if (stat.ne.0) write(*,*)'ERROR:could not allocate nested gridunc'
  !  endif
  ! allocate a dummy to avoid compilator complaints
    else if (.not.lroot.and.mpi_mode.gt.0) then
      allocate(wetgriduncn0(1,1,1,1,1,1),stat=stat)
      allocate(drygriduncn0(1,1,1,1,1,1),stat=stat)
    end if
  end if

  ! Compute surface area and volume of each grid cell: area, volume;
  ! and the areas of the northward and eastward facing walls: areaeast, areanorth
  !***********************************************************************

  do jy=0,numygridn-1
    ylat=outlat0n+(real(jy)+0.5)*dyoutn
    ylatp=ylat+0.5*dyoutn
    ylatm=ylat-0.5*dyoutn
    if ((ylatm.lt.0).and.(ylatp.gt.0.)) then
      hzone=dyoutn*r_earth*pi180
    else

  ! Calculate area of grid cell with formula M=2*pi*R*h*dx/360,
  ! see Netz, Formeln der Mathematik, 5. Auflage (1983), p.90
  !************************************************************

      cosfactp=cos(ylatp*pi180)
      cosfactm=cos(ylatm*pi180)
      if (cosfactp.lt.cosfactm) then
        hzone=sqrt(1-cosfactp**2)- &
             sqrt(1-cosfactm**2)
        hzone=hzone*r_earth
      else
        hzone=sqrt(1-cosfactm**2)- &
             sqrt(1-cosfactp**2)
        hzone=hzone*r_earth
      endif
    endif



  ! Surface are of a grid cell at a latitude ylat
  !**********************************************

    gridarea=2.*pi*r_earth*hzone*dxoutn/360.

    do ix=0,numxgridn-1
      arean(ix,jy)=gridarea

  ! Volume = area x box height
  !***************************

      volumen(ix,jy,1)=arean(ix,jy)*outheight(1)
      do kz=2,numzgrid
        volumen(ix,jy,kz)=arean(ix,jy)*(outheight(kz)-outheight(kz-1))
      end do
    end do
  end do


  !**************************************************************************
  ! Determine average height of model topography in nesteed output grid cells
  !**************************************************************************

  ! Loop over all output grid cells
  !********************************

  do jjy=0,numygridn-1
    do iix=0,numxgridn-1
      oroh=0.

  ! Take 100 samples of the topography in every grid cell
  !******************************************************

      do j1=1,10
        ylat=outlat0n+(real(jjy)+real(j1)/10.-0.05)*dyoutn
        yl=(ylat-ylat0)/dy
        do i1=1,10
          xlon=outlon0n+(real(iix)+real(i1)/10.-0.05)*dxoutn
          xl=(xlon-xlon0)/dx

  ! Determine the nest we are in
  !*****************************

          ngrid=0
          do j=numbnests,1,-1
            if ((xl.gt.xln(j)+eps).and.(xl.lt.xrn(j)-eps).and. &
                 (yl.gt.yln(j)+eps).and.(yl.lt.yrn(j)-eps)) then
              ngrid=j
              exit
            endif
          end do

  ! Determine (nested) grid coordinates and auxiliary parameters used for interpolation
  !*****************************************************************************

          if (ngrid.gt.0) then
            xtn=(xl-xln(ngrid))*xresoln(ngrid)
            ytn=(yl-yln(ngrid))*yresoln(ngrid)
            ix=int(xtn)
            jy=int(ytn)
            ddy=ytn-real(jy)
            ddx=xtn-real(ix)
          else
            ix=int(xl)
            jy=int(yl)
            ddy=yl-real(jy)
            ddx=xl-real(ix)
          endif
          ixp=ix+1
          jyp=jy+1
          rddx=1.-ddx
          rddy=1.-ddy
          p1=rddx*rddy
          p2=ddx*rddy
          p3=rddx*ddy
          p4=ddx*ddy

          if (ngrid.gt.0) then
            oroh=oroh+p1*oron(ix ,jy ,ngrid) &
                 + p2*oron(ixp,jy ,ngrid) &
                 + p3*oron(ix ,jyp,ngrid) &
                 + p4*oron(ixp,jyp,ngrid)
          else
            oroh=oroh+p1*oro(ix ,jy) &
                 + p2*oro(ixp,jy) &
                 + p3*oro(ix ,jyp) &
                 + p4*oro(ixp,jyp)
          endif
        end do
      end do

  ! Divide by the number of samples taken
  !**************************************

      orooutn(iix,jjy)=oroh/100.
    end do
  end do

  !*******************************
  ! Initialization of output grids
  !*******************************

  do kp=1,maxpointspec_act
    do ks=1,nspec
      do nage=1,nageclass
        do jy=0,numygridn-1
          do ix=0,numxgridn-1
            do l=1,nclassunc
  ! Deposition fields
              if (ldirect.gt.0) then
                wetgriduncn(ix,jy,ks,kp,l,nage)=0.
                drygriduncn(ix,jy,ks,kp,l,nage)=0.
#ifdef _OPENMP
                wetgriduncn_omp(ix,jy,ks,kp,l,nage,:)=0.
                drygriduncn_omp(ix,jy,ks,kp,l,nage,:)=0.
#endif
              endif
  ! Concentration fields
              do kz=1,numzgrid
                griduncn(ix,jy,kz,ks,kp,l,nage)=0.
              end do
            end do
          end do
        end do
      end do
    end do
  end do
end subroutine outgrid_init_nest

subroutine initial_cond_calc(itime,i,thread)
  !                               i   i
  !*****************************************************************************
  !                                                                            *
  !     Calculation of the sensitivity to initial conditions for BW runs       *
  !                                                                            *
  !     Author: A. Stohl                                                       *
  !                                                                            *
  !     15 January 2010                                                        *
  !                                                                            *
  !*****************************************************************************

  use par_mod
  use com_mod
  use interpol_mod, only: interpol_density,ix,jy,ixp,jyp
  use coordinates_ecmwf
  use particle_mod

  implicit none

  integer, intent(in) :: itime,i,thread
  integer :: kz,ks
  integer :: il,ind,indz,indzp,nrelpointer
  real :: rddx,rddy,p1,p2,p3,p4,dz1,dz2,dz
  real :: ddx,ddy
  real :: rhoprof(2),rhoi,xl,yl,wx,wy,w
  integer :: mind2
  ! mind2        eso: pointer to 2nd windfield in memory


  ! For forward simulations, make a loop over the number of species;
  ! for backward simulations, make an additional loop over the release points
  !**************************************************************************


  if (.not. part(i)%alive) return

  ! Depending on output option, calculate air density or set it to 1
  ! linit_cond: 1=mass unit, 2=mass mixing ratio unit
  !*****************************************************************


  if (linit_cond.eq.1) then     ! mass unit
    call update_zeta_to_z(itime,i)
    call interpol_density(itime,i,rhoi)
  elseif (linit_cond.eq.2) then    ! mass mixing ratio unit
    rhoi=1.
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
    if (real(outheight(kz),kind=dp).gt.part(i)%z) exit
  end do

  if (kz.le.numzgrid) then           ! inside output domain


    xl=(part(i)%xlon*dx+xoutshift)/dxout
    yl=(part(i)%ylat*dy+youtshift)/dyout
    ix=int(xl)
    if (xl.lt.0.) ix=ix-1
    jy=int(yl)
    if (yl.lt.0.) jy=jy-1


  ! If a particle is close to the domain boundary, do not use the kernel either
  !****************************************************************************

    if ((xl.lt.0.5).or.(yl.lt.0.5).or. &
         (xl.gt.real(numxgrid-1)-0.5).or. &
         (yl.gt.real(numygrid-1)-0.5)) then             ! no kernel, direct attribution to grid cell
      if ((ix.ge.0).and.(jy.ge.0).and.(ix.le.numxgrid-1).and. &
           (jy.le.numygrid-1)) then
        do ks=1,nspec
#ifdef _OPENMP
          init_cond_omp(ix,jy,kz,ks,nrelpointer,thread)= &
               init_cond_omp(ix,jy,kz,ks,nrelpointer,thread)+ &
               part(i)%mass(ks)/rhoi
#else
          init_cond(ix,jy,kz,ks,nrelpointer)= &
               init_cond(ix,jy,kz,ks,nrelpointer)+ &
               part(i)%mass(ks)/rhoi
#endif
        end do
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
          do ks=1,nspec
#ifdef _OPENMP
            init_cond_omp(ix,jy,kz,ks,nrelpointer,thread)= &
                 init_cond_omp(ix,jy,kz,ks,nrelpointer,thread)+part(i)%mass(ks)/rhoi*w
#else
            init_cond(ix,jy,kz,ks,nrelpointer)= &
                 init_cond(ix,jy,kz,ks,nrelpointer)+part(i)%mass(ks)/rhoi*w
#endif
          end do
        endif

        if ((jyp.ge.0).and.(jyp.le.numygrid-1)) then
          w=wx*(1.-wy)
          do ks=1,nspec
#ifdef _OPENMP
            init_cond_omp(ix,jyp,kz,ks,nrelpointer,thread)= &
                 init_cond_omp(ix,jyp,kz,ks,nrelpointer,thread)+part(i)%mass(ks)/rhoi*w
#else
            init_cond(ix,jyp,kz,ks,nrelpointer)= &
                 init_cond(ix,jyp,kz,ks,nrelpointer)+part(i)%mass(ks)/rhoi*w
#endif
          end do
        endif
      endif


      if ((ixp.ge.0).and.(ixp.le.numxgrid-1)) then
        if ((jyp.ge.0).and.(jyp.le.numygrid-1)) then
          w=(1.-wx)*(1.-wy)
          do ks=1,nspec
#ifdef _OPENMP
            init_cond_omp(ixp,jyp,kz,ks,nrelpointer,thread)= &
                 init_cond_omp(ixp,jyp,kz,ks,nrelpointer,thread)+part(i)%mass(ks)/rhoi*w
#else
            init_cond(ixp,jyp,kz,ks,nrelpointer)= &
                 init_cond(ixp,jyp,kz,ks,nrelpointer)+part(i)%mass(ks)/rhoi*w
#endif
          end do
        endif

        if ((jy.ge.0).and.(jy.le.numygrid-1)) then
          w=(1.-wx)*wy
          do ks=1,nspec
#ifdef _OPENMP
            init_cond_omp(ixp,jy,kz,ks,nrelpointer,thread)= &
                 init_cond_omp(ixp,jy,kz,ks,nrelpointer,thread)+part(i)%mass(ks)/rhoi*w
#else
            init_cond(ixp,jy,kz,ks,nrelpointer)= &
                 init_cond(ixp,jy,kz,ks,nrelpointer)+part(i)%mass(ks)/rhoi*w
#endif
          end do
        endif
      endif
    endif

  endif

end subroutine initial_cond_calc


end module outg_mod
