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
    if (iout.ne.0) then ! No gridded output
#ifdef USE_NCF
      if (lnetcdfout.eq.1) then 
        call writeheader_netcdf(lnest=.false.)
      else 
        call writeheader_binary
      end if

      if (nested_output.eq.1) then
        if (lnetcdfout.eq.1) then
          call writeheader_netcdf(lnest=.true.)
        else if ((nested_output.eq.1).and.(surf_only.ne.1)) then
          call writeheader_binary_nest
        else if ((nested_output.eq.1).and.(surf_only.eq.1)) then
          call writeheader_binary_nest_surf
        else if ((nested_output.ne.1).and.(surf_only.eq.1)) then 
          call writeheader_binary_surf
        endif
      endif
#else
      call writeheader_binary

      !if (nested_output.eq.1) call writeheader_nest
      if ((nested_output.eq.1).and.(surf_only.ne.1)) call writeheader_binary_nest
      if ((nested_output.eq.1).and.(surf_only.eq.1)) call writeheader_binary_nest_surf
      if ((nested_output.ne.1).and.(surf_only.eq.1)) call writeheader_binary_surf
#endif
    endif ! iout.ne.0
    ! FLEXPART 9.2 ticket ?? write header in ASCII format 
    call writeheader_txt

    ! NetCDF only: Create file for storing initial particle positions.
#ifdef USE_NCF
    if (mdomainfill.eq.0) then
      if (ldirect.eq.1) then
        call create_particles_initialoutput(ibtime,ibdate,ibtime,ibdate)
      else
        call create_particles_initialoutput(ietime,iedate,ietime,iedate)
      endif
    endif
    ! Create header files for files that store the particle dump output
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

  use interpol_mod
  use coordinates_ecmwf
  use particle_mod
#ifdef USE_NCF
  use netcdf
  use netcdf_output_mod, only: partoutput_netcdf,open_partoutput_file,close_partoutput_file
  use omp_lib, only: OMP_GET_THREAD_NUM
#endif

  implicit none

  real(kind=dp) :: jul
  integer :: itime,i,j,jjjjmmdd,ihmmss
  real :: tr(2),hm(2)
  character :: adate*8,atime*6

  real :: xlon(numpart),ylat(numpart),ztemp1,ztemp2
  real :: tti(numpart),rhoi(numpart),pvi(numpart),qvi(numpart),pri(numpart)
  real :: topo(numpart),hmixi(numpart),tri(numpart),ztemp(numpart)
  real :: masstemp(numpart,nspec)

#ifdef USE_NCF
  integer  :: ncid, mythread, thread_divide(12),mass_divide(nspec)
#endif

  ! Some variables needed for temporal interpolation
  !*************************************************
  call find_time_variables(itime)

!$OMP PARALLEL PRIVATE(i,tr,hm)
!$OMP DO
  do i=1,numpart
  ! Take only valid particles
  !**************************
    xlon(i)=-1.
    ylat(i)=-1.
    tti(i)=-1.
    rhoi(i)=-1.
    pri(i)=-1.
    pvi(i)=-1.
    qvi(i)=-1.
    topo(i)=-1.
    hmixi(i)=-1.
    tri(i)=-1.
    ztemp(i)=-1.
    do j=1,nspec
      masstemp(i,j)=-1.
    end do
    if (part(i)%alive) then
      xlon(i)=xlon0+part(i)%xlon*dx
      ylat(i)=ylat0+part(i)%ylat*dy

  !*****************************************************************************
  ! Interpolate several variables (PV, specific humidity, etc.) to particle position
  !*****************************************************************************
      call determine_grid_coordinates(real(part(i)%xlon),real(part(i)%ylat))
      call find_grid_distances(real(part(i)%xlon),real(part(i)%ylat))
  ! Topography
  !***********
      call bilinear_horizontal_interpolation_2dim(oro,topo(i))

      ! First set dz1out from interpol_mod to -1 so it only is calculated once per particle
      !************************************************************************************
      dz1out=-1
      ! Pressure
      call interpol_partoutput_value('PR',pri(i),i)
      ! Potential vorticity
      call interpol_partoutput_value('PV',pvi(i),i)
      ! Specific humidity
      call interpol_partoutput_value('QV',qvi(i),i)
      ! Temperature
      call interpol_partoutput_value('TT',tti(i),i)
      ! Density
      call interpol_partoutput_value('RH',rhoi(i),i)
      ! Reset dz1out
      !*************
      dz1out=-1

  ! Tropopause and PBL height
  !**************************
  ! Tropopause
      call bilinear_horizontal_interpolation(tropopause,tr,1,1)
      call temporal_interpolation(tr(1),tr(2),tri(i))
  ! PBL height
      call bilinear_horizontal_interpolation(hmix,hm,1,1)
      call temporal_interpolation(hm(1),hm(2),hmixi(i))


  ! Convert eta z coordinate to meters if necessary
  !************************************************
      call update_zeta_to_z(itime, i)
      ztemp(i)=part(i)%z

  ! Assign the masses
  !******************
      do j=1,nspec
        masstemp(i,j)=part(i)%mass(j)
      end do
    endif 
  end do

!$OMP END DO
!$OMP END PARALLEL
  if (numpart.gt.0) then
    write(*,*) 'topo: ', topo(1), 'z:', part(1)%zeta,part(1)%z
    write(*,*) 'xtra,xeta: ', part(1)%xlon
    write(*,*) 'ytra,yeta: ', part(1)%ylat
    write(*,*) 'mass: ', part(1)%mass(:)
    write(*,*) pvi(1),qvi(1),tti(1),rhoi(1),part(1)%alive,&
      count%alive,count%spawned,count%terminated
  endif

  ! Determine current calendar date, needed for the file name
  !**********************************************************

  jul=bdate+real(itime,kind=dp)/86400._dp
  call caldate(jul,jjjjmmdd,ihmmss)
  write(adate,'(i8.8)') jjjjmmdd
  write(atime,'(i6.6)') ihmmss

  if (lnetcdfout.eq.1) then
  ! open output file
    call open_partoutput_file(ncid)

    ! Dividing the openmp threads for writing
    j=0
    do i=1,11 !number of fields
      if (j.eq.numthreads) j = 0
      thread_divide(i) = j
      j = j + 1
    end do
    do i=1,nspec
      if (j.eq.numthreads) j = 0
      mass_divide(i) = j
      j = j + 1
    end do

    ! First allocate the time and particle dimention within the netcdf file
    call partoutput_netcdf(itime,xlon,'TI',j,ncid)
    call partoutput_netcdf(itime,xlon,'PA',j,ncid)

    ! Fill the fields in parallel
    if (numpart.gt.0) then
!$OMP PARALLEL PRIVATE(j,mythread)
#ifdef USE_NCF
      mythread = omp_get_thread_num()
      if (mythread.eq.thread_divide(1)) call partoutput_netcdf(itime,xlon,'LO',j,ncid)
      if (mythread.eq.thread_divide(2)) call partoutput_netcdf(itime,ylat,'LA',j,ncid)
      if (mythread.eq.thread_divide(3)) call partoutput_netcdf(itime,ztemp,'ZZ',j,ncid)
      if (mythread.eq.thread_divide(5)) call partoutput_netcdf(itime,pvi,'PV',j,ncid)
      if (mythread.eq.thread_divide(6)) call partoutput_netcdf(itime,qvi,'QV',j,ncid)
      if (mythread.eq.thread_divide(7)) call partoutput_netcdf(itime,rhoi,'RH',j,ncid)
      if (mythread.eq.thread_divide(10)) call partoutput_netcdf(itime,tti,'TT',j,ncid)
      if (mythread.eq.thread_divide(11)) call partoutput_netcdf(itime,pri,'PR',j,ncid)
      if (mythread.eq.thread_divide(4)) call partoutput_netcdf(itime,topo,'TO',j,ncid)
      if (mythread.eq.thread_divide(9)) call partoutput_netcdf(itime,tri,'TR',j,ncid)
      if (mythread.eq.thread_divide(8)) call partoutput_netcdf(itime,hmixi,'HM',j,ncid)
      if (mdomainfill.lt.1) then
        do j=1,nspec
          if (mythread.eq.mass_divide(j)) call partoutput_netcdf(itime,masstemp(:,j),'MA',j,ncid)
        end do
      endif
#endif
!$OMP END PARALLEL
    endif
    call close_partoutput_file(ncid)
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

      if (part(i)%alive) then
    ! Write the output
    !*****************      
        write(unitpartout) part(i)%npoint,xlon(i),ylat(i),part(i)%z, &
             part(i)%tstart,topo(i),pvi(i),qvi(i),rhoi(i),hmixi(i),tri(i),tti(i), &
             (part(i)%mass(j),j=1,nspec)
      endif
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
    if (iout.ne.0) call conccalc(itime,weight)
  endif

  ! If no grid is to be written to file, return
  !********************************************
  if (iout.eq.0) then 
    if (itime.ne.loutend) return
    loutnext=loutnext+loutstep
    loutstart=loutnext-loutaver/2
    loutend=loutnext+loutaver/2
    if (itime.eq.loutstart) then
      weight=0.5
      outnum=outnum+weight
    endif
    return
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
      ! call update_zeta_to_z(itime,i)
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
  if (numreceptor.eq.0) return

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

subroutine partpos_average(itime,j)

  !**********************************************************************
  ! This subroutine averages particle quantities, to be used for particle
  ! dump (in partoutput.f90). Averaging is done over output interval.
  !**********************************************************************

  use par_mod
  use com_mod
  use interpol_mod
  use coordinates_ecmwf

  implicit none

  integer :: itime,j
  real :: xlon,ylat,x,y,z
  real :: topo,hm(2),hmixi,pvi,qvi
  real :: tti,rhoi,ttemp
  real :: uui,vvi
  real :: tr(2),tri!,energy



 ! Some variables needed for temporal interpolation
  !*************************************************
  call find_time_variables(itime)

  xlon=xlon0+part(j)%xlon*dx
  ylat=ylat0+part(j)%ylat*dy

  !*****************************************************************************
  ! Interpolate several variables (PV, specific humidity, etc.) to particle position
  !*****************************************************************************

  call determine_grid_coordinates(real(part(j)%xlon),real(part(j)%ylat))
  call find_grid_distances(real(part(j)%xlon),real(part(j)%ylat))

  ! Topography
  !***********
  call bilinear_horizontal_interpolation_2dim(oro,topo)

  ! Potential vorticity, specific humidity, temperature, and density
  !*****************************************************************
  ! First set dz1out from interpol_mod to -1 so it only is calculated once per particle
  !************************************************************************************
  dz1out=-1
  ! Potential vorticity
  call interpol_partoutput_value('PV',pvi,j)
  ! Specific humidity
  call interpol_partoutput_value('QV',qvi,j)
  ! Temperature
  call interpol_partoutput_value('TT',tti,j)
  ! U wind
  call interpol_partoutput_value('UU',uui,j)
  ! V wind
  call interpol_partoutput_value('VV',vvi,j)
  ! Density
  call interpol_partoutput_value('RH',rhoi,j)
  ! Reset dz1out
  !*************
  dz1out=-1

  ! Convert eta z coordinate to meters if necessary. Can be moved to output only
  !************************************************
  call update_zeta_to_z(itime,j)
  
  ! energy=tti*cpa+(ztemp1+topo)*9.81+qvi*2501000.+(uui**2+vvi**2)/2.

  ! Add new values to sum and increase counter by one
  !**************************************************

  npart_av(j)=npart_av(j)+1

  ! Calculate Cartesian 3D coordinates suitable for averaging
  !**********************************************************

  xlon=xlon*pi180
  ylat=ylat*pi180
  x = cos(ylat)*sin(xlon)
  y = -1.*cos(ylat)*cos(xlon)
  z = sin(ylat)

  part_av_cartx(j)=part_av_cartx(j)+x
  part_av_carty(j)=part_av_carty(j)+y
  part_av_cartz(j)=part_av_cartz(j)+z
  part_av_z(j)=part_av_z(j)+part(j)%z
  part_av_pv(j)=part_av_pv(j)+pvi
  part_av_qv(j)=part_av_qv(j)+qvi
  part_av_tt(j)=part_av_tt(j)+tti
  part_av_uu(j)=part_av_uu(j)+uui
  part_av_vv(j)=part_av_vv(j)+vvi
  ! part_av_energy(j)=part_av_energy(j)+energy

  return
end subroutine partpos_average

end module output_mod