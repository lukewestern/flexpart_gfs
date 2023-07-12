! SPDX-FileCopyrightText: FLEXPART 1998-2019, see flexpart_license.txt
! SPDX-License-Identifier: GPL-3.0-or-later

!*****************************************************************************
!                                                                            *
!   L. Bakels 2022: This module contains most output related subroutines     *
!                                                                            *
!*****************************************************************************

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

subroutine init_output(itime,filesize)

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
  if (itime.eq.itime_init) then
    if (iout.ne.0) then ! No gridded output
#ifdef USE_NCF
      if (lnetcdfout.eq.1) then 
        call writeheader_netcdf(lnest=.false.)
      else 
        call writeheader_bin
      end if

      if (nested_output.eq.1) then
        if (lnetcdfout.eq.1) then
          call writeheader_netcdf(lnest=.true.)
        else if ((nested_output.eq.1).and.(surf_only.ne.1)) then
          call writeheader_bin_nest
        else if ((nested_output.eq.1).and.(surf_only.eq.1)) then
          call writeheader_bin_sfc_nest
        else if ((nested_output.ne.1).and.(surf_only.eq.1)) then 
          call writeheader_bin_sfc
        endif
      endif
#else
      call writeheader_bin

      !if (nested_output.eq.1) call writeheader_nest
      if ((nested_output.eq.1).and.(surf_only.ne.1)) call writeheader_bin_nest
      if ((nested_output.eq.1).and.(surf_only.eq.1)) call writeheader_bin_sfc_nest
      if ((nested_output.ne.1).and.(surf_only.eq.1)) call writeheader_bin_sfc
#endif
    endif ! iout.ne.0
    ! FLEXPART 9.2 ticket ?? write header in ASCII format 
    call writeheader_txt

    ! NetCDF only: Create file for storing initial particle positions.
#ifdef USE_NCF
    if (itime_init.ne.0) then
      jul=bdate+real(itime,kind=dp)/86400._dp
      call caldate(jul,jjjjmmdd,ihmmss)      
    endif
    if ((mdomainfill.eq.0).and.(ipout.ge.1).and.(ipin.le.1)) then
      if (itime_init.ne.0) then
        if (ldirect.eq.1) then
          call create_particles_initialoutput(ihmmss,jjjjmmdd,ibtime,ibdate)
        else
          call create_particles_initialoutput(ihmmss,jjjjmmdd,ietime,iedate)
        endif
      else if (ldirect.eq.1) then
        call create_particles_initialoutput(ibtime,ibdate,ibtime,ibdate)
      else
        call create_particles_initialoutput(ietime,iedate,ietime,iedate)
      endif
    endif
    ! Create header files for files that store the particle dump output
    if (ipout.ge.1) then
      if (itime_init.ne.0) then
        if (ldirect.eq.1) then
          call writeheader_partoutput(ihmmss,jjjjmmdd,ibtime,ibdate)
        else
          call writeheader_partoutput(ihmmss,jjjjmmdd,ietime,iedate)
        endif
      else if (ldirect.eq.1) then
        call writeheader_partoutput(ibtime,ibdate,ibtime,ibdate)
      else 
        call writeheader_partoutput(ietime,iedate,ietime,iedate)
      endif
    endif
#endif
  
  ! In case the particle output file is becoming larger than the maximum set
  ! in par_mod, create a new one while keeping track of the filesize.
  ! Also if a new restart file is created.
  else if ((mod(itime,ipoutfac*loutstep).eq.0).and.(ipout.ge.1)) then
#ifdef USE_NCF
    if ((filesize.ge.max_partoutput_filesize).or. &
      ((loutrestart.ne.-1).and.(mod(itime,loutrestart).eq.0))) then 
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
end subroutine init_output

subroutine finalise_output(itime)
  ! Complete the calculation of initial conditions for particles not yet terminated
  
  implicit none 

  integer, intent(in) :: itime
  integer :: j,ithread

  if (linit_cond.ge.1) then
    do j=1,numpart
      call initcond_calc(itime,j,1)
    end do
#ifdef _OPENMP
    do ithread=1,numthreads
      init_cond(:,:,:,:,:)=init_cond(:,:,:,:,:)+init_cond_omp(:,:,:,:,:,ithread)
    end do
#endif
  endif


  if (ipout.eq.2) call output_particles(itime)!,active_per_rel)     ! dump particle positions

  if (linit_cond.ge.1) then
    if(linversionout.eq.1) then
      call initcond_output_inv(itime)   ! dump initial cond. field
    else
      call initcond_output(itime)   ! dump initial cond. fielf
    endif
  endif
end subroutine finalise_output

subroutine output_particles(itime,initial_output)
  !                        i
  !*****************************************************************************
  !                                                                            *
  !     Dump all particle positions                                            *
  !     Author: A. Stohl                                                       *
  !                                                                            *
  !     12 March 1999                                                          *
  !                                                                            *
  !     Changes L. Bakels, 2021                                                *
  !     Output is chosen by the fields set in PARTOPTIONS                      *
  !     Binary output is no longer supported. If required, function can be     *
  !     added below at "Put binary function here"                              *
  !*****************************************************************************
  !                                                                            *
  ! Variables:                                                                 *
  !                                                                            *
  !*****************************************************************************

  use interpol_mod
  use coord_ecmwf_mod
  use particle_mod
#ifdef USE_NCF
  use netcdf
  use netcdf_output_mod, only: partoutput_netcdf,open_partoutput_file, &
                               close_partoutput_file,partinitpointer1
  use omp_lib, only: OMP_GET_THREAD_NUM
#endif

  implicit none

  integer,intent(in) :: itime
  logical,optional,intent(in) :: initial_output
  logical :: init_out
  integer :: i,j,m,jjjjmmdd,ihmmss,np,ns,i_av
  real(kind=dp) :: jul
  real :: tmp(2)
  character :: adate*8,atime*6

  real :: xlon(numpart),ylat(numpart),ztemp1,ztemp2,val_av(numpart,2),z_av(numpart)
  real :: tti(numpart),rhoi(numpart),pvi(numpart),qvi(numpart),pri(numpart)
  real :: topo(numpart),hmixi(numpart),tri(numpart),ztemp(numpart)
  real :: masstemp(numpart,nspec),masstemp_av(numpart,nspec)
  real :: wetdepotemp(numpart,nspec),drydepotemp(numpart,nspec)

  real :: output(num_partopt, numpart)

  ! For averaged output
  real :: xlon_av(numpart),ylat_av(numpart)

  real :: cartxyz(3)
  logical :: cartxyz_comp

#ifdef USE_NCF
  integer  :: ncid, mythread, thread_divide(12),mass_divide(nspec)
#else
  error stop 'NETCDF missing! Please compile with netcdf if you want the particle dump.'
#endif

#ifdef USE_NCF
  if (present(initial_output)) then
    init_out=initial_output
  else
    init_out=.false.
  endif

!$OMP PARALLEL PRIVATE(i,j,m,tmp,ns,i_av,cartxyz_comp,cartxyz,np)
  ! Some variables needed for temporal interpolation
  !*************************************************
  call find_time_vars(itime)

!$OMP DO
  do i=1,numpart
    if (((.not. part(i)%alive).and.(abs(part(i)%tend-itime).ge.ipoutfac*loutstep)) .or. &
      (init_out .and. (i.lt.partinitpointer1-1))) then ! Only freshly spawned particles need to be computed for init_out
      output(:,i) = -1
      masstemp(i,:) = -1
      masstemp_av(i,:) = -1
      wetdepotemp(i,:) = -1
      drydepotemp(i,:) = -1
      cycle
    endif
    !*****************************************************************************
    ! Interpolate several variables (PV, specific humidity, etc.) to particle position
    !*****************************************************************************
    ! Where in the grid? Stereographic (ngrid<0) or nested (ngrid>0)
    !***************************************************************
    call find_ngrid(real(part(i)%xlon),real(part(i)%ylat))
    call find_grid_indices(real(part(i)%xlon),real(part(i)%ylat))
    call find_grid_distances(real(part(i)%xlon),real(part(i)%ylat))
    ! First set dz1out from interpol_mod to -1 so it only is calculated once per particle
    !************************************************************************************
    dz1out=-1
    cartxyz_comp=.false.
    do np=1,num_partopt
      if (.not. partopt(np)%print) cycle ! Only compute when field should be printed
      i_av = partopt(np)%i_average
      if (init_out.and.(i_av.ne.0)) cycle ! no averages for initial particle output
      if ((i_av.ne.0).and.(part(i)%ntime.eq.0)) then
        if (partopt(np)%name.eq.'ma') then
          masstemp_av(i,1:nspec) = -1
        else
          output(np,i) = -1
        endif
        cycle ! no averages for freshly spawned particles
      endif
      select case (partopt(np)%name)
        case ('LO')
          output(np,i)=xlon0+part(i)%xlon*dx
          cycle
        case ('LA')
          output(np,i)=ylat0+part(i)%ylat*dy
          cycle
        case ('TO') ! Topography
          if (ngrid.le.0) then
            call hor_interpol(oro,output(np,i))
          else
            call hor_interpol_nest(oron,output(np,i))
          endif 
          cycle
        case ('TR') ! Tropopause
          if (ngrid.le.0) then
            do m=1,2
              call hor_interpol(tropopause,tmp(m),1,memind(m),1)
            end do
          else
            do m=1,2
              call hor_interpol_nest(tropopausen,tmp(m),1,memind(m),1)
            end do
          endif
          call temporal_interpolation(tmp(1),tmp(2),output(np,i))
          cycle
        case ('HM') ! PBL height
          if (ngrid.le.0) then
            do m=1,2
              call hor_interpol(hmix,tmp(m),1,memind(m),1)
            end do
          else
            do m=1,2
              call hor_interpol_nest(hmixn,tmp(m),1,memind(m),1)
            end do
          endif
          call temporal_interpolation(tmp(1),tmp(2),output(np,i))
          cycle
        case ('ZZ') ! Height
          call update_zeta_to_z(itime, i) ! Convert eta z coordinate to meters if necessary
          output(np,i)=part(i)%z
          cycle
        ! case ('UU') ! Longitudinal velocity
        !   output(np,i)=part(i)%vel%u !This would be preferred, but not implemented yet
        !   cycle
        case ('VS') ! Settling velocity
          output(np,i)=part(i)%settling
          cycle
        case ('MA') ! Mass
          do ns=1,nspec
            masstemp(i,ns)=part(i)%mass(ns)
          end do
          cycle
        case ('ma') ! Mass averaged
          do ns=1,nspec
            masstemp_av(i,ns)=part(i)%val_av(i_av+(ns-1))/part(i)%ntime
          end do
          cycle
        case ('WD') ! Wet deposition
          do ns=1,nspec
            wetdepotemp(i,ns)=part(i)%wetdepo(ns)
          end do
          cycle
        case ('DD') ! dry deposition
          do ns=1,nspec
            drydepotemp(i,ns)=part(i)%drydepo(ns)
          end do
          cycle
        case ('lo')
          if (.not. cartxyz_comp) then
            cartxyz(1) = part(i)%cartx_av/part(i)%ntime
            cartxyz(2) = part(i)%carty_av/part(i)%ntime
            cartxyz(3) = part(i)%cartz_av/part(i)%ntime
            cartxyz_comp=.true.
          endif
          output(np,i) = atan2(cartxyz(1),-1.*cartxyz(2))/pi180
          if (output(np,i).gt.360.) output(np,i)=output(np,i)-360.
          if (output(np,i).lt.0.) output(np,i)=output(np,i)+360.
          cycle
        case ('la')
          if (.not. cartxyz_comp) then
            cartxyz(1) = part(i)%cartx_av/part(i)%ntime
            cartxyz(2) = part(i)%carty_av/part(i)%ntime
            cartxyz(3) = part(i)%cartz_av/part(i)%ntime
            cartxyz_comp=.true.
          endif
          output(np,i) = atan2(cartxyz(3),sqrt(cartxyz(1)*cartxyz(1)+ &
            cartxyz(2)*cartxyz(2)))/pi180
        case default
          if (.not. partopt(np)%average) then
            call interpol_partoutput_val(partopt(np)%name,output(np,i),i)
          else
            output(np,i) = part(i)%val_av(i_av)/part(i)%ntime
          endif
      end select
    end do
    ! Reset dz1out
    !*************
    dz1out=-1
    cartxyz_comp=.false.

    if ((.not. init_out).and.(n_average.gt.0)) then
      part(i)%val_av = 0.
      part(i)%ntime = 0.
      part(i)%cartx_av = 0.
      part(i)%carty_av = 0.
      part(i)%cartz_av = 0.
    endif
  end do

!$OMP END DO
!$OMP END PARALLEL

  if ((.not. init_out).and.(numpart.gt.0)) then
    do np=1,num_partopt
      if (.not. partopt(np)%print) cycle
      if (partopt(np)%name.eq.'MA') then
        write(*,*) partopt(np)%long_name, masstemp(1,:)
      else if (partopt(np)%name.eq.'ma') then
        write(*,*) partopt(np)%long_name, masstemp_av(1,:)
      else if (partopt(np)%name.eq.'WD') then
        write(*,*) partopt(np)%long_name, wetdepotemp(1,:)
      else if (partopt(np)%name.eq.'DD') then
        write(*,*) partopt(np)%long_name, drydepotemp(1,:)
      else
        write(*,*) partopt(np)%long_name, output(np,1)
      endif
    end do
    write(*,*) part(1)%prob,part(1)%alive
    write(*,*) 'Alive: ', count%alive, 'Total spawned: ', count%spawned, 'Terminated: ', count%terminated
  endif

  ! Determine current calendar date, needed for the file name
  !**********************************************************

  jul=bdate+real(itime,kind=dp)/86400._dp
  call caldate(jul,jjjjmmdd,ihmmss)
  write(adate,'(i8.8)') jjjjmmdd
  write(atime,'(i6.6)') ihmmss
  j=1
  if (lnetcdfout.eq.1) then
  ! open output file
    if (init_out) then
      call open_partinit_file(ncid)
    else
      call open_partoutput_file(ncid)

      ! First allocate the time and particle dimensions within the netcdf file
      call partoutput_netcdf(itime,xlon,'TI',j,ncid)
      call partoutput_netcdf(itime,xlon,'PA',j,ncid)
    endif

    ! Fill the fields in parallel
    if (numpart.gt.0) then
    ! OpenMP output does not work on all systems depending on how they are set-up
! !$OMP PARALLEL PRIVATE(np,ns)
! !$OMP DO SCHEDULE(dynamic)
      do np=1,num_partopt
        !write(*,*) partopt(np)%name, output(np,1)
        if (.not. partopt(np)%print) cycle
        if (init_out.and.(partopt(np)%i_average.ne.0)) cycle ! no averages for initial particle output
        !write(*,*) partopt(np)%name
        if (partopt(np)%name.eq.'MA') then
          do ns=1,nspec
            if (init_out) then
              call partinit_netcdf(itime,masstemp(:,ns),'MA',ns,ncid)
            else
              call partoutput_netcdf(itime,masstemp(:,ns),'MA',ns,ncid)
            endif
          end do
        else if (partopt(np)%name.eq.'ma') then
          do ns=1,nspec
            call partoutput_netcdf(itime,masstemp_av(:,ns),'ma',ns,ncid)
          end do
        else if ((.not. init_out).and.(partopt(np)%name.eq.'WD')) then
          do ns=1,nspec
            call partoutput_netcdf(itime,wetdepotemp(:,ns),'WD',ns,ncid)
          end do
        else if ((.not. init_out).and.(partopt(np)%name.eq.'DD')) then
          do ns=1,nspec
            call partoutput_netcdf(itime,drydepotemp(:,ns),'DD',ns,ncid)
          end do
        else
          if (init_out) then
            call partinit_netcdf(itime,output(np,:),partopt(np)%name,j,ncid)
          else
            call partoutput_netcdf(itime,output(np,:),partopt(np)%name,j,ncid)
          endif
        endif
      end do
! !$OMP END DO
! !$OMP END PARALLEL
    endif
    call close_partoutput_file(ncid)
    if (.not. init_out) then
      mass_written=.true. ! needs to be reduced within openmp loop
      topo_written=.true. ! same
    endif
  else
    ! Put binary function here
  endif
#else
    ! Put binary function here
#endif
end subroutine output_particles

subroutine output_conc(itime,loutstart,loutend,loutnext,outnum)
  use unc_mod
  use outgrid_mod
  use par_mod
  use com_mod
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

  ! If no grid is to be written to file, return (LB)
  !*************************************************
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
      call concoutput_sfc_netcdf(itime,outnum,gridtotalunc,wetgridtotalunc,drygridtotalunc)
#else
      if (linversionout.eq.1) then
        call concoutput_inversion(itime,outnum,gridtotalunc,wetgridtotalunc,drygridtotalunc)
      else
        call concoutput_sfc(itime,outnum,gridtotalunc,wetgridtotalunc,drygridtotalunc)
      endif
#endif
    endif

    if (nested_output .eq. 1) then
#ifdef USE_NCF
      if (surf_only.ne.1) then
        call concoutput_nest_netcdf(itime,outnum)
      else 
        call concoutput_sfc_nest_netcdf(itime,outnum)
      endif
#else
      if (surf_only.ne.1) then
        call concoutput_nest(itime,outnum)
      else 
        if(linversionout.eq.1) then
          call concoutput_inversion_nest(itime,outnum)
        else 
          call concoutput_sfc_nest(itime,outnum)
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
end subroutine output_conc

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
  !     2021, LB: OpenMP parallelisation                                       *
  !                                                                            *
  !*****************************************************************************
  !                                                                            *
  ! Variables:                                                                 *
  ! nspeciesdim     = nspec for forward runs, 1 for backward runs              *
  !                                                                            *
  !*****************************************************************************

  use unc_mod
  use outgrid_mod
  use par_mod
  use com_mod
  use omp_lib, only: OMP_GET_THREAD_NUM
  use interpol_mod, only: interpol_density
  use coord_ecmwf_mod
  use particle_mod

  implicit none

  integer,intent(in) :: itime
  real,intent(in) :: weight
  integer :: itage,i,j,kz,ks,n,nage,inage,thread,ithread
  integer :: il,ind,indz,indzp,nrelpointer
  integer :: ix,jy,ixp,jyp
  real :: ddx,ddy
  real(kind=dp) :: mm3
  real :: hx,hy,hz,hxyz,xd,yd,zd,xkern,r2,c(maxspec)
  real :: rhoi
  real :: xl,yl,wx,wy,w
  real,parameter :: factor=.596831, hxmax=6.0, hymax=4.0, hzmax=150.
  !  integer xscav_count

  ! For forward simulations, make a loop over the number of species;
  ! for backward simulations, make an additional loop over the
  ! releasepoints
  !***************************************************************************
  !  xscav_count=0
#ifdef _OPENMP
  call omp_set_num_threads(numthreads_grid)
#endif
!$OMP PARALLEL PRIVATE(i,itage,nage,inage,rhoi,nrelpointer,kz,xl,yl,ks,wx,wy,w,thread,ddx,ddy, &
!$OMP ix,jy,ixp,jyp)
#if (defined _OPENMP)
    thread = OMP_GET_THREAD_NUM()+1 ! Starts with 1
#else
    thread = 1
#endif

!$OMP DO
  do j=1,count%alive

    i=count%ialive(j)

  ! Determine age class of the particle
    itage=abs(itime-part(i)%tstart)
    nage=1
    do inage=1,nageclass
      nage=inage
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
      call update_zeta_to_z(itime,i)
      call interpol_density(itime,i,rhoi)
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
#ifdef _OPENMP
              gridunc_omp(ix,jy,kz,ks,nrelpointer,part(i)%nclass,nage,thread)= &
                   gridunc_omp(ix,jy,kz,ks,nrelpointer,part(i)%nclass,nage,thread)+ &
                   part(i)%mass(ks)/rhoi*weight*max(xscav_frac1(i,ks),0.0)
#else
              gridunc(ix,jy,kz,ks,nrelpointer,part(i)%nclass,nage)= &
                   gridunc(ix,jy,kz,ks,nrelpointer,part(i)%nclass,nage)+ &
                   part(i)%mass(ks)/rhoi*weight*max(xscav_frac1(i,ks),0.0)
#endif
            end do
          else
            if (lparticlecountoutput) then
              do ks=1,nspec
#ifdef _OPENMP
                gridunc_omp(ix,jy,kz,ks,nrelpointer,part(i)%nclass,nage,thread)= &
                     gridunc_omp(ix,jy,kz,ks,nrelpointer,part(i)%nclass,nage,thread)+1
#else
                gridunc(ix,jy,kz,ks,nrelpointer,part(i)%nclass,nage)= &
                     gridunc(ix,jy,kz,ks,nrelpointer,part(i)%nclass,nage)+1
#endif
              end do
            else
              do ks=1,nspec
#ifdef _OPENMP
                gridunc_omp(ix,jy,kz,ks,nrelpointer,part(i)%nclass,nage,thread)= &
                     gridunc_omp(ix,jy,kz,ks,nrelpointer,part(i)%nclass,nage,thread)+ &
                     part(i)%mass(ks)/rhoi*weight
#else
                gridunc(ix,jy,kz,ks,nrelpointer,part(i)%nclass,nage)= &
                     gridunc(ix,jy,kz,ks,nrelpointer,part(i)%nclass,nage)+ &
                     part(i)%mass(ks)/rhoi*weight
#endif
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
#ifdef _OPENMP
                 gridunc_omp(ix,jy,kz,ks,nrelpointer,part(i)%nclass,nage,thread)= &
                   gridunc_omp(ix,jy,kz,ks,nrelpointer,part(i)%nclass,nage,thread)+ &
                   part(i)%mass(ks)/rhoi*w*weight*max(xscav_frac1(i,ks),0.0)
#else
                 gridunc(ix,jy,kz,ks,nrelpointer,part(i)%nclass,nage)= &
                   gridunc(ix,jy,kz,ks,nrelpointer,part(i)%nclass,nage)+ &
                   part(i)%mass(ks)/rhoi*w*weight*max(xscav_frac1(i,ks),0.0)
#endif
               end do
            else
               do ks=1,nspec
#ifdef _OPENMP
                 gridunc_omp(ix,jy,kz,ks,nrelpointer,part(i)%nclass,nage,thread)= &
                   gridunc_omp(ix,jy,kz,ks,nrelpointer,part(i)%nclass,nage,thread)+ &
                   part(i)%mass(ks)/rhoi*weight*w
#else
                 gridunc(ix,jy,kz,ks,nrelpointer,part(i)%nclass,nage)= &
                   gridunc(ix,jy,kz,ks,nrelpointer,part(i)%nclass,nage)+ &
                   part(i)%mass(ks)/rhoi*weight*w
#endif
               end do
            endif
          endif

          if ((jyp.ge.0).and.(jyp.le.numygrid-1)) then
            w=wx*(1.-wy)
            if (DRYBKDEP.or.WETBKDEP) then
              do ks=1,nspec
#ifdef _OPENMP
                 gridunc_omp(ix,jyp,kz,ks,nrelpointer,part(i)%nclass,nage,thread)= &
                   gridunc_omp(ix,jyp,kz,ks,nrelpointer,part(i)%nclass,nage,thread)+ &
                   part(i)%mass(ks)/rhoi*weight*w*max(xscav_frac1(i,ks),0.0)
#else
                 gridunc(ix,jyp,kz,ks,nrelpointer,part(i)%nclass,nage)= &
                   gridunc(ix,jyp,kz,ks,nrelpointer,part(i)%nclass,nage)+ &
                   part(i)%mass(ks)/rhoi*weight*w*max(xscav_frac1(i,ks),0.0)
#endif
               end do
             else
              do ks=1,nspec
#ifdef _OPENMP
                 gridunc_omp(ix,jyp,kz,ks,nrelpointer,part(i)%nclass,nage,thread)= &
                   gridunc_omp(ix,jyp,kz,ks,nrelpointer,part(i)%nclass,nage,thread)+ &
                   part(i)%mass(ks)/rhoi*weight*w
#else
                 gridunc(ix,jyp,kz,ks,nrelpointer,part(i)%nclass,nage)= &
                   gridunc(ix,jyp,kz,ks,nrelpointer,part(i)%nclass,nage)+ &
                   part(i)%mass(ks)/rhoi*weight*w
#endif
               end do
             endif
          endif
        endif !ix ge 0


        if ((ixp.ge.0).and.(ixp.le.numxgrid-1)) then
          if ((jyp.ge.0).and.(jyp.le.numygrid-1)) then
            w=(1.-wx)*(1.-wy)
            if (DRYBKDEP.or.WETBKDEP) then
               do ks=1,nspec
#ifdef _OPENMP
                 gridunc_omp(ixp,jyp,kz,ks,nrelpointer,part(i)%nclass,nage,thread)= &
                   gridunc_omp(ixp,jyp,kz,ks,nrelpointer,part(i)%nclass,nage,thread)+ &
                   part(i)%mass(ks)/rhoi*w*weight*max(xscav_frac1(i,ks),0.0)
#else
                 gridunc(ixp,jyp,kz,ks,nrelpointer,part(i)%nclass,nage)= &
                   gridunc(ixp,jyp,kz,ks,nrelpointer,part(i)%nclass,nage)+ &
                   part(i)%mass(ks)/rhoi*w*weight*max(xscav_frac1(i,ks),0.0)
#endif
               end do
            else
               do ks=1,nspec
#ifdef _OPENMP
                 gridunc_omp(ixp,jyp,kz,ks,nrelpointer,part(i)%nclass,nage,thread)= &
                   gridunc_omp(ixp,jyp,kz,ks,nrelpointer,part(i)%nclass,nage,thread)+ &
                   part(i)%mass(ks)/rhoi*weight*w
#else
                 gridunc(ixp,jyp,kz,ks,nrelpointer,part(i)%nclass,nage)= &
                   gridunc(ixp,jyp,kz,ks,nrelpointer,part(i)%nclass,nage)+ &
                   part(i)%mass(ks)/rhoi*weight*w
#endif
               end do
            endif
          endif

          if ((jy.ge.0).and.(jy.le.numygrid-1)) then
            w=(1.-wx)*wy
            if (DRYBKDEP.or.WETBKDEP) then
               do ks=1,nspec
#ifdef _OPENMP
                 gridunc_omp(ixp,jy,kz,ks,nrelpointer,part(i)%nclass,nage,thread)= &
                   gridunc_omp(ixp,jy,kz,ks,nrelpointer,part(i)%nclass,nage,thread)+ &
                   part(i)%mass(ks)/rhoi*weight*w*max(xscav_frac1(i,ks),0.0)
#else
                 gridunc(ixp,jy,kz,ks,nrelpointer,part(i)%nclass,nage)= &
                   gridunc(ixp,jy,kz,ks,nrelpointer,part(i)%nclass,nage)+ &
                   part(i)%mass(ks)/rhoi*weight*w*max(xscav_frac1(i,ks),0.0)
#endif
               end do
            else
               do ks=1,nspec
#ifdef _OPENMP
                 gridunc_omp(ixp,jy,kz,ks,nrelpointer,part(i)%nclass,nage,thread)= &
                   gridunc_omp(ixp,jy,kz,ks,nrelpointer,part(i)%nclass,nage,thread)+ &
                   part(i)%mass(ks)/rhoi*weight*w
#else
                 gridunc(ixp,jy,kz,ks,nrelpointer,part(i)%nclass,nage)= &
                   gridunc(ixp,jy,kz,ks,nrelpointer,part(i)%nclass,nage)+ &
                   part(i)%mass(ks)/rhoi*weight*w
#endif
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
#ifdef _OPENMP
                 griduncn_omp(ix,jy,kz,ks,nrelpointer,part(i)%nclass,nage,thread)= &
                   griduncn_omp(ix,jy,kz,ks,nrelpointer,part(i)%nclass,nage,thread)+ &
                   part(i)%mass(ks)/rhoi*weight*max(xscav_frac1(i,ks),0.0)
#else            
                 griduncn(ix,jy,kz,ks,nrelpointer,part(i)%nclass,nage)= &
                   griduncn(ix,jy,kz,ks,nrelpointer,part(i)%nclass,nage)+ &
                   part(i)%mass(ks)/rhoi*weight*max(xscav_frac1(i,ks),0.0)
#endif
               end do
            else
              if (lparticlecountoutput) then
                do ks=1,nspec
#ifdef _OPENMP
                  griduncn_omp(ix,jy,kz,ks,nrelpointer,part(i)%nclass,nage,thread)= &
                       griduncn_omp(ix,jy,kz,ks,nrelpointer,part(i)%nclass,nage,thread)+1
#else  
                  griduncn(ix,jy,kz,ks,nrelpointer,part(i)%nclass,nage)= &
                       griduncn(ix,jy,kz,ks,nrelpointer,part(i)%nclass,nage)+1
#endif
                end do
              else
                do ks=1,nspec
#ifdef _OPENMP
                  griduncn_omp(ix,jy,kz,ks,nrelpointer,part(i)%nclass,nage,thread)= &
                       griduncn_omp(ix,jy,kz,ks,nrelpointer,part(i)%nclass,nage,thread)+ &
                       part(i)%mass(ks)/rhoi*weight
#else            
                  griduncn(ix,jy,kz,ks,nrelpointer,part(i)%nclass,nage)= &
                       griduncn(ix,jy,kz,ks,nrelpointer,part(i)%nclass,nage)+ &
                       part(i)%mass(ks)/rhoi*weight
#endif
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
#ifdef _OPENMP
                   griduncn_omp(ix,jy,kz,ks,nrelpointer,part(i)%nclass,nage,thread)= &
                     griduncn_omp(ix,jy,kz,ks,nrelpointer,part(i)%nclass,nage,thread)+ &
                     part(i)%mass(ks)/rhoi*weight*w*max(xscav_frac1(i,ks),0.0)
#else              
                   griduncn(ix,jy,kz,ks,nrelpointer,part(i)%nclass,nage)= &
                     griduncn(ix,jy,kz,ks,nrelpointer,part(i)%nclass,nage)+ &
                     part(i)%mass(ks)/rhoi*weight*w*max(xscav_frac1(i,ks),0.0)
#endif
                 end do
              else
                do ks=1,nspec
#ifdef _OPENMP
                   griduncn_omp(ix,jy,kz,ks,nrelpointer,part(i)%nclass,nage,thread)= &
                     griduncn_omp(ix,jy,kz,ks,nrelpointer,part(i)%nclass,nage,thread)+ &
                     part(i)%mass(ks)/rhoi*weight*w
#else            
                   griduncn(ix,jy,kz,ks,nrelpointer,part(i)%nclass,nage)= &
                     griduncn(ix,jy,kz,ks,nrelpointer,part(i)%nclass,nage)+ &
                     part(i)%mass(ks)/rhoi*weight*w
#endif
                 end do
              endif
            endif

            if ((jyp.ge.0).and.(jyp.le.numygridn-1)) then
              w=wx*(1.-wy)
              if (DRYBKDEP.or.WETBKDEP) then
                 do ks=1,nspec
#ifdef _OPENMP
                   griduncn_omp(ix,jyp,kz,ks,nrelpointer,part(i)%nclass,nage,thread)= &
                     griduncn_omp(ix,jyp,kz,ks,nrelpointer,part(i)%nclass,nage,thread)+ &
                     part(i)%mass(ks)/rhoi*weight*w*max(xscav_frac1(i,ks),0.0)
#else              
                   griduncn(ix,jyp,kz,ks,nrelpointer,part(i)%nclass,nage)= &
                     griduncn(ix,jyp,kz,ks,nrelpointer,part(i)%nclass,nage)+ &
                     part(i)%mass(ks)/rhoi*weight*w*max(xscav_frac1(i,ks),0.0)
#endif
                 end do
              else
                 do ks=1,nspec
#ifdef _OPENMP
                   griduncn_omp(ix,jyp,kz,ks,nrelpointer,part(i)%nclass,nage,thread)= &
                     griduncn_omp(ix,jyp,kz,ks,nrelpointer,part(i)%nclass,nage,thread)+ &
                     part(i)%mass(ks)/rhoi*weight*w
#else              
                   griduncn(ix,jyp,kz,ks,nrelpointer,part(i)%nclass,nage)= &
                     griduncn(ix,jyp,kz,ks,nrelpointer,part(i)%nclass,nage)+ &
                     part(i)%mass(ks)/rhoi*weight*w
#endif
                 end do
              endif
            endif
          endif


          if ((ixp.ge.0).and.(ixp.le.numxgridn-1)) then
            if ((jyp.ge.0).and.(jyp.le.numygridn-1)) then
              w=(1.-wx)*(1.-wy)
              if (DRYBKDEP.or.WETBKDEP) then
                 do ks=1,nspec
#ifdef _OPENMP
                   griduncn_omp(ixp,jyp,kz,ks,nrelpointer,part(i)%nclass,nage,thread)= &
                     griduncn_omp(ixp,jyp,kz,ks,nrelpointer,part(i)%nclass,nage,thread)+ &
                     part(i)%mass(ks)/rhoi*weight*w*max(xscav_frac1(i,ks),0.0)
#else              
                   griduncn(ixp,jyp,kz,ks,nrelpointer,part(i)%nclass,nage)= &
                     griduncn(ixp,jyp,kz,ks,nrelpointer,part(i)%nclass,nage)+ &
                     part(i)%mass(ks)/rhoi*weight*w*max(xscav_frac1(i,ks),0.0)
#endif
                 end do
              else
                 do ks=1,nspec
#ifdef _OPENMP
                   griduncn_omp(ixp,jyp,kz,ks,nrelpointer,part(i)%nclass,nage,thread)= &
                     griduncn_omp(ixp,jyp,kz,ks,nrelpointer,part(i)%nclass,nage,thread)+ &
                     part(i)%mass(ks)/rhoi*weight*w
#else              
                   griduncn(ixp,jyp,kz,ks,nrelpointer,part(i)%nclass,nage)= &
                     griduncn(ixp,jyp,kz,ks,nrelpointer,part(i)%nclass,nage)+ &
                     part(i)%mass(ks)/rhoi*weight*w
#endif
                 end do
              endif
            endif

            if ((jy.ge.0).and.(jy.le.numygridn-1)) then
              w=(1.-wx)*wy
              if (DRYBKDEP.or.WETBKDEP) then
                 do ks=1,nspec
#ifdef _OPENMP
                   griduncn_omp(ixp,jy,kz,ks,nrelpointer,part(i)%nclass,nage,thread)= &
                     griduncn_omp(ixp,jy,kz,ks,nrelpointer,part(i)%nclass,nage,thread)+ &
                     part(i)%mass(ks)/rhoi*weight*w*max(xscav_frac1(i,ks),0.0)
#else              
                   griduncn(ixp,jy,kz,ks,nrelpointer,part(i)%nclass,nage)= &
                     griduncn(ixp,jy,kz,ks,nrelpointer,part(i)%nclass,nage)+ &
                     part(i)%mass(ks)/rhoi*weight*w*max(xscav_frac1(i,ks),0.0)
#endif
                 end do
              else
                 do ks=1,nspec
#ifdef _OPENMP
                    griduncn_omp(ixp,jy,kz,ks,nrelpointer,part(i)%nclass,nage,thread)= &
                     griduncn_omp(ixp,jy,kz,ks,nrelpointer,part(i)%nclass,nage,thread)+ &
                     part(i)%mass(ks)/rhoi*weight*w
#else              
                    griduncn(ixp,jy,kz,ks,nrelpointer,part(i)%nclass,nage)= &
                     griduncn(ixp,jy,kz,ks,nrelpointer,part(i)%nclass,nage)+ &
                     part(i)%mass(ks)/rhoi*weight*w
#endif
                 end do
              endif
            endif
          endif
        endif
      endif
    endif
  end do
!$OMP END DO
!$OMP END PARALLEL
#ifdef _OPENMP
  call omp_set_num_threads(numthreads)
#endif
  ! Reduction of gridunc and griduncn
#ifdef _OPENMP
  do ithread=1,numthreads_grid
    gridunc(:,:,:,:,:,:,:)=gridunc(:,:,:,:,:,:,:)+gridunc_omp(:,:,:,:,:,:,:,ithread)
    gridunc_omp(:,:,:,:,:,:,:,ithread)=0.
  end do
  if (nested_output.eq.1) then 
    do ithread=1,numthreads_grid
      griduncn(:,:,:,:,:,:,:)=griduncn(:,:,:,:,:,:,:)+griduncn_omp(:,:,:,:,:,:,:,ithread)
      griduncn_omp(:,:,:,:,:,:,:,ithread)=0.
    end do
  endif
#endif

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

    do j=1,count%alive

      i=count%ialive(j)

      itage=abs(itime-part(i)%tstart)

      hz=min(50.+0.3*sqrt(real(itage)),hzmax)
      zd=real(part(i)%z)/hz
      if (zd.gt.1.) cycle          ! save computing time, leave loop

      hx=min((0.29+2.222e-3*sqrt(real(itage)))*dx+ &
           real(itage)*1.2e-5,hxmax)                     ! 80 km/day
      xd=(real(part(i)%xlon)-xreceptor(n))/hx
      if (xd*xd.gt.1.) cycle       ! save computing time, leave loop

      hy=min((0.18+1.389e-3*sqrt(real(itage)))*dy+ &
           real(itage)*7.5e-6,hymax)                     ! 80 km/day
      yd=(real(part(i)%ylat)-yreceptor(n))/hy
      if (yd*yd.gt.1.) cycle       ! save computing time, leave loop
      hxyz=hx*hy*hz

      r2=xd*xd+yd*yd+zd*zd
      if (r2.lt.1.) then
        xkern=factor*(1.-r2)
        do ks=1,nspec
          c(ks)=c(ks)+part(i)%mass(ks)*xkern/hxyz
        end do
      endif
    end do

    do ks=1,nspec
      creceptor(n,ks)=creceptor(n,ks)+2.*weight*c(ks)/receptorarea(n)
    end do
  end do
end subroutine conccalc

subroutine partpos_avg(itime,j)

  !**********************************************************************
  ! This subroutine averages particle quantities, to be used for particle
  ! dump (in partoutput.f90). Averaging is done over output interval.
  ! Author: A. Stohl
  ! Changes L Bakels:
  !    - Computing fields defined in PARTOPTIONS
  !**********************************************************************

  use par_mod
  use com_mod
  use interpol_mod
  use coord_ecmwf_mod

  implicit none

  integer,intent(in) :: itime,j
  integer :: np,i_av,ns,m
  real :: xlon,ylat,x,y,z
  real :: topo,hm(2),hmixi,pvi,qvi
  real :: tti,rhoi,ttemp
  real :: uui,vvi,output
  real :: tr(2),tri!,energy

  logical :: cart_comp

  if (ipout.eq.0) return ! No need to compute averages since there is no particle output

  if (n_average.eq.0) return

  if (.not. part(j)%alive) return

  if (part(j)%nstop) return ! If particle is to be killed, averages cannot be computed

 ! Some variables needed for temporal interpolation
  !*************************************************
  call find_time_vars(itime)

  xlon=xlon0+real(part(j)%xlon)*dx
  ylat=ylat0+real(part(j)%ylat)*dy

  !*****************************************************************************
  ! Interpolate several variables (PV, specific humidity, etc.) to particle position
  !*****************************************************************************
  ! Where in the grid? Stereographic (ngrid<0) or nested (ngrid>0)
  !***************************************************************
  call find_ngrid(real(part(j)%xlon),real(part(j)%ylat))
  call find_grid_indices(real(part(j)%xlon),real(part(j)%ylat))
  call find_grid_distances(real(part(j)%xlon),real(part(j)%ylat))

  ! First set dz1out from interpol_mod to -1 so it only is calculated once per particle
  !************************************************************************************
  part(j)%ntime=part(j)%ntime + 1
  dz1out=-1
  cart_comp=.false.
  do np=1,num_partopt
    if ((.not. partopt(np)%print) .or. (.not. partopt(np)%average)) cycle
    i_av = partopt(np)%i_average
    select case (partopt(np)%name)
      case ('to')
        if (ngrid.le.0) then
          call hor_interpol(oro,output)
        else
          call hor_interpol_nest(oron,output)
        endif
        part(j)%val_av(i_av)=part(j)%val_av(i_av)+output
      case ('tr')
        if (ngrid.le.0) then
          do m=1,2
            call hor_interpol(tropopause,tr(m),1,memind(m),1)
          end do
        else
          do m=1,2
            call hor_interpol_nest(tropopausen,tr(m),1,memind(m),1)
          end do
        endif
        call temporal_interpolation(tr(1),tr(2),output)
        part(j)%val_av(i_av)=part(j)%val_av(i_av)+output
      case ('hm')
        if (ngrid.le.0) then
          do m=1,2
            call hor_interpol(hmix,hm(m),1,memind(m),1)
          end do
        else
          do m=1,2
            call hor_interpol_nest(hmixn,hm(m),1,memind(m),1)
          end do
        endif
        call temporal_interpolation(hm(1),hm(2),output)
        part(j)%val_av(i_av)=part(j)%val_av(i_av)+output
      case ('lo')
        if (.not. cart_comp) then
          ! Calculate Cartesian 3D coordinates suitable for averaging
          !**********************************************************

          xlon=xlon*pi180
          ylat=ylat*pi180
          x = cos(ylat)*sin(xlon)
          y = -1.*cos(ylat)*cos(xlon)
          z = sin(ylat)

          part(j)%cartx_av=part(j)%cartx_av+x
          part(j)%carty_av=part(j)%carty_av+y
          part(j)%cartz_av=part(j)%cartz_av+z
          cart_comp=.true.
        endif
      case ('la')
        if (.not. cart_comp) then
          ! Calculate Cartesian 3D coordinates suitable for averaging
          !**********************************************************

          xlon=xlon*pi180
          ylat=ylat*pi180
          x = cos(ylat)*sin(xlon)
          y = -1.*cos(ylat)*cos(xlon)
          z = sin(ylat)

          part(j)%cartx_av=part(j)%cartx_av+x
          part(j)%carty_av=part(j)%carty_av+y
          part(j)%cartz_av=part(j)%cartz_av+z
          cart_comp=.true.
        endif
      case ('zz')
        ! Convert eta z coordinate to meters if necessary. Can be moved to output only
        !************************************************
        call update_zeta_to_z(itime,j)
        part(j)%val_av(i_av)=part(j)%val_av(i_av)+part(j)%z
      case ('ma')
        do ns=1,nspec
          part(j)%val_av(i_av+(ns-1))=part(j)%val_av(i_av+(ns-1))+part(j)%mass(ns)
        end do
      case ('vs')
        part(j)%val_av(i_av)=part(j)%val_av(i_av)+part(j)%settling
      case default
        call interpol_partoutput_val(partopt(np)%name,output,j)
        part(j)%val_av(i_av)=part(j)%val_av(i_av)+output
    end select
  end do
  ! Reset dz1out
  !*************
  dz1out=-1
  cart_comp=.false.

  return
end subroutine partpos_avg

end module output_mod
