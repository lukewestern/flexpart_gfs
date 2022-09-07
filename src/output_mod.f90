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

  character(len=256) :: restart_filename1,restart_filename2,restart_filename3
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
  if (itime.eq.itime_init) then
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
    if (itime_init.ne.0) then
      jul=bdate+real(itime,kind=dp)/86400._dp
      call caldate(jul,jjjjmmdd,ihmmss)      
    endif
    if ((mdomainfill.eq.0).and.(ipout.ge.1)) then
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
    if ((filesize.ge.max_partoutput_filesize).or.(mod(itime,loutrestart).eq.0)) then 
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
  integer :: j,ithread

  if (linit_cond.ge.1) then
    do j=1,numpart
      call initial_cond_calc(itime,j,1)
    end do
    do ithread=1,numthreads
      init_cond(:,:,:,:,:)=init_cond(:,:,:,:,:)+init_cond_omp(:,:,:,:,:,ithread)
    end do
  endif


  if (ipout.eq.2) call output_particles(itime)!,active_per_rel)     ! dump particle positions

  if (linit_cond.ge.1) then
    if(linversionout.eq.1) then
      call initial_cond_output_inversion(itime)   ! dump initial cond. field
    else
      call initial_cond_output(itime)   ! dump initial cond. fielf
    endif
  endif
end subroutine finalise_output

subroutine output_restart(itime,loutnext,outnum)
  use particle_mod
  use coordinates_ecmwf
  use netcdf_output_mod
  use unc_mod

  implicit none

  integer, intent(in) :: itime,loutnext
  real, intent(in) :: outnum
  integer :: i,j,jjjjmmdd,ihmmss,stat
  integer :: ks,kp,kz,nage,jy,ix,l
  real(kind=dp) :: jul
  character :: adate*8,atime*6


  jul=bdate+real(itime,kind=dp)/86400._dp
  call caldate(jul,jjjjmmdd,ihmmss)
  write(adate,'(i8.8)') jjjjmmdd
  write(atime,'(i6.6)') ihmmss

  restart_filename3 = restart_filename2
  restart_filename2 = restart_filename1
  restart_filename1 = path(2)(1:length(2))//'restart_'//adate//atime

  write(*,*) 'Writing Restart file:', trim(restart_filename1)

  open(unitrestart,file=restart_filename1,form='unformatted')

  ! Write current time to file
  !***************************

  write(unitrestart) itime
  write(unitrestart) count%allocated
  write(unitrestart) loutnext
  write(unitrestart) outnum

  do i=1,count%allocated
    if (part(i)%alive) then
      call update_zeta_to_z(itime,i)
      call update_z_to_zeta(itime,i)
    endif
    write(unitrestart) part(i)%xlon,part(i)%ylat,part(i)%z,part(i)%zeta, &
      part(i)%npoint,part(i)%nclass,part(i)%idt,part(i)%tend, &
      part(i)%tstart,part(i)%alive,part(i)%turbvel%u, &
      part(i)%turbvel%v,part(i)%turbvel%w,part(i)%mesovel%u, &
      part(i)%mesovel%v,part(i)%mesovel%w,(part(i)%mass(j),j=1,nspec), &
      (part(i)%wetdepo(j),j=1,nspec),(part(i)%drydepo(j),j=1,nspec)
  end do
  if (iout.gt.0) then
    write(unitrestart) tpointer
    do ks=1,nspec
      do kp=1,maxpointspec_act
        do nage=1,nageclass
          do jy=0,numygrid-1
            do ix=0,numxgrid-1
              do l=1,nclassunc
                do kz=1,numzgrid
                  write(unitrestart) gridunc(ix,jy,kz,ks,kp,l,nage)
                end do
                if ((wetdep).and.(ldirect.gt.0)) then
                  write(unitrestart) wetgridunc(ix,jy,ks,kp,l,nage)
                endif
                if ((drydep).and.(ldirect.gt.0)) then
                  write(unitrestart) drygridunc(ix,jy,ks,kp,l,nage)
                endif
              end do
            end do
          end do
          if (nested_output.eq.1) then
            do jy=0,numygridn-1
              do ix=0,numxgridn-1
                do l=1,nclassunc
                  do kz=1,numzgrid
                    write(unitrestart) griduncn(ix,jy,kz,ks,kp,l,nage)
                  end do
                  if ((wetdep).and.(ldirect.gt.0)) then
                    write(unitrestart) wetgriduncn(ix,jy,ks,kp,l,nage)
                  endif
                  if ((drydep).and.(ldirect.gt.0)) then
                    write(unitrestart) drygriduncn(ix,jy,ks,kp,l,nage)
                  endif
                end do
              end do
            end do
          endif
        end do
      end do
      if ((drybkdep).or.(wetbkdep)) then
        do i=1,count%allocated
          write(unitrestart) xscav_frac1(i,ks)
        end do
      endif
    end do
  endif
  close(unitrestart)

  open(unit=1234, iostat=stat, file=restart_filename3, status='old')
  if(stat == 0) close(1234, status='delete')
end subroutine output_restart

subroutine output_heightlevels(height_tmp,nmixz_tmp)
  implicit none

  real,intent(in) :: height_tmp(nzmax)
  integer,intent(in) :: nmixz_tmp
  integer :: kz
  character(len=256) :: heightlevels_filename

  heightlevels_filename = path(2)(1:length(2))//'heightlevels.bin'

  write(*,*) 'Writing Initialised heightlevels to file:', trim(heightlevels_filename)
  
  open(unitheightlevels,file=trim(heightlevels_filename),form='unformatted')

  write(unitheightlevels) nmixz_tmp

  do kz=1,nz
    write(unitheightlevels) height_tmp(kz)
  end do
  close(unitheightlevels)
end subroutine output_heightlevels

subroutine output_particles(itime,initial_output)
  !                        i
  !*****************************************************************************
  !                                                                            *
  !     Dump all particle positions                                            *
  !     Binary output is no longer supported. If required, function can be     *
  !     added below at "Put binary function here"                              *
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
  write(*,*) 'NETCDF missing! Please compile with netcdf if you want the particle dump.'
  stop
#endif

  if (present(initial_output)) then
    init_out=initial_output
  else
    init_out=.false.
  endif

!$OMP PARALLEL PRIVATE(i,j,m,tmp,ns,i_av,cartxyz_comp,cartxyz,np)
  ! Some variables needed for temporal interpolation
  !*************************************************
  call find_time_variables(itime)

!$OMP DO
  do i=1,numpart
    if ((.not. part(i)%alive) .or. &
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
    call determine_grid_coordinates(real(part(i)%xlon),real(part(i)%ylat))
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
          call horizontal_interpolation(oro,output(np,i))
          cycle
        case ('TR') ! Tropopause
          do m=1,2
            call horizontal_interpolation(tropopause,tmp(m),1,memind(m),1)
          end do
          call temporal_interpolation(tmp(1),tmp(2),output(np,i))
          cycle
        case ('HM') ! PBL height
          do m=1,2
            call horizontal_interpolation(hmix,tmp(m),1,memind(m),1)
          end do
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
            call interpol_partoutput_value(partopt(np)%name,output(np,i),i)
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
    write(*,*) part(1)%prob,part(1)%alive,count%alive,count%spawned,count%terminated
  endif

  ! Determine current calendar date, needed for the file name
  !**********************************************************

  jul=bdate+real(itime,kind=dp)/86400._dp
  call caldate(jul,jjjjmmdd,ihmmss)
  write(adate,'(i8.8)') jjjjmmdd
  write(atime,'(i6.6)') ihmmss
  j=1
  if (lnetcdfout.eq.1) then
#ifdef USE_NCF
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
!$OMP PARALLEL PRIVATE(np,ns)
!$OMP DO SCHEDULE(dynamic)
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
!$OMP END DO
!$OMP END PARALLEL
    endif
    call close_partoutput_file(ncid)
    if (.not. init_out) then
      mass_written=.true. ! needs to be reduced within openmp loop
      topo_written=.true. ! same
    endif
#endif
  else
    ! Put binary function here
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
  use interpol_mod, only: interpol_density
  use coordinates_ecmwf
  use particle_mod

  implicit none

  integer,intent(in) :: itime
  real,intent(in) :: weight
  integer :: itage,i,kz,ks,n,nage,thread,ithread
  integer :: il,ind,indz,indzp,nrelpointer
  integer :: ix,jy,ixp,jyp
  real :: ddx,ddy
  real(kind=dp) :: mm3
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

!$OMP PARALLEL PRIVATE(i,itage,nage,rhoi,nrelpointer,kz,xl,yl,ks,wx,wy,w,thread,ddx,ddy, &
!$OMP ix,jy,ixp,jyp) num_threads(numthreads_grid)
#if (defined _OPENMP)
    thread = OMP_GET_THREAD_NUM()+1 ! Starts with 1
#else
    thread = 1
#endif

!$OMP DO
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

 ! Some variables needed for temporal interpolation
  !*************************************************
  call find_time_variables(itime)

  xlon=xlon0+real(part(j)%xlon)*dx
  ylat=ylat0+real(part(j)%ylat)*dy

  !*****************************************************************************
  ! Interpolate several variables (PV, specific humidity, etc.) to particle position
  !*****************************************************************************

  call determine_grid_coordinates(real(part(j)%xlon),real(part(j)%ylat))
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
        call horizontal_interpolation(oro,output)
        part(j)%val_av(i_av)=part(j)%val_av(i_av)+output
      case ('tr')
        do m=1,2
          call horizontal_interpolation(tropopause,tr(m),1,memind(m),1)
        end do
        call temporal_interpolation(tr(1),tr(2),output)
        part(j)%val_av(i_av)=part(j)%val_av(i_av)+output
      case ('hm')
        do m=1,2
          call horizontal_interpolation(hmix,hm(m),1,memind(m),1)
        end do
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
        call interpol_partoutput_value(partopt(np)%name,output,j)
        part(j)%val_av(i_av)=part(j)%val_av(i_av)+output
    end select
  end do
  ! Reset dz1out
  !*************
  dz1out=-1
  cart_comp=.false.

  return
end subroutine partpos_average

end module output_mod
