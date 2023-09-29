module restart_mod
  use particle_mod
#ifdef ETA
  use coord_ecmwf_mod
#endif
  use outgrid_mod
  use unc_mod
  use date_mod
#ifdef USE_NCF
  use netcdf_output_mod
#endif

  character(len=256) :: restart_filename1,restart_filename2,restart_filename3

contains

subroutine output_restart(itime,loutnext,outnum)

  implicit none

  integer, intent(in) :: itime,loutnext
  real, intent(in) :: outnum
  integer :: i,j,jjjjmmdd,ihmmss,stat
  integer :: ks,kp,kz,nage,jy,ix,l,n
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
  write(unitrestart) numreceptor


  do i=1,count%allocated
#ifdef ETA
    if (part(i)%alive) then
      call update_zeta_to_z(itime,i)
      call update_z_to_zeta(itime,i)
    endif
#endif
    write(unitrestart) part(i)%xlon,part(i)%ylat,part(i)%z, &
#ifdef ETA
      part(i)%zeta, &
#endif
      part(i)%npoint,part(i)%nclass,part(i)%idt,part(i)%tend, &
      part(i)%tstart,part(i)%alive,part(i)%turbvel%u, &
      part(i)%turbvel%v,part(i)%turbvel%w,part(i)%mesovel%u, &
      part(i)%mesovel%v,part(i)%mesovel%w,(part(i)%mass(j),j=1,nspec), &
      (part(i)%mass_init(j),j=1,nspec)
    if (wetdep) write(unitrestart) (part(i)%wetdepo(j),j=1,nspec)
    if (drydep) write(unitrestart) (part(i)%drydepo(j),j=1,nspec)
  end do
  if (iout.gt.0) then
#ifdef USE_NCF
    if (lnetcdfout.eq.1) write(unitrestart) tpointer
#endif
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
          if (iflux.eq.1) then
            do kz=1,numzgrid
              do jy=0,numygridn-1
                do ix=0,numxgridn-1
                  do i=1,5
                    write(unitrestart) flux(i,ix,jy,kz,ks,kp,nage)
                  end do
                end do
              end do
            end do        
          endif
        end do
        if (linit_cond.gt.0) then
          do kz=1,numzgrid
            do jy=0,numygridn-1
              do ix=0,numxgridn-1
                write(unitrestart) init_cond(ix,jy,kz,ks,kp)
              end do 
            end do
          end do 
        endif
      end do
      if ((drybkdep).or.(wetbkdep)) then
        do i=1,count%allocated
          write(unitrestart) xscav_frac1(i,ks)
        end do
      endif
      if (numreceptor.gt.0) then 
        do n=1,numreceptor
          write(unitrestart) creceptor(n,ks)
        end do
      endif
    end do
  endif
  close(unitrestart)

  ! open(unit=1234, iostat=stat, file=restart_filename3, status='old')
  ! if(stat == 0) close(1234, status='delete')
end subroutine output_restart

subroutine readrestart

  !*****************************************************************************
  !                                                                            *
  !   This routine opens the particle dump file and reads all the particle     *
  !   positions and gridded information from a previous run to initialize      *
  !   the current run.                                                         *
  !                                                                            *
  !     Author: L. Bakels 2022                                                 *
  !                                                                            *
  !*****************************************************************************

  implicit none

  integer :: i,j,ios
  integer :: id1,id2,it1,it2
  integer :: ks,kp,kz,nage,jy,ix,l,n
  real(kind=dp) :: julin

  numparticlecount=0


  open(unitpartin,file=path(2)(1:length(2))//'restart.bin', &
       form='unformatted',err=9989)

  write(*,*) 'Reading Restart file:', path(2)(1:length(2))//'restart.bin'
  
  read(unitpartin,iostat=ios) itime_init
  read(unitpartin) numpart ! count%allocated
  read(unitpartin) loutnext_init
  read(unitpartin) outnum_init
  read(unitpartin) numreceptor

  call spawn_particles(itime_init, numpart)
  do i=1,numpart
    read(unitpartin) part(i)%xlon,part(i)%ylat,part(i)%z, &
#ifdef ETA
      part(i)%zeta, &
#endif
      part(i)%npoint,part(i)%nclass,part(i)%idt,part(i)%tend, &
      part(i)%tstart,part(i)%alive,part(i)%turbvel%u, &
      part(i)%turbvel%v,part(i)%turbvel%w,part(i)%mesovel%u, &
      part(i)%mesovel%v,part(i)%mesovel%w,(part(i)%mass(j),j=1,nspec), &
      (part(i)%mass_init(j),j=1,nspec)
    if (wetdep) read(unitrestart) (part(i)%wetdepo(j),j=1,nspec)
    if (drydep) read(unitrestart) (part(i)%drydepo(j),j=1,nspec)
#ifdef ETA
    part(i)%etaupdate=.true.
    part(i)%meterupdate=.true.
#endif
    if (.not. part(i)%alive) then
      if (part(i)%tstart.le.itime_init) then
        call terminate_particle(i,part(i)%tend)
      else ! Particle is not spawned yet (original run with ipin=3)
        count%alive = count%alive - 1
        count%spawned = count%spawned -1
      endif
    endif
  end do
  if (iout.gt.0) then 
#ifdef USE_NCF
    if (lnetcdfout.eq.1) read(unitpartin) tpointer
#endif
    do ks=1,nspec
      do kp=1,maxpointspec_act
        do nage=1,nageclass
          do jy=0,numygrid-1
            do ix=0,numxgrid-1
              do l=1,nclassunc
                do kz=1,numzgrid
                  read(unitpartin) gridunc(ix,jy,kz,ks,kp,l,nage)
                end do
                if ((wetdep).and.(ldirect.gt.0)) then
                  read(unitpartin) wetgridunc(ix,jy,ks,kp,l,nage)
                endif
                if ((drydep).and.(ldirect.gt.0)) then
                  read(unitpartin) drygridunc(ix,jy,ks,kp,l,nage)
                endif
              end do
            end do
          end do
          if (nested_output.eq.1) then
            do jy=0,numygridn-1
              do ix=0,numxgridn-1
                do l=1,nclassunc
                  do kz=1,numzgrid
                    read(unitpartin) griduncn(ix,jy,kz,ks,kp,l,nage)
                  end do
                  if ((wetdep).and.(ldirect.gt.0)) then
                    read(unitpartin) wetgriduncn(ix,jy,ks,kp,l,nage)
                  endif
                  if ((drydep).and.(ldirect.gt.0)) then
                    read(unitpartin) drygriduncn(ix,jy,ks,kp,l,nage)
                  endif
                end do
              end do
            end do
          endif
          if (iflux.eq.1) then
            do kz=1,numzgrid
              do jy=0,numygridn-1
                do ix=0,numxgridn-1
                  do i=1,5
                    read(unitpartin) flux(i,ix,jy,kz,ks,kp,nage)
                  end do
                end do
              end do
            end do        
          endif
        end do
        if (linit_cond.gt.0) then
          do kz=1,numzgrid
            do jy=0,numygridn-1
              do ix=0,numxgridn-1
                read(unitpartin) init_cond(ix,jy,kz,ks,kp)
              end do 
            end do
          end do 
        endif
      end do
      if ((drybkdep).or.(wetbkdep)) then
        do i=1,numpart
          read(unitpartin) xscav_frac1(i,ks)
        end do
      endif
      if (numreceptor.gt.0) then 
        do n=1,numreceptor
          read(unitpartin) creceptor(n,ks)
        end do
      endif
    end do
  endif
  close(unitpartin)

  numpart=count%spawned
  
  julin=juldate(ibdate,ibtime)+real(itime_init,kind=dp)/86400._dp
  if (abs(julin-bdate).le.1.e-5) then
    write(*,*) ' #### FLEXPART ERROR: PLEASE KEEP IBDATE     #### '
    write(*,*) ' #### AND IBTIME INTACT FROM THE INITIAL RUN!#### '
    error stop
  endif
  call caldate(julin,id1,it1)
  call caldate(bdate,id2,it2)
  write(*,*) ' #### Restarting Flexpart from restart.bin.    #### '
  write(*,*) ' #### Original run started on                  #### '
  write(*,*) 'bdate: ',bdate,id2,it2
  write(*,*) ' #### Restarting run starts on                 #### '
  write(*,*) 'julin: ',julin,id1,it1

  return

9989   write(*,*) ' #### FLEXPART MODEL ERROR!   THE FILE             #### '
  write(*,*) ' #### '//path(2)(1:length(2))//'restart.bin'//'    #### '
  write(*,*) ' #### CANNOT BE OPENED. IF A FILE WITH THIS        #### '
  write(*,*) ' #### NAME DOES NOT EXISTS, RENAME THE APPROPRIATE #### '
  write(*,*) ' #### RESTART FILE TO restart.bin.                 #### '
end subroutine readrestart

end module restart_mod