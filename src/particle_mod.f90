
  !*****************************************************************************
  !                                                                            *
  ! Module that organises particle information in derived types                *
  ! Particles are terminated and spawned using routines from this module so    *
  ! that global particle counts (spawned, allocated, terminated) are kept      *
  ! up to date                                                                 *
  !                                                                            *
  ! Author: L. Bakels 2022                                                     *
  !                                                                            *
  !*****************************************************************************

module particle_mod
  use com_mod, only: maxspec,DRYBKDEP,WETBKDEP,iout,n_average
  use par_mod, only: dp

  implicit none
  
  type :: coordinates
    real(kind=dp) ::              &
      xlon,                       & ! longitude in grid coordinates
      ylat                          ! latitude in grid coordinates
    real          ::              &
      z                             ! height in meters
    real          ::              &
      zeta                          ! height in eta (ECMWF) coordinates
  end type coordinates

  type :: velocities
    real               ::         &
      u,                          & ! x velocity
      v,                          & ! y velocity
      w                             ! z velocity
    real               ::         &
      weta                          ! z velocity in eta (ECMWF) coordinates
  end type velocities

  type :: particle
    real(kind=dp)      ::         &
      xlon,                       & ! Longitude in grid coordinates
      ylat,                       & ! Latitude in grid coordinates
      xlon_prev, ylat_prev,       & ! Keeping the previous positions in memory
      z,                          & ! height in meters
      z_prev                        ! Previous position
    real(kind=dp)      ::         &
      zeta,                       & ! Height in eta (ECMWF) coordinates
      zeta_prev                     ! Previous position
    type(velocities)   ::         &
      vel,                        & ! Velocities from interpolated windfields
      turbvel,                    & ! Random turbulent velocities
      mesovel                       ! Mesoscale turbulent velocities
    real               ::         &
      settling                      ! Settling velocity for dry and wet(?) deposit
    logical            ::         &
      alive=.false.,              & ! Flag to show if the particle is still in the running
      etaupdate=.false.,          & ! If false, z(meter) is more up-to-date than z(eta)
      meterupdate=.false.,         & ! If false, z(eta) is more up-to-date than z(meter)
      nstop=.false.                 ! Flag to stop particle (used in advance, stopped in timemanager)
    integer(kind=2)    ::         &
      icbt                          ! Forbidden state flag   
    integer            ::         &
      tstart,                     & ! spawning time in seconds after start
      tend,                       & ! termination time in seconds after start
      npoint,                     & ! release point
      nclass,                     &
      !species(maxspec),           & ! the number of the corresponding species file of the particle
      idt                           ! internal time of the particle
    real               ::         &
      mass(maxspec),              & ! Particle mass for each particle species
      mass_init(maxspec),         & ! Initial mass of each particle
      wetdepo(maxspec)=0.,        & ! Wet deposition (cumulative)
      drydepo(maxspec)=0.,        & ! Dry deposition (cumulative)
      prob(maxspec)                 ! Probability of absorption at ground due to dry deposition
    
    real,allocatable   ::         &
      val_av(:)                     ! Averaged values; only used when average_output=.true.
    real               ::         &
      ntime=0.,                   & ! Number of timesteps to average over
      cartx_av=0.,                & ! Averaged x pos;
      carty_av=0.,                & ! Averaged y pos;
      cartz_av=0.                   ! Averaged z pos;

  end type particle

  type :: particlecount          
    integer              ::       &
      alive=0,                    & ! Number of particles that are alive
      spawned=0,                  & ! Total number of spawned particles
      terminated=0,               & ! Total number of particles that have been terminated
      allocated=0,                & ! Number of total allocated particle spaces
      ninmem=0                      ! Number of particles currently in memory
    logical,allocatable  ::       &
      inmem(:)
  end type

  type(particle), allocatable ::  &
    part(:)                         ! This is where all particles are being stored
  type(particlecount)         ::  &
    count                           ! Keeping track of global particle number within the simulation
  real,allocatable            ::  &
    xscav_frac1(:,:)                ! Only allocated when wet or dry deposit backward mode is switched on
  real,allocatable            ::  &
    xplum(:),yplum(:),zplum(:)      ! Only allocated for iout=4 or 5 (plumetraj)
  integer,allocatable         ::  &
    nclust(:)                       ! Only allocated for iout=4 or 5 (plumetraj)
  ! private ::                      &
  !   count             
  public ::                       &
    particle,                     &
    part,                         &
    alloc_particles,              &
    dealloc_particle_range,       &
    dealloc_particle,             &
    dealloc_all_particles,        &
    terminate_particle,           &
    spawn_particle,               &
    spawn_particles,              &
    get_totalpart_num,            &
    get_alivepart_num,            &
    get_newpart_index,            &
    particle_allocated,           &
    update_xlon,                  &
    update_ylat,                  &
    update_z,                     &
    count

  interface update_xlon
    procedure update_xlon_dp, update_xlon_sp, update_xlon_int
  end interface update_xlon

  interface set_xlon
    procedure set_xlon_dp, set_xlon_sp, set_xlon_int
  end interface set_xlon

  interface update_ylat
    procedure update_ylat_dp, update_ylat_sp, update_ylat_int
  end interface update_ylat

  interface set_ylat
    procedure set_ylat_dp, set_ylat_sp, set_ylat_int
  end interface set_ylat

  interface update_z
    procedure update_z_dp,update_z_sp
  end interface update_z

  interface set_z
    procedure set_z_dp,set_z_sp
  end interface set_z

  interface update_zeta
    procedure update_zeta_dp,update_zeta_sp
  end interface update_zeta

  interface set_zeta
    procedure set_zeta_dp,set_zeta_sp
  end interface set_zeta
   
contains

  logical function particle_allocated(ipart)
    !******************************************
    ! Checks if the memory of the particle is *
    ! still allocated                         *
    !******************************************

    implicit none 

    integer, intent(in)    :: ipart   ! Particle index
    !logical :: is_particle_allocated
    
    if (ipart.gt.count%allocated) then
      particle_allocated = .false.
    else
      particle_allocated = count%inmem(ipart)
    endif
  end function particle_allocated

  subroutine get_newpart_index(ipart)
    !**************************************************
    ! Returns the first free spot to put a new particle
    !**************************************************
    implicit none

    integer, intent(inout) :: ipart   ! First free index

    ipart = count%spawned + 1
  end subroutine get_newpart_index

  subroutine get_totalpart_num(npart)
    !********************************************
    ! Returns total number of particles spawned *
    !********************************************
    implicit none 

    integer, intent(inout) :: npart   ! Number of particles

    npart = count%spawned
  end subroutine get_totalpart_num

  subroutine get_alivepart_num(npart)
    !**********************************************
    ! Returns number of particles currently alive *
    !**********************************************
    implicit none 

    integer, intent(inout) :: npart   ! Number of particles

    npart = count%alive
  end subroutine get_alivepart_num

  subroutine spawn_particles(itime, nmpart)
    !******************************************************
    ! Spawning particles
    !
    ! This routine spawns new particles and allocates the 
    ! memory if necessary.
    !******************************************************
    implicit none 

    integer, intent(in) :: &
      itime,               &  ! spawning time
      nmpart                  ! number of particles that are being spawned

    ! Check if new memory needs to be allocated 
    !*******************************************
    if (nmpart+count%spawned.gt.count%allocated) then
      call alloc_particles( (nmpart+count%spawned) - count%allocated )
    endif
    ! Update the number of particles that are currently alive
    !********************************************************
    count%alive = count%alive + nmpart

    ! Set the spawning time for each new particle and mark it as alive
    !*****************************************************************
    part(count%spawned+1:count%spawned+nmpart)%tstart = itime
    part(count%spawned+1:count%spawned+nmpart)%alive = .true.

    ! Update the total number of spawned particles
    !*********************************************
    count%spawned = count%spawned + nmpart
  end subroutine spawn_particles

  subroutine spawn_particle(itime, ipart)
    !******************************************************
    ! Spawning particles
    !
    ! This routine spawns new particles and allocates the 
    ! memory if necessary.
    !******************************************************
    implicit none 

    integer, intent(in) :: &
      itime,               & ! spawning time
      ipart                  ! number of particles that are being spawned

    ! Check if new memory needs to be allocated 
    !*******************************************
    if (.not. particle_allocated(ipart)) call alloc_particle(ipart)

    if (part(ipart)%alive) error stop 'Attempting to overwrite existing particle'

    ! Update the number of particles that are currently alive
    !********************************************************
    count%alive = count%alive + 1

    ! Set the spawning time for each new particle and mark it as alive
    !*****************************************************************
    part(ipart)%tstart = itime
    part(ipart)%alive = .true.

    ! Update the total number of spawned particles
    !*********************************************
    count%spawned = count%spawned + 1
  end subroutine spawn_particle

  subroutine terminate_particle(ipart,itime)
    !*****************************************************
    ! Terminating specified particle
    !
    ! This routine terminates a selected particle
    !***************************************************** 
    implicit none

    integer, intent(in) :: &
      ipart,               & ! to be terminated particle index
      itime                  ! Time at which particle is terminated

    ! Flagging the particle as having been terminated
    !************************************************
    part(ipart)%alive=.false.
    part(ipart)%tend=itime

    ! Update the number of current particles that are alive
    !******************************************************
    count%alive = count%alive - 1

    ! Update the total number of terminated particles during the whole run
    !**********************************************************************
    count%terminated = count%terminated + 1
  end subroutine terminate_particle
 
  subroutine alloc_particles(nmpart)

    implicit none 

    integer, intent(in)        :: nmpart
    type(particle),allocatable :: tmppart(:)
    logical, allocatable       :: tmpcount(:)
    real, allocatable          :: tmpxscav(:,:)
    real, allocatable          :: tmpxl(:),tmpyl(:),tmpzl(:)
    integer, allocatable       :: tmpnclust(:)
    integer                    :: i

    if (nmpart.gt.100) write(*,*) 'Allocating ',nmpart,' particles'

    ! Keeping track of the allocated memory in case 
    ! there is a reason for deallocating some of it
    !**********************************************
    allocate( tmpcount(count%allocated+nmpart) )
    if (count%allocated.gt.0) tmpcount(1:count%allocated) = count%inmem
    call move_alloc(tmpcount,count%inmem)
    count%inmem(count%allocated+1:count%allocated+nmpart) = .true.

    ! Allocating new particle spaces
    !*******************************

    allocate( tmppart(count%allocated+nmpart) )
    if (n_average.gt.0) then 
      do i=1,count%allocated+nmpart
        allocate( tmppart(i)%val_av(n_average) )
        tmppart(i)%val_av = 0
      end do
    endif
    if (count%allocated.gt.0) tmppart(1:count%allocated) = part
    call move_alloc(tmppart,part)

    ! If wet or dry deposition backward mode is switched on, xscav_frac1
    ! needs to be allocated
    !*******************************************************************
    if (WETBKDEP.or.DRYBKDEP) then
      allocate( tmpxscav(count%allocated+nmpart,maxspec) )
      if (count%allocated.gt.0) tmpxscav(1:count%allocated,:) = xscav_frac1
      call move_alloc(tmpxscav,xscav_frac1)
    endif

    if ((iout.eq.4).or.(iout.eq.5)) then
      allocate( tmpxl(count%allocated+nmpart) )
      if (count%allocated.gt.0) tmpxl(1:count%allocated) = xplum
      call move_alloc(tmpxl,xplum)
      
      allocate( tmpyl(count%allocated+nmpart) )
      if (count%allocated.gt.0) tmpyl(1:count%allocated) = yplum
      call move_alloc(tmpyl,yplum)

      allocate( tmpzl(count%allocated+nmpart) )
      if (count%allocated.gt.0) tmpzl(1:count%allocated) = zplum
      call move_alloc(tmpzl,zplum)

      allocate( tmpnclust(count%allocated+nmpart) )
      if (count%allocated.gt.0) tmpnclust(1:count%allocated) = nclust
      call move_alloc(tmpnclust,nclust)
    endif

    count%allocated = count%allocated+nmpart
    if (nmpart.gt.100) write(*,*) 'Finished allocation'
  end subroutine alloc_particles

  subroutine alloc_particle(ipart)

    implicit none 

    integer, intent(in) :: ipart

    ! Keeping track of the allocated memory in case 
    ! there is a reason for deallocating some of it
    !**********************************************
    if (ipart.gt.count%allocated) then 
      call alloc_particles(ipart-count%allocated)
    else
      error stop 'Error: You are trying to allocate an already existing particle'
    endif

  end subroutine alloc_particle

  subroutine dealloc_particle_range(istart,iend)

    implicit none

    integer, intent(in) :: istart,iend

    !deallocate( part(istart:iend) )
    count%inmem(istart:iend) = .false.
  end subroutine dealloc_particle_range

  subroutine dealloc_particle(ipart)

    implicit none

    integer, intent(in) :: ipart ! particle index

    !deallocate( part(ipart) )
    part = part(1:ipart) ! FORTRAN 2008 only
    count%inmem(ipart+1:) = .false.
  end subroutine dealloc_particle

  subroutine dealloc_all_particles()

    implicit none

    integer :: i

    if (n_average.gt.0) then 
      do i=1,count%allocated
        deallocate( part(i)%val_av )
      end do
    endif
    deallocate( part )
    deallocate( count%inmem )

    if (WETBKDEP.or.DRYBKDEP) then
      deallocate( xscav_frac1 )
    endif

    if ((iout.eq.4).or.(iout.eq.5)) then
      deallocate( xplum )
      deallocate( yplum )
      deallocate( zplum )
      deallocate( nclust )
    endif
  end subroutine dealloc_all_particles

! Update_xlon
  subroutine update_xlon_dp(ipart,xchange)
    !**************************************
    ! Updates the longitude of the particle
    !**************************************
    implicit none

    integer, intent(in)       :: ipart   ! particle index
    real(kind=dp), intent(in) :: xchange

    part(ipart)%xlon = part(ipart)%xlon + xchange
  end subroutine update_xlon_dp

  subroutine update_xlon_sp(ipart,xchange)
    !**************************************
    ! Updates the longitude of the particle
    !**************************************
    implicit none

    integer, intent(in) :: ipart    ! particle index
    real, intent(in)    :: xchange

    part(ipart)%xlon = part(ipart)%xlon + real(xchange,kind=dp)
  end subroutine update_xlon_sp

  subroutine update_xlon_int(ipart,xchange)
    !**************************************
    ! Updates the longitude of the particle
    !**************************************
    implicit none

    integer, intent(in) :: ipart  ! particle index
    integer, intent(in) :: xchange

    part(ipart)%xlon = part(ipart)%xlon + real(xchange,kind=dp)
  end subroutine update_xlon_int
! End Update_xlon

! Set_xlon
  subroutine set_xlon_dp(ipart,xvalue)
    !**************************************
    ! Sets the longitude of the particle
    !**************************************
    implicit none

    integer, intent(in)       :: ipart  ! particle index
    real(kind=dp), intent(in) :: xvalue

    part(ipart)%xlon = xvalue
  end subroutine set_xlon_dp

  subroutine set_xlon_sp(ipart,xvalue)
    !**************************************
    ! Sets the longitude of the particle
    !**************************************
    implicit none

    integer, intent(in)    :: ipart   ! particle index
    real, intent(in)       :: xvalue

    part(ipart)%xlon = real(xvalue,kind=dp)
  end subroutine set_xlon_sp

  subroutine set_xlon_int(ipart,xvalue)
    !**************************************
    ! Sets the longitude of the particle
    !**************************************
    implicit none

    integer, intent(in)    :: ipart  ! particle index
    integer, intent(in)    :: xvalue

    part(ipart)%xlon = real(xvalue,kind=dp)
  end subroutine set_xlon_int
! End Set_xlon 

! Update_ylat
  subroutine update_ylat_dp(ipart,ychange)
    !**************************************
    ! Updates the latitude of the particle
    !**************************************
    implicit none

    integer, intent(in)       :: ipart  ! particle index
    real(kind=dp), intent(in) :: ychange

    part(ipart)%ylat = part(ipart)%ylat + ychange
  end subroutine update_ylat_dp

  subroutine update_ylat_sp(ipart,ychange)
    !**************************************
    ! Updates the latitude of the particle
    !**************************************
    implicit none

    integer, intent(in)    :: ipart  ! particle index
    real, intent(in)       :: ychange

    part(ipart)%ylat = part(ipart)%ylat + real(ychange,kind=dp)
  end subroutine update_ylat_sp

  subroutine update_ylat_int(ipart,ychange)
    !**************************************
    ! Updates the latitude of the particle
    !**************************************
    implicit none

    integer, intent(in)    :: ipart ! particle index
    integer, intent(in)    :: ychange

    part(ipart)%ylat = part(ipart)%ylat + real(ychange,kind=dp)
  end subroutine update_ylat_int
! End Update_ylat

! Set_ylat
  subroutine set_ylat_dp(ipart,yvalue)
    !**************************************
    ! Sets the latitude of the particle
    !**************************************
    implicit none

    integer, intent(in)       :: ipart  ! particle index
    real(kind=dp), intent(in) :: yvalue

    part(ipart)%ylat = yvalue
  end subroutine set_ylat_dp

  subroutine set_ylat_sp(ipart,yvalue)
    !**************************************
    ! Sets the latitude of the particle
    !**************************************
    implicit none

    integer, intent(in)    :: ipart  ! particle index
    real, intent(in)       :: yvalue

    part(ipart)%ylat = real(yvalue,kind=dp)
  end subroutine set_ylat_sp

  subroutine set_ylat_int(ipart,yvalue)
    !**************************************
    ! Sets the latitude of the particle
    !**************************************
    implicit none

    integer, intent(in)    :: ipart  ! particle index
    integer, intent(in) :: yvalue

    part(ipart)%ylat = real(yvalue,kind=dp)
  end subroutine set_ylat_int
! End Set_ylat

! Update z positions
  subroutine update_z_dp(ipart,zchange)
    !**************************************
    ! Updates the height of the particle
    !**************************************
    implicit none

    integer, intent(in)        :: ipart  ! particle index
    real(kind=dp), intent(in)  :: zchange

    part(ipart)%z = part(ipart)%z + zchange
    part(ipart)%meterupdate=.false.
    part(ipart)%etaupdate=.true.
  end subroutine update_z_dp

  subroutine update_z_sp(ipart,zchange)
    !**************************************
    ! Updates the height of the particle
    !**************************************
    implicit none

    integer, intent(in)    :: ipart  ! particle index
    real, intent(in)       :: zchange

    part(ipart)%z = part(ipart)%z + real(zchange,kind=dp)
    part(ipart)%meterupdate=.false.
    part(ipart)%etaupdate=.true.
  end subroutine update_z_sp  

  subroutine update_zeta_dp(ipart,zchange)
    !**************************************
    ! Updates the height of the particle
    !**************************************
    implicit none

    integer, intent(in)        :: ipart  ! particle index
    real(kind=dp), intent(in)  :: zchange

    part(ipart)%zeta = part(ipart)%zeta + zchange
    part(ipart)%etaupdate=.false.
    part(ipart)%meterupdate=.true.
  end subroutine update_zeta_dp

  subroutine update_zeta_sp(ipart,zchange)
    !**************************************
    ! Updates the height of the particle
    !**************************************
    implicit none

    integer, intent(in)    :: ipart  ! particle index
    real, intent(in)       :: zchange

    part(ipart)%zeta = part(ipart)%zeta + real(zchange,kind=dp)
    part(ipart)%etaupdate=.false.
    part(ipart)%meterupdate=.true.
  end subroutine update_zeta_sp
! End update z positions

! Update z positions
  subroutine set_z_dp(ipart,zvalue)
    !**************************************
    ! Updates the height of the particle
    !**************************************
    implicit none

    integer, intent(in)        :: ipart  ! particle index
    real(kind=dp), intent(in)  :: zvalue

    part(ipart)%z = zvalue
    part(ipart)%meterupdate=.false.
    part(ipart)%etaupdate=.true.
  end subroutine set_z_dp  

  subroutine set_z_sp(ipart,zvalue)
    !**************************************
    ! Updates the height of the particle
    !**************************************
    implicit none

    integer, intent(in)    :: ipart  ! particle index
    real, intent(in)       :: zvalue

    part(ipart)%z = real(zvalue,kind=dp)
    part(ipart)%meterupdate=.false.
    part(ipart)%etaupdate=.true.
  end subroutine set_z_sp

  subroutine set_zeta_dp(ipart,zvalue)
    !**************************************
    ! Updates the height of the particle
    !**************************************
    implicit none

    integer, intent(in)        :: ipart  ! particle index
    real(kind=dp), intent(in)  :: zvalue

    part(ipart)%zeta = zvalue
    part(ipart)%etaupdate=.false.
    part(ipart)%meterupdate=.true.
  end subroutine set_zeta_dp

  subroutine set_zeta_sp(ipart,zvalue)
    !**************************************
    ! Updates the height of the particle
    !**************************************
    implicit none

    integer, intent(in)    :: ipart  ! particle index
    real, intent(in)       :: zvalue

    part(ipart)%zeta = real(zvalue,kind=dp)
    part(ipart)%etaupdate=.false.
    part(ipart)%meterupdate=.true.
  end subroutine set_zeta_sp
! End update z positions

end module particle_mod
