module particle_mod
  use com_mod, only: maxspec
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
      xlon_prev, ylat_prev          ! Keeping the previous positions in memory
    real               ::         &
      z,                          & ! height in meters
      z_prev                        ! Previous position
    real               ::         &
      zeta                          ! Height in eta (ECMWF) coordinates
    type(velocities)   ::         &
      vel,                        & ! Velocities from interpolated windfields
      turbvel,                    & ! Random turbulent velocities
      mesovel                       ! Mesoscale turbulent velocities
    logical            ::         &
      alive=.false.,              & ! Flag to show if the particle is still in the running
      etaupdate=.false.,          & ! Set when the z(meter) is up-to-date with z(eta)
      nstop=.false.                 ! Flag to stop particle (used in advance, stopped in timemanager)
    integer(kind=2)    ::         &
      icbt                          ! Forbidden state flag   
    integer            ::         &
      tstart,                     & ! spawning time in seconds after start
      tend,                       & ! termination time in seconds after start
      npoint,                     & ! release point
      nclass,                     &
      idt
    real               ::         &
      mass(maxspec),              & ! Particle mass for each particle species
      prob(maxspec)                 ! Probability of absorption at ground due to dry deposition
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

  ! private ::                      &
  !   count             
  public ::                       &
    particle,                     &
    part,                         &
    allocate_particles,           &
    deallocate_particle_range,    &
    deallocate_particle,          &
    deallocate_all_particles,     &
    terminate_particle,           &
    spawn_particles,              &
    spawn_particle,               &
    get_total_part_num,           &
    get_alive_part_num,           &
    get_new_part_index,           &
    is_particle_allocated
    
contains

  logical function is_particle_allocated(ipart)
    !******************************************
    ! Checks if the memory of the particle is *
    ! still allocated                         *
    !******************************************

    implicit none 

    integer, intent(in)    :: &
      ipart                     ! Particle index
    !logical :: is_particle_allocated
    
    if (ipart.gt.count%allocated) then
      is_particle_allocated = .false.
    else
      is_particle_allocated = count%inmem(ipart)
    endif
  end function is_particle_allocated

  ! function is_particle_allocated(ipart) result(answer)
  !   !******************************************
  !   ! Checks if the memory of the particle is *
  !   ! still allocated                         *
  !   !******************************************

  !   implicit none 

  !   integer, intent(in)    :: &
  !     ipart                     ! Particle index
  !   logical                :: &
  !     answer 

  !   answer(ipart) = count%inmem(ipart)
  ! end function is_particle_allocated
  subroutine get_new_part_index(ipart)
    !**************************************************
    ! Returns the first free spot to put a new particle
    !**************************************************
    implicit none

    integer, intent(inout) :: &
      ipart                     ! First free index

    ipart = count%spawned + 1
  end subroutine get_new_part_index

  subroutine get_total_part_num(npart)
    !********************************************
    ! Returns total number of particles spawned *
    !********************************************
    implicit none 

    integer, intent(inout) :: &
      npart                     ! Number of particles

    npart = count%spawned
  end subroutine get_total_part_num

  subroutine get_alive_part_num(npart)
    !**********************************************
    ! Returns number of particles currently alive *
    !**********************************************
    implicit none 

    integer, intent(inout) :: &
      npart                     ! Number of particles

    npart = count%alive
  end subroutine get_alive_part_num

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
      call allocate_particles( count%allocated-(nmpart+count%spawned) )
    endif
    ! Update the number of particles that are currently alive
    !********************************************************
    count%alive = count%alive + nmpart

    ! Set the spawning time for each new particle and mark it as alive
    !*****************************************************************
    part(count%spawned:count%spawned+nmpart)%tstart = itime
    part(count%spawned:count%spawned+nmpart)%alive = .true.

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
    if (.not. is_particle_allocated(ipart)) then
      call allocate_particle(ipart)
    endif

    if (part(ipart)%alive) stop 'Attempting to overwrite existing particle'

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

  
  subroutine terminate_particle(ipart)
    !*****************************************************
    ! Terminating specified particle
    !
    ! This routine terminates a selected particle
    !***************************************************** 
    implicit none

    integer, intent(in) :: ipart ! to be terminated particle index

    ! Flagging the particle as having been terminated
    !************************************************
    part(ipart)%alive=.false.

    ! Update the number of current particles that are alive
    !******************************************************
    count%alive = count%alive - 1

    ! Update the total number of terminated particles during the whole run
    !**********************************************************************
    count%terminated = count%terminated + 1
  end subroutine terminate_particle
 
  subroutine allocate_particles(nmpart)

    implicit none 

    integer, intent(in) :: nmpart
    type(particle),allocatable :: tmppart(:)
    logical, allocatable :: tmpcount(:)

    if (nmpart.gt.100) write(*,*) 'Allocating ',nmpart,' particles'

    ! Keeping track of the allocated memory in case 
    ! there is a reason for deallocating some of it
    !**********************************************
    allocate( tmpcount(count%allocated+nmpart) )
    tmpcount(1:count%allocated) = count%inmem
    call move_alloc(tmpcount,count%inmem)
    count%inmem(count%allocated+1:count%allocated+nmpart) = .true.

    ! Allocating new particle spaces
    !*******************************

    allocate( tmppart(count%allocated+nmpart) )
    tmppart(1:count%allocated) = part
    call move_alloc(tmppart,part)
    count%allocated = count%allocated+nmpart
    if (nmpart.gt.100) write(*,*) 'Finished allocation'
  end subroutine allocate_particles

  subroutine allocate_particle(ipart)

    implicit none 

    integer, intent(in) :: ipart

    ! Keeping track of the allocated memory in case 
    ! there is a reason for deallocating some of it
    !**********************************************
    if (ipart.gt.count%allocated) then 
      call allocate_particles(1)
    else
      stop 'Error: You are trying to allocate an already existing particle'
    endif

  end subroutine allocate_particle

  subroutine deallocate_particle_range(istart,iend)

    implicit none

    integer, intent(in) :: istart,iend

    !deallocate( part(istart:iend) )
    count%inmem(istart:iend) = .false.
  end subroutine deallocate_particle_range

  subroutine deallocate_particle(ipart)

    implicit none

    integer, intent(in) :: ipart ! particle index

    !deallocate( part(ipart) )
    part = part(1:ipart)
    count%inmem(ipart:) = .false.
  end subroutine deallocate_particle

  subroutine deallocate_all_particles()

    implicit none

    deallocate( part )
    deallocate( count%inmem )
  end subroutine deallocate_all_particles

end module particle_mod