! SPDX-FileCopyrightText: FLEXPART 1998-2019, see flexpart_license.txt
! SPDX-License-Identifier: GPL-3.0-or-later

! DJM - 2017-05-09 - added #ifdef USE_MPIINPLACE cpp directive to     *
! enable declaration of a gridunc0 array if required by MPI code in   *
! mpi_mod.f90                                                         *
!                                                                     *
!**********************************************************************

module unc_mod

  use par_mod, only:dp,dep_prec,nclassunc

  implicit none

  real,allocatable, dimension (:,:,:,:,:,:,:) :: gridunc
#ifdef USE_MPIINPLACE
#else
  ! If MPI_IN_PLACE option is not used in mpi_mod.f90::mpif_tm_reduce_grid(),
  ! then an aux array is needed for parallel grid reduction
  real(kind=dp),allocatable, dimension (:,:,:,:,:,:,:) :: gridunc0
  real(kind=dp),allocatable, dimension (:,:,:,:,:,:,:) :: griduncn0
#endif
  real(kind=dp),allocatable, dimension (:,:,:,:,:,:,:) :: griduncn
  real(dep_prec),allocatable, dimension (:,:,:,:,:,:) :: drygridunc
  real(dep_prec),allocatable, dimension (:,:,:,:,:,:) :: drygriduncn
  real(dep_prec),allocatable, dimension (:,:,:,:,:,:) :: wetgridunc
  real(dep_prec),allocatable, dimension (:,:,:,:,:,:) :: wetgriduncn
#ifdef _OPENMP
  real(kind=dp),allocatable, dimension (:,:,:,:,:,:,:,:) :: gridunc_omp
  real(kind=dp),allocatable, dimension (:,:,:,:,:,:,:,:) :: griduncn_omp
  real(dep_prec),allocatable, dimension (:,:,:,:,:,:,:) :: drygridunc_omp
  real(dep_prec),allocatable, dimension (:,:,:,:,:,:,:) :: drygriduncn_omp
  real(dep_prec),allocatable, dimension (:,:,:,:,:,:,:) :: wetgridunc_omp
  real(dep_prec),allocatable, dimension (:,:,:,:,:,:,:) :: wetgriduncn_omp
#endif
! For sum of individual contributions, used for the MPI version
  real(dep_prec),allocatable, dimension (:,:,:,:,:,:) :: drygridunc0
  real(dep_prec),allocatable, dimension (:,:,:,:,:,:) :: drygriduncn0
  real(dep_prec),allocatable, dimension (:,:,:,:,:,:) :: wetgridunc0
  real(dep_prec),allocatable, dimension (:,:,:,:,:,:) :: wetgriduncn0

contains

subroutine radioactive_decay()
  ! Accumulated deposited mass radioactively decays
  use com_mod

  implicit none

  integer ::                &
    j,i,                    & ! loop variable over grid
    ks,                     & ! loop variable species
    kp,                     & ! loop variable for maxpointspec_act
    l,                      & ! loop variable over nclassunc
    nage,                   & ! loop variable over age classes
    n                         ! loop variable over particles

!$OMP PARALLEL PRIVATE(ks,kp,nage,l,j,i)
!$OMP DO COLLAPSE(2)
  do ks=1,nspec
  do kp=1,maxpointspec_act
    if (decay(ks).gt.0.) then
      do nage=1,nageclass
        do l=1,nclassunc
  ! Mother output grid
          do j=0,numygrid-1
            do i=0,numxgrid-1
              wetgridunc(i,j,ks,kp,l,nage)= &
                   wetgridunc(i,j,ks,kp,l,nage)* &
                   exp(-1.*outstep*decay(ks))
              drygridunc(i,j,ks,kp,l,nage)= &
                   drygridunc(i,j,ks,kp,l,nage)* &
                   exp(-1.*outstep*decay(ks))
            end do
          end do
  ! Nested output grid
          if (nested_output.eq.1) then
            do j=0,numygridn-1
              do i=0,numxgridn-1
                wetgriduncn(i,j,ks,kp,l,nage)= &
                     wetgriduncn(i,j,ks,kp,l,nage)* &
                     exp(-1.*outstep*decay(ks))
                drygriduncn(i,j,ks,kp,l,nage)= &
                     drygriduncn(i,j,ks,kp,l,nage)* &
                     exp(-1.*outstep*decay(ks))
              end do
            end do
          endif
        end do
      end do
    endif
  end do
  end do
!$OMP END DO
!$OMP END PARALLEL
end subroutine radioactive_decay

end module unc_mod
