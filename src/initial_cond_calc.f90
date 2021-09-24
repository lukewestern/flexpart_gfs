! SPDX-FileCopyrightText: FLEXPART 1998-2019, see flexpart_license.txt
! SPDX-License-Identifier: GPL-3.0-or-later

subroutine initial_cond_calc(itime,i)
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

  use unc_mod
  use outg_mod
  use par_mod
  use com_mod
  use interpol_mod, only: interpol_density,ix,jy,ixp,jyp
  use coordinates_ecmwf
  use particle_mod

  implicit none

  integer :: itime,i,kz,ks
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
    call update_zcoord(itime,i)
    call interpol_density(i,rhoi)
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
    if (outheight(kz).gt.part(i)%z) exit
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
          init_cond(ix,jy,kz,ks,nrelpointer)= &
               init_cond(ix,jy,kz,ks,nrelpointer)+ &
               part(i)%mass(ks)/rhoi
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
            init_cond(ix,jy,kz,ks,nrelpointer)= &
                 init_cond(ix,jy,kz,ks,nrelpointer)+part(i)%mass(ks)/rhoi*w
          end do
        endif

        if ((jyp.ge.0).and.(jyp.le.numygrid-1)) then
          w=wx*(1.-wy)
          do ks=1,nspec
            init_cond(ix,jyp,kz,ks,nrelpointer)= &
                 init_cond(ix,jyp,kz,ks,nrelpointer)+part(i)%mass(ks)/rhoi*w
          end do
        endif
      endif


      if ((ixp.ge.0).and.(ixp.le.numxgrid-1)) then
        if ((jyp.ge.0).and.(jyp.le.numygrid-1)) then
          w=(1.-wx)*(1.-wy)
          do ks=1,nspec
            init_cond(ixp,jyp,kz,ks,nrelpointer)= &
                 init_cond(ixp,jyp,kz,ks,nrelpointer)+part(i)%mass(ks)/rhoi*w
          end do
        endif

        if ((jy.ge.0).and.(jy.le.numygrid-1)) then
          w=(1.-wx)*wy
          do ks=1,nspec
            init_cond(ixp,jy,kz,ks,nrelpointer)= &
                 init_cond(ixp,jy,kz,ks,nrelpointer)+part(i)%mass(ks)/rhoi*w
          end do
        endif
      endif
    endif

  endif

end subroutine initial_cond_calc
