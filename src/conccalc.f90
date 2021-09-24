! SPDX-FileCopyrightText: FLEXPART 1998-2019, see flexpart_license.txt
! SPDX-License-Identifier: GPL-3.0-or-later

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

  integer :: itime,itage,i,kz,ks,n,nage
  integer :: il,ind,indz,indzp,nrelpointer
  real :: weight,hx,hy,hz,h,xd,yd,zd,xkern,r2,c(maxspec)
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
      call update_zcoord(itime,i)
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
