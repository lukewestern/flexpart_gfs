! SPDX-FileCopyrightText: FLEXPART 1998-2019, see flexpart_license.txt
! SPDX-License-Identifier: GPL-3.0-or-later

subroutine boundcond_domainfill(itime,loutend)
  !                                  i      i
  !*****************************************************************************
  !                                                                            *
  ! Particles are created by this subroutine continuously throughout the       *
  ! simulation at the boundaries of the domain-filling box.                    *
  ! All particles carry the same amount of mass which alltogether comprises the*
  ! mass of air within the box, which remains (more or less) constant.         *
  !                                                                            *
  !     Author: A. Stohl                                                       *
  !                                                                            *
  !     16 October 2002                                                        *
  !                                                                            *
  !*****************************************************************************
  !                                                                            *
  ! Variables:                                                                 *
  !                                                                            *
  ! nx_we(2)       grid indices for western and eastern boundary of domain-    *
  !                filling trajectory calculations                             *
  ! ny_sn(2)       grid indices for southern and northern boundary of domain-  *
  !                filling trajectory calculations                             *
  !                                                                            *
  !*****************************************************************************

  use point_mod
  use par_mod
  use com_mod
  use random_mod, only: ran1
  use particle_mod
  use coordinates_ecmwf

  implicit none

  real :: dz,dz1,dz2,dt1,dt2,dtt,ylat,xm,cosfact,accmasst
  integer :: itime,in,indz,indzp,i,loutend
  integer :: j,k,ix,jy,m,indzh,indexh,minpart,ipart,mmass
  integer :: numactiveparticles

  real :: windl(2),rhol(2)
  real :: windhl(2),rhohl(2)
  real :: windx,rhox
  real :: deltaz,boundarea,fluxofmass

  integer :: ixm,ixp,jym,jyp,indzm,mm
  real :: pvpart,ddx,ddy,rddx,rddy,p1,p2,p3,p4,y1(2),yh1(2)

  integer :: idummy = -11


  ! If domain-filling is global, no boundary conditions are needed
  !***************************************************************

  if (gdomainfill) return

  accmasst=0.
  numactiveparticles=0

  ! Terminate trajectories that have left the domain, if domain-filling
  ! trajectory calculation domain is not global
  !********************************************************************

  do i=1,numpart
    if (.not. part(i)%alive) cycle

    if ((part(i)%ylat.gt.real(ny_sn(2))).or. &
         (part(i)%ylat.lt.real(ny_sn(1)))) call terminate_particle(i)
    if (((.not.xglobal).or.(nx_we(2).ne.(nx-2))).and. &
         ((part(i)%xlon.lt.real(nx_we(1))).or. &
         (part(i)%xlon.gt.real(nx_we(2))))) call terminate_particle(i)
    if (part(i)%alive) numactiveparticles = numactiveparticles+1
  end do


  ! Determine auxiliary variables for time interpolation
  !*****************************************************

  dt1=real(itime-memtime(1))
  dt2=real(memtime(2)-itime)
  dtt=1./(dt1+dt2)

  !***************************************
  ! Western and eastern boundary condition
  !***************************************

  ! Loop from south to north
  !*************************

  do jy=ny_sn(1),ny_sn(2)

  ! Loop over western (index 1) and eastern (index 2) boundary
  !***********************************************************

    do k=1,2

  ! Loop over all release locations in a column
  !********************************************

      do j=1,numcolumn_we(k,jy)

  ! Determine, for each release location, the area of the corresponding boundary
  !*****************************************************************************

        if (j.eq.1) then
          deltaz=(zcolumn_we(k,jy,2)+zcolumn_we(k,jy,1))/2.
        else if (j.eq.numcolumn_we(k,jy)) then
  ! In order to avoid taking a very high column for very many particles,
  ! use the deltaz from one particle below instead
          deltaz=(zcolumn_we(k,jy,j)-zcolumn_we(k,jy,j-2))/2.
        else
          deltaz=(zcolumn_we(k,jy,j+1)-zcolumn_we(k,jy,j-1))/2.
        endif
        if ((jy.eq.ny_sn(1)).or.(jy.eq.ny_sn(2))) then
          boundarea=deltaz*111198.5/2.*dy
        else
          boundarea=deltaz*111198.5*dy
        endif


  ! Interpolate the wind velocity and density to the release location
  !******************************************************************

  ! Determine the model level below the release position
  !*****************************************************

        do i=2,nz
          if (height(i).gt.zcolumn_we(k,jy,j)) then
            indz=i-1
            indzp=i
            exit
          endif
        end do

  ! Vertical distance to the level below and above current position
  !****************************************************************

        dz1=zcolumn_we(k,jy,j)-height(indz)
        dz2=height(indzp)-zcolumn_we(k,jy,j)
        dz=1./(dz1+dz2)

  ! Vertical and temporal interpolation
  !************************************

        do m=1,2
          indexh=memind(m)
          do in=1,2
            indzh=indz+in-1
            windl(in)=uu(nx_we(k),jy,indzh,indexh)
            rhol(in)=rho(nx_we(k),jy,indzh,indexh)
          end do

          windhl(m)=(dz2*windl(1)+dz1*windl(2))*dz
          rhohl(m)=(dz2*rhol(1)+dz1*rhol(2))*dz
        end do

        windx=(windhl(1)*dt2+windhl(2)*dt1)*dtt
        rhox=(rhohl(1)*dt2+rhohl(2)*dt1)*dtt

  ! Calculate mass flux
  !********************

        fluxofmass=windx*rhox*boundarea*real(lsynctime)


  ! If the mass flux is directed into the domain, add it to previous mass fluxes;
  ! if it is out of the domain, set accumulated mass flux to zero
  !******************************************************************************

        if (k.eq.1) then
          if (fluxofmass.ge.0.) then
            acc_mass_we(k,jy,j)=acc_mass_we(k,jy,j)+fluxofmass
          else
            acc_mass_we(k,jy,j)=0.
          endif
        else
          if (fluxofmass.le.0.) then
            acc_mass_we(k,jy,j)=acc_mass_we(k,jy,j)+abs(fluxofmass)
          else
            acc_mass_we(k,jy,j)=0.
          endif
        endif
        accmasst=accmasst+acc_mass_we(k,jy,j)

  ! If the accumulated mass exceeds half the mass that each particle shall carry,
  ! one (or more) particle(s) is (are) released and the accumulated mass is
  ! reduced by the mass of this (these) particle(s)
  !******************************************************************************

        if (acc_mass_we(k,jy,j).ge.xmassperparticle/2.) then
          mmass=int((acc_mass_we(k,jy,j)+xmassperparticle/2.)/ &
               xmassperparticle)
          acc_mass_we(k,jy,j)=acc_mass_we(k,jy,j)- &
               real(mmass)*xmassperparticle
        else
          mmass=0
        endif

        do m=1,mmass
          call get_new_part_index(ipart)
          call spawn_particle(itime, ipart)

  ! Assign particle positions
  !**************************

          call set_xlon(ipart,real(nx_we(k)))
          if (jy.eq.ny_sn(1)) then
            call set_ylat(ipart,real(jy)+0.5*ran1(idummy))
          else if (jy.eq.ny_sn(2)) then
            call set_ylat(ipart,real(jy)-0.5*ran1(idummy))
          else
            call set_ylat(ipart,real(jy)+(ran1(idummy)-.5))
          endif
          if (j.eq.1) then
            part(ipart)%z=zcolumn_we(k,jy,1)+(zcolumn_we(k,jy,2)- &
                 zcolumn_we(k,jy,1))/4.
          else if (j.eq.numcolumn_we(k,jy)) then
            part(ipart)%z=(2.*zcolumn_we(k,jy,j)+ &
                 zcolumn_we(k,jy,j-1)+height(nz))/4.
          else
            part(ipart)%z=zcolumn_we(k,jy,j-1)+ran1(idummy)* &
                 (zcolumn_we(k,jy,j+1)-zcolumn_we(k,jy,j-1))
          endif


          if (wind_coord_type.eq.'ETA') then
            call z_to_zeta(itime,part(ipart)%xlon,part(ipart)%ylat, &
              part(ipart)%z,part(ipart)%zeta)
            part(ipart)%etaupdate = .true. ! The z(meter) coordinate is up to date
          endif

  ! Interpolate PV to the particle position
  !****************************************
          ixm=int(part(ipart)%xlon)
          jym=int(part(ipart)%ylat)
          ixp=ixm+1
          jyp=jym+1
          ddx=part(ipart)%xlon-real(ixm)
          ddy=part(ipart)%ylat-real(jym)
          rddx=1.-ddx
          rddy=1.-ddy
          p1=rddx*rddy
          p2=ddx*rddy
          p3=rddx*ddy
          p4=ddx*ddy
          do i=2,nz
            if (height(i).gt.part(ipart)%z) then
              indzm=i-1
              indzp=i
              exit
            endif
          end do
          dz1=part(ipart)%z-height(indzm)
          dz2=height(indzp)-part(ipart)%z
          dz=1./(dz1+dz2)
          do mm=1,2
            indexh=memind(mm)
            do in=1,2
              indzh=indzm+in-1
              y1(in)=p1*pv(ixm,jym,indzh,indexh) &
                   +p2*pv(ixp,jym,indzh,indexh) &
                   +p3*pv(ixm,jyp,indzh,indexh) &
                   +p4*pv(ixp,jyp,indzh,indexh)
            end do
            yh1(mm)=(dz2*y1(1)+dz1*y1(2))*dz
          end do
          pvpart=(yh1(1)*dt2+yh1(2)*dt1)*dtt
          ylat=ylat0+part(ipart)%ylat*dy
          if (ylat.lt.0.) pvpart=-1.*pvpart


  ! For domain-filling option 2 (stratospheric O3), do the rest only in the stratosphere
  !*****************************************************************************

          if (((part(ipart)%z.gt.3000.).and. &
               (pvpart.gt.pvcrit)).or.(mdomainfill.eq.1)) then
            part(ipart)%nclass=min(int(ran1(idummy)* &
                 real(nclassunc))+1,nclassunc)
            numactiveparticles=numactiveparticles+1
            numparticlecount=numparticlecount+1
            part(ipart)%npoint=numparticlecount
            part(ipart)%idt=mintime
            part(ipart)%tstart=itime
            part(ipart)%mass(1)=xmassperparticle
            if (mdomainfill.eq.2) part(ipart)%mass(1)= &
                 part(ipart)%mass(1)*pvpart*48./29.*ozonescale/10.**9
          else
            stop 'boundcond_domainfill error: look into original to understand what should happen here'
          endif
        end do ! particles
      end do ! release locations in column
    end do ! western and eastern boundary
  end do ! south to north


  !*****************************************
  ! Southern and northern boundary condition
  !*****************************************

  ! Loop from west to east
  !***********************

  do ix=nx_we(1),nx_we(2)

  ! Loop over southern (index 1) and northern (index 2) boundary
  !*************************************************************

    do k=1,2
      ylat=ylat0+real(ny_sn(k))*dy
      cosfact=cos(ylat*pi180)

  ! Loop over all release locations in a column
  !********************************************

      do j=1,numcolumn_sn(k,ix)

  ! Determine, for each release location, the area of the corresponding boundary
  !*****************************************************************************

        if (j.eq.1) then
          deltaz=(zcolumn_sn(k,ix,2)+zcolumn_sn(k,ix,1))/2.
        else if (j.eq.numcolumn_sn(k,ix)) then
  !        deltaz=height(nz)-(zcolumn_sn(k,ix,j-1)+
  !    +        zcolumn_sn(k,ix,j))/2.
  ! In order to avoid taking a very high column for very many particles,
  ! use the deltaz from one particle below instead
          deltaz=(zcolumn_sn(k,ix,j)-zcolumn_sn(k,ix,j-2))/2.
        else
          deltaz=(zcolumn_sn(k,ix,j+1)-zcolumn_sn(k,ix,j-1))/2.
        endif
        if ((ix.eq.nx_we(1)).or.(ix.eq.nx_we(2))) then
          boundarea=deltaz*111198.5/2.*cosfact*dx
        else
          boundarea=deltaz*111198.5*cosfact*dx
        endif


  ! Interpolate the wind velocity and density to the release location
  !******************************************************************

  ! Determine the model level below the release position
  !*****************************************************

        do i=2,nz
          if (height(i).gt.zcolumn_sn(k,ix,j)) then
            indz=i-1
            indzp=i
            exit
          endif
        end do

  ! Vertical distance to the level below and above current position
  !****************************************************************

        dz1=zcolumn_sn(k,ix,j)-height(indz)
        dz2=height(indzp)-zcolumn_sn(k,ix,j)
        dz=1./(dz1+dz2)

  ! Vertical and temporal interpolation
  !************************************

        do m=1,2
          indexh=memind(m)
          do in=1,2
            indzh=indz+in-1
            windl(in)=vv(ix,ny_sn(k),indzh,indexh)
            rhol(in)=rho(ix,ny_sn(k),indzh,indexh)
          end do

          windhl(m)=(dz2*windl(1)+dz1*windl(2))*dz
          rhohl(m)=(dz2*rhol(1)+dz1*rhol(2))*dz
        end do

        windx=(windhl(1)*dt2+windhl(2)*dt1)*dtt
        rhox=(rhohl(1)*dt2+rhohl(2)*dt1)*dtt

  ! Calculate mass flux
  !********************

        fluxofmass=windx*rhox*boundarea*real(lsynctime)

  ! If the mass flux is directed into the domain, add it to previous mass fluxes;
  ! if it is out of the domain, set accumulated mass flux to zero
  !******************************************************************************

        if (k.eq.1) then
          if (fluxofmass.ge.0.) then
            acc_mass_sn(k,ix,j)=acc_mass_sn(k,ix,j)+fluxofmass
          else
            acc_mass_sn(k,ix,j)=0.
          endif
        else
          if (fluxofmass.le.0.) then
            acc_mass_sn(k,ix,j)=acc_mass_sn(k,ix,j)+abs(fluxofmass)
          else
            acc_mass_sn(k,ix,j)=0.
          endif
        endif
        accmasst=accmasst+acc_mass_sn(k,ix,j)

  ! If the accumulated mass exceeds half the mass that each particle shall carry,
  ! one (or more) particle(s) is (are) released and the accumulated mass is
  ! reduced by the mass of this (these) particle(s)
  !******************************************************************************

        if (acc_mass_sn(k,ix,j).ge.xmassperparticle/2.) then
          mmass=int((acc_mass_sn(k,ix,j)+xmassperparticle/2.)/ &
               xmassperparticle)
          acc_mass_sn(k,ix,j)=acc_mass_sn(k,ix,j)- &
               real(mmass)*xmassperparticle
        else
          mmass=0
        endif

        do m=1,mmass
          call get_new_part_index(ipart)
          call spawn_particle(itime, ipart)
  
  ! Assign particle positions
  !**************************
          call set_ylat(ipart,real(ny_sn(k)))
          if (ix.eq.nx_we(1)) then
            call set_xlon(ipart,real(ix)+0.5*ran1(idummy))
          else if (ix.eq.nx_we(2)) then
            call set_xlon(ipart,real(ix)-0.5*ran1(idummy))
          else
            call set_xlon(ipart,real(ix)+(ran1(idummy)-.5))
          endif
          if (j.eq.1) then
            part(ipart)%z=zcolumn_sn(k,ix,1)+(zcolumn_sn(k,ix,2)- &
                 zcolumn_sn(k,ix,1))/4.
          else if (j.eq.numcolumn_sn(k,ix)) then
            part(ipart)%z=(2.*zcolumn_sn(k,ix,j)+ &
                 zcolumn_sn(k,ix,j-1)+height(nz))/4.
          else
            part(ipart)%z=zcolumn_sn(k,ix,j-1)+ran1(idummy)* &
                 (zcolumn_sn(k,ix,j+1)-zcolumn_sn(k,ix,j-1))
          endif

          if (wind_coord_type.eq.'ETA') then
            call z_to_zeta(itime,part(ipart)%xlon,part(ipart)%ylat, &
              part(ipart)%z,part(ipart)%zeta)
            part(ipart)%etaupdate = .true. ! The z(meter) coordinate is up to date
          endif

  ! Interpolate PV to the particle position
  !****************************************
          ixm=int(part(ipart)%xlon)
          jym=int(part(ipart)%ylat)
          ixp=ixm+1
          jyp=jym+1
          ddx=part(ipart)%xlon-real(ixm)
          ddy=part(ipart)%ylat-real(jym)
          rddx=1.-ddx
          rddy=1.-ddy
          p1=rddx*rddy
          p2=ddx*rddy
          p3=rddx*ddy
          p4=ddx*ddy
          do i=2,nz
            if (height(i).gt.part(ipart)%z) then
              indzm=i-1
              indzp=i
              exit
            endif
          end do
          dz1=part(ipart)%z-height(indzm)
          dz2=height(indzp)-part(ipart)%z
          dz=1./(dz1+dz2)
          do mm=1,2
            indexh=memind(mm)
            do in=1,2
              indzh=indzm+in-1
              y1(in)=p1*pv(ixm,jym,indzh,indexh) &
                   +p2*pv(ixp,jym,indzh,indexh) &
                   +p3*pv(ixm,jyp,indzh,indexh) &
                   +p4*pv(ixp,jyp,indzh,indexh)
            end do
            yh1(mm)=(dz2*y1(1)+dz1*y1(2))*dz
          end do
          pvpart=(yh1(1)*dt2+yh1(2)*dt1)*dtt
          if (ylat.lt.0.) pvpart=-1.*pvpart


  ! For domain-filling option 2 (stratospheric O3), do the rest only in the stratosphere
  !*****************************************************************************

          if (((part(ipart)%z.gt.3000.).and. &
               (pvpart.gt.pvcrit)).or.(mdomainfill.eq.1)) then
            part(ipart)%nclass=min(int(ran1(idummy)* &
                 real(nclassunc))+1,nclassunc)
            numactiveparticles=numactiveparticles+1
            numparticlecount=numparticlecount+1
            part(ipart)%npoint=numparticlecount
            part(ipart)%idt=mintime
            part(ipart)%mass(1)=xmassperparticle
            if (mdomainfill.eq.2) part(ipart)%mass(1)= &
                 part(ipart)%mass(1)*pvpart*48./29.*ozonescale/10.**9
          else
            stop 'boundcond_domainfill error: look into original to understand what should happen here'
          endif
        end do ! particles
      end do ! releases per column
    end do ! east west
  end do ! north south

  ! If particles shall be dumped, then accumulated masses at the domain boundaries
  ! must be dumped, too, to be used for later runs
  !*****************************************************************************

  if ((ipout.gt.0).and.(itime.eq.loutend)) then
    open(unitboundcond,file=path(2)(1:length(2))//'boundcond.bin', &
         form='unformatted')
    write(unitboundcond) numcolumn_we,numcolumn_sn, &
         zcolumn_we,zcolumn_sn,acc_mass_we,acc_mass_sn
    close(unitboundcond)
  endif

end subroutine boundcond_domainfill
