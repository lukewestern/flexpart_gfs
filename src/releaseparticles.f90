! SPDX-FileCopyrightText: FLEXPART 1998-2019, see flexpart_license.txt
! SPDX-License-Identifier: GPL-3.0-or-later

subroutine releaseparticles(itime)
  !                              o
  !*****************************************************************************
  !                                                                            *
  !     This subroutine releases particles from the release locations.         *
  !                                                                            *
  !     It searches for a "vacant" storage space and assigns all particle      *
  !     information to that space. A space is vacant either when no particle   *
  !     is yet assigned to it, or when it's particle is expired and, thus,     *
  !     the storage space is made available to a new particle.                 *
  !                                                                            *
  !     Author: A. Stohl                                                       *
  !                                                                            *
  !     29 June 2002                                                           *
  !                                                                            *
  !*****************************************************************************
  !                                                                            *
  ! Variables:                                                                 *
  ! itime [s]            current time                                          *
  ! ireleasestart, ireleaseend          start and end times of all releases    *
  ! npart(maxpoint)      number of particles to be released in total           *
  ! numrel               number of particles to be released during this time   *
  !                      step                                                  *
  !                                                                            *
  !*****************************************************************************

  use point_mod
  use xmass_mod
  use par_mod
  use com_mod
  use windfields_mod
  use random_mod, only: ran1
  use interpol_mod
  use coordinates_ecmwf

  implicit none

  !real xaux,yaux,zaux,ran1,rfraction,xmasssave(maxpoint)
  real :: xaux,yaux,zaux,rfraction
  real :: topo,rhoaux(2),r,t,rhoout
  real :: dz1,dz2,dz,xlonav,timecorrect(maxspec),press,pressold
  real :: presspart,average_timecorrect
  integer :: itime,numrel,i,j,k,n,ipart,minpart,ii
  integer :: kz
  integer :: nweeks,ndayofweek,nhour,jjjjmmdd,ihmmss,mm
  real(kind=dp) :: juldate,julmonday,jul,jullocal,juldiff
  real,parameter :: eps=nxmax/3.e5,eps2=1.e-6

  integer :: idummy = -7
  !save idummy,xmasssave
  !data idummy/-7/,xmasssave/maxpoint*0./

  real :: frac,psint,zzlev,zzlev2,ttemp


  ! Determine the actual date and time in Greenwich (i.e., UTC + correction for daylight savings time)
  !*****************************************************************************

  julmonday=juldate(19000101,0)          ! this is a Monday
  jul=bdate+real(itime,kind=dp)/86400._dp    ! this is the current day
  call caldate(jul,jjjjmmdd,ihmmss)
  mm=(jjjjmmdd-10000*(jjjjmmdd/10000))/100
  if ((mm.ge.4).and.(mm.le.9)) jul=jul+1._dp/24._dp   ! daylight savings time in summer


  ! For every release point, check whether we are in the release time interval
  !***************************************************************************
  ! First allocate all particles that are going to be in the simulation (this could be done differently...)
  if (itime.eq.0) then
    do i=1,numpoint 
      call allocate_particles(npart(i))
    end do 
  end if 

  minpart=1
  do i=1,numpoint
    if ((itime.ge.ireleasestart(i)).and. &! are we within release interval?
         (itime.le.ireleaseend(i))) then

  ! Determine the local day and time
  !*********************************

      xlonav=xlon0+(xpoint2(i)+xpoint1(i))/2.*dx  ! longitude needed to determine local time
      if (xlonav.lt.-180.) xlonav=xlonav+360.
      if (xlonav.gt.180.) xlonav=xlonav-360.
      jullocal=jul+real(xlonav,kind=dp)/360._dp   ! correct approximately for time zone to obtain local time

      juldiff=jullocal-julmonday
      nweeks=int(juldiff/7._dp)
      juldiff=juldiff-real(nweeks,kind=dp)*7._dp
      ndayofweek=int(juldiff)+1              ! this is the current day of week, starting with Monday
      nhour=nint((juldiff-real(ndayofweek-1,kind=dp))*24._dp)    ! this is the current hour
      if (nhour.eq.0) then
        nhour=24
        ndayofweek=ndayofweek-1
        if (ndayofweek.eq.0) ndayofweek=7
      endif

  ! Calculate a species- and time-dependent correction factor, distinguishing between
  ! area (those with release starting at surface) and point (release starting above surface) sources
  ! Also, calculate an average time correction factor (species independent)
  !*****************************************************************************
      average_timecorrect=0.
      do k=1,nspec
        if(abs(xpoint2(i)-xpoint1(i)).lt.1.E-4.and.abs(ypoint2(i)-ypoint1(i)).lt.1.E-4) then
!        if (zpoint1(i).gt.0.5) then      ! point source
          timecorrect(k)=point_hour(k,nhour)*point_dow(k,ndayofweek)
        else                             ! area source
          timecorrect(k)=area_hour(k,nhour)*area_dow(k,ndayofweek)
        endif
        average_timecorrect=average_timecorrect+timecorrect(k)
      end do
      average_timecorrect=average_timecorrect/real(nspec)

  ! Determine number of particles to be released this time; at start and at end of release,
  ! only half the particles are released
  !*****************************************************************************

      if (ireleasestart(i).ne.ireleaseend(i)) then
        rfraction=abs(real(npart(i))*real(lsynctime)/ &
             real(ireleaseend(i)-ireleasestart(i)))
        if ((itime.eq.ireleasestart(i)).or. &
             (itime.eq.ireleaseend(i))) rfraction=rfraction/2.

  ! Take the species-average time correction factor in order to scale the
  ! number of particles released this time
  !**********************************************************************
        rfraction=rfraction*average_timecorrect

        rfraction=rfraction+xmasssave(i)  ! number to be released at this time
        numrel=int(rfraction)
        xmasssave(i)=rfraction-real(numrel)
      else
        numrel=npart(i)
      endif

      xaux=xpoint2(i)-xpoint1(i)
      yaux=ypoint2(i)-ypoint1(i)
      zaux=zpoint2(i)-zpoint1(i)

      do j=1,numrel                       ! loop over particles to be released this time
        call get_new_part_index(ipart)
        call spawn_particle(itime, ipart)

  ! Particle coordinates are determined by using a random position within the release volume
  !*****************************************************************************

  ! Determine horizontal particle position
  !***************************************
        call set_xlon(ipart,real(xpoint1(i)+ran1(idummy)*xaux,kind=dp))
        if (xglobal) then
          if (part(ipart)%xlon.gt.real(nxmin1,kind=dp)) &
            call set_xlon(ipart,-real(nxmin1,kind=dp))
          if (part(ipart)%xlon.lt.0.) &
            call set_xlon(ipart,real(nxmin1,kind=dp))
        endif
        call set_ylat(ipart,real(ypoint1(i)+ran1(idummy)*yaux,kind=dp))

  ! Assign mass to particle: Total mass divided by total number of particles.
  ! Time variation has partly been taken into account already by a species-average
  ! correction factor, by which the number of particles released this time has been
  ! scaled. Adjust the mass per particle by the species-dependent time correction factor
  ! divided by the species-average one
  ! for the scavenging calculation the mass needs to be multiplied with rho of the particle layer and
  ! divided by the sum of rho of all particles.
  !*****************************************************************************
        do k=1,nspec
          part(ipart)%mass(k)=xmass(i,k)/real(npart(i)) &
                *timecorrect(k)/average_timecorrect
          if (DRYBKDEP.or.WETBKDEP) then ! if there is no scavenging in wetdepo it will be set to 0
!              if ( henry(k).gt.0 .or. &
!                   crain_aero(k).gt.0. .or. csnow_aero(k).gt.0. .or. &
!                   ccn_aero(k).gt.0. .or. in_aero(k).gt.0. )  then
            xscav_frac1(ipart,k)=-1.
          endif
  ! Assign certain properties to particle
  !**************************************
        end do
        part(ipart)%nclass=min(int(ran1(idummy)*real(nclassunc))+1, &
             nclassunc)
        numparticlecount=numparticlecount+1
        if (mquasilag.eq.0) then
          part(ipart)%npoint=i
        else
          part(ipart)%npoint=numparticlecount
        endif
        part(ipart)%idt=mintime               ! first time step

  ! Determine vertical particle position
  !*************************************
        call set_z(ipart,zpoint1(i)+ran1(idummy)*zaux)
  ! Interpolation of topography and density
  !****************************************

  ! Determine the nest we are in
  !*****************************

        ngrid=0
        do k=numbnests,1,-1
          if ((part(ipart)%xlon.gt.xln(k)+eps).and. &
               (part(ipart)%xlon.lt.xrn(k)-eps).and. &
               (part(ipart)%xlon.gt.yln(k)+eps).and. &
               (part(ipart)%xlon.lt.yrn(k)-eps)) then
            ngrid=k
            exit
          endif
        end do

  ! Determine (nested) grid coordinates and auxiliary parameters used for interpolation
  !*****************************************************************************

        if (ngrid.gt.0) then
          xtn=(part(ipart)%xlon-xln(ngrid))*xresoln(ngrid)
          ytn=(part(ipart)%ylat-yln(ngrid))*yresoln(ngrid)
          ix=int(xtn)
          jy=int(ytn)
          ddy=ytn-real(jy)
          ddx=xtn-real(ix)
        else
          ix=int(part(ipart)%xlon)
          jy=int(part(ipart)%ylat)
          ddy=part(ipart)%ylat-real(jy)
          ddx=part(ipart)%xlon-real(ix)
        endif
        ixp=ix+1
        jyp=jy+1
        rddx=1.-ddx
        rddy=1.-ddy
        p1=rddx*rddy
        p2=ddx*rddy
        p3=rddx*ddy
        p4=ddx*ddy

        if (ngrid.gt.0) then
          topo=p1*oron(ix ,jy ,ngrid) &
               + p2*oron(ixp,jy ,ngrid) &
               + p3*oron(ix ,jyp,ngrid) &
               + p4*oron(ixp,jyp,ngrid)
        else
          topo=p1*oro(ix ,jy) &
               + p2*oro(ixp,jy) &
               + p3*oro(ix ,jyp) &
               + p4*oro(ixp,jyp)
        endif

  ! If starting height is in pressure coordinates, retrieve pressure profile and convert zpart1 to meters
  !*****************************************************************************
        if (kindz(i).eq.3) then
          presspart=part(ipart)%z
          do kz=1,nz
            if (ngrid.gt.0) then
              r=p1*rhon(ix ,jy ,kz,2,ngrid) &
                   +p2*rhon(ixp,jy ,kz,2,ngrid) &
                   +p3*rhon(ix ,jyp,kz,2,ngrid) &
                   +p4*rhon(ixp,jyp,kz,2,ngrid)
              t=p1*ttn(ix ,jy ,kz,2,ngrid) &
                   +p2*ttn(ixp,jy ,kz,2,ngrid) &
                   +p3*ttn(ix ,jyp,kz,2,ngrid) &
                   +p4*ttn(ixp,jyp,kz,2,ngrid)
            else
              r=p1*rho(ix ,jy ,kz,2) &
                   +p2*rho(ixp,jy ,kz,2) &
                   +p3*rho(ix ,jyp,kz,2) &
                   +p4*rho(ixp,jyp,kz,2)
              t=p1*tt(ix ,jy ,kz,2) &
                   +p2*tt(ixp,jy ,kz,2) &
                   +p3*tt(ix ,jyp,kz,2) &
                   +p4*tt(ixp,jyp,kz,2)
            endif
            press=r*r_air*t/100.
            if (kz.eq.1) pressold=press

            if (press.lt.presspart) then
              if (kz.eq.1) then
                call set_z(ipart,height(1)/2.)
              else
                dz1=pressold-presspart
                dz2=presspart-press
                call set_z(ipart,(height(kz-1)*dz2+height(kz)*dz1) &
                     /(dz1+dz2))
              endif
              exit
            endif
            pressold=press
          end do
        endif


  ! If release positions are given in meters above sea level, subtract the
  ! topography from the starting height
  !***********************************************************************

        if (kindz(i).eq.2) call update_z(ipart,-topo)
        if (part(ipart)%z.lt.eps2) call set_z(ipart,eps2)   ! Minimum starting height is eps2
        if (part(ipart)%z.gt.height(nz)-0.5) &
          call set_z(ipart,height(nz)-0.5) ! Maximum starting height is uppermost level - 0.5 meters

        if (wind_coord_type.eq.'ETA') then
          call z_to_zeta(itime,part(ipart)%xlon,part(ipart)%ylat,part(ipart)%z,part(ipart)%zeta)
          part(ipart)%etaupdate = .true. ! The z(meter) coordinate is up to date
        end if

  ! For special simulations, multiply particle concentration air density;
  ! Simply take the 2nd field in memory to do this (accurate enough)
  !***********************************************************************
  !AF IND_SOURCE switches between different units for concentrations at the source
  !Af    NOTE that in backward simulations the release of particles takes place at the
  !Af         receptor and the sampling at the source.
  !Af          1="mass"
  !Af          2="mass mixing ratio"
  !Af IND_RECEPTOR switches between different units for concentrations at the receptor
  !Af          1="mass"
  !Af          2="mass mixing ratio"

  !Af switches for the releasefile:
  !Af IND_REL =  1 : xmass * rho
  !Af IND_REL =  0 : xmass * 1

  !Af ind_rel is defined in readcommand.f

        if ((ind_rel .eq. 1).or.(ind_rel .eq. 3).or.(ind_rel .eq. 4)) then

  ! Interpolate the air density
  !****************************

          do ii=2,nz
            if (height(ii).gt.part(ipart)%z) then
              indz=ii-1
              indzp=ii
              exit
            endif
          end do

          dz1=part(ipart)%z-height(indz)
          dz2=height(indzp)-part(ipart)%z
          dz=1./(dz1+dz2)

          if (ngrid.gt.0) then
            do n=1,2
              rhoaux(n)=p1*rhon(ix ,jy ,indz+n-1,2,ngrid) &
                   +p2*rhon(ixp,jy ,indz+n-1,2,ngrid) &
                   +p3*rhon(ix ,jyp,indz+n-1,2,ngrid) &
                   +p4*rhon(ixp,jyp,indz+n-1,2,ngrid)
            end do
          else
            do n=1,2
              rhoaux(n)=p1*rho(ix ,jy ,indz+n-1,2) &
                   +p2*rho(ixp,jy ,indz+n-1,2) &
                   +p3*rho(ix ,jyp,indz+n-1,2) &
                   +p4*rho(ixp,jyp,indz+n-1,2)
            end do
          endif
          rhoout=(dz2*rhoaux(1)+dz1*rhoaux(2))*dz
          rho_rel(i)=rhoout


  ! Multiply "mass" (i.e., mass mixing ratio in forward runs) with density
  !********************************************************************

          do k=1,nspec
            part(ipart)%mass(k)=part(ipart)%mass(k)*rhoout
          end do
        endif

        call get_total_part_num(numpart)

      end do  ! numrel 
    endif ! releasepoint
  end do ! numpoint


  return

996   continue
  write(*,*) '#####################################################'
  write(*,*) '#### FLEXPART MODEL SUBROUTINE RELEASEPARTICLES: ####'
  write(*,*) '####                                             ####'
  write(*,*) '#### ERROR - TOTAL NUMBER OF PARTICLES REQUIRED  ####'
  write(*,*) '#### EXCEEDS THE MAXIMUM ALLOWED NUMBER. REDUCE  ####'
  write(*,*) '#### EITHER NUMBER OF PARTICLES PER RELEASE POINT####'
  write(*,*) '#### OR REDUCE NUMBER OF RELEASE POINTS.         ####'
  write(*,*) '#####################################################'
  stop

end subroutine releaseparticles
