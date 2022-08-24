module wetdepo_mod
  use point_mod
  use par_mod
  use com_mod
  use particle_mod

  implicit none

contains

subroutine wetdepo(itime,ltsample,loutnext)
  !                  i      i        i
  !*****************************************************************************
  !                                                                            *
  ! Calculation of wet deposition using the concept of scavenging coefficients.*
  ! For lack of detailed information, washout and rainout are jointly treated. *
  ! It is assumed that precipitation does not occur uniformly within the whole *
  ! grid cell, but that only a fraction of the grid cell experiences rainfall. *
  ! This fraction is parameterized from total cloud cover and rates of large   *
  ! scale and convective precipitation.                                        *
  !                                                                            *
  !    Author: A. Stohl                                                        *
  !                                                                            *
  !    1 December 1996                                                         *
  !                                                                            *
  ! Correction by Petra Seibert, Sept 2002:                                    *
  ! use centred precipitation data for integration                             *
  ! Code may not be correct for decay of deposition!                           *
  !                                                                            *
  !*****************************************************************************
  !                                                                            *
  ! Variables:                                                                 *
  ! ix,jy              indices of output grid cell for each particle           *
  ! itime [s]          actual simulation time [s]                              *
  ! jpart              particle index                                          *
  ! ldeltat [s]        interval since radioactive decay was computed           *
  ! loutnext [s]       time for which gridded deposition is next output        *
  ! loutstep [s]       interval at which gridded deposition is output          *
  ! ltsample [s]       interval over which mass is deposited                   *
  ! wetdeposit         mass that is wet deposited                              *
  ! wetgrid            accumulated deposited mass on output grid               *
  ! wetscav            scavenging coefficient                                  *
  !                                                                            *
  ! Constants:                                                                 *
  !                                                                            *
  !*****************************************************************************
#ifdef _OPENMP
  use omp_lib
#endif
  use unc_mod

  implicit none

  integer :: jpart,itime,ltsample,loutnext,ldeltat
  integer :: itage,nage,ithread,thread
  integer :: ks, kp
  integer(selected_int_kind(16)), dimension(nspec) :: blc_count, inc_count
  real :: grfraction(3),wetscav
  real :: wetdeposit(maxspec),restmass
  real,parameter :: smallnum = tiny(0.0) ! smallest number that can be handled

  ! Compute interval since radioactive decay of deposited mass was computed
  !************************************************************************

  if (itime.le.loutnext) then
    ldeltat=itime-(loutnext-loutstep)
  else                                  ! first half of next interval
    ldeltat=itime-loutnext
  endif

  ! Loop over all particles
  !************************
  blc_count(:)=0
  inc_count(:)=0

  ! OMP doesn't work yet, a reduction is necessary for the kernel function
!$OMP PARALLEL PRIVATE(jpart,itage,nage,ks,kp,thread,wetscav,wetdeposit, &
!$OMP restmass, grfraction) REDUCTION(+:blc_count,inc_count)

#if (defined _OPENMP)
    thread = OMP_GET_THREAD_NUM() ! Starts with 0
#else
    thread = 1
#endif

!$OMP DO 
  do jpart=1,numpart

    ! Check if memory has been deallocated
    if (.not. is_particle_allocated(jpart)) cycle

    ! Check if particle is still allive
    if (.not. part(jpart)%alive) cycle

  ! Determine age class of the particle - nage is used for the kernel
  !******************************************************************
    itage=abs(itime-part(jpart)%tstart)
    do nage=1,nageclass
     if (itage.lt.lage(nage)) exit
    end do

    do ks=1,nspec      ! loop over species

      if (WETDEPSPEC(ks).eqv..false.) cycle 

  !**************************************************
  ! CALCULATE DEPOSITION 
  !**************************************************

      call get_wetscav(itime,ltsample,loutnext,jpart,ks,grfraction,inc_count,blc_count,wetscav) ! OMP carefully check

      if (WETBKDEP) then
        if ((xscav_frac1(jpart,ks).lt.-0.1)) then   ! particle marked as starting particle
          if (wetscav.gt.0.) then
             xscav_frac1(jpart,ks)=wetscav*(zpoint2(part(jpart)%npoint)-&
             zpoint1(part(jpart)%npoint))*grfraction(1)
          else
            part(jpart)%mass(ks)=0.
            xscav_frac1(jpart,ks)=0.
          endif
        endif
      endif

      if (wetscav.gt.0.) then
        wetdeposit(ks)=part(jpart)%mass(ks)* &
             (1.-exp(-wetscav*abs(ltsample)))*grfraction(1)  ! wet deposition
      else ! if no scavenging
        wetdeposit(ks)=0.
      endif
      part(jpart)%wetdepo(ks)=part(jpart)%wetdepo(ks)+wetdeposit(ks)
      restmass = part(jpart)%mass(ks)-wetdeposit(ks)
      if (ioutputforeachrelease.eq.1) then
        kp=part(jpart)%npoint
      else
        kp=1
      endif
      if (restmass .gt. smallnum) then
        part(jpart)%mass(ks)=restmass
  !   depostatistic
  !   wetdepo_sum(ks,kp)=wetdepo_sum(ks,kp)+wetdeposit(ks)
  !   depostatistic
      else
        part(jpart)%mass(ks)=0.
      endif
  !   Correct deposited mass to the last time step when radioactive decay of
  !   gridded deposited mass was calculated
      if (decay(ks).gt.0.) then
        wetdeposit(ks)=wetdeposit(ks)*exp(abs(ldeltat)*decay(ks))
      endif

  !    endif ! no deposition
    end do ! loop over species

  ! Sabine Eckhardt, June 2008 create deposition runs only for forward runs
  ! Add the wet deposition to accumulated amount on output grid and nested output grid
  !*****************************************************************************

    if ((ldirect.eq.1).and.(iout.ne.0)) then !OMP reduction necessary for wetgridunc
      call wetdepokernel(part(jpart)%nclass,wetdeposit,real(part(jpart)%xlon), &
           real(part(jpart)%ylat),nage,kp,thread+1)
      if (nested_output.eq.1) call wetdepokernel_nest(part(jpart)%nclass, &
           wetdeposit,real(part(jpart)%xlon),real(part(jpart)%ylat),nage,kp,thread+1)
    endif

  end do ! all particles

!$OMP END DO
!$OMP END PARALLEL

#ifdef _OPENMP
    if ((ldirect.eq.1).and.(iout.ne.0)) then
      do ithread=1,numthreads
        wetgridunc(:,:,:,:,:,:)=wetgridunc(:,:,:,:,:,:)+wetgridunc_omp(:,:,:,:,:,:,ithread)
        wetgridunc_omp(:,:,:,:,:,:,ithread)=0.
      end do
      if (nested_output.eq.1) then
        do ithread=1,numthreads
          wetgriduncn(:,:,:,:,:,:)=wetgriduncn(:,:,:,:,:,:)+wetgriduncn_omp(:,:,:,:,:,:,ithread)
          wetgriduncn_omp(:,:,:,:,:,:,ithread)=0.
        end do
      endif
    endif
#endif
  !write(*,*) 'WETGRIDUNC:',sum(wetgridunc),wetgridunc(20,270,1,1,1,1),wetgridunc(19,269,1,1,1,1)
  ! count the total number of below-cloud and in-cloud occurences:
  tot_blc_count(1:nspec)=tot_blc_count(1:nspec)+blc_count(1:nspec)
  tot_inc_count(1:nspec)=tot_inc_count(1:nspec)+inc_count(1:nspec)
end subroutine wetdepo

subroutine get_wetscav(itime,ltsample,loutnext,jpart,ks,grfraction,inc_count,blc_count,wetscav)
  !                          i      i        i     i   i    o           o          o       o
  !*****************************************************************************
  !                                                                            *
  ! Calculation of wet deposition using the concept of scavenging coefficients.*
  ! For lack of detailed information, washout and rainout are jointly treated. *
  ! It is assumed that precipitation does not occur uniformly within the whole *
  ! grid cell, but that only a fraction of the grid cell experiences rainfall. *
  ! This fraction is parameterized from total cloud cover and rates of large   *
  ! scale and convective precipitation.                                        *
  !                                                                            *
  !    Author: A. Stohl                                                        *
  !                                                                            *
  !    1 December 1996                                                         *
  !                                                                            *
  ! Correction by Petra Seibert, Sept 2002:                                    *
  ! use centred precipitation data for integration                             *
  ! Code may not be correct for decay of deposition!                           *
  !                                                                            *
  !*****************************************************************************
  !                                                                            *
  ! Variables:                                                                 *
  ! cc [0-1]           total cloud cover                                       *
  ! convp [mm/h]       convective precipitation rate                           *
  ! grfraction [0-1]   fraction of grid, for which precipitation occurs        *
  ! ix,jy              indices of output grid cell for each particle           *
  ! itime [s]          actual simulation time [s]                              *
  ! jpart              particle index                                          *
  ! lfr, cfr           area fraction covered by precipitation for large scale  *
  !                    and convective precipitation (dependent on prec. rate)  *
  ! loutnext [s]       time for which gridded deposition is next output        *
  ! loutstep [s]       interval at which gridded deposition is output          *
  ! lsp [mm/h]         large scale precipitation rate                          *
  ! ltsample [s]       interval over which mass is deposited                   *
  ! prec [mm/h]        precipitation rate in subgrid, where precipitation occurs*
  ! wetgrid            accumulated deposited mass on output grid               *
  ! wetscav            scavenging coefficient                                  *
  !                                                                            *
  ! Constants:                                                                 *
  !                                                                            *
  !*****************************************************************************

  use interpol_mod, only: interpol_rain, interpol_rain_nests
  use windfields_mod
  use coordinates_ecmwf

  implicit none

  integer :: jpart,itime,ltsample,loutnext,i,j,ix,jy
  integer :: ngrid,hz,il,interp_time, n
  integer(kind=1) :: clouds_v
  integer :: ks, kp
  integer(selected_int_kind(16)), dimension(nspec) :: blc_count, inc_count

  !  integer :: n1,n2, icbot,ictop, indcloud !TEST
  real :: S_i, act_temp, cl, cle ! in cloud scavenging
  real :: clouds_h ! cloud height for the specific grid point
  real :: xtn,ytn,lsp,convp,cc,grfraction(3),prec(3),wetscav,totprec
  real :: restmass
  real,parameter :: smallnum = tiny(0.0) ! smallest number that can be handled
  !save lfr,cfr

  real, parameter :: lfr(5) = (/ 0.5,0.65,0.8,0.9,0.95/)
  real, parameter :: cfr(5) = (/ 0.4,0.55,0.7,0.8,0.9 /)

  !ZHG aerosol below-cloud scavenging removal polynomial constants for rain and snow
  real, parameter :: bclr(6) = (/274.35758, 332839.59273, 226656.57259, 58005.91340, 6588.38582, 0.244984/) !rain (Laakso et al 2003)
  real, parameter :: bcls(6) = (/22.7, 0.0, 0.0, 1321.0, 381.0, 0.0/) !now (Kyro et al 2009)
  real :: frac_act, liq_frac, ice_frac, dquer_m

  real    :: Si_dummy, wetscav_dummy
  logical :: readclouds_this_nest


  wetscav=0.

  ! Determine which nesting level to be used
  !*****************************************
  ngrid=0
  do j=numbnests,1,-1
    if ((part(jpart)%xlon.gt.xln(j)).and.(part(jpart)%xlon.lt.xrn(j)).and. &
         (part(jpart)%ylat.gt.yln(j)).and.(part(jpart)%ylat.lt.yrn(j))) then
      ngrid=j
      exit
    endif
  end do

  ! Determine nested grid coordinates
  !**********************************
  readclouds_this_nest=.false.

  if (ngrid.gt.0) then
    xtn=(part(jpart)%xlon-xln(ngrid))*xresoln(ngrid)
    ytn=(part(jpart)%ylat-yln(ngrid))*yresoln(ngrid)
    ix=int(xtn)
    jy=int(ytn)
    if (readclouds_nest(ngrid)) readclouds_this_nest=.true.
  else
    ix=int(part(jpart)%xlon)
    jy=int(part(jpart)%ylat)
  endif

  ! Interpolate large scale precipitation, convective precipitation and
  ! total cloud cover
  ! Note that interpolated time refers to itime-0.5*ltsample [PS]
  !********************************************************************
  interp_time=nint(itime-0.5*ltsample) 

  n=memind(2)
  if (abs(memtime(1)-interp_time).lt.abs(memtime(2)-interp_time)) &
       n=memind(1)

  if (ngrid.eq.0) then
    call interpol_rain(lsprec,convprec,tcc,nxmax,nymax, &
         1,nx,ny,n,real(part(jpart)%xlon),real(part(jpart)%ylat),1, &
         memtime(1),memtime(2),interp_time,lsp,convp,cc)
  else
    call interpol_rain_nests(lsprecn,convprecn,tccn, &
         nxmaxn,nymaxn,1,maxnests,ngrid,nxn,nyn,n,xtn,ytn,1, &
         memtime(1),memtime(2),interp_time,lsp,convp,cc)
  endif

  !  If total precipitation is less than 0.01 mm/h - no scavenging occurs
  if ((lsp.lt.0.01).and.(convp.lt.0.01)) return

  ! get the level were the actual particle is in
  call update_zeta_to_z(itime,jpart)
  do il=2,nz
    if (height(il).gt.part(jpart)%z) then
      hz=il-1
      exit
    endif
  end do

  if (ngrid.eq.0) then
    clouds_v=clouds(ix,jy,hz,n)
    clouds_h=cloudsh(ix,jy,n)
  else
    clouds_v=cloudsn(ix,jy,hz,n,ngrid)
    clouds_h=cloudshn(ix,jy,n,ngrid)
  endif

  ! if there is no precipitation or the particle is above the clouds no
  ! scavenging is done

  if (clouds_v.le.1) return

  ! 1) Parameterization of the the area fraction of the grid cell where the
  !    precipitation occurs: the absolute limit is the total cloud cover, but
  !    for low precipitation rates, an even smaller fraction of the grid cell
  !    is used. Large scale precipitation occurs over larger areas than
  !    convective precipitation.
  !**************************************************************************

  if (lsp.gt.20.) then
    i=5
  else if (lsp.gt.8.) then
    i=4
  else if (lsp.gt.3.) then
    i=3
  else if (lsp.gt.1.) then
    i=2
  else
    i=1
  endif

  if (convp.gt.20.) then
    j=5
  else if (convp.gt.8.) then
    j=4
  else if (convp.gt.3.) then
    j=3
  else if (convp.gt.1.) then
    j=2
  else
    j=1
  endif


  !ZHG oct 2014 : Calculated for 1) both 2) lsp 3) convp - 2 and 3 not used removed by SE
  ! Tentatively differentiate the grfraction for lsp and convp for treating differently the two forms
  ! for now they are treated the same
  grfraction(1)=max(0.05,cc*(lsp*lfr(i)+convp*cfr(j))/(lsp+convp))

  ! 2) Computation of precipitation rate in sub-grid cell
  !******************************************************
  prec(1)=(lsp+convp)/grfraction(1)

  ! 3) Computation of scavenging coefficients for all species
  !    Computation of wet deposition
  !**********************************************************

  if (ngrid.gt.0) then
    act_temp=ttn(ix,jy,hz,n,ngrid)
  else
    act_temp=tt(ix,jy,hz,n)
  endif

  !***********************
  ! BELOW CLOUD SCAVENGING
  !***********************  
  if (clouds_v.ge.4) then !below cloud

  ! For gas: if positive below-cloud parameters (A or B), and dquer<=0
  !******************************************************************
    if ((dquer(ks).le.0.).and.(weta_gas(ks).gt.0..or.wetb_gas(ks).gt.0.)) then
      blc_count(ks)=blc_count(ks)+1
      wetscav=weta_gas(ks)*prec(1)**wetb_gas(ks)

  ! For aerosols: if positive below-cloud parameters (Crain/Csnow or B), and dquer>0
  !*********************************************************************************
    else if ((dquer(ks).gt.0.).and.(crain_aero(ks).gt.0..or.csnow_aero(ks).gt.0.)) then
      blc_count(ks)=blc_count(ks)+1

  !NIK 17.02.2015
  ! For the calculation here particle size needs to be in meter and not um as dquer is
  ! changed in readreleases
  ! For particles larger than 10 um use the largest size defined in the parameterizations (10um)
      dquer_m=min(10.,dquer(ks))/1000000. !conversion from um to m

  ! Rain:
      if (act_temp .ge. 273. .and. crain_aero(ks).gt.0.)  then

  ! ZHG 2014 : Particle RAIN scavenging coefficient based on Laakso et al 2003, 
  ! the below-cloud scavenging (rain efficienty) parameter Crain (=crain_aero) from SPECIES file
        wetscav=crain_aero(ks)*10**(bclr(1)+(bclr(2)*(log10(dquer_m))**(-4))+ &
             & (bclr(3)*(log10(dquer_m))**(-3))+ (bclr(4)*(log10(dquer_m))**(-2))+&
             &(bclr(5)*(log10(dquer_m))**(-1))+bclr(6)* (prec(1))**(0.5))

  ! Snow:
      elseif (act_temp .lt. 273. .and. csnow_aero(ks).gt.0.)  then 
  ! ZHG 2014 : Particle SNOW scavenging coefficient based on Kyro et al 2009, 
  ! the below-cloud scavenging (Snow efficiency) parameter Csnow (=csnow_aero) from SPECIES file
        wetscav=csnow_aero(ks)*10**(bcls(1)+(bcls(2)*(log10(dquer_m))**(-4))+&
             &(bcls(3)*(log10(dquer_m))**(-3))+ (bcls(4)*(log10(dquer_m))**(-2))+&
             &(bcls(5)*(log10(dquer_m))**(-1))+ bcls(6)* (prec(1))**(0.5))

      endif
            
    endif ! gas or particle
  !      endif ! positive below-cloud scavenging parameters given in Species file
  endif !end BELOW

  !********************
  ! IN CLOUD SCAVENGING
  !********************
  if (clouds_v.lt.4) then ! In-cloud
  ! NIK 13 may 2015: only do incloud if positive in-cloud scavenging parameters are
  ! given in species file, or if gas and positive Henry's constant
    if ((ccn_aero(ks).gt.0. .or. in_aero(ks).gt.0.).or.(henry(ks).gt.0.and.dquer(ks).le.0)) then 
      inc_count(ks)=inc_count(ks)+1
  ! if negative coefficients (turned off) set to zero for use in equation
      if (ccn_aero(ks).lt.0.) ccn_aero(ks)=0.
      if (in_aero(ks).lt.0.) in_aero(ks)=0.

  !ZHG 2015 Cloud liquid & ice water (CLWC+CIWC) from ECMWF
  ! nested fields
      if (ngrid.gt.0.and.readclouds_this_nest) then
        cl=ctwcn(ix,jy,n,ngrid)*(grfraction(1)/cc)
      else if (ngrid.eq.0.and.readclouds) then
        cl=ctwc(ix,jy,n)*(grfraction(1)/cc)
      else                                  !parameterize cloudwater m2/m3
  !ZHG updated parameterization of cloud water to better reproduce the values coming from ECMWF
  ! sec test
  !           cl=1E6*1E-7*prec(1)**0.3 !Sec GFS new
        cl=1E6*2E-7*prec(1)**0.36 !Sec ECMWF new, is also suitable for GFS
  !           cl=2E-7*prec(1)**0.36 !Andreas
  !           cl=1.6E-6*prec(1)**0.36 !Henrik
      endif

  !ZHG: Calculate the partition between liquid and water phase water. 
      if (act_temp .le. 253.) then
        liq_frac=0
        ice_frac=1
      else if (act_temp .ge. 273.) then
        liq_frac=1
        ice_frac=0
      else
  ! sec bugfix after FLEXPART paper review, liq_frac was 1-liq_frac
  ! IP bugfix v10.4, calculate ice_frac and liq_frac
        ice_frac= ((act_temp-273.)/(273.-253.))**2.
        !liq_frac = 1-ice_frac   !((act_temp-253.)/(273.-253.))**2.
        liq_frac=max(0.,1.-ice_frac)
      end if
  ! ZHG: Calculate the aerosol partition based on cloud phase and Ai and Bi
  !         frac_act = liq_frac*ccn_aero(ks) +(1-liq_frac)*in_aero(ks)
  ! IP, use ice_frac and liq_frac 
      frac_act = liq_frac*ccn_aero(ks) + ice_frac*in_aero(ks)

  !ZHG Use the activated fraction and the liqid water to calculate the washout

  ! AEROSOL
  !********
      if (dquer(ks).gt.0.) then
        S_i= frac_act/cl
  ! GAS
  !****
      else
        cle=(1-cl)/(henry(ks)*(r_air/3500.)*act_temp)+cl
        S_i=1/cle
      endif ! gas or particle

  ! scavenging coefficient based on Hertel et al 1995 - using the S_i for either gas or aerosol
  !SEC wetscav fix, the cloud height is no longer needed, it gives wrong results
      wetscav=incloud_ratio*S_i*(prec(1)/3.6E6)
    endif ! positive in-cloud scavenging parameters given in Species file
  endif !incloud
end subroutine get_wetscav

subroutine wetdepokernel(nunc,deposit,x,y,nage,kp,thread)
  !                          i      i    i i  i
  !*****************************************************************************
  !                                                                            *
  !     Attribution of the deposition from an individual particle to the       *
  !     deposition fields using a uniform kernel with bandwidths dxout and dyout.*
  !                                                                            *
  !     Author: A. Stohl                                                       *
  !                                                                            *
  !     26 December 1996                                                       *
  !                                                                            *
  !*****************************************************************************
  !                                                                            *
  ! Variables:                                                                 *
  !                                                                            *
  ! nunc             uncertainty class of the respective particle              *
  ! nage             age class of the respective particle                      *
  ! deposit          amount (kg) to be deposited                               *
  !                                                                            *
  !*****************************************************************************
  ! Changes:
  ! eso 10/2016: Added option to disregard kernel 
  ! 
  !*****************************************************************************

  use unc_mod

  implicit none

  integer,intent(in) :: thread
  real :: x,y,deposit(maxspec),ddx,ddy,xl,yl,wx,wy,w
  integer :: ix,jy,ixp,jyp,nunc,nage,ks,kp

  xl=(x*dx+xoutshift)/dxout
  yl=(y*dy+youtshift)/dyout
  ix=int(xl)
  jy=int(yl)
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

  ! If no kernel is used, direct attribution to grid cell
  !******************************************************

  if (.not.lusekerneloutput) then
    do ks=1,nspec
      if ((ix.ge.0).and.(jy.ge.0).and.(ix.le.numxgrid-1).and. &
           (jy.le.numygrid-1)) then
#ifdef _OPENMP
        wetgridunc_omp(ix,jy,ks,kp,nunc,nage,thread)= &
             wetgridunc_omp(ix,jy,ks,kp,nunc,nage,thread)+deposit(ks)
#else
        wetgridunc(ix,jy,ks,kp,nunc,nage)= &
             wetgridunc(ix,jy,ks,kp,nunc,nage)+deposit(ks)
#endif
      end if
    end do
  else ! use kernel 
    
  ! Determine mass fractions for four grid points
  !**********************************************

  do ks=1,nspec

    if ((ix.ge.0).and.(jy.ge.0).and.(ix.le.numxgrid-1).and. &
       (jy.le.numygrid-1)) then
      w=wx*wy
#ifdef _OPENMP
      wetgridunc_omp(ix,jy,ks,kp,nunc,nage,thread)= &
           wetgridunc_omp(ix,jy,ks,kp,nunc,nage,thread)+deposit(ks)*w
#else
      wetgridunc(ix,jy,ks,kp,nunc,nage)= &
           wetgridunc(ix,jy,ks,kp,nunc,nage)+deposit(ks)*w
#endif
    endif

    if ((ixp.ge.0).and.(jyp.ge.0).and.(ixp.le.numxgrid-1).and. &
       (jyp.le.numygrid-1)) then
      w=(1.-wx)*(1.-wy)
#ifdef _OPENMP
      wetgridunc_omp(ixp,jyp,ks,kp,nunc,nage,thread)= &
           wetgridunc_omp(ixp,jyp,ks,kp,nunc,nage,thread)+deposit(ks)*w
#else
      wetgridunc(ixp,jyp,ks,kp,nunc,nage)= &
           wetgridunc(ixp,jyp,ks,kp,nunc,nage)+deposit(ks)*w
#endif
    endif

    if ((ixp.ge.0).and.(jy.ge.0).and.(ixp.le.numxgrid-1).and. &
       (jy.le.numygrid-1)) then
      w=(1.-wx)*wy
#ifdef _OPENMP
      wetgridunc_omp(ixp,jy,ks,kp,nunc,nage,thread)= &
           wetgridunc_omp(ixp,jy,ks,kp,nunc,nage,thread)+deposit(ks)*w
#else
      wetgridunc(ixp,jy,ks,kp,nunc,nage)= &
           wetgridunc(ixp,jy,ks,kp,nunc,nage)+deposit(ks)*w
#endif
    endif

    if ((ix.ge.0).and.(jyp.ge.0).and.(ix.le.numxgrid-1).and. &
       (jyp.le.numygrid-1)) then
      w=wx*(1.-wy)
#ifdef _OPENMP
      wetgridunc_omp(ix,jyp,ks,kp,nunc,nage,thread)= &
           wetgridunc_omp(ix,jyp,ks,kp,nunc,nage,thread)+deposit(ks)*w
#else
      wetgridunc(ix,jyp,ks,kp,nunc,nage)= &
           wetgridunc(ix,jyp,ks,kp,nunc,nage)+deposit(ks)*w
#endif
    endif

  end do
  end if
end subroutine wetdepokernel

subroutine wetdepokernel_nest(nunc,deposit,x,y,nage,kp,thread)
  !                           i    i       i i i    i
  !*****************************************************************************
  !                                                                            *
  !     Attribution of the deposition from an individual particle to the       *
  !     nested deposition fields using a uniform kernel with bandwidths        *
  !     dxoutn and dyoutn.                                                     *
  !                                                                            *
  !     Author: A. Stohl                                                       *
  !                                                                            *
  !     26 December 1996                                                       *
  !                                                                            *
  !      2 September 2004: Adaptation from wetdepokernel.                      *
  !                                                                            *
  !                                                                            *
  !*****************************************************************************
  !                                                                            *
  ! Variables:                                                                 *
  !                                                                            *
  ! nunc             uncertainty class of the respective particle              *
  ! nage             age class of the respective particle                      *
  ! deposit          amount (kg) to be deposited                               *
  !                                                                            *
  !*****************************************************************************

  use unc_mod

  implicit none

  integer,intent(in) :: thread
  real :: x,y,deposit(maxspec),ddx,ddy,xl,yl,wx,wy,w
  integer :: ix,jy,ixp,jyp,ks,kp,nunc,nage

  xl=(x*dx+xoutshiftn)/dxoutn
  yl=(y*dy+youtshiftn)/dyoutn

  ! old:
  ! ix=int(xl) 
  ! jy=int(yl)

  ! ESO: for xl,yl in range <-.5,-1> we get ix,jy=0 and thus negative
  ! wx,wy as the int function rounds upwards for negative numbers.
  ! Either use the floor function, or (perhaps more correctly?) use "(xl.gt.-0.5)" 
  ! in place of "(ix.ge.0)" and similar for the upper boundary.

  ! new:
  ix=floor(xl)
  jy=floor(yl)

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

  do ks=1,nspec
    if ((ix.ge.0).and.(jy.ge.0).and.(ix.le.numxgridn-1).and. &
         (jy.le.numygridn-1)) then
      w=wx*wy
#ifdef _OPENMP
      wetgriduncn_omp(ix,jy,ks,kp,nunc,nage,thread)= &
           wetgriduncn_omp(ix,jy,ks,kp,nunc,nage,thread)+deposit(ks)*w
#else
      wetgriduncn(ix,jy,ks,kp,nunc,nage)= &
           wetgriduncn(ix,jy,ks,kp,nunc,nage)+deposit(ks)*w
#endif
    endif

    if ((ixp.ge.0).and.(jyp.ge.0).and.(ixp.le.numxgridn-1).and. &
         (jyp.le.numygridn-1)) then
      w=(1.-wx)*(1.-wy)
#ifdef _OPENMP
      wetgriduncn_omp(ixp,jyp,ks,kp,nunc,nage,thread)= &
           wetgriduncn_omp(ixp,jyp,ks,kp,nunc,nage,thread)+deposit(ks)*w
#else
      wetgriduncn(ixp,jyp,ks,kp,nunc,nage)= &
           wetgriduncn(ixp,jyp,ks,kp,nunc,nage)+deposit(ks)*w
#endif
    endif

    if ((ixp.ge.0).and.(jy.ge.0).and.(ixp.le.numxgridn-1).and. &
         (jy.le.numygridn-1)) then
      w=(1.-wx)*wy
#ifdef _OPENMP
      wetgriduncn_omp(ixp,jy,ks,kp,nunc,nage,thread)= &
           wetgriduncn_omp(ixp,jy,ks,kp,nunc,nage,thread)+deposit(ks)*w
#else
      wetgriduncn(ixp,jy,ks,kp,nunc,nage)= &
           wetgriduncn(ixp,jy,ks,kp,nunc,nage)+deposit(ks)*w
#endif
    endif

    if ((ix.ge.0).and.(jyp.ge.0).and.(ix.le.numxgridn-1).and. &
         (jyp.le.numygridn-1)) then
      w=wx*(1.-wy)
#ifdef _OPENMP
      wetgriduncn_omp(ix,jyp,ks,kp,nunc,nage,thread)= &
           wetgriduncn_omp(ix,jyp,ks,kp,nunc,nage,thread)+deposit(ks)*w
#else
      wetgriduncn(ix,jyp,ks,kp,nunc,nage)= &
           wetgriduncn(ix,jyp,ks,kp,nunc,nage)+deposit(ks)*w
#endif
    endif

  end do
end subroutine wetdepokernel_nest

subroutine writeprecip(itime,imem)

  !*****************************************************************************
  !                                                                            *
  !  This routine produces a file containing total precipitation for each      *
  !  releases point.                                                           *
  !                                                                            *
  !     Author: S. Eckhardt                                                    * 
  !     7 Mai 2017                                                             *
  !*****************************************************************************

  use point_mod
  use par_mod
  use com_mod
  use date_mod
  use windfields_mod

  implicit none

  integer :: jjjjmmdd,ihmmss,itime,i
  real(kind=dp) :: jul
  character :: adate*8,atime*6

  integer :: ix,jy,imem
  real :: xp1,yp1

  
  if (itime.eq.0) then
      open(unitprecip,file=path(2)(1:length(2))//'wetscav_precip.txt', &
       form='formatted',err=998)
  else
      open(unitprecip,file=path(2)(1:length(2))//'wetscav_precip.txt', &
       ACCESS='APPEND',form='formatted',err=998)
  endif

  jul=bdate+real(itime,kind=dp)/86400._dp
  call caldate(jul,jjjjmmdd,ihmmss)
  write(adate,'(i8.8)') jjjjmmdd
  write(atime,'(i6.6)') ihmmss

  do i=1,numpoint
    xp1=xpoint1(i)*dx+xlon0 !lat, long (real) coord
    yp1=ypoint1(i)*dy+ylat0 !lat, long (real) coord
    ix=int((xpoint1(i)+xpoint2(i))/2.)
    jy=int((ypoint1(i)+ypoint2(i))/2.)
    write(unitprecip,*)  jjjjmmdd, ihmmss, & 
           xp1,yp1,lsprec(ix,jy,1,imem),convprec(ix,jy,1,imem) !time is the same as in the ECMWF windfield
! units mm/h, valid for the time given in the windfield
  end do

  close(unitprecip)

  return


998   write(*,*) ' #### FLEXPART MODEL ERROR!   THE FILE         #### '
  write(*,*) ' #### '//path(2)(1:length(2))//'header_txt'//' #### '
  write(*,*) ' #### CANNOT BE OPENED. IF A FILE WITH THIS    #### '
  write(*,*) ' #### NAME ALREADY EXISTS, DELETE IT AND START #### '
  write(*,*) ' #### THE PROGRAM AGAIN.                       #### '
  stop
end subroutine writeprecip

end module wetdepo_mod