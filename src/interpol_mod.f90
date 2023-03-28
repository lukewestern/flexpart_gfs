  !*****************************************************************************
  !                                                                            *
  ! L. Bakels 2022: This module contains all interpolation subroutines         *
  !                 Code has been organised into subroutines                   *
  !                 Vertical logarithmic interpolation is optional (par_mod)   *
  !                                                                            *
  !*****************************************************************************

module interpol_mod
  use par_mod
  use com_mod
  use windfields_mod
  use particle_mod

  implicit none

  real,dimension(nzmax) ::          &
    uprof,vprof,wprof,wprofeta,             &
    usigprof,vsigprof,wsigprof,wsigprofeta, &
    rhoprof,rhogradprof
  logical,dimension(nzmax) ::       &
    indzindicator

  real :: u,v,w,usig,vsig,wsig,ueta,veta,weta,wsigeta

  real :: p1,p2,p3,p4,ddx,ddy,rddx,rddy,dtt,dt1,dt2
  real :: xtn,ytn
  real :: dz1out,dz2out
  integer :: nix,njy
  integer :: ix,jy,ixp,jyp,ngrid,indz,indzp,indzeta,indzpeta
  integer :: induv,indpuv
  logical :: depoindicator(maxspec)
  logical :: lbounds(2),lbounds_w(2),lbounds_uv(2) ! marking particles below or above bounds

  private :: interpol_wind_meter,interpol_wind_eta
  private :: standard_deviation_meter,standard_deviation_eta
  private :: interpol_partoutput_value_eta,interpol_partoutput_value_meter

  interface horizontal_interpolation
    procedure horizontal_interpolation_4d,horizontal_interpolation_2d
  end interface horizontal_interpolation
  
  interface horizontal_interpolation_nests
    procedure horizontal_interpolation_4d_nests,horizontal_interpolation_2d_nests
  end interface horizontal_interpolation_nests


  interface find_ngrid
    procedure find_ngrid_dp, find_ngrid_float
  end interface find_ngrid
!$OMP THREADPRIVATE(uprof,vprof,wprof,usigprof,vsigprof,wsigprof, &
!$OMP rhoprof,rhogradprof,u,v,w,usig,vsig,wsig, &
!$OMP p1,p2,p3,p4,ddx,ddy,rddx,rddy,dtt,dt1,dt2,ix,jy,ixp,jyp, &
!$OMP ngrid,indz,indzp,depoindicator,indzindicator, &
!$OMP wprofeta,wsigprofeta,induv,indpuv,lbounds,lbounds_w,lbounds_uv, &
!$OMP indzeta,indzpeta,ueta,veta,weta,wsigeta, &
!$OMP xtn,ytn,nix,njy,dz1out,dz2out)

contains

subroutine interpol_allocate
  ! allocate(uprof(nzmax),vprof(nzmax),wprof(nzmax),wprofeta(nzmax),      &
  !   usigprof(nzmax),vsigprof(nzmax),wsigprof(nzmax),wsigprofeta(nzmax), &
  !   rhoprof(nzmax),rhogradprof(nzmax),indzindicator(nzmax))
end subroutine interpol_allocate

subroutine interpol_deallocate
  ! deallocate(uprof,vprof,wprof,wprofeta,      &
  !   usigprof,vsigprof,wsigprof,wsigprofeta, &
  !   rhoprof,rhogradprof,indzindicator)
end subroutine interpol_deallocate

subroutine initialise_interpol_mod(itime,xt,yt,zt,zteta)
  ! This routine initialises all important values used in the interpol module
  ! This includes:
  ! - The current grid number in which the particle is positioned
  ! - The interpolation fractions of the grid (x,y,z) and of time

  implicit none

  integer, intent(in) :: itime             ! time step
  real, intent(in)    :: xt,yt             ! particle positions
  real, intent(in)    :: zt                ! height in meters
  real, intent(in)    :: zteta             ! height in eta coordinates

  call find_ngrid(xt,yt)
  call determine_grid_coordinates(xt,yt)
  call find_grid_distances(xt,yt)
  call find_time_variables(itime)
  call find_z_level(zt,zteta)
end subroutine initialise_interpol_mod

subroutine determine_grid_coordinates(xt,yt)
  implicit none 

  real, intent(in) :: xt,yt                 ! particle positions

  if (ngrid.gt.0) then
    xtn=(xt-xln(ngrid))*xresoln(ngrid)
    ytn=(yt-yln(ngrid))*yresoln(ngrid)
    ix=int(xtn)
    jy=int(ytn)
    nix=nint(xtn)
    njy=nint(ytn)
    if (ix+1 .ge. nxmaxn-1) then
      ixp=ix
    else
      ixp=ix+1
    endif
    if (jy+1 .ge. nymaxn-1) then
      jyp=jy
    else
      jyp=jy+1
    endif
    return
  else
    ix=int(xt)
    jy=int(yt)
    nix=nint(xt)
    njy=nint(yt)
    ixp=ix+1
    jyp=jy+1
  endif

  ! eso: Temporary fix for particle exactly at north pole
  if (jyp.ge.nymax) then
    write(*,*) 'WARNING: interpol_mod.f90 jyp >= nymax. xt,yt:',xt,yt
    jyp=jyp-1
  end if

  if (ixp.ge.nxmax) then
    write(*,*) 'WARNING: interpol_mod.f90 ixp >= nxmax. xt,yt:',xt,yt
    ixp=ixp-nxmax
  end if
end subroutine determine_grid_coordinates

subroutine find_grid_distances(xt,yt)

  implicit none 

  real, intent(in) :: xt,yt                 ! particle positions

  if (ngrid.le.0) then
    ddx=xt-real(ix)
    ddy=yt-real(jy)
  else
    ddx=xtn-real(ix)
    ddy=ytn-real(jy)
  endif
  rddx=1.-ddx
  rddy=1.-ddy
  p1=rddx*rddy
  p2=ddx*rddy
  p3=rddx*ddy
  p4=ddx*ddy
end subroutine find_grid_distances

subroutine find_time_variables(itime)
  
  implicit none  

  integer, intent(in) :: itime             ! time step
  
  dt1=real(itime-memtime(1))
  dt2=real(memtime(2)-itime)
  dtt=1./(dt1+dt2)  
end subroutine find_time_variables

subroutine find_z_level(zt,zteta)
  implicit none 
  real, intent(in)     :: &
    zt,                   & ! height in meters
    zteta                   ! height in eta

  select case (wind_coord_type)
    case('ETA')
      call find_z_level_meters(zt)
      call find_z_level_eta(zteta)
    case('METER')
      call find_z_level_meters(zt)
    case default
      call find_z_level_meters(zt)
  end select
end subroutine find_z_level

subroutine find_z_level_meters(zt)
  implicit none
  real, intent(in)     :: zt       ! height in meters
  integer              :: i

  indz=nz-1
  indzp=nz
  if (zt.le.height(1)) then
    lbounds(1)=.true.
    lbounds(2)=.false.
    indz=1
    indzp=2
  else if (zt.ge.height(nz)) then
    lbounds(1)=.false.
    lbounds(2)=.true.
  else
    lbounds(1)=.false.
    lbounds(2)=.false.
    do i=2,nz
      if (height(i).gt.zt) then
        indz=i-1
        indzp=i
        exit
      endif
    end do
  endif
end subroutine find_z_level_meters

subroutine find_z_level_eta(zteta)
  implicit none
  real, intent(in)       :: zteta    ! height in eta coordinates
  integer                :: i        ! loop variable

  call find_z_level_eta_w(zteta)

  call find_z_level_eta_uv(zteta)
end subroutine find_z_level_eta

subroutine find_z_level_eta_w(zteta)
  implicit none
  real, intent(in)       :: zteta    ! height in eta coordinates
  integer                :: i        ! loop variable

  indzeta=nz-1
  indzpeta=nz
  ! Flag particles that are above or below bounds
  if (zteta.ge.wheight(1)) then
    lbounds_w(1)=.true.
    lbounds_w(2)=.false.
    indzeta=1
    indzpeta=2
  else if (zteta.le.wheight(nz)) then
    lbounds_w(1)=.false.
    lbounds_w(2)=.true.
  else
    lbounds_w(1)=.false.
    lbounds_w(2)=.false.
    do i=2,nz
      if (wheight(i).lt.zteta) then
        indzeta=i-1
        indzpeta=i
        exit
      endif
    end do
  endif
end subroutine find_z_level_eta_w

subroutine find_z_level_eta_uv(zteta)
  implicit none
  real, intent(in)       :: zteta    ! height in eta coordinates
  integer                :: i        ! loop variable

  induv=nz-1
  indpuv=nz
  if (zteta.gt.uvheight(1)) then
    lbounds_uv(1)=.true.
    lbounds_uv(2)=.false.
    induv=1 
    indpuv=2
  else if (zteta.lt.uvheight(nz)) then
    lbounds_uv(1)=.false.
    lbounds_uv(2)=.true.
  else
    lbounds_uv(1)=.false.
    lbounds_uv(2)=.false.
    do i=2,nz
      if (uvheight(i).lt.zteta) then
        induv=i-1
        indpuv=i
        exit
      endif
    end do
  endif
end subroutine find_z_level_eta_uv

subroutine find_vertical_variables(vertlevels,zpos,zlevel,dz1,dz2,bounds,wlevel)
  !*****************************************************************************
  !                                                                            *
  ! This subroutine computes the vertical interpolation variables              *
  ! logarithmically, unless logarithmic_interpolation=.false. in the par_mod   *
  !                                                                            *
  ! Author: L. Bakels                                                          *
  !*****************************************************************************

  implicit none
  real, intent(in)    :: vertlevels(:)     ! vertical levels in coordinate system
  real, intent(in)    :: zpos              ! verticle particle position
  integer, intent(in) :: zlevel            ! vertical level of interest
  logical, intent(in) :: bounds(2),wlevel         ! flag marking if particles are outside bounds  
  real, intent(inout) :: dz1,dz2           ! fractional distance to point 1 (closer to ground) and 2
  real                :: dz,dh1,dh,pfact
  real                :: psint1(2),psint,pr1,pr2,pr_test       ! pressure of encompassing levels
  integer             :: m

  ! Only do logarithmic interpolation when using ETA coordinates, since the
  ! levels are following pressure, while METER levels are linear.
  !##############################################################
  if (.not. logarithmic_interpolation) then
    call find_vertical_variables_lin(vertlevels,zpos,zlevel,dz1,dz2,bounds,wlevel)
    return
  endif
  
  ! To check if taking the logarithm is safe
  if (wlevel) then
    pr_test=akm(zlevel+1)+bkm(zlevel+1)
  else
    pr_test=akz(zlevel+1)+bkz(zlevel+1)
  endif

  ! If the particle is below bounds (bounds(1)==.true.):
  if (bounds(1)) then
    dz1=0.
    dz2=1.
  ! If above bounds (bounds(2)==.true.):
  else if (bounds(2)) then
    dz1=1.
    dz2=0.

  ! Instead of the linear z variables, we need the ones that correspond to 
  ! the pressure of the height of the particle in relation to the model levels
  !***************************************************************************
  else if (pr_test.eq.0) then
    dz=1./(vertlevels(zlevel+1)-vertlevels(zlevel))
    dz1=(zpos-vertlevels(zlevel))*dz
    dz2=(vertlevels(zlevel+1)-zpos)*dz
  else
    if (ngrid.le.0) then
      do m=1,2
        call horizontal_interpolation(ps,psint1(m),1,memind(m),1)
      end do
    else
      do m=1,2
        call horizontal_interpolation_nests(psn,psint1(m),1,memind(m),1)
      end do
    endif
    call temporal_interpolation(psint1(1),psint1(2),psint)
    dh = vertlevels(zlevel+1)-vertlevels(zlevel)
    dh1 = zpos - vertlevels(zlevel)
    if (wlevel) then
      pr1=akm(zlevel) + bkm(zlevel)*psint
      pr2=akm(zlevel+1) + bkm(zlevel+1)*psint
    else
      pr1=akz(zlevel) + bkz(zlevel)*psint
      pr2=akz(zlevel+1) + bkz(zlevel+1)*psint
    endif
    pfact = log(pr2/pr1)*dh1/dh  
    dz = 1./(pr2-pr1)
    dz1 = pr1*(exp(pfact)-1.)*dz 
    dz2 = 1.-dz1
  endif
  ! else if ((vertlevels(zlevel).eq.0).or.(vertlevels(zlevel+1).eq.0)) then
  !   ! Linear interpolation for bottom or top layer is zero
  !   dz=1./(vertlevels(zlevel+1)-vertlevels(zlevel))
  !   dz1=(zpos-vertlevels(zlevel))*dz
  !   dz2=(vertlevels(zlevel+1)-zpos)*dz
  ! else 
  !   ! Logaritmic interpolation
  !   dz=1./(log(vertlevels(zlevel+1))-log(vertlevels(zlevel)))
  !   dz1=(log(zpos)-log(vertlevels(zlevel)))*dz
  !   dz2=(log(vertlevels(zlevel+1))-log(zpos))*dz
  ! endif
end subroutine find_vertical_variables

subroutine find_vertical_variables_lin(vertlevels,zpos,zlevel,dz1,dz2,bounds,wlevel)
  implicit none
  real, intent(in)    :: vertlevels(:)     ! vertical levels in coordinate system
  real, intent(in)    :: zpos              ! verticle particle position
  integer, intent(in) :: zlevel            ! vertical level of interest
  logical, intent(in) :: bounds(2),wlevel         ! flag marking if particles are outside bounds  
  real, intent(inout) :: dz1,dz2           ! fractional distance to point 1 (closer to ground) and 2
  real                :: dz,dh1,dh,pfact
  real                :: psint1(2),psint,pr1,pr2,temp       ! pressure of encompassing levels

  ! If the particle is below bounds (bounds(1)==.true.):
  if (bounds(1)) then
    dz1=0.
    dz2=1.
  ! If above bounds (bounds(2)==.true.):
  else if (bounds(2)) then
    dz1=1.
    dz2=0.
  else
    dz=1./(vertlevels(zlevel+1)-vertlevels(zlevel))
    dz1=(zpos-vertlevels(zlevel))*dz
    dz2=(vertlevels(zlevel+1)-zpos)*dz
  endif
end subroutine find_vertical_variables_lin

subroutine find_ngrid_dp(xt,yt)

  implicit none
  real ::                      &
    eps           
  real(kind=dp), intent(in) :: &
    xt,yt                           ! particle positions on grid
  integer ::                   &
    j

  eps=nxmax/3.e5
  if (nglobal.and.(real(yt).gt.switchnorthg)) then
    ngrid=-1
  else if (sglobal.and.(real(yt).lt.switchsouthg)) then
    ngrid=-2
  else
    ngrid=0
    do j=numbnests,1,-1
      if ((real(xt).gt.xln(j)+eps).and.(real(xt).lt.xrn(j)-eps).and. &
           (real(yt).gt.yln(j)+eps).and.(real(yt).lt.yrn(j)-eps)) then
        ngrid=j
        exit
      endif
    end do
  endif
end subroutine find_ngrid_dp

subroutine find_ngrid_float(xt,yt)

  implicit none
  real ::                      &
    eps           
  real, intent(in) :: &
    xt,yt                           ! particle positions on grid
  integer ::                   &
    j

  eps=nxmax/3.e5
  if (nglobal.and.(yt.gt.switchnorthg)) then
    ngrid=-1
  else if (sglobal.and.(yt.lt.switchsouthg)) then
    ngrid=-2
  else
    ngrid=0
    do j=numbnests,1,-1
      if ((xt.gt.xln(j)+eps).and.(xt.lt.xrn(j)-eps).and. &
           (yt.gt.yln(j)+eps).and.(yt.lt.yrn(j)-eps)) then
        ngrid=j
        exit
      endif
    end do
  endif
end subroutine find_ngrid_float

subroutine horizontal_interpolation_4d(field,output,zlevel,indexh,ztot)

  implicit none

  integer, intent(in) :: zlevel,ztot,indexh                       ! interpolation z level, z
  real, intent(in)    :: field(0:nxmax-1,0:nymax-1,ztot,numwfmem) ! input field to interpolate over
  real, intent(inout) :: output                                   ! interpolated values

  output=p1*field(ix ,jy ,zlevel,indexh) &
       + p2*field(ixp,jy ,zlevel,indexh) &
       + p3*field(ix ,jyp,zlevel,indexh) &
       + p4*field(ixp,jyp,zlevel,indexh)
end subroutine horizontal_interpolation_4d

subroutine horizontal_interpolation_2d(field,output)
  implicit none 
  real, intent(in)    :: field(0:nxmax-1,0:nymax-1)       ! 2D imput field
  real, intent(inout) :: output                           ! Interpolated value

  output=p1*field(ix ,jy) &
         + p2*field(ixp,jy) &
         + p3*field(ix ,jyp) &
         + p4*field(ixp,jyp)
end subroutine horizontal_interpolation_2d

subroutine horizontal_interpolation_4d_nests(field,output,zlevel,indexh,ztot)

  implicit none

  integer, intent(in) :: zlevel,ztot,indexh                       ! interpolation z level, z
  real, intent(in)    :: field(0:nxmaxn-1,0:nymaxn-1,ztot,numwfmem,numbnests) ! input field to interpolate over
  real, intent(inout) :: output                                   ! interpolated values

  output=p1*field(ix ,jy ,zlevel,indexh,ngrid) &
       + p2*field(ixp,jy ,zlevel,indexh,ngrid) &
       + p3*field(ix ,jyp,zlevel,indexh,ngrid) &
       + p4*field(ixp,jyp,zlevel,indexh,ngrid)
end subroutine horizontal_interpolation_4d_nests

subroutine horizontal_interpolation_2d_nests(field,output)

  implicit none

  real, intent(in)    :: field(0:nxmaxn-1,0:nymaxn-1,numbnests) ! input field to interpolate over
  real, intent(inout) :: output                                   ! interpolated values

  output=p1*field(ix ,jy ,ngrid) &
       + p2*field(ixp,jy ,ngrid) &
       + p3*field(ix ,jyp,ngrid) &
       + p4*field(ixp,jyp,ngrid)
end subroutine horizontal_interpolation_2d_nests

subroutine temporal_interpolation(time1,time2,output)

  implicit none

  real, intent(in)    :: time1,time2     ! input data at two timesteps 
  real, intent(inout) :: output          ! interpolated data

  output=(time1*dt2+time2*dt1)*dtt
end subroutine temporal_interpolation

subroutine vertical_interpolation(input1,input2,dz1,dz2,output)

  implicit none

  real, intent(in)    :: input1,input2   ! input data at two vertical levels, 1 being closer to ground
  real, intent(in)    :: dz1,dz2         ! logarithmic interpolation values
  real, intent(inout) :: output          ! interpolated data

  output = input1*dz2 + input2*dz1!input1**dz2 * input2**dz1
end subroutine vertical_interpolation

subroutine bilinear_spatial_interpolation(field,output,zlevel,dz1,dz2,ztot)
  implicit none
  integer, intent(in) :: zlevel,ztot                               ! interpolation z level
  real, intent(in)    :: field(0:nxmax-1,0:nymax-1,ztot,numwfmem)  ! input field to interpolate over
  real, intent(in)    :: dz1,dz2
  real, intent(inout) :: output(2)                                 ! interpolated values
  integer             :: m,n,indzh
  real                :: output1(2)

  do m=1,2
    do n=1,2
      indzh=zlevel+n-1
      call horizontal_interpolation_4d(field,output1(n),indzh,memind(m),ztot)
    end do
    !**********************************
    ! 2.) Linear vertical interpolation on logarithmic scale
    !**********************************
    call vertical_interpolation(output1(1),output1(2),dz1,dz2,output(m))
  end do
end subroutine bilinear_spatial_interpolation

subroutine bilinear_spatial_interpolation_nests(field,output,zlevel,dz1,dz2,ztot)
  implicit none
  integer, intent(in) :: zlevel,ztot                               ! interpolation z level
  real, intent(in)    :: field(0:nxmaxn-1,0:nymaxn-1,ztot,numwfmem,numbnests)  ! input field to interpolate over
  real, intent(in)    :: dz1,dz2
  real, intent(inout) :: output(2)                                 ! interpolated values
  integer             :: m,n,indzh
  real                :: output1(2)

  do m=1,2
    do n=1,2
      indzh=zlevel+n-1
      call horizontal_interpolation_4d_nests(field,output1(n),indzh,memind(m),ztot)
    end do
    !**********************************
    ! 2.) Linear vertical interpolation on logarithmic scale
    !**********************************
    call vertical_interpolation(output1(1),output1(2),dz1,dz2,output(m))
  end do
end subroutine bilinear_spatial_interpolation_nests

subroutine compute_sl_sq(field,sl,sq,zlevel,indexh,ztot)
  implicit none

  integer, intent(in) :: zlevel,ztot,indexh                       ! interpolation z levels
  real, intent(in)    :: field(0:nxmax-1,0:nymax-1,ztot,numwfmem) ! input field to interpolate over
  real, intent(inout) :: sl,sq                                   ! standard deviation


  sl=sl+field(ix ,jy ,zlevel,indexh)+field(ixp,jy ,zlevel,indexh) &
       +field(ix ,jyp,zlevel,indexh)+field(ixp,jyp,zlevel,indexh)
  sq=sq+field(ix ,jy ,zlevel,indexh)*field(ix ,jy ,zlevel,indexh)+ &
        field(ixp,jy ,zlevel,indexh)*field(ixp,jy ,zlevel,indexh)+ &
        field(ix ,jyp,zlevel,indexh)*field(ix ,jyp,zlevel,indexh)+ &
        field(ixp,jyp,zlevel,indexh)*field(ixp,jyp,zlevel,indexh)
end subroutine compute_sl_sq

subroutine compute_sl_sq_nests(field,sl,sq,zlevel,indexh,ztot)
  implicit none

  integer, intent(in) :: zlevel,ztot,indexh                       ! interpolation z levels
  real, intent(in)    :: field(0:nxmaxn-1,0:nymaxn-1,ztot,numwfmem,numbnests) ! input field to interpolate over
  real, intent(inout) :: sl,sq                                   ! standard deviation


  sl=sl+field(ix ,jy ,zlevel,indexh,ngrid)+field(ixp,jy ,zlevel,indexh,ngrid) &
       +field(ix ,jyp,zlevel,indexh,ngrid)+field(ixp,jyp,zlevel,indexh,ngrid)
  sq=sq+field(ix ,jy ,zlevel,indexh,ngrid)*field(ix ,jy ,zlevel,indexh,ngrid)+ &
        field(ixp,jy ,zlevel,indexh,ngrid)*field(ixp,jy ,zlevel,indexh,ngrid)+ &
        field(ix ,jyp,zlevel,indexh,ngrid)*field(ix ,jyp,zlevel,indexh,ngrid)+ &
        field(ixp,jyp,zlevel,indexh,ngrid)*field(ixp,jyp,zlevel,indexh,ngrid)
end subroutine compute_sl_sq_nests

subroutine standard_deviation(sl,sq,ndivide,output)
  implicit none

  real, intent(in) :: sl,sq,ndivide
  real, intent(out) :: output
  real :: xaux
  real,parameter      :: eps=1.0e-30

  xaux=sq-sl*sl/ndivide

  if (xaux.lt.eps) then
    output=0.
  else
    output=sqrt(xaux/(ndivide-1.))
  endif
end subroutine standard_deviation

! Interpolation functions
!************************
subroutine interpol_PBL(itime,xt,yt,zt,zteta)
  !                          i   i  i  i
  !*****************************************************************************
  !                                                                            *
  !  This subroutine interpolates everything that is needed for calculating the*
  !  dispersion.                                                               *
  !                                                                            *
  !    Author: A. Stohl                                                        *
  !                                                                            *
  !    16 December 1997                                                        *
  !                                                                            *
  !  Revision March 2005 by AST : all output variables in common block cal-    *
  !                               culation of standard deviation done in this  *
  !                               routine rather than subroutine call in order *
  !                               to save computation time                     *
  !                                                                            *
  !*****************************************************************************
  !                                                                            *
  ! Variables:                                                                 *
  ! itime [s]          current temporal position                               *
  ! memtime(3) [s]     times of the wind fields in memory                      *
  ! xt,yt,zt           coordinates position for which wind data shall be       *
  !                    culated                                                 *
  !                                                                            *
  ! Constants:                                                                 *
  !                                                                            *
  !*****************************************************************************

  use turbulence_mod

  implicit none

  integer, intent(in) :: itime
  real, intent(in)    :: xt,yt,zt,zteta
  integer             :: m,n,indexh
  integer             :: iw(2),iweta(2)
  real                :: uh1(2),vh1(2),wh1(2),wetah1(2),rho1(2),rhograd1(2)
  real                :: dz1weta,dz2weta
  real,parameter      :: eps=1.0e-30

  ! Auxiliary variables needed for interpolation
  real :: ust1(2),wst1(2),oli1(2),oliaux
  !********************************************
  ! Multilinear interpolation in time and space
  !********************************************

  ! Where in the grid? Stereographic (ngrid<0) or nested (ngrid>0)
  !***************************************************************
  call find_ngrid(xt,yt)

  ! Determine the lower left corner and its distance to the current position
  !*************************************************************************
  call find_grid_distances(xt,yt)

  ! Calculate variables for time interpolation
  !*******************************************
  call find_time_variables(itime)

  !********************************************************
  ! 1. Interpolate u*, w* and Obukhov length for turbulence
  !********************************************************

  ! a) Bilinear horizontal interpolation
  if (ngrid.le.0) then ! No nest
    do m=1,2
      indexh=memind(m)
      call horizontal_interpolation(ustar,ust1(m),1,memind(m),1)
      call horizontal_interpolation(wstar,wst1(m),1,memind(m),1)
      call horizontal_interpolation(oli,oli1(m),1,memind(m),1)
    end do
  else ! Nest
    do m=1,2
      indexh=memind(m)
      call horizontal_interpolation_nests(ustarn,ust1(m),1,memind(m),1)
      call horizontal_interpolation_nests(wstarn,wst1(m),1,memind(m),1)
      call horizontal_interpolation_nests(olin,oli1(m),1,memind(m),1)
    end do
  endif    
  ! b) Temporal interpolation
  call temporal_interpolation(ust1(1),ust1(2),ust)
  call temporal_interpolation(wst1(1),wst1(2),wst)
  call temporal_interpolation(oli1(1),oli1(2),oliaux)

  if (oliaux.ne.0.) then
    ol=1./oliaux
  else
    ol=99999.
  endif

  ! Within the PBL, only METER coordinates are used
  ! with the exception of mesoscale turbulence,
  ! which uses wsigeta computed in interpol_mesoscale
  !**************************************************

  ! Determine the level below the current position
  !***********************************************
  call find_z_level_meters(zt)

  iw(:)=(/ indz, indzp /)

  ! w(eta) velocities are necessary for the Petterssen correction
  !**************************************************************
  if (wind_coord_type.eq.'ETA') then
    call find_z_level_eta(zteta)
    iweta(:)=(/ indzeta, indzpeta /)
  endif

  !**************************************
  ! 1.) Bilinear horizontal interpolation
  ! 2.) Temporal interpolation (linear)
  !**************************************
  
  ! Loop over 2 time steps and indz levels
  !***************************************
  if (ngrid.le.0) then ! No nest
    do n=1,2
      do m=1,2
        call horizontal_interpolation(ww,wh1(m),iw(n),memind(m),nzmax)
        if (wind_coord_type.eq.'ETA') &
          call horizontal_interpolation(wweta,wetah1(m),iweta(n),memind(m),nzmax)
        call horizontal_interpolation(rho,rho1(m),iw(n),memind(m),nzmax)
        call horizontal_interpolation(drhodz,rhograd1(m),iw(n),memind(m),nzmax)
        if (ngrid.lt.0) then
          call horizontal_interpolation(uupol,uh1(m),iw(n),memind(m),nzmax)
          call horizontal_interpolation(vvpol,vh1(m),iw(n),memind(m),nzmax)
        else
          call horizontal_interpolation(uu,uh1(m),iw(n),memind(m),nzmax)
          call horizontal_interpolation(vv,vh1(m),iw(n),memind(m),nzmax)
        endif
      end do
      call temporal_interpolation(wh1(1),wh1(2),wprof(iw(n)))
      if (wind_coord_type.eq.'ETA') &
        call temporal_interpolation(wetah1(1),wetah1(2),wprofeta(iweta(n)))
      call temporal_interpolation(uh1(1),uh1(2),uprof(iw(n)))
      call temporal_interpolation(vh1(1),vh1(2),vprof(iw(n)))
      call temporal_interpolation(rho1(1),rho1(2),rhoprof(iw(n)))
      call temporal_interpolation(rhograd1(1),rhograd1(2),rhogradprof(iw(n)))
    end do
  else ! Nest
    do n=1,2
      do m=1,2
        call horizontal_interpolation_nests(wwn,wh1(m),iw(n),memind(m),nzmax)
        if (wind_coord_type.eq.'ETA') &
          call horizontal_interpolation_nests(wwetan,wetah1(m),iweta(n),memind(m),nzmax)
        call horizontal_interpolation_nests(uun,uh1(m),iw(n),memind(m),nzmax)
        call horizontal_interpolation_nests(vvn,vh1(m),iw(n),memind(m),nzmax)
        call horizontal_interpolation_nests(rhon,rho1(m),iw(n),memind(m),nzmax)
        call horizontal_interpolation_nests(drhodzn,rhograd1(m),iw(n),memind(m),nzmax)
      end do
      call temporal_interpolation(wh1(1),wh1(2),wprof(iw(n)))
      if (wind_coord_type.eq.'ETA') &
        call temporal_interpolation(wetah1(1),wetah1(2),wprofeta(iweta(n)))
      call temporal_interpolation(uh1(1),uh1(2),uprof(iw(n)))
      call temporal_interpolation(vh1(1),vh1(2),vprof(iw(n)))
      call temporal_interpolation(rho1(1),rho1(2),rhoprof(iw(n)))
      call temporal_interpolation(rhograd1(1),rhograd1(2),rhogradprof(iw(n)))

      indzindicator(iw(n))=.false.
    end do
  endif

  ! Only necessary for the Petterssen correction
  if (wind_coord_type.eq.'ETA') then
    call find_vertical_variables(wheight,zteta,indzeta,dz1weta,dz2weta,lbounds_w,.true.)
    call vertical_interpolation(wprofeta(indzeta),wprofeta(indzpeta),dz1weta,dz2weta,weta)
  endif
end subroutine interpol_PBL

subroutine interpol_PBL_misslev()
  !                            
  !*****************************************************************************
  !                                                                            *
  !  This subroutine interpolates u,v,w, density and density gradients.        *
  !                                                                            *
  !    Author: A. Stohl                                                        *
  !                                                                            *
  !    16 December 1997                                                        *
  !    Update: 2 March 1999                                                    *
  !                                                                            *
  !  Revision March 2005 by AST : all output variables in common block cal-    *
  !                               culation of standard deviation done in this  *
  !                               routine rather than subroutine call in order *
  !                               to save computation time                     *
  !                                                                            *
  !*****************************************************************************
  !                                                                            *
  ! Variables:                                                                 *
  ! n                  level                                                   *
  !                                                                            *
  ! Constants:                                                                 *
  !                                                                            *
  !*****************************************************************************
  implicit none

  integer             :: n,iw(2)
  real                :: uh1(2),vh1(2),wh1(2),rho1(2),rhograd1(2)
  integer             :: m


  ! Within the PBL, only METER coordinates are used
  ! with the exception of mesoscale turbulence,
  ! which uses wsigeta computed in interpol_mesoscale
  !**************************************************

  !********************************************
  ! Multilinear interpolation in time and space
  !********************************************
  iw(:)=(/ indz, indzp /)
  do n=1,2
    if (indzindicator(iw(n))) then
      if (ngrid.le.0) then ! No nest
        do m=1,2
          call horizontal_interpolation(ww,wh1(m),iw(n),memind(m),nzmax)
          call horizontal_interpolation(rho,rho1(m),iw(n),memind(m),nzmax)
          call horizontal_interpolation(drhodz,rhograd1(m),iw(n),memind(m),nzmax)
          if (ngrid.lt.0) then
            call horizontal_interpolation(uupol,uh1(m),iw(n),memind(m),nzmax)
            call horizontal_interpolation(vvpol,vh1(m),iw(n),memind(m),nzmax)
          else
            call horizontal_interpolation(uu,uh1(m),iw(n),memind(m),nzmax)
            call horizontal_interpolation(vv,vh1(m),iw(n),memind(m),nzmax)
          endif
        end do
      else ! Nest
        do m=1,2
          call horizontal_interpolation_nests(wwn,wh1(m),iw(n),memind(m),nzmax)
          call horizontal_interpolation_nests(uun,uh1(m),iw(n),memind(m),nzmax)
          call horizontal_interpolation_nests(vvn,vh1(m),iw(n),memind(m),nzmax)
          call horizontal_interpolation_nests(rhon,rho1(m),iw(n),memind(m),nzmax)
          call horizontal_interpolation_nests(drhodzn,rhograd1(m),iw(n),memind(m),nzmax)
        end do
      endif
      call temporal_interpolation(wh1(1),wh1(2),wprof(iw(n)))
      call temporal_interpolation(uh1(1),uh1(2),uprof(iw(n)))
      call temporal_interpolation(vh1(1),vh1(2),vprof(iw(n)))
      call temporal_interpolation(rho1(1),rho1(2),rhoprof(iw(n)))
      call temporal_interpolation(rhograd1(1),rhograd1(2),rhogradprof(iw(n)))

      indzindicator(iw(n))=.false.
    endif
  end do
end subroutine interpol_PBL_misslev

subroutine interpol_PBL_short(zt,rhoa,rhograd)
  implicit none 
  real, intent(in)    :: zt
  real, intent(inout) :: rhoa,rhograd
  real                :: dz1,dz2

  call find_vertical_variables(height,zt,indz,dz1,dz2,lbounds,.false.)

  call vertical_interpolation(wprof(indz),wprof(indzp),dz1,dz2,w)
  call vertical_interpolation(uprof(indz),uprof(indzp),dz1,dz2,u)
  call vertical_interpolation(vprof(indz),vprof(indzp),dz1,dz2,v)
  call vertical_interpolation(rhoprof(indz),rhoprof(indzp),dz1,dz2,rhoa)
  call vertical_interpolation(rhogradprof(indz),rhogradprof(indzp),dz1,dz2,rhograd)
end subroutine interpol_PBL_short

subroutine interpol_mesoscale(itime,xt,yt,zt,zteta)
  use turbulence_mod

  implicit none

  integer, intent(in) :: itime
  real, intent(in)    :: xt,yt,zt,zteta
  integer             :: m,indexh
  integer             :: iw(2),iuv(2),iweta(2)

  ! Where in the grid? Stereographic (ngrid<0) or nested (ngrid>0)
  !***************************************************************
  call find_ngrid(xt,yt)

  call determine_grid_coordinates(xt,yt)

  ! Determine the level below the current position
  !***********************************************
  call find_z_level_meters(zt)
  iw(:)=(/ indz, indzp /)
  
  select case (wind_coord_type)
    case ('ETA')  
      call find_z_level_eta(zteta)
      iuv(:)=(/ induv, indpuv /)
      iweta(:)=(/ indzeta, indzpeta /)
      call standard_deviation_eta(iw,iuv,iweta)
    case ('METER')
      iw(:)=(/ indz, indzp /)
      call standard_deviation_meter(iw)
    case default
      write(*,*) 'ERROR: wind_coord_type is not allowed ', wind_coord_type
      write(*,*) 'Choose ETA or METER.'
      stop
  end select
end subroutine interpol_mesoscale

subroutine interpol_wind(itime,xt,yt,zt,zteta,pp)
  !                           i   i  i  i
  !*****************************************************************************
  !                                                                            *
  !  This subroutine interpolates the wind data to current trajectory position.*
  !                                                                            *
  !    Author: A. Stohl                                                        *
  !                                                                            *
  !    16 December 1997                                                        *
  !                                                                            *
  !  Revision March 2005 by AST : all output variables in common block cal-    *
  !                               culation of standard deviation done in this  *
  !                               routine rather than subroutine call in order *
  !                               to save computation time                     *
  !                                                                            *
  !*****************************************************************************
  !                                                                            *
  ! Variables:                                                                 *
  ! u,v,w              wind components                                         *
  ! itime [s]          current temporal position                               *
  ! memtime(3) [s]     times of the wind fields in memory                      *
  ! xt,yt,zt           coordinates position for which wind data shall be       *
  !                    calculated                                              *
  !                                                                            *
  ! Constants:                                                                 *
  !                                                                            *
  !*****************************************************************************


  implicit none

  integer, intent(in) :: itime,pp
  real, intent(in)    :: xt,yt,zt
  real, intent(in)    :: zteta
  integer             :: iw(2),iuv(2),iweta(2)


  ! Where in the grid? Stereographic (ngrid<0) or nested (ngrid>0)
  !***************************************************************
  call find_ngrid(xt,yt)

  call determine_grid_coordinates(xt,yt)
  ! ! Multilinear interpolation in time and space
  ! !********************************************

  ! Determine the lower left corner and its distance to the current position
  !*************************************************************************
  call find_grid_distances(xt,yt)

  ! Calculate variables for time interpolation
  !*******************************************
  call find_time_variables(itime)

  ! Interpolate over the windfields depending on the prefered
  ! coordinate system
  !**********************************************************
  select case (wind_coord_type)
    case ('ETA')
      ! Same for eta coordinates
      !*************************
      call find_z_level_eta(zteta)

      iuv(:)  = (/ induv, indpuv /)
      iweta(:)= (/ indzeta, indzpeta /)
      call interpol_wind_eta(zteta,iuv,iweta)
      !call standard_deviation_wind_eta(iw,iuv,iweta)
    case ('METER')
      ! Determine the level below the current position for u,v
      !*******************************************************
      call find_z_level_meters(zt)

      iw(:)=(/ indz, indzp /)
      call interpol_wind_meter(zt,iw)
      !call standard_deviation_wind_meter(iw)

    case default
      write(*,*) 'ERROR: wind_coord_type is not allowed ', wind_coord_type
      write(*,*) 'Choose ETA or METER.'
      stop
  end select
end subroutine interpol_wind

subroutine interpol_wind_short(itime,xt,yt,zt,zteta)
  !                                 i   i  i  i
  !*****************************************************************************
  !                                                                            *
  !  This subroutine interpolates the wind data to current trajectory position.*
  !                                                                            *
  !    Author: A. Stohl                                                        *
  !                                                                            *
  !    16 December 1997                                                        *
  !                                                                            *
  !  Revision March 2005 by AST : all output variables in common block         *
  !                                                                            *
  !*****************************************************************************
  !                                                                            *
  ! Variables:                                                                 *
  ! u,v,w              wind components                                         *
  ! itime [s]          current temporal position                               *
  ! memtime(3) [s]     times of the wind fields in memory                      *
  ! xt,yt,zt           coordinates position for which wind data shall be       *
  !                    calculated                                              *
  !                                                                            *
  ! Constants:                                                                 *
  !                                                                            *
  !*****************************************************************************


  implicit none

  integer, intent(in) :: itime
  real, intent(in) :: xt,yt,zt
  real, intent(in) :: zteta
  integer             :: iw(2),iuv(2),iweta(2)

  !********************************************
  ! Multilinear interpolation in time and space
  !********************************************

  ! Where in the grid? Stereographic (ngrid<0) or nested (ngrid>0)
  !***************************************************************
  call find_ngrid(xt,yt)
  call determine_grid_coordinates(xt,yt)
  call find_grid_distances(xt,yt)

  ! Calculate variables for time interpolation
  !*******************************************
  call find_time_variables(itime)

  ! Interpolate over the windfields depending on the prefered
  ! coordinate system
  !**********************************************************
  select case (wind_coord_type)
    case ('ETA')
      ! Determine the level below the current position for eta coordinates
      !*******************************************************************
      call find_z_level_eta(zteta)

      iuv(:)=(/ induv, indpuv /)
      iweta(:)=(/ indzeta, indzpeta /)
      ! Interpolate the u, v, weta windfields
      !**************************************
      call interpol_wind_eta(zteta,iuv,iweta)
    case ('METER')

      ! Determine the level below the current position for u,v
      !*******************************************************
      call find_z_level_meters(zt)

      iw(:)=(/ indz, indzp /)
      call interpol_wind_meter(zt,iw)
    case default
      write(*,*) 'ERROR: wind_coord_type is not allowed ', wind_coord_type
      write(*,*) 'Choose ETA or METER.'
      stop
  end select
end subroutine interpol_wind_short

subroutine interpol_partoutput_value(fieldname,output,j)
  implicit none
  integer, intent(in)         :: j          ! particle number
  character(2), intent(in)    :: fieldname  ! input field to interpolate over
  real, intent(inout)         :: output
  ! Interpolate over the windfields depending on the prefered
  ! coordinate system
  !**********************************************************
  select case (wind_coord_type)
    case ('ETA')
      call interpol_partoutput_value_eta(fieldname,output,j)
    case ('METER')
      call interpol_partoutput_value_meter(fieldname,output,j)
    case default
      call interpol_partoutput_value_meter(fieldname,output,j)
  end select
end subroutine interpol_partoutput_value

subroutine interpol_htropo_hmix(tropop,h)
  implicit none 
  real, intent(inout) :: &
    tropop,              &  ! height of troposphere
    h                       ! mixing height
  real                :: &
    h1(2)                   ! mixing height of 2 timesteps
  integer             :: &
    mind,                &  ! windfield index
    i,j,k,m                 ! loop variables

  h=0.
  if (ngrid.le.0) then
    if (interpolhmix) then
      do m=1,2
        call horizontal_interpolation(hmix,h1(m),1,memind(m),1)
      end do
    else
      do k=1,2
        mind=memind(k) ! eso: compatibility with 3-field version
        do j=jy,jyp
          do i=ix,ixp
             if (hmix(i,j,1,mind).gt.h) h=hmix(i,j,1,mind)
          end do
        end do
      end do
    endif
    tropop=tropopause(nix,njy,1,memind(1))
  else
    do k=1,2
      mind=memind(k)
      do j=jy,jyp
        do i=ix,ixp
          if (hmixn(i,j,1,mind,ngrid).gt.h) h=hmixn(i,j,1,mind,ngrid)
        end do
      end do
    end do
    tropop=tropopausen(nix,njy,1,memind(1),ngrid)
  endif

  if (interpolhmix) h=(h1(1)*dt2+h1(2)*dt1)*dtt 
end subroutine interpol_htropo_hmix

subroutine interpol_density(itime,ipart,output)

  implicit none

  integer, intent(in) :: itime,ipart  ! time and particle index
  real, intent(inout) :: output ! output density (rhoi)
  integer :: ind
  real :: dz1,dz2
  real :: rhoprof(2)

  ! Where in the grid? Stereographic (ngrid<0) or nested (ngrid>0)
  !***************************************************************
  call find_ngrid(part(ipart)%xlon,part(ipart)%ylat)
  call determine_grid_coordinates(real(part(ipart)%xlon),real(part(ipart)%ylat))
  call find_grid_distances(real(part(ipart)%xlon),real(part(ipart)%ylat))
  call find_time_variables(itime)

  ! Take density from 2nd wind field in memory (accurate enough, no time interpolation needed)
  !*****************************************************************************
  select case (wind_coord_type)
    case ('ETA')
      call find_z_level_eta(real(part(ipart)%zeta))
      call find_vertical_variables(uvheight,real(part(ipart)%zeta),induv,dz1,dz2,lbounds_uv,.false.)
      if (ngrid.le.0) then
        do ind=induv,indpuv
          call horizontal_interpolation(rhoeta,rhoprof(ind-induv+1),ind,memind(2),nzmax)
        end do
      else
        do ind=induv,indpuv
          call horizontal_interpolation_nests(rhoetan,rhoprof(ind-induv+1),ind,memind(2),nzmax)
        end do
      endif
    case ('METER')
      call find_z_level_meters(real(part(ipart)%z))
      call find_vertical_variables(height,real(part(ipart)%z),indz,dz1,dz2,lbounds,.false.)
      if (ngrid.le.0) then
        do ind=indz,indzp
          call horizontal_interpolation(rho,rhoprof(ind-indz+1),ind,memind(2),nzmax)
        end do
      else
        do ind=indz,indzp
          call horizontal_interpolation_nests(rhon,rhoprof(ind-indz+1),ind,memind(2),nzmax)
        end do
      endif
    case default
      stop 'wind_coord_type not defined in conccalc.f90'
  end select
  call vertical_interpolation(rhoprof(1),rhoprof(2),dz1,dz2,output)
end subroutine interpol_density

!*********************
!* PRIVATE FUNCTIONS *
!*********************
! Interpolation of wind fields
!*****************************
subroutine interpol_wind_eta(zteta,iuv,iweta)
  implicit none

  real, intent(in)    :: zteta
  integer,intent(in)  :: iuv(2),iweta(2)
  integer             :: n,m
  real                :: uh(2),vh(2),wetah(2),uh1(2),vh1(2),wetah1(2)
  real                :: dz1uv,dz2uv,dz1weta,dz2weta
  !**********************************************************************
  ! 1.) Bilinear horizontal interpolation
  ! This has to be done separately for 6 fields (Temporal(2)*Vertical(3))
  !**********************************************************************

  ! Vertical distance to the level below and above current position
  !****************************************************************
  call find_vertical_variables(uvheight,zteta,induv,dz1uv,dz2uv,lbounds_uv,.false.)
  call find_vertical_variables(wheight,zteta,indzeta,dz1weta,dz2weta,lbounds_w,.true.)

  ! Loop over 2 time steps and 2 levels
  !************************************
  if (ngrid.le.0) then ! No nest
    do m=1,2
      do n=1,2
        call horizontal_interpolation(wweta,wetah1(n),iweta(n),memind(m),nzmax)
        if (ngrid.lt.0) then
          call horizontal_interpolation(uupoleta,uh1(n),iuv(n),memind(m),nzmax)
          call horizontal_interpolation(vvpoleta,vh1(n),iuv(n),memind(m),nzmax)
        else
          call horizontal_interpolation(uueta,uh1(n),iuv(n),memind(m),nzmax)
          call horizontal_interpolation(vveta,vh1(n),iuv(n),memind(m),nzmax)
        endif
      end do
      call vertical_interpolation(uh1(1),uh1(2),dz1uv,dz2uv,uh(m))
      call vertical_interpolation(vh1(1),vh1(2),dz1uv,dz2uv,vh(m))
      call vertical_interpolation(wetah1(1),wetah1(2),dz1weta,dz2weta,wetah(m))
    end do 
  else ! Nest
    do m=1,2
      do n=1,2
        
        ! wetah1(n) = p1*wwetan(ix ,jy ,iweta(n),memind(m),ngrid) &
        !           + p2*wwetan(ixp,jy ,iweta(n),memind(m),ngrid) &
        !           + p3*wwetan(ix ,jyp,iweta(n),memind(m),ngrid) &
        !           + p4*wwetan(ixp,jyp,iweta(n),memind(m),ngrid)
        call horizontal_interpolation_nests(wwetan,wetah1(n),iweta(n),memind(m),nzmax)
        call horizontal_interpolation_nests(uuetan,uh1(n),iuv(n),memind(m),nzmax)
        call horizontal_interpolation_nests(vvetan,vh1(n),iuv(n),memind(m),nzmax)
      end do
      call vertical_interpolation(uh1(1),uh1(2),dz1uv,dz2uv,uh(m))
      call vertical_interpolation(vh1(1),vh1(2),dz1uv,dz2uv,vh(m))
      call vertical_interpolation(wetah1(1),wetah1(2),dz1weta,dz2weta,wetah(m))
    end do    
  endif
  call temporal_interpolation(uh(1),uh(2),u)
  call temporal_interpolation(vh(1),vh(2),v)
  call temporal_interpolation(wetah(1),wetah(2),weta)
end subroutine interpol_wind_eta

subroutine interpol_wind_meter(zt,iw)
  implicit none

  real, intent(in)    :: zt
  integer,intent(in)  :: iw(2)
  integer             :: n,m
  real                :: uh(2),vh(2),wh(2),uh1(2),vh1(2),wh1(2)
  real                :: dz1w,dz2w

  !**********************************************************************
  ! 1.) Bilinear horizontal interpolation
  ! This has to be done separately for 6 fields (Temporal(2)*Vertical(3))
  !**********************************************************************

  ! Vertical distance to the level below and above current position
  !****************************************************************
  call find_vertical_variables(height,zt,indz,dz1w,dz2w,lbounds,.false.)

  ! Loop over 2 time steps and 2 levels
  !************************************
  if (ngrid.le.0) then ! No nest
    do m=1,2
      do n=1,2
        call horizontal_interpolation(ww,wh1(n),iw(n),memind(m),nzmax)
        if (ngrid.lt.0) then
          call horizontal_interpolation(uupol,uh1(n),iw(n),memind(m),nzmax)
          call horizontal_interpolation(vvpol,vh1(n),iw(n),memind(m),nzmax)
        else
          call horizontal_interpolation(uu,uh1(n),iw(n),memind(m),nzmax)
          call horizontal_interpolation(vv,vh1(n),iw(n),memind(m),nzmax)
        endif
      end do
      call vertical_interpolation(wh1(1),wh1(2),dz1w,dz2w,wh(m))
      call vertical_interpolation(uh1(1),uh1(2),dz1w,dz2w,uh(m))
      call vertical_interpolation(vh1(1),vh1(2),dz1w,dz2w,vh(m))
    end do 
  else ! Nest
    do m=1,2
      do n=1,2
        call horizontal_interpolation_nests(wwn,wh1(n),iw(n),memind(m),nzmax)
        call horizontal_interpolation_nests(uun,uh1(n),iw(n),memind(m),nzmax)
        call horizontal_interpolation_nests(vvn,vh1(n),iw(n),memind(m),nzmax)
      end do
      call vertical_interpolation(wh1(1),wh1(2),dz1w,dz2w,wh(m))
      call vertical_interpolation(uh1(1),uh1(2),dz1w,dz2w,uh(m))
      call vertical_interpolation(vh1(1),vh1(2),dz1w,dz2w,vh(m))
    end do    
  endif
  call temporal_interpolation(wh(1),wh(2),w)
  call temporal_interpolation(uh(1),uh(2),u)
  call temporal_interpolation(vh(1),vh(2),v)
end subroutine interpol_wind_meter

subroutine interpol_partoutput_value_eta(fieldname,output,j)
  implicit none
  integer, intent(in)         :: j          ! particle number
  character(2), intent(in)    :: fieldname  ! input field to interpolate over
  real, intent(inout)         :: output
  real                        :: field1(2)

  if (int(dz1out).eq.-1) then
    call find_z_level_eta(real(part(j)%zeta))
    call find_vertical_variables(uvheight,real(part(j)%zeta),induv,dz1out,dz2out,lbounds_uv,.false.)
  endif

  select case(fieldname)
    case('PR','pr')
      if (ngrid.le.0) then
        call bilinear_spatial_interpolation(prseta,field1,induv,dz1out,dz2out,nzmax)
      else
        call bilinear_spatial_interpolation_nests(prsetan,field1,induv,dz1out,dz2out,nzmax)
      endif
      call temporal_interpolation(field1(1),field1(2),output)
    case('PV','pv')
      if (ngrid.le.0) then
        call bilinear_spatial_interpolation(pveta,field1,induv,dz1out,dz2out,nzmax)
      else
        call bilinear_spatial_interpolation_nests(pvetan,field1,induv,dz1out,dz2out,nzmax)
      endif
      call temporal_interpolation(field1(1),field1(2),output)
    case('QV','qv')
      if (ngrid.le.0) then
        call bilinear_spatial_interpolation(qv,field1,induv,dz1out,dz2out,nzmax)
      else
        call bilinear_spatial_interpolation_nests(qvn,field1,induv,dz1out,dz2out,nzmax)
      endif
      call temporal_interpolation(field1(1),field1(2),output)
    case('TT','tt')
      if (ngrid.le.0) then
        call bilinear_spatial_interpolation(tteta,field1,induv,dz1out,dz2out,nzmax)
      else
        call bilinear_spatial_interpolation_nests(ttetan,field1,induv,dz1out,dz2out,nzmax)
      endif
      call temporal_interpolation(field1(1),field1(2),output)
    case('UU','uu')
      if (ngrid.le.0) then
        call bilinear_spatial_interpolation(uueta,field1,induv,dz1out,dz2out,nzmax)
      else
        call bilinear_spatial_interpolation_nests(uuetan,field1,induv,dz1out,dz2out,nzmax)
      endif
      call temporal_interpolation(field1(1),field1(2),output)
    case('VV','vv')
      if (ngrid.le.0) then
        call bilinear_spatial_interpolation(vveta,field1,induv,dz1out,dz2out,nzmax)
      else
        call bilinear_spatial_interpolation_nests(vvetan,field1,induv,dz1out,dz2out,nzmax)
      endif
      call temporal_interpolation(field1(1),field1(2),output)
    case('WW','ww')
      call find_z_level_meters(real(part(j)%z))
      call find_vertical_variables(height,real(part(j)%z),indz,dz1out,dz2out,lbounds,.false.)
      if (ngrid.le.0) then
        call bilinear_spatial_interpolation(ww,field1,induv,dz1out,dz2out,nzmax)
      else
        call bilinear_spatial_interpolation_nests(wwn,field1,induv,dz1out,dz2out,nzmax)
      endif
      call temporal_interpolation(field1(1),field1(2),output)
      dz1out = -1
    case('RH','rh')
      if (ngrid.le.0) then
        call bilinear_spatial_interpolation(rhoeta,field1,induv,dz1out,dz2out,nzmax)
      else
        call bilinear_spatial_interpolation_nests(rhoetan,field1,induv,dz1out,dz2out,nzmax)
      endif
      call temporal_interpolation(field1(1),field1(2),output)
  end select
end subroutine interpol_partoutput_value_eta

subroutine interpol_partoutput_value_meter(fieldname,output,j)
  implicit none
  integer, intent(in)         :: j          ! particle number
  character(2), intent(in)    :: fieldname  ! input field to interpolate over
  real, intent(inout)         :: output
  real                        :: field1(2)

  if (int(dz1out).eq.-1) then
    call find_z_level_meters(real(part(j)%z))
    call find_vertical_variables(height,real(part(j)%z),indz,dz1out,dz2out,lbounds,.false.)
  endif

  select case(fieldname)
    case('PR','pr')
      if (ngrid.le.0) then
        call bilinear_spatial_interpolation(prs,field1,indz,dz1out,dz2out,nzmax)
      else
        call bilinear_spatial_interpolation_nests(prsn,field1,indz,dz1out,dz2out,nzmax)
      endif
      call temporal_interpolation(field1(1),field1(2),output)
    case('PV','pv')
      if (ngrid.le.0) then
        call bilinear_spatial_interpolation(pv,field1,indz,dz1out,dz2out,nzmax)
      else
        call bilinear_spatial_interpolation_nests(pvn,field1,indz,dz1out,dz2out,nzmax)
      endif
      call temporal_interpolation(field1(1),field1(2),output)
    case('QV','qv')
      if (ngrid.le.0) then
        call bilinear_spatial_interpolation(qv,field1,indz,dz1out,dz2out,nzmax)
      else
        call bilinear_spatial_interpolation_nests(qvn,field1,indz,dz1out,dz2out,nzmax)
      endif
      call temporal_interpolation(field1(1),field1(2),output)
    case('TT','tt')
      if (ngrid.le.0) then
        call bilinear_spatial_interpolation(tt,field1,indz,dz1out,dz2out,nzmax)
      else
        call bilinear_spatial_interpolation_nests(ttn,field1,indz,dz1out,dz2out,nzmax)
      endif
      call temporal_interpolation(field1(1),field1(2),output)
    case('UU','uu')
      if (ngrid.le.0) then
        call bilinear_spatial_interpolation(uu,field1,indz,dz1out,dz2out,nzmax)
      else
        call bilinear_spatial_interpolation_nests(uun,field1,indz,dz1out,dz2out,nzmax)
      endif
      call temporal_interpolation(field1(1),field1(2),output)
    case('VV','vv')
      if (ngrid.le.0) then
        call bilinear_spatial_interpolation(vv,field1,indz,dz1out,dz2out,nzmax)
      else
        call bilinear_spatial_interpolation_nests(vvn,field1,indz,dz1out,dz2out,nzmax)
      endif
      call temporal_interpolation(field1(1),field1(2),output)
    case('WW','ww')
      if (ngrid.le.0) then
        call bilinear_spatial_interpolation(ww,field1,indz,dz1out,dz2out,nzmax)
      else
        call bilinear_spatial_interpolation_nests(wwn,field1,indz,dz1out,dz2out,nzmax)
      endif
      call temporal_interpolation(field1(1),field1(2),output)
    case('RH','rh')
      if (ngrid.le.0) then
        call bilinear_spatial_interpolation(rho,field1,indz,dz1out,dz2out,nzmax)
      else
        call bilinear_spatial_interpolation_nests(rhon,field1,indz,dz1out,dz2out,nzmax)
      endif
      call temporal_interpolation(field1(1),field1(2),output)
  end select
end subroutine interpol_partoutput_value_meter

subroutine interpol_mixinglayer_eta(zt,zteta,rhoa,rhograd)
  implicit none 
  real, intent(in)    :: zt,zteta
  real, intent(inout) :: rhoa,rhograd
  real                :: dz1w,dz2w,dz1uv,dz2uv,dz1weta,dz2weta

  call find_vertical_variables(height,zt,indz,dz1w,dz2w,lbounds,.false.)
  call find_vertical_variables(uvheight,zteta,induv,dz1uv,dz2uv,lbounds_uv,.false.)
  call find_vertical_variables(wheight,zteta,indzeta,dz1weta,dz2weta,lbounds_w,.true.)

  call vertical_interpolation(wprof(indz),wprof(indzp),dz1w,dz2w,w)
  call vertical_interpolation(uprof(induv),uprof(indpuv),dz1uv,dz2uv,u)
  call vertical_interpolation(vprof(induv),vprof(indpuv),dz1uv,dz2uv,v)
  call vertical_interpolation(rhoprof(induv),rhoprof(indpuv),dz1uv,dz2uv,rhoa)
  call vertical_interpolation(rhogradprof(induv),rhogradprof(indpuv),dz1uv,dz2uv,rhograd)
  call vertical_interpolation(wprofeta(indzeta),wprofeta(indzpeta),dz1weta,dz2weta,weta)
end subroutine interpol_mixinglayer_eta

subroutine standard_deviation_eta(iw,iuv,iweta)
  ! Standard deviation of surrounding grid points
  ! Only used in mesoscale turbulence calculations
  !***********************************************
  implicit none

  integer,intent(in)  :: iw(2),iuv(2),iweta(2)
  real                :: wsl,wsq,wxaux,usl,usq,uxaux,vsl,vsq,vxaux,wetasl,wetasq,wetaxaux
  integer             :: n,m
  real,parameter      :: eps=1.0e-30
  ! Standard deviations
  !********************
  wsl=0.
  wsq=0.
  usl=0.
  usq=0.
  vsl=0.
  vsq=0.
  wetasl=0.
  wetasq=0.

  if (ngrid.le.0) then ! No nest  
    do m=1,2
      do n=1,2
        call compute_sl_sq(ww,wsl,wsq,iw(n),memind(m),nzmax)
        call compute_sl_sq(wweta,wetasl,wetasq,iweta(n),memind(m),nzmax)
        if (ngrid.lt.0) then
          call compute_sl_sq(uupoleta,usl,usq,iuv(n),memind(m),nzmax)
          call compute_sl_sq(vvpoleta,vsl,vsq,iuv(n),memind(m),nzmax)
        else
          call compute_sl_sq(uueta,usl,usq,iuv(n),memind(m),nzmax)
          call compute_sl_sq(vveta,vsl,vsq,iuv(n),memind(m),nzmax)
        endif
      end do
    end do
  else ! Nest
    do m=1,2
      do n=1,2
        call compute_sl_sq_nests(wwn,wsl,wsq,iw(n),memind(m),nzmax)
        call compute_sl_sq_nests(wwetan,wetasl,wetasq,iweta(n),memind(m),nzmax)
        call compute_sl_sq_nests(uuetan,usl,usq,iuv(n),memind(m),nzmax)
        call compute_sl_sq_nests(vvetan,vsl,vsq,iuv(n),memind(m),nzmax)
      end do
    end do
  endif

  call standard_deviation(wsl,wsq,16.,wsig)
  call standard_deviation(usl,usq,16.,usig)
  call standard_deviation(vsl,vsq,16.,vsig)
  call standard_deviation(wetasl,wetasq,16.,wsigeta)
end subroutine standard_deviation_eta

subroutine standard_deviation_meter(iw)
  ! Standard deviation of surrounding grid points
  ! Only used in mesoscale turbulence calculations
  !***********************************************
  implicit none

  integer,intent(in)  :: iw(2)
  real                :: wsl,wsq,wxaux,usl,usq,uxaux,vsl,vsq,vxaux
  integer             :: n,m
  real,parameter      :: eps=1.0e-30

  ! Standard deviations
  !********************
  wsl=0.
  wsq=0.
  usl=0.
  usq=0.
  vsl=0.
  vsq=0.

  if (ngrid.le.0) then ! No nest  
    do m=1,2
      do n=1,2
        call compute_sl_sq(ww,wsl,wsq,iw(n),memind(m),nzmax)
        if (ngrid.lt.0) then
          call compute_sl_sq(uupol,usl,usq,iw(n),memind(m),nzmax)
          call compute_sl_sq(vvpol,vsl,vsq,iw(n),memind(m),nzmax)
        else
          call compute_sl_sq(uu,usl,usq,iw(n),memind(m),nzmax)
          call compute_sl_sq(vv,vsl,vsq,iw(n),memind(m),nzmax)
        endif
      end do
    end do
  else ! Nest
    do m=1,2
      do n=1,2
        call compute_sl_sq_nests(wwn,wsl,wsq,iw(n),memind(m),nzmax)
        call compute_sl_sq_nests(uun,usl,usq,iw(n),memind(m),nzmax)
        call compute_sl_sq_nests(vvn,vsl,vsq,iw(n),memind(m),nzmax)
      end do
    end do
  endif

  call standard_deviation(wsl,wsq,16.,wsig)
  call standard_deviation(usl,usq,16.,usig)
  call standard_deviation(vsl,vsq,16.,vsig)
end subroutine standard_deviation_meter

end module interpol_mod