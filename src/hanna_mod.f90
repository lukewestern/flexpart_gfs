module hanna_mod

  implicit none

  real :: ust,wst,ol,h,zeta,sigu,sigv,tlu,tlv,tlw
  real :: sigw,dsigwdz,dsigw2dz

  ! openmp change
!$OMP THREADPRIVATE(ust,wst,ol,h,zeta,sigu,sigv,tlu,tlv,tlw, &
!$OMP sigw,dsigwdz,dsigw2dz)
  ! openmp change end

end module hanna_mod
