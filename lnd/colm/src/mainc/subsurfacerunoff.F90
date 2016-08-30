
 subroutine subsurfacerunoff (nl_soil,dtime,pondmx,dzmm,&
                              wliq,eff_porosity,hk,dhkdw,dwat,rsubst) 
	           
!=======================================================================
! Original author : Yongjiu Dai, 07/29/2002
!=======================================================================

  use precision
  implicit none
    
!-----------------------Arguments---------------------------------------

  integer, INTENT(in) :: nl_soil   ! number of soil layers
  real(r8), INTENT(in) :: &
         dtime,           &! time step (s)
         pondmx,          &! ponding depth (mm)
         dzmm(1:nl_soil), &! layer thickness (mm)
 eff_porosity(1:nl_soil), &! effective porosity = porosity - vol_ice
           hk(1:nl_soil), &! hydraulic conductivity (mm h2o/s)
        dhkdw(1:nl_soil), &! d(hk)/d(vol_liq)
         dwat(1:nl_soil)   ! change in soil water

  real(r8), INTENT(inout) :: wliq(1:nl_soil)  ! liquid water (kg/m2)
  real(r8), INTENT(out) :: rsubst             ! subsurface runoff (mm h2o/s)

!-----------------------Local Variables---------------------------------

  real(r8) :: xs           ! excess soil water above saturation
  integer  :: i            ! loop counter

!-----------------------End Variable List-------------------------------

      rsubst = 0.

! determine water in excess of saturation and
! add excess water to the above soil layer like the bucket
! any left over put into runoff

      xs = 0.
      do i = nl_soil, 2, -1
         wliq(i) = wliq(i) + xs
         xs = max(0., wliq(i)-eff_porosity(i)*dzmm(i))    ! [mm]
         wliq(i) = min(eff_porosity(i)*dzmm(i), wliq(i))
      end do

      wliq(1) = wliq(1) + xs
      xs = max(0., wliq(1)-(pondmx+eff_porosity(1)*dzmm(1)))
      wliq(1) = min(pondmx+eff_porosity(1)*dzmm(1), wliq(1))

! sub-surface runoff and drainage 
      rsubst = rsubst + xs/dtime &
             + hk(nl_soil) + dhkdw(nl_soil)*dwat(nl_soil) ! [mm/s]
      
 end subroutine subsurfacerunoff
