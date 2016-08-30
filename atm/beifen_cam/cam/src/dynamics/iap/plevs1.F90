subroutine plevs00 (ncol    , ncold   ,nver    ,ps      ,pint    , &
                   pmid    ,pdel)

!-----------------------------------------------------------------------
! Purpose:
! Define the pressures of the interfaces and midpoints from the
! coordinate definitions and the surface pressure.
!
! Method:
!
! Author: B. Boville
!
! Revised by: ZhangHe, 2007.4.26
!           : Jiang Jinrong, 2012 October, for 2D version
!-----------------------------------------------------------------------
!
! $Id: plevs0.F90,v 1.1.2.2 2004/09/23 17:20:09 erik Exp $
! $Author: erik $
!
!-----------------------------------------------------------------------

  use shr_kind_mod, only: r8 => shr_kind_r8
  use Dyn_const,    only: PMTOP, SIG, SIGL
  use pmgrid !jjr     

  implicit none
!zhh #include <comhyb.h>

!----------------------------------Arguments------------------------------------
  integer , intent(in)  :: ncol               ! Longitude dimension
  integer , intent(in)  :: ncold              ! Declared longitude dimension
  integer , intent(in)  :: nver               ! vertical dimension
  real(r8), intent(in)  :: ps(ncold)          ! Surface pressure (pascals)
  real(r8), intent(out) :: pint(ncold,beglev:endlev+1) ! Pressure at model interfaces
  real(r8), intent(out) :: pmid(ncold,beglev:endlev)   ! Pressure at model levels
  real(r8), intent(out) :: pdel(ncold,beglev:endlev)   ! Layer thickness (pint(k+1) - pint(k))
!-----------------------------------------------------------------------

!---------------------------Local workspace-----------------------------
  integer i,k             ! Longitude, level indices
!-----------------------------------------------------------------------
!
! Set interface pressures
!
  do k=beglev,endlev+1
     do i=1,ncol
        pint(i,k) = SIG(k)*( ps(i)-PMTOP*100.0 ) + PMTOP*100.0
     end do
  end do
!
! Set midpoint pressures and layer thicknesses
!
  do k=beglev,endlev
     do i=1,ncol
        pmid(i,k) = SIGL(k)*( ps(i)-PMTOP*100.0 ) + PMTOP*100.0
        pdel(i,k) = pint(i,k+1) - pint(i,k)
     end do
  end do

  return
end subroutine plevs00

