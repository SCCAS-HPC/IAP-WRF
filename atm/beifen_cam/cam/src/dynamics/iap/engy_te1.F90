subroutine engy_te1(cwava   ,w       ,t       ,u      ,v        , &
                   phis    ,pdel    ,engy    , nlon  )

!----------------------------------------------------------------------- 
! 
! Purpose: 
! Calculate contribution of current latitude to total energy
! 
! Method: 
! 
! Author: J. Olson
! 
!-----------------------------------------------------------------------
!
! $Id: engy_te.F90,v 1.1.2.2 2004/09/23 17:20:08 erik Exp $
! $Author: erik $
!
!-----------------------------------------------------------------------
!
  use shr_kind_mod, only: r8 => shr_kind_r8
  use pmgrid,       only: plev, plon,beglev,endlev
  use physconst, only: cpair

  implicit none
!
!------------------------------Arguments--------------------------------
!
  integer , intent(in)  :: nlon                 ! longitude dimension  
  real(r8), intent(in)  :: cwava                ! normalization factor    l/(g*plon)
  real(r8), intent(in)  :: w                    ! gaussian weight this latitude
  real(r8), intent(in)  :: t   (plon,beglev:endlev)      ! temperature
  real(r8), intent(in)  :: u   (plon,beglev:endlev)      ! u-component
  real(r8), intent(in)  :: v   (plon,beglev:endlev)      ! v-component
  real(r8), intent(in)  :: phis(plon)           ! Geopotential
  real(r8), intent(in)  :: pdel(plon,beglev:endlev)      ! pressure diff between interfaces
  real(r8), intent(out) :: engy                 ! accumulator
!
!---------------------------Local variables-----------------------------
!
  integer i,k               ! longitude, level indices
  real(r8) const            ! temporary constant
!
!-----------------------------------------------------------------------
!
! Integration factor (the 0.5 factor arises because gaussian weights sum to 2)
!
  const = cwava*w*0.5
  engy = 0.
!
  do k=beglev,endlev
     do i=1,nlon
        engy = engy + ( cpair*t(i,k)                          + &
                        0.5*( u(i,k)*u(i,k) + v(i,k)*v(i,k) ) + &
                        phis(i) )*pdel(i,k)
     end do
  end do

  engy = engy*const

  return
end subroutine engy_te1
