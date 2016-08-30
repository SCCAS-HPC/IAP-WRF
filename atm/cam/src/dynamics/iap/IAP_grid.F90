module IAP_grid
!----------------------------------------------------------------------------------------------
! Purpose: Parameters and variables related to the dynamics grid
! Author : ZhangHe
! Completed : 2005.8.23
! Update : 2006.12.28, ZhangHe, use 'parame.inc' instead of statement in this module
!          2007.5.10, ZhangHe, change the module name 'Dyn_grid' ==> 'IAP_grid'
!          2008.05, WuJianping, delete sub. period3 
!----------------------------------------------------------------------------------------------
   use pmgrid, only : plon, plat, plev
   implicit none

   save
   public

   integer, parameter :: NLON = PLON
   integer, parameter :: NLAT = PLAT
   integer, parameter :: NLAY = PLEV
   integer, parameter :: IM    = NLON
   integer, parameter :: NY    = NLAT
   integer, parameter :: NL    = NLAY
   integer, parameter :: EX    = 3            ! extended grids due to leap-frog difference scheme 
   integer, parameter :: NX    = NLON + 2*EX  ! extended domain longitude
   integer, parameter :: NZ    = NL + 1       ! number of interface levels
   integer, parameter :: IB    = 1  + EX      ! begin index of longitude for computation
   integer, parameter :: IE    = NX - EX      ! end   index of longitude for computation
   integer, parameter :: JB    = 2            ! begin index of latitude for computation
   integer, parameter :: JE    = NY - 1       ! end   index of latitude for computation
   integer, parameter :: KE    = NL + 2       ! extended domain vertical levels
   integer, parameter :: NM    = NL - 1
   integer, parameter :: IB2   = IB + 1       ! uesd for periodic condition
   integer, parameter :: IB3   = IB + 2       ! uesd for periodic condition
   integer, parameter :: IE1   = IE - 2       ! uesd for periodic condition
   integer, parameter :: IE2   = IE - 1       ! uesd for periodic condition
   integer, parameter :: NX1   = IE + 1       ! uesd for periodic condition
   integer, parameter :: NX2   = IE + 2       ! uesd for periodic condition
   integer, parameter :: MGH   = NX * NY
    
!================================================================================================
CONTAINS
!================================================================================================

!================================================================================================
   subroutine period( F )
!-----------------------------------------------------------------------------------------------
!     Set the spherical cyclicity condition for variable F
!-----------------------------------------------------------------------------------------------
      use shr_kind_mod, only: r8 => shr_kind_r8
      implicit none 
!------------------------------------Arguments------------------------------------------------
      real(r8), intent(inout) :: F( NX )    ! input variable
!---------------------------------------------------------------------------------------------
      F(1  )  = F(IE1)                 
      F(2  )  = F(IE2)
      F(3  )  = F(IE )
      F(NX1)  = F(IB )
      F(NX2)  = F(IB2)
      F(NX )  = F(IB3) 
      return
   end subroutine

!================================================================================================
   subroutine periodp( F )
!-----------------------------------------------------------------------------------------------
!     Set the spherical cyclicity condition for variable F
!-----------------------------------------------------------------------------------------------
      use shr_kind_mod, only: r8 => shr_kind_r8
      implicit none 
!------------------------------------Arguments------------------------------------------------
      real(r8), intent(inout) :: F( NX )    ! input variable
!---------------------------------------------------------------------------------------------
      F(IE1)  = F(1  )                
      F(IE2)  = F(2  ) 
      F(IE )  = F(3  ) 
      F(NX1)  = F(IB )
      F(NX2)  = F(IB2)
      F(NX )  = F(IB3) 
      return
   end subroutine

!================================================================================================
end module
