module IAP_prog
!------------------------------------------------------------------------------------------------
! Purpose:    Set meteorological variables  
! Author :    ZhangHe
! Completed : 2005.8.26
! Update    : 2007.4.25, ZhangHe, corrected the unit of pressure from 'Pa' to 'hPa'
!             2007.5.10, ZhangHe, changed the module name 'prognostic' ==> 'IAP_prog'
!             2007.5.14, ZhangHe, changed 'GHS(IM,NY)' ==> 'GHS(NX,NY)'
!             2008.6.9, ZhangHe, parallel version: NY ==> beglatdynex:endlatdynex
!             2008.6.10, ZhangHe, added sub. initialize_IAPprog
!             2011.7.10 Juanxiong He, adapted to CESM 	
!             2012 July Jiang Jinrong, 2D parellel
!             2013-03-29, Zhang He, removed Q, Qliq, Qice, QT, QTliq, QTice
!------------------------------------------------------------------------------------------------
   use shr_kind_mod, only: r8 => shr_kind_r8
   use IAP_grid,     only: NX, NY, NL, NZ
   use pmgrid,       only: beglatdynex, endlatdynex, beglev, endlev
   use mathconst,    only: ZERO
   use infnan,       only: inf

   implicit none
   save
   public

   real(r8), allocatable :: P(:,:)        ! P = Pes = Ps - Pt  (unit: hPa )
   real(r8), allocatable :: Psa(:,:)      ! deviation from average surface pressure (P')  (unit: hPa  )
   real(r8), allocatable :: U(:,:,:)      ! latitudinal wind   (unit: m/s)
   real(r8), allocatable :: V(:,:,:)      ! meridional wind    (unit: m/s)
   real(r8), allocatable :: T(:,:,:)      ! temperature        (unit: K  )
   real(r8), allocatable :: Tsa(:,:,:)    ! deviation from standard temperature  (unit: K  )
   real(r8), allocatable :: Q(:,:,:)      ! specific humidity (water vapor) (unit:    )
   real(r8), allocatable :: PT(:,:)       ! IAP transformation parameter, PT = sqrt(Pes* / P0)
   real(r8), allocatable :: UT(:,:,:)     ! UT = PT*U
   real(r8), allocatable :: VT(:,:,:)     ! VT = PT*V
   real(r8), allocatable :: TT(:,:,:)     ! TT = PT*R*T' / b
!!   real(r8), allocatable :: QT(:,:,:)     ! QT = Pes*Q     ! 2007.4.17
   real(r8), allocatable :: WST(:,:,:)    ! WST = PT*WS
   real(r8), allocatable :: WS(:,:,:)     ! sigma-surface vertical velocity  (unit: s^-1)
   real(r8), allocatable :: Ustar(:,:,:)  ! flexible substitute of U 
   real(r8), allocatable :: Vstar(:,:,:)  ! flexible substitute of V
   real(r8), allocatable :: WSstar(:,:,:) ! flexible substitute of WS
   real(r8), allocatable :: Pstar1(:,:)   ! flexible substitute of Pes
   real(r8), allocatable :: Pstar2(:,:)   ! flexible substitute of Pes
   real(r8), allocatable :: PIN(:,:,:)    ! pressure at interface sigma layer   (unit: hPa)
   real(r8), allocatable :: PLY(:,:,:)    ! pressure at model sigma layer       (unit: hPa)
   real(r8), allocatable :: PS2(:,:)      ! surface pressure                    (unit: hPa)
   real(r8), allocatable :: GHI(:,:,:)    ! geopotential departure (= g*z')           (unit: m^2/s^2)
   real(r8), allocatable :: GZ(:,:,:)     ! geopotential above surface ( = g*(z-zs) ) (unit: m^2/s^2)
   real(r8), allocatable :: GHS(:,:)      ! surface geopotential                      (unit: m^2/s^2)
   real(r8), allocatable :: GZsm(:,:)     ! GZsm = g*Zsm                              (unit: m^2/s^2)
   real(r8), allocatable :: Zg (:,:,:)    ! geopotential height                       (unit: m )
   real(r8), allocatable :: WPA(:,:,:)    ! time averaged P-surface vertical velocity (unit: hPa/s)
   real(r8), allocatable :: WPV(:,:,:)    ! P-surface vertical velocity               (unit: hPa/s)
   real(r8), allocatable :: deltap(:,:,:) ! deltap = pt / p, if pt = 0, then deltap = 0
! ================================== zhh ====================================
   real(r8), allocatable :: Qliq(:,:,:)   ! cloud liquid water (unit: kg/kg)
   real(r8), allocatable :: Qice(:,:,:)   ! cloud ice water (unit: kg/kg)
!!   real(r8), allocatable :: QTliq(:,:,:)  ! QTliq = Pes*Qliq
!!   real(r8), allocatable :: QTice(:,:,:)  ! QTice = Pes*Qice
! ===============================  2007.12.20 ================================
   real(r8), allocatable :: UVW0(:,:,:,:) ! wind velocity used for Q advection
 
!================================================================================================
CONTAINS
!================================================================================================

!================================================================================================
   subroutine initialize_IAPprog
!-----------------------------------------------------------------------------------------------
! Purpose:  Allocate and initialize the prognostic arrays.
!-----------------------------------------------------------------------------------------------
   use tendency
!
      allocate ( P     (NX   ,beglatdynex:endlatdynex) )
      allocate ( Psa   (NX   ,beglatdynex:endlatdynex) )
      allocate ( U     (NX,beglev:endlev,beglatdynex:endlatdynex) )
      allocate ( V     (NX,beglev:endlev,beglatdynex:endlatdynex) )
      allocate ( T     (NX,beglev:endlev,beglatdynex:endlatdynex) )
      allocate ( Tsa   (NX,beglev:endlev,beglatdynex:endlatdynex) )
      allocate ( Q     (NX,beglev:endlev,beglatdynex:endlatdynex) )
      allocate ( PT    (NX   ,beglatdynex:endlatdynex) )
      allocate ( UT    (NX,beglev:endlev,beglatdynex:endlatdynex) )
      allocate ( VT    (NX,beglev:endlev,beglatdynex:endlatdynex) )
      allocate ( TT    (NX,beglev:endlev,beglatdynex:endlatdynex) )
!      allocate ( QT    (NX,beglev:endlev,beglatdynex:endlatdynex) )
      allocate ( WST   (NX,beglev:endlev+1,beglatdynex:endlatdynex) )
      allocate ( WS    (NX,beglev:endlev+1,beglatdynex:endlatdynex) )
      allocate ( Ustar (NX,beglev:endlev,beglatdynex:endlatdynex) )
      allocate ( Vstar (NX,beglev:endlev,beglatdynex:endlatdynex) )
      allocate ( WSstar(NX,beglev:endlev+1,beglatdynex:endlatdynex) ) !jjr
      allocate ( Pstar1(NX   ,beglatdynex:endlatdynex) )
      allocate ( Pstar2(NX   ,beglatdynex:endlatdynex) )
      allocate ( PIN   (NX,NZ,beglatdynex:endlatdynex) ) !jjr no need change
      allocate ( PLY   (NX,NZ,beglatdynex:endlatdynex) )!jjr NZ, in
      allocate ( PS2   (NX   ,beglatdynex:endlatdynex) )
      allocate ( GHI   (NX,beglev:endlev+1,beglatdynex:endlatdynex) ) 
      allocate ( GZ    (NX,beglev:endlev+1,beglatdynex:endlatdynex) )
      allocate ( GHS   (NX   ,beglatdynex:endlatdynex) )
      allocate ( GZsm  (NX   ,beglatdynex:endlatdynex) )
      allocate ( Zg    (NX,beglev:endlev+1,beglatdynex:endlatdynex) )
      allocate ( WPA   (NX,beglev:endlev+1,beglatdynex:endlatdynex) )
      allocate ( WPV   (NX,beglev:endlev,beglatdynex:endlatdynex) )
      allocate ( deltap(NX,NL,beglatdynex:endlatdynex) ) !jjr no need change,in
      allocate ( Qliq  (NX,NL,beglatdynex:endlatdynex) )
      allocate ( Qice  (NX,NL,beglatdynex:endlatdynex) )
!!      allocate ( QTliq (NX,NL,beglatdynex:endlatdynex) )
!!      allocate ( QTice (NX,NL,beglatdynex:endlatdynex) )
      allocate ( UVW0  (NX,beglev:endlev,beglatdynex:endlatdynex,3) )
!
! for tendency
      allocate ( DPsa  (NX   ,beglatdynex:endlatdynex) )
      allocate ( DPstar1(NX   ,beglatdynex:endlatdynex) )
      allocate ( DU    (NX,beglev:endlev,beglatdynex:endlatdynex) )
      allocate ( DV    (NX,beglev:endlev,beglatdynex:endlatdynex) )
      allocate ( DT    (NX,beglev:endlev,beglatdynex:endlatdynex) )
!      allocate ( DQ    (NX,beglev:endlev,beglatdynex:endlatdynex) )
      allocate ( SU    (NX,beglev:endlev,beglatdynex:endlatdynex) )
      allocate ( SV    (NX,beglev:endlev,beglatdynex:endlatdynex) )
      allocate ( ST    (NX,beglev:endlev,beglatdynex:endlatdynex) )
!
      P(:,:)        = inf
      Psa(:,:)      = inf
      U(:,:,:)      = inf
      V(:,:,:)      = inf
      T(:,:,:)      = inf
      Tsa(:,:,:)    = inf
!      Q(:,:,:)      = inf
      PT(:,:)       = inf
      UT(:,:,:)     = inf
      VT(:,:,:)     = inf
      TT(:,:,:)     = inf
!      QT(:,:,:)     = inf
      WST(:,:,:)    = ZERO
      WS(:,:,:)     = ZERO
      Ustar(:,:,:)  = inf
      Vstar(:,:,:)  = inf
      WSstar(:,:,:) = ZERO
      Pstar1(:,:)   = inf
      Pstar2(:,:)   = inf
      PIN(:,:,:)    = inf
      PLY(:,:,:)    = inf
      PS2(:,:)      = inf
      GHI(:,:,:)    = inf
      GZ(:,:,:)     = inf
      GHS(:,:)      = inf
      GZsm(:,:)     = ZERO
      Zg(:,:,:)     = inf
      WPA(:,:,:)    = ZERO
      WPV(:,:,:)    = ZERO
      deltap(:,:,:) = inf
      Qliq(:,:,:)   = inf
      Qice(:,:,:)   = inf
!!      QTliq(:,:,:)  = inf
!!      QTice(:,:,:)  = inf
      UVW0(:,:,:,:) = ZERO
!
! for tendency
      DPsa(:,:)     = inf
      DPstar1(:,:)  = inf
      DU(:,:,:)     = inf
      DV(:,:,:)     = inf
      DT(:,:,:)     = inf
!      DQ(:,:,:)     = inf
      SU(:,:,:)     = inf
      SV(:,:,:)     = inf
      ST(:,:,:)     = inf

      return
   end subroutine initialize_IAPprog

end module
