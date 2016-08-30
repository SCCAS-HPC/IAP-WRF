module tendency
!------------------------------------------------------------------------------------------------
! Purpose:    Set the tendency of prognostic variables  
! Author :    ZhangHe
! Completed : 2005.8.26
! Update: 2007.4.23, ZhangHe, changed 'LADDSS' to 'NADDSS'
!         2007.7.20, ZhangHe, added the data of NADDSS
!         2008.6.11, ZhangHe, make the arrys allocatable & allocated them at 
!                             sub. initialize_IAPprog in IAP_prog.F90
!         2013-03-29, ZhangHe, removed DQ, DQliq, and DQice
!------------------------------------------------------------------------------------------------
   use shr_kind_mod, only: r8 => shr_kind_r8

   implicit none
   save
   public

   real(r8), allocatable :: DPsa(:,:)     ! tendency of Psa
   real(r8), allocatable :: DPstar1(:,:)  ! tendency of Pes*
   real(r8), allocatable :: DU(:,:,:)     ! tendency of UT
   real(r8), allocatable :: DV(:,:,:)     ! tendency of VT
   real(r8), allocatable :: DT(:,:,:)     ! tendency of TT
   real(r8), allocatable :: SU(:,:,:)     ! soure & sink tendency of UT
   real(r8), allocatable :: SV(:,:,:)     ! soure & sink tendency of VT
   real(r8), allocatable :: ST(:,:,:)     ! soure & sink tendency of TT
   integer  :: NADDSS         ! type of adding source & sink
!                               NADDSS = 0 , add source & sink by physics timestep;  
!                               NADDSS = 1 , add source & sink by adaption timestep;  
!                               NADDSS = 2 , add source & sink by advective timestep;  
!  Set the data of NADDSS
   data NADDSS  / 1 /
!!   data NADDSS  / 4 /

end module
