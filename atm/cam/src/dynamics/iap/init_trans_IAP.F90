subroutine init_trans_IAP
!----------------------------------------------------------------------------------
! Purpose: Transform P, U, V, T to PT, UT, VT, TT at first time step.
! Author : ZhangHe
! Completed: 2007.5.8
! Update : ZhangHe, 2007.12.21
!        : Juanxiong He, 201008
!        : Jiang Jinrong, October 2012
! Reviewed: ZhangHe, 2011-11-18
!           ZhangHe, 2012-10-26
! Modified: Zhang He, 2013-03-12, removed transform of Q, Qliq, and Qice
!           Zhang He, 2013-04-17, moved DIAGPP and DIAGBB from MDINIT 
!----------------------------------------------------------------------------------

   use shr_kind_mod, only: r8 => shr_kind_r8
   use IAP_grid,   only: NX, NY, NL
   use pmgrid,     only: beglatdyn, endlatdyn,beglev,endlev
!zhh   use spmd_utils,   only: iam, masterproc  !juanxiong he

   use stdatm,     only: PSB, P00, TB, DIAGBB
   use Trans_coef, only: PTU, PTV, PTT, TRANSC 
   use IAP_prog,   only: U, V, T, P, PS2, Psa, Pstar1, PT, UT, VT, TT
   use Dyn_const,  only: PMTOP   

   implicit none
!------------------------------Local workspace--------------------------------
   real(r8) :: sqP0
   integer  :: I, J, K
!-----------------------------------------------------------------------------
!
!!   sqP0 = sqrt(P00)
! 
   DO J = beglatdyn, endlatdyn
      DO I = 1, NX
         P  (I,J) = PS2(I,J) - PMTOP
         Psa(I,J) = PS2(I,J) - PSB(I,J)
         Pstar1(I,J) = P(I,J)
!!         PT(I,J)     = sqrt( Pstar1(I,J) ) / sqP0
         PT(I,J)     = sqrt( Pstar1(I,J)/P00 )
      END DO
   END DO       
!
!------------------ COMPUTE PRESSURE AT MODEL & INTERFACE SIGMA LAYERS ----------------------
   CALL DIAGPP
!--------------- SET MSA TEMPERATURE, GEOPOTENTIAL & CHARACTERISTIC PHASE-SPEED -------------
   CALL DIAGBB( 3 )
!------------------- CALCULATE TRANSFORMATION COEFFECIENTS -------------------------
   CALL TRANSC   
!-----------------------------------------------------------------------------------
!
   DO J = beglatdyn, endlatdyn
      DO K = beglev,endlev 
         DO I = 1, NX
            UT(I,K,J) = U(I,K,J) * PTU(I,J)
            VT(I,K,J) = V(I,K,J) * PTV(I,J)
            TT(I,K,J) = ( T(I,K,J)-TB(I,K,J) ) * PTT(I,J)
         END DO
      END DO
   END DO

   return
end subroutine init_trans_IAP
