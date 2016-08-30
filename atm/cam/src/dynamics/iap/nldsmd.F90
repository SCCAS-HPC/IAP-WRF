SUBROUTINE NLDSMD(X,F,EPS)
!------------------------------------------------------------------------------------------------
! Purpose: DESCENT METHOD SOLVING NONLINEAR EQUATIONS FOR ONE POINT X
! Original version : NLDSMD.f (IAP 9L)
! Reconstructed  : ZhangHe
! Completed : 2005.9.8
!------------------------------------------------------------------------------------------------
   use shr_kind_mod, only: r8 => shr_kind_r8

   implicit none
!------------------------------------Arguments--------------------------------------------------
   real(r8), intent(in)  :: EPS
   real(r8), intent(out) :: X
!----------------------------------Local workspace-----------------------------------------------
   real(r8) :: R, DXI, DFDX, A, RLD
   real(r8) :: F
   EXTERNAL  F
!------------------------------------------------------------------------------------------------
!
10 R     = F( X )  
!   JUDGING IF THE ITERATION HAS FINISHED
   IF ( R.LT.EPS )  RETURN
!
!     COMPUTE PARTIAL DERIVATIVE OF F & LAMDA:DFDX & RLD
   DXI   = MAX(1.0D-6*X,1.0D-12)
   DFDX  = (F( X+DXI ) - R) / DXI
   A     = DFDX*DFDX + 1.0D-20
   RLD   = R / A
!
!     FINDING A NEW TRIAL SOLUTION  
   X     = X - RLD*DFDX
   GOTO 10
END
