SUBROUTINE DIAGPP
!------------------------------------------------------------------------------------------------
! Purpose: COMPUTE PRESSURE AT MODEL & INTERFACE SIGMA LAYERS
! Original version : DIAGPP.f (IAP 9L)
! Reconstructed & supplemented : ZhangHe, added the calculation of PIN
! Completed : 2005.9.1
! Update: October,2012, Jiang Jinrong, 2D parellel

!------------------------------------------------------------------------------------------------
   use IAP_grid,   only : NX, NY, NL, NZ
   use Dyn_const,  only : PMTOP, SIG, SIGL
   use IAP_prog,   only : P, PLY, PIN
   use pmgrid,     only : beglatdyn, endlatdyn,beglev,endlev

   implicit none
!----------------------------------Local workspace-----------------------------------------------
   integer  :: I,J,K             ! loop index
!------------------------------------------------------------------------------------------------
!
!---------------------------- compute model layer pressure ------------------------------
   DO J = beglatdyn, endlatdyn
      DO K = 1 ,NL
         DO I = 1 ,NX
            PLY(I,K,J) = P(I,J)*SIGL(K) + PMTOP
         END DO
      END DO
      DO I = 1 ,NX
         PLY(I,NZ,J)   = P(I,J)         + PMTOP
      END DO
   END DO
!---------------------------- compute interface layer pressure ------------------------------
   do J = beglatdyn, endlatdyn
      do K = 2 ,NL      
         do I = 1 ,NX
            PIN(I,K,J) = P(I,J)*SIG(K)  + PMTOP
         end do
      end do
	  
      do I = 1 ,NX
         PIN(I,1,J)   = PMTOP
         PIN(I,NZ,J)  = PMTOP + P(I,J)
      end do
   end do
      
   RETURN
END
