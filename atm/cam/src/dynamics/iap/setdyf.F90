!=================================================================================
SUBROUTINE SETDYF
!---------------------------------------------------------------------------------
! Purpose: SET CONSTANTS USED IN DYNAMIC FRAME INTEGRATION
! Original version: ENTRY SETDYF at DYFRAM.f (IAP 9L)
! Reconstructed & revised : ZhangHe
! Completed : 2005.10.13
! Update    : 2006.10.10, Zhanghe, revised the calling of sub. CONPDA     
!           : 2006.12.12, Zhanghe, revised the calling of sub. STFRAM     
!           : 2007.05.03, Zhanghe,
!           : 2007.08.25, Zhanghe, added the calling of sub. setfle
!           : 2007.11.11, Zhanghe, moved the calling of STFRAM to sub. stepon 
!                                  moved the calling of setfle to sub. STFRAM
!           : 2013-03-27, ZhangHe, removed CONPDA
!---------------------------------------------------------------------------------
   use shr_kind_mod, only: r8 => shr_kind_r8
   use flexib,     only: EPS, IDSOSS 
   use stdatm,     only: SETMSA     
   use Dyn_const,  only: DLAT, DLON        
   use Engy_const, only: STENGC
   use smoother,   only: STSMOT
   use sm9h,       only: STSM9C
   use Filt,       only: CONFIL
   use hdif,       only: STDFSC

   IMPLICIT NONE
!------------------------------------------------------------------------------------------------

!  FOR THE MODEL STANDARD ATMOSPHERE
!--------------- SET CONST OF THE MODEL STANDARD ATMOSPHERE ----------------------
   CALL SETMSA( EPS )  
!---------------------------------------------------------------------------------
!------------------ FOR THE ENERGY BUGET CALCULATION ----------------------------
   CALL STENGC    
!---------------------------------------------------------------------------------
!
!----------------------------- FOR THE SMOOTHERS/FILTERS ---------------------------------
   CALL STSMOT   
!-----------------------------------------------------------------------------------------
!--------------------- SET CONSTANTS USED IN SUB.SM9HAS ----------------------------------
   CALL STSM9C(1)
!-----------------------------------------------------------------------------------------
!--------------------- SET CONSTANTS USED IN SUB.FILTER ----------------------------------
   CALL CONFIL( IDSOSS,DLAT,DLON )  
!-----------------------------------------------------------------------------------------
!------------- Set constants for computation of the horizontal diffusion ----------------- 
   CALL STDFSC       
!-----------------------------------------------------------------------------------------
   RETURN
END

