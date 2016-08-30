SUBROUTINE FILT2D(CH,IV,IP,IC)
!------------------------------------------------------------------------------------------------
! Purpose: PERFORM 2-D FILTER
! Original version : FILT2D.f (IAP 9L)
! Reconstructed    : ZhangHe
! Completed : 2005.9.2
! Update: 2007.5.14, ZhangHe, specify which variables or subroutines to use only
!------------------------------------------------------------------------------------------------
   use shr_kind_mod, only: r8 => shr_kind_r8
   use IAP_grid, only: NX, NY, JB, JE, period
   use smoother, only: SOSS1D, SOSS2D
   use Filt,     only: FILTER
   use pmgrid,   only: beglatdyn, endlatdyn

   implicit none
!------------------------------------Arguments---------------------------------------------------
   real(r8), intent(inout) :: CH(NX,NY)     ! input variable needed to be filtered
!!   real(r8), intent(inout) :: CH(NX,beglatdyn:endlatdyn)    ! input variable needed to be filtered
   integer , intent(in)    :: IV, IP, IC    ! index
!----------------------------------Local workspace-----------------------------------------------
   integer  :: J, i             ! loop index
   integer  :: ID,JN,JS
!------------------------------------------------------------------------------------------------

!  APPLY THE SPHERICAL CYCLICITY
   DO J = beglatdyn, endlatdyn
      call period( CH(1,J) )
   END DO
   ID = IP*2 - 1
!------------------------------------------------------------------------------      
   CALL FILTER( CH,ID )
!------------------------------------------------------------------------------      
   IF ( IC.GE.1 ) THEN
      IF ( IV.EQ.1 ) THEN  ! do second-order Shapiro smoother of 1D (LAT) for u & v
         JN = JB - IP + 1
         JS = JE
!------------------------------------------------------------------------------      
         CALL SOSS1D( CH,JN,JS,IC )  ! at module smoother
!------------------------------------------------------------------------------      
      ELSE  ! do second-order Shapiro smoother of 2D (LAT & LON) for PS, T & q
!------------------------------------------------------------------------------      
         CALL SOSS2D( CH,IV,IP,IC )  
!------------------------------------------------------------------------------      
      ENDIF
   ENDIF

   RETURN
END
