SUBROUTINE nliter_uvt(CCU,CCV,CCT,CT)
!------------------------------------------------------------------------------------------------
! Purpose: THE START OF NONLINEAR ITERATIVE TIME INTEGERATION (for U, V & T)
! Original version : NLITER.f (IAP 21L)
! Reconstructed & supplemented : ZhangHe
! Completed : 2006.2.22
! Update : 2007.5.7, ZhangHe, change subroutine's name from NLITER3 to nliter_uvt
!          2008.6.10, WuJianping & ZhangHe, parallel version
!          October 2012, Jiang Jinrong, for 2D parallel
!------------------------------------------------------------------------------------------------
   use shr_kind_mod, only: r8 => shr_kind_r8
   use IAP_grid,  only : NX, NY, NL, IB, IE, JB, JE, period
   use IAP_prog,  only : UT, VT, TT 
   use tendency,  only : DU, DV, DT 
   use mathconst, only : ZERO
   use pmgrid,    only : beglatdyn, endlatdyn, beglev, endlev
   use perf_mod,  only : t_startf, t_stopf       ! zhh  2013-02-06

   implicit none
!------------------------------------Arguments---------------------------------------------------
   real(r8), intent(in) :: CCU(NX,beglev:endlev,beglatdyn:endlatdyn)   ! input prognostic variables for iteration,  UT
   real(r8), intent(in) :: CCV(NX,beglev:endlev,beglatdyn:endlatdyn)   ! input prognostic variables for iteration,  VT
   real(r8), intent(in) :: CCT(NX,beglev:endlev,beglatdyn:endlatdyn)   ! input prognostic variables for iteration,  TT
   real(r8), intent(in) :: CT              ! input dynamical integration timestep
!----------------------------------Local workspace-----------------------------------------------
   integer  :: I,J,K             ! loop index
!------------------------------------------------------------------------------------------------

   call t_startf('nliter_uvt')
   DO J = beglatdyn, endlatdyn
      IF (J.GE.JB.AND.J.LE.JE) THEN
         DO K = beglev ,endlev
            DO I = IB,IE
               IF (DU(I,K,J).GT.1E5.OR.CT.GT.4E3) THEN
                  PRINT*,'NLITER-----DU,DV,DT:', DU(I,K,J), DV(I,K,J), DT(I,K,J)
!!                  PAUSE 'NLITER---1'
               ENDIF
!**************************** integrate onwards one step *************************************
               UT(I,K,J) = CCU(I,K,J)+CT*DU(I,K,J)  ! CCU = UT
               VT(I,K,J) = CCV(I,K,J)+CT*DV(I,K,J)  ! CCV = VT
               TT(I,K,J) = CCT(I,K,J)+CT*DT(I,K,J)  ! CCT = TT
!********************************************************************************************
            ENDDO
            call period( UT(1,K,J) )      ! for spherical cyclicity
            call period( VT(1,K,J) )
            call period( TT(1,K,J) )
         ENDDO
!
      ELSE IF(J.EQ.1) THEN  ! AT the north polar
         DO K = beglev ,endlev
            DO I = IB, IE
               VT(I,K,J) = CCV(I,K,J) +CT*DV(I,K,J)
               TT(I,K,J) = CCT(IB,K,J)+CT*DT(IB,K,J)
               UT(I,K,J) = ZERO
            ENDDO
            call period( UT(1,K,J) )      ! for spherical cyclicity
            call period( VT(1,K,J) )
            call period( TT(1,K,J) )
         ENDDO
!
      ELSE                 ! AT the south polar
         DO K = beglev ,endlev
            DO I = 1 ,NX
               UT(I,K,J) = ZERO
               VT(I,K,J) = ZERO
               TT(I,K,J) = CCT(IB,K,J)+CT*DT(IB,K,J)
            ENDDO
         ENDDO
      ENDIF
   ENDDO
!
     call t_stopf('nliter_uvt')
   RETURN
END    
      
 
