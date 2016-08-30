SUBROUTINE nliter_uvtp(CCU,CCV,CCT,CCP,CCP2,CT, Istar)
!------------------------------------------------------------------------------------------------
! Purpose: THE START OF NONLINEAR ITERATIVE TIME INTEGERATION (for U, V, T & Psa)
! Original version : NLITER.f (IAP 21L)
! Reconstructed & supplemented : ZhangHe
! Completed : 2005.9.2
! Update : 2007.5.7, ZhangHe, change subroutine's name from NLITER2 to nliter_uvtp
!          2008.04.10, ZhangHe, add dumb parameter 'Istar' 
!          2008.6.10, WuJianping & ZhangHe, parallel version
!          October 2012, Jiang Jinrong, for 2D parallel
!------------------------------------------------------------------------------------------------
   use shr_kind_mod, only: r8 => shr_kind_r8
   use IAP_grid,   only : NX, NY, NL, IB, IE, JB, JE, period
   use Dyn_const,  only : PMTOP
   use IAP_prog,   only : P, PT, PS2, Psa, Pstar1, UT, VT, TT 
   use tendency,   only : DU, DV, DT, DPsa, DPstar1 
   use stdatm,     only : PSB, P00
   use mathconst,  only : ZERO
   use pmgrid,     only : beglatdyn, endlatdyn, beglev, endlev
   use perf_mod,   only : t_startf, t_stopf       ! zhh  2013-02-06

   implicit none
!------------------------------------Arguments---------------------------------------------------
   real(r8), intent(in) :: CCU(NX,beglev:endlev,beglatdyn:endlatdyn)   ! input prognostic variables for iteration,  UT
   real(r8), intent(in) :: CCV(NX,beglev:endlev,beglatdyn:endlatdyn)   ! input prognostic variables for iteration,  VT
   real(r8), intent(in) :: CCT(NX,beglev:endlev,beglatdyn:endlatdyn)   ! input prognostic variables for iteration,  TT
   real(r8), intent(in) :: CCP(NX,beglatdyn:endlatdyn   )   ! input prognostic variables for iteration,  Psa
   real(r8), intent(in) :: CCP2(NX,beglatdyn:endlatdyn  )   ! input prognostic variables for iteration,  Pstar1  
   real(r8), intent(in) :: CT              ! input dynamical integration timestep
   integer,  intent(in) :: Istar           ! index of the methods to defeine flexible substitute 
!----------------------------------Local workspace-----------------------------------------------
   real(r8) :: sqP0
   integer  :: I,J,K             ! loop index
!------------------------------------------------------------------------------------------------

   call t_startf('nliter_uvtp')
   DO J = beglatdyn, endlatdyn
      IF (J.GE.JB.AND.J.LE.JE) THEN
         DO K = beglev ,endlev
            DO I = IB,IE
               IF (DU(I,K,J).GT.1E5.OR.CT.GT.1E3) THEN
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
10 CONTINUE
   sqP0 = sqrt(P00)
   DO J = beglatdyn, endlatdyn
      IF (J.GE.JB.AND.J.LE.JE) THEN
         DO I = IB,IE
!**************************** integrate onwards one step *************************************
            Psa(I,J)    = CCP (I,J) + CT*DPsa(I,J)
            PS2(I,J)    = Psa(I,J)  + PSB(I,J)
            P(I,J)      = PS2(I,J)  - PMTOP 
!======================= zhh 2008.4.10 ===========================            
            if ( Istar == 1 ) then
               Pstar1(I,J) = P(I,J)
            else
			   Pstar1(I,J) = CCP2(I,J) + CT*DPstar1(I,J) 
            end if
!======================= zhh 2008.4.10 ===========================            
            PT(I,J) = sqrt( Pstar1(I,J) ) / sqP0
         ENDDO
         call period( Psa   (1,J) )
         call period( PS2   (1,J) )
         call period( P     (1,J) )
         call period( Pstar1(1,J) )
         call period( PT    (1,J) )
      ELSE         ! at the north and south polar
         DO I = 1 ,NX
            Psa(I,J)    = CCP (IB,J)   + CT*DPsa(IB,J)
            PS2(I,J)    = Psa(I,J)  + PSB(IB,J)
            P(I,J)      = PS2(I,J)  - PMTOP
!======================= zhh 2008.4.10 ===========================            
            if ( Istar == 1 ) then
               Pstar1(I,J) = P(I,J)
            else
			   Pstar1(I,J) = CCP2(IB,J) + CT*DPstar1(IB,J) 
            end if
!======================= zhh 2008.4.10 ===========================            
            PT(I,J) = sqrt( Pstar1(IB,J) ) / sqP0
         END DO
      END IF
   END DO
!
   call t_stopf('nliter_uvtp')

   RETURN
END    
      
 
