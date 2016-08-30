module smoother
!------------------------------------------------------------------------------------------------
! Purpose: Set the consts for SHAPIRO SMOOTHER
! Original version : SMOOTH.f & SOSS2D.f  (IAP 9L)
! Reconstructed to module : ZhangHe
! Completed : 2005.9.15
! Update: 2007.12.22, ZhangHe, added Qliq, Qice in sub. smooth
!         2008.4.7, ZhangHe, added caculation of Psa in sub. SMOTHP
!         2008.5, Wujianping, parallel version
!         2008.6.11, ZhangHe, available for both serial & parallel
!         2008.6.14, ZhangHe, added 'Istar' in sub. SMOTHP
! Reviewed: ZhangHe, 2011-11-19
!           Jiang Jinrong, 2012 October, for 2D parallel
! Modified: ZhangHe, 2013-01-24, Included dyn_state
!           ZhangHe, 2013-03-21, removed smoothing Q
!------------------------------------------------------------------------------------------------
   use shr_kind_mod, only: r8 => shr_kind_r8
   use IAP_grid,   only: NX, NY, NL, IB, IE, IM, JB, JE, periodp, period
   use Dyn_const,  only: PMTOP, SIGL
   use stdatm,     only: PSB, TSB   
   use pmgrid,     only: beglatdyn, endlatdyn, beglatdynex, endlatdynex, &
                         loc_JB, loc_JE,beglev,endlev,npr_y
   use spmd_utils, only: iam, masterproc
#if ( defined SPMD )
  use mod_comm, only: mp_send3d, mp_recv3d
#endif

   implicit none

   save
   public
 
   real(r8) :: SMONE(3,3)          ! coefficient of 1D 2-ORDER SHAPIRO SMOOTHER
   real(r8) :: SMTWO(6,3)          ! coefficient of 2D 2-ORDER SHAPIRO SMOOTHER
   real(r8) :: TRSM(NX,NL,NY)      ! temperature of standard atmosphere 

!================================================================================================
CONTAINS
!================================================================================================

!================================================================================================
   SUBROUTINE SMOOTH
!------------------------------------------------------------------------------------------------
!  Popurse: SMOOTHING T & Q  BY 2-ORDER SHAPIRO SMOOTHER(LON & LAT)
!                     U & V  BY 2-ORDER SHAPIRO SMOOTHER(LON  ONLY)
!------------------------------------------------------------------------------------------------
      use IAP_prog,   only: U, V, T, Q
      implicit none
!-------------------------------------Local workspace-------------------------------------------
      integer, parameter :: IC = 1   ! index controlling smoothing intensity
      real(r8) :: WW(NX,NY)
      real(r8) :: WG(NX,NL,NY)
      integer  :: I, J, K
!------------------------------------------------------------------------------------------------
!
!     SMOOTHING AIR TEMPERATURE BY USING ITS DEPARTURE FROM MSA
!
!     OTHER MODEL LAYERS
      DO J = beglatdyn, endlatdyn
         DO K = beglev ,endlev
            DO I = 1 ,NX
               WG(I,K,J) =  T(I,K,J) - TRSM(I,K,J)  
            END DO
		 END DO
      END DO
!
      DO K = beglev ,endlev
!-------------------------SECOND-ORDER SHAPIRO SMOOTHER OF TWO DIMENSION ------------------------
         WW(:,beglatdyn:endlatdyn)=WG(:,K,beglatdyn:endlatdyn)
         CALL SOSS2D( WW,0,1,IC )
         WG(:,K,beglatdyn:endlatdyn)=WW(:,beglatdyn:endlatdyn)
!------------------------------------------------------------------------------------------------
      END DO
      DO J = beglatdyn, endlatdyn
         DO K = beglev ,endlev
            DO I = 1 ,NX
               T(I,K,J) = WG(I,K,J) + TRSM(I,K,J)
            END DO
		 END DO
      END DO
!
!     SMOOTHING MIXING RATIO
!     
      DO K = beglev , endlev
!------------------------------------------------------------------------------------------------
!!         WW(:,beglatdyn:endlatdyn)=Q(:,K,beglatdyn:endlatdyn)
!!         CALL SOSS2D( WW,0,1,IC )
!!         Q(:,K,beglatdyn:endlatdyn)=WW(:,beglatdyn:endlatdyn)
!         CALL SOSS2D( Qliq(1,1,K),0,1,IC )    !zhh 2007.12.22
!        CALL SOSS2D( Qice(1,1,K),0,1,IC )    !zhh 2007.12.22
!------------------------------------------------------------------------------------------------
      END DO
!
!     SMOOTHING WIND
!     
      DO K = beglev ,endlev
!----------------SECOND-ORDER SHAPIRO SMOOTHER OF ONE DIMENSION  ( LON  ONLY )-------------------
         WW(:,beglatdyn:endlatdyn)=U(:,K,beglatdyn:endlatdyn)
         CALL SOSS1D( WW,JB,JE,IC )
         U(:,K,beglatdyn:endlatdyn)=WW(:,beglatdyn:endlatdyn)
         WW(:,beglatdyn:endlatdyn)=V(:,K,beglatdyn:endlatdyn)
         CALL SOSS1D( WW,01,JE,IC )
         V(:,K,beglatdyn:endlatdyn)=WW(:,beglatdyn:endlatdyn)
!------------------------------------------------------------------------------------------------
      END DO
      RETURN
   END SUBROUTINE

!================================================================================================
   SUBROUTINE SMOTHP
!------------------------------------------------------------------------------------------------
!   SMOOTHING P & PT BY 2-ORDER SHAPIRO SMOOTHER(LON & LAT)
!------------------------------------------------------------------------------------------------
      use IAP_prog, only: P, PT, Pstar1, Psa
	  use stdatm,   only: P00
      use flexib,   only: Istar        !zhh 2008.6.14
	  implicit none
!-------------------------------------Local workspace-------------------------------------------
      integer, parameter :: IC = 1   ! index controlling smoothing intensity
      real(r8) :: WK(NX,NY)
	  real(r8) :: WW(NX,NY)
	  real(r8) :: WP(NX,NY)  !zhh 2008.6.17
	  real(r8) :: PBI, WWI
	  integer  :: I, J
!------------------------------------------------------------------------------------------------
!
      DO J = beglatdyn, endlatdyn
         DO I = 1 ,NX
            WK(I,J) = PSB(I,J) - PMTOP
            WW(I,J) = P  (I,J) - WK(I,J)
         END DO
      END DO
!------------------------------------------------------------------------------------------------
      CALL SOSS2D( WW,0,1,IC )
	  if ( Istar /= 1 ) then
         WP(:,beglatdyn:endlatdyn)=Pstar1(:,beglatdyn:endlatdyn)
         call SOSS2D( WP,0,1,IC )     !zhh 2008.6.14 
         Pstar1(:,beglatdyn:endlatdyn)=WP(:,beglatdyn:endlatdyn)
         end if
!------------------------------------------------------------------------------------------------
      DO J = beglatdyn, endlatdyn
         DO I = 1 ,NX
            WWI      = WW  (I,J)
            P (I,J)  = WK  (I,J) + WWI
            if ( Istar == 1 ) Pstar1(I,J) = P(I,J)   !zhh 2008.6.14 
            PT(I,J)  = sqrt( Pstar1(I,J)/P00 )
            Psa(I,J) = P(I,J) + PMTOP - PSB(I,J)     !zhh 2008.3.19
         END DO
      END DO
      RETURN
   END SUBROUTINE

!================================================================================================
   SUBROUTINE STSMOT
!------------------------------------------------------------------------------------------------
!   SET CONSTANTS USED IN THE SMOOTHER
!------------------------------------------------------------------------------------------------
      use stdatm, only: TMSA
      implicit none
!-------------------------------------Local workspace-------------------------------------------
      real(r8) :: TW(NL)          ! temperature of standard atmosphere
      real(r8) :: PES, PWK, TWK
      integer  :: I, J, K
!------------------------------------------------------------------------------------------------
!
!--------------------------- SET CONSTANTS FOR SUB.SOSS2D & SOSS1D ------------------------------
      CALL STSOSS    
!------------------------------------------------------------------------------------------------
!      DO J = loc_JB,loc_JE
      DO J = beglatdyn, endlatdyn
         if ( j >= JB .and. j <= JE ) then 
            DO K = 1 ,NL
               DO I = IB,IE
                  PES         = PSB(I,J) - PMTOP
                  PWK         = PES*SIGL(K) + PMTOP  ! PWK: pressure of standard atmosphere
                  TRSM(I,K,J) = TMSA( PWK )   
               END DO
		    END DO
         else   ! at north and south pole         
            PES      = PSB(IB,J)   - PMTOP
            DO K = 1 ,NL
               PWK   = PES*SIGL(K) + PMTOP
               TW(K) = TMSA( PWK )
            END DO
            DO K = 1 ,NL
               TWK   = TW(K)
               DO I = IB,IE
                  TRSM(I,K,J) = TWK
               END DO
            END DO
         end if
!
         DO K = 1 ,NL
            call period( TRSM(1,K,J) )
         END DO
	  END DO
!
      RETURN
   END SUBROUTINE

!================================================================================================
   SUBROUTINE SOSS2D(A,IV,IP,IC)
!-----------------------------------------------------------------------------------------------
!  Purpose : SECOND-ORDER SHAPIRO SMOOTHER OF TWO DIMENSION  ( LON & LAT )
!  Reference : 
!  [1]. R. Shapiro,1970: Smoothing, filtering and boundary effects, 
!       Rev. Geophys. and Space Phys. (8), 359-387
!  [2]. X.-Z. Liang,1986: The Design of IAP GCM and the Simulation of Climate
!       and Its Interseasonal Variability, Ph.D. Thesis  250pp
!-----------------------------------------------------------------------------------------------
      use mathconst,  only: ZERO
      use dynamics_vars,   only: T_FVDYCORE_STATE
      use dyn_internal_state,   only : get_dyn_state      !zhh 2013-01-24
      
	  implicit none 
!------------------------------------Arguments------------------------------------------------
      integer, intent(in)     :: IV        ! IV = 0, A is a scalar; IV = 1, A is a vector
      integer, intent(in)     :: IP        ! index controlling boundary effect
      integer, intent(in)     :: IC        ! index controlling smoothing intensity
      real(r8), intent(inout) :: A(NX,NY)  ! the variable need to be smoothed 
!-------------------------------------Local workspace-------------------------------------------
      integer,  parameter :: IG  = IM + 4
      integer,  parameter :: JG  = NY + 4
      integer,  parameter :: IMH = IM / 2
      integer,  parameter :: NEX = IM + 2
      real(r8), parameter :: FIM = IM 
      real(r8) :: B(IG,JG), C(IG,JG)
      real(r8) :: D(IG)
      real(r8) :: F00, F01, F02, F11, F22, F12, BIJ01, BIJ02, BIJ11, BIJ22, BIJ12, SMA
      type (T_FVDYCORE_STATE), pointer :: dyn_state   !zhh 2013-01-24
      integer  :: JF, JF2, JQ, JJ, J2, I2, JN, JS, M, II, IX, J1, JJ1, JJ2, I1, II1, II2
      integer  :: I, J, J0              ! loop index
      integer  :: ibeg, iend, dest, src
      integer  :: commyz     !zhh 2013-01-24
!-----------------------------------------------------------------------------------------------
!
      dyn_state => get_dyn_state()
      commyz = dyn_state%grid%commyz
!
!     EXTEND THE FIELD DUE TO THE BOUNDARY EFFECT ON SHAPIRO BY THE SPHERICAL CYCLICITY
!         
      JF   = NY+ 1 - IP
      JF2  = JF+ 2
      JQ   = JF+ 4
      JJ   = 5 - IP
      DO J = beglatdynex ,min(JF,endlatdynex)        !zhh
         J2          = J + 2
         DO I = 1 ,IM
            I2       = I + 2
            B(I2,J2) = A(I,J)
         END DO
      END DO
!
!     NEAR NORTH/SOUTH POLES
      DO J = 1 ,2
         JN           = JJ  - J
         JS           = NY  - J
         M            = JF2 + J
         DO I = 1 ,IM
            I2         = I  + 2
            IF ( I.LE.IMH ) THEN
               II      = I  + IMH
            ELSE
               II      = I  - IMH
            ENDIF
            IF ( IV.EQ.1 ) THEN
!       FOR A VECTOR FIELD
               B(I2,J) = - A(II,JN)
               B(I2,M) = - A(II,JS)
            ELSE 
!       FOR A SCALAR FIELD
               B(I2,J) =   A(II,JN)
               B(I2,M) =   A(II,JS)
            ENDIF
         END DO
      END DO
!     LONGITUDINAL B.C. ( Boundary Condition )
      DO I = 1 ,2
         II         = I + IM
         IX         = I + NEX           !! zhh
         I2         = I + 2
         DO J = beglatdyn ,min(JQ,endlatdyn+2)
            B(I ,J) = B(II,J)
            B(IX,J) = B(I2,J)
         END DO
      END DO
!
      F00        = SMTWO(1,IC)
      F01        = SMTWO(2,IC)
      F02        = SMTWO(3,IC)
      F11        = SMTWO(4,IC)
      F22        = SMTWO(5,IC)
      F12        = SMTWO(6,IC)
!
!     SMOOTHER FOR INNER POINTS
!      
#if (defined SPMD)
!       iam=2               iam=1            iam=0                     !wjp 2007.05
!   01 02 03 04 05 06 | 07 08 09 10 11 | 12 13 14 15 16            (A) !wjp 2007.05
!   01 02 03 04 05 06   07 08(iam=2)                               (B) !wjp 2007.05
!                      (07 08)09 10 11   12 13(iam=1)              (B) !wjp 2007.05
!                                       (12 13)14 15 16 17 18(iam=0,B) !wjp 2007.05
!!jjr      call mpi_move_left(B(1,endlatdyn-1), B(1,beglatdyn-2), 4*IG)
      src = iam+1
      dest  = iam-1
      if ( mod(iam,npr_y) == 0 ) dest = -1
      if ( mod(iam+1,npr_y) == 0 ) src = -1
      call mp_send3d( commyz, dest, src, IG, JG,1,                &
                      1, IG, 1,JG ,1,1,       &
                      1, IG, endlatdyn-1, endlatdyn+2, 1, 1, B )
      call mp_recv3d( commyz, src, IG,  JG,1,                             &
                      1, IG, 1, JG,1,1,       &
                      1, IG, beglatdyn-2, beglatdyn+1,1,1, B )
!
!   01 02 03 04 05 06   07 08(iam=2)                               (B) !wjp 2007.05
!              <05 06   07 08>09 10 11   12 13(iam=1)              (B) !wjp 2007.05
!                               <10 11   12 13>14 15 16 17 18(iam=0,B) !wjp 2007.05
      ibeg = beglatdyn
      iend = endlatdyn
      if(beglatdyn.eq.1)  ibeg = 3
      if(endlatdyn.eq.NY) iend = JF2
!
      DO J0 = ibeg, iend
#else
      DO J0 = 3 ,JF2
#endif
         J1        = J0 - 1
         J2        = J0 - 2
         JJ1       = J0 + 1
         JJ2       = J0 + 2
         DO I = 3 ,NEX
            I1     = I  - 1
            I2     = I  - 2
            II1    = I  + 1
            II2    = I  + 2
            BIJ01  = B(I1,J0) + B(II1,J0) + B(I ,J1 ) + B(I  ,JJ1)
            BIJ02  = B(I2,J0) + B(II2,J0) + B(I ,J2 ) + B(I  ,JJ2)
            BIJ11  = B(I1,J1) + B(II1,J1) + B(I1,JJ1) + B(II1,JJ1)
            BIJ22  = B(I2,J2) + B(II2,J2) + B(I2,JJ2) + B(II2,JJ2)
            BIJ12  = B(I1,J2) + B(II1,J2) + B(I1,JJ2) + B(II1,JJ2)        &
                   + B(I2,J1)  +B(I2,JJ1) + B(II2,J1) + B(II2,JJ1)
!
            C(I,J0)= F00*B(I,J0) + F01*BIJ01 + F02*BIJ02 + F11*BIJ11 + F22*BIJ22 + F12*BIJ12          
         END DO
      END DO
!
#if (defined SPMD)
!   01 02 03 04 05 06  (07 08)(iam=2)                                  !wjp 2007.05
!                       07 08 09 10 11  (12 13)(iam=1)                 !wjp 2007.05
!                                        12 13 14 15 16 17 18(iam=0)   !wjp 2007.05
!!jjr      call mpi_move_right(C(1,beglatdyn), C(1,endlatdyn+1), 2*IG)
!   01 02 03 04 05 06   07 08(iam=2)                                   !wjp 2007.05
!                       07 08 09 10 11   12 13(iam=1)                  !wjp 2007.05
!                                        12 13 14 15 16 17 18(iam=0)   !wjp 2007.05
!       iam=2                      iam=1               iam=0           !wjp 2007.05
!   01 02 03 04 05 06   07 08 | 09 10 11   12 13 | 14 15 16 17 18(C)   !wjp 2007.05
!   01 02 03 04 05 06 | 07 08   09 10 11 | 12 13   14 15 16      (A)   !wjp 2007.05
      src = iam-1
      dest  = iam+1
      if ( mod(iam,npr_y) == 0 ) src = -1
      if ( mod(iam+1,npr_y) == 0 ) dest = -1
      call mp_send3d( commyz, dest, src,IG ,  JG,1,                  &
                      1, IG, 1, JG,1,1,       &
                      1, IG, beglatdyn, beglatdyn+1,1,1,  C )
      call mp_recv3d( commyz, src, IG,  JG,1,                         &
                      1, IG, 1, JG,1,1,       &
                      1, IG, endlatdyn+1, endlatdyn+2,1,1, C)
!
#endif
!
      DO J = beglatdyn ,min(JF,endlatdyn)
         DO I = 1 ,IM
            A(I,J) = C(I+2,J+2)
         END DO
      END DO
      DO J = beglatdyn ,min(JF,endlatdyn)
         call periodp( A(1,j) )
      END DO
!
!     RESET VALUES AT POLES FOR A SCALAR FIELD
      IF ( IV.NE.1 ) THEN
         DO J = 1, NY, JE
            DO I = 1 ,IM
               D(I) = A(I,J)
            END DO
            SMA     = ZERO
            DO I = 1 ,IM
               SMA  = SMA + D(I)
            END DO
            SMA     = SMA / FIM
            DO I = 1 ,NX
               A(I,J) = SMA
            END DO
         END DO
      ENDIF
!
      RETURN
   END SUBROUTINE 

!================================================================================================
   SUBROUTINE SOSS1D(A,JST,JED,IC)
!-----------------------------------------------------------------------------------------------
!  SECOND-ORDER SHAPIRO SMOOTHER OF ONE DIMENSION  ( LON  ONLY )
!-----------------------------------------------------------------------------------------------
      implicit none 
!------------------------------------Arguments------------------------------------------------
      integer, intent(in)     :: JST       ! start index of latitude
      integer, intent(in)     :: JED       ! end   index of latitude
      integer, intent(in)     :: IC        ! index controlling smoothing intensity
      real(r8), intent(inout) :: A(NX,NY)  ! the variable need to be smoothed 
!-------------------------------------Local workspace-------------------------------------------
      integer,  parameter :: IG  = IM + 4
      integer,  parameter :: JG  = NY + 4
      integer,  parameter :: NEX = IM + 2
      real(r8) :: B(IG,JG), C(IG,JG)
      real(r8) :: F00, F01, F02
      integer  :: I2, II, IX
      integer  :: I, J                     ! loop index
!-----------------------------------------------------------------------------------------------

      DO J = max(JST,beglatdyn), min(JED,endlatdyn)
         DO I = 1, IM
            B(I+2,J) = A(I,J)
         END DO
      END DO
      DO I = 1, 2
         II          = I + IM
         IX          = I + NEX
         I2          = I + 2
         DO J = max(JST,beglatdyn), min(JED,endlatdyn)
            B(I ,J)  = B(II,J)
            B(IX,J)  = B(I2,J)
         END DO
      END DO
!
      F00          = SMONE(1,IC)
      F01          = SMONE(2,IC)
      F02          = SMONE(3,IC)
      DO J = max(JST,beglatdyn), min(JED,endlatdyn)
         DO I = 3 ,NEX
            C(I,J) = F00*B(I,J) + F01*( B(I-1,J)+B(I+1,J) ) + F02*( B(I-2,J)+B(I+2,J) )          
         END DO
      END DO
      DO J = max(JST,beglatdyn), min(JED,endlatdyn)
         DO I = 1 ,IM
            A(I,J) = C(I+2,J)
         END DO
      END DO
      DO J = max(JST,beglatdyn), min(JED,endlatdyn)
         call periodp( A(1,j) )
      END DO
!
      RETURN
   END SUBROUTINE

!================================================================================================
   SUBROUTINE STSOSS
!-----------------------------------------------------------------------------------------------
!  SET CONSTANTS FOR SUB.SOSS2D & SOSS1D
!-----------------------------------------------------------------------------------------------
      use flexib, only: BETA     ! DATA BETA /0.02E0,0.05E0,0.1E0 /
      implicit none 
!-------------------------------------Local workspace-------------------------------------------
      real(r8) :: A38, A16, A332, A256, A964, BETK, BETK2
      integer  :: K        ! loop index
!-----------------------------------------------------------------------------------------------
!     
      A38  = 3.0E0  /  8.0D0
      A16  = 1.0E0  / 16.0D0
      A332 = 3.0E0  / 32.0D0
      A256 = A16 * A16
      A964 = A38 * A38
      DO K = 1 ,3
         BETK       = BETA(K)
!
!     FOR SUB.SOSS1D
         SMONE(1,K) = 1.00D0 - BETK*A38
         SMONE(2,K) = 0.25D0 * BETK
         SMONE(3,K) =        - BETK*A16
!
!     FOR SUB.SOSS2D
         BETK2      = BETK    * BETK
         SMTWO(1,K) = 1.00D0  - 0.75D0*BETK + A964*BETK2
         SMTWO(2,K) = 0.250D0 * BETK - A332*BETK2
         SMTWO(3,K) = -0.25D0 * SMTWO(2,K)
         SMTWO(4,K) = BETK2   * A16
         SMTWO(5,K) = BETK2   * A256
         SMTWO(6,K) = -4.00D0 * SMTWO(5,K)
      END DO
!
      RETURN
  END SUBROUTINE

end module
