module Filt
!------------------------------------------------------------------------------------------------
! Purpose: FILTER FOR ALL TENDENCES : LATITUDINAL
! Original version : FILTER.f (IAP 9L)
! Reconstructed to module & revised by : ZhangHe
! Completed : 2006.01.12
! Update : 2006.12.13, ZhangHe, set ADEF & KDEF by linear interpolation	  
!          2008.5, WuJianping, parallel version
!          2010.08, delete r16, Juanxiong he
! Review: ZhangHe, 2011-11-18
!------------------------------------------------------------------------------------------------
   use shr_kind_mod, only: r8 => shr_kind_r8  
   use IAP_grid,  only: NX, NY, IB, IE, JB, JE, IM, NX1, period, periodp
   use mathconst, only: ZERO, HALF, ONE
   use pmgrid,    only: beglatdyn, endlatdyn

   implicit none
!zhh 2008.6.9   include 'mpif.h'
   save
   public

   integer,  parameter :: N   = IM / 2
   integer,  parameter :: JH  = (NY+1) / 2        ! zhh
   integer,  parameter :: IS  = IB - 1        ! zhh
   integer  :: LB(4), LH(4), LM(4)
   integer  :: NCUT(JH,4)         ! critical wave number at latitude J    
   integer  :: KN(JH,4)
   real(r8) :: SF(N,JH,4)         ! damping factor of amplitude
   real(r8) :: ALPHA(JH,4)        ! smoothing factor of successive three point smoothing
   
!================================================================================================
CONTAINS
!================================================================================================

!================================================================================================
   SUBROUTINE FILTER(CH,ID)
!------------------------------------------------------------------------------------------------
! Purpose: FILTER FOR ALL TENDENCES : LATITUDINAL
!------------------------------------------------------------------------------------------------      
      use smoother, only: SOSS1D
      use iap_fft,  only: N0, TRIGS0, IFAX0, IAP_FFT99

      implicit none
!------------------------------------Arguments---------------------------------------------------
      real(r8), intent(inout) :: CH(NX,NY)      ! input variable needed to be filtered
!!      real(r8), intent(inout) :: CH(NX,beglatdyn:endlatdyn)   ! input variable needed to be filtered
      integer , intent(in)    :: ID             ! index
!----------------------------------Local workspace-----------------------------------------------
      integer, parameter :: INC  = 1
      integer, parameter :: JUMP = N0 + 2
      integer, parameter :: LOT  = 1
      integer  :: ISIGN 
      real(r8) :: WORK( (N0+1)*LOT )
!=============================== zhh =====================================
      real*8   :: A( (N0+2)*LOT )        !zhh 2007.11.10
      COMPLEX*16 ::  AC( N0/2+1 )        !zhh 2007.11.10
!============================ 2007.8.5 ====================================
      EQUIVALENCE ( AC(1),A(1) )

      real(r8), parameter :: FIM = IM
      real(r8) :: PUP(NX)
      real(r8) :: ALPHAJ, PUM, X0
!      integer  :: JBH    ! begin index of latitude to do FFT + Arakawa.damping 
      integer  :: JBH, JEH, JBM, JEM, JBL, JEL, JFP, NCT, MJ, IP, JBS, JES, KNJ, II
      integer  :: I, J, K, NP, KK          ! loop index
!------------------------------------------------------------------------------------------------

!     LU.FFT + ARAKAWA.DAMPING       NEAR POLES
      JBH  = LB(ID)
      JEH  = LH(ID) + 1
      JFP  = JE  + JBH
        
!!	  print*, 'JBH =', JBH
!!	  print*, 'JEH =', JEH
!!	  print*, 'JFP =', JFP
!!	  pause 'FILTER 1'

      DO J = JBH, JEH
         NCT        = NCUT(J,ID)
         MJ         = JFP - J       ! the value of MJ's southern latitude equals to J's northern latitude 

!zhh  do FFT + Arakawa.damping at north pole
         if (beglatdyn.le.J.and.J.le.endlatdyn) then
            DO I = 1, N0 + 2             !zhh  
               A(I)    = CH(I,J)
            END DO
!------------------------------------------------------------------------------------------------
            ISIGN = -1
            CALL IAP_FFT99(A,WORK,TRIGS0,IFAX0,INC,JUMP,N0,LOT,ISIGN)
!------------------------------------------------------------------------------------------------
!  do Arakawa.damping when wave number NP > NCT (critical wave number)
            DO NP  = NCT,N
               AC(NP)  = AC(NP) * SF(NP,J,ID)
            END DO
            A(N0+1) = ZERO
            A(N0+2) = ZERO
!------------------------------------------------------------------------------------------------
            ISIGN = 1
            CALL IAP_FFT99(A,WORK,TRIGS0,IFAX0,INC,JUMP,N0,LOT,ISIGN)
!------------------------------------------------------------------------------------------------
            DO I = 1, N0 + 2            !zhh
               CH(I,J) = A(I)         ! EQUIVALENCE ( AC(1),A(1) )
            END DO
            call periodp( CH(1,J) )
         endif
!****************************************

!zhh  do FFT + Arakawa.damping near south pole
         if (beglatdyn.le.MJ.and.MJ.le.endlatdyn) then
            DO I = 1, N0 + 2            !zhh
               A(I)    = CH(I,MJ)
            END DO
!------------------------------------------------------------------------------------------------
            ISIGN = -1
            CALL IAP_FFT99(A,WORK,TRIGS0,IFAX0,INC,JUMP,N0,LOT,ISIGN)
!------------------------------------------------------------------------------------------------
            DO NP = NCT,N
               AC(NP)  = AC(NP) * SF(NP,J,ID)
            END DO
            A(N0+1) = ZERO
            A(N0+2) = ZERO
!------------------------------------------------------------------------------------------------
            ISIGN = 1
            CALL IAP_FFT99(A,WORK,TRIGS0,IFAX0,INC,JUMP,N0,LOT,ISIGN)
!------------------------------------------------------------------------------------------------
            DO I = 1, N0 + 2                !zhh
               CH(I,MJ)= A(I)
            END DO
!------------------------------------------------------------------------------------------------
            call periodp( CH(1,MJ) )
!------------------------------------------------------------------------------------------------
         endif
      END DO
!
! **********  SECOND-ORDER SHAPIRO SMOOTHER IF REQUESTED in middle & low latitudes ************
!
      IF ( LM(ID).LE.-1 ) THEN
         IP     = - LM(ID)
         JBS    = JEH + 1
         JES    = JFP - JBS
         CALL SOSS1D( CH,JBS,JES,IP )
! --------------------------------------------
         RETURN
      ENDIF
!
!     SUCCESSIVE THREEE POINT SMOOTHING IN MIDDLE LATITUDES
!
      JBM = JEH + 1
      JEM = JEH + LM(ID) + 1   ! zhh  2006.12.13
!!	  print*, 'JBM =', JBM
!!	  print*, 'JEM =', JEM
      IF ( JEM.GE.JBM ) THEN
         DO J = JBM,JEM
            MJ     = JFP - J
            KNJ    = KN(J,ID)
            ALPHAJ = ALPHA(J,ID)

! zhh  do THREEE POINT SMOOTHING IN north MIDDLE LATITUDES
            DO KK  = 1  ,KNJ
               if (beglatdyn.le.J.and.J.le.endlatdyn) then
                  DO I = IB  ,NX1
                     PUP(I)  = CH(I-1,J)  - CH(I,J)
                  END DO
                  DO I = IB ,IE
                     PUM     = PUP(I)     - PUP(I+1)
                     CH(I,J) = CH(I ,J)   + ALPHAJ*PUM
                  END DO
                  call period( CH(1,J) )
               endif
               IF ( MJ.EQ.J ) GOTO 10
!
! zhh  do THREEE POINT SMOOTHING IN south MIDDLE LATITUDES
               if (beglatdyn.le.MJ.and.MJ.le.endlatdyn) then
                  DO I = IB  ,NX1
                     PUP(I)   = CH(I-1,MJ) - CH(I,MJ)
                  END DO
                  DO I = IB ,IE
                     PUM      = PUP(I)     - PUP(I+1)
                     CH(I,MJ) = CH(I ,MJ)  + ALPHAJ*PUM
                  END DO
                  call period( CH(1,MJ) )
               endif
10          END DO
         END DO
      ENDIF
!=============================== zhh =====================================
!!      print*, 'after THREEE POINT SMOOTHING  ', 'CH(4,2) =', CH(4,2)
!============================ 2007.8.8 ====================================
!*******************************************************************************************
!     FFT 2-GRID WAVE IN LOWER LATITUDES
!
      JBL = JEM + 1
      JEL = JFP - JBL
!!	  print*, 'JBL =', JBL
!!	  print*, 'JEL =', JEL

      IF ( JEL.GE.JBL ) THEN
         DO J = max(JBL,beglatdyn), min(JEL,endlatdyn)
            X0 = ZERO
            DO I = 1  ,IM ,2
               X0 = (CH(I,J) - CH(I+1,J)) + X0
            END DO
            X0 = X0 / FIM
            DO I = 1 ,IM ,2
               II      = I + 1
               CH(I,J) = CH(I ,J) - X0
               CH(II,J)= CH(II,J) + X0
            END DO
            call periodp( CH(1,J) )
         END DO
      END IF 
!********************************************************************************************
!=============================== zhh =====================================
!!      print*, 'at the end of FILTER  ', 'CH(4,2) =', CH(4,2)
!============================ 2007.8.8 ====================================
      RETURN
   END SUBROUTINE

!================================================================================================
   SUBROUTINE CONFIL(IDSOSS,DLAT,DLON)
!------------------------------------------------------------------------------------------------
!     SET CONSTANTS USED IN SUB.FILTER
!------------------------------------------------------------------------------------------------
      use iap_fft, only: IAP_SET99

      implicit none
!------------------------------------Arguments---------------------------------------------------
      integer,  intent(in)  :: IDSOSS          ! switch of filter
      real(r8), intent(in)  :: DLAT
      real(r8), intent(in)  :: DLON
!----------------------------------Local workspace----------------------------------------------
      real(r8), parameter :: FJM = JE
      integer  :: KDEF(JH)
      real(r8) :: ASN(N), DEFF(4), SNL(NY,2), ADEF(JH)
      real(r8) :: DEGL, WLCPP, WLCMP, WLCPV, WLCMV, SLCPP
      real(r8) :: SLCMP, SLCPV, SLCMV, DLNT, HALFDN, HALFDT, WEAK
      real(r8) :: ALTV, ALTP, DEFK, DLNTSA, SFNJ, SNL0, SNL1
!================================ added by zhh 2006.12.13 ==================================
      real(r8) :: lat(9), rlat(9), Ar(8), deglat(JH)
!===========================================================================================
      integer, parameter :: NS = 1 
      integer  :: JEH, JBM
      integer  :: J, K, L, jr          ! loop index
!================================ added by zhh 2006.12.13 ==================================
      data lat  / 70.0, 66.0, 62.0, 58.0, 54.0, 50.0, 46.0, 42.0, 38.0 /
      data rlat / 2.3955E-1, 1.433E-1, 1.076E-1, 8.441E-2, 6.632E-2, 4.913E-2,   &
	               3.292E-2, 1.777E-2, 3.747E-3 /
!===========================================================================================
!------------------------------------ for 4x5 ------------------------------------------------
!      DATA ( KDEF(J),J=06,JH )  / 5,7*1,10*1 /
!      DATA ( ADEF(J),J=06,JH )  / 4.791E-2, 1.433E-1, 1.076E-1, 8.441E-2, 6.632E-2,  &
!                                  4.913E-2, 3.292E-2, 1.777E-2, 10*3.747E-3 /        
!----------------------------------------------------------------------------------------------
!     
!-------------------------- SET CONSTANTS FOR FFT USAGE ---------------------------------------
      CALL IAP_SET99
!----------------------------------------------------------------------------------------------
!     SET DEFAULT VALUE
!=============================== zhh =====================================
      DO J = 1, JH                           
         ADEF(J) = ZERO
         KDEF(J) = 0
      END DO
!============================ 2007.8.6 ====================================
      DO K = 1 ,4
         DO J = 1 ,JH
            KN(J,K)    = 0
            ALPHA(J,K) = ZERO
            NCUT(J,K)  = N
            DO L = 1 ,N
               SF(L,J,K)  = ONE
            END DO
         END DO
      END DO
      DEGL       = 180.0E0 / FJM
!     LB         : INDEX  OF THE FIRST LATITUDE UNDER HIGH-FFT
!     LH         : NUMBER OF THE HIGH-LATITUDES UNDER HIGH-FFT
!     LM         : NUMBER OF THE MIDL-LATITUDES UNDER MIDL-3PS
      LB(1)      = 2    ! for P-grid (U,T,P,Q)
      LB(2)      = 2
      LB(3)      = 1    ! for V-grid (V)
      LB(4)      = 1
!     FOR THE STANDARD VERSION:      WEAK=1.0E0
!             ARAKAWA DAMPING        FOR POLAR AREA    [90,70)
!             RECURSIVE SMOOTHING    FOR MIDDLE LAT    [70,30)
!             FFT 2-GRID WAVE        FOR LOWER  LAT    [30,00]
!!      WLCPP      = 90.0E0 - 66.0E0
!!      WLCMP      = 66.0E0 - 30.0E0
!!      WLCPV      = 88.0E0 - 68.0E0
!!      WLCMV      = 68.0E0 - 32.0E0
      WLCPP      = 90.0E0 - 70.0E0           ! zhh   for P-grid
!!      WLCMP      = 70.0E0 - 30.0E0           ! zhh
      WLCMP      = 70.0E0 - 38.0E0           ! zhh  2006.12.13
      WLCPV      = WLCPP                     ! zhh   for V-grid
      WLCMV      = WLCMP                     ! zhh 
      LH(1)      = WLCPP / DEGL + 0.1E0
      LM(1)      = WLCMP / DEGL + 0.1E0
      LH(3)      = WLCPV / DEGL + 0.1E0 - 1  ! zhh
      LM(3)      = WLCMV / DEGL + 0.1E0
!     FOR THE STONG VERSION:
!             ARAKAWA DAMPING        FOR POLAR AREA    [90,32)
!             RECURSIVE SMOOTHING    FOR MIDDLE LAT    [32,00]
      SLCPP      = 90.0E0 - 34.0E0
      SLCMP      = 34.0E0 - 02.0E0
      SLCPV      = 88.0E0 - 32.0E0
      SLCMV      = 36.0E0 + 00.0E0
      LH(2)      = SLCPP  / DEGL   + 0.1E0
      LM(2)      = SLCMP  / DEGL   + 0.1E0
      LH(4)      = SLCPV  / DEGL   + 0.1E0
      LM(4)      = SLCMV  / DEGL   + 0.1E0
!
!     ARAKAWA DAMPING METHOD NEAR POLES
!
      DLNT       = DLON   / DLAT        ! DLNT   = d(lamda) / d(sita)
      HALFDN     = HALF   * DLON        ! HALFDN = d(lamda) / 2
      HALFDT     = HALF   * DLAT        ! HALFDT = d(sita) / 2
!     WEAK   HIGH-FFT
      WEAK       = 1.0E0
      DEFF(1)    = DLNT   * MAX( ONE , WEAK )
      DEFF(3)    = DEFF(1)
!     STRONG HIGH-FFT
      DEFF(2)    = DLNT
      DEFF(4)    = DEFF(2)
!
      DO L = 1 ,N                       ! N = IM / 2 
         ASN(L)  = ONE / SIN( HALFDN*DBLE( L ) )  ! ASN(L) = 1/sin(L*d(lamda)/2) = 1/sin(L*PI/IM)
      END DO
      ALTV       = HALFDT
      SNL(1,1)   = ZERO
      SNL(1,2)   = SIN( ALTV )
      DO J = JB ,JE
         ALTP    = DLAT   * DBLE(J-1)
         ALTV    = ALTP   + HALFDT
         SNL(J,1)= SIN( ALTP )
         SNL(J,2)= SIN( ALTV )
      END DO
      SNL(NY,1)  = ZERO

!     FOR P-GRID FIELDS
      DO K = 1 ,2
         DEFK          = DEFF(K)   
         DO J = 1 ,JH
            DLNTSA     = DEFK   * SNL(J,1)
            DO L = 1 ,N
               SFNJ    = ( DLNTSA * ASN(L) ) ** NS
               IF( SFNJ.LT.ONE ) SF(L,J,K) = SFNJ  ! else SF(L,J,K) = 1
!   SFNJ    = [ D(lamda)*sin(sita(j)) ] / [ D(sita)*sin(L*PI/IM) ]
!   SF(L,j) = Min[1, SFNJ] , L is the zonal wave number
            END DO
         END DO
      END DO

!     FOR V-GRID FIELDS
      DO K = 3 ,4
         DEFK          = DEFF(K)
         DO J = 1 ,JH
            DLNTSA     = DEFK   * SNL(J,2)
            DO L = 1 ,N
               SFNJ    = ( DLNTSA * ASN(L) ) ** NS
               IF( SFNJ.LT.ONE ) SF(L,J,K) = SFNJ
            END DO
         END DO
      END DO
      DO K = 1 ,4
         DO J = 1 ,JH
            DO L = 1 ,N
               IF ( SF(L,J,K).LT.ONE ) THEN    ! => SF(L+n,J,K) < 1, n >= 1
                  NCUT(J,K)= L     ! L is critical wave number at latitude j
                  GOTO 40
               ENDIF
            END DO
40       END DO
      END DO
!
!     RECURSIVE SMOOTHING OVER MIDDLE LATITUDES
!
      JEH           = LH(1) + 1                    !zhh
      JBM           = JEH + 1                      !zhh    
!
      do jr = 1, 8
         Ar(jr) = ( rlat(jr)-rlat(jr+1) ) / ( lat(jr)-lat(jr+1) )
      end do

!   set ADEF	  
      do j = JEH, JH
         deglat(j) = 90.0 - DBLE(j-1)*DEGL
         do jr = 1, 8		 
            if ( deglat(j) <= lat(jr) .and. deglat(j) >= lat(jr+1) ) then
               ADEF(j) = rlat(jr+1) + Ar(jr)*( deglat(j)-lat(jr+1) )
               exit
            else if ( deglat(j) < lat(9) .and. deglat(j) > 0.01 ) then
               ADEF(j) = 0.37470E-02			
               exit
            end if      
         end do
!   set KDEF	  
         if ( deglat(j) < lat(1) .and. deglat(j) >= lat(9) ) then
            KDEF(j) = 1	  
         else if ( deglat(j) < lat(9) .and. deglat(j) > 0.01 ) then
            KDEF(j) = 1	  
         else
            KDEF(j) = 0         
         end if
	  end do
      KDEF(JEH) = 5
      KDEF(JH)  = 1
      ADEF(JEH) = ADEF(JEH) / KDEF(JEH)
      ADEF(JH)  = 0.37470E-02


!     FOR P-GRID FIELDS
      DO J = JEH, JH
         KN(J,1)    = KDEF(J)
         ALPHA(J,1) = ADEF(J)
      END DO
!     FOR V-GRID FIELDS   [JEH, JH]
!!      J          = 06
      KN(JEH,2)     = MAX0( KN(JEH,1)-1, 2 )
      ALPHA(JEH,2)  = ALPHA(JEH,1)
      SNL0          = SNL(JEH,2)**2
      DO J = JBM, JH-1                            !zhh
         SNL1       = SNL(J,2)**2
         KN(J,2)    = KN(J,1)
         ALPHA(J,2) = (SNL0*ALPHA(J,1)+SNL1*ALPHA(J+1,1))/(SNL0+SNL1)
         SNL0       = SNL1
      END DO
!      J          = JH
      KN(JH,2)    = KN(JH,1)
      ALPHA(JH,2) = ALPHA(JH,1)
      DO K = 3 ,4
         DO J = 1 ,JH
            KN(J,K)    = KN(J,K-2)
            ALPHA(J,K) = ALPHA(J,K-2)
         END DO
      END DO
      IF ( IDSOSS.EQ.+1 ) THEN
!     USE FFT+SHAPIRO_SMOOTHER
         LM(1)      = -2
         LM(2)      = -3
         LM(3)      = -2
         LM(4)      = -3
!     USE PURE FFT IF REQUESTED
      ELSE IF ( IDSOSS.EQ.+2 ) THEN
         LM(1)      = 00
         LM(2)      = 00
         LM(3)      = 00
         LM(4)      = 00
         DO K = 1 ,4
            DO J = JH,1 ,-1
               IF ( NCUT(J,K).LT.N ) THEN
                  LH(K)    = J
                  GOTO 80
               ENDIF
            END DO
80       END DO
      ENDIF

      RETURN
   END SUBROUTINE
END module
