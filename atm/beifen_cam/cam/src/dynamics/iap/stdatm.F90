module stdatm 
!----------------------------------------------------------------------------------
! Purpose: set constants & variables related to standard atmosphere
! Reconstructed to module : ZhangHe
! Original version: STFRAM.f, DIAGBB.f, icaosa.f, SDGHI.f, HMSA.f & TMSA.f (IAP 9L)
! Completed: 2005.8.25
! Update: 2007.4.23, ZhangHe, removed the computation of 'DCB' in sub. DIAGBB
!         2008.4.30, ZhangHe, modified Itypy = 2 in sub. DIAGBB ( CB,DCB => TB )
!         2011.05.03, ZhangHe, let deltac = 0 if adiabatic run
! Reviewed: ZhangHe, 2011-11-19
!------------------------------------------------------------------------------------------------
   use shr_kind_mod, only: r8 => shr_kind_r8
   use IAP_grid,  only: NX, NY, NL, NZ, IB, IE, JB, JE, period
   use physconst, only: RD, GRAV, CAPA, b0
   use mathconst, only: ZERO, HALF, ONE, TWO, THREE, SIX, FOURTH
   use pmgrid,    only: beglatdyn, endlatdyn, loc_JB, loc_JE, beglev, endlev

   IMPLICIT NONE
   save
   public

   real(r8),parameter :: PEALIB = 1160.0       ! the pressure at bottom of standard atmosphere  (unit: hPa)
   real(r8),parameter :: PBALIB = 0.5          ! the pressure at top of standard atmosphere  (unit: hPa)
   real(r8),parameter :: DPALIB = 0.5          ! the difference pressure of two neighbour levels
   real(r8),parameter :: PGALIB = 1060.0       ! the pressure at bottom of actual atmosphere
   real(r8),parameter :: P00    = 1000.0
   integer, parameter :: NA = PEALIB / DPALIB  ! the level number of standard atmosphere
   integer, parameter :: NG = PGALIB / DPALIB  !
!
   real(r8) :: TBB(NA)         ! temperature of standard atmosphere
   real(r8) :: CBB(NA)         ! characteristic phase-speed of standard atmosphere
   real(r8) :: DCBB(NA)        ! DCBB = (1/CBB)*(dCBB/dln(p))
   real(r8) :: HBB(NA)         ! geopotential of standard atmosphere
   real(r8) :: P00SL           ! standard pressure at sea level
   real(r8) :: T00SL           ! standard temperature at sea level
   real(r8) :: TB(NX,NL,NY)    ! temperature of standard atmosphere at model grids
   real(r8) :: CB(NX,NL,NY)    ! characteristic phase-speed of standard atmosphere at model grids
   real(r8) :: DCB(NX,NL,NY)   ! DCB = (1/CB)*(dCB/dln(p))
   real(r8) :: GHI0(NX,NZ,NY)  ! geopotential of standard atmosphere at model grids
   real(r8) :: deltac(NX,NZ,NY)! CB^2 = (b0^2) * (1 + deltac)
   real(r8) :: PSB(NX,NY)      ! surface pressure of standard atmosphere
   real(r8) :: PSB2(NX,NY)     ! PSB2 = SQRT( PSB - Pt )
   real(r8) :: TSB(NX,NY)      ! surface temperature of standard atmosphere
   real(r8) :: H0B(NX,NY)      ! H0B = TSB * RD / PSB 
   real(r8) :: PS0             ! standard pressure at sea level
   real(r8) :: T00             ! standard temperature at sea leve

   private GHS0
   real(r8) :: GHS0 
    
!================================================================================================
CONTAINS
!================================================================================================

!================================================================================================
   SUBROUTINE crestdatm
!------------------------------------------------------------------------------------------------
! Popurse: CREATE LOOK-UP TABLES FOR THE MODEL STANDARD ATMOSPHERE
! Method : SPLINE-FITTING METHOD
!------------------------------------------------------------------------------------------------
      IMPLICIT NONE
!----------------------------------Local workspace-----------------------------------------------
      integer, parameter :: KS = 16       ! the level of reference atmosphere
      integer, parameter :: KC = KS - 1       
      REAL(r8) :: SK(KS)    ! SK = ln(P) 
      REAL(r8) :: TK(KS)    ! TK : temperature of standard atmosphere
      REAL(r8) :: CK(KC,3)  ! CK(KC,1) = d(TK)/d(lnP)
                            ! CK(KC,2) = 0.5 * [d^2(TK)/d(lnP)^2]
                            ! CK(KC,3) = (1/6) * [d^3(TK)/d(lnP)^3]
!
!WB   Add the data above 10hPa From CIRA climate data
! CIRA : Cooperative Institute for Research in the Atmosphere 
!========================================================================
!  Pressure (P) SK (Log(P))  TK (T)     CK(N,1)     CK(N,2)     CK(N,3)
!
!   0.2162828   -1.531168   245.2323    19.87025   0.2870293   -2.276788
!   0.4207161  -0.8657969   257.9098    18.09905   -4.314825   -4.329503
!   0.7959216  -0.2282546   266.5730    5.243720   -12.55026   0.8447471
!    1.492622   0.4005346   265.1181   -10.83125   -11.00837    4.035521
!    2.844819    1.045499   254.6358   -19.11973   -3.131378    2.037444
!    5.611613    1.724838   240.8407   -20.45873    1.029068    1.844211
!------------------------------------------------------------------------
!   10.0         2.302585   227.5813   -6.422815  -0.9810501    0.654469
!   95.042770    4.554327   215.6167  -0.8858135    3.440529    7.644918
!  174.847646    5.163915   218.0869    11.83129    17.42203    4.133282
!  288.471127    5.664595   228.8967    32.38542    23.63394   -22.63474
!  332.735207    5.807347   233.9356    37.74927    13.94737  -0.7294737
!  427.134521    6.057099   244.2222    44.57953    13.40285   -11.82933
!  526.326636    6.265922   254.0081    48.62964    5.924198  -1.9718339
!  663.976661    6.498247   265.6011    51.06303    4.503868   -5.541373
!  800.877499    6.685708   275.2952    52.16743    1.418971   -1.652569
! 1062.000848    6.967910   290.0925
!========================================================================

!WB
      DATA SK /-0.1531168D+01 ,-0.8657969D+00 ,-0.2282546D+00      & 
              , 0.4005346D+00 , 0.1045499D+01 , 0.1724838D+01      & 
!WB
              , 0.2302585D+01 , 0.4554327D+01 , 0.5163915D+01      &   
              , 0.5664595D+01 , 0.5807347D+01 , 0.6057099D+01      & 
              , 0.6265922D+01 , 0.6498247D+01 , 0.6685708D+01      &
              , 0.6967910D+01 /
!WB
      DATA TK / 0.2452323D+03 , 0.2579098D+03 , 0.2665730D+03      &
              , 0.2651181D+03 , 0.2546358D+03 , 0.2408407D+03      &
!WB
              , 0.2275813D+03 , 0.2156167D+03 , 0.2180869D+03      &
              , 0.2288967D+03 , 0.2339356D+03 , 0.2442222D+03      &
              , 0.2540081D+03 , 0.2656011D+03 , 0.2752952D+03      &
              , 0.2900925D+03 /
!WB ----------------------------- CK(N,1) --------------------------------------------------
      DATA CK / 19.870250D+00 , 18.099050D+00 , 5.2437200D+00      &
              ,-10.831250D+00 ,-19.119730D+00 ,-20.458730D+00      &
!WB  
              ,-6.4228150D+00 ,-0.8858135D+00 , 11.831290D+00      &
              , 32.385420D+00 , 37.749270D+00 , 44.579530D+00      &
              , 48.629640D+00 , 51.063030D+00 , 52.167430D+00      &
!WB ------------------------------- CK(N,2) ---------------------------------------------
              , 0.2870293D+00 ,-4.3148250D+00 ,-12.550260D+00      &
              ,-11.008370D+00 ,-3.1313780D+00 ,  1.029068D+00      &
!WB 
              ,-0.9810501D+00 ,  3.440529D+00 , 17.422030D+00      &
              , 23.633940D+00 , 13.947370D+00 , 13.402850D+00      &
              ,  5.924198D+00 ,  4.503868D+00 ,  1.418971D+00      &
!WB ------------------------------- CK(N,3) ----------------------------------------------
              , -2.276788D+00 , -4.329503D+00 , 0.8447471D+00      &
              ,  4.035521D+00 ,  2.037444D+00 ,  1.844211D+00      &
!WB 
              ,  0.654469D+00 ,  7.644918D+00 ,  4.133282D+00      &
              ,-22.634740D+00 ,-0.7294737D+00 ,-11.829330D+00      &
              ,-1.9718339D+00 ,-5.5413730D+00 ,-1.6525690D+00 /
!
      REAL(r8) :: THIRD ! THIRD = 1/3
      REAL(r8) :: SIXTH ! SIXTH = 1/6
      REAL(r8) :: HM0   ! THE GEOPOTENTIAL HEIGHT AT PGALIB (1062 hPa)
      REAL(r8) :: P     ! pressure
      integer  :: N0    ! index of first computation level
      REAL(r8) :: S     ! S   = ln(P)
      REAL(r8) :: DD,D  ! D   = ln(P) - ln(NP)
      REAL(r8) :: TW    ! TW  = TK(N)
      REAL(r8) :: CW    ! CW  = d(T)/d(lnP)
      REAL(r8) :: DW    ! DW  = d(CW)/d(lnP)
      REAL(r8) :: HW    ! HW  = (g/R) * [z(NP) - z(P)]
      REAL(r8) :: C02   ! C02 = C0 * C0
      REAL(r8) :: HMW   ! HMW = (g/R) * [z(IS) - z(IS+1)]
      REAL(r8) :: CKIII ! DS  = SK(I+1) - SK(I)
      REAL(r8) :: HUSA  ! geopotential height of 1976 US. standard atmosphere
      REAL(r8) :: TUSA  ! temperature of 1976 US. standard atmosphere
      REAL(r8) :: P0J   ! P0J = 1010.0 hPa
      REAL(r8) :: P0I   ! P0I = 1020.0 hPa
      REAL(r8) :: DP0   ! DP0 = P00SL / DPALIB
      REAL(r8) :: FI0, DP1, DS   
      integer  :: K, N, IS, NP, I, II, J0B, I0B, I00
!-----------------------------------------------------------------------------------------------

      THIRD = ONE / THREE
      SIXTH = ONE / SIX
!
      HM0   = GRAV * (-399.0E0)  
!     -399m is the lowest elevation of the world. (the Dead Sea)
      P     = PBALIB - DPALIB          ! PBALIB=0.5E0, DPALIB=0.5E0; P=0 hPa
      N0    = PBALIB / DPALIB + 0.001  ! N0 = 1.001  K = N0 = 1
!
      DO K  = N0,NA                    ! 1 -- 2320
         P  = P + DPALIB               ! P = ( 0.5,1160 ) hPa
         IF ( P.LT.PBALIB .OR. P.GT.PEALIB ) THEN  !PEALIB=1160 hPa
            PRINT*,'P IN STD.ATM IS OUT OF [0.5hPa,1160hPa] ==> STOP !'  
!!            STOP
         ELSE IF( P.LE.PGALIB ) THEN  !zhh PGALIB=1060.0E0
            S    = LOG( P )
            DD   = S - SK(1)  
            DO N = 1 ,KC   
               IS   = N
               D    = DD
               IF ( D.EQ.ZERO ) THEN     ! P = 10,95,175... 1060 mb
                  TW  = TK(N)
                  CW  = CK(N,1)          ! CW = d(T)/d(lnP)
                  DW  = TWO  * CK(N,2)   ! DW = d(CW)/d(lnP) 
                  HW  = ZERO                                           
                  GOTO 105
               ELSE                      ! P /= 0.2162828 mb
                  NP = N + 1
                  DD = S    - SK(NP)
                  IF ( DD.EQ.ZERO ) THEN  ! P = 95,175,288.47 ... 1062 mb
                     TW  = TK(NP)
                     CW  = CK(NP,1)
                     DW  = TWO  * CK(NP,2)
                     HW  = ZERO                                           
                     IF ( N.EQ.KC ) THEN
                        IS = KC
                     ELSE
                        IS = NP
                     END IF
                     GOTO 105
                  ELSE IF( DD.LT.ZERO ) THEN
                     TW  = ((CK(N,3)*D+CK(N,2))*D+CK(N,1))*D + TK(N)    !   (*1) 
                     CW  = (THREE*CK(N,3)*D + TWO*CK(N,2))*D + CK(N,1)  !   (*2)
                     DW  =    SIX*CK(N,3)*D + TWO*CK(N,2)               !   (*3)
! -------------------- INTEGRATE FROM S TO SK(IS=N) ------------------------------------------
                     HW   = ( ( (FOURTH*CK(N,3)*D + THIRD*CK(N,2) )*D   & 
                          +      HALF*CK(N,1) )*D + TK(N) )*D           !   (*4)
                     GOTO 105
!                 ELSE  ! DD.GT.ZERO              
                  END IF
               END IF
            END DO
!***********************************************************************************************
105    TBB (K)  = TW
!***********************************************************************************************
            C02      = RD * (CAPA*TW - CW)
            CBB (K)  = SQRT( C02 )
            DCBB(K)  = RD * (CAPA*CW - DW) / (C02+C02)
!           INTEGRATE FROM SK(IS) TO S(KC)
            HMW      = ZERO
            DO I = IS,KC
               II    = I  + 1
               IF ( I.EQ.KC ) THEN
                  CKIII  = CK(I ,2) + CK(I,2)
               ELSE
                  CKIII  = CK(II,2) + CK(I,2)
               END IF
               DS    = SK(II) - SK(I)
               HMW   = HMW + HALF*DS*( TK(II)+TK(I) - CKIII*DS*DS*SIXTH )  ! (*5)           
            END DO
            HBB (K)  = RD * ( HMW - HW ) + HM0   ! HM0 = GRAV * (-399.0E0)   (*6) 
         ELSE    !  1062 hPa < P < 1160 hPa
!--------------------------- Set the US Standard Atmosphere,1976 -------------------------------
            CALL USA76( P,HUSA,TUSA ) ! 
!-----------------------------------------------------------------------------------------------
            HBB(K)   = HUSA*GRAV      
            TBB(K)   = TUSA           
            CBB(K)   = CBB(NG)
            DCBB(K)  = 0.0D0          !  CBB is a constant.
         END IF
      END DO
!
      if ( N0.GT.1 ) then          
         DO K  = 1 ,N0-1         
            HBB(K)   = HBB(N0)
            TBB(K)   = TBB(N0)
            CBB(K)   = CBB(N0)
            DCBB(K)  = DCBB(N0)
         END DO
      end if
! ------------------ CALCULATE PRESSURE P00 & TEMPERATURE T00 AT THE SEA LEVEL ---------------
      P0J        = 1010.0D0
      P0I        = 1020.0D0
      J0B        = P0J/DPALIB + 1.D-4   ! INTEGER J0B
      I0B        = P0I/DPALIB + 1.D-4   ! INTEGER I0B
      P00SL      = (P0I*HBB(J0B)-P0J*HBB(I0B))/(HBB(J0B)-HBB(I0B))          !  (*7)
      DP0        = P00SL  / DPALIB      ! assumed P00SL=1013, then  DP0 = 202.6
      I00        = DP0    + 1.E-4       ! INTEGER I00 , I00 = 202    
      FI0        = I00                  ! REAL  FIO , FIO = 202.0
      DP1        = DP0    - FI0         ! REAL  DP1 , DP1 = 0.6
      T00SL      = (ONE - DP1)*TBB(I00) + DP1*TBB(I00+1)
      RETURN
   END SUBROUTINE

!================================================================================================
   SUBROUTINE SETMSA( EROR )
!------------------------------------------------------------------------------------------------
! SET CONST. & DICTIONARY OF THE MODEL STANDARD ATMOSPHERE              
!------------------------------------------------------------------------------------------------
      use IAP_prog,  only: PLY, PIN, GHS
      use Dyn_const, only: PMTOP
      IMPLICIT NONE
!------------------------------------Arguments--------------------------------------
      real(r8), intent(in) :: EROR          ! 
!----------------------------------Local workspace----------------------------------
      real(r8) :: RDT0,XP,XP2,PSB0,TSBJ,H0BJ
      INTEGER  :: I,J,ID
!------------------------------------------------------------------------------------------------

!----------------- CREATE LOOK-UP TABLES FOR THE MODEL STANDARD ATMOSPHERE -------------------
      CALL crestdatm	        
!-----------------------------------------------------------------------------------------------
      PS0        = P00SL  
      T00        = T00SL   
!
!    CALCULATE SURFACE PRESSURE BY NONLINEAR DESCENT METHOD
	
      RDT0       = - ONE / (RD*T00)
      DO J = beglatdyn, endlatdyn
         if ( j >= JB .and. j <= JE ) then 
            DO I = IB,IE
               GHS0 = GHS(I,J)  ! surface geopotential  zhh 2007.5.14
               XP   = PS0 * (ONE + GHS0*RDT0)  !    (*1)    see note P85
!---------------------- DESCENT METHOD SOLVING NONLINEAR EQUATIONS -----------------------------
               CALL NLDSMD( XP,SDGHI,EROR )  
!-----------------------------------------------------------------------------------------------
!   SDGHI=SDH*SDH;  SDH=HMSA(P)-GHS0;  HMSA: geopotential function of std atm     
               PSB (I,J)  = XP
               PSB2(I,J)  = SQRT( XP-PMTOP )  
            END DO
            call period( PSB (1,J) )
            call period( PSB2(1,J) )
!
!     CALCULATE H0B=RD*TMSA(PSB)/PSB
            DO I = IB, IE
               PSB0       = PSB(I,J)
               TSB(I,J)   = TMSA( PSB0 ) 
               H0B(I,J)   = TSB(I,J) * RD / PSB0
            END DO
            call period( TSB(1,J) )
            call period( H0B(1,J) )
!
         else   ! at north and south pole         
            GHS0 = GHS(IB,J)
            XP   = PS0 * (ONE + GHS0*RDT0)
!-----------------------------------------------------------------------------------
            CALL NLDSMD( XP,SDGHI,EROR )
!-----------------------------------------------------------------------------------
            XP2        = SQRT( XP-PMTOP )
            DO I = 1 ,NX
               PSB (I,J)  = XP
               PSB2(I,J)  = XP2
            END DO
!
            PSB0       = PSB(IB,J)
            TSBJ       = TMSA( PSB0 )
            H0BJ       = TSBJ * RD / PSB0
            DO I = 1 ,NX
               TSB(I,J)   = TSBJ
               H0B(I,J)   = H0BJ
            END DO
         end if
!
      END DO
!
      RETURN
   END SUBROUTINE

!================================================================================================
   SUBROUTINE DIAGBB(ITYPE)
!------------------------------------------------------------------------------------------------
! SET MSA TEMPERATURE, GEOPOTENTIAL & CHARACTERISTIC PHASE-SPEED BY LOOK-UP DICTIONARY              
!------------------------------------------------------------------------------------------------
      use IAP_prog, only: PLY, PIN
      use cam_control_mod, only: adiabatic        !zhh, 2011-11-19
!
      IMPLICIT NONE
!------------------------------------Arguments--------------------------------------
      integer, intent(in) :: ITYPE                        ! selective index
!----------------------------------Local workspace----------------------------------
      INTEGER  :: I, J, K, KPK, KPI
      REAL(r8) :: WPK, WPJ, WPI
!-----------------------------------------------------------------------------------------------
!
      GOTO (100,200,300),ITYPE 
100   CONTINUE
      DO J = beglatdyn, endlatdyn
         DO K = beglev ,endlev
            DO I = 1 ,NX
               WPK       = PLY(I,K,J) / DPALIB ! DPALIB = 0.5 hPa  REAL  WPK
               KPK       = WPK + 1.0E-4        ! INTEGER KPK
               if ( KPK == 0 ) then 
                  CB(I,K,J) = CBB(1)
               else
                  WPJ       = WPK - KPK           ! REAL  WPJ
                  WPI       = ONE - WPJ           ! REAL  WPI
                  KPI       = KPK + 1
                  CB(I,K,J) = WPI * CBB(KPK) + WPJ * CBB(KPI) ! interpolate CB with pressure
               end if 
            END DO
         END DO
      END DO
      RETURN
!
200   CONTINUE
      DO J = beglatdyn, endlatdyn
         DO K = beglev ,endlev
            DO I = 1 ,NX
               WPK        = PLY(I,K,J) / DPALIB 
               KPK        = WPK + 1.0E-4        
               if ( KPK == 0 ) then 
                  TB(I,K,J) = TBB(1)
               else
                  WPJ       = WPK - KPK           ! REAL  WPJ
                  WPI       = ONE - WPJ           ! REAL  WPI
                  KPI       = KPK + 1
                  TB(I,K,J) = WPI * TBB (KPK) + WPJ*TBB (KPI)    ! TB = T(P)
               end if 
            END DO
         END DO
      END DO
!
      DO J = beglatdyn, endlatdyn
         DO K = beglev,endlev+1
            do I = 1 ,NX
               WPK         = PIN(I,K,J) / DPALIB
               KPK         = WPK + 1.0E-4
               if ( KPK == 0 ) then 
                  GHI0(I,K,J) = HBB(1)
               else
                  WPJ      = WPK - KPK           ! REAL  WPJ
                  WPI      = ONE - WPJ           ! REAL  WPI
                  KPI      = KPK + 1
                  GHI0(I,K,J) = WPI*HBB (KPK) + WPJ*HBB (KPI)       ! GHI0(P) = phi(P)
               end if 
            end do
         end do
      end do
      RETURN
! 
300   CONTINUE
      DO J = beglatdyn, endlatdyn
         DO K = beglev ,endlev
            DO I = 1 ,NX
               WPK        = PLY(I,K,J) / DPALIB 
               KPK        = WPK + 1.0E-4        
               if ( KPK == 0 ) then 
                  TB(I,K,J) = TBB(1)
                  CB(I,K,J) = CBB(1)
                  if (adiabatic) then
                     deltac(I,K,J) = 0.0
                  else
                     deltac(I,K,J) = CBB(1) * CBB(1) / (b0 * b0) - 1.0
                  end if    
               else
                  WPJ       = WPK - KPK           ! REAL  WPJ
                  WPI       = ONE - WPJ           ! REAL  WPI
                  KPI       = KPK + 1
                  TB(I,K,J) = WPI * TBB (KPK) + WPJ*TBB (KPI)    ! TB = T(P)
                  CB(I,K,J) = WPI * CBB(KPK)  + WPJ * CBB(KPI)   ! CB(P) = C0(P)
                  if (adiabatic) then
                     deltac(I,K,J) = 0.0
                  else
		             deltac(I,K,J) = CB(I,K,J) * CB(I,K,J) / (b0 * b0) - 1.0
                  end if
               end if 
            END DO
         END DO
      END DO
!
      do J = beglatdyn, endlatdyn
         do K = beglev ,endlev+1
            do I = 1 ,NX
               WPK         = PIN(I,K,J) / DPALIB
               KPK         = WPK + 1.0E-4
               if ( KPK == 0 ) then 
                  GHI0(I,K,J) = HBB(1)
               else
                  WPJ      = WPK - KPK           ! REAL  WPJ
                  WPI      = ONE - WPJ           ! REAL  WPI
                  KPI      = KPK + 1
                  GHI0(I,K,J) = WPI*HBB (KPK) + WPJ*HBB (KPI)       ! GHI0(P) = phi(P)
               end if 
            end do
         end do
      end do
!
      RETURN
   END SUBROUTINE

!================================================================================================     
   SUBROUTINE USA76 (P,HUSA,TUSA)
!------------------------------------------------------------------------------------------------
! Purpose:
!    Set the US Standard Atmosphere,1976 which is identical with the earlier 1962 standard up to
! 51 km, and with the International Civil Aviation Organization (ICAO) standard up to 32 km.      
! Original version: ICAOSA.F (IAP 9L )
! Supplemented:  ZhangHe 2005.05.18  (add to layer 32000m -- 47000m, improvement the accuracy
!                                          of some constants )     
!----------------------------------------------------------------------------------------------
      IMPLICIT NONE
!------------------------------------ Arguments ----------------------------------------------
      real(r8), intent(in)  ::  P        ! pressure
      real(r8), intent(out) ::  HUSA     ! geopotential height of 1976 US. standard atmosphere
      real(r8), intent(out) ::  TUSA     ! temperature of 1976 US. standard atmosphere
!----------------------------------Local workspace----------------------------------
      real(r8) :: HTROP0,HSTRA0,HINVE1,HINVE2  ! consts for computing HUSA & TUSA
      real(r8) :: GAMA1,GAMA3,GAMA4            ! (- dT/dz) for different kinds of layers
!----------------------------------------------------------------------------------------------

!----------------------- Set constants for computing HUSA & TUSA ------------------------------
      HTROP0  = 44330.77D0 / (1013.25D0**0.19026D0)
      HSTRA0  = 6341.62D0  * LOG(226.32D0) + 11000.D0
      HINVE1  = 2.1665D+5  * (54.7487D0**0.029271D0)
      HINVE2  = 8.16607D+5 * ( 8.6802D0**0.081960D0)
      GAMA1   = 6.5D-3
      GAMA3   = -1.0D-3
      GAMA4   = -2.8D-3  

!------------------See 'Fundamental Atmospheric Physics'(P.X.Sheng) P59-63 -------------------

      IF ( P.GE.226.32D0 ) THEN        !troposphere,           h < 11000m ; gama1 = 0.0065K/m
         HUSA =  44330.77D0 - HTROP0*(P**0.19026D0)
         TUSA =  288.15D0 - GAMA1*(HUSA-0)
      ELSE IF ( P.GE.54.7489D0 ) THEN  ! Isothermal,  11000m < h < 20000m ; gama2 = 0
         HUSA =  HSTRA0 - 6341.62D0*LOG(P)
         TUSA =  216.65D0
      ELSE IF ( P.GE.8.68014D0 ) THEN  ! Inversion,   20000m < h < 32000m ; gama3 = -0.0010K/m
         HUSA =  HINVE1/(P**0.029271D0) - 1.9665D5
         TUSA =  216.65D0 - GAMA3*(HUSA-20000.D0)
      ELSE IF ( P.GE.1.1091D0 ) THEN   ! Inversion,   32000m < h < 47000m ; gama4 = -0.0028K/m
         HUSA =  HINVE2/(P**0.081960D0) - 4.96607D5
         TUSA =  228.65D0 - GAMA4*(HUSA-32000.D0)
      ELSE
         print '(A22)','No SA data available !'
      ENDIF   
      RETURN
   END SUBROUTINE

!================================================================================================
   FUNCTION SDGHI(P)
!------------------------------------------------------------------------------------------------
! VARIANCE FUNCTION (HMSA(P)-GHS0)^2              
!------------------------------------------------------------------------------------------------
      IMPLICIT NONE
!------------------------------------ Arguments -------------------------------------------------
      real(r8), intent(in)  ::  P        ! pressure
!----------------------------------Local workspace-----------------------------------------------
      real(r8) :: SDGHI
      real(r8) :: SDH
!------------------------------------------------------------------------------------------------
!     GHS0 : SURFACE GEOPOTENTIAL HIGHT AT A GIVEN POINT
      SDH    = HMSA( P ) - GHS0
      SDGHI  = SDH * SDH
      RETURN
   END FUNCTION

!================================================================================================
   FUNCTION HMSA( P )
!------------------------------------------------------------------------------------------------
! SET H OF STD.ATM BY LOOK-UP DICTIONARY              
!------------------------------------------------------------------------------------------------
      IMPLICIT NONE
!------------------------------------ Arguments -------------------------------------------------
      real(r8), intent(in)  ::  P        ! pressure
!----------------------------------Local workspace-----------------------------------------------
      real(r8) :: HMSA
      real(r8) :: WPK, WPJ, WPI
      INTEGER  :: KPK, KPI
!------------------------------------------------------------------------------------------------

      WPK     = P / DPALIB
      KPK     = WPK + 1.0E-4
      if ( KPK == 0 ) then 
         HMSA = HBB(1)
      else
         WPJ  = WPK - KPK           
         WPI  = ONE - WPJ           
         KPI  = KPK + 1
         HMSA = WPI*HBB(KPK) + WPJ*HBB(KPI)
      end if 
      RETURN
   END FUNCTION

!================================================================================================
   FUNCTION TMSA( P )
!------------------------------------------------------------------------------------------------
! SET T OF STD.ATM BY LOOK-UP DICTIONARY              
!------------------------------------------------------------------------------------------------
      IMPLICIT NONE
!------------------------------------ Arguments -------------------------------------------------
      real(r8), intent(in)  ::  P        ! pressure
!----------------------------------Local workspace-----------------------------------------------
      real(r8) :: TMSA
      real(r8) :: WPK, WPJ, WPI
      INTEGER  :: KPK, KPI
!------------------------------------------------------------------------------------------------

      WPK     = P / DPALIB
      KPK     = WPK + 1.0E-4
      if ( KPK == 0 ) then 
         TMSA = TBB(1)
      else
         WPJ  = WPK - KPK           
         WPI  = ONE - WPJ           
         KPI  = KPK + 1
         TMSA = WPI*TBB(KPK) + WPJ*TBB(KPI)
      end if 
      RETURN
   END FUNCTION

end module

