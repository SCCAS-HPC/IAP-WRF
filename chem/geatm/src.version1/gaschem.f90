module gaschem
contains
           SUBROUTINE HETERO_REAC(ASO4,BC,DUST01,DUST02,DUST03,DUST04,SSA01,SSA02,SSA03,SSA04,&
                               DUSTSO4,DUSTNO3,DSSASO4, DSSANO3,FSO4_DUST,FNO3_DUST,FSO4_SSA,FNO3_SSA)
            include 'chm1.inc'
            include 'gas1.inc' 
            REAL :: ASO4, BC  ! MASS CONCENTRATIONS IN UG/M3
!            REAL :: RH,TE ! IN % AND K
            REAL :: DENSI(2)   ! DENSITY OF AEROSOLS IN G/CM3 MEAN ! ASO4 AND BC
            REAL :: R(2)       ! MEAN RADIUS OF AEROSOLS IN UM
            REAL :: DENSIDUST(4), DENSISSA(4) ! DENSITY OF DUST AND SEA SALT FOR 4 BINS
            REAL :: RDUST(4), RSSA(4)         ! MEAN RADIUS OF DUST AND SEA SALT  IN UM 
            REAL :: AREA(2)    ! AEROSOL SURFACE AREA DENSITY UM2/CM3 FOR ASO4, BC
            REAL :: AREADUST(4),AREASSA(4) ! AEROSOL SURFACE AREA DENSITY UM2/CM3 FOR ASO4 AND BC
            REAL :: FRH,FRHSSA  ! FRH is the Hygroscopic growth for ASO4 and sea salt 
            REAL :: GAMMA1(28),A,B  ! 9 HETEROGENEOUS REACTIONS UPTAKE COFFICIENT 
            REAL :: DG(28)      ! gas diffusion in CM2/S  
            REAL :: C(28)       ! AVERAGE  molec speed [cm/s] 
!            REAL :: RK_HET(28)  ! 1/s HETEROGENEOUS rate 28 reactions  and have been defined in gas.inc
            REAL :: RK_HET_TMP(11:28, 4) ! 1/s HETEROGENEOUS rate FOR DUST AND SEA FOUR BINS
           !!!! 1: (NH4)2SO4 2: BC  3: DUST 4 : SEA SALT     
            REAL :: FSO4_DUST(4),FNO3_DUST(5),FSO4_SSA(4),FNO3_SSA(4) ! 4 BINS
            DATA DENSI/1.7, 1.0/
            DATA R/0.24, 0.04/
            DATA DENSIDUST / 2.5, 2.65, 2.65, 2.65 /
            DATA DENSISSA  / 2.2, 2.2 , 2.2 , 2.2  /
            DATA RDUST     / 0.15, 0.8,  1.75,  2.95  /
            DATA RSSA      / 0.15, 0.8,  1.75,  2.95  /

! 28 HETEROGENEOUS REACTIONS
!    1  N2O5 + ASO4 -> 2HNO3
!    2  NO2  + BC   -> 0.5HONO+0.5HNO3
!    3  NO3  + ASO4 -> HNO3
!    4  HO2  + ASO4 -> 0.5H2O2
!    5  HCHO + ASO4 -> PRODUCTS
!    6  OH   + ASO4 -> PRODUCTS
!    7  O3   + BC   -> PRODUCTS
!    8  NO2  + BC   -> HONO
!    9  HNO3 + BC   -> NO2
!   10  N2O5 + BC   -> 2HNO3
!   11  O3   + DUST -> PRODUCTS 
!   12  HNO3 + DUST -> ANO3 + PRODUCTS 
!   13  NO2  + DUST -> 0.5HONO + 0.5HNO3
!   14  NO3  + DUST -> HNO3
!   15  N2O5 + DUST -> 2HNO3
!   16  OH   + DUST -> PRODUCTS
!   17  HO2  + DUST -> 0.5H2O2
!   18  H2O2 + DUST -> PRODUCTS
!   19  SO2  + DUST -> ASO4 
!   20  CH3COOH + DUST -> PRODUCTS
!   21  CH3OH   + DUST -> PRODUCTS
!   22  HCHO    + DUST -> PRODUCTS
!   23  N2O5 + SSA  -> 2HNO3
!   24  NO3  + SSA  -> HNO3
!   25  HO2  + SSA  -> 0.5HONO
!   26  SO2  + SSA  -> ASO4
!   27  NO3  + SSA  -> ANO3
!   28  HNO3 + SSA  -> ANO3



            DATA GAMMA1 / 0.1    , 1.E-4,  3.E-3,  2.5E-1, 2.2E-2,&
                          0.2    , 3.3E-4, 3.3E-4, 2.1e-2, 5.E-3, &
                          2.7E-5 , 1.7E-1, 2.1E-6, 1.0E-3, 3.0E-2,&
                          0.1    , 0.2,    2.0E-3, 1.0E-4, 1.0E-3,&
                          1.0E-5 , 1.0E-5, 5.0E-3, 1.0E-3, 2.0E-1,&
                          5.0E-2 , 1.7E-2, 0.5/

            DATA DG / 0.1, 0.1, 0.1, 0.25, 0.1, 0.1, 0.1, 0.1,0.1, & 
                    ! N2O5,NO2, NO3, HO2, HCHO, OH, O3,  NO2, HNO3
                      0.1, 0.1, 0.1, 0.1,  0.1, 0.1, 0.1, 0.25, 0.1, &
                    ! N2O5, O3, HNO3,NO2, NO3, N2O5,OH, HO2,H2O2,
                      0.1, 0.1, 0.1, 0.1,  0.1, 0.1, 0.25, 0.1, 0.1,&
                    ! SO2, CH3COOH, CH3OH, HCHO, N2O5, NO3, HO2, SO2, NO3
                      0.1 /
                    ! HNO3

!      FRH IS CALCULATED FRH = 1+A*(RH/100)**B BY OBS BY PAN, X.(2009)IN  ACPD IN   BEIJING  
            FRH = 1. + 2.3 * (RH/100.)**6.27
           
        IF ( RH < 100 ) FRHSSA = 4.8
        IF ( RH < 099 ) FRHSSA = 2.9
        IF ( RH < 095 ) FRHSSA = 2.4
        IF ( RH < 090 ) FRHSSA = 2.0
        IF ( RH < 080 ) FRHSSA = 1.8
        IF ( RH < 070 ) FRHSSA = 1.6
        IF ( RH < 050 ) FRHSSA = 1.0

            AREA(1) = 3.0 * ASO4 / ( DENSI(1) * R (1)) * FRH**2
            AREA(2) = 3.0 * BC   / ( DENSI(2) * R (2)) 
            AREADUST(1) = 3.0 * DUST01 / (  DENSIDUST(1) * RDUST (1) )
            AREADUST(2) = 3.0 * DUST02 / (  DENSIDUST(2) * RDUST (2) )
            AREADUST(3) = 3.0 * DUST03 / (  DENSIDUST(3) * RDUST (3) )
            AREADUST(4) = 3.0 * DUST04 / (  DENSIDUST(4) * RDUST (4) )
            AREASSA(1) = 3.0 * SSA01 / (  DENSISSA(1) * RSSA (1) ) * FRHSSA**2
            AREASSA(2) = 3.0 * SSA02 / (  DENSISSA(2) * RSSA (2) ) * FRHSSA**2 
            AREASSA(3) = 3.0 * SSA03 / (  DENSISSA(3) * RSSA (3) ) * FRHSSA**2
            AREASSA(4) = 3.0 * SSA04 / (  DENSISSA(4) * RSSA (4) ) * FRHSSA**2

!    TO ADJUST GAMMA DEPENDING ON RH AND TEMPERATURE
            GAMMA1(7) = 1.8e-4 * EXP (-1000./TE)
            A = 2.79E-4 + 1.3E-4*RH - 3.43E-6 * RH**2.&
               + 7.52E-8 * RH**3. 
            IF(TE .GE. 282) B = 4.E-2*(TE-294)
            IF(TE .LT. 282) B = 0.48
            GAMMA1(1) = A * 10.**B

            IF (RH < 62  ) GAMMA1 (23) = 0.005
            IF (RH >= 62 ) GAMMA1 (23) = 0.03
            IF (RH < 50  ) GAMMA1 (26) = 0.005
            IF (RH >=50  ) GAMMA1 (26) = 0.05

            IF(RH <= 15) THEN
                GAMMA1(18) = 3.33E-4
            ELSE IF(RH <= 25) THEN
                GAMMA1(18) = 3.5E-4
            ELSE IF(RH <= 35) THEN
                GAMMA1(18) = 3.55E-4
            ELSE IF(RH <= 40) THEN
                GAMMA1(18) = 3.6E-4
            ELSE IF(RH <= 50) THEN
                GAMMA1(18) = 4.1E-4
            ELSE IF(RH <= 60) THEN
                GAMMA1(18) = 4.6E-4
            ELSE IF(RH <= 65) THEN
                GAMMA1(18) = 5.2E-4
            ELSE IF(RH <= 70) THEN
                GAMMA1(18) = 6.03E-4
            ELSE
                GAMMA1(18) = 6.03E-4
            ENDIF

                RH_TMP = RH / 100. - 0.15
                RH_TMP = MAX(RH_TMP,0.1)
                TMP1 = 8. * RH_TMP
                TMP2 = (1.-8.)* RH_TMP
                BET  = TMP1 / (1.-RH_TMP) / (1.-TMP2)
                GAMMA1(12) = BET * 0.033
                GAMMA1(12) = GAMMA1(12) * 0.5535 - 0.0058
! GAMMA1(12) dericed by Vlasenko (2006,ACP) and Wei (2010, Phd thesis Modelng the effect of Heteorogeneous reqction on atmospgeric chemistry and aerosol properties ) , using liear with latetr data with former datea
! ***   TEST IMPACTS of HETEOGENEOUS CHEMISTRY
!                 GAMMA1(12) = 0.0  ! upper 0.17 lower 0.001
                 GAMMA1(19) = 0.0 !5.0E-07  ! 2.6E-04 5.0E-07                       
!    TO ESTIMATE THE avg. molec speed [cm/s]

            C(1) = 1.455e4 * sqrt(te/108.) ! N2O5
            C(2) = 1.455e4 * sqrt(te/46.)  ! NO2
            C(3) = 1.455e4 * sqrt(te/62.)  ! NO3
            C(4) = 1.455e4 * sqrt(te/33.)  ! HO2
            C(5) = 1.455e4 * sqrt(te/30.)  ! HCHO
            C(6) = 1.455e4 * sqrt(te/17.)  ! OH
            C(7) = 1.455e4 * sqrt(te/48.)  ! O3
            C(8) = 1.455e4 * sqrt(te/46.)  ! NO2
            C(9) = 1.455e4 * sqrt(te/63.)  ! HNO3
            C(10) = 1.455e4 * sqrt(te/108.)! N2O5
            C(11) = 1.455e4 * sqrt(te/48.) ! O3
            C(12) = 1.455e4 * sqrt(te/63.) ! HNO3
            C(13) = 1.455e4 * sqrt(te/46.) ! NO2
            C(14) = 1.455e4 * sqrt(te/62.) ! NO3
            C(15) = 1.455e4 * sqrt(te/108.)! N2O5
            C(16) = 1.455e4 * sqrt(te/17.) ! OH
            C(17) = 1.455e4 * sqrt(te/33.) ! HO2
            C(18) = 1.455e4 * sqrt(te/24.) ! H2O2
            C(19) = 1.455e4 * sqrt(te/80.) ! SO2
            C(20) = 1.455e4 * sqrt(te/60.) ! CH3COOH
            C(21) = 1.455e4 * sqrt(te/32.) ! CH3OH
            C(22) = 1.455e4 * sqrt(te/30.) ! HCHO
            C(23) = 1.455e4 * sqrt(te/108.)! N2O5
            C(24) = 1.455e4 * sqrt(te/62.) ! NO3
            C(25) = 1.455e4 * sqrt(te/33.) ! HO2
            C(26) = 1.455e4 * sqrt(te/80.) ! SO2
            C(27) = 1.455e4 * sqrt(te/62.) ! NO3
            C(28) = 1.455e4 * sqrt(te/63.) ! HNO3


!    TO ESTIMATE THE HETEROGENEOUS RATE CONSTANHT 
!               K = [r/Dg + 4/ (c * gamma)]-1 A 

            
            RK_HET(1) = 1./( 1.E-4*R(1)/DG(1) + 4./C(1)/GAMMA1(1)) * AREA(1) * 1.E-8
            RK_HET(2) = 1./( 1.E-4*R(2)/DG(2) + 4./C(2)/GAMMA1(2)) * AREA(2) * 1.E-8             
            RK_HET(3) = 1./( 1.E-4*R(1)/DG(3) + 4./C(3)/GAMMA1(3)) * AREA(1) * 1.E-8
            RK_HET(4) = 1./( 1.E-4*R(1)/DG(4) + 4./C(4)/GAMMA1(4)) * AREA(1) * 1.E-8
            RK_HET(5) = 1./( 1.E-4*R(1)/DG(5) + 4./C(5)/GAMMA1(5)) * AREA(1) * 1.E-8
            RK_HET(6) = 1./( 1.E-4*R(1)/DG(6) + 4./C(6)/GAMMA1(6)) * AREA(1) * 1.E-8
            RK_HET(7) = 1./( 1.E-4*R(2)/DG(7) + 4./C(7)/GAMMA1(7)) * AREA(2) * 1.E-8
            RK_HET(8) = 1./( 1.E-4*R(2)/DG(8) + 4./C(8)/GAMMA1(8)) * AREA(2) * 1.E-8
            RK_HET(9) = 1./( 1.E-4*R(2)/DG(9) + 4./C(9)/GAMMA1(9)) * AREA(2) * 1.E-8
            RK_HET(10)= 1./( 1.E-4*R(2)/DG(10) + 4./C(10)/GAMMA1(10)) * AREA(2) * 1.E-8

           DO I = 11, 22   ! FOR DUST
              RK_HET_TMP (I, 1) = 1./( 1.E-4*RDUST(1)/DG(I) + 4./C(I)/GAMMA1(I)) * AREADUST(1) * 1.E-8
              RK_HET_TMP (I, 2) = 1./( 1.E-4*RDUST(2)/DG(I) + 4./C(I)/GAMMA1(I)) * AREADUST(2) * 1.E-8
              RK_HET_TMP (I, 3) = 1./( 1.E-4*RDUST(3)/DG(I) + 4./C(I)/GAMMA1(I)) * AREADUST(3) * 1.E-8
              RK_HET_TMP (I, 4) = 1./( 1.E-4*RDUST(4)/DG(I) + 4./C(I)/GAMMA1(I)) * AREADUST(4) * 1.E-8

              RK_HET(I) = RK_HET_TMP (I, 1) + RK_HET_TMP (I, 2) + RK_HET_TMP (I, 3) + RK_HET_TMP (I, 4)
           ENDDO 

           DO I = 23, 28   ! FOR SEA SALT (SSA)
              RK_HET_TMP (I, 1) = 1./( 1.E-4*RSSA(1)/DG(I) + 4./C(I)/GAMMA1(I)) * AREASSA(1) * 1.E-8
              RK_HET_TMP (I, 2) = 1./( 1.E-4*RSSA(2)/DG(I) + 4./C(I)/GAMMA1(I)) * AREASSA(2) * 1.E-8
              RK_HET_TMP (I, 3) = 1./( 1.E-4*RSSA(3)/DG(I) + 4./C(I)/GAMMA1(I)) * AREASSA(3) * 1.E-8
              RK_HET_TMP (I, 4) = 1./( 1.E-4*RSSA(4)/DG(I) + 4./C(I)/GAMMA1(I)) * AREASSA(4) * 1.E-8

              RK_HET(I) = RK_HET_TMP (I, 1) + RK_HET_TMP (I, 2) + RK_HET_TMP (I, 3) + RK_HET_TMP (I, 4)

           ENDDO
 
           DO I = 1,28
             RK_HET(I) = AMAX1(RK_HET(I),1.E-20)
           ENDDO
              
              RK_ASO4 = RK_HET(19) + RK_HET(26)
              RK_ANO3 = RK_HET(12) + RK_HET(27) + RK_HET(28)

              FSO4_DUST (1) = RK_HET_TMP(19,1) / AMAX1(RK_ASO4, 1.E-20)
              FSO4_DUST (2) = RK_HET_TMP(19,2) / AMAX1(RK_ASO4, 1.E-20)
              FSO4_DUST (3) = RK_HET_TMP(19,3) / AMAX1(RK_ASO4, 1.E-20)
              FSO4_DUST (4) = RK_HET_TMP(19,4) / AMAX1(RK_ASO4, 1.E-20)
              
              FNO3_DUST (1) = RK_HET_TMP(12,1) / AMAX1(RK_ANO3, 1.E-20)
              FNO3_DUST (2) = RK_HET_TMP(12,2) / AMAX1(RK_ANO3, 1.E-20)
              FNO3_DUST (3) = RK_HET_TMP(12,3) / AMAX1(RK_ANO3, 1.E-20)
              FNO3_DUST (4) = RK_HET_TMP(12,4) / AMAX1(RK_ANO3, 1.E-20)

              FSO4_SSA (1) = RK_HET_TMP(26,1) / AMAX1(RK_ASO4, 1.E-20)
              FSO4_SSA (2) = RK_HET_TMP(26,2) / AMAX1(RK_ASO4, 1.E-20)
              FSO4_SSA (3) = RK_HET_TMP(26,3) / AMAX1(RK_ASO4, 1.E-20)
              FSO4_SSA (4) = RK_HET_TMP(26,4) / AMAX1(RK_ASO4, 1.E-20)

              FNO3_SSA (1) = ( RK_HET_TMP(27,1) + RK_HET_TMP(28,1) )/ AMAX1(RK_ANO3, 1.E-20)
              FNO3_SSA (2) = ( RK_HET_TMP(27,2) + RK_HET_TMP(28,2) )/ AMAX1(RK_ANO3, 1.E-20)
              FNO3_SSA (3) = ( RK_HET_TMP(27,3) + RK_HET_TMP(28,3) )/ AMAX1(RK_ANO3, 1.E-20)
              FNO3_SSA (4) = ( RK_HET_TMP(27,4) + RK_HET_TMP(28,4) )/ AMAX1(RK_ANO3, 1.E-20)

           IF ( DUSTNO3.GE. 0.2*(DUST01+DUST02+DUST03+DUST04) ) RK_HET(12) = 0.0
           IF ( DSSASO4.GE. (SSA01+SSA02+SSA03+SSA04) )  RK_HET(26) = 0.0
           IF ( DSSANO3.GE. (SSA01+SSA02+SSA03+SSA04) )  RK_HET(28) = 0.0

              DO IR = 9, 9 
                RK_HET ( IR ) = 0.0           
              ENDDO 

           END subroutine
           
          SUBROUTINE FEEVOLUTION(MYID, FEIIIC, FEIII, FEII, SO2, HNO3,&
                         RK_HETSO2, RK_HETHNO3,FRCL,FRCM,FRCH,&
                         SWDOWN,RH1,DT,SX,EX,SY,EY,K,IS)
! ********  THIS ROUTINE IS TO CALCULATE THE FE(III) TO FE(II) IN DUST *********
! *************  PARTICLES BY LUO(2005) JGR, DOI:10.1029/2005JD006059  *********
! ****************    FAN ET AL (2006,GRL, DOI:10.1029/2005GL024852    *********
          INTEGER :: MYID, SX,EX,SY,EY
          REAL :: FAVG ! AVERAGE GLOBAL MEAN SHORTWAVE FLUX
          REAL :: TOBS ! ESTIMATED DAECAY LIFETIME
          REAL :: KSR  ! DECAY RATE BU SOLAR 
          REAL :: KFC  ! THE RATE COEFFICIENT FROM FRESH DUST TO COATED DUST BY HETEOROGENEOUS REACTIONS
          REAL :: KN, KS ! TEMPORARY VARAIABLE
          REAL :: RFE  ! RATE OF FE DISSOLUTION GRAMS OF FE IN DISSOLVED PER GRAM OF FE IN FE2O3 PER SECOND
          REAL :: RD   ! GRAMS OF FE DISSOLVED PER GRAM OF FE2O3
          REAL :: A    ! THE SPECIIFIC SURFACE AREA OF FE2O3
          REAL :: M,N,W
          REAL :: KCLD ! DECAY RATE BY CLOUD
          REAL :: FRCAVG ! GLOBAL AVERAGE CLOUD FRACTION
          REAL :: FRC   
          REAL,DIMENSION(SX-1:EX+1,SY-1:EY+1) :: FEIII,FEIIIC,FEII,SO2,HNO3 ! SO2 AND HNO3 IN PPBV
          REAL,DIMENSION(SX-1:EX+1,SY-1:EY+1) :: RK_HETSO2,RK_HETHNO3 ! HETEROGENOUS REACTIONS RATES OF SO2 AND HNO3 ON DUST
          REAL,DIMENSION(SX-1:EX+1,SY-1:EY+1) :: SWDOWN,RH1 ! in %
          REAL,DIMENSION(SX-1:EX+1,SY-1:EY+1) :: FRCL,FRCM,FRCH !LOW, MIDDLE AND HIGH CLOUD FRACTION
         

          DATA FAVG / 535.25 /  ! W/M2
          DATA TOBS / 2.59E07 / ! SECONDS          
          DATA RD   / 1.E-10 / ! MOL/M2/S
          DATA A    / 100./  ! M2/G
          DATA FRCAVG /0.05/ ! NO UNIT

          DO I = SX-1,EX+1
          DO J = SY-1,EY+1

!            IF(I==44.AND.J==48.AND.K==1.AND.IS==1) PRINT*, SWDOWN(I,J),RH1(I,J),FEIII(I,J),FEII(I,J),'before'
! --- TO ESTIMATE THE IMPACT OF SOLAR FROM LUO(JGR, 2005) 
           KSR = SWDOWN(I,J) / FAVG / TOBS
           FEII(I,J)  = FEII(I,J) + FEIII(I,J) * ( 1. -  EXP(-DT*KSR) )
           FEIII(I,J) = FEIII(I,J)* EXP( -DT * KSR )  
!            IF(I==44.AND.J==48.AND.K==1.AND.IS==1) PRINT*, SWDOWN(I,J),RH1(I,J),FEIII(I,J),FEII(I,J),'solar'

! --- TO ESTIMATE THE IMPACT OF  CLOUD  BY LUO  (2005,JGR)
           IF ( K<=10.and.K>=4) FRC = FRCL(I,J)
           IF ( K>10.AND.K<=13.) FRC = FRCM(I,J)
           IF ( K>13) FRC = FRCH(I,J)
           KCLD =  FRC/FRCAVG/TOBS
           FEII(I,J)  = FEII(I,J) + FEIII(I,J) * ( 1. -  EXP(-DT*KCLD) )
           FEIII(I,J) = FEIII(I,J)* EXP( -DT * KCLD )
  

! --- TO ESTIMATE THE IMPACT OF HETEOROGENEOUS REACTIONS BY FAN ! (2006,GRL)
           RH = RH1(I,J)
           IF(RH < 25. ) KN = 0
           IF(RH <= 35..AND.RH >= 25.) KN = 5.E-6 * ( RH - 35.) / (35.-25.)
           IF(RH > 35. ) KN = 5.E-6
           IF(RH < 50. ) KS = 0.
           IF(RH > 60. ) KS = 3.E-6
           IF(RH>=50..AND.RH<=.60) KS = 3.E-6 * (RH -50.)/(60.-50.)

           KN = RK_HETHNO3(I,J)! 5.E-6
           KS = RK_HETSO2(I,J) ! 1.E-3

           KFC = KN * HNO3(I,J) + KS * SO2(I,J)
           FEIIIC(I,J) = FEIIIC(I,J) + FEIII(I,J) * KFC
           FEIIIC(I,J) = AMAX1(FEIIIC(I,J),1.E-20) 

           N = 2. ! MOLES FE/MOLE FE2O3
           M = 55.8 ! 55.8 G/MOLE
           W = 0.7 ! THE MASS FRACTION OF FE IN FE2O3
           RFE = RD * A * N * M / W ! g(feII)/s/g(feIII)

           FEII(I,J) = FeII(I,J) +    RFE * FEIIIC(I,J) * DT
           FEIIIC(I,J) = FEIIIC(I,J) -  RFE * FEIIIC(I,J) * DT
           FEIII(I,J) = FEIII(I,J) -  RFE * FEIIIC(I,J) * DT
 

           FEIII(I,J) = AMAX1 ( FEIII(I,J), 1.E-20)              
           FEII(I,J) = AMAX1 ( FEIII(I,J)*0.005,FEII(I,J) )

!           IF(I==44.AND.J==48.AND.K==1.AND.IS==1) PRINT*,SWDOWN(I,J),RH1(I,J),FEIII(I,J),FEII(I,J),'chemistry'
! --- TO ESTIMATE THE IMPACT OF SO4  HETEROGENEOUS REACTIONS(FE(III)+H2O+SO4-->FE(II)+SO4) BY LUO  (2005,JGR)
           



          ENDDO ! J
          ENDDO ! I

          RETURN
          END SUBROUTINE
          
       SUBROUTINE CHEMOPE(MYID, OPE, DT, I, J, K, NE)
!C  TO CALCULATE THE NET OPE (OZONE PRODUCTION EFFICIENCY)FOR O3      
!C     OPE = (P(O3)-P(O3))/P(NOZ) 
   
      include 'chm1.inc'
      include 'gas1.inc'
      integer myid,i,j,k,ne
      real    DT,OPE,PO3,LO3,PNOZ

!ozone production and loss,and loss due to nox,radicals from radical-radical reactions
!P(o3)=k4[NO][HO2]+k5[NO][CH3O2]+k6[RO2][NO]
!L(o3)=k3[O1D][H2O]+k8[O3][HO2]+k7[O3][OH]+(k10[NO][O3]*k12[NO2][OH]/{k11[NO2]+k12[NO2][OH]})
      
      PO3  = 60.*dt*(rk_com(33)*cnn(kno)*cnn(kho2)& !k4[NO][HO2] inmolec/cc
            +rk_com(57)*cnn(kch3o2)*cnn(kno)& !k5[NO][CH3O2])
            +rk_urb(17)*cnn(kto2)*cnn(kno)&   !k[TO2][NO]
            +rk_com(58)*cnn(kethp)*cnn(kno)&  !K[ETHP][NO]
            +rk_urb(28)*cnn(kro2)*cnn(kno)&   !K[RO2][NO]
            +rk_com(71)*cnn(kc2o3)*cnn(kno)&  !K[c2o3][no]
            +rk_urb(29)*cnn(kano2)*cnn(kno)&  !k[ano2][no]
            +rk_urb(30)*cnn(knap)*cnn(kno)&   !k[nap][no]
            +rk_urb(31)*cnn(kxo2)*cnn(kno)&   !k[xo2][no]
            +rk_bio(8)*cnn(kisopp)*cnn(kno)&  !k[isopp][no]
            +rk_bio(9)*cnn(kisopn)*cnn(kno)&  !k[isopn][no]
            +rk_bio(10)*cnn(kisopo2)*cnn(kno)&!k[isopo2][no]
            )/cair_mlc*1.e+9

 
      LO3 = 60.*dt*(rk_com(12)*cnn(ko1d)*h2o& ! k3[O1D][H2O]
            +rk_com(21)*cnn(ko3)*cnn(kho2)&  ! k8[O3][HO2]
            +rk_com(20)*cnn(ko3)*cnn(koh)&   ! k7[O3][OH]
            +(rk_com(18)*cnn(kno)*cnn(ko3)*& !k10[NO][O3]
             rk_com(24)*cnn(kno2)*cnn(koh)/& !*k12[NO2][OH]
            (rk_com(1)*cnn(kno2)+rk_com(24)*cnn(kno2)*cnn(koh))& !k11[NO2]+k12[NO2][OH]
            )&      !(k10[NO][O3]*k12[NO2][OH]/{k11[NO2]+k12[NO2][OH]})
            )/cair_mlc*1.e+9

! CCCCCCCCCCC   NOZ PRODUCTION
! OH + NO2       = HNO3    
! NO3 + NO2      = N2O5      
! C2O3 + NO2     = PAN           
! CRO + NO2      = ONIT


      PNOZ =  60.*dt*(rk_com(24)*cnn(koh)*cnn(kno2) & ! OH +NO2
                  + rk_com(39)*cnn(kno3)*cnn(kno2) & ! NO3+NO2
                  + rk_com(69)*cnn(kc2o3)*cnn(kno2)& ! C2O3+NO2
                  + rk_urb(20)*cnn(kcro)*cnn(kno2) &  ! CRO +NO2 
                  ) / cair_mlc*1.e+9

      IF ( (PO3-LO3) .GE. 1.E-20 .AND. PNOZ.GE.1.E-20) THEN
       OPE = (PO3 - LO3) / PNOZ
      ELSE
       OPE = -1.E20
      ENDIF

      RETURN 
      END SUBROUTINE
      
         SUBROUTINE CALCLF(MYID,CLFLO,CLFMI,CLFHI,FCLD,RAINCV,RAINNCV,&
                       HE,TER,LON,LAT,TBEG_DD,TBEG_HH,TBEG_MM,TBEG_SS,&
                       I,J,K )
         INTEGER   :: MYID,SX,EX,SY,EY,NE
         REAL      :: CLFLO,CLFMI,CLFHI,CLF,FCLD,CTYPE,COPD !COPD:OPITICAL DEPTH
         REAL      :: H,TER,HE,LON,LAT
         INTEGER   :: ITOP  ! ITOP IS THE TYPE OF CLOUD
         REAL      :: cos_sa,tmid_sec1 !cosine of zenith angle 
         REAL      :: TBEG_DD,TBEG_HH,TBEG_MM,TNEG_SS
         REAL      :: HTOP(7) ! HEIGHT OF CLOUD M
         REAL      :: TR,AC,RAINCV,RAINNCV
         DATA HTOP / 2000., 4000., 6000., 4000., 6000., 6000., 6000./

         FCLD = 1.0
         tmid_sec1  = ((tbeg_dd*24.+tbeg_hh)*60.+tbeg_mm)*60.+tbeg_ss 
   
         CALL  ZenithAngle(tmid_sec1,LON,LAT,cos_sa)
 
         COS_SA = MAX(0.5,COS_SA)
         
         CLF = MAX(CLFLO,CLFMI,CLFHI)
         
         IF(CLFLO.GE.0.5.AND.CLFMI.GE.0.5.AND.CLFHI.GE.0.5) THEN
           CTYPE = 7.      
           ITOP  = 7
         ELSE IF(CLFLO.GE.0.5.AND.CLFMI.GE.0.5) THEN
           CTYPE = 6.
           ITOP  = 4
         ELSE IF(CLFLO.GE.0.5.AND.CLFHI.GE.0.5) THEN     
           CTYPE = 5.
           ITOP  = 5
         ELSE IF(CLFMI.GE.0.5.AND.CLFHI.GE.0.5) THEN
           CTYPE = 4.      
           ITOP  = 6
         ELSE IF(CLF == CLFLO.AND.CLF.GE.0.1) THEN
           CTYPE = 3.      
           ITOP  = 1
         ELSE IF(CLF == CLFMI.AND.CLF.GE.0.1) THEN 
           CTYPE = 2.
           ITOP  = 2 
         ELSE IF(CLF == CLFHI.AND.CLF.GE.0.1) THEN
           CTYPE = 1.
           ITOP  = 3 
         ELSE 
           CTYPE = 0.0      
         ENDIF

         IF(CTYPE == 7 ) COPD = 82.
         IF(CTYPE == 6 ) COPD = 80.
         IF(CTYPE == 5 ) COPD = 52.
         IF(CTYPE == 4 ) COPD = 32.
         IF(CTYPE == 3 ) COPD = 50.
         IF(CTYPE == 2 ) COPD = 30.
         IF(CTYPE == 1 ) COPD = 2.

         IF(COPD>=5.) THEN
          TR = (5.- 1./EXP(COPD))/(4.+3.*COPD*(1.-0.86)) ! THE ENEGERGY TRANSMISSION COFFICIENT
          H = HE - TER
          IF(H.GT.HTOP(ITOP)) THEN         
           AC = 1. +  (1.-TR)*COS_SA
          ELSE  
           AC = 1.6*TR*COS_SA 
          ENDIF
           FCLD = 1.+ CLF*(AC -1. )    
         ELSE
           FCLD = 1.0
         ENDIF       
!         IF(FCLD>1) PRINT*,FCLD,AC,TR,COPD,H
         IF(CLF<=0.1) FCLD = 1.0 
         IF(RAINCV>=0.1) COPD = 100.
        RETURN
        END  subroutine    

        subroutine ZenithAngle(tmid_sec1,lon,lat,cos_sa)
          real :: tmid_sec1,rlon,rlat,cos_sa,lon,lat
          real :: tlocal,tdec,codec,sidec,tloc,thou
          data deg2rad/0.017453293/ 
          rlon = lon*deg2rad
          rlat = lat*deg2rad
          tlocal=tmid_sec1
          tdec=0.4092797*sin(1.992385E-7*tlocal)
          sidec=sin(tdec)
          codec=cos(tdec)
          tloc=7.272205E-5*tlocal
          thou=cos(rlon+tloc)
          cos_sa=sin(rlat)*sidec+cos(rlat)*codec*thou 
        return
        end subroutine

      SUBROUTINE CALCLDOPD (CLWP, CLDOPD, RAINCV, RAINNCV) 
! ***    TO CALCULATE THE CLOUD OPTICAL DEPTH  ****                     
! **  USING LIU AND MAO(2008) IN ACTA SCIENTRIARUM NATURALIUM UNVERSITAR
                                                                        
                                                                        
      REAL ::CLWP, CLDOPD, RAINCV, RAINNCV 
                                                                        
                                                                        
!***  CLWP IS COLUMN CLOUD WATER IN G/M2 ASSUMING RATUS IS 11.85UM      
!***  CLDOPD IS CLOUD OPTICAL DEPTH                                     
!***  RAINCV IS CONVECTIVE RAINFALL                                     
!***  RAINNCV IS NON-CONVECTIVE RAINFALL                                
!***  WE ASSUME RADIUS OF CLOUD DROPETS IS 11.85 UM FROM LIU AND MAO, FO
!      35UM, NON-CONVECTIVE IS 22.5UM                                   
      IF (RAINCV.GE.0.05) THEN 
         R = 35. 
      ELSEIF (RAINNCV.GE.0.01) THEN 
         R = 22.5 
      ELSE 
         R = 11.85 
      ENDIF 
                                                                        
      CLDOPD = CLWP / (2. * R / 3.) 
                                                                        
      CLDOPD = MAX (CLDOPD, 1.0) 
                                                                        
      END SUBROUTINE CALCLDOPD                      

                                                                    
      subroutine JDAY (iday, imonth, iyear, julday, iseason)            
!-----------------------------------------------------------------------
!                                                                       
!                                                                       
! PURPOSE:  Compute Julian day and season index from date               
!                                                                       
! INPUTS:                                                               
!                                                                       
!     IDAY, IMONTH, IYEAR integer     day, month, year                  
!                                                                       
! OUTPUT:                                                               
!                                                                       
!     JULDAY      integer     Julian day                                
!     ISEASON     integer     season index (1 = spring, 2 = summer,     
!                                           3 = fall,   4 = winter)     
!                                                                       
! CALLS:                                                                
!                                                                       
!     none                                                              
!                                                                       
!-----------------------------------------------------------------------
                                                                        
      implicit none                                                     
                                                                        
      integer iday, imonth, iyear, julday, iseason                      
      integer ndaynorm(13), ndaybis(13)                                 
                                                                        
      data ndaynorm /0,31,59,90,120,151,181,212,243,273,304,334,365/    
      data ndaybis /0,31,60,91,121,152,182,213,244,274,305,335,366/     
                                                                        
                                                                        
      if (mod(iyear,4) .eq. 0) then                                     
        julday = ndaybis(imonth) + iday                                 
        if (julday .lt. 80 .or. julday .ge. 355) then                   
          iseason = 4                                                   
        else if (julday .lt. 172) then                                  
          iseason = 1                                                   
        else if (julday .lt. 264) then                                  
          iseason = 2                                                   
        else                                                            
          iseason = 3                                                   
        endif                                                           
      else                                                              
        julday = ndaynorm(imonth) + iday                                
        if (julday .lt. 81 .or. julday .ge. 356) then                   
          iseason = 4                                                   
        else if (julday .lt. 173) then                                  
          iseason = 1                                                   
        else if (julday .lt. 265) then                                  
          iseason = 2                                                   
        else                                                            
          iseason = 3                                                   
        endif                                                           
      endif                                                             
                                                                        
      return                                                            
      end  subroutine 
end module gaschem      
