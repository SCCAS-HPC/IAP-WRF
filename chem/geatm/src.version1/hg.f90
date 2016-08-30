 module hgchem
 contains
      SUBROUTINE HGAQSCHEM (CHG0, CHG2, CO3, CSO2, COH, CHO2, CPM, CPH, &
      LWC, TEMP, PRES, LAT, LON, HEIZ, LANDMASK, TMID_SEC, DELT)        
                                                                        
!                                                                       
!----CAMx v5.40 111010                                                  
!                                                                       
!   HGAQSCHEM performs the Hg aqueous-phase chemistry calculations.     
!                                                                       
!     Copyright 2003                                                    
!     AER, Inc.                                                         
!                                                                       
! Revision History:                                                     
! chenhs modify it to incorporate into GRACT, 201207                    
!***********************************************************************
!   Version 1.0 written March 2003 by Prakash Karamchandani, AER       *
!   Code updated August 2003 to use Henry's Law and SO2 dissociation   *
!   constants consistent with those used in the CAMx aqueous-phase     *
!   chemistry module                                                   *
!***********************************************************************
!                                                                       
!  Called by:  CHEMDRIV                                                 
!                                                                       
!***********************************************************************
                                                                        
      IMPLICIT NONE 
                                                                        
! Includes:                                                             
                                                                        
                              ! Hg aqueous-phase chemistry parameters   
      INCLUDE'hgaqschm.prm' 
                                                                        
! Arguments:                                                            
                                                                        
                      ! Hg(0) concentration, ppm (Input/Output)         
      REAL CHG0 
                      ! Hg(2) concentration, ppm (Input/Output)         
      REAL CHG2 
                                                                        
                      ! O3 concentration, ppm (Input only)              
      REAL CO3 
                      ! SO2 concentration, ppm (Input only)             
      REAL CSO2 
                      ! OH concentration, ppm (Input only)              
      REAL COH 
                      ! HO2 concentration, ppm (Input only)             
      REAL CHO2 
                      ! Total PM concentration, ug/m3 (Input only)      
      REAL CPM 
                                                                        
                      ! Cloudwater pH (Input only)                      
      REAL CPH 
                      ! Liquid water content, g/m3 (Input only)         
      REAL LWC 
                      ! Temperature, K (Input only)                     
      REAL TEMP 
                      ! Pressure, mb (Input only)                       
      REAL PRES 
                      ! Chemistry time step, seconds (Input only)       
      REAL DELT 
                                                                        
                      ! Latitude, Degree (Input only)                   
      REAL LAT 
                      ! Longitude, Degree (Input only)                  
      REAL LON 
                      ! Model layer height, m (Input only)              
      REAL HEIZ 
                      ! 1=land, 0=water (Input only)                    
      REAL LANDMASK 
                      ! time in seconds from Greenich Noon March 21 (Inp
      REAL TMID_SEC 
                                                                        
! Local variables                                                       
                      ! HCL concentration, ppm (Input only)             
      REAL CHCL 
                      ! CL2 concentration, ppm (Input only)             
      REAL CCL2 
                        ! 1=daytime, 0=nighttime                        
      INTEGER ifdaytime 
                      ! cosine of solar zenith angle                    
      REAL cos_sa 
                      ! cutoff solar zenith angle                       
      REAL sza_cut 
      PARAMETER (sza_cut = 89.0) 
                       !cos of sza_cut                                  
      REAL cos_sza_cut 
      PARAMETER (cos_sza_cut = 0.017452406) 
      INTEGER ihg, iht 
                                                                        
                      ! Henry's Law constant for SO2, M/atm (Input only)
      REAL SO2H 
      PARAMETER (SO2H = 1.22e+00) 
                      ! Henry's Law constant for O3, M/atm (Input only) 
      REAL O3H 
      PARAMETER (O3H = 1.10e-02) 
                      ! SO2 Henry's Law temperature dependence (Input on
      REAL SO2TMP 
      PARAMETER (SO2TMP = - 3156.) 
                      ! O3 Henry's Law temperature dependence (Input onl
      REAL O3TMP 
      PARAMETER (O3TMP = - 2415.) 
! Hg Molecular Weight                                                   
      REAL HG_MOLWT 
      PARAMETER (HG_MOLWT = 200.6) 
                                                                        
! Standard atmosphere in mb                                             
      REAL STDATMMB 
      PARAMETER (STDATMMB = 1013.25) 
                                                                        
! Pressure in atmospheres                                               
      REAL PRES_ATM 
                                                                        
! Conv. factor from ppm to atmospheres                                  
      REAL CFACT 
                                                                        
! Conversion factor (M to atm), liter-atm/mole                          
      REAL XL 
                                                                        
! Henry's Law constant for OH and HO2, M/atm                            
      REAL OHH, HO2H 
                                                                        
! Temperature in K and Centigrade and inverse temperature (K-1)         
      REAL TEMPK, TEMPC, TEMPI 
                                                                        
! Double precision local variables (for accuracy):                      
                                                                        
                              ! Initial total Hg(0) partial pressure, at
      DOUBLE PRECISION Y00 
                              ! Initial total Hg(2) partial pressure, at
      DOUBLE PRECISION Y20 
                                                                        
                              ! Final total Hg(0) partial pressure, atm 
      DOUBLE PRECISION Z0T 
                              ! Final total Hg(2) partial pressure, atm 
      DOUBLE PRECISION Z2T 
                                                                        
! --- Local concentration array                                         
      DOUBLE PRECISION Y (NUMSP) 
!     Y(1)         [Hg(0)]-Gaseous (variable), atm                      
!     Y(2)         [Hg(2)]-Gaseous (variable), atm                      
!     Y(3)         [O3]-Gaseous (constant), atm                         
!     Y(4)         [Cl2]-Gaseous (constant), atm                        
!     Y(5)         [HCl]-Gaseous (constant), atm                        
!     Y(6)         [SO2]-Gaseous (constant), atm                        
!     Y(7)         [OH]-Aqueous (constant), M                           
!     Y(8)         [HO2]-Aqueous (constant), M                          
                                                                        
                            ! PM10 conc. (g/L)                          
      DOUBLEPRECISION PM 
                                                                        
! --- Intermediate variables for calculating Hg(2) complexes with SO2   
! --- in solution and their reaction rates                              
      DOUBLEPRECISION CON1, CON2, CON3, CON4, CON51, CON52 
                                                                        
! --- Intermediate variables for calculating Hg(2) concentrations in    
! --- solution as Hg(2)2+, HgCl2 and Hg(OH)2 and total Hg(2)            
      DOUBLEPRECISION CON5, CON6, CON7, TOTHG 
                                                                        
! --- Aqueous concentrations (M) of dissociation products of Cl2,       
! --- HOCL and OCL-                                                     
      DOUBLEPRECISION CONHOCL, CONOCL 
                                                                        
! --- Aqueous concentrations (M) of HCl and Cl2                         
      DOUBLEPRECISION CONHCL, CONCL 
                                                                        
! --- Intermediate variables for calculating Cl species distribution    
! --- in solution                                                       
      DOUBLEPRECISION THCL, TCL2, CA, CB, CC, CD, AA, AB, AC, AD 
                                                                        
! --- Aqueous equilibrium constants for SO2 and HSO3                    
! --- SO2       <=>  H+ + HSO3-       (M)                               
! --- HSO3-     <=>  SO3-- + H+       (M)                               
      DOUBLEPRECISION EKSO2, EKHSO3 
                                                                        
! --- Henry's Law constants (M/atm)                                     
      DOUBLEPRECISION HNRHG0, HNRHGCL2, HNRHGOH2 
      DOUBLEPRECISION HNRSO2, HNRO3, HNRHCL, HNRCL2 
                                                                        
! --- H+ concentration, OH- concentration                               
      DOUBLEPRECISION PHCON, POHCON 
                                                                        
                                ! Conversion factor (M to atm)          
      DOUBLEPRECISION RLCONV 
                                                                        
! --- Chemistry parameters                                              
                             ! Chemistry time step, sec                 
      DOUBLEPRECISION DT 
                                               ! Intermediate variables 
      DOUBLEPRECISION ALPHA1, ALPHA2, ALPH12 
      DOUBLEPRECISION DENOM1, DENOM2 
                                                                        
! --- Temperature-dependent rate for HGSO3 reaction                     
      DOUBLEPRECISION RKHGSO3I 
                                                                        
! --- Vertical profiles of HCl and Cl2 (ppm) used by the Hg chemistry   
      INTEGER NHTCL 
      PARAMETER (NHTCL = 6) 
                                                                        
      REAL htcl (NHTCL) 
      REAL hclprof (NHTCL) 
      REAL cl2day (NHTCL) 
      REAL cl2nite (NHTCL) 
      DATA htcl / 0.06, 0.15, 0.45, 0.85, 2.00, 6.00 / 
      DATA hclprof / 1.0E-03, 1.0E-03, 1.0E-03, 1.0E-03, 1.0E-03,       &
      1.0E-03 /                                                         
      DATA cl2day / 1.0E-06, 1.0E-06, 1.0E-06, 1.0E-06, 1.0E-06,        &
      1.0E-06 /                                                         
      DATA cl2nite / 1.5E-04, 1.0E-04, 7.5E-05, 5.0E-05, 5.0E-05,       &
      5.0E-05 /                                                         
!---- to judge if it is daytime or night time                           
      CALL ZenithAngle00 (tmid_sec, lon, lat, cos_sa) 
                                      ! daytime                         
      IF (cos_sa.ge.cos_sza_cut) then 
         ifdaytime = 1 
            !nighttime                                                  
      ELSE 
         ifdaytime = 0 
      ENDIF 
!---- Assign HCl and Cl2 concentrations based on height(km), ocean and d
                                                                        
      DO ihg = 1, NHTCL 
      IF ( (heiz / 1000.) .le.htcl (ihg) ) then 
         iht = ihg 
         GOTO 115 
      ENDIF 
      enddo 
                           !! layer altitude higher than 6.0 km         
      iht = NHTCL 
  115 CONTINUE 
      CHCL = hclprof (iht) 
      IF (LANDMASK.eq.0) then 
         IF (ifdaytime.eq.0) then 
            CCL2 = cl2nite (iht) 
         ELSE 
            CCL2 = cl2day (iht) 
         ENDIF 
      ELSE 
         CCL2 = 0.0 
      ENDIF 
!-----------------------------------------------------------------------
! --- Assign values to local variables                                  
                                                                        
! --- Local copy of temperature (don't allow temperature to go below    
! --- freezing)                                                         
      TEMPK = MAX (TEMP, REAL (TZERO) ) 
                                                                        
! --- Inverse temperature                                               
      TEMPI = 1.0 / TEMPK 
                                                                        
! --- Calculate partial pressures in atm.                               
      PRES_ATM = PRES / STDATMMB 
      CFACT = 1.E-6 * PRES_ATM 
                                                                        
! --- convert conc unit, ng/m3 to ppm, by chenhs                        
      CHG0 = CHG0 * 1.0E-3 * (0.08206 * TEMP) / PRES_ATM / HG_MOLWT *   &
      1.0E-3                                                            
      CHG2 = CHG2 * 1.0E-3 * (0.08206 * TEMP) / PRES_ATM / HG_MOLWT *   &
      1.0E-3                                                            
! --- convert conc unit, ppb to ppm, by chenhs                          
      CO3 = CO3 * 1.0E-3 
      CSO2 = CSO2 * 1.0E-3 
      COH = COH * 1.0E-3 
      CHO2 = CHO2 * 1.0E-3 
                                                                        
! --- XL is the conversion factor from M (moles/liter of water) to atm  
      XL = (LWC / H2ODENS) * (MOLVOL * TEMPK / TZERO) 
                                                                        
      Y00 = CHG0 * CFACT 
      Y20 = CHG2 * CFACT 
                                                                        
      Y (KO3G) = CO3 * CFACT 
      Y (KSO2G) = CSO2 * CFACT 
      Y (KHCLG) = CHCL * CFACT 
      Y (KCL2G) = CCL2 * CFACT 
                                                                        
! --- Calculate liquid-phase concs of OH and HO2 radicals               
                                                                        
! --- Henry's Law constants for OH and HO2                              
! --- Values below are from Jacobson's book (1999)                      
      OHH = 25. * EXP (17.72 * (298.15 * TEMPI - 1.0) ) 
      HO2H = 2000. * EXP (22.28 * (298.15 * TEMPI - 1.0) ) 
                                                                        
! --- Since liquid-phase chemistry of OH and HO2 radicals is not explici
! --- simulated in CAMx, assume that in-cloud processes remove 50% of OH
! --- 90% of HO2                                                        
      Y (KOHAQ) = COH * CFACT * OHH * 0.5 / (1.0 + OHH * XL) 
      Y (KHO2AQ) = CHO2 * CFACT * HO2H * 0.1 / (1.0 + HO2H * XL) 
                                                                        
! --- H+ and OH- aqueous concentrations in moles/liter                  
      PHCON = 10.** ( - CPH) 
      POHCON = 10.** (CPH - 14) 
                                                                        
! --- concentration of PM in g/L of water                               
                                        ! No need to check if LWC = 0. s
      PM = CPM * H2ODENS * 1.E-9 / LWC 
                                        ! this routine is not called if 
                                        ! LWC < LWMIN                   
                                                                        
! --- Assign some double precision variables                            
                                                                        
! --- Conversion factor from M to atm.                                  
      RLCONV = XL 
                                                                        
! --- Chemistry time step                                               
      DT = DELT 
                                                                        
! --- Henry's Law constants for SO2 and O3                              
      HNRO3 = O3H * exp (O3TMP * (1. / 298. - 1. / tempk) ) 
      HNRSO2 = SO2H * exp (SO2TMP * (1. / 298. - 1. / tempk) ) 
      EKSO2 = 10.** (853. / tempk) / 54950. 
      EKHSO3 = 10.** (621.9 / tempk) / 1.897e+9 
                                                                        
! --- Henry's Law constants for Hg(0), HgCl2, Hg(OH)2, HCl and Cl2      
! --- (Not available from CAMx)                                         
! --- Special expression for Hg(0) Henry's Law constant                 
      HNRHG0 = 54.78775913129319 * EXP ( - 55.7339 + 9540.36 * TEMPI +  &
      16.0477 * ALOG (TEMPK / 100.) )                                   
      HNRHGCL2 = 1.4E6 
      HNRHGOH2 = 1.2E4 
                                                                        
      HNRHCL = 1.1 
                                                                        
! --- Special expression for Cl2 Henry's Law constant                   
      TEMPC = TEMPK - TZERO 
      HNRCL2 = 0.149 + TEMPC * ( - 3.59E-03 + TEMPC * (2.5E-05 + TEMPC *&
      5.67E-08) )                                                       
                                                                        
! --- Gas-liquid partitioning of reactant species                       
      Y (KO3G) = Y (KO3G) / (1.0 + RLCONV * HNRO3) 
                                                                        
! --- Calculate distribution of chlorine species                        
                                  ! Check for RLCONV not really required
      IF (RLCONV.GT.0.) THEN 
                                  ! since this routine is only called fo
                                  ! liquid water >= lwmin               
         THCL = Y (KHCLG) / RLCONV 
         TCL2 = 2. * Y (KCL2G) / RLCONV 
!                                                                       
! Calculate subcomponents                                               
         CA = 2. * (1.0 + PHCON / EKHOCL) 
         CB = 2. * (1.0 + 1.0 / (HNRCL2 * RLCONV) ) * PHCON**2. / EKCL2 &
         / EKHOCL                                                       
         CC = 1.0 + PHCON / EKHCL * (1.0 + 1.0 / (HNRHCL * RLCONV) ) 
         CD = TCL2 * (1.0 + PHCON / EKHOCL) 
!                                                                       
! Solve quadratic equation                                              
         AA = CB * CC 
         AB = CA * CC - CB * THCL 
         AC = - (CD+CA * THCL) 
         AD = AB * AB - 4. * AA * AC 
                                                                        
         CONCL = ( - AB + SQRT (AD) ) / (2. * AA) 
         CONHCL = CONCL * PHCON / EKHCL 
         Y (KHCLG) = CONHCL / HNRHCL 
         Y (KCL2G) = 0.5 * TCL2 / (1. / (RLCONV) + HNRCL2 * (1.0 +      &
         EKCL2 / (EKHCL * HNRHCL * Y (KHCLG) ) * (1.0 + EKHOCL / PHCON) &
         ) )                                                            
                                                                        
              ! ( RLCONV > 0. )                                         
      ENDIF 
                                                                        
! --- Calculate partitioning of S(IV), assuming that mercury sulfites   
! --- are negligible compared to other sulfite ions.                    
      Y (KSO2G) = Y (KSO2G) / (1.0 + RLCONV * HNRSO2 * (1.0 + EKSO2 /   &
      PHCON * (1.0 + EKHSO3 / PHCON) ) )                                
                                                                        
! --- Calculate constants                                               
      CON2 = EKHGSO3 * EKSO2 * EKHSO3 * HNRSO2 / (PHCON * PHCON) 
      CON3 = EKHGSO3 * EKHGSO32 * (EKSO2 * EKSO2) * (EKHSO3 * EKHSO3)   &
      * (HNRSO2 * HNRSO2) / (PHCON**4)                                  
      CON5 = (1.0 + Y (KSO2G) * (CON2 + CON3 * Y (KSO2G) ) ) 
      CON6 = ( (EKHCL * HNRHCL * Y (KHCLG) / PHCON) **2) / EKHGCL2 
      CON7 = (POHCON * POHCON) / EKHGOH2 
      CONHOCL = HNRCL2 * Y (KCL2G) * EKCL2 / (EKHCL * HNRHCL * Y (KHCLG)&
      )                                                                 
      CONOCL = CONHOCL * EKHOCL / PHCON 
                                                                        
      CON51 = CON2 * Y (KSO2G) 
      CON52 = CON3 * Y (KSO2G) * Y (KSO2G) 
      DENOM1 = (1.0 + RLCONV * HNRHG0) 
      TOTHG = CON5 + (CON51 + CON52) * EKP * PM + (CON6 + CON7) *       &
      (1.0 + EKP * PM)                                                  
      DENOM2 = CON6 / HNRHGCL2 + CON7 / HNRHGOH2 + RLCONV * TOTHG 
                                                                        
! --- Hg aqueous chemistry calculations for time-step DT                
                                                                        
! --- Calculate HGSO3 reaction rate (Van Loon et al., 2000)             
      RKHGSO3I = RKHGSO3 * EXP ( - 105000. / 8.3 * (TEMPI - TREFI) ) 
                                                                        
      CON1 = RKHGSO32 * CON3 
      CON4 = RKHGSO3I * CON2 
                                                                        
      ALPHA1 = (RLCONV * HNRHG0 * (RKO3 * HNRO3 * Y (KO3G) + RKOH * Y ( &
      KOHAQ) + RKHOCL * CONHOCL + RKOCL * CONOCL) ) / DENOM1            
                                                                        
      ALPHA2 = RLCONV / DENOM2 * (CON1 * (Y (KSO2G) * Y (KSO2G) )       &
      + CON4 * Y (KSO2G) + RKHO2 * Y (KHO2AQ) * (CON5 + CON6 + CON7) )  
                                                                        
      ALPH12 = ALPHA1 + ALPHA2 
                                                                        
! --- Calculate analytical solution                                     
      Z0T = ALPHA2 * (Y00 + Y20) / ALPH12 * (1.0 - EXP ( - ALPH12 * DT) &
      ) + Y00 * EXP ( - ALPH12 * DT)                                    
      Z2T = Y00 + Y20 - Z0T 
                                                                        
      CHG0 = Z0T / CFACT 
      CHG2 = Z2T / CFACT 
                                                                        
! The lower bound of HG0 and HG2 conc 1.0E-12 ppm                       
      CHG0 = max (CHG0, 1.0E-12) 
      CHG2 = max (CHG2, 1.0E-12) 
! convert conc unit, ppm to ng/m3, by chenhs                            
      CHG0 = CHG0 * 1.0E+3 * HG_MOLWT * PRES_ATM / (0.08206 * TEMP)     &
      * 1.0E+3                                                          
      CHG2 = CHG2 * 1.0E+3 * HG_MOLWT * PRES_ATM / (0.08206 * TEMP)     &
      * 1.0E+3                                                          
                                                                        
      RETURN 
      END SUBROUTINE HGAQSCHEM
      
      SUBROUTINE HGGASCHEM (CHG0, CHG2, CO3, CH2O2, COH, TEMP, PRES,    &
      LAT, LON, HEIZ, LANDMASK, TMID_SEC, DELT)                         
!                                                                       
!----CAMx v5.40 111010                                                  
!                                                                       
!   HGGASCHEM performs the Hg gas-phase chemistry calculations.         
!                                                                       
!     Copyright 2003                                                    
!     AER, Inc.                                                         
!                                                                       
! Revision History:                                                     
! chenhs modify it to incorporate into GRACT, 201207                    
!***********************************************************************
!   Version 1.0 written March 2003 by Prakash Karamchandani, AER       *
!***********************************************************************
!                                                                       
!  Called by:  CHEMDRIV                                                 
!                                                                       
!***********************************************************************
                                                                        
      IMPLICIT NONE 
                                                                        
! Includes:                                                             
                                                                        
                              ! Hg gas-phase chemistry parameters       
      INCLUDE'hggaschm.prm' 
                                                                        
! Arguments:                                                            
                                                                        
                      ! Hg(0) concentration, ppm (Input/Output)         
      REAL CHG0 
                      ! Hg(2) concentration, ppm (Input/Output)         
      REAL CHG2 
                                                                        
                      ! O3 concentration, ppm (Input only)              
      REAL CO3 
                      ! H2O2 concentration, ppm (Input only)            
      REAL CH2O2 
                      ! OH concentration, ppm (Input only)              
      REAL COH 
                                                                        
                      ! Temperature, K (Input only)                     
      REAL TEMP 
                      ! Pressure, mb (Input only)                       
      REAL PRES 
                      ! Chemistry time step, seconds (Input only)       
      REAL DELT 
                                                                        
                      ! Latitude, Degree (Input only)                   
      REAL LAT 
                      ! Longitude, Degree (Input only)                  
      REAL LON 
                      ! Model layer height, m (Input only)              
      REAL HEIZ 
                      ! 1=land, 0=water (Input only)                    
      REAL LANDMASK 
                      ! time in seconds from Greenich Noon March 21 (Inp
      REAL TMID_SEC 
                                                                        
! Local variables                                                       
                      ! HCL concentration, ppm (Input only)             
      REAL CHCL 
                      ! CL2 concentration, ppm (Input only)             
      REAL CCL2 
                        ! 1=daytime, 0=nighttime                        
      INTEGER ifdaytime 
                      ! cosine of solar zenith angle                    
      REAL cos_sa 
                      ! cutoff solar zenith angle                       
      REAL sza_cut 
      PARAMETER (sza_cut = 89.0) 
                       !cos of sza_cut                                  
      REAL cos_sza_cut 
      PARAMETER (cos_sza_cut = 0.017452406) 
      INTEGER ihg, iht 
! Standard atmosphere in mb                                             
      REAL STDATMMB 
      PARAMETER (STDATMMB = 1013.25) 
! Hg Molecular Weight                                                   
      REAL HG_MOLWT 
      PARAMETER (HG_MOLWT = 200.6) 
! Pressure in atmospheres                                               
      REAL PRES_ATM 
                                                                        
! Conv. factor for rate constants (from (mol/cc)-1 s-1 to ppm-1 s-1)    
      REAL CFACT 
                                                                        
! Local rate constants (ppm-1 s-1 units)                                
      REAL RKO3L, RKOHL, RKH2O2L, RKCL2L, RKHCLL 
                                                                        
! Effective 1st order rate constant for Hg(0) -> Hg(2)                  
      REAL RKHG0L 
                                                                        
! Decrease in Hg(0) conc. due to chemistry                              
      REAL DHG0 
!                                                                       
!-----Vertical profiles of HCl and Cl2 (ppm) used by the Hg chemistry   
!                                                                       
      INTEGER NHTCL 
!                                                                       
      PARAMETER (NHTCL = 6) 
!                                                                       
      REAL htcl (NHTCL) 
      REAL hclprof (NHTCL) 
      REAL cl2day (NHTCL) 
      REAL cl2nite (NHTCL) 
      DATA htcl / 0.06, 0.15, 0.45, 0.85, 2.00, 6.00 / 
      DATA hclprof / 1.0E-03, 1.0E-03, 1.0E-03, 1.0E-03, 1.0E-03,       &
      1.0E-03 /                                                         
      DATA cl2day / 1.0E-06, 1.0E-06, 1.0E-06, 1.0E-06, 1.0E-06,        &
      1.0E-06 /                                                         
      DATA cl2nite / 1.5E-04, 1.0E-04, 7.5E-05, 5.0E-05, 5.0E-05,       &
      5.0E-05 /                                                         
!     to judge if it is daytime or night time                           
      CALL ZenithAngle00 (tmid_sec, lon, lat, cos_sa) 
                                      ! daytime                         
      IF (cos_sa.ge.cos_sza_cut) then 
         ifdaytime = 1 
            !nighttime                                                  
      ELSE 
         ifdaytime = 0 
      ENDIF 
!                                                                       
!-----Assign HCl and Cl2 concentrations based on height(km), ocean      
!     and day/night                                                     
!                                                                       
      DO ihg = 1, NHTCL 
      IF ( (heiz / 1000.) .le.htcl (ihg) ) then 
         iht = ihg 
         GOTO 115 
      ENDIF 
      enddo 
                           !! layer altitude higher than 6.0 km         
      iht = NHTCL 
  115 CONTINUE 
      CHCL = hclprof (iht) 
      IF (LANDMASK.eq.0) then 
         IF (ifdaytime.eq.0) then 
            CCL2 = cl2nite (iht) 
         ELSE 
            CCL2 = cl2day (iht) 
         ENDIF 
      ELSE 
         CCL2 = 0.0 
      ENDIF 
!-----------------------------------------------------------------------
! Convert 2nd-order rate constants from mol/cc-s units to ppm-s units   
      PRES_ATM = PRES / STDATMMB 
      CFACT = COEF1 * PRES_ATM / TEMP 
                                                                        
      RKO3L = RKO3 * CFACT 
      RKOHL = RKOH * CFACT 
      RKH2O2L = RKH2O2 * CFACT 
      RKCL2L = RKCL2 * CFACT 
      RKHCLL = RKHCL * CFACT 
! convert conc unit, ng/m3 to ppm, by chenhs                            
      CHG0 = CHG0 * 1.0E-3 * (0.08206 * TEMP) / PRES_ATM / HG_MOLWT *   &
      1.0E-3                                                            
      CHG2 = CHG2 * 1.0E-3 * (0.08206 * TEMP) / PRES_ATM / HG_MOLWT *   &
      1.0E-3                                                            
! convert conc unit, ppb to ppm, by chenhs                              
      CO3 = CO3 * 1.0E-3 
      CH2O2 = CH2O2 * 1.0E-3 
      COH = COH * 1.0E-3 
                                                                        
! Calculate effective 1st-order rate constant for Hg(0) to Hg(2) oxidati
! assuming that oxidant species concentrations are constant             
      RKHG0L = (RKO3L * CO3 + RKOHL * COH + RKH2O2L * CH2O2 + RKCL2L *  &
      CCL2 + RKHCLL * CHCL)                                             
                                                                        
! Decrease in Hg(0) conc. due to gas-phase oxidation !! chenhs change   
      IF (lat.ge.15.and.lat.le.55.and.lon.ge. - 130.and.lon.le. - 55)   &
      then                                                              
         DHG0 = 0.2 * CHG0 * (1. - EXP ( - RKHG0L * DELT) ) 
      ELSE 
         DHG0 = 0.1 * CHG0 * (1. - EXP ( - RKHG0L * DELT) ) 
      ENDIF 
                                                                        
      CHG0 = CHG0 - DHG0 
      CHG2 = CHG2 + DHG0 
! The lower bound of HG0 and HG2 conc 1.0E-12 ppm                       
      CHG0 = max (CHG0, 1.0E-12) 
      CHG2 = max (CHG2, 1.0E-12) 
! convert conc unit, ppm to ng/m3, by chenhs                            
      CHG0 = CHG0 * 1.0E+3 * HG_MOLWT * PRES_ATM / (0.08206 * TEMP)     &
      * 1.0E+3                                                          
      CHG2 = CHG2 * 1.0E+3 * HG_MOLWT * PRES_ATM / (0.08206 * TEMP)     &
      * 1.0E+3                                                          
                                                                        
      RETURN 
      END SUBROUTINE HGGASCHEM                      
                                                                        
      SUBROUTINE ZenithAngle00 (tmid_sec, lon, lat, cos_sa) 
      REAL ::tmid_sec, rlon, rlat, cos_sa, lon, lat 
      REAL ::tlocal, tdec, codec, sidec, tloc, thou 
      DATA deg2rad / 0.017453293 / 
      rlon = lon * deg2rad 
      rlat = lat * deg2rad 
      tlocal = tmid_sec 
      tdec = 0.4092797 * sin (1.992385E-7 * tlocal) 
      sidec = sin (tdec) 
      codec = cos (tdec) 
      tloc = 7.272205E-5 * tlocal 
      thou = cos (rlon + tloc) 
      cos_sa = sin (rlat) * sidec + cos (rlat) * codec * thou 
      RETURN 
      END SUBROUTINE ZenithAngle00                        
end module hgchem      
