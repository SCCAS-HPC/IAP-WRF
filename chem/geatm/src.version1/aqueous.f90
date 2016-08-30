module aqueous_chem
contains
       SUBROUTINE AQUEOUS(T,PLEV,QVAPOR,CLW,CGAS, CAER, CPH, DT, SX,EX,SY,EY, I,J,K)
        INTEGER :: SX,EX,SY,EY,I,J,K
        REAL    :: TCELL,PCELL,TCHM,WCHM,DT
        REAL    :: CLDPH ! CLOUD PH VALUE
        REAL    :: H2  ! hydrogen concentration (ppm)
        REAL    :: O2  ! oxygen concentration (ppm)
        REAL    :: CH4 ! methane concentration (ppm)
        REAL    :: atm ! total gas concentration (M) in ppm
        REAL    :: CLIQ 
        REAL    :: convfac ! UNIT CHANGES FROM PPM to UMOL/M3
        REAL    :: TAMIN   ! Cloud water freezing threshold (K)
        REAL    :: CWMIN   ! Minimum cloud water threshold (g/m3)
        REAL    :: RH,PRES_PA, CW_KGM3
        REAL    :: R_GAS(11),R_AER(9)
        REAL    :: T,PLEV,CPH ! CPH IS CLOUD PH
        REAL    :: QVAPOR,WATER  ! WATER VAPOR IN KG/KG AND PPM
        REAL    :: CWC           ! cloud water content (g/m3):
        REAL    :: CLW           ! cloud water content IN KG/KG
        REAL    :: CGAS(11)
        REAL    :: CAER(9)     ! THE INPUT GAS AND AEROSOLS CONCENTRATIONS IN PPB AND UG/M3
        REAL    :: MAXBLD, MINBLD 

        DATA CWMIN / 0.05 /
        DATA TAMIN / 243./
! CGASNAME /SO2, HNO3, NXOY, CO2, NH3, HH2O2, O3, FOA, MHP, PAA, H2SO4/
! CAERNAME / ASO4, ANH4, ANO3, ACACO3, AMGCO3, ANACL, AFEIII, AMNII,POTCL/

          TCELL =    T 
          PCELL = PLEV ! HPA
          TCHM  = TCELL
          WATER = QVAPOR * 29. /18. * 1.E06 ! FROM KG/KG TO PPM
          WCHM  = WATER
                    
          O2  = 2.095e5
          CH4 = 1.75
          H2  = 0.60
          ATM = 1.E06
          convfac = 44.9 * (273./TCELL)*(PCELL/1013.)

! ****   FROM KG/KG TO G/M3
          
         !CWC  = CLW  * 29. * convfac 
         CWC  = CLW   !! by chenhs
          
         CLIQ = CWC 

         ! by chenhs
         !IF ( TCELL.LT.273) THEN
         ! CLIQ = AMAX1 (0., CLIQ* (TCELL - TAMIN) / (273. - TAMIN) )
         !ENDIF
 
           
!CCCCCCCCCCCCCC  DO RADM AQUEOUS CHEMISTRY IF CWC IS ABOVE THRESHOLD 
!CCCCCCCCCC    ALL CONC UNITS MYST BE MOL/MOL (MIXING RATIO)            
         IF ( CLIQ.GE.CWMIN .AND. TCELL.GE.TAMIN ) THEN
           PRES_PA = 100. * PCELL
           CW_KGM3 = CLIQ/1000.

           DO IG = 1 ,11
            R_GAS(IG) = CGAS(IG)*1.E-9 
           ENDDO
!            IF(I.EQ.16.AND.J.EQ.68) PRINT*, CGAS(16,68,4),R_GAS(4)
           R_AER(1) = (CAER(1)/96./CONVFAC)*1.E-6
           R_AER(2) = (CAER(2)/18./CONVFAC)*1.E-6
           R_AER(3) = (CAER(3)/62./CONVFAC)*1.E-6  
           R_AER(4) = (CAER(4)/100./CONVFAC)*1.E-6
           R_AER(5) = (CAER(5)/84./CONVFAC)*1.E-6
           R_AER(6) = (CAER(6)/58./CONVFAC)*1.E-6
           R_AER(7) = (CAER(7)/56./CONVFAC)*1.E-6
           R_AER(8) = (CAER(8)/55./CONVFAC)*1.E-6
           R_AER(9) = (CAER(9)/74./CONVFAC)*1.E-6
            
           CLDPH = CPH 

           IF(CLDPH.EQ.-1.E20) CLDPH = 5.0
              
            
           IF(R_AER(1).LT.5.0E-7.AND.R_AER(3).LT.1.5E-6) THEN

!             IF(I.EQ.46.AND.J.EQ.44.AND.K.EQ.10) PRINT*,CLDPH,'11111'
                                !R_AER(1),R_AER(2),R_AER(3),'111'
!            PRINT*, I,J,K,CLDPH,'111'
            CALL RAQCHEM(TCELL,PRES_PA,DT,CW_KGM3,R_GAS,R_AER,CLDPH,I,J,K) 
!             IF(I.EQ.46.AND.J.EQ.44.AND.K.EQ.10) PRINT*,CLDPH,R_AER(1),R_AER(2),R_AER(3),'222'     
!            PRINT*,I,J,K,CLDPH,'222'    
           ENDIF
          

           DO IG = 1, 11
            MINBLD =  CGAS(IG)*0.8
            MAXBLD =  CGAS(IG)*1.2

            IF(R_GAS (IG)* 1.E9 .GE. MAXBLD) THEN
             CGAS(IG) = MAXBLD
            ELSE  IF(R_GAS (IG)*1.E9 .LT.MINBLD) THEN
             CGAS(IG) = MINBLD
            ELSE
             CGAS(IG) = R_GAS (IG)* 1.E9
            ENDIF
             CGAS(IG) = AMAX1(CGAS(IG),1.E-20)

           ENDDO

          

!           CAER(I,J,1) = AMAX1(R_AER(1)* CONVFAC*96.*1.E6, 1.E-20) ! PSO4 UG/M3
!           CAER(I,J,2) = AMAX1(R_AER(2)* CONVFAC*18.*1.E6, 1.E-20) ! PNH4 UG/M3
!           CAER(I,J,3) = AMAX1(R_AER(3)* CONVFAC*62.*1.E6, 1.E-20) ! PNO3 UG/M3
!           CAER(I,J,6) = AMAX1(R_AER(6)* CONVFAC*58.*1.E6, 1.E-20) ! PCL  UG/M3
          
            CPH   = CLDPH
         ELSE
            CPH   = -1.E20  !! by chenhs
         ENDIF  ! CLIQ

         

      RETURN
      END SUBROUTINE

      SUBROUTINE RAQCHEM (TEMP, PRES_PA, TAUCLD, WCAVG, GAS, AEROSOL,   &
      BB, II, JJ, KK)                                                   
                                                                        
!amx                                                                    
!  This is a stand-alone version of the aqueous phase chemistry from    
!  CMAQ version 4 (November 2002). Deposition calculations are          
!  commented out. Aerosol concentrations by size mode have been replaced
!  by total concentrations.                                             
!amx                                                                    
!      SUBROUTINE AQCHEM ( JDATE, JTIME, TEMP, PRES_PA, TAUCLD, PRCRATE,
!     &                    WCAVG, WTAVG, AIRM, ALFA0, ALFA2, ALFA3, GAS,
!     &                    AEROSOL, GASWDEP, AERWDEP, HPWDEP )          
!-----------------------------------------------------------------------
!                                                                       
!  DESCRIPTION:                                                         
!    Compute concentration changes in cloud due to aqueous chemistry,   
!    scavenging and wet deposition amounts.                             
!                                                                       
!  Revision History:                                                    
!      No   Date   Who	What                                             
!      -- -------- ---  -----------------------------------------       
!      0  / /86    CW   BEGIN PROGRAM - Walceks's Original Code         
!      1  / /86    RB   INCORPORATE INTO RADM                           
!      2  03/23/87 DH   REFORMAT                                        
!      3  04/11/88 SJR  STREAMLINED CODE - ADDED COMMENTS               
!      4  08/27/88 SJR  COMMENTS, MODIFIED FOR RPM                      
!      4a 03/15/96 FSB  Scanned hard copy to develop Models3            
!                       Version.                                        
!      5  04/24/96 FSB  Made into Models3 Format                        
!      6  02/18/97 SJR  Revisions to link with Models3                  
!      7  08/12/97 SJR  Revised for new concentration units (moles/mole)
!                       and new treatment of nitrate and nitric acid    
!      8  01/15/98 sjr  revised to add new aitken mode scavenging       
!                       and aerosol number scavenging                   
!      9  12/15/98 David Wong at LM:                                    
!             -- change division of XL, TEMP to multiplication of XL, TE
!                reciprocal, respectively                               
!             -- change / TOTOX / TSIV to / ( TOTOX * TSIV )            
!     10  03/18/99 David Wong at LM:                                    
!             -- removed "* 1.0" redundant calculation at TEMP1 calculat
!     11  04/27/00 sjr  Added aerosol surface area as modeled species   
!                                                                       
!  Reference:                                                           
!     Walcek & Taylor, 1986, A theoretical Method for computing         
!      vertical distributions of acidity and sulfate within cumulus     
!      clouds, J. Atmos Sci.,  Vol. 43, no. 4 pp 339 - 355              
!                                                                       
!  Called by:  AQMAP                                                    
!                                                                       
!  Calls the following subroutines:  none                               
!                                                                       
!  Calls the following functions:  HLCONST                              
!                                                                       
!  ARGUMENTS     TYPE      I/O       DESCRIPTION                        
!  ---------     ----  ------------  --------------------------------   
!  GAS(ngas)     real  input&output  Concentration for species i=1,11   
!  GASWDEP(ngas) real     output     wet deposition for species         
!                                    (1) = SO2   conc (mol/mol of S02)  
!                                    (2) = HNO3  conc (mol/mol of HNO3) 
!                                    (3) = N2O5  conc (mol/mol of N2O5) 
!                                    (4) = CO2   conc (mol/mol of CO2)  
!                                    (5) = NH3   conc (mol/mol of NH3)  
!                                    (6) = H2O2  conc (mol/mol of H2O2) 
!                                    (7) = O3    conc (mol/mol of O3)   
!                                    (8) = FOA   conc (mol/mol of FOA)  
!                                    (9) = MHP   conc (mol/mol of MHP)  
!                                    (10)= PAA   conc (mol/mol of PAA)  
!                                    (11)= H2SO4 conc (mol/mol of H2SO4)
!                                                                       
!amx                                                                    
!                                                                       
!  AEROSOL concentrations for species i=1,9 in camx version             
!                                    (1) = SO4    conc (mol/mol)        
!                                    (2) = NH4    conc (mol/mol)        
!                                    (3) = NO3    conc (mol/mol)        
!                                    (4) = CAC    conc (mol/mol)        
!                                    (5) = MGC    conc (mol/mol)        
!                                    (6) = NACL   conc (mol/mol)        
!                                    (7) = A3FE   conc (mol/mol)        
!                                    (8) = B2MN   conc (mol/mol)        
!                                    (9) = KCL    conc (mol/mol)        
!amx                                                                    
!  AEROSOL(naer) real input&output   Concentration for species i=1,21   
!amx                                                                    
!  BB            real     output     cloud water pH                     
!amx                                                                    
!  AERWDEP(naer) real     output     wet deposition for species         
!                                    (1) = SO4AKN conc (mol/mol)        
!                                    (2) = SO4ACC conc (mol/mol)        
!                                    (3) = NH4AKN conc (mol/mol)        
!                                    (4) = NH4ACC conc (mol/mol)        
!                                    (5) = NO3AKN conc (mol/mol)        
!                                    (6) = NO3ACC conc (mol/mol)        
!                                    (7) = NO3COR conc (mol/mol)        
!                                    (8) = ORGAKN conc (mol/mol)        
!                                    (9) = ORGACC conc (mol/mol)        
!                                    (10)= PRIAKN conc (mol/mol)        
!                                    (11)= PRIACC conc (mol/mol)        
!                                    (12)= PRICOR conc (mol/mol)        
!                                    (13)= CACO3  conc (mol/mol)        
!                                    (14)= MGCO3  conc (mol/mol)        
!                                    (15)= NACL   conc (mol/mol)        
!                                    (16)= A3FE   conc (mol/mol)        
!                                    (17)= B2MN   conc (mol/mol)        
!                                    (18)= KCL    conc (mol/mol)        
!                                    (19)= NUMAKN conc (#/mol)          
!                                    (20)= NUMACC conc (#/mol)          
!                                    (21)= NUMCOR conc (#/mol)          
!                                    (22)= SRFAKN conc (m2/mol)         
!                                    (23)= SRFACC conc (m2/mol)         
!                                                                       
!-----------------------------------------------------------------------
                                                                        
      IMPLICIT NONE 
                                                                        
!      INCLUDE SUBST_CONST          ! constants                         
!      INCLUDE SUBST_XSTAT          ! M3EXIT status codes               
!      INCLUDE 'AQ_PARAMS.EXT'      ! aqueous chemistry shared parameter
!amx                                                                    
! Essential parameters from the include files                           
!                                                                       
                            ! local number of gasses for aqchem         
      INTEGER NGAS 
                            ! local number of aerosols for aqchem       
      INTEGER NAER 
                            ! Molar volume at STP [ L/mol ] Non MKS unit
      REAL MOLVOL 
                            ! standard atmosphere  [ Pa ]               
      REAL STDATMPA 
                            ! Standard Temperature [ K ]                
      REAL STDTEMP 
      PARAMETER (NGAS = 11) 
      PARAMETER (NAER = 9) 
      PARAMETER (MOLVOL = 22.41410) 
      PARAMETER (STDATMPA = 101325.0) 
      PARAMETER (STDTEMP = 273.15) 
!amx                                                                    
!      CHARACTER*120 XMSG           ! Exit status message               
!      DATA          XMSG / ' ' /                                       
                                                                        
!...........PARAMETERS and their descriptions:                          
                                                                        
                                   ! number of oxidizing reactions      
      INTEGER NUMOX 
      PARAMETER (NUMOX = 5) 
                                                                        
                                   ! density of water at 20 C and 1 ATM 
      REAL H2ODENS 
                                       ! (kg/m3)                        
      PARAMETER (H2ODENS = 1000.0) 
                                                                        
                                   ! number of liquid phase species     
      INTEGER NLIQS 
      PARAMETER (NLIQS = 33) 
                                                                        
                                  ! 1/3                                 
      REAL ONETHIRD 
      PARAMETER (ONETHIRD = 1.0 / 3.0) 
                                                                        
                                   ! 2/3                                
      REAL TWOTHIRDS 
      PARAMETER (TWOTHIRDS = 2.0 / 3.0) 
!amx                                                                    
!........... Gas species pointers                                       
      INTEGER lso2 
      INTEGER lhno3 
      INTEGER ln2o5 
      INTEGER lco2 
      INTEGER lnh3 
      INTEGER lh2o2 
      INTEGER lo3 
      INTEGER lfoa 
      INTEGER lmhp 
      INTEGER lpaa 
      INTEGER lh2so4 
      PARAMETER (lso2 = 1) 
      PARAMETER (lhno3 = 2) 
      PARAMETER (ln2o5 = 3) 
      PARAMETER (lco2 = 4) 
      PARAMETER (lnh3 = 5) 
      PARAMETER (lh2o2 = 6) 
      PARAMETER (lo3 = 7) 
      PARAMETER (lfoa = 8) 
      PARAMETER (lmhp = 9) 
      PARAMETER (lpaa = 10) 
      PARAMETER (lh2so4 = 11) 
!........... Aerosol species pointers                                   
      INTEGER lso4 
      INTEGER lnh4 
      INTEGER lno3 
      INTEGER lcaco3 
      INTEGER lmgco3 
      INTEGER lnacl 
      INTEGER la3fe 
      INTEGER lb2mn 
      INTEGER lkcl 
      PARAMETER (lso4 = 1) 
      PARAMETER (lnh4 = 2) 
      PARAMETER (lno3 = 3) 
      PARAMETER (lcaco3 = 4) 
      PARAMETER (lmgco3 = 5) 
      PARAMETER (lnacl = 6) 
      PARAMETER (la3fe = 7) 
      PARAMETER (lb2mn = 8) 
      PARAMETER (lkcl = 9) 
!amx                                                                    
!...........ARGUMENTS and their descriptions                            
                                                                        
!      INTEGER      JDATE           ! current model date, coded YYYYDDD 
!      INTEGER      JTIME           ! current model time, coded HHMMSS  
                                                                        
!      REAL         AIRM            ! total air mass in cloudy layers (m
!      REAL         ALFA0           ! scav coef for aitken aerosol numbe
!      REAL         ALFA2           ! scav coef for aitken aerosol sfc a
!      REAL         ALFA3           ! scav coef for aitken aerosol mass 
!      REAL         HPWDEP          ! hydrogen wet deposition (mm mol/li
!      REAL         PRCRATE         ! precip rate (mm/hr)               
                                   ! pressure (Pa)                      
      REAL PRES_PA 
                                   ! timestep for cloud (s)             
      REAL TAUCLD 
                                   ! temperature (K)                    
      REAL TEMP 
                                   ! liquid water content (kg/m3)       
      REAL WCAVG 
!      REAL         WTAVG           ! total water content (kg/m3)       
                                   ! gas phase concentrations (mol/molV)
      REAL GAS (NGAS) 
                                   ! aerosol concentrations (mol/molV)  
      REAL AEROSOL (NAER) 
!      REAL         GASWDEP( NGAS ) ! gas phase wet deposition array (mm
!      REAL         AERWDEP( NAER ) ! aerosol wet deposition array (mm m
                                   ! unit number for diagnostic output  
      INTEGER idiag 
                                   ! unit number for message output     
      INTEGER iout 
                                   ! grid number                        
      INTEGER igrd 
                                   ! i grid cell index (column)         
      INTEGER iaq 
                                   ! j grid cell index (row)            
      INTEGER jaq 
                                   ! k grid cell index (layer)          
      INTEGER kaq 
                                                                        
!...........LOCAL VARIABLES (scalars) and their descriptions:           
                                                                        
                                  ! driver program name                 
      CHARACTER(16) PNAME 
      SAVE PNAME 
      DATA PNAME / 'AQCHEM' / 
                                                                        
                                   ! loop counter for do loop 20        
      INTEGER I20C 
                                   ! loop counter for do loop 30        
      INTEGER I30C 
                                   ! # iterations of aqueaous chemistry 
      INTEGER ITERAT 
                                   ! aqueous chem iteration counter     
      INTEGER I7777C 
                                   ! aqueous chem iteration counter     
      INTEGER ICNTAQ 
!      INTEGER      LIQ             ! loop counter for liquid species   
                                   ! index over oxidation reactions     
      INTEGER IOX 
                                                                        
!      REAL         DEPSUM                                              
!      REAL         BETASO4                                             
                                   ! iron's anion concentration         
      REAL A 
                                   ! H+ concentration in cloudwater (mol
      REAL AC 
                                   ! activity corretion factor!single io
      REAL ACT1 
                                   ! activity factor correction!double i
      REAL ACT2 
                                   !                                    
      REAL ACTB 
                                   ! guess for H+ conc in cloudwater (mo
      REAL AE 
                                   ! manganese's anion concentration    
      REAL B 
                                   ! pressure (Atm)                     
      REAL PRES_ATM 
                                   ! lower limit guess of cloudwater pH 
      REAL BB 
                                   ! Calcium conc in cloudwater (mol/lit
      REAL CA 
                                   ! inital Calcium in cloudwater (mol/l
      REAL CAA 
!      REAL         NO3CORA         ! initial NO3COR in cloudwater (mol/
                                   ! Cl-  conc in cloudwater (mol/liter)
      REAL CL 
                                   ! initial Cl in cloudwater (mol/liter
      REAL CLA 
                                   ! Henry's Law constant for CO2       
      REAL CO2H 
                                   ! First dissociation constant for CO2
      REAL CO21 
                                   ! Second dissociation constant for CO
      REAL CO22 
                                   ! CO21*CO22                          
      REAL CO212 
                                   ! CO2H*CO21*CO22                     
      REAL CO212H 
                                   ! CO2H*CO21                          
      REAL CO21H 
                                   ! CO2 conc in cloudwater (mol/liter) 
      REAL CO2L 
                                   ! CO3= conc in cloudwater (mol/liter)
      REAL CO3 
                                   ! initial CO3 in cloudwater (mol/lite
      REAL CO3A 
!      REAL         CTHK1           ! cloud thickness (m)               
                                   !                                    
      REAL DTRMV 
                                   !                                    
      REAL DTS6 
                                   ! functional value ??                
      REAL FA 
                                   ! functional value ??                
      REAL FB 
                                   ! Fe+++ conc in cloudwater (mol/liter
      REAL FE 
                                   ! initial Fe in cloudwater (mol/liter
      REAL FEA 
                                   ! frac weight of NH3 to total ammonia
      REAL FNH3 
!      REAL         FNH4ACC         ! frac weight of NH4 acc to total am
                                   ! frac weight of HNO3 to total NO3   
      REAL FHNO3 
!      REAL         FNO3ACC         ! frac weight of NO3 acc to total NO
!      REAL         FNO3COR         ! frac weight of NO3 cor to total NO
!      REAL         FRACACC         ! frac ACC that was from accum mode 
!      REAL         FRACCOR         ! frac NO3 that was from coarse mode
                                   ! First dissociation constant for FOA
      REAL FOA1 
                                   ! Henry's Law constant for FOA       
      REAL FOAH 
                                   ! FOAH*FOA1                          
      REAL FOA1H 
                                   ! FOA conc in cloudwater (mol/liter) 
      REAL FOAL 
                                   !                                    
      REAL FTST 
                                   !                                    
      REAL GM 
                                   !                                    
      REAL GM1 
                                   !                                    
      REAL GM1LOG 
                                   ! activity correction factor         
      REAL GM2 
                                   !                                    
      REAL GM2LOG 
                                   !                                    
      REAL HA 
                                   !                                    
      REAL HB 
                                   !                                    
      REAL H2OW 
                                   ! Henry's Law Constant for H2O2      
      REAL H2O2H 
                                   ! H2O2 conc in cloudwater (mol/liter)
      REAL H2O2L 
                                   ! HCO2 conc in cloudwater (mol/liter)
      REAL HCO2 
                                   ! HCO3 conc in cloudwater (mol/liter)
      REAL HCO3 
                                   ! Henry's Law Constant for HNO3      
      REAL HNO3H 
                                   ! First dissociation constant for HNO
      REAL HNO31 
                                   !                                    
      REAL HNO31H 
                                   ! HNO3 conc in cloudwater (mol/liter)
      REAL HNO3L 
                                   ! HSO3 conc in cloudwater (mol/liter)
      REAL HSO3 
                                   ! HSO4 concn in cloudwater (mol/liter
      REAL HSO4 
                                   !                                    
      REAL HTST 
                                   ! K conc in cloudwater (mol/liter)   
      REAL K 
                                   ! initial K in cloudwater (mol/liter)
      REAL KA 
                                   ! log of TEMP                        
      REAL LGTEMP 
!      REAL         M3NEW           ! accumulation mode mass at time t  
!      REAL         M3OLD           ! accumulation mode mass at time 0  
                                   !                                    
      REAL MG 
                                   ! inital Mg in cloudwater (mol/liter)
      REAL MGA 
                                   ! Henry's Law Constant for MHP       
      REAL MHPH 
                                   ! MHP conc in cloudwater (mol/liter) 
      REAL MHPL 
                                   ! Mn++ conc in cloudwater (mol/liter)
      REAL MN 
                                   ! initial Mn in cloudwater (mol/liter
      REAL MNA 
                                   ! Na conc in cloudwater (mol/liter)  
      REAL NA 
                                   ! initial Na in cloudwater (mol/liter
      REAL NAA 
                                   ! First dissociation constant for NH3
      REAL NH31 
                                   ! Henry's Law Constant for NH3       
      REAL NH3H 
                                   !                                    
      REAL NH3DH20 
                                   !                                    
      REAL NH31HDH 
                                   ! NH3 conc in cloudwater (mol/liter) 
      REAL NH3L 
                                   ! NH4+ conc in cloudwater (mol/liter)
      REAL NH4 
!      REAL         NH4AKNA         ! init NH4 akn conc in cloudwater (m
                                   ! init NH4 acc conc in cloudwater (mo
      REAL NH4ACCA 
!      REAL         NITAER          ! total aerosol nitrate             
                                   ! NO3 conc in cloudwater (mol/liter) 
      REAL NO3 
                                   ! init NO3 acc conc in cloudwater (mo
      REAL NO3ACCA 
!      REAL         NO3AKNA         ! init NO3 akn conc in cloudwater (m
                                   ! Henry's Law Constant for O3        
      REAL O3H 
                                   ! O3 conc in cloudwater (mol/liter)  
      REAL O3L 
                                   ! OH conc in cloudwater (mol/liter)  
      REAL OH 
!      REAL         ORGN            ! ORGANIC aerosol in cloudwater (mol
!      REAL         ORGACCA         ! init ORG ACC aerosol in cloudwater
!      REAL         ORGAKNA         ! init ORG AKN aerosol in cloudwater
                                   ! Henry's Law Constant for PAA       
      REAL PAAH 
                                   ! PAA conc in cloudwater (mol/liter) 
      REAL PAAL 
                                   ! total CO2 partial pressure (atm)   
      REAL PCO20 
                                   ! gas only CO2 partial pressure (atm)
      REAL PCO2F 
                                   ! total ORGANIC acid partial pressure
      REAL PFOA0 
                                   ! gas only ORGANIC ACID partial press
      REAL PFOAF 
                                   ! total H2O2 partial pressure (atm)  
      REAL PH2O20 
                                   ! gas only H2O2 partial pressure (atm
      REAL PH2O2F 
                                   ! total HNO3 partial pressure (atm)  
      REAL PHNO30 
                                   ! gas only HNO3 partial pressure (atm
      REAL PHNO3F 
                                   ! total MHP partial pressure (atm)   
      REAL PMHP0 
                                   ! gas only MHP partial pressure (atm)
      REAL PMHPF 
                                   ! total NH3 partial pressure (atm)   
      REAL PNH30 
                                   ! gas only NH3 partial pressure (atm)
      REAL PNH3F 
                                   ! total O3 partial pressure (atm)    
      REAL PO30 
                                   ! gas only O3 partial pressure (atm) 
      REAL PO3F 
                                   ! total PAA partial pressure (atm)   
      REAL PPAA0 
                                   ! gas only PAA partial pressure (atm)
      REAL PPAAF 
!      REAL         PRIM            ! PRIMARY acc+akn aerosol in cloudwa
!      REAL         PRIMCOR         ! PRIMARY coarse aerosol in cloudwat
!      REAL         PRIACCA         ! init PRI ACC aerosol in cloudwater
!      REAL         PRIAKNA         ! init PRI AKN aerosol in cloudwater
!      REAL         PRICORA         ! init PRI COR aerosol in cloudwater
                                   ! total SO2 partial pressure (atm)   
      REAL PSO20 
                                   ! gas only SO2 partial pressure (atm)
      REAL PSO2F 
                                   !                                    
      REAL RATE 
                                   !                                    
      REAL RECIPA1 
                                   !                                    
      REAL RECIPA2 
                                   ! one over pressure (/atm)           
      REAL RECIPAP1 
                                   !                                    
      REAL RH2O2 
                                   !                                    
      REAL RMHP 
                                   !                                    
      REAL RPAA 
                                   ! gas const * temperature (liter atm/
      REAL RT 
                                   ! Scavenging efficiency (%)          
      REAL SCVEFF 
      SAVE SCVEFF 
                                    ! currently set to 100%             
      DATA SCVEFF / 100.0 / 
                                   ! dissolved so2 in cloudwater (mol/li
      REAL SIV 
                                   !                                    
      REAL SK6 
                                   !                                    
      REAL SK6TS6 
                                   ! First dissociation constant for SO2
      REAL SO21 
                                   ! Second dissociation constant for SO
      REAL SO22 
                                   ! Henry's Law Constant for SO2       
      REAL SO2H 
                                   ! SO21*SO22                          
      REAL SO212 
                                   ! SO21*SO22*SO2H                     
      REAL SO212H 
                                   ! SO21*SO2H                          
      REAL SO21H 
                                   ! SO2 conc in cloudwater (mol/liter) 
      REAL SO2L 
                                   ! SO3= conc in cloudwater (mol/liter)
      REAL SO3 
                                   ! SO4= conc in cloudwater (mol/liter)
      REAL SO4 
                                   ! ionic strength                     
      REAL STION 
                                   !                                    
      REAL TAC 
                                   !                                    
      REAL TEMP1 
                                   ! cloud chemistry clock (sec)        
      REAL TIMEW 
                                   !                                    
      REAL TOTOX 
                                   ! total ammonium                     
      REAL TOTAMM 
                                   ! total nitrate                      
      REAL TOTNIT 
                                   ! SO4 conc in cloudwater (mol/liter) 
      REAL TS6 
!      REAL         TS6AKNA         ! init SO4 akn conc in cloudwater (m
                                   ! init SO4 acc conc in cloudwater (mo
      REAL TS6ACCA 
                                   !                                    
      REAL TSIV 
                                   !                                    
      REAL TST 
!      REAL         XC1             ! (/mm)                             
!      REAL         XC2             ! (liter-atm/mol/mm)                
                                   ! conversion factor (liter-atm/mol)  
      REAL XL 
                                   ! 1.0 / XL                           
      REAL ONE_OVER_XL 
                                        ! PRES_ATM / XL                 
      REAL PRES_ATM_OVER_XL 
                                   !                                    
      REAL XLCO2 
                                   !                                    
      REAL XLH2O2 
                                   !                                    
      REAL XLHNO3 
                                   !                                    
      REAL XLMHP 
                                   !                                    
      REAL XLNH3 
                                   !                                    
      REAL XLO3 
                                   !                                    
      REAL XLPAA 
                                   !                                    
      REAL XLSO2 
                                                                        
!...........LOCAL VARIABLES (arrays) and their descriptions:            
                                                                        
                                   ! liquid concentration array (mol/lit
      REAL LIQUID (NLIQS) 
!      REAL         WETDEP( NLIQS ) ! wet deposition array (mm mol/liter
                                     ! rate of so2 oxid incloud (mol/lit
      REAL DSIVDT (0:NUMOX) 
                                     ! S(IV) oxidized over timestep DTW(
      REAL DS4 (0:NUMOX) 
                                     ! cloud chemistry timestep (sec)   
      REAL DTW (0:NUMOX) 
                                                                        
                                     ! 1.0 / TEMP                       
      REAL ONE_OVER_TEMP 
      INTEGER II, JJ, KK 
!...........EXTERNAL FUNCTIONS and their descriptions:                  
                                                                        
!      REAL HLCONST 
                                                                        
!*********************************************************************  
!     begin body of subroutine AQCHEM                                   
                                                                        
!      IF(II.EQ.65.AND.JJ.EQ.11.AND.KK.EQ.5) PRINT*, 'AQUE STEP 1'      
                                                                        
      ONE_OVER_TEMP = 1.0 / TEMP 
                                                                        
!...check for bad temperature, cloud air mass, or pressure              
                                                                        
!      IF ( TEMP .LE. 0.0 ) THEN                                        
!        IF ( AIRM .LE. 0.0 ) THEN                                      
!          IF ( PRES_PA .LE. 0.0 ) THEN                                 
      IF (TEMP.LE.0.0.OR.PRES_PA.LE.0.0) THEN 
!            XMSG = 'MET DATA ERROR'                                    
      WRITE ( * ,  * ) 'Error in RADM aqueous chemistry (RAQCHEM) , Inva&
      lid               met data: T, P', TEMP, PRES_PA                  
!          END IF                                                       
!        END IF                                                         
      ENDIF 
                                                                        
!...compute several conversion factors                                  
                                                                        
      ICNTAQ = 0 
      ITERAT = 0 
                                                   ! R * T (liter atm / 
      RT = (MOLVOL / STDTEMP) * TEMP 
                                                   ! pressure (atm)     
      PRES_ATM = PRES_PA / STDATMPA 
!      CTHK1 = AIRM * RT / ( PRES_ATM * 1000.0 )    ! cloud thickness (m
                                      ! conversion factor (l-atm/mol)   
      XL = WCAVG * RT / H2ODENS 
      ONE_OVER_XL = 1.0 / XL 
      PRES_ATM_OVER_XL = PRES_ATM / XL 
      TST = 0.999 
      GM = SCVEFF / 100.0 
      ACT1 = 1.0 
      ACT2 = 1.0 
      GM2 = 1.0 
      TIMEW = 0.0 
      RECIPAP1 = 1.0 / PRES_ATM 
!      XC1  = 1.0 / ( WCAVG * CTHK1 )                                   
!      XC2  = RT / ( 1000.0 * CTHK1 )                                   
                                                                        
!...set equilibrium constants as a function of temperature              
!...   Henry's law constants                                            
                                                                        
      SO2H = HLCONST ('SO2             ', TEMP) 
      CO2H = HLCONST ('CO2             ', TEMP) 
      NH3H = HLCONST ('NH3             ', TEMP) 
      H2O2H = HLCONST ('H2O2            ', TEMP) 
      O3H = HLCONST ('O3              ', TEMP) 
      HNO3H = HLCONST ('HNO3            ', TEMP) 
      MHPH = HLCONST ('METHYLHYDROPEROX', TEMP) 
      PAAH = HLCONST ('PEROXYACETIC_ACI', TEMP) 
      FOAH = HLCONST ('FORMIC_ACID     ', TEMP) 
                                                                        
!...dissociation constants                                              
                                                                        
      FOA1 = 1.71E-4 
                                                                        
!...From ERT book                                                       
                                                                        
      SK6 = 10.0** (1180.0 * ONE_OVER_TEMP - 5.95) 
                                                                        
!...From Maahs <1982>                                                   
                                                                        
      SO21 = 10.0** (853.0 * ONE_OVER_TEMP) / 5.495E4 
      SO22 = 10.0** (621.9 * ONE_OVER_TEMP) / 1.897E9 
                                                                        
!...From Edwards et al. <1978>                                          
                                                                        
      LGTEMP = ALOG (TEMP) 
                                                                        
      CO21 = 10.0** ( - 5251.5 * ONE_OVER_TEMP - 15.94 * LGTEMP +       &
      102.269)                                                          
      CO22 = 10.0** ( - 5401.4 * ONE_OVER_TEMP - 15.4096 * LGTEMP +     &
      95.574)                                                           
      H2OW = 10.0** ( - 5839.5 * ONE_OVER_TEMP - 9.7618 * LGTEMP +      &
      61.206)                                                           
                                                                        
!...From Morgan and Maas <1931>                                         
                                                                        
      NH31 = 10.0** ( - 189.1 * ONE_OVER_TEMP - 4.117) 
                                                                        
!...Schwartz and White <1981, 1983>                                     
                                                                        
      HNO31 = 15.4 
                                                                        
!...Kinetic oxidation rates                                             
!...   From Chamedies (1982)                                            
                                                                        
      TEMP1 = ONE_OVER_TEMP - 1.0 / 298.0 
                                                                        
      RH2O2 = 80000.0 * EXP ( - 3650.0 * TEMP1) 
                                                                        
!...From Kok                                                            
                                                                        
      RMHP = 1.75E7 * EXP ( - 3801.0 * TEMP1) 
      RPAA = 3.64E7 * EXP ( - 3994.0 * TEMP1) 
                                                                        
!...make initializations                                                
                                                                        
!      DO LIQ = 1, NLIQS                                                
!        WETDEP( LIQ ) = 0.0                                            
!      END DO                                                           
                                                                        
      DO IOX = 0, NUMOX 
      DSIVDT (IOX) = 0.0 
      DTW (IOX) = 0.0 
      DS4 (IOX) = 0.0 
      ENDDO 
                                                                        
!...compute the initial accumulation aerosol 3rd moment                 
!                                                                       
!      M3OLD = ( AEROSOL( LSO4ACC ) * SGRAERMW( LSO4ACC ) / 1.8e6       
!     &      +   AEROSOL( LNH4ACC ) * SGRAERMW( LNH4ACC ) / 1.8e6       
!     &      +   AEROSOL( LNO3ACC ) * SGRAERMW( LNO3ACC ) / 1.8e6       
!     &      +   AEROSOL( LORGACC ) * SGRAERMW( LORGACC ) / 2.0e6       
!     &      +   AEROSOL( LPRIACC ) * SGRAERMW( LPRIACC ) / 2.2e6 )     
!cc     &      * 6.0 / PI    ! cancels out in division at end of subrout
!                                                                       
!...compute fractional weights for several species                      
!amx                                                                    
      FHNO3 = 1.0 
      TOTNIT = GAS (LHNO3) + AEROSOL (LNO3) 
      IF (TOTNIT.GT.0.0) FHNO3 = GAS (LHNO3) / TOTNIT 
!                                                                       
      FNH3 = 1.0 
      TOTAMM = GAS (LNH3) + AEROSOL (LNH4) 
      IF (TOTAMM.GT.0.0) FNH3 = GAS (LNH3) / TOTAMM 
!amx                                                                    
!                                                                       
!      NITAER = AEROSOL( LNO3ACC ) + AEROSOL( LNO3COR )                 
!      IF ( NITAER .GT. 0.0 ) THEN                                      
!        FRACACC = AEROSOL( LNO3ACC ) / NITAER                          
!        FRACCOR = AEROSOL( LNO3COR ) / NITAER                          
!      ELSE                                                             
!        FRACACC = 1.0                                                  
!        FRACCOR = 0.0                                                  
!      END IF                                                           
!                                                                       
!      TOTNIT = GAS( LHNO3 ) + AEROSOL( LNO3ACC ) + AEROSOL( LNO3COR )  
!      IF ( TOTNIT .GT. 0.0 ) THEN                                      
!        FHNO3   = GAS( LHNO3 ) / TOTNIT                                
!        FNO3ACC = AEROSOL( LNO3ACC ) / TOTNIT                          
!        FNO3COR = AEROSOL( LNO3COR ) / TOTNIT                          
!      ELSE                                                             
!        FHNO3   = 1.0                                                  
!        FNO3ACC = 0.0                                                  
!        FNO3COR = 0.0                                                  
!      END IF                                                           
!                                                                       
!      TOTAMM = GAS( LNH3 ) + AEROSOL( LNH4ACC )                        
!      IF ( TOTAMM .GT. 0.0 ) THEN                                      
!        FNH3    = GAS( LNH3 ) / TOTAMM                                 
!        FNH4ACC = AEROSOL( LNH4ACC ) / TOTAMM                          
!      ELSE                                                             
!        FNH3    = 1.0                                                  
!        FNH4ACC = 0.0                                                  
!      END IF                                                           
!                                                                       
!...initial concentration from accumulation-mode aerosol loading (mol/li
!...  an assumption is made that all of the accumulation-mode           
!...  aerosol mass in incorporated into the cloud droplets              
!amx                                                                    
      TS6ACCA = (AEROSOL (LSO4) + GAS (LH2SO4) ) * PRES_ATM_OVER_XL 
      NO3ACCA = AEROSOL (LNO3) * PRES_ATM_OVER_XL 
      NH4ACCA = AEROSOL (LNH4) * PRES_ATM_OVER_XL 
!amx                                                                    
!      TS6ACCA = ( AEROSOL( LSO4ACC )                                   
!     &        +   GAS    ( LH2SO4  ) ) * PRES_ATM_OVER_XL              
!      NO3ACCA =   AEROSOL( LNO3ACC )   * PRES_ATM_OVER_XL              
!      NH4ACCA =   AEROSOL( LNH4ACC )   * PRES_ATM_OVER_XL              
!      ORGACCA =   AEROSOL( LORGACC )   * PRES_ATM_OVER_XL              
!      PRIACCA =   AEROSOL( LPRIACC )   * PRES_ATM_OVER_XL              
                                                                        
!...initial concentration from coarse-mode aerosol loading (mol/liter)  
!...  an assumption is made that all of the coarse-mode                 
!...  aerosol mass in incorporated into the cloud droplets              
                                                                        
      CLA = (AEROSOL (LNACL) * PRES_ATM_OVER_XL) 
!     &        +   AEROSOL( LKCL    ) ) * PRES_ATM_OVER_XL              
!      NO3CORA =   AEROSOL( LNO3COR )   * PRES_ATM_OVER_XL              
      CAA = AEROSOL (LCACO3) * PRES_ATM_OVER_XL 
      MGA = AEROSOL (LMGCO3) * PRES_ATM_OVER_XL 
      NAA = AEROSOL (LNACL) * PRES_ATM_OVER_XL 
      KA = AEROSOL (LKCL) * PRES_ATM_OVER_XL 
      FEA = AEROSOL (LA3FE) * PRES_ATM_OVER_XL 
      MNA = AEROSOL (LB2MN) * PRES_ATM_OVER_XL 
      CO3A = (AEROSOL (LCACO3) + AEROSOL (LMGCO3) ) * PRES_ATM_OVER_XL 
!      PRICORA =   AEROSOL( LPRICOR )   * PRES_ATM_OVER_XL              
                                                                        
!...set constant factors that will be used in later multiplications (mol
                                                                        
      XLH2O2 = H2O2H * XL 
      XLO3 = O3H * XL 
      XLMHP = MHPH * XL 
      XLPAA = PAAH * XL 
      XLSO2 = SO2H * XL 
      XLNH3 = NH3H * XL 
      XLHNO3 = HNO3H * XL 
      XLCO2 = CO2H * XL 
                                                                        
      SO212 = SO21 * SO22 
      SO21H = SO21 * SO2H 
      SO212H = SO212 * SO2H 
      CO212 = CO21 * CO22 
      CO21H = CO21 * CO2H 
      CO212H = CO22 * CO21H 
      NH3DH20 = NH31 / H2OW 
      NH31HDH = NH3H * NH3DH20 
      FOA1H = FOA1 * FOAH 
      HNO31H = HNO31 * HNO3H 
                                                                        
!...If kinetic calculations are made, return to this point              
!      IF(II.EQ.65.AND.JJ.EQ.11.AND.KK.EQ.5) PRINT*, 'AQUE STEP 1_1'    
                                                                        
      I20C = 0 
   20 CONTINUE 
!!***** LIJIE                                                           
!      IF(II.EQ.65.AND.JJ.EQ.11.AND.KK.EQ.5) PRINT*, I20C,'AQUE STEP 1_2
                                                                        
      I20C = I20C + 1 
      IF (I20C.GE.1000) THEN 
!        XMSG = 'EXCESSIVE LOOPING AT I20C'                             
!        CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT2 )              
         WRITE ( * , * ) 'Error in RADM aqueous chemistry (AQCHEM)' 
         WRITE ( * , * ) 'EXCESSIVE LOOPING AT I20C' 
      ENDIF 
                                                                        
!...set aitken-mode aerosol loading (mol/liter)                         
!                                                                       
!      NO3AKNA = AEROSOL( LNO3AKN ) * PRES_ATM_OVER_XL                  
!     &        * ( 1.0 - EXP( -ALFA3 * TIMEW ) )                        
!      NH4AKNA = AEROSOL( LNH4AKN ) * PRES_ATM_OVER_XL                  
!     &        * ( 1.0 - EXP( -ALFA3 * TIMEW ) )                        
!      TS6AKNA = AEROSOL( LSO4AKN ) * PRES_ATM_OVER_XL                  
!     &        * ( 1.0 - EXP( -ALFA3 * TIMEW ) )                        
!      ORGAKNA = AEROSOL( LORGAKN ) * PRES_ATM_OVER_XL                  
!     &        * ( 1.0 - EXP( -ALFA3 * TIMEW ) )                        
!      PRIAKNA = AEROSOL( LPRIAKN ) * PRES_ATM_OVER_XL                  
!     &        * ( 1.0 - EXP( -ALFA3 * TIMEW ) )                        
!                                                                       
!...Initial gas phase partial pressures (atm)                           
!...   = initial partial pressure - amount deposited partial pressure   
                                                                        
      PSO20 = GAS (LSO2) * PRES_ATM + DS4 (0) * XL 
!     &       - ( WETDEP(  8 ) + WETDEP(  9 ) + WETDEP( 10 ) ) * XC2    
      PNH30 = GAS (LNH3) * PRES_ATM + (NH4ACCA) * XL 
!     &       + ( NH4ACCA + NH4AKNA ) * XL                              
!     &       - ( WETDEP(  2 ) + WETDEP( 15 ) ) * XC2                   
      PHNO30 = (GAS (LHNO3) + 2.0 * GAS (LN2O5) ) * PRES_ATM + (NO3ACCA)&
      * XL                                                              
!      PHNO30 = ( GAS( LHNO3 ) + 2.0 * GAS( LN2O5 ) ) * PRES_ATM        
!     &       + ( NO3ACCA + NO3CORA + NO3AKNA ) * XL                    
!     &       - ( WETDEP( 14 ) + WETDEP( 32 ) ) * XC2                   
      PH2O20 = GAS (LH2O2) * PRES_ATM 
      PO30 = GAS (LO3) * PRES_ATM 
!      PH2O20 = GAS( LH2O2 ) * PRES_ATM - WETDEP( 17 ) * XC2            
!      PO30   = GAS( LO3   ) * PRES_ATM - WETDEP( 18 ) * XC2            
      PFOA0 = GAS (LFOA) * PRES_ATM 
!     &       - ( WETDEP( 22 ) + WETDEP( 23 ) ) * XC2                   
      PMHP0 = GAS (LMHP) * PRES_ATM 
      PPAA0 = GAS (LPAA) * PRES_ATM 
!      PMHP0  = GAS( LMHP  ) * PRES_ATM - WETDEP( 24 ) * XC2            
!      PPAA0  = GAS( LPAA  ) * PRES_ATM - WETDEP( 25 ) * XC2            
      PCO20 = GAS (LCO2) * PRES_ATM + CO3A * XL 
!     &       - ( WETDEP( 11 ) + WETDEP( 12 ) + WETDEP( 13 ) ) * XC2    
                                                                        
!...don't allow gas concentrations to go below zero                     
                                                                        
      PSO20 = MAX (PSO20, 0.0) 
      PNH30 = MAX (PNH30, 0.0) 
      PH2O20 = MAX (PH2O20, 0.0) 
      PO30 = MAX (PO30, 0.0) 
      PFOA0 = MAX (PFOA0, 0.0) 
      PMHP0 = MAX (PMHP0, 0.0) 
      PPAA0 = MAX (PPAA0, 0.0) 
      PCO20 = MAX (PCO20, 0.0) 
      PHNO30 = MAX (PHNO30, 0.0) 
                                                                        
!...Molar concentrations of soluble aerosols                            
!...   = Initial amount - amount deposited  (mol/liter)                 
                                                                        
      TS6 = TS6ACCA - DS4 (0) 
      CL = CLA 
      CA = CAA 
      MG = MGA 
      NA = NAA 
      K = KA 
      FE = FEA 
      MN = MNA 
      A = 3.0 * FE 
      B = 2.0 * MN 
!      TS6     = TS6ACCA  + TS6AKNA                                     
!     &        - ( WETDEP(  6 ) + WETDEP(  7 ) ) * XC1                  
!     &        - DS4( 0 )                                               
!      CL      = CLA      -   WETDEP( 16 )  * XC1                       
!      CA      = CAA      -   WETDEP(  3 )  * XC1                       
!      MG      = MGA      -   WETDEP( 29 )  * XC1                       
!      NA      = NAA      -   WETDEP(  4 )  * XC1                       
!      K       = KA       -   WETDEP( 30 )  * XC1                       
!      FE      = FEA      -   WETDEP( 19 )  * XC1                       
!      MN      = MNA      -   WETDEP( 20 )  * XC1                       
!      ORGN    = ORGACCA + ORGAKNA - WETDEP( 27 )  * XC1                
!      PRIM    = PRIACCA + PRIAKNA - WETDEP( 28 )  * XC1                
!      PRIMCOR = PRICORA  -   WETDEP( 33 )  * XC1                       
!      A       = 3.0 * FE                                               
!      B       = 2.0 * MN                                               
!                                                                       
!...don't allow aerosol concentrations to go below zero                 
                                                                        
      TS6 = MAX (TS6, 0.0) 
      CL = MAX (CL, 0.0) 
      CA = MAX (CA, 0.0) 
      MG = MAX (MG, 0.0) 
      NA = MAX (NA, 0.0) 
!      K       = MAX( K,       0.0 )                                    
      FE = MAX (FE, 0.0) 
      MN = MAX (MN, 0.0) 
!      ORGN    = MAX( ORGN,    0.0 )                                    
!      PRIM    = MAX( PRIM,    0.0 )                                    
!      PRIMCOR = MAX( PRIMCOR, 0.0 )                                    
      A = MAX (A, 0.0) 
      B = MAX (B, 0.0) 
                                                                        
      SK6TS6 = SK6 * TS6 
                                                                        
!...find solution of the equation using a method of reiterative         
!...  bisections Make initial guesses for pH:   between .01  to  10.    
                                                                        
      HA = 0.01 
      HB = 10.0 
                                                                        
      I7777C = 0 
 7777 CONTINUE 
                                                                        
!! *** LIJIE                                                            
!      IF(II.EQ.65.AND.JJ.EQ.11.AND.KK.EQ.5)                            
!     &             PRINT*, I7777C,'AQUE STEP 1_3'                      
                                                                        
      I7777C = I7777C + 1 
      IF (I7777C.GE.1000) THEN 
!        XMSG = 'EXCESSIVE LOOPING AT I7777C'                           
!        CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT2 )              
         WRITE ( * , * ) 'Error in RADM aqueous chemistry (AQCHEM)' 
         WRITE ( * , * ) 'EXCESSIVE LOOPING AT I7777C' 
                   ! LIJIE                                              
         GOTO 8888 
      ENDIF 
                                                                        
      HA = MAX (HA - 0.8, 0.1) 
      HB = MIN (HB + 0.8, 9.9) 
      AE = 10.0** ( - HA) 
                                                                        
      RECIPA1 = 1.0 / (AE * ACT1) 
      RECIPA2 = 1.0 / (AE * AE * ACT2) 
                                                                        
!...calculate final gas phase partial pressure of SO2, NH3, HNO3        
!...  HCOOH, and CO2 (atm)                                              
                                                                        
      PSO2F = PSO20 / (1.0 + XLSO2 * (1.0 + SO21 * RECIPA1 + SO212 *    &
      RECIPA2) )                                                        
                                                                        
      PNH3F = PNH30 / (1.0 + XLNH3 * (1.0 + NH3DH20 * AE) ) 
                                                                        
      PFOAF = PFOA0 / (1.0 + XL * (FOAH + FOA1H * RECIPA1) ) 
                                                                        
      PHNO3F = PHNO30 / (1.0 + XLHNO3 * (1.0 + HNO31 * RECIPA1) ) 
                                                                        
      PCO2F = PCO20 / (1.0 + XLCO2 * (1.0 + CO21 * RECIPA1 + CO212 *    &
      RECIPA2) )                                                        
                                                                        
!...calculate liquid phase concentrations (moles/liter)                 
                                                                        
      SO4 = SK6TS6 / (AE * GM2 + SK6) 
      HSO4 = TS6 - SO4 
      SO3 = SO212H * PSO2F * RECIPA2 
      HSO3 = SO21H * PSO2F * RECIPA1 
      CO3 = CO212H * PCO2F * RECIPA2 
      HCO3 = CO21H * PCO2F * RECIPA1 
      OH = H2OW * RECIPA1 
      NH4 = NH31HDH * PNH3F * AE 
      HCO2 = FOA1H * PFOAF * RECIPA1 
      NO3 = HNO31H * PHNO3F * RECIPA1 
                                                                        
!...compute functional value                                            
                                                                        
      FA = AE+NH4 + 2.0 * (CA + MG - CO3 - SO3 - SO4) - OH - HCO3 -     &
      HSO3 - NO3 - HSO4 - HCO2                                          
                                                                        
!...Start iteration and bisection ****************<<<<<<<               
                                                                        
      I30C = 0 
   30 CONTINUE 
!! *** LIJIE ***                                                        
!      IF(II.EQ.65.AND.JJ.EQ.11.AND.KK.EQ.5.AND.I30C.EQ.10)             
!     &      PRINT*, I30C,'AQUE STEP 1_4'                               
                                                                        
                                                                        
      I30C = I30C + 1 
      IF (I30C.GE.1000) THEN 
!        XMSG = 'EXCESSIVE LOOPING AT I30C'                             
!        CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT2 )              
         WRITE ( * , * ) 'Error in RADM aqueous chemistry (AQCHEM)' 
         WRITE ( * , * ) 'EXCESSIVE LOOPING AT I30C' 
         GOTO 40 
      ENDIF 
                                                                        
      BB = (HA + HB) / 2.0 
      AE = 10.0** ( - BB) 
                                                                        
      ICNTAQ = ICNTAQ + 1 
      IF (ICNTAQ.GE.3000) THEN 
!        XMSG = 'Maximum AQCHEM total iterations exceeded'              
!        CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT2 )              
         WRITE ( * , * ) 'Error in RADM aqueous chemistry (AQCHEM)' 
         WRITE ( * , * ) 'EXCESSIVE LOOPING AT I3000C' 
         GOTO 40 
      ENDIF 
                                                                        
!      IF(II.EQ.65.AND.JJ.EQ.11.AND.KK.EQ.5.AND.I30C.EQ.10)             
!     &      PRINT*, I30C,'AQUE STEP 1_5'                               
                                                                        
                                                                        
      RECIPA1 = 1.0 / (AE * ACT1) 
      RECIPA2 = 1.0 / (AE * AE * ACT2) 
                                                                        
!...calculate final gas phase partial pressure of SO2, NH3, HNO3        
!...  HCOOH, and CO2 (atm)                                              
                                                                        
      PSO2F = PSO20 / (1.0 + XLSO2	 * (1.0 + SO21 * RECIPA1 + SO212 *   &
      RECIPA2) )                                                        
                                                                        
      PNH3F = PNH30 / (1.0 + XLNH3 * (1.0 + NH3DH20 * AE) ) 
                                                                        
      PHNO3F = PHNO30 / (1.0 + XLHNO3 * (1.0 + HNO31 * RECIPA1) ) 
                                                                        
      PFOAF = PFOA0 / (1.0 + XL * (FOAH + FOA1H * RECIPA1) ) 
                                                                        
      PCO2F = PCO20 / (1.0 + XLCO2 * (1.0 + CO21 * RECIPA1 + CO212 *    &
      RECIPA2) )                                                        
                                                                        
!...calculate liquid phase concentrations (moles/liter)                 
                                                                        
      SO4 = SK6TS6 / (AE * GM2 + SK6) 
      HSO4 = TS6 - SO4 
      SO3 = SO212H * PSO2F * RECIPA2 
      HSO3 = SO21H * PSO2F * RECIPA1 
      CO3 = CO212H * PCO2F * RECIPA2 
      HCO3 = CO21H * PCO2F * RECIPA1 
      OH = H2OW * RECIPA1 
      NH4 = NH31HDH * PNH3F * AE 
      HCO2 = FOA1H * PFOAF * RECIPA1 
      NO3 = HNO31H * PHNO3F * RECIPA1 
                                                                        
!...compute functional value                                            
                                                                        
      FB = AE+NH4 + 2.0 * (CA + MG - CO3 - SO3 - SO4) - OH - HCO3 -     &
      HSO3 - NO3 - HSO4 - HCO2                                          
                                                                        
!...Calculate and check the sign of the product of the two functional va
                                                                        
      FTST = FA * FB 
      IF (FTST.LE.0.0) THEN 
         HB = BB 
      ELSE 
         HA = BB 
         FA = FB 
      ENDIF 
                                                                        
!...Check convergence of solutions                                      
                                                                        
      HTST = HA / HB 
                                                                        
!      IF(II.EQ.65.AND.JJ.EQ.11.AND.KK.EQ.5.AND.I30C.EQ.10)             
!     &      PRINT*, I30C,'AQUE STEP 1_6'                               
                                                                        
      IF (HTST.LE.TST) GOTO 30 
                                                                        
!      IF(II.EQ.65.AND.JJ.EQ.11.AND.KK.EQ.5) PRINT*,'AQUE STEP 1_7'     
!...end of zero-finding routine ****************<<<<<<<<<<<<            
   40 CONTINUE 
!...compute Ionic strength and activity coefficient by the Davies equati
                                                                        
      STION = 0.5 * (AE+NH4 + OH + HCO3 + HSO3 + 4.0 * (SO4 + CO3 + SO3 &
      + CA + MG + MN) + NO3 + HSO4 + 9.0 * FE+NA + CL + A + B + HCO2)   
!     &      + NO3 + HSO4 + 9.0 * FE + NA + K + CL + A + B + HCO2)      
      GM1LOG = - 0.509 * (SQRT (STION) / (1.0 + SQRT (STION) ) - 0.2 *  &
      STION)                                                            
      GM2LOG = GM1LOG * 4.0 
      GM1 = 10.0**GM1LOG 
      GM2 = MAX (10.0**GM2LOG, 1.0E-30) 
      ACTB = ACT1 
      ACT1 = MAX (GM1 * GM1, 1.0E-30) 
      ACT2 = MAX (GM1 * GM1 * GM2, 1.0E-30) 
                                                                        
!...check for convergence and possibly go to 7777, to recompute         
!...  Gas and liquid phase concentrations                               
                                                                        
      TAC = ABS (ACTB - ACT1) / ACTB 
      IF (TAC.GE.1.0E-2) GOTO 7777 
               ! LIJIE                                                  
 8888 CONTINUE 
!...return an error if the pH is not in range                           
                                                                        
!cc      IF ( ( HA .LT. 0.02 ) .OR. ( HA .GT. 9.49 ) ) THEN             
      IF ( (HA.LT.0.1) .OR. (HA.GT.9.9) ) THEN 
         WRITE ( * , * ) 'Error in RADM aqueous chemistry (AQCHEM)' 
         WRITE ( * , * ) 'PH VALUE OUT OF RANGE: ', ha 
      ENDIF 
                                                                        
!...Make those concentration calculations which can be made outside     
!...  of the function.                                                  
                                                                        
      SO2L = SO2H * PSO2F 
      AC = 10.0** ( - BB) 
      SIV = SO3 + HSO3 + SO2L 
                                                                        
!      IF(II.EQ.65.AND.JJ.EQ.11.AND.KK.EQ.5) PRINT*,'AQUE STEP 1_8'     
!...Calculate final gas phase concentrations of oxidants (atm)          
                                                                        
      PH2O2F = (PH2O20 + XL * DS4 (1) ) / (1.0 + XLH2O2) 
      PO3F = (PO30 + XL * DS4 (2) ) / (1.0 + XLO3) 
      PMHPF = (PMHP0 + XL * DS4 (4) ) / (1.0 + XLMHP) 
      PPAAF = (PPAA0 + XL * DS4 (5) ) / (1.0 + XLPAA) 
                                                                        
      PH2O2F = MAX (PH2O2F, 0.0) 
      PO3F = MAX (PO3F, 0.0) 
      PMHPF = MAX (PMHPF, 0.0) 
      PPAAF = MAX (PPAAF, 0.0) 
                                                                        
!...Calculate liquid phase concentrations of oxidants (moles/liter)     
                                                                        
      H2O2L = PH2O2F * H2O2H 
      O3L = PO3F * O3H 
      MHPL = PMHPF * MHPH 
      PAAL = PPAAF * PAAH 
      FOAL = PFOAF * FOAH 
      NH3L = PNH3F * NH3H 
      CO2L = PCO2F * CO2H 
      HNO3L = PHNO3F * HNO3H 
                                                                        
!...load the liquid concentration array with current values             
                                                                        
      LIQUID (1) = AC 
      LIQUID (2) = NH4 
      LIQUID (3) = CA 
      LIQUID (4) = NA 
      LIQUID (5) = OH 
      LIQUID (6) = SO4 
      LIQUID (7) = HSO4 
      LIQUID (8) = SO3 
      LIQUID (9) = HSO3 
      LIQUID (10) = SO2L 
      LIQUID (11) = CO3 
      LIQUID (12) = HCO3 
      LIQUID (13) = CO2L 
      LIQUID (14) = NO3 
      LIQUID (15) = NH3L 
      LIQUID (16) = CL 
      LIQUID (17) = H2O2L 
      LIQUID (18) = O3L 
      LIQUID (19) = FE 
      LIQUID (20) = MN 
      LIQUID (21) = A 
      LIQUID (22) = FOAL 
      LIQUID (23) = HCO2 
      LIQUID (24) = MHPL 
      LIQUID (25) = PAAL 
      LIQUID (26) = 0.0 
!      LIQUID( 27 ) = ORGN                                              
!      LIQUID( 28 ) = PRIM                                              
      LIQUID (29) = MG 
      LIQUID (30) = K 
      LIQUID (31) = B 
      LIQUID (32) = HNO3L 
!      LIQUID( 33 ) = PRIMCOR                                           
                                                                        
!...if the maximum cloud lifetime has not been reached, the compute     
!...  the next timestep.                                                
!      IF(II.EQ.65.AND.JJ.EQ.11.AND.KK.EQ.5) PRINT*,                    
!     &                        TIMEW,TAUCLD,'AQUE STEP 1_9'             
                                                                        
      IF (TIMEW.LT.TAUCLD) THEN 
                                                                        
!...make kinetics calculations                                          
!...  note: DS4(i) and DSIV(I) are negative numbers!                    
                                                                        
                       !! by chenhs                                     
         DTRMV = 600.0 
!        IF ( ( CTHK1 .GT. 1.0E-10 ) .AND. ( PRCRATE .GT. 1.0E-10 ) )   
!     &     DTRMV = 3.6 * WTAVG * 1000.0 * CTHK1 / PRCRATE  ! <<<uma fou
!        DTRMV = MIN( DTRMV, 300.0 )                                    
         ITERAT = ITERAT + 1 
                                                                        
!...Define the total S(iv) available for oxidation                      
                                                                        
         TSIV = PSO20 * ONE_OVER_XL 
                                                                        
!...Calculate sulfur iv oxidation rate due to H2O2                      
                                                                        
         DSIVDT (1) = - RH2O2 * H2O2L * SO2L / (0.1 + AC) 
         IF ( (DSIVDT (1) .EQ.0.0) .OR. (TSIV.LE.1.0E-30) ) THEN 
            DTW (1) = DTRMV 
         ELSE 
            TOTOX = PH2O20 * ONE_OVER_XL 
            RATE = - DSIVDT (1) / (TOTOX * TSIV) 
            DTW (1) = 0.05 / (RATE * MAX (TOTOX, TSIV) ) 
         ENDIF 
                                                                        
!...Calculate sulfur iv oxidation rate due to O3                        
                                                                        
         IF (BB.GE.2.7) THEN 
            DSIVDT (2) = - 4.19E5 * (1.0 + 2.39E-4 / AC) * O3L * SIV 
         ELSE 
            DSIVDT (2) = - 1.9E4 * SIV * O3L / SQRT (AC) 
         ENDIF 
         IF ( (DSIVDT (2) .EQ.0.0) .OR. (TSIV.LE.1.0E-30) ) THEN 
            DTW (2) = DTRMV 
         ELSE 
            TOTOX = PO30 * ONE_OVER_XL 
            RATE = - DSIVDT (2) / (TOTOX * TSIV) 
            DTW (2) = 0.01 / (RATE * MAX (TOTOX, TSIV) ) 
         ENDIF 
                                                                        
                                                                        
!...Calculate sulfur iv oxidation rate due to 02 catalyzed by Mn++      
!...  and Fe+++  See Table IV Walcek & Taylor ( 1986)                   
                                                                        
                                  ! 4.0  < pH                           
         IF (BB.GE.4.0) THEN 
                                                                        
            IF (SIV.LE.1.0E-5) THEN 
               DSIVDT (3) = - 5000.0 * MN * HSO3 
            ELSEIF (SIV.GT.1.0E-5) THEN 
               DSIVDT (3) = - (4.7 * MN * MN / AC + 1.0E7 * FE * SIV *  &
               SIV)                                                     
                  ! end of first pass through SIV conc.                 
            ENDIF 
                                                                        
                      ! pH , + 4.0                                      
         ELSE 
                                                                        
            IF (SIV.LE.1.0E-5) THEN 
               DSIVDT (3) = - 3.0 * (5000.0 * MN * HSO3 + 0.82 * FE *   &
               SIV / AC)                                                
            ELSE 
               DSIVDT (3) = - (4.7 * MN * MN / AC + (0.82 * FE * SIV /  &
               AC) * (1.0 + 1.7E3 * MN**1.5 / (6.3E-6 + FE) ) )         
                 ! end of second pass through SIV conc.                 
            ENDIF 
                                                                        
                ! end of pass through pH                                
         ENDIF 
                                                                        
!         IF(II.EQ.65.AND.JJ.EQ.11.AND.KK.EQ.5) PRINT*,'AQUE STEP 1_10' 
                                                                        
         IF ( (DSIVDT (3) .EQ.0.0) .OR. (TSIV.LE.1.0E-30) ) THEN 
            DTW (3) = DTRMV 
         ELSE 
            RATE = - DSIVDT (3) / TSIV 
            DTW (3) = 0.1 / RATE 
         ENDIF 
                                                                        
!...Calculate sulfur oxidation rate due to MHP                          
                                                                        
         DSIVDT (4) = - RMHP * AC * MHPL * HSO3 
         IF ( (DSIVDT (4) .EQ.0.0) .OR. (TSIV.LE.1.0E-30) ) THEN 
            DTW (4) = DTRMV 
         ELSE 
            TOTOX = PMHP0 * ONE_OVER_XL 
            RATE = - DSIVDT (4) / (TOTOX * TSIV) 
            DTW (4) = 0.1 / (RATE * MAX (TOTOX, TSIV) ) 
         ENDIF 
                                                                        
!...Calculate sulfur oxidation due to PAA                               
                                                                        
         DSIVDT (5) = - RPAA * HSO3 * PAAL * (AC + 1.65E-5) 
         IF ( (DSIVDT (5) .EQ.0.0) .OR. (TSIV.LE.1.0E-30) ) THEN 
            DTW (5) = DTRMV 
         ELSE 
            TOTOX = PPAA0 * ONE_OVER_XL 
            RATE = - DSIVDT (5) / (TOTOX * TSIV) 
            DTW (5) = 0.1 / (RATE * MAX (TOTOX, TSIV) ) 
         ENDIF 
                                                                        
!...Calculate total sulfur iv oxidation rate                            
!         IF(II.EQ.65.AND.JJ.EQ.11.AND.KK.EQ.5)                         
!     &                     PRINT*,NUMOX,'AQUE STEP 1_11'               
                                                                        
         DSIVDT (0) = 0.0 
         DO IOX = 1, NUMOX 
         DSIVDT (0) = DSIVDT (0) + DSIVDT (IOX) 
         ENDDO 
                                                                        
!...Calculate a minimum time step required                              
                                                                        
         DTW (0) = MIN (DTW (1), DTW (2), DTW (3), DTW (4), DTW (5) ) 
                                                                        
                                                                        
         IF (DTW (0) .EQ.0.) DTW (0) = 1. 
!...check for large time step                                           
                                                                        
         IF (DTW (0) .GT.8.0E+37) THEN 
      WRITE ( * ,  * ) 'Warning in RADM aqueous chemistry (AQCHEM)' 
!          WRITE(idiag,1001)                                            
!     &          PRCRATE, DSIVDT(0), TS6, DTW(0), CTHK1, WTAVG          
            WRITE ( *, * ) DSIVDT (0), TS6, DTW (0) 
         ELSE 
                                                                        
!...calculate the change in sulfur iv for this time step                
                                                                        
!60        DTS6 = ABS( DTW( 0 ) * ( -DSIVDT( 0 ) - TS6 * PRCRATE        
!     &         / ( 3.6 * CTHK1 * WTAVG * 1000.0 ) ) )                  
   60       DTS6 = ABS (DTW (0) * DSIVDT (0) ) 
                                                                        
!...If DSIV(0), sulfur iv oxidized during this time step would be       
!... less than 5% of sulfur oxidized since time 0, then double DT       
                                                                        
            IF (DTW (0) .LE.TAUCLD) THEN 
               IF (DTS6.LT.0.05 * TS6) THEN 
                  DTW (0) = DTW (0) * 2.0 
                  GOTO 60 
               ENDIF 
            ENDIF 
         ENDIF 
         DTW (0) = MIN (DTW (0), DTRMV) 
                                                                        
!        IF(II.EQ.65.AND.JJ.EQ.11.AND.KK.EQ.5) PRINT*,'AQUE STEP 1_13'  
                                                                        
                                                                        
!   Set DTW( 0 ) to avoid overshooting - bkoo (09/16/2005)              
         IF (DSIVDT (0) .LT.0.0) THEN 
            DTW (0) = MIN (DTW (0), - TSIV * 1.00001 / DSIVDT (0) ) 
         ENDIF 
                                                                        
!...If the total time after this time increment will be greater than    
!...  TAUCLD sec., then set DTW(0) so that total time will be TAUCLD    
                                                                        
         IF (TIMEW + DTW (0) .GT.TAUCLD) DTW (0) = TAUCLD-TIMEW 
         IF (TS6.LT.1.0E-11) DTW (0) = TAUCLD-TIMEW 
!        IF ( ITERAT .GT. 100 ) DTW( 0 ) = TAUCLD - TIMEW               
         IF (ITERAT.GT.100) then 
            WRITE ( * , * ) 'AQCHEM: iterat > 100:', taucld 
            DTW (0) = TAUCLD-TIMEW 
         ENDIF 
!...Set DSIV(I), I = 0,NUMOX, the amount of S(IV) oxidized by each      
!... individual oxidizing agent, as well as the total.                  
!        IF(II.EQ.65.AND.JJ.EQ.11.AND.KK.EQ.5) PRINT*,'AQUE STEP 1_14'  
                                                                        
         DO IOX = 0, NUMOX 
         DS4 (IOX) = DS4 (IOX) + DTW (0) * DSIVDT (IOX) 
         ENDDO 
                                                                        
!...Compute depositions and concentrations for each species             
!                                                                       
!        DO LIQ = 1, NLIQS                                              
!          WETDEP( LIQ ) = WETDEP( LIQ )                                
!     &                  + PRCRATE * LIQUID( LIQ ) * DTW( 0 ) * WCAVG   
!     &                  * 1000.0 / ( 3600.0 * WTAVG * 1000.0 )         
!        END DO                                                         
                                                                        
         TIMEW = TIMEW + DTW (0) 
                                                                        
!...Return to make additional calculations                              
!       IF(II.EQ.65.AND.JJ.EQ.11.AND.KK.EQ.5) PRINT*,'AQUE STEP 1_15'   
         GOTO 20 
      ENDIF 
!C *** LIJIE ***                                                        
!C      IF(II.EQ.71.AND.JJ.EQ.72.AND.KK.EQ.1) PRINT*, I30C,'4444'       
!      IF(II.EQ.65.AND.JJ.EQ.11.AND.KK.EQ.5) PRINT*, 'AQUE STEP 2'      
                                                                        
!...At this point, TIMEW=TAUCLD                                         
!...  compute the scavenging coefficient for SO4 which will be used for 
!...  scavenging aerosol number in the accumulation and coarse mode     
!                                                                       
!      DEPSUM = ( WETDEP( 6 ) + WETDEP( 7 ) ) * XC1                     
!                                                                       
!      IF ( ( TS6ACCA + TS6AKNA - DS4( 0 ) ) .NE. 0.0 ) THEN            
!        BETASO4 = DEPSUM / ( ( TS6ACCA + TS6AKNA - DS4( 0 ) ) * TAUCLD 
!      ELSE                                                             
!        BETASO4 = 0.0                                                  
!      END IF                                                           
!                                                                       
!...Compute the output concentrations and wetdeposition amounts         
                                                                        
      TOTAMM = (PNH3F + (NH4 + NH3L) * XL) * RECIPAP1 
      TOTNIT = (PHNO3F + (NO3 + HNO3L) * XL) * RECIPAP1 
                                                                        
!                                                                       
!...gas-phase species wet deposition (mm mol/lit)                       
!                                                                       
!      GASWDEP( LSO2   ) = WETDEP(  8 ) + WETDEP(  9 ) + WETDEP( 10 )   
!      GASWDEP( LNH3   ) = WETDEP( 15 )                                 
!      GASWDEP( LH2O2  ) = WETDEP( 17 )                                 
!      GASWDEP( LO3    ) = WETDEP( 18 )                                 
!      GASWDEP( LCO2   ) = WETDEP( 11 ) + WETDEP( 12 ) + WETDEP( 13 )   
!      GASWDEP( LFOA   ) = WETDEP( 22 ) + WETDEP( 23 )                  
!      GASWDEP( LMHP   ) = WETDEP( 24 )                                 
!      GASWDEP( LPAA   ) = WETDEP( 25 )                                 
!      GASWDEP( LHNO3  ) = WETDEP( 32 )                                 
!      GASWDEP( LN2O5  ) = 0.0                                          
!      GASWDEP( LH2SO4 ) = 0.0                                          
!                                                                       
!...gas concentrations (mol/molV)                                       
                                                                        
      GAS (LSO2) = (PSO2F + XL * SIV) * RECIPAP1 
      GAS (LNH3) = FNH3 * TOTAMM 
      GAS (LH2O2) = (PH2O2F + XL * H2O2L) * RECIPAP1 
      GAS (LO3) = (PO3F + XL * O3L) * RECIPAP1 
      GAS (LCO2) = (PCO2F + XL * CO2L) * RECIPAP1 
      GAS (LFOA) = (PFOAF + XL * (FOAL + HCO2) ) * RECIPAP1 
      GAS (LMHP) = (PMHPF + XL * MHPL) * RECIPAP1 
      GAS (LPAA) = (PPAAF + XL * PAAL) * RECIPAP1 
      GAS (LHNO3) = FHNO3 * TOTNIT 
                          ! assume all into aerosol                     
      GAS (LN2O5) = 0.0 
                          ! assume all into aerosol                     
      GAS (LH2SO4) = 0.0 
                                                                        
!...aerosol species wet deposition (mm mol/lit)                         
!...  there is no wet deposition of aitken particles, they attached     
!...  to the accumulation mode particles                                
!                                                                       
!      AERWDEP( LSO4AKN ) = 0.0                                         
!      AERWDEP( LSO4ACC ) = WETDEP(  6 ) + WETDEP(  7 )                 
!      AERWDEP( LNH4AKN ) = 0.0                                         
!      AERWDEP( LNH4ACC ) = WETDEP(  2 )                                
!      AERWDEP( LNO3AKN ) = 0.0                                         
!      AERWDEP( LNO3ACC ) = WETDEP( 14 ) * FRACACC                      
!      AERWDEP( LNO3COR ) = WETDEP( 14 ) * FRACCOR                      
!      AERWDEP( LORGAKN ) = 0.0                                         
!      AERWDEP( LORGACC ) = WETDEP( 27 )                                
!      AERWDEP( LPRIAKN ) = 0.0                                         
!      AERWDEP( LPRIACC ) = WETDEP( 28 )                                
!      AERWDEP( LPRICOR ) = WETDEP( 33 )                                
!      AERWDEP( LNACL   ) = WETDEP(  4 )                                
!      AERWDEP( LA3FE   ) = WETDEP( 19 )                                
!      AERWDEP( LB2MN   ) = WETDEP( 20 )                                
!      AERWDEP( LCACO3  ) = WETDEP(  3 )                                
!      AERWDEP( LMGCO3  ) = WETDEP( 29 )                                
!      AERWDEP( LKCL    ) = WETDEP( 30 )                                
!      AERWDEP( LNUMAKN ) = 0.0                                         
!      AERWDEP( LNUMACC ) = AEROSOL( LNUMACC ) * AIRM                   
!     &                   * ( 1.0 - EXP( -BETASO4 * TAUCLD ) )          
!      AERWDEP( LNUMCOR ) = AEROSOL( LNUMCOR ) * AIRM                   
!     &                   * ( 1.0 - EXP( -BETASO4 * TAUCLD ) )          
!      AERWDEP( LSRFAKN ) = 0.0                                         
!      AERWDEP( LSRFACC ) = 0.0                                         
!                                                                       
!C...aerosol concentrations (mol/molV)                                  
!                                                                       
!amx                                                                    
      AEROSOL (LSO4) = TS6 * XL * RECIPAP1 
      AEROSOL (LNH4) = TOTAMM * (1.0 - FNH3) 
      AEROSOL (LNO3) = TOTNIT * (1.0 - FHNO3) 
      AEROSOL (LNACL) = NA * XL * RECIPAP1 
      AEROSOL (LA3FE) = FE * XL * RECIPAP1 
      AEROSOL (LB2MN) = MN * XL * RECIPAP1 
      AEROSOL (LCACO3) = CA * XL * RECIPAP1 
      AEROSOL (LMGCO3) = MG * XL * RECIPAP1 
      AEROSOL (LKCL) = K * XL * RECIPAP1 
!amx                                                                    
!      AEROSOL( LSO4AKN ) = AEROSOL( LSO4AKN ) * EXP( -ALFA3 * TAUCLD ) 
!      AEROSOL( LSO4ACC ) = TS6    * XL * RECIPAP1                      
!      AEROSOL( LNH4AKN ) = AEROSOL( LNH4AKN ) * EXP( -ALFA3 * TAUCLD ) 
!      AEROSOL( LNH4ACC ) = FNH4ACC * TOTAMM                            
!      AEROSOL( LNO3AKN ) = AEROSOL( LNO3AKN ) * EXP( -ALFA3 * TAUCLD ) 
!      AEROSOL( LNO3ACC ) = FNO3ACC * TOTNIT                            
!      AEROSOL( LNO3COR ) = FNO3COR * TOTNIT                            
!      AEROSOL( LORGAKN ) = AEROSOL( LORGAKN ) * EXP( -ALFA3 * TAUCLD ) 
!      AEROSOL( LORGACC ) = ORGN   * XL * RECIPAP1                      
!      AEROSOL( LPRIAKN ) = AEROSOL( LPRIAKN ) * EXP( -ALFA3 * TAUCLD ) 
!      AEROSOL( LPRIACC ) = PRIM   * XL * RECIPAP1                      
!      AEROSOL( LPRICOR ) = PRIMCOR* XL * RECIPAP1                      
!      AEROSOL( LNACL   ) = NA     * XL * RECIPAP1                      
!      AEROSOL( LA3FE   ) = FE     * XL * RECIPAP1                      
!      AEROSOL( LB2MN   ) = MN     * XL * RECIPAP1                      
!      AEROSOL( LCACO3  ) = CA     * XL * RECIPAP1                      
!      AEROSOL( LMGCO3  ) = MG     * XL * RECIPAP1                      
!      AEROSOL( LKCL    ) = K      * XL * RECIPAP1                      
!      AEROSOL( LNUMAKN ) = AEROSOL( LNUMAKN ) * EXP( -ALFA0 * TAUCLD ) 
!      AEROSOL( LNUMACC ) = AEROSOL( LNUMACC ) * EXP( -BETASO4 * TAUCLD 
!      AEROSOL( LNUMCOR ) = AEROSOL( LNUMCOR ) * EXP( -BETASO4 * TAUCLD 
!                                                                       
!C...compute the final accumulation aerosol 3rd moment                  
!                                                                       
!      M3NEW = ( AEROSOL( LSO4ACC ) * SGRAERMW( LSO4ACC ) / 1.8e6       
!     &      +   AEROSOL( LNH4ACC ) * SGRAERMW( LNH4ACC ) / 1.8e6       
!     &      +   AEROSOL( LNO3ACC ) * SGRAERMW( LNO3ACC ) / 1.8e6       
!     &      +   AEROSOL( LORGACC ) * SGRAERMW( LORGACC ) / 2.0e6       
!     &      +   AEROSOL( LPRIACC ) * SGRAERMW( LPRIACC ) / 2.2e6 )     
!CCC     &      * 6.0 / PI      ! cancels out in division below         
!                                                                       
!      AEROSOL( LSRFAKN ) = AEROSOL( LSRFAKN ) * EXP( -ALFA2 * TAUCLD ) 
!      AEROSOL( LSRFACC ) = AEROSOL( LSRFACC )                          
!     &                   * ( EXP( -BETASO4 * TAUCLD * ONETHIRD ) )     
!     &                   * ( M3NEW / MAX( M3OLD, 1.0E-30) ) ** TWOTHIRD
!                                                                       
!C...store the amount of hydrogen deposition                            
!                                                                       
!      HPWDEP = WETDEP( 1 )                                             
!                                                                       
!      IF(II.EQ.65.AND.JJ.EQ.11.AND.KK.EQ.5) PRINT*, 'AQUE STEP 3'      
      RETURN 
                                                                        
!...formats                                                             
                                                                        
!1001  FORMAT (1X,'STORM RATE=', F6.3, 'DSIVDT(0) =', F10.5,            
!     &       'TS6=', F10.5, 'DTW(0)=', F10.5, 'CTHK1=', F10.5,         
!     &       'WTAVG=', F10.5)                                          
 1001 FORMAT (1X, 'DSIVDT(0) =', F10.5,                                 &
             'TS6=', F10.5, 'DTW(0)=', F10.5)                           
                                                                        
      END SUBROUTINE RAQCHEM                        
      
      FUNCTION HLCONST ( NAME, TEMP ) 
                                                                        
!-----------------------------------------------------------------------
!                                                                       
!  FUNCTION: return the Henry's law constant for the specified substance
!            at the given temperature                                   
!                                                                       
!  revision history:                                                    
!    who        when           what                                     
!  ---------  --------  -------------------------------------           
!  S.Roselle  08/15/97  code written for Models-3                       
!  J.Gipson   06/18/01  added Henry's Law constants 50-55 for saprc99   
!  W.Hutzell  07/03/01  added Henry's Law constants 56-57 for Atrazine  
!                       and the daughter products from Atrazine and OH  
!                       reactions.                                      
!-----------------------------------------------------------------------
                                                                        
      IMPLICIT NONE 
                                                                        
!...........INCLUDES and their descriptions                             
                                                                        
!      INCLUDE SUBST_IODECL    ! I/O definitions and declarations       
                                                                        
!...........PARAMETERS and their descriptions:                          
                                                                        
                                        ! Number of substances          
      INTEGER       MXSPCS 
      PARAMETER   ( MXSPCS = 58 ) 
                                                                        
!...........ARGUMENTS and their descriptions                            
                                                                        
      REAL  HLCONST                                                       
                                        ! name of substance             
      CHARACTER*(*) NAME 
                                        ! temperature (K)               
      REAL          TEMP 
                                                                        
!...........SCRATCH LOCAL VARIABLES and their descriptions:             
                                                                        
                                        ! list of substance names       
      CHARACTER*16  SUBNAME( MXSPCS ) 
      SAVE          SUBNAME 
                                                                        
                                        ! species index                 
      INTEGER       SPC 
                                                                        
                                        ! Henry's law constants at 298.1
      REAL          A( MXSPCS ) 
      SAVE          A 
                                        ! enthalpy (like activation ener
      REAL          E( MXSPCS ) 
      SAVE          E 
                                                                        
!...........EXTERNAL FUNCTIONS and their descriptions:                  
                                                                        
!      INTEGER       HLINDEX 
                                                                        
!***********************************************************************
      DATA SUBNAME(  1), A(  1), E(  1)                                 &
          / 'O3              ', 1.2E-02, 2.7E+03 /                     
                                                     ! Chameides 1984   
      DATA SUBNAME(  2), A(  2), E(  2)                                 &
          / 'HO2             ', 9.0E+03, 0.0E+00 /                     
                                                     ! Chameides 1984   
      DATA SUBNAME(  3), A(  3), E(  3)                                 &
          / 'H2O2            ', 9.7E+04, 6.6E+03 /                     
                                                     ! Chameides 1984   
      DATA SUBNAME(  4), A(  4), E(  4)                                 &
          / 'NH3             ', 5.8E+01, 4.1E+03 /                     
                                                     ! Chameides 1984   
      DATA SUBNAME(  5), A(  5), E(  5)                                 &
          / 'NO              ', 1.9E-03, 1.4E+03 /                     
                                                     ! Lide and Frederik
      DATA SUBNAME(  6), A(  6), E(  6)                                 &
          / 'NO2             ', 1.2E-02, 2.5E+03 /                     
                                                     ! Chameides 1984   
      DATA SUBNAME(  7), A(  7), E(  7)                                 &
          / 'NO3             ', 1.2E+01, 1.9E+03 /                     
                                                     ! Chameides 1984   
      DATA SUBNAME(  8), A(  8), E(  8)                                 &
          / 'N2O5            ', 1.0E+30, 0.0E+00 /                     
                                                     ! "inf" Sander and 
      DATA SUBNAME(  9), A(  9), E(  9)                                 &
          / 'HNO2            ', 4.9E+01, 4.8E+03 /                     
                                                     ! Chameides 1984   
      DATA SUBNAME( 10), A( 10), E( 10)                                 &
          / 'HNO3            ', 2.6E+06, 8.7E+03 /                     
                                                     ! Chameides 1984   
      DATA SUBNAME( 11), A( 11), E( 11)                                 &
          / 'HNO4            ', 2.0E+04, 0.0E+00 /                     
                                                     ! Jacob et al. 1989
      DATA SUBNAME( 12), A( 12), E( 12)                                 &
          / 'SO2             ', 1.2E+00, 3.1E+03 /                     
                                                     ! Chameides 1984   
      DATA SUBNAME( 13), A( 13), E( 13)                                 &
          / 'H2SO4           ', 1.0E+30, 0.0E+00 /                     
                                                     ! infinity         
      DATA SUBNAME( 14), A( 14), E( 14)                                 &
          / 'METHANE         ', 1.3E-03, 0.0E+00 /                     
                                                     ! Mackay and Shiu 1
      DATA SUBNAME( 15), A( 15), E( 15)                                 &
          / 'ETHANE          ', 2.0E-03, 0.0E+00 /                     
                                                     ! Mackay and Shiu 1
      DATA SUBNAME( 16), A( 16), E( 16)                                 &
          / 'PROPANE         ', 1.4E-03, 0.0E+00 /                     
                                                     ! Mackay and Shiu 1
      DATA SUBNAME( 17), A( 17), E( 17)                                 &
          / 'BUTANE          ', 1.1E-03, 0.0E+00 /                     
                                                     ! Mackay and Shiu 1
      DATA SUBNAME( 18), A( 18), E( 18)                                 &
          / 'PENTANE         ', 8.1E-04, 0.0E+00 /                     
                                                     ! Mackay and Shiu 1
      DATA SUBNAME( 19), A( 19), E( 19)                                 &
          / 'HEXANE          ', 6.0E-04, 0.0E+00 /                     
                                                     ! Mackay and Shiu 1
      DATA SUBNAME( 20), A( 20), E( 20)                                 &
          / 'OCTANE          ', 3.4E-04, 0.0E+00 /                     
                                                     ! Mackay and Shiu 1
      DATA SUBNAME( 21), A( 21), E( 21)                                 &
          / 'NONANE          ', 2.0E-04, 0.0E+00 /                     
                                                     ! Mackay and Shiu 1
      DATA SUBNAME( 22), A( 22), E( 22)                                 &
          / 'DECANE          ', 1.4E-04, 0.0E+00 /                     
                                                     ! Mackay and Shiu 1
      DATA SUBNAME( 23), A( 23), E( 23)                                 &
          / 'ETHENE          ', 4.7E-03, 0.0E+00 /                     
                                                     ! Mackay and Shiu 1
      DATA SUBNAME( 24), A( 24), E( 24)                                 &
          / 'PROPENE         ', 4.8E-03, 0.0E+00 /                     
                                                     ! Mackay and Shiu 1
      DATA SUBNAME( 25), A( 25), E( 25)                                 &
          / 'ISOPRENE        ', 1.3E-02, 0.0E+00 /                     
                                                     ! Mackay and Shiu 1
      DATA SUBNAME( 26), A( 26), E( 26)                                 &
          / 'ACETYLENE       ', 4.1E-02, 1.8E+03 /                     
                                                     ! Wilhelm et al. 19
      DATA SUBNAME( 27), A( 27), E( 27)                                 &
          / 'BENZENE         ', 1.8E-01, 0.0E+00 /                     
                                                     ! Mackay and Shiu 1
      DATA SUBNAME( 28), A( 28), E( 28)                                 &
          / 'TOLUENE         ', 1.5E-01, 0.0E+00 /                     
                                                     ! Mackay and Shiu 1
      DATA SUBNAME( 29), A( 29), E( 29)                                 &
          / 'O-XYLENE        ', 2.0E-01, 0.0E+00 /                     
                                                     ! Mackay and Shiu 1
      DATA SUBNAME( 30), A( 30), E( 30)                                 &
          / 'METHANOL        ', 2.2E+02, 0.0E+00 /                     
                                                     ! Snider and Dawson
      DATA SUBNAME( 31), A( 31), E( 31)                                 &
          / 'ETHANOL         ', 1.6E+02, 0.0E+00 /                     
                                                     ! Betterton 1992   
      DATA SUBNAME( 32), A( 32), E( 32)                                 &
          / '2-CRESOL        ', 8.2E+02, 0.0E+00 /                     
                                                     ! Betterton 1992   
      DATA SUBNAME( 33), A( 33), E( 33)                                 &
          / '4-CRESOL        ', 1.3E+02, 0.0E+00 /                     
                                                     ! Betterton 1992   
      DATA SUBNAME( 34), A( 34), E( 34)                                 &
          / 'METHYLHYDROPEROX', 3.1E+02, 5.2E+03 /                     
                                                     ! O'Sullivan et al.
      DATA SUBNAME( 35), A( 35), E( 35)                                 &
          / 'FORMALDEHYDE    ', 7.0E+03, 6.4E+03 /                     
                                                     ! Chameides 1984   
      DATA SUBNAME( 36), A( 36), E( 36)                                 &
          / 'ACETALDEHYDE    ', 1.3E+01, 5.7E+03 /                     
                                                     ! Benkelberg et al.
      DATA SUBNAME( 37), A( 37), E( 37)                                 &
          / 'GENERIC_ALDEHYDE', 4.2E+03, 0.0E+00 /                     
                                                     ! Graedel and Goldb
      DATA SUBNAME( 38), A( 38), E( 38)                                 &
          / 'GLYOXAL         ', 3.6E+05, 0.0E+00 /                     
                                                     ! Zhou and Mopper 1
      DATA SUBNAME( 39), A( 39), E( 39)                                 &
          / 'ACETONE         ', 3.2E+01, 5.8E+03 /                     
                                                     ! Betterton 1991   
      DATA SUBNAME( 40), A( 40), E( 40)                                 &
          / 'FORMIC_ACID     ', 3.7E+03, 5.7E+03 /                     
                                                     ! Chameides 1984   
      DATA SUBNAME( 41), A( 41), E( 41)                                 &
          / 'ACETIC_ACID     ', 5.2E+03, 0.0E+00 /                     
                                                     ! Johnson et al. 19
      DATA SUBNAME( 42), A( 42), E( 42)                                 &
          / 'METHYL_GLYOXAL  ', 3.2E+04, 0.0E+00 /                     
                                                     ! Zhou and Mopper 1
      DATA SUBNAME( 43), A( 43), E( 43)                                 &
          / 'CO              ', 9.5E-04, 1.3E+03 /                     
                                                     ! Wilhelm et al. 19
      DATA SUBNAME( 44), A( 44), E( 44)                                 &
          / 'CO2             ', 3.1E-02, 2.4E+03 /                     
                                                     ! Chameides 1984   
      DATA SUBNAME( 45), A( 45), E( 45)                                 &
          / 'PAN             ', 2.9E+00, 5.9E+03 /                     
                                                     ! Pandis and Seinfe
      DATA SUBNAME( 46), A( 46), E( 46)                                 &
          / 'MPAN            ', 1.7E+00, 0.0E+00 /                     
                                                     ! Kames and Schurat
      DATA SUBNAME( 47), A( 47), E( 47)                                 &
          / 'OH              ', 3.0E+01, 4.5E+03 /                     
                                                     ! Hanson et al. 199
      DATA SUBNAME( 48), A( 48), E( 48)                                 &
          / 'METHYLPEROXY_RAD', 2.0E+03, 6.6E+03 /                     
                                                     ! Lelieveld and Cru
      DATA SUBNAME( 49), A( 49), E( 49)                                 &
          / 'PEROXYACETIC_ACI', 8.4E+02, 5.3E+03 /                     
                                                     ! O'Sullivan et al.
      DATA SUBNAME( 50), A( 50), E( 50)                                 &
          / 'PROPANOIC_ACID  ', 5.7E+03, 0.0E+00 /                     
                                                     ! Kahn et al. 1995 
      DATA SUBNAME( 51), A( 51), E( 51)                                 &
          / '2-NITROPHENOL   ', 7.0E+01, 4.6E+03 /                     
                                                     ! USEPA 1982       
      DATA SUBNAME( 52), A( 52), E( 52)                                 &
          / 'PHENOL          ', 1.9E+03, 7.3E+03 /                     
                                                     ! USEPA 1982       
      DATA SUBNAME( 53), A( 53), E( 53)                                 &
          / 'BIACETYL        ', 7.4E+01, 5.7E+03 /                     
                                                     ! Betteron 1991    
      DATA SUBNAME( 54), A( 54), E( 54)                                 &
          / 'BENZALDEHYDE    ', 4.2E+01, 4.6E+03 /                     
                                                     ! Zhou and Mopper 1
      DATA SUBNAME( 55), A( 55), E( 55)                                 &
          / 'PINENE          ', 4.9E-02, 0.0E+00 /                     
                                                     ! Karl and Lindinge
      DATA SUBNAME( 56), A( 56), E( 56)                                 &
          / 'ATRA            ', 4.1E+05, 6.0E+03 /                     
                                                     ! CIBA Corp (1989) 
      DATA SUBNAME( 57), A( 57), E( 57)                                 &
          / 'DATRA           ', 4.1E+05, 6.0E+03 /                     
                                                     ! assumed same as A
      DATA SUBNAME( 58), A( 58), E( 58)                                 &
          / 'ADIPIC_ACID     ', 2.0E+08, 0.0E+00 /                     
                                                     ! Saxena and Hildem
                                                                        
                                                                        
!-----------------------------------------------------------------------
!  begin body of subroutine HLCONST                                     
                                                                        
      SPC = HLINDEX( NAME, MXSPCS, SUBNAME ) 
                                                                        
      IF ( SPC .GT. 0 ) THEN 
        HLCONST = A( SPC ) * EXP( E( SPC )                              &
               * ( ( 298.0 - TEMP) / ( 298.0 * TEMP) ) )               
      ELSE 
        HLCONST = 0.0 
      END IF 
                                                                        
      RETURN 
      END  FUNCTION  HLCONST
      
      FUNCTION HLINDEX (NAME, N, NLIST) 
                                                                        
!***********************************************************************
!                                                                       
!  FUNCTION:                                                            
!                                                                       
!    Searches for NAME in list NLIST and returns the subscript          
!    (1...N) at which it is found, or returns 0 when NAME not           
!    found in NLIST                                                     
!                                                                       
!  PRECONDITIONS REQUIRED:  none                                        
!                                                                       
!  SUBROUTINES AND FUNCTIONS CALLED:  none                              
!                                                                       
!  REVISION HISTORY:                                                    
!                                                                       
!    5/88   Modified for ROMNET                                         
!    9/94   Modified for Models-3 by CJC                                
!                                                                       
!***********************************************************************
                                                                        
      IMPLICIT NONE 
                                                                        
!.......   Arguments and their descriptions:                            
      INTEGER  HLINDEX                                                         
                                !  Character string being searched for  
      CHARACTER*(*) NAME 
                                !  Length of array to be searched       
      INTEGER       N 
                                !  array to be searched                 
      CHARACTER*(*) NLIST(*) 
                                                                        
!.......   Local variable:                                              
                                                                        
                        !  loop counter                                 
      INTEGER       I 
                                                                        
!.....................................................................  
!.......   begin body of INDEX1()                                       
                                                                        
      DO 100 I = 1, N 
                                                                        
                                              ! Found NAME in NLIST     
          IF ( NAME .EQ. NLIST( I ) ) THEN 
              HLINDEX = I 
              RETURN 
          ENDIF 
                                                                        
  100 END DO 
                                                                        
                         !  not found                                   
      HLINDEX = 0 
      RETURN 
                                                                        
      END FUNCTION HLINDEX
      
end module aqueous_chem      
