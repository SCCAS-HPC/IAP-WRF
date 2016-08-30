!=======================================================================
!
! *** ISORROPIA CODE
! *** INCLUDE FILE 'ISRPIA.INC'
! *** THIS FILE CONTAINS THE DECLARATIONS OF THE GLOBAL CONSTANTS
!     AND VARIABLES. 
!
! *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
! *** GEORGIA INSTITUTE OF TECHNOLOGY
! *** WRITTEN BY ATHANASIOS NENES
! *** UPDATED BY CHRISTOS FOUNTOUKIS
!
!=======================================================================
!
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NCOMP=5,NIONS=7,NGASAQ=3,NSLDS=9,NPAIR=13,NZSR=100,&
                NERRMX=25)
!
! *** INPUT VARIABLES **************************************************
!
      INTEGER METSTBL
      COMMON /INPT/ W(NCOMP), WAER(NCOMP), TEMP, RH, IPROB, METSTBL, &
                   NADJ
!
! *** WATER ACTIVITIES OF PURE SALT SOLUTIONS **************************
!
      COMMON /ZSR / AWAS(NZSR), AWSS(NZSR), AWAC(NZSR), AWSC(NZSR), &
                   AWAN(NZSR), AWSN(NZSR), AWSB(NZSR), AWAB(NZSR), &
                   AWSA(NZSR), AWLC(NZSR)
!
! *** DELIQUESCENCE RELATIVE HUMIDITIES ********************************
!
      INTEGER WFTYP
      COMMON /DRH / DRH2SO4,  DRNH42S4, DRNAHSO4, DRNACL,   DRNANO3,  &
                    DRNA2SO4, DRNH4HS4, DRLC,     DRNH4NO3, DRNH4CL
      COMMON /MDRH/ DRMLCAB,  DRMLCAS,  DRMASAN,  DRMG1,    DRMG2, &
                    DRMG3,    DRMH1,    DRMH2,    DRMI1,    DRMI2, &
                    DRMI3,    DRMQ1,    DRMR1,    DRMR2,    DRMR3, &
                    DRMR4,    DRMR5,    DRMR6,    DRMR7,    DRMR8, &
                    DRMR9,    DRMR10,   DRMR11,   DRMR12,   DRMR13,& 
                    WFTYP
!
! *** VARIABLES FOR LIQUID AEROSOL PHASE *******************************
!
      DOUBLE PRECISION MOLAL, MOLALR, M0
      REAL IONIC
      LOGICAL CALAOU, CALAIN, FRST, DRYF
      COMMON /IONS/ MOLAL(NIONS), MOLALR(NPAIR), GAMA(NPAIR), ZZ(NPAIR), &
                    Z(NIONS),     GAMOU(NPAIR),  GAMIN(NPAIR),M0(NPAIR), &
                    GASAQ(NGASAQ),                                       &  
                    EPSACT,       COH,           CHNO3,       CHCL,      &    
                    WATER,        IONIC,         IACALC,                 &
                    FRST,         CALAIN,        CALAOU,      DRYF
!
! *** VARIABLES FOR SOLID AEROSOL PHASE ********************************
!
      COMMON /SALT/ CH2SO4,  CNH42S4, CNH4HS4, CNACL,   CNA2SO4,& 
                    CNANO3,  CNH4NO3, CNH4CL,  CNAHSO4, CLC
!
! *** VARIABLES FOR GAS PHASE ******************************************
!
      COMMON /GAS / GNH3, GHNO3, GHCL 
!
! *** EQUILIBRIUM CONSTANTS ********************************************
!
      COMMON /EQUK/ XK1, XK2, XK3, XK4, XK5, XK6, XK7, XK8, XK9, XK10,&
                    XK11,XK12,XK13,XK14,XKW, XK21,XK22,XK31,XK32,XK41,&
                    XK42
!
! *** MOLECULAR WEIGHTS ************************************************
!
      DOUBLE PRECISION IMW
      COMMON /OTHR/ R, IMW(NIONS), WMW(NCOMP), SMW(NPAIR)
!
! *** SOLUTION/INFO VARIABLES ******************************************
!
      CHARACTER SCASE*15
      COMMON /CASE/ SULRATW, SULRAT, SODRAT, SCASE
!
      COMMON /SOLN/ EPS, MAXIT, NSWEEP, NDIV, ICLACT
!
! *** ERROR SYSTEM *****************************************************
!
      CHARACTER ERRMSG*40
      INTEGER   ERRSTK, NOFER   
      LOGICAL   STKOFL   
      COMMON /EROR/ STKOFL, NOFER, ERRSTK(NERRMX), ERRMSG(NERRMX)
!
! *** GENERIC VARIABLES ************************************************
!
      CHARACTER VERSION*15
      COMMON /CGEN/ GREAT, TINY, TINY2, ZERO, ONE, VERSION
!
! *** END OF INCLUDE FILE **********************************************
!
