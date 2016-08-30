module geatm_vartype
  parameter(ndomain=5) ! added by chenhs
  integer, dimension(ndomain) :: imother !,ichild  ! added by chenhs
  integer                     :: iichild
  integer, dimension(ndomain) :: sratio,tratio   !! by chenhs
  integer, dimension(ndomain) :: nx, ny ,nz, ntbeg,ntend,dtout,dtmet  !added by chenhs
  integer, dimension(ndomain) :: nxlo, nylo
  integer  ::  NEMIT,NEMIT_NOx,NEMIT_AIRC,NLAY_EM,NLAY_NOx,NLAY_AIRC,ifnox,ifairc,& 
               NEMIT_HgA,NLAY_HgA,ifHgA,ifHgN,NEMIT_HgN,NLAY_HgN  !! for Hg emission, by chenhs     
  !--------- make balance----------------------------
  REAL,ALLOCATABLE,DIMENSION(:,:) :: amount_bf,amount_af
  real :: total_af,total_bf,total_tmp(3),acdry,acwet
                                                      
  !-----------------------------------------------------
  ! added to read and transfer the contro lparameters by juanxiong he 
  !-----------------------------------------------------
  integer :: idifvert,ichemgas,idry, &
             iglobal,ifseacom,ifdustcom,imodis, & 
             ifprocess,ifglobal,ifHg,&
             iaer, isize, nzz, nest, ntt, &
             idmSet,iSrcDefined,ismMax,iHgtLMax,& 
             KDUSTTOP,ismMaxHere,igas,igasCBM,&
             iprecise,ikosaline,imasskeep,iprocess,&
             NSOA,ndustcom,nseacom,ifbalance,itotalspe
  integer :: ifsmt
  real :: hh

  !----------------------------------------------------
  integer, dimension(ndomain) :: comm2d, nbrleft, nbrright, nbrtop, nbrbottom, &
                                 sx, ex, sy, ey, stride,                       &
                                 irec, irecwet,irecdry,irecg,irecnox,irecairc, irec_as,&
                                 irec80,irecglobal,irecdust, irecsea,irec60,&
                                 irecMOZART,irecHgA,irecHgN
  integer, dimension(4,ndomain) :: bdysx, bdyex, bdysy, bdyey
                              !  1, south, 2, north, 3, west, 4, east
  integer, dimension(2,ndomain) :: dims, coords
  logical, dimension(2,ndomain) :: periods
  integer,dimension(:,:), allocatable :: procs ! for cam, added by Juanxiong He
  integer, dimension(ndomain) :: csx, cex, csy, cey, cnx, cny ! for cam, added by Juanxiong He
  integer myid, newid, numprocs, myid_x, myid_y, local_com ! for cam, added by Juanxiong He
  
  real,    dimension(ndomain) :: time,dtstep  !added by chenhs
  character*4 cdnum
  character*2 cdnum2
  integer iyear1,imonth1,idate1,ihour1
  integer   chem
  integer, allocatable, dimension(:,:) :: sxc, exc, syc, eyc
  integer, dimension(8,200) :: ISNDMRK  ! chenhs modified 200-->400
  real :: TTIME,TZ,ALONG,LAT
  
  !++++++++++++++++++ chenhs,polar transport ++++++++++++!
  INTEGER,DIMENSION(5,720) :: IPOLARMRK  !!! for 1x1 degree is 720
  INTEGER  IPOLARNUM
  !++++++++++++++++++ chenhs,polar transport ++++++++++++!
  !++++++++++++++++++ chenhs,from GEATM +++++++++++++++++!
  INTEGER MEMARK
  INTEGER,ALLOCATABLE,DIMENSION(:) :: STAMARK
  !++++++++++++++++++ chenhs,from GEATM +++++++++++++++++!

  ! for 2D variables 
  real, allocatable, dimension(:) :: TERRAIN, RAINCON,RAINNON,      &
                        LATITCRS, LONGICRS, LAND_USE, UST, U10, V10,&
                        T2, PSFC, RHSFC, PSEALVLC, SWDOWN,          &
                        PBL_HGT, REGIME, SHFLUX, LHFLUX, LWDOWN,    &
                        SWOUT, LWOUT, SOILT1, Q2,topo3,RMOL,clflo,&
                        clfmi,clfhi,HGT1,CLDOPD,landmask ! chenhs add landmask
  integer,allocatable,dimension(:) :: NPBL                         
                        !topo3 added by lijie 05-06-02
  integer, allocatable, dimension(:) :: ip2mem

  ! FOR DUST EMISSIONS 
  REAL,ALLOCATABLE,DIMENSION(:) :: Z0, UST0, FSOIL,FICE, FSNOW,FVEG
  REAL,ALLOCATABLE,DIMENSION(:) :: EMITFACT !DUSTEMISSION FACT 
  INTEGER,ALLOCATABLE,DIMENSION(:,:,:) :: IP3MEMAER
  REAL,ALLOCATABLE,DIMENSION(:) :: DUSTK,DUSTHGTF,TOTALDUST
  REAL,ALLOCATABLE,DIMENSION(:) :: DUSTEMISS,DUSTDRY,DUSTWET,DUSTGRAV ! DUST EMISS,DRY and WET and GRAVITY  DEPOSITION : kg/m2/hour

  ! FOR SEA SALT EMISSIONS
  REAL,ALLOCATABLE,DIMENSION(:) :: SEAEMISS

  ! FOR HETEROGENEOUS CHEMISTRY
  REAL,ALLOCATABLE,DIMENSION(:) :: DUSTCOMP ! DUST COMPOSITIONS 
  ! 1: CACO3, 2 MAGCO3; 3. K2CO3, 4:NA2CO3, 5: SIO2, 6: AL2O3,  7: SO4, 8: NO3, 9: Fe(II),10: Total Fe(III) 11: Coated Fe(III)
  REAL,ALLOCATABLE,DIMENSION(:) :: SEACOMP  ! SEA SALT COMPOSITIONS
                                               ! 1: CL, 2: NA 3 SS-SO4, 4: MG, 5: Ca, 6: K 7: SO4 8: NO3           
  REAL,ALLOCATABLE,DIMENSION(:) ::  DUSTDRYSO4,DUSTDRYNO3,DUSTDRYFeII, DUSTDRYFeIII ! kg/m2/hour
  REAL,ALLOCATABLE,DIMENSION(:) ::  DUSTWETSO4,DUSTWETNO3,DUSTWETFeII, DUSTWETFeIII ! kg/m2/hour
  REAL,ALLOCATABLE,DIMENSION(:) ::  DUSTGRAVSO4,DUSTGRAVNO3,DUSTGRAVFeII,DUSTGRAVFeIII ! kg/m2/hour
  ! FRACTION OF DIFFERENT BINS OF SSA AND DUST FOR HETEROGENEOUS PRODUCED SO4  AND NO3
  REAL,ALLOCATABLE,DIMENSION(:) :: FSO4_DUST,FNO3_DUST, FSO4_SSA,FNO3_SSA 
 
  ! for 2D variables such as deposition, emission, wet deposition
  real, allocatable, dimension(:) :: EmtAntGas,EmtShpGas,EmtBBGas,EmtBioGas, &
                                     EmtOceGas,EmtSoiGas,EmtLigGas,          &
                                     DryDGas,DryVelGas,        &
                                     WetDGas,                  &
                                     EmtBBHg,EmtOceHg,EmtGeoHg,EmtReeHg  !! for Hg natural emission, by chenhs
  integer, allocatable, dimension(:,:) :: ip2memGas


  ! for 3D variables 
  real, allocatable, dimension(:) :: u,v,t,clw,rnw,h,rh1,w,kh,kv, &
                                        dx,dy,dz,heiz,KOSASFC,Plev,QVAPOR,&
                                        jo1d,jno2,TAUCLDI,TAUCLDC
  integer, allocatable, dimension(:,:) :: ip3mem

  real,allocatable,dimension(:)  :: globalno2,globalo3,globalco

  ! additional 3D variable to keep mass balance
  ! Zifa 2004 09 02 at Acadimic Sinica
  real, allocatable, dimension(:) ::  kpmass_m1,RatioMass,kpmass_m2
                                      ! kpmass_m1 to mark the copy tracer
                                      ! kpmass_m2 to remember the real tracer

  ! for 4D variables 
  real, allocatable, dimension(:) :: gas, EmtAircGas, EmtHgAGas  !! added by chenhs
  !+++++++++++++++++++++++++for process by lijie+++++++++++++++++++++++
  real, allocatable,dimension(:)  ::  gasOLD
  !+++++++++++++++++++++++++for process by lijie+++++++++++++++++++++++
  integer, allocatable, dimension(:,:,:) :: ip4mem !,ib4mem,jb4mem

  ! for 5D variables 
  real, allocatable, dimension(:) :: aer!!, aer_src  by chenhs
  integer, allocatable, dimension(:,:,:,:) :: ip5mem !,ib5mem,jb5mem
  integer, allocatable, dimension(:,:,:,:) :: ip5memc,ip5memcs  ! DUST and Sea SALT COMPOSITIONS

  ! FOR DUST AND SEA SALT
  REAL,ALLOCATABLE,DIMENSION(:) :: MSIZDIS,MSIZDID ! sea salt and dust PARTICLES SIZE
  REAL :: SFT(4)   !! SFT IS THE MASS FRACTION OF DIFFERENT DUST PARTICLES SIZE IN TOTAL 
                !! DUST FOUR SIZE : 0.1-1.0um, 1-2.0um, 2.0-5.0um, 5.0-10.0um           
  real, allocatable, dimension(:) :: DryVeldust
  real, allocatable, dimension(:) :: SOILRH,SOILT
  real, allocatable, dimension(:) :: RK_HETSO2_DUST, RK_HETHNO3_DUST
 
  ! For Source Mark
  integer, dimension(ndomain) :: irecSM,ifsm
  real, allocatable, dimension(:,:,:) :: sm,smb,smp,sm1 !sm1 is for ACM2
  real,allocatable,dimension(:,:,:,:) :: smconv  !convective diffusion
  real,allocatable,dimension(:,:)  :: smconv1 !moist convection 
  real,allocatable,dimension(:) :: c00,smthis,smother
  real,allocatable, dimension(:) :: MapSource
  real,allocatable, dimension(:) :: tmpMarkCon
  real, allocatable, dimension(:) :: contem0  
  real, allocatable, dimension(:) :: TmpSM    
  integer, allocatable, dimension(:) :: igMark  ! to mark gas   (idmSet)
  integer, allocatable, dimension(:) :: iaMarkAer ! to mark aer type (idmSet)
  integer, allocatable, dimension(:) :: iaMarkSiz ! to mark aer size (idmSet)
  
  ! for 5D variables (for O3/SO2/SO4/et al.) <idm,ism,i,j,k>
  real, allocatable, dimension(:) :: SourceMark          !source mark of c 
  integer, allocatable, dimension(:,:,:,:) :: ipSMmem    !nz,ism,idm,ne
  real :: delta,delta1,delta2,delta3,delta4
 
  !++++++++++++++++++++++for  process by lijie+++++++++++++++++++++++
  real, allocatable, dimension(:) :: GasTermBal
  integer, allocatable, dimension(:) :: IGGPOS,IGOPos
  integer, allocatable, dimension(:,:,:,:) :: ipGasTermBal
  !++++++++++++++++++++++++for process by lijie+++++++++++++++++++++++

  !+++++++++++++++++==for ACM2++++++++++++++++++++
   real,allocatable,dimension(:) :: fconv,xmol ! xmol is the mon-ob length,kpbl is layer of pbl,
   integer,allocatable,dimension(:) :: kpbl    ! fconv is the fraction of  convection in pbl
   real :: nstep0(1000,1000)
   real, allocatable, dimension(:) :: gastmp
  !+++++++++++++++++++for ACM2++++++++++++++++++++

   integer, allocatable, dimension(:) :: IPIG,IPIG_NOx,IPIG_AIRC,IPIG_HgA,IPIG_HgN
   real, allocatable, dimension(:,:,:,:) :: ratioemitAnt,ratioemitShp,ratioemitBB,ratioemitBio,ratioemitOce,&
                                           ratioemitSoi,ratioemitLig,ratioemitAirc
   real, allocatable, dimension(:,:) :: contem
  
   real, allocatable, dimension(:,:)   :: atestR,atestS
   real, allocatable, dimension(:)   :: atestR0,atestS0
   !++++++++++++++++++ chenhs, polar transport ++++++++++++++++++++++++!
   REAL,ALLOCATABLE,DIMENSION(:,:)   :: ATESTRU,ATESTSU
   REAL,ALLOCATABLE,DIMENSION(:)   :: ATESTRU0,ATESTSU0
   !++++++++++++++++++ chenhs, polar transport ++++++++++++++++++++++++!
   real, allocatable, dimension(:,:) :: gravel
   integer, dimension(105) :: PrintTermGas  ! process by lijie

   integer, dimension(105) :: PrintGas,InitPrintGas ! chenhs, 103:HG0, 104:HG2, 105:HGP
   real, dimension(105) ::  dryvelinit,GC_MOLWT     ! chenhs, 103:HG0, 104:HG2, 105:HGP
   character*40, dimension(105) :: GC_Unit,GC_NAME  ! chenhs, 103:HG0, 104:HG2, 105:HGP
   REAL, dimension(1:74) ::  cppb,cnnzifa
   REAL, PARAMETER ::  STDATMPA = 101325.0 ! standard atmosphere  [ Pa ]
   REAL, PARAMETER ::  PA2ATM = 1.0 / STDATMPA
   LOGICAL LSUNLIGHT                ! Flag for daytime
   INTEGER ISTEPI
   INTEGER JDATE           ! Current date (YYYYDDD)
   INTEGER JTIME           ! Current time (HHMMSS)
   REAL   :: FCLD          ! THE COFFI OF CLOUD TO PHOTO
   REAL   :: CLFLO1,CLFMI1,CLFHI1,TER
   character*4 cyear
   character*2 cmonth,cday
   character*7 CHJDATE !character Julian date
!
!-------------------------------cloud and convection---------------
!the tempory variables
    real,allocatable,dimension(:,:,:) ::ppp,ttn,rkv,dzz,atm,ffn,conc,concmark,kp_tmp 
                                       !conc is the temporary gas concentrations
    real,allocatable,dimension(:,:,:,:) :: sea,dust ! tmporary sea and dust compositions
    real,allocatable,dimension(:) :: T1,Q1,QS1,U1,V1,P1,PH1 
    real,allocatable,dimension(:,:) :: TRA1,FTRA1,FTRA1D,FTRA1U,FTRA1O,FTRA1E
    real,allocatable,dimension(:)   :: CBMF1
!****** the cbmf(cloud base mass flux,must be remerbered in convect)
    
    integer, allocatable, dimension(:) :: ip2mem2dconv
    integer, allocatable, dimension(:,:) ::  ip3mem3dconv
    integer :: nxx,nyy ! the maximum grids in all domains

!----------------------------------dry------------
!the tempory variables
    real,allocatable,dimension(:,:,:) ::  height1,uu,vv,ww,ddx,ddy,ddz,water,QQ
    real,allocatable,dimension(:,:) :: tsurf,xlat,xlon,SWDOWN1,pblhgt,vdep
    integer,allocatable,dimension(:,:) :: land,icddrt,icdsno
    REAL,ALLOCATABLE,DIMENSION(:) :: DRYDEP2  ! THE DRYDEP AMOUNT UG/M2  !! by chenhs 
!***  ISORROPIA PARTS **************************************
!
    DOUBLE PRECISION,DIMENSION(5)  ::  WI  ! TOTAL ( GAS + AER0SOL)PRECURSORS OF ISOPPOPIA
!   WI(1) : TOTAL SODIUM   AS EQUIVALENT NA
!   WI(2) : TOTAL SULFATE  AS EQUIVALENT H2SO4
!   WI(3) : TOTAL AMMONIUM AS EQUIVALENT NH3
!   WI(4) : TOTAL NITRATE  AS EQUIVALENT HNO3
!   WI(5) : TOTAL CHLORIDE AS EQUIVALENT HCL
    DOUBLE PRECISION, DIMENSION(17)  :: WO
!   WO :  18 AEROSOLS SPECIES
!   WOG:  4  GASEOUS  SPECIES
    DOUBLE PRECISION ::  AWATER ! AEROSOL
    DOUBLE PRECISION,DIMENSION(4)   :: WOG
    REAL,ALLOCATABLE,DIMENSION(:) :: ANA,ASO4,ANH4,ANO3,ACL
    DOUBLE PRECISION              :: TEMPI,RHI 
!   TEMPI: TEMPERATURE IN K    
!   RHI  : RELATIVE IN 0-100%
!
!***  SOA MODULE *********************************************
!
    REAL,DIMENSION(6)  ::  SOA, SVG
    REAL               ::  POA , TEMP, PRESS0  
!    6: 6 SOA SPECIES AND SVG GASES  SV1,....SV6

!***   EXTINCTION AND AOD  ***********************************
!
    REAL :: NH4,NO3,SO4,BC,OC,PM10,PM25,NA,DUST01,DUST02,DUST03,DUST04,SEA01,SEA02,SEA03,SEA04
    REAL,ALLOCATABLE,DIMENSION(:) :: EXT  ! EXTINCTIONS IN KM-1
    REAL,ALLOCATABLE,DIMENSION(:) :: DUSTEXT  ! DUST EXTINCTIONS IN KM-1
    REAL,ALLOCATABLE,DIMENSION(:) :: VISIB ! VISIBITY in KM
    REAL,ALLOCATABLE,DIMENSION(:) :: EXT0,DUSTEXT0 ! TEMPORY 
    REAL,ALLOCATABLE,DIMENSION(:) :: AOD  ! AEROSOL OPITICAL DEPTH
    REAL,ALLOCATABLE,DIMENSION(:) :: DUSTAOD  ! DUST AEROSOL OPITICAL DEPTH 
    REAL,ALLOCATABLE,DIMENSION(:) :: PBLAOD  ! PBL AEROSOL OPITICAL DEPTH
    REAL,ALLOCATABLE,DIMENSION(:) :: SSA  ! Single scatter albedo
    REAL,ALLOCATABLE,DIMENSION(:) :: EXTS ! TEMPORY THE AOD DUE TO SCATTER
    REAL                          :: AODS ! AOD DUE TO SCATTER
    !
!***  GET THE COLUMN SO2 NO2 O3 IN DU   **********************
    REAL,ALLOCATABLE,DIMENSION(:) :: DSO2,DO3,DNO2 ! IN MOLES/M3
    REAL,ALLOCATABLE,DIMENSION(:) :: DUSO2,DUO3,DUNO2 ! IN DU

!***  FOR UV IRRADITION BY TUV  ********************************
!1 UV-B, 280-315 nm
!2 UV-B*, 280-320 nm
!3 UV-A, 315-400 nm
!4 vis+, > 400 nm 
    REAL, ALLOCATABLE,DIMENSION(:) :: UVA, UVB ,UVBS, VIS
    REAL                           :: CLWP,RAINNCV,RAINCV,CLD 
!----------------to define the layer of top conditions from global model
    real :: press
    real,allocatable,dimension(:) :: ktop !the layer
    real,allocatable,dimension(:,:) :: kktop,tp !tempory variable
    real,allocatable,dimension(:) :: tropp  ! the pressure of tropause
    real                          :: plimu,pliml,plimlex
    integer                       :: tperr
    LOGICAL                       :: dofill
!---------------------------------------                                    


!CCCC   AQUEOUS  CCCCCC
    REAL   :: CAER(9),CGAS (11)
    REAL,ALLOCATABLE,DIMENSION(:)     :: CPH  ! cloud water PH value
    REAL   :: CLW_TMP   !! tmp variable, by chenhs
!CCC

!CCC   WET DEPOSITION CCCCCCCCC
    LOGICAL :: LPREC  ! IF CONTAINING PRECITATION
    INTEGER    :: KTOPC, KBOTC ! THE BOTTOM AND TOP LAYER CONTAINING PRECITATION
    REAL,ALLOCATABLE,DIMENSION(:) :: CLWC,RNWC,TWET,PWET,RR,VOLRAT ! RR mm/hr drop volume/air volume
    REAL,ALLOCATABLE,DIMENSION(:) :: WETDEP, WETDEP2               ! 2-D array of wet deposited mass (mol/ha,g/ha)    
                                                                   ! and surface liquid concentrations (mol/l,g/l)
    REAL  :: TMASS
    REAL,ALLOCATABLE,DIMENSION(:,:,:) :: DEPFLD,DEPFLD2   ! TEMPORY 
    REAL,ALLOCATABLE,DIMENSION(:) :: CWC_C,PWR_C,CON ! TEMPORY

!CCC  net OPE  CCC
    REAL,ALLOCATABLE,DIMENSION(:) :: OPE           ! NET OZONE PRODUCTION EFFICIENCY

    logical :: lrddrt, lrdsno
    integer :: kkbb, nstep, nstep00, ndt      !! chenhs, timestep
    real :: dttmp,dttmp0                 !! by chenhs 
 
    integer, dimension(4000) :: IISCPU,IIRCPU,IsLocX,IsLocY,&
                              IsSNWE,IsNest,IrNest,IrLocX,IrLocY    ! chenhs modifid 4000-->8000,add IrNest
   !++++++++++++ Hg by chenhs +++++++++++++++++++++++++++!
    real :: CHG0, CHG2, CO3, CH2O2, COH, CSO2, CHO2, CPM10, tmid_sec00 
   !++++++++++++    end Hg    +++++++++++++++++++++++++++!
    DATA GC_MOLWT/98.,  63.,  36.,  17.,  30.,  46.,  &
              62.,  108., 47.,  79.,  48.,  16.,  &
              16.,  17.,  33.,  34.,  28.,  64.,  &
              16.,  30.,  47.,  28.,  30.,  32., &
              45.,  48.,  61.,  44.,  46.,  1.,   &
              72.,  121., 14.,  75.3,  72.,  28., &
              27.,  27.,  92.,  106., 108., 109., & 
              139., 84.,  14.,  1.,   1.,   1.,   &
              14.,  1.,   14.,   68.,  68.,  68., &
              68.,  68.,   1.,   1.,   1.,   1.,  &
              80.,  93.,  79.,  96.,  111., 125., &
              132., 160., 136., 136., 168., 168., &
              130., 130.,  1.,   1.,   12., 220., &
              1.0,  23.0, 18.0, 35.5, 96.0, 97.0, &
              62.0, 58.5, 142.0,85.0, 132.0,80.0, &
              53.5, 98.0, 115.0,120.0,247.0,136., &
              136., 168., 168., 130., 130., 18.,  &
              200.6,200.6,200.6/
              
    DATA GC_Unit/'ppb' , 'ppb','ppb','ppb','ppb','ppb' ,&
             'ppb' , 'ppb','ppb','ppb','ppb','ppb' ,&
             'ppb' , 'ppb','ppb','ppb','ppb','ppb' ,&
             'ppb' , 'ppb','ppb','ppb','ppb','ppb' ,&
             'ppb' , 'ppb','ppb','ppb','ppb','ppb' ,&
             'ppb' , 'ppb','ppb','ppb','ppb','ppb' ,&
             'ppb' , 'ppb','ppb','ppb','ppb','ppb' ,&
             'ppb' , 'ppb','ppb','ppb','ppb','ppb' ,&
             'ppb' , 'ppb','ppb','ppb','ppb','ppb' ,&
             'ppb' , 'ppb','ppb','ppb','ppb','ppb' ,&
             'ppb' , 'ppb','ppb','ppb','ppb','ppb' ,&
             'ppb' , 'ppb','ppb','ppb','ppb','ppb' ,&
             'ppb' , 'ppb',                         &
             'ugm-3','ugm-3','ugm-3','ugm-3','ugm-3',&
             'ugm-3','ugm-3','ugm-3','ugm-3','ugm-3',&
             'ugm-3','ugm-3','ugm-3','ugm-3','ugm-3',&
             'ugm-3','ugm-3','ugm-3','ugm-3','ugm-3',& 
             'ugm-3','ugm-3','ugm-3','ugm-3','ugm-3',&
             'ugm-3','ugm-3','ugm-3',                &
             'ngm-3','ngm-3','ngm-3'/
              
    DATA GC_NAME / 'H2SO4',   'HNO3'  , 'HCL'   , 'NH3' ,   'NO'   ,  'NO2'  ,&
               'NO3'  ,   'N2O5'  , 'HONO'  , 'HNO4',   'O3'   ,  'O1D'  ,&
               'O3P'  ,   'OH'    , 'HO2'   , 'H2O2',   'CO'   ,  'SO2'  ,&
               'CH4'  ,   'C2H6'  , 'CH3O2' , 'ETHP',   'HCHO' ,  'CH3OH',&
               'ANOL' ,   'CH3OOH', 'ETHOOH', 'ALD2',   'HCOOH',  'RCOOH',&
               'C2O3' ,   'PAN'   , 'PAR'   , 'AONE',   'MGLY' ,  'ETH'  ,&
               'OLET' ,   'OLEI'  , 'TOL'   , 'XYL' ,   'CRES' ,  'TO2'  ,&  
               'CRO'  ,   'OPEN'  , 'ONIT'  , 'ROOH',   'RO2'  ,  'ANO2' ,&
               'NAP'  ,   'XO2'   , 'XPAR'  , 'ISOP',   'ISOPRD', 'ISOPP',&
               'ISOPN',   'ISOPO2', 'DMS'   , 'MSA' ,   'DMSO' ,  'DMSO2',&
             'CH3SO2H', 'CH3SCH2OO','CH3SO2','CH3SO3','CH3SO2OO','CH3SO2CH2OO',&
             'SULFHOX',  'TERP'   , 'SV1'   , 'SV2' ,   'SV3'  ,  'SV4'  ,&
             'SV5'    ,  'SV6'    ,                                       &
             'PM25'   , 'PM10',   'BC'  ,   'OC'   , 'H+(AQ)',&
             'NA+AQ',  'NH4+AQ', 'CL-AQ','SO4--AQ','HSO4-AQ','NO3-AQ',&
             'NACLS',  'NA2SO4','NANO3','NH42SO4','NH4NO3','NH4CL',&
             'H2SO4AQ','NH4HSO4S','NAHSO4S','NH44HSO42', 'SOA1',&
             'SOA2'  ,   'SOA3'   ,   'SOA4'   ,  'SOA5'  ,  'SOA6', 'AH2O',&
             'HG0', 'HG2', 'HGP'/

! for gas (CBM-Z) and ISORROPIA AEROSOLS--------------------------------
!                 1        2      3      4      5      6       !
!                 H2SO4(g) HNO3   HCL    NH3    NO     NO2
!                 7        8      9      10     11     12      +
!                 NO3     N2O5   HONO   HNO4   O3     O1D
!                 13       14     15     16     17     18      +
!                 O3P     OH     HO2    H2O2   CO     SO2
!                 19       20     21     22     23     24      +
!                 CH4     C2H6    CH3O2  ETHP   HCHO  CH3OH
!                 25       26     27     28     29     30
!                 ANOL    CH3OOH  ETHOOH ALD2   HCOOH RCOOH
!                 31       32   , 33   , 34     35     36
!                 C2O3    PAN     PAR   AONE   MGLY   ETH
!                 37       38     39     40     41     42
!                 OLET  , OLEI ,  TOL   XYL    CRES   TO2    
!                 43       44     45     46     47     48
!                 CRO     OPEN    ONIT  ROOH   RO2    ANO2
!                 49       50     51    52      53     54
!                 NAP     XO2     XPAR  ISOP   ISOPRD ISOPP
!                 55       56     57    58      59     60 
!                 ISOPN   ISOPO2  DMS  MSA     DMSO   DMSO2
!                 61       62         63     64     65     66
!                 CH3SO2H CH3SCH2OO CH3SO2  CH3SO3 CH3SO2OO CH3SO2CH2OO
!                  67     68      69     70    71    72
!                 SULFH   TERP    SV1    SV2   SV3   SV4
!                  73     74 
!                  SV5    SV6
!                  75     76     77     78   79      80
!                   PM25  PM10   BC     OC   H+(AQ) NA+(AQ)
!                   81         82       83   84          85
!                 NH4+(AQ) CL-(AQ) SO4--(AQ) HSO4-(AQ)  NO3-(AQ)
!                  86      87        88         89        90
!                 NACL(S)  NA2SO4(S) NANO3(S)  NH42SO4(S) NH4NO3(S)
!                  91       92        93          94         95
!                  NH4CL(S) H2SO4(AQ) NH4HSO4(S)  NAHSO4(S) (NH4)4H(SO4)2(S)
!                  96      97        98      99    100   101, 102
!                  SOA1    SOA2     SOA3    SOA4   SOA5  SOA6 AH2O'
!                  103    104   105
!                  HG0    HG2   HGP
    DATA dryvelinit /0.02 ,  0.02 , 0.01, 0.65, 0.1,  0.5 ,&
                0.2,    0.3 ,  0.05, 0.05, 0.2,  0.001,&
                0.001 , 0.01 , 0.1,  0.1,  0.1,  0.48 ,&
                0.1 ,   0.1 ,  0.02, 0.02, 0.1,  0.1  ,&
                0.1,    0.1 ,  0.1,  0.1,  0.1,  0.1  ,&
                0.02,   0.1 ,  0.1,  0.1,  0.1,  0.1  ,&
                0.1,    0.1 ,  0.1,  0.1,  0.1,  0.1  ,&
                0.1,    0.1 ,  0.1,  0.1,  0.02, 0.02 ,&
                0.1,    0.1 ,  0.1,  0.1,  0.1,  0.1  ,&
                0.1,    0.1 ,  0.01, 0.01, 0.01, 0.01 ,&
                0.01,   0.01,  0.01, 0.01, 0.01, 0.01 ,&
                0.01,   0.02,  0.03, 0.01, 0.01, 0.01 ,&
                0.01,   0.01,  0.01, 0.01, 0.01, 0.01 ,&
                0.01,   0.01,  0.01, 0.01, 0.01, 0.01 ,&
                0.01,   0.01,  0.01, 0.01, 0.01, 0.01 ,&
                0.01,   0.01,  0.01, 0.01, 0.01, 0.01 ,&
                0.01,   0.01,  0.01, 0.01, 0.01, 0.01 ,&
                0.01,   0.01,  0.01/                    
    DATA PrintGas/0    , 0    , 0    , 1    , 1    , 1   , &
              0    , 0    , 1    , 0    , 1    , 0   , &
              0    , 1    , 1    , 1    , 1    , 1   , &
              0    , 0    , 0    , 0    , 0    , 0   , &
              0    , 0    , 0    , 0    , 0    , 0   , &
              0    , 0    , 0    , 0    , 0    , 0   , &
              0    , 0    , 0    , 0    , 0    , 0   , &
              0    , 0    , 0    , 0    , 0    , 0   , &
              0    , 0    , 0    , 1    , 0    , 0   , &
              0    , 0    , 0    , 0    , 0    , 0   , &
              0    , 0    , 0    , 0    , 0    , 0   , &
              0    , 1    , 0    , 0    , 0    , 0   , &
              0    , 0    , 1    , 1    , 1    , 1   , &
              0    , 0    , 0    , 0    , 0    , 0   , &
              0    , 0    , 0    , 0    , 0    , 0   , &
              0    , 0    , 0    , 0    , 0    , 1   , &
              1    , 1    , 1    , 1    , 1    , 0   , &
              1    , 1    , 1/
!!!!!!!!!!!!!!!!!!!!! for init conditions, chenhs !!!!!!!!!!!!!!!
    DATA InitPrintGas/1    , 1    , 1    , 1    , 1    , 1   , &
                  1    , 1    , 1    , 1    , 1    , 1   , &
                  1    , 1    , 1    , 1    , 1    , 1   , &
                  1    , 1    , 1    , 1    , 1    , 1   , &
                  1    , 1    , 1    , 1    , 1    , 1   , &
                  1    , 1    , 1    , 1    , 1    , 1   , &
                  1    , 1    , 1    , 1    , 1    , 1   , &
                  1    , 1    , 1    , 1    , 1    , 1   , &
                  1    , 1    , 1    , 1    , 1    , 1   , &
                  1    , 1    , 1    , 1    , 1    , 1   , &
                  1    , 1    , 1    , 1    , 1    , 1   , &
                  1    , 1    , 1    , 1    , 1    , 1   , &
                  1    , 1    , 1    , 1    , 1    , 1   , &
                  1    , 1    , 1    , 1    , 1    , 1   , &
                  1    , 1    , 1    , 1    , 1    , 1   , &
                  1    , 1    , 1    , 1    , 1    , 1   , &
                  1    , 1    , 1    , 1    , 1    , 1   , &
                  1    , 1    , 1/
!++++++++++++++++++++++for  process by lijie+++++++++++++++++++++++
    DATA PrintTermGas/0    , 0    , 0    , 0    , 0    , 0   , &
                  0    , 0    , 0    , 0    , 1    , 0   , &
                  0    , 0    , 0    , 0    , 0    , 0   , &
                  0    , 0    , 0    , 0    , 0    , 0   , &
                  0    , 0    , 0    , 0    , 0    , 0   , &
                  0    , 0    , 0    , 0    , 0    , 0   , &
                  0    , 0    , 0    , 0    , 0    , 0   , &
                  0    , 0    , 0    , 0    , 0    , 0   , &
                  0    , 0    , 0    , 0    , 0    , 0   , &
                  0    , 0    , 0    , 0    , 0    , 0   , &
                  0    , 0    , 0    , 0    , 0    , 0   , &
                  0    , 0    , 0    , 0    , 0    , 0   , &
                  0    , 0    , 0    , 0    , 0    , 0   , &
                  0    , 0    , 0    , 0    , 0    , 0   , &
                  0    , 0    , 0    , 0    , 0    , 0   , &
                  0    , 0    , 0    , 0    , 0    , 0   , &
                  0    , 0    , 0    , 0    , 0    , 0   , &
                  0    , 0    , 0/  
!++++++++++++++++++++++for process by lijie+++++++++++++++++++++++
    data periods/10*.false./

    DATA SFT/0.03,0.10,0.75,0.12/  !! LiJie

  contains
  
  subroutine initial_geatm_var
 
   IPOLARNUM=720 
   isrnum = itotalspe
   allocate(atestR(nzz,itotalspe),atestS(nzz,itotalspe))
   allocate(atestR0(itotalspe),atestS0(itotalspe))
 !++++++++++++++++++ chenhs,polar transport ++++++++++++!
   if(ifglobal.eq.-1) then
    ALLOCATE(ATESTRU(NZZ,2),ATESTSU(NZZ,2))
    ALLOCATE(ATESTRU0(2),ATESTSU0(2))
   endif
 !++++++++++++++++++ chenhs,polar transport ++++++++++++! 
!-------------------------------------------------------------------
  mem2d=0         ! number of memory for 2d variables
  do ne=1,nest
  mem2d=mem2d+(nx(ne)+2)*(ny(ne)+2)
  enddo
  allocate(TERRAIN(mem2d), RAINCON(mem2d), RAINNON(mem2d), &
           LATITCRS(mem2d), LONGICRS(mem2d), LAND_USE(mem2d), &
           UST(mem2d), U10(mem2d), V10(mem2d),T2(mem2d), &
           PSEALVLC(mem2d), SWDOWN(mem2d),PSFC(mem2d), &
           PBL_HGT(mem2d), REGIME(mem2d),RHSFC(mem2d), &
           SHFLUX(mem2d),LWDOWN(mem2d),LHFLUX(mem2d),  &
           SWOUT(mem2d),LWOUT(mem2d),SOILT1(mem2d),Q2(mem2d),&
           topo3(mem2d),RMOL(mem2d),NPBL(mem2d),CLDOPD(mem2d),&
           clflo(mem2d),clfmi(mem2d),clfhi(mem2d),HGT1(mem2d),landmask(mem2d))  !added by lijie 050602
  allocate(fconv(mem2d),xmol(mem2d),kpbl(mem2d))  ! ACM2 lijie
  allocate(ktop(mem2d),tropp(mem2d),CBMF1(mem2d))  ! the top layer of global model
  ALLOCATE(AOD(MEM2D),DUSTAOD(MEM2D),PBLAOD(MEM2D))        
  ALLOCATE(DUSO2(MEM2D),DUO3(MEM2D),DUNO2(MEM2D))
  ALLOCATE(FSOIL(MEM2D),FICE(MEM2D),FSNOW(MEM2D),FVEG(MEM2D),UST0(MEM2D),Z0(MEM2D)) ! DUST and SEA SALT
  ALLOCATE(EMITFACT(MEM2D),TOTALDUST(MEM2D),DUSTEMISS(MEM2D)) ! DUST AND SEA SALT
  ALLOCATE(SEAEMISS(MEM2D)) ! SEA SALT
  ALLOCATE(DUSTDRY(MEM2D), DUSTWET(MEM2D), DUSTGRAV(MEM2D) ) ! DUST AND SEA SALT
  if(ifdustcom.eq.1) then
   ALLOCATE(DUSTDRYSO4(MEM2D),DUSTDRYNO3(MEM2D),DUSTDRYFeII(MEM2D),DUSTDRYFeIII(MEM2D)) ! DUST AND SEA SALT
   ALLOCATE(DUSTWETSO4(MEM2D),DUSTWETNO3(MEM2D),DUSTWETFeII(MEM2D),DUSTWETFeIII(MEM2D)) ! DUST AND SEA SALT
   ALLOCATE(DUSTGRAVSO4(MEM2D),DUSTGRAVNO3(MEM2D),DUSTGRAVFeII(MEM2D),DUSTGRAVFeIII(MEM2D)) ! DUST AND SEA SALT
  endif

  !!!!!!!!!!!!!!!!!!!!!
  if(ifsmt>0)then ! For Source Mark
     allocate(MapSource(mem2d))
     allocate(tmpMarkCon(mem2d))
  endif        
  !!!!!!!!!!!!!!!!!!!!!

  allocate(ip2mem(nest))
  ii=1
  do ne=1,nest
  ip2mem(ne)=ii
  ii=ii+(nx(ne)+2)*(ny(ne)+2)
  enddo

  mem2dgas=0   !number of memory for 2d variables such emit,dep
  do ne=1,nest
  do ig=1,igas
  mem2dgas=mem2dgas+(nx(ne)+2)*(ny(ne)+2)
  enddo
  enddo
  allocate(EmtAntGas(mem2dgas),EmtShpGas(mem2dgas),EmtBBGas(mem2dgas),&
           EmtBioGas(mem2dgas),EmtOceGas(mem2dgas),EmtSoiGas(mem2dgas),EmtLigGas(mem2dgas))  
  allocate(EmtBBHg(mem2dgas),EmtOceHg(mem2dgas),EmtGeoHg(mem2dgas),EmtReeHg(mem2dgas)) !! for Hg natural emission, by chenhs         
  allocate(DryDGas(mem2dgas),DryVelGas(mem2dgas),DRYDEP2(mem2dgas))
  allocate(WetDGas(mem2dgas))
  ALLOCATE(WETDEP(MEM2DGAS),WETDEP2(MEM2DGAS))


  allocate(ip2memGas(igas,nest))

  ii=1
  do ne=1,nest
  do ig=1,igas
  ip2memGas(ig,ne)=ii
  ii=ii+(nx(ne)+2)*(ny(ne)+2)
  enddo
  enddo

!----- FOR DUST AND SEA SALT -----
 MEM3DAER =0
 DO NE = 1, NEST
 DO IA = 1, IAER
 DO IS = 1, ISIZE
  MEM3DAER = MEM3DAER + (NX(NE)+2)*(NY(NE)+2) 
 ENDDO
 ENDDO
 ENDDO 

 ALLOCATE(IP3MEMAER(ISIZE, IAER, NEST))
 II = 1
 DO NE = 1, NEST
 DO IA = 1, IAER
 DO IS  =1,ISIZE
 IP3MEMAER(IS,IA,NE) = II
 II = II + (NX(NE)+2)*(NY(NE)+2)
 ENDDO
 ENDDO
 ENDDO
 ALLOCATE(DryVeldust(MEM3DAER))

!------------------------for the cloud variable ------------
!for 2d variable 
      mem2dconv=0
     do ne=1,nest
        mem2dconv=mem2dconv+(nx(ne)+2)*(ny(ne)+2)
     enddo
     allocate(ip2mem2dconv(nest))
     ii=1
     do ne=1,nest
      ip2mem2dconv(ne)=ii
      ii=ii+(nx(ne)+2)*(ny(ne)+2)
     enddo
! for 3d varia
      mem3dconv=0
     do ne=1,nest
       do k=1,nz(ne)
        mem3dconv=mem3dconv+(nx(ne)+2)*(ny(ne)+2)
       enddo
     enddo
     allocate(ip3mem3dconv(nzz,nest))
     ii=1
     do ne=1,nest
       do k=1,nzz
        ip3mem3dconv(k,ne)=ii
        ii=ii+(nx(ne)+2)*(ny(ne)+2)
       enddo
     enddo
!-------------------------finished--------------------------                 
      
  mem3d=0        ! number of memory for 3d variables
  do ne=1,nest
  do k=1,nz(ne)
  mem3d=mem3d+(nx(ne)+2)*(ny(ne)+2)
  enddo
  enddo
  allocate(u(mem3d),v(mem3d),t(mem3d),h(mem3d),w(mem3d))
  allocate(clw(mem3d),rnw(mem3d),rh1(mem3d),Plev(mem3d),QVAPOR(mem3d))
  allocate(TAUCLDI(mem3d),TAUCLDC(mem3d))
  allocate(kh(mem3d),kv(mem3d))
  allocate(dx(mem3d),dy(mem3d),dz(mem3d),heiz(mem3d))
  allocate(KOSASFC(mem3d))
  allocate(globalno2(mem3d),globalo3(mem3d),globalco(mem3d))
  ALLOCATE(EXT(MEM3D),VISIB(MEM3D))
!  ALLOCATE(UVA(MEM3D),UVB(MEM3D),UVBS(MEM3D),VIS(MEM3D))  ! by chenhs
  ALLOCATE(DUSTEXT(MEM3D))
  allocate(jo1d(mem3d),jno2(mem3d),SSA(MEM3D))
  ALLOCATE(ANA(MEM3D),ASO4(MEM3D),ANO3(MEM3D),ANH4(MEM3D),ACL(MEM3D))
  ALLOCATE(CPH(MEM3D) )  ! AQUEOUS CHEMISTRY 
  ALLOCATE(OPE(MEM3D) )  ! OPE 
  ALLOCATE(SOILT(MEM3D),SOILRH(MEM3D))
  ALLOCATE(DUSTK(MEM3D),DUSTHGTF(MEM3D)) ! DUST EMISSIONS
  if(ifdustcom.eq.1) then
    ALLOCATE(RK_HETSO2_DUST(MEM3D),RK_HETHNO3_DUST(MEM3D)) 
  endif
  ! Zifa add 2004.09.02 Acedimic Sinica
  !
  if(imasskeep==1)then
  allocate(kpmass_m1(mem3d))
  allocate(kpmass_m2(mem3d))
  allocate(RatioMass(mem3d))
  endif

  allocate(ip3mem(nzz,nest))
  ii=1
  do ne=1,nest
  do k=1,nzz
  ip3mem(k,ne)=ii
  ii=ii+(nx(ne)+2)*(ny(ne)+2)
  enddo
  enddo

  mem4d=0        ! number of memory for 4d variables
  do ne=1,nest
  do ig=1,igas
  do k=1,nz(ne)
  mem4d=mem4d+(nx(ne)+2)*(ny(ne)+2)
  enddo
  enddo
  enddo
  allocate(gas(mem4d))
  allocate(EmtAircGas(mem4d),EmtHgAGas(mem4d))
 !----------------------acm2----
 ! allocate(gastmp(mem4d))   !by chenhs
 !--------------------------
!++++++++++++++++++++++for  process by lijie+++++++++++++++++++++++
 if(ifprocess.eq.1) then
   allocate(gasOLD(mem4d))
 endif
!++++++++++++++++++++++for  process by lijie+++++++++++++++++++++++
  
  allocate(ip4mem(nzz,igas,nest))
  ii=1
  do ne=1,nest
  do ig=1,igas
  do k=1,nzz
  ip4mem(k,ig,ne)=ii
  ii=ii+(nx(ne)+2)*(ny(ne)+2)
  enddo
  enddo
  enddo
!++++++++++++++++++++++for  process by lijie+++++++++++++++++++++++
if(ifprocess.eq.1) then
   iPrintTermGas = 0
     do ig=1,igas
        iPrintTermGas = iPrintTermGas + PrintTermGas(ig)
     enddo
  allocate(IGGPOS( iPrintTermGas ))
  allocate(IGOPos( igas ))
   iPrintTermGas = 0
   IGOPos=0
 do ig=1,igas
   if( PrintTermGas(ig) == 1 )then
     iPrintTermGas = iPrintTermGas + 1
     IGGPOS(iPrintTermGas) = ig
     IGOPos(ig) = iPrintTermGas
   endif
 enddo
                                               
   mem5d=0        ! number of memory for Track -GAS -Balance
     do ne=1,nest
       do ig=1,iPrintTermGas
         do ip=1,iprocess ! this term is very large
           do k=1,nzz
             mem5d=mem5d+(nx(ne)+2)*(ny(ne)+2)
            enddo
         enddo
       enddo
      enddo
  allocate(GasTermBal(mem5d))
  allocate(ipGasTermBal(nzz,iprocess,iPrintTermGas,nest))

  ii=1
    do ne=1,nest
      do ig=1,iPrintTermGas
        do ip=1,iprocess
          do k=1,nzz
            ipGasTermBal(k,ip,ig,ne)=ii
              ii=ii+(nx(ne)+2)*(ny(ne)+2)
          enddo
         enddo
      enddo
     enddo
endif
!++++++++++++++++++++++for  process by lijie+++++++++++++++++++++++

! FOR DUST AND SEA SALT COMPOSITIONS
if(ifdustcom.eq.1) then
  allocate(ip5memc(nzz,isize,ndustcom,nest))

  mem5dc=0
  do ne = 1, nest
  do iduc = 1, ndustcom
  do is = 1, isize
  do k = 1, nzz
   mem5dc = mem5dc + (nx(ne)+2)*(ny(ne)+2)
  enddo
  enddo
  enddo
  enddo
   allocate (DUSTCOMP(mem5dc))
   
  ii = 1
  do ne = 1,nest
  do iduc = 1, ndustcom
  do is = 1, isize
  do k = 1, nzz
    ip5memc (k, is, iduc, ne) = ii
    ii = ii + (nx(ne) + 2) * ( ny(ne) + 2 )
  enddo
  enddo
  enddo
  enddo
endif
if(ifseacom.eq.1) then
  allocate(ip5memcs(nzz,isize,nseacom,nest))
  mem5dc = 0
  do ne = 1, nest
  do iduc = 1, nseacom
  do is = 1, isize
  do k = 1, nzz
      mem5dc = mem5dc + (nx(ne)+2)*(ny(ne)+2)
  enddo
  enddo
  enddo
  enddo

    allocate ( SEACOMP(mem5dc) )
  ii = 1
  do ne = 1, nest
  do iduc = 1, nseacom
  do is = 1, isize
  do k = 1, nzz
   ip5memcs(k,is,iduc, ne) = ii
   ii = ii + ( nx(ne) + 2) * ( ny(ne) + 2 )
  enddo
  enddo
  enddo
  enddo
endif
! END DUST AND SEA SALT COMPOSITIONS

! FOR DUST AND SEA SALT TOTAL CONCENTRATIONS--------------------------------------------------------
  mem5d=0        ! number of memory for 5D variables
  do ne=1,nest
  do ia=1,iaer
  do is=1,isize
  do k=1,nzz
  mem5d=mem5d+(nx(ne)+2)*(ny(ne)+2)
  enddo
  enddo
  enddo
  enddo

  allocate(aer(mem5d))  !!! by chenhs
  allocate(ip5mem(nzz,isize,iaer,nest))

  ii=1
  do ne=1,nest
  do ia=1,iaer
  do is=1,isize
  do k=1,nzz
  ip5mem(k,is,ia,ne)=ii
  ii=ii+(nx(ne)+2)*(ny(ne)+2)
  enddo
  enddo
  enddo
  enddo

 !!!!!!!!!!!!!!!!!!!!!
 ! For Source Mark
 if(ifsmt>0)then  ! For Source Mark
  mem5d=0        ! number of memory for 5d variables
  do ne=1,nest
  if(ifsm(ne)==1)then
     do idm=1,idmSet
     do ism=1,ismMax
     do k=1,nz(ne)
        mem5d=mem5d+(nx(ne)+2)*(ny(ne)+2)
        enddo
        enddo
        enddo
    endif
  enddo
  allocate(SourceMark(mem5d))

  allocate(ipSMmem(nzz,ismMax,idmSet,nest))
  ii=1
  do ne=1,nest
    if(ifsm(ne)==1)then
      do idm=1,idmSet
      do ism=1,ismMax
      do k=1,nzz
         ipSMmem(k,ism,idm,ne)=ii
         ii=ii+(nx(ne)+2)*(ny(ne)+2)
         enddo
         enddo
         enddo
     endif
  enddo
 endif

 !+++++++++++++++ chenhs,from GEATM ++++++++++++++!
   if(ifglobal.eq.1) then
     MEMARK=4*NUMPROCS
     ALLOCATE(STAMARK(MEMARK))
   endif

  allocate(sxc(0:numprocs-1,nest))
  allocate(exc(0:numprocs-1,nest))
  allocate(syc(0:numprocs-1,nest))
  allocate(eyc(0:numprocs-1,nest))
  
  if(ifbalance.eq.1) then
   ALLOCATE(amount_bf(NUMPROCS,igas),amount_af(NUMPROCS,igas))
  endif
  
  end subroutine initial_geatm_var

  subroutine final_geatm_var
  
  deallocate(ip2mem)
  deallocate(TERRAIN, RAINCON, RAINNON, T2, PSEALVLC,SWDOWN, &
             LATITCRS, LONGICRS, LAND_USE, UST, U10, V10,RMOL,&
             PBL_HGT,NPBL,clflo,clfmi,clfhi,HGT1,CLDOPD,landmask)
  deallocate(ip2memGas)
  deallocate(EmtAntGas,EmtShpGas,EmtBBGas,EmtBioGas,EmtOceGas,EmtSoiGas,EmtLigGas)
  deallocate(EmtBBHg,EmtOceHg,EmtGeoHg,EmtReeHg)  !! for Hg natural emission, by chenhs 
  deallocate(DryDGas,DryVelGas,WetDGas,DRYDEP2)
  DEALLOCATE(WETDEP,WETDEP2)

  deallocate(ip3mem)
  deallocate(u,v,t,h,w,dx,dy,dz,rnw,clw,rh1,Plev,TAUCLDI,TAUCLDC)
  !--------------lijie modify--------------------
  deallocate(globalno2,globalo3,globalco)
  !----------------finish------------------------
  ! Zifa 2004/09/02
  if(imasskeep==1)then
  deallocate(kpmass_m1, RatioMass ,kpmass_m2)
  endif

  deallocate(ip4mem)
  deallocate(gas)
  deallocate(EmtAircGas,EmtHgAGas)
  !deallocate(gastmp) !acm2 ! by chenhs
!++++++++++++++++++++++++for process by lijie+++++++++++++++++++++++
  if(ifprocess.eq.1) then
  deallocate(gasOLD)
  deallocate(GasTermBal,ipGasTermBal)
  deallocate(IGGPOS,IGOPos)
  endif
!++++++++++++++++++++++++for process by lijie+++++++++++++++++++++++
!!++++++++++++++ chenhs,from GEATM +++++++++++++++++!!
  if(ifglobal.eq.1) then
   DEALLOCATE(STAMARK)
  endif
!!++++++++++++++ chenhs,from GEATM +++++++++++++++++!!
! *** DUST and SEA SALT
  DEALLOCATE(GRAVEL,DUSTK,DUSTHGTF,TOTALDUST)
  DEALLOCATE(FICE,FSNOW,SOILT,SOILRH,FVEG,EMITFACT)
  DEALLOCATE(DUSTEMISS,DUSTDRY,DUSTWET,DUSTGRAV)
  if(ifdustcom.eq.1) then
   DEALLOCATE(DUSTDRYSO4,DUSTDRYNO3,DUSTDRYFEII,DUSTDRYFEIII)
   DEALLOCATE(DUSTWETSO4,DUSTWETNO3,DUSTWETFEII,DUSTWETFEIII)
   DEALLOCATE(DUSTGRAVSO4,DUSTGRAVNO3,DUSTGRAVFEII,DUSTGRAVFEIII)
   DEALLOCATE(RK_HETSO2_DUST,RK_HETHNO3_DUST)
  endif
  DEALLOCATE(SEAEMISS)
! ***  AQUEOUS CHEMISTRY ******
  deallocate(CPH)
  DEALLOCATE (OPE) ! OPE
!=========================cloud and convection=========================
!  deallocate(ppp,ttn,ffn,conc,ip2mem2dconv,ip3mem3dconv)                              
!========================finished======================================
  deallocate(ktop,tropp) 
  deallocate(ASO4,ANO3,ACL,ANA,ANH4)
  deallocate(ip5mem)
  if(ifdustcom.eq.1) deallocate (ip5memc) ! DUST and SEA SALT
  if(ifseacom.eq.1)  deallocate (ip5memcs)  ! DUST and SEA SALT
  !deallocate(aer,aer_src)
  deallocate(aer)  !! by chenhs
  deallocate(jo1d,jno2)
  deallocate(DUSTEXT,EXT,VISIB,SSA,AOD,PBLAOD,DUSTAOD)
  deallocate(DUSO2,DUNO2,DUO3)
  !deallocate(UVB,UVBS,UVA,VIS)
  deallocate(syc,sxc,exc,eyc)
  deallocate(atestR,atestS)
  deallocate(atestR0,atestS0)
  !!!!!!! make balance !!!!!!!!!!!!
  if(ifbalance.eq.1) then
  DEALLOCATE(amount_af,amount_bf)
  endif
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!++++++++++++++++++++++ chenhs,polar transport +++++++++++++++++++++++++!!
  if(ifglobal.eq.-1) then
  DEALLOCATE(ATESTRU,ATESTSU)
  DEALLOCATE(ATESTRU0,ATESTSU0)
  endif
!!++++++++++++++++++++++ chenhs,polar transport +++++++++++++++++++++++++!!
  !!!!!!!!!!!!!!!!!!!!!!
  ! for Source Mark
  if(ifsmt>0)then
  deallocate(igMark,iaMarkAer,iaMarkSiz)
  deallocate(MapSource,tmpMarkCon)
  deallocate(SourceMark)
  endif

  end subroutine final_geatm_var
end module geatm_vartype  
