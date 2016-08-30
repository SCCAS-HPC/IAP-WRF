#include <define.h>

 subroutine CLMDRIVER (dolai,doalb,dosst) 

!=======================================================================
!
! CLM MODEL DRIVER
!
! Original author : Yongjiu Dai, 09/30/1999; 08/30/2002
!
!=======================================================================

 use precision
 use phycon_module, only : tfrz, rgas, vonkar
 use paramodel
 use spmd
 use colm_varMod, only : numpatch, numcolumn, numgrid, ftune, forc, &
                         fcon_col, fvar_col, fldv_col, &
                         fcon_pft, fvar_pft, fldv_pft, &
                         wxy_patch, wxy_column, itypwat, oro, &
                         fLitterSoil, fLitterAtmos
 use timemgr, only : idate, dtime
 use spmd_decomp, only: cgmap, gxmap, gymap
#ifdef VEGDATA
 use vegdata, only : interp_vegdata
#endif
 use debug, only : c_bug, mybug
 implicit none

! ------------------- arguments  ----------------------------------

  logical, INTENT(in) :: dolai    ! true if time for time-varying vegetation paramter
  logical, INTENT(in) :: doalb    ! true if time for surface albedo calculation
  logical, INTENT(in) :: dosst    ! true if time for update sst/ice/snow

  real(r8),pointer  :: ivt(:)     ! land cover type number  add by zhq 06/02/2009

! ----------------------------------------------------------------
! I. Time invariant model variables
! ----------------------------------------------------------------

  integer  :: itypwat_c                 ! land water type
  real(r8), pointer :: dlat             ! latitude in radians
  real(r8), pointer :: dlon             ! longitude in radians
  real(r8), pointer :: wt_col
  real(r8), pointer :: wt_pft(:)

                       ! Soil physical parameters
  real(r8), pointer :: albsol           ! soil albedo for different coloured soils [-]
  real(r8), pointer :: csol  (:)        ! heat capacity of soil solids [J/(m3 K)]
  real(r8), pointer :: porsl (:)        ! fraction of soil that is voids [-]
  real(r8), pointer :: phi0  (:)        ! minimum soil suction [mm]
  real(r8), pointer :: bsw   (:)        ! clapp and hornbereger "b" parameter [-]
  real(r8), pointer :: dkmg  (:)        ! thermal conductivity of soil minerals [W/m-K]
  real(r8), pointer :: dksatu(:)        ! thermal conductivity of saturated soil [W/m-K]
  real(r8), pointer :: dkdry (:)        ! thermal conductivity for dry soil  [W/(m-K)]
  real(r8), pointer :: hksati(:)        ! hydraulic conductivity at saturation [mm h2o/s]
   
                       ! Vegetation static parameters
  real(r8), pointer :: z0m(:)              ! aerodynamic roughness length [m]
  real(r8), pointer :: displa(:)           ! displacement height [m]
  real(r8), pointer :: sqrtdi(:)           ! inverse sqrt of leaf dimension [m**-0.5]
  real(r8), pointer :: effcon(:)           ! quantum efficiency of RuBP regeneration (molCO2/molquanta)
  real(r8), pointer :: vmax25(:)           ! maximum carboxylation rate at 25 C at canopy top
  real(r8), pointer :: slti(:)             ! s3: slope of low temperature inhibition function     
  real(r8), pointer :: hlti(:)             ! s4: 1/2 point of low temperature inhibition function
  real(r8), pointer :: shti(:)             ! s1: slope of high temperature inhibition function  
  real(r8), pointer :: hhti(:)             ! s2: 1/2 point of high temperature inhibition function 
  real(r8), pointer :: trda(:)             ! s5: temperature coefficient in gs-a model            
  real(r8), pointer :: trdm(:)             ! s6: temperature coefficient in gs-a model           
  real(r8), pointer :: trop(:)             ! temperature coefficient in gs-a model          
  real(r8), pointer :: gradm(:)            ! conductance-photosynthesis slope parameter
  real(r8), pointer :: binter(:)           ! conductance-photosynthesis intercep
  real(r8), pointer :: extkn(:)            ! coefficient of leaf nitrogen allocation
  real(r8), pointer :: chil(:)             ! leaf angle distribution factor
  real(r8), pointer :: ref   (:,:)        ! leaf reflectance (iw=iband, il=life and dead)
  real(r8), pointer :: tran  (:,:)        ! leaf transmittance (iw=iband, il=life and dead)
  real(r8), pointer :: rootfr(:,:)        ! fraction of roots in each soil layer
#if(defined DGVM)
  real(r8), pointer :: pftpar(:,:)           !32 parameters of PFTs
  real(r8), pointer :: vegclass(:)            !1.tree 2.shrub 3.grass 4.crop -1.others
  real(r8), pointer :: summergreen(:)         !1. for summergreen; otherwise -1.
  real(r8), pointer :: raingreen(:)           !1. for raingreen; otherwise -1.
  real(r8), pointer :: sla(:)                 !sla
  real(r8), pointer :: lm_sapl(:)             !leafmass
  real(r8), pointer :: sm_sapl(:)             !sapwood mass
  real(r8), pointer :: hm_sapl(:)             !heartwood mass
  real(r8), pointer :: rm_sapl(:)             !rootmass
#endif
#if(defined DyN)
  real(r8), pointer :: cton_soil(:)          ! soil C:N mass ratio
  real(r8), pointer :: cton_pro (:)          ! C:N mass ratio in production
#endif

              ! CLM time step and TUNABLE constants
  real(r8) :: zlnd             ! roughness length for soil [m]
  real(r8) :: zsno             ! roughness length for snow [m]
  real(r8) :: csoilc           ! drag coefficient for soil under canopy [-]
  real(r8) :: dewmx            ! maximum dew
  real(r8) :: wtfact           ! fraction of model area with high water table
  real(r8) :: capr             ! tuning factor to turn first layer T into surface T
  real(r8) :: cnfac            ! Crank Nicholson factor between 0 and 1 
  real(r8) :: ssi              ! irreducible water saturation of snow
  real(r8) :: wimp             ! water impremeable if porosity less than wimp
  real(r8) :: pondmx           ! ponding depth (mm)
  real(r8) :: smpmax           ! wilting point potential in mm
  real(r8) :: smpmin           ! restriction for min of soil poten. (mm)
  real(r8) :: trsmx0           ! max transpiration for moist soil+100% veg.  [mm/s]
  real(r8) :: tcrit            ! critical temp. to determine rain or snow

! -----------------------------------------------------------------
! II. Time-varying state variables which reaquired by restart run
! -----------------------------------------------------------------
                       ! Currrnt cal2ar 
  integer :: year             ! current year of model run
  integer :: jday             ! current julian day of model run 
  integer :: msec             ! current seconds of model run (0 - 86400)

                       ! Main land surface variables 
  real(r8), pointer :: z   (:)          ! node depth [m]
  real(r8), pointer :: dz  (:)          ! interface depth [m]
  real(r8), pointer :: tss (:)          ! soil temperature [K]
  real(r8), pointer :: wliq(:)          ! liquid water in layers [kg/m2]
  real(r8), pointer :: wice(:)          ! ice lens in layers [kg/m2]
  real(r8), pointer :: tg               ! ground surface temperature [K]
  real(r8), pointer :: tlsun(:)         ! sunlit leaf temperature [K]
  real(r8), pointer :: tlsha(:)         ! shaded leaf temperature [K]
  real(r8), pointer :: ldew(:)          ! depth of water on foliage [mm]
  real(r8), pointer :: sag              ! non dimensional snow age [-]
  real(r8), pointer :: scv              ! snow cover, water equivalent [mm]
  real(r8), pointer :: snowdp           ! snow depth [meter]

                       ! Vegetation dynamic parameters 
  real(r8), pointer :: fveg(:)             ! fraction of vegetation cover
  real(r8), pointer :: fsno                ! fraction of snow cover on ground
  real(r8), pointer :: sigf(:)             ! fraction of veg cover, excluding snow-covered veg [-]
  real(r8), pointer :: green(:)            ! leaf greenness
  real(r8), pointer :: lai(:)              ! leaf area index
  real(r8), pointer :: sai(:)              ! stem area index

#if (defined DGVM)
  real(r8), pointer :: t10min(:)              !annual minimum of 10-day running mean (K)
  real(r8), pointer :: lai_ind(:)             !LAI per individual
  real(r8), pointer :: dphen(:)               !phenology [0 to 1]
  real(r8), pointer :: leafon(:)              !leafon days
  real(r8), pointer :: leafof(:)              !leafoff days
  real(r8), pointer :: firelength(:)          !fire season in days
  real(r8), pointer :: litterag(:)            !above ground litter
  real(r8), pointer :: litterbg(:)            !below ground litter
  real(r8), pointer :: cpool_fast(:)          !fast carbon pool
  real(r8), pointer :: cpool_slow(:)          !slow carbon pool
  real(r8), pointer :: k_fast_ave(:)          !decomposition rate
  real(r8), pointer :: k_slow_ave(:)          !decomposition rate
  real(r8), pointer :: litter_decom_ave(:)    !decomposition rate
  real(r8), pointer :: fmicr(:)               !microbial respiration (mol CO2 /m**2 /s)
  real(r8), pointer :: nind(:)                !number of individuals (#/m**2)
  real(r8), pointer :: lm_ind(:)              !individual leaf mass
  real(r8), pointer :: sm_ind(:)              !individual sapwood mass
  real(r8), pointer :: hm_ind(:)              !individual heartwood mass
  real(r8), pointer :: rm_ind(:)              !individual root mass
  real(r8), pointer :: tmomin20(:)            !20-yr running mean of tmomin
  real(r8), pointer :: agdd0(:)               !growing dgree days above 0
  real(r8), pointer :: agdd(:)                !growing dgree days above 5
  real(r8), pointer :: agddtw(:)              !growing dgree days above twmax
  real(r8), pointer :: agdd20(:)              !20-yr running mean of agdd
  real(r8), pointer :: t_mo(:)                !30-day mean temperature of 2m (K)
  real(r8), pointer :: t_mo_sum(:)            !30-day accumulated temperature of 2m (K)
  real(r8), pointer :: t_mo_min(:)            !annual min of t_mo (Kelvin)
  real(r8), pointer :: crownarea(:)           !area that each individual tree takes up (m^2)
  real(r8), pointer :: htop(:)                !canopy top
  real(r8), pointer :: tsai(:)                !one-sided stem area index, no burying by snow
  real(r8), pointer :: fpcgrid(:)             !foliar projective cover on gridcell (fraction)
  real(r8), pointer :: bm_inc(:)              !biomass increment
  real(r8), pointer :: afmicr(:)              !annual microbial respiration
  real(r8), pointer :: annpsn(:)              !annual photosynthesis (umol CO2 /m**2)
  real(r8), pointer :: annpsnpot(:)           !annual potential photosynthesis (same units)
  real(r8), pointer :: tref10(:)              !10-day averaged temperature at 2m
  real(r8), pointer :: tref_sum(:)            !sum of temperature in current day
  real(r8), pointer :: t10(:,:)               !array to record the 10 day temperature
  real(r8), pointer :: assimn10(:)            !10-day averaged assimilation rate
  real(r8), pointer :: assimn_sum(:)          !sum of assimn of current day
  real(r8), pointer :: an10(:,:)              !array to record 10 day assimn
  real(r8), pointer :: anngpp(:)              !annual gpp
  real(r8), pointer :: annfrmf(:)             !annual frmf
  real(r8), pointer :: annfrms(:)             !annual frms
  real(r8), pointer :: annfrmr(:)             !annual frmr
  real(r8), pointer :: annfrg(:)              !annual frg
  real(r8), pointer :: cflux_litter_soil(:)   !litter->soil
  real(r8), pointer :: cflux_litter_atmos(:)  !litter->atmos
  real(r8), pointer :: nday                   !counting the model days
  real(r8), pointer :: nyr                    !counting the model years
  real(r8), pointer :: turnover_ind(:)        !individual turnover biomass
  real(r8), pointer :: fpc_inc(:)             !fpc increase
  real(r8), pointer :: prec365                !yearly running mean of precipitation(mm/s)
  real(r8), pointer :: ifpre(:)               !-1=no PFT present;1=PFT present in this grid
#if(defined DyN)
  real(r8), pointer :: litter_leaf(:)      ! leaf-derived litter for PFT on modelled area basis (gC/m2)
  real(r8), pointer :: litter_wood(:)      ! heart&sapwood-derived litter for PFT on modelled area basis(gC/m2)
  real(r8), pointer :: litter_root(:)      ! fine root-derived litter for PFT on modelled area basis(gC/m2)
  real(r8), pointer :: litter_repr(:)      ! litter derived from allocation to reproduction for PFT on modelled

  real(r8), pointer :: litter_leaf_n(:) ! leaf-derived N litter for PFT on modelled area basis (gN/m2)
  real(r8), pointer :: litter_wood_n(:) ! heart&sapwood-derived N litter for PFT on modelled area basis(gN/m2)
  real(r8), pointer :: litter_root_n(:) ! fine root-derived N litter for PFT on modelled area basis (gN/m2)
  real(r8), pointer :: litter_repr_n(:) ! litter derived from allocation to reproduction N for PFT on modelled
                                        ! area basis (gN/m2)
  real(r8), pointer :: afcton_leaf(:)   ! annual floating leaf C:N ratio
  real(r8), pointer :: afcton_root(:)   ! annual floating root C:N ratio
  real(r8), pointer :: afcton_sap(:)    ! annual floating sapwood C:N ratio
  real(r8), pointer :: lm_ind_n(:)      ! individual leaf nitrogen mass
  real(r8), pointer :: sm_ind_n(:)      ! individual sapwood nitrogen mass
  real(r8), pointer :: hm_ind_n(:)      ! individual heartwood nitrogen mass
  real(r8), pointer :: rm_ind_n(:)      ! individual root nitrogen mass
                                        ! gN/m2 veget'd area for each pft
  real(r8), pointer :: an_up(:)         ! annual plant nitrogen uptake(gN/m2 vegt'd area)
  real(r8), pointer :: an_stress(:)     ! annual plant nitrogen stress(-)

  real(r8), pointer :: soil_no3        
  real(r8), pointer :: soil_no2       
  real(r8), pointer :: soil_no       
  real(r8), pointer :: soil_n2o     
  real(r8), pointer :: soil_n2     
  real(r8), pointer :: soil_nh4    
#endif
#endif

#if(defined VEGDATA)
  real(r8) :: lai_r  (numpft+1)           ! read-in leaf area index
  real(r8) :: sai_r  (numpft+1)           ! read-in stem+dead leaf area index
  real(r8) :: green_r(numpft+1)           ! read-in leaf greenness
  real(r8) :: fveg_r (numpft+1)           ! read-in vegetation fraction cover
#endif

                       ! Radiation  related (albedoes)
  real(r8), pointer :: coszen             ! cosine of solar zenith angle
  real(r8), pointer :: albg (:,:)         ! albedo, ground [-]
  real(r8), pointer :: albv (:,:)         ! albedo, vegetation [-]
  real(r8), pointer :: alb  (:,:)         ! averaged albedo [-]
  real(r8), pointer :: ssun (:,:)         ! sunlit canopy absorption for solar radiation (0-1)
  real(r8), pointer :: ssha (:,:)         ! shaded canopy absorption for solar radiation (0-1)
  real(r8), pointer :: thermk(:)           ! canopy gap fraction for tir radiation
  real(r8), pointer :: extkb(:)            ! (k, g(mu)/mu) direct solar extinction coefficient
  real(r8), pointer :: extkd(:)            ! diffuse and scattered diffuse PAR extinction coefficient

                       ! Additional variables required by reginal model (WRF & RSM) 
  real(r8), pointer :: trad(:)             ! radiative temperature of surface [K]
  real(r8), pointer :: tref(:)             ! 2 m height air temperature [kelvin]
  real(r8), pointer :: qref(:)             ! 2 m height air specific humidity
  real(r8), pointer :: rst(:)              ! canopy stomatal resistance (s/m)

  real(r8), pointer :: emis(:)             ! averaged bulk surface emissivity
  real(r8), pointer :: z0ma(:)             ! effective roughness [m]
  real(r8), pointer :: zol(:)              ! dimensionless height (z/L) used in Monin-Obukhov theory
  real(r8), pointer :: rib(:)              ! bulk Richardson number in surface layer
  real(r8), pointer :: ustar(:)            ! u* in similarity theory [m/s]
  real(r8), pointer :: qstar(:)            ! q* in similarity theory [kg/kg]
  real(r8), pointer :: tstar(:)            ! t* in similarity theory [K]
  real(r8), pointer :: fm(:)               ! integral of profile function for momentum
  real(r8), pointer :: fh(:)               ! integral of profile function for heat
  real(r8), pointer :: fq(:)               ! integral of profile function for moisture

! -----------------------------------------------------------------
! III. Forcing 
! -----------------------------------------------------------------

  real(r8), pointer :: pco2m            ! CO2 concentration in atmos. (35 pa)
  real(r8), pointer :: po2m             ! O2 concentration in atmos. (20900 pa)
  real(r8), pointer :: us               ! wind in eastward direction [m/s]
  real(r8), pointer :: vs               ! wind in northward direction [m/s]
  real(r8), pointer :: tm               ! temperature at reference height [kelvin]
  real(r8), pointer :: qm               ! specific humidity at reference height [kg/kg]
  real(r8), pointer :: prc              ! convective precipitation [mm/s]
  real(r8), pointer :: prl              ! large scale precipitation [mm/s]
  real(r8), pointer :: pbot             ! atm bottom level pressure (or reference height) (pa)
  real(r8), pointer :: psrf             ! atmospheric pressure at the surface [pa]
  real(r8), pointer :: sols             ! atm vis direct beam solar rad onto srf [W/m2]
  real(r8), pointer :: soll             ! atm nir direct beam solar rad onto srf [W/m2]
  real(r8), pointer :: solsd            ! atm vis diffuse solar rad onto srf [W/m2]
  real(r8), pointer :: solld            ! atm nir diffuse solar rad onto srf [W/m2]
  real(r8), pointer :: frl              ! atmospheric infrared (longwave) radiation [W/m2]
  real(r8), pointer :: hu               ! observational height of wind [m]
  real(r8), pointer :: ht               ! observational height of temperature [m]
  real(r8), pointer :: hq               ! observational height of humidity [m]

  real(r8) :: rhoair                    ! air density [kg/m3]

! -----------------------------------------------------------------
! IV. Fluxes
! -----------------------------------------------------------------

  real(r8), pointer :: taux(:)             ! wind stress: E-W [kg/m/s2]
  real(r8), pointer :: tauy(:)             ! wind stress: N-S [kg/m/s2]
  real(r8), pointer :: fsena(:)            ! sensible heat from canopy height to atmosphere [W/m2]
  real(r8), pointer :: lfevpa(:)           ! latent heat flux from canopy height to atmosphere [W/m2]
  real(r8), pointer :: fevpa(:)            ! evapotranspiration from canopy height to atmosphere [mm/s]
  real(r8), pointer :: fsenl(:)            ! sensible heat from leaves [W/m2]
  real(r8), pointer :: fevpl(:)            ! evaporation+transpiration from leaves [mm/s]
  real(r8), pointer :: etr(:)              ! transpiration rate [mm/s]
  real(r8), pointer :: fseng(:)            ! sensible heat flux from ground [W/m2]
  real(r8), pointer :: fevpg(:)            ! evaporation heat flux from ground [mm/s]
  real(r8), pointer :: fgrnd(:)            ! ground heat flux [W/m2]
  real(r8), pointer :: sabvsun(:)          ! solar absorbed by sunlit vegetation [W/m2]
  real(r8), pointer :: sabvsha(:)          ! solar absorbed by shaded vegetation [W/m2]
  real(r8), pointer :: sabg(:)             ! solar absorbed by ground  [W/m2]
  real(r8), pointer :: olrg(:)             ! outgoing long-wave radiation from ground+canopy [W/m2]
! real(r8), pointer :: rnet(:)             ! net radiation by surface [W/m2]
  real(r8), pointer :: xerr                ! the error of water banace [mm/s]
  real(r8), pointer :: zerr(:)             ! the error of energy balance [W/m2]

  real(r8), pointer :: rsur                ! surface runoff (mm h2o/s)
  real(r8), pointer :: rnof                ! total runoff (mm h2o/s)
  real(r8), pointer :: assim(:)            ! canopy assimilation rate (mol m-2 s-1)
  real(r8), pointer :: respc(:)            ! canopy respiration (mol m-2 s-1)

  real(r8), pointer :: u10m(:)             ! 10m u-velocity 
  real(r8), pointer :: v10m(:)             ! 10m v-velocity 
  real(r8), pointer :: f10m(:)             ! integral of profile function for momentum at 10m
! -----------------------------------------------------------------
! V. Local declaration
! -----------------------------------------------------------------
  integer, parameter :: npftpar=32         ! 32 parameter of each PFT 
  integer, parameter :: maxpft=17          ! number of pfts in natural vegetation land unit, including bare soil 

  real(r8) :: parsun(maxpft)               ! PAR by sunlit leaves [W/m2] 
  real(r8) :: parsha(maxpft)               ! PAR by shaded leaves [W/m2]
  real(r8) :: sabvg(maxpft)                ! solar absorbed by ground + vegetation [W/m2]

  real(r8) :: sm                           ! rate of snowmelt [kg/(m2 s)]
  real(r8) :: qsubl(maxpft)                ! sublimation rate from snow pack (mm h2o /s) [+]

  integer  :: snl                          ! number of snow layers
  real(r8) :: lwsnl                        ! mass of liquid water of snow layers [kg/m2]
  real(r8) :: nsnow                        ! number of snow events [-]
  real(r8) :: tsn                          ! snow internal temperature [K]
  real(r8) :: mrsos                        ! mass of water of all phases in the upper 0.1 meters of soil [kg/m2]
  real(r8) :: mrso                         ! mass of water of all phases over all soil layers [kg/m2]
  real(r8) :: mrfso                        ! mass of frozen water over all soil layers [kg/m2]
  real(r8) :: mrlsl(1:nl_soil)             ! mass of water of all phases in each soil layer [kg/m2]

  integer i,j,lc,uc,lb,ub,jm               ! loop/array indices
  integer c,p,p1,p2,n_pft 

! -----------------------------------------------------------------

      fLitterSoil = 0._r8
      fLitterAtmos = 0._r8

      jm = nl_soil+abs(maxsnl) 

      ! CLM time step and TUNABLE constants

      zlnd   = ftune(1)
      zsno   = ftune(2)
      csoilc = ftune(3)
      dewmx  = ftune(4)
      wtfact = ftune(5)
      capr   = ftune(6)
      cnfac  = ftune(7)
      ssi    = ftune(8)
      wimp   = ftune(9)
      pondmx = ftune(10) 
      smpmax = ftune(11) 
      smpmin = ftune(12) 
      trsmx0 = ftune(13) 
      tcrit  = ftune(14)  

      year   = idate(1)      
      jday   = idate(2)      
      msec   = idate(3)      

      p1 = 1
      p2 = 0

      DO c = 1, numcolumn

! ======================================================================
!  [1] Transfer the time invariant and time-varying variables
! ======================================================================
         ! Time invariant model variables

         uc = 1
         dlat      => fcon_col(uc,c)             ; uc = uc + 1                     !1
         dlon      => fcon_col(uc,c)             ; uc = uc + 1                     !2
         itypwat_c =  nint(fcon_col(uc,c))       ; uc = uc + 1                     !3
         albsol    => fcon_col(uc,c)             ; lc = uc + 1; uc = uc + nl_soil  !4
         csol      => fcon_col(lc:uc,c)          ; lc = uc + 1; uc = uc + nl_soil  !1_
         porsl     => fcon_col(lc:uc,c)          ; lc = uc + 1; uc = uc + nl_soil  !2_
         phi0      => fcon_col(lc:uc,c)          ; lc = uc + 1; uc = uc + nl_soil  !3_
         bsw       => fcon_col(lc:uc,c)          ; lc = uc + 1; uc = uc + nl_soil  !4_
         dkmg      => fcon_col(lc:uc,c)          ; lc = uc + 1; uc = uc + nl_soil  !5_
         dksatu    => fcon_col(lc:uc,c)          ; lc = uc + 1; uc = uc + nl_soil  !6_
         dkdry     => fcon_col(lc:uc,c)          ; lc = uc + 1; uc = uc + nl_soil  !7_
         hksati    => fcon_col(lc:uc,c)                                            !8_

#ifdef MYBUG
         write(6,"('grid:',I4,2(F10.3),I4)") p_iam, dlat*180./3.14, dlon*180./3.14, itypwat_c
#endif

         if(itypwat_c.ne.itypwat(c)) then
            write(6,*), 'fatal error on itypwat checking [CLMDRIVER]', itypwat_c, itypwat(c)
            call abort
         end if

#if(defined SINGLE)
         if(itypwat_c==0)      n_pft = 1       ! natural vegetation + bare soil
#else
         if(itypwat_c==0)      n_pft = numpft + 1  ! natural vegetation + bare soil
#endif
         if(itypwat_c==1)      n_pft = 1       ! urban and built-up
         if(itypwat_c==2)      n_pft = 1       ! wetland
         if(itypwat_c==3)      n_pft = 1       ! land ice
         if(itypwat_c==4)      n_pft = 1       ! river or deep lake
         if(itypwat_c==99)then                 ! ocean
            write(6,*) 'ocean column',c
            stop
         endif

         p2 = p2 + n_pft

         ub = 1
#if(!defined DGVM)
         ivt       => fcon_pft(ub,p1:p2)             ; ub = ub + 1                     !*
#endif
         z0m       => fcon_pft(ub,p1:p2)             ; ub = ub + 1                     !1
         displa    => fcon_pft(ub,p1:p2)             ; ub = ub + 1                     !2
         sqrtdi    => fcon_pft(ub,p1:p2)             ; ub = ub + 1                     !3
         effcon    => fcon_pft(ub,p1:p2)             ; ub = ub + 1                     !4
         vmax25    => fcon_pft(ub,p1:p2)             ; ub = ub + 1                     !5
         slti      => fcon_pft(ub,p1:p2)             ; ub = ub + 1                     !6
         hlti      => fcon_pft(ub,p1:p2)             ; ub = ub + 1                     !7
         shti      => fcon_pft(ub,p1:p2)             ; ub = ub + 1                     !8
         hhti      => fcon_pft(ub,p1:p2)             ; ub = ub + 1                     !9
         trda      => fcon_pft(ub,p1:p2)             ; ub = ub + 1                     !10
         trdm      => fcon_pft(ub,p1:p2)             ; ub = ub + 1                     !11
         trop      => fcon_pft(ub,p1:p2)             ; ub = ub + 1                     !12
         gradm     => fcon_pft(ub,p1:p2)             ; ub = ub + 1                     !13
         binter    => fcon_pft(ub,p1:p2)             ; ub = ub + 1                     !14
         extkn     => fcon_pft(ub,p1:p2)             ; ub = ub + 1                     !15
         chil      => fcon_pft(ub,p1:p2)             ; lb = ub + 1; ub = ub + 4        !16
         ref       => fcon_pft(lb:ub,p1:p2)          ; lb = ub + 1; ub = ub + 4        !17-20
         tran      => fcon_pft(lb:ub,p1:p2)          ; lb = ub + 1; ub = ub + nl_soil  !21-24
         rootfr    => fcon_pft(lb:ub,p1:p2)                                            !25-34
#if(defined DGVM)
         ub = ub + 1
         pftpar    => fcon_pft(ub:ub+npftpar-1,p1:p2)    ; ub = ub+npftpar       !35-66
         vegclass  => fcon_pft(ub,p1:p2)                 ; ub = ub+1             !67
         summergreen=>fcon_pft(ub,p1:p2)                 ; ub = ub+1             !68
         raingreen => fcon_pft(ub,p1:p2)                 ; ub = ub+1             !69
         sla       => fcon_pft(ub,p1:p2)                 ; ub = ub+1             !70
         lm_sapl   => fcon_pft(ub,p1:p2)                 ; ub = ub+1             !71
         sm_sapl   => fcon_pft(ub,p1:p2)                 ; ub = ub+1             !72
         hm_sapl   => fcon_pft(ub,p1:p2)                 ; ub = ub+1             !73
         rm_sapl   => fcon_pft(ub,p1:p2)                                         !74
#endif
#if(defined DyN)
         ub = ub + 1
         cton_soil => fcon_pft(ub,p1:p2)                 ; ub = ub+1             !75
         cton_pro  => fcon_pft(ub,p1:p2)                                         !76
#endif
                             ! Time-varying variables 
         lc = 1
         uc = jm
         z         => fvar_col(lc:uc,c)                  ; lc = uc+1; uc = uc+jm !1_
         dz        => fvar_col(lc:uc,c)                  ; lc = uc+1; uc = uc+jm !2_
         tss       => fvar_col(lc:uc,c)                  ; lc = uc+1; uc = uc+jm !3_
         wliq      => fvar_col(lc:uc,c)                  ; lc = uc+1; uc = uc+jm !4_
         wice      => fvar_col(lc:uc,c)                  ; uc = uc+1             !5_
         tg        => fvar_col(uc,c)                     ; uc = uc+1             !1
         sag       => fvar_col(uc,c)                     ; uc = uc+1             !2
         scv       => fvar_col(uc,c)                     ; uc = uc+1             !3
         snowdp    => fvar_col(uc,c)                     ; uc = uc+1             !4
         fsno      => fvar_col(uc,c)                     ; uc = uc+1             !5
         coszen    => fvar_col(uc,c)                                             !6

         if(scv<0.)then
            write(6,*),'driversnow',itypwat_c,tg,scv,snowdp
         endif

#if(defined DGVM)
         uc = uc+1
         nday      => fvar_col(uc,c)                     ; uc = uc+1             !1
         nyr       => fvar_col(uc,c)                     ; uc = uc+1             !2
         prec365   => fvar_col(uc,c)                                             !3
#endif
#if(defined DyN)
         uc = uc+1
         soil_no3  => fvar_col(uc,c)                   ; uc = uc+1             !1
         soil_no2  => fvar_col(uc,c)                   ; uc = uc+1             !2
         soil_no   => fvar_col(uc,c)                   ; uc = uc+1             !3
         soil_n2o  => fvar_col(uc,c)                   ; uc = uc+1             !4
         soil_n2   => fvar_col(uc,c)                   ; uc = uc+1             !5
         soil_nh4  => fvar_col(uc,c)                                           !6
#endif
         ub = 1
         tlsun     => fvar_pft(ub,p1:p2)                 ; ub = ub+1             !1
         tlsha     => fvar_pft(ub,p1:p2)                 ; ub = ub+1             !2
         ldew      => fvar_pft(ub,p1:p2)                 ; ub = ub+1             !3
         fveg      => fvar_pft(ub,p1:p2)                 ; ub = ub+1             !4
         sigf      => fvar_pft(ub,p1:p2)                 ; ub = ub+1             !5
         green     => fvar_pft(ub,p1:p2)                 ; ub = ub+1             !6
         lai       => fvar_pft(ub,p1:p2)                 ; ub = ub+1             !7
         sai       => fvar_pft(ub,p1:p2)                 ; lb = ub+1; ub = ub+4  !8
         albg      => fvar_pft(lb:ub,p1:p2)              ; lb = ub+1; ub = ub+4  !9-12
         albv      => fvar_pft(lb:ub,p1:p2)              ; lb = ub+1; ub = ub+4  !13-16
         alb       => fvar_pft(lb:ub,p1:p2)              ; lb = ub+1; ub = ub+4  !17-20
         ssun      => fvar_pft(lb:ub,p1:p2)              ; lb = ub+1; ub = ub+4  !21-24
         ssha      => fvar_pft(lb:ub,p1:p2)              ; ub = ub+1             !25-28
         thermk    => fvar_pft(ub,p1:p2)                 ; ub = ub+1             !29
         extkb     => fvar_pft(ub,p1:p2)                 ; ub = ub+1             !30
         extkd     => fvar_pft(ub,p1:p2)                 ; ub = ub+1             !31

                    ! Additional variables required by reginal model (WRF & RSM) 
         trad      => fvar_pft(ub,p1:p2)                 ; ub = ub+1             !32
         tref      => fvar_pft(ub,p1:p2)                 ; ub = ub+1             !33
         qref      => fvar_pft(ub,p1:p2)                 ; ub = ub+1             !34
         rst       => fvar_pft(ub,p1:p2)                 ; ub = ub+1             !35
         emis      => fvar_pft(ub,p1:p2)                 ; ub = ub+1             !36
         z0ma      => fvar_pft(ub,p1:p2)                 ; ub = ub+1             !37
         zol       => fvar_pft(ub,p1:p2)                 ; ub = ub+1             !38
         rib       => fvar_pft(ub,p1:p2)                 ; ub = ub+1             !39
         ustar     => fvar_pft(ub,p1:p2)                 ; ub = ub+1             !40
         qstar     => fvar_pft(ub,p1:p2)                 ; ub = ub+1             !41
         tstar     => fvar_pft(ub,p1:p2)                 ; ub = ub+1             !42
         fm        => fvar_pft(ub,p1:p2)                 ; ub = ub+1             !43
         fh        => fvar_pft(ub,p1:p2)                 ; ub = ub+1             !44
         fq        => fvar_pft(ub,p1:p2)                                         !45
#if(defined DGVM)
         ub = ub+1
         t10min    => fvar_pft(ub,p1:p2)                 ; ub = ub+1             !1
         lai_ind   => fvar_pft(ub,p1:p2)                 ; ub = ub+1             !2
         dphen     => fvar_pft(ub,p1:p2)                 ; ub = ub+1             !3
         leafon    => fvar_pft(ub,p1:p2)                 ; ub = ub+1             !4
         leafof    => fvar_pft(ub,p1:p2)                 ; ub = ub+1             !5
         firelength=> fvar_pft(ub,p1:p2)                 ; ub = ub+1             !6
         litterag  => fvar_pft(ub,p1:p2)                 ; ub = ub+1             !7
         litterbg  => fvar_pft(ub,p1:p2)                 ; ub = ub+1             !8
         cpool_fast=> fvar_pft(ub,p1:p2)                 ; ub = ub+1             !9
         cpool_slow=> fvar_pft(ub,p1:p2)                 ; ub = ub+1             !10
         k_fast_ave=> fvar_pft(ub,p1:p2)                 ; ub = ub+1             !11
         k_slow_ave=> fvar_pft(ub,p1:p2)                 ; ub = ub+1             !12
         litter_decom_ave=> fvar_pft(ub,p1:p2)           ; ub = ub+1             !13
         fmicr     => fvar_pft(ub,p1:p2)                 ; ub = ub+1             !14
         nind      => fvar_pft(ub,p1:p2)                 ; ub = ub+1             !15
         lm_ind    => fvar_pft(ub,p1:p2)                 ; ub = ub+1             !16
         sm_ind    => fvar_pft(ub,p1:p2)                 ; ub = ub+1             !17
         hm_ind    => fvar_pft(ub,p1:p2)                 ; ub = ub+1             !18
         rm_ind    => fvar_pft(ub,p1:p2)                 ; ub = ub+1             !19
         tmomin20  => fvar_pft(ub,p1:p2)                 ; ub = ub+1             !20
         agdd0     => fvar_pft(ub,p1:p2)                 ; ub = ub+1             !21
         agdd      => fvar_pft(ub,p1:p2)                 ; ub = ub+1             !22
         agdd20    => fvar_pft(ub,p1:p2)                 ; ub = ub+1             !23
         t_mo_min  => fvar_pft(ub,p1:p2)                 ; ub = ub+1             !24
         crownarea => fvar_pft(ub,p1:p2)                 ; ub = ub+1             !25
         htop      => fvar_pft(ub,p1:p2)                 ; ub = ub+1             !26
         tsai      => fvar_pft(ub,p1:p2)                 ; ub = ub+1             !27
         fpcgrid   => fvar_pft(ub,p1:p2)                 ; ub = ub+1             !28
         bm_inc    => fvar_pft(ub,p1:p2)                 ; ub = ub+1             !29
         afmicr    => fvar_pft(ub,p1:p2)                 ; ub = ub+1             !30
         annpsn    => fvar_pft(ub,p1:p2)                 ; ub = ub+1             !31
         annpsnpot => fvar_pft(ub,p1:p2)                 ; ub = ub+1             !32
         tref10    => fvar_pft(ub,p1:p2)                 ; ub = ub+1             !33
         tref_sum  => fvar_pft(ub,p1:p2)                 ; ub = ub+1             !34
         t10       => fvar_pft(ub:ub+9,p1:p2)            ; ub = ub+10            !35-44
         assimn10  => fvar_pft(ub,p1:p2)                 ; ub = ub+1             !45
         assimn_sum=> fvar_pft(ub,p1:p2)                 ; ub = ub+1             !46
         an10      => fvar_pft(ub:ub+9,p1:p2)            ; ub = ub+10            !47-56
         turnover_ind => fvar_pft(ub,p1:p2)              ; ub = ub+1             !57
         fpc_inc   => fvar_pft(ub,p1:p2)                 ; ub = ub+1             !58
         ivt       => fvar_pft(ub,p1:p2)                 ; ub = ub+1             !59
         agddtw    => fvar_pft(ub,p1:p2)                 ; ub = ub+1             !60
         ifpre     => fvar_pft(ub,p1:p2)                 ; ub = ub+1             !61 
         t_mo      => fvar_pft(ub,p1:p2)                 ; ub = ub+1             !62
         t_mo_sum  => fvar_pft(ub,p1:p2)                 ; ub = ub+1             !63 
         anngpp    => fvar_pft(ub,p1:p2)                 ; ub = ub+1             !64
         annfrmf   => fvar_pft(ub,p1:p2)                 ; ub = ub+1             !65
         annfrms   => fvar_pft(ub,p1:p2)                 ; ub = ub+1             !66
         annfrmr   => fvar_pft(ub,p1:p2)                 ; ub = ub+1             !67
         annfrg    => fvar_pft(ub,p1:p2)                                         !68
#if(defined DyN)
         ub=ub+1
         litter_leaf       =>fvar_pft(ub,p1:p2)          ; ub = ub+1             !1
         litter_wood       =>fvar_pft(ub,p1:p2)          ; ub = ub+1             !2
         litter_root       =>fvar_pft(ub,p1:p2)          ; ub = ub+1             !3
         litter_repr       =>fvar_pft(ub,p1:p2)          ; ub = ub+1             !4
         litter_leaf_n     =>fvar_pft(ub,p1:p2)          ; ub = ub+1             !5
         litter_wood_n     =>fvar_pft(ub,p1:p2)          ; ub = ub+1             !6
         litter_root_n     =>fvar_pft(ub,p1:p2)          ; ub = ub+1             !7
         litter_repr_n     =>fvar_pft(ub,p1:p2)          ; ub = ub+1             !8
         afcton_leaf       =>fvar_pft(ub,p1:p2)          ; ub = ub+1             !9
         afcton_sap        =>fvar_pft(ub,p1:p2)          ; ub = ub+1             !10
         afcton_root       =>fvar_pft(ub,p1:p2)          ; ub = ub+1             !11
         lm_ind_n          =>fvar_pft(ub,p1:p2)          ; ub = ub+1             !12
         sm_ind_n          =>fvar_pft(ub,p1:p2)          ; ub = ub+1             !13
         hm_ind_n          =>fvar_pft(ub,p1:p2)          ; ub = ub+1             !14
         rm_ind_n          =>fvar_pft(ub,p1:p2)          ; ub = ub+1             !15
         an_up             =>fvar_pft(ub,p1:p2)          ; ub = ub+1             !16
         an_stress         =>fvar_pft(ub,p1:p2)                                  !17
#endif
         cflux_litter_soil =>fLitterSoil (p1:p2)
         cflux_litter_atmos=>fLitterAtmos(p1:p2)
#endif

#if(defined VEGDATA)
         call interp_vegdata(itypwat_c,p1,p2,lai_r,sai_r)
#endif

! ======================================================================
!  [2] atmospheric fields to force clm
! ======================================================================

         pco2m     => forc( 1,c)       !CO2 concentration in atmos. (35 pa)
         po2m      => forc( 2,c)       !O2 concentration in atmos. (20900 pa)
         us        => forc( 3,c)       !wind in eastward direction [m/s]
         vs        => forc( 4,c)       !wind in northward direction [m/s]
         tm        => forc( 5,c)       !temperature at reference height [kelvin]
         qm        => forc( 6,c)       !specific humidity at reference height [kg/kg]
         prc       => forc( 7,c)       !convective precipitation [mm/s]
         prl       => forc( 8,c)       !large scale precipitation [mm/s]
         pbot      => forc( 9,c)       !atm bottom level pressure (or reference height) (pa)
         psrf      => forc(10,c)       !atmospheric pressure at the surface [pa]
         sols      => forc(11,c)       !atm vis direct beam solar rad onto srf [W/m2]
         soll      => forc(12,c)       !atm nir direct beam solar rad onto srf [W/m2]
         solsd     => forc(13,c)       !atm vis diffuse solar rad onto srf [W/m2]
         solld     => forc(14,c)       !atm nir diffuse solar rad onto srf [W/m2]
         frl       => forc(15,c)       !atmospheric infrared (longwave) radiation [W/m2]
         hu        => forc(16,c)       !observational height of wind [m]
         ht        => forc(17,c)       !observational height of temperature [m]
         hq        => forc(18,c)       !observational height of humidity [m]

         rhoair    =  (pbot-0.378*qm*pbot/(0.622+0.378*qm))/(rgas*tm)
!print*,'force data', tm, prl

! ======================================================================
!  [3] surface fluxes    add by zhq 03/06/2009
! ======================================================================

         xerr     => fldv_col(1,c)              !the error of water banace [mm/s]
         rsur     => fldv_col(2,c)              !surface runoff [mm/s]
         rnof     => fldv_col(3,c)              !total runoff [mm/s]

         taux     => fldv_pft(1,p1:p2)          !wind stress: E-W [kg/m/s2]
         tauy     => fldv_pft(2,p1:p2)          !wind stress: N-S [kg/m/s2]
         fsena    => fldv_pft(3,p1:p2)          !sensible heat from canopy height to atmosphere [W/m2]
         lfevpa   => fldv_pft(4,p1:p2)          !latent heat flux from canopy height to atmosphere [W/m2]
         fevpa    => fldv_pft(5,p1:p2)          !evapotranspiration from canopy height to atmosphere [mm/s]
         fsenl    => fldv_pft(6,p1:p2)          !sensible heat from leaves [W/m2]
         fevpl    => fldv_pft(7,p1:p2)          !evaporation+transpiration from leaves [mm/s]
         etr      => fldv_pft(8,p1:p2)          !transpiration rate [mm/s]
         fseng    => fldv_pft(9,p1:p2)          !sensible heat flux from ground [W/m2]
         fevpg    => fldv_pft(10,p1:p2)         !evaporation heat flux from ground [mm/s]
         fgrnd    => fldv_pft(11,p1:p2)         !ground heat flux [W/m2]
         sabvsun  => fldv_pft(12,p1:p2)         !solar absorbed by sunlit canopy [W/m2]
         sabvsha  => fldv_pft(13,p1:p2)         !solar absorbed by shaded [W/m2]
         sabg     => fldv_pft(14,p1:p2)         !solar absorbed by ground  [W/m2]
         olrg     => fldv_pft(15,p1:p2)         !outgoing long-wave radiation from ground+canopy [W/m2]
         zerr     => fldv_pft(17,p1:p2)         !the error of energy balance [W/m2]
         assim    => fldv_pft(18,p1:p2)         !canopy assimilation rate [mol m-2 s-1]
         respc    => fldv_pft(19,p1:p2)         !respiration (plant+soil) [mol m-2 s-1]
         u10m     => fldv_pft(45,p1:p2)         !
         v10m     => fldv_pft(46,p1:p2)         !
         f10m     => fldv_pft(47,p1:p2)         !

         wt_pft   => wxy_patch(p1:p2)
         wt_col   => wxy_column(c)

#ifdef MYBUG
i=gxmap(cgmap(c))
j=gymap(cgmap(c))
if(i.eq.108 .and. j.eq.32) then
   mybug = .true.
else
   mybug = .false.
endif

if(p_master .and. c.eq.c_bug) then
   mybug = .true.
else
   mybug = .false.
endif

if(sum(lai(1:n_pft)).gt.1.0E-6 .and. dlat/3.1415926535897932*180..lt.-60.) then
   write(6,*) "lat,lon", dlat/3.1416*180., dlon/3.1416*180.
   write(6,*) "col", c, itypwat_c, wt_col
   write(6,*) "wt_pft", wt_pft(1:n_pft)
   write(6,*) "lai", lai(1:n_pft)
   write(6,*) "ivt", ivt(1:n_pft)
endif
if(p_master .and. itypwat_c.eq.0) then
   write(6,*) 'wt_pft', wt_pft
   write(6,*) 'lai_r', lai_r
end if
#endif

! ======================================================================
!  [2] Main driver for CLM
! ======================================================================

         CALL CLMMAIN (dtime, doalb, dolai, dosst, &
         nl_soil, maxsnl, dlon, dlat, itypwat_c, n_pft, wt_col, wt_pft, ivt, oro(c), &
   
       ! soil information
         albsol, csol, porsl, phi0, &
         bsw, dkmg, dksatu, dkdry, hksati, &
   
       ! vegetation information
         z0m, displa, sqrtdi, &
         effcon, vmax25, slti, hlti, shti, hhti, &
         trda, trdm, trop, gradm, binter, extkn, &
         chil, ref, tran, rootfr, &
#if (defined DGVM)
         pftpar, vegclass, summergreen, raingreen, &
         sla, lm_sapl, sm_sapl, hm_sapl, rm_sapl,& 
         t10min, lai_ind,&
         dphen,leafon,leafof,firelength,litterag, &
         litterbg,cpool_fast,cpool_slow,k_fast_ave, &
         k_slow_ave,litter_decom_ave,fmicr, &
         ifpre,prec365,nind,lm_ind,sm_ind, agddtw, &
         hm_ind,rm_ind,tmomin20,agdd0,agdd,agdd20, &
         t_mo,t_mo_sum,t_mo_min,crownarea,htop,tsai,fpcgrid,&
         bm_inc,afmicr,annpsn,annpsnpot,tref10,&
         tref_sum,t10,assimn10,assimn_sum,an10,&
         anngpp,annfrmf,annfrms,annfrmr,annfrg,&
         cflux_litter_soil,cflux_litter_atmos,&
         nday,nyr,&
#endif   
#if (defined DyN)
         cton_soil, cton_pro, afcton_leaf,afcton_sap,afcton_root,&
         litter_leaf, litter_wood, litter_root, litter_repr,&
         litter_leaf_n, litter_wood_n, litter_root_n, litter_repr_n,&
         an_up, an_stress, soil_no3, soil_no2, soil_no, soil_n2o,&
         soil_n2, soil_nh4,&
#endif

       ! atmospheric forcing
         frl, sols, soll, solsd, solld, &
         pco2m, po2m, us, vs, tm, qm, &
         prc, prl, psrf, &
         rhoair, &
         hu, ht, hq, &
   
       ! model variables needed by restart run
         year, jday, msec, &
         z, dz, tss, &
         wliq, wice, &
         tg, tlsun, tlsha, ldew, &
         sag, scv, snowdp, &
         fveg, fsno, sigf, green, lai, sai, &
         coszen, albg, albv, alb, &
         ssun, ssha, thermk, extkb, extkd, &
   
       ! fluxes
         taux, tauy, &
         fsena, fevpa, lfevpa, fsenl, fevpl, etr, &
         fseng, fevpg, olrg, fgrnd, trad, tref, qref, &
         rsur, rnof, rst, assim, respc, &
         parsun(1:n_pft), parsha(1:n_pft), sabvsun, sabvsha, sabg, sabvg(1:n_pft), &
         xerr, zerr, &
         sm, qsubl(1:n_pft), &
   
       ! TUNABLE modle constants
         zlnd, zsno, csoilc, dewmx, wtfact, &
         capr, cnfac, ssi, wimp, pondmx, &
         smpmax, smpmin, trsmx0, tcrit, &
   
       ! time-varying vegetation from read-in file
#if (defined VEGDATA)
         lai_r(1:n_pft), sai_r(1:n_pft), green_r(1:n_pft), fveg_r(1:n_pft), &
#endif
   
       ! additional variables required by coupling with regional model (WRF & RSM) 
         emis, z0ma, zol, rib, ustar, qstar, tstar, &
         u10m, v10m, f10m, fm, fh, fq )

#ifdef MYBUG
if(p_master .and. itypwat_c.eq.0) write(6,*) 'lai', lai(:)
#endif

! ======================================================================
!  [4] Return required surface flux fields 
! ======================================================================

                                         lc =  3+1; uc =  3+nl_soil 
         fldv_col(lc:uc,c) = tss (abs(maxsnl)+1:abs(maxsnl)+nl_soil)    !4_13
                                         lc = uc+1; uc = uc+nl_soil
         fldv_col(lc:uc,c) = wliq(abs(maxsnl)+1:abs(maxsnl)+nl_soil)    !14_23
                                         lc = uc+1; uc = uc+nl_soil
         fldv_col(lc:uc,c) = wice(abs(maxsnl)+1:abs(maxsnl)+nl_soil)    !24_33
                                                    uc = uc+1       
         fldv_col(uc,c) = tg                      ; uc = uc+1           !34
         fldv_col(uc,c) = scv                     ; uc = uc+1           !35
         fldv_col(uc,c) = snowdp                  ; uc = uc+1           !36
         fldv_col(uc,c) = fsno                    ; uc = uc+1           !37

       ! forcing
         fldv_col(uc,c) = us                      ; uc = uc+1           !38
         fldv_col(uc,c) = vs                      ; uc = uc+1           !39
         fldv_col(uc,c) = tm                      ; uc = uc+1           !40
         fldv_col(uc,c) = qm                      ; uc = uc+1           !41
         fldv_col(uc,c) = prc                     ; uc = uc+1           !42
         fldv_col(uc,c) = prl                     ; uc = uc+1           !43
         fldv_col(uc,c) = pbot                    ; uc = uc+1           !44
         fldv_col(uc,c) = frl                     ; uc = uc+1           !45
         fldv_col(uc,c) = sols + soll + solsd + solld                   !46

         fldv_pft(16,p1:p2) = sabvg(1:n_pft) + frl - olrg(1:n_pft)      !16

                                                ; ub = 19+1 
         fldv_pft(ub,p1:p2) = fmicr(:)          ; ub = ub+1             !20
         fldv_pft(ub,p1:p2) = tlsun(:)          ; ub = ub+1             !21
         fldv_pft(ub,p1:p2) = tlsha(:)          ; ub = ub+1             !22
         fldv_pft(ub,p1:p2) = ldew(:)           ; ub = ub+1             !23
         fldv_pft(ub,p1:p2) = sigf(:)           ; ub = ub+1             !24
         fldv_pft(ub,p1:p2) = green(:)          ; ub = ub+1             !25
         fldv_pft(ub,p1:p2) = lai(:)            ; ub = ub+1             !26
         fldv_pft(ub,p1:p2) = sai(:)            ; ub = ub+1             !27
         fldv_pft(ub,p1:p2) = alb(1,:)          ; ub = ub+1             !28
         fldv_pft(ub,p1:p2) = alb(3,:)          ; ub = ub+1             !29
         fldv_pft(ub,p1:p2) = alb(2,:)          ; ub = ub+1             !30
         fldv_pft(ub,p1:p2) = alb(4,:)          ; ub = ub+1             !31
         fldv_pft(ub,p1:p2) = emis(:)           ; ub = ub+1             !32
         fldv_pft(ub,p1:p2) = z0ma(:)           ; ub = ub+1             !33

         fldv_pft(ub,p1:p2) = trad(:)           ; ub = ub+1             !34
         fldv_pft(ub,p1:p2) = ustar(:)          ; ub = ub+1             !35
         fldv_pft(ub,p1:p2) = tstar(:)          ; ub = ub+1             !36
         fldv_pft(ub,p1:p2) = qstar(:)          ; ub = ub+1             !37

         fldv_pft(ub,p1:p2) = zol(:)            ; ub = ub+1             !38
         fldv_pft(ub,p1:p2) = rib(:)            ; ub = ub+1             !39
         fldv_pft(ub,p1:p2) = fm(:)             ; ub = ub+1             !40
         fldv_pft(ub,p1:p2) = fh(:)             ; ub = ub+1             !41
         fldv_pft(ub,p1:p2) = fq(:)             ; ub = ub+1             !42

       ! diagnostic variables
         fldv_pft(ub,p1:p2) = tref(:)           ; ub = ub+1             !43
         fldv_pft(ub,p1:p2) = qref(:)           ; ub = ub+1             !44

#ifdef CMIP

         snl = 0
         do j = 1, abs(maxsnl)
            if((wliq(j)+wice(j))>0.) snl = snl + 1
         end do

         lb = abs(maxsnl)-snl+1
         ub = abs(maxsnl)

       ! liquid water content of snow layer [lwsnl:kg/m2]
         if (snl>0) then
            lwsnl = sum(wliq(lb:ub))
         else
            lwsnl = 0._r8
         end if

       ! snow age [agesno:day]
       ! time samples weighted by snow mass and accumulate
       ! finally divided by ave(scv)

       ! agesno = dtime*scv

       ! snow internal temperature [tsn:K]

         if (snl>0) then
            tsn = sum(tss(lb:ub)*dz(lb:ub))/sum(dz(lb:ub))
            nsnow = 1._r8
         else
            tsn = 0._r8
            nsnow = 0._r8
         end if

         lb = abs(maxsnl)+1
         ub = abs(maxsnl)+3

         mrsos = sum(wliq(lb:ub)) + sum(wice(lb:ub))

         lb = abs(maxsnl)+1
         ub = abs(maxsnl)+nl_soil

         mrso  = sum(wliq(lb:ub)) + sum(wice(lb:ub))

         mrfso = sum(wice(lb:ub))

         mrlsl = wice(lb:ub) + wliq(lb:ub)

       ! --------------------------------PFT LEVEL-------------------------------- !

                                                ; ub = 47+1
         fldv_pft(ub,p1:p2) = qsubl(1:n_pft)    ; ub = ub+1         !48


       ! -----------------------------COLUMN LEVEL-------------------------------- !

         lc = 46+1                              ; uc = 46+nl_soil

         fldv_col(lc:uc,c) = mrlsl(1:nl_soil)   ; uc = uc+1         !47-56
         fldv_col(uc,c) = mrsos                 ; uc = uc+1         !57
         fldv_col(uc,c) = mrso                  ; uc = uc+1         !58
         fldv_col(uc,c) = mrfso                 ; uc = uc+1         !59
         fldv_col(uc,c) = lwsnl                 ; uc = uc+1         !60
         fldv_col(uc,c) = sm                    ; uc = uc+1         !61
         fldv_col(uc,c) = tsn                   ; uc = uc+1         !62
         fldv_col(uc,c) = nsnow                 ; uc = uc+1         !63

#endif
         p1 = p2 + 1

      ENDDO

 end subroutine CLMDRIVER
