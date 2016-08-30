#include <define.h>

 subroutine CLMMAIN (dtime, doalb, dolai, dosst, &
            nl_soil, maxsnl, dlon, dlat, itypwat, n_pft, wt_column, wt_patch, ivt_r, oro, &

          ! soil information
            albsol, csol, porsl, phi0, bsw, &
            dkmg, dksatu, dkdry, hksati, &

          ! vegetation information
            z0m_r, displa_r, sqrtdi, &
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
            ifpre,prec365,nind,lm_ind,sm_ind,agddtw, &
            hm_ind,rm_ind,tmomin20,agdd0,agdd,agdd20, &
            t_mo,t_mo_sum,t_mo_min,crownarea,htop,tsai,&
            fpcgrid,bm_inc,afmicr,annpsn,annpsnpot,tref10,&
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

          ! land surface variables required for restart
            year, jday, msec, &
            z, dz, tss, wliq, wice, &
            tg, tlsun, tlsha, ldew, sag, scv, snowdp, &
            fveg, fsno, sigf, green, lai, sai, &
            coszen, albg, albv, alb, ssun, ssha, thermk, extkb, extkd, &

          ! fluxes
            taux,  tauy, &
            fsena, fevpa, lfevpa, fsenl, fevpl, etr, &
            fseng, fevpg, olrg, fgrnd, trad, tref, qref, &
            rsur, rnof, rst, assim, respc, &
            parsun, parsha, sabvsun, sabvsha, sabg, sabvg, xerr, zerr, &
            sm, qsubl, &

          ! TUNABLE modle constants
            zlnd,   zsno,   csoilc, dewmx,  wtfact, & 
            capr,   cnfac,  ssi,    wimp,   pondmx, &  
            smpmax, smpmin, trsmx0, tcrit, & 

          ! time-varying vegetation parameters from read-in file
#if(defined VEGDATA)
            lai_r, sai_r, green_r, fveg_r, &
#endif

          ! additional variables required by coupling with WRF model 
            emis, z0ma, zol, rib, ustar, qstar, tstar, &
            u10m, v10m, f10m, fm, fh, fq )

!=======================================================================
!
! Main subroutine, advance time information
!
! Original author : Yongjiu Dai, 09/30/1999; 08/30/2002
! New structure modificated by: zhq 06/06/2009
!
!    FLOW DIAGRAM FOR CLMMAIN
!
!    CLMMAIN ===> netsolar                       |> all surface
!
!                 leafinterception               |]  
!                 newsnow                        |] itypwat = 0 (soil ground) 
!                 THERMAL                        |]         = 1 (urban & built-up) 
!                 WATER                          |]         = 2 (wetland) 
!                 snowcompaction                 |]         = 3 (land ice) 
!                 snowlayerscombine              |]         = 4 (deep lake)
!                 snowlayersdivide               |]         = 5 (shallow lake)
!                 snowage                        |]         = 99(ocean)
!
!                 LAKE                           |> lake 
!                 SOCEAN                         |> ocean and sea ice
!                
!                 orb_coszen                     |> all surface
!                 EcoModel (DGVM/lai_empirical)  |> land 
!                 snowfraction                   |> land
!                 albland                        |> land 
!                 albocean                       |> ocean & sea ice
!
!=======================================================================

  use precision
  use paramodel, only : numpft
  use phycon_module, only : tfrz 
  use DGVMtimevar
  use SubgridMod, only : p2c, BuildFilter
  use debug
 
  implicit none
 
! ------------------------ Dummy Argument ------------------------------
  integer, INTENT(in) :: &
        nl_soil    , &! number of soil layers
        maxsnl     , &! maximum number of snow layers
        n_pft      , &! plant function type in a single column add by zhq. 06/02/2009
        itypwat       ! land water type (0=soil, 1=urban and built-up, 
                      ! 2=wetland, 3=land ice, 4=deep lake, 5=shallow lake, 99 = ocean)

  logical, INTENT(in) :: doalb   !true if time for surface albedo calculation
  logical, INTENT(in) :: dolai   !true if time for leaf area index calculation
  logical, INTENT(in) :: dosst   !true to update sst/ice/snow before calculation

! Parameters
! ----------------------
  real(r8), INTENT(in) :: &
        dtime      , &!model time step [second]
        zlnd       , &!roughness length for soil [m]
        zsno       , &!roughness length for snow [m]
        csoilc     , &!drag coefficient for soil under canopy [-]
        dewmx      , &!maximum dew
        wtfact     , &!fraction of model area with high water table
        capr       , &!tuning factor to turn first layer T into surface T
        cnfac      , &!Crank Nicholson factor between 0 and 1
        ssi        , &!irreducible water saturation of snow
        wimp       , &!water impremeable if porosity less than wimp
        pondmx     , &!ponding depth (mm)
        smpmax     , &!wilting point potential in mm
        smpmin     , &!restriction for min of soil poten.  (mm)
        trsmx0     , &!max transpiration for moist soil+100% veg.  [mm/s]
        tcrit      , &!critical temp. to determine rain or snow

        dlon       , &! logitude in radians
        dlat       , &! latitude in radians

        wt_column  , &! weight of column relative to grid area
        !wt_patch(n_pft), &! weight of pft relative to grid area
        
        ! soil physical parameters
        albsol         , &! soil albedo for different coloured soils [-]
        csol(nl_soil)  , &! heat capacity of soil solids [J/(m3 K)]
        porsl(nl_soil) , &! fraction of soil that is voids [-]
        phi0(nl_soil)  , &! minimum soil suction [mm]
        bsw(nl_soil)   , &! clapp and hornbereger "b" parameter [-]
        dkmg(nl_soil)  , &! thermal conductivity of soil minerals [W/m-K]
        dksatu(nl_soil), &! thermal conductivity of saturated soil [W/m-K]
        dkdry(nl_soil) , &! thermal conductivity for dry soil  [J/(K s m)]
        hksati(nl_soil), &! hydraulic conductivity at saturation [mm h2o/s]

        ! vegetation static, dynamic, derived parameters
        ivt_r(n_pft)      , &! land cover type add by zhq. 06/02/2009
        z0m_r(n_pft)      , &! aerodynamic roughness length [m] zhq: different in DGVM [-]
        displa_r(n_pft)   , &! displacement height [m] zhq: different in DGVM [-]
        sqrtdi(n_pft)     , &! inverse sqrt of leaf dimension [m**-0.5]

        effcon(n_pft)     , &! quantum efficiency of RuBP regeneration (mol CO2/mol quanta)
        vmax25(n_pft)     , &! maximum carboxylation rate at 25 C at canopy top
        slti(n_pft)       , &! slope of low temperature inhibition function      [s3] 
        hlti(n_pft)       , &! 1/2 point of low temperature inhibition function  [s4]
        shti(n_pft)       , &! slope of high temperature inhibition function     [s1]
        hhti(n_pft)       , &! 1/2 point of high temperature inhibition function [s2]
        trda(n_pft)       , &! temperature coefficient in gs-a model             [s5]
        trdm(n_pft)       , &! temperature coefficient in gs-a model             [s6]
        trop(n_pft)       , &! temperature coefficient in gs-a model          
        gradm(n_pft)      , &! conductance-photosynthesis slope parameter
        binter(n_pft)     , &! conductance-photosynthesis intercep
        extkn(n_pft)      , &! coefficient of leaf nitrogen allocation

        chil(n_pft)       , &! leaf angle distribution factor
        ref(2,2,n_pft)    , &! leaf reflectance (iw=iband, il=life and dead)
        tran(2,2,n_pft)   , &! leaf transmittance (iw=iband, il=life and dead)
        rootfr(nl_soil,n_pft) ! fraction of roots in each soil layer

! Forcing
! ----------------------
  real(r8), INTENT(in) :: &
        frl        , &! atmospheric infrared (longwave) radiation [W/m2]
        sols       , &! atm vis direct beam solar rad onto srf [W/m2]
        soll       , &! atm nir direct beam solar rad onto srf [W/m2]
        solsd      , &! atm vis diffuse solar rad onto srf [W/m2]
        solld      , &! atm nir diffuse solar rad onto srf [W/m2]

        pco2m      , &! partial pressure of CO2 at observational height [pa]
        po2m       , &! partial pressure of O2 at observational height [pa]
        us         , &! wind speed in eastward direction [m/s]
        vs         , &! wind speed in northward direction [m/s]
        tm         , &! temperature at agcm reference height [kelvin]
        qm         , &! specific humidity at agcm reference height [kg/kg]
        prc        , &! convective precipitation [mm/s]
        prl        , &! large scale precipitation [mm/s]
        psrf       , &! atmosphere pressure at the surface [pa]
        rhoair     , &! density air [kg/m3]
        hu         , &! observational height of wind [m]
        ht         , &! observational height of temperature [m]
        hq            ! observational height of humidity [m]
#if(defined DGVM)
  real(r8),INTENT(in)  :: &
        pftpar(32,n_pft)      , &!PFT biological and bioclimate limitation parameters
        vegclass(n_pft)       , &!1.tree 2.shrub 3.grass 4.crop -1.others
        summergreen(n_pft)    , &!1. for summergreen tree and shrub
        raingreen(n_pft)      , &!1. for raingreen tree and shrub
        sla(n_pft)            , &!specialized leaf area
        lm_sapl(n_pft)        , &!sapling leafmass
        sm_sapl(n_pft)        , &!sapling sapwood mass
        hm_sapl(n_pft)        , &!sapling heartwood mass
        rm_sapl(n_pft)           !sapling rootmass
#endif        

#if(defined DyN)
  real(r8),INTENT(in) :: cton_soil(n_pft)        ! soil C:N mass ratio
  real(r8),INTENT(in) :: cton_pro (n_pft)        ! C:N mass ratio in production
#endif

#if(defined VEGDATA)
  real(r8),INTENT(in) :: lai_r(n_pft), sai_r(n_pft), green_r(n_pft), fveg_r(n_pft)
#endif

! Variables required for restart run
! ----------------------------------------------------------------------
  integer, INTENT(in) :: &
        year,jday,msec ! next time-step /year/julian day/second in a day/

  real(r8), INTENT(inout) :: oro  ! ocean(0)/seaice(2)/ flag
  real(r8), INTENT(inout) :: &
        wt_patch(n_pft)       , &! weight of pft relative to grid area 
        z(maxsnl+1:nl_soil)   , &! layer depth (m)
        dz(maxsnl+1:nl_soil)  , &! layer thickness (m)
        tss(maxsnl+1:nl_soil) , &! soil + snow layer temperature [K]
        wliq(maxsnl+1:nl_soil), &! liquid water (kg/m2)
        wice(maxsnl+1:nl_soil), &! ice lens (kg/m2)

        tg               , &! ground surface temperature [k]
        tlsun    (n_pft) , &! sunlit leaf temperature [K]
        tlsha    (n_pft) , &! shaded leaf temperature [K]
        ldew     (n_pft) , &! depth of water on foliage [kg/m2/s]
        sag              , &! non dimensional snow age [-]
        scv              , &! snow mass (kg/m2)
        snowdp           , &! snow depth (m)

        fveg     (n_pft) , &! fraction of vegetation cover
        fsno             , &! fractional snow cover
        sigf     (n_pft) , &! fraction of veg cover, excluding snow-covered veg [-]
        green    (n_pft) , &! greenness
        lai      (n_pft) , &! leaf area index
        sai      (n_pft) , &! stem area index
 
        coszen           , &! cosine of solar zenith angle
        albg (2,2,n_pft) , &! albedo, ground [-]
        albv (2,2,n_pft) , &! albedo, vegetation [-]
        alb  (2,2,n_pft) , &! averaged albedo [-]
        ssun (2,2,n_pft) , &! sunlit canopy absorption for solar radiation
        ssha (2,2,n_pft) , &! shaded canopy absorption for solar radiation
        thermk   (n_pft) , &! canopy gap fraction for tir radiation
        extkb    (n_pft) , &! (k, g(mu)/mu) direct solar extinction coefficient
        extkd    (n_pft)    ! diffuse and scattered diffuse PAR extinction coefficient

#if (defined DGVM)
  real(r8),INTENT(inout) :: &
        t10min    (n_pft) , &! annual minimum of 10-day running mean (K)
        lai_ind   (n_pft) , &! LAI per individual
        dphen     (n_pft) , &! phenology [0 to 1]
        leafon    (n_pft) , &! leafon days
        leafof    (n_pft) , &! leafoff days
        firelength(n_pft) , &! fire season in days
        litterag (n_pft)  , &! above ground litter
        litterbg (n_pft)  , &! below ground litter
        cpool_fast(n_pft) , &! fast carbon pool
        cpool_slow(n_pft) , &! slow carbon pool
        k_fast_ave(n_pft) , &! decomposition rate
        k_slow_ave(n_pft) , &! decomposition rate
        litter_decom_ave(n_pft) , &! decomposition rate
        fmicr     (n_pft) , &! microbial respiration (umol CO2 /m**2 /s)
        nind      (n_pft) , &! number of individuals (#/m**2)
        lm_ind    (n_pft) , &! individual leaf mass
        sm_ind    (n_pft) , &! individual sapwood mass
        hm_ind    (n_pft) , &! individual heartwood mass
        rm_ind    (n_pft) , &! individual root mass
        tmomin20  (n_pft) , &! 20-yr running mean of tmomin
        agdd0     (n_pft) , &! growing dgree days above 0
        agdd      (n_pft) , &! growing dgree days above 5
        agddtw    (n_pft) , &! growing dgree days above twmax
        agdd20    (n_pft) , &! 20-yr running mean of agdd
        t_mo      (n_pft) , &! 30-day mean temperature of 2m (K)
        t_mo_sum  (n_pft) , &! 30-day accumulated temperature of 2m (K)
        t_mo_min  (n_pft) , &! annual min of t_mo (Kelvin)
        crownarea (n_pft) , &! area that each individual tree takes up (m^2)
        htop      (n_pft) , &! canopy top (m)
        tsai      (n_pft) , &! one-sided stem area index, no burying by snow
        fpcgrid   (n_pft) , &! foliar projective cover on gridcell (fraction)
        bm_inc    (n_pft) , &! biomass increment
        afmicr    (n_pft) , &! annual microbial respiration
        annpsn    (n_pft) , &! annual photosynthesis (umol CO2 /m**2)
        annpsnpot (n_pft) , &! annual potential photosynthesis (same units)
        tref10    (n_pft) , &! 10-day averaged temperature at 2m
        tref_sum  (n_pft) , &! sum of temperature in current day
        t10    (10,n_pft) , &! array to record 10 day temperature
        assimn10  (n_pft) , &! 10-day averaged assimilation rate
        assimn_sum(n_pft) , &! sum of assimn of current day
        an10   (10,n_pft) , &! aarry to record 10 day assimn
        anngpp    (n_pft) , &! annual gpp
        annfrmf   (n_pft) , &! annual frmf
        annfrms   (n_pft) , &! annual frms
        annfrmr   (n_pft) , &! annual frmr
        annfrg    (n_pft) , &! annual frg
        cflux_litter_soil(n_pft)  , &! 
        cflux_litter_atmos(n_pft) , &! 
        nday              , &! counting the model days
        nyr               , &! counting the model days
        prec365           , &! yearly running mean of precipitation(mm/s)
        ifpre     (n_pft)    ! whether PFT present in patch. 1. for present

! Maximum fractional cover of vegetation [-]   ! added by zhq. feb.10.09
  real(r8), dimension(21), parameter :: &
  vegc=(/1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, &
         1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, &
         1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0 /)
#endif

#if(defined DyN)
  real(r8),INTENT(inout) :: &
         litter_leaf(n_pft)   , &! leaf-derived litter for PFT on modelled area basis (gC/m2)
         litter_wood(n_pft)   , &! heart&sapwood-derived litter for PFT on modelled area basis(gC/m2)
         litter_root(n_pft)   , &! fine root-derived litter for PFT on modelled area basis (gC/m2)
         litter_repr(n_pft)   , &! litter derived from allocation to reproduction for PFT on modelled area (gC/m2)

         litter_leaf_n(n_pft) , &! leaf-derived N litter for PFT on modelled area basis (gN/m2)
         litter_wood_n(n_pft) , &! heart&sapwood-derived N litter for PFT on modelled area basis(gN/m2)
         litter_root_n(n_pft) , &! fine root-derived N litter for PFT on modelled area basis (gN/m2)
         litter_repr_n(n_pft) , &! litter derived from allocation to reproduction N for PFT on modelled

         afcton_leaf  (n_pft) , &! annual floating leaf C:N ratio
         afcton_root  (n_pft) , &! annual floating root C:N ratio
         afcton_sap   (n_pft) , &! annual floating sapwood C:N ratio
         an_up        (n_pft) , &! annual plant nitrogen uptake gN/m2 veget'd area for each pft
         an_stress    (n_pft) , &! annual nitrogen stress for production

         soil_no3             , &! soil NO3 pool(gN/m2)
         soil_no2             , &! soil NO2 pool(gN/m2)
         soil_no              , &! soil NO  pool(gN/m2)
         soil_n2o             , &! soil N2O pool(gN/m2)
         soil_n2              , &! soil N2  pool(gN/m2)
         soil_nh4                ! soil NH4 pool(gN/m2)
#endif       

! Fluxes
! ----------------------------------------------------------------------
  real(r8), INTENT(out) :: &
        taux     (n_pft), &! wind stress: E-W [kg/m/s**2]
        tauy     (n_pft), &! wind stress: N-S [kg/m/s**2]
        fsena    (n_pft), &! sensible heat from canopy height to atmosphere [W/m2]
        fevpa    (n_pft), &! evapotranspiration from canopy height to atmosphere [mm/s]
        lfevpa   (n_pft), &! latent heat flux from canopy height to atmosphere [W/2]
        fsenl    (n_pft), &! ensible heat from leaves [W/m2]
        fevpl    (n_pft), &! evaporation+transpiration from leaves [mm/s]
        etr      (n_pft), &! transpiration rate [mm/s]
        fseng    (n_pft), &! sensible heat flux from ground [W/m2]
        fevpg    (n_pft), &! evaporation heat flux from ground [mm/s]
        olrg     (n_pft), &! outgoing long-wave radiation from ground+canopy
        fgrnd    (n_pft), &! ground heat flux [W/m2]
        xerr            , &! water balance error at current time-step [mm/s]
        zerr     (n_pft), &! energy balnce errore at current time-step [W/m2]
        sm              , &! rate of snowmelt [kg/(m2 s)]
        qsubl    (n_pft), &! sublimation rate from snow pack (mm h2o /s) [+]

        tref     (n_pft), &! 2 m height air temperature [K]
        qref     (n_pft), &! 2 m height air specific humidity
        trad     (n_pft), &! radiative temperature [K]
        rsur            , &! surface runoff (mm h2o/s)
        rnof            , &! total runoff (mm h2o/s)
       
        rst      (n_pft), &! canopy stomatal resistance 
        assim    (n_pft), &! canopy assimilation
        respc    (n_pft), &! canopy respiration

        parsun   (n_pft), &! PAR by sunlit leaves [W/m2]
        parsha   (n_pft), &! PAR by shaded leaves [W/m2]
        sabvsun  (n_pft), &! solar absorbed by sunlit vegetation [W/m2]
        sabvsha  (n_pft), &! solar absorbed by shaded vegetation [W/m2]
        sabg     (n_pft), &! solar absorbed by ground  [W/m2]
        sabvg    (n_pft), &! solar absorbed by ground + vegetation [W/m2]

        emis     (n_pft), &! averaged bulk surface emissivity
        z0ma     (n_pft), &! effective roughness [m]
        zol      (n_pft), &! dimensionless height (z/L) used in Monin-Obukhov theory
        rib      (n_pft), &! bulk Richardson number in surface layer
        ustar    (n_pft), &! u* in similarity theory [m/s]
        qstar    (n_pft), &! q* in similarity theory [kg/kg]
        tstar    (n_pft), &! t* in similarity theory [K]
        u10m     (n_pft), &! 10m u-velocity
        v10m     (n_pft), &! 10m v-velocity
        f10m     (n_pft), &! integral of profile function for momentum at 10m
        fm       (n_pft), &! integral of profile function for momentum
        fh       (n_pft), &! integral of profile function for heat
        fq       (n_pft)   ! integral of profile function for moisture

! ----------------------- Local  Variables -----------------------------
   real(r8) :: &
        calday            , &! Julian cal day (1.xx to 365.xx)
        endwb             , &! water mass at the end of time step
        errore     (n_pft), &! energy balnce errore (Wm-2)
        errorw            , &! water balnce errore (mm)
        fiold      (maxsnl+1:nl_soil), &! fraction of ice relative to the total water
        orb_coszen        , &! cosine of the solar zenith angle
        pg         (n_pft), &! water onto ground including canopy runoff [kg/(m2 s)]
        pg_rain    (n_pft), &! liquid water onto ground [kg/(m2 s)]
        pg_snow    (n_pft), &! ice onto ground [kg/(m2 s)]
        qseva      (n_pft), &! ground surface evaporation rate (mm h2o/s)
        qsdew      (n_pft), &! ground surface dew formation (mm h2o /s) [+]
        qfros      (n_pft), &! surface dew added to snow pack (mm h2o /s) [+]
        rootr      (1:nl_soil,n_pft), &! root resistance of a layer(n_pft), all layers add to 1.0
        scvold            , &! snow cover for previous time step [mm]
        ssw               , &! water volumetric content of soil surface layer [m3/m3]
        tssub      (7)    , &! surface/sub-surface temperatures [K]
        tssea             , &! sea surface temperature [K]
        totwb             , &! water mass at the begining of time step
        wt         (n_pft), &! fraction of vegetation buried (covered) by snow [-]
        zi(maxsnl:nl_soil), &! interface level below a "z" level (m)
        ldew_col          , &! column average of canopy water [mm/s], add by zhq. 07/03/2009
        fevpa_col         , &! column average of evapotranspiration from canopy height to atmosphere [mm/s]
        z0m        (n_pft), &! aerodynamic roughness length [m] add by zhq. 07/27/2010
        displa     (n_pft)   ! displacement height [m] add by zhq. 07/27/2010

  integer :: &
        ivt       (n_pft), &! add by zhq. 06/02/2009
        filterp   (n_pft), &! array stores filtered index, add by zhq. 07/29/2009
        num_filterp         ! actual number of patches in a column after filtered.
  integer snl            , &! number of snow layers
        imelt(maxsnl+1:nl_soil), &! flag for: melting=1, freezing=2, Nothing happended=0
        lb               , &! lower bound of arrays
        j                , &! do looping index
        c,p,fp              ! indices of different sturcture add by zhq. 06/06/2009
#if (defined DGVM)
  real(r8),parameter:: T0=273.16         ! freezing temperature[k]
  real(r8),parameter:: T5=278.16         ! 5 degree above freezing temperature
  real(r8) twmax(n_pft)                  ! upper limit of temperature of the warmest month
  real(r8) dnpp                          ! NPP per step (gC/m2 patch area/step) 
  real(r8) dgpp                          ! GPP per step (gC/m2 patch area/step)
  real(r8) dfrmf                         ! frmf per step (gC/m2 patch area/step)
  real(r8) dfrms                         ! frms per step (gC/m2 patch area/step)
  real(r8) dfrmr                         ! frmr per step (gC/m2 patch area/step)
  real(r8) dfrg                          ! frg per step (gC/m2 patch area/step)
  real(r8) frmf(n_pft)                   ! leaf maintenance respiration  (mol CO2 /m**2 /s)
  real(r8) tl(n_pft)                     ! leaf temperature[k]
  real(r8) rstfac(n_pft)                 ! factor of soil water stress 
  real(r8) wf                            ! soil water as frac. of whc for top 0.5 m ! only for DGVM
  real(r8) tsoi25                        ! soil temperature to 0.25 m (Kelvin)
  real(r8) l_cton, r_cton, s_cton        ! C:N ratio for leaf, root, sapwood
  real(r8) soil_doc                      ! especially for DyN
#endif

#ifdef MYBUG
if(mybug) write(6,*) 'wt_patch', wt_patch
#endif

!jidy: turning off the DGVM to temporarily run LAI forcing cases
#if(defined VEGDATA)
  z0m = z0m_r
  displa = displa_r
#else
#if(defined DGVM)
  if(itypwat==0) then
  ! zhq: z0m and displa can not be 0. 07/27/2010
    z0m = max(z0m_r * htop,0.01)
    displa = max(displa_r * htop,0.1)
  ! zhq: avoid incredible restart. 07/29/2010
    wt_patch = fpcgrid * wt_column
  else
    z0m = z0m_r
    displa = displa_r
  endif
#else
  z0m = z0m_r
  displa = displa_r
#endif
#endif
!print*,'z0m',z0m
!print*,'displa',displa
  pg = 0.
  pg_snow = 0.
  pg_rain = 0.

      ivt = nint(ivt_r)

  ! Build filter of current patches. add by zhq. 07/30/2009 

      call BuildFilter(n_pft, itypwat, wt_patch, num_filterp, filterp)

!======================================================================
!  [1] Solar absorbed by vegetation and ground
!======================================================================
      do fp = 1, num_filterp
          p = filterp(fp)

      call netsolar (itypwat,sigf(p),alb(:,:,p),ssun(:,:,p),ssha(:,:,p),&
                     sols,soll,solsd,solld,&
                     parsun(p),parsha(p),sabvsun(p),sabvsha(p),sabg(p),sabvg(p))
      end do
!======================================================================

if(itypwat<=3)then     ! <=== not lake and ocean (itypwat = 0, 1, 2, 3)

!======================================================================
                          !initial set
      scvold = scv        !snow mass at previous time step

      snl = 0
      do j=maxsnl+1,0
         if(wliq(j)+wice(j)>0.) snl=snl-1
      enddo

      zi(0)=0.
      if(snl<0)then
      do j = -1, snl, -1
         zi(j)=zi(j+1)-dz(j+1)
      enddo
      endif
      do j = 1,nl_soil
         zi(j)=zi(j-1)+dz(j)
      enddo

      ! pft-level canopy water averaged to column

      call p2c(n_pft,num_filterp,filterp(:),wt_patch(:),wt_column,ldew(:),ldew_col) 

      totwb = ldew_col + scv + sum(wice(1:)+wliq(1:))
      if(snl<0) fiold(snl+1:0)=wice(snl+1:0)/(wliq(snl+1:0)+wice(snl+1:0))

!----------------------------------------------------------------------
! [2] Canopy interception and precipitation onto ground surface
!----------------------------------------------------------------------
      do fp = 1, num_filterp
         p = filterp(fp)
         call leafinterception (dtime,dewmx,chil(p), &
                                prc,prl,tm,scv,sigf(p),lai(p),sai(p),ldew(p),pg(p))
      enddo

!----------------------------------------------------------------------
! [3] Initialize new snow nodes for snowfall / sleet
!----------------------------------------------------------------------
      call newsnow (itypwat,maxsnl,n_pft,num_filterp,filterp,dtime,wt_patch,wt_column,tm,tg,pg,tcrit, &
                    zi(:0),z(:0),dz(:0),tss(:0),wliq(:0),wice(:0),fiold(:0),&
                    snl,sag,scv,snowdp,pg_rain,pg_snow)
!----------------------------------------------------------------------
! [4] Energy AND Water balance 
!----------------------------------------------------------------------
!print*,'here-1'
      lb  = snl + 1           ! lower bound of array 
      CALL THERMAL  &
           ( itypwat ,lb        ,nl_soil  ,num_filterp,filterp,&
             wt_patch,wt_column ,dtime    ,n_pft    ,trsmx0   ,&
             zlnd    ,zsno      ,csoilc   ,dewmx    ,capr     ,&
             cnfac   ,csol      ,porsl    ,phi0     ,bsw      ,&
             dkmg    ,dkdry     ,dksatu   ,lai      ,sai      ,&
             z0m     ,displa    ,sqrtdi   ,rootfr(:,:),effcon ,&
             vmax25  ,slti      ,hlti     ,shti     ,hhti     ,&
             trda    ,trdm      ,trop     ,gradm    ,binter   ,&
             extkn   ,hu        ,ht       ,hq       ,us       ,&
             vs      ,tm        ,qm       ,rhoair   ,psrf     ,&
             pco2m   ,po2m      ,coszen   ,parsun   ,parsha   ,&
             sabvsun ,sabvsha   ,sabg     ,frl      ,extkb    ,&
             extkd   ,thermk    ,fsno     ,sigf     ,dz(lb:)  ,&
             z(lb:)  ,zi(lb-1:) ,tlsun    ,tlsha    ,tss(lb:) ,&
             wice(lb:),wliq(lb:),ldew     ,scv      ,snowdp   ,&
             imelt(lb:),taux    ,tauy     ,fsena    ,fevpa    ,&
             lfevpa  ,fsenl     ,fevpl    ,etr      ,fseng    ,&
             fevpg   ,olrg      ,fgrnd    ,rootr(:,:),qseva   ,&
             qsdew   ,qsubl     ,qfros    ,sm       ,tref     ,&
             qref    ,trad      ,rst      ,assim    ,respc    ,&
             errore  ,emis      ,z0ma     ,zol      ,rib      ,&
             ustar   ,qstar     ,tstar    ,u10m     ,v10m     ,&
#if (defined DGVM)
             annpsn  ,annpsnpot ,tl       ,rstfac   ,&    
#endif
             f10m    ,fm        ,fh       ,fq       ,ivt     ) 
      
!print*,'here-2'
       CALL WATER  (   itypwat      ,lb         ,nl_soil  ,dtime     ,&
             n_pft    ,num_filterp  ,filterp    ,wt_patch ,wt_column ,&
             z(lb:)   ,dz(lb:)      ,zi(lb-1:)  ,bsw      ,porsl     ,&
             phi0     ,hksati       ,rootr(:,:) ,tss(lb:) ,wliq(lb:) ,&
             wice(lb:),pg_rain      ,sm         ,etr      ,qseva     ,&
             qsdew    ,qsubl        ,qfros      ,rsur     ,rnof      ,&
             wtfact   ,pondmx       ,ssi        ,wimp     ,smpmin   )

!print*,'here-3'
      if(snl<0)then
         ! Compaction rate for snow 
         ! Natural compaction and metamorphosis. The compaction rate
         ! is recalculated for every new timestep
         lb  = snl + 1   ! lower bound of array 
         call snowcompaction (lb,dtime,&
                         imelt(lb:0),fiold(lb:0),tss(lb:0),&
                         wliq(lb:0),wice(lb:0),dz(lb:0))

         ! Combine thin snow elements
         lb = maxsnl + 1
         call snowlayerscombine (lb,snl,&
                         z(lb:1),dz(lb:1),zi(lb-1:1),&
                         wliq(lb:1),wice(lb:1),tss(lb:1),scv,snowdp)

         ! Divide thick snow elements
         if(snl<0) &
         call snowlayersdivide (lb,snl,&
                         z(lb:0),dz(lb:0),zi(lb-1:0),&
                         wliq(lb:0),wice(lb:0),tss(lb:0))
      endif

      ! Set zero to the empty node
      if(snl==0) sag=0.

      if(snl>maxsnl)then
         wice(maxsnl+1:snl)=0.
         wliq(maxsnl+1:snl)=0.
         tss (maxsnl+1:snl)=0.
         z   (maxsnl+1:snl)=0.
         dz  (maxsnl+1:snl)=0.
      endif

      lb = snl + 1
      tg = tss(lb)
      ssw = min(1.,1.e-3*wliq(1)/dz(1))

      ! ----------------------------------------
      ! Update the snow age 
      ! ----------------------------------------
      call snowage (dtime, tg, scv, scvold, sag)

      ! ----------------------------------------
      ! energy balance
      ! ----------------------------------------
      do fp = 1, num_filterp
         p = filterp(fp)
         zerr(p)=errore(p)
#ifdef MYBUG
         if(abs(errore(p))>.2) then
            write(6,*) 'Warning: energy balance violation ',errore(p),ivt(p)
         end if
#endif
      end do
      ! ----------------------------------------
      ! water balance
      ! ----------------------------------------
      call p2c(n_pft,num_filterp,filterp,wt_patch,wt_column,fevpa,fevpa_col) 
      call p2c(n_pft,num_filterp,filterp,wt_patch,wt_column,ldew,ldew_col) 

      endwb=sum(wice(1:)+wliq(1:))+ldew_col+scv
      errorw=(endwb-totwb)-(prc+prl-fevpa_col-rnof)*dtime

      if(itypwat>1) errorw=0.      !wetland, glacier

      xerr=errorw/dtime

#ifdef MYBUG
      if(abs(errorw)>1.e-3) then
         write(6,*) 'Warning: water balance violation', errorw, itypwat
      end if
#endif

!======================================================================

else if(itypwat <= 5)then   ! <=== is lake (itypwat = 4 or 5) 

!======================================================================

      scvold = scv          ! snow mass at previous time step
      zi(0)=0.
      do j = 1,nl_soil
      zi(j)=zi(j-1)+dz(j)
      enddo

      do fp = 1, num_filterp
          p = filterp(fp)

         CALL LAKE (nl_soil,itypwat,dlat,dtime,&
           z(1:),dz(1:),zi(0:),hu,ht,hq,us,vs,tm,qm,prc,prl,rhoair,psrf,&
           sabg(p),frl,tg,tss(1:),wliq(1:),wice(1:),scv,snowdp,trad(p),tref(p),qref(p),&
           taux(p),tauy(p),fsena(p),fevpa(p),lfevpa(p),fseng(p),fevpg(p),olrg(p),fgrnd(p),tcrit,&
           emis(p),z0ma(p),zol(p),rib(p),ustar(p),qstar(p),tstar(p),u10m(p),v10m(p),f10m(p),fm(p),fh(p),fq(p))
   
         ! null data for lake component
           tlsun(p) = tm
           tlsha(p) = tm
           ldew(p)  = 0.0

           fsenl(p) = 0.0
           fevpl(p) = 0.0
           etr(p)   = 0.0
           rst(p)   = -9999.
           assim(p) = 0.0
           respc(p) = 0.0

           zerr(p)=0.
      end do 

           snl = 0
           z   (:0) = 0.0
           dz  (:0) = 0.0
           tss (:0) = 0.0
           wliq(:) = 0.0
           wice(:) = 0.0
           rsur  = 0.0
           rnof  = 0.0
           xerr=0.

           ! Update the snow age 
           call snowage (dtime, tg, scv, scvold, sag)
           ssw = 0.0

!======================================================================

else                     ! <=== is ocean (itypwat >= 99) 

!======================================================================
! simple ocean-sea ice model

    tssea = tg
    tssub (1:7) = tss (1:7) 

    do fp = 1, num_filterp
        p = filterp(fp)

    CALL SOCEAN (dosst,dtime,oro,hu,ht,hq,&
                 us,vs,tm,qm,rhoair,psrf,sabg,frl,tssea,tssub(1:7),scv,&
                 taux(p),tauy(p),fsena(p),fevpa(p),lfevpa(p),fseng(p),fevpg(p),tref(p),qref(p),&
                 z0ma(p),zol(p),rib(p),ustar(p),qstar(p),tstar(p),u10m(p),v10m(p),f10m(p),fm(p),&
                 fh(p),fq(p),emis(p),olrg(p))
                 
               ! null data for sea component
                 tlsun  (p) = tm
                 tlsha  (p) = tm
                 ldew   (p) = 0.0

                 trad   (p) = tssea
                 fsenl  (p) = 0.0
                 fevpl  (p) = 0.0  
                 etr    (p) = 0.0  
                 fgrnd  (p) = 0.0
                 rst    (p) = -9999.
                 assim  (p) = 0.0
                 respc  (p) = 0.0

                 zerr   (p) = 0.0
    end do
                 z    (:) = 0.0
                 dz   (:) = 0.0
                 tss  (:) = 0.0; tss(1:7) = tssub(1:7)
                 wliq (:) = 0.0
                 wice (:) = 0.0
                 tg      = tssea
                 sag     = 0.0
                 snowdp  = scv/1000.*20.
                 rsur    = 0.0
                 rnof    = 0.0
                 xerr    = 0.0

!======================================================================

endif

!======================================================================
! Preparation for the next time step
! 1) time-varying parameters for vegatation
! 2) fraction of snow cover 
! 3) solar zenith angle and
! 4) albedos 
!======================================================================

  ! cosine of solar zenith angle 
    calday = float(jday) + float(msec)/86400.
    coszen = orb_coszen(calday,dlon,dlat)

if(itypwat == 0)then   ! column of natural vegetation

#if (defined DGVM)
     
    ! calculate water fraction of top 50cm soil layer
      call cal_wf(nl_soil,porsl,wliq(1:),wice(1:),z(1:),dz(1:),wf)
     
    ! calculate soil temperature of top 25cm soil layer
      call cal_tsoi25(nl_soil,z(1:),dz(1:),tss(1:),tsoi25)

    ! calculate the yearly averaged precipitation rate
    ! including the convective and large scale precipitation
      call cal_prec365(prec365,prc,prl)

#if (defined DyN)
    !*PLACEHOLDER (wait further implementation)
    ! nitrate leaching
    ! call nleaching(dtime,nl_soil,porsl,wice(1:),dz(1:),rnof,soil_no3)
    ! zhq: Accumulate dissolved soil carbon per step on vegt'd area 
    ! soil_doc = 0.
#endif
    do p = 1, numpft 
    ! calculate the fire season around a year
      call Fireseason(dtime, ivt(p), tref(p), litterag(p), &
                      pftpar(6,p),firelength(p),wf)

    ! litter and soil decomposition; returns litterag, litterbg and fmicr
#if(!defined DyN)
      call LitterSOM(dtime,nyr,wf,tsoi25,litterag(p),litterbg(p),&
                     cpool_fast(p),cpool_slow(p),k_fast_ave(p),&
                     k_slow_ave(p),litter_decom_ave(p),fmicr(p),afmicr(p),&
                     cflux_litter_soil(p),cflux_litter_atmos(p))
#else
    !*PLACEHOLDER (wait further implementation)
    ! call DGVMSomDynam(dtime, nyr, wf, tsoi25, litterag(p), litterbg(p), &
    !                cton_soil(p), afcton_leaf(p), afcton_sap(p), afcton_root(p), &
    !                litter_leaf(p), litter_wood(p), litter_root(p), litter_repr(p), &
    !                litter_leaf_n(p), litter_wood_n(p), litter_root_n(p), litter_repr_n(p),&
    !                soil_doc, soil_nh4, cpool_fast(p), cpool_slow(p), &
    !                k_fast_ave(p), k_slow_ave(p), litter_decom_ave(p), fmicr(p), afmicr(p))
#endif

    end do

#if(defined DyN)
  !*PLACEHOLDER (wait further implementation)
  ! call DGVMntransform(dtime,wf,tsoi25,soil_nh4,soil_no3,&
  !                     soil_no2,soil_no,soil_n2o,soil_n2,soil_doc)
#endif

    do fp = 1, num_filterp
        p = filterp(fp)
    ! Calculate plant respiration 
      dnpp = 0.
#if(defined DyN)
      respc(p) = 0.
      l_cton = afcton_leaf(p)
      s_cton = afcton_sap(p)
      r_cton = afcton_root(p)
#else
      respc(p) = 0.
      l_cton = pftpar(13,p)
      s_cton = pftpar(14,p)
      r_cton = pftpar(15,p)
#endif
    ! print*,'clmmain-npp1',assim(p),respc(p),dnpp,fmicr(p)
      call DGVMRespiration(dtime, ivt(p), fpcgrid(p), nind(p), dphen(p), lm_ind(p), sm_ind(p),&
                           rm_ind(p), pftpar(5,p), l_cton, s_cton, r_cton,&
                           tl(p), tsoi25, assim(p), respc(p), frmf(p), dnpp, &
                           dgpp, dfrmf, dfrms, dfrmr, dfrg) 

!print*,'clmmain-npp2',p,assim(p),respc(p),dnpp
#if(defined DyN)
    !*PLACEHOLDER (wait further implementation)
    ! Calculate N uptake by plant
    ! call n_uptake_by_plant(ivt(p),dtime,soil_no3,soil_nh4,wf,tsoi25,cton_pro(p),dnpp,&
    !                        an_up(p),an_stress(p),assim(p),respc(p))
    ! print*,p,'N4', soil_no3, soil_no2, soil_no, soil_n2o, soil_n2,soil_nh4
    ! patch.fluxes.dcflux_veg+=(indiv.resp-indiv.assim)*indiv.n_stress;
#endif

    ! bm_inc=[gC/m2 vegt'd area] 

      bm_inc(p)  = bm_inc(p)  + dnpp

      anngpp(p)  = anngpp(p)  + dgpp
      annfrmf(p) = annfrmf(p) + dfrmf
      annfrms(p) = annfrms(p) + dfrms
      annfrmr(p) = annfrmr(p) + dfrmr
      annfrg(p)  = annfrg(p)  + dfrg

    ! calculate the 10-day averaged temperature at 2m
      call cal_T10(msec,tref(p),tref10(p),tref_sum(p),t10(1:,p),nday)

#ifdef MYBUG
      if (msec.eq.0) print *, 'tref10', p, tref10(p)
#endif
 
    ! calculate the 30-day averaged temperature at 2m
      call cal_Tmomin(msec,tref(p),t_mo(p),t_mo_sum(p),t_mo_min(p),nday)

    ! calculate the 10-day averaged net assimilation rate
      call cal_AN10(msec,assim(p),frmf(p),assimn_sum(p),assimn10(p),an10(1:,p),nday)

      if (msec .eq. 0) then    ! daily time step
       ! growing degree days above 0
         call cal_GDD(agdd0(p),tref10(p),T0)

       ! growing degree days above 5              
         call cal_GDD(agdd(p),tref10(p),T5)

       ! growing degree days above twmax(PFT parameter)
         twmax(p)= pftpar(31,p) + T0           !convert to K
         call cal_GDD(agddtw(p),tref10(p),twmax(p))

       ! need to update lai and sai, fveg, green, they are done once in a day only
         call DGVMphenology(ivt(p),lai_ind(p),agdd0(p),tref10(p),assimn10(p),tmomin20(p),&
                            vegclass(p),raingreen(p),summergreen(p),t10min(p),&
                            leafon(p),leafof(p),dphen(p),htop(p),lai(p),sai(p),rstfac(p))

       ! print*,'lai1',ivt(p),dphen(p),lai(p),rstfac(p)
      endif ! end daily time step

    ! update fveg and green, for temporary usage. added by zhq. feb.11.09

      fveg(p) = vegc(ivt(p))
      if(fveg(p) > 0.)then
         green(p) = 1.
      else
         green(p) = 0.
      endif
    enddo ! end pft loop of dgvm 

  ! print*,'N3', soil_no3, soil_nh4

  ! update nday on natural vegetation column 
    if (msec .eq. 0) nday = nday +1    

#endif
  endif  ! end natural vegetation
  
  if(itypwat<=5)then    !land grid
#if(defined VEGDATA)
     lai(:) = 0.
     sai(:) = 0.
#endif
     do fp = 1, num_filterp
         p = filterp(fp)
#if(defined EcoDynamics)
        if(dolai)then
           call lai_empirical(ivt(p),nl_soil,rootfr(:,p),tss(1:),lai(p),sai(p),fveg(p),green(p))
        endif
#endif
#if(defined VEGDATA)
      ! lai(p)=lai_r(p); sai(p)=sai_r(p); green(p)=green_r(p); fveg(p)=fveg_r(p)
        lai(p)=lai_r(p); sai(p)=sai_r(p);
        fveg(p) = vegc(ivt(p))
        if(fveg(p) > 0.)then
           green(p) = 1.
        else
           green(p) = 0.
        endif
#endif

      ! fraction of snow cover.
      ! call snowfraction (fveg(p),z0m(p),snowdp,wt(p),sigf(p),fsno)
        call snowfraction (itypwat,fveg(p),z0m(p),snowdp,scv,wt(p),sigf(p),fsno)
      ! albedos 
      ! we supposed call it every time-step, just because 
      ! other vegeation related parameters are needed to create
        if(doalb)then
           call albland (itypwat,albsol,chil(p),ref(:,:,p),tran(:,:,p),&
                         fveg(p),green(p),lai(p),sai(p),coszen,wt(p),fsno,scv,sag,ssw,tg,&
                         alb(:,:,p),albg(:,:,p),albv(:,:,p),ssun(:,:,p),ssha(:,:,p),thermk(p),extkb(p),extkd(p))
        endif
     enddo ! pft loop

  else ! ocean grid

     do fp = 1, num_filterp
         p = filterp(fp)

      !*if(doalb)then
        call albocean (oro,scv,coszen,alb(:,:,p))
      !*endif
 
      ! null data for sea component
        lai(p) = 0.0
        sai(p) = 0.0
        green(p) = 0.0
        fveg(p) = 0.0
        sigf(p) = 0.0
 
        albg(:,:,p) = alb(:,:,p)
        albv(:,:,p) = 0.0
        ssun(:,:,p) = 0.0
        ssha(:,:,p) = 0.0
        thermk(p) = 0.0
        extkb(p) = 0.0
        extkd(p) = 0.0
     enddo                 ! pft loop

     fsno = 0.0
  endif

!----------------------------------------------------------------------

 end subroutine CLMMAIN
