
#include <define.h>

 SUBROUTINE THERMAL (itypwat  ,lb       ,nl_soil  ,num_filterp,filterp,&
                     wt_patch ,wt_column,dtime    ,n_pft    ,trsmx0   ,&
                     zlnd     ,zsno     ,csoilc   ,dewmx    ,capr     ,&
                     cnfac    ,csol     ,porsl    ,phi0     ,bsw      ,&
                     dkmg     ,dkdry    ,dksatu   ,lai      ,sai      ,&
                     z0m      ,displa   ,sqrtdi   ,rootfr   ,effcon   ,&
                     vmax25   ,slti     ,hlti     ,shti     ,hhti     ,&
                     trda     ,trdm     ,trop     ,gradm    ,binter   ,&
                     extkn    ,hu       ,ht       ,hq       ,us       ,&
                     vs       ,tm       ,qm       ,rhoair   ,psrf     ,&
                     pco2m    ,po2m     ,coszen   ,parsun   ,parsha   ,&
                     sabvsun  ,sabvsha  ,sabg     ,frl      ,extkb    ,&
                     extkd    ,thermk   ,fsno     ,sigf     ,dz       ,&
                     z        ,zi       ,tlsun    ,tlsha    ,tss      ,&
                     wice     ,wliq     ,ldew     ,scv      ,snowdp   ,&
                     imelt    ,taux     ,tauy     ,fsena    ,fevpa    ,&
                     lfevpa   ,fsenl    ,fevpl    ,etr      ,fseng    ,&
                     fevpg    ,olrg     ,fgrnd    ,rootr    ,qseva    ,&
                     qsdew    ,qsubl    ,qfros    ,sm       ,tref     ,&
                     qref     ,trad     ,rst      ,assim    ,respc    ,&
                     errore   ,emis     ,z0ma     ,zol      ,rib      ,&
                     ustar    ,qstar    ,tstar    ,u10m     ,v10m     ,&
#if (defined DGVM)
                     annpsn  ,annpsnpot ,tl       ,rstfac_r ,&
#endif
                     f10m    ,fm      ,fh     ,fq     ,ivt   ) 

!=======================================================================
! this is the main subroutine to execute the calculation 
! of thermal processes and surface fluxes
!
! Original author : Yongjiu Dai, 09/15/1999; 08/30/2002
! 
! FLOW DIAGRAM FOR THERMAL.F90
! 
! THERMAL ===> qsadv
!              groundfluxes
!              eroot                      |dewfraction
!              leaftemone |               |qsadv
!              leaftemtwo |  ---------->  |moninobukini
!                                         !moninobuk
!                                         |stomata
!
!              groundTem     ---------->   meltf
!                               
!=======================================================================

  use precision
  use phycon_module, only : denh2o,roverg,hvap,hsub,rgas,cpair,&
                            stefnc,denice,tfrz,vonkar,grav 
  implicit none
 
!---------------------Argument------------------------------------------

  integer, INTENT(in) ::   &
        lb                ,&! lower bound of array 
        nl_soil           ,&! upper bound of array
        n_pft             ,&! number of pfts in a column
        ivt(n_pft)        ,&! land cover type
        num_filterp       ,&!
        filterp(n_pft)    ,&!
        itypwat             ! land water type (0=soil, 1=urban or built-up, 2=wetland,
                            ! 3=land ice, 4=deep lake, 5=shallow lake)
  real(r8), INTENT(in) ::  &
        wt_patch(n_pft),   &! weight of pft relative to grid
        wt_column,         &! weight of column relative to grid
        dtime,             &! model time step [second]
        trsmx0,            &!max transpiration for moist soil+100% veg.  [mm/s]
        zlnd,              &!roughness length for soil [m]
        zsno,              &!roughness length for snow [m]
        csoilc,            &!drag coefficient for soil under canopy [-]
        dewmx,             &!maximum dew
        capr,              &!tuning factor to turn first layer T into surface T
        cnfac,             &!Crank Nicholson factor between 0 and 1

        ! soil physical parameters
        csol(1:nl_soil)   ,&! heat capacity of soil solids [J/(m3 K)]
        porsl(1:nl_soil)  ,&! soil porosity [-]
        phi0(1:nl_soil)   ,&! soil water suction, negative potential [m]
        bsw(1:nl_soil)    ,&! clapp and hornbereger "b" parameter [-]
        dkmg(1:nl_soil)   ,&! thermal conductivity of soil minerals [W/m-K]
        dkdry(1:nl_soil)  ,&! thermal conductivity of dry soil [W/m-K]
        dksatu(1:nl_soil) ,&! thermal conductivity of saturated soil [W/m-K]

        ! vegetation parameters
        lai(n_pft)        ,&! adjusted leaf area index for seasonal variation [-]
        sai(n_pft)        ,&! stem area index  [-]
        z0m(n_pft)        ,&! roughness length(n_pft), momentum [m]
        displa(n_pft)     ,&! displacement height [m]
        sqrtdi(n_pft)     ,&! inverse sqrt of leaf dimension [m**-0.5]
        rootfr(1:nl_soil,n_pft),&! root fraction 
       
        effcon(n_pft)     ,&! quantum efficiency of RuBP regeneration (mol CO2/mol quanta)
        vmax25(n_pft)     ,&! maximum carboxylation rate at 25 C at canopy top
        slti(n_pft)       ,&! slope of low temperature inhibition function      [s3]
        hlti(n_pft)       ,&! 1/2 point of low temperature inhibition function  [s4]
        shti(n_pft)       ,&! slope of high temperature inhibition function     [s1]
        hhti(n_pft)       ,&! 1/2 point of high temperature inhibition function [s2]
        trda(n_pft)       ,&! temperature coefficient in gs-a model             [s5]
        trdm(n_pft)       ,&! temperature coefficient in gs-a model             [s6]
        trop(n_pft)       ,&! temperature coefficient in gs-a model          
        gradm(n_pft)      ,&! conductance-photosynthesis slope parameter
        binter(n_pft)     ,&! conductance-photosynthesis intercept
        extkn(n_pft)      ,&! coefficient of leaf nitrogen allocation

        ! atmospherical variables and observational height
        hu,          &! observational height of wind [m]
        ht,          &! observational height of temperature [m]
        hq,          &! observational height of humidity [m]
        us,          &! wind component in eastward direction [m/s]
        vs,          &! wind component in northward direction [m/s]
        tm,          &! temperature at agcm reference height [kelvin]
        qm,          &! specific humidity at agcm reference height [kg/kg]
        rhoair,      &! density air [kg/m3]
        psrf,        &! atmosphere pressure at the surface [pa]
        pco2m,       &! CO2 concentration in atmos. (35 pa)
        po2m,        &! O2 concentration in atmos. (20900 pa)
        frl,         &! atmospheric infrared (longwave) radiation [W/m2]

        ! radiative fluxes
        coszen            ,&! cosine of the solar zenith angle
        parsun(n_pft)     ,&! photosynthetic active radiation by sunlit leaves (W m-2)
        parsha(n_pft)     ,&! photosynthetic active radiation by shaded leaves (W m-2)
        sabvsun(n_pft)    ,&! solar radiation absorbed by vegetation [W/m2]
        sabvsha(n_pft)    ,&! solar radiation absorbed by vegetation [W/m2]
        sabg(n_pft)       ,&! solar radiation absorbed by ground [W/m2]
        extkb(n_pft)      ,&! (k(n_pft), g(mu)/mu) direct solar extinction coefficient
        extkd(n_pft)      ,&! diffuse and scattered diffuse PAR extinction coefficient
        thermk(n_pft)     ,&! canopy gap fraction for tir radiation

        ! state variable (1)
        fsno              ,&! fraction of ground covered by snow
        sigf(n_pft)       ,&! fraction of veg cover, excluding snow-covered veg [-]
        dz(lb:nl_soil)    ,&! layer thickiness [m]
        z (lb:nl_soil)    ,&! node depth [m]
        zi(lb-1:nl_soil)    ! interface depth [m]

        ! state variables (2)
  real(r8), INTENT(inout) :: &
        tlsun(n_pft)      ,&! sunlit leaf temperature [K]
        tlsha(n_pft)      ,&! shaded leaf temperature [K]
        tss (lb:nl_soil)  ,&! soil temperature [K]
        wice(lb:nl_soil)  ,&! ice lens [kg/m2]
        wliq(lb:nl_soil)  ,&! liqui water [kg/m2]
        ldew(n_pft)       ,&! depth of water on foliage [kg/(m2 s)] 
        scv               ,&! snow cover, water equivalent [mm, kg/m2]
        snowdp              ! snow depth [m]
#if (defined DGVM)
  real(r8), INTENT(inout) :: &
        annpsn(n_pft)     ,&! annual photosynthesis (umol CO2 /m**2)
        annpsnpot(n_pft)  ,&! annual potential photosynthesis (umol CO2/ m**2)
        tl(n_pft)         ,&! leaf temperature
        rstfac_r(n_pft)     ! factor of soil water stress 
#endif

  integer, INTENT(out) ::  & 
       imelt(lb:nl_soil)    ! flag for melting or freezing [-]
 
        ! Output fluxes
  real(r8), INTENT(out) :: &
        taux(n_pft),       &! wind stress: E-W [kg/m/s**2]
        tauy(n_pft),       &! wind stress: N-S [kg/m/s**2]
        fsena(n_pft),      &! sensible heat from canopy height to atmosphere [W/m2]
        fevpa(n_pft),      &! evapotranspiration from canopy height to atmosphere [mm/s]
        lfevpa(n_pft),     &! latent heat flux from canopy height to atmosphere [W/m2]
        fsenl(n_pft),      &! sensible heat from leaves [W/m2]
        fevpl(n_pft),      &! evaporation+transpiration from leaves [mm/s]
        etr(n_pft),        &! transpiration rate [mm/s]
        fseng(n_pft),      &! sensible heat flux from ground [W/m2]
        fevpg(n_pft),      &! evaporation heat flux from ground [mm/s]
        olrg(n_pft),       &! outgoing long-wave radiation from ground+canopy
        fgrnd(n_pft),      &! ground heat flux [W/m2]
        errore(n_pft),     &! energy balnce error [w/m2]

        rootr(1:nl_soil,n_pft),&! root resistance of a layer(n_pft), all layers add to 1

        qseva(n_pft),      &! ground surface evaporation rate (mm h2o/s)
        qsdew(n_pft),      &! ground surface dew formation (mm h2o /s) [+]
        qsubl(n_pft),      &! sublimation rate from snow pack (mm h2o /s) [+]
        qfros(n_pft),      &! surface dew added to snow pack (mm h2o /s) [+]

        sm,                &! rate of snowmelt [kg/(m2 s)]
        tref(n_pft),       &! 2 m height air temperature [kelvin]
        qref(n_pft),       &! 2 m height air specific humidity
        trad(n_pft),       &! radiative temperature [K]

        rst(n_pft),        &! stomatal resistance (s m-1)
        assim(n_pft),      &! assimilation
        respc(n_pft),      &! respiration

       ! additional variables required by coupling with WRF or RSM model
        emis(n_pft),       &! averaged bulk surface emissivity
        z0ma(n_pft),       &! effective roughness [m]
        zol(n_pft),        &! dimensionless height (z/L) used in Monin-Obukhov theory
        rib(n_pft),        &! bulk Richardson number in surface layer
        ustar(n_pft),      &! u* in similarity theory [m/s]
        qstar(n_pft),      &! q* in similarity theory [kg/kg]
        tstar(n_pft),      &! t* in similarity theory [K]
        u10m(n_pft),       &! 10m u-velocity
        v10m(n_pft),       &! 10m v-velocity
        f10m(n_pft),       &! integral of profile function for momentum at 10m
        fm(n_pft),         &! integral of profile function for momentum
        fh(n_pft),         &! integral of profile function for heat
        fq(n_pft)           ! integral of profile function for moisture

!---------------------Local Variables-----------------------------------

  integer i,j,p,fp

  real(r8) :: &
       cgrnd(n_pft),       &! deriv. of soil energy flux wrt to soil temp [w/m2/k]
       cgrndl(n_pft),      &! deriv, of soil sensible heat flux wrt soil temp [w/m2/k]
       cgrnds(n_pft),      &! deriv of soil latent heat flux wrt soil temp [w/m**2/k]
       degdT,              &! d(eg)/dT
       dqgdT,              &! d(qg)/dT
       dlrad(n_pft),       &! downward longwave radiation blow the canopy [W/m2]
       eg,                 &! water vapor pressure at temperature T [pa]
       egsmax,             &! max. evaporation which soil can provide at one time step
       egidif(n_pft),      &! the excess of evaporation over "egsmax"
       emg,                &! ground emissivity (0.97 for snow, 
                            ! glaciers and water surface; 0.96 for soil and wetland)
       errsoi_pft(n_pft),  &! soil energy balnce error [w/m2]
       errsoi_col,         &!    
       etrc(n_pft),        &! maximum possible transpiration rate [mm/s]
       fac,                &! soil wetness of surface layer
       fact(lb:nl_soil),   &! used in computing tridiagonal matrix
       fsun,               &! fraction of sunlit canopy
       hr,                 &! relative humidity
       htvp,               &! latent heat of vapor of water (or sublimation) [j/kg]
       olru(n_pft),        &! olrg excluding dwonwelling reflection [W/m2]
       olrb(n_pft),        &! olrg assuming blackbody emission [W/m2]
       psit,               &! negative potential of soil
       par,                &! PAR absorbed by canopy [W/m2]
       qg,                 &! ground specific humidity [kg/kg]
       qsatg,              &! saturated humidity [kg/kg]
       qsatgdT,            &! d(qsatg)/dT
       qred,               &! soil surface relative humidity
       rstfac(n_pft),      &! factor of soil water stress 
       sabv(n_pft),        &! solar absorbed by canopy [W/m2]
       thm,                &! intermediate variable (tm+0.0098*ht)
       th,                 &! potential temperature (kelvin)
       thv,                &! virtual potential temperature (kelvin)
       tg,                 &! ground surface temperature [K]
       tssbef(lb:nl_soil), &! soil/snow temperature before update
       tinc,               &! temperature difference of two time step
       ur,                 &! wind speed at reference height [m/s]
       ulrad(n_pft),       &! upward longwave radiation above the canopy [W/m2]
       wice0(lb:nl_soil),  &! ice mass from previous time-step
       wliq0(lb:nl_soil),  &! liquid mass from previous time-step
       wx,                 &! patitial volume of ice and water of surface layer
       xmf                  ! total latent heat of phase change of ground water

  real(r8) :: &
       z0ma_g(n_pft),      &! 
       zol_g(n_pft),       &!
       rib_g(n_pft),       &!
       ustar_g(n_pft),     &!
       qstar_g(n_pft),     &!
       tstar_g(n_pft),     &!
       f10m_g(n_pft),      &!
       fm_g(n_pft),        &!
       fh_g(n_pft),        &!
       fq_g(n_pft)
  real(r8) :: temp1,temp2,temp12m,temp22m,um,obu,tssum
  real(r8) extkb_t

!=======================================================================
! [1] Initial set and propositional variables
!=======================================================================

      ! fluxes 
      taux   = 0.;  tauy   = 0.    
      fsena  = 0.;  fevpa  = 0.  
      lfevpa = 0.;  fsenl  = 0.    
      fevpl  = 0.;  etr    = 0.  
      fseng  = 0.;  fevpg  = 0.    
      dlrad  = 0.;  ulrad  = 0. 
      cgrnds = 0.;  cgrndl = 0.    
      cgrnd  = 0.;  tref   = 0. 
      qref   = 0.;  rst    = 2.0e4
      assim  = 0.;  respc  = 0. 

      emis   = 0.;  z0ma   = 0.
      zol    = 0.;  rib    = 0.
      ustar  = 0.;  qstar  = 0.
      tstar  = 0.;  rootr  = 0.
      rstfac = 0.;  errsoi_col = 0.

      ! temperature and water mass from previous time step
      tg = tss(lb)
      tssbef(lb:) = tss(lb:)
      wice0(lb:) = wice(lb:)
      wliq0(lb:) = wliq(lb:)

      ! emissivity
      emg = 0.96
      if(scv>0. .OR. itypwat==3) emg = 0.97

      ! latent heat, assumed that the sublimation occured only as wliq=0
      htvp = hvap
      if(wliq(lb)<=0. .AND. wice(lb)>0.) htvp = hsub

      ! potential temperatur at the reference height
      thm = tm + 0.0098*ht              ! intermediate variable equivalent to
                                        ! tm*(pgcm/psrf)**(rgas/cpair)
      th = tm*(100000./psrf)**(rgas/cpair) ! potential T
      thv = th*(1.+0.61*qm)             ! virtual potential T
      ur = max(0.1,sqrt(us*us+vs*vs))   ! limit set to 0.1

!=======================================================================
! [2] specific humidity and its derivative at ground surface
!=======================================================================

      qred = 1.
      call qsadv(tg,psrf,eg,degdT,qsatg,qsatgdT)

      if(itypwat<=1)then            ! soil ground
         wx   = (wliq(1)/denh2o + wice(1)/denice)/dz(1)
         if(porsl(1)<1.e-6)then     ! bed rock
            fac  = 0.001
         else 
            fac  = min(1.,wx/porsl(1))
            fac  = max( fac, 0.001 )
         endif

         psit = -phi0(1) * fac ** (- bsw(1) )   ! psit = max(smpmin, psit)
         hr   = exp(psit/roverg/tg)
         qred = (1. - fsno)*hr + fsno
      endif

      qg = qred*qsatg  
      dqgdT = qred*qsatgdT

      if(qsatg > qm .AND. qm > qred*qsatg)then
        qg = qm; dqgdT = 0.
      endif

!=======================================================================
! [3] Compute sensible and latent fluxes and their derivatives with respect 
!     to ground temperature using ground temperatures from previous time step.
!=======================================================================
!pft level looping, zhq. 07/15/2009
!print*,'thermal-1'

      do fp = 1, num_filterp
          p = filterp(fp)

      if(sigf(p) <= 0.999) then
         call groundfluxes (zlnd,zsno,hu,ht,hq, &
                            us,vs,tm,qm,rhoair,psrf, &
                            ur,thm,th,thv,tg,qg,dqgdT,htvp, &
                            fsno,sigf(p),cgrnd(p),cgrndl(p),cgrnds(p), &
                            taux(p),tauy(p),fsena(p),fevpa(p),fseng(p),fevpg(p),tref(p),qref(p), &
              z0ma_g(p),zol_g(p),rib_g(p),ustar_g(p),qstar_g(p),tstar_g(p),f10m_g(p),fm_g(p),fh_g(p),fq_g(p))
      end if
!print*,p,'thermal-2-1'
!=======================================================================
! [4] Canopy temperature, fluxes from the canopy
!=======================================================================

      par = parsun(p) + parsha(p)
      sabv(p) = sabvsun(p) + sabvsha(p)
      if(sigf(p) >= 0.001) then

         ! soil water strees factor on stomatal resistance
         call eroot(nl_soil,trsmx0,porsl,bsw,phi0,rootfr(:,p),&
                             dz(1:),tss(1:),wliq(1:),rootr(:,p),etrc(p),rstfac(p))

         ! fraction of sunlit and shaded leaves of canopy
         fsun = ( 1. - exp(-extkb(p)*lai(p)) ) / max( extkb(p)*lai(p), 1.e-6 )
         if(coszen<=0.0 .OR. sabv(p)<1.) fsun = 0.

         if(fsun.le.0.1)then

            fsun = 0.
           !tl(p) = tlsha(p)
            call leaftemone (dtime ,csoilc ,dewmx  ,htvp   ,lai(p)    ,&
                 sai (p) ,displa(p) ,sqrtdi(p) ,z0m (p) ,effcon(p) ,vmax25(p) ,&
                 slti(p) ,hlti  (p) ,shti  (p) ,hhti(p) ,trda  (p) ,trdm  (p) ,&
                 trop (p),gradm (p) ,binter(p) ,extkn(p),extkb (p) ,extkd (p) ,&
                 hu    ,ht    ,hq     ,us   ,vs     ,thm ,&
                 th ,thv ,qm ,psrf ,rhoair ,par ,&
                 sabv(p) ,frl ,thermk(p) ,rstfac(p) ,po2m ,pco2m  ,&
                 sigf (p) ,etrc  (p) ,tg   ,qg   ,dqgdT  ,emg  ,&
                 tlsha(p) ,ldew  (p) ,taux(p) ,tauy(p) ,fseng(p)  ,fevpg(p)  ,&
                 cgrnd(p) ,cgrndl(p) ,cgrnds(p) ,tref(p) ,qref(p) ,rst (p) ,&
                 assim(p) ,respc (p) ,fsenl(p)  ,fevpl(p)  ,etr (p) ,dlrad(p)  ,&
                 ulrad(p) ,z0ma  (p) ,zol (p) ,rib (p) ,ustar(p)  ,qstar(p)  ,&
#if (defined DGVM)
                 annpsn(p)  ,annpsnpot(p),&
#endif
                 tstar(p) ,f10m(p) ,fm(p) ,fh(p) ,fq(p) ,ivt(p) )
                   

           !tlsun(p) = tl(p)
           !tlsha(p) = tl(p)
            tlsun(p) = tlsha(p)
         else
            call leaftemtwo (dtime ,csoilc  ,dewmx  ,htvp   ,lai(p)    ,&
                 sai   (p),displa (p),sqrtdi(p),z0m  (p),effcon(p) ,vmax25(p) ,&
                 slti  (p),hlti   (p),shti  (p),hhti (p),trda (p),trdm (p),&
                 trop  (p),gradm  (p),binter(p),extkn(p),extkb(p),extkd(p),&
                 hu      ,ht       ,hq      ,us     ,vs     ,thm    ,&
                 th      ,thv      ,qm      ,psrf   ,rhoair ,parsun(p) ,&
                 parsha(p)  ,sabvsun(p)  ,sabvsha(p) ,frl    ,fsun   ,thermk(p) ,&
                 rstfac(p),po2m    ,pco2m   ,sigf(p),etrc(p),tg     ,&
                 qg      ,dqgdT    ,emg     ,tlsun(p)  ,tlsha(p)  ,ldew(p)   ,&
                 taux   (p) ,tauy    (p) ,fseng  (p) ,fevpg (p) ,cgrnd (p) ,cgrndl(p) ,&
                 cgrnds (p) ,tref    (p) ,qref   (p) ,rst   (p) ,assim (p) ,respc (p) ,&
                 fsenl  (p) ,fevpl   (p) ,etr    (p) ,dlrad (p) ,ulrad (p) ,z0ma  (p) ,&
                 zol    (p) ,rib     (p) ,ustar  (p) ,qstar (p) ,tstar (p) ,f10m  (p) ,&
#if (defined DGVM)
                 annpsn (p) ,annpsnpot(p),&
#endif
                 fm(p) ,fh(p) ,fq(p) ,ivt(p) )
                 
!print*,p,'thermal-2-4'
         endif

      endif

    ! equate canopy temperature to air over bareland.
    ! required as sigf=0 carried over to next time step
      if(sigf(p) < 0.001)then
         fsun = 1.
         tlsun(p) = tm
         tlsha(p) = tm
         ldew(p) = 0.
      endif

#if (defined DGVM)
    if(itypwat == 0)then
!     tl = (tlsha + tlsun)/2.
     tl(p) = tlsha(p)*(1.-fsun) + tlsun(p)*fsun  !revised by zhq Nov.2008
     rstfac_r(p)=rstfac(p)
    endif
#endif
     
     end do ! end pft looping zhq. 07/15/2009
!=======================================================================
! [5] Gound temperature
!=======================================================================
      call groundtem (itypwat,lb,nl_soil,dtime,n_pft,num_filterp,filterp,wt_patch,&
                      wt_column,capr,cnfac,csol,porsl,dkmg,dkdry,dksatu, &
                      sigf,dz,z,zi,tss,wice,wliq,scv,snowdp, &
                      frl,dlrad,sabg,fseng,fevpg,cgrnd,htvp,emg, &
                      imelt,sm,xmf,fact)

!print*,'thermal-3'
!=======================================================================
! [6] Correct fluxes to present soil temperature
!=======================================================================

      tg = tss(lb)
      tinc = tss(lb) - tssbef(lb)
      fseng = fseng + tinc*cgrnds 
      fevpg = fevpg + tinc*cgrndl

! calculation of evaporative potential; flux in kg m-2 s-1.  
! egidif holds the excess energy if all water is evaporated
! during the timestep.  this energy is later added to the sensible heat flux.
! a pft loop is added to determine if corrections are required for excess evaporation
! from the top of soil by zhq. 07/20/2009

      egsmax = (wice(lb)+wliq(lb)) / dtime
!--------------------------------------------------------------------
! pft loop begin
      do fp = 1, num_filterp
         p = filterp(fp)

      egidif(p) = max( 0., fevpg(p) - egsmax )
      fevpg(p) = min ( fevpg(p), egsmax )
      fseng(p) = fseng(p) + htvp*egidif(p)

! total fluxes to atmosphere
      fsena(p) = fsenl(p) + fseng(p)
      fevpa(p) = fevpl(p) + fevpg(p)
      lfevpa(p)= hvap*fevpl(p) + htvp*fevpg(p)   ! w/m2 (accouting for sublimation)
      
      qseva(p) = 0.
      qsubl(p) = 0.
      qfros(p) = 0.
      qsdew(p) = 0.

      if(fevpg(p) >= 0.)then
! not allow for sublimation in melting (melting ==> evap. ==> sublimation)
         qseva(p) = min(wliq(lb)/dtime, fevpg(p))
         qsubl(p) = fevpg(p) - qseva(p)
      else
         if(tg < tfrz)then
            qfros(p) = abs(fevpg(p))
         else
            qsdew(p) = abs(fevpg(p))
         endif
      endif

! ground heat flux
      fgrnd(p) = sabg(p) + dlrad(p) + (1.-sigf(p))*emg*frl &
            - emg*stefnc*tssbef(lb)**3*(tssbef(lb) + 4.*tinc) &
            - (fseng(p)+fevpg(p)*htvp)

! outgoing long-wave radiation from canopy + ground
      olrg(p) = ulrad(p) &
           + (1.-sigf(p))*(1.-emg)*frl &
           + (1.-sigf(p))*emg*stefnc * tssbef(lb)**4 &
! for conservation we put the increase of ground longwave to outgoing
           + 4.*emg*stefnc*tssbef(lb)**3*tinc

! averaged bulk surface emissivity 
      olrb(p) = stefnc*tssbef(lb)**3*((1.-sigf(p))*tssbef(lb) + 4.*tinc)
      olru(p) = ulrad(p) + emg*olrb(p)
      olrb(p) = ulrad(p) + olrb(p)
      emis(p) = olru(p) / olrb(p)

! radiative temperature
      trad(p) = (olrg(p)/stefnc)**0.25

! additonal variables required by WRF and RSM model
      if(sigf(p) < 0.001)then
         ustar(p) = ustar_g(p)
         tstar(p) = tstar_g(p)
         qstar(p) = qstar_g(p)
         rib(p)   = rib_g(p)
         zol(p)   = zol_g(p)
	 z0ma(p)  = z0ma_g(p)
         f10m(p)  = f10m_g(p)
         fm(p)    = fm_g(p)
         fh(p)    = fh_g(p)
         fq(p)    = fq_g(p)

      else if(sigf(p) <= 0.7)then
         z0ma(p)  = sigf(p)*z0ma(p)  + (1.-sigf(p))*z0ma_g(p)

       ! assumed um ~= ur here
         um = ur
         ustar(p) =   sqrt(max(1.e-6,sqrt(taux(p)*taux(p)+tauy(p)*tauy(p)))/rhoair)
         tstar(p) = - fsena(p)/(cpair*ustar(p)*rhoair)
         qstar(p) = - fevpa(p)/(ustar(p)*rhoair)
         zol(p) = (hu-displa(p))*vonkar*grav*(tstar(p)+0.61*th*qstar(p))/(ustar(p)**2*thv)
         if(zol(p) .ge. 0.)then   !stable
            zol(p) = min(2.,max(zol(p),1.e-6))
         else                  !unstable
            zol(p) = max(-100.,min(zol(p),-1.e-6))
         endif

         obu = (hu-displa(p))/zol(p)
         call moninobuk(hu,ht,hq,displa(p),z0ma(p),z0ma(p),z0ma(p),obu,um,&
              ustar(p),temp1,temp2,temp12m,temp22m,f10m(p),fm(p),fh(p),fq(p))
         rib(p) = min(5.,zol(p)*ustar(p)**2/(vonkar*temp1*um**2))

      else
         ustar(p) = ustar(p)
         tstar(p) = tstar(p)
         qstar(p) = qstar(p)
         rib(p)   = rib(p)
         zol(p)   = zol(p)
         z0ma(p)  = z0ma(p)
         f10m(p)  = f10m(p)
         fm(p)    = fm(p)
         fh(p)    = fh(p)
         fq(p)    = fq(p)

      endif

      u10m(p)  = us/ur * ustar(p)/vonkar * f10m(p)
      v10m(p)  = vs/ur * ustar(p)/vonkar * f10m(p)

      end do ! end pft loop

!=======================================================================
! [7] energy balance error
!=======================================================================
! Surface energy balance

      do fp = 1, num_filterp
         p = filterp(fp)
         errore(p) = sabv(p) + sabg(p) + frl - olrg(p) - fsena(p) - lfevpa(p) - fgrnd(p) 

#ifdef MYBUG
         if(abs(errore(p))>.1)then 
            write(6,*) 'THERMAL.F90 : surface flux energy  balance violation'
            write(6,100) errore(p),sabv(p),sabg(p),frl,olrg(p),fsenl(p),fseng(p),hvap*fevpl(p),htvp*fevpg(p),fgrnd(p)
         endif
#endif
      end do 
100 format(10(f15.3))

! soil energy balance
      do fp = 1, num_filterp
         p = filterp(fp)
         errsoi_pft(p) =  fgrnd(p) - xmf
         tssum = 0.
         do j = lb, nl_soil
            errsoi_pft(p) = errsoi_pft(p) - (tss(j)-tssbef(j))/fact(j)
            tssum = tssum + (tss(j)-tssbef(j))/fact(j)
         enddo
         errsoi_col = errsoi_col + errsoi_pft(p) * wt_patch(p)/wt_column
      end do

#ifdef MYBUG
      if(abs(errsoi_col)>.1)then
         write(6,*) 'THERMAL.F90 : soil energy  balance violation',itypwat
         write(6,101) errsoi_col
      endif
#endif

101 format(1(f15.3))

!      do fp = 1, num_filterp
!         p = filterp(fp)
!        errore(p) = sabv(p) + sabg(p) + frl - olrg(p) - fsena(p) - lfevpa(p) - xmf
!        do j = lb, nl_soil
!           errore(p) = errore(p) - (tss(j)-tssbef(j))/fact(j)
!        enddo
!
!       if(abs(errore(p))>.2)then
!        write(6,*) 'THERMAL.F90 : energy  balance violation'
!        write(6,100) errore(p),sabv(p),sabg(p),frl,olrg(p),fsenl(p),fseng(p),hvap*fevpl(p),htvp*fevpg(p),xmf
!       endif
!
!100    format(10(f15.3))
!      end do

 END SUBROUTINE THERMAL
