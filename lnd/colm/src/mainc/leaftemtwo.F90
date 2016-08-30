
#include <define.h>

 subroutine leaftemtwo (dtime ,csoilc  ,dewmx  ,htvp   ,lai    ,&
            sai     ,displa   ,sqrtdi  ,z0m    ,effcon ,vmax25 ,&
            slti    ,hlti     ,shti    ,hhti   ,trda   ,trdm   ,&
            trop    ,gradm    ,binter  ,extkn  ,extkb  ,extkd  ,&
            hu      ,ht       ,hq      ,us     ,vs     ,thm    ,&
            th      ,thv      ,qm      ,psrf   ,rhoair ,parsun ,&
            parsha  ,sabvsun  ,sabvsha ,frl    ,fsun   ,thermk ,&
            rstfac  ,po2m     ,pco2m   ,sigf   ,etrc   ,tg     ,&
            qg      ,dqgdT    ,emg     ,tlsun  ,tlsha  ,ldew   ,&
            taux    ,tauy     ,fseng   ,fevpg  ,cgrnd  ,cgrndl ,&
            cgrnds  ,tref     ,qref    ,rst    ,assim  ,respc  ,&
            fsenl   ,fevpl    ,etr     ,dlrad  ,ulrad  ,z0ma   ,&
            zol     ,rib      ,ustar   ,qstar  ,tstar  ,f10m   ,&
#if (defined DGVM)
            annpsn     ,annpsnpot,&
#endif
            fm      ,fh       ,fq      ,ivt) 
          

!=======================================================================
! Original author : Yongjiu Dai, August 15, 2001
!
! Foliage energy conservation is given by foliage energy budget equation
!                   sunlit:  [Rnet - Hf - LEf] = 0
!                   shaded:  [Rnet - Hf - LEf] = 0
! The equation is solved by Newton-Raphson iteration, in which this iteration
! includes the calculation of the photosynthesis and stomatal resistance, and the
! integration of turbulent flux profiles. The sensible and latent heat
! transfer between foliage and atmosphere and ground is linked by the equations:
!                   Ha = [Hf]sunlit + [Hf]shaded + Hg 
!                   Ea = [Ef]sunlit + [Ef]shaded + Eg
!
!=======================================================================

  use precision
  use phycon_module, only : vonkar, grav, hvap, cpair, stefnc
  implicit none
 
!-----------------------Arguments---------------------------------------

  integer , INTENT(in) :: &
        ivt          ! land cover type
!        nl_soil      ! upper bound of array  ! deleted by zhq

  real(r8), INTENT(in) :: &
        dtime,      &! time step [second]
        csoilc,     &! drag coefficient for soil under canopy [-]
        dewmx,      &! maximum dew
!        tss (1:nl_soil),&! soil temperature [K]    !deleted by zhq                     
        htvp         ! latent heat of evaporation (/sublimation) [J/kg]

! vegetation parameters
  real(r8), INTENT(in) :: &
        lai,        &! adjusted leaf area index for seasonal variation [-]
        sai,        &! stem area index  [-]
        displa,     &! displacement height [m]
        sqrtdi,     &! inverse sqrt of leaf dimension [m**-0.5]
        z0m,        &! roughness length, momentum [m]
!        rootr(1:nl_soil),&!effective fraction of roots in each soil layer !deleted by zhq

        effcon,     &! quantum efficiency of RuBP regeneration (mol CO2 / mol quanta)
        vmax25,     &! maximum carboxylation rate at 25 C at canopy top
                     ! the range : 30.e-6 <-> 100.e-6 (mol co2 m-2 s-1)
        shti,       &! slope of high temperature inhibition function     (s1)
        hhti,       &! 1/2 point of high temperature inhibition function (s2)
        slti,       &! slope of low temperature inhibition function      (s3)
        hlti,       &! 1/2 point of low temperature inhibition function  (s4)
        trda,       &! temperature coefficient in gs-a model             (s5)
        trdm,       &! temperature coefficient in gs-a model             (s6)
        trop,       &! temperature coefficient in gs-a model             (273+25)
        gradm,      &! conductance-photosynthesis slope parameter
        binter,     &! conductance-photosynthesis intercept
        extkn        ! coefficient of leaf nitrogen allocation

! input variables
  real(r8), INTENT(in) :: &
        hu,         &! observational height of wind [m]
        ht,         &! observational height of temperature [m]
        hq,         &! observational height of humidity [m]
        us,         &! wind component in eastward direction [m/s]
        vs,         &! wind component in northward direction [m/s]
        thm,        &! intermediate variable (tm+0.0098*ht)
        th,         &! potential temperature (kelvin)
        thv,        &! virtual potential temperature (kelvin)
        qm,         &! specific humidity at reference height [kg/kg]
        psrf,       &! pressure at reference height [pa]
        rhoair,     &! density air [kg/m**3]

        parsun,     &! par absorbed per unit sunlit lai [w/m**2]
        parsha,     &! par absorbed per unit shaded lai [w/m**2]
        sabvsun,    &! solar radiation absorbed by vegetation [W/m2]
        sabvsha,    &! solar radiation absorbed by vegetation [W/m2]
        frl,        &! atmospheric infrared (longwave) radiation [W/m2]

        fsun,       &! sunlit fraction of canopy
        extkb,      &! (k, g(mu)/mu) direct solar extinction coefficient
        extkd,      &! diffuse and scattered diffuse PAR extinction coefficient
        thermk,     &! canopy gap fraction for tir radiation
        rstfac,     &! factor of soil water stress to transpiration  

        po2m,       &! atmospheric partial pressure  o2 (pa)
        pco2m,      &! atmospheric partial pressure co2 (pa)

        sigf,       &! fraction of veg cover, excluding snow-covered veg [-]
        etrc,       &! maximum possible transpiration rate (mm/s)
        tg,         &! ground surface temperature [K]
        qg,         &! specific humidity at ground surface [kg/kg]
        dqgdT,      &! temperature derivative of "qg"
        emg          ! vegetation emissivity

  real(r8), INTENT(inout) :: &
        tlsun,      &! sunlit leaf temperature [K]
        tlsha,      &! shaded leaf temperature [K]
        ldew,       &! depth of water on foliage [mm]
        taux,       &! wind stress: E-W [kg/m/s**2]
        tauy,       &! wind stress: N-S [kg/m/s**2]
        fseng,      &! sensible heat flux from ground [W/m2]
        fevpg,      &! evaporation heat flux from ground [mm/s]
        cgrnd,      &! deriv. of soil energy flux wrt to soil temp [w/m2/k]
        cgrndl,     &! deriv, of soil sensible heat flux wrt soil temp [w/m2/k]
        cgrnds,     &! deriv of soil latent heat flux wrt soil temp [w/m**2/k]
        tref,       &! 2 m height air temperature (kelvin)
#if (defined DGVM)
        annpsn,     &! annual photosynthesis (umol CO2 /m**2)
        annpsnpot,  &! annual potential photosynthesis (umol CO2/ m**2)
#endif
        qref         ! 2 m height air specific humidity


  real(r8), INTENT(out) :: &
        rst,        &! total stomatal resistance (s m-1)
        assim,      &! total assimilation (mol m-2 s-1)
        respc,      &! total respiration (mol m-2 s-1)
        fsenl,      &! sensible heat from leaves [W/m2]
        fevpl,      &! evaporation+transpiration from leaves [mm/s]
        etr,        &! transpiration rate [mm/s]
        dlrad,      &! downward longwave radiation blow the canopy [W/m2]
        ulrad,      &! upward longwave radiation above the canopy [W/m2]

        z0ma,       &! effective roughness [m]
        zol,        &! dimensionless height (z/L) used in Monin-Obukhov theory
        rib,        &! bulk Richardson number in surface layer
        ustar,      &! friction velocity [m/s]
        tstar,      &! temperature scaling parameter
        qstar,      &! moisture scaling parameter
        f10m,       &! integral of profile function for momentum at 10m
        fm,         &! integral of profile function for momentum
        fh,         &! integral of profile function for heat
        fq           ! integral of profile function for moisture

!-----------------------Local Variables---------------------------------

! assign iteration parameters
   integer,  parameter :: itmax  = 40   ! maximum number of iteration
   integer,  parameter :: itmin  = 6    ! minimum number of iteration
   real(r8), parameter :: delmax = 1.0  ! maximum change in leaf temperature [K]
   real(r8), parameter :: dtmin  = 0.01 ! max limit for temperature convergence [K]
   real(r8), parameter :: dlemin = 0.1  ! max limit for energy flux convergence [w/m2]

   real(r8) dtlsun(0:itmax+1)     ! difference of tlsun between two iterative step
   real(r8) dtlsha(0:itmax+1)     ! difference of tlsha between two iterative step

   real(r8) :: &
        zldis,      &! reference height "minus" zero displacement heght [m]
        zii,        &! convective boundary layer height [m]
        z0mv,       &! roughness length, momentum [m]
        z0hv,       &! roughness length, sensible heat [m]
        z0qv,       &! roughness length, latent heat [m]
        zeta,       &! dimensionless height used in Monin-Obukhov theory
        beta,       &! coefficient of conective velocity [-]
        wc,         &! convective velocity [m/s]
        wc2,        &! wc**2
        dth,        &! diff of virtual temp. between ref. height and surface 
        dthv,       &! diff of vir. poten. temp. between ref. height and surface
        dqh,        &! diff of humidity between ref. height and surface
        obu,        &! monin-obukhov length (m)
        um,         &! wind speed including the stablity effect [m/s]
        ur,         &! wind speed at reference height [m/s]
        uaf,        &! velocity of air within foliage [m/s]
        temp1,      &! relation for potential temperature profile
        temp2,      &! relation for specific humidity profile
        temp12m,    &! relation for temperature at 2m
        temp22m,    &! relation for specific humidity at 2m
        thvstar,    &! virtual potential temperature scaling parameter
        taf,        &! air temperature within canopy space [K]
        qaf,        &! humidity of canopy air [kg/kg]
        eah,        &! canopy air vapor pressure (pa)
        pco2a,      &! canopy air co2 pressure (pa)

        fdry,       &! fraction of foliage that is green and dry [-]
        fwet,       &! fraction of foliage covered by water [-]
        cf,         &! heat transfer coefficient from leaves [-]
        rb,         &! leaf boundary layer resistance [s/m]
        rbsun,      &! bulk boundary layer resistance of sunlit fraction of canopy
        rbsha,      &! bulk boundary layer resistance of shaded fraction of canopy
        rd,         &! aerodynamical resistance between ground and canopy air
        ram,        &! aerodynamical resistance [s/m]
        rah,        &! thermal resistance [s/m]
        raw,        &! moisture resistance [s/m]
        cah,        &! heat conduactance for air [m/s]
        cgh,        &! heat conduactance for ground [m/s]
        cfsunh,     &! heat conduactance for sunlit leaf [m/s]
        cfshah,     &! heat conduactance for shaded leaf [m/s]
        caw,        &! latent heat conduactance for air [m/s]
        cgw,        &! latent heat conduactance for ground [m/s]
        cfsunw,     &! latent heat conduactance for sunlit leaf [m/s]
        cfshaw,     &! latent heat conduactance for shaded leaf [m/s]
        wtshi,      &! sensible heat resistance for air, grd and leaf [-]
        wtsqi,      &! latent heat resistance for air, grd and leaf [-]
        wta0,       &! normalized heat conduactance for air [-]
        wtg0,       &! normalized heat conduactance for ground [-]
        wtlsun0,    &! normalized heat conductance for air and sunlit leaf [-]
        wtlsha0,    &! normalized heat conductance for air and shaded leaf [-]
        wtaq0,      &! normalized latent heat conduactance for air [-]
        wtgq0,      &! normalized heat conduactance for ground [-]
        wtlsunq0,   &! normalized latent heat cond. for air and sunlit leaf [-]
        wtlshaq0,   &! normalized latent heat cond. for air and shaded leaf [-]

        del,        &! absolute change in leaf temp in current iteration [K]
        del2,       &! change in leaf temperature in previous iteration [K]
        dele,       &! change in heat fluxes from leaf [K]
        dele2,      &! change in heat fluxes from leaf [K]
        det,        &! maximum leaf temp. change in two consecutive iter [K]
 
        obuold,     &! monin-obukhov length from previous iteration
        tlsunbef,   &! sunlit leaf temperature from previous iteration [K]
        tlshabef,   &! shaded leaf temperature from previous iteration [K]
        ecidif,     &! excess energies [W/m2]
        err,        &! balance error

        fsha,       &! shaded fraction of canopy
        laisun,     &! sunlit leaf area index, one-sided
        laisha,     &! shaded leaf area index, one-sided
        rssun,      &! sunlit leaf stomatal resistance [s/m]
        rssha,      &! shaded leaf stomatal resistance [s/m]
        assimsun,   &! sunlit leaf assimilation rate [mol co2 /m**2/ s] [+]
        assimsha,   &! shaded leaf assimilation rate [mol co2 /m**2/ s] [+]
#if (defined DGVM)
        assimpot,   &! total potential assimilation (mol m-2 s-1), added by zhq. 
        assimsun_pot, &! potential sunlit leaf assimilation rate, added by zhq.
        assimsha_pot, &! potential shaded leaf assimilation rate, added by zhq. 
#endif
        respcsun,   &! sunlit leaf respiration rate [mol co2 /m**2/ s] [+]
        respcsha,   &! shaded leaf respiration rate [mol co2 /m**2/ s] [+]
        rsoil,      &! soil respiration
        gah2o,      &! 
        tprcor       !

   integer  it, nmozsgn
   real(r8) delta1, delta2, fac, fac1, fac2 
   real(r8) cintsun(3), cintsha(3)
   real(r8) etrsun, etrsha, evplwetsun, evplwetsha, evplwet
   real(r8) etrsun_dtlsun, etrsha_dtlsun, evplwetsun_dtlsun, evplwetsha_dtlsun
   real(r8) etrsun_dtlsha, etrsha_dtlsha, evplwetsun_dtlsha, evplwetsha_dtlsha

   real(r8) irabsun, senlsun, evplsun, ftsun
   real(r8) irabsha, senlsha, evplsha, ftsha
   real(r8) dirabsun_dtlsun, senlsun_dtlsun, dirabsun_dtlsha, senlsun_dtlsha  
   real(r8) dirabsha_dtlsun, senlsha_dtlsun, dirabsha_dtlsha, senlsha_dtlsha
   real(r8) evplsun_dtlsun, evplsun_dtlsha, dftsunDTlsun, dftsunDTlsha
   real(r8) evplsha_dtlsun, evplsha_dtlsha, dftshaDTlsun, dftshaDTlsha     

   real(r8) eisun, eisha, deisundT, deishadT
   real(r8) qsatlsun, qsatlsha, qsatlsunDT, qsatlshaDT
   real(r8) z0mg, w, csoilcn
   real(r8) clai, elwmax, elwdif
   real(r8) taf0 ! zhq for test @09/02/2010

!-----------------------------------------------------------------------
! initialization of errors and  iteration parameters
!-----------------------------------------------------------------------
 
       it     = 1    ! counter for leaf temperature iteration
       del    = 0.0  ! change in leaf temperature from previous iteration
       dele   = 0.0  ! latent head flux from leaf for previous iteration
 
       dtlsun(0) = 0.
       dtlsha(0) = 0.

!-----------------------------------------------------------------------
! leaf area index
!-----------------------------------------------------------------------
! partion visible canopy absorption to sunlit and shaded fractions
! to get average absorbed par for sunlit and shaded leaves
       fsha = 1. - fsun

       laisun = lai*fsun
       laisha = lai*fsha

! scaling-up coefficients from leaf to canopy
       cintsun(1) = (1.-exp(-(extkn+extkb)*lai))/(extkn+extkb)
       cintsun(2) = (1.-exp(-(extkb+extkd)*lai))/(extkb+extkd)
       cintsun(3) = (1.-exp(-extkb*lai))/extkb

       cintsha(1) = (1.-exp(-extkn*lai))/extkn - cintsun(1)
       cintsha(2) = (1.-exp(-extkd*lai))/extkd - cintsun(2)
       cintsha(3) = lai - cintsun(3)

!-----------------------------------------------------------------------
! get fraction of wet and dry canopy surface (fwet & fdry)
! initial saturated vapor pressure and humidity and their derivation
!-----------------------------------------------------------------------

      !clai = 4.2 * 1000. * 0.2
       clai = 0.0
       call dewfraction (sigf,lai,sai,dewmx,ldew,fwet,fdry)

       call qsadv(tlsun,psrf,eisun,deisunDT,qsatlsun,qsatlsunDT)
       call qsadv(tlsha,psrf,eisha,deishaDT,qsatlsha,qsatlshaDT)

!-----------------------------------------------------------------------
! initial for fluxes profile
!-----------------------------------------------------------------------
                          
       assimsun = 0.     ! add by zhq  
       assimsha = 0.     ! add by zhq
#if (defined DGVM)
       assimsun_pot = 0. ! add by zhq
       assimsha_pot = 0. ! add by zhq
       assimpot = 0.
#endif
       respcsun = 0.     ! add by zhq
       respcsha = 0.     ! add by zhq

       nmozsgn = 0    ! number of times moz changes sign
       obuold = 0.    ! monin-obukhov length from previous iteration
       zii = 1000.    ! m  (pbl height)
       beta = 1.      ! -  (in computing W_*)

       z0mv = z0m; z0hv = z0m; z0qv = z0m

       taf = 0.5 * (tg + thm)
       qaf = 0.5 * (qm + qg)
       taf0 = taf     ! add by zhq for convergence test

       pco2a = pco2m
       tprcor = 44.6*273.16*psrf/1.013e5
!       rsoil = 0. !soil respiration (mol m-2 s-1)
!      rsoil = 1.22e-6*exp(308.56*(1./56.02-1./(tg-227.13))) 
!      rsoil = rstfac * 0.23 * 15. * 2.**((tg-273.16-10.)/10.) * 1.e-6
!      rsoil = 5.22 * 1.e-6
       rsoil = 0.22 * 1.e-6

       ur = max(0.1, sqrt(us*us+vs*vs))    ! limit set to 0.1
       dth = thm - taf
       dqh = qm - qaf
       dthv = dth*(1.+0.61*qm) + 0.61*th*dqh
       zldis = hu - displa

    !* DGVM may crash the model here
       if(zldis.lt.0.) then
          write(6,*), 'negative zldis', hu, displa, zldis
          call abort
       end if

       call moninobukini(ur,th,thm,thv,dth,dqh,dthv,zldis,z0mv,um,obu)

! ======================================================================
!     BEGIN stability iteration 
! ======================================================================

      do while (it .le. itmax) 

         tlsunbef = tlsun
         tlshabef = tlsha

         del2 = del
         dele2 = dele
 
!-----------------------------------------------------------------------
! Aerodynamical resistances
!-----------------------------------------------------------------------
! Evaluate stability-dependent variables using moz from prior iteration
        call moninobuk(hu,ht,hq,displa,z0mv,z0hv,z0qv,obu,um,&
                       ustar,temp1,temp2,temp12m,temp22m,f10m,fm,fh,fq)
 
! Aerodynamic resistance
        ram = 1./(ustar*ustar/um) 
        rah = 1./(temp1*ustar)
        raw = 1./(temp2*ustar)
 
! Bulk boundary layer resistance of leaves
        uaf = ustar
        cf = 0.01*sqrtdi/sqrt(uaf)
        rb = 0.5/(cf*uaf) 

!<->    rd = 1./(csoilc*uaf)          ! legacy from BATS
! modified by Xubin Zeng's suggestion at 08-07-2002
        z0mg = 0.01
        w = exp(-(lai+sai))
        csoilcn = (vonkar/(0.13*(z0mg*uaf/1.5e-5)**0.45))*w + csoilc*(1.-w)
        rd = 1./(csoilcn*uaf)

!-----------------------------------------------------------------------
! stomatal resistances for sunlit and shaded fractions of canopy.
! should do each iteration to account for differences in eah, leaf temperatures.
!-----------------------------------------------------------------------
        if(lai .gt. 0.001) then
           rbsun = rb / laisun
           rbsha = rb / laisha

           eah = qaf * psrf / ( 0.622 + 0.378 * qaf )    ! pa

           !-------- Sunlit leaves

           CALL stomata (   vmax25 ,effcon  ,slti     ,hlti     ,shti   ,&
               hhti  ,trda ,trdm   ,trop    ,gradm    ,binter   ,thm    ,&
               psrf  ,po2m ,pco2m  ,pco2a   ,eah      ,eisun    ,tlsun  ,parsun ,&
               rbsun ,raw  ,rstfac ,cintsun ,assimsun ,respcsun ,rssun  ,& 
#if (defined DGVM)                                   
               assimsun_pot,&                             ! added by zhq.
#endif                                                 
               ivt )

           !-------- Shaded leaves

           CALL stomata (   vmax25 ,effcon  ,slti     ,hlti     ,shti   ,&
               hhti  ,trda ,trdm   ,trop    ,gradm    ,binter   ,thm    ,&
               psrf  ,po2m ,pco2m  ,pco2a   ,eah      ,eisha    ,tlsha  ,parsha ,&
               rbsha ,raw  ,rstfac ,cintsha ,assimsha ,respcsha ,rssha  ,&
#if (defined DGVM)                                   
               assimsha_pot,&                              ! added by zhq.  
#endif
               ivt )

        else
           rssun = 2.0e4 ; rssha = 2.0e4
           assimsun = 0. ; assimsha = 0.
           respcsun = 0. ; respcsha = 0.
#if (defined DGVM)
           assimsun_pot =0. ; assimsha_pot = 0.      
#endif
        endif

! above stomatal resistances are for the canopy, the stomatal resistance 
! the "rb" in the following calculations are the averages for single leaf. thus,
        rssun = rssun * laisun
        rssha = rssha * laisha

!-----------------------------------------------------------------------
! dimensional and non-dimensional sensible and latent heat conductances
! for canopy and soil flux calculations.
!-----------------------------------------------------------------------
        delta1 = 0.0
        delta2 = 0.0
        if( (fwet .lt. .99) .AND. (sigf*etrc .gt. 1.e-12) )then
           if(qsatlsun-qaf .gt. 0.) delta1 = 1.0
           if(qsatlsha-qaf .gt. 0.) delta2 = 1.0
        endif
 
        cah = sigf / rah
        cgh = sigf / rd
        cfsunh = sigf * laisun / rb
        cfshah = sigf * (laisha + sai) / rb

        caw = sigf / raw
        cgw = sigf / rd
        cfsunw = sigf * ( (1.-delta1*(1.-fwet)) * laisun / rb &
	                           + (1. - fwet) * delta1 * laisun / (rb + rssun) )
        cfshaw = sigf * ( (1.-delta2*(1.-fwet)) * (laisha + sai) / rb &
				   + (1. - fwet) * delta2 * laisha / (rb + rssha) )

        wtshi = 1. / ( cah + cgh + cfsunh + cfshah )
        wtsqi = 1. / ( caw + cgw + cfsunw + cfshaw )

        wta0 = cah * wtshi
        wtg0 = cgh * wtshi
        wtlsun0 = cfsunh * wtshi
        wtlsha0 = cfshah * wtshi

        wtaq0 = caw * wtsqi
        wtgq0 = cgw * wtsqi
        wtlsunq0 = cfsunw * wtsqi
        wtlshaq0 = cfshaw * wtsqi
 
!-----------------------------------------------------------------------
! IR radiation, sensible and latent heat fluxes and their derivatives
!-----------------------------------------------------------------------
! the partial derivatives of areodynamical resistance are ignored 
! which cannot be determined analtically
        fac = sigf * (1. - thermk)
        fac1 = fac * fsun
        fac2 = fac * fsha

! longwave absorption and their derivatives 
        irabsun = (frl - 2. * stefnc * tlsun**4 + emg*stefnc*tg**4 ) * fac1 
        dirabsun_dtlsun = - 8.* stefnc * tlsun**3                    * fac1
        dirabsun_dtlsha = 0.

        irabsha = (frl - 2. * stefnc * tlsha**4 + emg*stefnc*tg**4 ) * fac2 
        dirabsha_dtlsha = - 8.* stefnc * tlsha**3                    * fac2 
        dirabsha_dtlsun = 0.

! sensible heat fluxes and their derivatives
        senlsun = rhoair * cpair * cfsunh &
                * ( (wta0 + wtg0 + wtlsha0)*tlsun &
                   - wta0*thm - wtg0*tg - wtlsha0*tlsha )
        senlsun_dtlsun = rhoair * cpair * cfsunh * (wta0 + wtg0 + wtlsha0)
        senlsun_dtlsha = rhoair * cpair * cfsunh * ( - wtlsha0 )

        senlsha = rhoair * cpair * cfshah &
                * ( (wta0 + wtg0 + wtlsun0)*tlsha &
                   - wta0*thm - wtg0*tg - wtlsun0*tlsun )
        senlsha_dtlsun = rhoair * cpair * cfshah * ( - wtlsun0 ) 
        senlsha_dtlsha = rhoair * cpair * cfshah * ( wta0 + wtg0 + wtlsun0 )

! latent heat fluxes and their derivatives
        etrsun = sigf * rhoair * (1.-fwet) * delta1 * laisun / (rb + rssun) &
	              * ( (wtaq0 + wtgq0 + wtlshaq0)*qsatlsun &
                      - wtaq0*qm - wtgq0*qg - wtlshaq0*qsatlsha )
        etrsun_dtlsun = sigf * rhoair * (1.-fwet) * delta1 * laisun / (rb + rssun) &
		      * (wtaq0 + wtgq0 + wtlshaq0)*qsatlsunDT 
        etrsun_dtlsha = sigf * rhoair * (1.-fwet) * delta1 * laisun / (rb + rssun) &
		      * ( - wtlshaq0*qsatlshaDT )

        if(etrsun .ge.  sigf*etrc*laisun/lai)then
           etrsun =  sigf*etrc*laisun/lai
           etrsun_dtlsun = 0.
           etrsun_dtlsha = 0.
        end if

        evplwetsun = sigf * rhoair * (1.-delta1*(1.-fwet)) * laisun / rb &
		          * ( (wtaq0 + wtgq0 + wtlshaq0)*qsatlsun &
                          - wtaq0*qm - wtgq0*qg - wtlshaq0*qsatlsha )
        evplwetsun_dtlsun = sigf * rhoair * (1.-delta1*(1.-fwet)) * laisun / rb &
			  * (wtaq0 + wtgq0 + wtlshaq0)*qsatlsunDT 
        evplwetsun_dtlsha = sigf * rhoair * (1.-delta1*(1.-fwet)) * laisun / rb &
			  * ( - wtlshaq0*qsatlshaDT )

        if(evplwetsun .ge. ldew/dtime*laisun/lai)then
           evplwetsun = ldew/dtime*laisun/lai
           evplwetsun_dtlsun = 0.
           evplwetsun_dtlsha = 0.
        endif

        evplsun = etrsun + evplwetsun
        evplsun_dtlsun = etrsun_dtlsun + evplwetsun_dtlsun
        evplsun_dtlsha = etrsun_dtlsha + evplwetsun_dtlsha

        ! ---------------

        etrsha = sigf * rhoair * (1.-fwet) * delta2 * laisha / (rb + rssha) &
	       * ( (wtaq0 + wtgq0 + wtlsunq0)*qsatlsha &
               - wtaq0*qm - wtgq0*qg - wtlsunq0*qsatlsun )
        etrsha_dtlsun = sigf * rhoair * (1.-fwet) * delta2 * laisha / (rb + rssha) &
		      * ( - wtlsunq0*qsatlsunDT )
        etrsha_dtlsha = sigf * rhoair * (1.-fwet) * delta2 * laisha / (rb + rssha) &
		      * ( (wtaq0 + wtgq0 + wtlsunq0)*qsatlshaDT )

        if(etrsha .ge. sigf*etrc*laisha/lai)then
           etrsha = sigf*etrc*laisha/lai
           etrsha_dtlsun = 0.
           etrsha_dtlsha = 0.
        endif

        evplwetsha = sigf * rhoair * (1.-delta2*(1.-fwet)) * (laisha+sai) / (rb) &
		   * ( (wtaq0 + wtgq0 + wtlsunq0)*qsatlsha &
                      - wtaq0*qm - wtgq0*qg - wtlsunq0*qsatlsun )
        evplwetsha_dtlsun = sigf * rhoair * (1.-delta2*(1.-fwet)) * (laisha+sai) / rb &
			  * ( - wtlsunq0*qsatlsunDT )
        evplwetsha_dtlsha = sigf * rhoair * (1.-delta2*(1.-fwet)) * (laisha+sai) / rb &
			  * ( (wtaq0 + wtgq0 + wtlsunq0)*qsatlshaDT )

        if(evplwetsha .ge. ldew/dtime*(laisha+sai)/(lai+sai))then
           evplwetsha = ldew/dtime*(laisha+sai)/(lai+sai) 
           evplwetsha_dtlsun = 0.
           evplwetsha_dtlsha = 0
        endif

        evplsha = etrsha + evplwetsha

        evplsha_dtlsun = etrsha_dtlsun + evplwetsha_dtlsun
        evplsha_dtlsha = etrsha_dtlsha + evplwetsha_dtlsha

! functions and their derivatives with respect to temperatures
        ftsun = sabvsun + irabsun - senlsun - hvap*evplsun 
        ftsha = sabvsha + irabsha - senlsha - hvap*evplsha

        dftsunDTlsun = dirabsun_dtlsun - senlsun_dtlsun - hvap*evplsun_dtlsun
        dftsunDTlsha = dirabsun_dtlsha - senlsun_dtlsha - hvap*evplsun_dtlsha

        dftshaDTlsun = dirabsha_dtlsun - senlsha_dtlsun - hvap*evplsha_dtlsun
        dftshaDTlsha = dirabsha_dtlsha - senlsha_dtlsha - hvap*evplsha_dtlsha

!-----------------------------------------------------------------------
! difference of temperatures by quasi-newton-raphson method for the non-linear system equations
!-----------------------------------------------------------------------
        dtlsun(it) = - (ftsun * (dftshaDTlsha-(laisha+sai)*clai/dtime) - ftsha * dftsunDTlsha) &
               / ((dftsunDTlsun-laisun*clai/dtime) * (dftshaDTlsha-(laisha+sai)*clai/dtime) &
		  - dftsunDTlsha * dftshaDTlsun)

        dtlsha(it) = - (ftsun * dftshaDTlsun - ftsha * (dftsunDTlsun-laisun*clai/dtime)) &
               / (dftsunDTlsha * dftshaDTlsun &
               - (dftsunDTlsun-laisun*clai/dtime) * (dftshaDTlsha-(laisha+sai)*clai/dtime))

        if(it .lt. itmax)then

      ! put brakes on large temperature excursions
        if(abs(dtlsun(it)).gt.delmax) dtlsun(it) = delmax*dtlsun(it)/abs(dtlsun(it))
        if(abs(dtlsha(it)).gt.delmax) dtlsha(it) = delmax*dtlsha(it)/abs(dtlsha(it))

        if(it.ge.2)then
        if(dtlsun(it-1)*dtlsun(it) .lt. 0.) dtlsun(it) = 0.5*(dtlsun(it-1) + dtlsun(it))
        if(dtlsha(it-1)*dtlsha(it) .lt. 0.) dtlsha(it) = 0.5*(dtlsha(it-1) + dtlsha(it))
        endif

        endif

        tlsun = tlsunbef + dtlsun(it)
        tlsha = tlshabef + dtlsha(it)

!       test add by zhq @ 09/03/2010
        if(abs(dtlsun(it))>10.)then 
          tlsun = taf0
          dtlsun(it) = 10. * dtlsun(it)/abs(dtlsun(it))
        endif
        if(abs(dtlsha(it))>10.)then 
          tlsha = taf0
          dtlsha(it) = 10. * dtlsha(it)/abs(dtlsha(it))
        endif
!       end of test 

!-----------------------------------------------------------------------
! square roots differences of temperatures and fluxes for use as the condition of convergences
!-----------------------------------------------------------------------
        del  = sqrt( dtlsun(it)*dtlsun(it) + dtlsha(it)*dtlsha(it) )
        dele = (    dirabsun_dtlsun * dtlsun(it) +     dirabsun_dtlsha * dtlsha(it))**2 &
             + (    dirabsha_dtlsun * dtlsun(it) +     dirabsha_dtlsha * dtlsha(it))**2 &
             + (     senlsun_dtlsun * dtlsun(it) +      senlsun_dtlsha * dtlsha(it))**2 &
             + (     senlsha_dtlsun * dtlsun(it) +      senlsha_dtlsha * dtlsha(it))**2 &
             + ( hvap*evplsun_dtlsun * dtlsun(it) +  hvap*evplsun_dtlsha * dtlsha(it))**2 &
             + ( hvap*evplsha_dtlsun * dtlsun(it) +  hvap*evplsha_dtlsha * dtlsha(it))**2  
        dele = sqrt(dele)

!-----------------------------------------------------------------------
!  saturated vapor pressures and canopy air temperature, canopy air humidity
!-----------------------------------------------------------------------
! Recalculate leaf saturated vapor pressure (ei_)for updated leaf temperature
! and adjust specific humidity (qsatl_) proportionately
        call qsadv(tlsun,psrf,eisun,deisunDT,qsatlsun,qsatlsunDT)

        call qsadv(tlsha,psrf,eisha,deishaDT,qsatlsha,qsatlshaDT)
 
! update vegetation/ground surface temperature, canopy air temperature, 
! canopy air humidity
        taf = wta0*thm + wtg0*tg + wtlsun0*tlsun +  wtlsha0*tlsha

        qaf = wtaq0*qm + wtgq0*qg + wtlsunq0*qsatlsun +  wtlshaq0*qsatlsha

! update co2 partial pressure within canopy air
        gah2o  = 1.0/raw * tprcor/thm                     ! mol m-2 s-1
        pco2a = pco2m - 1.37*psrf/max(0.446,gah2o) &
              * (assimsun + assimsha - respcsun - respcsha - rsoil)

!-----------------------------------------------------------------------
! Update monin-obukhov length and wind speed including the stability effect
!-----------------------------------------------------------------------
        dth = thm - taf       
        dqh = qm - qaf

        tstar = temp1*dth
        qstar = temp2*dqh

        thvstar = tstar + 0.61*th*qstar
        zeta = zldis*vonkar*grav*thvstar / (ustar**2*thv)

        if(zeta .ge. 0.)then                             !stable
           zeta = min(2.,max(zeta,1.e-6))
        else                                             !unstable
           zeta = max(-100.,min(zeta,-1.e-6))
        endif
        obu = zldis/zeta

        if(zeta .ge. 0.)then
          um = max(ur,0.1)
        else
          wc = (-grav*ustar*thvstar*zii/thv)**(1./3.)
          wc2 = beta*beta*wc*wc
          um = sqrt(ur*ur+wc2)
        endif

        if(obuold*obu .lt. 0.) nmozsgn = nmozsgn+1
        if(nmozsgn .ge. 4) obu = zldis/(-0.01)
        obuold = obu
 
!-----------------------------------------------------------------------
! Test for convergence
!-----------------------------------------------------------------------
      it = it+1

      if(it .gt. itmin) then
         det = max(del,del2)
         del = max(dele,dele2)
         if(det .lt. dtmin .AND. del .lt. dlemin) exit 
      endif
 
      end do                ! ITERATION

! ======================================================================
!     END stability iteration 
! ======================================================================
!if(it>40)then
! print*,ivt,'leaftemtwo->convergen',it-1,det,del,tlsun,tlsha
!endif

      z0ma = z0mv
      zol = zeta
      rib = min(5.,zol*ustar**2/(vonkar*temp1*um**2))

! canopy fluxes and total assimilation amd respiration

      if(lai .gt. 0.001) then
         rst = 1./(laisun/rssun + laisha/rssha)
      else
         rssun = 2.0e4 ; rssha = 2.0e4
         assimsun = 0. ; assimsha = 0.
#if (defined DGVM)                            
         assimsun_pot =0. ; assimsha_pot = 0. ! added by zhq. dec27,08
#endif 
         respcsun = 0. ; respcsha = 0.
         rst = 2.0e4
      endif
      assim = assimsun + assimsha
      respc = respcsun + respcsha + rsoil
#if (defined DGVM)                             
! In DGVM mode, respc represents only canopy respiration. 
! I think rsoil must be deducted here. zhq. 02/20/2010
      respc = respc - rsoil

      assimpot = assimsun_pot + assimsha_pot  !  added by zhq. dec27,08       
      annpsn = annpsn + assim * dtime
      annpsnpot = annpsnpot + assimpot * dtime
#endif

      fsenl = senlsun + senlsha + senlsun_dtlsun*dtlsun(it-1) + senlsun_dtlsha*dtlsha(it-1) &
                                + senlsha_dtlsun*dtlsun(it-1) + senlsha_dtlsha*dtlsha(it-1)
      etr = etrsun + etrsha     +  etrsun_dtlsun*dtlsun(it-1) +  etrsun_dtlsha*dtlsha(it-1) &
                                +  etrsha_dtlsun*dtlsun(it-1) +  etrsha_dtlsha*dtlsha(it-1)
      evplwet = evplwetsun + evplwetsun_dtlsun*dtlsun(it-1) + evplwetsun_dtlsha*dtlsha(it-1) &
              + evplwetsha + evplwetsha_dtlsun*dtlsun(it-1) + evplwetsha_dtlsha*dtlsha(it-1) 
      fevpl = evplsun + evplsha + evplsun_dtlsun*dtlsun(it-1) + evplsun_dtlsha*dtlsha(it-1) &
                                + evplsha_dtlsun*dtlsun(it-1) + evplsha_dtlsha*dtlsha(it-1)

      elwmax = ldew/dtime
      elwdif = max(0., evplwet-elwmax)
      evplwet = min (evplwet, elwmax)

      fevpl = fevpl - elwdif
      fsenl = fsenl + hvap*elwdif

      ! wind stresses 
      taux  = taux - sigf*rhoair*us/ram
      tauy  = tauy - sigf*rhoair*vs/ram

!-----------------------------------------------------------------------
! fluxes from ground to canopy space
!-----------------------------------------------------------------------

      fseng = fseng + cpair*rhoair*cgh*(tg-taf)
      fevpg = fevpg +    rhoair*cgw*(qg-qaf)

!-----------------------------------------------------------------------
! downward (upward) longwave radiation below (above) the canopy
!-----------------------------------------------------------------------

      dlrad = sigf * thermk * frl  &
            + stefnc * ( fac1*tlsunbef**3*(tlsunbef+4.*dtlsun(it-1)) &
		       + fac2*tlshabef**3*(tlshabef+4.*dtlsha(it-1)) )

      ulrad = stefnc * ( fac1*tlsunbef**3*(tlsunbef+4.*dtlsun(it-1)) &
		       + fac2*tlshabef**3*(tlshabef+4.*dtlsha(it-1)) &
		       + sigf*emg*thermk*tg**4 )  

!-----------------------------------------------------------------------
! Derivative of soil energy flux with respect to soil temperature (cgrnd)
!-----------------------------------------------------------------------
      cgrnds = cgrnds + cpair*rhoair*cgh*(1.-wtg0)
      cgrndl = cgrndl + rhoair*cgw*(1.-wtgq0)*dqgdT
      cgrnd  = cgrnds + cgrndl*htvp

!-----------------------------------------------------------------------
! balance check
!-----------------------------------------------------------------------
      err = sabvsun+sabvsha &
          + irabsun+dirabsun_dtlsun*dtlsun(it-1)+dirabsun_dtlsha*dtlsha(it-1) &
          + irabsha+dirabsha_dtlsun*dtlsun(it-1)+dirabsha_dtlsha*dtlsha(it-1) &
          - fsenl-hvap*fevpl
 
#ifdef MYBUG
      if(abs(err) .gt. .2) write(6,*) 'leaftemtwo.F90 : energy imbalance',err, it-1
#endif
 
!-----------------------------------------------------------------------
! Update dew accumulation (kg/m2)
!-----------------------------------------------------------------------

      ldew = max(0.,ldew-evplwet*dtime)

!-----------------------------------------------------------------------
! 2 m height air temperature
!-----------------------------------------------------------------------
      tref = tref + sigf*(thm + temp1*dth * (1./temp12m - 1./temp1))
      qref = qref + sigf*( qm + temp2*dqh * (1./temp22m - 1./temp2))

 end subroutine leaftemtwo
