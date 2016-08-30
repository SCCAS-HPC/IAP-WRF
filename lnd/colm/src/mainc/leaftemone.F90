
#include <define.h>

 subroutine  leaftemone (dtime ,csoilc ,dewmx  ,htvp   ,lai    ,&
             sai     ,displa   ,sqrtdi ,z0m    ,effcon ,vmax25 ,&
             slti    ,hlti     ,shti   ,hhti   ,trda   ,trdm   ,&
             trop    ,gradm    ,binter ,extkn  ,extkb  ,extkd  ,&
             hu      ,ht       ,hq     ,us     ,vs     ,thm    ,&
             th      ,thv      ,qm     ,psrf   ,rhoair ,par    ,&
             sabv    ,frl      ,thermk ,rstfac ,po2m   ,pco2m  ,&
             sigf    ,etrc     ,tg     ,qg     ,dqgdT  ,emg    ,&
             tl      ,ldew     ,taux   ,tauy   ,fseng  ,fevpg  ,&
             cgrnd   ,cgrndl   ,cgrnds ,tref   ,qref   ,rst    ,&
             assim   ,respc    ,fsenl  ,fevpl  ,etr    ,dlrad  ,&
             ulrad   ,z0ma     ,zol    ,rib    ,ustar  ,qstar  ,&
#if (defined DGVM)
             annpsn  ,annpsnpot,&
#endif
             tstar   ,f10m     ,fm     ,fh     ,fq     ,ivt )
              
 
!=======================================================================
! Original author : Yongjiu Dai, August 15, 2001
!
! Foliage energy conservation is given by foliage energy budget equation
!                      Rnet - Hf - LEf = 0
! The equation is solved by Newton-Raphson iteration, in which this iteration
! includes the calculation of the photosynthesis and stomatal resistance, and the
! integration of turbulent flux profiles. The sensible and latent heat
! transfer between foliage and atmosphere and ground is linked by the equations:
!                      Ha = Hf + Hg and Ea = Ef + Eg
!
!=======================================================================

  use precision
  use phycon_module, only : vonkar, grav, hvap, cpair, stefnc
  implicit none
 
!-----------------------Arguments---------------------------------------

  integer , INTENT(in) :: &
        ivt        ! land cover type
  real(r8), INTENT(in) :: &
        dtime,      &! time step [second]
        csoilc,     &! drag coefficient for soil under canopy [-]
        dewmx,      &! maximum dew
        htvp         ! latent heat of evaporation (/sublimation) [J/kg]

! vegetation parameters
  real(r8), INTENT(in) :: &
        lai,        &! adjusted leaf area index for seasonal variation [-]
        sai,        &! stem area index  [-]
        displa,     &! displacement height [m]
        sqrtdi,     &! inverse sqrt of leaf dimension [m**-0.5]
        z0m,        &! roughness length, momentum [m]

        effcon,     &! quantum efficiency of RuBP regeneration (mol CO2 / mol quanta)
        vmax25,     &! maximum carboxylation rate at 25 C at canopy top
                     ! the range : 30.e-6 <-> 100.e-6 (mol co2 m-2 s-1)
        shti,       &! slope of high temperature inhibition function     (s1)
        hhti,       &! 1/2 point of high temperature inhibition function (s2)
        slti,       &! slope of low temperature inhibition function      (s3)
        hlti,       &! 1/2 point of low temperature inhibition function  (s4)
        trda,       &! temperature coefficient in gs-a model             (s5)
        trdm,       &! temperature coefficient in gs-a model             (s6)
        trop,       &! temperature coefficient in gs-a model         (273+25)
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

        par,        &! par absorbed per unit lai [w/m**2]
        sabv,       &! solar radiation absorbed by vegetation [W/m2]
        frl,        &! atmospheric infrared (longwave) radiation [W/m2]

        extkb,      &! (k, g(mu)/mu) direct solar extinction coefficient
        extkd,      &! diffuse and scattered diffuse PAR extinction coefficient
        thermk,     &! canopy gap fraction for tir radiation

        po2m,       &! atmospheric partial pressure  o2 (pa)
        pco2m,      &! atmospheric partial pressure co2 (pa)

        sigf,       &! fraction of veg cover, excluding snow-covered veg [-]
        etrc,       &! maximum possible transpiration rate (mm/s)
        rstfac,     &!  
        tg,         &! ground surface temperature [K]
        qg,         &! specific humidity at ground surface [kg/kg]
        dqgdT,      &! temperature derivative of "qg"
        emg          ! vegetation emissivity

  real(r8), INTENT(inout) :: &
        tl,         &! leaf temperature [K]
        ldew,       &! depth of water on foliage [mm]
        taux,       &! wind stress: E-W [kg/m/s**2]
        tauy,       &! wind stress: N-S [kg/m/s**2]
        fseng,      &! sensible heat flux from ground [W/m2]
        fevpg,      &! evaporation heat flux from ground [mm/s]
        cgrnd,      &! deriv. of soil energy flux wrt to soil temp [w/m2/k]
        cgrndl,     &! deriv, of soil sensible heat flux wrt soil temp [w/m2/k]
        cgrnds,     &! deriv of soil latent heat flux wrt soil temp [w/m**2/k]
        tref,       &! 2 m height air temperature (kelvin)
        qref         ! 2 m height air specific humidity
#if (defined DGVM)
  real(r8), INTENT(inout) :: &
        annpsn,     &! annual photosynthesis (umol CO2 /m**2)
        annpsnpot    ! annual potential photosynthesis (umol CO2/ m**2)
#endif

  real(r8), INTENT(out) :: &
        rst,        &! stomatal resistance
        assim,      &! rate of assimilation
        respc,      &! rate of respiration
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
   integer, parameter :: itmax  = 40   ! maximum number of iteration
   integer, parameter :: itmin  = 6    ! minimum number of iteration
   real(r8),parameter :: delmax = 3.0  ! maximum change in leaf temperature [K]
   real(r8),parameter :: dtmin  = 0.01 ! max limit for temperature convergence [K]
   real(r8),parameter :: dlemin = 0.1  ! max limit for energy flux convergence [w/m2]

   real(r8) dtl(0:itmax+1)     ! difference of tl between two iterative step

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
        pco2g,      &! co2 pressure (pa) at ground surface (pa)
        pco2a,      &! canopy air co2 pressure (pa)
#if (defined DGVM)
        assimpot,   &! potential canopy assimilation rate(mol/m2/s), added by zhq.             
#endif
        fdry,       &! fraction of foliage that is green and dry [-]
        fwet,       &! fraction of foliage covered by water [-]
        cf,         &! heat transfer coefficient from leaves [-]
        rb,         &! leaf boundary layer resistance [s/m]
        rbone,      &! canopy bulk boundary layer resistance 
        rd,         &! aerodynamical resistance between ground and canopy air
        ram,        &! aerodynamical resistance [s/m]
        rah,        &! thermal resistance [s/m]
        raw,        &! moisture resistance [s/m]
        clai,       &! canopy heat capacity [Jm-2K-1]
        cah,        &! heat conduactance for air [m/s]
        cgh,        &! heat conduactance for ground [m/s]
        cfh,        &! heat conduactance for leaf [m/s]
        caw,        &! latent heat conduactance for air [m/s]
        cgw,        &! latent heat conduactance for ground [m/s]
        cfw,        &! latent heat conduactance for leaf [m/s]
        wtshi,      &! sensible heat resistance for air, grd and leaf [-]
        wtsqi,      &! latent heat resistance for air, grd and leaf [-]
        wta0,       &! normalized heat conduactance for air [-]
        wtg0,       &! normalized heat conduactance for ground [-]
        wtl0,       &! normalized heat conductance for air and leaf [-]
        wtaq0,      &! normalized latent heat conduactance for air [-]
        wtgq0,      &! normalized heat conduactance for ground [-]
        wtlq0,      &! normalized latent heat cond. for air and leaf [-]

        ei,         &! vapor pressure on leaf surface [pa]
        deidT,      &! derivative of "ei" on "tl" [pa/K]
        qsatl,      &! leaf specific humidity [kg/kg]
        qsatldT,    &! derivative of "qsatl" on "tlef"

        del,        &! absolute change in leaf temp in current iteration [K]
        del2,       &! change in leaf temperature in previous iteration [K]
        dele,       &! change in heat fluxes from leaf [K]
        dele2,      &! change in heat fluxes from leaf [K]
        det,        &! maximum leaf temp. change in two consecutive iter [K]
 
        obuold,     &! monin-obukhov length from previous iteration
        tlbef,      &! leaf temperature from previous iteration [K]
        ecidif,     &! excess energies [W/m2]
        err,        &! balance error

        rs,         &! leaf stomatal resistance [s/m]
        rsoil,      &! soil respiration
        gah2o,      &! conductance between canopy and atmosphere
        gdh2o,      &! conductance between canopy and ground
        tprcor       ! tf*psur*100./1.013e5

   integer it, nmozsgn 

   real(r8) delta, fac
   real(r8) evplwet, evplwet_dtl, etr_dtl, elwmax, elwdif
   real(r8) irab, dirab_dtl, fsenl_dtl, fevpl_dtl  
   real(r8) w, csoilcn, z0mg, cint(3)

!-----------------------End Variable List-------------------------------

!print *, 'leaftemone---', ivt, lai

! initialization of errors and  iteration parameters
       it     = 1    ! counter for leaf temperature iteration
       del    = 0.0  ! change in leaf temperature from previous iteration
       dele   = 0.0  ! latent head flux from leaf for previous iteration

       dtl(0) = 0.

!-----------------------------------------------------------------------
! scaling-up coefficients from leaf to canopy
!-----------------------------------------------------------------------

!print*,'leaf-1',z0m
       cint(1) = (1.-exp(-extkn*lai))/extkn
       cint(2) = (1.-exp(-extkd*lai))/extkd
       cint(3) = lai

!-----------------------------------------------------------------------
! get fraction of wet and dry canopy surface (fwet & fdry)
! initial saturated vapor pressure and humidity and their derivation
!-----------------------------------------------------------------------

       !clai = 4.2 * 1000. * 0.2
       clai = 0.0

       call dewfraction (sigf,lai,sai,dewmx,ldew,fwet,fdry)

       call qsadv(tl,psrf,ei,deiDT,qsatl,qsatlDT)

!-----------------------------------------------------------------------
! initial for fluxes profile
!-----------------------------------------------------------------------
       assim = 0.  
#if (defined DGVM)
       assimpot = 0.         
#endif  
       nmozsgn = 0    ! number of times moz changes sign
       obuold = 0.    ! monin-obukhov length from previous iteration
       zii = 1000.    ! m  (pbl height)
       beta = 1.      ! -  (in computing W_*)

       z0mv = z0m; z0hv = z0m; z0qv = z0m

       taf = 0.5 * (tg + thm)
       qaf = 0.5 * (qm + qg)

       pco2a = pco2m
       tprcor = 44.6*273.16*psrf/1.013e5
!      rsoil = 1.22e-6*exp(308.56*(1./56.02-1./(tg-227.13)))
!      rsoil = rstfac * 0.23 * 15. * 2.**((tg-273.16-10.)/10.) * 1.e-6
!      rsoil = 5.22 * 1.e-6
!      rsoil = 0.                          !respiration (mol m-2 s-1)
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

         tlbef = tl

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

!print*,'leaf-2',uaf
      ! rd = 1./(csoilc*uaf)                 ! BATS legacy
      ! w = exp(-0.5*(lai+sai))              ! Dickinson's modification :
      ! csoilc = ( 1.-w + w*um/uaf)/rah      ! "rah" here is the resistance over
      ! rd = 1./(csoilc*uaf)                 ! bare ground fraction

! modified by Xubin Zeng's suggestion at 08-07-2002
        z0mg = 0.01
        w = exp(-(lai+sai))
        csoilcn = (vonkar/(0.13*(z0mg*uaf/1.5e-5)**0.45))*w + csoilc*(1.-w)
        rd = 1./(csoilcn*uaf)

!-----------------------------------------------------------------------
! stomatal resistances 
!-----------------------------------------------------------------------

        if(lai .gt. 0.001) then
           rbone = rb / lai
           eah = qaf * psrf / ( 0.622 + 0.378 * qaf )    ! pa

           CALL stomata (   vmax25 ,effcon ,slti  ,hlti   ,shti ,&
               hhti  ,trda ,trdm   ,trop   ,gradm ,binter ,thm  ,&
               psrf  ,po2m ,pco2m  ,pco2a  ,eah   ,ei     ,tl   ,par  ,&
               rbone ,raw  ,rstfac ,cint   ,assim ,respc  ,rs   ,&
#if (defined DGVM)   
               assimpot,&                                               ! added by zhq. 
#endif                                                                   
               ivt)
        else
           rs = 2.e4; assim = 0.; respc = 0.
#if (defined DGVM)                                                      
           assimpot = 0.
#endif           
        endif

! above stomatal resistances are for the canopy, the stomatal rsistances 
! and the "rb" in the following calculations are the average for single leaf. thus,
        rs = rs * lai

!-----------------------------------------------------------------------
! dimensional and non-dimensional sensible and latent heat conductances
! for canopy and soil flux calculations.
!-----------------------------------------------------------------------

        delta = 0.0
        if( (fwet .lt. 0.99) .AND. (sigf*etrc .gt. 1.e-12) )then
           if(qsatl-qaf .gt. 0.) delta = 1.0
        endif
 
        cah = sigf / rah
        cgh = sigf / rd
        cfh = sigf * (lai + sai) / rb

        caw = sigf / raw
        cgw = sigf / rd
        cfw = sigf * ( (1.-delta*(1.-fwet))*(lai+sai)/rb + (1.-fwet)*delta*lai/(rb+rs) )

        wtshi = 1. / ( cah + cgh + cfh )
        wtsqi = 1. / ( caw + cgw + cfw )

        wta0 = cah * wtshi
        wtg0 = cgh * wtshi
        wtl0 = cfh * wtshi

        wtaq0 = caw * wtsqi
        wtgq0 = cgw * wtsqi
        wtlq0 = cfw * wtsqi
 
!-----------------------------------------------------------------------
! IR radiation, sensible and latent heat fluxes and their derivatives
!-----------------------------------------------------------------------
! the partial derivatives of areodynamical resistance are ignored 
! which cannot be determined analtically
        fac = sigf * (1. - thermk)

! longwave absorption and their derivatives 
        irab = (frl - 2. * stefnc * tl**4 + emg*stefnc*tg**4 ) * fac 
        dirab_dtl = - 8. * stefnc * tl**3                      * fac

! sensible heat fluxes and their derivatives
        fsenl = rhoair * cpair * cfh * ( (wta0 + wtg0)*tl - wta0*thm - wtg0*tg )
        fsenl_dtl = rhoair * cpair * cfh * (wta0 + wtg0)

! latent heat fluxes and their derivatives
        etr = sigf * rhoair * (1.-fwet) * delta * lai / (rb + rs) &
            * ( (wtaq0 + wtgq0)*qsatl - wtaq0*qm - wtgq0*qg )
        etr_dtl = sigf * rhoair * (1.-fwet) * delta * lai / (rb + rs) &
		* (wtaq0 + wtgq0)*qsatlDT 
     
        if(etr.ge.sigf*etrc)then
           etr = sigf*etrc
           etr_dtl = 0.
        endif

        evplwet = sigf * rhoair * (1.-delta*(1.-fwet)) * (lai+sai) / rb &
                 * ( (wtaq0 + wtgq0)*qsatl - wtaq0*qm - wtgq0*qg )
        evplwet_dtl = sigf * rhoair * (1.-delta*(1.-fwet)) * (lai+sai) / rb &
		 * (wtaq0 + wtgq0)*qsatlDT 
        if(evplwet.ge.ldew/dtime)then
           evplwet = ldew/dtime
           evplwet_dtl = 0.
        endif
 
        fevpl = etr + evplwet
        fevpl_dtl = etr_dtl + evplwet_dtl

!-----------------------------------------------------------------------
! difference of temperatures by quasi-newton-raphson method for the non-linear system equations
!-----------------------------------------------------------------------

        dtl(it) = (sabv + irab - fsenl - hvap*fevpl) &
            / ((lai+sai)*clai/dtime - dirab_dtl + fsenl_dtl + hvap*fevpl_dtl)

        ! check magnitude of change in leaf temperature limit to maximum allowed value

        if(it .lt. itmax) then

        ! put brakes on large temperature excursions
          if(abs(dtl(it)).gt.delmax)then
              dtl(it) = delmax*dtl(it)/abs(dtl(it))
          endif

          if((it.ge.2) .and. (dtl(it-1)*dtl(it).lt.0.))then
              dtl(it) = 0.5*(dtl(it-1) + dtl(it))
          endif

        endif

        tl = tlbef + dtl(it)

!-----------------------------------------------------------------------
! square roots differences of temperatures and fluxes for use as the condition of convergences
!-----------------------------------------------------------------------

        del  = sqrt( dtl(it)*dtl(it) )
        dele = dtl(it) * dtl(it) * ( dirab_dtl**2 + fsenl_dtl**2 + hvap*fevpl_dtl**2 ) 
        dele = sqrt(dele)

!-----------------------------------------------------------------------
!  saturated vapor pressures and canopy air temperature, canopy air humidity
!-----------------------------------------------------------------------
! Recalculate leaf saturated vapor pressure (ei_)for updated leaf temperature
! and adjust specific humidity (qsatl_) proportionately
        call qsadv(tl,psrf,ei,deiDT,qsatl,qsatlDT)

! update vegetation/ground surface temperature, canopy air temperature, 
! canopy air humidity
        taf = wta0*thm + wtg0*tg + wtl0*tl 

        qaf = wtaq0*qm + wtgq0*qg + wtlq0*qsatl

! update co2 partial pressure within canopy air
        gah2o = 1.0/raw * tprcor/thm                     ! mol m-2 s-1
        gdh2o = 1.0/rd  * tprcor/thm                     ! mol m-2 s-1
        pco2a = pco2m - 1.37*psrf/max(0.446,gah2o) * (assim - respc - rsoil)

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
          um = max(ur,.1)
        else
          wc = (-grav*ustar*thvstar*zii/thv)**(1./3.)
         wc2 = beta*beta*(wc*wc)
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
 
      end do 

! ======================================================================
!     END stability iteration 
! ======================================================================
      z0ma = z0mv
      zol = zeta
      rib = min(5.,zol*ustar**2/(vonkar*temp1*um**2))

! canopy fluxes and total assimilation amd respiration

      if(lai .gt. 0.001) then
         rst = rs / lai
      else
        rst = 2.0e4; assim = 0.; respc = 0. 
#if (defined DGVM)                                                      
        assimpot = 0.
#endif  
      endif
      respc = respc + rsoil

#if (defined DGVM)                            
! In DGVM mode, respc represents only canopy respiration. 
! I think rsoil must be deducted here. zhq. 02/20/2010
      respc = respc - rsoil                        
      annpsn = annpsn + assim * dtime               ! added by zhq. dec27,08
      annpsnpot = annpsnpot + assimpot * dtime
#endif

! canopy fluxes and total assimilation amd respiration
      fsenl = fsenl + fsenl_dtl*dtl(it-1)

      etr     = etr     +     etr_dtl*dtl(it-1)
      evplwet = evplwet + evplwet_dtl*dtl(it-1)
      fevpl   = fevpl   +   fevpl_dtl*dtl(it-1)

      elwmax = ldew/dtime
      elwdif = max(0., evplwet-elwmax)
      evplwet = min(evplwet, elwmax)

      fevpl = fevpl - elwdif
      fsenl = fsenl + hvap*elwdif

      ! wind stresses 
      taux = taux - sigf*rhoair*us/ram
      tauy = tauy - sigf*rhoair*vs/ram

!-----------------------------------------------------------------------
! fluxes from ground to canopy space
!-----------------------------------------------------------------------

      fseng = fseng + cpair*rhoair*cgh*(tg-taf)
      fevpg = fevpg + rhoair*cgw*(qg-qaf)
!print*,'leaftemone->taf',wta0,thm,wtg0,tg,wtl0,tl
!print*,'leaftemone->fseng',ivt,fseng,cpair*rhoair*cgh, tg-taf  
!print*,'leaftemone->fseng',ivt,fevpg 

!-----------------------------------------------------------------------
! downward (upward) longwave radiation below (above) the canopy
!-----------------------------------------------------------------------

      dlrad = sigf * thermk * frl &
              + stefnc * fac * tlbef**3 * (tlbef + 4.*dtl(it-1)) 
      ulrad = stefnc * ( fac * tlbef**3 * (tlbef + 4.*dtl(it-1)) &
              + sigf*thermk*emg*tg**4 )  

!-----------------------------------------------------------------------
! Derivative of soil energy flux with respect to soil temperature (cgrnd)
!-----------------------------------------------------------------------

      cgrnds = cgrnds + cpair*rhoair*cgh*(1.-wtg0)
      cgrndl = cgrndl + rhoair*cgw*(1.-wtgq0)*dqgdT
      cgrnd  = cgrnds + cgrndl*htvp

!-----------------------------------------------------------------------
! balance check
! (the computational error was created by the assumed 'dtl' in line 406-408) 
!-----------------------------------------------------------------------

      err = sabv + irab + dirab_dtl*dtl(it-1) - fsenl - hvap*fevpl
#ifdef MYBUG
      if(abs(err) .gt. .2) &
      write(6,*) 'energy imbalance in leaftemone.F90',it-1,err,sabv,irab,fsenl,hvap*fevpl
#endif

!-----------------------------------------------------------------------
! Update dew accumulation (kg/m2)
!-----------------------------------------------------------------------

      ldew = max(0., ldew-evplwet*dtime)

!-----------------------------------------------------------------------
! 2 m height air temperature
!-----------------------------------------------------------------------
!print*,'temp1',ivt,temp1
      tref = tref + sigf*(thm + temp1*dth * (1./temp12m - 1./temp1)) 
      qref = qref + sigf*( qm + temp2*dqh * (1./temp22m - 1./temp2))

 end subroutine leaftemone
