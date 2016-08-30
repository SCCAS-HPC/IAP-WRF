
 subroutine lake (nl_lake ,itypwat ,dlat    ,dtime  ,&
                  zlak    ,dzlak   ,zilak   ,hu     ,ht    ,&
                  hq      ,us      ,vs      ,tm     ,qm    ,&
                  prc     ,prl     ,rhoair  ,psrf   ,sabg  ,&
                  frl     ,tg      ,tlak    ,wliq   ,wice  ,&
                  scv     ,snowdp  ,trad    ,tref   ,qref  ,&
                  taux    ,tauy    ,fsena   ,fevpa  ,lfevpa,&
                  fseng   ,fevpg   ,olrg    ,fgrnd  ,tcrit ,&
                  emis    ,z0ma    ,zol     ,rib    ,ustar ,&
                  qstar   ,tstar   ,u10m    ,v10m   ,f10m  ,&
                  fm      ,fh      ,fq)

! ------------------------ code history ---------------------------
! purpose:           lake temperature and snow on frozen lake
! initial author:    Gordon Bonan
! revised:           Yongjiu Dai, September 15, 1999
!
! ------------------------ notes ----------------------------------
! calculate lake temperatures from one-dimensional thermal
! stratification model based on eddy diffusion concepts to 
! represent vertical mixing of heat
!
! d ts    d            d ts     1 ds
! ---- = -- [(km + ke) ----] + -- --
!  dt    dz             dz     cw dz   
 
! where: ts = temperature (kelvin)
!         t = time (s)
!         z = depth (m)
!        km = molecular diffusion coefficient (m**2/s)
!        ke = eddy diffusion coefficient (m**2/s)
!        cw = heat capacity (j/m**3/kelvin)
!         s = heat source term (w/m**2)
 
! there are two types of lakes: 
!    deep lakes are 50 m. shallow lakes are 10 m deep.
!    for unfrozen deep lakes:    ke > 0 and    convective mixing
!    for unfrozen shallow lakes: ke = 0 and no convective mixing
 
! use crank-nicholson method to set up tridiagonal system of equations to
! solve for ts at time n+1, where the temperature equation for layer i is
! r_i = a_i [ts_i-1] n+1 + b_i [ts_i] n+1 + c_i [ts_i+1] n+1
 
! the solution conserves energy as
 
! cw*([ts(  1)] n+1 - [ts(  1)] n)*dz(  1)/dt + ... +
! cw*([ts(nl_lake)] n+1 - [ts(nl_lake)] n)*dz(nl_lake)/dt = fin
 
! where 
! [ts] n   = old temperature (kelvin)
! [ts] n+1 = new temperature (kelvin)
! fin      = heat flux into lake (w/m**2)
!          = beta*sabg+frl-olrg-fsena-lfevpa-hm + phi(1) + ... + phi(nl_lake) 
 
! -----------------------------------------------------------------

      use precision
      use phycon_module, only : tfrz,hvap,hfus,hsub,tkwat,tkice,stefnc,&
                                vonkar,grav,cpliq,cpair,denh2o,rgas
      implicit none
      
! ------------------------ input/output variables -----------------

  integer, INTENT(in) :: &
        nl_lake,   &!number of soil layers
        itypwat     !land water type (4=deep lake, 5=shallow lake)

  real(r8), INTENT(in) :: &
        dlat,      &!latitude (radians)
        dtime,     &!time step (s)
        tcrit,     &!critical temp. to determine rain or snow
        dzlak(nl_lake),&!soil layer thickness (m)
        zlak(nl_lake), &!depth (m)
        zilak(0:nl_lake), &!
        hu,        &!observational height of wind [m]
        ht,        &!observational height of temperature [m]
        hq,        &!observational height of humidity [m]
        us,        &!wind component in eastward direction [m/s]
        vs,        &!wind component in northward direction [m/s]
        tm,        &!temperature at agcm reference height [kelvin]
        qm,        &!specific humidity at agcm reference height [kg/kg]
        prc,       &!convective precipitation [mm/s]
        prl,       &!large scale precipitation [mm/s]
        rhoair,    &!density air [kg/m3]
        psrf,      &!atmosphere pressure at the surface [pa]
        sabg,      &!solar radiation absorbed by ground [W/m2]
        frl         !atmospheric infrared (longwave) radiation [W/m2]

  real(r8), INTENT(inout) :: &
        tg,        &!surface temperature (kelvin)
        tlak(nl_lake), &!lake temperature (kelvin)
        wliq(nl_lake), &!
        wice(nl_lake), &!
        scv,       &!snow water equivalent [mm]
        snowdp      !snow depth [m]

  real(r8), INTENT(out) :: &
        taux,      &!wind stress: E-W [kg/m/s**2]
        tauy,      &!wind stress: N-S [kg/m/s**2]
        fsena,     &!sensible heat from canopy height to atmosphere [W/m2]
        fevpa,     &!evapotranspiration from canopy height to atmosphere [mm/s]
        lfevpa,    &!latent heat flux from canopy height to atmosphere [W/m2]
        fseng,     &!sensible heat flux from ground [W/m2]
        fevpg,     &!evaporation heat flux from ground [mm/s]
        olrg,      &!outgoing long-wave radiation from ground+canopy
        fgrnd,     &!ground heat flux [W/m2]

        tref,      &!2 m height air temperature [kelvin]
        qref,      &!2 m height air specific humidity
        trad,      &!radiative temperature [K]

        emis,      &!averaged bulk surface emissivity
        z0ma,      &!effective roughness [m]
        zol,       &!dimensionless height (z/L) used in Monin-Obukhov theory
        rib,       &!bulk Richardson number in surface layer
        ustar,     &!u* in similarity theory [m/s]
        qstar,     &!q* in similarity theory [kg/kg]
        tstar,     &!t* in similarity theory [K]
        u10m,      &!10m u-velocity
        v10m,      &!10m v-velocity
        f10m,      &!integral of profile function for momentum at 10m
        fm,        &!integral of profile function for momentum
        fh,        &!integral of profile function for heat
        fq          !integral of profile function for moisture

! ------------------------ local variables ------------------------
  integer &
        niters,    &!maximum number of iterations for surface temperature
        iter,      &!iteration index
        nmozsgn,   &!number of times moz changes sign
        idlak       !index of lake, 1 = deep lake, 2 = shallow lake

  real(r8)  ax,    &!
        bx,        &!
        beta1,     &!coefficient of conective velocity [-]
        degdT,     &!d(eg)/dT
        displax,   &!zero- displacement height [m]
        dqh,       &!diff of humidity between ref. height and surface
        dth,       &!diff of virtual temp. between ref. height and surface
        dthv,      &!diff of vir. poten. temp. between ref. height and surface
        dzsur,     &!
        eg,        &!water vapor pressure at temperature T [pa]
        emg,       &!ground emissivity (0.97 for snow,
        errore,    &!lake temperature energy conservation error (w/m**2)
        hm,        &!energy residual [W/m2]
        htvp,      &!latent heat of vapor of water (or sublimation) [j/kg]
        obu,       &!monin-obukhov length (m)
        obuold,    &!monin-obukhov length of previous iteration
        qsatg,     &!saturated humidity [kg/kg]
        qsatgdT,   &!d(qsatg)/dT
        qseva,     &!ground surface evaporation rate (mm h2o/s)
        qsdew,     &!ground surface dew formation (mm h2o /s) [+]
        qsubl,     &!sublimation rate from snow pack (mm h2o /s) [+]
        qfros,     &!surface dew added to snow pack (mm h2o /s) [+]
        qmelt,     &!snow melt [mm/s]
        ram,       &!aerodynamical resistance [s/m]
        rah,       &!thermal resistance [s/m]
        raw,       &!moisture resistance [s/m]
        snowrate,  &!rate of snowfall [mm/s]
        stftg3,    &!
        temp1,     &!relation for potential temperature profile
        temp2,     &!relation for specific humidity profile
        temp12m,   &! relation for temperature at 2m
        temp22m,   &! relation for specific humidity at 2m
        tgbef,     &!
        thm,       &!intermediate variable (tm+0.0098*ht)
        th,        &!potential temperature (kelvin)
        thv,       &!virtual potential temperature (kelvin)
        thvstar,   &!virtual potential temperature scaling parameter
        tksur,     &!thermal conductivity of snow/soil (w/m/kelvin)

        um,        &!wind speed including the stablity effect [m/s]
        ur,        &!wind speed at reference height [m/s]
        visa,      &! kinematic viscosity of dry air [m2/s]
        wc,        &!convective velocity [m/s]
        wc2,       &!wc*wc
        xt,        &!
        xq,        &!
        zeta,      &!dimensionless height used in Monin-Obukhov theory
        zii,       &!convective boundary height [m]
        zldis,     &!reference height "minus" zero displacement heght [m]
        z0mg,      &!roughness length over ground, momentum [m]
        z0hg,      &!roughness length over ground, sensible heat [m]
        z0qg,      &!roughness length over ground, latent heat [m]

        beta(2),   &!fraction solar rad absorbed at surface: depends on lake type
        za(2),     &!base of surface absorption layer (m): depends on lake type
        eta(2),    &!light extinction coefficient (/m): depends on lake type
        p0,        &!neutral value of turbulent prandtl number
     
        a(nl_lake),    &!"a" vector for tridiagonal matrix
        b(nl_lake),    &!"b" vector for tridiagonal matrix
        c(nl_lake),    &!"c" vector for tridiagonal matrix
        r(nl_lake),    &!"r" vector for tridiagonal solution
        rhow(nl_lake), &!density of water (kg/m**3)
        phi(nl_lake),  &!solar radiation absorbed by layer (w/m**2)
        kme(nl_lake),  &!molecular + eddy diffusion coefficient (m**2/s)

        cwat,      &!specific heat capacity of water (j/m**3/kelvin)
        ws,        &!surface friction velocity (m/s)
        ks,        &!coefficient
        in,        &!relative flux of solar radiation into layer
        out,       &!relative flux of solar radiation out of layer
        ri,        &!richardson number
        fin,       &!heat flux into lake - flux out of lake (w/m**2)
        ocvts,     &!(cwat*(tlak[n  ])*dzlak
        ncvts,     &!(cwat*(tlak[n+1])*dzlak

        m1,        &!intermediate variable for calculating r, a, b, c
        m2,        &!intermediate variable for calculating r, a, b, c
        m3,        &!intermediate variable for calculating r, a, b, c
        ke,        &!eddy diffusion coefficient (m**2/s)
        km,        &!molecular diffusion coefficient (m**2/s)
        zin,       &!depth at top of layer (m)
        zout,      &!depth at bottom of layer (m)
        drhodz,    &!d [rhow] /dz (kg/m**4)
        n2,        &!brunt-vaisala frequency (/s**2)
        num,       &!used in calculating ri
        den,       &!used in calculating ri
        tav,       &!used in aver temp for convectively mixed layers
        nav,       &!used in aver temp for convectively mixed layers
        phidum,    &!temporary value of phi
        u2m         !2 m wind speed (m/s)

  integer i,j       !do loop or array index

! -----------------------------------------------------------------
!*[1] constants and model parameters
! -----------------------------------------------------------------
! constants for lake temperature model
      beta = (/0.4, 0.4/)                              ! (deep lake, shallow lake)
      za   = (/0.6, 0.5/)    
      eta  = (/0.1, 0.5/)  
      p0   = 1.  

! latent heat 
      if (tm > tfrz) then
         htvp = hvap
      else
         htvp = hsub
      end if

! surface emissivity
      emg = 0.97

! deep lake or shallow
      idlak = 1
      if(itypwat==5) idlak = 2

      snowrate = 0
      if((prc+prl)>0. .and. (tm<=tfrz+tcrit)) snowrate = prc+prl

! ----------------------------------------------------------------------
!*[2] surface temperature and fluxes
! ----------------------------------------------------------------------

      dzsur = dzlak(1) + snowdp

      call qsadv(tg,psrf,eg,degdT,qsatg,qsatgdT)

! potential temperatur at the reference height

      beta1=1.       ! -  (in computing W_*)
      zii = 1000.    ! m  (pbl height)

      thm = tm + 0.0098*ht              ! intermediate variable equivalent to
                                        ! tm*(pgcm/psrf)**(rgas/cpair)
      th = tm*(100000./psrf)**(rgas/cpair) ! potential T
      thv = th*(1.+0.61*qm)             ! virtual potential T
      ur = max(0.1,sqrt(us*us+vs*vs))   ! limit set to 1

! Initialization variables

      nmozsgn = 0
      obuold = 0.

      dth   = thm-tg
      dqh   = qm-qsatg
      dthv  = dth*(1.+0.61*qm)+0.61*th*dqh
      zldis = hu-0.

! aerodynamical roughness 
! Kinematic viscosity of dry air (m2/s)- Andreas (1989) CRREL Rep. 89-11

    ! 2011.07.28 (Yuqing Wang)
    ! visa=1.326e-5*(1.+6.542e-3*tm + 8.301e-6*tm**2 - 4.84e-9*tm**3)
      visa=1.326e-5*(1.+6.542e-3*(tm-tfrz) + 8.301e-6*(tm-tfrz)**2 - 4.84e-9*(tm-tfrz)**3)
      ustar=0.06
      wc=0.5

      if(tg.ge.tfrz)then          ! unfrozen lake
       ! loop to obtain initial and good ustar and zo
         if(dthv.ge.0.) then
            um=max(ur,0.1)
         else
            um=sqrt(ur*ur+wc*wc)
         endif
         do i=1,5
            z0mg=0.013*ustar*ustar/grav+0.11*visa/ustar
            ustar=vonkar*um/log(zldis/z0mg)
         enddo
      else                        ! frozen lake
         z0mg = 0.04
      endif
      z0qg = z0mg
      z0hg = z0mg

      call moninobukini(ur,th,thm,thv,dth,dqh,dthv,zldis,z0mg,um,obu)

! ----------------------------------------------------------------------
      niters = 3

      do iter = 1, niters         ! begin stability iteration
         tgbef = tg
         if(tg.ge.tfrz) then
            tksur = tkwat
         else
            tksur = tkice
         end if

         if(tg.ge.tfrz)then       ! unfrozen lake
            z0mg=0.013*ustar*ustar/grav + 0.11*visa/ustar
            xq=2.67*(ustar*z0mg/visa)**0.25 - 2.57
            xt= xq
            z0qg=z0mg/exp(xq)
            z0hg=z0mg/exp(xt)
         endif

! Evaluated stability-dependent variables using moz from prior iteration
         displax = 0.
         call moninobuk(hu,ht,hq,displax,z0mg,z0hg,z0qg,obu,um,&
                        ustar,temp1,temp2,temp12m,temp22m,f10m,fm,fh,fq)

         obuold = obu
!
! Get derivative of fluxes with repect to ground temperature

         ram    = 1./(ustar*ustar/um)
         rah    = 1./(temp1*ustar)
         raw    = 1./(temp2*ustar)

         stftg3 = emg*stefnc*tgbef*tgbef*tgbef

         ax  = sabg + emg*frl + 3.*stftg3*tgbef &
             + rhoair*cpair/rah*thm &
             - htvp*rhoair/raw*(qsatg-qsatgdT*tgbef - qm) &
             + tksur*tlak(1)/dzsur
 
         bx  = 4.*stftg3 + rhoair*cpair/rah &
             + htvp*rhoair/raw*qsatgdT + tksur/dzsur
 
         tg = ax/bx

! surface fluxes of momentum, sensible and latent
! using ground temperatures from previous time step

         fseng = rhoair*cpair*(tg-thm)/rah
         fevpg = rhoair*(qsatg+qsatgdT*(tg-tgbef)-qm)/raw
 
         call qsadv(tg,psrf,eg,degdT,qsatg,qsatgdT)
         dth=thm-tg
         dqh=qm-qsatg

         tstar = temp1*dth
         qstar = temp2*dqh

       ! thvstar=tstar+0.61*th*qstar
         thvstar=tstar*(1+0.61*qm)+0.61*th*qstar
         zeta=zldis*vonkar*grav*thvstar/(ustar**2*thv)
         if(zeta >= 0.) then     !stable
           zeta = min(2.,max(zeta,1.e-6))
         else                    !unstable
           zeta = max(-100.,min(zeta,-1.e-6))
         endif
         obu = zldis/zeta

         if(zeta >= 0.)then
           um = max(ur,0.1)
         else
           wc = (-grav*ustar*thvstar*zii/thv)**(1./3.)
          wc2 = beta1*beta1*(wc*wc)
           um = sqrt(ur*ur+wc2)
         endif

         if (obuold*obu < 0.) nmozsgn = nmozsgn+1
         if(nmozsgn >= 4) EXIT

      enddo
! ----------------------------------------------------------------------

! if snow on ground and tg > tfrz: reset tg = tfrz. reevaluate ground fluxes.
! energy inbalance used to melt snow. scv > 0.5 prevents spurious fluxes

      if (scv > 0.5 .AND. tg > tfrz) then
         tg = tfrz
         fseng = rhoair*cpair*(tg-thm)/rah
         fevpg = rhoair*(qsatg+qsatgdT*(tg-tgbef)-qm)/raw !*qsatg and qsatgdT
                                                          !*should be f(tgbef)
      end if

! net longwave from ground to atmosphere
      olrg = (1.-emg)*frl + stftg3*(-3.*tgbef+4.*tg)

! additional variables for WRF model
      emis = emg
      z0ma = z0mg
      zol  = zeta
      rib  = min(5.,zol*ustar**2/(vonkar*temp1*um**2))

! radiative temperature
      trad = (olrg/stefnc)**0.25

! ground heat flux
      fgrnd = sabg + frl - olrg - fseng - htvp*fevpg

      taux   = -rhoair*us/ram
      tauy   = -rhoair*vs/ram

      fsena  = fseng
      fevpa  = fevpg
      lfevpa = htvp*fevpg

! 2 m height air temperature
      tref   = thm + temp1*dth * (1./temp12m - 1./temp1)
      qref   =  qm + temp2*dqh * (1./temp22m - 1./temp2)

! 10 m wind
      u10m = us/max(0.1,ur) * ustar/vonkar * f10m
      v10m = vs/max(0.1,ur) * ustar/vonkar * f10m

! energy residual for snow melting
      if (scv > 0. .AND. tg >= tfrz) then
         hm = min( scv*hfus/dtime, max(fgrnd,0.) )
      else
         hm = 0.
      end if
      qmelt = hm/hfus             ! snow melt (mm/s)

! ----------------------------------------------------------------------
!*[3] lake layer temperature
! ----------------------------------------------------------------------

! lake density

      do j = 1, nl_lake
         rhow(j) = 1000.*( 1.0 - 1.9549e-05*(abs(tlak(j)-277.))**1.68 )
      end do

! eddy diffusion +  molecular diffusion coefficient:
! eddy diffusion coefficient used for unfrozen deep lakes only

      cwat = cpliq*denh2o
      km = tkwat/cwat

      fin = beta(idlak)*sabg + frl - (olrg+fsena+lfevpa+hm)
      u2m = max(1.0,ustar/vonkar*log(2./z0mg))

      ws = 1.2e-03 * u2m
      ks = 6.6 * sqrt( abs(sin(dlat)) ) * (u2m**(-1.84))

      do j = 1, nl_lake-1
         drhodz = (rhow(j+1)-rhow(j)) / (zlak(j+1)-zlak(j))
       ! 2011.07.28 (Xin-Zhong Liang)
       ! n2  = -grav / rhow(j) * drhodz
         n2  = max(7.5e-5, grav / rhow(j) * drhodz)
         num = 40. * n2 * (vonkar*zlak(j))**2
         den = max( (ws**2) * exp(-2.*ks*zlak(j)), 1.e-10 )
         ri = ( -1. + sqrt( max(1.+num/den, 0.) ) ) / 20.
         if (idlak == 1 .AND. tg > tfrz) then
            ke = vonkar*ws*zlak(j)/p0 * exp(-ks*zlak(j)) / (1.+37.*ri*ri)
         else
            ke = 0.
         end if
         kme(j) = km + ke 
      end do

      kme(nl_lake) = kme(nl_lake-1)

! heat source term: unfrozen lakes only

      do j = 1, nl_lake
         zin  = zlak(j) - 0.5*dzlak(j)
         zout = zlak(j) + 0.5*dzlak(j)
         in  = exp( -eta(idlak)*max(  zin-za(idlak),0. ) )
         out = exp( -eta(idlak)*max( zout-za(idlak),0. ) )
         !assumed solar absorption is only in the considered depth
         if(j == nl_lake) out = 0.  
         if (tg > tfrz) then
            phidum = (in-out) * sabg * (1.-beta(idlak))
         else if (j == 1) then
            phidum= sabg * (1.-beta(idlak))
         else
            phidum = 0.
         end if
         phi(j) = phidum
      end do

! sum cwat*tlak*dzlak for energy check

      ocvts = 0.
      do j = 1, nl_lake
         ocvts = ocvts + cwat*tlak(j)*dzlak(j) 
      end do

! set up vector r and vectors a, b, c that define tridiagonal matrix

      j = 1
      m2 = dzlak(j)/kme(j) + dzlak(j+1)/kme(j+1)
      m3 = dtime/dzlak(j)
      r(j) = tlak(j) + (fin+phi(j))*m3/cwat - (tlak(j)-tlak(j+1))*m3/m2
      a(j) = 0.
      b(j) = 1. + m3/m2
      c(j) = -m3/m2

      j = nl_lake
      m1 = dzlak(j-1)/kme(j-1) + dzlak(j)/kme(j)
      m3 = dtime/dzlak(j)
      r(j) = tlak(j) + phi(j)*m3/cwat + (tlak(j-1)-tlak(j))*m3/m1
      a(j) = -m3/m1
      b(j) = 1. + m3/m1
      c(j) = 0.

      do j = 2, nl_lake-1
         m1 = dzlak(j-1)/kme(j-1) + dzlak(j  )/kme(j  )
         m2 = dzlak(j  )/kme(j  ) + dzlak(j+1)/kme(j+1)
         m3 = dtime/dzlak(j)
         r(j) = tlak(j) + phi(j)*m3/cwat &
              +(tlak(j-1)-tlak(j))*m3/m1 - (tlak(j)-tlak(j+1))*m3/m2

         a(j) = -m3/m1
         b(j) = 1. + m3/m1 + m3/m2
         c(j) = -m3/m2
      end do

! solve for tlak: a, b, c, r, u go from 1 to nsoi. tlak = 1 to npt

      call tridia (nl_lake ,a ,b ,c ,r ,tlak) 

! convective mixing: make sure cwat*dzlak*ts is conserved. mixing
! is only allowed for unfrozen deep lakes. mix every 3 time steps

      if(idlak == 1 .AND. tg > tfrz) then
         do j = 1, nl_lake-1
            if(rhow(j) > rhow(j+1)) then
               tav = 0.
               nav = 0.
               do i = 1, j+1
                  tav = tav + tlak(i)*dzlak(i)
                  nav = nav + dzlak(i)
               end do
               tav = tav/nav
  
               do i = 1, j+1
                  tlak(i) = tav
                  rhow(i) = 1000.*( 1.0 &
                          - 1.9549e-05*(abs(tlak(i)-277.))**1.68 )
               end do

            end if
         end do
      end if

! sum cwat*tlak*dzlak and total energy into lake for energy check

      ncvts = 0.
      do j = 1, nl_lake
         ncvts = ncvts + cwat*tlak(j)*dzlak(j) 
         fin = fin + phi(j)
      end do

      errore = (ncvts-ocvts) / dtime - fin

! ----------------------------------------------------------------------
! [4] snow on the lake ice 
! ----------------------------------------------------------------------

      qseva = 0.
      qsubl = 0.
      qfros = 0.
      qsdew = 0.

      if(fevpg >= 0.)then
! sublimation. do not allow for more sublimation than there is snow
! after melt. remaining surface evaporation used for infiltration
         qsubl = min( fevpg, scv/dtime-qmelt )
         qseva = fevpg - qsubl
      else
         if(tg < tfrz-0.1)then
            qfros = abs(fevpg)
         else
            qsdew = abs(fevpg)
         endif
      endif

! update snow pack
      scv = scv + (snowrate-qmelt-qsubl+qfros)*dtime
      scv = max( scv, 0. )

! no snow if lake unfrozen
      if (tg > tfrz) scv = 0.

! snow height and fractional coverage
      snowdp = scv/250.       !assumed a constant snow bulk density = 250.

! null water mass
      do j = 1, nl_lake
         wice(j) = 0.
         wliq(j) = 0.
      enddo

 end subroutine lake
