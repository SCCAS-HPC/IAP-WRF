module scavrat_c
contains
     SUBROUTINE scavrat (laero, lcloud, lgraupl, tamin, pmm, temp,     &
      cwat, dz, rhoair, conc, hlaw, difrat, rscale, prtdia, rhoprt,     &
      gscav, ascav)                                                     
!                                                                       
!----CAMx v4.42 070603                                                  
!                                                                       
!     SCAVRAT calculates wet scavenging rates for gases and aerosols.   
!     Rates are determined for:                                         
!       1) Uptake of cloud water with dissolved gasses                  
!       2) Uptake of ambient gasses into precip                         
!       3) Uptake of cloud water with PM (all PM in cloudy layers is ass
!          to reside in cloud water)                                    
!       4) Uptake of ambient PM into precip, dependent on particle size 
!          ice form                                                     
!     Super-cooled liquid cloud water is assumed to exist in the tempera
!     range tamin < T < 273K using a ramp function (100% liquid at 273 t
!     0% liquid at tamin).                                              
!                                                                       
!     Copyright 1996-2007                                               
!     ENVIRON International Corporation                                 
!                                                                       
!     Modifications:                                                    
!        02/11/04   Updated empirical rainfall relationships            
!        05/03/05   Revised gas scavenging calculations                 
!        06/20/05   Revised to handle liquid/frozen cloud and precip wat
!        10/21/05   Henry's Law calculation moved out to WETDEP         
!                                                                       
!     Input arguments:                                                  
!        laero               aerosol flag                               
!        lcloud              in-cloud flag (F=no cloud water)           
!        lgraupl             graupel frozen precip flag                 
!        tamin               minimum temperature for liquid cloud water 
!        pmm                 precip rate (mm/hr)                        
!        temp                temperature (K)                            
!        cwat                cloud water content (g/m3)                 
!        dz                  cell depth (m)                             
!        rhoair              atmospheric density (kg/m3)                
!        conc                cell gas concentration (umol/m3)           
!        hlaw                Henry's Law constant (M/M)                 
!        difrat              Ratio of H2O to gas diffusivity            
!        prtdia              mean aerosol size (m)                      
!        rhoprt              aerosol density (g/m3)                     
!                                                                       
!     Output arguments:                                                 
!        gscav               Gas scavenging rate (1/s)                  
!        ascav               Aerosol scavenging rate (1/s)              
!                                                                       
!     Routines called:                                                  
!        none                                                           
!                                                                       
!     Called by:                                                        
!        WETDEP                                                         
!                                                                       
      IMPLICIT none 
      REAL tamin, pmm, temp, cwat, dz, rhoair, conc, difrat, rscale,    &
      prtdia, rhoprt, hlaw, gscav, ascav                                
      LOGICAL laero, lcloud, lgraupl 
!                                                                       
      REAL pi, rhoh2o, difh2o, boltz, xmfp, cldeff, dscav, gdiff,       &
      drpdia, drpvel, cscav, cgas, caq, diff, term1, term2, expo,       &
      reynold, power, scf, difbrwn, schmidt, stoke, top, bot, star,     &
      terma, termb, phi, term3, eff, drpmas                             
      REAL nuair, muair, muh2o, kc 
!                                                                       
!-----Constants                                                         
!                                                                       
      DATA pi / 3.1415927 / 
                             !g/m3                                      
      DATA rhoh2o / 1.e6 / 
                             !m2/s                                      
      DATA difh2o / 2.3e-5 / 
                             !kg/ms                                     
      DATA muair / 1.8e-5 / 
                             !kg/ms                                     
      DATA muh2o / 1.e-3 / 
                             !J/K                                       
      DATA boltz / 1.38e-23 / 
                             !m                                         
      DATA xmfp / 6.5e-8 / 
      DATA cldeff / 0.9 / 
!                                                                       
!-----Entry point                                                       
!                                                                       
      dscav = 0. 
      gscav = 0. 
      ascav = 0. 
      gdiff = 1. 
      eff = 1. 
!                                                                       
!-----Calculate environmental parameters                                
!                                                                       
                                           !air molecular diffusivity (m
      nuair = muair / rhoair 
                                           !rain drop diameter (m)      
      drpdia = 9.0e-4 * (pmm**0.21) 
                                           !rain drop mass (g)          
      drpmas = 1.e6 * (pi / 6.) * (drpdia) **3 
      IF (temp.gt.273.) then 
                                           !rain drop fall speed (m/s)  
         drpvel = 3.1e3 * drpdia 
      ELSEIF (lgraupl) then 
                                           !graupel diameter (mm)       
         drpdia = (17. * 1000. * drpmas) **0.38 
                                           !graupel fall speed (m/s)    
         drpvel = 1.10 * drpdia**0.61 
                                           !(m)                         
         drpdia = 1.e-3 * drpdia 
      ELSE 
                                           !snow diameter (mm)          
         drpdia = (29. * 1000. * drpmas) **0.56 
                                           !snow fall speed (m/s)       
         drpvel = 0.83 * drpdia**0.20 
                                           !(m)                         
         drpdia = 1.e-3 * drpdia 
      ENDIF 
                                           !cloud scavenging rate (1/s) 
      cscav = 4.2e-7 * pmm * cldeff / drpdia 
                                                                        
      IF (laero) goto 1000 
!                                                                       
!-----Gas scavenging                                                    
!                                                                       
      IF (temp.le.tamin) return 
      cgas = conc 
      caq = 0. 
!                                                                       
!-----If in cloud, partition total gas into equilibrium aqueous and gas 
!     and calculate scavenging rate for dissolved gasses in cloud water 
!                                                                       
      IF (lcloud) then 
         cgas = conc / (1. + hlaw * cwat / rhoh2o) 
         caq = conc - cgas 
         dscav = cscav * caq / conc 
      ENDIF 
!                                                                       
!-----Calculate scavenging rate for ambient gas dissolving into precip  
!                                                                       
      IF (temp.gt.273..or. (temp.le.273..and.rscale.eq.0.) ) then 
         diff = difh2o / difrat 
         term1 = (drpvel * drpdia / nuair) **0.5 
         term2 = (nuair / diff) **0.333 
         kc = diff / drpdia * (2. + 0.6 * term1 * term2) 
         expo = - 6. * kc * dz / (drpdia * drpvel * hlaw) 
         gdiff = 1. - exp (expo) 
         gscav = 2.8e-7 * pmm * cgas * hlaw * gdiff / (dz * conc) 
      ENDIF 
      gscav = gscav + dscav 
!                                                                       
      RETURN 
!                                                                       
!-----Aerosol scavenging                                                
!                                                                       
 1000 CONTINUE 
!                                                                       
!-----Scavenging rate for aerosol in cloud water is set equal to        
!     cloud water scavenging rate                                       
!                                                                       
      IF (lcloud) then 
         ascav = cscav 
      ELSE 
!                                                                       
!-----Calculate scavenging rate of dry aerosols below cloud as f(size)  
!                                                                       
         reynold = drpdia * drpvel / (2. * nuair) 
         power = amin1 (7.6, 0.55 * prtdia / xmfp) 
         scf = 1. + (2.514 + 0.8 * exp ( - power) ) * xmfp / prtdia 
         difbrwn = boltz * temp * scf / (3. * pi * muair * prtdia) 
         schmidt = nuair / difbrwn 
         stoke = drpvel * prtdia * prtdia * rhoprt * scf / (9000. *     &
         muair * drpdia)                                                
                                                                        
         top = 1.2 + alog (1. + reynold) / 12. 
         bot = 1.0 + alog (1. + reynold) 
         star = amin1 (top / bot, stoke) 
                                                                        
         terma = reynold**0.5 * schmidt**0.333 
         termb = reynold**0.5 * schmidt**0.5 
         term1 = 4. / (reynold * schmidt) * (1. + 0.4 * terma + 0.16 *  &
         termb)                                                         
         phi = prtdia / drpdia 
         term2 = 4. * phi * (muair / muh2o + (1. + 2. * reynold**0.5)   &
         * phi)                                                         
         term3 = (stoke-star) / (stoke-star + 2. / 3.) 
         term3 = term3**1.5 
         eff = term1 + term2 + term3 
         IF (temp.le.273..and..not.lgraupl) eff = amax1 (eff, 1.e-3) 
                                                                        
         ascav = cscav * eff / cldeff 
      ENDIF 
                                                                        
      RETURN 
      END SUBROUTINE scavrat                        
end module scavrat_c
