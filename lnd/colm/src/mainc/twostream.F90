
 subroutine twostream ( chil, ref,  tran, green, lai, sai, &
            coszen, albg, albv, tranc, thermk, extkb, extkd, ssun, ssha )
                                                                              
!-----------------------------------------------------------------------
!                                                                               
!     calculation of canopy albedos via two stream approximation (direct
!     and diffuse ) and partition of incident solar
!                                                                               
! Original author: Yongjiu Dai, June 11, 2001
!
!-----------------------------------------------------------------------        
                                                                               
  use precision
  implicit none

! parameters
  real(r8), intent(in) :: &
          ! static parameters associated with vegetation type
            chil,          &! leaf angle distribution factor
            ref(2,2),      &! leaf reflectance (iw=iband, il=life and dead)
            tran(2,2),     &! leaf transmittance (iw=iband, il=life and dead)

          ! time-space varying vegetation parameters
            green,         &! green leaf fraction
            lai,           &! leaf area index of exposed canopy (snow-free)
            sai             ! stem area index

! environmental variables
  real(r8), intent(in) :: &
            coszen,        &! consine of solar zenith angle
            albg(2,2)       ! albedos of ground

! output
  real(r8), intent(out) :: &
            albv(2,2),     &! albedo, vegetation [-]
            tranc(2,2),    &! canopy transmittances for solar radiation
            thermk,        &! canopy gap fraction for tir radiation                    
            extkb,         &! (k, g(mu)/mu) direct solar extinction coefficient  
            extkd,         &! diffuse and scattered diffuse PAR extinction coefficient
            ssun(2,2),     &! sunlit canopy absorption for solar radiation
            ssha(2,2)       ! shaded canopy absorption for solar radiation,
                            ! normalized by the incident flux 
                                                                               
!-------------------------- local -----------------------------------           
  real(r8) :: &
           phi1,           &! (phi-1)
           phi2,           &! (phi-2)
           scat,           &! (omega)        
           proj,           &! (g(mu))        
           zmu,            &! (int(mu/g(mu)) 
           zmu2,           &! (zmu * zmu)    
           as,             &! (a-s(mu))      
           upscat,         &! (omega-beta)   
           betao,          &! (beta-0)       
           psi,            &! (h)            

           be,             &! (b)            
           ce,             &! (c)            
           de,             &! (d)            
           fe,             &! (f)            

           power1,         &! (h*lai)
           power2,         &! (k*lai)
           power3,         &!  

           sigma,          &! 
           s1,             &! 
           s2,             &! 
           p1,             &! 
           p2,             &! 
           p3,             &! 
           p4,             &! 
           f1,             &! 
           f2,             &!
           h1,             &!
           h4,             &!
           m1,             &!
           m2,             &!
           m3,             &!
           n1,             &!
           n2,             &!
           n3,             &!

           hh1,            &! (h1/sigma)     
           hh2,            &! (h2)           
           hh3,            &! (h3)           
           hh4,            &! (h4/sigma)     
           hh5,            &! (h5)           
           hh6,            &! (h6)           
           hh7,            &! (h7)           
           hh8,            &! (h8)           
           hh9,            &! (h9)           
           hh10,           &! (h10)          

           eup(2,2),       &! (integral of i_up*exp(-kx) )
           edown(2,2)       ! (integral of i_down*exp(-kx) )
  
  integer iw                !
                                                                                
!-----------------------------------------------------------------------        
! projected area of phytoelements in direction of mu and 
! average inverse diffuse optical depth per unit leaf area

      phi1 = 0.5 - 0.633 * chil - 0.33 * chil * chil
      phi2 = 0.877 * ( 1. - 2. * phi1 )
         
      proj = phi1 + phi2 * coszen
      extkb = (phi1 + phi2 * coszen) / coszen

      extkd = 0.719

      if (abs(phi1).gt.1.e-6 .and. abs(phi2).gt.1.e-6) then
         zmu = 1. / phi2 * ( 1. - phi1 / phi2 * log ( ( phi1 + phi2 ) / phi1 ) )
      else if (abs(phi1).le.1.e-6) then
         zmu = 1./0.877
      else if (abs(phi2).le.1.e-6) then
         zmu = 1./(2.*phi1)
      endif
      zmu2 = zmu * zmu

      power3 = (lai+sai) / zmu
      power3 = min( 50., power3 )
      power3 = max( 1.e-5, power3 )
      thermk = exp(-power3)

      do iw = 1, 2    ! WAVE_BAND_LOOP                                               

!-----------------------------------------------------------------------        
!     calculate average scattering coefficient, leaf projection and             
!     other coefficients for two-stream model.                                  
!-----------------------------------------------------------------------        
                                                                               
      scat = green * ( tran(iw,1) + ref(iw,1) ) + &
             ( 1. - green ) * ( tran(iw,2) + ref(iw,2) ) 

      as = scat / 2. * proj / ( proj + coszen * phi2 )                               
      as = as * ( 1. - coszen * phi1 / ( proj + coszen * phi2 ) * &
               log ( ( proj + coszen * phi2 + coszen * phi1 ) / ( coszen * phi1 ) ) ) 

      upscat = green * tran(iw,1) + ( 1. - green ) * tran(iw,2)                           
      upscat = 0.5 * ( scat + ( scat - 2. * upscat ) * &         
               (( 1. - chil ) / 2. ) ** 2 )                                     
      betao = ( 1. + zmu * extkb ) / ( scat * zmu * extkb ) * as            
                                                                                
!-----------------------------------------------------------------------        
!     intermediate variables identified in appendix of SE-85.                   
!-----------------------------------------------------------------------        
                                                                               
      be = 1. - scat + upscat                                                   
      ce = upscat                                                               
      de = scat * zmu * extkb * betao                                          
      fe = scat * zmu * extkb * ( 1. - betao )                                 

      psi = sqrt(be**2 - ce**2)/zmu                                            
      power1 = min( psi*lai, 50. )                                            
      power2 = min( extkb*lai, 50. )                                          
      s1 = exp( - power1 )                                                    
      s2 = exp ( - power2 )                                                     

!-----------------------------------------------------------------------        
!     calculation of direct albedos and canopy transmittances.                  
!     albv (iw,1)     ( i-up ) 
!     tranc(iw,irad)  ( i-down ) 
!-----------------------------------------------------------------------        

      p1 = be + zmu * psi
      p2 = be - zmu * psi
      p3 = be + zmu * extkb
      p4 = be - zmu * extkb

      f1 = 1. - albg(iw,2)*p1/ce
      f2 = 1. - albg(iw,2)*p2/ce
      
      h1 = - ( de * p4 + ce * fe )
      h4 = - ( fe * p3 + ce * de )

      sigma = ( zmu * extkb ) ** 2 + ( ce**2 - be**2 )                           

      if (abs(sigma) .gt. 1.e-10) then     !<======

         hh1 = h1 / sigma
         hh4 = h4 / sigma
                                                                                
         m1 = f1 * s1
         m2 = f2 / s1
         m3 = ( albg(iw,1) - ( hh1 - albg(iw,2) * hh4 ) ) * s2
 
         n1 = p1 / ce
         n2 = p2 / ce
         n3 = - hh4

         hh2 = (m3*n2 - m2*n3) / (m1*n2 - m2*n1)
         hh3 = (m3*n1 - m1*n3) / (m2*n1 - m1*n2)

         hh5 = hh2 * p1 / ce
         hh6 = hh3 * p2 / ce

         albv(iw,1) = hh1 + hh2 + hh3                                   
         tranc(iw,1) = hh4 * s2 + hh5 * s1 + hh6 / s1                    

         eup(iw,1) = hh1 * (1. - s2*s2) / (2.*extkb) &
                   + hh2 * (1. - s1*s2) / (extkb + psi) &
                   + hh3 * (1. - s2/s1) / (extkb - psi)

         edown(iw,1) = hh4 * (1. - s2*s2) / (2.*extkb) &
                     + hh5 * (1. - s1*s2) / (extkb + psi) &
                     + hh6 * (1. - s2/s1) / (extkb - psi)

      else                               !<======

         m1 = f1 * s1
         m2 = f2 / s1
         m3 = h1 / zmu2 * ( lai + 1. / (2.*extkb) ) * s2 &
            + albg(iw,2) / ce * ( - h1 / (2.*extkb) / zmu2 * &
              ( p3*lai + p4 / (2.*extkb) ) - de ) * s2 &
            + albg(iw,1) * s2
 
         n1 = p1 / ce
         n2 = p2 / ce
         n3 = 1./ce * ( h1*p4 / (4.*extkb*extkb) / zmu2 + de)

         hh2 = (m3*n2 - m2*n3) / (m1*n2 - m2*n1)
         hh3 = (m3*n1 - m1*n3) / (m2*n1 - m1*n2)

         hh5 = hh2 * p1 / ce
         hh6 = hh3 * p2 / ce

         albv(iw,1) =  - h1 / (2.*extkb*zmu2) + hh2 + hh3
         tranc(iw,1) = 1./ce * ( -h1 / (2.*extkb*zmu2) * &
                                ( p3*lai + p4 / (2.*extkb) ) - de ) * s2 &
                     + hh5 * s1 + hh6 / s1

         eup(iw,1) = (hh2 - h1/(2.*extkb*zmu2)) * (1. - s2*s2) / (2.*extkb) &
                   + hh3 * (lai - 0.) &
                   + h1/(2.*extkb*zmu2) * ( lai*s2*s2 - (1. - s2*s2)/(2.*extkb) )

         edown(iw,1) = (hh5 - (h1*p4/(4.*extkb*extkb*zmu) + de)/ce) * &
                             (1. - s2*s2) / (2.*extkb) &
                     + hh6 * (lai - 0.) &
                     + h1*p3/(ce*4.*extkb*extkb*zmu2) * &
                                         ( lai*s2*s2 - (1. - s2*s2)/(2.*extkb) )
      
      endif                              !<======

      ssun(iw,1) = (1.-scat) * ( 1.-s2 + 1. / zmu * (eup(iw,1) + edown(iw,1)) ) 
      ssha(iw,1) = scat * (1.-s2) &
               + ( albg(iw,2)*tranc(iw,1) + albg(iw,1)*s2 - tranc(iw,1) ) - albv(iw,1) &
               - ( 1. - scat ) / zmu * ( eup(iw,1) + edown(iw,1) ) 

!-----------------------------------------------------------------------        
!     calculation of diffuse albedos and canopy transmittances
!     albv (iw,2) ( i-up ) 
!     tranc(iw,2) ( i-down ) 
!-----------------------------------------------------------------------        
                                                                               
      m1 = f1 * s1
      m2 = f2 / s1
      m3 = 0.

      n1 = p1 / ce
      n2 = p2 / ce
      n3 = 1.

      hh7 = -m2 / (m1*n2 - m2*n1)
      hh8 = -m1 / (m2*n1 - m1*n2)

      hh9 = hh7 * p1 / ce
      hh10 = hh8 * p2 / ce

      albv(iw,2) =  hh7 + hh8                                            
      tranc(iw,2) = hh9 * s1 + hh10 / s1

      if (abs(sigma) .gt. 1.e-10) then 
         eup(iw,2)   = hh7 * (1. - s1*s2) / (extkb + psi) &
                     + hh8 * (1. - s2/s1) / (extkb - psi)
         edown(iw,2) = hh9 * (1. - s1*s2) / (extkb + psi) &
                     + hh10 * (1. - s2/s1) / (extkb - psi)
      else
         eup(iw,2)   = hh7 * (1. - s1*s2) / ( extkb + psi) + hh8 * (lai - 0.)
         edown(iw,2) = hh9 * (1. - s1*s2) / ( extkb + psi) + hh10 * (lai - 0.)
      endif

      ssun(iw,2) = (1.-scat) / zmu * (eup(iw,2) + edown(iw,2))
      ssha(iw,2) = tranc(iw,2) * ( albg(iw,2) -1. ) - ( albv(iw,2) - 1. ) &
                 - ( 1. - scat ) / zmu * ( eup(iw,2) + edown(iw,2) ) 

      enddo           ! WAVE_BAND_LOOP

 end subroutine twostream
