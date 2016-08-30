  subroutine albland (itypwat,albsol,&
                      chil,ref,tran,fveg,green,lai,sai,coszen,&
                      wt,fsno,scv,sag,ssw,tg,&
                      alb,albg,albv,ssun,ssha,thermk,extkb,extkd) 

!=======================================================================
! Calculates fragmented albedos (direct and diffuse) in
! wavelength regions split at 0.7um.
! 
! (1) soil albedos: as in BATS formulations, which are the function of
!     soil color and moisture in the surface soil layer
! (2) snow albedos: as in BATS formulations, which are inferred from
!     the calculations of Wiscombe and Warren (1980) and the snow model
!     and data of Anderson(1976), and the function of snow age, grain size,
!     solar zenith angle, pollution, the amount of the fresh snow
! (3) canopy albedo: two-stream approximation model 
! (4) glacier albedos: as in BATS, which are set to constants (0.8 for visible beam,
!     0.55 for near-infrared)
! (5) lake and wetland albedos: as in BATS, which depend on cosine solar zenith angle,
!     based on data in Henderson-Sellers (1986). The frozen lake and wetland albedos
!     are set to constants (0.6 for visible beam, 0.4 for near-infrared)
! (6) over the snow covered tile, the surface albedo is estimated by a linear
!     combination of albedos for snow, canopy and bare soil (or lake, wetland, glacier).
!
! Original author : Yongjiu Dai, 09/15/1999; 08/30/2002
!=======================================================================

  use precision
  use phycon_module, only : tfrz
  implicit none

!------------------------- Dummy Arguments -----------------------------
! ground cover index
 integer, INTENT(in) :: &
      itypwat     ! land water type (0=soil, 1=urban or built-up, 2=wetland,
                  ! 3=land ice, 4=deep lake, 5=shallow lake)

                  ! parameters
 real(r8), INTENT(in) :: & 
      albsol,    &! soil albedo for different coloured soils [-]
      chil,      &! leaf angle distribution factor
      ref(2,2),  &! leaf reflectance (iw=iband, il=life and dead)
      tran(2,2), &! leaf transmittance (iw=iband, il=life and dead)
      fveg,      &! fractional vegetation cover [-]
      green,     &! green leaf fraction
      lai,       &! leaf area index (LAI+SAI) [m2/m2]
      sai,       &! stem area index (LAI+SAI) [m2/m2]

                  ! variables
      coszen,    &! cosine of solar zenith angle [-]
      wt,        &! fraction of vegetation covered by snow [-]
      fsno,      &! fraction of soil covered by snow [-]
      ssw,       &! water volumetric content of soil surface layer [m3/m3]
      scv,       &! snow cover, water equivalent [mm]
      sag,       &! non dimensional snow age [-]
      tg          ! ground surface temperature [K]

 real(r8), INTENT(out) :: &
      alb(2,2),  &! averaged albedo [-]
      albg(2,2), &! albedo, ground
      albv(2,2), &! albedo, vegetation [-]
      ssun(2,2), &! sunlit canopy absorption for solar radiation
      ssha(2,2), &! shaded canopy absorption for solar radiation,
                  ! normalized by the incident flux
      thermk,    &! canopy gap fraction for tir radiation
      extkb,     &! (k, g(mu)/mu) direct solar extinction coefficient
      extkd       ! diffuse and scattered diffuse PAR extinction coefficient

!-------------------------- Local variables ----------------------------
 integer         &!
      iw,        &! wavelength (1=visible, 2=near-infrared)
      id,        &! 1=direct, 2=diffuse
      k           ! looping indx

 real(r8) age,       &! factor to reduce visible snow alb due to snow age [-]
      albg0,     &! temporary varaiable [-]
      albsno(2,2),&! snow albedo [-]
      albv0(2),  &! vegetation albedo [-]
      alwet,     &! decrease in soil albedo due to wetness [-]
      beta0,     &! upscattering parameter for direct beam [-]
      cff,       &! snow alb correction factor for zenith angle > 60 [-]
      conn,      &! constant (=0.5) for visible snow alb calculation [-]
      cons,      &! constant (=0.2) for nir snow albedo calculation [-]
      czen,      &! cosine of solar zenith angle > 0 [-]
      czf,       &! solar zenith correction for new snow albedo [-]
      dfalbl,    &! snow albedo for diffuse nir radiation [-]
      dfalbs,    &! snow albedo for diffuse visible solar radiation [-]
      dralbl,    &! snow albedo for visible radiation [-]
      dralbs,    &! snow albedo for near infrared radiation [-]
      fsol1,     &! solar flux fraction for wavelength < 0.7 micron [-]
      fsol2,     &! solar flux fraction for wavelength > 0.7 micron [-]
      lsai,      &! leaf and stem area index (LAI+SAI) [m2/m2]
      scat(2),   &! single scattering albedo for vir/nir beam [-]
      sl,        &! factor that helps control alb zenith dependence [-]
      snal0,     &! alb for visible,incident on new snow (zen ang<60) [-]
      snal1,     &! alb for NIR, incident on new snow (zen angle<60) [-]
      tdiffs,    &! difference of air temperature and freezing temp [K]
      tff,       &! exp(-LSAI)
      tffd,      &! exp(-0.5*LSAI/czen)
      ti,        &! correction due to scattering
      upscat,    &! upward scattered fraction for direct beam [-]
      tranc(2,2),&! canopy transmittances for solar radiation
      zkat(2),   &! temporary
      zkatd(2)    ! temporary
     
! ----------------------------------------------------------------------
! 1. Initial set
! ----------------------------------------------------------------------
! division of solar flux for wavelength less or greater than 0.7 micron
      fsol1 = 0.5      ! shortwave
      fsol2 = 0.5      ! longwave

! short and long wave albedo for new snow
      snal0 = 0.85     ! shortwave
      snal1 = 0.65     ! long wave

! set initial leaf scattering reflectance. Note: "scat" may use different
! value for different vegetation latter
      beta0 = 0.5
      scat(1) = 0.15
      scat(2) = 0.85

! ----------------------------------------------------------------------
! set default soil and vegetation albedos and solar absorption
      alb (:,:) = 0. ! averaged
      albg(:,:) = 0. ! ground
      albv(:,:) = 0. ! vegetation
      ssun(:,:) = 0.
      ssha(:,:) = 0.
     tranc(:,:) = 0.
      thermk = 1.e-3
      extkb = 1.e-6
      extkd = 0.718

      lsai=lai+sai
      if(coszen<=0.) RETURN  !only do albedo when coszen > 0

      czen=max(coszen,0.001) 
      albsno(:,:)=0.         !set initial snow albedo

! ----------------------------------------------------------------------
! 2. albedo for snow cover.
!    snow albedo depends on snow-age, zenith angle, and thickness
!    of snow age gives reduction of visible radiation
! ----------------------------------------------------------------------
      if(scv>0.)then
         cons = 0.2
         conn = 0.5
         sl  = 2.0           !sl helps control albedo zenith dependence

         ! correction for snow age
         age = 1.-1./(1.+sag) !correction for snow age
         dfalbs = snal0*(1.-cons*age)

         ! czf corrects albedo of new snow for solar zenith
         cff    = ((1.+1./sl)/(1.+czen*2.*sl )- 1./sl)
         cff    = max(cff,0.)
         czf    = 0.4*cff*(1.-dfalbs)
         dralbs = dfalbs+czf
         dfalbl = snal1*(1.-conn*age)
         czf    = 0.4*cff*(1.-dfalbl)
         dralbl = dfalbl+czf
   
         albsno(1,1) = dralbs
         albsno(2,1) = dralbl
         albsno(1,2) = dfalbs
         albsno(2,2) = dfalbl
      endif

! ----------------------------------------------------------------------
! 3. get albedo over land
! ----------------------------------------------------------------------
! 3.1 bare soil albedos, depends on moisture
      if(itypwat<=1)then    ! not wetland, permanent ice and water
         alwet = max((11.-40.0*ssw),0.)*0.01
         alwet = min(alwet,albsol)
         albg0 = albsol+alwet
         albg(1,1) = albg0
         albg(2,1) = 2.*albg0
         albg(:,2) = albg(:,1)         !diffused albedos for bare soil

! 3.2 albedos for permanent ice sheet. 
      else if(itypwat==3) then         !permanent ice sheet
         albg(1,:) = 0.8
         albg(2,:) = 0.55

! 3.3 albedo for wetland (swamps, rice paddies etc) and inland water
      else if(itypwat==2 .OR. itypwat>=4) then             
         albg0 = 0.05/(czen+0.15)
         albg(:,:) = albg0

       ! if(tg<tfrz)then               !frozen lake and wetland
         if(tg<tfrz .OR. scv>0.5)then  !frozen lake and wetland
            albg(1,:) = 0.6
            albg(2,:) = 0.4
         endif
      end if

! 3.4 correction due to snow cover
      albg(:,:) = (1.-fsno)*albg(:,:) + fsno*albsno(:,:)
      alb(:,:) = albg(:,:)

! ----------------------------------------------------------------------
! 4. canopy albedos : two stream approximation  
! ----------------------------------------------------------------------
      if(fveg>0.001)then
         call twostream (chil,ref,tran,green,lai,sai,&
                         czen,albg,albv,tranc,thermk,extkb,extkd,ssun,ssha) 

         albv(:,:) = (1.-wt)*albv(:,:) + wt*albsno(:,:)
         alb(:,:) = (1.-fveg)*albg(:,:) + fveg*albv(:,:)
      end if

!-----------------------------------------------------------------------

 end subroutine albland
