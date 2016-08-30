#include <define.h>

subroutine lpjreset(itypwat,albsol,nfcon_pft,nfvar_col,nfvar_pft,n_pft,fcon_pft,fvar_col,fvar_pft,&
                    ivt_r,ivt_r_old,ifpre,ifpre_old,wt_patch,wt_patch_old,prec365)

   use precision
   implicit none

   integer,  intent(in)    :: nfcon_pft
   integer,  intent(in)    :: nfvar_col
   integer,  intent(in)    :: nfvar_pft
   integer,  intent(in)    :: n_pft
   real(r8), target, intent(in) :: fcon_pft(nfcon_pft,n_pft)
   real(r8), target, intent(inout) :: fvar_col(nfvar_col)
   real(r8), target, intent(inout) :: fvar_pft(nfvar_pft,n_pft)
   real(r8), intent(in)    :: ivt_r(n_pft)
   real(r8), intent(in)    :: ivt_r_old(n_pft)
   real(r8), intent(in)    :: ifpre(n_pft)
   real(r8), intent(in)    :: ifpre_old(n_pft)
   real(r8), intent(in)    :: wt_patch(n_pft)
   real(r8), intent(in)    :: wt_patch_old(n_pft)
   real(r8), intent(inout) :: prec365   !

! ----------------------------------------------------------------
! I. Time invariant model variables
! ----------------------------------------------------------------

  integer,   intent(in)    :: itypwat   ! land water type
  real(r8),  intent(in)    :: albsol    ! soil albedo for different coloured soils [-]
  integer  :: ivt, ivt_old              ! land cover type  

! -----------------------------------------------------------------
! II. Time-varying state variables which reaquired by restart run
! -----------------------------------------------------------------

                       ! Main land surface variables 
  real(r8), pointer :: z   (:)          ! node depth [m]
  real(r8), pointer :: dz  (:)          ! interface depth [m]
  real(r8), pointer :: tss (:)          ! soil temperature [K]
  real(r8), pointer :: wliq(:)          ! liquid water in layers [kg/m2]
  real(r8), pointer :: wice(:)          ! ice lens in layers [kg/m2]
  real(r8), pointer :: tg               ! ground surface temperature [K]
  real(r8), pointer :: sag              ! non dimensional snow age [-]
  real(r8), pointer :: scv              ! snow cover, water equivalent [mm]
  real(r8), pointer :: snowdp           ! snow depth [meter]
  real(r8), pointer :: fsno             ! fraction of snow cover on ground
  real(r8), pointer :: coszen           ! cosine of solar zenith angle

  real(r8), pointer :: tlsun            ! sunlit leaf temperature [K]
  real(r8), pointer :: tlsha            ! shaded leaf temperature [K]
  real(r8), pointer :: ldew             ! depth of water on foliage [mm]

                       ! Vegetation dynamic parameters 
  real(r8), pointer :: fveg             ! fraction of vegetation cover
  real(r8), pointer :: sigf             ! fraction of veg cover, excluding snow-covered veg [-]
  real(r8), pointer :: green            ! leaf greenness
  real(r8), pointer :: lai              ! leaf area index
  real(r8), pointer :: sai              ! stem area index

  real(r8), pointer :: t10min              !annual minimum of 10-day running mean (K)
  real(r8), pointer :: lai_ind             !LAI per individual
  real(r8), pointer :: dphen               !phenology [0 to 1]
  real(r8), pointer :: leafon              !leafon days
  real(r8), pointer :: leafof              !leafoff days
  real(r8), pointer :: firelength          !fire season in days
  real(r8), pointer :: litterag            !above ground litter
  real(r8), pointer :: litterbg            !below ground litter
  real(r8), pointer :: cpool_fast          !fast carbon pool
  real(r8), pointer :: cpool_slow          !slow carbon pool
  real(r8), pointer :: k_fast_ave          !decomposition rate
  real(r8), pointer :: k_slow_ave          !decomposition rate
  real(r8), pointer :: litter_decom_ave    !decomposition rate
  real(r8), pointer :: fmicr               !microbial respiration (umol CO2 /m**2 /s)
  real(r8), pointer :: nind                !number of individuals (#/m**2)
  real(r8), pointer :: lm_ind              !individual leaf mass
  real(r8), pointer :: sm_ind              !individual sapwood mass
  real(r8), pointer :: hm_ind              !individual heartwood mass
  real(r8), pointer :: rm_ind              !individual root mass
  real(r8), pointer :: tmomin20            !20-yr running mean of tmomin
  real(r8), pointer :: agdd0               !growing dgree days above 0
  real(r8), pointer :: agdd                !growing dgree days above 5
  real(r8), pointer :: agddtw              !growing dgree days above twmax
  real(r8), pointer :: agdd20              !20-yr running mean of agdd
  real(r8), pointer :: t_mo                !30-day mean temperature of 2m
  real(r8), pointer :: t_mo_sum            !30-day accumulated temperature of 2m
  real(r8), pointer :: t_mo_min            !annual min of t_mo (Kelvin)
  real(r8), pointer :: crownarea           !area that each individual tree takes up (m^2)
  real(r8), pointer :: htop                !canopy top
  real(r8), pointer :: tsai                !one-sided stem area index, no burying by snow
  real(r8), pointer :: fpcgrid             !foliar projective cover on gridcell (fraction)
  real(r8), pointer :: bm_inc              !biomass increment
  real(r8), pointer :: afmicr              !annual microbial respiration
  real(r8), pointer :: annpsn              !annual photosynthesis (umol CO2 /m**2)
  real(r8), pointer :: annpsnpot           !annual potential photosynthesis (same units)
  real(r8), pointer :: tref10              !10-day averaged temperature at 2m
  real(r8), pointer :: tref_sum            !sum of temperature in current day
  real(r8), pointer :: t10(:)              !array ro record the 10 day temperature
  real(r8), pointer :: assimn10            !10-day averaged assimilation rate
  real(r8), pointer :: assimn_sum          !sum of assimn of current day
  real(r8), pointer :: an10(:)             !arry to record 10 day assimn
  real(r8), pointer :: anngpp              !annual gpp
  real(r8), pointer :: annfrmf             !annual frmf
  real(r8), pointer :: annfrms             !annual frms
  real(r8), pointer :: annfrmr             !annual frmr
  real(r8), pointer :: annfrg              !annual frg
  real(r8), pointer :: turnover_ind        !individual turnover biomass
! real(r8), pointer :: litter_ag           !above ground litter mass
! real(r8), pointer :: litter_bg           !below ground litter mass
  real(r8), pointer :: fpc_inc             !fpc increase
  real(r8), pointer :: ivt_x               !ivt
  real(r8), pointer :: pftpar(:)           !32 parameters of PFTs
  real(r8), pointer :: vegclass            !1.tree 2.shrub 3.grass 4.crop -1.others
  real(r8), pointer :: summergreen         !1. for summergreen; otherwise -1.
  real(r8), pointer :: raingreen           !1. for raingreen; otherwise -1.
  real(r8), pointer :: sla                 !sla
  real(r8), pointer :: lm_sapl             !leafmass
  real(r8), pointer :: sm_sapl             !sapwood mass
  real(r8), pointer :: hm_sapl             !heartwood mass
  real(r8), pointer :: rm_sapl             !rootmass
! real(r8), pointer :: ifpre               !-1=no PFT present;1=PFT present in this grid
#if(defined DyN)
   real(r8), pointer :: litter_leaf      ! leaf-derived litter for PFT on modelled area basis (gC/m2)
   real(r8), pointer :: litter_wood      ! heart&sapwood-derived litter for PFT on modelled area basis(gC/m2)
   real(r8), pointer :: litter_root      ! fine root-derived litter for PFT on modelled area basis(gC/m2)
   real(r8), pointer :: litter_repr      ! litter derived from allocation to reproduction for PFT on modelled

   real(r8), pointer :: litter_leaf_n ! leaf-derived N litter for PFT on modelled area basis (gN/m2)
   real(r8), pointer :: litter_wood_n ! heart&sapwood-derived N litter for PFT on modelled area basis(gN/m2)
   real(r8), pointer :: litter_root_n ! fine root-derived N litter for PFT on modelled area basis (gN/m2)
   real(r8), pointer :: litter_repr_n ! litter derived from allocation to reproduction N for PFT on modelled
                                         ! area basis (gN/m2)
   real(r8), pointer :: afcton_leaf   ! annual floating leaf C:N ratio
   real(r8), pointer :: afcton_root   ! annual floating root C:N ratio
   real(r8), pointer :: afcton_sap    ! annual floating sapwood C:N ratio
   real(r8), pointer :: lm_ind_n      ! individual leaf nitrogen mass
   real(r8), pointer :: sm_ind_n      ! individual sapwood nitrogen mass
   real(r8), pointer :: hm_ind_n      ! individual heartwood nitrogen mass
   real(r8), pointer :: rm_ind_n      ! individual root nitrogen mass
                                         ! gN/m2 veget'd area for each pft
   real(r8), pointer :: an_up         ! annual plant nitrogen uptake(gN/m2 vegt'd area)
   real(r8), pointer :: an_stress     ! annual nitrogen stress for production 
#endif
                       ! Radiation  related (albedoes)
  real(r8), pointer :: albg (:)         ! albedo, ground [-]
  real(r8), pointer :: albv (:)         ! albedo, vegetation [-]
  real(r8), pointer :: alb  (:)         ! averaged albedo [-]
  real(r8), pointer :: ssun (:)         ! sunlit canopy absorption for solar radiation (0-1)
  real(r8), pointer :: ssha (:)         ! shaded canopy absorption for solar radiation (0-1)
  real(r8), pointer :: thermk           ! canopy gap fraction for tir radiation
  real(r8), pointer :: extkb            ! (k, g(mu)/mu) direct solar extinction coefficient
  real(r8), pointer :: extkd            ! diffuse and scattered diffuse PAR extinction coefficient

                       ! Additional variables required by reginal model (WRF & RSM) 
  real(r8), pointer :: trad             ! radiative temperature of surface [K]
  real(r8), pointer :: tref             ! 2 m height air temperature [kelvin]
  real(r8), pointer :: qref             ! 2 m height air specific humidity
  real(r8), pointer :: rst              ! canopy stomatal resistance (s/m)

  real(r8), pointer :: emis             ! averaged bulk surface emissivity
  real(r8), pointer :: z0ma             ! effective roughness [m]
  real(r8), pointer :: zol              ! dimensionless height (z/L) used in Monin-Obukhov theory
  real(r8), pointer :: rib              ! bulk Richardson number in surface layer
  real(r8), pointer :: ustar            ! u* in similarity theory [m/s]
  real(r8), pointer :: qstar            ! q* in similarity theory [kg/kg]
  real(r8), pointer :: tstar            ! t* in similarity theory [K]
  real(r8), pointer :: fm               ! integral of profile function for momentum
  real(r8), pointer :: fh               ! integral of profile function for heat
  real(r8), pointer :: fq               ! integral of profile function for moisture

                       ! vegetation time invariant variables
  real(r8), pointer :: z0m_r            ! aerodynamic roughness length [m]
  real(r8), pointer :: chil             ! leaf angle distribution factor
  real(r8), pointer :: ref(:)           ! leaf reflectance (iw=iband, il=life and dead)
  real(r8), pointer :: tran(:)          ! leaf transmittance (iw=iband, il=life and dead)

! -----------------------------------------------------------------
! Local declaration
! -----------------------------------------------------------------

  real(r8) :: ssw                       ! water volumetric content of soil surface layer [m3/m3]
  real(r8) :: wt                        ! fraction of vegetation buried (covered) by snow [-]
  real(r8) :: z0m                       ! aerodynamic roughness length [m] zhq: 07/27/2010. z0m vary with htop
  integer lc,uc,lb,ub,jm                ! column/pft indices

  integer pfrom(1), p

  real(r8) wt_tmp(n_pft), vegwt

      jm = 15            ! nl_soil+abs(maxsnl)       

!      vegwt = sum(wt_patch(1:n_pft))
!      if(vegwt<1.0E-6) return

      wt_tmp = wt_patch_old-wt_patch

      pfrom = maxloc(wt_tmp)
      if(.not.(pfrom(1)>=1 .and. pfrom(1)<=n_pft)) then
          print *, 'pfrom', pfrom
          print *, wt_patch_old
          print *, wt_patch
          print *, 'error pfrom'
          return
      end if

      lc = 1
      uc = jm
      z         => fvar_col(lc:uc)                    ; lc = uc+1; uc = uc+jm !1_
      dz        => fvar_col(lc:uc)                    ; lc = uc+1; uc = uc+jm !2_
      tss       => fvar_col(lc:uc)                    ; lc = uc+1; uc = uc+jm !3_
      wliq      => fvar_col(lc:uc)                    ; lc = uc+1; uc = uc+jm !4_
      wice      => fvar_col(lc:uc)                    ; uc = uc+1             !5_
      tg        => fvar_col(uc)                       ; uc = uc+1             !1
      sag       => fvar_col(uc)                       ; uc = uc+1             !2
      scv       => fvar_col(uc)                       ; uc = uc+1             !3
      snowdp    => fvar_col(uc)                       ; uc = uc+1             !4
      fsno      => fvar_col(uc)                       ; uc = uc+1             !5
      coszen    => fvar_col(uc)                                               !6

      ssw = min(1.,1.e-3*wliq(6)/dz(6))

   do p = 1, n_pft

      ivt     = nint(ivt_r(p))
      ivt_old = nint(ivt_r_old(p))

      !if(ifpre(p).lt.0.) cycle       zhq: 07/27/2010. move below.   

      z0m_r     => fcon_pft(1,p)
      chil      => fcon_pft(16,p)        
      ref       => fcon_pft(17:20,p)        
      tran      => fcon_pft(21:24,p)        

      ub = 1
      tlsun     => fvar_pft(ub,p)                 ; ub = ub+1             !1
      tlsha     => fvar_pft(ub,p)                 ; ub = ub+1             !2
      ldew      => fvar_pft(ub,p)                 ; ub = ub+1             !3
      fveg      => fvar_pft(ub,p)                 ; ub = ub+1             !4
      sigf      => fvar_pft(ub,p)                 ; ub = ub+1             !5
      green     => fvar_pft(ub,p)                 ; ub = ub+1             !6
      lai       => fvar_pft(ub,p)                 ; ub = ub+1             !7
      sai       => fvar_pft(ub,p)                 ; lb = ub+1; ub = ub+4  !8
      albg      => fvar_pft(lb:ub,p)              ; lb = ub+1; ub = ub+4  !9-12
      albv      => fvar_pft(lb:ub,p)              ; lb = ub+1; ub = ub+4  !13-16
      alb       => fvar_pft(lb:ub,p)              ; lb = ub+1; ub = ub+4  !17-20
      ssun      => fvar_pft(lb:ub,p)              ; lb = ub+1; ub = ub+4  !21-24
      ssha      => fvar_pft(lb:ub,p)              ; ub = ub+1             !25-28
      thermk    => fvar_pft(ub,p)                 ; ub = ub+1             !29
      extkb     => fvar_pft(ub,p)                 ; ub = ub+1             !30
      extkd     => fvar_pft(ub,p)                 ; ub = ub+1             !31

                 ! Additional variables required by reginal model (WRF & RSM) 
      trad      => fvar_pft(ub,p)                 ; ub = ub+1             !32
      tref      => fvar_pft(ub,p)                 ; ub = ub+1             !33
      qref      => fvar_pft(ub,p)                 ; ub = ub+1             !34
      rst       => fvar_pft(ub,p)                 ; ub = ub+1             !35
      emis      => fvar_pft(ub,p)                 ; ub = ub+1             !36
      z0ma      => fvar_pft(ub,p)                 ; ub = ub+1             !37
      zol       => fvar_pft(ub,p)                 ; ub = ub+1             !38
      rib       => fvar_pft(ub,p)                 ; ub = ub+1             !39
      ustar     => fvar_pft(ub,p)                 ; ub = ub+1             !40
      qstar     => fvar_pft(ub,p)                 ; ub = ub+1             !41
      tstar     => fvar_pft(ub,p)                 ; ub = ub+1             !42
      fm        => fvar_pft(ub,p)                 ; ub = ub+1             !43
      fh        => fvar_pft(ub,p)                 ; ub = ub+1             !44
      fq        => fvar_pft(ub,p)                 ; ub = ub+1             !45

      t10min    => fvar_pft(ub,p)                 ; ub = ub+1             !1
      lai_ind   => fvar_pft(ub,p)                 ; ub = ub+1             !2
      dphen     => fvar_pft(ub,p)                 ; ub = ub+1             !3
      leafon    => fvar_pft(ub,p)                 ; ub = ub+1             !4
      leafof    => fvar_pft(ub,p)                 ; ub = ub+1             !5
      firelength=> fvar_pft(ub,p)                 ; ub = ub+1             !6
      litterag  => fvar_pft(ub,p)                 ; ub = ub+1             !7
      litterbg  => fvar_pft(ub,p)                 ; ub = ub+1             !8
      cpool_fast=> fvar_pft(ub,p)                 ; ub = ub+1             !9
      cpool_slow=> fvar_pft(ub,p)                 ; ub = ub+1             !10
      k_fast_ave=> fvar_pft(ub,p)                 ; ub = ub+1             !11
      k_slow_ave=> fvar_pft(ub,p)                 ; ub = ub+1             !12
      litter_decom_ave=> fvar_pft(ub,p)           ; ub = ub+1             !13
      fmicr     => fvar_pft(ub,p)                 ; ub = ub+1             !14
      nind      => fvar_pft(ub,p)                 ; ub = ub+1             !15
      lm_ind    => fvar_pft(ub,p)                 ; ub = ub+1             !16
      sm_ind    => fvar_pft(ub,p)                 ; ub = ub+1             !17
      hm_ind    => fvar_pft(ub,p)                 ; ub = ub+1             !18
      rm_ind    => fvar_pft(ub,p)                 ; ub = ub+1             !19
      tmomin20  => fvar_pft(ub,p)                 ; ub = ub+1             !20
      agdd0     => fvar_pft(ub,p)                 ; ub = ub+1             !21
      agdd      => fvar_pft(ub,p)                 ; ub = ub+1             !22
      agdd20    => fvar_pft(ub,p)                 ; ub = ub+1             !23
      t_mo_min  => fvar_pft(ub,p)                 ; ub = ub+1             !24
      crownarea => fvar_pft(ub,p)                 ; ub = ub+1             !25
      htop      => fvar_pft(ub,p)                 ; ub = ub+1             !26
      tsai      => fvar_pft(ub,p)                 ; ub = ub+1             !27
      fpcgrid   => fvar_pft(ub,p)                 ; ub = ub+1             !28
      bm_inc    => fvar_pft(ub,p)                 ; ub = ub+1             !29
      afmicr    => fvar_pft(ub,p)                 ; ub = ub+1             !30
      annpsn    => fvar_pft(ub,p)                 ; ub = ub+1             !31
      annpsnpot => fvar_pft(ub,p)                 ; ub = ub+1             !32
      tref10    => fvar_pft(ub,p)                 ; ub = ub+1             !33
      tref_sum  => fvar_pft(ub,p)                 ; ub = ub+1             !34
      t10       => fvar_pft(ub:ub+9,p)            ; ub = ub+10            !35-44
      assimn10  => fvar_pft(ub,p)                 ; ub = ub+1             !45
      assimn_sum=> fvar_pft(ub,p)                 ; ub = ub+1             !46
      an10      => fvar_pft(ub:ub+9,p)            ; ub = ub+10            !47-56
      turnover_ind => fvar_pft(ub,p)              ; ub = ub+1             !57
!     litter_ag => fvar_pft(ub,p)                 ; ub = ub+1             !58
!     litter_bg => fvar_pft(ub,p)                 ; ub = ub+1             !59
      fpc_inc   => fvar_pft(ub,p)                 ; ub = ub+1             !58
      ivt_x     => fvar_pft(ub,p)                 ; ub = ub+1             !59
      agddtw    => fvar_pft(ub,p)                 ; ub = ub+2             !60
   !  ifpre     => fvar_pft(ub,p)                                         !* 
      t_mo      => fvar_pft(ub,p)                 ; ub = ub+1             !62
      t_mo_sum  => fvar_pft(ub,p)                 ; ub = ub+1             !63
      anngpp    => fvar_pft(ub,p)                 ; ub = ub+1             !64
      annfrmf   => fvar_pft(ub,p)                 ; ub = ub+1             !65
      annfrms   => fvar_pft(ub,p)                 ; ub = ub+1             !66
      annfrmr   => fvar_pft(ub,p)                 ; ub = ub+1             !67
      annfrg    => fvar_pft(ub,p)                                         !68
#if(defined DyN)
      ub=ub+1
      litter_leaf       =>fvar_pft(ub,p)                  ; ub = ub+1             !1
      litter_wood       =>fvar_pft(ub,p)                  ; ub = ub+1             !2
      litter_root       =>fvar_pft(ub,p)                  ; ub = ub+1             !3
      litter_repr       =>fvar_pft(ub,p)                  ; ub = ub+1             !4
      litter_leaf_n     =>fvar_pft(ub,p)                  ; ub = ub+1             !5
      litter_wood_n     =>fvar_pft(ub,p)                  ; ub = ub+1             !6
      litter_root_n     =>fvar_pft(ub,p)                  ; ub = ub+1             !7
      litter_repr_n     =>fvar_pft(ub,p)                  ; ub = ub+1             !8
      afcton_leaf       =>fvar_pft(ub,p)                  ; ub = ub+1             !9
      afcton_sap        =>fvar_pft(ub,p)                  ; ub = ub+1             !10
      afcton_root       =>fvar_pft(ub,p)                  ; ub = ub+1             !11
      lm_ind_n          =>fvar_pft(ub,p)                  ; ub = ub+1             !12
      lm_ind_n          =>fvar_pft(ub,p)                  ; ub = ub+1             !13
      lm_ind_n          =>fvar_pft(ub,p)                  ; ub = ub+1             !14
      lm_ind_n          =>fvar_pft(ub,p)                  ; ub = ub+1             !15
      an_up             =>fvar_pft(ub,p)                  ; ub = ub+1             !16
      an_stress         =>fvar_pft(ub,p)                                          !17
#endif
!     print*,'lpjreset->',p,tlsun,tlsha 
! ======================================================================

! reset accumulated variables on pft level
      annpsn = 0.
      annpsnpot = 0.
      firelength = 0.
      bm_inc = 0.
      afmicr = 0.
      agdd0 = 0.
      agdd = 0.
      t_mo_min = 1.0E36
      anngpp  = 0.
      annfrmf = 0.
      annfrms = 0.
      annfrmr = 0.
      annfrg  = 0.
#if(defined DyN)
      an_up=0.
      an_stress=0.
#endif

! for new pfts established, zhq. 10/09/2009
      if(ifpre(p).lt.0.) cycle       

      if(ivt.ne.nint(ivt_x)) then
         print *, ivt, ivt_x
         stop 'lpjreset, ivt check failed'
      end if


      if(ifpre(p).gt.0. .and. ifpre_old(p).lt.0.) then

        !fvar_pft(1:125,p) = fvar_pft(1:125,pfrom(1))
         fvar_pft(65,p) = fvar_pft(65,pfrom(1))            ! tmomin20 added by zhq. dec24,08  
         fvar_pft(68,p) = fvar_pft(68,pfrom(1))            ! agdd20 added by zhq. dec24,08  
!!!xxx***fvar_pft(110,p) = fvar_pft(110,pfrom(1))          ! t_mo_sum added by zhq. dec24,08  
         fvar_pft(108,p) = fvar_pft(108,pfrom(1))          ! t_mo_sum added by zhq. dec24,08  

   ! update albedo for new pfts for the first timestep of next year, zhq. 11/11/2009

         lai = max(lai, 0.05)     
         sai = lai * 0.25

         z0m = max(z0m_r*htop, 0.01)
 
       ! call snowfraction (fveg,z0m,snowdp,wt,sigf,fsno)
         call snowfraction (itypwat,fveg,z0m,snowdp,scv,wt,sigf,fsno)

         call albland (itypwat,albsol,chil,ref,tran,&
                       fveg,green,lai,sai,coszen,wt,fsno,scv,sag,ssw,tg,&
                       alb,albg,albv,ssun,ssha,thermk,extkb,extkd)
      end if

!! reset accumulated variables on pft level
!      annpsn = 0.
!      annpsnpot = 0.
!      firelength = 0.
!      bm_inc = 0.
!      afmicr = 0.
!      agdd0 = 0.
!      agdd = 0.
!      t_mo_min = 1.0E36  
!#if(defined DyN)
!      an_up=0.
!      an_stress=0.
!#endif
   end do

! reset accumulated variables on column level
      prec365 = 0.  

end subroutine lpjreset
