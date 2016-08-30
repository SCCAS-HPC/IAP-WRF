#include <define.h>
 
  subroutine LPJDRIVER
 
 !=======================================================================
 !
 ! LPJ MODEL DRIVER
 !
 ! Adopted from NCAR-CLM3.0 : Ming Chen 05/10/2008
 !
 !=======================================================================
 
  use precision
  use paramodel
  use spmd_decomp, only: cgmap, ggmap, gmask, numgrid_proc, gxmap, gymap
  use colm_varMod, only: lon_points, lat_points, numpatch, numcolumn, numgrid, &
                         fcon_col, fcon_pft, fvar_col, fvar_pft, fldv_dgvm, forc, &
                         wxy_patch, wxy_column, itypwat, nep_residual, lnep_adjust, &
                         fldv, idx_assim, idx_respc, idx_fmicr, fLitterSoil, fLitterAtmos
  use timemgr, only: idate_p, idate, dtime
  use landuse, only: landuse_cflux
  use debug, only: c_bug
  implicit none
 
 ! ----------------------------------------------------------------
 ! I. Time invariant model variables
 ! ----------------------------------------------------------------
 
   integer :: itypwat_c                        ! land water type
   real(r8), pointer :: dlat                   ! latitude in radians
   real(r8), pointer :: dlon                   ! longitude in radians

   real(r8), pointer :: albsol                 ! soil albedo for different coloured soils [-]
!  real(r8), pointer :: dz  (:)                ! interface depth [m]
!  real(r8), pointer :: wliq(:)                ! liquid water in layers [kg/m2]

   real(r8), pointer :: ivt(:)                 !land cover type  
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
   real(r8), pointer :: fmicr(:)               !microbial respiration (umol CO2 /m**2 /s)
   real(r8), pointer :: nind(:)                !number of individuals (*/m**2)
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
   real(r8), pointer :: crownarea(:)           !area that each individual tree takes up (m2)
   real(r8), pointer :: htop(:)                !canopy top
   real(r8), pointer :: tsai(:)                !one-sided stem area index, no burying by snow
   real(r8), pointer :: fpcgrid(:)             !foliar projective cover on gridcell (fraction)
   real(r8), pointer :: bm_inc(:)              !biomass increment
   real(r8), pointer :: afmicr(:)              !annual microbial respiration
   real(r8), pointer :: annpsn(:)              !annual photosynthesis (umol CO2 /m**2)
   real(r8), pointer :: annpsnpot(:)           !annual potential photosynthesis (same units)
   real(r8), pointer :: tref10(:)              !10-day averaged temperature at 2m
   real(r8), pointer :: tref_sum(:)            !sum of temperature in current day
   real(r8), pointer :: t10(:,:)               !array ro record the 10 day temperature
   real(r8), pointer :: assimn10(:)            !10-day averaged assimilation rate
   real(r8), pointer :: assimn_sum(:)          !sum of assimn of current day
   real(r8), pointer :: an10(:,:)              !arry to record 10 day assimn
   real(r8), pointer :: anngpp(:)              !annual gpp
   real(r8), pointer :: annfrmf(:)             !annual frmf
   real(r8), pointer :: annfrms(:)             !annual frms
   real(r8), pointer :: annfrmr(:)             !annual frmr
   real(r8), pointer :: annfrg(:)              !annual frg
   real(r8), pointer :: turnover_ind(:)        !individual turnover biomass
!  real(r8), pointer :: litter_ag(:)           !above ground litter mass
!  real(r8), pointer :: litter_bg(:)           !below ground litter mass
!  real(r8), pointer :: litter_leaf(:)         !leaf-derived litter for PFT on modelled area basis (gC/m2)
!  real(r8), pointer :: litter_wood(:)         !heart&sapwood-derived litter for PFT on modelled area basis(gC/m2)
!  real(r8), pointer :: litter_root(:)         !fine root-derived litter for PFT on modelled area basis(gC/m2)
!  real(r8), pointer :: litter_repr(:)         !litter derived from allocation to reproduction for PFT on modelled
   real(r8), pointer :: fpc_inc(:)             !fpc increase
   real(r8), pointer :: pftpar(:,:)            !32 parameters of PFTs
   real(r8), pointer :: vegclass(:)            !1.tree 2.shrub 3.grass 4.crop -1.others
   real(r8), pointer :: summergreen(:)         !1. for summergreen; otherwise -1.
   real(r8), pointer :: raingreen(:)           !1. for raingreen; otherwise -1.
   real(r8), pointer :: sla(:)                 !sla
   real(r8), pointer :: lm_sapl(:)             !leafmass
   real(r8), pointer :: sm_sapl(:)             !sapwood mass
   real(r8), pointer :: hm_sapl(:)             !heartwood mass
   real(r8), pointer :: rm_sapl(:)             !rootmass
   real(r8), pointer :: ifpre(:)               !-1=no PFT present;1=PFT present in this grid

   real(r8), pointer :: nday                   !counting the model days
   real(r8), pointer :: nyr                    !counting the model years
   real(r8), pointer :: prec365                !yearly running mean of precipitation(mm/s)

   real(r8), pointer :: wt_pft(:)              !weight of each patch in a grid
   real(r8), pointer :: wt_col                 !weight of each patch in a grid

#if(defined DyN)
   real(r8), pointer :: cton_soil(:)        ! soil C:N mass ratio

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

 ! -----------------------------------------------------------------
 ! V. Local declaration
 ! -----------------------------------------------------------------

   integer, parameter :: npftpar=32 ! 32 parameter of each PFT 
 
   integer i,j,lb,ub,lc,uc,jm,n           ! loop/array indices
 
   integer n_pft 
   integer g, c, p
   integer p1, p2 

   logical newyear
 
   real(r8) :: ifpre_old(17), ivt_old(17), wt_pft_old(17)

   real(r8) :: leafc, woodc, rootc, vegc
   real(r8) :: litc, litc_ag, litc_bg
   real(r8) :: soic, soic_fast, soic_slow
   real(r8) :: fveg2litter, flitter2soil, flitter2atmos
   real(r8) :: gpp, npp, nep, nbp
   real(r8) :: ra, rh, ffirec
   real(r8), pointer :: litter2soil_ind(:), litter2atmos_ind(:)

   real(r8) :: lucflux(lon_points,lat_points)

 !-------------------------------------------------------------------
 ! variables for histDGVM
 !-------------------------------------------------------------------
   real(r8) afirec, afiref, avegc, aestabc, anpp, amrh, alitc_ag, alitc_bg, asoic_fast, asoic_slow
   real(r8) npp_ind(17)
   real(r8) gpp_ind(17)
   real(r8) frmf_ind(17)
   real(r8) frms_ind(17)
   real(r8) frmr_ind(17)
   real(r8) frg_ind(17)
#if(defined DyN) 
   real(r8) an_up_total, an_stress_total, avegn, alitn_ag, alitn_bg, asoin
#endif
 !---------------------------------------------------------------------------------------
 !LPJ is called on the gridcell level, because there are communicates between
 !patches in the anual processes.(light competetion, for example). Here we assume
 !that there are 21 patches in each gridcell. If one patch do not exist, then the
 !weight of the patch is set to zero. This part can be improved in further development.
 !---------------------------------------------------------------------------------------

    if(idate_p(1).ne.idate(1)) then
       newyear = .true.
       lnep_adjust = .true.
    else
       newyear = .false.
    endif

    fldv_dgvm(:,:) = 0.

    nep_residual(:) = 0.

    call landuse_cflux(idate(1),lucflux)

  ! g(C)/m2/s -> kg(C)/m2/s
    lucflux(:,:) = lucflux(:,:)*1.e-3_r8

    p1 = 1
    p2 = 0
 
    DO c = 1, numcolumn

       itypwat_c =  nint(fcon_col(3,c))                                      !3
       albsol    => fcon_col(4,c)                                            !4

       uc = 82                      
       nday      => fvar_col(uc,c)                   ; uc = uc+1             !1
       nyr       => fvar_col(uc,c)                   ; uc = uc+1             !2
       prec365   => fvar_col(uc,c)                                           !3

#if(defined DyN)
       uc = uc+1
       soil_no3  => fvar_col(uc,c)                   ; uc = uc+1             !1
       soil_no2  => fvar_col(uc,c)                   ; uc = uc+1             !2
       soil_no   => fvar_col(uc,c)                   ; uc = uc+1             !3
       soil_n2o  => fvar_col(uc,c)                   ; uc = uc+1             !4
       soil_n2   => fvar_col(uc,c)                   ; uc = uc+1             !5
       soil_nh4  => fvar_col(uc,c)                                           !6
#endif

       if(itypwat_c.ne.itypwat(c)) then
          write(6,*) 'fatal error on itypwat checking [LPJDRIVER]'
          call abort
       end if

       if(itypwat_c==0)      n_pft = 16 + 1  ! natural vegetation + bare soil
       if(itypwat_c==1)      n_pft = 1       ! urban and built-up
       if(itypwat_c==2)      n_pft = 1       ! wetland
       if(itypwat_c==3)      n_pft = 1       ! land ice
       if(itypwat_c==4)      n_pft = 1       ! river or deep lake
       if(itypwat_c==99)then                 ! ocean
          write(6,*) 'ocean column', c
          call abort
       endif

       p2 = p2 + n_pft
 
       ub = 34+1

       pftpar            => fcon_pft(ub:ub+npftpar-1,p1:p2)    ; ub = ub+npftpar       !1-32
       vegclass          => fcon_pft(ub,p1:p2)                 ; ub = ub+1             !33
       summergreen       => fcon_pft(ub,p1:p2)                 ; ub = ub+1             !34
       raingreen         => fcon_pft(ub,p1:p2)                 ; ub = ub+1             !35
       sla               => fcon_pft(ub,p1:p2)                 ; ub = ub+1             !36
       lm_sapl           => fcon_pft(ub,p1:p2)                 ; ub = ub+1             !37
       sm_sapl           => fcon_pft(ub,p1:p2)                 ; ub = ub+1             !38
       hm_sapl           => fcon_pft(ub,p1:p2)                 ; ub = ub+1             !39
       rm_sapl           => fcon_pft(ub,p1:p2)                                         !40

#if(defined DyN)
       ub = ub + 1
       cton_soil         => fcon_pft(ub,p1:p2)                                         !75
#endif
 
       ub = 45+1
 
       t10min            => fvar_pft(ub,p1:p2)                 ; ub = ub+1             !1
       lai_ind           => fvar_pft(ub,p1:p2)                 ; ub = ub+1             !2
       dphen             => fvar_pft(ub,p1:p2)                 ; ub = ub+1             !3
       leafon            => fvar_pft(ub,p1:p2)                 ; ub = ub+1             !4
       leafof            => fvar_pft(ub,p1:p2)                 ; ub = ub+1             !5
       firelength        => fvar_pft(ub,p1:p2)                 ; ub = ub+1             !6
       litterag          => fvar_pft(ub,p1:p2)                 ; ub = ub+1             !7
       litterbg          => fvar_pft(ub,p1:p2)                 ; ub = ub+1             !8
       cpool_fast        => fvar_pft(ub,p1:p2)                 ; ub = ub+1             !9
       cpool_slow        => fvar_pft(ub,p1:p2)                 ; ub = ub+1             !10
       k_fast_ave        => fvar_pft(ub,p1:p2)                 ; ub = ub+1             !11
       k_slow_ave        => fvar_pft(ub,p1:p2)                 ; ub = ub+1             !12
       litter_decom_ave  => fvar_pft(ub,p1:p2)                 ; ub = ub+1             !13
       fmicr             => fvar_pft(ub,p1:p2)                 ; ub = ub+1             !14
       nind              => fvar_pft(ub,p1:p2)                 ; ub = ub+1             !15
       lm_ind            => fvar_pft(ub,p1:p2)                 ; ub = ub+1             !16
       sm_ind            => fvar_pft(ub,p1:p2)                 ; ub = ub+1             !17
       hm_ind            => fvar_pft(ub,p1:p2)                 ; ub = ub+1             !18
       rm_ind            => fvar_pft(ub,p1:p2)                 ; ub = ub+1             !19
       tmomin20          => fvar_pft(ub,p1:p2)                 ; ub = ub+1             !20
       agdd0             => fvar_pft(ub,p1:p2)                 ; ub = ub+1             !21
       agdd              => fvar_pft(ub,p1:p2)                 ; ub = ub+1             !22
       agdd20            => fvar_pft(ub,p1:p2)                 ; ub = ub+1             !23
       t_mo_min          => fvar_pft(ub,p1:p2)                 ; ub = ub+1             !24
       crownarea         => fvar_pft(ub,p1:p2)                 ; ub = ub+1             !25
       htop              => fvar_pft(ub,p1:p2)                 ; ub = ub+1             !26
       tsai              => fvar_pft(ub,p1:p2)                 ; ub = ub+1             !27
       fpcgrid           => fvar_pft(ub,p1:p2)                 ; ub = ub+1             !28
       bm_inc            => fvar_pft(ub,p1:p2)                 ; ub = ub+1             !29
       afmicr            => fvar_pft(ub,p1:p2)                 ; ub = ub+1             !30
       annpsn            => fvar_pft(ub,p1:p2)                 ; ub = ub+1             !31
       annpsnpot         => fvar_pft(ub,p1:p2)                 ; ub = ub+1             !32
       tref10            => fvar_pft(ub,p1:p2)                 ; ub = ub+1             !33
       tref_sum          => fvar_pft(ub,p1:p2)                 ; ub = ub+1             !34
       t10               => fvar_pft(ub:ub+9,p1:p2)            ; ub = ub+10            !35-44
       assimn10          => fvar_pft(ub,p1:p2)                 ; ub = ub+1             !45
       assimn_sum        => fvar_pft(ub,p1:p2)                 ; ub = ub+1             !46
       an10              => fvar_pft(ub:ub+9,p1:p2)            ; ub = ub+10            !47-56
       turnover_ind      => fvar_pft(ub,p1:p2)                 ; ub = ub+1             !57
       fpc_inc           => fvar_pft(ub,p1:p2)                 ; ub = ub+1             !58
       ivt               => fvar_pft(ub,p1:p2)                 ; ub = ub+1             !59
       agddtw            => fvar_pft(ub,p1:p2)                 ; ub = ub+1             !60
       ifpre             => fvar_pft(ub,p1:p2)                 ; ub = ub+1             !61 
       t_mo              => fvar_pft(ub,p1:p2)                 ; ub = ub+1             !62
       t_mo_sum          => fvar_pft(ub,p1:p2)                 ; ub = ub+1             !63  
       anngpp            => fvar_pft(ub,p1:p2)                 ; ub = ub+1             !64
       annfrmf           => fvar_pft(ub,p1:p2)                 ; ub = ub+1             !65
       annfrms           => fvar_pft(ub,p1:p2)                 ; ub = ub+1             !66
       annfrmr           => fvar_pft(ub,p1:p2)                 ; ub = ub+1             !67
       annfrg            => fvar_pft(ub,p1:p2)                                         !68

#if(defined DyN)
       ub=ub+1
       litter_leaf       => fvar_pft(ub,p1:p2)                 ; ub = ub+1             !1
       litter_wood       => fvar_pft(ub,p1:p2)                 ; ub = ub+1             !2
       litter_root       => fvar_pft(ub,p1:p2)                 ; ub = ub+1             !3
       litter_repr       => fvar_pft(ub,p1:p2)                 ; ub = ub+1             !4
       litter_leaf_n     => fvar_pft(ub,p1:p2)                 ; ub = ub+1             !5
       litter_wood_n     => fvar_pft(ub,p1:p2)                 ; ub = ub+1             !6
       litter_root_n     => fvar_pft(ub,p1:p2)                 ; ub = ub+1             !7
       litter_repr_n     => fvar_pft(ub,p1:p2)                 ; ub = ub+1             !8
       afcton_leaf       => fvar_pft(ub,p1:p2)                 ; ub = ub+1             !9
       afcton_sap        => fvar_pft(ub,p1:p2)                 ; ub = ub+1             !10
       afcton_root       => fvar_pft(ub,p1:p2)                 ; ub = ub+1             !11
       lm_ind_n          => fvar_pft(ub,p1:p2)                 ; ub = ub+1             !12
       sm_ind_n          => fvar_pft(ub,p1:p2)                 ; ub = ub+1             !13
       hm_ind_n          => fvar_pft(ub,p1:p2)                 ; ub = ub+1             !14
       rm_ind_n          => fvar_pft(ub,p1:p2)                 ; ub = ub+1             !15
       an_up             => fvar_pft(ub,p1:p2)                 ; ub = ub+1             !16
       an_stress         => fvar_pft(ub,p1:p2)                                         !17
#endif

       litter2soil_ind   => fLitterSoil(p1:p2)
       litter2atmos_ind  => fLitterAtmos(p1:p2)

       wt_pft            => wxy_patch(p1:p2)
       wt_col            => wxy_column(c)
 
       IF (itypwat_c==0) THEN

          g = cgmap(c)
          i = gxmap(g)
          j = gymap(g)

          leafc        = 0._r8
          woodc        = 0._r8
          rootc        = 0._r8
          vegc         = 0._r8
          litc         = 0._r8
          litc_ag      = 0._r8
          litc_bg      = 0._r8
          soic         = 0._r8
          soic_fast    = 0._r8
          soic_slow    = 0._r8
          fveg2litter  = 0._r8
          flitter2soil = 0._r8
          flitter2atmos= 0._r8
          gpp          = 0._r8
          npp          = 0._r8
          nep          = 0._r8
          nbp          = 0._r8
          ra           = 0._r8
          rh           = 0._r8
          ffirec       = 0._r8

          CALL lpjave(nyr         ,n_pft       ,npftpar       ,pftpar           ,&
                      vegclass    ,bm_inc      ,fpcgrid       ,agdd             ,&
                      t_mo_min    ,tmomin20    ,agdd20        ,afmicr           ,&
                      ifpre       ,lm_ind      ,hm_ind        ,sm_ind           ,&
                      rm_ind      ,litterag    ,litterbg      ,turnover_ind     ,&
                      lai_ind     ,fpc_inc     ,crownarea     ,htop             ,&
                      nind        ,ivt         ,annpsn        ,annpsnpot        ,&
                      sla         ,agddtw      ,firelength    ,lm_sapl          ,&
                      sm_sapl     ,hm_sapl     ,rm_sapl       ,wt_pft           ,&
                      wt_col      ,prec365     ,cpool_fast    ,cpool_slow       ,&
                      vegc        ,litc        ,litc_ag       ,litc_bg          ,&
                      soic        ,soic_fast   ,soic_slow     ,&
                      litter2soil_ind, flitter2soil           ,&
                      litter2atmos_ind,flitter2atmos          ,&
                      leafc       ,woodc       ,rootc          )

        ! Units of kgC/m2

          leafc     = leafc*wxy_column(c)*0.001
          woodc     = woodc*wxy_column(c)*0.001
          rootc     = rootc*wxy_column(c)*0.001
          vegc      = vegc*wxy_column(c)*0.001
          litc      = litc*wxy_column(c)*0.001
          litc_ag   = litc_ag*wxy_column(c)*0.001
          litc_bg   = litc_bg*wxy_column(c)*0.001
          soic      = soic*wxy_column(c)*0.001
          soic_fast = soic_fast*wxy_column(c)*0.001
          soic_slow = soic_slow*wxy_column(c)*0.001

        ! Units of kgC/m2/s

          flitter2soil  = flitter2soil*wxy_column(c)*0.001
          flitter2atmos = flitter2atmos*wxy_column(c)*0.001

        ! GPP, NPP, NEP, NBP have already been weighted.

          gpp = fldv(idx_assim,g)*0.012
          ra  = fldv(idx_respc,g)*0.012
          rh  = fldv(idx_fmicr,g)*0.012
          npp = gpp-ra
          nep = gpp-ra-rh
          nbp = gpp-ra-rh-lucflux(i,j)

          ivt_old = ivt
          ifpre_old = ifpre
          wt_pft_old = wt_pft

#ifndef VEGDATA
          IF (newyear) THEN
   
             gpp_ind  = anngpp
             frmf_ind = annfrmf
             frms_ind = annfrms
             frmr_ind = annfrmr
             frg_ind  = annfrg
   
             CALL LPJ(nyr         ,n_pft       ,npftpar       ,pftpar        ,&
                      vegclass    ,bm_inc      ,fpcgrid       ,agdd          ,&
                      t_mo_min    ,tmomin20    ,agdd20        ,afmicr        ,&
                      ifpre       ,lm_ind      ,hm_ind        ,sm_ind        ,&
                      rm_ind      ,litterag    ,litterbg      ,turnover_ind  ,&
                      lai_ind     ,fpc_inc     ,crownarea     ,htop          ,&
                      nind        ,ivt         ,annpsn        ,annpsnpot     ,&
                      sla         ,agddtw      ,firelength    ,lm_sapl       ,&
                      sm_sapl     ,hm_sapl     ,rm_sapl       ,wt_pft        ,&
                      wt_col      ,prec365     ,afirec        ,afiref        ,&
                      avegc       ,anpp        ,amrh          ,alitc_ag      ,&
                      alitc_bg    ,asoic_fast  ,asoic_slow    ,cpool_fast    ,&
                      cpool_slow  ,aestabc     ,npp_ind)
   
             CALL lpjreset(itypwat_c,albsol,nfcon_pft,nfvar_col,nfvar_pft,n_pft,&
                           fcon_pft(1:nfcon_pft,p1:p2),fvar_col(1:nfvar_col,c),&
                           fvar_pft(1:nfvar_pft,p1:p2),ivt,&
                           ivt_old,ifpre,ifpre_old,wt_pft,wt_pft_old,prec365)
    
             nyr = nyr + 1.

           ! gC/m2/yr => molC/m2/s
             nep_residual(g) = (afirec-aestabc)*wxy_column(c)/(12.*dtime)

             nep = nep - (afirec-aestabc)*wxy_column(c)*0.001/dtime

             nbp = nbp - (afirec-aestabc)*wxy_column(c)*0.001/dtime

           ! gC/m2/yr => kgC/m2/s
             ffirec = afirec*wxy_column(c)*0.001/dtime
             fveg2litter = ((alitc_ag+alitc_bg)*wxy_column(c)*0.001-litc_ag-litc_bg)/dtime

             litc_ag = alitc_ag*wxy_column(c)*0.001
             litc_bg = alitc_bg*wxy_column(c)*0.001

          END IF
#endif

        !*Monthly variables (weighted)

          uc = 1

          fldv_dgvm(uc,g)    = leafc                      ; uc = uc+1                     !1
          fldv_dgvm(uc,g)    = woodc                      ; uc = uc+1                     !2
          fldv_dgvm(uc,g)    = rootc                      ; uc = uc+1                     !3
          fldv_dgvm(uc,g)    = vegc                       ; uc = uc+1                     !4
          fldv_dgvm(uc,g)    = litc_ag                    ; uc = uc+1                     !5
          fldv_dgvm(uc,g)    = litc_bg                    ; uc = uc+1                     !6
          fldv_dgvm(uc,g)    = litc                       ; uc = uc+1                     !7
          fldv_dgvm(uc,g)    = soic_fast                  ; uc = uc+1                     !8
          fldv_dgvm(uc,g)    = soic_slow                  ; uc = uc+1                     !9
          fldv_dgvm(uc,g)    = soic                       ; uc = uc+1                     !10
          fldv_dgvm(uc,g)    = fveg2litter                ; uc = uc+1                     !11
          fldv_dgvm(uc,g)    = flitter2soil               ; uc = uc+1                     !12
          fldv_dgvm(uc,g)    = flitter2atmos              ; uc = uc+1                     !13
          fldv_dgvm(uc,g)    = gpp                        ; uc = uc+1                     !14
          fldv_dgvm(uc,g)    = npp                        ; uc = uc+1                     !15
          fldv_dgvm(uc,g)    = nep                        ; uc = uc+1                     !16
          fldv_dgvm(uc,g)    = nbp                        ; uc = uc+1                     !17
          fldv_dgvm(uc,g)    = ra                         ; uc = uc+1                     !18
          fldv_dgvm(uc,g)    = rh                         ; uc = uc+1                     !19
          fldv_dgvm(uc,g)    = ffirec                                                     !20

#ifdef MYBUG
if(p_master) then
   if(newyear) write(6,*) 'newyear'
   if(c_bug.eq.c) write(6,*) 'wt_patc2', wt_pft_old(1:numpft_nat)
endif
#endif

          lc = uc+1                                       ; uc = uc+numpft_nat
          fldv_dgvm(lc:uc,g) = wt_pft_old(1:numpft_nat)                                   !21-34
   
        !*Yearly variables (unweighted)

          lc = nfldv_dgvm_accum +1                        ; uc = nfldv_dgvm_accum + numpft_nat
   
          fldv_dgvm(lc:uc,g) = fpcgrid(1:numpft_nat)      ; uc = uc+1                     !35-48
          fldv_dgvm(uc,g)    = fpcgrid(17)                ; uc = uc+1                     !49
          fldv_dgvm(uc,g)    = afirec                     ; uc = uc+1                     !50
          fldv_dgvm(uc,g)    = afiref                     ; uc = uc+1                     !51
          fldv_dgvm(uc,g)    = avegc                      ; uc = uc+1                     !52
          fldv_dgvm(uc,g)    = aestabc                    ; uc = uc+1                     !53
          fldv_dgvm(uc,g)    = anpp                       ; uc = uc+1                     !54
          fldv_dgvm(uc,g)    = amrh                       ; uc = uc+1                     !55
          fldv_dgvm(uc,g)    = alitc_ag                   ; uc = uc+1                     !56
          fldv_dgvm(uc,g)    = alitc_bg                   ; uc = uc+1                     !57
          fldv_dgvm(uc,g)    = asoic_fast                 ; uc = uc+1                     !58
          fldv_dgvm(uc,g)    = asoic_slow                 ; lc = uc+1; uc = uc+numpft_nat !59
          fldv_dgvm(lc:uc,g) = npp_ind(1:numpft_nat)      ; lc = uc+1; uc = uc+numpft_nat !60-73
          fldv_dgvm(lc:uc,g) = lm_ind(1:numpft_nat)       ; lc = uc+1; uc = uc+numpft_nat !74-87
          fldv_dgvm(lc:uc,g) = sm_ind(1:numpft_nat)       ; lc = uc+1; uc = uc+numpft_nat !88-101
          fldv_dgvm(lc:uc,g) = hm_ind(1:numpft_nat)       ; lc = uc+1; uc = uc+numpft_nat !102-115
          fldv_dgvm(lc:uc,g) = rm_ind(1:numpft_nat)       ; lc = uc+1; uc = uc+numpft_nat !116-129
          fldv_dgvm(lc:uc,g) = crownarea(1:numpft_nat)    ; lc = uc+1; uc = uc+numpft_nat !130-143
          fldv_dgvm(lc:uc,g) = htop(1:numpft_nat)         ; lc = uc+1; uc = uc+numpft_nat !144-157
          fldv_dgvm(lc:uc,g) = nind(1:numpft_nat)         ; lc = uc+1; uc = uc+numpft_nat !158-171
          fldv_dgvm(lc:uc,g) = lai_ind(1:numpft_nat)      ; lc = uc+1; uc = uc+numpft_nat !172-185
          fldv_dgvm(lc:uc,g) = gpp_ind(1:numpft_nat)      ; lc = uc+1; uc = uc+numpft_nat !186-199
          fldv_dgvm(lc:uc,g) = frmf_ind(1:numpft_nat)     ; lc = uc+1; uc = uc+numpft_nat !200-213
          fldv_dgvm(lc:uc,g) = frms_ind(1:numpft_nat)     ; lc = uc+1; uc = uc+numpft_nat !214-227
          fldv_dgvm(lc:uc,g) = frmr_ind(1:numpft_nat)     ; lc = uc+1; uc = uc+numpft_nat !228-241
          fldv_dgvm(lc:uc,g) = frg_ind(1:numpft_nat)      ;                               !242-255

#ifdef MYBUG
          do p=1,numpft_nat
             if(fpcgrid(p)>0.)then
                print*,'bm',p,npp_ind(p)
                print*,'lai',p,lai_ind(p)    
                print*,'lm',p,lm_ind(p)
                print*,'sm',p,sm_ind(p)
                print*,'hm',p,hm_ind(p)
                print*,'rm',p,rm_ind(p)
             endif
          enddo
#endif
   
#ifdef DyN
          lc = uc+1 ; uc = uc+numpft_nat
          fldv_dgvm(lc:uc,g) = afcton_leaf(1:numpft_nat)  ; lc = uc+1; uc = uc+numpft_nat !1-14
          fldv_dgvm(lc:uc,g) = afcton_sap(1:numpft_nat)   ; lc = uc+1; uc = uc+numpft_nat !15-28
          fldv_dgvm(lc:uc,g) = afcton_root(1:numpft_nat)  ; uc = uc+1                     !29-42
          fldv_dgvm(uc,g) = an_up_total                   ; uc = uc+1                     !43
          fldv_dgvm(uc,g) = an_stress_total               ; uc = uc+1                     !44
          fldv_dgvm(uc,g) = avegn                         ; uc = uc+1                     !45
          fldv_dgvm(uc,g) = alitn_ag                      ; uc = uc+1                     !46
          fldv_dgvm(uc,g) = alitn_bg                      ; uc = uc+1                     !47
          fldv_dgvm(uc,g) = asoin                         ; uc = uc+1                     !48
          fldv_dgvm(uc,g) = soil_no3                      ; uc = uc+1                     !49
          fldv_dgvm(uc,g) = soil_nh4                                                      !50
   
#ifdef MYBUG
          print*, 'npp_ind',npp_ind
          print*, 'inorganic N',soil_no3,soil_nh4
          print*, 'nitrogen',an_up_total,an_stress_total,avegn,alitn_ag,alitn_bg,asoin,soil_no3,soil_nh4
#endif

#endif
       END IF

       p1 = p2 + 1
 
    END DO  !column
 
 end subroutine LPJDRIVER
