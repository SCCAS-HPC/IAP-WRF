  subroutine lpjave(nyr_r       ,n_pft         ,npftpar       ,pftpar          ,&
                    vegclass    ,bm_inc        ,fpcgrid       ,agdd            ,&
                    t_mo_min    ,tmomin20      ,agdd20        ,afmicr          ,&
                    ifpre       ,lm_ind        ,hm_ind        ,sm_ind          ,&
                    rm_ind      ,litter_ag     ,litter_bg     ,turnover_ind    ,&
                    lai_ind     ,fpc_inc       ,crownarea     ,htop            ,&
                    nind        ,ivt_r         ,annpsn        ,annpsnpot       ,&
                    sla         ,agddtw        ,firelength    ,lm_sapl         ,&
                    sm_sapl     ,hm_sapl       ,rm_sapl       ,wt_patch        ,&
                    wt_column   ,prec365       ,cpool_fast    ,cpool_slow      ,&
                    vegc        ,litc          ,litc_ag       ,litc_bg         ,&
                    soic        ,soic_fast     ,soic_slow     ,&
                    litter2soil_ind, flitter2soil,&
                    litter2atmos_ind, flitter2atmos,&
                    leafc       ,woodc         ,rootc          )

!---------------------------------------------------------------------

   use precision
   use timemgr, only: dtime
   implicit none

! INTENT IN VARIAVLES:

   integer, INTENT(in) :: n_pft                      ! total patch number in this grid
   integer, INTENT(in) :: npftpar                    ! number of pft parameter
   real(r8),INTENT(in) :: nyr_r                      ! counting the model years
   real(r8),INTENT(in) :: pftpar(npftpar,n_pft)      ! 32 pft parameters
   real(r8),INTENT(in) :: vegclass(n_pft)            ! 1.tree 2.shrub 3.grass 4.crop -1.others
   real(r8),INTENT(in) :: agdd(n_pft)                ! accumulated growing degree days above 5
   real(r8),INTENT(in) :: t_mo_min(n_pft)            ! annual min of t_mo (Kelvin)
   real(r8),INTENT(in) :: annpsn(n_pft)              ! annual photosynthesis (umol CO2 /m**2)
   real(r8),INTENT(in) :: annpsnpot(n_pft)           ! annual potential photosynthesis (..)
   real(r8),INTENT(in) :: sla(n_pft)                 ! specific leaf area [m2 leaf g-1 carbon]
   real(r8),INTENT(in) :: agddtw(n_pft)              ! accumulated growing degree days above twmax
   real(r8),INTENT(in) :: firelength(n_pft)          ! fire season in days
   real(r8),INTENT(in) :: lm_sapl(n_pft)             ! ecophys const - leaf mass of sapling
   real(r8),INTENT(in) :: sm_sapl(n_pft)             ! ecophys const - stem mass of sapling
   real(r8),INTENT(in) :: hm_sapl(n_pft)             ! ecophys const - heartwood mass of sapling
   real(r8),INTENT(in) :: rm_sapl(n_pft)             ! ecophys const - root mass of saping
   real(r8),INTENT(in) :: cpool_fast(n_pft)          ! fast soil C pool (gC/m2 veget'd area)
   real(r8),INTENT(in) :: cpool_slow(n_pft)          ! slow soil C pool (gC/m2 veget'd area)
   real(r8),INTENT(in) :: prec365                    ! yearly running mean of precipitation [mm/s]
   real(r8),INTENT(in) :: wt_column                  ! relative weight of natural vegetation to grid

! INTENT INOUT VARIABLES:
   real(r8),INTENT(in) :: wt_patch(n_pft)         ! pft weight relative to grid cell
   real(r8),INTENT(in) :: tmomin20(n_pft)         ! 20-yr running mean of tmomin
   real(r8),INTENT(in) :: agdd20(n_pft)           ! 20-yr running mean of agdd
   real(r8),INTENT(in) :: bm_inc(n_pft)           ! biomass increment
   real(r8),INTENT(in) :: fpcgrid(n_pft)          ! foliar projective cover on gridcell
   real(r8),INTENT(in) :: afmicr(n_pft)           ! annual microbial respiration
   real(r8),INTENT(in) :: ifpre(n_pft)            ! whether pft present in this patch
   real(r8),INTENT(in) :: lm_ind(n_pft)           ! individual leaf mass
   real(r8),INTENT(in) :: sm_ind(n_pft)           ! individual sapwood mass
   real(r8),INTENT(in) :: hm_ind(n_pft)           ! individual heartwood mass
   real(r8),INTENT(in) :: rm_ind(n_pft)           ! individual root mass
   real(r8),INTENT(in) :: turnover_ind(n_pft)     ! individual turnover biomass
   real(r8),INTENT(in) :: litter_ag(n_pft)        ! above ground litter mass
   real(r8),INTENT(in) :: litter_bg(n_pft)        ! below ground litter mass
   real(r8),INTENT(in) :: lai_ind(n_pft)          ! max lai for individual
   real(r8),INTENT(in) :: fpc_inc(n_pft)          ! fpc increase
   real(r8),INTENT(in) :: crownarea(n_pft)        ! area each individual tree takes up (m^2)
   real(r8),INTENT(in) :: htop(n_pft)             ! canopy top 
   real(r8),INTENT(in) :: nind(n_pft)             ! population density over gridcell 
   real(r8),INTENT(in) :: ivt_r(n_pft)            ! land cover type 
   real(r8),INTENT(in) :: litter2soil_ind(n_pft)
   real(r8),INTENT(in) :: litter2atmos_ind(n_pft)
  
! INTENT OUT VARIABLES:
   real(r8),INTENT(out) :: vegc                   ! gridcell vegetaion biomass
   real(r8),INTENT(out) :: litc                   ! gridcell litter
   real(r8),INTENT(out) :: litc_ag                ! gridcell above ground litter
   real(r8),INTENT(out) :: litc_bg                ! gridcell below ground litter
   real(r8),INTENT(out) :: soic                   ! gridcell soil cpool
   real(r8),INTENT(out) :: soic_fast              ! gridcell fast soil cpool
   real(r8),INTENT(out) :: soic_slow              ! gridcell slow soil cpool(gC/m2 veget'd area)
   real(r8),INTENT(out) :: flitter2soil           ! 
   real(r8),INTENT(out) :: flitter2atmos          ! 
   real(r8),INTENT(out) :: leafc                  ! 
   real(r8),INTENT(out) :: woodc                  ! 
   real(r8),INTENT(out) :: rootc                  ! 

! LOCAL VARIABLES:
   integer  :: g,p,fp                  ! indices
   integer  :: ivt(n_pft)              ! integer form of ivt_r(n_pft)
   integer  :: nyr

  ! Transform ivt_r from real to integer
    ivt = nint(ivt_r)

    leafc = 0._r8
    woodc = 0._r8
    rootc = 0._r8
    vegc = 0._r8
    litc = 0._r8
    litc_ag = 0._r8
    litc_bg = 0._r8
    soic = 0._r8
    soic_fast = 0._r8
    soic_slow = 0._r8
    flitter2soil = 0._r8
    flitter2atmos= 0._r8

    do p = 1, n_pft
       leafc = leafc + lm_ind(p)*nind(p)
       woodc = woodc + (sm_ind(p) + hm_ind(p))*nind(p)
       rootc = rootc + rm_ind(p)*nind(p)
       litc_ag = litc_ag + litter_ag(p)
       litc_bg = litc_bg + litter_bg(p)
       soic_fast = soic_fast + cpool_fast(p)
       soic_slow = soic_slow + cpool_slow(p)
       flitter2soil = flitter2soil + litter2soil_ind(p)
       flitter2atmos = flitter2atmos + litter2atmos_ind(p)
    end do

    vegc = leafc + woodc + rootc

    litc = litc_ag + litc_bg

    soic = soic_fast + soic_slow

    flitter2soil = flitter2soil/dtime
    flitter2atmos = flitter2atmos/dtime

  end subroutine lpjave
