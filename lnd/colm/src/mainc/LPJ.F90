  subroutine LPJ(nyr_r       ,n_pft         ,npftpar       ,pftpar        ,&
                 vegclass    ,bm_inc        ,fpcgrid       ,agdd          ,&
                 t_mo_min    ,tmomin20      ,agdd20        ,afmicr        ,&
                 ifpre       ,lm_ind        ,hm_ind        ,sm_ind        ,&
                 rm_ind      ,litter_ag     ,litter_bg     ,turnover_ind  ,&
                 lai_ind     ,fpc_inc       ,crownarea     ,htop          ,&
                 nind        ,ivt_r         ,annpsn        ,annpsnpot     ,&
                 sla         ,agddtw        ,firelength    ,lm_sapl       ,&
                 sm_sapl     ,hm_sapl       ,rm_sapl       ,wt_patch      ,&
                 wt_column   ,prec365       ,afirec_gcell  ,afiref_gcell  ,&
                 avegc_gcell ,anpp_gcell    ,amrh_gcell    ,alitc_ag      ,&
                 alitc_bg    ,asoic_fast    ,asoic_slow    ,cpool_fast    ,&
                 cpool_slow  ,aestabc_gcell ,npp_ind)

!---------------------------------------------------------------------
! Drives the annual portion of lpj, called once per year from LPJDRIVER
! Shrub is included in adapted from Xiaodong Zeng's code.
! subroutine with * is based on new consideration.
!
!   |- BuildNatVegFilter
!   |
!   |- Reproduction
!   |
!   |- Turnover
!   |
!   |- Kill
!   |
!   |- BuildNatVegFilter
!   |
!   |- * Allocation
!   |
!   |- Light       
!   |
!   |- Mortality
!   |
!   |- Fire
!   |
!   |- * Establishment
!   |
!   |- BuildNatVegFilter
!   |
!   |- Light        
!   |
!   |- * ResetBareFPC
!---------------------------------------------------------------------

   use precision
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
   real(r8),INTENT(inout) :: wt_patch(n_pft)         ! pft weight relative to grid cell
   real(r8),INTENT(inout) :: tmomin20(n_pft)         ! 20-yr running mean of tmomin
   real(r8),INTENT(inout) :: agdd20(n_pft)           ! 20-yr running mean of agdd
   real(r8),INTENT(inout) :: bm_inc(n_pft)           ! biomass increment
   real(r8),INTENT(inout) :: fpcgrid(n_pft)          ! foliar projective cover on gridcell
   real(r8),INTENT(inout) :: afmicr(n_pft)           ! annual microbial respiration
   real(r8),INTENT(inout) :: ifpre(n_pft)            ! whether pft present in this patch
   real(r8),INTENT(inout) :: lm_ind(n_pft)           ! individual leaf mass
   real(r8),INTENT(inout) :: sm_ind(n_pft)           ! individual sapwood mass
   real(r8),INTENT(inout) :: hm_ind(n_pft)           ! individual heartwood mass
   real(r8),INTENT(inout) :: rm_ind(n_pft)           ! individual root mass
   real(r8),INTENT(inout) :: turnover_ind(n_pft)     ! individual turnover biomass
   real(r8),INTENT(inout) :: litter_ag(n_pft)        ! above ground litter mass
   real(r8),INTENT(inout) :: litter_bg(n_pft)        ! below ground litter mass
   real(r8),INTENT(inout) :: lai_ind(n_pft)          ! max lai for individual
   real(r8),INTENT(inout) :: fpc_inc(n_pft)          ! fpc increase
   real(r8),INTENT(inout) :: crownarea(n_pft)        ! area each individual tree takes up (m^2)
   real(r8),INTENT(inout) :: htop(n_pft)             ! canopy top 
   real(r8),INTENT(inout) :: nind(n_pft)             ! population density over gridcell 
   real(r8),INTENT(inout) :: ivt_r(n_pft)            ! land cover type 
  
! INTENT OUT VARIABLES:
   real(r8),INTENT(out) :: afiref_gcell            ! fraction of gridcell affected by fire
   real(r8),INTENT(out) :: afirec_gcell            ! gridcell C flux to atmosphere from burning
   real(r8),INTENT(out) :: anpp_gcell              ! gridcell annual net primary production
   real(r8),INTENT(out) :: amrh_gcell              ! gridcell microbial respiration
   real(r8),INTENT(out) :: avegc_gcell             ! gridcell vegetaion biomass
   real(r8),INTENT(out) :: aestabc_gcell           ! gridcell established vegetation biomass
   real(r8),INTENT(out) :: alitc_ag                ! gridcell above ground litter
   real(r8),INTENT(out) :: alitc_bg                ! gridcell below ground litter
   real(r8),INTENT(out) :: asoic_fast              ! gridcell fast soil cpool
   real(r8),INTENT(out) :: asoic_slow              ! gridcell slow soil cpool(gC/m2 veget'd area)
   real(r8),INTENT(out) :: npp_ind(n_pft)          ! gridcell slow soil cpool(gC/m2 veget'd area)

! LOCAL VARIABLES:
   integer  :: g,p,fp                  ! indices
   integer  :: num_natvegp             ! number of naturally-vegetated pfts in filter
   integer  :: filter_natvegp(n_pft)   ! filter for naturally-vegetated pfts
   real(r8) :: afirefrac(n_pft)        ! for history write
   real(r8) :: acfluxfire(n_pft)       ! for history write
   real(r8) :: prec                    ! for history write
   real(r4) :: fpcout(n_pft)           ! for history write  ! add by zhq 
   integer  :: ivt(n_pft)              ! integer form of ivt_r(n_pft)
   integer  :: nyr

! Avoid silly gfortran intrinsic subroutine name collision problem
#ifdef __GFORTRAN__
   EXTERNAL Kill
#endif

  ! Transform ivt_r from real to integer
    ivt = nint(ivt_r)
    nyr = nint(nyr_r)

#ifdef MYBUG
    print *, 'ivt-1', ivt
    print *, 'ifpre-1', ifpre
    print *, 'bm_inc-1', bm_inc
    print *, 'fpcgrid-1', fpcgrid
#endif

    do p = 1,n_pft
       if (nyr == 1) then
          tmomin20(p) = t_mo_min(p)
          agdd20(p) = agdd(p)
       end if
       tmomin20(p) = (19.0 * tmomin20(p) + t_mo_min(p)) / 20.0
       agdd20(p)   = (19.0 * agdd20(p)   + agdd(p)    ) / 20.0
    end do

   ! Determine grid values of npp and microbial respiration

    anpp_gcell = 0._r8
    amrh_gcell = 0._r8
    npp_ind = 0._r8

    do p = 1,n_pft
       if(fpcgrid(p)>0.0.and.ivt(p).ne.17)then
         npp_ind(p) = bm_inc(p)/fpcgrid(p)    ![gC/m2 Patch area]
         anpp_gcell = anpp_gcell + bm_inc(p)  ![gC/m2 vegetated area] for output
       endif

       amrh_gcell = amrh_gcell + afmicr(p)    ![gC/m2 cell veg'd area]
    end do


    ! Build filter of present natually-vegetated pfts
    do p = 1, n_pft
        if(wt_patch(p).gt.1.0E-16 .and. ivt(p).ne.17) ifpre(p)=1.           !revised by zhq. dec30.08
    end do

    call BuildNatVegFilter(n_pft, ivt, ifpre, num_natvegp, filter_natvegp)

#ifdef MYBUG
    print *, 'num_natvegp', num_natvegp
    print *, 'filter', filter_natvegp
#endif

    ! Returns updated bm_inc, litterag

    call Reproduction(n_pft, num_natvegp, filter_natvegp, litter_ag, bm_inc)

#ifdef MYBUG
    print *, 'bm_inc<-reproduction', bm_inc
    print *, 'litter_ag<-reproduction', litter_ag
#endif

    ! Returns turnover_ind and updated litterag,bg, l,s,h,rm_ind

    call Turnover(nind(1:), pftpar(9,1:), pftpar(11,1:), pftpar(12,1:),&
                  n_pft, num_natvegp, filter_natvegp, turnover_ind,&
                  litter_ag(1:), litter_bg(1:), lm_ind(1:), sm_ind(1:),&
                  hm_ind(1:), rm_ind(1:) )

#ifdef MYBUG
    print *, 'turnover_ind<-turnover', turnover_ind
    print *, 'lm_ind<-turnover', lm_ind
#endif

    ! Returns updated litterag, bg, and ifpre

    call Kill(n_pft, num_natvegp, filter_natvegp,&
              litter_ag(1:), litter_bg(1:), ifpre(1:),&
              nind(1:), lm_ind(1:), sm_ind(1:), hm_ind(1:),&
              rm_ind(1:), bm_inc(1:), vegclass(1:), ivt(1:))

#ifdef MYBUG
    print *, 'ifpre<-Kill', ifpre
    print *, 'litter_ag<-Kill', litter_ag
#endif

    ! Rebuild filter of present natually-vegetated pfts after Kill()

    call BuildNatVegFilter(n_pft, ivt, ifpre, num_natvegp, filter_natvegp)

    call Allocation(n_pft, num_natvegp, filter_natvegp,&
                    lai_ind(1:), fpc_inc(1:), crownarea(1:), htop(1:),&
                    lm_ind(1:), sm_ind(1:), hm_ind(1:), rm_ind(1:),&
                    litter_ag(1:), litter_bg(1:), fpcgrid(1:), sla(1:),&
                    vegclass(1:), pftpar(20,1:), pftpar(18,1:),&
                    bm_inc(1:), nind(1:), annpsn(1:), annpsnpot(1:),ivt(1:))

#ifdef MYBUG
    print *, 'crownarea<-alloc', crownarea
    print *, 'lai_ind<-alloc', lai_ind
    print *, 'htop<-alloc', htop
    print *, 'fpc_inc<-alloc', fpc_inc
#endif

    call Light(n_pft, num_natvegp, filter_natvegp,ivt, &
               fpc_inc, sm_ind, hm_ind, lai_ind, crownarea, &
               sla, vegclass, fpcgrid, nind, litter_ag, &
               litter_bg, lm_ind, rm_ind, pftpar(20,1:))

    ! Obtain updated ifpre, nind, litterag and bg

#ifdef MYBUG
    print *, 'fpcgrid->mort', fpcgrid
    print *, 'lai_ind->mort', lai_ind
    print *, 'lm_ind->mort', lm_ind
    print *, 'nind->mort', nind
    print *, 'litter_ag->mort',litter_ag
#endif

    call Mortality(n_pft, num_natvegp, filter_natvegp(1:) ,&
                   ivt(1:), bm_inc(1:), lm_ind(1:), sm_ind(1:), hm_ind(1:),&
                   rm_ind(1:), agddtw(1:), turnover_ind(1:), vegclass(1:),&
                   sla(1:), ifpre(1:), nind(1:), litter_ag(1:), litter_bg(1:))

#ifdef MYBUG
    print *, 'ifpre<-mort', ifpre
    print *, 'lm_ind<-mort', lm_ind
#endif

    ! Returns updated litterag and nind

    call Fire(n_pft, afirefrac(1:), acfluxfire(1:),&
              ivt(1:), lm_ind(1:), sm_ind(1:),&
              hm_ind(1:), rm_ind(1:), fpcgrid(1:),&
              ifpre(1:), firelength(1:), vegclass(1:), &
              pftpar(8,1:), litter_ag(1:), nind(1:))

#ifdef MYBUG
    print *, 'ifpre<-fire', ifpre
    print *, 'fpcgrid<-fire', fpcgrid
    print *, 'lm_ind<-Fire', lm_ind
    print*,'litter_ag<-Fire',litter_ag
#endif

    afiref_gcell = 0.0
    afirec_gcell = 0.0
    do p = 1,n_pft
       afiref_gcell = afiref_gcell + afirefrac(p)*fpcgrid(p)
       afirec_gcell = afirec_gcell + acfluxfire(p)
    end do

    ! Returns updated present, nind, *m_ind, crownarea, fpcgrid, htop,
    ! litter*g, and vegetation type
    ! reason for different set up (ie, no external patch loop):
    ! in this routine sub-grid patches (k) communicate at the grid cell level (i,j)

    call Establishment(n_pft,tmomin20,agdd20,agddtw,&
                       prec365,sla,pftpar(20,1:),lm_sapl,&
                       sm_sapl,hm_sapl,rm_sapl,pftpar(28,1:),pftpar(29,1:),&
                       pftpar(30,1:),ivt,ifpre,nind,lm_ind,sm_ind,&
                       hm_ind,rm_ind,litter_ag,litter_bg,&
                       vegclass,fpcgrid,htop,lai_ind,crownarea,aestabc_gcell)

#ifdef MYBUG
    print *, 'fpcgrid<-estab', fpcgrid
    print *, 'litter_ag<-estab', litter_ag
    print *, 'lai_ind<-estab', lai_ind
    print *, 'crownarea<-estab', crownarea
    print *, 'ifpre<-estab', ifpre
    print *, 'ivt<-estab', ivt
    print *, 'lm_ind<-estab', lm_ind
    print *, 'nind<-estab', nind
#endif

    ! Rebuild filter of present natually-vegetated pfts after Kill()

    call BuildNatVegFilter(n_pft, ivt, ifpre, num_natvegp, filter_natvegp)
    
    ! reset fpc_inc following X.D.Z, before recalling the Light module.
    fpc_inc = 0.0

    call Light(n_pft, num_natvegp, filter_natvegp, ivt, &
               fpc_inc, sm_ind, hm_ind, lai_ind, crownarea, &
               sla, vegclass, fpcgrid, nind, litter_ag, &
               litter_bg, lm_ind, rm_ind, pftpar(20,1:))

#ifdef MYBUG
    print *, 'lai_ind<-light2', lai_ind
    print *, 'fpcgrid<-Light2', fpcgrid
    print *, 'nind<-Light2', nind
    print *, 'litter_ag<-Light2', litter_ag
#endif

    ! Returns lm,rm_ind, fpcgrid, nind, litterag,bg via modules
    ! reason for different set up (ie, no external patch loop):
    ! in this routine sub-grid patches (k) communicate at the grid cell level (i,j)

    alitc_ag = 0.
    alitc_bg = 0.
    asoic_fast = 0.
    asoic_slow = 0.
    avegc_gcell = 0.

    do p = 1, n_pft  
       alitc_ag = alitc_ag + litter_ag(p)
       alitc_bg = alitc_bg + litter_bg(p)
       asoic_fast = asoic_fast + cpool_fast(p)
       asoic_slow = asoic_slow + cpool_slow(p)
       avegc_gcell = avegc_gcell + (lm_ind(p) + sm_ind(p) + hm_ind(p) + rm_ind(p))*nind(p)
    end do

    call ResetBareFPC(n_pft, ivt, fpcgrid, lai_ind, ifpre)

#ifdef MYBUG
    print *, 'fpcgrid<-reset', fpcgrid
#endif

    do p = 1, n_pft
       wt_patch(p) = wt_column * fpcgrid(p)
    enddo

!   Transform ivt from integer to real
    ivt_r = ivt

#ifdef MYBUG
    print *, 'ivt-2', ivt
    print *, 'ifpre-2', ifpre
    print *, 'lai-2',lai_ind
    print *, 'fpcgrid-2',fpcgrid
    print *, 'wt_patch-2', wt_column
#endif

  end subroutine LPJ


!====================================================================================

  subroutine BuildNatVegFilter(n_pft, ivt, ifpre, num_natvegp, filter_natvegp)

!------------------------------------------------------------------------
! Reconstruct a filter of naturally-vegetated PFTs for use in DGVM
! CALLED FROM: subroutine lpj
!------------------------------------------------------------------------

    use precision
    implicit none
    integer , INTENT(in) :: n_pft                     ! total patches
    integer , INTENT(in) :: ivt(n_pft)                ! pft vegetation (pft level)
    real(r8), INTENT(in) :: ifpre(n_pft)              ! whether this pft present in patch
    integer , INTENT(out):: num_natvegp               ! number of pfts in naturally-vegetated filter
    integer , INTENT(out):: filter_natvegp(n_pft)     ! pft filter for naturally-vegetated points
!
! LOCAL VARIABLES:
    integer :: p
!-----------------------------------------------------------------------

    num_natvegp = 0
    do p = 1,n_pft
       if (ivt(p) > 0 .and. ivt(p) .lt. 17) then
          if(ifpre(p) .gt. 0.)then
             num_natvegp = num_natvegp + 1
             filter_natvegp(num_natvegp) = p
          endif
       end if
    end do

  end subroutine BuildNatVegFilter

!=================================================================================

  subroutine Reproduction(n_pft, num_natvegp, filter_natvegp,&
                          litter_ag, bm_inc)
!-------------------------------------------------------------------
! 10% of biomass increase is assumed to reproduct new generation
! CALLED FROM: subroutine lpj
!--------------------------------------------------------------------

    use precision
    implicit none
    integer , INTENT(in) :: n_pft             ! total patches
    integer , INTENT(in) :: num_natvegp        ! number of vegetated pfts in this grid
    integer , INTENT(in) :: filter_natvegp(n_pft) ! pft filter for naturally-vegetated points
    real(r8), INTENT(inout) :: litter_ag(n_pft)    ! above ground litter
    real(r8), INTENT(inout) :: bm_inc(n_pft)       ! biomass increment

! LOCAL VARIABLES:
    real(r8), parameter :: reprod_cost = 0.1 ! proportion of NPP lost to reproduction (Harper 1977)
    integer  :: p, fp                        ! pft index
    real(r8) :: reprod                       ! temporary
!-----------------------------------------------------------------------

    ! Compute reproduction costs

    do fp = 1,num_natvegp
       p = filter_natvegp(fp)

       ! Calculate allocation to reproduction
       ! Reproduction costs taken simply as a constant fraction of annual NPP

       reprod = max(bm_inc(p) * reprod_cost, 0.0)

       ! assume the costs go to reproductive structures which will
       ! eventually enter the litter pool

       litter_ag(p) = litter_ag(p) + reprod

       ! Reduce biomass increment by reproductive cost

       bm_inc(p) = bm_inc(p) - reprod

    end do

  end subroutine Reproduction

!==============================================================================


  subroutine Turnover(nind, l_turn, s_turn, r_turn,&
                      n_pft, num_natvegp, filter_natvegp, turnover_ind,&
                      litter_ag, litter_bg, lm_ind, sm_ind, hm_ind, rm_ind)
!-----------------------------------------------------------------------------
! Turnover of PFT-specific fraction from each living C pool
! Leaf and root C transferred to litter, sapwood C to heartwood
! Called once per year
! CALLED FROM: subroutine lpj
!----------------------------------------------------------------------------
    use precision
    implicit none
    integer, INTENT(in) :: n_pft                 ! total patches
    integer, INTENT(in) :: num_natvegp            ! number of vegetated patches in grid
    integer, intent(in) :: filter_natvegp(n_pft) ! pft filter for naturally-vegetated points

! INTENT IN VARIABLES:
    real(r8), INTENT(in) :: nind(n_pft)            ! number of individuals (#/m**2)
    real(r8), INTENT(in) :: l_turn(n_pft)          ! ecophys const - leaf turnover period [years]
    real(r8), INTENT(in) :: s_turn(n_pft)          ! ecophys const - sapwood turnover period [years]
    real(r8), INTENT(in) :: r_turn(n_pft)          ! ecophys const - root turnover period [years]

! INTENT INOUT VARIABLES:
    real(r8), INTENT(inout) :: litter_ag(n_pft)       ! above ground litter
    real(r8), INTENT(inout) :: litter_bg(n_pft)       ! below ground litter
    real(r8), INTENT(inout) :: lm_ind(n_pft)          ! individual leaf mass
    real(r8), INTENT(inout) :: sm_ind(n_pft)          ! individual sapwood mass
    real(r8), INTENT(inout) :: hm_ind(n_pft)          ! individual heartwood mass
    real(r8), INTENT(inout) :: rm_ind(n_pft)          ! individual root mass

! INTENT OUT VARIABLES:
    real(r8), INTENT(out) :: turnover_ind(n_pft)    !

! LOCAL VARIABLES:
    integer  :: p, fp
    real(r8) :: l_torate
    real(r8) :: s_torate
    real(r8) :: r_torate
    real(r8) :: lm_turn
    real(r8) :: sm_turn
    real(r8) :: rm_turn

    ! Determine turnover of pft-specific fraction from each living
    ! C pool

    do fp = 1,num_natvegp
       p = filter_natvegp(fp)

       ! Turnover rates are reciprocals of tissue longevity

       l_torate = 1.0 / l_turn(p)
       s_torate = 1.0 / s_turn(p)
       r_torate = 1.0 / r_turn(p)
!print *, 'l_turn', l_turn(p), 'l_torate', l_torate
!print *, 's_turn', s_turn(p), 's_torate', s_torate
!print *, 'r_turn', r_turn(p), 'r_torate', r_torate

       ! Calculate the biomass turnover in this year

       lm_turn = lm_ind(p) * l_torate
       sm_turn = sm_ind(p) * s_torate
       rm_turn = rm_ind(p) * r_torate

       ! Update the pools

       lm_ind(p) = lm_ind(p) - lm_turn
       sm_ind(p) = sm_ind(p) - sm_turn
       rm_ind(p) = rm_ind(p) - rm_turn

       ! Convert sapwood to heartwood

       hm_ind(p) = hm_ind(p) + sm_turn

       ! Transfer to litter pools

       litter_ag(p) = litter_ag(p) + lm_turn * nind(p)
       litter_bg(p) = litter_bg(p) + rm_turn * nind(p)

       ! Record total turnover

       turnover_ind(p) = lm_turn + sm_turn + rm_turn

    end do

  end subroutine Turnover

!==========================================================================

  subroutine Kill(n_pft, num_natvegp, filter_natvegp,&
                  litter_ag, litter_bg, ifpre, &
                  nind, lm_ind, sm_ind, hm_ind,&
                  rm_ind, bm_inc, vegclass, ivt)

!---------------------------------------------------------------------
! Removal of PFTs with negative annual C increment
! NB: PFTs newly beyond their bioclimatic limits are removed in
! subroutine establishment
! Called once per year
! CALLED FROM: subroutine lpj
!---------------------------------------------------------------------
    use precision
    implicit none

    integer, INTENT(in) :: n_pft            ! total patches
    integer, INTENT(in) :: num_natvegp       ! number of vegetated patches in grid
    integer, INTENT(in) :: filter_natvegp(n_pft) ! pft filter for naturally-vegetated points

! INTENT IN VARIABLES:
    real(r8), INTENT(in) :: nind(n_pft)            ! number of individuals (#/m**2)
    real(r8), INTENT(in) :: lm_ind(n_pft)          ! individual leaf mass
    real(r8), INTENT(in) :: sm_ind(n_pft)          ! individual sapwood mass
    real(r8), INTENT(in) :: hm_ind(n_pft)          ! individual heartwood mass
    real(r8), INTENT(in) :: rm_ind(n_pft)          ! individual root mass
    real(r8), INTENT(in) :: bm_inc(n_pft)          ! biomass increment
    real(r8), INTENT(in) :: vegclass(n_pft)        ! ecophys const 

! INTENT INOUT VARIABLES:
!
    integer , INTENT(inout) :: ivt(n_pft)             ! land cover type
    real(r8), INTENT(inout) :: ifpre(n_pft)           ! whether PFT present in patch
    real(r8), INTENT(inout) :: litter_ag(n_pft)       ! above ground litter
    real(r8), INTENT(inout) :: litter_bg(n_pft)       ! below ground litter
!EOP
!
! !LOCAL VARIABLES:
    integer :: p,fp

    do fp = 1,num_natvegp
       p = filter_natvegp(fp)

       if (bm_inc(p) < 0.0) then !negative C increment this year

          ifpre(p) = -1.   !remove PFT
          ivt(p) = 17

          ! Transfer killed biomass to litter

          if (vegclass(p) .le. 2.) then  ! woody PFTs
             litter_ag(p) = litter_ag(p) + (lm_ind(p) + sm_ind(p) + hm_ind(p)) * nind(p)
          else         ! herbaceous PFTs
             litter_ag(p) = litter_ag(p) + lm_ind(p) * nind(p)
          end if

          litter_bg(p) = litter_bg(p) + rm_ind(p) * nind(p)

       end if
    end do

  end subroutine Kill

!===================================================================================

  subroutine Allocation (n_pft, num_natvegp, filter_natvegp,&
                         lai_ind, fpc_inc, crownarea, htop,&
                         lm_ind, sm_ind, hm_ind, rm_ind,&
                         litter_ag, litter_bg, fpcgrid, sla,&
                         vegclass, crownarea_max, init_lmtorm,&
!                        bm_inc,nind, annpsn, annpsnpot,ivt,fveg)
                         bm_inc,nind, annpsn, annpsnpot,ivt)

!------------------------------------------------------------------
! Performs yearly allocation calculation
! Allocation of this year's biomass increment (bm_inc_ind) to the
! three living carbon pools, such that the basic allometric
! relationships (A-C below) are always satisfied.
! CALLED FROM: subroutine lpj
! -------------------------------------------------------------------
! TREE ALLOCATION
! (A) (leaf area) = latosa * (sapwood xs area)
!       (Pipe Model, Shinozaki et al. 1964a,b; Waring et al 1982)
! (B) (leaf mass) = lmtorm * (root mass)
! (C) height = allom2 * (stem diameter)**allom3 (source?)
! (D) (crown area) = min (allom1 * (stem diameter)**reinickerp, crownarea_max)
!---------------------------------------------------------------------
! Adapted from Xiaodong Zeng's version of DGVMAllocationa, 
! Different parameterization for shrub is included in.

    use precision
    implicit none

    integer, intent(in) :: n_pft                    ! total patches
    integer, intent(in) :: num_natvegp               ! number of naturally-vegetated pfts in filter
    integer, intent(in) :: filter_natvegp(n_pft)    ! pft filter for naturally-vegetated points

! INTENT IN VARIABLES:

    integer,  INTENT(in) :: ivt(n_pft)             ! land cover types
    real(r8), INTENT(in) :: sla(n_pft)             ! specific leaf area [m2 leaf g-1 carbon]
    real(r8), INTENT(in) :: vegclass(n_pft)        ! 1.tree 2.shrub 3.grass 4.crop 
    real(r8), INTENT(in) :: crownarea_max(n_pft)   ! tree maximum crown area [m2]
    real(r8), INTENT(in) :: init_lmtorm(n_pft)     ! non-water stressed leaf:root ratio
    real(r8), INTENT(in) :: bm_inc(n_pft)          ! biomass increment
    real(r8), INTENT(in) :: nind(n_pft)            ! number of individuals (#/m**2)
    real(r8), INTENT(in) :: annpsn(n_pft)          ! annual photosynthesis (umol CO2 /m**2)
    real(r8), INTENT(in) :: annpsnpot(n_pft)       ! annual potential photosynthesis (..) 

! INTENT INOUT VARIABLES:
    real(r8), INTENT(inout) :: fpcgrid(n_pft)      ! foliar projective cover on gridcell
    real(r8), INTENT(inout) :: crownarea(n_pft)    ! area each individual takes up (m^2)
    real(r8), INTENT(inout) :: htop(n_pft)         ! canopy top (m)
    real(r8), INTENT(inout) :: lm_ind(n_pft)       ! individual leaf mass
    real(r8), INTENT(inout) :: sm_ind(n_pft)       ! individual sapwood mass
    real(r8), INTENT(inout) :: hm_ind(n_pft)       ! individual heartwood mass
    real(r8), INTENT(inout) :: rm_ind(n_pft)       ! individual root mass
    real(r8), INTENT(inout) :: litter_ag(n_pft)    ! above ground litter
    real(r8), INTENT(inout) :: litter_bg(n_pft)    ! below ground litter

! INTENT OUT VARIABLES:
    real(r8), INTENT(out) :: lai_ind(n_pft)          ! LAI per individual
    real(r8), INTENT(out) :: fpc_inc(n_pft)          ! foliar projective cover increment 

! LOCAL VARIABLES:
    integer , parameter :: nseg = 20
    real(r8), parameter :: xacc = 0.1           ! threshold x-axis and threshold
    real(r8), parameter :: yacc = 1.0e-10       ! y-axis precision of allocation soln
    real(r8), parameter :: reinickerp = 1.6     ! parameter in allometric equation
    real(r8), parameter :: wooddens = 2.0e5     ! wood density (gC/m3)
    real(r8), parameter :: PI = 3.14159265358979323846 
    real(r8) :: latosa(21)                      ! leafarea:sapwood cross-sectional area
    real(r8) :: allom1(21)                      ! parameters in allometric
    real(r8) :: allom2(21)                      ! parameters in allometric
    real(r8) :: allom3(21)                      ! parameters in allometric
    integer  :: p, fp                           ! index
    integer  :: fn,fnold                        ! number of elements in local filter
    real(r8) :: wscal
    real(r8) :: lmtorm(n_pft)
    real(r8) :: bm_inc_ind(n_pft)
    real(r8) :: lminc_ind_min(n_pft)
    real(r8) :: rminc_ind_min(n_pft)
    real(r8) :: lminc_ind(n_pft)
    real(r8) :: rminc_ind(n_pft)
    real(r8) :: sminc_ind(n_pft)
    real(r8) :: fpc_grid_old(n_pft)
    real(r8) :: x1(n_pft), x2(n_pft), dx(n_pft)
    real(r8) :: sap_xsa(n_pft)
    real(r8) :: stemdiam(n_pft)
    real(r8) :: fx1(n_pft)
    real(r8) :: fmid(n_pft)
    real(r8) :: xmid(n_pft)
    real(r8) :: sign(n_pft)
    real(r8) :: rtbis(n_pft)
    real(r8) :: fpc_ind(n_pft)                 ! individual fpc

    real(r8) :: lm, rm, sm, hm                 ! local biomass; added by X.D.Z
    real(r8) :: lm0, lm1, err_hm, csum         ! temporal variables; added by X.D.Z
    real(r8) :: sa, hd, dd, vd                 ! temporal variables; added by X.D.Z
    real(r8) :: k_rl, ksla, k_lasa, rou, k_a2, k_a3 ! temporal variables; added by X.D.Z
    real(r8) :: eps_alloc                      ! interation eps; added by X.D.Z
    integer  :: n_step                         ! indices; added by X.D.Z

    latosa(1:21) = 8.0e3   ! leafarea:sapwood cross-sectional area
    allom1(1:21) = 100.0   ! parameters in allometric
    allom2(1:21) =  40.0   ! parameters in allometric
    allom3(1:21) =   0.5   ! parameters in allometric

    latosa(9:11) = 4000.0                        ! 2000.0 in last test
    allom1(9:11) = 200.0
    allom2(9:11) = 10.0

    eps_alloc = 1.0E-6

    do fp = 1,num_natvegp
       p = filter_natvegp(fp)
       bm_inc_ind(p) = bm_inc(p) / nind(p)
!print*, p, 'alloc-1',bm_inc_ind(p),fpcgrid(p),nind(p)
    end do

    do fp = 1,num_natvegp
       p = filter_natvegp(fp)

       if (vegclass(p) .le. 2) then               ! tree allocation

          ! calculate this year's leaf to fine root mass ratio from mean annual
          ! water scalar and pft specific parameter
          ! slevis: in lpj wscal=awscal(pft)/aleafdays(pft), awscal=SUM(dwscal) (dphen>0)
          !         dwscal=min(supply(pft)/demandpot(pft),1) or =1 when demand(pft)=0 etc
          !         here wscal=annpsn/annpsnpot
          if (annpsnpot(p) > 0.0) then
             wscal = annpsn(p)/annpsnpot(p)
          else
             wscal = 1.0
          end if

          hm = hm_ind(p)
          csum = bm_inc_ind(p) + lm_ind(p) + rm_ind(p) + sm_ind(p)
  
          k_lasa = latosa(ivt(p))
          ksla = sla(p)
          rou = wooddens
          k_a2 = allom2(ivt(p))
          k_a3 = allom3(ivt(p))
          lmtorm(p) = init_lmtorm(p) * wscal
          k_rl = 1.0 / lmtorm(p)
  
          lm0 = 0.0
          lm1 = csum / (1.0 + k_rl)    !give a guessed value of above ground cpool.
          err_hm = eps_alloc + eps_alloc                    !
          n_step = 0

#ifdef MYBUG
          print *, 'wscal', wscal, annpsn(p), annpsnpot(p), k_rl, lm1, csum
#endif
  
        ! Avoid unbelievable value for very small lm. Modified by zhq @09/14/2010
        ! do while ((lm1 - lm0) > eps_alloc .and. (err_hm > eps_alloc .or. err_hm < -eps_alloc) .and. n_step < 25)
          do while ((lm1 - lm0) > eps_alloc .or. (err_hm > eps_alloc .or. err_hm < -eps_alloc) .and. n_step < 25)
             lm = (lm0 + lm1) * 0.5
             rm = k_rl * lm
             sm = csum - lm - rm
             sa = lm * ksla / k_lasa
             hd = sm / (rou * sa)
             dd = (hd/k_a2) ** (1.0/k_a3)
             vd = PI * .25 * dd * dd * hd;
             err_hm = rou * vd - sm - hm; 
             if (err_hm > 0) then
                lm0 = lm
             else
                lm1 = lm
             end if
             n_step = n_step + 1
          end do

#ifdef MYBUG
          print *,'alloc-2', lm0, lm1, lm, n_step,dd
#endif

          lm_ind(p) = lm
          rm_ind(p) = rm
          sm_ind(p) = sm

          ! Calculate new height, diameter and crown area

          sap_xsa(p) = lm_ind(p) * ksla / k_lasa  !eqn (5)

          !BUGFIX
          if (lm_ind(p) == 0) then
             htop(p) = 0.
          else
             htop(p) = sm_ind(p) / sap_xsa(p) / rou
          end if
          !BUGFIX

        ! stemdiam = (4.0*(sm_ind(p) + hm_ind(p))/wooddens/PI/allom2(ivt(p)))**(1.0/(2.0+allom3(ivt(p))))  ! Eqn 9
          stemdiam(p) = dd 
          crownarea(p) = min(allom1(ivt(p))*stemdiam(p)**reinickerp, crownarea_max(p)) !eqn (D)
        ! crownarea(p) = allom1(ivt(p))*stemdiam(p)**reinickerp !eqn (D)

#ifdef MYBUG
          print *, 'Allocate-CA', ivt(p),dd, crownarea(p),lm,sm,rm
#endif
       else !grasses

          lmtorm(p) = init_lmtorm(p) !slevis: eliminate influence from soil H2O

          ! GRASS ALLOCATION
          ! Distribute this year's production among leaves and fine roots
          ! according to leaf to rootmass ratio [eqn (33)]
          ! Relocation of C from one compartment to the other not allowed:
          ! negative increment in either compartment transferred to litter

          lminc_ind(p) = (bm_inc_ind(p) - lm_ind(p)/lmtorm(p) + rm_ind(p)) / (1.0 + 1.0/lmtorm(p))
          rminc_ind(p) = bm_inc_ind(p) - lminc_ind(p)

          if (lminc_ind(p) >= 0.0) then

             ! Add killed roots (if any) to below-ground litter

             ! CHECK: take out if statement because if rminc is negative than
             ! root mass has been translocated to the leaves, therefore mass balance
             ! problem since this carbon stays in the vegetation but is in addition
             ! added to the litter pool. ALLOW translocation from roots to leaves
             ! i.e. assume carbon stores in the roots which can be delivered
             ! to the leaves under times of stress.

             ! if (rminc_ind < 0.0) litter_bg = litter_bg -rminc_ind * nind

          else

             ! Negative allocation to leaf mass

             rminc_ind(p) = bm_inc_ind(p)
             lminc_ind(p) = (rm_ind(p) + rminc_ind(p))*lmtorm(p) - lm_ind(p) !from eqn (9)

             ! Add killed leaves to litter

             litter_ag(p) = litter_ag(p) - lminc_ind(p) * nind(p)

          end if

          ! Increment C compartments

          lm_ind(p) = lm_ind(p) + lminc_ind(p)
          rm_ind(p) = rm_ind(p) + rminc_ind(p)
          ! crownarea(p) = 1.0                
          crownarea(p) = 1.0 - exp(-0.5 * lm_ind(p) * sla(p)) !zhq:follow X.D.Z's revision.07/22/2010                

       end if   ! end of tree or grass

       ! Update LAI and FPC

       if (crownarea(p) > 0.0) then
          lai_ind(p) = lm_ind(p) * sla(p) / crownarea(p)
       else
          lai_ind(p) = 0.0
       end if

       fpc_grid_old(p) = fpcgrid(p)
     ! fpc_ind(p) = 1. - exp (-0.5 * lai_ind(p))
     ! fpcgrid(p) = crownarea(p) * nind(p) * fpc_ind(p)
       fpcgrid(p) = crownarea(p) * nind(p) 
       fpc_inc(p) = max(0.0, fpcgrid(p) - fpc_grid_old(p))

    end do

  end subroutine Allocation

!=======================================================================

  subroutine Mortality(n_pft, num_natvegp, filter_natvegp,&
                       ivt, bm_inc, lm_ind, sm_ind, hm_ind,&
                       rm_ind, agddtw, turnover_ind, vegclass,&
                       sla, ifpre, nind, litter_ag, litter_bg)
!---------------------------------------------------------------------
! Tree background and stress mortality
! CALLED FROM:subroutine lpj
!---------------------------------------------------------------------
    use precision
    implicit none

! INTENT IN VARIABLES:
    integer ,INTENT(in) :: n_pft               ! total patches
    integer ,INTENT(in) :: num_natvegp         !number of naturally-vegetated pfts in filter
    integer ,INTENT(in) :: filter_natvegp(n_pft)!pft filter for naturally-vegetated points
    real(r8),INTENT(in) :: bm_inc(n_pft)       ! biomass increment
    real(r8),INTENT(in) :: lm_ind(n_pft)       ! individual leaf mass
    real(r8),INTENT(in) :: sm_ind(n_pft)       ! individual sapwood mass
    real(r8),INTENT(in) :: hm_ind(n_pft)       ! individual heartwood mass
    real(r8),INTENT(in) :: rm_ind(n_pft)       ! individual root mass
    real(r8),INTENT(in) :: agddtw(n_pft)       ! accumulated growing degree days above twmax
    real(r8),INTENT(in) :: turnover_ind(n_pft) !
    real(r8),INTENT(in) :: vegclass(n_pft)     ! pft class
    real(r8),INTENT(in) :: sla(n_pft)          ! specific leaf area [m2 leaf g-1 carbon]

! INTENT INOUT VARIABLES: 
    integer ,INTENT(inout) :: ivt(n_pft)          ! pft vegetation type
    real(r8),INTENT(inout) :: ifpre(n_pft)        ! whether PFT present in patch
    real(r8),INTENT(inout) :: nind(n_pft)         ! number of individuals
    real(r8),INTENT(inout) :: litter_ag(n_pft)    ! above ground litter
    real(r8),INTENT(inout) :: litter_bg(n_pft)     ! below ground litter

! LOCAL VARIABLES:
    real(r8), parameter :: k_mort = 0.3 !coefficient of growth efficiency in mortality equation
    real(r8), parameter :: ramp_agddtw = 300.0
    integer  :: p,fp      !index
    real(r8) :: mort_max  ! asymptotic maximum mortality rate (/yr)
    real(r8) :: bm_delta  ! net individual living biomass increment
                          ! (incorporating loss through leaf, root and sapwood turnover) (gC)
    real(r8) :: mort      ! tree mortality rate
    real(r8) :: nind_kill ! reduction in individual density due to mortality (indiv/m2)
    real(r8) :: greffic
    real(r8) :: heatstress

    ! Compute stress mortality

!   print *, 'agddtw->mort', agddtw 
    do fp = 1,num_natvegp
       p = filter_natvegp(fp)
       if (vegclass(p) .le. 2.) then

          if (ivt(p)==3 .or. ivt(p) == 8) then
             mort_max = 0.03 !testing diff values
          else
             mort_max = 0.01 !original value for all pfts
          end if

          heatstress = min(1.0, agddtw(p) / ramp_agddtw)

          ! Calculate net individual living biomass increment
          bm_delta = max(0.0, bm_inc(p) / nind(p) - turnover_ind(p))

          ! Calculate growth efficiency (net biomass increment per unit leaf area)

          !BUGFIX
          if (lm_ind(p) == 0) then
             greffic = 0.
          else
             greffic = bm_delta / lm_ind(p) / sla(p)
          end if
!print *, ivt(p),'greffic', greffic
          !BUGFIX

          ! Mortality rate inversely related to growth efficiency (Prentice et al 1993)
          mort = mort_max / (1.0 + k_mort * greffic)

          ! Reduce individual density (=> gridcell-level biomass) by mortality rate
          mort = min(1.0, mort + heatstress)
          nind_kill = nind(p) * mort
          nind(p) = nind(p) - nind_kill

          ! Transfer lost biomass to litter
          litter_ag(p) = litter_ag(p) + nind_kill * (lm_ind(p) + sm_ind(p) + hm_ind(p))
          litter_bg(p) = litter_bg(p) + nind_kill * rm_ind(p)
       end if

       if (nind(p) < 1.0E-6) then
          ifpre(p) = -1.
          ivt(p) = 17
       end if
    end do
  end subroutine Mortality

!=================================================================================

  subroutine Fire(n_pft, afire_frac, acflux_fire, &
                  ivt, lm_ind, sm_ind, hm_ind, rm_ind, &
                  fpcgrid, ifpre, firelength, vegclass, &
                  resist, litter_ag, nind)
                    
!---------------------------------------------------------------------
! Effect of the fire on vegetation structure and litter
! CALLED FROM: subroutine lpj
!---------------------------------------------------------------------
    use precision
    implicit none

    integer , intent(in)  :: n_pft             ! total patches
    real(r8), intent(out) :: afire_frac(n_pft)
    real(r8), intent(out) :: acflux_fire(n_pft)

! INTENT IN VARIABLES: 
    integer , INTENT(in) :: ivt(n_pft)           ! pft vegetation type
    real(r8), INTENT(in) :: lm_ind(n_pft)        ! individual leaf mass
    real(r8), INTENT(in) :: sm_ind(n_pft)        ! individual sapwood mass
    real(r8), INTENT(in) :: hm_ind(n_pft)        ! individual heartwood mass
    real(r8), INTENT(in) :: rm_ind(n_pft)        ! individual root mass
    real(r8), INTENT(in) :: fpcgrid(n_pft)       ! foliar projective cover on gridcell
    real(r8), INTENT(in) :: ifpre(n_pft)         ! whether this pft present in patch
    real(r8), INTENT(in) :: firelength(n_pft)    ! fire season in days
    real(r8), INTENT(in) :: vegclass(n_pft)      ! 1.tree 2.shrub 3.grass 4.crop
    real(r8), INTENT(in) :: resist(n_pft)        ! fire resistance index 

! INTENT INOUT VARIABLES:
    real(r8), INTENT(inout) :: litter_ag(n_pft)  ! above ground litter
    real(r8), INTENT(inout) :: nind(n_pft)       ! number of individuals (/m**2)

! LOCAL VARIABLES:
    real(r8), parameter :: minfuel = 200.0  ! fuel threshold to carry a fire (gC/m2)
    real(r8) :: fire_index
    real(r8) :: fire_term
    real(r8) :: disturb
    integer  :: p                           ! index
!-----------------------------------------------------------------------
 
    acflux_fire = 0.
    do p = 1, n_pft
       if (ivt(p) .le. 16) then

          !fire_length comes from subroutine FireSeason, which is called once per step

          ! Calculate annual fire index

          fire_index = firelength(p) / 365.0

          ! Calculate the fraction of the grid cell affected by fire

          fire_term = fire_index - 1.0
          afire_frac(p) = max(fire_index * &
               exp(fire_term / (-0.13*fire_term**3 + 0.6*fire_term**2 + 0.8*fire_term + 0.45)), &
               0.001)
!         print*,'fire',p,firelength(p), afire_frac(p)

          ! Reduce fraction of patch affected by fire when fuel
          ! becomes limiting (reduced carrying capacity)

          if (litter_ag(p) < minfuel * fpcgrid(p)) afire_frac(p) = 0.001

          ! Implement the effect of the fire on vegetation structure and litter
          ! in the disturbed fraction.

          ! Each PFT is assigned a resistance to fire, representing the fraction of
          ! the PFT which survives a fire. Grasses assumed already to have completed
          ! their life cycle and thus are not affected by fire, giving them
          ! a competitive advantage against woody PFTs.

          if (ifpre(p) .gt. 0.) then
             if ( vegclass(ivt(p)) .le. 2.) then
                ! Calculate the fraction of individuals in grid cell which die
                ! (slevis: 'in grid cell' because nind is grid average)

                disturb = (1.0-resist(ivt(p))) * afire_frac(p)
!print *, 'disturb', disturb

                ! Calculate carbon flux to atmosphere (gC/m2) due to burnt biomass

                acflux_fire(p) = disturb * nind(p) * (lm_ind(p) + sm_ind(p) + hm_ind(p) + rm_ind(p))

                ! Update the individual density

                nind(p) = nind(p) * (1.0-disturb)

                ! Add combusted litter to carbon flux to atmosphere term

                acflux_fire(p) = acflux_fire(p) + afire_frac(p) * litter_ag(p)
             end if
          else
             acflux_fire(p) = afire_frac(p) * litter_ag(p)
          end if

          ! Update the above ground litter term

          litter_ag(p) = (1.0 - afire_frac(p)) * litter_ag(p)

       else

          afire_frac(p) = 0.0
          acflux_fire(p) = 0.0

       end if
    end do

  end subroutine Fire

!==================================================================================

  subroutine Establishment(n_pft,tmomin20,agdd20,agddtw,&
                           prec365,sla,crownarea_max,lm_sapl,&
                           sm_sapl,hm_sapl,rm_sapl,tcmin,tcmax,&
                           gddmin,ivt,ifpre,nind,lm_ind,sm_ind,&
                           hm_ind,rm_ind,litter_ag,litter_bg,&
                           vegclass,fpcgrid,htop,lai_ind,crownarea,&
                           aestabc_gcell)

!---------------------------------------------------------------------------------
! Calculates establishment of new pfts
! Called once per year
! CALLED FROM: subroutine lpj
!--------------------------------------------------------------------------------
    use precision
    use paramodel, only: numpft
    use timemgr, only: dtime
    implicit none

! INTENT IN VARIABLES:
    integer , INTENT(in) :: n_pft                ! total patches
!   real(r8), INTENT(in) :: wt_patch(n_pft)      ! pft weight relative to grid cell
    real(r8), INTENT(in) :: tmomin20(n_pft)      ! 20-yr running mean of tmomin
    real(r8), INTENT(in) :: agdd20(n_pft)        ! 20-yr running mean of agdd
    real(r8), INTENT(in) :: agddtw(n_pft)        ! accumulated growing degree days above twmax
    real(r8), INTENT(in) :: sla(n_pft)           ! ecophys const - sp. leaf area [m2 leaf g-1 carbon]
    real(r8), INTENT(in) :: crownarea_max(n_pft) ! ecophys const - tree maximum crown area [m2]
    real(r8), INTENT(in) :: lm_sapl(n_pft)       ! ecophys const - leaf mass of sapling
    real(r8), INTENT(in) :: sm_sapl(n_pft)       ! ecophys const - stem mass of sapling
    real(r8), INTENT(in) :: hm_sapl(n_pft)       ! ecophys const - heartwood mass of sapling
    real(r8), INTENT(in) :: rm_sapl(n_pft)       ! ecophys const - root mass of saping
    real(r8), INTENT(in) :: tcmin(n_pft)         ! ecophys const - minimum coldest monthly mean temperature
    real(r8), INTENT(in) :: tcmax(n_pft)         ! ecophys const - maximum coldest monthly mean temperature
    real(r8), INTENT(in) :: gddmin(n_pft)        ! ecophys const - minimum growing degree days (at or above 5 C)
    real(r8), INTENT(in) :: prec365              ! yearly running precipitation [mm/s]
    real(r8), INTENT(in) :: vegclass(n_pft)      ! 1.tree 2.shrub 3.grass 4.crop

! INTENT INOUT VARIABLES:
    integer , INTENT(inout) :: ivt(n_pft)           ! vegetation type for this pft
    real(r8), INTENT(inout) :: ifpre(n_pft)         ! true=> PFT present in patch
    real(r8), INTENT(inout) :: nind(n_pft)          ! number of individuals (#/m**2)
    real(r8), INTENT(inout) :: lm_ind(n_pft)        ! individual leaf mass
    real(r8), INTENT(inout) :: sm_ind(n_pft)        ! individual sapwood mass
    real(r8), INTENT(inout) :: hm_ind(n_pft)        ! individual heartwood mass
    real(r8), INTENT(inout) :: rm_ind(n_pft)        ! individual root mass
    real(r8), INTENT(inout) :: litter_ag(n_pft)     ! above ground litter
    real(r8), INTENT(inout) :: litter_bg(n_pft)     ! below ground litter
    real(r8), INTENT(inout) :: fpcgrid(n_pft)       ! foliar projective cover on gridcell zhq.

! INTENT OUT VARIABLES:
    real(r8), INTENT(out) :: htop(n_pft)          ! canopy top (m)
    real(r8), INTENT(out) :: lai_ind(n_pft)       ! LAI per individual
    real(r8), INTENT(out) :: crownarea(n_pft)     ! area each individual tree takes up (m^2)
    real(r8), INTENT(out) :: aestabc_gcell        ! biomass increment due to establishment (gC/m2 veget'd area)

! OTHER LOCAL VARIABLES:
    integer  :: g,l,p,m                         ! indices
    integer  :: fn, filterg                     ! local gridcell filter for error check
    real(r8), parameter :: reinickerp = 1.6     ! parameter in allometric equation
    real(r8), parameter :: wooddens = 2.0e5     ! wood density (gC/m3)
    real(r8), parameter :: T0 = 273.16          ! temperature(K)
    real(r8), parameter :: PI = 3.14159265358979323846 
    real(r8) :: latosa(numpft)       ! leafarea:sapwood cross-sectional area
    real(r8) :: allom1(numpft)       ! parameters in allometric
    real(r8) :: allom2(numpft)       ! parameters in allometric
    real(r8) :: allom3(numpft)       ! parameters in allometric
    logical  :: grid_present(n_pft)  !true=>pft is present in gridcell
    logical  :: grid_survive(n_pft)  !true=>pft survives in gridcell
    logical  :: grid_estab(n_pft)    !true=>pft is established in grirdcell
    integer  :: ngrass               !counter
    integer  :: npft_estab           !counter
    real(r8) :: fpc_tree_total       !total fractional cover of trees in vegetated portion of gridcell
    real(r8) :: fpc_grass_total      !total fractional cover of grass in vegetated portion of gridcell
    real(r8) :: fpc_total            !old-total fractional vegetated portion of gridcell (without bare ground)
    real(r8) :: fpcsum               !old-total fractional vegetated portion of gridcell (include bare ground)
    real(r8) :: grid_tmomin20        !20 year running mean of minimum monthly temperature
    real(r8) :: grid_agdd20          !20 year running mean of growing degree days
    real(r8) :: grid_agddtw          !growing degree base tw
    real(r8) :: grid_prec365         !yearly running precipitation [mm/s]   !added by zhq. dec23
    real(r8) :: fpc_ind(n_pft)       !individual fpc
    real(r8) :: estab_rate(n_pft)    !establishment rate
    real(r8) :: estab_grid(n_pft)    !establishment rate on grid cell
    real(r8) :: bare_max             !maximum bare soil
    real(r8) :: bare                 !fractional cover of bare soil
    real(r8) :: nind_old             !old number of individuals
    real(r8) :: fpcgridtemp          !temporary
    real(r8) :: sm_ind_temp          !temporary
    real(r8) :: stemdiam             !stem diameter
!
! minimum individual density for persistence of PFT (indiv/m2)
!
    real(r8), parameter :: nind_min = 1.0e-10
!
! minimum precip. for establishment (mm/s)
!
    real(r8), parameter :: prec_min_estab = 100./(365.*86400)
!
! maximum sapling establishment rate (indiv/m2)
!
    real(r8), parameter :: estab_max = 0.24
!   real(r8), parameter :: estab_max = 0.12
!-----------------------------------------------------------------------
    ! **********************************************************************
    ! Zeng's version of LPJ's subr. bioclim
    ! Limits based on 20-year running averages of coldest-month mean
    ! temperature and growing degree days (5 degree base).
    ! For SURVIVAL, coldest month temperature and GDD should be
    ! For REGENERATION, PFT must be able to survive AND coldest month
    ! temperature should be no higher than a PFT-specific limit.
    ! **********************************************************************

    latosa(1:numpft) = 8.0e3   ! leafarea:sapwood cross-sectional area
    allom1(1:numpft) = 100.0   ! parameters in allometric
    allom2(1:numpft) =  40.0   ! parameters in allometric
    allom3(1:numpft) =   0.5   ! parameters in allometric

    latosa(9:11) = 4000.0        ! 2000.0 in last test
    allom1(9:11) = 200.0
    allom2(9:11) = 10.0

    ! Initialize gridcell-level metrics

    do p = 1,n_pft        !17,including bare soil
       grid_present(p) = .false.
       grid_survive(p) = .false.
       grid_estab(p)   = .false.
    end do

    ngrass = 0
    npft_estab = 0
    fpc_tree_total = 0._r8
    fpc_grass_total = 0._r8
    fpc_total = 0._r8
    fpcsum = 0._r8
    grid_tmomin20 = 0._r8
    grid_agdd20 = 0._r8
    grid_agddtw = 0._r8
    grid_prec365 = 0._r8

    ! Calculate total woody FPC, FPC increment and grass cover (= crown area)

  ! grid_prec365 = prec365 / 365. / 48.
    grid_prec365 = prec365 / 365. / (86400/dtime)
    fpcsum = sum(fpcgrid(1:n_pft))

    do p = 1,n_pft

       ! Calculate the grid-average bioclimate variables for
       ! survival and establishment

       if(fpcgrid(p) <= 1.0E-6) cycle                             !add by zhq. 0812 

       grid_tmomin20 = grid_tmomin20 + tmomin20(p) * fpcgrid(p) / fpcsum
       grid_agdd20   = grid_agdd20   + agdd20(p)   * fpcgrid(p) / fpcsum
       grid_agddtw   = grid_agddtw   + agddtw(p)   * fpcgrid(p) / fpcsum

    end do

!print *, 'estab', grid_tmomin20,grid_agdd20,grid_agddtw,fpcsum

    ! Must go thru all 16 pfts and decide which can/cannot establish or survive
    ! Determine present, survive, estab.  Note - if tmomin20 > tcmax then  crops,
    ! shrubs and 2nd boreal summergreen tree cannot exist yet (see EcosystemDynini)
    ! because they often coexist using up all pft patches and allowing for no bare
    ! ground. Thus they make fpc_grid_total < 1. Solve by allowing one more pft
    ! per grid cell than the number of pfts.

    do p = 1,numpft 
       ! Set the presence of pft for this gridcell
       ! Note: this modifies the pft-level vegetation type if present is false

       if (ifpre(p) .gt. 0.) then
          grid_present(ivt(p)) = .true.
       else
          ivt(p) = 17
       end if

       if (grid_tmomin20 >= tcmin(p) + T0 ) then
          if (grid_tmomin20 <= tcmax(p) + T0  .and. &
              grid_agdd20 >= gddmin(p) .and. nint(grid_agddtw) == 0) then
             grid_estab(p) = .true.          
          end if

          grid_survive(p) = .true.       
       end if
    end do

    do p = 1, numpft

       ! Case 1 -- pft ceases to exist -kill pfts not adapted to current climate

       if (ifpre(p) .gt. 0. .and. (.not. grid_survive(p) .or. nind(p)<nind_min)) then
          ifpre(p) = -1. !PFT do not exist in this patch
          grid_present(ivt(p)) = .false.
          
          ivt(p) = 17

          ! Add killed biomass to litter
          if (vegclass(p) .le. 2.) then !woody
             litter_ag(p) = litter_ag(p) + nind(p) * (lm_ind(p) + sm_ind(p) + hm_ind(p))
          else              !herbaceous
             litter_ag(p) = litter_ag(p) + nind(p) * lm_ind(p)
          end if
          litter_bg(p) = litter_bg(p) + nind(p) * rm_ind(p)
          lm_ind(p) = 0.0
          sm_ind(p) = 0.0
          rm_ind(p) = 0.0
          hm_ind(p) = 0.0
          fpcgrid(p) = 0.0
       end if

       ! Case 2 -- pft begins to exist - introduce newly "adapted" pfts
!       if (ivt(p) == 17) then  !bare soil
!          if (ifpre(p) .lt. 0. .and. grid_prec365 >= prec_min_estab) then
!             do m = 1, numpft  !excluding bare soil
!               ! if (m /= 17 .and. (.not. grid_present(m)) .and. grid_estab(m)) then
!                if (ivt(p).ne.m .and.ifpre(p).lt.0.) then 
!                  if (.not. grid_present(m) .and. grid_estab(m)) then     ! revised by zhq.
!                    ifpre(p) = 1.
!                    grid_present(m) = .true.
!                    ivt(p) = m
!                    if (vegclass(p) .le. 2.) then
!                       nind(p) = 0.0
!                    else
!                       nind(p) = 1.0 !each grass PFT = 1 "individual"
!                    end if
!                    lm_ind(p) = 0.0
!                    sm_ind(p) = 0.0
!                    rm_ind(p) = 0.0
!                    hm_ind(p) = 0.0
!                    fpcgrid(p) = 0.0
!                    if (vegclass(ivt(p)) .gt. 2) crownarea(p) = 1 !add by zhq
!
!                  end if                                                   ! revised by zhq.
!                end if   ! conditions suitable for establishment
!             end do   ! numpft
!          end if   !  patch conditions for establishment
!       end if   ! if soil


       ! Case 2 -- pft begins to exist - introduce newly "adapted" pfts. revised by zhq.
          if (ifpre(p) .lt. 0. .and. grid_prec365 >= prec_min_estab) then
             if (.not. grid_present(p) .and. grid_estab(p)) then     
                    ifpre(p) = 1.
                    grid_present(p) = .true.
                    ivt(p) =p 
                    if (vegclass(p) .le. 2.) then
                       nind(p) = 0.0
                    else
                       nind(p) = 1.0 !each grass PFT = 1 "individual"
                    end if
                    lm_ind(p) = 0.0
                    sm_ind(p) = 0.0
                    rm_ind(p) = 0.0
                    hm_ind(p) = 0.0
                    fpcgrid(p) = 0.0
                  ! if (vegclass(ivt(p)) .gt. 2) crownarea(p) = 1    ! attention!

             end if                                                   
          end if   !  patch conditions for establishment

       ! Case 3 -- some pfts continue to exist (no change) and some pfts
       ! continue to not exist (no change). Do nothing for this case.

    end do

    ! Sapling and grass establishment
    ! Calculate total woody FPC and number of woody PFTs present and able to establish

    do p = 1, numpft
       if (ifpre(p).gt. 0.) then
          if (vegclass(p) .le. 2.) then
             fpc_tree_total = fpc_tree_total + fpcgrid(p)
             if (grid_estab(ivt(p))) npft_estab = npft_estab + 1
          else if (vegclass(p) .gt. 2.) then !grass
             ngrass = ngrass + 1
             fpc_grass_total = fpc_grass_total + fpcgrid(p)
          end if
          fpc_total = fpc_total + fpcgrid(p)
       end if
    end do

    ! Above establishment counters at the grid level required for the next steps.
    ! Note that ngrass, npft_estab, fpc_tree_total, fpc_grass_total, and fpc_total
    ! complete for vegetated area.

    aestabc_gcell = 0.

    do p = 1, n_pft

       ! Prohibit establishment under extreme temperature or water stress.

       if (grid_prec365 >= prec_min_estab .and. npft_estab > 0) then

          ! Calculate establishment rate over available space, per tree PFT
          ! Maximum establishment rate reduced by shading as tree FPC approaches 1
          ! Total establishment rate partitioned equally among regenerating woody PFTs

          if (fpc_tree_total < 1.0) then    ! new test at Aug 5, 2008, 2:43pm, try to avoid negetive nind

            estab_rate(p) = estab_max * (1.0-exp(5.0*(fpc_tree_total-1.0)))/real(npft_estab)

          ! Calculate grid-level establishment rate per woody PFT
          ! Space available for woody PFT establishment is proportion of grid cell
          ! not currently occupied by woody PFTs

            estab_grid(p) = estab_rate(p) * (1.0-fpc_tree_total)
          else
            estab_grid(p) = 0.0
          endif

       else ! if unsuitable climate for establishment

          estab_grid(p) = 0.0

       end if

       if (ifpre(p) .gt. 0. .and. grid_estab(ivt(p))) then

          if(vegclass(p) .le. 2.) then
             ! Add new saplings to current population

             nind_old = nind(p)
             nind(p) = nind_old + estab_grid(p)

             ! Avoid unbelievable value when nind<nind_min, add by zhq @ 09/14/2010
             if(nind(p) < nind_min) cycle
             
             lm_ind(p)   = (lm_ind(p) * nind_old + lm_sapl(ivt(p)) * estab_grid(p)) / nind(p)
             sm_ind_temp = (sm_ind(p) * nind_old + sm_sapl(ivt(p)) * estab_grid(p)) / nind(p)
             hm_ind(p)   = (hm_ind(p) * nind_old + hm_sapl(ivt(p)) * estab_grid(p)) / nind(p)
             rm_ind(p)   = (rm_ind(p) * nind_old + rm_sapl(ivt(p)) * estab_grid(p)) / nind(p)

             ! Calculate height, diameter and crown area for new average
             ! individual such that the basic allometric relationships (A-C below)
             ! are satisfied.
             ! (A) (leaf area) = latosa * (sapwood xs area)
             !        (Pipe Model, Shinozaki et al. 1964a,b; Waring et al 1982)
             ! (B) (leaf mass) = lmtorm * (root mass)
             ! (C) height = allom2 * (stem diameter)**allom3  (source?)
             ! (D) (crown area) = min (allom1 * (stem diameter)**reinickerp,
             !                                 crownarea_max)
             ! From (A),
             !  (1) sap_xsa = lm_ind * sla / latosa
             !  (2) wooddens = (sm_ind + hm_ind) / stemvolume
             !  (3) stemvolume = stem_xsa * height
             ! From (1), (2) & (3),
             !  (4) stem_xsa = (sm_ind + hm_ind) / wooddens / height
             !  (5) stem_xsa = PI * (stemdiam**2) / 4
             ! From (5),
             !  (6) stemdiam = ( 4 * stem_xsa / PI )**0.5
             ! From (4) & (6),
             !  (7) stemdiam = ( 4 * (sm_ind + hm_ind) / wooddens / height /
             !        PI )**0.5
             ! From (C) & (7),
             !  (8) stemdiam = ( 4 * (sm_ind + hm_ind) / wooddens /
             !        ( allom2 * stemdiam**allom3 ) / PI )**0.5
             ! From (8),
             !  (9) stemdiam = ( 4 * (sm_ind + hm_ind ) / wooddens / PI /
             !        allom2 )**( 1 / (2 + allom3) )

             stemdiam = (4.0*(sm_ind_temp + hm_ind(p))/wooddens/PI/allom2(ivt(p)))**(1.0/(2.0+allom3(ivt(p)))) ! Eqn 9
             htop(p) = allom2(ivt(p)) * stemdiam**allom3(ivt(p))                       ! Eqn C
             crownarea(p) = min(crownarea_max(p), allom1(ivt(p))*stemdiam**reinickerp) ! Eqn D

             ! Recalculate sapwood mass, transferring excess sapwood to heartwood
             ! compartment, if necessary to satisfy Eqn A

             sm_ind(p) = lm_ind(p) * htop(p) * wooddens * sla(p) / latosa(ivt(p))
             hm_ind(p) = hm_ind(p) + (sm_ind_temp-sm_ind(p))

             ! Accumulate biomass increment due to sapling establishment

             aestabc_gcell=aestabc_gcell+(lm_sapl(p)+sm_sapl(p)+hm_sapl(p)+rm_sapl(p))*estab_grid(p)

             ! Update LAI and FPC
             if (crownarea(p) > 0.0) then
                lai_ind(p) = lm_ind(p) * sla(p) / crownarea(p)
             else
                lai_ind(p) = 0.0
             end if
          !  fpc_ind(p) = 1.0 - exp(-0.5*lai_ind(p))
             fpcgrid(p) = crownarea(p) * nind(p)            ! zhq. Adapting X.D.Z's definition of FPC
          end if  
       end if   ! add new saplings block

    end do   ! close loop to update fpc_total_new

!    ! Adjustments- don't allow trees to exceed 95% of vegetated landunit 
!    ! added by zhq. dec24.08
! 
!    do p = 1, n_pft
!       if (fpc_tree_total > 0.95) then
!          if (vegclass(p).le.2. .and. ifpre(p).gt.0.) then
!             nind_old = nind(p)
!             nind(p) = nind(p) / (fpc_tree_total/0.95)
!             fpcgrid(p) = fpcgrid(p) / (fpc_tree_total/0.95)
!             litter_ag(p) = litter_ag(p) + (nind_old - nind(p)) * (lm_ind(p) + sm_ind(p) + hm_ind(p))
!             litter_bg(p) = litter_bg(p) + (nind_old - nind(p)) * rm_ind(p)
!          end if
!          fpc_total=fpc_total-(fpc_tree_total-0.95)
!          fpc_tree_total=0.95
!       end if
!    end do                                   
!print*, 'litter_ag-1',litter_ag(1:n_pft)
    ! Section for grasses. Grasses can establish in non-vegetated areas

    do p = 1, n_pft
       if (ifpre(p).gt. 0. .and. vegclass(p).gt.2.) then  !herbaceous
          bare = 0.

          if (grid_estab(p)) then
             if (ngrass > 0) then
                bare = (1.0 - fpc_total) / real(ngrass)
             else
                bare = 0.0
             end if
             bare_max = ((-2.0 * log(max(1.0 - bare - fpcgrid(p), 0.000001_r8))) &
                  / sla(ivt(p)) - lm_ind(p)) / lm_sapl(ivt(p))
             bare = max(0.0_r8, min(bare, bare_max))
             lm_ind(p) = lm_ind(p) + bare * lm_sapl(ivt(p))
             rm_ind(p) = rm_ind(p) + bare * rm_sapl(ivt(p))
             crownarea(p) = 1.0 - exp(-0.5 * lm_ind(p) * sla(p))        ! added by X.D.Z
             ! Update LAI and FPC
             if (crownarea(p) > 0.0) then
                lai_ind(p) = lm_ind(p) * sla(p) / crownarea(p)
             else
                lai_ind(p) = 0.0
             end if

!            lai_ind(p) = lm_ind(p) * sla(p) / crownarea(p)             ! added by X.D.Z
             fpcgrid(p) = crownarea(p) * nind(p)                        ! added by X.D.Z
!if(p==14) print*,'estab2',crownarea(p),lai_ind(p)
          end if

        ! Accumulate biomass increment due to grass establishment

          aestabc_gcell=aestabc_gcell+bare*(lm_sapl(p)+rm_sapl(p))  ! zhq.

          if (lm_ind(p) <= 0.0) then
             ifpre(p) = -1.
             litter_bg(p) = litter_bg(p) + rm_ind(p) * nind(p)
          end if
       end if
    end do   ! end of pft-loop

!print*, 'litter_ag-2',litter_ag(1:n_pft)
    ! Recalculate fpc_total after grass establishment
    ! added by zhq. jan7,09

       fpc_tree_total=0.0
       fpc_grass_total=0.0
       fpc_total=0.0

    do p = 1, n_pft
       if (ifpre(p).gt. 0.) then
          if (vegclass(p) .le. 2.) then
             fpc_tree_total = fpc_tree_total + fpcgrid(p)
          else if (vegclass(p) .gt. 2.) then !grass
             fpc_grass_total = fpc_grass_total + fpcgrid(p)
          end if
          fpc_total = fpc_total + fpcgrid(p)
       end if
    enddo

  end subroutine Establishment

!=======================================================================================

  subroutine Light(n_pft, num_natvegp, filter_natvegp, ivt,&
                   fpc_inc, sm_ind, hm_ind, lai_ind, crownarea,&
                   sla, vegclass, fpcgrid, nind, litter_ag,&
                   litter_bg, lm_ind, rm_ind, crownarea_max)

!-------------------------------------------------------------------------
! Calculate light competition
! Update fpc 
! Called once per year
! CALLED FROM: subroutine lpj
!-------------------------------------------------------------------------

    use precision
    implicit none

! INTENT IN VARIABLES:
    integer, INTENT(in) :: n_pft             ! total patches
    integer, INTENT(in) :: num_natvegp       ! number of vegetated pfts in grid
    integer, INTENT(in) :: filter_natvegp(n_pft) ! naturally-vegetated points
    integer, INTENT(in) :: ivt(n_pft)        ! pft vegetation type
    real(r8), INTENT(in) :: fpc_inc(n_pft)   ! foliar projective cover increment (fraction)
    real(r8), INTENT(in) :: sm_ind(n_pft)    ! individual stem mass
    real(r8), INTENT(in) :: hm_ind(n_pft)    ! individual heartwood mass
    real(r8), INTENT(in) :: sla(n_pft)       ! specific leaf area [m2 leaf g-1 carbon]
    real(r8), INTENT(in) :: vegclass(n_pft)  ! 1.tree 2.shrub 3.grass 4.crop
    real(r8), INTENT(in) :: crownarea_max(n_pft) ! ecophys const - tree maximum crown area [m2]

! INTENT INOUT VARIABLES:
    real(r8), INTENT(inout) :: fpcgrid(n_pft)   ! foliar projective cover on gridcell 
    real(r8), INTENT(inout) :: nind(n_pft)      ! number of individuals
    real(r8), INTENT(inout) :: litter_ag(n_pft) ! above ground litter
    real(r8), INTENT(inout) :: litter_bg(n_pft) ! below ground litter
    real(r8), INTENT(inout) :: lm_ind(n_pft)    ! individual leaf mass
    real(r8), INTENT(inout) :: rm_ind(n_pft)    ! individual root mass
    real(r8), INTENT(inout) :: lai_ind(n_pft)   ! LAI per individual, added by X.D.Z
    real(r8), INTENT(inout) :: crownarea(n_pft) ! area each individual tree takes up (m^2)

! OTHER LOCAL VARIABLES:
    real(r8), parameter :: fpc_tree_max = 0.95  ! maximum total tree FPC

    real(r8) :: fpc_tree_total
    real(r8) :: fpc_inc_tree
    real(r8) :: fpc_grass_total
    real(r8) :: fpc_shrub_total        ! added by X.D.Z
    real(r8) :: fpc_grass_max          ! added by X.D.Z
    real(r8) :: fpc_shrub_max          ! added by X.D.Z
    integer  :: p, fp, vc              ! indices
    integer  :: ntree
    real(r8) :: excess
    real(r8) :: nind_kill
    real(r8) :: lm_old
    real(r8) :: lm_kill
    real(r8) :: rm_kill

    fpc_tree_total = 0.
    fpc_inc_tree = 0.
    fpc_grass_total = 0.
    fpc_shrub_total = 0.         ! added by X.D.Z
    ntree = 0

    do fp = 1,num_natvegp
       p = filter_natvegp(fp)
       vc = nint(vegclass(p))

       if (vc == 1) then            ! tree
          ntree = ntree + 1
          fpc_tree_total = fpc_tree_total + fpcgrid(p)
          fpc_inc_tree = fpc_inc_tree + fpc_inc(p)
       else if (vc == 2) then       ! shrub
          fpc_shrub_total = fpc_shrub_total + fpcgrid(p)
       else                         ! grass
          fpc_grass_total = fpc_grass_total + fpcgrid(p)
       end if
    end do
    fpc_grass_max = 1.0 - min(fpc_tree_total, fpc_tree_max)
    fpc_shrub_max = max(0.0, fpc_grass_max - fpc_grass_total) 

    do fp = 1,num_natvegp
       p = filter_natvegp(fp)
       vc = nint(vegclass(p))

       ! light competition
       if (vc == 1) then                ! tree

          if (fpc_tree_total > fpc_tree_max) then

             if (fpc_inc_tree > 0.0) then
                excess = (fpc_tree_total - fpc_tree_max) * &
                          fpc_inc(p) / fpc_inc_tree
             else
              ! excess = (fpc_tree_total - fpc_tree_max) / &
              !           real(ntree)
                excess = (fpc_tree_total - fpc_tree_max) * &      ! change by X.D.Z to avoid negative fpc
                          fpcgrid(p) / fpc_tree_total
             end if

             ! Reduce individual density (and thereby gridcell-level biomass)
             ! so that total tree FPC reduced to 'fpc_tree_max'

             nind_kill = nind(p) * excess / fpcgrid(p)
             nind(p) = nind(p) - nind_kill

             ! Transfer lost biomass to litter

             litter_ag(p) = litter_ag(p) + nind_kill * (lm_ind(p) + sm_ind(p) + hm_ind(p))
             litter_bg(p) = litter_bg(p) + nind_kill * rm_ind(p)
             fpcgrid(p) = crownarea(p) * nind(p)        

          end if   ! tree
          
       else if (vc > 2) then        

          if (fpc_grass_total > fpc_grass_max) then    ! added by X.D.Z

             ! grass competes with itself if total fpc exceeds 1 (**add comments**)

             excess = (fpc_grass_total - fpc_grass_max) * fpcgrid(p) / fpc_grass_total  
             lm_old = lm_ind(p)
             lm_ind(p) = -2.0 * log(1.0-(fpcgrid(p) - excess)) / sla(p)
             lm_kill = lm_old - lm_ind(p)
             rm_kill = rm_ind(p) * lm_kill/lm_old
             rm_ind(p) = rm_ind(p) - rm_kill

             ! Transfer lost biomass to litter

             litter_ag(p) = litter_ag(p) + lm_kill 
             litter_bg(p) = litter_bg(p) + rm_kill

             crownarea(p) = 1.0 - exp(-0.5 * lm_ind(p) * sla(p))   ! added by X.D.Z
             lai_ind(p) = lm_ind(p) * sla(p) / crownarea(p)        ! added by X.D.Z
             fpcgrid(p) = crownarea(p)                             ! added by X.D.Z 

          end if  ! grass
          
       else                              ! shrub
          if (fpc_shrub_total > fpc_shrub_max) then

             excess = 1.0 - fpc_shrub_max / fpc_shrub_total

             ! Reduce individual density (and thereby gridcell-level biomass)
             ! so that total shrub FPC reduced to 'fpc_shrub_max'

             nind_kill = nind(p) * excess
             nind(p) = nind(p) - nind_kill

             ! Transfer lost biomass to litter

             litter_ag(p) = litter_ag(p) + nind_kill * (lm_ind(p) + sm_ind(p) + hm_ind(p))
             litter_bg(p) = litter_bg(p) + nind_kill * rm_ind(p)
             fpcgrid(p) = crownarea(p) * nind(p) 

          end if ! shrub

       end if   ! end of tree, grass, shrub

    end do

  end subroutine Light

!===========================================================================

  subroutine ResetBareFPC(n_pft, ivt, fpcgrid, lai_ind, ifpre)
!----------------------------------------------------------
! Calculates baresoil FPC
! CALLED FROM: subroutine lpj
!----------------------------------------------------------
    use precision
    implicit none

! INTENT IN VARIABLES:
    integer , INTENT(in) :: n_pft           ! total patches

! INTENT INOUT VARIABLES:
    integer , INTENT(inout) :: ivt(n_pft)        ! vegetation type for this pft
    real(r8), INTENT(inout) :: fpcgrid(n_pft)    ! foliar projective cover on gridcell 
    real(r8), INTENT(inout) :: lai_ind(n_pft)    ! LAI per individual
    real(r8), INTENT(inout) :: ifpre(n_pft)      ! true=> PFT present in patch

! OTHER LOCAL VARIABLES:
    integer  :: p                          ! indices
    real(r8) :: fpc_total                  ! total fractional vegetated portion of gridcell

    ! Initialize gridcell-level metrics
    fpc_total = 0._r8

    do p = 1, n_pft
       ! Avoid fpcgrid=0 with ifpre=1.

       if (fpcgrid(p) == 0.) then
          ifpre(p) = -1.
          ivt(p) = 17
       endif

       ! Avoid ivt(p)>16 with lai>0.

       if (ifpre(p) .lt. 0. .and. ivt(p) .gt. 16) then
          fpcgrid(p) = 0.0
          lai_ind(p) = 0.0
       end if
    end do

    ! confine fpc_total .le. 1. added by zhq. dec25,08
 
    fpc_total=sum(fpcgrid(1:16))    ! except bare soil

    if ((fpc_total-1.0).gt.1.0e-6) then
        print*,'error: fpc_total>1',fpc_total,fpcgrid
        stop
    end if 

    ! reset bare soil

    fpcgrid(n_pft) = max(0.0, 1.0-fpc_total)

    do p =1, n_pft
       if (fpcgrid(p)>=1.0E-6) ifpre(p)=1.   ! include bare soil pft
    end do

  end subroutine ResetBareFPC
