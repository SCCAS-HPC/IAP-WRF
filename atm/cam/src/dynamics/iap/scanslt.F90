module scanslt
!-----------------------------------------------------------------------
!
! Module to handle Semi-Lagrangian transport in the context of
! Eulerian Spectral dynamics.
! Modified: zhanghe, 2008.4
! Update: Juanxiong he, 2010.08	
! Reviewed: ZhangHe, 2011-11-19 
! Modified: ZhangHe, 2012-01-15, copy qminus to adv_state%qminus          
!         : Jiang Jinrong, October 2012, for 2D parallel
!         : Zhang He, 2013-01-25, new mp_send3d, dyn_state was added
!         : Zhang He, 2013-03-12, removed Q, Qliq and Qice
!         : Zhang He, 2013-03-15
!         : Zhang He, 2013-03-21, n3m2 --> n3m1
!         : Zhang He, 2013-04-12, modified ad_coupling and da_coupling
!-----------------------------------------------------------------------
!
! $Id$
!
!-----------------------------------------------------------------------
   use shr_kind_mod, only: r8 => shr_kind_r8
   use pmgrid,       only: plon, plat, plev, beglat, endlat, plevp, beglatdyn, endlatdyn, &
                           beglev, endlev, endlevp, npr_y, beglatdynex
   use Dyn_const,  only: SIG
   use infnan,       only: bigint
   use constituents, only: pcnst
   use abortutils,   only: endrun
   use scamMod,      only: single_column
   use perf_mod
   use dynamics_vars,   only: T_FVDYCORE_GRID, T_FVDYCORE_STATE      !zhh 2013-01-25
   use dyn_internal_state,   only : get_dyn_state      !zhh 2013-01-25
#if (defined SPMD)
   use mod_comm, only: mp_send3d, mp_recv3d, mp_sendirr, mp_recvirr 
   use spmd_dyn
#endif
!-----------------------------------------------------------------------
   implicit none
!-----------------------------------------------------------------------
   private
!
! Public interfaces
!
   public scanslt_initial ! Advection initialization method
   public scanslt_run     ! Advection run method
   public scanslt_final   ! Advection finalization method
   public scanslt_alloc   ! Allocate some slt data needed for restarting
!
! Public extended grid parameters
!
   integer, public, parameter :: nxpt   = 1                        ! no. of pts outside active domain of interpolant
   integer, public, parameter :: nxptj  = 1                        ! zhh 2008.4.15   
   integer, public, parameter :: jintmx = 2                        ! number of extra latitudes in polar region
   integer, public, parameter :: i1     = 1 + nxpt                 ! model starting longitude index
   integer, public, parameter :: j1     = 1 + nxptj + jintmx       ! zhh 2008.4.15   
   integer, public, parameter :: plond  = plon + 1 + 2*nxpt        ! slt extended domain longitude
   integer, public, parameter :: plond1 = plond - i1 +1            ! slt extended domain longitude starting at i1
   integer, public, parameter :: platd  = plat + 2*nxptj + 2*jintmx ! zhh 2008.4.15
   integer, public, parameter :: numbnd = nxptj + jintmx           ! zhh 2008.4.15
   integer, public, parameter :: plndlv = plond*plev               ! Length of multilevel 3-d field slice

   integer, public :: beglatex   ! extended grid beglat
   integer, public :: endlatex   ! extended grid endlat

#if ( ! defined SPMD )
   parameter (beglatex = 1)
   parameter (endlatex = platd)
#endif

   public engy1lat    ! For calculation of total energy
   public hw1lat      ! For calculation of total moisture
!
! Public data structures
!
   public advection_state

   ! advection data structure of data that will be on the extended grid for SLT
   type advection_state
     real(r8), pointer :: u3(:,:,:)        ! u-wind
     real(r8), pointer :: v3(:,:,:)        ! v-wind
     real(r8), pointer :: qminus(:,:,:,:)  ! constituents on previous step
   end type advection_state

   public lammp, phimp, sigmp, qfcst       ! Needed for restart
!
   integer, public :: nlonex(platd) = bigint ! num longitudes per lat (extended grid)
   real(r8) :: hw1lat (pcnst,plat)           ! lat contribution to const. mass integral
   real(r8) :: engy1lat(plat)                ! lat contribution to total energy integral
   real(r8), allocatable :: lammp(:,:,:)   ! Lamda midpoint coordinate
   real(r8), allocatable :: phimp(:,:,:)   ! Phi midpoint coordinate
   real(r8), allocatable :: sigmp(:,:,:)   ! Sigma midpoint coordinate
   real(r8), allocatable :: qfcst(:,:,:,:) ! slt forecast of moisture and constituents
!
! Private data
!
   integer, parameter  :: pmap = 20000 
!                    ! max dimension of evenly spaced vert. 
!                    ! grid used by SLT code to map the departure pts into true
!                    ! model levels.
!
   real(r8) :: etaint(plevp)               ! vertical coords at interfaces
   real(r8) :: dlam(platd)                 ! longitudinal grid interval (radians)
   real(r8) :: lam(plond,platd)            ! longitude coords of extended grid
   real(r8) :: phi(platd)                  ! latitude  coords of extended grid
   real(r8) :: dphi(platd)                 ! latitude intervals (radians)
   real(r8) :: sinlam(plond,platd)         ! sin(lam) model domain only
   real(r8) :: coslam(plond,platd)         ! cos(lam) model domain only
   real(r8) :: lbasdy(4,2,platd)           ! latitude derivative weights
   real(r8) :: lbasdz(4,2,plev)            ! vert (full levels) deriv wghts
   real(r8) :: lbassd(4,2,plevp)           ! vert (half levels) deriv wghts
   real(r8) :: lbasiy(4,2,platd)           ! Lagrange cubic interp wghts (lat.)
   real(r8) :: detai(plevp)                ! intervals between vert half levs.
   integer  :: kdpmpf(pmap)                ! artificial full vert grid indices
   integer  :: kdpmph(pmap)                ! artificial half vert grid indices
   real(r8) :: gravit                      ! gravitational constant

!-----------------------------------------------------------------------
contains
!-----------------------------------------------------------------------

!
!-----------------------------------------------------------------------
!

subroutine scanslt_alloc()
!----------------------------------------------------------------------- 
! 
! Purpose: 
! Allocate some scanslt data
! 
! Author: 
!
! Erik Kluzek
! Jiang Jinrong
!----------------------------------------------------------------------- 
   use infnan,       only: nan

   allocate (lammp(plon,beglev:endlev,beglat:endlat))
   allocate (phimp(plon,beglev:endlev,beglat:endlat))
   allocate (sigmp(plon,beglev:endlev,beglat:endlat))
   allocate (qfcst(plon,beglev:endlev,pcnst,beglat:endlat))

   lammp (:,:,:)   = nan
   phimp (:,:,:)   = nan
   sigmp (:,:,:)   = nan
   qfcst (:,:,:,:) = nan
end subroutine scanslt_alloc

!
!-----------------------------------------------------------------------
!
subroutine scanslt_initial( adv_state, etamid, gravit_in, gw, detam, cwava )
!----------------------------------------------------------------------- 
! 
! Purpose: 
! SLT initialization for Eulerian dynamics
! 
! Author: 
!
! Erik Kluzek
!
!----------------------------------------------------------------------- 
   use commap,       only: clat
   use prognostics,  only: ps, n3
   use rgrid,        only: nlon
   use time_manager, only: is_first_step
   use eul_control_mod, only : pdela
!
! Input arguments
!
   real(r8), intent(inout) :: etamid(plev)  ! vertical coords at midpoints
   real(r8), intent(in) :: gravit_in        ! Gravitational constant
!
! Output arguments
!
   real(r8), intent(out) :: gw(plat)               ! Gaussian weights
   real(r8), intent(out) :: detam(plev)            ! intervals between vert full levs.
   real(r8), intent(out) :: cwava(plat)            ! weight applied to global integrals
   type(advection_state), intent(out) :: adv_state ! Advection state data

!
!  Local variables
!
   integer :: i, j, k, lat        ! indices
   real(r8) :: hyad (plev)        ! del (A)
   real(r8) :: pmid(plon,plev)    ! pressure at model levels
   real(r8) :: pint(plon,plevp)   ! pressure at interfaces
   real(r8) :: pdel(plon,plev)    ! pressure difference between
!-------------------------------------------------------------------------------
!
! Allocate memory for scanslt variables
   beglatex=beglat
   endlatex=endlat+2*numbnd
!
   call adv_state_alloc( adv_state )
!
! Eta at interfaces
!
   do k=1,plevp
      etaint(k) = SIG(k)
   end do
!
! For SCAM compute pressure levels to use for eta interface
!
   if (single_column) then
      lat = beglat
      call plevs0(nlon(lat), plon, plev, ps(1,lat,n3), pint, pmid, pdel)
      etamid(:) = pmid(lat,:)
      etaint(:) = pint(lat,:)
      if ( any(etamid == 0.0_r8) ) call endrun('etamid == 0')
      if ( any(etaint == 0.0_r8) ) call endrun('etaint == 0')
   endif
!
! Set slt module variables
!
   gravit = gravit_in
   call grdini(pmap    ,etamid  ,etaint  ,gravit  ,dlam    , &
               lam     ,phi     ,dphi    ,gw      ,sinlam  , &
               coslam  ,lbasdy  ,lbasdz  ,lbassd  ,lbasiy  , &
               detam   ,detai   ,kdpmpf  ,kdpmph  ,cwava   )
!
! Initial guess for trajectory midpoints in spherical coords.
! nstep = 0:  use arrival points as initial guess for trajectory midpoints.
! nstep > 0:  use calculated trajectory midpoints from previous time
! step as first guess.
! NOTE:  reduce number of iters necessary for convergence after nstep = 1.
!
   if (is_first_step()) then
      do lat=beglat,endlat
         j = j1 - 1 + lat
!
! Set current time pressure arrays for model levels etc.
!
         call plevs0(nlon(lat), plon, plev, ps(1,lat,n3), pint, pmid, pdel)

         do k=beglev,endlev !jjr
            do i=1,nlon(lat)
               if (single_column) then
                  sigmp(i,k,lat) = pmid(i,k)
               else
                  lammp(i,k,lat) = real(i-1,r8)*dlam(j1-1+lat)
                  phimp(i,k,lat) = clat(lat)
                  sigmp(i,k,lat) = etamid(k)
               endif
            end do
         end do
      end do
   end if
!
end subroutine scanslt_initial

!
!-----------------------------------------------------------------------
!
subroutine scanslt_run(adv_state, ztodt, detam, etamid, cwava )   !zhh 2008.1.8
!----------------------------------------------------------------------- 
! 
! Purpose: 
! Driving routine for semi-lagrangian transport.
! 
! Method: 
! The latitude loop in this routine is multitasked.
! 
! Author: 
! Original version:  J. Rosinski
! Standardized:      J. Rosinski, June 1992
! Reviewed:          D. Williamson, P. Rasch, August 1992
! Reviewed:          D. Williamson, P. Rasch,  March 1996
!
!-----------------------------------------------------------------------
   use rgrid,       only: nlon
   use physconst,   only: ra
   use prognostics, only: hadv, qminus, q3, n3, n3m1
   use time_manager, only: get_nstep
   use pmgrid,       only: plon, plat
   use IAP_prog,    only: Q, Qliq, Qice   !debug
#if (defined SPMD)
   use mpishorthand, only: mpicom
#endif
   use spmd_utils, only: masterproc, iam    ! zhh for debug, 2013-02-18
!------------------------------Parameters-------------------------------
   integer itermx  ! number of iterations to be used in departure
!                     ! point calculation for nstep = 0 and 1
   integer itermn  ! number of iterations to be used in departure
!                     ! point calculation for nstep > 1
   parameter(itermx=4,itermn=1)
!-----------------------------------------------------------------------
!
! Input arguments
!
   real(r8), intent(in) :: ztodt              ! twice the time step unless nstep = 0
   real(r8), intent(in) :: etamid(plev)       ! eta at levels
!
! In/Output arguments
!
   real(r8), intent(inout) :: detam(plev)     ! delta eta at levels 
                                              ! needs intent(out) because of SCAM
   real(r8), intent(inout) :: cwava(plat)     ! weight for global water vapor int.
                                              ! needs intent(out) because of SCAM
   type(advection_state), intent(inout) :: adv_state ! Advection state data
!
!---------------------------Local workspace-----------------------------
!
   integer iter                ! number of iterations for
!                                 ! departure point calculation
   integer m
   integer lat                 ! latitude index
   integer irow                ! N/S latitude pair index
   integer jcen                ! lat index (extended grid) of forecast
   integer :: nstep            ! current timestep number
   real(r8) :: pmid(plon,plev) ! pressure at model levels
   real(r8) :: pint(plon,plevp)! pressure at interfaces
   real(r8) :: pdel(plon,plev) ! pressure difference between
!
! Dynamic (SPMD) vs stack (shared memory)
!
   real(r8) uxl(plond,beglev:endlev,beglatex:endlatex)     ! left  x-deriv of u/v
   real(r8) uxr(plond,beglev:endlev,beglatex:endlatex)     ! left  x-deriv of u/v
   real(r8) vxl(plond,beglev:endlev,beglatex:endlatex)     ! left  x-deriv of u/v
   real(r8) vxr(plond,beglev:endlev,beglatex:endlatex)     ! left  x-deriv of u/v
   real(r8) qxl(plond,beglev:endlev,pcnst,beglatex:endlatex) ! left  x-deriv of constituents
   real(r8) qxr(plond,beglev:endlev,pcnst,beglatex:endlatex) ! right  x-deriv of constituents
   real(r8) :: gw(plat)                             ! Gaussian weights needed for SCAM grdini call
   real(r8) :: etadot(plon,beglev:endlev+1,beglat:endlat)! vertical motion (slt)
   real(r8) :: etadotxy(beglonxy:endlonxy,plevp,beglatxy:endlatxy)
   real(r8):: fhr(plon,beglev:endlev,pcnst,beglat:endlat)
   real(r8):: dummy3(plon,beglev:endlev,beglat:endlat)
   real(r8):: fhrxy(beglonxy:endlonxy,plev,pcnst,beglatxy:endlatxy)
   real(r8):: dummy3xy(beglonxy:endlonxy,plev,beglatxy:endlatxy)
   real(r8):: sigmpxy(beglonxy:endlonxy,plev,beglatxy:endlatxy)
   real(r8):: qfcstxy(beglonxy:endlonxy,plev,pcnst,beglatxy:endlatxy)
   integer :: k,dest,src,l ,i,j                                    ! Vertical index needed for SCAM
!
   type (T_FVDYCORE_STATE), pointer :: dyn_state   !zhh 2013-01-25
   type (T_FVDYCORE_GRID) , pointer :: grid        ! For convenience
   integer  :: commyz, commxy     !zhh 2013-01-25
!-----------------------------------------------------------------------
!
   dyn_state => get_dyn_state()
   commyz = dyn_state%grid%commyz
   grid => dyn_state%grid     ! For convenience
!
! Copy dynamics data into SLT advection structure
!
   call t_startf ('scanslt_da_coup')
   call da_coupling( cwava, etadot, adv_state ) !zhang he  
   call t_stopf ('scanslt_da_coup')
!
! For SCAM reset vertical grid
!
   if (single_column) then
!
!     IF surface pressure changes with time we need to remap the vertical
!     coordinate for the slt advection process.  It has been empirically
!     determined that we can get away with 500 for pmap (instead of 20000)
!     This is necessary to make the procedure computationally feasible
!
      call grdini(pmap    ,etamid  ,etaint  ,gravit  ,dlam    , &
               lam     ,phi     ,dphi    ,gw      ,sinlam  , &
               coslam  ,lbasdy  ,lbasdz  ,lbassd  ,lbasiy  , &
               detam   ,detai   ,kdpmpf  ,kdpmph  ,cwava   )
!     
! Initial guess for trajectory midpoints in spherical coords.
! nstep = 0:  use arrival points as initial guess for trajectory midpoints.
! nstep > 0:  use calculated trajectory midpoints from previous time
! step as first guess.
! NOTE:  reduce number of iters necessary for convergence after nstep = 1.
!
      do k=beglev, endlev
         sigmp(1,k,beglat) = etamid(k)
      end do

   else
!
! Mpi barrier
!
#if ( defined SPMD )
!
      call t_startf ('sync_bndexch')
      call mpibarrier (mpicom)
      call t_stopf ('sync_bndexch')
!
! Communicate boundary information 
!u3,v3,qminus,need test!!!!!
      src = iam+1
      dest  = iam-1
      if ( mod(iam,npr_y) == 0 ) dest = -1
      if ( mod(iam+1,npr_y) == 0 ) src = -1
      call mp_send3d(commyz, dest,src,plond,plev,platd,1,plond,beglev,endlev,beglatex,endlatex,&
                      1,plond,beglev,endlev,beglatex+3,beglatex+5,adv_state%u3)     
      call mp_recv3d(commyz, src,plond,plev,platd,1,plond,beglev,endlev,beglatex,endlatex,&
                      1,plond,beglev,endlev,endlatex-2,endlatex,adv_state%u3)
!
      call mp_send3d(commyz, dest,src,plond,plev,platd,1,plond,beglev,endlev,beglatex,endlatex,&
                      1,plond,beglev,endlev,beglatex+3,beglatex+5,adv_state%v3)
      call mp_recv3d(commyz, src,plond,plev,platd,1,plond,beglev,endlev,beglatex,endlatex,&
                      1,plond,beglev,endlev,endlatex-2,endlatex,adv_state%v3)
!
      call mp_send3d( commyz, dest, src, plond,plev*pcnst , platd,                      &
                      1, plond,beglev,(endlev-beglev+1)*pcnst+beglev-1, beglatex, endlatex,       &
                      1, plond,beglev,(endlev-beglev+1)*pcnst+beglev-1,beglatex+3 , beglatex+5, adv_state%qminus )
      call mp_recv3d( commyz, src, plond,plev*pcnst,platd,                             &
                      1, plond, beglev,(endlev-beglev+1)*pcnst+beglev-1,beglatex, endlatex,       &
                      1, plond, beglev,(endlev-beglev+1)*pcnst+beglev-1,endlatex-2, endlatex,  adv_state%qminus )
!
      src = iam-1
      dest  = iam+1
      if ( mod(iam,npr_y) == 0 ) src = -1
      if ( mod(iam+1,npr_y) == 0 ) dest = -1
      call mp_send3d(commyz, dest,src,plond,plev,platd,1,plond,beglev,endlev,beglatex,endlatex,&
                      1,plond,beglev,endlev,endlatex-5,endlatex-3,adv_state%u3)
      call mp_recv3d(commyz, src,plond,plev,platd,1,plond,beglev,endlev,beglatex,endlatex, &
                      1,plond,beglev,endlev,beglatex,beglatex+2,adv_state%u3)
!
      call mp_send3d(commyz, dest,src,plond,plev,platd,1,plond,beglev,endlev,beglatex,endlatex,&
                      1,plond,beglev,endlev,endlatex-5,endlatex-3,adv_state%v3)
      call mp_recv3d(commyz, src,plond,plev,platd,1,plond,beglev,endlev,beglatex,endlatex, &
                      1,plond,beglev,endlev,beglatex,beglatex+2,adv_state%v3)
!
      call mp_send3d( commyz, dest, src, plond,plev*pcnst , platd,                      &
                      1, plond,beglev,(endlev-beglev+1)*pcnst+beglev-1, beglatex, endlatex,       &
                      1, plond,beglev,(endlev-beglev+1)*pcnst+beglev-1,endlatex-5 , endlatex-3, adv_state%qminus )
      call mp_recv3d( commyz, src, plond,plev*pcnst,platd,                             &
                      1, plond, beglev,(endlev-beglev+1)*pcnst+beglev-1,beglatex, endlatex,       &
                      1, plond, beglev,(endlev-beglev+1)*pcnst+beglev-1,beglatex, beglatex+2,  adv_state%qminus )
!
#endif

      nstep = get_nstep()
!
! Initialize extended arrays
!
      call t_startf('sltini')
      call sltini (dlam,    sinlam,  coslam,  uxl,     uxr, &
           vxl,     vxr,     qxl,     qxr,     adv_state )
      call t_stopf('sltini')
   endif
   nstep = get_nstep()
   if (nstep .le. 1) then
      iter = itermx
   else
      iter = itermn
   end if
!
! Loop through latitudes producing forecast
!
   call t_startf ('sltb1')
#ifdef OUTER_OMP
!wjp 2011.04.16 !$OMP PARALLEL DO PRIVATE (LAT, IROW, JCEN)
#endif
   do lat=beglat,endlat
      jcen = j1 - 1 + lat
!
! Call slt interface routine.
!
      call sltb1 (jcen    ,lat     ,ztodt   ,ra      , &
                  iter    ,uxl     ,uxr     ,vxl     ,vxr     , &
                  qxl     ,qxr     ,lam     , &
                  phi     ,dphi   , &
                  lbasdy  ,lbasiy  , &
                  lammp(:,:,lat), phimp(:,:,lat),  &
                  adv_state, nlon(lat), hadv, nlonex,fhr(:,:,:,lat)  )
   end do

   if( twod_decomp.eq.1) then
#if (defined SPMD)
      do l=1,pcnst
  
         do lat=beglat,endlat
            do k=beglev,endlev
               do i=1,plon
                  dummy3(i,k,lat)=fhr(i,k,l,lat)
               end do
            end do
         end do

         call mp_sendirr(grid%commxy, grid%ikj_yz_to_xy%SendDesc,         &
                         grid%ikj_yz_to_xy%RecvDesc, dummy3, dummy3xy,    &
                         modc=grid%modc_dynrun )
         call mp_recvirr(grid%commxy, grid%ikj_yz_to_xy%SendDesc,         &
                         grid%ikj_yz_to_xy%RecvDesc, dummy3, dummy3xy,    &
                         modc=grid%modc_dynrun )
         do lat=beglatxy,endlatxy
            do k=1,plev
               do i=beglonxy,endlonxy
                  fhrxy(i,k,l,lat)=dummy3xy(i,k,lat)
               end do
            end do
         end do

      end do
      call mp_sendirr(grid%commxy, grid%pe_to_pexy%SendDesc,         &   !zhh 2013-01-25
                      grid%pe_to_pexy%RecvDesc, etadot, etadotxy,    &
                      modc=grid%modc_dynrun )
      call mp_recvirr(grid%commxy, grid%pe_to_pexy%SendDesc,         &   !zhh 2013-01-25
                      grid%pe_to_pexy%RecvDesc, etadot, etadotxy,    &
                      modc=grid%modc_dynrun )
!
      call mp_sendirr(grid%commxy, grid%ikj_yz_to_xy%SendDesc,         &
                      grid%ikj_yz_to_xy%RecvDesc, sigmp, sigmpxy,      &
                      modc=grid%modc_dynrun )
      call mp_recvirr(grid%commxy, grid%ikj_yz_to_xy%SendDesc,         &
                      grid%ikj_yz_to_xy%RecvDesc, sigmp, sigmpxy,      &
                      modc=grid%modc_dynrun )
#endif
   else    
      fhrxy=fhr
      sigmpxy=sigmp
      etadotxy=etadot   
          !jjr need define xy   
   endif
   call t_stopf ('sltb1')
!
   call t_startf ('sltb2')
   do lat=beglatxy,endlatxy
      jcen = j1 - 1 + lat
!
! Call slt interface routine.
!
      call sltb2 (lat,pmap  ,ztodt        , &
                  iter    , &
                  etadotxy(:,:,lat)  ,  &
                  etamid  ,etaint  ,detam   , &
                  detai    ,lbasdz  ,lbassd  , &
                  kdpmpf  ,kdpmph  ,  sigmpxy(:,:,lat), &
                  qfcstxy(:,:,:,lat) , nlon(lat),fhrxy(:,:,:,lat)  )
   end do

   if( twod_decomp.eq.1) then

#if (defined SPMD)
      call mp_sendirr(grid%commxy, grid%ikj_xy_to_yz%SendDesc,         &
                      grid%ikj_xy_to_yz%RecvDesc, sigmpxy, sigmp,      &
                      modc=grid%modc_dynrun )
      call mp_recvirr(grid%commxy, grid%ikj_xy_to_yz%SendDesc,         &
                      grid%ikj_xy_to_yz%RecvDesc, sigmpxy, sigmp,      &
                      modc=grid%modc_dynrun )

      do l=1,pcnst

         do lat=beglatxy,endlatxy
            do k=1,plev
               do i=beglonxy,endlonxy
                  dummy3xy(i,k,lat)=qfcstxy(i,k,l,lat)
               end do
            end do
         end do
!!
         call mp_sendirr(grid%commxy, grid%ikj_xy_to_yz%SendDesc,         &
                         grid%ikj_xy_to_yz%RecvDesc, dummy3xy, dummy3,    &
                         modc=grid%modc_dynrun )
         call mp_recvirr(grid%commxy, grid%ikj_xy_to_yz%SendDesc,         &
                         grid%ikj_xy_to_yz%RecvDesc, dummy3xy, dummy3,    &
                         modc=grid%modc_dynrun )
         do lat=beglat,endlat
            do k=beglev,endlev
               do i=1,plon
                  qfcst(i,k,l,lat)=dummy3(i,k,lat)
               end do
            end do
         end do
      end do
!!
#endif
   else  !jjr nned define _yz
      sigmp=sigmpxy
      qfcst=qfcstxy
   endif
   call t_stopf ('sltb2')
!
! Copy SLT advection structure data back into dynamics data
!
   call t_startf ('scanslt_ad_coup')
   call ad_coupling( adv_state )
   call t_stopf ('scanslt_ad_coup')

   return
end subroutine scanslt_run

!
!-----------------------------------------------------------------------
!
subroutine scanslt_final( adv_state )
!----------------------------------------------------------------------- 
! 
! Purpose: 
! SLT finalization for Eulerian dynamics
! 
! Author: 
!
! Erik Kluzek
!
!----------------------------------------------------------------------- 
!
! Arguments
!
   type(advection_state), intent(inout) :: adv_state    ! Advection state data

   call adv_state_dealloc( adv_state )
end subroutine scanslt_final

!
!-----------------------------------------------------------------------
!

subroutine ad_coupling( adv_state )
!----------------------------------------------------------------------- 
! 
! Purpose: 
! Copy advection data into dynamics state.
! 
! Author: 
!
! Erik Kluzek
!
!----------------------------------------------------------------------- 
   use prognostics, only: u3, v3, qminus, q3, n3m1, n3
! ====================================== zhh ==================================
   use IAP_prog,    only: Q, Qliq, Qice
   use IAP_grid,    only: EX
! ==================================== 2008.1.8 ===============================
!
! Arguments
!
   type(advection_state), intent(in) :: adv_state    ! Advection state data

   integer :: i, j, k, c, jd   ! Indices

#ifdef OUTER_OMP
!wjp 2011.04.16 !$OMP PARALLEL DO PRIVATE (I,K,J,C)
#endif
   do j = beglat, endlat
      jd = plat+1-j   ! zhang he
!$OMP PARALLEL DO PRIVATE (K,I,C)
      do k=beglev,endlev
         do i = 1, plon
!!            u3(i,k,j,n3m1) = adv_state%u3(i+i1-1,k,j+beglatex+numbnd-beglat)
!!            v3(i,k,j,n3m1) = adv_state%v3(i+i1-1,k,j+beglatex+numbnd-beglat)
!!            Q(i+EX,k,jd)    = qfcst(i,k,1,j)
!!            Qliq(i+EX,k,jd) = qfcst(i,k,2,j)
!!            Qice(i+EX,k,jd) = qfcst(i,k,3,j)
!!            do c = 4, pcnst   !zhh 2012-01-14
!!            do c = 1, pcnst       !zhh 2013-03-12
!!               qminus(i,k,c,j) = adv_state%qminus(i+i1-1,k,c,j+beglatex+numbnd-beglat)
!!            end do
!!            do c = 1, 3      !zhh 2013-04-15
            do c = 1, pcnst      !zhh 2013-04-15
               q3(i,k,c,j,n3) = qfcst(i,k,c,j)
            end do
         end do
      end do
   end do

end subroutine ad_coupling

!
!-----------------------------------------------------------------------
!

subroutine da_coupling( cwava, etadot, adv_state )
!----------------------------------------------------------------------- 
! 
! Purpose: 
! Copy dynamics data into advection state
! Also find the total moisture mass before SLT.
! 
! Author: 
!
! Erik Kluzek
!
!----------------------------------------------------------------------- 
   use prognostics, only: u3, v3, qminus, n3m1, ps, q3, pdeld, n3
   use commap,      only: w
   use rgrid,       only: nlon
   use qmassa,      only: qmassarun
! ====================================== zhh ==================================
   use IAP_grid,  only : NX, NY, NL, NZ, JB, JE, IM, EX
   use IAP_prog,   only: U, V, WS, UVW0, Q, Qliq, Qice
   use Dyn_const,  only: DXVPN, DXVPS   
   use mathconst,  only: ZERO, HALF
! ==================================== 2008.1.7 ===============================
	
!
! Arguments
!
   real(r8), intent(in) :: cwava(plat)             ! weight for global water vapor int.
   type(advection_state), intent(out) :: adv_state ! Advection state data
   real(r8), intent(out) :: etadot(plon,beglev:endlev+1,beglat:endlat) ! vertical motion (slt)
!
! Local variables
!
   integer :: i, j, k, c, irow, lat ! Indices

   real(r8) :: pmid(plon,beglev:endlev)      ! pressure at model levels
   real(r8) :: pint(plon,beglev:endlev+1)     ! pressure at interfaces
   real(r8) :: pdel(plon,beglev:endlev)      ! pressure difference between
! ====================================== zhh ==================================
   real(r8) :: UQ(NX,beglev:endlev,beglatdyn:endlatdyn)
   real(r8) :: VQ(NX,beglev:endlev,beglatdynex:endlatdyn)
   real(r8) :: WQ(NX,beglev:endlevp,beglatdyn:endlatdyn)  
   real(r8) :: u3mean, v3mean
   integer  :: jd, jp,src,dest                  
! ==================================== 2008.1.8 ===============================
   type (T_FVDYCORE_STATE), pointer :: dyn_state   !zhh 2013-01-25
   integer  :: commyz, commxy     !zhh 2013-01-25
!-----------------------------------------------------------------------
!
   dyn_state => get_dyn_state()
   commyz = dyn_state%grid%commyz
!
! Initialize moisture mass integrals.
!
   hw1lat = 0.0_r8
!
! Find moisture mass before SLT
!
#ifdef OUTER_OMP
!wjp 2011.04.16 !$OMP PARALLEL DO PRIVATE (LAT, IROW, PINT, PMID, PDEL)
#endif
   do lat=beglat,endlat
      if(lat.le.plat/2) then
         irow = lat
      else
         irow = plat + 1 - lat
      end if
!
! Only pdel is needed inside SLT.  pint and pmid are not.
!
      call plevs00 (nlon(lat),plon,plev,ps(1,lat,n3m1), pint, pmid, pdel)
!
! Calculate mass of moisture in field being advected by slt. (hw1lat)
!
!  q3     is plon,plev,pcnst,beglat:endlat,ptimelevs
!  qminus is plon,plev,pcnst,beglat:endlat
!
      call qmassarun (cwava(lat),w(irow) ,qminus(1,beglev,1,lat),pdel    , &
                   hw1lat(1,lat),nlon(lat), q3(1,beglev,1,lat,n3m1), lat, pdeld(1,beglev,lat,n3m1 ))
   end do

#ifdef OUTER_OMP
!wjp 2011.04.16 !$OMP PARALLEL DO PRIVATE (I,K,J,C)
#endif
! ====================================== zhh ==================================
!!   do j = beglat, endlat
   do j = beglatdyn, endlatdyn
      do k = beglev, endlev
         do i = 1, NX
            UQ(I,K,J)     = HALF*(U (I,K,J)+UVW0(I,K,J,1))
            VQ(I,K,J)     = HALF*(V (I,K,J)+UVW0(I,K,J,2))
            WQ(I,K,J)     = HALF*(WS(I,K,J)+UVW0(I,K,J,3))
         end do
      end do
   end do
   do j = beglatdyn, endlatdyn
      do i = 1, NX
         if(beglev.eq.1)     WQ(I,1,J)  = ZERO      !zhh 2008.4.30
         if(endlev.eq.plev)  WQ(I,NZ,J) = ZERO !!!jjr need change
      end do
   end do
!
#if (defined SPMD)
!JJR 
   src = iam+1
   dest  = iam-1
   if ( mod(iam,npr_y) == 0 ) dest = -1
   if ( mod(iam+1,npr_y) == 0 ) src = -1
   call mp_send3d( commyz, dest, src, NX, NL, NY,                      &
                   1, NX,beglev,endlev, beglatdynex, endlatdyn,        &
                   1, NX,beglev,endlev, endlatdyn, endlatdyn, vq )
   call mp_recv3d( commyz, src, NX, NL, NY,                            &
                   1, NX, beglev,endlev,beglatdynex, endlatdyn,        &
                   1, NX, beglev,endlev,beglatdynex, beglatdynex, vq )

#endif
!
   do jd = beglatdyn, endlatdyn
      jp = plat+1-jd
! for etadot
      do k = beglev,endlevp
         do i = 1, plon
            etadot(i,k,jp) = WQ(i+EX,k,jd)
         end do
      end do
! for u3 & v3
      if ( jd >= JB .and. jd <= JE ) then
         do k = beglev, endlev
            do i = 1, plon
               u3mean = HALF * ( UQ(i+EX,k,jd) + UQ(i+EX+1,k,jd) )
               v3mean = - (DXVPN(jd)*VQ(i+EX,k,jd-1) + DXVPS(jd)*VQ(i+EX,k,jd))
               adv_state%u3(i+i1-1,k,jp+beglatex+numbnd-beglat) = u3mean
               adv_state%v3(i+i1-1,k,jp+beglatex+numbnd-beglat) = v3mean
!!            adv_state%u3(i+i1-1,k,j+beglatex+numbnd-beglat) = u3(i,k,j,n3m1)
!!            adv_state%v3(i+i1-1,k,j+beglatex+numbnd-beglat) = v3(i,k,j,n3m1)
! ===============================================================================
!!               do c = 4, pcnst       !zhh 2012-01-14
!!               do c = 1, pcnst       !zhh 2013-03-12
!!                  adv_state%qminus(i+i1-1,k,c,jp+beglatex+numbnd-beglat) = qminus(i,k,c,jp)
!!               end do
!!               adv_state%qminus(i+i1-1,k,1,jp+beglatex+numbnd-beglat) = Q(i+EX,k,jd)
!!               adv_state%qminus(i+i1-1,k,2,jp+beglatex+numbnd-beglat) = Qliq(i+EX,k,jd)
!!               adv_state%qminus(i+i1-1,k,3,jp+beglatex+numbnd-beglat) = Qice(i+EX,k,jd)
! ===============================================================================
            end do
         end do
      else   ! at north and south pole
         do k = beglev, endlev
            do i = 1, plon
               adv_state%u3(i+i1-1,k,jp+beglatex+numbnd-beglat) = ZERO
               adv_state%v3(i+i1-1,k,jp+beglatex+numbnd-beglat) = ZERO
! ===============================================================================
!!               do c = 4, pcnst     !zhh 2012-01-14
!!               do c = 1, pcnst       !zhh 2013-03-12
!!                  adv_state%qminus(i+i1-1,k,c,jp+beglatex+numbnd-beglat) = qminus(i,k,c,jp)
!!               end do
!!               adv_state%qminus(i+i1-1,k,1,jp+beglatex+numbnd-beglat) = Q(i+EX,k,jd)
!!               adv_state%qminus(i+i1-1,k,2,jp+beglatex+numbnd-beglat) = Qliq(i+EX,k,jd)
!!               adv_state%qminus(i+i1-1,k,3,jp+beglatex+numbnd-beglat) = Qice(i+EX,k,jd)
! ===============================================================================
            end do
         end do
      end if
!
      do k = beglev,endlevp
         do i = 1, plon
            do c = 1, pcnst     !zhh 2012-01-14
!!            do c = 1, 3     !zhh 2012-01-14
               adv_state%qminus(i+i1-1,k,c,jp+beglatex+numbnd-beglat) = q3(i,k,c,jp,n3)
            end do
         end do
      end do
   end do
! ==================================== 2008.6.5 ===============================

end subroutine da_coupling

!
!-----------------------------------------------------------------------
!

subroutine adv_state_alloc( adv_state )
!----------------------------------------------------------------------- 
! 
! Purpose: 
! Allocate advection state data
! 
! Author: 
!
! Erik Kluzek
!
!----------------------------------------------------------------------- 
   use infnan,       only: inf
!
! Arguments
!
   type(advection_state), intent(out) :: adv_state    ! Advection state data

   allocate (adv_state%u3    (plond,beglev:endlev            ,beglatex:endlatex) )
   allocate (adv_state%v3    (plond,beglev:endlev            ,beglatex:endlatex) )
   allocate (adv_state%qminus(plond,beglev:endlev,pcnst      ,beglatex:endlatex) )
   adv_state%u3    (:,:,  beglatex:endlatex) = inf
   adv_state%v3    (:,:,  beglatex:endlatex) = inf
   adv_state%qminus(:,:,:,beglatex:endlatex) = inf

end subroutine adv_state_alloc

!
!-----------------------------------------------------------------------
!

subroutine adv_state_dealloc( adv_state )
!----------------------------------------------------------------------- 
! 
! Purpose: 
! De-allocate advection state data
! 
! Author: 
!
! Erik Kluzek
!
!----------------------------------------------------------------------- 
!
! Arguments
!
   type(advection_state), intent(inout) :: adv_state    ! Advection state data

   deallocate (adv_state%u3    )
   deallocate (adv_state%v3    )
   deallocate (adv_state%qminus)

end subroutine adv_state_dealloc

!
!-----------------------------------------------------------------------
!

subroutine grdini(pmap    ,etamid  ,etaint  ,gravit  ,dlam    , &
                  lam     ,phi     ,dphi    ,gw      ,sinlam  , &
                  coslam  ,lbasdy  ,lbasdz  ,lbassd  ,lbasiy  , &
                  detam   ,detai   ,kdpmpf  ,kdpmph  ,cwava   )

!----------------------------------------------------------------------- 
! 
! Purpose: 
! Initialize model and extended grid parameters
! Initialize weights for Lagrange cubic derivative estimates
! Initialize weights for Lagrange cubic interpolant
! 
! Method: 
! 
! Author: 
! Original version:  J. Olson
! Standardized:      J. Rosinski, June 1992
! Reviewed:          D. Williamson, P. Rasch, August 1992
! Reviewed:          D. Williamson, P. Rasch, March 1996
!
!-----------------------------------------------------------------------
   use rgrid,   only: nlon
!------------------------------Parameters-------------------------------
!
! Input arguments
!
   integer, intent(in) :: pmap              ! dimension of artificial vert. grid
!
   real(r8), intent(in) :: etamid(plev)         ! full-level model vertical grid
   real(r8), intent(in) :: etaint(plevp)        ! half-level model vertical grid
   real(r8), intent(in) :: gravit               ! gravitational constant
!
! Output arguments
!
   real(r8), intent(out) :: dlam(platd)          ! longitudinal grid interval (radians)
   real(r8), intent(out) :: lam   (plond,platd)  ! longitudinal coords of extended grid
   real(r8), intent(out) :: phi   (platd)        ! latitudinal  coords of extended grid
   real(r8), intent(out) :: dphi  (platd)        ! latitude intervals (radians)
   real(r8), intent(out) :: gw    (plat)         ! Gaussian weights
   real(r8), intent(out) :: sinlam(plond,platd)  ! sin(lam) model domain only
   real(r8), intent(out) :: coslam(plond,platd)  ! cos(lam) model domain only
   real(r8), intent(out) :: lbasdy(4,2,platd)    ! latitude derivative weights
   real(r8), intent(out) :: lbasdz(4,2,plev)     ! vertical (full levels) deriv weights
   real(r8), intent(out) :: lbassd(4,2,plevp)    ! vertical (half levels) deriv weights
   real(r8), intent(out) :: lbasiy(4,2,platd)    ! Lagrange cubic interp weights (lat.)
   real(r8), intent(out) :: detam (plev)         ! intervals between vertical full levs.
   real(r8), intent(out) :: detai (plevp)        ! intervals between vertical half levs.
!
   integer, intent(out) :: kdpmpf(pmap)      ! artificial full vertical grid indices
   integer, intent(out) :: kdpmph(pmap)      ! artificial half vertical grid indices
!
   real(r8), intent(out) :: cwava(plat)          ! weight applied to global integrals
!
!-----------------------------------------------------------------------
!
!  pmap    Dimension of artificial evenly spaced vertical grid arrays
!  etamid  Full-index hybrid-levels in vertical grid.
!  etaint  Half-index hybrid-levels from sig(1/2) = etaint(1) = 0. to
!          sig(plev+1/2) = etaint(plevp) = 1.
!  gravit  Gravitational constant.
!  dlam    Length of increment in longitude grid.
!  lam     Longitude values in the extended grid.
!  phi     Latitude values in the extended grid.
!  dphi    Interval between latitudes in the extended grid
!  gw      Gauss weights for latitudes in the global grid.  (These sum
!          to 2.0.)
!  sinlam  Sine of longitudes in global grid (no extension points).
!  coslam  Cosine of longitudes in global grid (no extension points).
!  lbasdy  Weights for Lagrange cubic derivative estimates on the
!          unequally spaced latitude grid
!  lbasdz  Weights for Lagrange cubic derivative estimates on the
!          unequally spaced vertical grid (corresponding to model
!          full levels).
!  lbassd  Weights for Lagrange cubic derivative estimates on the
!          unequally spaced vertical grid (corresponding to model
!          half levels).
!  lbasiy  Weights for Lagrange cubic interpolation on the
!          unequally spaced latitude grid
!  detam   Increment between model mid-levels ("full" levels)
!  detai   Increment between model interfaces ("half" levels).
!  kdpmpf  Array of indicies of the model full levels which are mapped
!          into an artificial evenly spaced vertical grid.  Used to aid
!          in search for vertical position of departure point 
!  kdpmph  Array of indicies of the model half levels which are mapped
!          into an artificial evenly spaced vertical grid.  Used to aid
!          in search for vertical position of departure point 
!  cwava   1./(plon*gravit)
!
!---------------------------Local variables-----------------------------
!
   integer j                 ! index
   integer k                 ! index
!
   real(r8) etamln(plev)         ! log(etamid)
   real(r8) etailn(plevp)        ! log(etaint)
   real(r8) detamln(plev)        ! dlog(etamid)
   real(r8) detailn(plevp)       ! dlog(etaint)
!
!-----------------------------------------------------------------------
   if (single_column) then

      dlam(:)=0._r8
      lam(:,:)=0._r8
      phi(:)=0._r8
      dphi(:)=0._r8
      sinlam(:,:)=0._r8
      coslam(:,:)=0._r8
      detai(:)=0._r8
      kdpmpf(:)=0._r8
      kdpmph(:)=0._r8
      gw(:)=1._r8
      call basdz(plev    ,etamid  ,lbasdz  )
      call basdz(plevp   ,etaint  ,lbassd  )

   else
!
! Initialize extended horizontal grid coordinates.
!
      call grdxy(dlam    ,lam     ,phi     ,gw      ,sinlam  , &
                 coslam  )
!
! Basis functions for computing Lagrangian cubic derivatives
! on unequally spaced latitude and vertical grids.
!
      call basdy(phi     ,lbasdy  )

      call basdz(plev    ,etamid  ,lbasdz  )
      call basdz(plevp   ,etaint  ,lbassd  )


!
! Basis functions for computing weights for Lagrangian cubic
! interpolation on unequally spaced latitude grids.
!
      call basiy(phi     ,lbasiy  )
!
! Compute interval lengths in latitudinal grid
!
      do j = 1,platd-1
         dphi(j) = phi(j+1) - phi(j)
      end do

   endif
!
! Compute interval lengths in vertical grids.
!
   do k = 1,plev
      etamln(k) = log(etamid(k))
   end do
   do k = 1,plevp
      etailn(k) = log(etaint(k))
   end do
   do k = 1,plev-1
      detam  (k) = etamid(k+1) - etamid(k)
      detamln(k) = etamln(k+1) - etamln(k)
   end do
   do k = 1,plev
      detai  (k) = etaint(k+1) - etaint(k)
      detailn(k) = etailn(k+1) - etailn(k)
   end do
!
! Build artificial evenly spaced vertical grid for use in determining
! vertical position of departure point.
! Build one grid for full model levels and one for half levels.
!
   call vrtmap(plev    ,pmap    ,etamln  ,detamln ,kdpmpf  )
   call vrtmap(plevp   ,pmap    ,etailn  ,detailn ,kdpmph  )
!
! Compute moisture integration constant
!
if (single_column) then
   cwava = 1._r8
else
   do j=1,plat
      cwava(j) = 1._r8/(nlon(j)*gravit)
   end do
endif
!
   return
end subroutine grdini

!
!-----------------------------------------------------------------------
!

subroutine grdxy(dlam    ,lam     ,phi     ,w       ,sinlam  , &
                 coslam  )

!----------------------------------------------------------------------- 
! 
! Purpose: 
! Define the "extended" grid used in the semi-Lagrangian transport
! scheme.  The longitudes are equally spaced and the latitudes are
! Gaussian.  The global grid is extended to include "wraparound" points
! on all sides.
! 
! Method: 
! 
! Author: J. Olson
! 
! Modify: ZhangHe, 2008.4.18
!-----------------------------------------------------------------------
  use Dyn_const, only: wlat      !zhh 2008.4.15
!------------------------------Parameters-------------------------------
  integer, parameter :: istart = nxpt+1         ! index for first model long.
  integer, parameter :: jstart = nxptj+jintmx+1  ! index for first model lat.
  integer, parameter :: jstop  = jstart-1+plat  ! index for last  model lat.
!-----------------------------------------------------------------------

!------------------------------Arguments--------------------------------
  real(r8), intent(out) :: dlam(platd)          ! longitudinal increment
  real(r8), intent(out) :: lam   (plond,platd)  ! long. coords. in extended grid
  real(r8), intent(out) :: phi   (platd)        ! lat.  coords. in extended grid
  real(r8), intent(out) :: w     (plat)         ! Gaussian weights
  real(r8), intent(out) :: sinlam(plond,platd)  ! sin(lam)
  real(r8), intent(out) :: coslam(plond,platd)  ! cos(lam)
!
! dlam    Length of increment in longitude grid.
! lam     Longitude values in the extended grid.
! phi     Latitude values in the extended grid.
! w       Gauss weights for latitudes in the global grid.  (These sum
!         to 2.0 like the ones in CCM1.)
! sinlam  Sine of longitudes in global grid (no extension points).
! coslam  Cosine of longitudes in global grid (no extension points).
!-----------------------------------------------------------------------

!---------------------------Local variables-----------------------------
  integer i,j,ig            ! indices
  integer nlond             ! extended long dim
  real(r8) lam0             ! lamda = 0
  real(r8) pi               ! 3.14...
  real(r8) wrk(platd)       ! work space
!-----------------------------------------------------------------------
!
  lam0 = 0.0_r8
  pi = 4._r8*atan(1._r8)
!
! Interval length in equally spaced longitude grid.
!
  do j=1,platd
     dlam(j) = 2._r8*pi/real(nlonex(j),r8)
!
! Longitude values on extended grid.
!
     nlond = nlonex(j) + 1 + 2*nxpt
     do i = 1,nlond
        lam(i,j) = real(i-istart,r8)*dlam(j) + lam0
     end do
  end do
!
  do j=1,plat
     w(j) = wlat(j)  !zhang he
  end do
! Reorder and compute latitude values.
!
  do j = jstart,jstop
!!     wrk(j) = asin( phi(jstop-j+1) ) !zhang he
     phi(j) = -pi/2.0 + (j-jstart)*pi/dble(plat-1)  !zhang he
  end do
!  phi(jstart:jstop) = wrk(jstart:jstop)  !zhang he
!
! North and south poles.
!
!!  phi(jstart-1) = -pi/2.0  !zhang he
!!  phi(jstop +1) =  pi/2.0  !zhang he
  phi(jstart) = -pi/2.0
  phi(jstop ) =  pi/2.0
!
! Extend Gauss latitudes below south pole so that the spacing above
! the pole is symmetric, and phi is decreasing, i.e., phi < -pi/2
!
!!  if( jstart > 2 )then   !zhanghe
!!     do j = 1,jstart-2
!!        phi(j) = -pi - phi(2*jstart-2-j)
!!     end do
!!  end if
  if( jstart > 1 )then
     do j = 1,jstart-1
        phi(j) = -pi - phi(2*jstart-j)
     end do
  end if
!
! Analogously for Northern Hemisphere
!
!!  if( platd > jstop+1 )then  !zhanghe
!!     do j = jstop+2,platd
!!        phi(j) = pi - phi(2*jstop+2-j)
!!     end do
!!  end if
  if( platd > jstop )then
     do j = jstop+1,platd
        phi(j) = pi - phi(2*jstop-j)
     end do
  end if
!
! Sine and cosine of longitude.
!
  do j=1,platd
     ig = 0
     do i = istart,nlonex(j)+nxpt
        ig = ig + 1
        sinlam(ig,j) = sin( lam(i,j) )
        coslam(ig,j) = cos( lam(i,j) )
     end do
  end do

  return
end subroutine grdxy

!
!-----------------------------------------------------------------------
!

subroutine sltb1(jcen    ,jgc     ,dt      ,ra      , &
                 iterdp    ,uxl     ,uxr     ,vxl     ,vxr     , &
                 fxl     ,fxr     ,lam     ,phib    , &
                 dphib       , &
                 lbasdy    ,lbasiy    , &
                 lammp   ,phimp    , &
                 adv_state ,nlon    ,hadv    ,nlonex,fhr  )

!----------------------------------------------------------------------- 
! 
! Purpose: 
! Drive the slt algorithm on a given latitude slice in the extended
! data arrays using information from the entire latitudinal extent
! of the arrays.
! 
! Method: 
! Compute departure points and corresponding indices.
! Poleward of latitude phigs (radians), perform the computation in
! local geodesic coordinates.
! Equatorward of latitude phigs, perform the computation in global
! spherical coordinates
! 
! Author: J. Olson
! Modified: Jiang Jinrong, for 2D parallel
! 
!-----------------------------------------------------------------------
  use prognostics, only: ptimelevels

#include <parslt.h>

!------------------------------Parameters-------------------------------
  real(r8), parameter :: phigs = 1.221730_r8 ! cut-off latitude: about 70 degrees
!-----------------------------------------------------------------------

!------------------------------Arguments--------------------------------
  integer , intent(in) :: nlon                ! longitude dimension
  integer , intent(in) :: nlonex(platd)       ! extended longitude dimension
  integer , intent(in) :: jcen                ! index of lat slice(extend)
  integer , intent(in) :: jgc                 ! index of lat slice (model)
  real(r8), intent(in) :: dt                  ! time step (seconds)
  real(r8), intent(in) :: ra                  ! 1./(radius of earth)
  integer , intent(in) :: iterdp              ! iteration count
  real(r8), intent(in) :: uxl(plond,beglev:endlev,beglatex:endlatex) ! left  x-deriv of ub
  real(r8), intent(in) :: uxr(plond,beglev:endlev,beglatex:endlatex) ! right x-deriv of ub
  real(r8), intent(in) :: vxl(plond,beglev:endlev,beglatex:endlatex) ! left  x-deriv of vb
  real(r8), intent(in) :: vxr(plond,beglev:endlev,beglatex:endlatex) ! right x-deriv of vb
  real(r8), intent(in) :: fxl(plond,beglev:endlev,  pcnst,beglatex:endlatex) ! left  fb x-deriv
  real(r8), intent(in) :: fxr(plond,beglev:endlev,  pcnst,beglatex:endlatex) ! right fb x-deriv
  real(r8), intent(in) :: lam  (plond,platd)  ! long. coord of model grid
  real(r8), intent(in) :: phib (platd)        ! lat.  coord of model grid
  real(r8), intent(in) :: dphib(platd)        ! increment between lats.
  real(r8), intent(in) :: lbasdy(4,2,platd)   ! lat deriv weights
  real(r8), intent(in) :: lbasiy(4,2,platd)   ! lat interp wts(lagrng)
  real(r8), intent(inout) ::  hadv  (plon, beglev:endlev, pcnst, beglat:endlat) ! horizontal advection tendency
  real(r8), intent(inout) ::  lammp(plon,beglev:endlev)       ! long coord of mid-point
  real(r8), intent(inout) ::  phimp(plon,beglev:endlev)       ! lat  coord of mid-point
  type(advection_state), intent(in) :: adv_state    ! Advection state
!
!  pmap    Dimension of kdpmpX arrays
!  jcen    Latitude index in extended grid corresponding to lat slice
!          being forecasted.
!  jgc     Latitude index in model    grid corresponding to lat slice 
!          being forecasted.
!  dt      Time interval that parameterizes the parcel trajectory.
!  ra      Reciprocal of radius of earth.
!  iterdp  Number of iterations used for departure point calculation.
!  uxl     x-derivatives of u at the left  (west) edge of given interval
!  vxl     x-derivatives of v at the left  (west) edge of given interval
!  uxr     x-derivatives of u at the right (east) edge of given interval
!  vxr     x-derivatives of v at the right (east) edge of given interval
!  wb      z-velocity component (eta-dot).
!  fxl     x-derivatives at the left  edge of each interval containing 
!          the departure point.
!  fxr     x-derivatives at the right edge of each interval containing 
!          the departure point.
!  lam     Longitude values for the extended grid.
!  phib    Latitude  values for the extended grid.
!  dphib   Interval between latitudes in the extended grid.
!  sig     Hybrid eta values at the "full-index" levels.
!  sigh    Half-index eta-levels including sigh(i,1) = eta(1/2) = 0.0
!          and sigh(i,plev+1) = eta(plev+1/2) = 1.  Note that in general
!          sigh(i,k) .lt. sig(i,k)  where sig(i,k) is the hybrid value
!          at the k_th full-index level.
!  dsig    Interval lengths in full-index hybrid level grid.
!  dsigh   Interval lengths in half-index hybrid level grid.
!  lbasdy  Weights for Lagrange cubic derivative estimates on the
!          unequally spaced latitude grid.
!  lbasdz  Weights for Lagrange cubic derivative estimates on the
!          unequally spaced vertical grid (full levels).
!  lbassd  Weights for Lagrange cubic derivative estimates on the
!          unequally spaced vertical grid (half levels).
!  lbasiy  Weights for Lagrange cubic interpolation on the unequally
!          spaced latitude grid.
!  kdpmpf  indices of artificial grid mapped into the full level grid
!  kdpmph  indices of artificial grid mapped into the half level grid
!  lammp   Longitude coordinates of the trajectory mid-points of the
!          parcels that correspond to the global grid points contained
!          in the latitude slice being forecasted.  On entry lammp
!          is an initial guess.
!  phimp   Latitude coordinates of the trajectory mid-points of the
!          parcels that correspond to the global grid points contained
!          in the latitude slice being forecasted.  On entry phimp
!          is an initial guess.
!  sigmp   Hybrid value at the trajectory midpoint for each gridpoint
!          in a vertical slice from the global grid.  On entry sigmp is
!          an initial guess.
!  fbout   Extended array only one latitude of which, however, is filled
!          with forecasted (transported) values.  This routine must be
!          called multiple times to fill the entire array.  This is
!          done to facilitate multi-tasking.
!-----------------------------------------------------------------------

!---------------------------Local variables-----------------------------
  integer m                            ! constituent index
  integer idp(plon,beglev:endlev,4)             ! zonal      dep point index
  integer jdp(plon,beglev:endlev)               ! meridional dep point index
  real(r8),intent(out):: fhr(plon,beglev:endlev,pcnst)        ! horizontal interpolants
  real(r8) lamdp(plon,beglev:endlev)            ! zonal      departure pt. coord.
  real(r8) phidp(plon,beglev:endlev)            ! meridional departure pt. coord.
  real(r8) fint(plon,beglev:endlev,ppdy,pcnst)  ! work space
  real(r8) fyb(plon,beglev:endlev,pcnst)        ! work space
  real(r8) fyt(plon,beglev:endlev,pcnst)        ! work space
  logical locgeo                       ! flag indicating coordinate sys
  integer :: k,i                       ! indices (needed for SCAM)
!-----------------------------------------------------------------------
  if (.not. single_column) then

!
! Horizontal interpolation
!
     locgeo = abs(phib(jcen))>=phigs
!
     call sphdep(jcen    ,jgc     ,dt      ,ra      ,iterdp  ,                     &
              locgeo  ,adv_state%u3(1,beglev,beglatex)       ,uxl     ,uxr     ,lam     ,   &
              phib    ,lbasiy  ,lammp   ,phimp   ,lamdp   ,                     &
              phidp   ,idp     ,jdp     ,adv_state%v3(1,beglev,beglatex),                   &
              vxl     ,vxr     ,nlon    ,nlonex  )
!
! Interpolate scalar fields to the departure points.
!
     call hrintp(pcnst   ,pcnst   ,adv_state%qminus(1,beglev,1,beglatex), fxl     ,fxr     , &
              lam     ,phib    ,dphib   ,lbasdy  ,lamdp   ,                    &
              phidp   ,idp     ,jdp     ,jcen    ,plimdr  ,                    &
              fint    ,fyb     ,fyt     ,fhr     ,nlon    ,                    &   
              nlonex  )

     do m = 1,pcnst
!$OMP PARALLEL DO PRIVATE (K, I)
     do k = beglev,endlev
           do i = 1,nlon
              hadv(i,k,m,jgc) = (fhr(i,k,m) - adv_state%qminus(i1-1+i,k,m,jcen))/dt
           end do
        end do
     end do
else
!
! fill in fhr in leiu of horizontal interpolation
!
   do m = 1,pcnst
     do k = beglev,endlev
         do i = 1,nlon
            fhr(i,k,m) = adv_state%qminus(i1+i-1,k,m,jcen)
         end do
      end do
   end do
endif
!
  return
end subroutine sltb1
!
!============================================================================================
subroutine sltb2(lat,pmap             ,dt           , &
                 iterdp         , &
                 wb              , &
                 sig     ,sigh    ,dsig    ,dsigh   , &
                 lbasdz  ,lbassd  ,kdpmpf  , &
                 kdpmph      ,sigmp   ,fbout   , &
                 nlon,fhr       )

!-----------------------------------------------------------------------
!
! Purpose:
! Drive the slt algorithm on a given latitude slice in the extended
! data arrays using information from the entire latitudinal extent
! of the arrays.
!
! Method:
! Compute departure points and corresponding indices.
! Poleward of latitude phigs (radians), perform the computation in
! local geodesic coordinates.
! Equatorward of latitude phigs, perform the computation in global
! spherical coordinates
!
! Author: J. Olson
!
!-----------------------------------------------------------------------
  use prognostics, only: ptimelevels
  implicit none

#include <parslt.h>

!------------------------------Parameters-------------------------------
  real(r8), parameter :: phigs = 1.221730 ! cut-off latitude: about 70 degrees
!-----------------------------------------------------------------------

!------------------------------Arguments--------------------------------
  integer ,intent(in) :: lat
  integer , intent(in) :: nlon                ! longitude dimension
  integer , intent(in) :: pmap                ! artificial vert grid dim.
  real(r8), intent(in) :: dt                  ! time step (seconds)
  integer , intent(in) :: iterdp              ! iteration count
  real(r8), intent(in) :: wb(beglonxy:endlonxy,plevp)                    ! eta-dot
  real(r8), intent(in) :: sig  (plev)         ! vertical full levels
  real(r8), intent(in) :: sigh (plevp)        ! vertical half levels
  real(r8), intent(in) :: dsig (plev)         ! inc. between full levs
  real(r8), intent(in) :: dsigh(plevp)        ! inc. between half levs
  real(r8), intent(in) :: lbasdz(4,2,plev)    ! vert full level deriv wts
  real(r8), intent(in) :: lbassd(4,2,plevp)   ! vert half level deriv wts
  integer , intent(in) :: kdpmpf(pmap)        ! artificial vert grid index
  integer , intent(in) :: kdpmph(pmap)        ! artificial vert grid index
  real(r8), intent(in) :: fhr(beglonxy:endlonxy,plev,pcnst)        ! horizontal interpolants
  real(r8), intent(inout) ::  sigmp(beglonxy:endlonxy,plev)       ! vert coord of mid-point
  real(r8), intent(out) :: fbout(beglonxy:endlonxy,plev,pcnst)   ! advected constituents
!---------------------------Local variables-----------------------------
  integer m                            ! constituent index
  integer kdp(beglonxy:endlonxy,plev)               ! vertical   dep point index
  real(r8) sigdp(beglonxy:endlonxy,plev)            ! vertical   departure pt. coord.
  real(r8) fhst(beglonxy:endlonxy,plev,pcnst)       ! derivative at top of interval
  real(r8) fhsb(beglonxy:endlonxy,plev,pcnst)       ! derivative at bot of interval
  real(r8) wst(beglonxy:endlonxy,plevp)             ! w derivative at top of interval
  real(r8) wsb(beglonxy:endlonxy,plevp)             ! w derivative at bot of interval
  integer :: k,i                       ! indices (needed for SCAM)
!-----------------------------------------------------------------------
!
! Vertical interpolation. end jjr 2012.5.24
!!JJR need change MPP deomposition zy to xy
! Compute vertical derivatives of vertical wind

!
  call cubzdr(nlon    ,plevp   ,wb      ,lbassd  ,wst     , &
              wsb     )
!

! Compute departure points and corresponding indices.
!
  call vrtdep(pmap    ,dt      ,iterdp  ,wb      ,wst     , &
              wsb     ,sig     ,sigh    ,dsigh   ,kdpmpf  , &
              kdpmph  ,sigmp   ,sigdp   ,kdp     ,nlon    )
!

! Vertical derivatives of scalar fields.
! Loop over constituents.
!
!CSD$ PARALLEL DO 
  do m = 1,pcnst
     call cubzdr(nlon    ,plev    ,fhr(beglonxy,1,m), lbasdz  ,fhst(beglonxy,1,m), &
                 fhsb(beglonxy,1,m) )
  end do

!CSD$ END PARALLEL DO
  if( plimdr )then
     call limdz(fhr     ,dsig    ,fhst    ,fhsb    ,nlon    )
  end if

!
! Vertical interpolation of scalar fields.
!
  call herzin(plev    ,pcnst   ,fhr     ,fhst    ,fhsb    , &
              sig     ,dsig    ,sigdp   ,kdp     ,fbout   , &
              nlon    )


  return
end subroutine sltb2

!
!============================================================================================
!

subroutine vrtdep(pmap    ,dt      ,iterdp  ,wb      ,wst     , &
                  wsb     ,sig     ,sigh    ,dsigh   ,kdpmpf  , & 
                  kdpmph  ,sigmp   ,sigdp   ,kdp     ,nlon    )

!----------------------------------------------------------------------- 
! 
! Purpose: 
! Compute vertical departure point and departure point index.
! 
! Method: 
! 
! Author: J. Olson
! 
!-----------------------------------------------------------------------
!------------------------------Arguments--------------------------------
  integer , intent(in) :: nlon                ! longitude dimension
  integer , intent(in) :: pmap                ! dimension of artificial vert grid
  real(r8), intent(in) :: dt                  ! time step (seconds)
  integer , intent(in) :: iterdp              ! number of iterations
  real(r8), intent(in) :: wb (beglonxy:endlonxy,plevp)     ! vertical velocity
  real(r8), intent(in) :: wst(beglonxy:endlonxy,plevp)     ! z-derivative of wb at top of interval
  real(r8), intent(in) :: wsb(beglonxy:endlonxy,plevp)     ! z-derivative of wb at bot of interval
  real(r8), intent(in) :: sig  (plev )        ! sigma values of model full levels
  real(r8), intent(in) :: sigh (plevp)        ! sigma values of model half levels
  real(r8), intent(in) :: dsigh(plevp)        ! increment between half levels
  integer , intent(in) :: kdpmpf(pmap)        ! artificial grid indices
  integer , intent(in) :: kdpmph(pmap)        ! artificial grid indices
  real(r8), intent(inout) :: sigmp(beglonxy:endlonxy,plev) ! vert coords of traj mid-points
  real(r8), intent(out) :: sigdp(beglonxy:endlonxy,plev)   ! vert coords of traj departure points
  integer , intent(out) :: kdp(beglonxy:endlonxy,plev)     ! vertical departure point indices
!
!  pmap    Dimension of kdpmap arrays
!  dt      Time interval that parameterizes the parcel trajectory.
!  iterdp  Number of iterations used for departure point calculation.
!  wb      Vertical velocity component (sigma dot).
!  wst     z-derivs at the top edge of each interval contained in wb
!  wsb     z-derivs at the bot edge of each interval contained in wb
!  sig     Sigma values at the full-index levels.
!  sigh    Half-index sigma levels including sigh(1) = sigma(1/2) = 0.0
!          sigh(plev+1) = sigma(plev+1/2) = 1.0 .  Note that in general
!          sigh(k) .lt. sig(k)  where sig(k) is the sigma value at the
!          k_th full-index level.
!  dsigh   Increment in half-index sigma levels.
!  kdpmpf  Array of indices of the model full levels which are mapped
!          into an artificial evenly spaced vertical grid.  Used to aid
!          in search for vertical position of departure point 
!  kdpmph  Array of indices of the model half levels which are mapped
!          into an artificial evenly spaced vertical grid.  Used to aid
!          in search for vertical position of departure point 
!  sigmp   Sigma value at the trajectory midpoint for each gridpoint
!          in a vertical slice from the global grid.  On entry sigmp is
!          an initial guess.
!  sigdp   Sigma value at the trajectory endpoint for each gridpoint
!          in a vertical slice from the global grid.
!  kdp     Vertical index for each gridpoint.  This index points into a
!          vertical slice array whose vertical grid is given by sig.
!          E.g.,   sig(kdp(i,k)) .le. sigdp(i,k) .lt. sig(kdp(i,k)+1).
!-----------------------------------------------------------------------

!---------------------------Local variables-----------------------------
  integer i                 ! |
  integer iter              ! |-- indices
  integer k                 ! |
  real(r8) wmp(beglonxy:endlonxy,plev)   ! vert vel. at midpoint
!-----------------------------------------------------------------------
!
! Loop over departure point iterates.
!
  do iter = 1,iterdp
!
! Compute midpoint indices in half-index sigma-level arrays (use kdp
! as temporary storage).
!
     call kdpfnd1(plevp   ,pmap    ,sigh    ,sigmp   ,kdpmph  , &
                 kdp     ,nlon    )
!
! Interpolate sigma dot field to trajectory midpoints using Hermite
! cubic interpolant.
!
     call herzin(plevp   ,1       ,wb      ,wst     ,wsb     , &
                 sigh    ,dsigh   ,sigmp   ,kdp     ,wmp     , &
                 nlon    )
!
! Update estimate of trajectory midpoint.
!
!$OMP PARALLEL DO PRIVATE (K, I)
     do k = 1,plev
        do i = beglonxy,endlonxy
           sigmp(i,k) = sig(k) - .5_r8*dt*wmp(i,k)
        end do
     end do
!
! Restrict vertical midpoints to be between the top and bottom half-
! index sigma levels.
!
     call vdplim(plevp   ,sigh    ,sigmp   ,nlon)
  end do
!
! Compute trajectory endpoints.
!
!$OMP PARALLEL DO PRIVATE (K, I)
  do k = 1,plev
     do i = beglonxy,endlonxy
        sigdp(i,k) = sig(k) - dt*wmp(i,k)
     end do
  end do
!
! Restrict vertical departure points to be between the top and bottom
! full-index sigma levels.
!
  call vdplim(plev    ,sig     ,sigdp   ,nlon)
!
! Vertical indices for trajectory endpoints that point into full-index
! sigma level arrays.
!
  call kdpfnd1(plev    ,pmap    ,sig     ,sigdp   ,kdpmpf  , &
              kdp     ,nlon    )
!
  return
end subroutine vrtdep

!
!============================================================================================
!

subroutine vdplim(pkdim   ,sig     ,sigdp   ,nlon    )

!----------------------------------------------------------------------- 
! 
! Purpose: 
! Restrict vertical departure points to be between the top and bottom
! sigma levels of the "full-" or "half-" level grid
! 
! Method: 
! 
! Author: J. Olson
! 
!-----------------------------------------------------------------------
!---------------------- Arguments --------------------------------------
  integer , intent(in)    :: nlon               ! longitude dimension
  integer , intent(in)    :: pkdim              ! vertical dimension
  real(r8), intent(in)    :: sig(pkdim)         ! vertical coordinate of model grid
  real(r8), intent(inout) :: sigdp(beglonxy:endlonxy,plev)   ! vertical coords. of departure points.
! pkdim   Vertical dimension of "sig"
! sig     Sigma values at the "full" or "half" model levels
! sigdp   Sigma value at the trajectory endpoint or midpoint for each
!         gridpoint in a vertical slice from the global grid.  This
!         routine restricts those departure points to within the
!         model's vertical grid.
!-----------------------------------------------------------------------

!---------------------------Local variables-----------------------------
  integer i,k                 ! index
!-----------------------------------------------------------------------
!
!$OMP PARALLEL DO PRIVATE (K, I)
  do k=1,plev
     do i = beglonxy,endlonxy
        if (sigdp(i,k) < sig(1)) then
           sigdp(i,k) = sig(1)
        end if
        if (sigdp(i,k) >= sig(pkdim)) then
           sigdp(i,k) = sig(pkdim)*(1._r8 - 10._r8*epsilon(sigdp))
        end if
     end do
  end do

  return
end subroutine vdplim

!
!-----------------------------------------------------------------------
!

subroutine sltini(dlam,    sinlam,  coslam,  uxl,     uxr, &
                  vxl,     vxr,     qxl,     qxr,     adv_state )

!----------------------------------------------------------------------- 
! 
! Purpose: 
! Prepare the extended arrays for use in the SLT routines
!
!   1)  Fill latitude extensions.
!   2)  Fill longitude extensions.
!   3)  Compute x-derivatives
! 
! Method: 
! Computational note: The latitude loop in this routine is multitasked
! 
! Author: 
! Original version:  J. Olson
! Standardized:      J. Rosinski, June 1992
! Reviewed:          D. Williamson, P. Rasch, August 1992
! Reviewed:          D. Williamson, P. Rasch, March 1996
!
!-----------------------------------------------------------------------
   use prognostics,  only: ptimelevels
!-----------------------------------------------------------------------
#include <parslt.h>
!---------------------------Local parameters----------------------------
!
   integer puvpts            ! number of u/v pts in lat slice
   integer pqpts             ! number of constituent pts in lat slice
!
   parameter(puvpts = plond*plev, pqpts  = plond*plev*pcnst) 
!-----------------------------------------------------------------------
!
! Input arguments
!
   real(r8), intent(in) :: dlam(platd)          ! increment in x-direction
   real(r8), intent(in) :: sinlam(plond,platd)  ! sin(lamda)
   real(r8), intent(in) :: coslam(plond,platd)  ! cos(lamda)
   real(r8), intent(inout) :: uxl (plond,beglev:endlev,      beglatex:endlatex) 
   real(r8), intent(inout) :: uxr (plond,beglev:endlev,      beglatex:endlatex)  
   real(r8), intent(inout) :: vxl (plond,beglev:endlev,      beglatex:endlatex) 
   real(r8), intent(inout) :: vxr (plond,beglev:endlev,      beglatex:endlatex)  
   real(r8), intent(inout) :: qxl (plond,beglev:endlev,pcnst,beglatex:endlatex)  
   real(r8), intent(inout) :: qxr (plond,beglev:endlev,pcnst,beglatex:endlatex)
   type(advection_state), intent(inout) :: adv_state   ! Advection data state
!                                    
!
!-----------------------------------------------------------------------
!
!  dlam    Length of increment in longitude grid.
!  sinlam  Sin of longitudes in global grid (model grid pts only).
!  coslam  Cos of longitudes in global grid (model grid pts only).
!  uxl     x-derivatives of u at the left  (west) edge of given interval
!  vxl     x-derivatives of v at the left  (west) edge of given interval
!  uxr     x-derivatives of u at the right (east) edge of given interval
!  vxr     x-derivatives of v at the right (east) edge of given interval
!  qxl     x-derivatives of scalar species at the left  (west) edge
!          of given interval
!  qxr     x-derivatives of scalar species at the right (east) edge
!          of given interval
!
!---------------------------Local variables-----------------------------
!
   integer m,j,k             ! index
   integer nlond
!
!------------------------------Externals--------------------------------
!
   external cubxdr,extxx,extyss,extyvv,limdx
!
!-----------------------------------------------------------------------
!
! Fill latitude extensions beyond the southern- and northern-most
! latitudes in the global grid
!
   call t_startf ('slt_single')
   if (beglatex .le. endlatex) then
      call extyvv(1, plev, coslam, sinlam, adv_state%u3(1,beglev,beglatex), &
                                          adv_state%v3(1,beglev,beglatex))
      call extyss(pcnst, plev    ,adv_state%qminus, pcnst)
!
! Fill longitude extensions
!
      call extxx(1 ,plev    ,adv_state%u3(1,beglev,beglatex), 1)
      call extxx(1 ,plev    ,adv_state%v3(1,beglev,beglatex), 1)
      call extxx(pcnst, plev    ,adv_state%qminus, pcnst)
   endif
   call t_stopf ('slt_single')
!
! Compute x-derivatives.
!
#ifdef OUTER_OMP
!$OMP  PARALLEL DO PRIVATE (J, NLOND, K, M)
#endif
   do j = beglatex, endlatex
      nlond = 1 + 2*nxpt + nlonex(j)
!$OMP  PARALLEL DO PRIVATE (K, M)
      do k=beglev,endlev
         call cubxdr (nlond, 2, nlond-3, dlam(j), adv_state%u3(1,k,j), &
            uxl(1,k,j), uxr(1,k,j))
         call cubxdr (nlond, 2, nlond-3, dlam(j), adv_state%v3(1,k,j), &
            vxl(1,k,j), vxr(1,k,j))
         do m=1,pcnst
            call cubxdr (nlond, 2, nlond-3, dlam(j), adv_state%qminus(1,k,m,j), &
               qxl(1,k,m,j), qxr(1,k,m,j))
            if( plimdr )then
               call limdx (nlond, 2, nlond-3, dlam(j), adv_state%qminus(1,k,m,j), &
                  qxl(1,k,m,j), qxr(1,k,m,j))
            end if
         end do
      end do
   end do

   return
end subroutine sltini

!
!-----------------------------------------------------------------------
!

end module scanslt
