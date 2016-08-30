#include <define.h>

#ifdef CPL6

module colm_csmMod

#if (defined COUP_CSM)
!-----------------------------------------------------------------------
!BOP
!
! !MODULE: colm_csmMod
!
! !DESCRIPTION:
! Set of routines that define communication between the
! land model and flux coupler. The order of sends/receives is:
! 1) receive orbital data from coupler
! 2) send control data (grids and masks) to coupler
!    land grid does not have valid data, runoff grid does
! 3) receive valid land grid from flux coupler
! 4) send compressed runoff information to flux coupler
! 5) start normal send/recv communication patterm
!
! !USES:
  use precision , only : r8
  use nanMod
  use spmd      , only : masterproc=>p_master, mpicom=>p_comm 
  use mpiinc
  use cpl_fields_mod
  use cpl_contract_mod
  use cpl_interface_mod
  use RunoffMod        , only : runoff
  use shr_sys_mod      , only : shr_sys_irtc, shr_sys_flush ! csm_share system utility routines
  use system_messages  , only : allocation_err              ! allocation error output
  use abortutils, only : endrun
!
! !PUBLIC TYPES:
  implicit none
  save
!
! !PUBLIC MEMBER FUNCTIONS:
  public :: csm_setup          ! Setup, mpi_init
  public :: csm_shutdown       ! Shutdown, mpi_finalize
  public :: csm_initialize     ! Initialize contracts, etc
!*public :: csm_recvgrid       ! Receive grid and land mask (CCSM3.5 doesn't support -jidy-)
  public :: csm_dosndrcv       ! Logic for determining if send/recv
  public :: csm_recv           ! Receive data from flux coupler
  public :: csm_send           ! Send data to flux coupler
  public :: csm_sendalb        ! Send initial albedos, surface temp and snow data
  public :: csm_flxave         ! Flux averaging rougine
  public :: csm_restart        ! Restart code
  public :: compat_check_spval ! Checks that data sent from the coupler is valid
  public :: csm_compat         ! Checks compatibility of messages send/received

! !REVISION HISTORY:
!  02.09.17 Mariana Vertenstein Updated to clm2_1 data structures
!  03.01.15 T. Craig Update for cpl6
!  03.04.27 M. Vertenstein, added qref_2m to communication and
!           generalized global sums to include all fields
!
!EOP
!
  private
!
! PRIVATE MEMBER FUNCTIONS:
!
! PRIVATE TYPES:
!
  integer   :: nsend, nrecv, nroff            ! Buffer sizes
  integer   :: ibuffr(cpl_fields_ibuf_total)  ! Integer buffer from cpl
  integer   :: ibuffs(cpl_fields_ibuf_total)  ! Integer buffer to   cpl
  real(r8)  :: rbuffr(cpl_fields_rbuf_total)  ! Real    buffer from cpl
  real(r8)  :: rbuffs(cpl_fields_rbuf_total)  ! Real    buffer to   cpl
  type(cpl_contract)    :: contractRg         ! Contract for grid recvs from cpl
  type(cpl_contract)    :: contractR          ! Contract for recvs from cpl
  type(cpl_contract)    :: contractS          ! Contract for sends to   cpl
  type(cpl_contract)    :: contractSr         ! Contract for runoff sends to cpl
  real(r8), allocatable :: Gbuf(:,:)          ! Temporary generic buffer
  real(r8), allocatable :: bufS(:,:)          ! Send buffer for land
  real(r8), allocatable :: bufR(:,:)          ! Recv buffer for land
  real(r8), allocatable :: bufSr(:,:)         ! Send buffer for runoff
  real(r8), pointer     :: bufSglob(:,:)      ! Send global sum buffer for land
  real(r8), pointer     :: bufRglob(:,:)      ! Recv global sum buffer for land
  real(r8), pointer     :: bufSloc(:,:)       ! Send local sum buffer for land
  real(r8), pointer     :: bufRloc(:,:)       ! Recv local sum buffer for land
  real(r8), pointer     :: fieldS(:)          ! Global sum send field
  real(r8), pointer     :: fieldR(:)          ! Global sum receive field

  integer :: csm_nptg                         ! Loc sizes, grid coupling buffers
  integer :: csm_nptr                         ! Loc sizes, roff coupling buffers
  integer :: beg_lnd_rof,end_lnd_rof          ! beginning,ending landrunoff points
  integer :: beg_ocn_rof,end_ocn_rof          ! beginning,ending oceaan
!
!
! Flux averaging arrays and counters
!
  integer  :: icnt                         ! step counter for flux averager
  integer  :: ncnt                         ! number of steps over which to average output fluxes
  real(r8) :: rncnt                        ! reciprocal of ncnt

  real(r8), allocatable :: taux_ave(:)     ! averaged array
  real(r8), allocatable :: tauy_ave(:)     ! averaged array
  real(r8), allocatable :: lhflx_ave(:)    ! averaged array
  real(r8), allocatable :: shflx_ave(:)    ! averaged array
  real(r8), allocatable :: lwup_ave(:)     ! averaged array
  real(r8), allocatable :: qflx_ave(:)     ! averaged array
  real(r8), allocatable :: swabs_ave(:)    ! averaged array
  real(r8), allocatable :: nee_ave(:)      ! averaged array
!
! When to send/receive messages to coupler and when to make restart and stop
!
  integer, private:: ncpday         ! number of send/recv calls per day
  logical, public :: dorecv         ! receive data from coupler this step
  logical, public :: dosend         ! send data to coupler this step
  logical, public :: csmstop_next   ! received stop at eod signal and will stop on next ts
  logical, public :: csmstop_now    ! received stop now signal from coupler
  logical, public :: csmrstrt       ! restart write signal received from coupler
!
! Indices for send/recv fields
!
! Indices for send/recv fields

  integer :: index_l2c_Sl_t        ! temperature
  integer :: index_l2c_Sl_tref     ! 2m reference temperature
  integer :: index_l2c_Sl_qref     ! 2m reference specific humidity
!fengjm
  integer :: index_l2c_Sl_u10m     ! 10m u-wind
  integer :: index_l2c_Sl_v10m     ! 10m v-wind
  integer :: index_l2c_Fall_sublim ! sublimation
!fengjm
  integer :: index_l2c_Sl_avsdr    ! albedo: direct , visible
  integer :: index_l2c_Sl_avsdf    ! albedo: diffuse, visible
  integer :: index_l2c_Sl_anidr    ! albedo: direct , near-ir
  integer :: index_l2c_Sl_anidf    ! albedo: diffuse, near-ir
  integer :: index_l2c_Sl_snowh    ! snow height
  integer :: index_l2c_Fall_taux   ! wind stress, zonal
  integer :: index_l2c_Fall_tauy   ! wind stress, meridional
  integer :: index_l2c_Fall_lat    ! latent          heat flux
  integer :: index_l2c_Fall_sen    ! sensible        heat flux
  integer :: index_l2c_Fall_lwup   ! upward longwave heat flux
  integer :: index_l2c_Fall_evap   ! evaporation    water flux
  integer :: index_l2c_Fall_swnet  ! solar radiation absorbed (total)
  integer :: index_l2c_Fall_nee    ! co2 flux

  integer :: index_c2l_Sa_co2prog  ! bottom atm prognostic co2  ***
  integer :: index_c2l_Sa_co2diag  ! bottom atm diagnostic co2  ***
  integer :: index_c2l_Sa_z        ! bottom atm level height
  integer :: index_c2l_Sa_u        ! bottom atm level zon wind
  integer :: index_c2l_Sa_v        ! bottom atm level mer wind
  integer :: index_c2l_Sa_tbot     ! bottom atm level temp
  integer :: index_c2l_Sa_ptem     ! bottom atm level pot temp
  integer :: index_c2l_Sa_shum     ! bottom atm level spec hum
  integer :: index_c2l_Sa_dens     ! bottom atm level air dens  ***
  integer :: index_c2l_Sa_pbot     ! bottom atm level pressure
  integer :: index_c2l_Sa_pslv     ! sea level atm pressure     ***
  integer :: index_c2l_Faxa_lwdn   ! downward longwave heat flux
  integer :: index_c2l_Faxa_rainc  ! precip: liquid, convective
  integer :: index_c2l_Faxa_rainl  ! precip: liquid, large-scale
  integer :: index_c2l_Faxa_snowc  ! precip: frozen, convective
  integer :: index_c2l_Faxa_snowl  ! precip: frozen, large-scale
  integer :: index_c2l_Faxa_swndr  ! shortwave: nir direct  down
  integer :: index_c2l_Faxa_swvdr  ! shortwave: vis direct  down
  integer :: index_c2l_Faxa_swndf  ! shortwave: nir diffuse down
  integer :: index_c2l_Faxa_swvdf  ! shortwave: vis diffuse down

  integer :: index_r2c_Forr_roff   ! runoff to ocean

!
! CCSM timers
!
  logical  :: timer_lnd_sendrecv = .false. ! true => timer is on
  logical  :: timer_lnd_recvsend = .false. ! true => timer is on
!
! Input dtime
!
  integer, public :: csm_dtime   ! value passed to coupler on initialization set to namelist input

! !PRIVATE MEMBER FUNCTIONS:
  private :: global_sum_fld2d    ! global sum of 2d grid fields
  private :: global_sum_fld1d    ! global sum of 1d fluxes
!-----------------------------------------------------------------------

contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: csm_setup
!
! !INTERFACE:
  subroutine csm_setup(mpicom)
!
! !DESCRIPTION:
!  Initialize csm coupling, return the communicator group to the
!  application.
!
! !ARGUMENTS:
    implicit none
    integer, intent(out) :: mpicom  !MPI group communicator
!
! !REVISION HISTORY:
!  03.01.15 T. Craig: initial version
!
!EOP
!
! !LOCAL VARIABLES:
!------------------------------------------------------------------------

   call cpl_interface_init(cpl_fields_lndname,mpicom)

 end subroutine csm_setup

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: csm_shutdown
!
! !INTERFACE:
  subroutine csm_shutdown
!
! !DESCRIPTION:
!  Finalize csm coupling
!
! !ARGUMENTS:
    implicit none
!
! !REVISION HISTORY:
!  03.01.15 T. Craig: initial version
!
!EOP
!
! !LOCAL VARIABLES:
!------------------------------------------------------------------------

   call cpl_interface_finalize(cpl_fields_lndname)

 end subroutine csm_shutdown

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: csm_initialize
!
! !INTERFACE:
 subroutine csm_initialize(irad, eccen, obliqr, lambm0, mvelpp)
!
! !DESCRIPTION:
!  Send first control data to flux coupler and "invalid" grid
!  containing special value data.
!  The coupler treats points where the mask is nonzero as points where
!  you could possibly do a calculation (in the case of the runoff, this
!  corresponds to all rtm ocean points). The coupler then defines a "key"
!  as points where the model can give you valid data (in the case of runoff,
!  this corresponds to points where the land model will give you valid
!  compressed data points). The key can be 0 where the mask is 1. However,
!  the key cannot be 1 where the mask is 0 unless the data is also zero.
!  In the case of runoff, the key the coupler builds is time invariant.
!  Send first control data to flux coupler and "invalid" grid
!  containing special value data
!
! !USES:
    use phycon_module, only : re
    use shr_const_mod, only : SHR_CONST_CDAY
    use RunoffMod    , only : get_proc_rof_bounds, runoff
    use RtmMod       , only : area_r, longxy_r, latixy_r, mask_r
    use colm_rtmVar  , only : rtmlon, rtmlat
    use colm_varMod  , only : lsmlon=>lon_points, lsmlat=>lat_points, numgrid, numgrid_glob, &
                              area,latixy, longxy, landmask, landfrac
    use colm_varctl  , only : csm_doflxave
    use spmd_decomp  , only : gxmap, gymap
    use timemgr      , only : get_step_size
    use colm_varctl  , only : nsrest
!
! !ARGUMENTS:
    implicit none
    integer , intent(in)  :: irad   ! frequency of radiation computation
    real(r8), intent(out) :: eccen  ! Earth's eccentricity of orbit
    real(r8), intent(out) :: obliqr ! Earth's obliquity in radians
    real(r8), intent(out) :: lambm0 ! Mean longitude of perihelion at the vernal equinox (radians)
    real(r8), intent(out) :: mvelpp ! Earth's moving vernal equinox long of perihelion plus pi (radians)
!
! !REVISION HISTORY:
!  02.09.17 Mariana Vertenstein Updated to clm2_1 data structures
!  03.01.15 T.Craig Update for cpl6
!
!EOP
!
! !LOCAL VARIABLES:
    integer :: i,j,g,gi,n,ni              ! indices
    real(r8):: dtime                      ! step size
    real(r8):: spval                      ! special value
!------------------------------------------------------------------------

    ! Determine processor bounds

!   call get_proc_bounds(begg, endg, begl, endl, begc, endc, begp, endp)
!   call get_proc_global(numg, numl, numc, nump)
    call get_proc_rof_bounds(beg_lnd_rof, end_lnd_rof, beg_ocn_rof, end_ocn_rof)

    ! Determine number of send/recv calls steps per day to flux coupler

    if (nsrest == 0) then
       dtime = get_step_size()
    else
       dtime = csm_dtime
    endif

    if (csm_doflxave) then
       ncpday = nint(SHR_CONST_CDAY/dtime)/irad
    else
       ncpday = nint(SHR_CONST_CDAY/dtime)
    endif

    ! Setup contracts for grid communication

    ibuffs(:) = 0                                    ! initialize ibuffs
    ibuffs(cpl_fields_ibuf_gsize  ) = lsmlon*lsmlat  ! global array size
    ibuffs(cpl_fields_ibuf_gisize ) = lsmlon         ! global number of lons
    ibuffs(cpl_fields_ibuf_gjsize ) = lsmlat         ! global number of lats
    ibuffs(cpl_fields_ibuf_nfields) = cpl_fields_grid_total
    ibuffs(cpl_fields_ibuf_ncpl   ) = ncpday         ! number of land send/recv calls per day
  ! ibuffs(cpl_fields_ibuf_inimask) = 1              ! T(or 1) => requests cpl to send valid land domain info back

  ! ibuffs(cpl_fields_ibuf_lsize  ) = endg-begg+1
  ! ibuffs(cpl_fields_ibuf_lisize ) = endg-begg+1
    ibuffs(cpl_fields_ibuf_lsize  ) = numgrid
    ibuffs(cpl_fields_ibuf_lisize ) = numgrid
    ibuffs(cpl_fields_ibuf_ljsize ) = 1

  ! allocate(Gbuf(begg:endg,cpl_fields_grid_total))
    allocate(Gbuf(numgrid,cpl_fields_grid_total))

  ! do g = begg,endg
    do g = 1, numgrid
       i = gxmap(g)
       j = gymap(g)
       Gbuf(g,cpl_fields_grid_lon  ) = longxy(i,j)
       Gbuf(g,cpl_fields_grid_lat  ) = latixy(i,j)
       Gbuf(g,cpl_fields_grid_area ) = area(i,j)/(re*re)
       Gbuf(g,cpl_fields_grid_frac ) = landfrac(i,j)
       Gbuf(g,cpl_fields_grid_mask ) = float(landmask(i,j))
     ! gi = (gptr%jxy(g)-1)*lsmlon + gptr%ixy(g)
       Gbuf(g,cpl_fields_grid_index) = (j-1)*lsmlon+i
    end do

    call cpl_interface_contractInit(contractS, cpl_fields_lndname,cpl_fields_cplname,cpl_fields_l2c_fields, ibuffs,Gbuf)
    call cpl_interface_contractInit(contractR, cpl_fields_lndname,cpl_fields_cplname,cpl_fields_c2l_fields, ibuffs,Gbuf)

    deallocate(Gbuf)

    ! Setup contracts for  runoff communication

    ibuffs(:) = 0
    ibuffs(cpl_fields_ibuf_gsize  ) = rtmlon*rtmlat                 ! global array size
    ibuffs(cpl_fields_ibuf_gisize ) = rtmlon                        ! global number of lons
    ibuffs(cpl_fields_ibuf_gjsize ) = rtmlat                        ! global number of lats
    ibuffs(cpl_fields_ibuf_ncpl   ) = ncpday                        ! number of land send/recv calls per day
    ibuffs(cpl_fields_ibuf_lsize  ) = end_ocn_rof - beg_ocn_rof + 1 ! local array size
    ibuffs(cpl_fields_ibuf_lisize ) = end_ocn_rof - beg_ocn_rof + 1 ! local array size
    ibuffs(cpl_fields_ibuf_ljsize ) = 1                             ! local array size

    csm_nptr = ibuffs(cpl_fields_ibuf_lsize)

    allocate(Gbuf(csm_nptr,cpl_fields_grid_total))

    do n = beg_ocn_rof,end_ocn_rof
       ni = n - beg_ocn_rof+ 1
       gi = (runoff%ocn_jxy(n)-1)*rtmlon + runoff%ocn_ixy(n)
       j = (gi-1) / rtmlon + 1
       i = mod(gi-1,rtmlon) + 1
       Gbuf(ni,cpl_fields_grid_lon  ) = longxy_r(i,j)
       Gbuf(ni,cpl_fields_grid_lat  ) = latixy_r(i,j)
       Gbuf(ni,cpl_fields_grid_area ) = area_r(i,j)/(re*re)
       Gbuf(ni,cpl_fields_grid_mask ) = 1.0 - float(mask_r(i,j))
       Gbuf(ni,cpl_fields_grid_index) = gi
    end do

    call cpl_interface_contractInit(contractSr,cpl_fields_lndname,cpl_fields_cplname,cpl_fields_r2c_fields, ibuffs,Gbuf)

    deallocate(Gbuf)

    ! Receive initial ibuf message

    call cpl_interface_ibufRecv(cpl_fields_cplname,ibuffr,rbuffr)

    spval  = rbuffr(cpl_fields_rbuf_spval)
    eccen  = rbuffr(cpl_fields_rbuf_eccen)
    obliqr = rbuffr(cpl_fields_rbuf_obliqr)
    lambm0 = rbuffr(cpl_fields_rbuf_lambm0)
    mvelpp = rbuffr(cpl_fields_rbuf_mvelpp)

    ! Check that data is good data and not the special value

    if (masterproc) then
       call compat_check_spval(spval, eccen ,'Eccentricity'     )
       call compat_check_spval(spval, obliqr,'Obliquity'        )
       call compat_check_spval(spval, lambm0,'Long of perhelion')
       call compat_check_spval(spval, mvelpp,'Move long of perh')

       write(6,*)'(CSM_INITIALIZE): eccen:  ', eccen
       write(6,*)'(CSM_INITIALIZE): obliqr: ', obliqr
       write(6,*)'(CSM_INITIALIZE): lambm0: ', lambm0
       write(6,*)'(CSM_INITIALIZE): mvelpp: ', mvelpp
    end if

    write(6,*)'(CSM_INITIALIZE): there will be ',ncpday, &
         ' send/recv calls per day from the land model to the flux coupler'
    write(6,*)'(CSM_INITIALIZE):sent l->d control data '

    ! Allocate memory for grid and runoff communication

    nsend = cpl_interface_contractNumatt(contractS)
    allocate(bufS(numgrid,nsend))

    nrecv = cpl_interface_contractNumatt(contractR)
    allocate(bufR(numgrid,nrecv))

    nroff = cpl_interface_contractNumatt(contractSr)
    allocate(bufSr(csm_nptr,nroff))

    !---------------------------------------------------------------
    ! Determine indices
    !---------------------------------------------------------------

    ! Determine send indices

    index_l2c_Sl_t       = cpl_interface_contractIndex(contractS,'Sl_t')
    index_l2c_Sl_tref    = cpl_interface_contractIndex(contractS,'Sl_tref')
    index_l2c_Sl_qref    = cpl_interface_contractIndex(contractS,'Sl_qref')
!fengjm
    index_l2c_Sl_u10m    = cpl_interface_contractIndex(contractS,'Sl_u10m')
    index_l2c_Sl_v10m    = cpl_interface_contractIndex(contractS,'Sl_v10m')
    index_l2c_Fall_sublim= cpl_interface_contractIndex(contractS,'Fall_sublim')
!fengjm
    index_l2c_Sl_avsdr   = cpl_interface_contractIndex(contractS,'Sl_avsdr')
    index_l2c_Sl_avsdf   = cpl_interface_contractIndex(contractS,'Sl_avsdf')
    index_l2c_Sl_anidr   = cpl_interface_contractIndex(contractS,'Sl_anidr')
    index_l2c_Sl_anidf   = cpl_interface_contractIndex(contractS,'Sl_anidf')
    index_l2c_Sl_snowh   = cpl_interface_contractIndex(contractS,'Sl_snowh')
    index_l2c_Fall_taux  = cpl_interface_contractIndex(contractS,'Fall_taux')
    index_l2c_Fall_tauy  = cpl_interface_contractIndex(contractS,'Fall_tauy')
    index_l2c_Fall_lat   = cpl_interface_contractIndex(contractS,'Fall_lat')
    index_l2c_Fall_sen   = cpl_interface_contractIndex(contractS,'Fall_sen')
    index_l2c_Fall_lwup  = cpl_interface_contractIndex(contractS,'Fall_lwup')
    index_l2c_Fall_evap  = cpl_interface_contractIndex(contractS,'Fall_evap')
    index_l2c_Fall_swnet = cpl_interface_contractIndex(contractS,'Fall_swnet')
    index_l2c_Fall_nee   = cpl_interface_contractIndex(contractS,'Fall_nee', perrwith='quiet')

    ! Determine receive indices

    index_c2l_Sa_z       = cpl_interface_contractIndex(contractR,'Sa_z')
    index_c2l_Sa_u       = cpl_interface_contractIndex(contractR,'Sa_u')
    index_c2l_Sa_v       = cpl_interface_contractIndex(contractR,'Sa_v')
    index_c2l_Sa_tbot    = cpl_interface_contractIndex(contractR,'Sa_tbot')
    index_c2l_Sa_ptem    = cpl_interface_contractIndex(contractR,'Sa_ptem')
    index_c2l_Sa_shum    = cpl_interface_contractIndex(contractR,'Sa_shum')
    index_c2l_Sa_dens    = cpl_interface_contractIndex(contractR,'Sa_dens')
    index_c2l_Sa_pbot    = cpl_interface_contractIndex(contractR,'Sa_pbot')
    index_c2l_Sa_pslv    = cpl_interface_contractIndex(contractR,'Sa_pslv')
    index_c2l_Faxa_lwdn  = cpl_interface_contractIndex(contractR,'Faxa_lwdn')
    index_c2l_Faxa_rainc = cpl_interface_contractIndex(contractR,'Faxa_rainc')
    index_c2l_Faxa_rainl = cpl_interface_contractIndex(contractR,'Faxa_rainl')
    index_c2l_Faxa_snowc = cpl_interface_contractIndex(contractR,'Faxa_snowc')
    index_c2l_Faxa_snowl = cpl_interface_contractIndex(contractR,'Faxa_snowl')
    index_c2l_Faxa_swndr = cpl_interface_contractIndex(contractR,'Faxa_swndr')
    index_c2l_Faxa_swvdr = cpl_interface_contractIndex(contractR,'Faxa_swvdr')
    index_c2l_Faxa_swndf = cpl_interface_contractIndex(contractR,'Faxa_swndf')
    index_c2l_Faxa_swvdf = cpl_interface_contractIndex(contractR,'Faxa_swvdf')
    index_c2l_Sa_co2prog = cpl_interface_contractIndex(contractR,'Sa_co2prog', perrwith='quiet')
    index_c2l_Sa_co2diag = cpl_interface_contractIndex(contractR,'Sa_co2diag', perrwith='quiet')

    ! Determine runoff (send) indices

    index_r2c_Forr_roff  = cpl_interface_contractIndex(contractSR,'Forr_roff')

  end subroutine csm_initialize

!*!-----------------------------------------------------------------------
!*!BOP
!*!
!*! !IROUTINE: csm_recvgrid
!*!
!*! !INTERFACE:
!*  subroutine csm_recvgrid(cam_longxy, cam_latixy, cam_numlon, &
!*                          cam_landfrac, cam_landmask)
!*!
!*! !DESCRIPTION:
!*!  Receive valid land grid and land mask from coupler
!*!
!*! !USES:
!*    use colm_varMod, only : lsmlon=>lon_points, lsmlat=>lat_points
!*
!*! !ARGUMENTS:
!*    implicit none
!*    integer , intent(out) :: cam_numlon(lsmlat)           !cam number of longitudes
!*    real(r8), intent(out) :: cam_longxy(lsmlon,lsmlat)    !cam lon values
!*    real(r8), intent(out) :: cam_latixy(lsmlon,lsmlat)    !cam lat values
!*    real(r8), intent(out) :: cam_landfrac(lsmlon,lsmlat)  !cam fractional land
!*    integer , intent(out) :: cam_landmask(lsmlon,lsmlat)  !cam land mask
!*!
!*! !REVISION HISTORY:
!*!  02.09.17 Mariana Vertenstein Updated to clm2_1 data structures
!*!
!*!EOP
!*!
!*! !LOCAL VARIABLES:
!*    integer  :: i,j,n                    ! indices
!*    real(r8) :: area_a(lsmlon,lsmlat)    ! coupler atm grid areas
!*    integer  :: mask_a(lsmlon,lsmlat)    ! coupler atm valid grid mask
!*    real(r8) :: tmp(lsmlon,lsmlat)       ! temporary
!*!------------------------------------------------------------------------
!*
!*    ! Set integer control information and setup contracts
!*
!*    ibuffs(:)  = 0                                   !initialize ibuffs
!*    ibuffs(cpl_fields_ibuf_gsize  ) = lsmlon*lsmlat  !global array size
!*    ibuffs(cpl_fields_ibuf_gisize ) = lsmlon         !global number of lons
!*    ibuffs(cpl_fields_ibuf_gjsize ) = lsmlat         !global number of lats
!*    ibuffs(cpl_fields_ibuf_nfields) = cpl_fields_grid_total
!*    ibuffs(cpl_fields_ibuf_inimask) = 1              !T(or 1) => requests cpl to send valid land domain info back
!*    ibuffs(cpl_fields_ibuf_lsize  ) = lsmlon*lsmlat  !local array size
!*    ibuffs(cpl_fields_ibuf_lisize ) = lsmlon         !local array size
!*    ibuffs(cpl_fields_ibuf_ljsize ) = lsmlat         !local array size
!*
!*    csm_nptg = ibuffs(cpl_fields_ibuf_lsize)
!*
!*    allocate(Gbuf(csm_nptg,cpl_fields_grid_total))
!*    do n = 1,csm_nptg
!*       Gbuf(n,cpl_fields_grid_lon  ) = 1.0e30
!*       Gbuf(n,cpl_fields_grid_lat  ) = 1.0e30
!*       Gbuf(n,cpl_fields_grid_area ) = 1.0e30
!*       Gbuf(n,cpl_fields_grid_mask ) = 1.0e30
!*       Gbuf(n,cpl_fields_grid_index) = n
!*    enddo
!*    call cpl_interface_contractInit(contractRg,cpl_fields_lndname,cpl_fields_cplname,cpl_fields_c2lg_fields,ibuffs,Gbuf)
!*    deallocate(Gbuf)
!*
!*    allocate(Gbuf(csm_nptg,cpl_fields_c2lg_total))
!*    call cpl_interface_contractRecv(cpl_fields_cplname,contractRg,ibuffr,Gbuf)
!*    cam_numlon(:) = 0
!*    do n = 1,csm_nptg
!*       j = (n-1) / lsmlon + 1
!*       i = mod(n-1,lsmlon) + 1
!*       cam_longxy(i,j)   = Gbuf(n,cpl_fields_c2lg_alon)
!*       cam_latixy(i,j)   = Gbuf(n,cpl_fields_c2lg_alat)
!*       area_a(i,j)       = Gbuf(n,cpl_fields_c2lg_aarea)
!*       cam_landfrac(i,j) = Gbuf(n,cpl_fields_c2lg_lfrac)
!*       cam_landmask(i,j) = nint(Gbuf(n,cpl_fields_c2lg_lmask))
!*       mask_a(i,j)       = nint(Gbuf(n,cpl_fields_c2lg_amask))
!*       if (mask_a(i,j) /= 0) cam_numlon(j) = cam_numlon(j)+1
!*    end do
!*    deallocate(Gbuf)
!*
!*    if (masterproc) then
!*       write(6,*)'(CSM_SENDGRID):recd d->l land grid '
!*       write(6,100) 'lnd','recv', cpl_fields_c2lg_alon, global_sum_fld2d(cam_longxy  , 1.e30), ' lon'
!*       write(6,100) 'lnd','recv', cpl_fields_c2lg_alat, global_sum_fld2d(cam_latixy  , 1.e30), ' lat'
!*       write(6,100) 'lnd','recv', cpl_fields_c2lg_aarea,global_sum_fld2d(area_a      , 1.e30), ' aarea'
!*       write(6,100) 'lnd','recv', cpl_fields_c2lg_lfrac,global_sum_fld2d(cam_landfrac, 1.e30), ' lfrac'
!*       tmp = float(cam_landmask)
!*       write(6,100) 'lnd','recv', cpl_fields_c2lg_lmask,global_sum_fld2d(tmp, 1.e30), ' lmask'
!*       tmp = float(mask_a)
!*       write(6,100) 'lnd','recv', cpl_fields_c2lg_amask,global_sum_fld2d(tmp, 1.e30), ' amask'
!*100    format('comm_diag ',a3,1x,a4,1x,i3,es26.19,a)
!*    endif
!*
!*  end subroutine csm_recvgrid

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: csm_sendalb
!
! !INTERFACE:
  subroutine csm_sendalb()
!
! !DESCRIPTION:
! Send initial albedos, surface temperature and snow data to the
! flux coupler
!
! !USES:
!   use clm_varctl  , only : csm_doflxave, nsrest
!   use clm_varcon  , only : sb
!   use lnd2atmMod  , only : lnd2atm

    use colm_cplMod  , only : lnd_avsdr, lnd_avsdf, lnd_anidr, lnd_anidf, lnd_trad, lnd_scv
    use phycon_module, only : sb
    use timemgr      , only : get_curr_date, get_prev_date
    use colm_varctl  , only : nsrest
    use colm_varMod  , only : numgrid
!
! !ARGUMENTS:
    implicit none
!
! !REVISION HISTORY:
!  02.09.17 Mariana Vertenstein Updated to clm2_1 data structures
!
!EOP
!
! !LOCAL VARIABLES:
    integer :: i,j,g,gi,n,ni ! generic indices
    integer :: yr            ! current year
    integer :: mon           ! current month
    integer :: day           ! current day (0, 1, ...)
    integer :: ncsec         ! current seconds of current date (0, ..., 86400)
    integer :: ncdate        ! current date (yymmdd format) (e.g., 021105)
  ! type(gridcell_type), pointer :: gptr  ! pointer to gridcell derived subtype
! -----------------------------------------------------------------

    ! Fill send buffer for grid data

    bufS(:,:) = 0.0

    if (nsrest == 0) then   !initial run

       ! On initial timestep ONLY: determine 1d vector of states that will be sent
       ! to coupler and map fields from 1d subgrid vector to 2d [lsmlon]x[lsmlat] grid.

!      call lnd2atm(init=.true.)

!      do g = begg,endg
!         bufS(g,index_l2c_Sl_t ) = sqrt(sqrt(gptr%l2af%eflx_lwrad_out(g)/sb))
!         bufS(g,index_l2c_Sl_snowh) = gptr%l2as%h2osno(g)
!         bufS(g,index_l2c_Sl_avsdr) = gptr%l2as%albd(g,1)
!         bufS(g,index_l2c_Sl_anidr) = gptr%l2as%albd(g,2)
!         bufS(g,index_l2c_Sl_avsdf) = gptr%l2as%albi(g,1)
!         bufS(g,index_l2c_Sl_anidf) = gptr%l2as%albi(g,2)
!      end do

       do g = 1, numgrid
          bufS(g,index_l2c_Sl_t    ) = lnd_trad(g)
          bufS(g,index_l2c_Sl_snowh) = lnd_scv(g)
          bufS(g,index_l2c_Sl_avsdr) = lnd_avsdr(g)  ! avsdr
          bufS(g,index_l2c_Sl_avsdf) = lnd_avsdf(g)  ! avsdf
          bufS(g,index_l2c_Sl_anidr) = lnd_anidr(g)  ! anidr
          bufS(g,index_l2c_Sl_anidf) = lnd_anidf(g)  ! anidf
       end do

    else  ! restart run

       ! On a restart run, no meaningful data is sent to the flux coupler -
       ! this includes ocean runoff (which should only contain zero values)
       ! since the runoff code (riverfluxrtm) has not been called yet

       bufS(:,:) = 1.e30

    endif

    ! Determine time index to send to coupler. Note that for a restart run,
    ! the next time step is nstep+1. But must send current time step to flux couper here.

    if (nsrest == 0) then
       call get_curr_date (yr, mon, day, ncsec)
    else
       call get_prev_date (yr, mon, day, ncsec)
    endif

    ncdate = yr*10000 + mon*100 + day

    ibuffs(:)  = 0
    ibuffs(cpl_fields_ibuf_cdate)  = ncdate      !model date (yyyymmdd)
    ibuffs(cpl_fields_ibuf_sec  )  = ncsec       !elapsed seconds in current date

    ! Send grid data to coupler

    call cpl_interface_contractSend(cpl_fields_cplname,contractS,ibuffs,bufS)

    ! Send runoff data to coupler
    ! Must convert runoff to units of kg/m^2/s from m^3/s

    bufSr(:,:) = 0.
    do n = beg_ocn_rof, end_ocn_rof
       ni = n - beg_ocn_rof+ 1
       bufSr(ni,index_r2c_Forr_roff) = runoff%ocn(n) / (runoff%ocn_area(n)*1000.)
    end do
    call cpl_interface_contractSend(cpl_fields_cplname,contractSr,ibuffs,bufSr)

  end subroutine csm_sendalb

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: csm_dosndrcv
!
! !INTERFACE:
  subroutine csm_dosndrcv(doalb)
!
! !DESCRIPTION:
! Determine when to send and receive messages to/from the
! flux coupler on this time-step.
! Determine if send/receive information to/from flux coupler
! Send msgs (land state and fluxes) to the flux coupler only when
! doalb is true (i.e. on time steps before the atm does a solar
! radiation computation). Receive msgs (atm state) from the
! flux coupler only when dorad is true (i.e. on time steps
! when the atm does a solar radiation computation).
! The fluxes are then averaged between the send and receive calls.
!
! !USES:
    use colm_varctl  , only : csm_doflxave
    use shr_const_mod, only : SHR_CONST_CDAY
    use timemgr      , only : get_nstep, get_step_size
!
! !ARGUMENTS:
    implicit none
    logical, intent(in) :: doalb  !true=>next timestep a radiation time step
!
! !REVISION HISTORY:
!  02.09.17 Mariana Vertenstein Updated to clm2_1 data structures
!
!EOP
!
! !LOCAL VARIABLES:
    integer  :: ntspday           !model steps per day
    real(r8) :: dtime             !step size (seconds)
    integer  :: nstep             !time step
!-----------------------------------------------------------------------

    ! Determine if send/receive information to/from flux coupler

    nstep = get_nstep()
    if (csm_doflxave) then
       if (nstep == 0) then
          dorecv = .true.
          dosend = .false.
       else if (nstep == 1) then
          dorecv = .false.
          dosend = doalb
       else
          dorecv = dosend
          dosend = doalb
       endif
    else
       if (nstep == 0) then
          dorecv = .true.
          dosend = .false.
       else if (nstep == 1) then
          dorecv = .false.
          dosend = .true.
       else
          dorecv = .true.
          dosend = .true.
       endif
    endif

    ! If at end of day: check if should write restart file or stop
    ! at next time step. Note, these statements must appear here since
    ! ibuffr is not received at every time step when flux averaging occurs.

    csmstop_next = .false.
    csmrstrt     = .false.
    dtime        = get_step_size()
    ntspday      = nint(SHR_CONST_CDAY/dtime)
    if (mod(nstep,ntspday) == 0) then
       if (ibuffr(cpl_fields_ibuf_stopeod) /= 0) then  !stop at end of day
          csmstop_next = .true.  !will stop on next time step
          write(6,*)'(CSM_DOSNDRCV) output restart and history files at nstep = ',nstep
       endif
       if (ibuffr(cpl_fields_ibuf_resteod) /= 0) then !write restart at end of day
          csmrstrt = .true.      !will write restart now
          write(6,*)'(CSM_DOSNDRCV) output restart and history files at nstep = ',nstep
       endif
    endif

  end subroutine csm_dosndrcv

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: csm_recv
!
! !INTERFACE:
  subroutine csm_recv()
!
! !DESCRIPTION:
!  Receive and map data from flux coupler
!
! !USES:
    use phycon_module, only : rair
    use colm_varMod, only : numgrid, numgrid_glob
    use colm_cplMod, only : forcg
    use colm_varctl, only : co2_option, po2, pco2
    use timemgr, only : idate
    use spmdGathScatMod
!
! !ARGUMENTS:
    implicit none
!
! !REVISION HISTORY:
!  02.09.17 Mariana Vertenstein Updated to clm2_1 data structures
!
!EOP
!
! !LOCAL VARIABLES:
    integer :: i,j,g,gi,n,ni ! generic indices
    real(r8):: forc_rainc    ! rainxy Atm flux mm/s
    real(r8):: forc_rainl    ! rainxy Atm flux mm/s
    real(r8):: forc_snowc    ! snowfxy Atm flux  mm/s
    real(r8):: forc_snowl    ! snowfxl Atm flux  mm/s
    real(r8):: co2_ppmv_diag ! temporary
    real(r8):: co2_ppmv_prog ! temporary
    real(r8):: co2_ppmv_val  ! temporary
    integer :: ier           ! return error code
  ! type(gridcell_type), pointer :: gptr  ! pointer to gridcell derived subtype
!-----------------------------------------------------------------------

    ! Start timers

   ! if (timer_lnd_sendrecv) then
   !    call t_stopf ('lnd_sendrecv') ; timer_lnd_sendrecv = .false.
   ! endif

   ! call t_startf('lnd_recv')

     ibuffr(:) = 0

     call cpl_interface_contractRecv(cpl_fields_cplname,contractR,ibuffr,bufR)

     ! Do global integrals of fluxes if flagged

     if (ibuffr(cpl_fields_ibuf_infobug) >= 2) then
        if (masterproc) then
           if (.not. associated(bufRglob)) then
            ! allocate(bufRglob(nrecv,numg), stat=ier)
              allocate(bufRglob(nrecv,numgrid_glob), stat=ier)
              if (ier /= 0) then
                 write(6,*)'clm_csmMod: allocation error for bufRglob'; call endrun
              end if
              if (.not. associated(fieldR)) then
               ! allocate(fieldR(numg), stat=ier)
                 allocate(fieldR(numgrid_glob), stat=ier)
                 if (ier /= 0) then
                    write(6,*)'clm_csmMod: allocation error for fieldR'; call endrun
                 end if
              endif
           end if
        end if
        if (.not. associated(bufRloc)) then
         ! allocate(bufRloc(nrecv,begg:endg), stat=ier)
           allocate(bufRloc(nrecv,numgrid), stat=ier)
           if (ier /= 0) then
              write(6,*)'clm_csmMod: allocation error for bufRloc'; call endrun
           end if
        end if
      ! do g = begg,endg
        do g = 1, numgrid
           do n = 1,nrecv
              bufRloc(n,g) = bufR(g,n)
           end do
        end do
#if (defined SPMD)
        call gather_data_to_master(bufRloc, bufRglob, clmlevel='gridcell')
#else
        bufRglob(:,:) = bufRloc(:,:)
#endif
        if (masterproc) then
           write(6,*)

           if (index_c2l_Sa_co2prog /= 0) then
              fieldR(:) = bufRglob(index_c2l_Sa_co2prog,:)
              write(6,100) 'lnd','recv', index_c2l_Sa_co2prog, global_sum_fld1d(fieldR), ' co2prog'
           end if

           if (index_c2l_Sa_co2diag /= 0) then
              fieldR(:) = bufRglob(index_c2l_Sa_co2diag,:)
              write(6,100) 'lnd','recv', index_c2l_Sa_co2diag, global_sum_fld1d(fieldR), ' co2diag'
           end if

           fieldR(:) = bufRglob(index_c2l_Sa_z,:)
           write(6,100) 'lnd','recv', index_c2l_Sa_z      , global_sum_fld1d(fieldR), ' hgt'
           fieldR(:) = bufRglob(index_c2l_Sa_u,:)
           write(6,100) 'lnd','recv', index_c2l_Sa_u      , global_sum_fld1d(fieldR), ' u'
           fieldR(:) = bufRglob(index_c2l_Sa_v,:)
           write(6,100) 'lnd','recv', index_c2l_Sa_v      , global_sum_fld1d(fieldR), ' v'
           fieldR(:) = bufRglob(index_c2l_Sa_ptem,:)
           write(6,100) 'lnd','recv', index_c2l_Sa_ptem   , global_sum_fld1d(fieldR), ' th'
           fieldR(:) = bufRglob(index_c2l_Sa_shum,:)
           write(6,100) 'lnd','recv', index_c2l_Sa_shum   , global_sum_fld1d(fieldR), ' q'
           fieldR(:) = bufRglob(index_c2l_Sa_pbot,:)
           write(6,100) 'lnd','recv', index_c2l_Sa_pbot   , global_sum_fld1d(fieldR), ' pbot'
           fieldR(:) = bufRglob(index_c2l_Sa_tbot,:)
           write(6,100) 'lnd','recv', index_c2l_Sa_tbot   , global_sum_fld1d(fieldR), ' t'
           fieldR(:) = bufRglob(index_c2l_Faxa_lwdn,:)
           write(6,100) 'lnd','recv', index_c2l_Faxa_lwdn , global_sum_fld1d(fieldR), ' lwrad'
           fieldR(:) = bufRglob(index_c2l_Faxa_rainc,:)
           write(6,100) 'lnd','recv', index_c2l_Faxa_rainc, global_sum_fld1d(fieldR), ' rainc'
           fieldR(:) = bufRglob(index_c2l_Faxa_rainl,:)
           write(6,100) 'lnd','recv', index_c2l_Faxa_rainl, global_sum_fld1d(fieldR), ' rainl'
           fieldR(:) = bufRglob(index_c2l_Faxa_snowc,:)
           write(6,100) 'lnd','recv', index_c2l_Faxa_snowc, global_sum_fld1d(fieldR), ' snowc'
           fieldR(:) = bufRglob(index_c2l_Faxa_snowl,:)
           write(6,100) 'lnd','recv', index_c2l_Faxa_snowl, global_sum_fld1d(fieldR), ' snowl'
           fieldR(:) = bufRglob(index_c2l_Faxa_swndr,:)
           write(6,100) 'lnd','recv', index_c2l_Faxa_swndr, global_sum_fld1d(fieldR), ' soll '
           fieldR(:) = bufRglob(index_c2l_Faxa_swvdr,:)
           write(6,100) 'lnd','recv', index_c2l_Faxa_swvdr, global_sum_fld1d(fieldR), ' sols '
           fieldR(:) = bufRglob(index_c2l_Faxa_swndf,:)
           write(6,100) 'lnd','recv', index_c2l_Faxa_swndf, global_sum_fld1d(fieldR), ' solld'
           fieldR(:) = bufRglob(index_c2l_Faxa_swvdf,:)
           write(6,100) 'lnd','recv', index_c2l_Faxa_swvdf, global_sum_fld1d(fieldR), ' solsd'
100        format('comm_diag ',a3,1x,a4,1x,i3,es26.19,a)
           write(6,*)
        endif
     endif

     ! Stop timer

   ! call t_stopf('lnd_recv')

     ! Check if end of run now, if so stop (each processor does this)

     csmstop_now = .false.
     if (ibuffr(cpl_fields_ibuf_stopnow) /= 0) then
        csmstop_now = .true.
      ! if (timer_lnd_recvsend) call t_stopf('lnd_recvsend')
      ! if (timer_lnd_sendrecv) call t_stopf('lnd_sendrecv')
        write(6,*)'(CSM_RECV) stop now signal from flux coupler'
        write(6,*)'(CSM_RECV) ibuffr(cpl_fields_ibuf_stopnow) = ',ibuffr(cpl_fields_ibuf_stopnow)
        if (masterproc) then
           write(6,9001)
           write(6,9002) ibuffr(cpl_fields_ibuf_cdate)
           write(6,9003)
9001       format(/////' ===========> Terminating CLM Model')
9002       format(     '      Date: ',i8)
9003       format(/////' <=========== CLM Model Terminated')
        endif
        RETURN
     endif

     ! More timer logic

   ! if (.not. timer_lnd_recvsend) then
   !    call t_startf('lnd_recvsend') ; timer_lnd_recvsend = .true.
   ! endif

    ! Split data from coupler into component arrays.
    ! Note that the precipitation fluxes received  from the coupler
    ! are in units of kg/s/m^2. To convert these precipitation rates
    ! in units of mm/sec, one must divide by 1000 kg/m^3 and multiply
    ! by 1000 mm/m resulting in an overall factor of unity.
    ! Below the units are therefore given in mm/s.

!   do g = begg,endg
!      gptr%a2ls%forc_hgt(g)     = bufR(g,index_c2l_Sa_z)        !zgcmxy  Atm state m
!      gptr%a2ls%forc_u(g)       = bufR(g,index_c2l_Sa_u)        !forc_uxy  Atm state m/s
!      gptr%a2ls%forc_v(g)       = bufR(g,index_c2l_Sa_v)        !forc_vxy  Atm state m/s
!      gptr%a2ls%forc_th(g)      = bufR(g,index_c2l_Sa_ptem)     !forc_thxy Atm state K
!      gptr%a2ls%forc_q(g)       = bufR(g,index_c2l_Sa_shum)     !forc_qxy  Atm state kg/kg
!      gptr%a2ls%forc_pbot(g)    = bufR(g,index_c2l_Sa_pbot)     !ptcmxy  Atm state Pa
!      gptr%a2ls%forc_t(g)       = bufR(g,index_c2l_Sa_tbot)     !forc_txy  Atm state K
!      gptr%a2lf%forc_lwrad(g)   = bufR(g,index_c2l_Faxa_lwdn)   !flwdsxy Atm flux  W/m^2
!      forc_rainc                = bufR(g,index_c2l_Faxa_rainc)  !mm/s
!      forc_rainl                = bufR(g,index_c2l_Faxa_rainl)  !mm/s
!      forc_snowc                = bufR(g,index_c2l_Faxa_snowc)  !mm/s
!      forc_snowl                = bufR(g,index_c2l_Faxa_snowl)  !mm/s
!      gptr%a2lf%forc_solad(g,2) = bufR(g,index_c2l_Faxa_swndr)  !forc_sollxy  Atm flux  W/m^2
!      gptr%a2lf%forc_solad(g,1) = bufR(g,index_c2l_Faxa_swvdr)  !forc_solsxy  Atm flux  W/m^2
!      gptr%a2lf%forc_solai(g,2) = bufR(g,index_c2l_Faxa_swndf)  !forc_solldxy Atm flux  W/m^2
!      gptr%a2lf%forc_solai(g,1) = bufR(g,index_c2l_Faxa_swvdf)  !forc_solsdxy Atm flux  W/m^2

!      ! Determine derived quantities

!      gptr%a2ls%forc_hgt_u(g) = gptr%a2ls%forc_hgt(g)    !observational height of wind [m]
!      gptr%a2ls%forc_hgt_t(g) = gptr%a2ls%forc_hgt(g)    !observational height of temperature [m]
!      gptr%a2ls%forc_hgt_q(g) = gptr%a2ls%forc_hgt(g)    !observational height of humidity [m]
!      gptr%a2ls%forc_vp(g)    = gptr%a2ls%forc_q(g) * gptr%a2ls%forc_pbot(g) &
!                                / (0.622 + 0.378 * gptr%a2ls%forc_q(g))
!      gptr%a2ls%forc_rho(g)   = (gptr%a2ls%forc_pbot(g) - 0.378 * gptr%a2ls%forc_vp(g)) &
!                                / (rair * gptr%a2ls%forc_t(g))
!      gptr%a2ls%forc_co2(g)   = pco2 * gptr%a2ls%forc_pbot(g)
!      gptr%a2ls%forc_o2(g)    = po2 * gptr%a2ls%forc_pbot(g)
!      gptr%a2ls%forc_wind(g)  = sqrt(gptr%a2ls%forc_u(g)**2 + gptr%a2ls%forc_v(g)**2)
!      gptr%a2lf%forc_solar(g) = gptr%a2lf%forc_solad(g,1) + gptr%a2lf%forc_solai(g,1) + &
!                                gptr%a2lf%forc_solad(g,2) + gptr%a2lf%forc_solai(g,2)

!      ! Determine precipitation needed by clm

!      gptr%a2lf%forc_rain(g) = forc_rainc + forc_rainl
!      gptr%a2lf%forc_snow(g) = forc_snowc + forc_snowl

!   end do

    if (co2_option == 'co2prog' .and. index_c2l_Sa_co2prog == 0) then
       write(6,*)' must have nonzero index_c2l_Sa_co2prog for co2_option equal to co2prog'
       call endrun()
    else if (co2_option == 'co2diag' .and. index_c2l_Sa_co2diag == 0) then
       write(6,*)' must have nonzero index_c2l_Sa_co2diag for co2_option equal to co2diag'
       call endrun()
    end if

    do g = 1,numgrid

!      pco2m     => forc( 1,i)       !CO2 concentration in atmos. (35 pa)
!      po2m      => forc( 2,i)       !O2 concentration in atmos. (20900 pa)
!      us        => forc( 3,i)       !wind in eastward direction [m/s]
!      vs        => forc( 4,i)       !wind in northward direction [m/s]
!      tm        => forc( 5,i)       !temperature at reference height [kelvin]
!      qm        => forc( 6,i)       !specific humidity at reference height [kg/kg]
!      prc       => forc( 7,i)       !convective precipitation [mm/s]
!      prl       => forc( 8,i)       !large scale precipitation [mm/s]
!      pbot      => forc( 9,i)       !atm bottom level pressure (or reference height) (pa)
!      psrf      => forc(10,i)       !atmospheric pressure at the surface [pa]
!      sols      => forc(11,i)       !atm vis direct beam solar rad onto srf [W/m2]
!      soll      => forc(12,i)       !atm nir direct beam solar rad onto srf [W/m2]
!      solsd     => forc(13,i)       !atm vis diffuse solar rad onto srf [W/m2]
!      solld     => forc(14,i)       !atm nir diffuse solar rad onto srf [W/m2]
!      frl       => forc(15,i)       !atmospheric infrared (longwave) radiation [W/m2]
!      hu        => forc(16,i)       !observational height of wind [m]
!      ht        => forc(17,i)       !observational height of temperature [m]
!      hq        => forc(18,i)       !observational height of humidity [m]

!      rhoair    =  (pbot-0.378*qm*pbot/(0.622+0.378*qm))/(rgas*tm)

       ! Determine optional receive fields

       if (index_c2l_Sa_co2prog /= 0) then
          co2_ppmv_prog = bufR(g,index_c2l_Sa_co2prog)   ! co2 atm state prognostic
       else
          co2_ppmv_prog = pco2*1.e6
       end if

       if (index_c2l_Sa_co2diag /= 0) then
          co2_ppmv_diag = bufR(g,index_c2l_Sa_co2diag)   ! co2 atm state diagnostic
       else
          co2_ppmv_diag = pco2*1.e6
       end if

       ! Determine derived quantities for optional fields
       ! Note that the following does unit conversions from ppmv to partial pressures (Pa)
       ! Note that forc_pbot is in Pa

       if (co2_option == 'co2prog') then
          co2_ppmv_val = co2_ppmv_prog
       else if (co2_option == 'co2diag') then
          co2_ppmv_val = co2_ppmv_diag
       else
          co2_ppmv_val = pco2*1.e6
       end if

       forcg( 1,g)           = bufR(g,index_c2l_Sa_pbot)*co2_ppmv_val*1.e-6

       forcg( 2,g)           = bufR(g,index_c2l_Sa_pbot)*po2
       forcg( 3,g)           = bufR(g,index_c2l_Sa_u   )  !forc_uxy  Atm state m/s
       forcg( 4,g)           = bufR(g,index_c2l_Sa_v   )  !forc_vxy  Atm state m/s
       forcg( 5,g)           = bufR(g,index_c2l_Sa_tbot)  !forc_txy  Atm state K
       forcg( 6,g)           = bufR(g,index_c2l_Sa_shum)  !forc_qxy  Atm state kg/kg

       forc_rainc            = bufR(g,index_c2l_Faxa_rainc)  !mm/s
       forc_rainl            = bufR(g,index_c2l_Faxa_rainl)  !mm/s
       forc_snowc            = bufR(g,index_c2l_Faxa_snowc)  !mm/s
       forc_snowl            = bufR(g,index_c2l_Faxa_snowl)  !mm/s

       forcg( 7,g)           = forc_rainc + forc_snowc
       forcg( 8,g)           = forc_rainl + forc_snowl

       forcg( 9,g)           = bufR(g,index_c2l_Sa_pbot)  !ptcmxy  Atm state Pa
       forcg(10,g)           = bufR(g,index_c2l_Sa_pbot)  !ptcmxy  Atm state Pa

       forcg(11,g)           = bufR(g,index_c2l_Faxa_swvdr)  !forc_solsxy  Atm flux  W/m^2
       forcg(12,g)           = bufR(g,index_c2l_Faxa_swndr)  !forc_sollxy  Atm flux  W/m^2
       forcg(13,g)           = bufR(g,index_c2l_Faxa_swvdf)  !forc_solsdxy Atm flux  W/m^2
       forcg(14,g)           = bufR(g,index_c2l_Faxa_swndf)  !forc_solldxy Atm flux  W/m^2

       forcg(15,g)           = bufR(g,index_c2l_Faxa_lwdn)  !flwdsxy Atm flux  W/m^2
       forcg(16,g)           = bufR(g,index_c2l_Sa_z     )  !zgcmxy  Atm state m
       forcg(17,g)           = bufR(g,index_c2l_Sa_z     )  !zgcmxy  Atm state m
       forcg(18,g)           = bufR(g,index_c2l_Sa_z     )  !zgcmxy  Atm state m

    end do

#ifdef MYBUG
    write(6,*), 'index_c2l_Sa_pbot'   , minval(bufR(:,index_c2l_Sa_pbot))   , maxval(bufR(:,index_c2l_Sa_pbot))
    write(6,*), 'index_c2l_Sa_u'      , minval(bufR(:,index_c2l_Sa_u))      , maxval(bufR(:,index_c2l_Sa_u))
    write(6,*), 'index_c2l_Sa_v'      , minval(bufR(:,index_c2l_Sa_v))      , maxval(bufR(:,index_c2l_Sa_v))
    write(6,*), 'index_c2l_Sa_tbot'   , minval(bufR(:,index_c2l_Sa_tbot))   , maxval(bufR(:,index_c2l_Sa_tbot))
    write(6,*), 'index_c2l_Sa_shum'   , minval(bufR(:,index_c2l_Sa_shum))   , maxval(bufR(:,index_c2l_Sa_shum))
    write(6,*), 'index_c2l_Faxa_rainc', minval(bufR(:,index_c2l_Faxa_rainc)), maxval(bufR(:,index_c2l_Faxa_rainc))
    write(6,*), 'index_c2l_Faxa_rainl', minval(bufR(:,index_c2l_Faxa_rainl)), maxval(bufR(:,index_c2l_Faxa_rainl))
    write(6,*), 'index_c2l_Faxa_snowc', minval(bufR(:,index_c2l_Faxa_snowc)), maxval(bufR(:,index_c2l_Faxa_snowc))
    write(6,*), 'index_c2l_Faxa_snowl', minval(bufR(:,index_c2l_Faxa_snowl)), maxval(bufR(:,index_c2l_Faxa_snowl))
    write(6,*), 'index_c2l_Faxa_swvdr', minval(bufR(:,index_c2l_Faxa_swvdr)), maxval(bufR(:,index_c2l_Faxa_swvdr))
    write(6,*), 'index_c2l_Faxa_swndr', minval(bufR(:,index_c2l_Faxa_swndr)), maxval(bufR(:,index_c2l_Faxa_swndr))
    write(6,*), 'index_c2l_Faxa_swvdf', minval(bufR(:,index_c2l_Faxa_swvdf)), maxval(bufR(:,index_c2l_Faxa_swvdf))
    write(6,*), 'index_c2l_Faxa_swndf', minval(bufR(:,index_c2l_Faxa_swndf)), maxval(bufR(:,index_c2l_Faxa_swndf))
    write(6,*), 'index_c2l_Faxa_lwdn' , minval(bufR(:,index_c2l_Faxa_lwdn)) , maxval(bufR(:,index_c2l_Faxa_lwdn))
    write(6,*), 'index_c2l_Sa_z'      , minval(bufR(:,index_c2l_Sa_z))      , maxval(bufR(:,index_c2l_Sa_z))

    if (index_c2l_Sa_co2prog /= 0) then
       write(6,*) 'colm_csmMod: co2_ppmv_prog', minval(bufR(:,index_c2l_Sa_co2prog)), maxval(bufR(:,index_c2l_Sa_co2prog))
    end if
    if (index_c2l_Sa_co2diag /= 0) then
       write(6,*) 'colm_csmMod: co2_ppmv_diag', minval(bufR(:,index_c2l_Sa_co2diag)), maxval(bufR(:,index_c2l_Sa_co2diag))
    end if

    write(6,*) 'colm_csmMod: co2_ppmv_val', co2_option=='co2prog', co2_option=='co2diag', minval(forcg(1,:)), maxval(forcg(1,:))
#endif

  end subroutine csm_recv

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: csm_send
!
! !INTERFACE:
  subroutine csm_send()
!
! !DESCRIPTION:
! Send data to the flux coupler
!
! !USES:
  ! use clm_varctl  , only : csm_doflxave
  ! use clm_varsur  , only : landmask
  ! use clm_varcon  , only : sb
  ! use time_manager, only : get_curr_date, get_nstep
  ! use lnd2atmMod  , only : lnd2atm

    use colm_cplMod
    use colm_varMod  , only : landmask, numgrid, numgrid_glob
    use timemgr      , only : get_curr_date, get_nstep
    use phycon_module, only : sb
    use spmdGathScatMod, only : gather_data_to_master
!
! !ARGUMENTS:
    implicit none
!
! !REVISION HISTORY:
!  02.09.17 Mariana Vertenstein Updated to clm2_1 data structures
!
!EOP
!
! !LOCAL VARIABLES:
    integer :: i,j,n,ni,g,c,p  ! indices
    integer :: yr              ! current year
    integer :: mon             ! current month
    integer :: day             ! current day (0, 1, ...)
    integer :: ncsec           ! current seconds of current date (0, ..., 86400)
    integer :: ncdate          ! current date (yymmdd format) (e.g., 021105)
    integer :: ier             ! error status
    real(r8):: wt              ! weight
 !  type(gridcell_type), pointer :: gptr  ! pointer to gridcell derived type
 !  type(pft_type)     , pointer :: pptr  ! pointer to column   derived type
!-----------------------------------------------------------------------

    ! Set pointers into derived type

  ! gptr => clm3%g
  ! pptr => clm3%g%l%c%p

    ! Send data to the flux coupler

  ! if (timer_lnd_recvsend) then
  !    call t_stopf ('lnd_recvsend') ; timer_lnd_recvsend = .false.
  ! endif

    ! Start timer

  ! call t_startf('lnd_send')

    ! Determine 1d vector of fields that will be sent to coupler.
    ! Coupler has convention that fluxes are positive downward.

    bufS(:,:) = 0.0

    ! Calculate fluxes and states to send to atm

  ! call lnd2atm()

    ! Recalculate fluxes to send to atm in flux averaged case
    ! Also recalculate a new t_rad based on flux averged lwup

!   if (csm_doflxave) then
!      do g = begg,endg
!         gptr%l2af%taux(g) = 0.
!         gptr%l2af%tauy(g) = 0.
!         gptr%l2af%eflx_lh_tot(g) = 0.
!         gptr%l2af%eflx_sh_tot(g) = 0.
!         gptr%l2af%eflx_lwrad_out(g) = 0.
!         gptr%l2af%qflx_evap_tot(g) = 0.
!         gptr%l2af%fsa(g) = 0.
!      end do
!      do p = begp,endp
!         g  = pptr%gridcell(p)
!         wt = pptr%wtgcell(p)
!         gptr%l2af%taux(g)           = gptr%l2af%taux(g)           + taux_ave(p)  * wt
!         gptr%l2af%tauy(g)           = gptr%l2af%tauy(g)           + tauy_ave(p)  * wt
!         gptr%l2af%eflx_lh_tot(g)    = gptr%l2af%eflx_lh_tot(g)    + lhflx_ave(p) * wt
!         gptr%l2af%eflx_sh_tot(g)    = gptr%l2af%eflx_sh_tot(g)    + shflx_ave(p) * wt
!         gptr%l2af%eflx_lwrad_out(g) = gptr%l2af%eflx_lwrad_out(g) + lwup_ave(p)  * wt
!         gptr%l2af%qflx_evap_tot(g)  = gptr%l2af%qflx_evap_tot(g)  + qflx_ave(p)  * wt
!         gptr%l2af%fsa(g)            = gptr%l2af%fsa(g)            + swabs_ave(p) * wt
!      end do
!      do g = begg,endg
!         gptr%l2as%t_rad(g) = (abs(gptr%l2af%eflx_lwrad_out(g)/sb))**0.25
!      end do
!   endif

!   do g = begg,endg
!      bufS(g,index_l2c_Sl_t      ) =  gptr%l2as%t_rad(g)
!      bufS(g,index_l2c_Sl_snowh  ) =  gptr%l2as%h2osno(g)
!      bufS(g,index_l2c_Sl_avsdr  ) =  gptr%l2as%albd(g,1)
!      bufS(g,index_l2c_Sl_anidr  ) =  gptr%l2as%albd(g,2)
!      bufS(g,index_l2c_Sl_avsdf  ) =  gptr%l2as%albi(g,1)
!      bufS(g,index_l2c_Sl_anidf  ) =  gptr%l2as%albi(g,2)
!      bufS(g,index_l2c_Sl_tref   ) =  gptr%l2as%t_ref2m(g)
!      bufS(g,index_l2c_Sl_qref   ) =  gptr%l2as%q_ref2m(g)
!      bufS(g,index_l2c_Fall_taux ) = -gptr%l2af%taux(g)
!      bufS(g,index_l2c_Fall_tauy ) = -gptr%l2af%tauy(g)
!      bufS(g,index_l2c_Fall_lat  ) = -gptr%l2af%eflx_lh_tot(g)
!      bufS(g,index_l2c_Fall_sen  ) = -gptr%l2af%eflx_sh_tot(g)
!      bufS(g,index_l2c_Fall_lwup ) = -gptr%l2af%eflx_lwrad_out(g)
!      bufS(g,index_l2c_Fall_evap ) = -gptr%l2af%qflx_evap_tot(g)
!      bufS(g,index_l2c_Fall_swnet) = -gptr%l2af%fsa(g)
!   end do

    if (csm_doflxave) then

       do g = 1,numgrid
          lnd_trad(g) = abs(lwup_ave(g)/sb)**0.25
       end do

       do g = 1,numgrid
          bufS(g,index_l2c_Sl_t       ) =  lnd_trad(g)
          bufS(g,index_l2c_Sl_snowh   ) =  lnd_scv(g)
          bufS(g,index_l2c_Sl_avsdr   ) =  lnd_avsdr(g) ! avsdr
          bufS(g,index_l2c_Sl_avsdf   ) =  lnd_avsdf(g) ! avsdf
          bufS(g,index_l2c_Sl_anidr   ) =  lnd_anidr(g) ! anidr
          bufS(g,index_l2c_Sl_anidf   ) =  lnd_anidf(g) ! anidf
          bufS(g,index_l2c_Sl_tref    ) =  lnd_tref(g)
          bufS(g,index_l2c_Sl_qref    ) =  lnd_qref(g)
          bufS(g,index_l2c_Sl_u10m    ) =  lnd_u10m(g)
          bufS(g,index_l2c_Sl_v10m    ) =  lnd_v10m(g)
          bufS(g,index_l2c_Fall_sublim) =  lnd_sublim(g)
          bufS(g,index_l2c_Fall_taux  ) = -taux_ave(g)
          bufS(g,index_l2c_Fall_tauy  ) = -tauy_ave(g)
          bufS(g,index_l2c_Fall_lat   ) = -lhflx_ave(g)
          bufS(g,index_l2c_Fall_sen   ) = -shflx_ave(g)
          bufS(g,index_l2c_Fall_lwup  ) = -lwup_ave(g)
          bufS(g,index_l2c_Fall_evap  ) = -qflx_ave(g)
          bufS(g,index_l2c_Fall_swnet ) = -swabs_ave(g)
          if (index_l2c_Fall_nee /= 0) then
             bufS(g,index_l2c_Fall_nee) =  nee_ave(g)
          end if
       end do

    else

       do g = 1,numgrid
          lnd_trad(g) = abs(lnd_lwup(g)/sb)**0.25
       end do

       do g = 1,numgrid
          bufS(g,index_l2c_Sl_t       ) =  lnd_trad(g)
          bufS(g,index_l2c_Sl_snowh   ) =  lnd_scv(g)
          bufS(g,index_l2c_Sl_avsdr   ) =  lnd_avsdr(g) ! avsdr
          bufS(g,index_l2c_Sl_avsdf   ) =  lnd_avsdf(g) ! avsdf
          bufS(g,index_l2c_Sl_anidr   ) =  lnd_anidr(g) ! anidr
          bufS(g,index_l2c_Sl_anidf   ) =  lnd_anidf(g) ! anidf
          bufS(g,index_l2c_Sl_tref    ) =  lnd_tref(g)
          bufS(g,index_l2c_Sl_qref    ) =  lnd_qref(g)
          bufS(g,index_l2c_Sl_u10m    ) =  lnd_u10m(g)
          bufS(g,index_l2c_Sl_v10m    ) =  lnd_v10m(g)
          bufS(g,index_l2c_Fall_sublim) =  lnd_sublim(g)
          bufS(g,index_l2c_Fall_taux  ) = -lnd_taux(g)
          bufS(g,index_l2c_Fall_tauy  ) = -lnd_tauy(g)
          bufS(g,index_l2c_Fall_lat   ) = -lnd_lhflx(g)
          bufS(g,index_l2c_Fall_sen   ) = -lnd_shflx(g)
          bufS(g,index_l2c_Fall_lwup  ) = -lnd_lwup(g)
          bufS(g,index_l2c_Fall_evap  ) = -lnd_qflx(g)
          bufS(g,index_l2c_Fall_swnet ) = -lnd_swabs(g)
          if (index_l2c_Fall_nee /= 0) then
             bufS(g,index_l2c_Fall_nee) =  lnd_nee(g)
          end if
       end do

    end if

   !DEBUG
   !write(6,*), 'trad', minval(bufS(:,index_l2c_Sl_t)), maxval(bufS(:,index_l2c_Sl_t))
   !write(6,*), 'scv', minval(bufS(:,index_l2c_Sl_snowh)), maxval(bufS(:,index_l2c_Sl_snowh))
   !write(6,*), 'avsdr', minval(bufS(:,index_l2c_Sl_avsdr)), maxval(bufS(:,index_l2c_Sl_avsdr))
   !write(6,*), 'avsdf', minval(bufS(:,index_l2c_Sl_avsdf)), maxval(bufS(:,index_l2c_Sl_avsdf))
   !write(6,*), 'anidr', minval(bufS(:,index_l2c_Sl_anidr)), maxval(bufS(:,index_l2c_Sl_anidr))
   !write(6,*), 'anidf', minval(bufS(:,index_l2c_Sl_anidf)), maxval(bufS(:,index_l2c_Sl_anidf))
   !write(6,*), 'tref', minval(bufS(:,index_l2c_Sl_tref)), maxval(bufS(:,index_l2c_Sl_tref))
   !write(6,*), 'qref', minval(bufS(:,index_l2c_Sl_qref)), maxval(bufS(:,index_l2c_Sl_qref))
   !write(6,*), 'taux', minval(bufS(:,index_l2c_Fall_taux)), maxval(bufS(:,index_l2c_Fall_taux))
   !write(6,*), 'tauy', minval(bufS(:,index_l2c_Fall_tauy)), maxval(bufS(:,index_l2c_Fall_tauy))
   !write(6,*), 'lat', minval(bufS(:,index_l2c_Fall_lat)), maxval(bufS(:,index_l2c_Fall_lat))
   !write(6,*), 'sen', minval(bufS(:,index_l2c_Fall_sen)), maxval(bufS(:,index_l2c_Fall_sen))
   !write(6,*), 'lwup', minval(bufS(:,index_l2c_Fall_lwup)), maxval(bufS(:,index_l2c_Fall_lwup))
   !write(6,*), 'evap', minval(bufS(:,index_l2c_Fall_evap)), maxval(bufS(:,index_l2c_Fall_evap))
   !write(6,*), 'swnet', minval(bufS(:,index_l2c_Fall_swnet)), maxval(bufS(:,index_l2c_Fall_swnet))
   !write(6,*), 'nee', minval(bufS(:,index_l2c_Fall_nee)), maxval(bufS(:,index_l2c_Fall_nee))
   !DEBUG

    ! Map fields from 1d-grid vector to 2d-grid
    ! NOTE: snow is sent as zero over non-land because currently
    ! the ocn and sea-ice send no snow cover to cpl and so the cpl
    ! sends back zero snow over non-land to  the atm (the atm and
    ! land grid are currently assumed to be identical)

    call get_curr_date (yr, mon, day, ncsec)
    ncdate = yr*10000 + mon*100 + day

    ibuffs(:)  = 0
    ibuffs(cpl_fields_ibuf_cdate)  = ncdate      ! model date (yyyymmdd)
    ibuffs(cpl_fields_ibuf_sec  )  = ncsec       ! elapsed seconds in current date

    call cpl_interface_contractSend(cpl_fields_cplname,contractS ,ibuffs,bufS)

    ! Must convert runoff to units of kg/m^2/s from m^3/s

    bufSr(:,:) = 0.
    do n = beg_ocn_rof, end_ocn_rof
       ni = n - beg_ocn_rof+ 1
       bufSr(ni,index_r2c_Forr_roff) = runoff%ocn(n) / (runoff%ocn_area(n)*1000.)

       if(bufSr(ni,index_r2c_Forr_roff).lt.0) bufSr(ni,index_r2c_Forr_roff) = 0.
    end do

#ifdef MYBUG
    write(6,*) 'runoff%ocn', minval(bufSr(:,index_r2c_Forr_roff)), maxval(bufSr(:,index_r2c_Forr_roff))
#endif

    call cpl_interface_contractSend(cpl_fields_cplname,contractSr,ibuffs,bufSr)

    ! Do global integrals if flag is set

    if (ibuffr(cpl_fields_ibuf_infobug) >= 2) then
       if (masterproc) then
          if (.not. associated(bufSglob)) then
           ! allocate(bufSglob(nsend,numg), stat=ier)
             allocate(bufSglob(nsend,numgrid_glob), stat=ier)
             if (ier /= 0) then
                write(6,*)'clm_csmMod: allocation error for bufSglob'; call endrun
             end if
          end if
          if (.not. associated(fieldS)) then
           ! allocate(fieldS(numg), stat=ier)
             allocate(fieldS(numgrid_glob), stat=ier)
             if (ier /= 0) then
                write(6,*)'clm_csmMod: allocation error for fieldS'; call endrun
             end if
          endif
       end if
       if (.not. associated(bufSloc)) then
        ! allocate(bufSloc(nsend,begg:endg), stat=ier)
          allocate(bufSloc(nsend,numgrid), stat=ier)
          if (ier /= 0) then
             write(6,*)'clm_csmMod: allocation error for bufSloc'; call endrun
          end if
       end if
     ! do g = begg,endg
       do g = 1,numgrid
          do n = 1,nsend
             bufSloc(n,g) = bufS(g,n)
          end do
       end do
#if (defined SPMD)
       call gather_data_to_master(bufSloc, bufSglob, clmlevel='gridcell')
#else
       bufSglob(:,:) = bufSloc(:,:)
#endif
       if (masterproc) then
          write(6,*)
          fieldS(:) = bufSglob(index_l2c_Sl_t,:)
          write(6,100) 'lnd','send', index_l2c_Sl_t      , global_sum_fld1d(fieldS), ' trad'
          fieldS(:) = bufSglob(index_l2c_Sl_avsdr,:)
          write(6,100) 'lnd','send', index_l2c_Sl_avsdr  , global_sum_fld1d(fieldS),' asdir'
          fieldS(:) = bufSglob(index_l2c_Sl_anidr,:)
          write(6,100) 'lnd','send', index_l2c_Sl_anidr  , global_sum_fld1d(fieldS),' aldir'
          fieldS(:) = bufSglob(index_l2c_Sl_avsdf,:)
          write(6,100) 'lnd','send', index_l2c_Sl_avsdf  , global_sum_fld1d(fieldS),' asdif'
          fieldS(:) = bufSglob(index_l2c_Sl_anidf,:)
          write(6,100) 'lnd','send', index_l2c_Sl_anidf  , global_sum_fld1d(fieldS),' aldif'
          fieldS(:) = bufSglob(index_l2c_Fall_taux,:)
          write(6,100) 'lnd','send', index_l2c_Fall_taux , global_sum_fld1d(fieldS), ' taux'
          fieldS(:) = bufSglob(index_l2c_Fall_tauy,:)
          write(6,100) 'lnd','send', index_l2c_Fall_tauy , global_sum_fld1d(fieldS), ' tauy'
          fieldS(:) = bufSglob(index_l2c_Fall_lat,:)
          write(6,100) 'lnd','send', index_l2c_Fall_lat  , global_sum_fld1d(fieldS), ' lhflx'
          fieldS(:) = bufSglob(index_l2c_Fall_sen,:)
          write(6,100) 'lnd','send', index_l2c_Fall_sen  , global_sum_fld1d(fieldS), ' shflx'
          fieldS(:) = bufSglob(index_l2c_Fall_lwup,:)
          write(6,100) 'lnd','send', index_l2c_Fall_lwup , global_sum_fld1d(fieldS), ' lwup'
          fieldS(:) = bufSglob(index_l2c_Fall_evap,:)
          write(6,100) 'lnd','send', index_l2c_Fall_evap , global_sum_fld1d(fieldS), ' qflx'
          fieldS(:) = bufSglob(index_l2c_Fall_swnet,:)
          write(6,100) 'lnd','send', index_l2c_Fall_swnet, global_sum_fld1d(fieldS), ' swabs'

          if (index_l2c_Fall_nee /= 0) then
             fieldS(:) = bufSglob(index_l2c_Fall_nee,:)
             write(6,100) 'lnd','send', index_l2c_Fall_nee  , global_sum_fld1d(fieldS), ' nee'
          end if

          write(6,*)
100       format('comm_diag ',a3,1x,a4,1x,i3,es26.19,a)
       endif
    endif

    ! Stop timers

  ! call t_stopf('lnd_send')

  ! if (.not. timer_lnd_recvsend) then
  !    call t_startf('lnd_sendrecv') ; timer_lnd_sendrecv = .true.
  ! endif

  end subroutine csm_send

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: csm_flxave
!
! !INTERFACE:
  subroutine csm_flxave()
!
! !DESCRIPTION:
! Average output fluxes for flux coupler
! Add land surface model output fluxes to accumulators every time step.
! When icnt==ncnt, compute the average flux over the time interval.
!
! !USES:
!
  ! use clmtype
  ! use clm_varctl  , only : irad
  ! use time_manager, only : get_nstep
    use timemgr     , only : get_nstep
    use colm_varMod , only : numgrid
    use colm_cplMod
!
! !ARGUMENTS:
    implicit none
!
! !REVISION HISTORY:
!  02.09.17 Mariana Vertenstein Updated to clm2_1 data structures
!
!EOP
!
! !LOCAL VARIABLES:
    integer :: p            ! indices
    integer :: nstep        ! model time step
    integer :: g
  ! type(pft_type), pointer :: pptr  ! pointer to pft derived subtype
!-----------------------------------------------------------------------

    ! Set pointers into derived type

  ! pptr => clm3%g%l%c%p

    ! Allocate dynamic memory if necessary
    if (.not. allocated(taux_ave)) then
       allocate (taux_ave(numgrid)) ; taux_ave(:) = nan
    endif
    if (.not. allocated(tauy_ave)) then
       allocate (tauy_ave(numgrid)) ; tauy_ave(:) = nan
    endif
    if (.not. allocated(lhflx_ave)) then
       allocate (lhflx_ave(numgrid)); lhflx_ave(:) = nan
    endif
    if (.not. allocated(shflx_ave)) then
       allocate (shflx_ave(numgrid)); shflx_ave(:) = nan
    endif
    if (.not. allocated(lwup_ave)) then
       allocate (lwup_ave(numgrid)) ; lwup_ave(:) = nan
    endif
    if (.not. allocated(qflx_ave)) then
       allocate (qflx_ave(numgrid)) ; qflx_ave(:) = nan
    endif
    if (.not. allocated(swabs_ave)) then
       allocate (swabs_ave(numgrid)) ; swabs_ave(:) = nan
    endif
    if (index_l2c_Fall_nee /= 0) then
       if (.not. allocated(nee_ave)) then
          allocate (nee_ave(numgrid)) ; nee_ave(:) = nan
       endif
    end if

    ! Determine output flux averaging interval

    nstep = get_nstep()
    if (dorecv) then
       icnt = 1
       if ( nstep==0 ) then
          ncnt = irad + 1
       else
          ncnt = irad
       endif
       rncnt = 1./ncnt
    endif

    if (icnt == 1) then

       ! Initial call of averaging interval, copy data to accumulators

!      do p = begp,endp
!         taux_ave(p)  = pptr%pmf%taux(p)
!         tauy_ave(p)  = pptr%pmf%tauy(p)
!         lhflx_ave(p) = pptr%pef%eflx_lh_tot(p)
!         shflx_ave(p) = pptr%pef%eflx_sh_tot(p)
!         lwup_ave(p)  = pptr%pef%eflx_lwrad_out(p)
!         qflx_ave(p)  = pptr%pwf%qflx_evap_tot(p)
!         swabs_ave(p) = pptr%pef%fsa(p)
!      end do

       do g = 1, numgrid
          taux_ave(g)  = lnd_taux(g)
          tauy_ave(g)  = lnd_tauy(g)
          lhflx_ave(g) = lnd_lhflx(g)
          shflx_ave(g) = lnd_shflx(g)
          lwup_ave(g)  = lnd_lwup(g)
          qflx_ave(g)  = lnd_qflx(g)
          swabs_ave(g) = lnd_swabs(g)
          if (index_l2c_Fall_nee /= 0) then
             nee_ave(g)= lnd_nee(g)
          end if
       end do

    else if (icnt == ncnt) then

       ! Final call of averaging interval, complete averaging

!      do p = begp,endp
!         taux_ave(p)  = rncnt * (taux_ave(p)  + pptr%pmf%taux(p))
!         tauy_ave(p)  = rncnt * (tauy_ave(p)  + pptr%pmf%tauy(p))
!         lhflx_ave(p) = rncnt * (lhflx_ave(p) + pptr%pef%eflx_lh_tot(p))
!         shflx_ave(p) = rncnt * (shflx_ave(p) + pptr%pef%eflx_sh_tot(p))
!         lwup_ave(p)  = rncnt * (lwup_ave(p)  + pptr%pef%eflx_lwrad_out(p))
!         qflx_ave(p)  = rncnt * (qflx_ave(p)  + pptr%pwf%qflx_evap_tot(p))
!         swabs_ave(p) = rncnt * (swabs_ave(p) + pptr%pef%fsa(p))
!      end do

       do g = 1,numgrid
          taux_ave(g)  = rncnt * (taux_ave(g)  + lnd_taux(g))
          tauy_ave(g)  = rncnt * (tauy_ave(g)  + lnd_tauy(g))
          lhflx_ave(g) = rncnt * (lhflx_ave(g) + lnd_lhflx(g))
          shflx_ave(g) = rncnt * (shflx_ave(g) + lnd_shflx(g))
          lwup_ave(g)  = rncnt * (lwup_ave(g)  + lnd_lwup(g))
          qflx_ave(g)  = rncnt * (qflx_ave(g)  + lnd_qflx(g))
          swabs_ave(g) = rncnt * (swabs_ave(g) + lnd_swabs(g))
          if (index_l2c_Fall_nee /= 0) then
             nee_ave(g) = rncnt * (nee_ave(g) + lnd_nee(g))
          end if
       end do

    else

       ! Intermediate call, add data to accumulators

!      do p = begp,endp
!         taux_ave(p)  = (taux_ave(p)  + pptr%pmf%taux(p))
!         tauy_ave(p)  = (tauy_ave(p)  + pptr%pmf%tauy(p))
!         lhflx_ave(p) = (lhflx_ave(p) + pptr%pef%eflx_lh_tot(p))
!         shflx_ave(p) = (shflx_ave(p) + pptr%pef%eflx_sh_tot(p))
!         lwup_ave(p)  = (lwup_ave(p)  + pptr%pef%eflx_lwrad_out(p))
!         qflx_ave(p)  = (qflx_ave(p)  + pptr%pwf%qflx_evap_tot(p))
!         swabs_ave(p) = (swabs_ave(p) + pptr%pef%fsa(p))
!      end do

       do g = 1,numgrid
          taux_ave(g)  = (taux_ave(g)  + lnd_taux(g))
          tauy_ave(g)  = (tauy_ave(g)  + lnd_tauy(g))
          lhflx_ave(g) = (lhflx_ave(g) + lnd_lhflx(g))
          shflx_ave(g) = (shflx_ave(g) + lnd_shflx(g))
          lwup_ave(g)  = (lwup_ave(g)  + lnd_lwup(g))
          qflx_ave(g)  = (qflx_ave(g)  + lnd_qflx(g))
          swabs_ave(g) = (swabs_ave(g) + lnd_swabs(g))
          if (index_l2c_Fall_nee /= 0) then
             nee_ave(g) = (nee_ave(g)  + lnd_nee(g))
          end if
       end do

    end if

    ! Increment counter

    icnt = icnt + 1

  end subroutine csm_flxave

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: compat_check
!
! !INTERFACE:
  subroutine compat_check_spval(spval, data, string)
!
! !DESCRIPTION:
! Check that the given piece of real data sent from the coupler is valid
! data and not the couplers special data flag.  This ensures that the data
! you expect is actually being sent by the coupler.
!
! !ARGUMENTS:
    implicit none
    real(r8), intent(in) :: spval
    real(r8), intent(in) :: data
    character(len=*), intent(in) :: string
!
! !REVISION HISTORY:
!  02.09.17 Mariana Vertenstein Updated to clm2_1 data structures
!
!EOP
! -----------------------------------------------------------------

    if ( spval == data )then
       write(6,*)'ERROR:(compat_check_spval) msg incompatibility'
       write(6,*)'ERROR: I expect to recieve the data type: ',string
       write(6,*)'from CPL, but all I got was the special data flag'
       write(6,*)'coupler must not be sending this data, you are'
       write(6,*)'running with an incompatable version of the coupler'
       call endrun
    end if

  end subroutine compat_check_spval

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: csm_compat
!
! !INTERFACE:
  subroutine csm_compat(cpl_maj_vers, cpl_min_vers, expect_maj_vers, &
                        expect_min_vers)
!
! !DESCRIPTION:
! Checks that the message recieved from the coupler is compatable
! with the type of message that I expect to recieve.  If the minor
! version numbers differ I print a warning message.  If the major
! numbers differ I abort since that means that the change is
! drastic enough that I can't run with the differences.
! Original Author: Erik Kluzek Dec/97
!
! !PARAMETERS:
    implicit none
    integer, intent(in) :: cpl_maj_vers    ! major version from coupler initial ibuffr array
    integer, intent(in) :: cpl_min_vers    ! minor version from coupler initial ibuffr array
    integer, intent(in) :: expect_maj_vers ! major version of the coupler I'm expecting
    integer, intent(in) :: expect_min_vers ! minor version of the coupler I'm expecting
!
! !REVISION HISTORY:
!  02.09.17 Mariana Vertenstein Updated to clm2_1 data structures
!
!EOP
! -----------------------------------------------------------------

    write(6,*)'(cpl_COMPAT): This is revision: $Revision: 1.11.4.37 $'
    write(6,*)'              Tag: $Name: ccsm3_0_rel04 $'
    write(6,*)'              of the message compatability interface:'

    if ( cpl_min_vers /= expect_min_vers )then
       write(6,*) 'WARNING(cpl_compat):: Minor version of coupler ', &
            'messages different than expected: '
       write(6,*) 'The version of the coupler being used is: ',&
            cpl_min_vers
       write(6,*) 'The version I expect is:                  ',&
            expect_min_vers
    end if

    if ( cpl_maj_vers /= expect_maj_vers )then
       write(6,*) 'ERROR(cpl_compat):: Major version of coupler ', &
            'messages different than expected: '
       write(6,*) 'The version of the coupler being used is: ',&
            cpl_maj_vers
       write(6,*) 'The version I expect is:                  ',&
            expect_maj_vers
       call endrun
    end if

  end subroutine csm_compat

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: csm_restart
!
! !INTERFACE:
  subroutine csm_restart(nio, flag)
!
! !DESCRIPTION:
!  Read/write restart data needed for running in flux coupled mode
!
! !USES:
  ! use clm_varctl, only : csm_doflxave
    use colm_varctl, only : csm_doflxave
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: nio             !restart unit
    character(len=*), intent(in) :: flag   !"read" or "write"
!
! !REVISION HISTORY:
!  02.09.17  Mariana Vertenstein: moved code to be part of ccsm module
!
!EOP
!
! !LOCAL VARIABLES:
    logical :: flxave_res   !flux averaging flag read from restart file
    integer :: ier          !mpi return error code
!-----------------------------------------------------------------------

    if (flag == 'read') then
       if (masterproc) then
          read(nio) flxave_res
          read(nio) dosend
          if ((flxave_res .and. .not.csm_doflxave).or.(.not.flxave_res .and. csm_doflxave)) then
             write(6,*)'(CSM_RESTART): flxave value from namelist ',csm_doflxave, &
                  ' must be the same as flxave value from restart dataset ',flxave_res
             call endrun
          endif
          if (flxave_res .and. .not. dosend) then
             write(6,*)'(CSM_RESTART): assume that current flux coupled model ', &
                  'with flux averaging must stop on a time step where dosend (doalb) is true'
             call endrun
          end if
       endif
#if (defined SPMD)
       call mpi_bcast(dosend    , 1, MPI_LOGICAL, 0, mpicom, ier)
       if (ier /= MPI_SUCCESS) then
          write(6,*) 'MPI BCAST ERROR: for dosend in csm_restart'
          call endrun
       end if
       call mpi_bcast(flxave_res, 1, MPI_INTEGER, 0, mpicom, ier)
       if (ier /= MPI_SUCCESS) then
          write(6,*) 'MPI BCAST ERROR: for flxave_res in csm_restart'
          call endrun
       end if
#endif
    endif

    if (flag == 'write') then
       if (masterproc) then
          write(nio) csm_doflxave
          write(nio) dosend
       endif
    end if

  end subroutine csm_restart

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: global_sum_fld2d
!
! !INTERFACE:
  real(r8) function global_sum_fld2d(array, spval)
!
! !DESCRIPTION:
! Performs a global sum on an input 2d grid array
!
! !USES:
  ! use clm_varsur, only : area                 !km^2
    use colm_varMod, only : area                 !km^2
    use colm_varMod, only : lsmlon=>lon_points, lsmlat=>lat_points
!
! !ARGUMENTS:
    implicit none
    real(r8), intent(in) :: array(lsmlon,lsmlat) !W/m2, Kg/m2-s or N/m2
    real(r8), intent(in) :: spval                !points to not include in global sum
!
! !REVISION HISTORY:
!
!EOP
!
! !LOCAL VARIABLES:
    integer :: i,j         !indices
!------------------------------------------------------------------------

    global_sum_fld2d = 0.
    do j = 1,lsmlat
       do i = 1,lsmlon
          if (array(i,j) /= spval) then
             global_sum_fld2d = global_sum_fld2d + array(i,j) * area(i,j) * 1.e6
          endif
       end do
    end do

  end function global_sum_fld2d

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: global_sum_fld1d
!
! !INTERFACE:
  real(r8) function global_sum_fld1d(array)
!
! !DESCRIPTION:
! Performs a global sum on an input flux array
!
! !USES:
  ! use clmtype
  ! use clm_varsur, only : area
    use colm_varMod, only : area, numgrid_glob
    use spmd_decomp, only : gxmap_glob, gymap_glob
!
! !ARGUMENTS:
    implicit none
    real(r8), intent(in) :: array(:) !W/m2, Kg/m2-s or N/m2
!
! !REVISION HISTORY:
!
!EOP
!
! !LOCAL VARIABLES:
    integer :: g,i,j  ! indices
  ! type(gridcell_type), pointer :: gptr  ! pointer to gridcell derived type
!------------------------------------------------------------------------

    ! Set pointers into derived type

  ! gptr => clm3%g

    ! Note: area is in km^2

    global_sum_fld1d = 0.
!   do g = 1,numg
!      i = gptr%ixy(g)
!      j = gptr%jxy(g)
!      global_sum_fld1d = global_sum_fld1d + array(g) * area(i,j) * 1.e6
!   end do

    do g = 1,numgrid_glob
       i = gxmap_glob(g)
       j = gymap_glob(g)
       global_sum_fld1d = global_sum_fld1d + array(g) * area(i,j) * 1.e6
    end do

  end function global_sum_fld1d

#endif

end module colm_csmMod

#endif
