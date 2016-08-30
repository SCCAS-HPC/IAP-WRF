!BOP
!
! !MODULE: dyn_comp --- Dynamical Core Component
!
! !INTERFACE:

   Module dyn_comp

! !USES:

   use shr_kind_mod,       only: r8 => shr_kind_r8, r4 => shr_kind_r4
   use dynamics_vars,      only: T_FVDYCORE_GRID,            &
                                T_FVDYCORE_STATE, T_FVDYCORE_CONSTANTS
   use abortutils,         only: endrun

#if defined(SPMD)
   use mpishorthand,       only: mpicom, mpir8
   use mod_comm, only: mp_sendirr, mp_recvirr   !zhh 2013-01-30
                      
#endif
   use perf_mod
   use cam_logfile,        only: iulog
!---------------------- zhh 2013-01-17 ------------------------
   use prognostics,  only: phisxy, GHSxy, psxy, omgaxy, lammpxy, phimpxy, sigmpxy,    &
                           WSxy, u3xy, v3xy, t3xy, pdeldxy, qfcstxy, q3xy, ptimelevels,  &
                           ps, u3, v3, t3, q3, qminus, omga, phis, pdeld, Uxy, Vxy
   use IAP_prog, only: GHS, WS, U, V
   use scanslt, only: lammp, phimp, sigmp, qfcst
!---------------------- zhh 2013-01-17 ------------------------

  implicit none
  private

! !PUBLIC MEMBER FUNCTIONS:
  public dyn_init
  public yz2D_xy2D, xy2D_yz2D

! !PUBLIC DATA MEMBERS:
  public dyn_import_t, dyn_export_t, dyn_state

  type (T_FVDYCORE_STATE), save, target :: dyn_state ! to be moved up later

  type dyn_import_t
       real(r8), dimension(:,: ),    pointer     :: phis   ! Surface geopotential
       real(r8), dimension(:,: ),    pointer     :: ps     ! Surface pressure
       real(r8), dimension(:,:,:  ), pointer     :: u3s    ! U-winds (staggered)
       real(r8), dimension(:,:,:  ), pointer     :: v3s    ! V-winds (staggered)
       real(r8), dimension(:,:,:  ), pointer     :: pe     ! Pressure
       real(r8), dimension(:,:,:  ), pointer     :: pt     ! Potential temperature
       real(r8), dimension(:,:,:  ), pointer     :: t3     ! Temperatures
       real(r8), dimension(:,:,:  ), pointer     :: pk     ! Pressure to the kappa
       real(r8), dimension(:,:,:  ), pointer     :: pkz    ! Pressure to the kappa offset
       real(r8), dimension(:,:,:  ), pointer     :: delp   ! Delta pressure
       real(r8), dimension(:,:,:,:), pointer     :: tracer ! Tracers
  end type dyn_import_t

  type dyn_export_t
       real(r8), dimension(:,: ),    pointer     :: phis   ! Surface geopotential
       real(r8), dimension(:,: ),    pointer     :: ps     ! Surface pressure
       real(r8), dimension(:,:,:  ), pointer     :: u3s    ! U-winds (staggered)
       real(r8), dimension(:,:,:  ), pointer     :: v3s    ! V-winds (staggered)
       real(r8), dimension(:,:,:  ), pointer     :: pe     ! Pressure
       real(r8), dimension(:,:,:  ), pointer     :: pt     ! Potential temperature
       real(r8), dimension(:,:,:  ), pointer     :: t3     ! Temperatures
       real(r8), dimension(:,:,:  ), pointer     :: pk     ! Pressure to the kappa
       real(r8), dimension(:,:,:  ), pointer     :: pkz    ! Pressure to the kappa offset
       real(r8), dimension(:,:,:  ), pointer     :: delp   ! Delta pressure
       real(r8), dimension(:,:,:,:), pointer     :: tracer ! Tracers
       real(r8), dimension(:,:,:  ), pointer     :: peln   !
       real(r8), dimension(:,:,:  ), pointer     :: omga   ! Vertical velocity
       real(r8), dimension(:,:,:  ), pointer     :: mfx    ! Mass flux in X
       real(r8), dimension(:,:,:  ), pointer     :: mfy    ! Mass flux in Y
  end type dyn_export_t



! !DESCRIPTION: This module implements the FVCAM Dynamical Core as
!               an ESMF gridded component.  It is specific to FVCAM
!               and does not use ESMF.
!
! \paragraph{Overview}
!
!   This module contains an ESMF wrapper for the Finite-Volume
!   Dynamical Core used in the Community Atmospheric Model
!   (FVCAM). This component will hereafter be referred
!   to as the ``FVdycore'' ESMF gridded component.  FVdycore
!   consists of four sub-components,
!
!   \begin{itemize}
!      \item {\tt cd\_core:}  The C/D-grid dycore component
!      \item {\tt te\_map:}   Vertical remapping algorithm
!      \item {\tt trac2d:}    Tracer advection
!      \item {\tt benergy:}   Energy balance
!   \end{itemize}
!
!   Subsequently the ESMF component design for FV dycore
!   will be described.   
!
! \paragraph{Internal State}
!
!  FVdycore maintains an internal state consisting of the
!  following fields:  control variables
!
!   \begin{itemize}
!     \item {\tt U}:    U winds on a D-grid (m/s)
!     \item {\tt V}:    V winds on a D-grid (m/s)
!     \item {\tt PT}:   Scaled Virtual Potential Temperature (T_v/PKZ)
!     \item {\tt PE}:   Edge pressures
!     \item {\tt Q}:    Tracers
!     \item {\tt PKZ}:  Consistent mean for p^kappa
!   \end{itemize}
!
!  as well as a GRID (to be mentioned later) 
!  and same additional run-specific variables 
!  (dt, iord, jord, nsplit, nspltrac, nspltvrm -- to be mentioned later)
!
! Note: {\tt PT} is not updated if the flag {\tt CONVT} is true.
!
! The internal state is updated each time FVdycore is called.
!
! !REVISION HISTORY:
!
!   WS  05.06.10:  Adapted from FVdycore_GridCompMod
!   WS  05.09.20:  Renamed dyn_comp
!   WS  05.11.10:  Now using dyn_import/export_t containers
!   WS  06.03.01:  Removed tracertrans-related variables
!   WS  06.04.13:  dyn_state moved here from prognostics (temporary?)
!   CC  07.01.29:  Corrected calculation of OMGA
!   AM  07.10.31:  Supports overlap of trac2d and cd_core subcycles
!   Zhang, He 2013-01-16: For IAP dyncore
!   Zhang, He 2013-03-21, reviewed
!   Zhang, He 2013-04-15, added Uxy, Vxy
!
!EOP
!----------------------------------------------------------------------
!BOC

! Enumeration of DYNAMICS_IN_COUPLINGS


  logical, parameter         :: DEBUG = .true.

! The FV core is always called in its "full physics" mode.  We don't want
! the dycore to know what physics package is responsible for the forcing.
  logical, parameter         :: convt = .true.

  real(r8), parameter        :: ZERO                 = 0.0_r8
  real(r8), parameter        :: HALF                 = 0.5_r8
  real(r8), parameter        :: THREE_QUARTERS       = 0.75_r8
  real(r8), parameter        :: ONE                  = 1.0_r8
  real(r4), parameter        :: ONE_R4               = 1.0_r4
  real(r8), parameter        :: SECS_IN_SIX_HOURS    = 21600.0_r8

  character(*), parameter, public :: MODULE_NAME = "dyn_comp"
  character(*), parameter, public :: VERSION     = "$Id$" 

CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!-------------------------------------------------------------------------
!BOP
! !ROUTINE:  dyn_init --- Initialize the IAP finite-difference dynamical core
!
! !INTERFACE:

!zhh subroutine dyn_init(dyn_state, dyn_in, dyn_out, NLFileName )
 subroutine dyn_init(dyn_state, NLFileName)

! !USES:
   use constituents,         only : pcnst
   use pmgrid,               only : plon, plat, plev, plevp,                   &
                                    beglonxy, endlonxy, beglatxy, endlatxy,    &
                                    beglat,   endlat,   beglev,   endlev,      &
                                    npr_y, npr_z, nprxy_x, nprxy_y,            &
                                    twod_decomp, mod_geopk, mod_transpose,     &
                                    mod_gatscat
   use time_manager,         only : get_step_size
   use pmgrid, only : dyndecomp_set
   use dynamics_vars, only : dynamics_init
   use dycore,      only: get_resolution
#if ( defined OFFLINE_DYN )
   use metdata, only:  metdata_dyn_init
#endif
   use spmd_utils, only: npes, masterproc
#if defined(SPMD)
   use parutilitiesmodule, only : gid, parcollective, maxop
   use spmd_dyn, only : geopkdist, geopkblocks, geopk16byte,                 &
                        npes_xy, npes_yz, mpicom_xy, mpicom_yz, mpicom_nyz,  &
                        modc_sw_dynrun, modc_hs_dynrun,                      &
                        modc_send_dynrun, modc_mxreq_dynrun,                 &
                        modc_sw_cdcore, modc_hs_cdcore,                      &
                        modc_send_cdcore, modc_mxreq_cdcore,                 &
                        modc_sw_gather, modc_hs_gather,                      &
                        modc_send_gather, modc_mxreq_gather,                 &
                        modc_sw_scatter, modc_hs_scatter,                    &
                        modc_send_scatter, modc_mxreq_scatter,               &
                        modc_sw_tracer, modc_hs_tracer,                      &
                        modc_send_tracer, modc_mxreq_tracer,                 &
                        modc_onetwo, modc_tracers
   use spmd_dyn, only: spmd_readnl,spmdinit_dyn
   use mpishorthand, only: mpicom
#endif
   use hycoef,          only : hyai, hybi
   use eul_control_mod, only : dyn_eul_readnl    !zhh 3013-01-16
   use physconst,       only : omega,             & 
                               rearth,            &
                               rair,              &
                               cpair,             &
                               zvir,              &
                               pi

  implicit none

#if !defined( SPMD )
   integer :: npes_xy=1
   integer :: npes_yz=1
   integer :: mpicom=0
   integer :: mpicom_xy=0
   integer :: mpicom_yz=0
   integer :: mpicom_nyz=0
#endif

!
! !PARAMETERS:
   character(len=*)   , intent(in) :: NLFileName ! namelist file
   real(r8), parameter ::  D0_0                  =   0.0_r8
   real(r8), parameter ::  D1E5                  =   1.0e5_r8


  type (T_FVDYCORE_STATE), target    :: dyn_state

! !DESCRIPTION: Initialize the FV dynamical core
!
! !REVISION HISTORY:
!   05.06.18   Sawyer  Creation
!   06.03.03   Sawyer  Added dyn_state as argument (for reentrancy)
!   06.05.09   Sawyer  Added dyn_conservative to conserve total energy
!   2013-01-16   ZhangHe  Removed ak, bk
!   
!EOP
!==================================================================================
!BOC

! Local variables

  integer, parameter    :: MAXPES = 256

  type (T_FVDYCORE_GRID)      , pointer :: GRID      ! For convenience
  type (T_FVDYCORE_CONSTANTS) , pointer :: CONSTANTS ! For convenience
  integer              :: unit

#if defined(SPMD)
  integer :: tmp(npes)
#endif
  integer, allocatable :: jmyz(:),kmyz(:),imxy(:),jmxy(:)                         ! used for nonblocking receive

  integer :: nstep, nymd, nhms
  integer :: yr, mm, dd, h, m, s, itmp
  integer :: INT_PACK(6)

  integer             :: k
  integer             :: TE_METHOD = 0
  integer             :: NTOTQ                !  Method for total energy remapping
  integer             :: NQ                   !  No. advected tracers
  integer             :: IFIRSTXY             !  No. total tracers
  integer             :: ILASTXY
  integer             :: JFIRSTXY
  integer             :: JLASTXY
  integer             :: JFIRST
  integer             :: JLAST
  integer             :: KFIRST
  integer             :: KLAST
  integer             :: im, jm, km
  real(r8)            :: dt
  real(r8)            :: cp       ! heat capacity of air at constant pressure
  real(r8)            :: ae       ! radius of the earth (m)
!zhh  real(r8), allocatable :: ak(:), bk(:)     !  Vertical coordinates
  integer             :: ks               !  True # press. levs

! BEGIN
!
!!  allocate( ak(plev+1) )
!!  allocate( bk(plev+1) )
!!  do k = 1, plev+1
!!     ak(k) = hyai(k) * D1E5
!!     bk(k) = hybi(k)
!!     if( bk(k) == D0_0 ) ks = k-1
!!  end do  
!
! Get the layout and store directly in the GRID data structure
!
  GRID => DYN_STATE%GRID     ! For convenience
  CONSTANTS => DYN_STATE%CONSTANTS

  dt = get_step_size()

!!  call dyn_readnl(nlfilename)
  call dyn_eul_readnl(nlfilename)

#if defined(SPMD)
  call spmd_readnl(nlfilename)
  call spmdinit_dyn()
#endif

  IFIRSTXY = beglonxy
  ILASTXY  = endlonxy
  JFIRSTXY = beglatxy
  JLASTXY  = endlatxy
  JFIRST   = beglat
  JLAST    = endlat
  KFIRST   = beglev
  KLAST    = endlev
  NTOTQ  = pcnst
  NQ     = pcnst
  IM     = plon
  JM     = plat
  KM     = plev
  cp     = cpair
  ae     = rearth

! Set constants
  constants%pi    = pi
  constants%omega = omega
  constants%ae    = ae
  constants%rair  = rair
  constants%cp    = cp
  constants%cappa = rair/cpair
  constants%zvir  = zvir

  allocate (jmyz(npr_y))
  allocate (kmyz(npr_z))
  allocate (imxy(nprxy_x))
  allocate (jmxy(nprxy_y))
    
!
! SPMD-related stuff
!
#if defined(SPMD)
    grid%twod_decomp = twod_decomp
    grid%geopkdist= geopkdist
    grid%geopk16byte = geopk16byte
    grid%geopkblocks = geopkblocks
    grid%mod_method = mod_transpose
    grid%mod_geopk  = mod_geopk
    grid%mod_gatscat  = mod_gatscat

    grid%modc_dynrun(1)   = modc_sw_dynrun
    if (modc_hs_dynrun) then
       grid%modc_dynrun(2) = 1
    else
       grid%modc_dynrun(2) = 0
    endif
    if (modc_send_dynrun) then
       grid%modc_dynrun(3) = 1
    else
       grid%modc_dynrun(3) = 0
    endif
    grid%modc_dynrun(4) = modc_mxreq_dynrun

    grid%modc_cdcore(1)   = modc_sw_cdcore
    if (modc_hs_cdcore) then
       grid%modc_cdcore(2) = 1
    else
       grid%modc_cdcore(2) = 0
    endif
    if (modc_send_cdcore) then
       grid%modc_cdcore(3) = 1
    else
       grid%modc_cdcore(3) = 0
    endif
    grid%modc_cdcore(4) = modc_mxreq_cdcore

    grid%modc_gather(1)   = modc_sw_gather
    if (modc_hs_gather) then
       grid%modc_gather(2) = 1
    else
       grid%modc_gather(2) = 0
    endif
    if (modc_send_gather) then
       grid%modc_gather(3) = 1
    else
       grid%modc_gather(3) = 0
    endif
    grid%modc_gather(4) = modc_mxreq_gather

    grid%modc_scatter(1)   = modc_sw_scatter
    if (modc_hs_scatter) then
       grid%modc_scatter(2) = 1
    else
       grid%modc_scatter(2) = 0
    endif
    if (modc_send_scatter) then
       grid%modc_scatter(3) = 1
    else
       grid%modc_scatter(3) = 0
    endif
    grid%modc_scatter(4) = modc_mxreq_scatter

    grid%modc_tracer(1)   = modc_sw_tracer
    if (modc_hs_tracer) then
       grid%modc_tracer(2) = 1
    else
       grid%modc_tracer(2) = 0
    endif
    if (modc_send_tracer) then
       grid%modc_tracer(3) = 1
    else
       grid%modc_tracer(3) = 0
    endif
    grid%modc_tracer(4) = modc_mxreq_tracer

    grid%modc_onetwo       = modc_onetwo
    grid%modc_tracers      = modc_tracers

!
!  Define imxy, jmxy, jmyz, kmyz from ifirstxy, ilastxy, etc.
!
     tmp = 0
     tmp(gid+1) = ilastxy-ifirstxy+1
     call parcollective( mpicom, maxop, npes, tmp )
     imxy(1:nprxy_x) = tmp(1:nprxy_x)

     tmp = 0
     tmp(gid+1) = jlastxy-jfirstxy+1
     call parcollective( mpicom, maxop, npes, tmp )
     do k=1,nprxy_y
       jmxy(k) = tmp((k-1)*nprxy_x+1)
     enddo

     tmp = 0
     tmp(gid+1)   = jlast-jfirst+1
     call parcollective( mpicom, maxop, npes, tmp )
     jmyz(1:npr_y) = tmp(1:npr_y)

     tmp = 0
     tmp(gid+1)   = klast-kfirst+1
     call parcollective( mpicom, maxop, npes, tmp )
     do k=1,npr_z
        kmyz(k) = tmp((k-1)*npr_y+1)
     enddo

#else
!
! Sensible initializations for OMP-only  (hopefully none of these variables are used...)
!
    grid%twod_decomp = 0
    grid%geopkdist   = .false.
    grid%geopk16byte = .false.
    grid%geopkblocks = 1
    grid%mod_method  = 0
    grid%mod_geopk   = 0
    grid%mod_gatscat = 0

    grid%modc_dynrun(1)  = 0
    grid%modc_dynrun(2)  = 1
    grid%modc_dynrun(3)  = 1
    grid%modc_dynrun(4)  = -1

    grid%modc_cdcore(1)  = 0
    grid%modc_cdcore(2)  = 1
    grid%modc_cdcore(3)  = 1
    grid%modc_cdcore(4)  = -1

    grid%modc_gather(1)  = 0
    grid%modc_gather(2)  = 1
    grid%modc_gather(3)  = 1
    grid%modc_gather(4)  = -1

    grid%modc_scatter(1) = 0
    grid%modc_scatter(2) = 1
    grid%modc_scatter(3) = 1
    grid%modc_scatter(4) = -1

    grid%modc_tracer(1)  = 0
    grid%modc_tracer(2)  = 1
    grid%modc_tracer(3)  = 1
    grid%modc_tracer(4)  = -1

    grid%modc_onetwo     = 1
    grid%modc_tracers    = 0

#endif

! These are run-specific variables:  
!     DT              Time step
!     IORD            Order (mode) of X interpolation (1,..,6)
!     JORD            Order (mode) of Y interpolation (1,..,6)
!     NSPLIT          Ratio of big to small timestep (set to zero if in doubt)
!     NSPLTRAC        Ratio of big to tracer timestep
!     NSPLTVRM        Ratio of big to vertical re-mapping timestep
!

  DYN_STATE%DOTIME    = .TRUE.
  DYN_STATE%CHECK_DT  = SECS_IN_SIX_HOURS ! Check max and min every 6 hours.
  DYN_STATE%DT        = DT         ! Should this be part of state??
!!  DYN_STATE%NSPLIT    = NSPLIT
!!  DYN_STATE%NSPLTRAC  = NSPLTRAC
!!  DYN_STATE%NSPLTVRM  = NSPLTVRM
!!  DYN_STATE%IORD      = IORD
!!  DYN_STATE%JORD      = JORD
!!  DYN_STATE%KORD      = KORD
  DYN_STATE%TE_METHOD = TE_METHOD
!!  DYN_STATE%CONSV     = dyn_conservative
!!  DYN_STATE%FILTCW    = filtcw
!!  if (filtcw .gt. 0) then
!!     if (masterproc) then
!!        write (iulog,*) ' '
!!        write (iulog,*) 'Filtering of c-grid winds turned on'
!!        write (iulog,*) ' '
!!     endif
!!  endif

!
! Calculation of orders for the C grid is fixed by D-grid IORD, JORD
!!  if( iord <= 2 ) then
!!    DYN_STATE%ICD =  1
!!  else
!!    DYN_STATE%ICD = -2
!!  endif

!!  if( jord <= 2 ) then
!!    DYN_STATE%JCD =  1
!!  else
!!    DYN_STATE%JCD =  -2
!!  endif

!
! Calculate NSPLIT if it was specified as 0
!!  if ( NSPLIT <= 0 ) DYN_STATE%NSPLIT= INIT_NSPLIT(DYN_STATE%DT,IM,JM)

! Calculate NSPLTRAC if it was specified as 0
!!  if (NSPLTRAC <= 0) then
!!     if (get_resolution() == '0.23x0.31') then
!!        DYN_STATE%NSPLTRAC = max ( 1, DYN_STATE%NSPLIT/2 )
!!     else
!!        DYN_STATE%NSPLTRAC = max ( 1, DYN_STATE%NSPLIT/4 )
!!     endif
!!  endif



! Set NSPLTVRM to 1 if it was specified as 0
!!  if (NSPLTVRM <= 0) then
!!     DYN_STATE%NSPLTVRM = 1
!!  endif
!
!
! Create the dynamics interface
!
!!  call dyn_create_interface( ifirstxy, ilastxy, jfirstxy, jlastxy, &
!!                             1, km, ntotq, dyn_in, dyn_out )

!
! Now there is sufficient information to perform the dynamics initialization
! from FVCAM.  Gradually this will be removed, and all the initialization will
! be performed in this module.
!

!
! Initialize the FVDYCORE variables, which are now all in the GRID
!
! ------------------- revised by zhanghe, 2013-01-16 -------------------
  call dynamics_init( im, jm, km, nq,                         &
                      ifirstxy, ilastxy, jfirstxy, jlastxy,   &
                      jfirst, jlast, kfirst, klast,           &
                      npes_xy, npes_yz, mpicom, mpicom_xy,    &
                      mpicom_yz, mpicom_nyz,                  &
                      nprxy_x, nprxy_y, npr_y, npr_z,         &
                      imxy,    jmxy,    jmyz,    kmyz,        &
                      grid )
! ----------------------------------------------------------------------

!  Clear wall clock time clocks and global budgets

  DYN_STATE%RUN_TIMES = 0
  DYN_STATE%NUM_CALLS = 0

#if ( defined OFFLINE_DYN )
   call metdata_dyn_init(grid)
#endif
   dyndecomp_set=.true.
   call history_defaults()

!zhh   call ctem_init( NLFileName )

   deallocate (jmyz)
   deallocate (kmyz)
   deallocate (imxy)
   deallocate (jmxy)
!!   deallocate (ak)
!!   deallocate (bk)

  return

contains

!-----------------------------------------------------------------------
!BOP
! !ROUTINE:  dyn_create_interface --- create the dynamics import and export
!
! !INTERFACE:
subroutine dyn_create_interface ( I1, IN, J1, JN, K1, KN, LM, &
                                  dyn_in, dyn_out )
   use infnan, only : inf
!
! !USES:
  implicit none

! !PARAMETERS:
   integer, intent(in)                 :: I1, IN, J1, JN, K1, KN, LM
   type (dyn_import_t), intent(out)    :: dyn_in
   type (dyn_export_t), intent(out)    :: dyn_out

!EOP
!-----------------------------------------------------------------------

   integer :: l
   integer :: ierror

   allocate( dyn_in%phis( I1:IN, J1:JN  ), stat=ierror  )
   if ( ierror /= 0 ) call endrun('DYN_COMP ALLOC error: array PHIS')
   allocate( dyn_in%ps(   I1:IN, J1:JN  ), stat=ierror  )
   if ( ierror /= 0 ) call endrun('DYN_COMP ALLOC error: array PS')
   allocate( dyn_in%u3s(  I1:IN,J1:JN,K1:KN  ), stat=ierror )
   if ( ierror /= 0 ) call endrun('DYN_COMP ALLOC error: array U3S')
   allocate( dyn_in%v3s(  I1:IN,J1:JN,K1:KN  ), stat=ierror )
   if ( ierror /= 0 ) call endrun('DYN_COMP ALLOC error: array V3S')
   allocate( dyn_in%pe(   I1:IN,K1:KN+1,J1:JN  ), stat=ierror )
   if ( ierror /= 0 ) call endrun('DYN_COMP ALLOC error: array PE')
   allocate( dyn_in%pt(   I1:IN,J1:JN,K1:KN  ), stat=ierror )
   if ( ierror /= 0 ) call endrun('DYN_COMP ALLOC error: array PT')
   allocate( dyn_in%t3(   I1:IN,J1:JN,K1:KN  ), stat=ierror )
   if ( ierror /= 0 ) call endrun('DYN_COMP ALLOC error: array T3')
   allocate( dyn_in%pk(   I1:IN,J1:JN,K1:KN+1  ), stat=ierror )
   if ( ierror /= 0 ) call endrun('DYN_COMP ALLOC error: array PK')
   allocate( dyn_in%pkz(  I1:IN,J1:JN,K1:KN  ), stat=ierror )
   if ( ierror /= 0 ) call endrun('DYN_COMP ALLOC error: array PKZ')
   allocate( dyn_in%delp( I1:IN,J1:JN,K1:KN  ), stat=ierror )
   if ( ierror /= 0 ) call endrun('DYN_COMP ALLOC error: array DELP')
!
! allocate tracer contents
!
   allocate( dyn_in%tracer(I1:IN,J1:JN,K1:KN,LM), stat=ierror )
   if ( ierror /= 0 ) then
      write(iulog,*) "Allocation error", ierror, "for tracer"
      call endrun('DYN_COMP ALLOC error: array TRACER')
   endif

   dyn_in%tracer = inf     
   dyn_in%phis = inf
   dyn_in%ps = inf
   dyn_in%u3s = inf
   dyn_in%v3s = inf
   dyn_in%pe = inf
   dyn_in%pt = inf
   dyn_in%t3 = inf
   dyn_in%pk = inf
   dyn_in%pkz = inf
   dyn_in%delp = inf

!
! Output has all of these except phis
!
   dyn_out%phis => dyn_in%phis
   dyn_out%ps   => dyn_in%ps
   dyn_out%u3s  => dyn_in%u3s
   dyn_out%v3s  => dyn_in%v3s
   dyn_out%pe   => dyn_in%pe
   dyn_out%pt   => dyn_in%pt
   dyn_out%t3   => dyn_in%t3
   dyn_out%pk   => dyn_in%pk
   dyn_out%pkz  => dyn_in%pkz
   dyn_out%delp => dyn_in%delp
   dyn_out%tracer => dyn_in%tracer

!
! And several more which are not in the import container
!
   allocate( dyn_out%peln( I1:IN,K1:KN+1,J1:JN  ), stat=ierror )
   if ( ierror /= 0 ) call endrun('DYN_COMP ALLOC error: array PELN')
   allocate( dyn_out%omga( I1:IN,K1:KN,J1:JN    ), stat=ierror )
   if ( ierror /= 0 ) call endrun('DYN_COMP ALLOC error: array OMGA')
   allocate( dyn_out%mfx( I1:IN,J1:JN,K1:KN    ), stat=ierror )
   if ( ierror /= 0 ) call endrun('DYN_COMP ALLOC error: array MFX')
   allocate( dyn_out%mfy( I1:IN,J1:JN,K1:KN    ), stat=ierror )
   if ( ierror /= 0 ) call endrun('DYN_COMP ALLOC error: array MFY')

   dyn_out%peln = inf
   dyn_out%omga = inf
   dyn_out%mfx  = inf
   dyn_out%mfy  = inf



end subroutine dyn_create_interface
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
! !ROUTINE:  init_nsplit --- find proper value for nsplit if not specified
!
! !INTERFACE:
      integer function INIT_NSPLIT(dtime,im,jm) 
!
! !USES:
      implicit none

! !INPUT PARAMETERS:
      real (r8), intent(in) :: dtime      !  time step
      integer, intent(in)   :: im, jm     !  Global horizontal resolution

! !DESCRIPTION:
!
!    If nsplit=0 (module variable) then determine a good value 
!    for ns (used in dynpkg) based on resolution and the large-time-step 
!    (pdt). The user may have to set this manually if instability occurs.
!
! !REVISION HISTORY:
!   00.10.19   Lin     Creation
!   01.03.26   Sawyer  ProTeX documentation
!   01.06.10   Sawyer  Modified for dynamics_init framework
!   03.12.04   Sawyer  Moved here from dynamics_vars.  Now a function
!
!EOP
!-----------------------------------------------------------------------
!BOC
! !LOCAL VARIABLES:
      real (r8)   pdt                       ! Time-step in seconds
                                            ! Negative dt (backward in time
                                            ! integration) is allowed
      real (r8)   dim
      real (r8)   dim0                      ! base dimension
      real (r8)   dt0                       ! base time step
      real (r8)   ns0                       ! base nsplit for base dimension
      real (r8)   ns                        ! final value to be returned
      real (r8)   one                       ! equal to unity

      parameter ( dim0 = 191._r8  )
      parameter ( dt0  = 1800._r8 )
      parameter ( ns0  = 4._r8    )
      parameter ( one  = 1.0_r8   )

      pdt = int(dtime)   ! dtime is a variable internal to this module
      dim = max ( im, 2*(jm-1) )
      ns  = int ( ns0*abs(pdt)*dim/(dt0*dim0) + THREE_QUARTERS )
      ns  = max ( one, ns )   ! for cases in which dt or dim is too small

      init_nsplit = ns

      return
!EOC
      end function INIT_NSPLIT

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine history_defaults
        !----------------------------------------------------------------------- 
        ! 
        ! Purpose: 
        !
        ! Build Master Field List of all possible fields in a history file.  Each field has 
        ! associated with it a "long_name" netcdf attribute that describes what the field is, 
        ! and a "units" attribute.
        ! 
        ! Method: Call a subroutine to add each field
        ! 
        ! Author: CCM Core Group
        ! Update: Zhang He, 2013-01-16
        !         Zhang He, 2013-01-30
        !-----------------------------------------------------------------------

        use shr_kind_mod, only: r8 => shr_kind_r8, r4 => shr_kind_r4
        use constituents, only: pcnst, cnst_name, cnst_longname   !zhh
        use constituents, only: sflxnam, tendnam, fixcnam, tottnam, hadvnam, vadvnam, cnst_get_ind  !zhh
        use ppgrid,       only: pver, pverp
        use pmgrid,       only: plev, plevp
        use cam_history,  only: dyn_stagger_decomp, dyn_decomp, addfld, add_default
        use phys_control, only: phys_getopts

        implicit none

        !-----------------------------------------------------------------------
        !
        ! Local workspace
        !
        integer m                      ! Index
        integer :: ixcldice, ixcldliq  ! constituent indices for cloud liquid and ice water.
        logical :: history_budget      ! output tendencies and state variables for CAM4
                                       ! temperature, water vapor, cloud ice and cloud
                                       ! liquid budgets.
        integer :: history_budget_histfile_num  ! output history file number for budget fields

        !
        ! Call addfld to add each field to the Master Field List.
        !

        !----------------------------------------------------------------------------
        ! Dynamics variables which belong in dynamics specific initialization modules
        !----------------------------------------------------------------------------


    call addfld ('ETADOT  ','1/s ',plevp,'A','Vertical (eta) velocity',dyn_decomp)
    call addfld ('U&IC    ','m/s ',plev, 'I','Zonal wind'                                    ,dyn_decomp )
    call addfld ('V&IC    ','m/s ',plev, 'I','Meridional wind'                               ,dyn_decomp )
    call add_default ('U&IC       ',0, 'I')
    call add_default ('V&IC       ',0, 'I')

    call addfld ('PS&IC      ','Pa      ',1,    'I','Surface pressure'                              ,dyn_decomp )
    call addfld ('T&IC       ','K       ',plev, 'I','Temperature'                                   ,dyn_decomp )
    call add_default ('PS&IC      ',0, 'I')
    call add_default ('T&IC       ',0, 'I')
    do m = 1,pcnst
       call addfld (trim(cnst_name(m))//'&IC','kg/kg   ',plev, 'I',cnst_longname(m)                 ,dyn_decomp )
    end do

    do m = 1,pcnst
       call add_default(trim(cnst_name(m))//'&IC',0, 'I')
    end do
!
! Constituent tracers
!
    do m=1,pcnst
       call addfld (hadvnam(m), 'kg/kg/s ',pver, 'A',trim(cnst_name(m))//' horizontal advection tendency ',dyn_decomp)
       call addfld (vadvnam(m), 'kg/kg/s ',pver, 'A',trim(cnst_name(m))//' vertical advection tendency ',dyn_decomp)
       call addfld (tendnam(m), 'kg/kg/s ',pver, 'A',trim(cnst_name(m))//' total tendency ',dyn_decomp)
       call addfld (tottnam(m), 'kg/kg/s ',pver, 'A',trim(cnst_name(m))//' horz + vert + fixer tendency ',dyn_decomp)
       call addfld (fixcnam(m), 'kg/kg/s ',pver, 'A',trim(cnst_name(m))//' tendency due to slt fixer',dyn_decomp)
    end do
    call addfld ('DUH     ','K/s     ',plev, 'A','U horizontal diffusive heating',dyn_decomp)
    call addfld ('DVH     ','K/s     ',plev, 'A','V horizontal diffusive heating',dyn_decomp)
    call addfld ('DTH     ','K/s     ',plev, 'A','T horizontal diffusive heating',dyn_decomp)

    call addfld ('ENGYCORR','W/m2    ',plev, 'A','Energy correction for over-all conservation',dyn_decomp)
    call addfld ('TFIX    ','K/s     ',1,    'A','T fixer (T equivalent of Energy correction)',dyn_decomp)

    call addfld ('FU      ','m/s2    ',plev, 'A','Zonal wind forcing term',dyn_decomp)
    call addfld ('FV      ','m/s2    ',plev, 'A','Meridional wind forcing term',dyn_decomp)
    call addfld ('UTEND   ','m/s2    ',plev, 'A','U tendency',dyn_decomp)
    call addfld ('VTEND   ','m/s2    ',plev, 'A','V tendency',dyn_decomp)
    call addfld ('TTEND   ','K/s     ',plev, 'A','T tendency',dyn_decomp)
    call addfld ('LPSTEN  ','Pa/s    ',1,    'A','Surface pressure tendency',dyn_decomp)
    call addfld ('VAT     ','K/s     ',plev, 'A','Vertical advective tendency of T',dyn_decomp)
    call addfld ('KTOOP   ','K/s     ',plev, 'A','(Kappa*T)*(omega/P)',dyn_decomp)

    call add_default ('DTH     ', 1, ' ')
    call phys_getopts(history_budget_out = history_budget, history_budget_histfile_num_out = history_budget_histfile_num)
    if ( history_budget ) then
       call cnst_get_ind('CLDLIQ', ixcldliq)
       call cnst_get_ind('CLDICE', ixcldice)
       call add_default(hadvnam(       1), history_budget_histfile_num, ' ')
       call add_default(hadvnam(ixcldliq), history_budget_histfile_num, ' ')
       call add_default(hadvnam(ixcldice), history_budget_histfile_num, ' ')
       call add_default(vadvnam(       1), history_budget_histfile_num, ' ')
       call add_default(vadvnam(ixcldliq), history_budget_histfile_num, ' ')
       call add_default(vadvnam(ixcldice), history_budget_histfile_num, ' ')
       call add_default(fixcnam(       1), history_budget_histfile_num, ' ')
       call add_default(fixcnam(ixcldliq), history_budget_histfile_num, ' ')
       call add_default(fixcnam(ixcldice), history_budget_histfile_num, ' ')
       call add_default(tottnam(       1), history_budget_histfile_num, ' ')
       call add_default(tottnam(ixcldliq), history_budget_histfile_num, ' ')
       call add_default(tottnam(ixcldice), history_budget_histfile_num, ' ')
       call add_default(tendnam(       1), history_budget_histfile_num, ' ')
       call add_default(tendnam(ixcldliq), history_budget_histfile_num, ' ')
       call add_default(tendnam(ixcldice), history_budget_histfile_num, ' ')
       call add_default('TTEND   '       , history_budget_histfile_num, ' ')
       call add_default('TFIX    '       , history_budget_histfile_num, ' ')
       call add_default('KTOOP   '       , history_budget_histfile_num, ' ')
       call add_default('VAT     '       , history_budget_histfile_num, ' ')
       call add_default('DTH     '       , history_budget_histfile_num, ' ')
    end if

        !-----------------------------------------------------------------------
        ! End of dynamics variables
        !-----------------------------------------------------------------------

      end subroutine history_defaults

end subroutine dyn_init
!---------------------------------------------------------------------

!======================================================================================
  subroutine yz2D_xy2D

! Purpose:  Transpose yz 2D decomposition to xy 2D decomposition.
! Author: Zhang He, 2012-11-22

    use constituents,         only : pcnst
    use pmgrid,   only: beglat, endlat, plon, plat, plev, beglev, endlev, &
                        beglonxy, endlonxy, beglatxy, endlatxy, npr_z
    use IAP_grid, only: EX
    use dyn_internal_state,   only : get_dyn_state      !zhh 2013-01-24

!---------------------------Local variables-----------------------------
    real(r8), allocatable :: tmpxy3(:,:,:)
    real(r8), allocatable :: tmp2d(:,:)
    real(r8), allocatable :: tmp3d(:,:,:)
    real(r8), allocatable :: tmpxy(:,:,:)
    type (T_FVDYCORE_STATE), pointer :: dyn_state   !zhh 2013-01-24
    type (T_FVDYCORE_GRID) , pointer :: grid        ! For convenience
    integer :: i, j, jd, k, nt, m
!-----------------------------------------------------------------------
!
   dyn_state => get_dyn_state()
   grid => dyn_state%grid     ! For convenience
!
   allocate(tmpxy3(beglonxy:endlonxy, beglatxy:endlatxy,npr_z))
   allocate(tmp2d(plon, beglat:endlat))
   allocate(tmp3d(plon, beglat:endlat, beglev:endlev))
   allocate(tmpxy(beglonxy:endlonxy, beglatxy:endlatxy,plev))
!
   
   if (grid%twod_decomp .eq. 1) then
!
#if defined (SPMD)
! Embed in 3D array since transpose machinery cannot handle 2D arrays
! ---------- phis -----------
      call mp_sendirr( grid%commxy, grid%yz2d_to_xy2d%SendDesc,                   &
                       grid%yz2d_to_xy2d%RecvDesc, phis, tmpxy3,                  &
                       modc=grid%modc_dynrun )
      call mp_recvirr( grid%commxy, grid%yz2d_to_xy2d%SendDesc,                   &
                       grid%yz2d_to_xy2d%RecvDesc, phis, tmpxy3,                  &
                       modc=grid%modc_dynrun )

!$omp parallel do private(i,j)
      do j = beglatxy, endlatxy
         do i = beglonxy, endlonxy
            phisxy(i,j) = tmpxy3(i,j,1)
         enddo
      enddo
!
! ---------- GHS -----------
      do j = beglat, endlat
         jd = plat+1-j
         do i = 1, plon      
            tmp2d(i,j) = GHS(i+EX,jd)  
         end do
      end do
!
      call mp_sendirr( grid%commxy, grid%yz2d_to_xy2d%SendDesc,                   &
                       grid%yz2d_to_xy2d%RecvDesc, tmp2d, tmpxy3,                 &
                       modc=grid%modc_dynrun )
      call mp_recvirr( grid%commxy, grid%yz2d_to_xy2d%SendDesc,                   &
                       grid%yz2d_to_xy2d%RecvDesc, tmp2d, tmpxy3,                 &
                       modc=grid%modc_dynrun )

      do j = beglatxy, endlatxy
         do i = beglonxy, endlonxy
            GHSxy(i,j) = tmpxy3(i,j,1)
         enddo
      enddo
!
! ---------- ps -----------
      do nt = 1, ptimelevels
         tmp2d(:,:) = ps(:,:,nt)
!
         call mp_sendirr( grid%commxy, grid%yz2d_to_xy2d%SendDesc,                   &
                          grid%yz2d_to_xy2d%RecvDesc, tmp2d, tmpxy3,                 &
                          modc=grid%modc_dynrun )
         call mp_recvirr( grid%commxy, grid%yz2d_to_xy2d%SendDesc,                   &
                          grid%yz2d_to_xy2d%RecvDesc, tmp2d, tmpxy3,                 &
                          modc=grid%modc_dynrun )

         do j = beglatxy, endlatxy
            do i = beglonxy, endlonxy
               psxy(i,j,nt) = tmpxy3(i,j,1)
            enddo
         enddo
!
      end do    ! nt = 1, ptimelevels
!
! ---------- omga -----------
      do k = beglev, endlev
         do j = beglat, endlat
            do i = 1, plon
               tmp3d(i,j,k) = omga(i,k,j)
            end do
         end do
      end do
!
      call mp_sendirr( grid%commxy, grid%ijk_yz_to_xy%SendDesc,                    &
                       grid%ijk_yz_to_xy%RecvDesc, tmp3d, omgaxy,                  &
                       modc=grid%modc_dynrun )
      call mp_recvirr( grid%commxy, grid%ijk_yz_to_xy%SendDesc,                    &
                       grid%ijk_yz_to_xy%RecvDesc, tmp3d, omgaxy,                  &
                       modc=grid%modc_dynrun )
!
! ---------- WS -----------
      do j = beglat, endlat
         jd = plat+1-j
         do k = beglev, endlev
            do i = 1, plon      
               tmp3d(i,j,k) = WS(i+EX,k,jd)  
            end do
         end do
      end do
!
      call mp_sendirr( grid%commxy, grid%ijk_yz_to_xy%SendDesc,                    &
                       grid%ijk_yz_to_xy%RecvDesc, tmp3d, WSxy,                    &
                       modc=grid%modc_dynrun )
      call mp_recvirr( grid%commxy, grid%ijk_yz_to_xy%SendDesc,                    &
                       grid%ijk_yz_to_xy%RecvDesc, tmp3d, WSxy,                    &
                       modc=grid%modc_dynrun )
!
! ---------- U -----------
      do j = beglat, endlat
         jd = plat+1-j
         do k = beglev, endlev
            do i = 1, plon      
               tmp3d(i,j,k) = U(i+EX,k,jd)  
            end do
         end do
      end do
!
      call mp_sendirr( grid%commxy, grid%ijk_yz_to_xy%SendDesc,                    &
                       grid%ijk_yz_to_xy%RecvDesc, tmp3d, Uxy,                     &
                       modc=grid%modc_dynrun )
      call mp_recvirr( grid%commxy, grid%ijk_yz_to_xy%SendDesc,                    &
                       grid%ijk_yz_to_xy%RecvDesc, tmp3d, Uxy,                     &
                       modc=grid%modc_dynrun )
!
! ---------- U -----------
      do j = beglat, endlat
         jd = plat+1-j
         do k = beglev, endlev
            do i = 1, plon      
               tmp3d(i,j,k) = V(i+EX,k,jd)  
            end do
         end do
      end do
!
      call mp_sendirr( grid%commxy, grid%ijk_yz_to_xy%SendDesc,                    &
                       grid%ijk_yz_to_xy%RecvDesc, tmp3d, Vxy,                     &
                       modc=grid%modc_dynrun )
      call mp_recvirr( grid%commxy, grid%ijk_yz_to_xy%SendDesc,                    &
                       grid%ijk_yz_to_xy%RecvDesc, tmp3d, Vxy,                     &
                       modc=grid%modc_dynrun )
!
! ---------- lammp -----------
      do k = beglev, endlev
         do j = beglat, endlat
            do i = 1, plon
               tmp3d(i,j,k) = lammp(i,k,j)
            end do
         end do
      end do
!
      call mp_sendirr( grid%commxy, grid%ijk_yz_to_xy%SendDesc,                    &
                       grid%ijk_yz_to_xy%RecvDesc, tmp3d, lammpxy,                 &
                       modc=grid%modc_dynrun )
      call mp_recvirr( grid%commxy, grid%ijk_yz_to_xy%SendDesc,                    &
                       grid%ijk_yz_to_xy%RecvDesc, tmp3d, lammpxy,                 &
                       modc=grid%modc_dynrun )
!
! ---------- phimp -----------
      do k = beglev, endlev
         do j = beglat, endlat
            do i = 1, plon
               tmp3d(i,j,k) = phimp(i,k,j)
            end do
         end do
      end do
!
      call mp_sendirr( grid%commxy, grid%ijk_yz_to_xy%SendDesc,                    &
                       grid%ijk_yz_to_xy%RecvDesc, tmp3d, phimpxy,                 &
                       modc=grid%modc_dynrun )
      call mp_recvirr( grid%commxy, grid%ijk_yz_to_xy%SendDesc,                    &
                       grid%ijk_yz_to_xy%RecvDesc, tmp3d, phimpxy,                 &
                       modc=grid%modc_dynrun )
!
! ---------- sigmp -----------
      do k = beglev, endlev
         do j = beglat, endlat
            do i = 1, plon
               tmp3d(i,j,k) = sigmp(i,k,j)
            end do
         end do
      end do
!
      call mp_sendirr( grid%commxy, grid%ijk_yz_to_xy%SendDesc,                    &
                       grid%ijk_yz_to_xy%RecvDesc, tmp3d, sigmpxy,                 &
                       modc=grid%modc_dynrun )
      call mp_recvirr( grid%commxy, grid%ijk_yz_to_xy%SendDesc,                    &
                       grid%ijk_yz_to_xy%RecvDesc, tmp3d, sigmpxy,                 &
                       modc=grid%modc_dynrun )
!
! ---------- u3 -----------
      do nt = 1, ptimelevels
         do k = beglev, endlev
            do j = beglat, endlat
               do i = 1, plon
                  tmp3d(i,j,k) = u3(i,k,j,nt)
               end do
            end do
         end do
!      
         call mp_sendirr( grid%commxy, grid%ijk_yz_to_xy%SendDesc,                  &
                          grid%ijk_yz_to_xy%RecvDesc, tmp3d, tmpxy,                 &
                          modc=grid%modc_dynrun )
         call mp_recvirr( grid%commxy, grid%ijk_yz_to_xy%SendDesc,                  &
                          grid%ijk_yz_to_xy%RecvDesc, tmp3d, tmpxy,                 &
                          modc=grid%modc_dynrun )
!
         do k = 1, plev
            do j = beglatxy, endlatxy
               do i = beglonxy, endlonxy
                  u3xy(i,j,k,nt) = tmpxy(i,j,k)
               enddo
            enddo
         enddo
!
      end do   
!
! ---------- v3 -----------
      do nt = 1, ptimelevels
         do k = beglev, endlev
            do j = beglat, endlat
               do i = 1, plon
                  tmp3d(i,j,k) = v3(i,k,j,nt)
               end do
            end do
         end do
!      
         call mp_sendirr( grid%commxy, grid%ijk_yz_to_xy%SendDesc,                  &
                          grid%ijk_yz_to_xy%RecvDesc, tmp3d, tmpxy,                 &
                          modc=grid%modc_dynrun )
         call mp_recvirr( grid%commxy, grid%ijk_yz_to_xy%SendDesc,                  &
                          grid%ijk_yz_to_xy%RecvDesc, tmp3d, tmpxy,                 &
                          modc=grid%modc_dynrun )
!
         do k = 1, plev
            do j = beglatxy, endlatxy
               do i = beglonxy, endlonxy
                  v3xy(i,j,k,nt) = tmpxy(i,j,k)
               enddo
            enddo
         enddo
!
      end do   
!
! ---------- t3 -----------
      do nt = 1, ptimelevels
         do k = beglev, endlev
            do j = beglat, endlat
               do i = 1, plon
                  tmp3d(i,j,k) = t3(i,k,j,nt)
               end do
            end do
         end do
!      
         call mp_sendirr( grid%commxy, grid%ijk_yz_to_xy%SendDesc,                  &
                          grid%ijk_yz_to_xy%RecvDesc, tmp3d, tmpxy,                 &
                          modc=grid%modc_dynrun )
         call mp_recvirr( grid%commxy, grid%ijk_yz_to_xy%SendDesc,                  &
                          grid%ijk_yz_to_xy%RecvDesc, tmp3d, tmpxy,                 &
                          modc=grid%modc_dynrun )
!
         do k = 1, plev
            do j = beglatxy, endlatxy
               do i = beglonxy, endlonxy
                  t3xy(i,j,k,nt) = tmpxy(i,j,k)
               enddo
            enddo
         enddo
!
      end do   
!
! ---------- pdeld -----------
      do nt = 1, ptimelevels
         do k = beglev, endlev
            do j = beglat, endlat
               do i = 1, plon
                  tmp3d(i,j,k) = pdeld(i,k,j,nt)
               end do
            end do
         end do
!      
         call mp_sendirr( grid%commxy, grid%ijk_yz_to_xy%SendDesc,                  &
                          grid%ijk_yz_to_xy%RecvDesc, tmp3d, tmpxy,                 &
                          modc=grid%modc_dynrun )
         call mp_recvirr( grid%commxy, grid%ijk_yz_to_xy%SendDesc,                  &
                          grid%ijk_yz_to_xy%RecvDesc, tmp3d, tmpxy,                 &
                          modc=grid%modc_dynrun )
!
         do k = 1, plev
            do j = beglatxy, endlatxy
               do i = beglonxy, endlonxy
                  pdeldxy(i,j,k,nt) = tmpxy(i,j,k)
               enddo
            enddo
         enddo
!
      end do   
!
! ---------- qfcst -----------
      do m = 1, pcnst
         do k = beglev, endlev
            do j = beglat, endlat
               do i = 1, plon
                  tmp3d(i,j,k) = qfcst(i,k,m,j)
               end do
            end do
         end do
!      
         call mp_sendirr( grid%commxy, grid%ijk_yz_to_xy%SendDesc,                  &
                          grid%ijk_yz_to_xy%RecvDesc, tmp3d, tmpxy,                 &
                          modc=grid%modc_dynrun )
         call mp_recvirr( grid%commxy, grid%ijk_yz_to_xy%SendDesc,                  &
                          grid%ijk_yz_to_xy%RecvDesc, tmp3d, tmpxy,                 &
                          modc=grid%modc_dynrun )
!
         do k = 1, plev
            do j = beglatxy, endlatxy
               do i = beglonxy, endlonxy
                  qfcstxy(i,j,k,m) = tmpxy(i,j,k)
               enddo
            enddo
         enddo
!
      end do   
!
! ---------- q3 -----------
      do nt = 1, ptimelevels
         do m = 1, pcnst
            do k = beglev, endlev
               do j = beglat, endlat
                  do i = 1, plon
                     tmp3d(i,j,k) = q3(i,k,m,j,nt)
                  end do
               end do
            end do
!      
            call mp_sendirr( grid%commxy, grid%ijk_yz_to_xy%SendDesc,                  &
                             grid%ijk_yz_to_xy%RecvDesc, tmp3d, tmpxy,                 &
                             modc=grid%modc_dynrun )
            call mp_recvirr( grid%commxy, grid%ijk_yz_to_xy%SendDesc,                  &
                             grid%ijk_yz_to_xy%RecvDesc, tmp3d, tmpxy,                 &
                             modc=grid%modc_dynrun )
!
            do k = 1, plev
               do j = beglatxy, endlatxy
                  do i = beglonxy, endlonxy
                     q3xy(i,j,k,nt,m) = tmpxy(i,j,k)
                  enddo
               enddo
            enddo
!
         end do
      end do   
#endif
!
   else  ! grid%twod_decomp .eq. 0
!
      do j = beglat, endlat
         do i = 1, plon
            phisxy(i,j) = phis(i,j)
            do nt = 1, ptimelevels
               psxy(i,j,nt) = ps(i,j,nt)
            end do
         end do
      end do
!
      do j = beglat, endlat
         jd = plat+1-j
         do i = 1, plon      
            GHSxy(i,j) = GHS(i+EX,jd)  
         end do
      end do
!
      do j = beglat, endlat
         jd = plat+1-j
         do k = 1, plev
            do i = 1, plon      
               WSxy(i,j,k) = WS(i+EX,k,jd)  
               Uxy(i,j,k)  = U(i+EX,k,jd)  
               Vxy(i,j,k)  = V(i+EX,k,jd)  
            end do
         end do
      end do
!
      do k = 1, plev
         do j = beglat, endlat
            do i = 1, plon
               omgaxy(i,j,k) = omga(i,k,j)
               lammpxy(i,j,k) = lammp(i,k,j)
               phimpxy(i,j,k) = phimp(i,k,j)
               sigmpxy(i,j,k) = sigmp(i,k,j)
!
               u3xy(i,j,k,:) = u3(i,k,j,:)           
               v3xy(i,j,k,:) = v3(i,k,j,:)           
               t3xy(i,j,k,:) = t3(i,k,j,:)           
               pdeldxy(i,j,k,:) = pdeld(i,k,j,:)           
            end do
         end do
      end do
!
      do k = 1, plev
         do j = beglat, endlat
            do i = 1, plon
               do m = 1, pcnst
                  qfcstxy(i,j,k,m) = qfcst(i,k,m,j)
                  do nt = 1, ptimelevels
                     q3xy(i,j,k,nt,m) = q3(i,k,m,j,nt)
                  end do
			   end do
            end do
         end do
      end do
!
   endif  !  (grid%twod_decomp .eq. 1)

   deallocate(tmpxy3)
   deallocate(tmp2d)
   deallocate(tmp3d)
   deallocate(tmpxy)

  end subroutine yz2D_xy2D

!======================================================================================
  subroutine xy2D_yz2D

! Purpose:  Transpose xy 2D decomposition to yz 2D decomposition.
! Author: Zhang He, 2012-11-23

    use constituents,         only : pcnst
    use pmgrid,   only: beglat, endlat, plon, plat, plev, beglev, endlev, &
                        beglonxy, endlonxy, beglatxy, endlatxy, npr_z
    use IAP_grid, only: EX, NL, NZ, period
    use dyn_internal_state,   only : get_dyn_state      !zhh 2013-01-24
   use spmd_utils, only: masterproc

!---------------------------Local variables-----------------------------
    real(r8), allocatable :: tmpxy3(:,:,:)
    real(r8), allocatable :: tmp2d(:,:)
    real(r8), allocatable :: tmp3d(:,:,:)
    real(r8), allocatable :: tmpxy(:,:,:)
    type (T_FVDYCORE_STATE), pointer :: dyn_state   !zhh 2013-01-24
    type (T_FVDYCORE_GRID) , pointer :: grid        ! For convenience
    integer :: i, j, jd, k, nt, m
!------------------------------------------------------------------------------------
!
   dyn_state => get_dyn_state()
   grid => dyn_state%grid     ! For convenience
!
   allocate(tmpxy3(beglonxy:endlonxy, beglatxy:endlatxy,npr_z))
   allocate(tmp2d(plon, beglat:endlat))
   allocate(tmp3d(plon, beglat:endlat, beglev:endlev))
   allocate(tmpxy(beglonxy:endlonxy, beglatxy:endlatxy,plev))
!
    if (grid%twod_decomp .eq. 1) then
!
#if defined (SPMD)
! ---------- phis -----------
       do k = 1, npr_z
          do j = beglatxy, endlatxy
             do i = beglonxy, endlonxy
                tmpxy3(i,j,k) = phisxy(i,j)
             end do
          end do
       end do
!
       call mp_sendirr(grid%commxy, grid%xy2d_to_yz2d%SendDesc,                   &
                       grid%xy2d_to_yz2d%RecvDesc, tmpxy3, phis,                  &
                       modc=grid%modc_dynrun )
       call mp_recvirr(grid%commxy, grid%xy2d_to_yz2d%SendDesc,                   &
                       grid%xy2d_to_yz2d%RecvDesc, tmpxy3, phis,                  &
                       modc=grid%modc_dynrun )
!
! ---------- GHS -----------
       do k = 1, npr_z
          do j = beglatxy, endlatxy
             do i = beglonxy, endlonxy
                tmpxy3(i,j,k) = GHSxy(i,j)
             end do
          end do
       end do
!
       call mp_sendirr(grid%commxy, grid%xy2d_to_yz2d%SendDesc,                   &
                       grid%xy2d_to_yz2d%RecvDesc, tmpxy3, tmp2d,                 &
                       modc=grid%modc_dynrun )
       call mp_recvirr(grid%commxy, grid%xy2d_to_yz2d%SendDesc,                   &
                       grid%xy2d_to_yz2d%RecvDesc, tmpxy3, tmp2d,                 &
                       modc=grid%modc_dynrun )
!
       do j = beglat, endlat
          jd = plat+1-j
          do i = 1, plon
             GHS(i+EX,jd) = tmp2d(i,j)
          end do
          call period( GHS(1,jd) )
       end do
!
! ---------- ps -----------
       do nt = 1, ptimelevels
          do k = 1, npr_z
             do j = beglatxy, endlatxy
                do i = beglonxy, endlonxy
                   tmpxy3(i,j,k) = psxy(i,j,nt)
                end do
             end do
          end do
!
          call mp_sendirr(grid%commxy, grid%xy2d_to_yz2d%SendDesc,                   &
                          grid%xy2d_to_yz2d%RecvDesc, tmpxy3, tmp2d,                 &
                          modc=grid%modc_dynrun )
          call mp_recvirr(grid%commxy, grid%xy2d_to_yz2d%SendDesc,                   &
                          grid%xy2d_to_yz2d%RecvDesc, tmpxy3, tmp2d,                 &
                          modc=grid%modc_dynrun )
!
          ps(:,:,nt) = tmp2d(:,:)
       end do    ! nt = 1, ptimelevels
!
! ---------- omga -----------
       call mp_sendirr( grid%commxy, grid%ijk_xy_to_yz%SendDesc,                    &
                        grid%ijk_xy_to_yz%RecvDesc, omgaxy, tmp3d,                   &
                        modc=grid%modc_dynrun )
       call mp_recvirr( grid%commxy, grid%ijk_xy_to_yz%SendDesc,                    &
                        grid%ijk_xy_to_yz%RecvDesc, omgaxy, tmp3d,                  &
                        modc=grid%modc_dynrun )
!
       do k = beglev, endlev
          do j = beglat, endlat
             do i = 1, plon
                omga(i,k,j) = tmp3d(i,j,k)
             end do
          end do
       end do
!
! ---------- WS -----------
       call mp_sendirr( grid%commxy, grid%ijk_xy_to_yz%SendDesc,                    &
                        grid%ijk_xy_to_yz%RecvDesc, WSxy, tmp3d,                   &
                        modc=grid%modc_dynrun )
       call mp_recvirr( grid%commxy, grid%ijk_xy_to_yz%SendDesc,                    &
                        grid%ijk_xy_to_yz%RecvDesc, WSxy, tmp3d,                  &
                        modc=grid%modc_dynrun )
!
      do j = beglat, endlat
         jd = plat+1-j
         do k = beglev, endlev
            do i = 1, plon      
               WS(i+EX,k,jd) = tmp3d(i,j,k) 
            end do
            call period( WS(1,k,jd) )
         end do
      end do
! set WS to zero at top and bottom
      if(beglev==1) WS(:,1,:) = 0.0
      if(endlev==NL) WS(:,NZ,:) = 0.0
!
! ---------- U -----------
       call mp_sendirr( grid%commxy, grid%ijk_xy_to_yz%SendDesc,                    &
                        grid%ijk_xy_to_yz%RecvDesc, Uxy, tmp3d,                   &
                        modc=grid%modc_dynrun )
       call mp_recvirr( grid%commxy, grid%ijk_xy_to_yz%SendDesc,                    &
                        grid%ijk_xy_to_yz%RecvDesc, Uxy, tmp3d,                  &
                        modc=grid%modc_dynrun )
!
      do j = beglat, endlat
         jd = plat+1-j
         do k = beglev, endlev
            do i = 1, plon      
               U(i+EX,k,jd) = tmp3d(i,j,k) 
            end do
            call period( U(1,k,jd) )
         end do
      end do
!
! ---------- V -----------
       call mp_sendirr( grid%commxy, grid%ijk_xy_to_yz%SendDesc,                    &
                        grid%ijk_xy_to_yz%RecvDesc, Vxy, tmp3d,                   &
                        modc=grid%modc_dynrun )
       call mp_recvirr( grid%commxy, grid%ijk_xy_to_yz%SendDesc,                    &
                        grid%ijk_xy_to_yz%RecvDesc, Vxy, tmp3d,                  &
                        modc=grid%modc_dynrun )
!
      do j = beglat, endlat
         jd = plat+1-j
         do k = beglev, endlev
            do i = 1, plon      
               V(i+EX,k,jd) = tmp3d(i,j,k) 
            end do
            call period( V(1,k,jd) )
         end do
      end do
!
! ---------- lammp -----------
       call mp_sendirr( grid%commxy, grid%ijk_xy_to_yz%SendDesc,                    &
                        grid%ijk_xy_to_yz%RecvDesc, lammpxy, tmp3d,                   &
                        modc=grid%modc_dynrun )
       call mp_recvirr( grid%commxy, grid%ijk_xy_to_yz%SendDesc,                    &
                        grid%ijk_xy_to_yz%RecvDesc, lammpxy, tmp3d,                  &
                        modc=grid%modc_dynrun )
!
       do k = beglev, endlev
          do j = beglat, endlat
             do i = 1, plon
                lammp(i,k,j) = tmp3d(i,j,k)
             end do
          end do
       end do
!
! ---------- phimp -----------
       call mp_sendirr( grid%commxy, grid%ijk_xy_to_yz%SendDesc,                    &
                        grid%ijk_xy_to_yz%RecvDesc, phimpxy, tmp3d,                   &
                        modc=grid%modc_dynrun )
       call mp_recvirr( grid%commxy, grid%ijk_xy_to_yz%SendDesc,                    &
                        grid%ijk_xy_to_yz%RecvDesc, phimpxy, tmp3d,                  &
                        modc=grid%modc_dynrun )
!
       do k = beglev, endlev
          do j = beglat, endlat
             do i = 1, plon
                phimp(i,k,j) = tmp3d(i,j,k)
             end do
          end do
       end do
!
! ---------- sigmp -----------
       call mp_sendirr( grid%commxy, grid%ijk_xy_to_yz%SendDesc,                    &
                        grid%ijk_xy_to_yz%RecvDesc, sigmpxy, tmp3d,                   &
                        modc=grid%modc_dynrun )
       call mp_recvirr( grid%commxy, grid%ijk_xy_to_yz%SendDesc,                    &
                        grid%ijk_xy_to_yz%RecvDesc, sigmpxy, tmp3d,                  &
                        modc=grid%modc_dynrun )
!
       do k = beglev, endlev
          do j = beglat, endlat
             do i = 1, plon
                sigmp(i,k,j) = tmp3d(i,j,k)
             end do
          end do
       end do
!
! ---------- u3 -----------
       do nt = 1, ptimelevels
          do k = 1, plev
             do j = beglatxy, endlatxy
                do i = beglonxy, endlonxy
                   tmpxy(i,j,k) = u3xy(i,j,k,nt)
                enddo
             enddo
          enddo
!
          call mp_sendirr( grid%commxy, grid%ijk_xy_to_yz%SendDesc,                    &
                           grid%ijk_xy_to_yz%RecvDesc, tmpxy, tmp3d,                   &
                           modc=grid%modc_dynrun )
          call mp_recvirr( grid%commxy, grid%ijk_xy_to_yz%SendDesc,                    &
                           grid%ijk_xy_to_yz%RecvDesc, tmpxy, tmp3d,                  &
                           modc=grid%modc_dynrun )
!
          do k = beglev, endlev
             do j = beglat, endlat
                do i = 1, plon
                   u3(i,k,j,nt) = tmp3d(i,j,k)
                end do
             end do
          end do
!
       end do   
!
! ---------- v3 -----------
       do nt = 1, ptimelevels
          do k = 1, plev
             do j = beglatxy, endlatxy
                do i = beglonxy, endlonxy
                   tmpxy(i,j,k) = v3xy(i,j,k,nt)
                enddo
             enddo
          enddo
!
          call mp_sendirr( grid%commxy, grid%ijk_xy_to_yz%SendDesc,                    &
                           grid%ijk_xy_to_yz%RecvDesc, tmpxy, tmp3d,                   &
                           modc=grid%modc_dynrun )
          call mp_recvirr( grid%commxy, grid%ijk_xy_to_yz%SendDesc,                    &
                           grid%ijk_xy_to_yz%RecvDesc, tmpxy, tmp3d,                  &
                           modc=grid%modc_dynrun )
!
          do k = beglev, endlev
             do j = beglat, endlat
                do i = 1, plon
                   v3(i,k,j,nt) = tmp3d(i,j,k)
                end do
             end do
          end do
!
       end do   
!
! ---------- t3 -----------
       do nt = 1, ptimelevels
          do k = 1, plev
             do j = beglatxy, endlatxy
                do i = beglonxy, endlonxy
                   tmpxy(i,j,k) = t3xy(i,j,k,nt)
                enddo
             enddo
          enddo
!
          call mp_sendirr( grid%commxy, grid%ijk_xy_to_yz%SendDesc,                    &
                           grid%ijk_xy_to_yz%RecvDesc, tmpxy, tmp3d,                   &
                           modc=grid%modc_dynrun )
          call mp_recvirr( grid%commxy, grid%ijk_xy_to_yz%SendDesc,                    &
                           grid%ijk_xy_to_yz%RecvDesc, tmpxy, tmp3d,                  &
                           modc=grid%modc_dynrun )
!
          do k = beglev, endlev
             do j = beglat, endlat
                do i = 1, plon
                   t3(i,k,j,nt) = tmp3d(i,j,k)
                end do
             end do
          end do
!
       end do   
!
! ---------- pdeld -----------
       do nt = 1, ptimelevels
          do k = 1, plev
             do j = beglatxy, endlatxy
                do i = beglonxy, endlonxy
                   tmpxy(i,j,k) = pdeldxy(i,j,k,nt)
                enddo
             enddo
          enddo
!
          call mp_sendirr( grid%commxy, grid%ijk_xy_to_yz%SendDesc,                    &
                           grid%ijk_xy_to_yz%RecvDesc, tmpxy, tmp3d,                   &
                           modc=grid%modc_dynrun )
          call mp_recvirr( grid%commxy, grid%ijk_xy_to_yz%SendDesc,                    &
                           grid%ijk_xy_to_yz%RecvDesc, tmpxy, tmp3d,                  &
                           modc=grid%modc_dynrun )
!
          do k = beglev, endlev
             do j = beglat, endlat
                do i = 1, plon
                   pdeld(i,k,j,nt) = tmp3d(i,j,k)
                end do
             end do
          end do
!
       end do   
!
! ---------- qfcst -----------
       do m = 1, pcnst
          do k = 1, plev
             do j = beglatxy, endlatxy
                do i = beglonxy, endlonxy
                   tmpxy(i,j,k) = qfcstxy(i,j,k,m)
                enddo
             enddo
          enddo
!
          call mp_sendirr( grid%commxy, grid%ijk_xy_to_yz%SendDesc,                    &
                           grid%ijk_xy_to_yz%RecvDesc, tmpxy, tmp3d,                   &
                           modc=grid%modc_dynrun )
          call mp_recvirr( grid%commxy, grid%ijk_xy_to_yz%SendDesc,                    &
                           grid%ijk_xy_to_yz%RecvDesc, tmpxy, tmp3d,                  &
                           modc=grid%modc_dynrun )
!
          do k = beglev, endlev
             do j = beglat, endlat
                do i = 1, plon
                   qfcst(i,k,m,j) = tmp3d(i,j,k)
                end do
             end do
          end do
!
       end do   
!
! ---------- q3 -----------
       do nt = 1, ptimelevels
          do m = 1, pcnst
             do k = 1, plev
                do j = beglatxy, endlatxy
                   do i = beglonxy, endlonxy
                      tmpxy(i,j,k) = q3xy(i,j,k,nt,m)
                   enddo
                enddo
             enddo
!
             call mp_sendirr( grid%commxy, grid%ijk_xy_to_yz%SendDesc,                    &
                              grid%ijk_xy_to_yz%RecvDesc, tmpxy, tmp3d,                   &
                              modc=grid%modc_dynrun )
             call mp_recvirr( grid%commxy, grid%ijk_xy_to_yz%SendDesc,                    &
                              grid%ijk_xy_to_yz%RecvDesc, tmpxy, tmp3d,                  &
                              modc=grid%modc_dynrun )
!
             do k = beglev, endlev
                do j = beglat, endlat
                   do i = 1, plon
                      q3(i,k,m,j,nt) = tmp3d(i,j,k)
                   end do
                end do
             end do
!
          end do
       end do   
!
#endif
! 
    else  ! grid%twod_decomp .eq. 0
!
      do j = beglat, endlat
         do i = 1, plon
            phis(i,j) = phisxy(i,j)
            do nt = 1, ptimelevels
               ps(i,j,nt) = psxy(i,j,nt)
            end do
         end do
      end do
!
       do j = beglat, endlat
          jd = plat+1-j
          do i = 1, plon      
             GHS(i+EX,jd) = GHSxy(i,j) 
          end do
          call period( GHS(1,jd) )
       end do
!
       do j = beglat, endlat
          jd = plat+1-j
          do k = 1, plev
             do i = 1, plon      
                WS(i+EX,k,jd) = WSxy(i,j,k) 
                U(i+EX,k,jd)  = Uxy(i,j,k) 
                V(i+EX,k,jd)  = Vxy(i,j,k) 
             end do
             call period( WS(1,k,jd) )
             call period( U(1,k,jd) )
             call period( V(1,k,jd) )
          end do
       end do
!      set WS to zero at top and bottom
       WS(:,1,:) = 0.0
       WS(:,NZ,:) = 0.0
!
       do k = 1, plev
          do j = beglat, endlat
             do i = 1, plon
                omga(i,k,j) = omgaxy(i,j,k)
                lammp(i,k,j) = lammpxy(i,j,k)
                phimp(i,k,j) = phimpxy(i,j,k)
                sigmp(i,k,j) = sigmpxy(i,j,k)
!
                u3(i,k,j,:) = u3xy(i,j,k,:)           
                v3(i,k,j,:) = v3xy(i,j,k,:)           
                t3(i,k,j,:) = t3xy(i,j,k,:)           
                pdeld(i,k,j,:) = pdeldxy(i,j,k,:)           
             end do
          end do
       end do
!
       do k = 1, plev
          do j = beglat, endlat
             do i = 1, plon
                do m = 1, pcnst
                   qfcst(i,k,m,j) = qfcstxy(i,j,k,m)
                   do nt = 1, ptimelevels
                      q3(i,k,m,j,nt) = q3xy(i,j,k,nt,m)
                   end do
                end do
             end do
          end do
       end do
!
    endif  !  (grid%twod_decomp .eq. 1)

   deallocate(tmpxy3)
   deallocate(tmp2d)
   deallocate(tmp3d)
   deallocate(tmpxy)

  end subroutine xy2D_yz2D

!-----------------------------------------------------------------------


end module dyn_comp



