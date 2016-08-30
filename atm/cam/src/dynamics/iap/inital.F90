module inital
!-----------------------------------------------------------------------
! 
! Purpose:  CAM startup initial conditions module
! Update:  ZhangHe, 2008.4.18, add calling sub. STFRAM
!          ZhangHe, 2008.6.10, add calling sub. initialize_IAPprog
!          ZhangHe, 2008,6.12, move calling sub. STFRAM after calling sub. readinitial
! update to cesm: HeJuanxiong, 2010.08
!                 ZhangHe, 2011-11-14
!                 ZhangHe, 2011-12-17, move calling sub. STFRAM to sub. initcom
!                 Jiang Jinrong and Zhang He, 2012-10-30, for 2D parallel
!                 ZhangHe, 2013-01-16, update dyn_init
!-----------------------------------------------------------------------

   implicit none

   private   ! By default everything private to this module
!
! Public methods
!
   public cam_initial   ! Cam initialization (formally inital)

contains

!
!-----------------------------------------------------------------------
!

subroutine cam_initial( dyn_in, dyn_out, nlfilename )

!-----------------------------------------------------------------------
!
! Purpose:
! Define initial conditions for first run of case
!
! Method:
!
! Author:
! Original version:  CCM1
! Standardized:      L. Bath, June 1992
!                    T. Acker, March 1996
! Reviewed:          B. Boville, April 1996
!
!-----------------------------------------------------------------------
   use shr_kind_mod, only: r8 => shr_kind_r8
   use dyn_comp,             only: dyn_import_t, dyn_export_t
   use prognostics,          only: initialize_prognostics
   use phys_grid,            only: phys_grid_init
   use chem_surfvals,        only: chem_surfvals_init
   use camsrfexch_types,     only: atm2hub_alloc
   use scanslt,              only: scanslt_alloc
   use startup_initialconds, only: setup_initial, initial_conds
   use dyn_comp,             only: dyn_init
   use dynamics_vars,        only: dynamics_init  !jjr
   use constituents,         only: pcnst                   
   use pmgrid, only: plon, plat, plev, beglonxy, endlonxy, beglatxy, endlatxy, &
                     beglat, endlat, beglev, endlev 
! use runtime_opts   iord, jord, nsplit  !
! ========================== zhh =============================
!!   use Dyn_const,    only: STFRAM
   use IAP_prog,     only: initialize_IAPprog  !zhh 2008.6.10
   use spmd_utils,   only: masterproc
   use dynamics_vars,        only : T_FVDYCORE_STATE   !zhh 2013-01-16
   use dyn_internal_state,   only : get_dyn_state      !zhh 2013-01-16
! ============================================================

#if (defined SPMD)
!jjr   use spmd_dyn,             only: spmdbuf
#endif
!-----------------------------------------------------------------------
!
! Arguments
!
!  Arguments are not used in this dycore, included for compatibility
   type(dyn_import_t) :: dyn_in
   type(dyn_export_t) :: dyn_out
   character(len=*), intent(in) :: nlfilename

!---------------------------Local variables-----------------------------
   real(r8) :: dtime          ! timestep size
   type (T_FVDYCORE_STATE), pointer :: dyn_state   !zhh 2013-01-16
!-----------------------------------------------------------------------
!
   call setup_initial()
   !
!---------------------------- test zhh ------------------------
!!   stop
!---------------------------- test zhh ------------------------
   ! Initialize ghg surface values before default initial distributions
   ! are set in inidat.
   call chem_surfvals_init()
   !
!!   call dyn_init(nlfilename)
!
   ! Initialize dynamics
!
   dyn_state => get_dyn_state()

   call dyn_init(dyn_state, NLFileName )    !zhh 2013-01-16
!
!------------------------------jjr-----------------------
!!   dtime = get_step_size()
!!   if (masterproc) print*, 'dtime =', dtime
!!!!   call dynamics_init( dtime, iord, jord, nsplit, &
!!   call dynamics_init( plon, plat, plev, pcnst,   &   !ppcnst ==> pcnst
!!                       beglonxy, endlonxy,        &
!!                       beglatxy, endlatxy,        &
!!                       beglat,   endlat,          &
!!                       beglev,   endlev )
!!
!------------------------------jjr-----------------------
   ! Initialize prognostics variables
   !
   call initialize_prognostics
   call initialize_IAPprog     !zhh 2008.6.10
   call scanslt_alloc()
   !
   ! Set commons
   !
   call initcom
   !
!---------------------------- test zhh ------------------------
   if (masterproc) print*, 'success calling initcom'
!---------------------------- test zhh ------------------------
   ! Define physics data structures
   !
   call phys_grid_init
!---------------------------- test zhh ------------------------
   if (masterproc) print*, 'success calling phys_grid_init'
!---------------------------- test zhh ------------------------
#if (defined SPMD)
   ! Allocate communication buffers for
   ! collective communications in realloc
   ! routines and in dp_coupling
!jjr   call spmdbuf ()
!---------------------------- test zhh ------------------------
   if (masterproc) print*, 'success calling spmdbuf'
!---------------------------- test zhh ------------------------
#endif
!
   call initial_conds( dyn_in )
!---------------------------- test zhh ------------------------
   if (masterproc) print*, 'success calling initial_conds'
!---------------------------- test zhh ------------------------
   !
end subroutine cam_initial

!
!-----------------------------------------------------------------------
!

end module inital
