module stepon
!----------------------------------------------------------------------- 
! 
! Purpose: 
! Module for time-stepping of the IAP finite-difference dynamics.
! Modified:          ZhangHe, 2007.5.25
! Update:            ZhangHe, 2008.4.18, move calling of sub. STFRAM to 
!                                        sub. inital & sub. read_restart
! Update to cesm1:  HeJuanxion, 2010.10
!                   ZhangHe, 2011-11-14
!                   ZhangHe, 2011-12-16
!                   ZhangHe, 2011-12-27, removed calling sub. apply_fq
!                   ZhangHe, 2012-01-16, add get_curr_date
! Modified: Jiang Jinrong and Zhang He, 2012-11-13, for 2D parallel
!           Zhang He, 2013-01-29, new mp_sendirr, dyn_state is added
!           Zhang He, 2013-03-21, update diag_dynvar_ic
! Reviewed: Zhang He, 2013-03-20
! Modified: He Juanxiong, 2013-09-26, for nudging and two-way coupling
!-----------------------------------------------------------------------
  use shr_kind_mod,     only: r8 => shr_kind_r8
  use shr_sys_mod,      only: shr_sys_flush
  use pmgrid,           only: plev, plat, plevp, plon, beglat, endlat, &
                              twod_decomp, beglev, endlev,  &
                              beglonxy,endlonxy,beglatxy,endlatxy
  use spmd_utils,       only: masterproc
  use scanslt,          only: advection_state
  use prognostics,      only: ps, u3, v3, t3, q3, qminus, ptimelevels,   &
                              omga, phis, n3, n3m1, shift_time_indices
  use camsrfexch_types, only: cam_out_t     
#if (defined SPMD)
  use mpishorthand, only : mpicom
  use parutilitiesmodule, only : sumop, parcollective
  use mod_comm, only: mp_sendirr, mp_recvirr
#endif
  use ppgrid,           only: begchunk, endchunk
  use physics_types,    only: physics_state, physics_tend
  use time_manager,     only: is_first_step, is_first_restart_step, get_curr_date
  use infnan,           only: nan
  use iop,              only: setiopupdate, readiopdata
  use scamMod,          only: use_iop,doiopupdate,use_pert_frc,wfld,wfldh,single_column
  use perf_mod
  use times, only: time_dyn, time_phy, time_all, times_set, times_out

  implicit none

  private   ! By default make all data and methods private to this module
!
! Public methods
!
  public stepon_init     ! Initialization
  public stepon_run1     ! Run method phase 1
  public stepon_run2     ! Run method phase 2
  public stepon_run3     ! Run method phase 3
  public stepon_final    ! Finalization
  public phys_state      ! by Wang Yuzhu
!
! Private module data
!
  save
  type(physics_state), pointer :: phys_state(:)   ! Physics state data
  type(physics_tend ), pointer :: phys_tend(:)    ! Physics tendency data

  real(r8) :: detam(plev)               ! intervals between vert full levs.
  real(r8) :: cwava(plat)               ! weight applied to global integrals
  real(r8), allocatable :: t2(:,:,:)    ! temp tendency
  real(r8), allocatable :: fu(:,:,:)    ! u wind tendency
  real(r8), allocatable :: fv(:,:,:)    ! v wind tendency
  real(r8), allocatable :: flx_net(:,:) ! net flux from physics
!!  real(r8), allocatable :: fq(:,:,:,:)  ! Q tendencies,for eul_nsplit>1
!!  real(r8), allocatable :: t2_save(:,:,:)    ! temp tendency
!!  real(r8), allocatable :: fu_save(:,:,:)    ! u wind tendency
!!  real(r8), allocatable :: fv_save(:,:,:)    ! v wind tendency
  real(r8) :: coslat(plon)              ! cosine of latitude
  real(r8) :: rcoslat(plon)             ! Inverse of coseine of latitude
  real(r8) :: rpmid(plon,plev)          ! inverse of midpoint pressure
  real(r8) :: pdel(plon,plev)           ! Pressure depth of layer
  real(r8) :: pint(plon,plevp)          ! Pressure at interfaces
  real(r8) :: pmid(plon,plev)           ! Pressure at midpoint
  real(r8) :: dtime = nan               ! timestep size
  type(advection_state) :: adv_state    ! Advection state data
!-----------------------------------------------------------------------
! The following arrays are for secondary 2D x-y decomposition
!-----------------------------------------------------------------------
   real(r8), allocatable :: dummy3(:,:,:)
   real(r8), allocatable :: phisxy(:,:)       ! Surface geopotential
   real(r8), allocatable :: psxy(:,:,:)       ! Surface pressure
   real(r8), allocatable :: omgaxy(:,:,:)     ! vertical pressure velocity
   real(r8), allocatable :: u3xy(:,:,:,:)     ! Staggered grid winds, latitude
   real(r8), allocatable :: v3xy(:,:,:,:)     ! Satggered grid winds, longitude
   real(r8), allocatable :: pdeldxy(:,:,:,:)  ! delta pressure
   real(r8), allocatable :: t3xy(:,:,:,:)     ! virtual potential temperature
   real(r8), allocatable :: q3xy(:,:,:,:,:)   ! Moisture and constituents
   real(r8), allocatable :: dummy3xy(:,:,:)

!======================================================================= 
contains
!======================================================================= 

!
!======================================================================= 
!

subroutine stepon_init( gw, etamid, dyn_in, dyn_out )
!----------------------------------------------------------------------- 
! 
! Purpose:  Initialization, primarily of dynamics.
!
!----------------------------------------------------------------------- 
   use dyn_comp,       only: dyn_import_t, dyn_export_t
   use scanslt,        only: scanslt_initial
   use commap,         only: clat
   use constituents,   only: pcnst
   use physconst,      only: gravit
   use rgrid,          only: nlon
   use hycoef,         only: hyam, hybm
   use time_manager,   only: get_step_size
!!   use eul_control_mod,only: eul_nsplit
   use Dyn_const,      only: SIGL   !zhh 2011-11-14         
#if ( defined BFB_CAM_SCAM_IOP )
   use iop,            only:init_iop_fields
#endif
#if (defined SPMD)
#include <mpif.h>
#endif
!-----------------------------------------------------------------------
! Arguments
!
  real(r8), intent(out) :: gw(plat)                  ! Gaussian weights
  real(r8), intent(out) :: etamid(plev)              ! vertical coords at midpoints
  type(dyn_import_t) :: dyn_in                       ! included for compatibility
  type(dyn_export_t) :: dyn_out                      ! included for compatibility
!-----------------------------------------------------------------------
!  Local variables
!
   integer :: k, lat, i
   real(r8)  :: time_begin, time_end   ! begin & end of the CPU time
   real(r8)  :: time_begin2, time_end2 ! begin & end of the CPU time
!-----------------------------------------------------------------------

   call t_startf ('stepon_startup')

   dtime = get_step_size()
   !
   ! Define eta coordinates: Used for calculation etadot vertical velocity 
   ! for slt.
   !
   do k=1,plev
      etamid(k) = SIGL(k)
   end do

   call scanslt_initial( adv_state, etamid, gravit, gw, detam, cwava )
   !
   ! Initial guess for trajectory midpoints in spherical coords.
   ! nstep = 0:  use arrival points as initial guess for trajectory midpoints.
   ! nstep > 0:  use calculated trajectory midpoints from previous time 
   ! step as first guess.
   ! NOTE:  reduce number of iters necessary for convergence after nstep = 1.
   !
   if (is_first_step()) then
      do lat=beglat,endlat

	     if (.not. single_column) then
         !
         ! Calculate vertical motion field
         !
            omga(:,:,lat)=0.0
!
         else
         
            omga(1,:,lat)=wfld(:)
         endif
      end do

!------------------------------------
#if (defined SPMD)
      time_begin2=mpi_wtime()!wjp 2011.04.16
      call times_set         !wjp 2011.04.16  
#endif  
      call init_trans_pd     !! 
      call initial_dyn       ! zhh 2007.8.23
	  call init_trans_IAP    !! added by zhh  
      CALL DIAGHI( 0 )       !zhh 2013-01-21

      if(masterproc) write(6,*) 'this is first timestep'
!
   else if (is_first_restart_step()) then
      
#if (defined SPMD)
      time_begin2=mpi_wtime()!wjp 2011.04.16
      call times_set         !wjp 2011.04.16
#endif    
      call init_trans_pd     !! zhh 2012.11.08 
      call initial_dyn       !zhh 2007.5.25
	  call init_trans_IAP    !! zhh 2012.11.08  
      CALL DIAGHI( 0 )       !zhh 2013-01-21

	  if(masterproc) write(6,*) ' this is first restart timestep '

   end if

   allocate(t2(plon,plev,beglat:endlat))
   allocate(fu(plon,plev,beglat:endlat))
   allocate(fv(plon,plev,beglat:endlat))
   allocate( flx_net(plon,beglat:endlat))
   !
!----------------------------------------------------------
! Allocate variables for secondary 2D xy decomposition
!----------------------------------------------------------
   allocate (dummy3(plon,beglat:endlat,beglev:endlev))
   allocate (phisxy(beglonxy:endlonxy      , beglatxy:endlatxy))
   allocate (  psxy(beglonxy:endlonxy      , beglatxy:endlatxy, ptimelevels))
   allocate (omgaxy(beglonxy:endlonxy, plev, beglatxy:endlatxy))
   allocate ( u3xy(beglonxy:endlonxy,plev, beglatxy:endlatxy, ptimelevels))  ! now unghosted, was N1
   allocate ( v3xy(beglonxy:endlonxy,plev, beglatxy:endlatxy, ptimelevels))
   allocate (pdeldxy(beglonxy:endlonxy,plev, beglatxy:endlatxy, ptimelevels))
   allocate (  t3xy(beglonxy:endlonxy,plev, beglatxy:endlatxy, ptimelevels))
   allocate (  q3xy(beglonxy:endlonxy,plev,pcnst, beglatxy:endlatxy, ptimelevels)) !zhh
   allocate (dummy3xy(beglonxy:endlonxy,beglatxy:endlatxy,plev))
!
   ! Beginning of basic time step loop
   !
   call t_stopf ('stepon_startup')


#if ( defined BFB_CAM_SCAM_IOP )
   if (is_first_step()) then
      call init_iop_fields()
   endif
#endif
end subroutine stepon_init

!
!======================================================================= 
!
subroutine stepon_run1( ztodt, phys_state, phys_tend, wrf_state, wrf_tend,&
                        cam_state_ac, cam_tend_ac, dyn_in, dyn_out, &
                        nudging_state, nudging_tend, fdda, &
                        twoway_coupling, twoway_nudging)  ! juanxiong he
!----------------------------------------------------------------------- 
! 
! Purpose:  Phase 1 run method of dynamics. Set the time-step
!           to use for physics. And couple from dynamics to physics.
!
!----------------------------------------------------------------------- 
  use dyn_comp,       only: dyn_import_t, dyn_export_t
  use time_manager,   only: get_nstep
  use prognostics,    only: pdeld
  use constituents,   only: pcnst
  use phys_buffer,    only: pbuf
  use dp_coupling,    only: d_p_coupling, wrf_to_dynamics ! juanxiong he
  use eul_control_mod,only: eul_nsplit
  use dynamics_vars,   only: T_FVDYCORE_GRID, T_FVDYCORE_STATE    !zhh 2013-01-29
  use dyn_internal_state,   only : get_dyn_state      !zhh 2013-01-29

!---------------------------------Arguments------------------------------
  real(r8), intent(out) :: ztodt            ! twice time step unless nstep=0
  type(physics_state), intent(inout) :: phys_state(begchunk:endchunk)
  type(physics_tend),  intent(out)   :: phys_tend(begchunk:endchunk)
  type(physics_state), intent(inout) :: wrf_state(begchunk:endchunk) !juanxiong he
  type(physics_tend),  intent(inout)   :: wrf_tend(begchunk:endchunk)  !juanxiong he
  type(physics_state), intent(inout) :: nudging_state(begchunk:endchunk) !juanxiong he
  type(physics_tend),  intent(inout)   :: nudging_tend(begchunk:endchunk) !juanxiong he
  type(physics_state), intent(inout) :: cam_state_ac(begchunk:endchunk) !juanxiong he
  type(physics_tend),  intent(inout)   :: cam_tend_ac(begchunk:endchunk) !juanxiong he
  type(dyn_import_t) :: dyn_in                       ! included for compatibility
  type(dyn_export_t) :: dyn_out                      ! included for compatibility
  logical, intent(inout) :: twoway_coupling  !juanxiong he
  integer, intent(in) :: twoway_nudging  !juanxiong he
  integer, intent(in) :: fdda  !juanxiong he
!-----------------------------Local workspace-----------------------------
  integer :: i, j, k, l1, lchnk   ! for debug
  type (T_FVDYCORE_STATE), pointer :: dyn_state   !zhh 2013-01-29
  type (T_FVDYCORE_GRID) , pointer :: grid        ! For convenience
!----------------------------------------------------------------------- 
!
   dyn_state => get_dyn_state()
   grid => dyn_state%grid     ! For convenience

  !------------------------------------------------------------
  ! juanxiong he 
  !------------------------------------------------------------ 
  ! two way coupling
  if(twoway_coupling.and.twoway_nudging.eq.0) then
    call wrf_to_dynamics(wrf_state, ps(:,:,n3m1), t3(:,:,:,n3m1), u3(:,:,:,n3m1), &
                         v3(:,:,:,n3m1), q3(:,:,1,:,n3m1))
  endif

   if (grid%twod_decomp .eq. 1) then
#if defined( SPMD )
      do k = beglev,endlev
         do j = beglat,endlat
            do i = 1,plon
               dummy3(i,j,k) = phis(i,j)
            enddo
         enddo
      enddo
      call mp_sendirr(grid%commxy, grid%ijk_yz_to_xy%SendDesc,         &
                      grid%ijk_yz_to_xy%RecvDesc, dummy3, dummy3xy,    &
                      modc=grid%modc_dynrun )
      call mp_recvirr(grid%commxy, grid%ijk_yz_to_xy%SendDesc,         &
                      grid%ijk_yz_to_xy%RecvDesc, dummy3, dummy3xy,    &
                      modc=grid%modc_dynrun )
      do j = beglatxy,endlatxy
         do i = beglonxy,endlonxy
            phisxy(i,j) = dummy3xy(i,j,1)
         enddo
      enddo

      call mp_sendirr(grid%commxy, grid%ikj_yz_to_xy%SendDesc,         &
                      grid%ikj_yz_to_xy%RecvDesc, omga, omgaxy,        &
                      modc=grid%modc_dynrun )
      call mp_recvirr(grid%commxy, grid%ikj_yz_to_xy%SendDesc,         &
                      grid%ikj_yz_to_xy%RecvDesc, omga, omgaxy,        &
                      modc=grid%modc_dynrun )
!
      call mp_sendirr(grid%commxy, grid%ikj_yz_to_xy%SendDesc,         &
                      grid%ikj_yz_to_xy%RecvDesc, u3(:,:,:,n3m1),      &
                      u3xy(:,:,:,n3m1), modc=grid%modc_dynrun )
      call mp_recvirr(grid%commxy, grid%ikj_yz_to_xy%SendDesc,         &
                      grid%ikj_yz_to_xy%RecvDesc, u3(:,:,:,n3m1),      &
                      u3xy(:,:,:,n3m1), modc=grid%modc_dynrun )
!
      call mp_sendirr(grid%commxy, grid%ikj_yz_to_xy%SendDesc,         &
                      grid%ikj_yz_to_xy%RecvDesc, v3(:,:,:,n3m1),      &
                      v3xy(:,:,:,n3m1), modc=grid%modc_dynrun )
      call mp_recvirr(grid%commxy, grid%ikj_yz_to_xy%SendDesc,         &
                      grid%ikj_yz_to_xy%RecvDesc, v3(:,:,:,n3m1),      &
                      v3xy(:,:,:,n3m1), modc=grid%modc_dynrun )
!
      call mp_sendirr(grid%commxy, grid%ikj_yz_to_xy%SendDesc,         &
                      grid%ikj_yz_to_xy%RecvDesc, t3(:,:,:,n3m1),      &
                      t3xy(:,:,:,n3m1), modc=grid%modc_dynrun )
      call mp_recvirr(grid%commxy, grid%ikj_yz_to_xy%SendDesc,         &
                      grid%ikj_yz_to_xy%RecvDesc, t3(:,:,:,n3m1),      &
                      t3xy(:,:,:,n3m1), modc=grid%modc_dynrun )
!
      call mp_sendirr(grid%commxy, grid%ikj_yz_to_xy%SendDesc,         &
                      grid%ikj_yz_to_xy%RecvDesc, pdeld(:,:,:,n3m1),   &
                      pdeldxy(:,:,:,n3m1), modc=grid%modc_dynrun )
      call mp_recvirr(grid%commxy, grid%ikj_yz_to_xy%SendDesc,         &
                      grid%ikj_yz_to_xy%RecvDesc, pdeld(:,:,:,n3m1),   &
                      pdeldxy(:,:,:,n3m1), modc=grid%modc_dynrun )
!
      do k = beglev,endlev
         do j = beglat,endlat
            do i = 1,plon
               dummy3(i,j,k) = ps(i,j,n3m1)
            enddo
         enddo
      enddo
      call mp_sendirr(grid%commxy, grid%ijk_yz_to_xy%SendDesc,         &
                      grid%ijk_yz_to_xy%RecvDesc, dummy3, dummy3xy,    &
                      modc=grid%modc_dynrun )
      call mp_recvirr(grid%commxy, grid%ijk_yz_to_xy%SendDesc,         &
                      grid%ijk_yz_to_xy%RecvDesc, dummy3, dummy3xy,    &
                      modc=grid%modc_dynrun )
      do j = beglatxy,endlatxy
         do i = beglonxy,endlonxy
            psxy(i,j,n3m1) = dummy3xy(i,j,1)
         enddo
      end do
!
      do l1=1,pcnst           !zhh 
         do k = beglev,endlev
            do j = beglat,endlat
               do i =1,plon
                  dummy3(i,j,k) = q3(i,k,l1,j,n3m1)
               enddo
            enddo
         enddo
         call mp_sendirr(grid%commxy, grid%ijk_yz_to_xy%SendDesc,         &
                         grid%ijk_yz_to_xy%RecvDesc, dummy3, dummy3xy,    &
                         modc=grid%modc_dynrun )
         call mp_recvirr(grid%commxy, grid%ijk_yz_to_xy%SendDesc,         &
                         grid%ijk_yz_to_xy%RecvDesc, dummy3, dummy3xy,    &
                         modc=grid%modc_dynrun )
         do k=1,plev
            do j = beglatxy,endlatxy
               do i = beglonxy,endlonxy
                  q3xy(i,k,l1,j,n3m1) = dummy3xy(i,j,k)
               enddo
            end do
         enddo
      enddo !jjr end l1
!
#endif
   else !jjr need define _xy
      psxy(:,:,n3m1)=ps(:,:,n3m1)
      t3xy(:,:,:,n3m1)=t3(:,:,:,n3m1)
      u3xy(:,:,:,n3m1)=u3(:,:,:,n3m1)
      v3xy(:,:,:,n3m1)= v3(:,:,:,n3m1)
      q3xy(:,:,:,:,n3m1)=q3(:,:,:,:,n3m1)
      omgaxy=omga
      phisxy=phis
      pdeldxy(:,:,:,n3m1)=pdeld(:,:,:,n3m1)
    end if
!
      ztodt = dtime         ! ZhangHe
  !
  ! Dump state variables to IC file
  !
  call t_startf ('diag_dynvar_ic')
  call diag_dynvar_ic (phisxy, psxy(:,:,n3m1), t3xy(:,:,:,n3m1), u3xy(:,:,:,n3m1), &
                       v3xy(:,:,:,n3m1), q3xy(:,:,:,:,n3m1) )
  call t_stopf ('diag_dynvar_ic')

  !------------------------------------------------------------
  ! juanxiong he 
  !------------------------------------------------------------ 
  ! nudging
  if(fdda.gt.0) then
   nudging_state=phys_state
   nudging_tend=phys_tend
  end if
  
  !
  !----------------------------------------------------------
  ! Couple from dynamics to physics
  !----------------------------------------------------------
  !
  call t_startf ('d_p_coupling')
  call d_p_coupling (psxy(:,:,n3m1), t3xy(:,:,:,n3m1), u3xy(:,:,:,n3m1), &
                     v3xy(:,:,:,n3m1), q3xy(:,:,:,:,n3m1), &
                     omgaxy, phisxy, phys_state, phys_tend, pbuf, pdeldxy(:,:,:,n3m1))
  call t_stopf  ('d_p_coupling')
!

end subroutine stepon_run1

!
!======================================================================= 
!

subroutine stepon_run2( phys_state, phys_tend,  wrf_state, wrf_tend,&
                        cam_state_ac, cam_tend_ac, &
                        nudging_state, nudging_tend, dyn_in, dyn_out, &
                        twoway_coupling, twoway_nudging )
!----------------------------------------------------------------------- 
! 
! Purpose:  Phase 2 run method of dynamics. Couple from physics
!           to dynamics.
!
!----------------------------------------------------------------------- 
  use dyn_comp,       only: dyn_import_t, dyn_export_t
  use dp_coupling,    only: p_d_coupling, wrf_to_dynamics_tend
  type(physics_state), intent(in):: phys_state(begchunk:endchunk)
  type(physics_tend), intent(in):: phys_tend(begchunk:endchunk)
  type(physics_state), intent(in):: nudging_state(begchunk:endchunk)  !juanxiong he
  type(physics_tend), intent(in):: nudging_tend(begchunk:endchunk) !juanxiong he
  type(physics_state), intent(inout) :: cam_state_ac(begchunk:endchunk) !juanxiong he
  type(physics_tend),  intent(inout)   :: cam_tend_ac(begchunk:endchunk) !juanxiong he
  type(physics_state), intent(in):: wrf_state(begchunk:endchunk)  !juanxiong he
  type(physics_tend), intent(in):: wrf_tend(begchunk:endchunk) !juanxiong he
  type(dyn_import_t) :: dyn_in                       ! included for compatibility
  type(dyn_export_t) :: dyn_out                      ! included for compatibility
  logical, intent(in) :: twoway_coupling  !juanxiong he
  integer, intent(in) :: twoway_nudging  !juanxiong he

  call t_startf ('p_d_coupling')
  call p_d_coupling (phys_state, phys_tend, t2, fu, fv, flx_net, &
                     qminus(:,:,:,:) )     !zhh 2012-11-13
  if(twoway_nudging.eq.0.and.twoway_coupling) then
  call wrf_to_dynamics_tend(wrf_state, wrf_tend, t2, fu, fv) ! juanxiong he
  end if
  call t_stopf  ('p_d_coupling')
end subroutine stepon_run2

!
!======================================================================= 
!
subroutine stepon_run3( ztodt, etamid, cam_out, phys_state, phys_tend, &
                        dyn_in, dyn_out, cam_state, cam_tend, &
                        twoway_coupling, twoway_nudging )
!----------------------------------------------------------------------- 
! 
! Purpose:  Final phase of dynamics run method. Run the actual dynamics.
!
!----------------------------------------------------------------------- 
  use dyn_comp,       only: dyn_import_t, dyn_export_t
  use eul_control_mod,only: eul_nsplit
  use dp_coupling,    only: dynamics_to_wrf ! juanxiong he
  use prognostics,    only: pdeld  ! juanxiong he
  use phys_buffer,    only: pbuf   ! juanxiong he

  real(r8), intent(in) :: ztodt            ! twice time step unless nstep=0
  type(cam_out_t), intent(inout) :: cam_out(begchunk:endchunk)
  real(r8), intent(in) :: etamid(plev)     ! vertical coords at midpoints

  type(physics_tend), intent(inout):: phys_tend(begchunk:endchunk) ! juanxiong he
  type(physics_state), intent(inout):: cam_state(begchunk:endchunk) ! juanxiong he
  type(physics_tend), intent(inout):: cam_tend(begchunk:endchunk) ! juanxiong he
  logical, intent(in) :: twoway_coupling ! juanxiong he
  integer :: twoway_nudging ! juanxiong he

  type(physics_state), intent(in):: phys_state(begchunk:endchunk)
  type(dyn_import_t) :: dyn_in                       ! included for compatibility
  type(dyn_export_t) :: dyn_out                      ! included for compatibility
  real(r8) :: dt_dyn0,dt_dyn 
  integer :: stage
  integer :: yr, mon, day      ! year, month, and day components of a date
  integer :: ncsec             ! current time of day [seconds]
  real(r8)  :: time_begin, time_end   ! begin & end of the CPU time
#if (defined SPMD)
#include <mpif.h>
#endif
!
  if (single_column) then
     
     ! Determine whether it is time for an IOP update;
     ! doiopupdate set to true if model time step > next available IOP
     if (use_iop) then
        call setiopupdate
     end if
     
     ! Update IOP properties e.g. omega, divT, divQ
     
     if (doiopupdate) call readiopdata()
     
  endif

!----------------------------------------------------------
! DYNPKG Call the Dynamics Package
!----------------------------------------------------------
!
#if (defined SPMD)
      time_begin=mpi_wtime()
#endif
!
!============================ zhh ===============================
  if(masterproc) then
     call get_curr_date(yr, mon, day, ncsec)
     print*, 'The run time is', yr, mon, day, ncsec
  end if
!======================== 2012.01.16 =============================

  call t_startf ('dynpkg')
!====================== revised by zhh 08.04.28 ==============================
  call dynpkg(adv_state, t2, fu, fv, qminus, flx_net, t3, u3, v3, q3, ps,  &
              omga, phis, etamid, cwava, detam, dtime)


  if(twoway_coupling.and.twoway_nudging.eq.0) then
  call dynamics_to_wrf (cam_state, ps(:,:,n3), t3(:,:,:,n3), u3(:,:,:,n3),&
                        v3(:,:,:,n3), q3(:,:,1,:,n3), pdeldxy(:,:,:,n3)) ! added by Juanxiong He
  end if

  call shift_time_indices ()                  
!=============================================================================
  call t_stopf  ('dynpkg')

#if (defined SPMD)     
      time_end=mpi_wtime()
#endif
      time_dyn=time_dyn+time_end-time_begin
!
end subroutine stepon_run3


subroutine apply_fq(qminus,q3,fq,dt)
   use shr_kind_mod, only: r8 => shr_kind_r8
   use pmgrid,       only: plon, plat, plev, plevp, beglat, endlat
   use rgrid,        only: nlon
   use constituents,   only: pcnst

   real(r8), intent(in) :: q3(plon,beglev:endlev,beglat:endlat,pcnst)
   real(r8), intent(in) :: fq(plon,beglev:endlev,beglat:endlat,pcnst)
   real(r8), intent(out) :: qminus(plon,beglev:endlev,beglat:endlat,pcnst)
   real(r8), intent(in) :: dt

   !local 	
   real(r8) :: q_tmp,fq_tmp
   integer :: q,c,k,i

   do q=1,pcnst 
   do c=beglat,endlat
   do k=beglev,endlev
   do i=1,nlon(c)
      fq_tmp = dt*fq(i,k,c,q)
      q_tmp  = q3(i,k,c,q)
      ! if forcing is > 0, do nothing (it makes q less negative)
      if (fq_tmp<0 .and. q_tmp+fq_tmp<0 ) then
         ! reduce magnitude of forcing so it wont drive q negative 
         ! but we only reduce the magnitude of the forcing, dont increase
         ! its magnitude or change the sign
         
         ! if q<=0, then this will set fq=0  (q already negative)
         ! if q>0, then we know from above that fq < -q < 0, so we 
         ! can reduce the magnitive of fq by setting fq = -q:
         fq_tmp = min(-q_tmp,0d0)
      endif
      qminus(i,k,c,q) = q_tmp + fq_tmp
   enddo
   enddo
   enddo
   enddo
   	
end subroutine


!
!======================================================================= 
!

subroutine stepon_final(dyn_in, dyn_out)
!----------------------------------------------------------------------- 
! 
! Purpose:  Stepon finalization.
!
!----------------------------------------------------------------------- 
   use dyn_comp,       only: dyn_import_t, dyn_export_t
   use scanslt, only: scanslt_final
   type(dyn_import_t) :: dyn_in                       ! included for compatibility
   type(dyn_export_t) :: dyn_out                      ! included for compatibility

   call print_memusage ('End stepon')
   call scanslt_final( adv_state )
   deallocate(t2)
   deallocate(fu)
   deallocate(fv)
   deallocate(flx_net)
!
! deallocate 2D xy decomposition
   deallocate(dummy3)
   deallocate(phisxy)
   deallocate(psxy)
   deallocate(omgaxy)
   deallocate(u3xy)
   deallocate(v3xy)
   deallocate(pdeldxy)
   deallocate(t3xy)
   deallocate(q3xy)
   deallocate(dummy3xy)

end subroutine stepon_final
!
!======================================================================= 
!

End module stepon
