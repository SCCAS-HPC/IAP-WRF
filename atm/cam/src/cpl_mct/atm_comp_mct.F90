module atm_comp_mct

  use pio              , only: file_desc_t, io_desc_t, var_desc_t, pio_double, pio_def_dim, &
                               pio_put_att, pio_enddef, pio_initdecomp, pio_read_darray, pio_freedecomp, &
                               pio_closefile, pio_write_darray, pio_def_var, pio_inq_varid, &
	                       pio_noerr, pio_bcast_error, pio_internal_error, pio_seterrorhandling 
  use mct_mod
  use esmf_mod
  use seq_flds_mod
  use seq_cdata_mod
  use seq_infodata_mod
  use seq_timemgr_mod

  use shr_kind_mod     , only: r8 => shr_kind_r8, cl=>shr_kind_cl
  use shr_file_mod     , only: shr_file_getunit, shr_file_freeunit, &
                               shr_file_setLogUnit, shr_file_setLogLevel, &
                               shr_file_getLogUnit, shr_file_getLogLevel, &
		               shr_file_setIO
  use shr_sys_mod      , only: shr_sys_flush, shr_sys_abort

  use cam_cpl_indices
  use cam_comp
  use cam_control_mod  , only: nsrest, adiabatic, ideal_phys, aqua_planet, eccen, obliqr, lambm0, mvelpp
  use radiation        , only: radiation_get, radiation_do, radiation_nextsw_cday
  use phys_grid        , only: get_ncols_p, get_gcol_all_p, & 
                               ngcols, get_gcol_p, get_rlat_all_p, &
	                       get_rlon_all_p, get_area_all_p, &
                               get_rlon_all_p, get_area_all_p, get_lat_all_p, get_lon_all_p  ! juanxiong he
  use ppgrid           , only: pcols, begchunk, endchunk, pverp, pver ! for geopotential computation, juanxiong he 
  use dyn_grid         , only: get_horiz_grid_dim_d
  use camsrfexch_types , only: cam_out_t, cam_in_t     
  use cam_restart      , only: get_restcase, get_restartdir
  use cam_history      , only: outfld, ctitle
  use abortutils       , only: endrun
  use filenames        , only: interpret_filename_spec, caseid, brnch_retain_casename
#ifdef SPMD
  use spmd_utils       , only: spmdinit, masterproc, iam
  use mpishorthand     , only: mpicom
#else
  use spmd_utils       , only: spmdinit, masterproc, mpicom, iam
#endif
  use time_manager     , only: get_curr_calday, advance_timestep, get_curr_date, get_nstep, &
                               is_first_step, get_step_size, timemgr_init, timemgr_check_restart
  use ioFileMod             
  use perf_mod
  use cam_logfile      , only: iulog
  use co2_cycle        , only: c_i, co2_readFlux_ocn, co2_readFlux_fuel, co2_transport, &
                               co2_time_interp_ocn, co2_time_interp_fuel, data_flux_ocn, data_flux_fuel
  use physconst       ,  only: mwco2
  use runtime_opts     , only: read_namelist
  use phys_control     , only: cam_chempkg_is
  use physics_types    , only: physics_state, physics_tend, physics_ptend  ! juanxiong he
  use coupling_utils   !by Huiqun Hao
!
! !PUBLIC TYPES:
  implicit none
  save
  private ! except   !by Wang Yuzhu
  public :: cam_state, cam_in   !by Huiqun Hao
  public :: mgrid1, mgrid2   !by Huiqun Hao
!--------------------------------------------------------------------------
! Public interfaces
!--------------------------------------------------------------------------

  public :: atm_init_mct
  public :: atm_run_mct
  public :: atm_final_mct
  public :: wrf_read_srfrest_mct  ! juanxiong he
  public :: wrf_write_srfrest_mct ! juanxiong he
  public :: wrf_read_srfrest_mct_app ! juanxiong he
  public :: wrf_write_srfrest_mct_app ! juanxiong he
  public :: gea_read_srfrest_mct  ! juanxiong he
  public :: gea_write_srfrest_mct ! juanxiong he
  public :: gea_read_srfrest_mct_app ! juanxiong he
  public :: gea_write_srfrest_mct_app ! juanxiong he
  public :: pack_mgrid !huiqun hao
!--------------------------------------------------------------------------
! Private interfaces
!--------------------------------------------------------------------------

  private :: atm_SetgsMap_mct
  private :: atm_import_mct
  private :: atm_export_mct
  private :: atm_domain_mct
  private :: atm_read_srfrest_mct
  private :: atm_write_srfrest_mct
  private :: atm_read_wrfrest_mct  ! jaunxiong he
  private :: atm_write_wrfrest_mct ! juanxiong he
  
!--------------------------------------------------------------------------------------
! state and flux variables for wrf/cam coupling, added juanxiong he 
!--------------------------------------------------------------------------------------  
  private :: atm_SetgsMap_wrf ! juanxiong he for wrf/cam coupling, 3d
  private :: atm_domain_wrf ! juanxiong he for wrf/cam coupling, 3d
  private :: atm_import_wrf  ! juanxiong he for wrf/cam coupling, 3d
  private :: atm_export_wrf ! juanxiong he for wrf/cam coupling, 3d
!--------------------------------------------------------------------------------------
! state and flux variables for wrf/cam coupling, added juanxiong he 
!--------------------------------------------------------------------------------------  
!--------------------------------------------------------------------------------------
! state and flux variables for geatm/cam coupling, added juanxiong he 
!--------------------------------------------------------------------------------------  
  private :: atm_SetgsMap_geatm ! juanxiong he for geatm/cam coupling, 3d
  private :: atm_domain_geatm ! juanxiong he for geatm/cam coupling, 3d
  private :: atm_import_geatm  ! juanxiong he for geatm/cam coupling, 3d
  private :: atm_export_geatm ! juanxiong he for geatm/cam coupling, 3d
!--------------------------------------------------------------------------------------
! state and flux variables for geatm/cam coupling, added juanxiong he 
!--------------------------------------------------------------------------------------  

!--------------------------------------------------------------------------
! Private data
!--------------------------------------------------------------------------

  type(cam_in_t) , pointer :: cam_in(:)
  type(cam_out_t), pointer :: cam_out(:)
!--------------------------------------------------------------------------------------
! state and flux variables for wrf/cam coupling, added juanxiong he 
!--------------------------------------------------------------------------------------  
  type(physics_state), pointer :: cam_state(:) ! juaxniong he, impoart cam state variables
  type(physics_tend), pointer :: cam_tend(:)  ! juanxiong he, import cam tend variables
  type(physics_state), pointer :: wrf_state(:) ! juaxniong he, import wrf state variables
  type(physics_tend), pointer :: wrf_tend(:)  ! juanxiong he, import wrf tend variables
!--------------------------------------------------------------------------------------
! state and flux variables for wrf/cam coupling, added juanxiong he 
!--------------------------------------------------------------------------------------

  type(mct_aVect)   :: a2x_a_SNAP
  type(mct_aVect)   :: a2x_a_SUM

  integer, parameter  :: nlen = 256     ! Length of character strings
  character(len=nlen) :: fname_srf_cam  ! surface restart filename
  character(len=nlen) :: pname_srf_cam  ! surface restart full pathname
  character(len=nlen) :: fname_srf_wrf  ! surface restart filename, juanxiong he
  character(len=nlen) :: pname_srf_wrf  ! surface restart full pathname, juaxniong he

  character(len=nlen) :: fname_srf_camcam  ! surface restart filename, juanxiong he
  character(len=nlen) :: pname_srf_camcam  ! surface restart full pathname, juanxiong he
  character(len=nlen) :: fname_srf_gea  ! surface restart filename, juanxiong he
  character(len=nlen) :: pname_srf_gea  ! surface restart full pathname, juaxniong he
  
!
! Filename specifier for restart surface file
! (%c = caseid, $y = year, $m = month, $d = day, $s = seconds in day, %t = tape number)
!
  character(len=*), parameter :: rsfilename_spec_cam = '%c.cam2.rs.%y-%m-%d-%s.nc' ! cam srf restarts
  character(len=*), parameter :: wrsfilename_spec_cam = '%c.cam2.rw.%y-%m-%d-%s.nc' ! cam restart from wrf, juanxiong
  character(len=*), parameter :: wrsfilename_spec_wrf = '%c.wrf.rw.%y-%m-%d-%s.nc' ! wrf restarts, juanxiong
  character(len=*), parameter :: wrsfilename_spec_wrfold = '%c.wrfold.rw.%y-%m-%d-%s.nc' ! wrf restarts, juanxiong
  character(len=*), parameter :: wrsfilename_spec_camcam = '%c.cam2.rge.%y-%m-%d-%s.nc' ! cam restart from wrf, juanxiong
  character(len=*), parameter :: wrsfilename_spec_gea = '%c.gea.rge.%y-%m-%d-%s.nc' ! geatm restarts, juanxiong
  character(len=*), parameter :: wrsfilename_spec_geaold = '%c.geaold.rge.%y-%m-%d-%s.nc' ! geatm restarts, juanxiong

  integer :: nrg = -1   ! Logical unit number for cam srf restart dataset
!
! Time averaged counter for flux fields
!
  integer :: avg_count
!
! Time averaged flux fields
!  
  character(*), parameter :: a2x_avg_flds = "Faxa_rainc:Faxa_rainl:Faxa_snowc:Faxa_snowl"  
!
! Are all surface types present   
!
  logical :: lnd_present ! if true => land is present
  logical :: ocn_present ! if true => ocean is present

  integer  :: id_isop, id_c10h16
  type(metgrid_info_atm) ::mgrid1, mgrid2   !by Huiqun Hao  for wrf/cam coupling directly
!
!================================================================================
CONTAINS
!================================================================================

  subroutine atm_init_mct( EClock, cdata_a, cdata_c, cdata_ca, x2a_a, a2x_a, &
                           x2c_c1, x2c_c2, c2x_c1, c2x_c2, &
                           x2ca_caca1, x2ca_caca2, ca2x_caca1, ca2x_caca2, & ! juanxiong he for geatm 
                           twoway_coupling, twoway_nudging, NLFilename ) !juanxiong he
    use spmd_dyn,only:iam
    use constituents,  only: cnst_get_ind
    use pmgrid,          only: plev, plevp
    !-----------------------------------------------------------------------
    !
    ! Arguments
    !
    type(ESMF_Clock),intent(in)                 :: EClock
    type(seq_cdata), intent(inout)              :: cdata_a
    type(seq_cdata), intent(inout)              :: cdata_c !juanxiong he, for wrf/cam
    type(seq_cdata), intent(inout)              :: cdata_ca !juanxiong he, for geatm/cam
    type(mct_aVect), intent(inout)              :: x2a_a
    type(mct_aVect), intent(inout)              :: a2x_a   
    type(mct_aVect), intent(inout)              :: x2c_c1 !juanxiong he, for wrf/cam
    type(mct_aVect), intent(inout)              :: x2c_c2 !juanxiong he, for wrf/cam
    type(mct_aVect), intent(inout)              :: c2x_c1 !juanxiong he, for wrf/cam
    type(mct_aVect), intent(inout)              :: c2x_c2 !juanxiong he, for wrf/cam
    type(mct_aVect), intent(inout)              :: x2ca_caca1 !juanxiong he, for geatm/cam
    type(mct_aVect), intent(inout)              :: x2ca_caca2 !juanxiong he, for geatm/cam
    type(mct_aVect), intent(inout)              :: ca2x_caca1 !juanxiong he, for geatm/cam
    type(mct_aVect), intent(inout)              :: ca2x_caca2 !juanxiong he, for geatm/cam 
    character(len=*), optional,   intent(IN)    :: NLFilename ! Namelist filename
    logical, intent(inout) :: twoway_coupling        !jaunxiong he
    integer, intent(in) :: twoway_nudging        !jaunxiong he

    !
    ! Locals
    !
    type(mct_gsMap), pointer   :: gsMap_atm
    type(mct_gsMap), pointer   :: gsMap_cc  !juanxiong he for wrf/cam coupling
    type(mct_gsMap), pointer   :: gsMap_caca  !juanxiong he for geatm/cam coupling
    type(mct_gGrid), pointer   :: dom_a
    type(mct_gGrid), pointer   :: dom_c   !juanxiong he for wrf/cam coupling
    type(mct_gGrid), pointer   :: dom_ca   !juanxiong he for geatm/cam coupling
    type(seq_infodata_type),pointer :: infodata
    integer :: ATMID
    integer :: mpicom_atm
    integer :: lsize
    integer :: iradsw
    logical :: exists           ! true if file exists
    real(r8):: nextsw_cday      ! calendar of next atm shortwave
    integer :: stepno           ! time step			 
    integer :: dtime_sync       ! integer timestep size
    integer :: currentymd       ! current year-month-day
    integer :: dtime            ! time step increment (sec)
    integer :: atm_cpl_dt       ! driver atm coupling time step 
    integer :: nstep            ! CAM nstep
    real(r8):: caldayp1         ! CAM calendar day for for next cam time step
    integer :: dtime_cam        ! Time-step increment (sec)
    integer :: ymd              ! CAM current date (YYYYMMDD)
    integer :: yr               ! CAM current year
    integer :: mon              ! CAM current month
    integer :: day              ! CAM current day
    integer :: tod              ! CAM current time of day (sec)
    integer :: start_ymd        ! Start date (YYYYMMDD)
    integer :: start_tod        ! Start time of day (sec)
    integer :: ref_ymd          ! Reference date (YYYYMMDD)
    integer :: ref_tod          ! Reference time of day (sec)
    integer :: stop_ymd         ! Stop date (YYYYMMDD)
    integer :: stop_tod         ! Stop time of day (sec)
    logical :: perpetual_run    ! If in perpetual mode or not
    integer :: perpetual_ymd    ! Perpetual date (YYYYMMDD)
    logical :: single_column
    real(r8):: scmlat,scmlon
    integer :: shrlogunit,shrloglev ! old values
    logical :: first_time = .true.
    character(len=SHR_KIND_CS) :: calendar  ! Calendar type
    character(len=SHR_KIND_CS) :: starttype ! infodata start type
    integer :: lbnum
    integer :: hdim1_d, hdim2_d ! dimensions of rectangular horizontal grid
                                ! data structure, If 1D data structure, then
                                ! hdim2_d == 1.
    !-----------------------------------------------------------------------
    !
    ! Determine cdata points
    !
#if (defined _MEMTRACE)
    if(masterproc) then
      lbnum=1
      call memmon_dump_fort('memmon.out','atm_init_mct:start::',lbnum)
    endif                      
#endif                         
    call seq_cdata_setptrs(cdata_a, ID=ATMID, mpicom=mpicom_atm, &
         gsMap=gsMap_atm, dom=dom_a, infodata=infodata)

    call seq_cdata_setptrs(cdata_c, gsMap=gsMap_cc, dom=dom_c)  !juanxiong he

    call seq_cdata_setptrs(cdata_ca, gsMap=gsMap_caca, dom=dom_ca)  !juanxiong he
    if (first_time) then
       
       ! Determine attribute vector indices

       call cam_cpl_indices_set()

       ! Redirect share output to cam log
       
       call spmdinit(mpicom_atm)
       
       if (masterproc) then
          inquire(file='atm_modelio.nml',exist=exists)
          if (exists) then
             iulog = shr_file_getUnit()
             call shr_file_setIO('atm_modelio.nml',iulog)
          endif
          write(iulog,*) "CAM atmosphere model initialization"
       endif
       
       call shr_file_getLogUnit (shrlogunit)
       call shr_file_getLogLevel(shrloglev)
       call shr_file_setLogUnit (iulog)
       ! 
       ! Consistency check                              
       !
       if (co2_readFlux_ocn .and. index_x2a_Faxx_fco2_ocn /= 0) then
          write(iulog,*)'error co2_readFlux_ocn and index_x2a_Faxx_fco2_ocn cannot both be active'
          call shr_sys_abort()
       end if
       ! 
       ! Get data from infodata object
       !
       call seq_infodata_GetData( infodata,                                           &
            case_name=caseid, case_desc=ctitle,                                       &
            start_type=starttype,                                                     &
            atm_adiabatic=adiabatic,                                                  &
            atm_ideal_phys=ideal_phys,                                                &
            aqua_planet=aqua_planet,                                                  &
            brnch_retain_casename=brnch_retain_casename,                              &
            single_column=single_column, scmlat=scmlat, scmlon=scmlon,                &
            orb_eccen=eccen, orb_mvelpp=mvelpp, orb_lambm0=lambm0, orb_obliqr=obliqr, &
            lnd_present=lnd_present, ocn_present=ocn_present,                         & 
            perpetual=perpetual_run, perpetual_ymd=perpetual_ymd)
       !
       ! Get nsrest from startup type methods
       !
       if (     trim(starttype) == trim(seq_infodata_start_type_start)) then
          nsrest = 0
       else if (trim(starttype) == trim(seq_infodata_start_type_cont) ) then
          nsrest = 1
       else if (trim(starttype) == trim(seq_infodata_start_type_brnch)) then
          nsrest = 3
       else
          write(iulog,*) 'atm_comp_mct: ERROR: unknown starttype'
          call shr_sys_abort()
       end if
       !
       ! Initialize time manager.
       !
       call seq_timemgr_EClockGetData(EClock, &
                                      start_ymd=start_ymd, start_tod=start_tod, &
                                      ref_ymd=ref_ymd, ref_tod=ref_tod,         &
                                      stop_ymd=stop_ymd, stop_tod=stop_tod,     &
                                      calendar=calendar )
       !
       ! Read namelist
       !
       call read_namelist(single_column_in=single_column, scmlat_in=scmlat, scmlon_in=scmlon)
       !
       ! Initialize cam time manager
       !
       if ( nsrest == 0 )then
          call timemgr_init( calendar_in=calendar, start_ymd=start_ymd, &
                             start_tod=start_tod, ref_ymd=ref_ymd,      &
                             ref_tod=ref_tod, stop_ymd=stop_ymd,        &
                             stop_tod=stop_tod,                         &
                             perpetual_run=perpetual_run,               &
                             perpetual_ymd=perpetual_ymd )
       end if
       !
       ! First phase of cam initialization 
       ! Initialize mpicom_atm, allocate cam_in and cam_out and determine 
       ! atm decomposition (needed to initialize gsmap) 
       ! for an initial run, cam_in and cam_out are allocated in cam_initial
       ! for a restart/branch run, cam_in and cam_out are allocated in restart 
       ! Set defaults then override with user-specified input and initialize time manager
       ! Note that the following arguments are needed to cam_init for timemgr_restart only
       !
       call cam_init( cam_out, cam_in, cam_state, cam_tend, wrf_state, wrf_tend, mpicom_atm, &  !juanxiong he
                      start_ymd, start_tod, ref_ymd, ref_tod, stop_ymd, stop_tod, &
                      perpetual_run, perpetual_ymd, calendar)

       !
       ! Check consistency of restart time information with input clock
       !
       if (nsrest /= 0) then
          dtime_cam = get_step_size()
          call timemgr_check_restart( calendar, start_ymd, start_tod, ref_ymd, &
                                      ref_tod, dtime_cam, perpetual_run, perpetual_ymd)
       end if
       !
       ! Initialize MCT gsMap, domain and attribute vectors
       !
       call atm_SetgsMap_mct( mpicom_atm, ATMID, gsMap_atm )
       lsize = mct_gsMap_lsize(gsMap_atm, mpicom_atm)
       !
       ! Initialize MCT domain 
       !
       call atm_domain_mct( lsize, gsMap_atm, dom_a )
       !
       ! Initialize MCT attribute vectors
       !
       call mct_aVect_init(a2x_a, rList=seq_flds_a2x_fields, lsize=lsize)
       call mct_aVect_zero(a2x_a)
       
       call mct_aVect_init(x2a_a, rList=seq_flds_x2a_fields, lsize=lsize) 
       call mct_aVect_zero(x2a_a)
       
       call mct_aVect_init(a2x_a_SNAP, rList=a2x_avg_flds, lsize=lsize)
       call mct_aVect_zero(a2x_a_SNAP)
       
       call mct_aVect_init(a2x_a_SUM , rList=a2x_avg_flds, lsize=lsize)
       call mct_aVect_zero(a2x_a_SUM )

       call atm_SetgsMap_wrf( mpicom_atm, ATMID, gsMap_cc )  !juanxiong he
       lsize = mct_gsMap_lsize(gsMap_cc, mpicom_atm) !juanxiong he

       ! for wrf coupling
       call atm_domain_wrf( lsize, gsMap_cc, dom_c )  !juanxiong he
              
       call mct_aVect_init(c2x_c1, rList=seq_flds_c2x_fields, lsize=lsize) !juaxniong he
       call mct_aVect_zero(c2x_c1) !juanxiong he

       call mct_aVect_init(c2x_c2, rList=seq_flds_c2x_fields, lsize=lsize) !juaxniong he
       call mct_aVect_zero(c2x_c2) !juanxiong he

       call mct_aVect_init(x2c_c1, rList=seq_flds_x2c_fields, lsize=lsize) !juanxiong he
       call mct_aVect_zero(x2c_c1) !juanxiong he      

       call mct_aVect_init(x2c_c2, rList=seq_flds_x2c_fields, lsize=lsize) !juanxiong he
       call mct_aVect_zero(x2c_c2) !juanxiong he   

       call atm_SetgsMap_geatm( mpicom_atm, ATMID, gsMap_caca )  !juanxiong he
       lsize = mct_gsMap_lsize(gsMap_caca, mpicom_atm) !juanxiong he
       ! for geatm coupling
       call atm_domain_geatm( lsize, gsMap_caca, dom_ca )  !juanxiong he
              
       call mct_aVect_init(ca2x_caca1, rList=seq_flds_ca2x_fields, lsize=lsize) !juaxniong he
       call mct_aVect_zero(ca2x_caca1) !juanxiong he

       call mct_aVect_init(ca2x_caca2, rList=seq_flds_ca2x_fields, lsize=lsize) !juaxniong he
       call mct_aVect_zero(ca2x_caca2) !juanxiong he

       call mct_aVect_init(x2ca_caca1, rList=seq_flds_x2ca_fields, lsize=lsize) !juanxiong he
       call mct_aVect_zero(x2ca_caca1) !juanxiong he      

       call mct_aVect_init(x2ca_caca2, rList=seq_flds_x2ca_fields, lsize=lsize) !juanxiong he
       call mct_aVect_zero(x2ca_caca2) !juanxiong he   
       !
	   ! Initialize mgrid structure   by Huiqun Hao
	   !
	   call t_startf ('DC_init_mgrid')
       call set_metgrid(mgrid1, num_soil_layers, plev+1)
       call initial_metgrid(mgrid1, mgrid1%ids, mgrid1%ide, mgrid1%jds, mgrid1%jde, &
                            mgrid1%num_metgrid_levels, mgrid1%num_metgrid_soil_levels )

       call set_metgrid(mgrid2, num_soil_layers, plev+1)
       call initial_metgrid(mgrid2, mgrid2%ids, mgrid2%ide, mgrid2%jds, mgrid2%jde, &
                            mgrid2%num_metgrid_levels, mgrid2%num_metgrid_soil_levels )
	   call t_stopf ('DC_init_mgrid')
       !
       ! Initialize averaging counter
       !
       avg_count = 0
       !
       ! Create initial atm export state
       !
       call atm_export_mct( cam_out, a2x_a )
       !
       ! Set flag to specify that an extra albedo calculation is to be done (i.e. specify active)
       !
       call seq_infodata_PutData(infodata, atm_prognostic=.true.)
       call get_horiz_grid_dim_d(hdim1_d, hdim2_d)
       call seq_infodata_PutData(infodata, atm_nx=hdim1_d, atm_ny=hdim2_d)

       ! Set flag to indicate that CAM will provide carbon and dust deposition fluxes.
       ! This is now hardcoded to .true. since the ability of CICE to read these
       ! fluxes from a file has been removed.
       call seq_infodata_PutData(infodata, atm_aero=.true.)

       !
       ! Set time step of radiation computation as the current calday
       ! This will only be used on the first timestep of an initial run
       !
       if (nsrest == 0) then
          nextsw_cday = get_curr_calday()
          call seq_infodata_PutData( infodata, nextsw_cday=nextsw_cday )
       end if
       
       ! End redirection of share output to cam log
       
       call shr_file_setLogUnit (shrlogunit)
       call shr_file_setLogLevel(shrloglev)

       first_time = .false.

    else
       
       ! For initial run, run cam radiation/clouds and return
       ! For restart run, read restart x2a_a
       ! Note - a2x_a is computed upon the completion of the previous run - cam_run1 is called
       ! only for the purposes of finishing the flux averaged calculation to compute a2x_a
       ! Note - cam_run1 is called on restart only to have cam internal state consistent with the 
       ! a2x_a state sent to the coupler

       ! Redirect share output to cam log

       call shr_file_getLogUnit (shrlogunit)
       call shr_file_getLogLevel(shrloglev)
       call shr_file_setLogUnit (iulog)

       call seq_timemgr_EClockGetData(EClock,curr_ymd=CurrentYMD, StepNo=StepNo, dtime=DTime_Sync )
       if (StepNo == 0) then
          call atm_import_mct( x2a_a, cam_in )
           twoway_coupling=.False.
          call cam_run1 ( cam_in, cam_out, cam_state, cam_tend, wrf_state, wrf_tend, twoway_coupling, twoway_nudging )  !juanxiong he
          call atm_export_mct( cam_out, a2x_a )
          call atm_export_geatm( ca2x_caca1, cam_in, cam_out, cam_state, cam_tend ) !juanxiong he
          call atm_export_geatm( ca2x_caca2, cam_in, cam_out, cam_state, cam_tend ) !juanxiong he
       else
          call atm_read_srfrest_mct( EClock, cdata_a, x2a_a, a2x_a )
          call atm_read_wrfrest_mct( EClock, cdata_c, x2c_c1, c2x_c2 )
          call atm_read_gearest_mct( EClock, cdata_ca, x2ca_caca1, ca2x_caca2 )
          call atm_import_mct( x2a_a, cam_in )
          if(twoway_coupling) then
          call atm_import_wrf( x2c_c1, c2x_c2, wrf_state, wrf_tend, twoway_nudging ) !juanxiong he, wrf state restart
          call atm_import_geatm( x2ca_caca1, ca2x_caca2, wrf_state, wrf_tend, twoway_nudging ) !juanxiong he, wrf state restart          
          call atm_import_self( c2x_c2, cam_state, cam_tend ) ! juanxiong he, cam state before two way coupling
          endif
          call cam_run1 ( cam_in, cam_out, cam_state, cam_tend, wrf_state, wrf_tend, twoway_coupling, twoway_nudging )  !juanxiong he
       end if

       ! Compute time of next radiation computation, like in run method for exact restart

! tcx was
!       nextsw_cday = radiation_nextsw_cday() 

       call seq_timemgr_EClockGetData(Eclock,dtime=atm_cpl_dt)
       dtime = get_step_size()          
       nstep = get_nstep()
       if (nstep < 1 .or. dtime < atm_cpl_dt) then
          nextsw_cday = radiation_nextsw_cday() 
       else if (dtime == atm_cpl_dt) then
          caldayp1 = get_curr_calday(offset=int(dtime))
          nextsw_cday = radiation_nextsw_cday() 
          if (caldayp1 /= nextsw_cday) nextsw_cday = -1._r8
       else
          call shr_sys_abort('dtime must be less than or equal to atm_cpl_dt')
       end if
       call seq_infodata_PutData( infodata, nextsw_cday=nextsw_cday ) 

       ! End redirection of share output to cam log
       
       call shr_file_setLogUnit (shrlogunit)
       call shr_file_setLogLevel(shrloglev)
       
    end if

#if (defined _MEMTRACE )
    if(masterproc) then
      lbnum=1
      call memmon_dump_fort('memmon.out','atm_init_mct:end::',lbnum)
      call memmon_reset_addr()
    endif
#endif

    call cnst_get_ind ( 'ISOP', id_isop, abort=.false.)
    call cnst_get_ind ( 'C10H16', id_c10h16, abort=.false.)

    call shr_sys_flush(iulog)

 end subroutine atm_init_mct

!================================================================================

  subroutine atm_run_mct( EClock, EClock_w, cdata_a, cdata_c, cdata_ca, &
                          x2a_a, a2x_a, x2c_c1, x2c_c2, c2x_c1, c2x_c2, &
                          x2ca_caca1, x2ca_caca2, ca2x_caca1, ca2x_caca2, & ! juanxiong he for geatm 
                          twoway_coupling, twoway_nudging, integration_phase) !juanxiong he
    !-----------------------------------------------------------------------
    !
    ! Uses
    !
    use time_manager,    only: advance_timestep, get_curr_date, get_curr_calday, &
	                       get_nstep, get_step_size
    use scamMod,         only: single_column
!   use iop,             only: scam_use_iop_srf
    use pmgrid,          only: plev, plevp
    use constituents,    only: pcnst
    use shr_sys_mod, only: shr_sys_flush
    ! 
    ! Arguments
    !
    type(ESMF_Clock)            ,intent(in)    :: EClock
    type(ESMF_Clock)            ,intent(in)    :: EClock_w ! juanxiong he
    type(seq_cdata)             ,intent(inout) :: cdata_a
    type(seq_cdata)             ,intent(inout) :: cdata_c !juanxiong he
    type(seq_cdata)             ,intent(inout) :: cdata_ca !juanxiong he
    type(mct_aVect)             ,intent(inout) :: x2a_a
    type(mct_aVect)             ,intent(inout) :: a2x_a
    type(mct_aVect)             ,intent(inout) :: x2c_c1  !juanxiong he
    type(mct_aVect)             ,intent(inout) :: x2c_c2  !juanxiong he
    type(mct_aVect)             ,intent(inout) :: c2x_c1  !juanxiong he
    type(mct_aVect)             ,intent(inout) :: c2x_c2  !juanxiong he
    type(mct_aVect), intent(inout)              :: x2ca_caca1 !juanxiong he
    type(mct_aVect), intent(inout)              :: x2ca_caca2 !juanxiong he
    type(mct_aVect), intent(inout)              :: ca2x_caca1 !juanxiong he  
    type(mct_aVect), intent(inout)              :: ca2x_caca2 !juanxiong he     
    logical, intent(inout) :: twoway_coupling        !jaunxiong he
    integer, intent(in) :: twoway_nudging        !jaunxiong he
    integer, intent(inout) :: integration_phase !juanxiong he

    !
    ! Local variables
    !
    type(seq_infodata_type),pointer :: infodata
    integer :: lsize           ! size of attribute vector
    integer :: StepNo          ! time step			 
    integer :: DTime_Sync      ! integer timestep size
    integer :: CurrentYMD      ! current year-month-day
    integer :: iradsw          ! shortwave radation frequency (time steps) 
    logical :: dosend          ! true => send data back to driver
    integer :: dtime           ! time step increment (sec)
    integer :: atm_cpl_dt      ! driver atm coupling time step 
    integer :: wrf_cpl_dt      ! driver wrf coupling time step, juanxiong he
    integer :: ymd_sync        ! Sync date (YYYYMMDD)
    integer :: yr_sync         ! Sync current year
    integer :: mon_sync        ! Sync current month
    integer :: day_sync        ! Sync current day
    integer :: tod_sync        ! Sync current time of day (sec)
    integer :: ymd             ! CAM current date (YYYYMMDD)
    integer :: yr              ! CAM current year
    integer :: mon             ! CAM current month
    integer :: day             ! CAM current day
    integer :: tod             ! CAM current time of day (sec)
    integer :: nstep           ! CAM nstep
    integer :: shrlogunit,shrloglev ! old values
    real(r8):: caldayp1        ! CAM calendar day for for next cam time step
    real(r8):: nextsw_cday     ! calendar of next atm shortwave
    logical :: rstwr           ! .true. ==> write restart file before returning
    logical :: nlend           ! Flag signaling last time-step
    logical :: rstwr_sync      ! .true. ==> write restart file before returning
    logical :: nlend_sync      ! Flag signaling last time-step
    logical :: first_time = .true.
    character(len=*), parameter :: subname="atm_run_mct"
    !-----------------------------------------------------------------------
    integer :: lbnum
    integer :: loop  !by Wang Yuzhu

#if (defined _MEMTRACE)
    if(masterproc) then
      lbnum=1
      call memmon_dump_fort('memmon.out',SubName //':start::',lbnum)
    endif
#endif

    ! Redirect share output to cam log
    
    call shr_file_getLogUnit (shrlogunit)
    call shr_file_getLogLevel(shrloglev)
    call shr_file_setLogUnit (iulog)
    
    !----------------------------------------------------------------------- 
    ! phase =1, don't send to coupler, send to wrf only
    ! phase =2, send to coupler only
    ! juanxiong he 
    !-----------------------------------------------------------------------

    !----------------------------------------------------------------------- 
    ! before the integration during the coupling interval, juanxiong he 
    !----------------------------------------------------------------------- 
    if(integration_phase.eq.1) then
    ! Note that sync clock time should match cam time at end of time step/loop not beginning
    
    call seq_cdata_setptrs(cdata_a, infodata=infodata)
    call seq_timemgr_EClockGetData(EClock,curr_ymd=ymd_sync,curr_tod=tod_sync, &
       curr_yr=yr_sync,curr_mon=mon_sync,curr_day=day_sync)

    !load orbital parameters
    call seq_infodata_GetData( infodata,                                           &
       orb_eccen=eccen, orb_mvelpp=mvelpp, orb_lambm0=lambm0, orb_obliqr=obliqr)

    nlend_sync = seq_timemgr_StopAlarmIsOn(EClock)
    rstwr_sync = seq_timemgr_RestartAlarmIsOn(EClock)
    
    ! Map input from mct to cam data structure

    call t_startf ('CAM_import')
    call atm_import_mct( x2a_a, cam_in )
    call t_stopf  ('CAM_import')
    
    endif

    ! Cycle over all time steps in the atm coupling interval
    
    if(integration_phase.eq.1) then  ! juanxiong he
     call seq_timemgr_EClockGetData(Eclock,dtime=atm_cpl_dt)
     call seq_timemgr_EClockGetData(Eclock_w,dtime=wrf_cpl_dt)
     dosend = .false.
    else
     dosend = .true.
    endif

    do while (.not. dosend)
       
       ! Determine if dosend
       ! When time is not updated at the beginning of the loop - then return only if
       ! are in sync with clock before time is updated
       
       call get_curr_date( yr, mon, day, tod )
       ymd = yr*10000 + mon*100 + day
       tod = tod
       dosend = (seq_timemgr_EClockDateInSync( EClock, ymd, tod))
       
       ! Determine if time to write cam restart and stop
       
       rstwr = .false.
       if (rstwr_sync .and. dosend) rstwr = .true.
       nlend = .false.
       if (nlend_sync .and. dosend) nlend = .true.
       
       ! Single column specific input 
       
       if (single_column) then
          call scam_use_iop_srf( cam_in )
       endif
       ! by Wang Yuzhu
              
       loop = 1
       ! Run CAM (run2, run3, run4)
       ! juanxiong he
       nstep = get_nstep()

       if(.not.twoway_coupling) then
       if (nstep.eq.0) then
!         call atm_export_wrf( c2x_c1, cam_in, cam_state, cam_tend ) !juanxiong he !jjr
         call t_startf ('cam_export_wrf1_dc')
		 call wrf_import_cam_test( mgrid1, loop, cam_in, cam_state )   !by Huiqun Hao
		 mgrid2%sstold = mgrid1%sstold
		 mgrid2%xiceold = mgrid1%xiceold
		 call t_stopf ('cam_export_wrf1_dc')
       end if
       if(atm_cpl_dt/wrf_cpl_dt.gt.1) then
        if (nstep.gt.2) then
		 call t_startf ('cam_export_wrf2_dc')
!         call mct_aVect_copy( c2x_c1, c2x_c2 ) ! juanxiong he  !jjr
		 mgrid2=mgrid1
!         call atm_export_wrf( c2x_c1, cam_in, cam_state, cam_tend ) !juanxiong he !jjr
		 call wrf_import_cam_test( mgrid1, loop, cam_in, cam_state )   !by Huiqun Hao
		 call t_stopf ('cam_export_wrf2_dc')
        end if
       else
        if (nstep.gt.1) then
!         call mct_aVect_copy( c2x_c2, c2x_c1 ) ! juanxiong he !jjr
!                 mgrid1=mgrid2   !by Huiqun Hao
        end if
!         call atm_export_wrf( c2x_c2, cam_in, cam_state, cam_tend ) !juanxiong he !jjr
         call t_startf ('cam_export_wrf3_dc')
		 call wrf_import_cam_test( mgrid2, loop, cam_in, cam_state )   !by Huiqun Hao
		 call t_stopf ('cam_export_wrf3_dc')
       end if
       end if

       call t_startf ('CAM_run2')
       call cam_run2( cam_out, cam_in, cam_state, cam_tend, &    !juanxiong he 
                      wrf_state, wrf_tend, twoway_coupling, twoway_nudging ) 
       call t_stopf  ('CAM_run2')
		
         call atm_export_geatm( ca2x_caca1, cam_in, cam_out, cam_state, cam_tend ) !juanxiong he
         call atm_export_geatm( ca2x_caca2, cam_in, cam_out, cam_state, cam_tend ) !juanxiong he

       call t_startf ('CAM_run3')
       call cam_run3( cam_out, cam_state, cam_tend, wrf_state, wrf_tend, &   !juanxiong he
                      twoway_coupling, twoway_nudging)
       call t_stopf  ('CAM_run3')
       
       call t_startf ('CAM_run4')
       call cam_run4( cam_out, cam_in, rstwr, nlend, &
            yr_spec=yr_sync, mon_spec=mon_sync, day_spec=day_sync, sec_spec=tod_sync)
       call t_stopf  ('CAM_run4')
       
       if(twoway_coupling) then
       if(nstep.eq.0) then
!       call atm_export_wrf( c2x_c1, cam_in, cam_state, cam_tend ) !juanxiong he
	   call wrf_import_cam_test( mgrid1, loop, cam_in, cam_state )   !by Huiqun Hao
       call atm_export_geatm( ca2x_caca1, cam_in, cam_out, cam_state, cam_tend ) !juanxiong he
       endif
!       call atm_export_wrf( c2x_c2, cam_in, cam_state, cam_tend ) !juanxiong he
	   call wrf_import_cam_test( mgrid2, loop, cam_in, cam_state )   !by Huiqun Hao
       call atm_export_geatm( ca2x_caca2, cam_in, cam_out, cam_state, cam_tend ) !juanxiong he
       endif

       ! Advance cam time step 
       
       call t_startf ('CAM_adv_timestep')
       call advance_timestep()
       call t_stopf  ('CAM_adv_timestep')
       
       ! Run cam radiation/clouds (run1)
          
       if(.not.dosend) then  ! juanxiong he

       call t_startf ('CAM_run1')
       call cam_run1 ( cam_in, cam_out, cam_state, cam_tend, wrf_state, wrf_tend, twoway_coupling, twoway_nudging )  !juanxiong he
       call t_stopf  ('CAM_run1')
       
       ! Map output from cam to mct data structures
       
       call t_startf ('CAM_export')
       call atm_export_mct( cam_out, a2x_a )
       call t_stopf ('CAM_export')
       
       ! Compute snapshot attribute vector for accumulation
       
! don't accumulate on first coupling freq ts1 and ts2
! for consistency with ccsm3 when flxave is off
       nstep = get_nstep()
       if (nstep <= 2) then
          call mct_aVect_copy( a2x_a, a2x_a_SUM )
          avg_count = 1
       else
          call mct_aVect_copy( a2x_a, a2x_a_SNAP )
          call mct_aVect_accum( aVin=a2x_a_SNAP, aVout=a2x_a_SUM )
          avg_count = avg_count + 1
       endif

      end if ! juanxiong he
       
    end do

    !----------------------------------------------------------------------- 
    ! Finish the final step of the coupling interval, juanxiong he
    !----------------------------------------------------------------------- 

    if (integration_phase.eq.2) then  ! juanxiong he

       call t_startf ('CAM_import')
        if(twoway_coupling) then
           call atm_import_wrf( x2c_c1, c2x_c2, wrf_state, wrf_tend, twoway_nudging) !juanxiong he
           call atm_import_geatm( x2ca_caca1, ca2x_caca2, wrf_state, wrf_tend, twoway_nudging) !juanxiong he
        endif
       call t_stopf  ('CAM_import')

       !----------------------------------------------------------------------- 
       ! Finish the cam_run1 at the end of the coupling interval, juanxiong he
       !-----------------------------------------------------------------------

       call t_startf ('CAM_run1')
       call cam_run1 ( cam_in, cam_out, cam_state, cam_tend, wrf_state, wrf_tend, twoway_coupling, twoway_nudging )  !juanxiong he
       call t_stopf  ('CAM_run1')

       ! Map output from cam to mct data structures
        
       call t_startf ('CAM_export')
       call atm_export_mct( cam_out, a2x_a )
       call t_stopf ('CAM_export')
        
       ! Compute snapshot attribute vector for accumulation
        
! don't accumulate on first coupling freq ts1 and ts2
! for consistency with ccsm3 when flxave is off
       nstep = get_nstep()
       if (nstep <= 2) then
          call mct_aVect_copy( a2x_a, a2x_a_SUM )
          avg_count = 1
       else 
          call mct_aVect_copy( a2x_a, a2x_a_SNAP )
          call mct_aVect_accum( aVin=a2x_a_SNAP, aVout=a2x_a_SUM )
          avg_count = avg_count + 1
       endif

       !----------------------------------------------------------------------- 
       ! Finish cam_run1 at the end of the coupling interval, juanxiong he 
       !-----------------------------------------------------------------------

    ! Finish accumulation of attribute vector and average and copy accumulation 
    ! field into output attribute vector
    
    call mct_aVect_avg ( a2x_a_SUM, avg_count)
    call mct_aVect_copy( a2x_a_SUM, a2x_a )
    call mct_aVect_zero( a2x_a_SUM) 
    avg_count = 0                   
    
    ! Get time of next radiation calculation - albedos will need to be 
    ! calculated by each surface model at this time
    
    call seq_timemgr_EClockGetData(Eclock,dtime=atm_cpl_dt)
    dtime = get_step_size()          
    if (dtime < atm_cpl_dt) then
       nextsw_cday = radiation_nextsw_cday() 
    else if (dtime == atm_cpl_dt) then
       caldayp1 = get_curr_calday(offset=int(dtime))
       nextsw_cday = radiation_nextsw_cday() 
       if (caldayp1 /= nextsw_cday) nextsw_cday = -1._r8
    else
       call shr_sys_abort('dtime must be less than or equal to atm_cpl_dt')
    end if
    call seq_infodata_PutData( infodata, nextsw_cday=nextsw_cday ) 
    
    ! Write merged surface data restart file if appropriate
    rstwr_sync = seq_timemgr_RestartAlarmIsOn(EClock) ! rstwr_sync doesn't have the save attribute,juanxiong he
    if (rstwr_sync) then
       call seq_timemgr_EClockGetData(EClock,curr_ymd=ymd_sync,curr_tod=tod_sync, &
            curr_yr=yr_sync,curr_mon=mon_sync,curr_day=day_sync) ! get the time again,since it's seperated phase. addedd by juanxiong he
       call atm_write_srfrest_mct( cdata_a, x2a_a, a2x_a, &
            yr_spec=yr_sync, mon_spec=mon_sync, day_spec=day_sync, sec_spec=tod_sync)
       call atm_write_wrfrest_mct( cdata_c, x2c_c1, c2x_c2, & 
            yr_spec=yr_sync, mon_spec=mon_sync, day_spec=day_sync, sec_spec=tod_sync) ! write wrf restart, added by juanxiong he
       call atm_write_gearest_mct( cdata_c, x2c_c1, c2x_c2, & 
            yr_spec=yr_sync, mon_spec=mon_sync, day_spec=day_sync, sec_spec=tod_sync) ! write wrf restart, added by juanxiong he             
    end if
    
    ! Check for consistency of internal cam clock with master sync clock 
    
    dtime = get_step_size()
    call get_curr_date( yr, mon, day, tod, offset=-dtime )
    ymd = yr*10000 + mon*100 + day
    tod = tod
    if ( .not. seq_timemgr_EClockDateInSync( EClock, ymd, tod ) )then
       call seq_timemgr_EClockGetData(EClock, curr_ymd=ymd_sync, curr_tod=tod_sync )
       write(iulog,*)' cam ymd=',ymd     ,'  cam tod= ',tod
       write(iulog,*)'sync ymd=',ymd_sync,' sync tod= ',tod_sync
       call shr_sys_abort( subname//': CAM clock is not in sync with master Sync Clock' )
    end if
    
    endif
    !----------------------------------------------------------------------- 
    ! Finish the final step of the coupling interval, juanxiong he 
    !----------------------------------------------------------------------

    ! End redirection of share output to cam log

    call shr_file_setLogUnit (shrlogunit)
    call shr_file_setLogLevel(shrloglev)

#if (defined _MEMTRACE)
    if(masterproc) then
      lbnum=1
      call memmon_dump_fort('memmon.out',SubName //':end::',lbnum)
      call memmon_reset_addr()
    endif
#endif

  end subroutine atm_run_mct

!================================================================================

  subroutine atm_final_mct( )

     call cam_final( cam_out, cam_in, cam_state, cam_tend, wrf_state, wrf_tend )  ! added by juanxiong he

  end subroutine atm_final_mct

!================================================================================

  subroutine atm_SetgsMap_mct( mpicom_atm, ATMID, GSMap_atm )
    use phys_grid, only : get_nlcols_p
    !-------------------------------------------------------------------
    !
    ! Arguments
    !
    integer        , intent(in)  :: mpicom_atm
    integer        , intent(in)  :: ATMID
    type(mct_gsMap), intent(out) :: GSMap_atm
    !
    ! Local variables
    !
    integer, allocatable :: gindex(:)
    integer :: i, n, c, ncols, sizebuf, nlcols
    integer :: ier            ! error status
    !-------------------------------------------------------------------

    ! Build the atmosphere grid numbering for MCT
    ! NOTE:  Numbering scheme is: West to East and South to North
    ! starting at south pole.  Should be the same as what's used in SCRIP
    
    ! Determine global seg map

    sizebuf=0
    do c = begchunk, endchunk
       ncols = get_ncols_p(c)
       do i = 1,ncols
          sizebuf = sizebuf+1
       end do
    end do

    allocate(gindex(sizebuf))

    n=0
    do c = begchunk, endchunk
       ncols = get_ncols_p(c)
       do i = 1,ncols
          n=n+1
          gindex(n) = get_gcol_p(c,i)
       end do
    end do

    nlcols = get_nlcols_p()
    call mct_gsMap_init( gsMap_atm, gindex, mpicom_atm, ATMID, nlcols, ngcols)

    deallocate(gindex)

  end subroutine atm_SetgsMap_mct

!===============================================================================

  subroutine atm_import_mct( x2a_a, cam_in )

    !-----------------------------------------------------------------------
    !
    ! Uses	
    !
    use dust_intr,     only: dust_idx1
#if (defined MODAL_AERO)
    use mo_chem_utls,  only: get_spc_ndx
#endif
    use shr_const_mod, only: shr_const_stebol
    use seq_drydep_mod,only: n_drydep
    use cam_nudging         !juanxiong he
    use dycore,           only: dycore_is  !juanxiong he
    !
    ! Arguments
    !
    type(mct_aVect),    intent(inout) :: x2a_a
    type(cam_in_t),     intent(inout) :: cam_in(begchunk:endchunk)
    !
    ! Local variables
    !		
    integer  :: i,lat,n,c,ig,k  ! indices, juanxiong he
    integer  :: ncols         ! number of columns
    integer  :: dust_ndx
    integer  :: nstep ! juanxiong he
    logical, save :: first_time = .true.

    ! factors used to convert MEGAN source units [ug C/m2/hr] to CAM srf emis units [kg/m2/sec]
    real(r8), parameter :: isop_factor   = 1._r8/3600._r8/1.e9_r8* 68._r8/ 60._r8
    real(r8), parameter :: c10h16_factor = 1._r8/3600._r8/1.e9_r8*136._r8/120._r8

#if (defined MODAL_AERO)
    integer, parameter:: ndst =2
    integer, target   :: spc_ndx(ndst)
#if (defined MODAL_AERO_7MODE)
    integer, pointer  :: dst_a5_ndx, dst_a7_ndx
#elif (defined MODAL_AERO_3MODE)
    integer, pointer  :: dst_a1_ndx, dst_a3_ndx
#endif
#endif
    !-----------------------------------------------------------------------
    !
#if (defined MODAL_AERO)
#if (defined MODAL_AERO_7MODE)
    dst_a5_ndx => spc_ndx(1)
    dst_a7_ndx => spc_ndx(2)
    dst_a5_ndx = get_spc_ndx( 'dst_a5' )
    dst_a7_ndx = get_spc_ndx( 'dst_a7' )
#elif (defined MODAL_AERO_3MODE)
    dst_a1_ndx => spc_ndx(1)
    dst_a3_ndx => spc_ndx(2)
    dst_a1_ndx = get_spc_ndx( 'dst_a1' )
    dst_a3_ndx = get_spc_ndx( 'dst_a3' )
#endif
#endif

    ! ccsm sign convention is that fluxes are positive downward

    ig=1
    do c=begchunk,endchunk
       ncols = get_ncols_p(c)                                                 
       do i =1,ncols                                                               
          cam_in(c)%wsx(i)       = -x2a_a%rAttr(index_x2a_Faxx_taux,ig)     
          cam_in(c)%wsy(i)       = -x2a_a%rAttr(index_x2a_Faxx_tauy,ig)     
          cam_in(c)%lhf(i)       = -x2a_a%rAttr(index_x2a_Faxx_lat, ig)     
          cam_in(c)%shf(i)       = -x2a_a%rAttr(index_x2a_Faxx_sen, ig)     
          cam_in(c)%lwup(i)      = -x2a_a%rAttr(index_x2a_Faxx_lwup,ig)    
          cam_in(c)%cflx(i,1)    = -x2a_a%rAttr(index_x2a_Faxx_evap,ig)                
          cam_in(c)%asdir(i)     =  x2a_a%rAttr(index_x2a_Sx_avsdr, ig)  
          cam_in(c)%aldir(i)     =  x2a_a%rAttr(index_x2a_Sx_anidr, ig)  
          cam_in(c)%asdif(i)     =  x2a_a%rAttr(index_x2a_Sx_avsdf, ig)  
          cam_in(c)%aldif(i)     =  x2a_a%rAttr(index_x2a_Sx_anidf, ig)
          cam_in(c)%ts(i)        =  x2a_a%rAttr(index_x2a_Sx_t,     ig)  
          cam_in(c)%sst(i)       =  x2a_a%rAttr(index_x2a_So_t,     ig)             
          cam_in(c)%snowhland(i) =  x2a_a%rAttr(index_x2a_Sl_snowh, ig)  
          cam_in(c)%snowhice(i)  =  x2a_a%rAttr(index_x2a_Si_snowh, ig)  
          cam_in(c)%tref(i)      =  x2a_a%rAttr(index_x2a_Sx_tref,  ig)  
          cam_in(c)%qref(i)      =  x2a_a%rAttr(index_x2a_Sx_qref,  ig)
          cam_in(c)%u10(i)       =  x2a_a%rAttr(index_x2a_Sx_u10,   ig)
          cam_in(c)%icefrac(i)   =  x2a_a%rAttr(index_x2a_Sa_ifrac, ig)  
          cam_in(c)%ocnfrac(i)   =  x2a_a%rAttr(index_x2a_Sa_ofrac, ig)
	  cam_in(c)%landfrac(i)  =  x2a_a%rAttr(index_x2a_Sa_lfrac, ig)
          if ( associated(cam_in(c)%ram1) ) &
               cam_in(c)%ram1(i) =  x2a_a%rAttr(index_x2a_Sl_ram1 , ig)
          if ( associated(cam_in(c)%fv) ) &
               cam_in(c)%fv(i)   =  x2a_a%rAttr(index_x2a_Sl_fv   , ig)
          dust_ndx = dust_idx1()
          ! check that dust constituents are actually in the simulation
          if (dust_ndx>0) then
#if (defined MODAL_AERO)
#if (defined MODAL_AERO_7MODE)
            cam_in(c)%cflx(i,dust_ndx   )  = 0.13_r8  &  ! 1st mode, based on Zender et al (2003) Table 1
#elif (defined MODAL_AERO_3MODE)
            cam_in(c)%cflx(i,dust_ndx   )  = 0.032_r8  &  ! 1st mode, based on Zender et al (2003) Table 1
#endif
                                           * (-x2a_a%rAttr(index_x2a_Fall_flxdst1, ig) &
                                              -x2a_a%rAttr(index_x2a_Fall_flxdst2, ig) &
                                              -x2a_a%rAttr(index_x2a_Fall_flxdst3, ig) &
                                              -x2a_a%rAttr(index_x2a_Fall_flxdst4, ig))
#if (defined MODAL_AERO_7MODE)
            cam_in(c)%cflx(i,dust_ndx-spc_ndx(1)+spc_ndx(2))  = 0.87_r8 &  ! 2nd mode
#elif (defined MODAL_AERO_3MODE)
            cam_in(c)%cflx(i,dust_ndx-spc_ndx(1)+spc_ndx(2))  = 0.968_r8 &  ! 2nd mode
#endif
                                           * (-x2a_a%rAttr(index_x2a_Fall_flxdst1, ig) &
                                              -x2a_a%rAttr(index_x2a_Fall_flxdst2, ig) &
                                              -x2a_a%rAttr(index_x2a_Fall_flxdst3, ig) &
                                              -x2a_a%rAttr(index_x2a_Fall_flxdst4, ig))
#else
	    cam_in(c)%cflx(i,dust_ndx   )  = -x2a_a%rAttr(index_x2a_Fall_flxdst1, ig)
	    cam_in(c)%cflx(i,dust_ndx +1)  = -x2a_a%rAttr(index_x2a_Fall_flxdst2, ig)
	    cam_in(c)%cflx(i,dust_ndx +2)  = -x2a_a%rAttr(index_x2a_Fall_flxdst3, ig)
	    cam_in(c)%cflx(i,dust_ndx +3)  = -x2a_a%rAttr(index_x2a_Fall_flxdst4, ig)
#endif
          endif
          ! heald: add MEGAN VOC fluxes (convert from units C to units X)
          if ( id_isop > 0 .and. index_x2a_Faxx_flxvoc1 > 0 ) then
            cam_in(c)%cflx(i,id_isop)   = -x2a_a%rAttr(index_x2a_Faxx_flxvoc1, ig)*isop_factor
          endif
          if ( id_c10h16 > 0 .and. index_x2a_Faxx_flxvoc2 > 0 ) then
            cam_in(c)%cflx(i,id_c10h16) = -x2a_a%rAttr(index_x2a_Faxx_flxvoc2, ig)*c10h16_factor
          endif
          if ( index_x2a_Sx_ddvel/=0 .and. n_drydep>0 ) then
             cam_in(c)%depvel(i,:n_drydep) = &
                  x2a_a%rAttr(index_x2a_Sx_ddvel:index_x2a_Sx_ddvel+n_drydep-1, ig)
          endif
          !
          ! fields needed to calculate water isotopes to ocean evaporation processes
          !
          cam_in(c)%ustar(i) = x2a_a%rAttr(index_x2a_So_ustar,ig)
          cam_in(c)%re(i)    = x2a_a%rAttr(index_x2a_So_re   ,ig)
          cam_in(c)%ssq(i)   = x2a_a%rAttr(index_x2a_So_ssq  ,ig)
          !
          ! bgc scenarios
          !
          if (index_x2a_Faxx_fco2_lnd /= 0) then
             cam_in(c)%fco2_lnd(i) = -x2a_a%rAttr(index_x2a_Faxx_fco2_lnd,ig)
          end if
          if (index_x2a_Faxx_fco2_ocn /= 0) then
             cam_in(c)%fco2_ocn(i) = -x2a_a%rAttr(index_x2a_Faxx_fco2_ocn,ig)
          end if
          if (index_x2a_Faxx_fdms_ocn /= 0) then
             cam_in(c)%fdms(i)     = -x2a_a%rAttr(index_x2a_Faxx_fdms_ocn,ig)
          end if

!--------------------------------------------------------------------------------------
! soil depth/height and soil temperature/moisture for wrf/cam coupling, added juanxiong he 
!-------------------------------------------------------------------------------------- 
       do k = 1, 4
       cam_in(c)%soildepth(i,k) = x2a_a%rAttr(index_x2a_Sx_soildepth(k),ig)
       cam_in(c)%soilthick(i,k) = x2a_a%rAttr(index_x2a_Sx_soilthick(k),ig)
       cam_in(c)%soilt(i,k) = x2a_a%rAttr(index_x2a_Sx_soilt(k),ig)
       cam_in(c)%soilm(i,k) = x2a_a%rAttr(index_x2a_Sx_soilm(k),ig)
       enddo
!--------------------------------------------------------------------------------------
! soil depth/height and soil temperature/moisture for wrf/cam coupling, added juanxiong he 
!--------------------------------------------------------------------------------------  
          ig=ig+1

       end do
    end do

    ! Get total co2 flux from components,
    ! Note - co2_transport determines if cam_in(c)%cflx(i,c_i(1:4)) is allocated

    if (co2_transport()) then

       ! Interpolate in time for flux data read in
       if (co2_readFlux_ocn) then
          call co2_time_interp_ocn
       end if
       if (co2_readFlux_fuel) then
          call co2_time_interp_fuel
       end if
       
       ! from ocn : data read in or from coupler or zero
       ! from fuel: data read in or zero
       ! from lnd : through coupler or zero
       do c=begchunk,endchunk
          ncols = get_ncols_p(c)                                                 
          do i=1,ncols                                                               
             
             ! all co2 fluxes in unit kgCO2/m2/s ! co2 flux from ocn 
             if (index_x2a_Faxx_fco2_ocn /= 0) then
                cam_in(c)%cflx(i,c_i(1)) = cam_in(c)%fco2_ocn(i)
             else if (co2_readFlux_ocn) then 
                ! convert from molesCO2/m2/s to kgCO2/m2/s
                cam_in(c)%cflx(i,c_i(1)) = &
                     -data_flux_ocn%co2flx(i,c)*(1._r8- cam_in(c)%landfrac(i)) &
                     *mwco2*1.0e-3_r8
             else
                cam_in(c)%cflx(i,c_i(1)) = 0._r8
             end if
             
             ! co2 flux from fossil fuel
             if (co2_readFlux_fuel) then
                cam_in(c)%cflx(i,c_i(2)) = data_flux_fuel%co2flx(i,c)
             else
                cam_in(c)%cflx(i,c_i(2)) = 0._r8
             end if
             
             ! co2 flux from land (cpl already multiplies flux by land fraction)
             if (index_x2a_Faxx_fco2_lnd /= 0) then
                cam_in(c)%cflx(i,c_i(3)) = cam_in(c)%fco2_lnd(i)
             else
                cam_in(c)%cflx(i,c_i(3)) = 0._r8
             end if
             
             ! merged co2 flux
             cam_in(c)%cflx(i,c_i(4)) = cam_in(c)%cflx(i,c_i(1)) + &
                                        cam_in(c)%cflx(i,c_i(2)) + &
                                        cam_in(c)%cflx(i,c_i(3))
          end do
       end do
    end if
    !
    ! if first step, determine longwave up flux from the surface temperature 
    !
    if (first_time) then
       if (is_first_step()) then
          do c=begchunk, endchunk
             ncols = get_ncols_p(c)
             do i=1,ncols
                cam_in(c)%lwup(i) = shr_const_stebol*(cam_in(c)%ts(i)**4)
             end do
          end do
       end if
       first_time = .false.
    end if

    if (fdda_sfc.gt.0) then
    select case(data_sfcsrc)
     case (2)
       nstep = get_nstep()
       call cam_read_ncep2_sfc( cam_in, cam_state, nstep )
     case default
    end select
    end if

  end subroutine atm_import_mct

!===============================================================================

  subroutine atm_export_mct( cam_out, a2x_a )

    !-------------------------------------------------------------------
    !
    ! Arguments
    !
    type(cam_out_t),     intent(in)  :: cam_out(begchunk:endchunk) 
    type(mct_aVect),     intent(out) :: a2x_a
    !
    ! Local variables
    !
    integer :: avsize, avnat
    integer :: i,m,c,n,ig       ! indices
    integer :: ncols            ! Number of columns
    !-----------------------------------------------------------------------

    ! Copy from component arrays into chunk array data structure
    ! Rearrange data from chunk structure into lat-lon buffer and subsequently
    ! create attribute vector

    ig=1
    do c=begchunk, endchunk
       ncols = get_ncols_p(c)
       do i=1,ncols
          a2x_a%rAttr(index_a2x_Sa_pslv   ,ig) = cam_out(c)%psl(i)
          a2x_a%rAttr(index_a2x_Sa_z      ,ig) = cam_out(c)%zbot(i)   
          a2x_a%rAttr(index_a2x_Sa_u      ,ig) = cam_out(c)%ubot(i)   
          a2x_a%rAttr(index_a2x_Sa_v      ,ig) = cam_out(c)%vbot(i)   
          a2x_a%rAttr(index_a2x_Sa_tbot   ,ig) = cam_out(c)%tbot(i)   
          a2x_a%rAttr(index_a2x_Sa_ptem   ,ig) = cam_out(c)%thbot(i)  
          a2x_a%rAttr(index_a2x_Sa_pbot   ,ig) = cam_out(c)%pbot(i)   
          a2x_a%rAttr(index_a2x_Sa_shum   ,ig) = cam_out(c)%qbot(i,1) 
	  a2x_a%rAttr(index_a2x_Sa_dens   ,ig) = cam_out(c)%rho(i)
          a2x_a%rAttr(index_a2x_Faxa_swnet,ig) = cam_out(c)%netsw(i)      
          a2x_a%rAttr(index_a2x_Faxa_lwdn ,ig) = cam_out(c)%flwds(i)  
          a2x_a%rAttr(index_a2x_Faxa_rainc,ig) = (cam_out(c)%precc(i)-cam_out(c)%precsc(i))*1000._r8
          a2x_a%rAttr(index_a2x_Faxa_rainl,ig) = (cam_out(c)%precl(i)-cam_out(c)%precsl(i))*1000._r8
          a2x_a%rAttr(index_a2x_Faxa_snowc,ig) = cam_out(c)%precsc(i)*1000._r8
          a2x_a%rAttr(index_a2x_Faxa_snowl,ig) = cam_out(c)%precsl(i)*1000._r8
          a2x_a%rAttr(index_a2x_Faxa_swndr,ig) = cam_out(c)%soll(i)   
          a2x_a%rAttr(index_a2x_Faxa_swvdr,ig) = cam_out(c)%sols(i)   
          a2x_a%rAttr(index_a2x_Faxa_swndf,ig) = cam_out(c)%solld(i)  
          a2x_a%rAttr(index_a2x_Faxa_swvdf,ig) = cam_out(c)%solsd(i)  

          ! aerosol deposition fluxes
          a2x_a%rAttr(index_a2x_Faxa_bcphidry,ig) = cam_out(c)%bcphidry(i)
          a2x_a%rAttr(index_a2x_Faxa_bcphodry,ig) = cam_out(c)%bcphodry(i)
          a2x_a%rAttr(index_a2x_Faxa_bcphiwet,ig) = cam_out(c)%bcphiwet(i)
          a2x_a%rAttr(index_a2x_Faxa_ocphidry,ig) = cam_out(c)%ocphidry(i)
          a2x_a%rAttr(index_a2x_Faxa_ocphodry,ig) = cam_out(c)%ocphodry(i)
          a2x_a%rAttr(index_a2x_Faxa_ocphiwet,ig) = cam_out(c)%ocphiwet(i)
          a2x_a%rAttr(index_a2x_Faxa_dstwet1,ig)  = cam_out(c)%dstwet1(i)
          a2x_a%rAttr(index_a2x_Faxa_dstdry1,ig)  = cam_out(c)%dstdry1(i)
          a2x_a%rAttr(index_a2x_Faxa_dstwet2,ig)  = cam_out(c)%dstwet2(i)
          a2x_a%rAttr(index_a2x_Faxa_dstdry2,ig)  = cam_out(c)%dstdry2(i)
          a2x_a%rAttr(index_a2x_Faxa_dstwet3,ig)  = cam_out(c)%dstwet3(i)
          a2x_a%rAttr(index_a2x_Faxa_dstdry3,ig)  = cam_out(c)%dstdry3(i)
          a2x_a%rAttr(index_a2x_Faxa_dstwet4,ig)  = cam_out(c)%dstwet4(i)
          a2x_a%rAttr(index_a2x_Faxa_dstdry4,ig)  = cam_out(c)%dstdry4(i)

          if (index_a2x_Sa_co2prog /= 0) then
             a2x_a%rAttr(index_a2x_Sa_co2prog,ig) = cam_out(c)%co2prog(i) ! atm prognostic co2
          end if
          if (index_a2x_Sa_co2diag /= 0) then
             a2x_a%rAttr(index_a2x_Sa_co2diag,ig) = cam_out(c)%co2diag(i) ! atm diagnostic co2
          end if

          ig=ig+1
       end do
    end do
    
  end subroutine atm_export_mct

!===============================================================================

  subroutine atm_domain_mct( lsize, gsMap_a, dom_a )

    !-------------------------------------------------------------------
    !
    ! Arguments
    !
    integer        , intent(in)   :: lsize
    type(mct_gsMap), intent(in)   :: gsMap_a
    type(mct_ggrid), intent(inout):: dom_a  
    !
    ! Local Variables
    !
    integer  :: n,i,c,ncols           ! indices	
    real(r8) :: lats(pcols)           ! array of chunk latitudes
    real(r8) :: lons(pcols)           ! array of chunk longitude
    real(r8) :: area(pcols)           ! area in radians squared for each grid point
    real(r8), pointer  :: data(:)     ! temporary
    integer , pointer  :: idata(:)    ! temporary
    real(r8), parameter:: radtodeg = 180.0_r8/SHR_CONST_PI
    !-------------------------------------------------------------------
    !
    ! Initialize mct atm domain
    !
    call mct_gGrid_init( GGrid=dom_a, CoordChars=trim(seq_flds_dom_coord), OtherChars=trim(seq_flds_dom_other), lsize=lsize )
    !
    ! Allocate memory
    !
    allocate(data(lsize))
    !
    ! Initialize attribute vector with special value
    !
    call mct_gsMap_orderedPoints(gsMap_a, iam, idata)
    call mct_gGrid_importIAttr(dom_a,'GlobGridNum',idata,lsize)
    !
    ! Determine domain (numbering scheme is: West to East and South to North to South pole)
    ! Initialize attribute vector with special value
    !
    data(:) = -9999.0_R8 
    call mct_gGrid_importRAttr(dom_a,"lat"  ,data,lsize) 
    call mct_gGrid_importRAttr(dom_a,"lon"  ,data,lsize) 
    call mct_gGrid_importRAttr(dom_a,"area" ,data,lsize) 
    call mct_gGrid_importRAttr(dom_a,"aream",data,lsize) 
    data(:) = 0.0_R8     
    call mct_gGrid_importRAttr(dom_a,"mask" ,data,lsize) 
    data(:) = 1.0_R8
    call mct_gGrid_importRAttr(dom_a,"frac" ,data,lsize)
    !
    ! Fill in correct values for domain components
    !
    n=0
    do c = begchunk, endchunk
       ncols = get_ncols_p(c)
       call get_rlat_all_p(c, ncols, lats)
       do i=1,ncols
          n = n+1
          data(n) = lats(i)*radtodeg
       end do
    end do
    call mct_gGrid_importRAttr(dom_a,"lat",data,lsize) 

    n=0
    do c = begchunk, endchunk
       ncols = get_ncols_p(c)
       call get_rlon_all_p(c, ncols, lons)
       do i=1,ncols
          n = n+1
          data(n) = lons(i)*radtodeg
       end do
    end do
    call mct_gGrid_importRAttr(dom_a,"lon",data,lsize) 

    n=0
    do c = begchunk, endchunk
       ncols = get_ncols_p(c)
       call get_area_all_p(c, ncols, area)
       do i=1,ncols
          n = n+1
          data(n) = area(i) 
       end do
    end do
    call mct_gGrid_importRAttr(dom_a,"area",data,lsize) 

    n=0
    do c = begchunk, endchunk
       ncols = get_ncols_p(c)
       do i=1,ncols
          n = n+1
          data(n) = 1._r8 ! mask
       end do
    end do
    call mct_gGrid_importRAttr(dom_a,"mask"   ,data,lsize) 
    deallocate(data)

  end subroutine atm_domain_mct

!===========================================================================================
!
  subroutine atm_read_srfrest_mct( EClock, cdata_a, x2a_a, a2x_a)
    use cam_pio_utils
    !-----------------------------------------------------------------------
    !
    ! Arguments
    !
    type(ESMF_Clock),intent(in)    :: EClock
    type(seq_cdata), intent(inout) :: cdata_a
    type(mct_aVect), intent(inout) :: x2a_a
    type(mct_aVect), intent(inout) :: a2x_a
    ! 
    ! Local variables
    !
    integer         :: npts         ! array size
    integer         :: rcode        ! return error code
    type(mct_aVect) :: gData        ! global/gathered bundle data
    integer         :: yr_spec      ! Current year
    integer         :: mon_spec     ! Current month
    integer         :: day_spec     ! Current day
    integer         :: sec_spec     ! Current time of day (sec)
    !-----------------------------------------------------------------------
    !
    ! Determine and open surface restart dataset
    !
    integer, pointer :: dof(:)
    integer :: lnx, nf_x2a, nf_a2x, k
    real(r8), allocatable :: tmp(:)
    type(file_desc_t) :: file
    type(io_desc_t) :: iodesc
    type(var_desc_t) :: varid
    character(CL)    :: itemc       ! string converted to char
    type(mct_string) :: mstring     ! mct char type



    call seq_timemgr_EClockGetData( EClock, curr_yr=yr_spec,curr_mon=mon_spec, &
         curr_day=day_spec, curr_tod=sec_spec ) 
    fname_srf_cam = interpret_filename_spec( rsfilename_spec_cam, case=get_restcase(), &
         yr_spec=yr_spec, mon_spec=mon_spec, day_spec=day_spec, sec_spec= sec_spec )
    pname_srf_cam = trim(get_restartdir() )//fname_srf_cam
    call getfil(pname_srf_cam, fname_srf_cam)
    
    call cam_pio_openfile(File, fname_srf_cam, 0)
    call mct_gsmap_OrderedPoints(cdata_a%gsmap, iam, Dof)
    lnx = mct_gsmap_gsize(cdata_a%gsmap)
    call pio_initdecomp(pio_subsystem, pio_double, (/lnx/), dof, iodesc)
    allocate(tmp(size(dof)))
    deallocate(dof)
    
    nf_x2a = mct_aVect_nRattr(x2a_a)

    do k=1,nf_x2a
       call mct_aVect_getRList(mstring,k,x2a_a)
       itemc = mct_string_toChar(mstring)
       call mct_string_clean(mstring)

       call pio_seterrorhandling(File, pio_bcast_error)
       rcode = pio_inq_varid(File,'x2a_'//trim(itemc) ,varid)
       if (rcode == pio_noerr) then
          call pio_read_darray(File, varid, iodesc, tmp, rcode)
          x2a_a%rattr(k,:) = tmp(:)
       else
	  if (masterproc) then
             write(iulog,*)'srfrest warning: field ',trim(itemc),' is not on restart file'
             write(iulog,*)'for backwards compatibility will set it to 0'
          end if
          x2a_a%rattr(k,:) = 0._r8
       end if
       call pio_seterrorhandling(File, pio_internal_error)
    end do

    nf_a2x = mct_aVect_nRattr(a2x_a)

    do k=1,nf_a2x
       call mct_aVect_getRList(mstring,k,a2x_a)
       itemc = mct_string_toChar(mstring)
       call mct_string_clean(mstring)

       rcode = pio_inq_varid(File,'a2x_'//trim(itemc) ,varid)
       call pio_read_darray(File, varid, iodesc, tmp, rcode)
       a2x_a%rattr(k,:) = tmp(:)
    end do

    call pio_freedecomp(File,iodesc)
    call pio_closefile(File)
    deallocate(tmp)

  end subroutine atm_read_srfrest_mct
!
!===========================================================================================
!
  subroutine atm_write_srfrest_mct( cdata_a, x2a_a, a2x_a, & 
       yr_spec, mon_spec, day_spec, sec_spec)
    use cam_pio_utils
    !-----------------------------------------------------------------------
    !
    ! Arguments
    !
    type(seq_cdata), intent(in) :: cdata_a
    type(mct_aVect), intent(in) :: x2a_a
    type(mct_aVect), intent(in) :: a2x_a
    integer        , intent(in) :: yr_spec         ! Simulation year
    integer        , intent(in) :: mon_spec        ! Simulation month
    integer        , intent(in) :: day_spec        ! Simulation day
    integer        , intent(in) :: sec_spec        ! Seconds into current simulation day
    !
    ! Local variables
    !
    integer         :: rcode        ! return error code
    type(mct_aVect) :: gData        ! global/gathered bundle data
    !-----------------------------------------------------------------------
    !
    ! Determine and open surface restart dataset
    !

    integer, pointer :: dof(:)
    integer :: nf_x2a, nf_a2x, lnx, dimid(1), k
    type(file_desc_t) :: file
    type(var_desc_t), pointer :: varid_x2a(:), varid_a2x(:)
    type(io_desc_t)  :: iodesc
    character(CL)    :: itemc       ! string converted to char
    type(mct_string) :: mstring     ! mct char type


    fname_srf_cam = interpret_filename_spec( rsfilename_spec_cam, &
         yr_spec=yr_spec, mon_spec=mon_spec, day_spec=day_spec, sec_spec= sec_spec )
    call cam_pio_createfile(File, fname_srf_cam, 0)

    call mct_gsmap_OrderedPoints(cdata_a%gsmap, iam, Dof)
    lnx = mct_gsmap_gsize(cdata_a%gsmap)
    call pio_initdecomp(pio_subsystem, pio_double, (/lnx/), dof, iodesc)

    deallocate(dof)
    
    nf_x2a = mct_aVect_nRattr(x2a_a)
    allocate(varid_x2a(nf_x2a))
    
    rcode = pio_def_dim(File,'x2a_nx',lnx,dimid(1))
    do k = 1,nf_x2a
       call mct_aVect_getRList(mstring,k,x2a_a)
       itemc = mct_string_toChar(mstring)
       call mct_string_clean(mstring)
       rcode = pio_def_var(File,'x2a_'//trim(itemc),PIO_DOUBLE,dimid,varid_x2a(k))
       rcode = pio_put_att(File,varid_x2a(k),"_fillvalue",fillvalue)
    enddo

    nf_a2x = mct_aVect_nRattr(a2x_a)
    allocate(varid_a2x(nf_a2x))
    
    rcode = pio_def_dim(File,'a2x_nx',lnx,dimid(1))
    do k = 1,nf_a2x
       call mct_aVect_getRList(mstring,k,a2x_a)
       itemc = mct_string_toChar(mstring)
       call mct_string_clean(mstring)
       rcode = PIO_def_var(File,'a2x_'//trim(itemc),PIO_DOUBLE,dimid,varid_a2x(k))
       rcode = PIO_put_att(File,varid_a2x(k),"_fillvalue",fillvalue)
    enddo

    rcode = pio_enddef(File)  ! don't check return code, might be enddef already


    do k=1,nf_x2a
       call pio_write_darray(File, varid_x2a(k), iodesc, x2a_a%rattr(k,:), rcode)
    end do

    do k=1,nf_a2x
       call pio_write_darray(File, varid_a2x(k), iodesc, a2x_a%rattr(k,:), rcode)       
    end do

    deallocate(varid_x2a, varid_a2x)

    call pio_freedecomp(File,iodesc)
    call pio_closefile(file)


  end subroutine atm_write_srfrest_mct

  
!===============================================================================
!
!  The following subroutines are for wrf/cam coupling
!	
!===============================================================================	

  subroutine atm_import_self( c2x_c,cam_state, cam_tend)  ! juanxiong he

    !
    ! Arguments
    !
    type(mct_aVect)   , intent(inout) :: c2x_c
    type(physics_state), intent(inout) :: cam_state(begchunk:endchunk)
    type(physics_tend), intent(inout) :: cam_tend(begchunk:endchunk)  
    !
    ! Local variables
    !		
    integer  :: m             ! indices
    integer  :: i,j,k,lat,n,c,ig  ! indices
    integer  :: ncols         ! number of columns
    integer  :: dust_ndx

    ig=1
    do j=begchunk, endchunk
       ncols = get_ncols_p(j)   
       do i=1, ncols
          do k = 1, num_cam_levs 
             cam_state(j)%u(i,k) = c2x_c%rAttr(index_c2x_Sc_u3d(k)   ,ig)
             cam_state(j)%v(i,k) = c2x_c%rAttr(index_c2x_Sc_v3d(k)   ,ig)
             cam_state(j)%t(i,k) = c2x_c%rAttr(index_c2x_Sc_t3d(k)   ,ig)
             cam_state(j)%q(i,k,1) = c2x_c%rAttr(index_c2x_Sc_q3d(k)   ,ig)   ! note: q have several species and 1 is water vapor	  	  
             cam_state(j)%pdel(i,k) = c2x_c%rAttr(index_c2x_Fcxc_qtend(k)   ,ig)   
             cam_state(j)%pmid(i,k) = c2x_c%rAttr(index_c2x_Sc_p3d(k)   ,ig)
             cam_state(j)%zm(i,k) = c2x_c%rAttr(index_c2x_Sc_z3d(k)   ,ig) 
          end do
             cam_state(j)%ps(i) = c2x_c%rAttr(index_c2x_Sc_ps    ,ig) 
             cam_state(j)%phis(i) = c2x_c%rAttr(index_c2x_Sc_phis  ,ig) 
          ig=ig+1
      end do 
    end do

  end subroutine atm_import_self

!===============================================================================

  subroutine atm_import_wrf( x2c_c, c2x_c, wrf_state, wrf_tend, twoway_nudging)  ! juanxiong he

     use dycore,           only: dycore_is  !juanxiong he
    !
    ! Arguments
    !
    type(mct_aVect)   , intent(inout) :: x2c_c, c2x_c
    type(physics_state), intent(inout) :: wrf_state(begchunk:endchunk)
    type(physics_tend), intent(inout) :: wrf_tend(begchunk:endchunk)  
    integer, intent(in) :: twoway_nudging  
    !
    ! Local variables
    !		
    integer  :: m             ! indices
    integer  :: i,j,k,lat,n,c,ig  ! indices
    integer  :: ncols         ! number of columns
    integer  :: dust_ndx

    ig=1
    do j=begchunk, endchunk
       ncols = get_ncols_p(j)   
       do i=1, ncols
          
          wrf_state(j)%ps(i) = x2c_c%rAttr(index_x2c_Sx_ps   ,ig)
          do k = 1, num_cam_levs 

             wrf_state(j)%t(i,k) = x2c_c%rAttr(index_x2c_Sx_t3d(k)   ,ig)
           if(twoway_nudging.eq.0) then
             wrf_state(j)%u(i,k) = x2c_c%rAttr(index_x2c_Sx_u3d(k)   ,ig)
             wrf_state(j)%v(i,k) = x2c_c%rAttr(index_x2c_Sx_v3d(k)   ,ig)
             wrf_state(j)%q(i,k,1) = x2c_c%rAttr(index_x2c_Sx_q3d(k)   ,ig)    ! note: q have several species and 1 is water vapor	  	  
           else
             if(wrf_state(j)%t(i,k).gt.0) then
              wrf_tend(j)%dudt(i,k) = x2c_c%rAttr(index_x2c_Fcxx_dudt(k), ig) 
              wrf_tend(j)%dvdt(i,k) = x2c_c%rAttr(index_x2c_Fcxx_dvdt(k), ig)
              ! use pdel and pmid to store tendency 
              wrf_state(j)%pdel(i,k) = x2c_c%rAttr(index_x2c_Fcxx_dtdt(k), ig)
              wrf_state(j)%pmid(i,k) = x2c_c%rAttr(index_x2c_Fcxx_dqdt(k), ig)   
             else
              wrf_tend(j)%dudt(i,k) = 0.0_8
              wrf_tend(j)%dvdt(i,k) = 0.0_8
              ! use t and q to store tendency 
              wrf_state(j)%t(i,k) = 0.0_8
              wrf_state(j)%q(i,k,1) = 0.0_8
             endif
           endif

          end do
          ig=ig+1
      end do 
    end do

  end subroutine atm_import_wrf

!===============================================================================

  subroutine atm_export_wrf( c2x_c, cam_in, cam_state, cam_tend ) !juanxiong he
    !-------------------------------------------------------------------
    ! Arguments

    type(mct_aVect)    , intent(out) :: c2x_c
    type(cam_in_t), intent(in) :: cam_in(begchunk:endchunk)
    type(physics_state), intent(in) :: cam_state(begchunk:endchunk)
    type(physics_tend), intent(in) :: cam_tend(begchunk:endchunk)
    !
    ! Local variables
    !
    integer :: avsize, avnat
    integer :: i,j,k,n,ig       ! indices
    integer :: ncols            ! Number of columns
    !-----------------------------------------------------------------------

    ! Copy from component arrays into chunk array data structure
    ! Rearrange data from chunk structure into lat-lon buffer and subsequently
    ! create attribute vector
    	    
    ig=1
    do j=begchunk, endchunk
       ncols = get_ncols_p(j)   
       do i=1, ncols
          c2x_c%rAttr(index_c2x_Sc_lat   ,ig) = cam_state(j)%lat(i)  ! wrf has its own lat/lon.Only need horiztonal remappingg. 
          c2x_c%rAttr(index_c2x_Sc_lon   ,ig) = cam_state(j)%lon(i)  
          
          do k = 1, num_cam_levs 
          c2x_c%rAttr(index_c2x_Sc_p3d(k)   ,ig) = cam_state(j)%pmid(i,k)
          c2x_c%rAttr(index_c2x_Sc_z3d(k)   ,ig) = cam_state(j)%zm(i,k)          
          c2x_c%rAttr(index_c2x_Sc_u3d(k)   ,ig) = cam_state(j)%u(i,k)   
          c2x_c%rAttr(index_c2x_Sc_v3d(k)   ,ig) = cam_state(j)%v(i,k)   
          c2x_c%rAttr(index_c2x_Sc_t3d(k)   ,ig) = cam_state(j)%t(i,k)    
	  c2x_c%rAttr(index_c2x_Sc_w3d(k)   ,ig) = cam_state(j)%rh(i,k) ! use omega as the media of rh 	  
          c2x_c%rAttr(index_c2x_Sc_q3d(k)   ,ig) = cam_state(j)%q(i,k,1) ! for nudging 
          c2x_c%rAttr(index_c2x_Fcxc_utend(k)   ,ig) = cam_tend(j)%dudt(i,k) ! the difference  not the tendency
          c2x_c%rAttr(index_c2x_Fcxc_vtend(k)   ,ig) = cam_tend(j)%dvdt(i,k) ! the difference  not the tendency
          c2x_c%rAttr(index_c2x_Fcxc_qtend(k)   ,ig) = cam_state(j)%pdel(i,k) ! used for twoway coupling 
          c2x_c%rAttr(index_c2x_Fcxc_ttend(k)   ,ig) = cam_tend(j)%dtdt(i,k) ! the difference  not the tendency
          end do
          c2x_c%rAttr(index_c2x_Sc_q3d(1)   ,ig) = cam_state(j)%psl(i)  ! use q(1) as the PSL, since q in the uppermost atmospheric layer is almost zero and not used in WRF
          c2x_c%rAttr(index_c2x_Sc_q3d(2)   ,ig) = cam_in(j)%tref(i)  ! use q(2) as the 2m temperature, since q in the uppermost atmospheric layer is almost zero and not used in WRF 

          c2x_c%rAttr(index_c2x_Sc_ps    ,ig) = cam_state(j)%ps(i) 
          c2x_c%rAttr(index_c2x_Sc_phis  ,ig) = cam_state(j)%phis(i)          
          c2x_c%rAttr(index_c2x_Sc_ts  ,ig) = cam_in(j)%ts(i)                ! for surface layer information           
          c2x_c%rAttr(index_c2x_Sc_sst  ,ig) = cam_in(j)%sst(i)   
          c2x_c%rAttr(index_c2x_Sc_snowhland  ,ig) = cam_in(j)%snowhland(i)  
          c2x_c%rAttr(index_c2x_Sc_snowhice   ,ig) = cam_in(j)%snowhice(i)    
          c2x_c%rAttr(index_c2x_Sc_seaice  ,ig) = cam_in(j)%icefrac(i)  
          c2x_c%rAttr(index_c2x_Sc_ocnfrac  ,ig) = cam_in(j)%landfrac(i)   ! for land surface 

          do k = 1, num_soil_layers
          c2x_c%rAttr(index_c2x_Sc_soildepth(k),ig) = cam_in(j)%soildepth(i,k)   ! for soil layer information
          c2x_c%rAttr(index_c2x_Sc_soilthick(k),ig) = cam_in(j)%soilthick(i,k)
          c2x_c%rAttr(index_c2x_Sc_soilt(k),ig) = cam_in(j)%soilt(i,k)
          c2x_c%rAttr(index_c2x_Sc_soilm(k),ig) = cam_in(j)%soilm(i,k)
          enddo
       ig=ig+1
!       print *, c2x_c%rAttr(index_c2x_Sc_soildepth(1),ig),cam_state(j)%lat(i)
      end do 
    end do

  end subroutine atm_export_wrf

!===============================================================================
	
  subroutine atm_SetgsMap_wrf( mpicom_atm, ATMID, GSMap_cc )  !juanxiong he
    use phys_grid, only : get_nlcols_p
    !-------------------------------------------------------------------
    !
    ! Arguments
    !
    integer        , intent(in)  :: mpicom_atm
    integer        , intent(in)  :: ATMID
    type(mct_gsMap), intent(out) :: GSMap_cc
    !
    ! Local variables
    !
    integer, allocatable :: gindex(:)
    integer :: i, k, n, c, ncols, sizebuf, nlcols
    integer :: ier            ! error status
    !-------------------------------------------------------------------

    ! Build the atmosphere grid numbering for MCT
    ! NOTE:  Numbering scheme is: West to East, bottom to top, and South to North
    ! starting at south pole.  Should be the same as what's used in SCRIP
    
    ! Determine global seg map

    sizebuf=0
    do c = begchunk, endchunk
       ncols = get_ncols_p(c)
        do i = 1,ncols
          sizebuf = sizebuf+1
        end do
    end do

    allocate(gindex(sizebuf))

    n=0
    do c = begchunk, endchunk
       ncols = get_ncols_p(c)
        do i = 1,ncols
          n=n+1
          gindex(n) = get_gcol_p(c,i)  !global index
        end do
    end do

    nlcols = get_nlcols_p()
    call mct_gsMap_init( gsMap_cc, gindex, mpicom_atm, ATMID, nlcols, ngcols)

    deallocate(gindex)

  end subroutine atm_SetgsMap_wrf

!===============================================================================

  subroutine atm_domain_wrf( lsize, gsMap_c, dom_c )  !juanxiong he
    !-------------------------------------------------------------------
    !
    ! Arguments
    !
    integer        , intent(in)   :: lsize
    type(mct_gsMap), intent(in)   :: gsMap_c
    type(mct_ggrid), intent(inout):: dom_c  
    !
    ! Local Variables
    !
    integer  :: n,i,k,c,ncols           ! indices	
    real(r8) :: lats(pcols)           ! array of chunk latitudes
    real(r8) :: lons(pcols)           ! array of chunk longitude
    real(r8) :: area(pcols)           ! area in radians squared for each grid point
    real(r8), pointer  :: data(:)     ! temporary
    integer , pointer  :: idata(:)    ! temporary
    real(r8), parameter:: radtodeg = 180.0_r8/SHR_CONST_PI
    !-------------------------------------------------------------------
    !
    ! Initialize mct atm domain
    !
    call mct_gGrid_init( GGrid=dom_c, CoordChars=trim(seq_flds_dom_coord), OtherChars=trim(seq_flds_dom_other), lsize=lsize )
    !
    ! Allocate memory
    !
    allocate(data(lsize))
    !
    ! Initialize attribute vector with special value
    !
    call mct_gsMap_orderedPoints(gsMap_c, iam, idata)
    call mct_gGrid_importIAttr(dom_c,'GlobGridNum',idata,lsize)
    !
    ! Determine domain (numbering scheme is: West to East and South to North to South pole)
    ! Initialize attribute vector with special value
    !
    data(:) = -9999.0_R8 
    call mct_gGrid_importRAttr(dom_c,"lat"  ,data,lsize) 
    call mct_gGrid_importRAttr(dom_c,"lon"  ,data,lsize) 
    call mct_gGrid_importRAttr(dom_c,"area" ,data,lsize) 
    call mct_gGrid_importRAttr(dom_c,"aream",data,lsize) 
    data(:) = 0.0_R8     
    call mct_gGrid_importRAttr(dom_c,"mask" ,data,lsize) 
    data(:) = 1.0_R8
    call mct_gGrid_importRAttr(dom_c,"frac" ,data,lsize)
    !
    ! Fill in correct values for domain components
    !
    n=0
    do c = begchunk, endchunk
       ncols = get_ncols_p(c)
       call get_rlat_all_p(c, ncols, lats)
       do i=1,ncols
          n = n+1
          data(n) = lats(i)*radtodeg
       end do
    end do
    call mct_gGrid_importRAttr(dom_c,"lat",data,lsize) 

    n=0
    do c = begchunk, endchunk
       ncols = get_ncols_p(c)
       call get_rlon_all_p(c, ncols, lons)
       do i=1,ncols
          n = n+1
          data(n) = lons(i)*radtodeg
       end do
    end do
    call mct_gGrid_importRAttr(dom_c,"lon",data,lsize) 

    n=0
    do c = begchunk, endchunk
       ncols = get_ncols_p(c)
       call get_area_all_p(c, ncols, area)
       do i=1,ncols
          n = n+1
          data(n) = area(i) 
       end do
    end do
    call mct_gGrid_importRAttr(dom_c,"area",data,lsize) 

    n=0
    do c = begchunk, endchunk
       ncols = get_ncols_p(c)
       do i=1,ncols
          n = n+1
          data(n) = 1._r8 ! mask
       end do
    end do
    call mct_gGrid_importRAttr(dom_c,"mask"   ,data,lsize) 
    deallocate(data)

  end subroutine atm_domain_wrf  
!
!===========================================================================================
!
  subroutine atm_read_wrfrest_mct( EClock, cdata_c, x2c_c, c2x_c)
    use cam_pio_utils
    !-----------------------------------------------------------------------
    !
    ! Arguments
    !
    type(ESMF_Clock),intent(in)    :: EClock
    type(seq_cdata), intent(inout) :: cdata_c
    type(mct_aVect), intent(inout) :: x2c_c
    type(mct_aVect), intent(inout) :: c2x_c
    ! 
    ! Local variables
    !
    integer         :: npts         ! array size
    integer         :: rcode        ! return error code
    type(mct_aVect) :: gData        ! global/gathered bundle data
    integer         :: yr_spec      ! Current year
    integer         :: mon_spec     ! Current month
    integer         :: day_spec     ! Current day
    integer         :: sec_spec     ! Current time of day (sec)
    !-----------------------------------------------------------------------
    !
    ! Determine and open surface restart dataset
    !
    integer, pointer :: dof(:)
    integer :: lnx, nf_x2c, nf_c2x, k
    real(r8), allocatable :: tmp(:)
    type(file_desc_t) :: file
    type(io_desc_t) :: iodesc
    type(var_desc_t) :: varid
    character(CL)    :: itemc       ! string converted to char
    type(mct_string) :: mstring     ! mct char type

    call seq_timemgr_EClockGetData( EClock, curr_yr=yr_spec,curr_mon=mon_spec, &
         curr_day=day_spec, curr_tod=sec_spec ) 
    fname_srf_cam = interpret_filename_spec( wrsfilename_spec_cam, case=get_restcase(), &
         yr_spec=yr_spec, mon_spec=mon_spec, day_spec=day_spec, sec_spec= sec_spec )
    pname_srf_cam = trim(get_restartdir() )//fname_srf_cam
    call getfil(pname_srf_cam, fname_srf_cam)
    
    call cam_pio_openfile(File, fname_srf_cam, 0)
    call mct_gsmap_OrderedPoints(cdata_c%gsmap, iam, Dof)
    lnx = mct_gsmap_gsize(cdata_c%gsmap)
    call pio_initdecomp(pio_subsystem, pio_double, (/lnx/), dof, iodesc)
    allocate(tmp(size(dof)))
    deallocate(dof)
    
    nf_x2c = mct_aVect_nRattr(x2c_c)

    do k=1,nf_x2c
       call mct_aVect_getRList(mstring,k,x2c_c)
       itemc = mct_string_toChar(mstring)
       call mct_string_clean(mstring)

       call pio_seterrorhandling(File, pio_bcast_error)
       rcode = pio_inq_varid(File,'x2c_'//trim(itemc) ,varid)
       if (rcode == pio_noerr) then
          call pio_read_darray(File, varid, iodesc, tmp, rcode)
          x2c_c%rattr(k,:) = tmp(:)
       else
	  if (masterproc) then
             write(iulog,*)'srfrest warning: field ',trim(itemc),' is not on restart file'
             write(iulog,*)'for backwards compatibility will set it to 0'
          end if
          x2c_c%rattr(k,:) = 0._r8
       end if
       call pio_seterrorhandling(File, pio_internal_error)
    end do

    nf_c2x = mct_aVect_nRattr(c2x_c)

    do k=1,nf_c2x
       call mct_aVect_getRList(mstring,k,c2x_c)
       itemc = mct_string_toChar(mstring)
       call mct_string_clean(mstring)

       rcode = pio_inq_varid(File,'c2x_'//trim(itemc) ,varid)
       call pio_read_darray(File, varid, iodesc, tmp, rcode)
       c2x_c%rattr(k,:) = tmp(:)
    end do

    call pio_freedecomp(File,iodesc)
    call pio_closefile(File)
    deallocate(tmp)

  end subroutine atm_read_wrfrest_mct
!
!===========================================================================================
!
  subroutine atm_write_wrfrest_mct( cdata_c, x2c_c, c2x_c, & 
       yr_spec, mon_spec, day_spec, sec_spec)
    use cam_pio_utils
    !-----------------------------------------------------------------------
    !
    ! Arguments
    !
    type(seq_cdata), intent(in) :: cdata_c
    type(mct_aVect), intent(in) :: x2c_c
    type(mct_aVect), intent(in) :: c2x_c
    integer        , intent(in) :: yr_spec         ! Simulation year
    integer        , intent(in) :: mon_spec        ! Simulation month
    integer        , intent(in) :: day_spec        ! Simulation day
    integer        , intent(in) :: sec_spec        ! Seconds into current simulation day
    !
    ! Local variables
    !
    integer         :: rcode        ! return error code
    type(mct_aVect) :: gData        ! global/gathered bundle data
    !-----------------------------------------------------------------------
    !
    ! Determine and open surface restart dataset
    !

    integer, pointer :: dof(:)
    integer :: nf_x2c, nf_c2x, lnx, dimid(1), k
    type(file_desc_t) :: file
    type(var_desc_t), pointer :: varid_x2c(:), varid_c2x(:)
    type(io_desc_t)  :: iodesc
    character(CL)    :: itemc       ! string converted to char
    type(mct_string) :: mstring     ! mct char type


    fname_srf_cam = interpret_filename_spec( wrsfilename_spec_cam, &
         yr_spec=yr_spec, mon_spec=mon_spec, day_spec=day_spec, sec_spec= sec_spec )
    call cam_pio_createfile(File, fname_srf_cam, 0)

    call mct_gsmap_OrderedPoints(cdata_c%gsmap, iam, Dof)
    lnx = mct_gsmap_gsize(cdata_c%gsmap)
    call pio_initdecomp(pio_subsystem, pio_double, (/lnx/), dof, iodesc)

    deallocate(dof)
    
    nf_x2c = mct_aVect_nRattr(x2c_c)
    allocate(varid_x2c(nf_x2c))
    
    rcode = pio_def_dim(File,'x2c_nx',lnx,dimid(1))
    do k = 1,nf_x2c
       call mct_aVect_getRList(mstring,k,x2c_c)
       itemc = mct_string_toChar(mstring)
       call mct_string_clean(mstring)
       rcode = pio_def_var(File,'x2c_'//trim(itemc),PIO_DOUBLE,dimid,varid_x2c(k))
       rcode = pio_put_att(File,varid_x2c(k),"_fillvalue",fillvalue)
    enddo

    nf_c2x = mct_aVect_nRattr(c2x_c)
    allocate(varid_c2x(nf_c2x))
    
    rcode = pio_def_dim(File,'c2x_nx',lnx,dimid(1))
    do k = 1,nf_c2x
       call mct_aVect_getRList(mstring,k,c2x_c)
       itemc = mct_string_toChar(mstring)
       call mct_string_clean(mstring)
       rcode = PIO_def_var(File,'c2x_'//trim(itemc),PIO_DOUBLE,dimid,varid_c2x(k))
       rcode = PIO_put_att(File,varid_c2x(k),"_fillvalue",fillvalue)
    enddo

    rcode = pio_enddef(File)  ! don't check return code, might be enddef already

    do k=1,nf_x2c
       call pio_write_darray(File, varid_x2c(k), iodesc, x2c_c%rattr(k,:), rcode)
    end do

    do k=1,nf_c2x
       call pio_write_darray(File, varid_c2x(k), iodesc, c2x_c%rattr(k,:), rcode)       
    end do

    deallocate(varid_x2c, varid_c2x)

    call pio_freedecomp(File,iodesc)
    call pio_closefile(file)


  end subroutine atm_write_wrfrest_mct
!
!===========================================================================================
!
  subroutine wrf_read_srfrest_mct( EClock, cdata_c, x2c_c, c2x_c)
    use cam_pio_utils
    !-----------------------------------------------------------------------
    !
    ! Arguments
    !
    type(ESMF_Clock),intent(in)    :: EClock
    type(seq_cdata), intent(inout) :: cdata_c
    type(mct_aVect), intent(inout) :: x2c_c
    type(mct_aVect), intent(inout) :: c2x_c
    ! 
    ! Local variables
    !
    integer         :: npts         ! array size
    integer         :: rcode        ! return error code
    type(mct_aVect) :: gData        ! global/gathered bundle data
    integer         :: yr_spec      ! Current year
    integer         :: mon_spec     ! Current month
    integer         :: day_spec     ! Current day
    integer         :: sec_spec     ! Current time of day (sec)
    !-----------------------------------------------------------------------
    !
    ! Determine and open surface restart dataset
    !
    integer, pointer :: dof(:)
    integer :: lnx, nf_x2c, nf_c2x, k
    real(r8), allocatable :: tmp(:)
    type(file_desc_t) :: file
    type(io_desc_t) :: iodesc
    type(var_desc_t) :: varid
    character(CL)    :: itemc       ! string converted to char
    type(mct_string) :: mstring     ! mct char type

    call seq_timemgr_EClockGetData( EClock, curr_yr=yr_spec,curr_mon=mon_spec, &
         curr_day=day_spec, curr_tod=sec_spec ) 
    fname_srf_wrf = interpret_filename_spec( wrsfilename_spec_wrf, case=get_restcase(), &
         yr_spec=yr_spec, mon_spec=mon_spec, day_spec=day_spec, sec_spec= sec_spec )
    pname_srf_wrf = trim(get_restartdir() )//fname_srf_wrf
    call getfil(pname_srf_wrf, fname_srf_wrf)
    
    call cam_pio_openfile(File, fname_srf_wrf, 0)
    call mct_gsmap_OrderedPoints(cdata_c%gsmap, iam, Dof)
    lnx = mct_gsmap_gsize(cdata_c%gsmap)
    call pio_initdecomp(pio_subsystem, pio_double, (/lnx/), dof, iodesc)
    allocate(tmp(size(dof)))
    deallocate(dof)
    
    nf_x2c = mct_aVect_nRattr(x2c_c)

    do k=1,nf_x2c
       call mct_aVect_getRList(mstring,k,x2c_c)
       itemc = mct_string_toChar(mstring)
       call mct_string_clean(mstring)

       call pio_seterrorhandling(File, pio_bcast_error)
       rcode = pio_inq_varid(File,'x2c_'//trim(itemc) ,varid)
       if (rcode == pio_noerr) then
          call pio_read_darray(File, varid, iodesc, tmp, rcode)
          x2c_c%rattr(k,:) = tmp(:)
       else
	  if (masterproc) then
             write(iulog,*)'srfrest warning: field ',trim(itemc),' is not on restart file'
             write(iulog,*)'for backwards compatibility will set it to 0'
          end if
          x2c_c%rattr(k,:) = 0._r8
       end if
       call pio_seterrorhandling(File, pio_internal_error)
    end do

    nf_c2x = mct_aVect_nRattr(c2x_c)

    do k=1,nf_c2x
       call mct_aVect_getRList(mstring,k,c2x_c)
       itemc = mct_string_toChar(mstring)
       call mct_string_clean(mstring)

       rcode = pio_inq_varid(File,'c2x_'//trim(itemc) ,varid)
       call pio_read_darray(File, varid, iodesc, tmp, rcode)
       c2x_c%rattr(k,:) = tmp(:)
    end do

    call pio_freedecomp(File,iodesc)
    call pio_closefile(File)
    deallocate(tmp)

  end subroutine wrf_read_srfrest_mct
!
!===========================================================================================
!
  subroutine wrf_write_srfrest_mct( cdata_c, x2c_c, c2x_c, & 
       yr_spec, mon_spec, day_spec, sec_spec)
    use cam_pio_utils
    !-----------------------------------------------------------------------
    !
    ! Arguments
    !
    type(seq_cdata), intent(in) :: cdata_c
    type(mct_aVect), intent(in) :: x2c_c
    type(mct_aVect), intent(in) :: c2x_c
    integer        , intent(in) :: yr_spec         ! Simulation year
    integer        , intent(in) :: mon_spec        ! Simulation month
    integer        , intent(in) :: day_spec        ! Simulation day
    integer        , intent(in) :: sec_spec        ! Seconds into current simulation day
    !
    ! Local variables
    !
    integer         :: rcode        ! return error code
    type(mct_aVect) :: gData        ! global/gathered bundle data
    !-----------------------------------------------------------------------
    !
    ! Determine and open surface restart dataset
    !

    integer, pointer :: dof(:)
    integer :: nf_x2c, nf_c2x, lnx, dimid(1), k
    type(file_desc_t) :: file
    type(var_desc_t), pointer :: varid_x2c(:), varid_c2x(:)
    type(io_desc_t)  :: iodesc
    character(CL)    :: itemc       ! string converted to char
    type(mct_string) :: mstring     ! mct char type


    fname_srf_wrf = interpret_filename_spec( wrsfilename_spec_wrf, &
         yr_spec=yr_spec, mon_spec=mon_spec, day_spec=day_spec, sec_spec= sec_spec )
    call cam_pio_createfile(File, fname_srf_wrf, 0)

    call mct_gsmap_OrderedPoints(cdata_c%gsmap, iam, Dof)
    lnx = mct_gsmap_gsize(cdata_c%gsmap)
    call pio_initdecomp(pio_subsystem, pio_double, (/lnx/), dof, iodesc)

    deallocate(dof)
    
    nf_x2c = mct_aVect_nRattr(x2c_c)
    allocate(varid_x2c(nf_x2c))
    
    rcode = pio_def_dim(File,'x2c_nx',lnx,dimid(1))
    do k = 1,nf_x2c
       call mct_aVect_getRList(mstring,k,x2c_c)
       itemc = mct_string_toChar(mstring)
       call mct_string_clean(mstring)
       rcode = pio_def_var(File,'x2c_'//trim(itemc),PIO_DOUBLE,dimid,varid_x2c(k))
       rcode = pio_put_att(File,varid_x2c(k),"_fillvalue",fillvalue)
    enddo

    nf_c2x = mct_aVect_nRattr(c2x_c)
    allocate(varid_c2x(nf_c2x))
    
    rcode = pio_def_dim(File,'c2x_nx',lnx,dimid(1))
    do k = 1,nf_c2x
       call mct_aVect_getRList(mstring,k,c2x_c)
       itemc = mct_string_toChar(mstring)
       call mct_string_clean(mstring)
       rcode = PIO_def_var(File,'c2x_'//trim(itemc),PIO_DOUBLE,dimid,varid_c2x(k))
       rcode = PIO_put_att(File,varid_c2x(k),"_fillvalue",fillvalue)
    enddo

    rcode = pio_enddef(File)  ! don't check return code, might be enddef already

    do k=1,nf_x2c
       call pio_write_darray(File, varid_x2c(k), iodesc, x2c_c%rattr(k,:), rcode)
    end do

    do k=1,nf_c2x
       call pio_write_darray(File, varid_c2x(k), iodesc, c2x_c%rattr(k,:), rcode)       
    end do

    deallocate(varid_x2c, varid_c2x)

    call pio_freedecomp(File,iodesc)
    call pio_closefile(file)


  end subroutine wrf_write_srfrest_mct
!===========================================================================================
!
  subroutine wrf_read_srfrest_mct_app( EClock, cdata_c, x2c_c, c2x_c)
    use cam_pio_utils
    !-----------------------------------------------------------------------
    !
    ! Arguments
    !
    type(ESMF_Clock),intent(in)    :: EClock
    type(seq_cdata), intent(inout) :: cdata_c
    type(mct_aVect), intent(inout) :: x2c_c
    type(mct_aVect), intent(inout) :: c2x_c
    ! 
    ! Local variables
    !
    integer         :: npts         ! array size
    integer         :: rcode        ! return error code
    type(mct_aVect) :: gData        ! global/gathered bundle data
    integer         :: yr_spec      ! Current year
    integer         :: mon_spec     ! Current month
    integer         :: day_spec     ! Current day
    integer         :: sec_spec     ! Current time of day (sec)
    !-----------------------------------------------------------------------
    !
    ! Determine and open surface restart dataset
    !
    integer, pointer :: dof(:)
    integer :: lnx, nf_x2c, nf_c2x, k
    real(r8), allocatable :: tmp(:)
    type(file_desc_t) :: file
    type(io_desc_t) :: iodesc
    type(var_desc_t) :: varid
    character(CL)    :: itemc       ! string converted to char
    type(mct_string) :: mstring     ! mct char type

    call seq_timemgr_EClockGetData( EClock, curr_yr=yr_spec,curr_mon=mon_spec, &
         curr_day=day_spec, curr_tod=sec_spec )
    fname_srf_wrf = interpret_filename_spec( wrsfilename_spec_wrfold, case=get_restcase(), &
         yr_spec=yr_spec, mon_spec=mon_spec, day_spec=day_spec, sec_spec= sec_spec )
    pname_srf_wrf = trim(get_restartdir() )//fname_srf_wrf
    call getfil(pname_srf_wrf, fname_srf_wrf)

    call cam_pio_openfile(File, fname_srf_wrf, 0)
    call mct_gsmap_OrderedPoints(cdata_c%gsmap, iam, Dof)
    lnx = mct_gsmap_gsize(cdata_c%gsmap)
    call pio_initdecomp(pio_subsystem, pio_double, (/lnx/), dof, iodesc)
    allocate(tmp(size(dof)))
    deallocate(dof)

    nf_x2c = mct_aVect_nRattr(x2c_c)

    do k=1,nf_x2c
       call mct_aVect_getRList(mstring,k,x2c_c)
       itemc = mct_string_toChar(mstring)
       call mct_string_clean(mstring)

       call pio_seterrorhandling(File, pio_bcast_error)
       rcode = pio_inq_varid(File,'x2c_'//trim(itemc) ,varid)
       if (rcode == pio_noerr) then
          call pio_read_darray(File, varid, iodesc, tmp, rcode)
          x2c_c%rattr(k,:) = tmp(:)
       else
          if (masterproc) then
             write(iulog,*)'srfrest warning: field ',trim(itemc),' is not on restart file'
             write(iulog,*)'for backwards compatibility will set it to 0'
          end if
          x2c_c%rattr(k,:) = 0._r8
       end if
       call pio_seterrorhandling(File, pio_internal_error)
    end do

    nf_c2x = mct_aVect_nRattr(c2x_c)

    do k=1,nf_c2x
       call mct_aVect_getRList(mstring,k,c2x_c)
       itemc = mct_string_toChar(mstring)
       call mct_string_clean(mstring)

       rcode = pio_inq_varid(File,'c2x_'//trim(itemc) ,varid)
       call pio_read_darray(File, varid, iodesc, tmp, rcode)
       c2x_c%rattr(k,:) = tmp(:)
    end do

    call pio_freedecomp(File,iodesc)
    call pio_closefile(File)
    deallocate(tmp)

  end subroutine wrf_read_srfrest_mct_app
!
!===========================================================================================
!
  subroutine wrf_write_srfrest_mct_app( cdata_c, x2c_c, c2x_c, &
       yr_spec, mon_spec, day_spec, sec_spec)
    use cam_pio_utils
    !-----------------------------------------------------------------------
    !
    ! Arguments
    !
    type(seq_cdata), intent(in) :: cdata_c
    type(mct_aVect), intent(in) :: x2c_c
    type(mct_aVect), intent(in) :: c2x_c
    integer        , intent(in) :: yr_spec         ! Simulation year
    integer        , intent(in) :: mon_spec        ! Simulation month
    integer        , intent(in) :: day_spec        ! Simulation day
    integer        , intent(in) :: sec_spec        ! Seconds into current simulation day
    !
    ! Local variables
    !
    integer         :: rcode        ! return error code
    type(mct_aVect) :: gData        ! global/gathered bundle data
    !-----------------------------------------------------------------------
    !
    ! Determine and open surface restart dataset
    !

    integer, pointer :: dof(:)
    integer :: nf_x2c, nf_c2x, lnx, dimid(1), k
    type(file_desc_t) :: file
    type(var_desc_t), pointer :: varid_x2c(:), varid_c2x(:)
    type(io_desc_t)  :: iodesc
    character(CL)    :: itemc       ! string converted to char
    type(mct_string) :: mstring     ! mct char type

    fname_srf_wrf = interpret_filename_spec( wrsfilename_spec_wrfold, &
         yr_spec=yr_spec, mon_spec=mon_spec, day_spec=day_spec, sec_spec= sec_spec )
    call cam_pio_createfile(File, fname_srf_wrf, 0)

    call mct_gsmap_OrderedPoints(cdata_c%gsmap, iam, Dof)
    lnx = mct_gsmap_gsize(cdata_c%gsmap)
    call pio_initdecomp(pio_subsystem, pio_double, (/lnx/), dof, iodesc)

    deallocate(dof)

    nf_x2c = mct_aVect_nRattr(x2c_c)
    allocate(varid_x2c(nf_x2c))

    rcode = pio_def_dim(File,'x2c_nx',lnx,dimid(1))
    do k = 1,nf_x2c
       call mct_aVect_getRList(mstring,k,x2c_c)
       itemc = mct_string_toChar(mstring)
       call mct_string_clean(mstring)
       rcode = pio_def_var(File,'x2c_'//trim(itemc),PIO_DOUBLE,dimid,varid_x2c(k))
       rcode = pio_put_att(File,varid_x2c(k),"_fillvalue",fillvalue)
    enddo

    nf_c2x = mct_aVect_nRattr(c2x_c)
    allocate(varid_c2x(nf_c2x))

    rcode = pio_def_dim(File,'c2x_nx',lnx,dimid(1))
    do k = 1,nf_c2x
       call mct_aVect_getRList(mstring,k,c2x_c)
       itemc = mct_string_toChar(mstring)
       call mct_string_clean(mstring)
       rcode = PIO_def_var(File,'c2x_'//trim(itemc),PIO_DOUBLE,dimid,varid_c2x(k))
       rcode = PIO_put_att(File,varid_c2x(k),"_fillvalue",fillvalue)
    enddo

    rcode = pio_enddef(File)  ! don't check return code, might be enddef already

    do k=1,nf_x2c
       call pio_write_darray(File, varid_x2c(k), iodesc, x2c_c%rattr(k,:), rcode)
    end do

    do k=1,nf_c2x
       call pio_write_darray(File, varid_c2x(k), iodesc, c2x_c%rattr(k,:), rcode)
    end do

    deallocate(varid_x2c, varid_c2x)

    call pio_freedecomp(File,iodesc)
    call pio_closefile(file)


  end subroutine wrf_write_srfrest_mct_app
!==========================================================================================
!
!  The following subroutines are for geatm/cam coupling
!	
!==========================================================================================	
  subroutine atm_import_geatm( x2c_c, c2x_c, wrf_state, wrf_tend, twoway_nudging)  ! juanxiong he

     use dycore,           only: dycore_is  !juanxiong he
    !
    ! Arguments
    !
    type(mct_aVect)   , intent(inout) :: x2c_c, c2x_c
    type(physics_state), intent(inout) :: wrf_state(begchunk:endchunk)
    type(physics_tend), intent(inout) :: wrf_tend(begchunk:endchunk)  
    integer, intent(in) :: twoway_nudging  
    !
    ! Local variables
    !		
    integer  :: m             ! indices
    integer  :: i,j,k,lat,n,c,ig  ! indices
    integer  :: ncols         ! number of columns
    integer  :: dust_ndx


  end subroutine atm_import_geatm

!===============================================================================

  subroutine atm_export_geatm( c2x_c, cam_in, cam_out, cam_state, cam_tend ) !juanxiong he
    !-------------------------------------------------------------------
    ! Arguments

    type(mct_aVect)    , intent(out) :: c2x_c
    type(cam_in_t), intent(in) :: cam_in(begchunk:endchunk)
    type(cam_out_t), intent(in) :: cam_out(begchunk:endchunk)
    type(physics_state), intent(in) :: cam_state(begchunk:endchunk)
    type(physics_tend), intent(in) :: cam_tend(begchunk:endchunk)
    !
    ! Local variables
    !
    integer :: avsize, avnat
    integer :: i,j,k,n,ig       ! indices
    integer :: ncols            ! Number of columns
    real :: vv
    !-----------------------------------------------------------------------

    ! Copy from component arrays into chunk array data structure
    ! Rearrange data from chunk structure into lat-lon buffer and subsequently
    ! create attribute vector
    	    
    ig=1
    do j=begchunk, endchunk
       ncols = get_ncols_p(j)   
       do i=1, ncols
          c2x_c%rAttr(index_ca2x_Sca_lat   ,ig) = cam_state(j)%lat(i)  ! wrf has its own lat/lon.Only need horiztonal remappingg. 
          c2x_c%rAttr(index_ca2x_Sca_lon   ,ig) = cam_state(j)%lon(i)  
          c2x_c%rAttr(index_ca2x_Sca_phis  ,ig) = cam_state(j)%phis(i)
          do k = 1, num_cam_levs 
          c2x_c%rAttr(index_ca2x_Sca_p3d(k)   ,ig) = cam_state(j)%pmid(i,k)
          c2x_c%rAttr(index_ca2x_Sca_z3d(k)   ,ig) = cam_state(j)%zm(i,k)          
          c2x_c%rAttr(index_ca2x_Sca_u3d(k)   ,ig) = cam_state(j)%u(i,k)   
          c2x_c%rAttr(index_ca2x_Sca_v3d(k)   ,ig) = cam_state(j)%v(i,k)   
          c2x_c%rAttr(index_ca2x_Sca_t3d(k)   ,ig) = cam_state(j)%t(i,k)    
	  c2x_c%rAttr(index_ca2x_Sca_rh3d(k)   ,ig) = cam_state(j)%rh(i,k)  
          c2x_c%rAttr(index_ca2x_Sca_qv3d(k)   ,ig) = cam_state(j)%q(i,k,1) 
          c2x_c%rAttr(index_ca2x_Sca_qc3d(k)   ,ig) = cam_state(j)%q(i,k,2) 
          c2x_c%rAttr(index_ca2x_Sca_qi3d(k)   ,ig) = cam_state(j)%q(i,k,3) 
          c2x_c%rAttr(index_ca2x_Sca_taucldi3d(k)   ,ig) = cam_state(j)%taucldi3d(i,k) 
          c2x_c%rAttr(index_ca2x_Sca_taucldv3d(k)   ,ig) = cam_state(j)%taucldv3d(i,k) 
          end do

          c2x_c%rAttr(index_ca2x_Sca_q2   ,ig) = cam_in(j)%qref(i)  ! use q(2) as the 2m temperature, since q in the uppermost atmospheric layer is almost zero and not used in WRF
          c2x_c%rAttr(index_ca2x_Sca_rh2   ,ig) = cam_in(j)%rhref(i)  ! use q(2) as the 2m temperature, since q in the uppermost atmospheric layer is almost zero and not used in WRF
          c2x_c%rAttr(index_ca2x_Sca_t2   ,ig) = cam_in(j)%tref(i)  ! use q(2) as the 2m temperature, since q in the uppermost atmospheric layer is almost zero and not used in WRF 
          c2x_c%rAttr(index_ca2x_Sca_ps    ,ig) = cam_state(j)%ps(i)          
          c2x_c%rAttr(index_ca2x_Sca_ts  ,ig) = cam_in(j)%ts(i)                ! for surface layer information           
          c2x_c%rAttr(index_ca2x_Sca_sst  ,ig) = cam_in(j)%sst(i)   
          c2x_c%rAttr(index_ca2x_Sca_snowhland  ,ig) = cam_in(j)%snowhland(i)  
          c2x_c%rAttr(index_ca2x_Sca_snowhice   ,ig) = cam_in(j)%snowhice(i)    
          c2x_c%rAttr(index_ca2x_Sca_seaice  ,ig) = cam_in(j)%icefrac(i)  
          c2x_c%rAttr(index_ca2x_Sca_ocnfrac  ,ig) = cam_in(j)%landfrac(i)   ! for land surface      
          vv = sqrt(cam_in(j)%wsx(i)*cam_in(j)%wsx(i)+cam_in(j)%wsy(i)*cam_in(j)%wsy(i))
          if(vv.ne.0) then
          c2x_c%rAttr(index_ca2x_Sca_u10  ,ig) =  -1.0*cam_in(j)%u10(i)*cam_in(j)%wsx(i)/vv    ! u10          
          c2x_c%rAttr(index_ca2x_Sca_v10  ,ig) =  -1.0*cam_in(j)%u10(i)*cam_in(j)%wsy(i)/vv    ! v10
          else
          c2x_c%rAttr(index_ca2x_Sca_u10  ,ig) = 0.0    ! u10          
          c2x_c%rAttr(index_ca2x_Sca_v10  ,ig) =  0.0    ! v10         	  
          end if
          if(cam_in(j)%landfrac(i).gt.0) then
          c2x_c%rAttr(index_ca2x_Sca_ust  ,ig) =  cam_in(j)%fv(i)! ust          
          else
          c2x_c%rAttr(index_ca2x_Sca_ust  ,ig) = cam_in(j)%ustar(i)
          endif
          c2x_c%rAttr(index_ca2x_Sca_rmol  ,ig) =  cam_in(j)%ram1(i)! rmol          
          c2x_c%rAttr(index_ca2x_Sca_pblh  ,ig) =  cam_out(j)%pblh(i)! pblh          
          c2x_c%rAttr(index_ca2x_Sca_raincv  ,ig) =  cam_out(j)%precc(i)! raincv 
          c2x_c%rAttr(index_ca2x_Sca_rainncv  ,ig) =  cam_out(j)%precl(i)! rainncv          
          c2x_c%rAttr(index_ca2x_Sca_swdown  ,ig) =  cam_out(j)%soll(i)+cam_out(j)%sols(i)+cam_out(j)%solld(i)+cam_out(j)%solsd(i)!swdown          
          c2x_c%rAttr(index_ca2x_Sca_clflo  ,ig) =  cam_out(j)%clflo(i)!clflo
          c2x_c%rAttr(index_ca2x_Sca_clfmi  ,ig) =  cam_out(j)%clfmi(i)!clfmi
          c2x_c%rAttr(index_ca2x_Sca_clfhi  ,ig) =  cam_out(j)%clfhi(i)!clfhi
          
          do k = 1, num_soil_layers
          c2x_c%rAttr(index_ca2x_Sca_soildepth(k),ig) = cam_in(j)%soildepth(i,k)   ! for soil layer information
          c2x_c%rAttr(index_ca2x_Sca_soilthick(k),ig) = cam_in(j)%soilthick(i,k)
          c2x_c%rAttr(index_ca2x_Sca_soilt(k),ig) = cam_in(j)%soilt(i,k)
          c2x_c%rAttr(index_ca2x_Sca_soilm(k),ig) = cam_in(j)%soilm(i,k)
          enddo
       ig=ig+1
      end do 
    end do

  end subroutine atm_export_geatm

!===============================================================================
	
  subroutine atm_SetgsMap_geatm( mpicom_atm, ATMID, GSMap_cc )  !juanxiong he
    use phys_grid, only : get_nlcols_p
    !-------------------------------------------------------------------
    !
    ! Arguments
    !
    integer        , intent(in)  :: mpicom_atm
    integer        , intent(in)  :: ATMID
    type(mct_gsMap), intent(out) :: GSMap_cc
    !
    ! Local variables
    !
    integer, allocatable :: gindex(:)
    integer :: i, k, n, c, ncols, sizebuf, nlcols
    integer :: ier            ! error status
    !-------------------------------------------------------------------

    ! Build the atmosphere grid numbering for MCT
    ! NOTE:  Numbering scheme is: West to East, bottom to top, and South to North
    ! starting at south pole.  Should be the same as what's used in SCRIP
    
    ! Determine global seg map

    sizebuf=0
    do c = begchunk, endchunk
       ncols = get_ncols_p(c)
        do i = 1,ncols
          sizebuf = sizebuf+1
        end do
    end do

    allocate(gindex(sizebuf))

    n=0
    do c = begchunk, endchunk
       ncols = get_ncols_p(c)
        do i = 1,ncols
          n=n+1
          gindex(n) = get_gcol_p(c,i)  !global index
        end do
    end do

    nlcols = get_nlcols_p()
    call mct_gsMap_init( gsMap_cc, gindex, mpicom_atm, ATMID, nlcols, ngcols)

    deallocate(gindex)

  end subroutine atm_SetgsMap_geatm

!===============================================================================

  subroutine atm_domain_geatm( lsize, gsMap_c, dom_c )  !juanxiong he
    !-------------------------------------------------------------------
    !
    ! Arguments
    !
    integer        , intent(in)   :: lsize
    type(mct_gsMap), intent(in)   :: gsMap_c
    type(mct_ggrid), intent(inout):: dom_c  
    !
    ! Local Variables
    !
    integer  :: n,i,k,c,ncols           ! indices	
    real(r8) :: lats(pcols)           ! array of chunk latitudes
    real(r8) :: lons(pcols)           ! array of chunk longitude
    real(r8) :: area(pcols)           ! area in radians squared for each grid point
    real(r8), pointer  :: data(:)     ! temporary
    integer , pointer  :: idata(:)    ! temporary
    real(r8), parameter:: radtodeg = 180.0_r8/SHR_CONST_PI
    !-------------------------------------------------------------------
    !
    ! Initialize mct atm domain
    !
    call mct_gGrid_init( GGrid=dom_c, CoordChars=trim(seq_flds_dom_coord), OtherChars=trim(seq_flds_dom_other), lsize=lsize )
    !
    ! Allocate memory
    !
    allocate(data(lsize))
    !
    ! Initialize attribute vector with special value
    !
    call mct_gsMap_orderedPoints(gsMap_c, iam, idata)
    call mct_gGrid_importIAttr(dom_c,'GlobGridNum',idata,lsize)
    !
    ! Determine domain (numbering scheme is: West to East and South to North to South pole)
    ! Initialize attribute vector with special value
    !
    data(:) = -9999.0_R8 
    call mct_gGrid_importRAttr(dom_c,"lat"  ,data,lsize) 
    call mct_gGrid_importRAttr(dom_c,"lon"  ,data,lsize) 
    call mct_gGrid_importRAttr(dom_c,"area" ,data,lsize) 
    call mct_gGrid_importRAttr(dom_c,"aream",data,lsize) 
    data(:) = 0.0_R8     
    call mct_gGrid_importRAttr(dom_c,"mask" ,data,lsize) 
    data(:) = 1.0_R8
    call mct_gGrid_importRAttr(dom_c,"frac" ,data,lsize)
    !
    ! Fill in correct values for domain components
    !
    n=0
    do c = begchunk, endchunk
       ncols = get_ncols_p(c)
       call get_rlat_all_p(c, ncols, lats)
       do i=1,ncols
          n = n+1
          data(n) = lats(i)*radtodeg
       end do
    end do
    call mct_gGrid_importRAttr(dom_c,"lat",data,lsize) 

    n=0
    do c = begchunk, endchunk
       ncols = get_ncols_p(c)
       call get_rlon_all_p(c, ncols, lons)
       do i=1,ncols
          n = n+1
          data(n) = lons(i)*radtodeg
       end do
    end do
    call mct_gGrid_importRAttr(dom_c,"lon",data,lsize) 

    n=0
    do c = begchunk, endchunk
       ncols = get_ncols_p(c)
       call get_area_all_p(c, ncols, area)
       do i=1,ncols
          n = n+1
          data(n) = area(i) 
       end do
    end do
    call mct_gGrid_importRAttr(dom_c,"area",data,lsize) 

    n=0
    do c = begchunk, endchunk
       ncols = get_ncols_p(c)
       do i=1,ncols
          n = n+1
          data(n) = 1._r8 ! mask
       end do
    end do
    call mct_gGrid_importRAttr(dom_c,"mask"   ,data,lsize) 
    deallocate(data)

  end subroutine atm_domain_geatm 
!
!===========================================================================================
!
  subroutine atm_read_gearest_mct( EClock, cdata_c, x2c_c, c2x_c)
    use cam_pio_utils
    !-----------------------------------------------------------------------
    !
    ! Arguments
    !
    type(ESMF_Clock),intent(in)    :: EClock
    type(seq_cdata), intent(inout) :: cdata_c
    type(mct_aVect), intent(inout) :: x2c_c
    type(mct_aVect), intent(inout) :: c2x_c
    ! 
    ! Local variables
    !
    integer         :: npts         ! array size
    integer         :: rcode        ! return error code
    type(mct_aVect) :: gData        ! global/gathered bundle data
    integer         :: yr_spec      ! Current year
    integer         :: mon_spec     ! Current month
    integer         :: day_spec     ! Current day
    integer         :: sec_spec     ! Current time of day (sec)
    !-----------------------------------------------------------------------
    !
    ! Determine and open surface restart dataset
    !
    integer, pointer :: dof(:)
    integer :: lnx, nf_x2c, nf_c2x, k
    real(r8), allocatable :: tmp(:)
    type(file_desc_t) :: file
    type(io_desc_t) :: iodesc
    type(var_desc_t) :: varid
    character(CL)    :: itemc       ! string converted to char
    type(mct_string) :: mstring     ! mct char type

    call seq_timemgr_EClockGetData( EClock, curr_yr=yr_spec,curr_mon=mon_spec, &
         curr_day=day_spec, curr_tod=sec_spec ) 
    fname_srf_camcam = interpret_filename_spec( wrsfilename_spec_camcam, case=get_restcase(), &
         yr_spec=yr_spec, mon_spec=mon_spec, day_spec=day_spec, sec_spec= sec_spec )
    pname_srf_camcam = trim(get_restartdir() )//fname_srf_camcam
    call getfil(pname_srf_camcam, fname_srf_camcam)
    
    call cam_pio_openfile(File, fname_srf_camcam, 0)
    call mct_gsmap_OrderedPoints(cdata_c%gsmap, iam, Dof)
    lnx = mct_gsmap_gsize(cdata_c%gsmap)
    call pio_initdecomp(pio_subsystem, pio_double, (/lnx/), dof, iodesc)
    allocate(tmp(size(dof)))
    deallocate(dof)
    
    nf_x2c = mct_aVect_nRattr(x2c_c)

    do k=1,nf_x2c
       call mct_aVect_getRList(mstring,k,x2c_c)
       itemc = mct_string_toChar(mstring)
       call mct_string_clean(mstring)

       call pio_seterrorhandling(File, pio_bcast_error)
       rcode = pio_inq_varid(File,'x2c_'//trim(itemc) ,varid)
       if (rcode == pio_noerr) then
          call pio_read_darray(File, varid, iodesc, tmp, rcode)
          x2c_c%rattr(k,:) = tmp(:)
       else
	  if (masterproc) then
             write(iulog,*)'srfrest warning: field ',trim(itemc),' is not on restart file'
             write(iulog,*)'for backwards compatibility will set it to 0'
          end if
          x2c_c%rattr(k,:) = 0._r8
       end if
       call pio_seterrorhandling(File, pio_internal_error)
    end do

    nf_c2x = mct_aVect_nRattr(c2x_c)

    do k=1,nf_c2x
       call mct_aVect_getRList(mstring,k,c2x_c)
       itemc = mct_string_toChar(mstring)
       call mct_string_clean(mstring)

       rcode = pio_inq_varid(File,'c2x_'//trim(itemc) ,varid)
       call pio_read_darray(File, varid, iodesc, tmp, rcode)
       c2x_c%rattr(k,:) = tmp(:)
    end do

    call pio_freedecomp(File,iodesc)
    call pio_closefile(File)
    deallocate(tmp)

  end subroutine atm_read_gearest_mct
!
!===========================================================================================
!
  subroutine atm_write_gearest_mct( cdata_c, x2c_c, c2x_c, & 
       yr_spec, mon_spec, day_spec, sec_spec)
    use cam_pio_utils
    !-----------------------------------------------------------------------
    !
    ! Arguments
    !
    type(seq_cdata), intent(in) :: cdata_c
    type(mct_aVect), intent(in) :: x2c_c
    type(mct_aVect), intent(in) :: c2x_c
    integer        , intent(in) :: yr_spec         ! Simulation year
    integer        , intent(in) :: mon_spec        ! Simulation month
    integer        , intent(in) :: day_spec        ! Simulation day
    integer        , intent(in) :: sec_spec        ! Seconds into current simulation day
    !
    ! Local variables
    !
    integer         :: rcode        ! return error code
    type(mct_aVect) :: gData        ! global/gathered bundle data
    !-----------------------------------------------------------------------
    !
    ! Determine and open surface restart dataset
    !

    integer, pointer :: dof(:)
    integer :: nf_x2c, nf_c2x, lnx, dimid(1), k
    type(file_desc_t) :: file
    type(var_desc_t), pointer :: varid_x2c(:), varid_c2x(:)
    type(io_desc_t)  :: iodesc
    character(CL)    :: itemc       ! string converted to char
    type(mct_string) :: mstring     ! mct char type


    fname_srf_camcam = interpret_filename_spec( wrsfilename_spec_camcam, &
         yr_spec=yr_spec, mon_spec=mon_spec, day_spec=day_spec, sec_spec= sec_spec )
    call cam_pio_createfile(File, fname_srf_camcam, 0)

    call mct_gsmap_OrderedPoints(cdata_c%gsmap, iam, Dof)
    lnx = mct_gsmap_gsize(cdata_c%gsmap)
    call pio_initdecomp(pio_subsystem, pio_double, (/lnx/), dof, iodesc)

    deallocate(dof)
    
    nf_x2c = mct_aVect_nRattr(x2c_c)
    allocate(varid_x2c(nf_x2c))
    
    rcode = pio_def_dim(File,'x2c_nx',lnx,dimid(1))
    do k = 1,nf_x2c
       call mct_aVect_getRList(mstring,k,x2c_c)
       itemc = mct_string_toChar(mstring)
       call mct_string_clean(mstring)
       rcode = pio_def_var(File,'x2c_'//trim(itemc),PIO_DOUBLE,dimid,varid_x2c(k))
       rcode = pio_put_att(File,varid_x2c(k),"_fillvalue",fillvalue)
    enddo

    nf_c2x = mct_aVect_nRattr(c2x_c)
    allocate(varid_c2x(nf_c2x))
    
    rcode = pio_def_dim(File,'c2x_nx',lnx,dimid(1))
    do k = 1,nf_c2x
       call mct_aVect_getRList(mstring,k,c2x_c)
       itemc = mct_string_toChar(mstring)
       call mct_string_clean(mstring)
       rcode = PIO_def_var(File,'c2x_'//trim(itemc),PIO_DOUBLE,dimid,varid_c2x(k))
       rcode = PIO_put_att(File,varid_c2x(k),"_fillvalue",fillvalue)
    enddo

    rcode = pio_enddef(File)  ! don't check return code, might be enddef already

    do k=1,nf_x2c
       call pio_write_darray(File, varid_x2c(k), iodesc, x2c_c%rattr(k,:), rcode)
    end do

    do k=1,nf_c2x
       call pio_write_darray(File, varid_c2x(k), iodesc, c2x_c%rattr(k,:), rcode)       
    end do

    deallocate(varid_x2c, varid_c2x)

    call pio_freedecomp(File,iodesc)
    call pio_closefile(file)


  end subroutine atm_write_gearest_mct
!
!===========================================================================================
!
  subroutine gea_read_srfrest_mct( EClock, cdata_c, x2c_c, c2x_c)
    use cam_pio_utils
    !-----------------------------------------------------------------------
    !
    ! Arguments
    !
    type(ESMF_Clock),intent(in)    :: EClock
    type(seq_cdata), intent(inout) :: cdata_c
    type(mct_aVect), intent(inout) :: x2c_c
    type(mct_aVect), intent(inout) :: c2x_c
    ! 
    ! Local variables
    !
    integer         :: npts         ! array size
    integer         :: rcode        ! return error code
    type(mct_aVect) :: gData        ! global/gathered bundle data
    integer         :: yr_spec      ! Current year
    integer         :: mon_spec     ! Current month
    integer         :: day_spec     ! Current day
    integer         :: sec_spec     ! Current time of day (sec)
    !-----------------------------------------------------------------------
    !
    ! Determine and open surface restart dataset
    !
    integer, pointer :: dof(:)
    integer :: lnx, nf_x2c, nf_c2x, k
    real(r8), allocatable :: tmp(:)
    type(file_desc_t) :: file
    type(io_desc_t) :: iodesc
    type(var_desc_t) :: varid
    character(CL)    :: itemc       ! string converted to char
    type(mct_string) :: mstring     ! mct char type

    call seq_timemgr_EClockGetData( EClock, curr_yr=yr_spec,curr_mon=mon_spec, &
         curr_day=day_spec, curr_tod=sec_spec ) 
    fname_srf_gea = interpret_filename_spec( wrsfilename_spec_gea, case=get_restcase(), &
         yr_spec=yr_spec, mon_spec=mon_spec, day_spec=day_spec, sec_spec= sec_spec )
    pname_srf_gea = trim(get_restartdir() )//fname_srf_gea
    call getfil(pname_srf_gea, fname_srf_gea)
    
    call cam_pio_openfile(File, fname_srf_gea, 0)
    call mct_gsmap_OrderedPoints(cdata_c%gsmap, iam, Dof)
    lnx = mct_gsmap_gsize(cdata_c%gsmap)
    call pio_initdecomp(pio_subsystem, pio_double, (/lnx/), dof, iodesc)
    allocate(tmp(size(dof)))
    deallocate(dof)
    
    nf_x2c = mct_aVect_nRattr(x2c_c)

    do k=1,nf_x2c
       call mct_aVect_getRList(mstring,k,x2c_c)
       itemc = mct_string_toChar(mstring)
       call mct_string_clean(mstring)

       call pio_seterrorhandling(File, pio_bcast_error)
       rcode = pio_inq_varid(File,'x2c_'//trim(itemc) ,varid)
       if (rcode == pio_noerr) then
          call pio_read_darray(File, varid, iodesc, tmp, rcode)
          x2c_c%rattr(k,:) = tmp(:)
       else
	  if (masterproc) then
             write(iulog,*)'srfrest warning: field ',trim(itemc),' is not on restart file'
             write(iulog,*)'for backwards compatibility will set it to 0'
          end if
          x2c_c%rattr(k,:) = 0._r8
       end if
       call pio_seterrorhandling(File, pio_internal_error)
    end do

    nf_c2x = mct_aVect_nRattr(c2x_c)

    do k=1,nf_c2x
       call mct_aVect_getRList(mstring,k,c2x_c)
       itemc = mct_string_toChar(mstring)
       call mct_string_clean(mstring)

       rcode = pio_inq_varid(File,'c2x_'//trim(itemc) ,varid)
       call pio_read_darray(File, varid, iodesc, tmp, rcode)
       c2x_c%rattr(k,:) = tmp(:)
    end do

    call pio_freedecomp(File,iodesc)
    call pio_closefile(File)
    deallocate(tmp)

  end subroutine gea_read_srfrest_mct
!
!===========================================================================================
!
  subroutine gea_write_srfrest_mct( cdata_c, x2c_c, c2x_c, & 
       yr_spec, mon_spec, day_spec, sec_spec)
    use cam_pio_utils
    !-----------------------------------------------------------------------
    !
    ! Arguments
    !
    type(seq_cdata), intent(in) :: cdata_c
    type(mct_aVect), intent(in) :: x2c_c
    type(mct_aVect), intent(in) :: c2x_c
    integer        , intent(in) :: yr_spec         ! Simulation year
    integer        , intent(in) :: mon_spec        ! Simulation month
    integer        , intent(in) :: day_spec        ! Simulation day
    integer        , intent(in) :: sec_spec        ! Seconds into current simulation day
    !
    ! Local variables
    !
    integer         :: rcode        ! return error code
    type(mct_aVect) :: gData        ! global/gathered bundle data
    !-----------------------------------------------------------------------
    !
    ! Determine and open surface restart dataset
    !

    integer, pointer :: dof(:)
    integer :: nf_x2c, nf_c2x, lnx, dimid(1), k
    type(file_desc_t) :: file
    type(var_desc_t), pointer :: varid_x2c(:), varid_c2x(:)
    type(io_desc_t)  :: iodesc
    character(CL)    :: itemc       ! string converted to char
    type(mct_string) :: mstring     ! mct char type


    fname_srf_gea = interpret_filename_spec( wrsfilename_spec_gea, &
         yr_spec=yr_spec, mon_spec=mon_spec, day_spec=day_spec, sec_spec= sec_spec )
    call cam_pio_createfile(File, fname_srf_gea, 0)

    call mct_gsmap_OrderedPoints(cdata_c%gsmap, iam, Dof)
    lnx = mct_gsmap_gsize(cdata_c%gsmap)
    call pio_initdecomp(pio_subsystem, pio_double, (/lnx/), dof, iodesc)

    deallocate(dof)
    
    nf_x2c = mct_aVect_nRattr(x2c_c)
    allocate(varid_x2c(nf_x2c))
    
    rcode = pio_def_dim(File,'x2c_nx',lnx,dimid(1))
    do k = 1,nf_x2c
       call mct_aVect_getRList(mstring,k,x2c_c)
       itemc = mct_string_toChar(mstring)
       call mct_string_clean(mstring)
       rcode = pio_def_var(File,'x2c_'//trim(itemc),PIO_DOUBLE,dimid,varid_x2c(k))
       rcode = pio_put_att(File,varid_x2c(k),"_fillvalue",fillvalue)
    enddo

    nf_c2x = mct_aVect_nRattr(c2x_c)
    allocate(varid_c2x(nf_c2x))
    
    rcode = pio_def_dim(File,'c2x_nx',lnx,dimid(1))
    do k = 1,nf_c2x
       call mct_aVect_getRList(mstring,k,c2x_c)
       itemc = mct_string_toChar(mstring)
       call mct_string_clean(mstring)
       rcode = PIO_def_var(File,'c2x_'//trim(itemc),PIO_DOUBLE,dimid,varid_c2x(k))
       rcode = PIO_put_att(File,varid_c2x(k),"_fillvalue",fillvalue)
    enddo

    rcode = pio_enddef(File)  ! don't check return code, might be enddef already

    do k=1,nf_x2c
       call pio_write_darray(File, varid_x2c(k), iodesc, x2c_c%rattr(k,:), rcode)
    end do

    do k=1,nf_c2x
       call pio_write_darray(File, varid_c2x(k), iodesc, c2x_c%rattr(k,:), rcode)       
    end do

    deallocate(varid_x2c, varid_c2x)

    call pio_freedecomp(File,iodesc)
    call pio_closefile(file)


  end subroutine gea_write_srfrest_mct
!===========================================================================================
!
  subroutine gea_read_srfrest_mct_app( EClock, cdata_c, x2c_c, c2x_c)
    use cam_pio_utils
    !-----------------------------------------------------------------------
    !
    ! Arguments
    !
    type(ESMF_Clock),intent(in)    :: EClock
    type(seq_cdata), intent(inout) :: cdata_c
    type(mct_aVect), intent(inout) :: x2c_c
    type(mct_aVect), intent(inout) :: c2x_c
    ! 
    ! Local variables
    !
    integer         :: npts         ! array size
    integer         :: rcode        ! return error code
    type(mct_aVect) :: gData        ! global/gathered bundle data
    integer         :: yr_spec      ! Current year
    integer         :: mon_spec     ! Current month
    integer         :: day_spec     ! Current day
    integer         :: sec_spec     ! Current time of day (sec)
    !-----------------------------------------------------------------------
    !
    ! Determine and open surface restart dataset
    !
    integer, pointer :: dof(:)
    integer :: lnx, nf_x2c, nf_c2x, k
    real(r8), allocatable :: tmp(:)
    type(file_desc_t) :: file
    type(io_desc_t) :: iodesc
    type(var_desc_t) :: varid
    character(CL)    :: itemc       ! string converted to char
    type(mct_string) :: mstring     ! mct char type

    call seq_timemgr_EClockGetData( EClock, curr_yr=yr_spec,curr_mon=mon_spec, &
         curr_day=day_spec, curr_tod=sec_spec )
    fname_srf_gea = interpret_filename_spec( wrsfilename_spec_geaold, case=get_restcase(), &
         yr_spec=yr_spec, mon_spec=mon_spec, day_spec=day_spec, sec_spec= sec_spec )
    pname_srf_gea = trim(get_restartdir() )//fname_srf_gea
    call getfil(pname_srf_gea, fname_srf_gea)

    call cam_pio_openfile(File, fname_srf_gea, 0)
    call mct_gsmap_OrderedPoints(cdata_c%gsmap, iam, Dof)
    lnx = mct_gsmap_gsize(cdata_c%gsmap)
    call pio_initdecomp(pio_subsystem, pio_double, (/lnx/), dof, iodesc)
    allocate(tmp(size(dof)))
    deallocate(dof)

    nf_x2c = mct_aVect_nRattr(x2c_c)

    do k=1,nf_x2c
       call mct_aVect_getRList(mstring,k,x2c_c)
       itemc = mct_string_toChar(mstring)
       call mct_string_clean(mstring)

       call pio_seterrorhandling(File, pio_bcast_error)
       rcode = pio_inq_varid(File,'x2c_'//trim(itemc) ,varid)
       if (rcode == pio_noerr) then
          call pio_read_darray(File, varid, iodesc, tmp, rcode)
          x2c_c%rattr(k,:) = tmp(:)
       else
          if (masterproc) then
             write(iulog,*)'srfrest warning: field ',trim(itemc),' is not on restart file'
             write(iulog,*)'for backwards compatibility will set it to 0'
          end if
          x2c_c%rattr(k,:) = 0._r8
       end if
       call pio_seterrorhandling(File, pio_internal_error)
    end do

    nf_c2x = mct_aVect_nRattr(c2x_c)

    do k=1,nf_c2x
       call mct_aVect_getRList(mstring,k,c2x_c)
       itemc = mct_string_toChar(mstring)
       call mct_string_clean(mstring)

       rcode = pio_inq_varid(File,'c2x_'//trim(itemc) ,varid)
       call pio_read_darray(File, varid, iodesc, tmp, rcode)
       c2x_c%rattr(k,:) = tmp(:)
    end do

    call pio_freedecomp(File,iodesc)
    call pio_closefile(File)
    deallocate(tmp)

  end subroutine gea_read_srfrest_mct_app
!
!===========================================================================================
!
  subroutine gea_write_srfrest_mct_app( cdata_c, x2c_c, c2x_c, &
       yr_spec, mon_spec, day_spec, sec_spec)
    use cam_pio_utils
    !-----------------------------------------------------------------------
    !
    ! Arguments
    !
    type(seq_cdata), intent(in) :: cdata_c
    type(mct_aVect), intent(in) :: x2c_c
    type(mct_aVect), intent(in) :: c2x_c
    integer        , intent(in) :: yr_spec         ! Simulation year
    integer        , intent(in) :: mon_spec        ! Simulation month
    integer        , intent(in) :: day_spec        ! Simulation day
    integer        , intent(in) :: sec_spec        ! Seconds into current simulation day
    !
    ! Local variables
    !
    integer         :: rcode        ! return error code
    type(mct_aVect) :: gData        ! global/gathered bundle data
    !-----------------------------------------------------------------------
    !
    ! Determine and open surface restart dataset
    !

    integer, pointer :: dof(:)
    integer :: nf_x2c, nf_c2x, lnx, dimid(1), k
    type(file_desc_t) :: file
    type(var_desc_t), pointer :: varid_x2c(:), varid_c2x(:)
    type(io_desc_t)  :: iodesc
    character(CL)    :: itemc       ! string converted to char
    type(mct_string) :: mstring     ! mct char type

    fname_srf_gea = interpret_filename_spec( wrsfilename_spec_geaold, &
         yr_spec=yr_spec, mon_spec=mon_spec, day_spec=day_spec, sec_spec= sec_spec )
    call cam_pio_createfile(File, fname_srf_gea, 0)

    call mct_gsmap_OrderedPoints(cdata_c%gsmap, iam, Dof)
    lnx = mct_gsmap_gsize(cdata_c%gsmap)
    call pio_initdecomp(pio_subsystem, pio_double, (/lnx/), dof, iodesc)

    deallocate(dof)

    nf_x2c = mct_aVect_nRattr(x2c_c)
    allocate(varid_x2c(nf_x2c))

    rcode = pio_def_dim(File,'x2c_nx',lnx,dimid(1))
    do k = 1,nf_x2c
       call mct_aVect_getRList(mstring,k,x2c_c)
       itemc = mct_string_toChar(mstring)
       call mct_string_clean(mstring)
       rcode = pio_def_var(File,'x2c_'//trim(itemc),PIO_DOUBLE,dimid,varid_x2c(k))
       rcode = pio_put_att(File,varid_x2c(k),"_fillvalue",fillvalue)
    enddo

    nf_c2x = mct_aVect_nRattr(c2x_c)
    allocate(varid_c2x(nf_c2x))

    rcode = pio_def_dim(File,'c2x_nx',lnx,dimid(1))
    do k = 1,nf_c2x
       call mct_aVect_getRList(mstring,k,c2x_c)
       itemc = mct_string_toChar(mstring)
       call mct_string_clean(mstring)
       rcode = PIO_def_var(File,'c2x_'//trim(itemc),PIO_DOUBLE,dimid,varid_c2x(k))
       rcode = PIO_put_att(File,varid_c2x(k),"_fillvalue",fillvalue)
    enddo

    rcode = pio_enddef(File)  ! don't check return code, might be enddef already

    do k=1,nf_x2c
       call pio_write_darray(File, varid_x2c(k), iodesc, x2c_c%rattr(k,:), rcode)
    end do

    do k=1,nf_c2x
       call pio_write_darray(File, varid_c2x(k), iodesc, c2x_c%rattr(k,:), rcode)
    end do

    deallocate(varid_x2c, varid_c2x)

    call pio_freedecomp(File,iodesc)
    call pio_closefile(file)


  end subroutine gea_write_srfrest_mct_app
  
!===============================================================================
    subroutine pack_mgrid(mgrid, pack_buf, comm)
      !for packing mgrid data
      implicit none
	  include "mpif.h"
	  !----arguments-------
      type(metgrid_info_atm), intent(in) :: mgrid
!	  character(len=*)                   :: pack_buf   ! message
	  integer                            :: pack_buf   ! message
	  integer                            :: comm
	  integer                            :: ierr
	  integer                            :: size_tmp
	  integer                            :: d1, d2, d3
	  integer                            :: mgrid_attr2, mgrid_attr3, mgrid_attr_tmp_2
	  !--------------------
      integer            :: rank
	  integer*8          :: packsize
	  integer            :: pos
	  !--------------------
	  packsize = 727379960   ! can pack till 4th 2D array
!     packsize = 9727379960
!     packsize = 10000000000 ! out of bound
	  pos      = 1
	  !--------------------Pack the attributes-----
	  
	  !----integers--------------------------------
	  call MPI_Pack(mgrid%num_metgrid_levels,      1, MPI_INTEGER, pack_buf, packsize, pos, comm, ierr)
	  call MPI_Pack(mgrid%num_metgrid_soil_levels, 1, MPI_INTEGER, pack_buf, packsize, pos, comm, ierr)
	  call MPI_Pack(mgrid%iproj,                   1, MPI_INTEGER, pack_buf, packsize, pos, comm, ierr)
	  call MPI_Pack(mgrid%ids,  1, MPI_INTEGER, pack_buf, packsize, pos, comm, ierr)
	  call MPI_Pack(mgrid%ide,  1, MPI_INTEGER, pack_buf, packsize, pos, comm, ierr)
	  call MPI_Pack(mgrid%jds,  1, MPI_INTEGER, pack_buf, packsize, pos, comm, ierr)
	  call MPI_Pack(mgrid%jde,  1, MPI_INTEGER, pack_buf, packsize, pos, comm, ierr)
	  call MPI_Pack(mgrid%kds,  1, MPI_INTEGER, pack_buf, packsize, pos, comm, ierr)
	  call MPI_Pack(mgrid%kde,  1, MPI_INTEGER, pack_buf, packsize, pos, comm, ierr)
	  call MPI_Pack(mgrid%ims,  1, MPI_INTEGER, pack_buf, packsize, pos, comm, ierr)
	  call MPI_Pack(mgrid%ime,  1, MPI_INTEGER, pack_buf, packsize, pos, comm, ierr)
	  call MPI_Pack(mgrid%jms,  1, MPI_INTEGER, pack_buf, packsize, pos, comm, ierr)
	  call MPI_Pack(mgrid%jme,  1, MPI_INTEGER, pack_buf, packsize, pos, comm, ierr)
	  call MPI_Pack(mgrid%kms,  1, MPI_INTEGER, pack_buf, packsize, pos, comm, ierr)
	  call MPI_Pack(mgrid%kme,  1, MPI_INTEGER, pack_buf, packsize, pos, comm, ierr)
	  call MPI_Pack(mgrid%imsx, 1, MPI_INTEGER, pack_buf, packsize, pos, comm, ierr)
	  call MPI_Pack(mgrid%imex, 1, MPI_INTEGER, pack_buf, packsize, pos, comm, ierr)
	  call MPI_Pack(mgrid%jmsx, 1, MPI_INTEGER, pack_buf, packsize, pos, comm, ierr)
	  call MPI_Pack(mgrid%jmex, 1, MPI_INTEGER, pack_buf, packsize, pos, comm, ierr)
	  call MPI_Pack(mgrid%kmsx, 1, MPI_INTEGER, pack_buf, packsize, pos, comm, ierr)
	  call MPI_Pack(mgrid%kmex, 1, MPI_INTEGER, pack_buf, packsize, pos, comm, ierr)
	  call MPI_Pack(mgrid%imsy, 1, MPI_INTEGER, pack_buf, packsize, pos, comm, ierr)
	  call MPI_Pack(mgrid%imey, 1, MPI_INTEGER, pack_buf, packsize, pos, comm, ierr)
	  call MPI_Pack(mgrid%jmsy, 1, MPI_INTEGER, pack_buf, packsize, pos, comm, ierr)
	  call MPI_Pack(mgrid%jmey, 1, MPI_INTEGER, pack_buf, packsize, pos, comm, ierr)
	  call MPI_Pack(mgrid%kmsy, 1, MPI_INTEGER, pack_buf, packsize, pos, comm, ierr)
	  call MPI_Pack(mgrid%kmey, 1, MPI_INTEGER, pack_buf, packsize, pos, comm, ierr)
	  call MPI_Pack(mgrid%ips,  1, MPI_INTEGER, pack_buf, packsize, pos, comm, ierr)
	  call MPI_Pack(mgrid%ipe,  1, MPI_INTEGER, pack_buf, packsize, pos, comm, ierr)
	  call MPI_Pack(mgrid%jps,  1, MPI_INTEGER, pack_buf, packsize, pos, comm, ierr)
	  call MPI_Pack(mgrid%jpe,  1, MPI_INTEGER, pack_buf, packsize, pos, comm, ierr)
	  call MPI_Pack(mgrid%kps,  1, MPI_INTEGER, pack_buf, packsize, pos, comm, ierr)
	  call MPI_Pack(mgrid%kpe,  1, MPI_INTEGER, pack_buf, packsize, pos, comm, ierr)
	  call MPI_Pack(mgrid%ipsx, 1, MPI_INTEGER, pack_buf, packsize, pos, comm, ierr)
	  call MPI_Pack(mgrid%ipex, 1, MPI_INTEGER, pack_buf, packsize, pos, comm, ierr)
	  call MPI_Pack(mgrid%jpsx, 1, MPI_INTEGER, pack_buf, packsize, pos, comm, ierr)
	  call MPI_Pack(mgrid%jpex, 1, MPI_INTEGER, pack_buf, packsize, pos, comm, ierr)
	  call MPI_Pack(mgrid%kpsx, 1, MPI_INTEGER, pack_buf, packsize, pos, comm, ierr)
	  call MPI_Pack(mgrid%kpex, 1, MPI_INTEGER, pack_buf, packsize, pos, comm, ierr)
	  call MPI_Pack(mgrid%ipsy, 1, MPI_INTEGER, pack_buf, packsize, pos, comm, ierr)
	  call MPI_Pack(mgrid%ipey, 1, MPI_INTEGER, pack_buf, packsize, pos, comm, ierr)
	  call MPI_Pack(mgrid%jpsy, 1, MPI_INTEGER, pack_buf, packsize, pos, comm, ierr)
	  call MPI_Pack(mgrid%jpey, 1, MPI_INTEGER, pack_buf, packsize, pos, comm, ierr)
	  call MPI_Pack(mgrid%kpsy, 1, MPI_INTEGER, pack_buf, packsize, pos, comm, ierr)
	  call MPI_Pack(mgrid%kpey, 1, MPI_INTEGER, pack_buf, packsize, pos, comm, ierr)
	  !----2D arrays--------------------------------
	  size_tmp = 0
      d1 = size(mgrid%xlat,1)
      d2 = size(mgrid%xlat,2)
	  call MPI_Type_vector(d1, d2, d2, MPI_REAL, mgrid_attr2, ierr)
	  call MPI_Type_commit(mgrid_attr2, ierr)
!	  print *, 'before 1st pack'
	  call MPI_Pack(mgrid%xlat,    1, mgrid_attr2, pack_buf, packsize, pos, comm, ierr)
!	  print *, '1st pack'
	  call MPI_Pack(mgrid%xlon,    1, mgrid_attr2, pack_buf, packsize, pos, comm, ierr)
!	  print *, '2nd pack'
	  call MPI_Pack(mgrid%tsk,     1, mgrid_attr2, pack_buf, packsize, pos, comm, ierr)
!	  print *, '3rd pack'
!	  call MPI_Pack(mgrid%sst,     1, mgrid_attr2, pack_buf, packsize, pos, comm, ierr)
!	  print *, '4th pack'
!	  call MPI_Pack(mgrid%snow,    1, mgrid_attr2, pack_buf, packsize, pos, comm, ierr)
!	  print *, '5th pack'
!	  call MPI_Pack(mgrid%xice,    1, mgrid_attr2, pack_buf, packsize, pos, comm, ierr)
!	  print *, '6th pack'
!	  call MPI_Pack(mgrid%pslv,    1, mgrid_attr2, pack_buf, packsize, pos, comm, ierr)
!	  print *, '7th pack'
!	  call MPI_Pack(mgrid%psfc,    1, mgrid_attr2, pack_buf, packsize, pos, comm, ierr)
!	  print *, '8th pack'
!	  call MPI_Pack(mgrid%xland,   1, mgrid_attr2, pack_buf, packsize, pos, comm, ierr)
!	  print *, '9th pack'
!	  call MPI_Pack(mgrid%ht,      1, mgrid_attr2, pack_buf, packsize, pos, comm, ierr)
!	  print *, '10th pack'
!	  call MPI_Pack(mgrid%sstold,  1, mgrid_attr2, pack_buf, packsize, pos, comm, ierr)
!	  print *, '11th pack'
!	  call MPI_Pack(mgrid%xiceold, 1, mgrid_attr2, pack_buf, packsize, pos, comm, ierr)
!	  print *, '12th pack'
	  !----3D arrays--------------------------------
	  size_tmp = 0
      d1 = size(mgrid%u3d,1)
      d2 = size(mgrid%u3d,2)
	  d3 = size(mgrid%u3d,3)
	  call MPI_Type_vector(d2, d3, d3, MPI_REAL, mgrid_attr_tmp_2, ierr)
	  call MPI_Type_commit(mgrid_attr_tmp_2, ierr)
	  call MPI_Type_vector(d1, d2, d2, mgrid_attr_tmp_2, mgrid_attr3, ierr)
	  call MPI_Type_commit(mgrid_attr3, ierr)
	  
!	  call MPI_Pack(mgrid%u3d, 1, mgrid_attr3, pack_buf, packsize, pos, comm, ierr)
!	  call MPI_Pack(mgrid%v3d, 1, mgrid_attr3, pack_buf, packsize, pos, comm, ierr)
!	  call MPI_Pack(mgrid%w3d, 1, mgrid_attr3, pack_buf, packsize, pos, comm, ierr)
!	  call MPI_Pack(mgrid%q3d, 1, mgrid_attr3, pack_buf, packsize, pos, comm, ierr)
!	  call MPI_Pack(mgrid%z3d, 1, mgrid_attr3, pack_buf, packsize, pos, comm, ierr)
!	  call MPI_Pack(mgrid%t3d, 1, mgrid_attr3, pack_buf, packsize, pos, comm, ierr)
!	  call MPI_Pack(mgrid%p3d, 1, mgrid_attr3, pack_buf, packsize, pos, comm, ierr)
!	  call MPI_Pack(mgrid%rh3d,      1, mgrid_attr3, pack_buf, packsize, pos, comm, ierr)
!	  call MPI_Pack(mgrid%soilt,     1, mgrid_attr3, pack_buf, packsize, pos, comm, ierr)
!	  call MPI_Pack(mgrid%soilm,     1, mgrid_attr3, pack_buf, packsize, pos, comm, ierr)
!	  call MPI_Pack(mgrid%soildepth, 1, mgrid_attr3, pack_buf, packsize, pos, comm, ierr)
!	  call MPI_Pack(mgrid%soilthick, 1, mgrid_attr3, pack_buf, packsize, pos, comm, ierr)
	  
	  
    end subroutine pack_mgrid
!==========================================================================================
	
end module atm_comp_mct
