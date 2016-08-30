module mrg_x2ge_mct

  use shr_kind_mod, only: r8 => shr_kind_r8
  use mct_mod
  use seq_flds_mod
  use seq_flds_indices
  use seq_comm_mct
  use seq_cdata_mod

  implicit none
  save
  private ! except

!--------------------------------------------------------------------------
! Public interfaces
!--------------------------------------------------------------------------

  public :: mrg_x2ge_init_mct
  public :: mrg_x2ge_run_mct
  public :: mrg_x2ge_final_mct 

!--------------------------------------------------------------------------
! Private interfaces
!--------------------------------------------------------------------------
!
!--------------------------------------------------------------------------
! Private data
!--------------------------------------------------------------------------
!
!===========================================================================================
contains
!===========================================================================================
!
  subroutine mrg_x2ge_init_mct( cdata_w, c2x_w )

    !----------------------------------------------------------------------- 
    !
    ! Arguments
    !
    type(seq_cdata), intent(in)    :: cdata_w
    type(mct_aVect), intent(inout) :: c2x_w
    !
    ! Locals
    !
    type(mct_gsmap), pointer :: gsMap_w
    integer                  :: lsize
    integer                  :: mpicom
    !----------------------------------------------------------------------- 

    ! Get atmosphere gsMap

    call seq_cdata_setptrs(cdata_w, gsMap=gsMap_w, mpicom=mpicom)
    lsize = mct_GSMap_lsize(GSMap_w, mpicom)

    ! Initialize av for cam export state on atmosphere grid/ decomp
    call mct_aVect_init(c2x_w, rList=seq_flds_ca2x_fields, lsize=lsize)
    call mct_aVect_zero(c2x_w)
    
  end subroutine mrg_x2ge_init_mct

!===========================================================================================

  subroutine mrg_x2ge_run_mct( cdata_w, c2x_w, x2w_w )

    !----------------------------------------------------------------------- 
    !
    ! Arguments
    !
    type(seq_cdata), intent(in)     :: cdata_w
    type(mct_aVect), intent(in)     :: c2x_w
    type(mct_aVect), intent(inout)  :: x2w_w
    !----------------------------------------------------------------------- 
    !
    ! Local workspace
    !
    logical  :: usevector    ! use vector-friendly mct_copy
    integer  :: n, k, ki, kl, ko  ! indices
    real(r8) :: frac         ! temporary
    integer  :: lsize        ! temporary
    !-----------------------------------------------------------------------
    !
    ! Zero attribute vector
    !
    call mct_avect_zero(x2w_w)
    !
    ! Copy attributes that do not need to be merged
    ! These are assumed to have the same name in 
    ! (o2x_a and x2a_a) and in (l2x_a and x2a_a), etc.
    !
#ifdef CPP_VECTOR
    usevector = .true.
#else
    usevector = .false.
#endif
    call mct_aVect_copy(aVin=c2x_w, aVout=x2w_w, vector=usevector)
    ! 
    ! Merge based on fractional cell coverage
    !	
    lsize = mct_avect_lsize(x2w_w)
    do n = 1,lsize

        do k=1,num_cam_levs
        x2w_w%rAttr(index_x2ge_Sx_z3d(k),n)  = c2x_w%rAttr(index_ca2x_Sca_z3d(k),n)             ! atm level height
        x2w_w%rAttr(index_x2ge_Sx_u3d(k),n)  = c2x_w%rAttr(index_ca2x_Sca_u3d(k),n)             ! atm level zon wind
        x2w_w%rAttr(index_x2ge_Sx_v3d(k),n)  = c2x_w%rAttr(index_ca2x_Sca_v3d(k),n)             ! atm level mer wind
        x2w_w%rAttr(index_x2ge_Sx_t3d(k),n)  = c2x_w%rAttr(index_ca2x_Sca_t3d(k),n)             ! atm level temp
        x2w_w%rAttr(index_x2ge_Sx_rh3d(k),n)  = c2x_w%rAttr(index_ca2x_Sca_rh3d(k),n)          
        x2w_w%rAttr(index_x2ge_Sx_qv3d(k),n)  = c2x_w%rAttr(index_ca2x_Sca_qv3d(k),n)          
        x2w_w%rAttr(index_x2ge_Sx_qi3d(k),n)  = c2x_w%rAttr(index_ca2x_Sca_qi3d(k),n)          
        x2w_w%rAttr(index_x2ge_Sx_qc3d(k),n)  = c2x_w%rAttr(index_ca2x_Sca_qc3d(k),n)          
        x2w_w%rAttr(index_x2ge_Sx_taucldi3d(k),n)  = c2x_w%rAttr(index_ca2x_Sca_taucldi3d(k),n)         
        x2w_w%rAttr(index_x2ge_Sx_taucldv3d(k),n)  = c2x_w%rAttr(index_ca2x_Sca_taucldv3d(k),n)          
        x2w_w%rAttr(index_x2ge_Sx_p3d(k),n)  = c2x_w%rAttr(index_ca2x_Sca_p3d(k),n)            ! atm level pressure 
        enddo
        x2w_w%rAttr(index_x2ge_Sx_ps,n)  = c2x_w%rAttr(index_ca2x_Sca_ps,n)          ! atm surface pressure
        x2w_w%rAttr(index_x2ge_Sx_phis,n)  = c2x_w%rAttr(index_ca2x_Sca_phis,n)          ! atm surface geopotential height
        x2w_w%rAttr(index_x2ge_Sx_lat,n)  = c2x_w%rAttr(index_ca2x_Sca_lat,n)          ! atm latitude
        x2w_w%rAttr(index_x2ge_Sx_lon,n)  = c2x_w%rAttr(index_ca2x_Sca_lon,n)          ! atm longittude
        x2w_w%rAttr(index_x2ge_Sx_ts,n)  = c2x_w%rAttr(index_ca2x_Sca_ts,n)          !surface temperature
        x2w_w%rAttr(index_x2ge_Sx_sst,n)  = c2x_w%rAttr(index_ca2x_Sca_sst,n)         ! sst   
        x2w_w%rAttr(index_x2ge_Sx_snowhland,n)  = c2x_w%rAttr(index_ca2x_Sca_snowhland,n)   ! snow height over land  
        x2w_w%rAttr(index_x2ge_Sx_snowhice,n)  = c2x_w%rAttr(index_ca2x_Sca_snowhice,n)  !snow height over ice    
        x2w_w%rAttr(index_x2ge_Sx_seaice,n)  = c2x_w%rAttr(index_ca2x_Sca_seaice,n)   ! seaice  
        x2w_w%rAttr(index_x2ge_Sx_ocnfrac,n)  = c2x_w%rAttr(index_ca2x_Sca_ocnfrac,n)  !ocn fraction          
        x2w_w%rAttr(index_x2ge_Sx_t2,n)  = c2x_w%rAttr(index_ca2x_Sca_t2,n) 
        x2w_w%rAttr(index_x2ge_Sx_q2,n)  = c2x_w%rAttr(index_ca2x_Sca_q2,n) 
        x2w_w%rAttr(index_x2ge_Sx_rh2,n)  = c2x_w%rAttr(index_ca2x_Sca_rh2,n) 
        x2w_w%rAttr(index_x2ge_Sx_u10,n)  = c2x_w%rAttr(index_ca2x_Sca_u10,n) 
        x2w_w%rAttr(index_x2ge_Sx_v10,n)  = c2x_w%rAttr(index_ca2x_Sca_v10,n)  
        x2w_w%rAttr(index_x2ge_Sx_ust,n)  = c2x_w%rAttr(index_ca2x_Sca_ust,n)      
        x2w_w%rAttr(index_x2ge_Sx_rmol,n)  = c2x_w%rAttr(index_ca2x_Sca_rmol,n)         
        x2w_w%rAttr(index_x2ge_Sx_pblh,n)  = c2x_w%rAttr(index_ca2x_Sca_pblh,n)         
        x2w_w%rAttr(index_x2ge_Sx_rainncv,n)  = c2x_w%rAttr(index_ca2x_Sca_rainncv,n) 
        x2w_w%rAttr(index_x2ge_Sx_raincv,n)  = c2x_w%rAttr(index_ca2x_Sca_raincv,n)     
        x2w_w%rAttr(index_x2ge_Sx_clflo,n)  = c2x_w%rAttr(index_ca2x_Sca_clflo,n)    
        x2w_w%rAttr(index_x2ge_Sx_clfmi,n)  = c2x_w%rAttr(index_ca2x_Sca_clfmi,n)      
        x2w_w%rAttr(index_x2ge_Sx_clfhi,n)  = c2x_w%rAttr(index_ca2x_Sca_clfhi,n)     
        x2w_w%rAttr(index_x2ge_Sx_swdown,n)  = c2x_w%rAttr(index_ca2x_Sca_swdown,n)   
        do k=1,num_soil_layers
        x2w_w%rAttr(index_x2ge_Sx_soildepth(k),n)  = c2x_w%rAttr(index_ca2x_Sca_soildepth(k),n)  ! soil layer depth
        x2w_w%rAttr(index_x2ge_Sx_soilthick(k),n)  = c2x_w%rAttr(index_ca2x_Sca_soilthick(k),n)  ! soil layer thickness
        x2w_w%rAttr(index_x2ge_Sx_soilt(k),n)  = c2x_w%rAttr(index_ca2x_Sca_soilt(k),n)  ! soil layer temperature
        x2w_w%rAttr(index_x2ge_Sx_soilm(k),n)  = c2x_w%rAttr(index_ca2x_Sca_soilm(k),n) ! soil layer moisture
        enddo

    end do

  end subroutine mrg_x2ge_run_mct
!
!===========================================================================================
!
  subroutine mrg_x2ge_final_mct
    ! ******************
    ! Do nothing for now
    ! ******************
  end subroutine mrg_x2ge_final_mct

end module mrg_x2ge_mct
