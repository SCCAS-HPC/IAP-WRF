module mrg_x2w_mct

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

  public :: mrg_x2w_init_mct
  public :: mrg_x2w_run_mct
  public :: mrg_x2w_final_mct 

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
  subroutine mrg_x2w_init_mct( cdata_w, c2x_w, l2x_w, o2x_w, i2x_w, xao_w )

    !----------------------------------------------------------------------- 
    !
    ! Arguments
    !
    type(seq_cdata), intent(in)    :: cdata_w
    type(mct_aVect), intent(inout) :: c2x_w
    type(mct_aVect), intent(inout), optional :: l2x_w
    type(mct_aVect), intent(inout), optional :: o2x_w
    type(mct_aVect), intent(inout), optional :: i2x_w
    type(mct_aVect), intent(inout), optional :: xao_w
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
    call mct_aVect_init(c2x_w, rList=seq_flds_c2x_fields, lsize=lsize)
    call mct_aVect_zero(c2x_w)
    
    ! Initialize av for land export state on atmosphere grid/ decomp
    if(present(l2x_w)) then	      
    call mct_aVect_init(l2x_w, rList=seq_flds_l2x_fields, lsize=lsize)
    call mct_aVect_zero(l2x_w)
    endif
    
    ! Initialize av for ocn export state on atmosphere grid/decomp
    if(present(o2x_w)) then
    call mct_aVect_init(o2x_w, rList=seq_flds_o2x_fields, lsize=lsize)
    call mct_aVect_zero(o2x_w)
    endif
    
    ! Initialize av for ice export state on atmosphere grid/decomp
    if(present(i2x_w)) then
    call mct_aVect_init(i2x_w, rList=seq_flds_i2x_fields, lsize=lsize)
    call mct_aVect_zero(i2x_w)
    endif
    
    ! Initialize av for atm/ocn flux calculation on atmosphere grid/decomp
    if(present(xao_w)) then
    call mct_aVect_init(xao_w, rList=seq_flds_xao_fields, lsize=lsize)
    call mct_aVect_zero(xao_w)
    endif
  end subroutine mrg_x2w_init_mct

!===========================================================================================

  subroutine mrg_x2w_run_mct( cdata_w, c2x_w, x2w_w, l2x_w, o2x_w, xao_w, i2x_w, fractions_w )

    !----------------------------------------------------------------------- 
    !
    ! Arguments
    !
    type(seq_cdata), intent(in)     :: cdata_w
    type(mct_aVect), intent(in)     :: c2x_w
    type(mct_aVect), intent(inout)  :: x2w_w
    type(mct_aVect), intent(in)   , optional  :: l2x_w
    type(mct_aVect), intent(in)   , optional  :: o2x_w
    type(mct_aVect), intent(in)   , optional  :: xao_w
    type(mct_aVect), intent(in)   , optional  :: i2x_w
    type(mct_aVect), intent(in)   , optional  :: fractions_w
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

        do k=1,26
        x2w_w%rAttr(index_x2w_Sx_z3d(k),n)  = c2x_w%rAttr(index_c2x_Sc_z3d(k),n)             ! atm level height
        x2w_w%rAttr(index_x2w_Sx_u3d(k),n)  = c2x_w%rAttr(index_c2x_Sc_u3d(k),n)             ! atm level zon wind
        x2w_w%rAttr(index_x2w_Sx_v3d(k),n)  = c2x_w%rAttr(index_c2x_Sc_v3d(k),n)             ! atm level mer wind
        x2w_w%rAttr(index_x2w_Sx_t3d(k),n)  = c2x_w%rAttr(index_c2x_Sc_t3d(k),n)             ! atm level temp
        x2w_w%rAttr(index_x2w_Sx_w3d(k),n)  = c2x_w%rAttr(index_c2x_Sc_w3d(k),n)          ! atm level vert wind
        x2w_w%rAttr(index_x2w_Sx_q3d(k),n)  = c2x_w%rAttr(index_c2x_Sc_q3d(k),n)          ! atm level spec hum
        x2w_w%rAttr(index_x2w_Sx_p3d(k),n)  = c2x_w%rAttr(index_c2x_Sc_p3d(k),n)            ! atm levelpressure 
        enddo
        x2w_w%rAttr(index_x2w_Sx_ps,n)  = c2x_w%rAttr(index_c2x_Sc_ps,n)          ! atm surface pressure
        x2w_w%rAttr(index_x2w_Sx_phis,n)  = c2x_w%rAttr(index_c2x_Sc_phis,n)          ! atm surface geopotential height
        x2w_w%rAttr(index_x2w_Sx_lat,n)  = c2x_w%rAttr(index_c2x_Sc_lat,n)          ! atm latitude
        x2w_w%rAttr(index_x2w_Sx_lon,n)  = c2x_w%rAttr(index_c2x_Sc_lon,n)          ! atm longittude
        x2w_w%rAttr(index_x2w_Sx_ts,n)  = c2x_w%rAttr(index_c2x_Sc_ts,n)          !surface temperature
        x2w_w%rAttr(index_x2w_Sx_sst,n)  = c2x_w%rAttr(index_c2x_Sc_sst,n)         ! sst   
        x2w_w%rAttr(index_x2w_Sx_snowhland,n)  = c2x_w%rAttr(index_c2x_Sc_snowhland,n)   ! snow height over land  
        x2w_w%rAttr(index_x2w_Sx_snowhice,n)  = c2x_w%rAttr(index_c2x_Sc_snowhice,n)  !snow height over ice    
        x2w_w%rAttr(index_x2w_Sx_seaice,n)  = c2x_w%rAttr(index_c2x_Sc_seaice,n)   ! seaice  
        x2w_w%rAttr(index_x2w_Sx_ocnfrac,n)  = c2x_w%rAttr(index_c2x_Sc_ocnfrac,n)  !ocn fraction          

        do k=1,4
        x2w_w%rAttr(index_x2w_Sx_soildepth(k),n)  = c2x_w%rAttr(index_c2x_Sc_soildepth(k),n)  ! soil layer depth
        x2w_w%rAttr(index_x2w_Sx_soilthick(k),n)  = c2x_w%rAttr(index_c2x_Sc_soilthick(k),n)  ! soil layer thickness
        x2w_w%rAttr(index_x2w_Sx_soilt(k),n)  = c2x_w%rAttr(index_c2x_Sc_soilt(k),n)  ! soil layer temperature
        x2w_w%rAttr(index_x2w_Sx_soilm(k),n)  = c2x_w%rAttr(index_c2x_Sc_soilm(k),n) ! soil layer moisture
        enddo

    end do

  end subroutine mrg_x2w_run_mct
!
!===========================================================================================
!
  subroutine mrg_x2w_final_mct
    ! ******************
    ! Do nothing for now
    ! ******************
  end subroutine mrg_x2w_final_mct

end module mrg_x2w_mct
