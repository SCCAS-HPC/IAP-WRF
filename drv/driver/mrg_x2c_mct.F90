module mrg_x2c_mct

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

  public :: mrg_x2c_init_mct
  public :: mrg_x2c_run_mct
  public :: mrg_x2c_final_mct 

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
  subroutine mrg_x2c_init_mct( cdata_c, w2x_c )

    !----------------------------------------------------------------------- 
    !
    ! Arguments
    !
    type(seq_cdata), intent(in)    :: cdata_c
    type(mct_aVect), intent(inout) :: w2x_c
    !
    ! Locals
    !
    type(mct_gsmap), pointer :: gsMap_a
    integer                  :: lsize
    integer                  :: mpicom
    !----------------------------------------------------------------------- 

    ! Get atmosphere gsMap

    call seq_cdata_setptrs(cdata_c, gsMap=gsMap_a, mpicom=mpicom)
    lsize = mct_GSMap_lsize(GSMap_a, mpicom)

    ! Initialize av for land export state on atmosphere grid/ decomp

    call mct_aVect_init(w2x_c, rList=seq_flds_w2x_fields, lsize=lsize)
    call mct_aVect_zero(w2x_c)

  end subroutine mrg_x2c_init_mct

!===========================================================================================

  subroutine mrg_x2c_run_mct( cdata_c, w2x_c, x2c_c )

    !----------------------------------------------------------------------- 
    !
    ! Arguments
    !
    type(seq_cdata), intent(in)     :: cdata_c
    type(mct_aVect), intent(in)     :: w2x_c
    type(mct_aVect), intent(inout)  :: x2c_c
    !----------------------------------------------------------------------- 
    !
    ! Local workspace
    !
    logical  :: usevector    ! use vector-friendly mct_copy
    integer  :: n, ki, kl, ko, k  ! indices
    real(r8) :: frac         ! temporary
    integer  :: lsize        ! temporary
    !-----------------------------------------------------------------------
    !
    ! Zero attribute vector
    !
    call mct_avect_zero(x2c_c)
    
    ! Copy attributes that do not need to be merged
    ! These are assumed to have the same name in 
    ! (o2x_a and x2a_a) and in (l2x_a and x2a_a), etc.
    !
#ifdef CPP_VECTOR
    usevector = .true.
#else
    usevector = .false.
#endif
    call mct_aVect_copy(aVin=w2x_c, aVout=x2c_c, vector=usevector)
    ! 
    ! Merge based on fractional cell coverage
    !	
    lsize = mct_avect_lsize(x2c_c)
    do n = 1,lsize
    
           x2c_c%rAttr(index_x2c_Sx_ps,n)  = w2x_c%rAttr(index_w2x_Sw_ps   ,n)
!           if(x2c_c%rAttr(index_x2c_Sx_ps,n).ne.0) print *,x2c_c%rAttr(index_x2c_Sx_ps,n),'mrg_x2c'
      	  do k=1,26
       	   x2c_c%rAttr(index_x2c_Sx_u3d(k),n)  = w2x_c%rAttr(index_w2x_Sw_u3d(k),n)
       	   x2c_c%rAttr(index_x2c_Sx_v3d(k),n)  = w2x_c%rAttr(index_w2x_Sw_v3d(k),n)
       	   x2c_c%rAttr(index_x2c_Sx_t3d(k),n)  = w2x_c%rAttr(index_w2x_Sw_t3d(k),n)
       	   x2c_c%rAttr(index_x2c_Sx_q3d(k),n)  = w2x_c%rAttr(index_w2x_Sw_q3d(k),n)
       	   x2c_c%rAttr(index_x2c_Fcxx_dudt(k),n)  = w2x_c%rAttr(index_w2x_Fwxw_dudt(k),n)
       	   x2c_c%rAttr(index_x2c_Fcxx_dvdt(k),n)  = w2x_c%rAttr(index_w2x_Fwxw_dvdt(k),n)
       	   x2c_c%rAttr(index_x2c_Fcxx_dtdt(k),n)  = w2x_c%rAttr(index_w2x_Fwxw_dtdt(k),n)
       	   x2c_c%rAttr(index_x2c_Fcxx_dqdt(k),n)  = w2x_c%rAttr(index_w2x_Fwxw_dqdt(k),n)
       	  enddo
       
    end do

  end subroutine mrg_x2c_run_mct
!
!===========================================================================================
!
  subroutine mrg_x2c_final_mct
    ! ******************
    ! Do nothing for now
    ! ******************
  end subroutine mrg_x2c_final_mct

end module mrg_x2c_mct
