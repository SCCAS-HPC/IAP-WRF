module map_wrfcam_mct

!---------------------------------------------------------------------
!
! Purpose:
!
! Collect coupling routines for sequential coupling of WRF-CAM.      
!
! Author: Juanxiong He
! Revision History:
! 05/06/2010 - first version	
!
!---------------------------------------------------------------------

  use shr_kind_mod      ,only: R8 => SHR_KIND_R8, IN=>SHR_KIND_IN
  use shr_sys_mod
  use shr_const_mod
  use shr_mct_mod, only: shr_mct_sMatPInitnc
  use mct_mod

  use seq_comm_mct, only : logunit, loglevel
  use seq_cdata_mod
  use seq_flds_indices
  use seq_infodata_mod
  use m_die

  implicit none
  save
  private  ! except

!--------------------------------------------------------------------------
! Public interfaces
!--------------------------------------------------------------------------

  public :: map_cam2wrf_init_mct
  public :: map_wrf2cam_init_mct
  public :: map_cam2wrf_mct
  public :: map_wrf2cam_mct

!--------------------------------------------------------------------------
! Public data
!--------------------------------------------------------------------------

!--------------------------------------------------------------------------
! Private data
!--------------------------------------------------------------------------

  type(mct_rearr), private :: Re_cam2wrf
  type(mct_rearr), private :: Re_wrf2cam
  type(mct_sMatp), private :: sMatp_Fw2c
  type(mct_sMatp), private :: sMatp_Sw2c
  type(mct_sMatp), private :: sMatp_Fc2w
  type(mct_sMatp), private :: sMatp_Sc2w

#ifdef CPP_VECTOR
    logical :: usevector = .true.
#else
    logical :: usevector = .false.
#endif

#ifdef SYSUNICOS
    logical :: usealltoall = .true.
#else
    logical :: usealltoall = .false.
#endif
  logical, private :: samegrid_mapw2c

!=======================================================================
   contains
!=======================================================================

  subroutine map_wrf2cam_init_mct( cdata_w, cdata_c)

    !--------------------------------------------------
    ! 
    ! Arguments
    !
    type(seq_cdata),intent(in) :: cdata_w
    type(seq_cdata),intent(in) :: cdata_c
    ! 
    ! Local variables
    !
    type(seq_infodata_type), pointer :: infodata
    type(mct_gsMap), pointer :: gsMap_w           ! atm gsMap
    type(mct_gsMap), pointer :: gsMap_c           ! lnd gsMap
    type(mct_gGrid), pointer :: dom_w             ! atm domain
    type(mct_gGrid), pointer :: dom_c             ! lnd domain
    integer                  :: mpicom            ! communicator spanning atm and lnd
    integer                  :: ka, km            ! indices
    integer                  :: lsize             ! size of attribute vector
    type(mct_aVect)          :: areasrc           ! atm areas from mapping file
    type(mct_aVect)          :: areadst           ! lnd areas set to atm areas
    character(*),parameter :: subName = '(map_wrf2cam_init_mct) '
    !--------------------------------------------------

    call seq_cdata_setptrs(cdata_w, gsMap=gsMap_w, dom=dom_w)
    call seq_cdata_setptrs(cdata_c, gsMap=gsMap_c, dom=dom_c)
    call seq_cdata_setptrs(cdata_c, mpicom=mpicom, infodata=infodata)

    call seq_infodata_GetData( infodata, samegrid_wa=samegrid_mapw2c)

    if (samegrid_mapw2c) then

       call mct_rearr_init(gsMap_w, gsMap_c, mpicom, Re_wrf2cam)

    else

       call shr_mct_sMatPInitnc(sMatp_Fw2c,gsMap_w,gsMap_c,"seq_maps.rc", &
          "wrf2camFmapname:","wrf2camFmaptype:",mpicom)

       call shr_mct_sMatPInitnc(sMatp_Sw2c,gsMap_w,gsMap_c,"seq_maps.rc", &
          "wrf2camSmapname:","wrf2camSmaptype:",mpicom)

    endif

  end subroutine map_wrf2cam_init_mct

!=======================================================================

  subroutine map_cam2wrf_init_mct( cdata_c, cdata_w )

    !--------------------------------------------------
    ! 
    ! Arguments
    !
    type(seq_cdata),intent(in) :: cdata_c
    type(seq_cdata),intent(in) :: cdata_w
    ! 
    ! Local variables
    !
    type(seq_infodata_type), pointer :: infodata
    type(mct_gsMap), pointer :: gsMap_w           ! atm gsMap
    type(mct_gsMap), pointer :: gsMap_c           ! lnd gsMap
    type(mct_gGrid), pointer :: dom_c             ! lnd domain
    type(mct_gGrid), pointer :: dom_w             ! atm domain
    integer                  :: kf,iam,ierr,lsize
    integer                  :: mpicom            ! communicator spanning atm and lnd
    character(*),parameter :: subName = '(map_cam2wrf_init_mct) '
    !--------------------------------------------------

    call seq_cdata_setptrs(cdata_c, gsMap=gsMap_c, dom=dom_c)
    call seq_cdata_setptrs(cdata_w, gsMap=gsMap_w, dom=dom_w)
    call seq_cdata_setptrs(cdata_w, mpicom=mpicom,infodata=infodata)
    call mpi_comm_rank(mpicom,iam,ierr)
    
    call seq_infodata_GetData( infodata, samegrid_wa=samegrid_mapw2c)

    ! Initialize lnd -> atm mapping or rearranger

    if (samegrid_mapw2c) then

       call mct_rearr_init(gsMap_c, gsMap_w, mpicom, Re_cam2wrf)

    else

       call shr_mct_sMatPInitnc(sMatp_Fc2w, gsMap_c, gsMap_w, "seq_maps.rc", &
            "cam2wrfFmapname:", "cam2wrfFmaptype:", mpicom)

       call shr_mct_sMatPInitnc(sMatp_Sc2w, gsMap_c, gsMap_w, "seq_maps.rc", &
            "cam2wrfSmapname:", "cam2wrfSmaptype:", mpicom)

    endif

  end subroutine map_cam2wrf_init_mct

!=======================================================================

  subroutine map_wrf2cam_mct( cdata_w, av_w, cdata_c, av_c, fluxlist, statelist )

    !--------------------------------------------------
    ! 
    ! Arguments
    !
    type(seq_cdata), intent(in)  :: cdata_w
    type(mct_aVect), intent(in)  :: av_w
    type(seq_cdata), intent(in)  :: cdata_c
    type(mct_aVect), intent(out) :: av_c
    character(len=*),intent(in), optional :: statelist
    character(len=*),intent(in), optional :: fluxlist
    !
    ! Local
    ! 
    integer                :: lsize,n
    type(mct_aVect)        :: av_c_f     ! temporary flux attribute vector
    type(mct_aVect)        :: av_c_s     ! temporary state attribute vector
    character(*),parameter :: subName = '(map_wrf2cam_mct) '
    !--------------------------------------------------

    if (samegrid_mapw2c) then

       if (present(fluxlist) .or. present(statelist)) then
          if (present(fluxlist)) then
             call mct_rearr_rearrange_fldlist(av_w, av_c, Re_wrf2cam, VECTOR=usevector, &
                  ALLTOALL=usealltoall, fldlist=fluxlist)
          endif
          if (present(statelist)) then
              call mct_rearr_rearrange_fldlist(av_w, av_c, Re_wrf2cam, VECTOR=usevector, &
                  ALLTOALL=usealltoall, fldlist=statelist)
          endif
       else
          call mct_rearr_rearrange(av_w, av_c, Re_wrf2cam, VECTOR=usevector, ALLTOALL=usealltoall)
       end if

    else

       if (present(fluxlist) .or. present(statelist)) then
          if (present(fluxlist)) then
             lsize = mct_aVect_lsize(av_c)
             call mct_aVect_init (av_c_f, rlist=fluxlist , lsize=lsize)
             call mct_sMat_avMult(av_w, sMatp_Fw2c, av_c_f, VECTOR=usevector, rList=fluxlist)
             call mct_aVect_copy (aVin=av_c_f, aVout=av_c, vector=usevector)
             call mct_aVect_clean(av_c_f)
          end if
          if (present(statelist)) then
             lsize = mct_aVect_lsize(av_c)
             call mct_aVect_init (av_c_s, rlist=statelist, lsize=lsize)
             call mct_sMat_avMult(av_w, sMatp_Sw2c, av_c_s, VECTOR=usevector, rList=statelist)
             call mct_aVect_copy (aVin=av_c_s, aVout=av_c, vector=usevector)
             call mct_aVect_clean(av_c_s)
          end if
       else
          !--- default is flux mapping
          call mct_sMat_avMult(av_w, sMatp_Fw2c, av_c, VECTOR=usevector)
       endif
    endif

!    lsize = mct_avect_lsize(av_c)  ! juanxiong he, debug
!    do n = 1,lsize
!      if(av_w%rAttr(index_w2x_Sw_ps,n).ne.0) print *, av_w%rAttr(index_w2x_Sw_ps,n), samegrid_mapw2c 
!    end do
       
  end subroutine map_wrf2cam_mct

!=======================================================================

  subroutine map_cam2wrf_mct( cdata_c, av_c, cdata_w, av_w, &
                              fractions_c, fractions_w, &
                              fluxlist, statelist )

    !--------------------------------------------------
    ! 
    ! Arguments
    !
    type(seq_cdata) ,intent(in)  :: cdata_c
    type(mct_aVect) ,intent(in)  :: av_c
    type(seq_cdata) ,intent(in)  :: cdata_w
    type(mct_aVect) ,intent(out) :: av_w
    type(mct_aVect) ,intent(in), optional :: fractions_c
    type(mct_aVect) ,intent(in), optional :: fractions_w
    character(len=*),intent(in), optional :: fluxlist
    character(len=*),intent(in), optional :: statelist
    !
    ! Local
    !
    integer  :: i,j,kl,lsize,numats,ier
    real(r8) :: recip
    type(mct_aVect)          :: av_w_f     ! temporary flux attribute vector
    type(mct_aVect)          :: av_w_s     ! temporary state attribute vector
    type(mct_aVect)          :: temp       ! temporary attribute vector
    character(*),parameter :: subName = '(map_cam2wrf_mct) '
    !--------------------------------------------------

    if (samegrid_mapw2c) then

       if (present(fluxlist) .or. present(statelist)) then
          if (present(fluxlist)) then
             call mct_rearr_rearrange_fldlist(av_c, av_w, Re_cam2wrf, VECTOR=usevector, &
                  ALLTOALL=usealltoall, fldlist=fluxlist)
          endif
          if (present(statelist)) then
             call mct_rearr_rearrange_fldlist(av_c, av_w, Re_cam2wrf, VECTOR=usevector, &
                  ALLTOALL=usealltoall, fldlist=statelist)
          endif
       else
          call mct_rearr_rearrange(av_c, av_w, Re_cam2wrf, VECTOR=usevector, ALLTOALL=usealltoall)
       end if

    else
       
          if (present(fluxlist) .or. present(statelist)) then
             if (present (fluxlist)) then
                lsize = mct_aVect_lsize(av_w)
                call mct_aVect_init (av_w_f, rlist=fluxlist , lsize=lsize)
                call mct_sMat_avMult(av_c, sMatp_Fc2w, av_w_f, VECTOR=usevector, rList=fluxlist )
                call mct_aVect_copy (aVin=av_w_f, aVout=av_w, vector=usevector)
                call mct_aVect_clean(av_w_f)
             end if
             if (present(statelist)) then
                lsize = mct_aVect_lsize(av_w)
                call mct_aVect_init (av_w_s, rlist=statelist, lsize=lsize)
                call mct_sMat_avMult(av_c, sMatp_Sc2w, av_w_s, VECTOR=usevector, rList=statelist)
                call mct_aVect_copy (aVin=av_w_s, aVout=av_w, vector=usevector)
                call mct_aVect_clean(av_w_s)
             end if
          else
             ! --- default is flux mapping
             call mct_sMat_avMult(av_c, sMatp_Fc2w, av_w, VECTOR=usevector)
          endif

    endif

  end subroutine map_cam2wrf_mct

end module map_wrfcam_mct
