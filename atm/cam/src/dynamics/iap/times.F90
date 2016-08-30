module times

!----------------------------------------------------------------------- 
! 
! Purpose: Time information for dynamic subroutines.
!          WuJianping,  April 2011
! 
!-----------------------------------------------------------------------

   use spmd_utils, only: masterproc          !zhh 2011-12-15
   use time_manager, only: get_curr_date
   implicit none
   double precision, public :: time_all
   double precision, public :: time_dyn
   double precision, public :: time_phy
   double precision, public :: time_trans_grid_pd
   double precision, public :: time_hdifus
   double precision, public :: time_dyfram
   double precision, public :: time_trans_grid_dp
   double precision, public :: time_mass_engy
   double precision, public :: time_trans_af_mass
   double precision, public :: time_nlitti
   integer :: yr, mon, day      ! year, month, and day components of a date
   integer :: ncsec             ! current time of day [seconds]

CONTAINS

!========================================================================

  subroutine times_set
    time_dyn=0.0D0
    time_phy=0.0D0
    time_trans_grid_pd=0.0D0
    time_hdifus=0.0D0
    time_dyfram=0.0D0
    time_trans_grid_dp=0.0D0
    time_mass_engy=0.0D0
    time_trans_af_mass=0.0D0
  end subroutine times_set

!========================================================================

  subroutine times_out
    if(masterproc) then
      call get_curr_date(yr, mon, day, ncsec)
       write(6,10) yr, mon, day, ncsec, time_all, time_dyn, time_phy
    endif
10  format('End of run time: ', 4I6, ' Time eclapsed(s). Tol: ', F12.1, ' Dyn: ', F12.1, &
           ' Phy: ', F12.1)
  end subroutine times_out

end module times
