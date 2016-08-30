
#include <define.h>

module colm_rtmVar

#ifdef RTM

   use precision, only : r8
   implicit none

!
! Define parameters for RTM river routing model
!
   integer, parameter :: rtmlon = 720       !number of rtm longitudes
   integer, parameter :: rtmlat = 360       !number of rtm latitudes

   character(len=255) :: frivinp_rtm        !RTM input data file name

!
! Rtm control variables
!
   integer :: rtm_nsteps                  !if > 1, average rtm over rtm_nsteps time steps

   real(r8), pointer :: fevpa(:)          !gridlevel(local proc) evap [mm/s]
   real(r8), pointer :: rnof(:)           !gridlevel(local proc) total runoff [mm/s]
   real(r8), pointer :: prc(:)            !gridlevel(local proc) rainc [mm/s]
   real(r8), pointer :: prl(:)            !gridlevel(local proc) rainl [mm/s]

 ! move these functions into nchistMod
 ! real(r8), pointer :: lnd_ave(:)        !average of lnd runoff(global proc)
 ! real(r8), pointer :: ocn_ave(:)        !average of ocn runoff(global proc)

 ! integer :: lnd_nac                     !number of lnd runoff accum times
 ! integer :: ocn_nac                     !number of ocn runoff accum times

#endif

end module colm_rtmVar
