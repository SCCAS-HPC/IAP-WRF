#include <define.h>

module timemgr

   use precision
   use colm_varctl
   use abortutils, only: endrun
   implicit none

   real(r8):: dtime                           ! time step (senconds)
   integer :: mstep                           ! model step for simulation [-]
   integer :: istep                           ! current model step

   integer :: idate(3)                        ! calendar (year, julian day, seconds)
   integer :: idate_p(3)                      ! current model calendar 

contains

   subroutine ticktime 

!=======================================================================
! Original author: Yongjiu Dai, September 15, 1999
!
!     Notes: Greenwich time
!
!=======================================================================

      integer maxday                     ! days of one year (365 or 366)

!-----------------------------------------------------------------------

      idate_p(:) = idate(:)

      idate(3) = idate(3) + nint(dtime)

      if(idate(3)>=86400)then
         idate(3) = idate(3) - 86400
#ifndef COUP_CSM
         if((mod(idate(1),4)==0 .AND. mod(idate(1),100)/=0) .OR. &
                                      mod(idate(1),400)==0)then
             maxday = 366
         else
             maxday = 365
         endif
#else
         maxday = 365
#endif
         idate(2) = idate(2) + 1

         if(idate(2)>maxday) then
            idate(1) = idate(1) + 1
            idate(2) = 1
         endif
      endif

#ifdef SPINUP
      if(spinup_curryr.lt.spinup_begdat) then
         if(idate_p(2).eq.1.and.idate_p(3).lt.dtime) then
            spinup_curryr = spinup_curryr + 1

            if(idate_p(1).gt.spinup_enddat) then
               idate_p(1) = spinup_begdat
               idate(1) = spinup_begdat
            end if
         endif
      endif
#endif

   end subroutine ticktime

   subroutine set_curr_date(idate_c)

      integer, intent(in) :: idate_c(3)

      idate = idate_c

   end subroutine set_curr_date

   subroutine get_curr_date(yr, mon, day, sec)

      integer, intent(out) :: yr
      integer, intent(out) :: mon
      integer, intent(out) :: day
      integer, intent(out) :: sec

      integer months(12)

      logical leapyear

      integer k

      leapyear = (mod(idate(1),4)==0.and.mod(idate(1),100)/=0).or.mod(idate(1),400)==0

#ifndef COUP_CSM
      if(leapyear)then
         months = (/31,60,91,121,152,182,213,244,274,305,335,366/)
      else
         months = (/31,59,90,120,151,181,212,243,273,304,334,365/)
      endif
#else
      months = (/31,59,90,120,151,181,212,243,273,304,334,365/)
#endif

      do k = 1, 12
         if(idate(2) .le. months(k)) then
            mon = k

            if(k.gt.1) then
               day = idate(2) - months(k-1)
            else
               day = idate(2)
            end if

            exit
         end if
      end do

      yr = idate(1)
      sec = idate(3)

   end subroutine get_curr_date

   subroutine get_mon_ndays(mon, ndays)

      integer, intent(in)  :: mon
      integer, intent(out) :: ndays

      integer months(12)
      logical leapyear
      integer k

      leapyear = (mod(idate(1),4)==0.and.mod(idate(1),100)/=0).or.mod(idate(1),400)==0

#ifndef COUP_CSM
      if(leapyear)then
         months = (/31,60,91,121,152,182,213,244,274,305,335,366/)
      else
         months = (/31,59,90,120,151,181,212,243,273,304,334,365/)
      endif
#else
      months = (/31,59,90,120,151,181,212,243,273,304,334,365/)
#endif

      do k = 12, 2, -1
         months(k) = months(k)-months(k-1)
      end do

      ndays = months(mon)

   end subroutine get_mon_ndays

   function get_curr_calday(offset)

! Return calendar day at end of current timestep with optional offset.
! Calendar day 1.0 = 0Z on Jan 1.

! Arguments
      integer, optional, intent(in) :: offset  ! Offset from current time in seconds.
                                            ! Positive for future times, negative 
                                            ! for previous times.
! Return value
      real(r8) :: get_curr_calday

      if(present(offset)) then
         get_curr_calday = idate(2) + float(idate(3))/86400 + float(offset)/86400
      else
         get_curr_calday = idate(2) + float(idate(3))/86400
      end if

      if (get_curr_calday.gt.366) then
          get_curr_calday = get_curr_calday-365
      end if

      if ( (get_curr_calday < 1.0) .or. (get_curr_calday > 366.0) )then
         call endrun( 'get_curr_calday'//': error get_curr_calday out of bounds' )
      end if
   
   end function get_curr_calday

   subroutine get_prev_date(yr, mon, day, sec)

      integer, intent(out) :: yr
      integer, intent(out) :: mon
      integer, intent(out) :: day
      integer, intent(out) :: sec

      integer months(12)

      logical leapyear

      integer k

      leapyear = (mod(idate_p(1),4)==0.and.mod(idate_p(1),100)/=0).or.mod(idate_p(1),400)==0

#ifndef COUP_CSM
      if(leapyear)then
         months = (/31,60,91,121,152,182,213,244,274,305,335,366/)
      else
         months = (/31,59,90,120,151,181,212,243,273,304,334,365/)
      endif
#else
      months = (/31,59,90,120,151,181,212,243,273,304,334,365/)
#endif

      do k = 1, 12
         if(idate_p(2) .le. months(k)) then
            mon = k

            if(k.gt.1) then
               day = idate_p(2) - months(k-1)
            else
               day = idate_p(2)
            end if

            exit
         end if
      end do

      yr = idate_p(1)
      sec = idate_p(3)

   end subroutine get_prev_date

   function get_curr_year()

      integer get_curr_year

      get_curr_year = idate(1) 

   end function get_curr_year

   function get_step_size()

      real(r8) get_step_size

      get_step_size = dtime

   end function get_step_size

   function get_nstep()

      integer get_nstep

      get_nstep = istep

   end function get_nstep

end module timemgr
