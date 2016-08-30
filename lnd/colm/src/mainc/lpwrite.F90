#include <define.h>
 
subroutine lpwrite(idate_p,idate,fout,frestart,cdate)

  use precision
  use colm_varctl
  use nchistMod
#ifdef SPINUP
  use timemgr, only : dtime, spinup_curryr
#else
  use timemgr, only : dtime
#endif

  implicit none

  integer, intent(in) :: idate_p(3)
  integer, intent(in) :: idate(3)
  character(LEN=255), intent(in)  :: fout 
  character(LEN=255), intent(inout) :: frestart
  character(LEN=255), intent(out) :: cdate

  logical :: leapyear
  integer :: months(0:12)
  integer year  , month  , day  , hour  , minute  , second
  integer yearp, monthp, dayp, hourp, minutep, secondp
  integer m, n, nyear


  leapyear = (mod(idate_p(1),4)==0.and.mod(idate_p(1),100)/=0).or.&
                                       mod(idate_p(1),400)==0

#ifdef COUP_CSM
  months = (/0,31,59,90,120,151,181,212,243,273,304,334,365/)
#else
  if(leapyear)then
     months = (/0,31,60,91,121,152,182,213,244,274,305,335,366/)
  else
     months = (/0,31,59,90,120,151,181,212,243,273,304,334,365/)
  end if
#endif

  monthp = 0
  do m = 1, 12
     if(idate_p(2).le.months(m))then
        monthp = m
        exit
     end if
  end do

  yearp   = idate_p(1)
  dayp    = idate_p(2) - months(monthp-1)
  hourp   = idate_p(3)/3600
  minutep = (idate_p(3)-3600*hourp)/60
  secondp = mod(idate_p(3),60)

  month = 0
  do m = 1, 12
     if(idate(2).le.months(m))then
        month = m
        exit
     end if
  end do

  year   = idate(1)
  day    = idate(2) - months(month-1)
  hour   = idate(3)/3600
  minute = (idate(3)-3600*hour)/60
  second = mod(idate(3),60)

#ifdef SPINUP
  if(spinup_curryr.lt.spinup_begdat) then
     nyear = yearp - spinup_curryr
     yearp = yearp - nyear
     year  = year  - nyear
  end if
#endif

  do n = 1, nMaxHist

     if(histArray(n)%is_valid) then
        histArray(n)%is_ready = .false.
        histArray(n)%is_newyear = .false.

        if(yearp.ne.year) histArray(n)%is_newyear = .true.

        if(histArray(n)%y_interval.gt.0 .and. &
           mod(yearp,histArray(n)%y_interval).eq.0 .and. yearp.ne.year) then
           write(cdate,'(i4.4)') yearp
           histArray(n)%fhistory = trim(fout)//'-'//trim(cdate)
           histArray(n)%is_ready = .true.
        end if

        if(histArray(n)%m_interval.gt.0 .and. &
           mod(monthp,histArray(n)%m_interval).eq.0 .and. monthp.ne.month) then
           write(cdate,'(i4.4,"-",i2.2)') yearp,monthp
           histArray(n)%fhistory = trim(fout)//'-'//trim(cdate)
           histArray(n)%is_ready = .true.
        end if

        if(histArray(n)%d_interval.gt.0 .and. &
           mod(dayp,histArray(n)%d_interval).eq.0 .and. dayp.ne.day) then
           write(cdate,'(i4.4,"-",i2.2,"-",i2.2)') yearp,monthp,dayp
           histArray(n)%fhistory = trim(fout)//'-'//trim(cdate)
           histArray(n)%is_ready = .true.
        end if

        if(histArray(n)%h_interval.gt.0 .and. &
           mod(hourp,histArray(n)%h_interval).eq.0 .and. hourp.ne.hour) then
           write(cdate,'(i4.4,"-",i2.2,"-",i2.2,"-",i2.2)') yearp,monthp,dayp,hourp
           histArray(n)%fhistory = trim(fout)//'-'//trim(cdate)
           histArray(n)%is_ready = .true.
        end if
     end if

  end do

  if(yearp.ne.year) then
     write(cdate,'(i4.4)') yearp
     frestart = trim(fout)//'-restart-'//trim(cdate)
  endif

#ifdef DAILYRST
  if(dayp.ne.day) then
     write(cdate,'(i4.4,"-",i2.2,"-",i2.2)') yearp,monthp,dayp
     frestart = trim(fout)//'-restart-'//trim(cdate)
  end if
#endif

! Used as metadata in netcdf file
  write(cdate,'(i4.4,"-",i2.2,"-",i2.2," ",i2.2,":",i2.2,":",i2.2)') &
                yearp,   monthp,  dayp,    hourp,   minutep, secondp

end subroutine lpwrite
! ------------------------------------------------------------------------
! EOP
