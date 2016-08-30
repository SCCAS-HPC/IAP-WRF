#include <define.h>

module vegdata

#ifdef VEGDATA

   use precision
   use spmd, only : p_master
   use paramodel, only : numpft
   use timemgr, only : get_curr_date, get_mon_ndays
   use colm_varMod, only : mlai, msai

   implicit none

contains 

   subroutine interp_vegdata(ctype,p1,p2,lai_r,sai_r)

      integer,  intent(in)  :: ctype
      integer,  intent(in)  :: p1
      integer,  intent(in)  :: p2
      real(r8), intent(out) :: lai_r(numpft+1)
      real(r8), intent(out) :: sai_r(numpft+1)

      integer p, m1, m2, yr, mon, day, sec, ndays, ndays1, ndays2
      real(r8) r

      if(ctype.ne.0) then ! non-vegetation
         lai_r(:) = 0.
         sai_r(:) = 0.
         return
      end if

      if((p2-p1).ne.numpft) stop 'Error on checking numpft in interp_vegdata'

      call get_curr_date(yr,mon,day,sec)

      call get_mon_ndays(mon,ndays)

      if(day.le.ndays/2.) then
         m1 = mon-1
         m2 = mon
         if(m1.lt.1) m1 = 12
      else
         m1 = mon
         m2 = mon+1
         if(m2.gt.12) m2 = 1
      end if

      call get_mon_ndays(m1, ndays1)
      call get_mon_ndays(m2, ndays2)

      if(m1.eq.mon) then
         r = (day-ndays1/2.)*2/(ndays1+ndays2)
      else
         r = (day+ndays1/2.)*2/(ndays1+ndays2)
      end if

      do p = p1, p2
         lai_r(p-p1+1) = mlai(m1,p)*(1-r)+mlai(m2,p)*r
         sai_r(p-p1+1) = msai(m1,p)*(1-r)+msai(m2,p)*r
      end do

   end subroutine interp_vegdata

#endif

end module vegdata

