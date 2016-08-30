 module DGVMtimevar
!---------------------------------------------------------------------------------
!Description:
!All of the time-averaged variables and accumulated variables used in
!the DGVM are calculated in this module.
!---------------------------------------------------------------------------------

use precision
use phycon_module, only : denice, denh2o
implicit none

save


contains

subroutine cal_T10(msec,tref,tref10,tref_sum,t10,nday)

 use timemgr, only : dtime

!intent in variables:
 integer,  intent(in) :: msec        ! current model second
 real(r8), INTENT(in) :: tref        ! model simulated temperature at 2m

!intent inout variables:
 real(r8), intent(inout):: t10(10)   ! an array to keep 10-day temperature
 real(r8), intent(in):: nday         ! counting the model day


!intent out variables:
 real(r8), INTENT(out):: tref10      ! 10-day running mean temperature at 2m (K)
 real(r8), INTENT(out):: tref_sum    ! sum tref of current day


!local variables:
 integer   :: n

 tref_sum = tref_sum + tref

 if(msec.eq.0) then
         n = mod(nday,10.)
         if(n.eq.0) n = 10
       ! t10(n) = tref_sum/48
         t10(n) = tref_sum/(86400/dtime)
         tref_sum = 0.
         tref10 = sum(t10)/10
 end if

end subroutine cal_T10




subroutine cal_Tmomin(msec,tref,t_mo,t_mo_sum,t_mo_min,nday)  !added by zhq. 08.11

 use timemgr, only : dtime

!intent in variables:
 integer,  intent(in) :: msec    ! current model second (0~86400)
 real(r8), intent(in) :: tref    ! model simulated temperature at 2m
 real(r8), intent(in) :: nday    ! counting the model day
 
!intent inout variables:
 real(r8), intent(inout) :: t_mo      ! 30-day mean temperature at 2m
 real(r8), intent(inout) :: t_mo_sum  ! 30-day accumulated temperature at 2m
 real(r8), intent(inout) :: t_mo_min  ! annual minimun t_mo

 t_mo_sum = t_mo_sum + tref
 
 if((msec.eq.0).and.(mod(nday,30.).eq.0)) then
  ! t_mo = t_mo_sum/30./48.
    t_mo = t_mo_sum/30./(86400/dtime)
    t_mo_sum = 0.
    t_mo_min = min(t_mo_min,t_mo)
 endif

end subroutine cal_Tmomin    




subroutine cal_AN10(msec,assim,frmf,assimn_sum,assimn10,an10,nday)

 use timemgr, only : dtime

!intent in variables:
 integer,  intent(in) :: msec        ! current model second
 real(r8), INTENT(in) :: assim       ! model simulated assimilation rate
 real(r8), INTENT(in) :: frmf        ! model simulated leaf respiration rate

!intent inout variables:
 real(r8), intent(inout):: an10(10)  ! an array to keep 10-day temperature
 real(r8), intent(in):: nday         ! counting the model day


!intent out variables:
 real(r8), INTENT(out):: assimn10      ! 10-day running mean temperature at 2m (K)
 real(r8), INTENT(out):: assimn_sum    ! sum tref of current day


!local variables:
 integer   :: n
 real(r8)  :: assimn      ! net assimilation rate(assim-frmf)

 assimn = assim - frmf
 assimn_sum = assimn_sum + assimn

 if(msec.eq.0) then
         n = mod(nday,10.)
         if(n.eq.0) n = 10
      !  an10(n) = assimn_sum/48
         an10(n) = assimn_sum/(86400/dtime)
         assimn_sum = 0.
         assimn10 = sum(an10)/10
 end if

end subroutine cal_AN10
 



subroutine cal_prec365(prec365,prc,prl)

!intent in variables:
 real(r8), INTENT(in) :: prc
 real(r8), INTENT(in) :: prl

!intent out variables:
 real(r8), INTENT(inout):: prec365

 prec365 = prec365 + (prc + prl)

end subroutine cal_prec365




subroutine cal_GDD(agdd,tref10,T0)

!intent in variables:
 real(r8), INTENT(in):: tref10      ! 10-day running mean temperature at 2m [K]
 real(r8), INTENT(in):: T0          ! 10-day running mean temperature at 2m [K]

!intent inout variables:
 real(r8), INTENT(inout):: agdd      !the annual growing dgree days

!local variables:
 real(r8) :: dgd                     !the daily growing dgree
 
 dgd = tref10 - T0
 if(dgd .gt. 0.) then
    agdd = agdd + dgd
 else                                !reset agdd0 when tref is less then T0
!    agdd = 0.                       !deliminated by zhq. dec30,08
    dgd = 0.                         !added by zhq. dec30,08
 endif

end subroutine cal_GDD




subroutine cal_wf(nl_soil,porsl,wliq,wice,z,dz,wf)

!intent in variables:
 integer,  INTENT(in) :: nl_soil      ! number of soil layers
 real(r8), INTENT(in) :: porsl(nl_soil) ! fraction of soil that is voids [-]
 real(r8), INTENT(in) :: wliq(nl_soil)! liquid water (kg/m2 = mm/m2)
 real(r8), INTENT(in) :: wice(nl_soil)! liquid ice (kg/m2)
 real(r8), INTENT(in) :: z(nl_soil)   ! layer depth (m)
 real(r8), INTENT(in) :: dz(nl_soil)  ! layer thickness (m)

!intent inout variables:
 real(r8), INTENT(out) :: wf          ! soil water as frac. of whc for top 0.5 m

!local variables:
 real(r8) porslsum      ! accumulated soil por. volume
 real(r8) wliqsum       ! accumulated liquid water volume
 integer :: j           ! index

!fraction of liquid water/soil por. volume to a depth of 0.5 m.
 porslsum = 0.
 wliqsum  = 0.

 do j = 1,nl_soil
   if(z(j)+0.5*dz(j) <= 0.5) then
      porslsum = porslsum + porsl(j)*dz(j)
      wliqsum  = wliqsum  + wliq(j)/denh2o + wice(j)/denice
   end if
 end do

 if(porslsum .gt. 0.)then
    wf = min(wliqsum/porslsum,1.0)
 else
    wf = 0
 end if

end subroutine cal_wf



subroutine cal_tsoi25(nl_soil,z,dz,tss,tsoi25)

!intent in variables:
 integer , INTENT(in) :: nl_soil     ! number of soil layers
 real(r8), INTENT(in) :: z(1:nl_soil)! layer depth (m)
 real(r8), INTENT(in) :: dz(1:nl_soil)! layer thickness (m)
 real(r8), INTENT(in) :: tss(1:nl_soil) ! soil temperature (Kelvin)

!intent inout variables:
 real(r8), INTENT(out) :: tsoi25     ! soil temperature to 0.25 m (Kelvin)

!local variables:
 real(r8) :: tsoi       ! temporary
 real(r8) :: dep        ! temporary
 integer :: j

! Soil temperature to a depth of 0.25 m.
 tsoi = 0.
 dep  = 0.
 do j = 1, nl_soil
    if (z(j)+0.5*dz(j) <= 0.50) then
       tsoi = tsoi + tss(j)*dz(j)
       dep  = dep + dz(j)
    end if
 end do

 if (dep /= 0.) then
    tsoi25 = tsoi/dep
 else
    tsoi25 = tss(1)
 end if

end subroutine cal_tsoi25



End module DGVMtimevar
