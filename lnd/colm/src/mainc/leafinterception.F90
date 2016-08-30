
#include <define.h>

 subroutine leafinterception (dtime,dewmx,chil, &
                              prc,prl,tm,scv,sigf,lai,sai,ldew,pg)

!=======================================================================
!
! calculation of  interception and drainage of precipitation
! the treatment are based on Sellers et al. (1996)
!
! modified by Yongjiu Dai, 08/31/2002
!----------------------------------------------------------------------

  use precision
  use phycon_module, only : tfrz
  implicit none

!-----------------------Arguments---------------------------------------

  real(r8), INTENT(in) :: dtime   ! time step [second]
  real(r8), INTENT(in) :: dewmx   ! maximum dew [mm]
  real(r8), INTENT(in) :: chil    ! leaf angle distribution factor
  real(r8), INTENT(in) :: prc     ! convective precipitation rate [mm/s]
  real(r8), INTENT(in) :: prl     ! large-scale precipitation rate [mm/s]
  real(r8), INTENT(in) :: tm      ! air temperature at reference height [K]
  real(r8), INTENT(in) :: scv     ! snow mass (mm)
  real(r8), INTENT(in) :: sigf    ! fraction of veg cover, excluding snow-covered veg [-]
  real(r8), INTENT(in) :: lai     ! leaf area index [-]
  real(r8), INTENT(in) :: sai     ! stem area index [-]

  real(r8), INTENT(inout) :: ldew ! depth of water on foliage [mm]
  real(r8), INTENT(out) :: pg     ! water onto ground including canopy runoff [kg/(m2 s)]

!-----------------------Local Variables---------------------------------

  real(r8) :: satcap   ! maximum allowed water on canopy [mm]
  real(r8) :: lsai     ! sum of leaf area index and stem area index [-]
  real(r8) :: chiv     ! leaf angle distribution factor
  real(r8) :: ppc      ! convective precipitation in time-step [mm]
  real(r8) :: ppl      ! large-scale precipitation in time-step [mm]
  real(r8) :: p0       ! precipitation in time-step [mm]
  real(r8) :: fpi      ! coefficient of interception
  real(r8) :: pinf     ! interception of precipitation in time step [mm]
  real(r8) :: tti      ! direct throughfall in time step [mm]
  real(r8) :: tex      ! canopy drainage in time step [mm]
  real(r8) :: vegt     ! sigf*lsai
  real(r8) :: xs       ! proportion of the grid area where the intercepted rainfall 
		       ! plus the preexisting canopy water storage

  real(r8) :: ap, cp, bp, aa, bb, exrain, arg, thru, xsc, w
  real pcoefs (2,2)

!-----------------------End Variable List-------------------------------

 if(sigf>=0.001)then

    pcoefs(1,1) = 20.
    pcoefs(1,2) = 0.206e-8
    pcoefs(2,1) = 0.0001 
    pcoefs(2,2) = 0.9999 
    bp = 20. 

    lsai = lai + sai
    vegt = sigf*lsai
    satcap = dewmx*vegt

    p0 = (prc + prl)*dtime
    ppc = prc*dtime
    ppl = prl*dtime

    w = ldew+p0

    if(scv>0. .or. tm<tfrz) ppc = 0.
    ppl = p0 - ppc

    xsc = max(0., ldew-satcap)
    ldew = ldew - xsc

    ap = pcoefs(2,1)
    cp = pcoefs(2,2)

    if(p0>1.e-8)then
       ap = ppc/p0 * pcoefs(1,1) + ppl/p0 * pcoefs(2,1)
       cp = ppc/p0 * pcoefs(1,2) + ppl/p0 * pcoefs(2,2)

!----------------------------------------------------------------------
!      proportional saturated area (xs) and leaf drainage(tex)
!-----------------------------------------------------------------------

       chiv = chil
       if ( abs(chiv) .le. 0.01 ) chiv = 0.01
       aa = 0.5 - 0.633 * chiv - 0.33 * chiv * chiv
       bb = 0.877 * ( 1. - 2. * aa )
       exrain = aa + bb

       fpi = ( 1.-exp(-exrain*lsai) ) * sigf
       tti = p0 * ( 1.-fpi )

       xs = 1.
       if(p0*fpi>1.e-9)then
          arg = (satcap-ldew)/(p0*fpi*ap) - cp/ap
          if(arg>1.e-9)then
             xs = -1./bp * log( arg )
             xs = min( xs, 1. )
             xs = max( xs, 0. )
          endif
       endif

       tex = p0 * fpi * ( ap/bp*(1.-exp(-bp*xs)) + cp*xs ) - ( satcap - ldew ) * xs
       tex = max( tex, 0. )

     ! if(tex+tti > p0) then
       if(tex+tti-p0 > 1.e-10) then
         write(6,*) tex, tti, p0,fpi,(ap/bp*(1.-exp(-bp*xs))+cp*xs), (satcap-ldew)*xs, xs, ldew
         print *, 'tex + tti > p0 in interception code : '
			call abort
       endif

    else
       tti = 0.
       tex = 0.
    endif

!----------------------------------------------------------------------
!   total throughfall (thru) and store augmentation
!----------------------------------------------------------------------

    thru = tti + tex
    pinf = p0 - thru
    ldew = ldew + pinf

    pg = (xsc + thru) / dtime

    w = w - ldew - pg*dtime
    if(abs(w)>1.e-6)then
       write(6,*) w, ldew, pg*dtime, satcap
       print *, 'something wrong in interception code : '
	    call abort
    endif

 else

    ldew=0.
    pg = prc + prl

 endif

 end subroutine leafinterception
