
 subroutine netsolar (itypwat,sigf,alb,ssun,ssha,&
                      sols,soll,solsd,solld,&
                      parsun, parsha,sabvsun,sabvsha,sabg,sabvg)
!=======================================================================
! Net solar absorbed by surface
! Original author : Yongjiu Dai, 09/15/1999; 09/11/2001
!=======================================================================

  use precision
  implicit none

! Dummy argument
  integer, INTENT(in) :: itypwat  ! land water type (99-sea)
  real(r8), dimension(1:2,1:2), INTENT(in) :: &
!       albg,   &! albedo, ground [-]
!       albv,   &! albedo, vegetation [-]
        alb,    &! averaged albedo [-]
        ssun,   &! sunlit canopy absorption for solar radiation
        ssha     ! shaded canopy absorption for solar radiation

  real(r8), INTENT(in) :: &
        sigf,   &! fraction of veg cover, excluding snow-buried veg [-]
        sols,   &! atm vis direct beam solar rad onto srf [W/m2]
        soll,   &! atm nir direct beam solar rad onto srf [W/m2]
        solsd,  &! atm vis diffuse solar rad onto srf [W/m2]
        solld    ! atm nir diffuse solar rad onto srf [W/m2]

  real(r8), INTENT(out) :: &
        parsun, &! PAR absorbed by sunlit vegetation [W/m2]
        parsha, &! PAR absorbed by shaded vegetation [W/m2]
        sabvsun,&! solar absorbed by sunlit vegetation [W/m2]
        sabvsha,&! solar absorbed by shaded vegetation [W/m2]
        sabg,   &! solar absorbed by ground  [W/m2]
        sabvg    ! solar absorbed by ground + vegetation [W/m2]

!=======================================================================

        sabvsun = 0.
        sabvsha = 0.
        parsun = 0.
        parsha = 0.

        sabg  = 0.
        sabvg = 0.

        if(sols+soll+solsd+solld>0.)then
           if(itypwat<4)then        !non lake and ocean
           ! Radiative fluxes onto surface
              parsun  = ssun(1,1)*sols + ssun(1,2)*solsd
              parsha  = ssha(1,1)*sols + ssha(1,2)*solsd
              sabvsun = sols*ssun(1,1) + soll*ssun(2,1) &
                      + solsd*ssun(1,2) + solld*ssun(2,2)
              sabvsha = sols*ssha(1,1) + soll*ssha(2,1) &
                      + solsd*ssha(1,2) + solld*ssha(2,2)
              sabvsun = sigf*sabvsun
              sabvsha = sigf*sabvsha
              sabvg = sols *(1.-alb(1,1)) + soll *(1.-alb(2,1)) &
                       + solsd*(1.-alb(1,2)) + solld*(1.-alb(2,2))
              sabg = sabvg - sabvsun - sabvsha
           else                     !lake or ocean
              sabvg = sols *(1.-alb(1,1)) + soll *(1.-alb(2,1)) &
                       + solsd*(1.-alb(1,2)) + solld*(1.-alb(2,2))
              sabg = sabvg
           endif
        endif

 end subroutine netsolar
