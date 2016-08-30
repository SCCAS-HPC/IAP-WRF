
  subroutine FireSeason (dtime,  ivt, tref, litterag,  &
                         flam,  firelength, wf)
                        
!-----------------------------------------------------------------------
! Calculate length of fire season in a year
! called every tstep.
! CALLED FROM: CLMMAIN
!-----------------------------------------------------------------------
    use precision
    implicit none

! INTENT IN VARIABLES:
    real(r8), INTENT(in) :: dtime        ! model step size
    integer,  INTENT(in) :: ivt          ! vegetation type for this pft
    real(r8), INTENT(in) :: tref         ! 2 m height surface air temperature (Kelvin)
    real(r8), INTENT(in) :: litterag    ! above ground litter (gC/m2 veget'd area)
    real(r8), INTENT(in) :: flam         ! flammability threshold
    real(r8), INTENT(in) :: wf          ! soil water as frac. of whc for top 0.5 m
!    integer,  INTENT(in) :: nl_soil      ! number of soil layers
!    real(r8), INTENT(in) :: porsl(nl_soil) ! fraction of soil that is voids [-]
!    real(r8), INTENT(in) :: wliq(nl_soil)! liquid water (kg/m2 = mm/m2)
!    real(r8), INTENT(in) :: wice(nl_soil)! liquid ice (kg/m2)
!    real(r8), INTENT(in) :: z(nl_soil)   ! layer depth (m)
!    real(r8), INTENT(in) :: dz(nl_soil)  ! layer thickness (m)


! INTENT INOUT VARIABLES:
    real(r8), INTENT(inout):: firelength ! fire season in days

! OTHER LOCAL VARIABLES:
!    integer  :: j   ! indices
    real(r8), parameter :: PI = 3.14159265358979323846
    real(r8) :: fire_prob  ! fire probability (tsteps)
!    real(r8) porslsum      ! accumulated soil por. volume
!    real(r8) wliqsum       ! accumulated liquid water volume    
    ! Calculate the length of the fire season (in days)
    ! Calculate today's fire probability, fire_prob
    ! Assume fire is only possible when temperature is above 273.16K
    ! slevis: *wf is top 0.5 m soil water as a fraction of the whc
    !         *divide fire_prob (days) by tsteps/day to get fire_prob (tsteps)
    !         *else need daily avg tref and wf to calculate fire_prob

    ! fraction of liquid water/soil por. volume to a depth of 0.5 m.
!    porslsum = 0.
!    wliqsum  = 0.
!
!    do j = 1,nl_soil
!          if (z(j)+0.5*dz(j) <= 0.5) then
!             porslsum = porslsum + porsl(j)*dz(j)
!             wliqsum  = wliqsum  + wliq(j)/denh2o + wice(j)/denice
!          end if
!    end do
!
!       if(porslsum .gt. 0.)then
!          wf = min(wliqsum/porslsum,1.0)
!       else
!          wf = 0
!       end if

       if (tref > 273.16 .and. litterag > 0.0) then
          fire_prob = exp((-PI/4.0) * (max(0.0,wf)/flam)**2) * dtime / 86400
       else
          fire_prob = 0.0
       end if
       firelength = firelength + fire_prob  ! reset 1/yr in subroutine lpj
!print*,wf, fire_prob
  end subroutine FireSeason

