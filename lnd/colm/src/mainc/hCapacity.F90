
  subroutine hCapacity (itypwat,lb,nl_soil,csol,porsl,wice,wliq,scv,dz,cv)

!-----------------------------------------------------------------------
! Original author : Yongjiu Dai, September 15, 1999
!
! calculation of heat capacities of snow / soil layers
! the volumetric heat capacity is calculated as a linear combination
! in terms of the volumetric fraction of the constituent phases.
!-----------------------------------------------------------------------

  use precision
  use phycon_module, only : cpice,cpliq
  implicit none

  integer, INTENT(in) :: lb           ! lower bound of array
  integer, INTENT(in) :: nl_soil      ! upper bound of array
  integer, INTENT(in) :: itypwat      ! land water type (0=soil, 1=urban, 2=wetland, 
  real(r8), INTENT(in) :: csol(1:nl_soil)  ! heat capacity of soil soilds [J/(m3 K)]
  real(r8), INTENT(in) :: porsl(1:nl_soil) ! soil porosity
  real(r8), INTENT(in) :: wice(lb:nl_soil) ! ice lens [kg/m2]
  real(r8), INTENT(in) :: wliq(lb:nl_soil) ! liqui water [kg/m2]
  real(r8), INTENT(in) :: dz(lb:nl_soil)   ! layer thickiness [m]
  real(r8), INTENT(in) :: scv              ! snow water equivalent [mm]
  real(r8), INTENT(out) :: cv(lb:nl_soil)  ! heat capacity [J/(m2 K)]

!-----------------------------------------------------------------------
! Soil heat capacity, which from de Vires (1963)

      if(itypwat<=1)then ! soil ground
         cv(1:) = csol(1:)*(1.-porsl(1:))*dz(1:) + wice(1:)*cpice + wliq(1:)*cpliq
      else               ! wet land or glacier
         cv(1:) = wice(1:)*cpice + wliq(1:)*cpliq
      endif
      if(lb==1 .AND. scv>0.) cv(1) = cv(1) + cpice*scv

! Snow heat capacity
      if(lb<1)then
        cv(:0) = cpliq*wliq(:0) + cpice*wice(:0)
      endif

 end subroutine hCapacity 
