
 subroutine dewfraction (sigf,lai,sai,dewmx,ldew,fwet,fdry)
       
!=======================================================================
! Original author: Yongjiu Dai, September 15, 1999
!
! determine fraction of foliage covered by water and
! fraction of foliage that is dry and transpiring
!
!=======================================================================

 use precision
 implicit none

  real(r8), INTENT(in) :: sigf   ! fraction of veg cover, excluding snow-covered veg [-]
  real(r8), INTENT(in) :: lai    ! leaf area index  [-]
  real(r8), INTENT(in) :: sai    ! stem area index  [-]
  real(r8), INTENT(in) :: dewmx  ! maximum allowed dew [0.1 mm]
  real(r8), INTENT(in) :: ldew   ! depth of water on foliage [kg/m2/s]

  real(r8), INTENT(out) :: fwet  ! fraction of foliage covered by water [-]
  real(r8), INTENT(out) :: fdry  ! fraction of foliage that is green and dry [-]

  real(r8) lsai                  ! lai + sai
  real(r8) dewmxi                ! inverse of maximum allowed dew [1/mm]
  real(r8) vegt                  ! sigf*lsai
!
!-----------------------------------------------------------------------
! Fwet is the fraction of all vegetation surfaces which are wet 
! including stem area which contribute to evaporation
      lsai = lai + sai
      dewmxi  = 1.0/dewmx
      vegt   =  sigf*lsai

      fwet = 0
      if(ldew > 0.) then
         fwet = ((dewmxi/vegt)*ldew)**.666666666666

! Check for maximum limit of fwet
         fwet = min(fwet,1.0)

      end if

! fdry is the fraction of lai which is dry because only leaves can 
! transpire. Adjusted for stem area which does not transpire
      fdry = (1.-fwet)*lai/lsai
!
!-----------------------------------------------------------------------
 end subroutine dewfraction
