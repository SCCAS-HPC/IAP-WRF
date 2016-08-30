
 subroutine snowfraction (itypwat,fveg,z0m,snowdp,scv,wt,sigf,fsno)

!=======================================================================
! Original author : Yongjiu Dai, September 15, 1999
! Provide snow cover fraction
!=======================================================================

  use precision
  implicit none

! dummy arguments
  integer , INTENT(in) :: itypwat! land water type (0=soil, 1=urban and built-up, 
                                 ! 2=wetland, 3=land ice, 4=deep lake, 5=shallow lake, 99 = ocean)
  real(r8), INTENT(in) :: fveg   ! fractional vegetation cover [-]
  real(r8), INTENT(in) :: z0m    ! aerodynamic roughness length [m]
  real(r8), INTENT(in) :: snowdp ! snow depth [m]
  real(r8), INTENT(in) :: scv    ! snow mass [kg/m2]

  real(r8), INTENT(out) :: wt    ! fraction of vegetation covered with snow [-]
  real(r8), INTENT(out) :: sigf  ! fraction of veg cover, excluding snow-covered veg [-]
  real(r8), INTENT(out) :: fsno  ! fraction of soil covered by snow [-]

  real(r8), parameter   :: zlnd = 0.01      !Roughness length for soil [m]

  real(r8) snowbd, fmelt

!-----------------------------------------------------------------------
      if(fveg > 0.001) then
! Fraction of vegetation buried (covered) by snow
         wt = 0.1*snowdp/z0m
         wt = wt/(1.+wt)

! Fraction of vegetation cover free of snow
         sigf = (1.-wt)*fveg

! Fraction of soil covered by snow
         fsno = snowdp/(0.1+snowdp)
      else
         wt = 0.
         sigf = 0.
         fsno = snowdp/(0.1+snowdp)
      endif

      if(snowdp>=1.E-6) then
         if(itypwat==1) then
            fsno = min(snowdp/0.05, 1.)
         else
            snowbd = min(800., scv/snowdp)
            fmelt  = (snowbd/100.)**1.
            ! 100 is the assumed fresh snow density; 1 is a melting factor that could be
            ! reconsidered, optimal value of 1.5 in Niu et al., 2007
            fsno = tanh(snowdp/(2.5*zlnd*fmelt))
         end if
      else
         fsno = 0.
      end if

      if(sigf < 0.001) sigf = 0.
      if(sigf > 0.999) sigf = 1.

 end subroutine snowfraction
