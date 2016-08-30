
 subroutine albocean (oro, scv, coszrs, alb)

!-----------------------------------------------------------------------
!
! Compute surface albedos
!
! Computes surface albedos for direct/diffuse incident radiation for
! two spectral intervals:
!   s = 0.2-0.7 micro-meters
!   l = 0.7-5.0 micro-meters
!
! Albedos specified as follows:
!
! Ocean           Uses solar zenith angle to compute albedo for direct
!                 radiation; diffuse radiation values constant; albedo
!                 independent of spectral interval and other physical
!                 factors such as ocean surface wind speed.
!
! Ocean with      Surface albs specified; combined with overlying snow
!   sea ice       
!
! For more details , see Briegleb, Bruce P., 1992: Delta-Eddington
! Approximation for Solar Radiation in the NCAR Community Climate Model,
! Journal of Geophysical Research, Vol 97, D7, pp7603-7612).
!
! yongjiu dai and xin-zhong liang (08/01/2001)
!-----------------------------------------------------------------------

   use precision
   implicit none

!------------------------------Arguments--------------------------------

  real(r8), INTENT(in) :: oro       ! /ocean(0)/seaice(2) flag
  real(r8), INTENT(in) :: scv       ! snow water equivalent) [mm]
  real(r8), INTENT(in) :: coszrs    ! Cosine solar zenith angle

  real(r8), INTENT(out) :: alb(2,2) ! srf alb for direct (diffuse) rad 0.2-0.7 micro-ms
                                    ! Srf alb for direct (diffuse) rad 0.7-5.0 micro-ms

!---------------------------Local variables-----------------------------

  real(r8) frsnow       ! horizontal fraction of snow cover
  real(r8) snwhgt       ! physical snow height
  real(r8) rghsnw       ! roughness for horizontal snow cover fractn

  real(r8) sasdir       ! snow alb for direct rad  0.2-0.7 micro-ms
  real(r8) saldir       ! snow alb for direct rad  0.7-5.0 micro-ms
  real(r8) sasdif       ! snow alb for diffuse rad  0.2-0.7 micro-ms
  real(r8) saldif       ! snow alb for diffuse rad  0.7-5.0 micro-ms

  real(r8), parameter :: asices = 0.70 ! sea ice albedo for 0.2-0.7 micro-meters [-]
  real(r8), parameter :: asicel = 0.50 ! sea ice albedo for 0.7-5.0 micro-meters [-]
  real(r8), parameter :: asnows = 0.95 ! snow    albedo for 0.2-0.7 micro-meters [-]
  real(r8), parameter :: asnowl = 0.70 ! snow    albedo for 0.7-5.0 micro-meters 

!-----------------------------------------------------------------------
! initialize all ocean/sea ice surface albedos to zero

      alb(:,:) = 0.
      if(coszrs<=0.0) return

      if(nint(oro)==2)then
        alb(1,1) = asices
        alb(2,1) = asicel
        alb(1,2) = alb(1,1) 
        alb(2,2) = alb(2,1)
        sasdif = asnows
        saldif = asnowl

        if(scv>0.)then
          if (coszrs<0.5) then
          ! zenith angle regime 1 ( coszrs < 0.5 ).
          ! set direct snow albedos (limit to 0.98 max)
            sasdir = min(0.98,sasdif+(1.-sasdif)*0.5*(3./(1.+4.*coszrs)-1.))
            saldir = min(0.98,saldif+(1.-saldif)*0.5*(3./(1.+4.*coszrs)-1.))
          else
          ! zenith angle regime 2 ( coszrs >= 0.5 )
            sasdir = asnows
            saldir = asnowl
          end if

        ! compute both diffuse and direct total albedos
          snwhgt = 20.*scv / 1000.
          rghsnw = 0.25
          frsnow = snwhgt/(rghsnw+snwhgt)
          alb(1,1) = alb(1,1)*(1.-frsnow) + sasdir*frsnow
          alb(2,1) = alb(2,1)*(1.-frsnow) + saldir*frsnow
          alb(1,2) = alb(1,2)*(1.-frsnow) + sasdif*frsnow
          alb(2,2) = alb(2,2)*(1.-frsnow) + saldif*frsnow
        end if
      end if 

! ice-free ocean albedos function of solar zenith angle only, and
! independent of spectral interval:

      if(nint(oro)==0)then
        alb(2,1) = .026/(coszrs**1.7+.065) & 
                 + .15*(coszrs-0.1)*(coszrs-0.5)*(coszrs-1.) 
        alb(1,1) = alb(2,1) 
        alb(1,2) = 0.06
        alb(2,2) = 0.06
      end if
   
 end subroutine albocean
