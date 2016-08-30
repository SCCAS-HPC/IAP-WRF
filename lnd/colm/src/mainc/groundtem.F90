
 subroutine groundtem (itypwat,lb,nl_soil,dtime,n_pft,num_filterp,filterp,wt_patch, &
                       wt_column,capr,cnfac,csol,porsl,dkmg,dkdry,dksatu, &
                       sigf,dz,z,zi,tss,wice,wliq,scv,snowdp, &
                       frl,dlrad,sabg,fseng,fevpg,cgrnd,htvp,emg, &
                       imelt,sm,xmf,fact)

!=======================================================================
! Snow and soil temperatures
! o The volumetric heat capacity is calculated as a linear combination
!   in terms of the volumetric fraction of the constituent phases.
! o The thermal conductivity of soil is computed from
!   the algorithm of Johansen (as reported by Farouki 1981), and of snow is from
!   the formulation used in SNTHERM (Jordan 1991).
! o Boundary conditions:
!   F = Rnet - Hg - LEg (top),  F= 0 (base of the soil column).
! o Soil / snow temperature is predicted from heat conduction
!   in 10 soil layers and up to 5 snow layers.
!   The thermal conductivities at the interfaces between two neighbor layers
!   (j, j+1) are derived from an assumption that the flux across the interface
!   is equal to that from the node j to the interface and the flux from the
!   interface to the node j+1. The equation is solved using the Crank-Nicholson
!   method and resulted in a tridiagonal system equation.
!
! Phase change (see meltf.F90)
! 
! Original author : Yongjiu Dai, 09/15/1999; 08/30/2002
!=======================================================================

  use precision
  use phycon_module, only : stefnc
  implicit none

  integer, INTENT(in) :: lb             !lower bound of array
  integer, INTENT(in) :: nl_soil        !upper bound of array
  integer, INTENT(in) :: n_pft          !number of pfts in a single column
  integer, INTENT(in) :: itypwat        !land water type (0=soil,1=urban or built-up,2=wetland,
                                        !3=land ice, 4=deep lake, 5=shallow lake)
  integer, INTENT(in) :: num_filterp    !
  integer, INTENT(in) :: filterp(n_pft) !
  real(r8), INTENT(in) :: wt_patch(n_pft) ! relative weight of pfts to grid area
  real(r8), INTENT(in) :: wt_column     ! relative weight of column to grid area

  real(r8), INTENT(in) :: dtime    !model time step [second]
  real(r8), INTENT(in) :: capr     !tuning factor to turn first layer T into surface T
  real(r8), INTENT(in) :: cnfac    !Crank Nicholson factor between 0 and 1

  real(r8), INTENT(in) :: csol(1:nl_soil)  !heat capacity of soil solids [J/(m3 K)]
  real(r8), INTENT(in) :: porsl(1:nl_soil) !soil porosity [-]
  real(r8), INTENT(in) :: dkmg(1:nl_soil)  !thermal conductivity of soil minerals [W/m-K]
  real(r8), INTENT(in) :: dkdry(1:nl_soil) !thermal conductivity of dry soil [W/m-K]
  real(r8), INTENT(in) :: dksatu(1:nl_soil)!thermal conductivity of saturated soil [W/m-K]

  real(r8), INTENT(in) :: sigf(n_pft)      !fraction of veg cover, excluding snow-covered veg [-]
  real(r8), INTENT(in) :: dz(lb:nl_soil)   !layer thickiness [m]
  real(r8), INTENT(in) :: z (lb:nl_soil)   !node depth [m]
  real(r8), INTENT(in) :: zi(lb-1:nl_soil) !interface depth [m]

  real(r8), INTENT(in) :: sabg(n_pft)     !solar radiation absorbed by ground [W/m2]
  real(r8), INTENT(in) :: frl             !atmospheric infrared (longwave) radiation [W/m2]
  real(r8), INTENT(in) :: dlrad(n_pft)    !downward longwave radiation blow the canopy [W/m2]
  real(r8), INTENT(in) :: fseng(n_pft)    !sensible heat flux from ground [W/m2]
  real(r8), INTENT(in) :: fevpg(n_pft)    !evaporation heat flux from ground [mm/s]
  real(r8), INTENT(in) :: cgrnd(n_pft)    !deriv. of soil energy flux wrt to soil temp [w/m2/k]
  real(r8), INTENT(in) :: htvp            !latent heat of vapor of water (or sublimation) [j/kg]
  real(r8), INTENT(in) :: emg             !ground emissivity (0.97 for snow,

  real(r8), INTENT(inout) :: tss (lb:nl_soil) !soil temperature [K]
  real(r8), INTENT(inout) :: wice(lb:nl_soil) !ice lens [kg/m2]
  real(r8), INTENT(inout) :: wliq(lb:nl_soil) !liqui water [kg/m2]
  real(r8), INTENT(inout) :: scv              !snow cover, water equivalent [mm, kg/m2]
  real(r8), INTENT(inout) :: snowdp           !snow depth [m]

  real(r8), INTENT(out) :: sm               !rate of snowmelt [kg/(m2 s)]
  real(r8), INTENT(out) :: xmf              !total latent heat of phase change of ground water
  real(r8), INTENT(out) :: fact(lb:nl_soil) !used in computing tridiagonal matrix
  integer, INTENT(out) :: imelt(lb:nl_soil)    !flag for melting or freezing [-]

!------------------------ local variables ------------------------------
  real(r8) cv(lb:nl_soil)     ! heat capacity [J/(m2 K)]
  real(r8) tk(lb:nl_soil)     ! thermal conductivity [W/(m K)]

  real(r8) at(lb:nl_soil)     !"a" vector for tridiagonal matrix
  real(r8) bt(lb:nl_soil)     !"b" vector for tridiagonal matrix
  real(r8) ct(lb:nl_soil)     !"c" vector for tridiagonal matrix
  real(r8) rt(lb:nl_soil)     !"r" vector for tridiagonal solution

  real(r8) fn  (lb:nl_soil)   ! heat diffusion through the layer interface [W/m2]
  real(r8) fn1 (lb:nl_soil)   ! heat diffusion through the layer interface [W/m2]
  real(r8) dzm                ! used in computing tridiagonal matrix
  real(r8) dzp                ! used in computing tridiagonal matrix

  real(r8) tssbef(lb:nl_soil) ! soil/snow temperature before update
  real(r8) hs_pft(n_pft)      ! net energy flux into the surface on pfts(w/m2)
  real(r8) hs                 ! net energy flux into the surface on column(w/m2)
  real(r8) dhsdT_pft(n_pft)   ! d(hs)/dT on pft
  real(r8) dhsdT              ! d(hs)/dT on column
  real(r8) brr(lb:nl_soil)    ! temporay set

  integer i,j,p,fp

!=======================================================================

! heat capacity 
      call hCapacity (itypwat,lb,nl_soil,csol,porsl,wice,wliq,scv,dz,cv)

! thermal conductivity
      call hConductivity (itypwat,lb,nl_soil,&
                          dkmg,dkdry,dksatu,porsl,dz,z,zi,tss,wice,wliq,tk)


! net ground heat flux into the surface and its temperature derivative
! Added a pfts loop here to get the average of hs and dhsdT over all 
! PFTs on the column  add by zhq. 07/20/2009
      hs = 0.
      dhsdT = 0.

      do fp = 1, num_filterp
          p = filterp(fp)

        hs_pft(p) = sabg(p) + dlrad(p) &
           + (1.-sigf(p))*emg*frl - emg*stefnc*tss(lb)**4 &
           - (fseng(p)+fevpg(p)*htvp) 

        dhsdT_pft(p) = - cgrnd(p) - 4.*emg * stefnc * tss(lb)**3
!print*,p,'cgrnd',cgrnd(p)
!print*,p,'dhsdT',dhsdT_pft(p)
        hs = hs + hs_pft(p) * wt_patch(p) / wt_column
        dhsdT = dhsdT + dhsdT_pft(p) *wt_patch(p) / wt_column

      end do
      tssbef(lb:) = tss(lb:)
      j       = lb
      fact(j) = dtime / cv(j) &
              * dz(j) / (0.5*(z(j)-zi(j-1)+capr*(z(j+1)-zi(j-1))))

      do j = lb + 1, nl_soil
        fact(j) = dtime/cv(j)
      enddo

      do j = lb, nl_soil - 1
        fn(j) = tk(j)*(tss(j+1)-tss(j))/(z(j+1)-z(j))
      enddo
      fn(nl_soil) = 0.

! set up vector r and vectors a, b, c that define tridiagonal matrix
      j     = lb
      dzp   = z(j+1)-z(j)
      at(j) = 0.
      bt(j) = 1+(1.-cnfac)*fact(j)*tk(j)/dzp-fact(j)*dhsdT
      ct(j) =  -(1.-cnfac)*fact(j)*tk(j)/dzp
      rt(j) = tss(j) + fact(j)*( hs - dhsdT*tss(j) + cnfac*fn(j) )
!print*,'bt-1',bt
!print*,'rt-1',rt


      do j = lb + 1, nl_soil - 1
         dzm   = (z(j)-z(j-1))
         dzp   = (z(j+1)-z(j))
         at(j) =   - (1.-cnfac)*fact(j)* tk(j-1)/dzm
         bt(j) = 1.+ (1.-cnfac)*fact(j)*(tk(j)/dzp + tk(j-1)/dzm)
         ct(j) =   - (1.-cnfac)*fact(j)* tk(j)/dzp
         rt(j) = tss(j) + cnfac*fact(j)*( fn(j) - fn(j-1) )
      end do

      j     =  nl_soil
      dzm   = (z(j)-z(j-1))
      at(j) =   - (1.-cnfac)*fact(j)*tk(j-1)/dzm
      bt(j) = 1.+ (1.-cnfac)*fact(j)*tk(j-1)/dzm
      ct(j) = 0.
      rt(j) = tss(j) - cnfac*fact(j)*fn(j-1)

! solve for tss
      i = size(at)
      call tridia (i ,at ,bt ,ct ,rt ,tss) 
!print*,'tss-3',tss
!=======================================================================
! melting or freezing 
!=======================================================================

      do j = lb, nl_soil - 1
         fn1(j) = tk(j)*(tss(j+1)-tss(j))/(z(j+1)-z(j))
      enddo
      fn1(nl_soil) = 0.

      j = lb
      brr(j) = cnfac*fn(j) + (1.-cnfac)*fn1(j)

      do j = lb + 1, nl_soil
         brr(j) = cnfac*(fn(j)-fn(j-1)) + (1.-cnfac)*(fn1(j)-fn1(j-1))
      enddo

      call meltf ( lb, nl_soil, dtime, &
                   fact(lb:), brr(lb:), hs, dhsdT, &
                   tssbef(lb:), tss(lb:), wliq(lb:), wice(lb:), imelt(lb:), &
                   scv, snowdp, sm, xmf )


!-----------------------------------------------------------------------

 end subroutine groundtem
