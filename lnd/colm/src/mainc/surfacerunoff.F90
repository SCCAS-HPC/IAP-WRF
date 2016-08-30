
 subroutine surfacerunoff (nl_soil,wtfact,wimp,bsw,porsl,phi0,hksati,&
                           z,dz,zi,vol_liq,vol_ice,eff_porosity,gwat,rsur) 
	           
!=======================================================================
! the original code was provide by Robert E. Dickinson based on following clues:
! a water table level determination level added including highland and 
! lowland levels and fractional area of wetland (water table above the surface. 
! Runoff is parametrized from the lowlands in terms of precip incident on 
! wet areas and a base flow, where these are estimated using ideas from TOPMODEL.
!
! Author : Yongjiu Dai, 07/29/2002
!=======================================================================

  use precision
  implicit none
    
!-----------------------Arguments---------------------------------------

  integer, INTENT(in) :: nl_soil   ! number of soil layers
  real(r8), INTENT(in) :: &
        wtfact,        &! fraction of model area with high water table
        wimp,          &! water impremeable if porosity less than wimp
       bsw(1:nl_soil), &! Clapp-Hornberger "B"
     porsl(1:nl_soil), &! saturated volumetric soil water content(porosity)
      phi0(1:nl_soil), &! saturated soil suction (mm)
    hksati(1:nl_soil), &! hydraulic conductivity at saturation (mm h2o/s)
         z(1:nl_soil), &! layer depth (m)
        dz(1:nl_soil), &! layer thickness (m)
        zi(0:nl_soil), &! interface level below a "z" level (m)
   vol_liq(1:nl_soil), &! partial volume of liquid water in layer [-]
   vol_ice(1:nl_soil), &! partial volume of ice lens in layer [-]
 eff_porosity(1:nl_soil), &! effective porosity = porosity - vol_ice
        gwat            ! net water input from top

  real(r8), INTENT(out) :: rsur    ! surface runoff (mm h2o/s)
                    
!-----------------------Local Variables---------------------------------

  real(r8) fcov         ! fractional area with water table at surface
  real(r8) fz           ! coefficient for water table depth
  real(r8) zwt          ! water table depth
  real(r8) deficit1     ! water deficit from 10-layer soil moisture
  real(r8) deficit2     ! water deficit from 100-layer
  real(r8) dzfine       ! layer thickness for fine mesh soil layers
  real(r8) zfine(1:100) ! depth for fine mesh soil layers
  real(r8) s1,su,v      ! local variables to calculate qinmax
  real(r8) qinmax       ! maximum infiltration capability
  real(r8) watsatm      ! averaged soil porosity (-)
  real(r8) sucsatm      ! averaged minimum soil suction (mm)
  real(r8) bswm         ! averaged Clapp and Hornberger "b"
  real(r8) a
  real(r8) belta

  integer i             ! loop counter
                      
!-----------------------End Variable List-------------------------------

! average soil properties

      watsatm = 0.
      sucsatm = 0.
      bswm    = 0.

      do i=1,nl_soil
         watsatm = watsatm + porsl(i)
         sucsatm = sucsatm + phi0(i)
         bswm    = bswm    + bsw(i)
      end do

      watsatm = watsatm / nl_soil
      sucsatm = sucsatm / nl_soil
      bswm    = bswm    / nl_soil

! Determine water table
      deficit1 = 0.
      do i = 1,nl_soil
      deficit1 = deficit1+ (porsl(i)-(vol_ice(i)+vol_liq(i)))*dz(i) 
      enddo

!     dzfine = zi(nl_soil)/100. 
! Correction 04/06/2007
      dzfine = 2. * zi(nl_soil)/100. 
      do i =1,100
        zfine(i) = i*dzfine
      enddo

      deficit2 = 0.
!     zwt = zi(nl_soil) - 0.001   
! Correction 04/06/2007
      zwt = 2. * zi(nl_soil) - 0.001   

      do i = 1, 100
!        a = (zwt-zfine(i))*1000./sucsatm
!        a = max(1.e-3, a)
!        belta = a**(-1./bswm)
!        deficit2 = deficit2 + belta*watsatm*dzfine
! Correction 04/06/2007
         deficit2 = deficit2 + watsatm * &
                   (1.-(1.+(zwt-zfine(i))*1000./sucsatm)**(-1./bswm))*dzfine

         if(abs(deficit2-deficit1).le.0.1) then
           zwt = zfine(i)
           exit
         endif
      enddo

! Saturation fraction
      fcov = wtfact
      if(eff_porosity(1)>=wimp) fcov = wtfact*min(1.,exp(-zwt))

! Maximum infiltration capacity
      s1 = max(0.01,vol_liq(1)/max(wimp, eff_porosity(1)))  
      su = max(0.,(s1-fcov)/(1.-fcov))
      v = -bsw(1)*phi0(1)/(0.5*dz(1)*1000.)
      qinmax = (1.+v*(su-1.))*hksati(1)
      if(eff_porosity(1)<wimp) qinmax = 0.

! Surface runoff
      rsur = fcov*gwat + (1.-fcov)*max(0.,gwat-qinmax)

 end subroutine surfacerunoff
