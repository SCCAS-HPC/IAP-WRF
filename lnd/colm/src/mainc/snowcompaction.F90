
 subroutine snowcompaction ( lb, dtime, &
                      imelt, fiold, tss, wliq, wice, dz )
 
!=======================================================================
! Original author: Yongjiu Dai, September 15, 1999
!
! three of metamorphisms of changing snow characteristics are implemented,
! i.e., destructive, overburden, and melt. The treatments of the former two
! are from SNTHERM.89 and SNTHERM.99 (1991, 1999). The contribution due to
! melt metamorphism is simply taken as a ratio of snow ice fraction after
! the melting versus before the melting.
!
!=======================================================================

  use precision
  use phycon_module, only : denice, denh2o, tfrz
  implicit none
!
!-------------------------- Dummy argument -----------------------------
!
  integer, INTENT(in) :: &
        lb             ! lower bound of array

  integer, INTENT(in) :: & 
        imelt(lb : 0)  ! signifies if node in melting (imelt = 1)

  real(r8), INTENT(in) :: &
        dtime,        &! time step [second]
        fiold(lb : 0),&! fraction of ice relative to 
                       ! the total water content at the previous time step 
        tss (lb : 0), &! nodal temperature [K]
        wice(lb : 0), &! ice lens [kg/m2]
        wliq(lb : 0)   ! liquid water [kg/m2]

  real(r8), INTENT(inout) :: &
        dz  (lb : 0)   ! layer thickness [m]
!
!----------------------- local variables ------------------------------
  integer i            ! Numeber of doing loop
  
  real(r8) :: &
       c1,            &! = 2.777e-7 [m2/(kg s)]
       c2,            &! = 21e-3 [m3/kg]
       c3,            &! = 2.777e-6 [1/s]
       c4,            &! = 0.04 [1/K]
       c5,            &! = 2.0 
       c6,            &! = 5.15e-7.
       c7,            &! = 4.
       dm,            &! Upper Limit on Destructive Metamorphism Compaction [kg/m3]
       eta0            ! The Viscosity Coefficient Eta0 [kg-s/m2]

  real(r8) :: &
       burden,        &! pressure of overlying snow [kg/m2]
       ddz1,          &! Rate of settling of snowpack due to destructive metamorphism.
       ddz2,          &! Rate of compaction of snowpack due to overburden.
       ddz3,          &! Rate of compaction of snowpack due to melt [1/s]
       dexpf,         &! expf=exp(-c4*(273.15-tss)).
       fi,            &! Fraction of ice relative to the total water content 
                       ! at the current time step
       td,            &! tss - tfrz [K]
       pdzdtc,        &! Nodal rate of change in fractional-thickness 
                       ! due to compaction [fraction/s]
       void,          &! void (1 - vol_ice - vol_liq)
       wx,            &! water mass (ice+liquid) [kg/m2]
       bi              ! partitial density of ice [kg/m3]

  data c2,c3,c4,c5/23.e-3, 2.777e-6, 0.04, 2.0/
  data c6/5.15e-7/, c7/4./
  data dm/100./         
  data eta0/9.e5/      

!=======================================================================

      burden = 0.0

      do i = lb, 0
         wx = wice(i) + wliq(i)
         void = 1.- (wice(i)/denice + wliq(i)/denh2o)/dz(i)

! Disallow compaction for water saturated node and lower ice lens node.
         if(void <= 0.001 .or. wice(i) <= .1)then
            burden = burden+wx
            CYCLE
         endif

         bi = wice(i) / dz(i)
         fi = wice(i) / wx
         td = tfrz-tss(i)
         
         dexpf = exp(-c4*td)

! Settling as a result of destructive metamorphism
         ddz1 = -c3*dexpf
         if(bi > dm) ddz1 = ddz1*exp(-46.0e-3*(bi-dm))

! Liquid water term
         if(wliq(i) > 0.01*dz(i)) ddz1=ddz1*c5

! Compaction due to overburden
         ddz2 = -burden*exp(-0.08*td - c2*bi)/eta0

! Compaction occurring during melt
         if(imelt(i) == 1)then
            ddz3 = - 1./dtime * max(0.,(fiold(i) - fi)/fiold(i))
         else
            ddz3 = 0.
         endif

! Time rate of fractional change in dz (units of s-1)
         pdzdtc = ddz1+ddz2+ddz3

! The change in dz due to compaction
         dz(i) = dz(i)*(1.+pdzdtc*dtime)

! Pressure of overlying snow
         burden = burden+wx

      end do

 end subroutine snowcompaction
