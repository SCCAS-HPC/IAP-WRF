
 subroutine GETMET (lumet,lon_points,lat_points, &
                    tair2d,qair2d,pres2d,rainc2d,rainl2d, &
                    windu2d,windv2d,dswrf2d,dlwrf2d,  &
                    tair_z,qair_z,wind_z)
! ======================================================================
!
! Access meteorological data
!
! Original author : Yongjiu Dai
! ======================================================================

   use precision
   implicit none

   integer, INTENT(in)  :: lumet       ! logical unit number of atmospheric forcing data 
   integer, INTENT(in)  :: lon_points  ! number of longitude points on model grid
   integer, INTENT(in)  :: lat_points  ! number of latitude points on model grid
   real(r8),intent(out) :: tair2d  (lon_points,lat_points)
   real(r8),intent(out) :: qair2d  (lon_points,lat_points)
   real(r8),intent(out) :: pres2d  (lon_points,lat_points)
   real(r8),intent(out) :: rainc2d (lon_points,lat_points)
   real(r8),intent(out) :: rainl2d (lon_points,lat_points)
   real(r8),intent(out) :: windu2d (lon_points,lat_points)
   real(r8),intent(out) :: windv2d (lon_points,lat_points)
   real(r8),intent(out) :: dswrf2d (lon_points,lat_points)
   real(r8),intent(out) :: dlwrf2d (lon_points,lat_points)
   real(r8),intent(out) :: tair_z  (lon_points,lat_points)
   real(r8),intent(out) :: qair_z  (lon_points,lat_points)
   real(r8),intent(out) :: wind_z  (lon_points,lat_points)

! local
   real(r8) solar                    ! incident solar radiation [W/m2]
   real(r8) frl                      ! atmospheric infrared (longwave) radiation [W/m2]
   real(r8) prcp                     ! precipitation [mm/s]
   real(r8) tm                       ! temperature at reference height [kelvin]
   real(r8) us                       ! wind component in eastward direction [m/s]
   real(r8) vs                       ! wind component in northward direction [m/s]
   real(r8) pres                     ! atmospheric pressure [pa]
   real(r8) qm                       ! specific humidity at reference height [kg/kg]
       
   integer i, j                      ! looping index

   integer IYEAR,IMONTH,IDAY,IHOUR
   real(r8) rnet,twet,esat,esatdT,qsat,qsatdT,us_d 
   real(r8) rrref, ddd, ggg
   character(len=2) site
   real(r8) time, VPA
   integer flag, Rflag

! ------------------ note ------------------
! the model required the same longitudinal resolution for
! all latitude strip. For the variable longitudinal resolution
! cases, please assign the meteorological values to 
! -999 to the land grids which are not included in the calcultion.
! ----------------------------------------------------------------------

![1] PILPS's Valdai (obs height: wind 10 m, tem & hum 2 m)
! & CABAUW (obs height = 20 m for all) & HAPEX (Obs heigh = 2 m for all)
 
      read (lumet,10) solar, frl, prcp, tm, us, vs, pres, qm
10    format (2f7.1, e14.3, 3f10.3, f10.1, e12.3)

!-----------------------------------------
!  ABRACOS site Reserva Jaru & Fazenda N.S. 
! (obs height forest: zwind=52.5)
! (obs height pasture: zwind=5.4, ztem=5.1)
 
!     read(lumet,30) site,iyear,iday,ihour,&
!                      solar,rrref,rnet,twet,tm,us,ddd,ggg,prcp
30    format(1x,a2,1x,i4,1x,i3,1x,i4,1x,9(f9.2))
! 
! for missing data
!     if (prcp.lt.0.0) prcp=0.0
!     if (solar.lt.0.0) solar=0.0
!     if (us.lt.0.0) us = 1.0
!     vs=0.0
!
!     if(rrref.eq.-99.99) rrref=solar*0.122
!     if(rrref.lt.0.0) rrref=0.0
!
!     if(rnet.eq.-99.99) then
!       if(site.eq.'RD' .or. site.eq.'rd') rnet = 0.79*solar - 26.08
!       if(site.eq.'FD' .or. site.eq.'fd') rnet = 0.79*solar - 16.26
!     endif
!
!     pres =1000.0*100.0 
! specific humidity
!     tm = tm + 273.16 - 50.
!     if(tm.lt.274.) stop 'meterological data error 2'
!     twet  = twet + 273.16
!     call qsadv(twet,pres,esat,esatdT,qsat,qsatdT)
!     qm = qsat - 3.5*2.8704E2*(tm-twet)/2.50036e6
!
!     !for Arme, Abracos, Tucson, FIFE 88 & 89:
!     !if(istep.gt.11760)then  ! only used for FIFE Forcing data
!     frl = rnet-(solar-rrref)+5.67e-8*tm**4
!     !endif
!
!     prcp=prcp/3600
! ----------------------------------------------------------------------
      do j = 1, lat_points
         do i = 1, lon_points
            tair2d(i,j) =  tm
            qair2d(i,j) =  qm
            pres2d(i,j) =  pres
           rainc2d(i,j) =  0.
           rainl2d(i,j) =  prcp
           windu2d(i,j) =  us
           windv2d(i,j) =  vs
           dswrf2d(i,j) =  solar
           dlwrf2d(i,j) =  frl
            tair_z(i,j) =  50.
            qair_z(i,j) =  50.
            wind_z(i,j) =  50.
         enddo 
      enddo 

 end subroutine GETMET

