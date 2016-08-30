module drydep
contains
       subroutine getland(landuse0,ilanduse0)
!cccccccc  the routine is convert the landuse in MM5  to RADMII
          real :: landuse0 ! the 25 catagories in MM5
          integer :: ilanduse0 ! the 11 catagories in RADM
!--------------------------------------Table----------------------
!C        Change MM5 24 landuse categories to RADM 11 categories
!C        Also, change percentages to fractions
!C        -----------------------------------------------
!C        MM5        CATEGORIES                 CMAQ/RADM
!C        -----------------------------------------------
!C         1      URBAN LAND                          1
!C         2      DRYLAND CROPLAND & PASTURE          2
!C         3      IRRIGATED CROPLAND & PASTURE        2
!C         4      MIXED DRYLAND & IRRIGATED CROPLAND  2
!C         5      CROPLAND/GRASSLAND MOSAIC           2
!C         6      CROPLAND/WOODLAND MOSAIC           10
!C         7      GRASSLAND                           3
!C         8      SHRUBLAND                           3
!C         9      MIXED SHRUBLAND/GRASSLAND           3
!C        10      SAVANNAH                           10
!C        11      DECIDUOUS BROADLEAF FOREST          4
!C        12      DECIDUOUS NEEDLELEAF FOREST         4
!C        13      EVERGREEN BROADLEAF FOREST          5
!C        14      EVERGREEN NEEDLELEAF FOREST         5
!C        15      MIXED FOREST                        6
!C        16      WATER                               7
!C        17      HERBACEOUS WETLAND                  9
!C        18      WOODED WETLAND                      5
!C        19      BARREN OR SPARSELY VEGETATED        8
!C        20      HERBACEOUS TUNDRA                  11
!C        21      WOODED TUNDRA                       5
!C        22      MIXED TUNDRA                       10
!C        23      BARE GROUND TUNDRA                  8
!C        24      SNOW OR ICE                        11
!-----------------------------------------------------------
         select case(int(landuse0))
             case(0)
                 ilanduse0=0 ! no data
             case(1)
                 ilanduse0=1
             case(2)
                 ilanduse0=2
             case(3)
                 ilanduse0=2
             case(4)
                 ilanduse0=2
             case(5)
                 ilanduse0=2
             case(6)
                 ilanduse0=10
             case(7)
                 ilanduse0=3
             case(8)
                 ilanduse0=3
             case(9)
                 ilanduse0=3
             case(10)
                 ilanduse0=10
             case(11)
                 ilanduse0=4
             case(12)
                 ilanduse0=4
             case(13)
                 ilanduse0=5
             case(14)
                 ilanduse0=5
             case(15)
                 ilanduse0=6
             case(16)
                 ilanduse0=7
             case(17)
                 ilanduse0=9
             case(18)
                 ilanduse0=5
             case(19)
                 ilanduse0=8
             case(20)
                 ilanduse0=11
             case(21)
                 ilanduse0=5
             case(22)
                 ilanduse0=10
             case(23)
                 ilanduse0=8
             case(24)
                 ilanduse0=11
          end select
        return

        end subroutine

    subroutine drydep_aer(myid,month,tsurf,cellat,cellon,height,&
                            press,windu,windv,landuse,tempk,&
                    pbl_hgt,diam,vdep_aer,sx,ex,sy,ey,nzz,ig)
      integer :: sx,ex,sy,ey
      real tsurf(sx-1:ex+1,sy-1:ey+1),cellat(sx-1:ex+1,sy-1:ey+1),&
          cellon(sx-1:ex+1,sy-1:ey+1)
      real height(sx-1:ex+1,sy-1:ey+1,nzz),&     
          press(sx-1:ex+1,sy-1:ey+1,nzz),&
          windu(sx-1:ex+1,sy-1:ey+1,nzz),&
          windv(sx-1:ex+1,sy-1:ey+1,nzz),&
          tempk(sx-1:ex+1,sy-1:ey+1,nzz)
      real vdep_aer(sx-1:ex+1,sy-1:ey+1),pbl_hgt(sx-1:ex+1,sy-1:ey+1)
      logical lstable
      REAL :: iseason(5,12)
      REAL :: diam  
! diam : log-mean sectional aerosol diameter (m)
      REAL :: rhop ! AEROSOL DENISITY in g/m3
      INTEGER ::  landuse(sx-1:ex+1,sy-1:ey+1)
      REAL,DIMENSION(11,5) :: z0lu
!
      data pi/3.1415927/
      data cwmin /3.8E-5/
!
!c-----Season indices by month and latitude band
!c     Season Indices            Latitude Bands
!c     1 = summer                1 = <20    Tropical
!c     2 = autumn                2 = 20-35  Sub-tropical
!c     3 = winter w/o snow       3 = 35-50  Temperate
!c     4 = winter w/ snow        4 = 50-75  Cool
!c     5 = spring                5 = >75    Polar
!c                    Latitude Band
      data iseason / 1, 3, 3, 3, 3,& ! Jan
                    1, 5, 3, 3, 3,& ! Feb
                    1, 5, 5, 3, 3,& ! Mar
                    1, 5, 5, 5, 3,& ! Apr
                    1, 1, 5, 5, 3,& ! May
                    1, 1, 1, 1, 5,& ! Jun
                    1, 1, 1, 1, 1,& ! Jul
                    1, 1, 1, 1, 2,& ! Aug
                    1, 1, 2, 2, 3,& ! Sep
                    1, 2, 2, 2, 3,& ! Oct
                    1, 2, 2, 3, 3,& ! Nov
                    1, 2, 3, 3, 3/ ! Dec
!c
!c-----Surface roughness (m) as a function of 11 landuse categories
!c     and 5 seasons; based on AERMET model (ref EPA SCRAM website)
!c
      data z0lu&
      /1.0,0.20,0.100,1.3,1.3,1.30,0.0001,0.002,0.20,0.150,0.30,&
       1.0,0.05,0.010,0.8,1.3,1.05,0.0001,0.002,0.20,0.030,0.30,&
       1.0,0.05,0.010,0.8,1.3,1.05,0.0001,0.002,0.20,0.030,0.30,&
       1.0,0.01,0.001,0.5,1.3,0.90,0.0001,0.002,0.05,0.006,0.15,&
       1.0,0.03,0.050,1.0,1.3,1.15,0.0001,0.002,0.20,0.040,0.30/

      IF(ig == 75 ) rhop = 1.0E06 ! Primary PM25
      IF(ig == 76 ) rhop = 1.0E06 ! Primary PM10
      IF(ig == 77 ) rhop = 2.0E06 ! BC
      IF(ig == 78 ) rhop = 1.0E06 ! Primary OC
      IF(ig == 79 ) rhop = 0.0    ! H+
      IF(ig == 80 ) rhop = 2.0E06 ! Na+
      IF(ig == 81 ) rhop = 1.5E06 ! NH4+
      IF(ig == 82 ) rhop = 2.0E06 ! CL-
      IF(ig == 83 ) rhop = 1.5E06 ! SO42-
      IF(ig == 84 ) rhop = 1.5E06 ! HSO4-
      IF(ig == 85 ) rhop = 1.5E06 ! NO3-
      IF(ig == 86 ) rhop = 2.0E06 ! NACL
      IF(ig == 87 ) rhop = 1.5E06 ! NA2SO4
      IF(ig == 88 ) rhop = 1.5E06 ! NANO3
      IF(ig == 89 ) rhop = 1.5E06 ! NH42SO4
      IF(ig == 90 ) rhop = 1.5E06 ! NH4NO3
      IF(ig == 91 ) rhop = 1.5E06 ! NH4CL
      IF(ig == 92 ) rhop = 1.5E06 ! H2SO4
      IF(ig == 93 ) rhop = 1.5E06 ! NH4HSO4
      IF(ig == 94)  rhop = 1.5E06 ! NAHSO4
      IF(ig == 95)  rhop = 1.5E06 ! (NH4)4H(SO4)2(S)
      IF(ig == 96)  rhop = 1.0E06 ! SOA1
      IF(ig == 97)  rhop = 1.0E06 ! SOA2
      IF(ig == 98)  rhop = 1.0E06 ! SOA3
      IF(ig == 99)  rhop = 1.0E06 ! SOA4
      IF(ig == 100) rhop = 1.0E06 ! SOA5
      IF(ig == 101) rhop = 1.0E06 ! SOA6
      IF(ig == 102) rhop = 1.0E06 ! AH2O       
      IF(ig == 103) rhop = 2.0E06 ! HGP, by chenhs       

      do 30 j = sy,ey
        do 20 i = sx,ex
!c
!c-----Determine season
!c
          mbin = month
          if (cellat(i,j).lt.0.) then
            mbin = mod(month+6,12)
            if (mbin.eq.0) mbin = 12
          endif
          latbin = 1
          if (abs(cellat(i,j)).gt.20.) then
          latbin = 2
          elseif (abs(cellat(i,j)).gt.35.) then
          latbin = 3
          elseif (abs(cellat(i,j)).gt.50.) then
          latbin = 4
          elseif (abs(cellat(i,j)).gt.75.) then
          latbin = 5
          endif
          if ((cellat(i,j).gt.50. .and. cellat(i,j).lt.75.) .and.&
             (cellon(i,j).gt.-15. .and. cellon(i,j).lt.15.)) latbin = 3
           isesn = iseason(latbin,mbin)
!c
!c-----Use input snow cover to set season, if specified
!c
!c
!c-----Load local met variables
          if( height(i,j,1).gt.0 ) then  !! by chenhs
            deltaz = height(i,j,1)/2.
          else
            deltaz=0.
          endif
          temp0 = tsurf(i,j) - 273.15
          prss0 = press(i,j,1) -&
                 2.*deltaz*(press(i,j,2) - press(i,j,1))/height(i,j,2)
          ucomp = (windu(i,j,1) + windu(i-1,j,1))/2.
          vcomp = (windv(i,j,1) + windv(i,j-1,1))/2.
          wind = sqrt(ucomp**2 + vcomp**2)
          wind = amax1(0.1,wind)
!c
!c
!c-----Loop over land use; surface roughness for water is dependent on        
!c     wind speed
         vdep_aer(i,j) = 0.
         totland = 1.
                   
        m=landuse(i,j)
        z0 = z0lu(m,isesn)
        if (m.eq.7) z0 = amax1(z0,2.0e-6*wind**2.5)
!c
!c-----Use input surface roughness, if specified
!c
!c
!c-----Get surface layer micrometeorological parameters for this cell
!and
!c     landuse type
!c
         if (prss0.lt.0) then
           write(iout,'(//,a)') 'ERROR in DRYDEP:'
           write(iout,*) 'Invalid pressure value'
           write(iout,*) 'Cell   Height  Deltaz'
           write(iout,*) i,j,height(i,j,1),deltaz,press(i,j,1),press(i,j,2),prss0
         endif

         pbl = pbl_hgt(i,j)
         call micromet(tempk(i,j,1),tsurf(i,j),press(i,j,1),prss0,&
                    deltaz,wind,z0,pbl,ustar,psih,wstar,lstable)
!cc
!cc         call vd_aer(z0,deltaz,psih,ustar,diam,rhop,temp0,vd)      
         call vd_aer(z0,deltaz,psih,ustar,diam,rhop,tsurf(i,j),vd)      

         vdep_aer(i,j) =  vd
 
   20     continue
   30   continue
    
                            
      return
      end  subroutine    

      subroutine vd_aer(z0,deltaz,psih,ustar,diam,rhop,ts,vd)
!c
!c----CAMx v4.42 070603
!c
!c     VD_AER calculates a deposition velocity for a specific aerosol size
!c     bin, grid cell, and land use category.  The parallel resistance approach
!c     of Slinn and Slinn (1980) is used, as implemented in UAM-AERO
!c     (STI, 1996).
!c
!c     Copyright 1996-2007
!c     ENVIRON International Corporation
!c
!c     Modifications:
!c        3/21/03       Modified Ra calculation to use layer 1 midpoint height
!c                      rather than default 10 m.
!c        9/18/03       Removed small bug in final Vd equation
!c
!c     Input arguments:
!c        z0                  surface roughness length (m)
!c        deltaz              Layer 1 midpoint height (m)
!c        psih                similarity stability correction term
!c        ustar               friction velocity (m/s)
!c        diam                log-mean sectional aerosol diameter (m)
!c        rhop                aerosol density (g/m3)
!c        ts                  surface temperature (K)
!c
!c     Output arguments:
!c        vd                  deposition velocity (m/s)
!c
!c     Routines called:
!c        none
!c
!c     Called by:
!c        DRYDEP

      data vk/0.4/, rmin/1.0/, xmfp/6.5e-8/, g/9.8/, vabs/1.81e-2/
      data boltz/1.38e-20/, pi/3.1415927/, vair/1.5e-5/
!c
!c-----Entry point
!c
!c-----Speed correction factor and sedimendation velocity
!c
      power = amin1(7.6,0.55*diam/xmfp)
      scf = 1. + (2.514 + 0.8*exp(-power))*xmfp/diam
      vsed = rhop*g*(diam*diam)*scf/(18.*vabs)
!c
!c-----Brownian diffusivity and Schmidt number
!c
      difbrwn = boltz*ts*scf/(3.*pi*vabs*diam)
      schmidt = vair/difbrwn
!c
!c-----Stokes number
!c
      stokes = vsed*(ustar*ustar)/(vair*g)
!c
!c-----Compute atmospheric resistance, RA
!c
      ra = (alog(deltaz/z0) - psih)/(vk*ustar)
      ra = amax1(ra,rmin)
!c
!c-----Compute the deposition layer resistance, RD
!c
      sc23 = schmidt**(-2./3.)
      power = -3./stokes
      if (power.lt.-37.) then
        xinert = 10.**(-37.)
      else
        xinert = 10.**(power)
      endif
      rd = 1./(ustar*(sc23 + xinert))
      rd = amax1(rd,rmin)
!c-----Final deposition velocity for this cell, land use, and aerosol size
!c
      vd = vsed + 1./(ra + rd + ra*rd*vsed)
!c
      return
      end subroutine

       subroutine drydep_gas(myid,month,tsurf,cellat,cellon,pwc,cwc,height,&
                       press,windu,windv,solflux,landuse,water,&
                       tempk,pbl_hgt,lrddrt,icddrt,lrdsno,icdsno,vdep,&
                       sx,ex,sy,ey,nzz,ig)
!c
!c-----CAMx v4.42 070603
!c
!c
!c     DRYDEP is the driver for the calculation of gridded dry deposition
!c     velocities for a given grid. Deposition velocities are calculated for
!c     each gas species and for each aerosol size bin, weighted by the
!c     fractional land use specified for each cell.
!c
!c     Copyright 1996-2007
!c     ENVIRON International Corporation
!c
!c     Modifications:
!c        4/4/00    Added aerosol deposition as f(size)
!c        4/4/01    Fixed a few bugs in the call for VD_AER
!c        1/9/02    Aerosol size and density now species-dependent
!c        3/26/03   Added scaling factor to surface resistance (provided
!c                  on chemparam file), and zero dry deposition for
!c                  insoluble gases
!c        4/9/03    Removed rain, added precip and cloud water contents:
!c                  surfaces can now be rain-wetted, fog-wetted, or dew-wetted
!c        6/6/03    Protect against divide by zero with totland
!c        6/11/03   Use optional surface roughness length, if available
!c        7/21/03   Use optional drought stress and snow cover, if available;
!c                  Introduced latitude-dependent specification of season;
!c                  Revised solar flux calculation to use RADM cloud adjustment
!c        8/18/03   Relate drought stress to Palmer Drought Index
!c        4/21/04   Incorporated sectional PM
!c        11/19/04  Incorporated season-dependent roughness length
!c
!c     Input arguments:
!c        ncol                number of columns
!c        nrow                number of rows
!c        nlay                number of layers
!c        tsurf               surface temperature field (K)
!c        cellat              cell centroid latitude (deg)
!c        cellon              cell centroid longitude (deg)
!c        pwc                 precipitation water content (g/m3)
!c        cwc                 cloud water content (g/m3)
!c        height              layer interface height field (m)
!c        press               layer pressure field (mb)
!c        windu               layer U-component wind field (m/s)
!c        windv               layer V-component wind field (m/s)
!c        solflux             solar flux(w/m2)
!c        landuse             landuse
!c        water               layer water vapor field (ppm)
!c        tempk               layer temperature field (K)
!c        lrddrt              flag that gridded drought stress is available
!c        icddrt              optional drought index
!c        lrdsno              flag that gridded snow cover is available
!c        icdsno              optional snow index
!c
!c     Output arguments:
!c        vdep                species-dependent deposition velocity field (m/s)
!c
!c     Routines called:
!c        CALDATE
!c        GETZNTH
!c        MICROMET
!c        VD_GAS
!c        VD_AER
!c        HENRYFNC
!c
!c     Called by:
!c        CAMx
!c
!c
      integer :: sx,ex,sy,ey                 
      real tsurf(sx-1:ex+1,sy-1:ey+1),cellat(sx-1:ex+1,sy-1:ey+1),&
         cellon(sx-1:ex+1,sy-1:ey+1)
      integer icddrt(sx-1:ex+1,sy-1:ey+1),&
         icdsno(sx-1:ex+1,sy-1:ey+1)
      real height(sx-1:ex+1,sy-1:ey+1,nzz),&
          press(sx-1:ex+1,sy-1:ey+1,nzz),&
          windu(sx-1:ex+1,sy-1:ey+1,nzz),&
          windv(sx-1:ex+1,sy-1:ey+1,nzz),&
          water(sx-1:ex+1,sy-1:ey+1,nzz),&
          tempk(sx-1:ex+1,sy-1:ey+1,nzz),pwc(sx-1:ex+1,sy-1:ey+1,nzz),&
          cwc(sx-1:ex+1,sy-1:ey+1,nzz)
      real vdep(sx-1:ex+1,sy-1:ey+1),pbl_hgt(sx-1:ex+1,sy-1:ey+1)
      logical lstable
      logical ldark,lrdruf,lrddrt,lrdsno
      REAL :: iseason(5,12),solflux(sx-1:ex+1,sy-1:ey+1)
      INTEGER ::  landuse(sx-1:ex+1,sy-1:ey+1)
      REAL :: henry0(76),tfact(76),f0(76),rscale(76) !! 75=HG0, 76=HG2
      REAL :: diffrat(76),henry,dstress(0:5)         !! 75=HG0, 76=HG2
      REAL,DIMENSION(11,5) :: rj,rlu,rac,rgss,rgso,rlcs,rlco,z0lu
!c
      data eps/0.622/, e0/6.11/, lv/2.5e6/, rv/461./, pi/3.1415927/   ! juanxiong he
      data cwmin /3.8E-5/
!c
!c-----Season indices by month and latitude band
!c     Season Indices            Latitude Bands
!c     1 = summer                1 = <20    Tropical
!c     2 = autumn                2 = 20-35  Sub-tropical
!c     3 = winter w/o snow       3 = 35-50  Temperate
!c     4 = winter w/ snow        4 = 50-75  Cool
!c     5 = spring                5 = >75    Polar
!c                    Latitude Band
      data iseason / 1, 3, 3, 3, 3,& ! Jan
                    1, 5, 3, 3, 3,& ! Feb
                    1, 5, 5, 3, 3,& ! Mar
                    1, 5, 5, 5, 3,& ! Apr
                    1, 1, 5, 5, 3,& ! May
                    1, 1, 1, 1, 5,& ! Jun
                    1, 1, 1, 1, 1,& ! Jul
                    1, 1, 1, 1, 2,& ! Aug
                    1, 1, 2, 2, 3,& ! Sep
                    1, 2, 2, 2, 3,& ! Oct
                    1, 2, 2, 3, 3,& ! Nov
                    1, 2, 3, 3, 3/  ! Dec
!c
!c-----Surface roughness (m) as a function of 11 landuse categories
!c     and 5 seasons; based on AERMET model (ref EPA SCRAM website)
!c
      data z0lu&
      /1.0,0.20,0.100,1.3,1.3,1.30,0.0001,0.002,0.20,0.150,0.30,&
       1.0,0.05,0.010,0.8,1.3,1.05,0.0001,0.002,0.20,0.030,0.30,&
       1.0,0.05,0.010,0.8,1.3,1.05,0.0001,0.002,0.20,0.030,0.30,&
       1.0,0.01,0.001,0.5,1.3,0.90,0.0001,0.002,0.05,0.006,0.15,&
       1.0,0.03,0.050,1.0,1.3,1.15,0.0001,0.002,0.20,0.040,0.30/

      data henry0&
      /1.00e+10,  2.10e+05, 1.00e+05, 5.76e+01, 1.90e-03, 1.00e-02,&
       0.0,       3.20e+04, 5.90e+01, 0.0     , 1.10e-02, 0.0,&
       0.0,       0.0     , 0.0,      7.40e+04, 1.00e-10, 1.22e+00,&
       0.0,       1.00e-03, 0.0,      0.0,      6.30e+03, 0.0,&
       0.0,       0.0,      2.20e+02, 6.30e+03, 0.0,      0.0,&
       0.0,       3.60e+00, 1.00e-03, 0.0,      2.70e+03, 1.00e-02,&
       5.00e-03,  5.00e-03, 1.20e+00, 1.40e+00, 2.70e+03, 0.0,&
       0.0,       2.70e+03, 9.40e+03, 0.0,      0.0,      0.0,&
       0.0,       0.0,      1.00e-03, 1.00e-02, 6.30e+03, 6.30e+03,&
       6.30e+03,  6.30e+03, 0.0,      0.0,      0.0,      0.0,&
       0.0,       0.0,      0.0,      0.0,      0.0,      0.0,&
       1.00e+10, 2.7e+03,  2.7e+03,  2.7e+03, 2.7e+03,   2.7e+03,&
       2.7e+03,  2.7e+03,   0.11, 6.00e+05/
     data tfact&
      /0.0,    -8707., 0.0,    -4100.,  -1480.,  -2516.,&
       0.0,    -8706., -4781., 0.0,     -2415.,  0.0,&
       0.0,    0.0,    0.0,    -6643.,  0.0,     -3156.,&
       0.0,    0.0,    0.0,     0.0,    -6492.,  0.0,&
       0.0,    0.0,    -4932., -6492.,  0.0,     0.0,&
       0.0,    -5910., 0.0,    0.0,     -6492.,  0.0,&
       0.0,    0.0,    0.0,    0.0,     -6492.,  0.0,&
       0.0,    -6492., -8706., 0.0,     0.0,     0.0,&
       0.0,    0.0,    0.0,    0.0,     -6492.,  -6492.,&
       -6492., -6492., 0.0,    0.0,      0.0,      0.0,&
       0.0,    0.0,    0.0,    0.0,      0.0,      0.0,&
       0.0,    -6492., -6492., -6492.,  -6492.,  -6492.,&
       -6492., -6492., -4970., -4000./

        data f0&
      /0.0,    0.0,    0.0,    0.0,      0.0,      0.1,&
       0.0,    0.1,    0.1,    0.0,      1.0,      0.0,&
       0.0,    0.0,    0.0,    1.0,      0.0,      0.0,&
       0.0,    0.0,    0.0,    0.0,      0.0,      0.0,&
       0.0,    0.0,    0.0,    0.0,      0.0,      0.0,&
       0.0,    0.1,    0.0,    0.0,      0.0,      0.0,&
       0.0,    0.0,    0.0,    0.0,      0.0,      0.0,&
       0.0,    0.0,    0.0,    0.0,      0.0,      0.0,&
       0.0,    0.0,    0.0,    0.0,      0.0,      0.0,&
       0.0,    0.0,    0.0,    0.0,      0.0,      0.0,&
       0.0,    0.0,    0.0,    0.0,      0.0,      0.0,&
       0.0,    0.0,    0.0,    0.0,      0.0,      0.0,&
       0.0,    0.0,    0.0,    0.0/
        data rscale&
      /0.0,    0.0,    0.0,    0.0,      1.0,      1.0,&
       1.0,    0.0,    1.0,    1.0,      1.0,      1.0,&
       1.0,    1.0,    1.0,    1.0,      1.0,      1.0,&
       1.0,    1.0,    1.0,    1.0,      1.0,      1.0,&
       1.0,    1.0,    1.0,    1.0,      1.0,      1.0,&
       1.0,    1.0,    1.0,    1.0,      1.0,      1.0,&
       1.0,    1.0,    1.0,    1.0,      1.0,      1.0,&
       1.0,    1.0,    1.0,    1.0,      1.0,      1.0,&
       1.0,    1.0,    1.0,    1.0,      1.0,      1.0,&
       1.0,    1.0,    1.0,    1.0,      1.0,      1.0,&
       1.0,    1.0,    1.0,    1.0,      1.0,      1.0,&
       1.0,    1.0,    1.0,    1.0,      1.0,      1.0,&
       1.0,    1.0,    1.0,    0.0 /

       data diffrat&
      /1.00,   1.87,   1.42,   0.97,     1.29,    1.29,&
       0.0,    2.45,   1.62,   0.0,      1.63,    0.0,&
       0.0,    0.0,    0.0,    1.37,     1.25,    1.89,&
       0.0,    2.00,   0.0,    0.0 ,     1.29,    0.0,&
       0.0,    0.0,    1.60,   1.56,     0.0,     0.0,&
       0.0,    2.59,   2.00,   0.0,      2.00,    1.25,&
       1.80,   1.80,   2.26,   2.43,     2.45,    0.0,&
       0.0,    2.47,   2.72,   0.0,      0.0,     0.0,&
       0.0,    0.0,    2.00,   1.94,     1.97,    1.97,&
       1.97,   1.97,   0.0,    0.0,      0.0,     0.0,&
       0.0,    0.0,    0.0,    0.0,      0.0,     0.0,&
       1.00,   2.50,   2.50,   2.50,      2.50,   2.50,&
       2.50,   2.50,   3.34,   3.76/
!c
!c-----Baseline resistances are from Wesely (1989)
!c
      data rj  /9999.,  60., 120.,  70., 130., 100.,9999.,9999.,&
                 80., 100., 150.,9999.,9999.,9999.,9999., 250.,&
                500.,9999.,9999.,9999.,9999.,9999.,9999.,9999.,&
               9999.,9999., 250., 500.,9999.,9999.,9999.,9999.,&
               9999.,9999.,9999.,9999.,9999., 400., 800.,9999.,&
               9999.,9999.,9999.,9999.,9999., 120., 240., 140.,&
                250., 190.,9999.,9999., 160., 200., 300./
!c
      data rlu /9999.,2000.,2000.,2000.,2000.,2000.,9999.,9999.,&
               2500.,2000.,4000.,9999.,9000.,9000.,9000.,4000.,&
               8000.,9999.,9999.,9000.,9000.,9000.,9999.,9999.,&
               9000.,9000.,4000.,8000.,9999.,9999.,9000.,9000.,&
               9000.,9999.,9999.,9999.,9999.,6000.,9000.,9999.,&
               9999.,9000.,9000.,9000.,9999.,4000.,4000.,4000.,&
               2000.,3000.,9999.,9999.,4000.,4000.,8000./
     data rac / 100., 200., 100.,2000.,2000.,2000.,0.001,0.001,&
                300., 150., 200., 100., 150., 100.,1500.,2000.,&
               1700.,0.001,0.001, 200., 120., 140., 100.,  10.,&
                100.,1000.,2000.,1500.,0.001,0.001, 100.,  50.,&
                120., 100.,  10.,  10.,1000.,2000.,1500.,0.001,&
               0.001,  50.,  10.,  50., 100.,  50.,  80.,1200.,&
               2000.,1500.,0.001,0.001, 200.,  60., 120./
!c
      data rgss/ 400., 150., 350., 500., 500., 100.,0.001,1000.,&
               0.001, 220., 400., 400., 200., 350., 500., 500.,&
                100.,0.001,1000.,0.001, 300., 400., 400., 150.,&
                350., 500., 500., 200.,0.001,1000.,0.001, 200.,&
                400., 100., 100., 100., 100., 100., 100.,0.001,&
               1000., 100., 100.,  50., 500., 150., 350., 500.,&
                500., 200.,0.001,1000.,0.001, 250., 400./
!c
!c
      data rgso/ 300., 150., 200., 200., 200., 300.,2000., 400.,&
               1000., 180., 200., 300., 150., 200., 200., 200.,&
                300.,2000., 400., 800., 180., 200., 300., 150.,&
                200., 200., 200., 300.,2000., 400.,1000., 180.,&
                200., 600.,3500.,3500.,3500.,3500.,3500.,2000.,&
                400.,3500.,3500.,3500., 300., 150., 200., 200.,&
                200., 300.,2000., 400.,1000., 180., 200./
!c
     data rlcs/9999.,2000.,2000.,2000.,2000.,2000.,9999.,9999.,&
               2500.,2000.,4000.,9999.,9000.,9000.,9000.,2000.,&
               4000.,9999.,9999.,9000.,9000.,9000.,9999.,9999.,&
               9000.,9000.,3000.,6000.,9999.,9999.,9000.,9000.,&
               9000.,9999.,9999.,9999.,9000., 200., 400.,9999.,&
               9999.,9000.,9999.,9000.,9999.,4000.,4000.,4000.,&
               2000.,3000.,9999.,9999.,4000.,4000.,8000./
!c
      data rlco/9999.,1000.,1000.,1000.,1000.,1000.,9999.,9999.,&
               1000.,1000.,1000.,9999., 400., 400., 400.,1000.,&
                600.,9999.,9999., 400., 400., 400.,9999.,1000.,&
                400., 400.,1000., 600.,9999.,9999., 800., 600.,&
                600.,9999.,1000.,1000., 400.,1500., 600.,9999.,&
               9999., 800.,1000., 800.,9999.,1000., 500., 500.,&
               1500., 700.,9999.,9999., 600., 800., 800./
!c
      data dstress/ 1.0, 1.1, 1.5, 2.0, 3.5, 10.0 /
!c
!c-----Entry point
!c
!c
!c-----Loop over rows and columns
!c
      do 30 j = sy,ey
        do 20 i = sx,ex
!c-----To adjust Hg0 Henry constant, by chenhs
          if(ig.eq.75) then
            if(cellat(i,j).ge.0.) then !! North Hemisphere
              henry0(ig)=0.055
            else
              henry0(ig)=0.033
            endif
          endif
!c
!c
!c-----Determine season
!c
          mbin = month
          if (cellat(i,j).lt.0.) then
            mbin = mod(month+6,12)
            if (mbin.eq.0) mbin = 12
          endif
          latbin = 1
          if (abs(cellat(i,j)).gt.20.) then
            latbin = 2
          elseif (abs(cellat(i,j)).gt.35.) then
            latbin = 3
          elseif (abs(cellat(i,j)).gt.50.) then
            latbin = 4
          elseif (abs(cellat(i,j)).gt.75.) then
            latbin = 5
          endif
          if ((cellat(i,j).gt.50. .and. cellat(i,j).lt.75.) .and.&
             (cellon(i,j).gt.-15. .and. cellon(i,j).lt.15.)) latbin = 3
          isesn = iseason(latbin,mbin)
!c
!c-----Use input snow cover to set season, if specified
!c
          if (lrdsno .and. icdsno(i,j).eq.1) isesn = 4
!c
!c-----Calculate solar flux
!c
          solflux(i,j) = amax1(0.,solflux(i,j))
!c
!c-----Load local met variables
!c
          if( height(i,j,1).gt.0 ) then  !! by chenhs
            deltaz = height(i,j,1)/2.
          else
            deltaz=0.
          endif
          temp0 = tsurf(i,j) - 273.15
          prss0 = press(i,j,1) -&
                 2.*deltaz*(press(i,j,2) - press(i,j,1))/height(i,j,2)
          ucomp = (windu(i,j,1) + windu(i-1,j,1))/2.
          vcomp = (windv(i,j,1) + windv(i,j-1,1))/2.
          wind = sqrt(ucomp**2 + vcomp**2)
          wind = amax1(0.1,wind)
!       if(i==53.and.j==65.and.ig==11)print*,press(i,j,2),press(i,j,1),height(i,j,2)
!c
!c-----Determine surface wetness
!c
          iwet = 0
          if (pwc(i,j,1).gt.cwmin) then
            iwet = 2
          elseif (cwc(i,j,1).gt.cwmin) then
            iwet = 1
          else
            qwatr = 1.e-6*water(i,j,1)*18./28.8
            ev = qwatr*prss0/(qwatr + eps)
            es = e0*exp((lv/rv)*(1./273. - 1./tsurf(i,j)))
            rh = amin1(1.,ev/es)
            dew = (1. - rh)*(wind + 0.6)
            if (dew.lt.0.19) iwet = 1
          endif
!c
!c-----Use input drought stress, if specified
!c
          istress = 0
          if (lrddrt .and. icddrt(i,j).gt.0) istress = icddrt(i,j)
!c
!c-----Loop over land use; surface roughness for water is dependent on
!c     wind speed
         vdep(i,j) = 0.

          totland = 1.
!c
!            m=nint(landuse)
            m=landuse(i,j)
            z0 = z0lu(m,isesn)
            if (m.eq.7) z0 = amax1(z0,2.0e-6*wind**2.5)
!c
!c-----Use input surface roughness, if specified
!c
!c
!c-----Get surface layer micrometeorological parameters for this cell and
!c     landuse type
!c
            if (prss0.lt.0) then
              write(iout,'(//,a)') 'ERROR in DRYDEP:'
              write(iout,*) 'Invalid pressure value'
              write(iout,*) 'Cell   Height  Deltaz'
              write(iout,*) i,j,height(i,j,1),deltaz,press(i,j,1),press(i,j,2),prss0
            endif
            
          pbl = pbl_hgt(i,j)
            
            call micromet(tempk(i,j,1),tsurf(i,j),press(i,j,1),prss0,&
                         deltaz,wind,z0,pbl,ustar,psih,wstar,lstable)
             
!c
!c-----Loop over GAS species, and calculate deposition velocity for this cell,
!c     landuse, and current species.
!c     Use input drought stress code, if specified
!c
!           if(i==53.and.j==65.and.ig==11)print*,tempk(i,j,1),tsurf(i,j),press(i,j,1),prss0
               knh3  =  4
               khno3 =  2
               kso2  = 18
               ko3   = 11
!               IF(ig==18) THEN   !SO2
                henso20  = 1.0e+5
                tfactso2 = -3156
              call henryfnc(0,henso20,tfactso2,tsurf(i,j),7.,knh3,khno3,&
                         kso2,henso2)
!               ENDIF
              if (henry0(ig).gt.1.e-6) then
                call henryfnc(ig,henry0(ig),tfact(ig),tsurf(i,j),7.,knh3,&
                             khno3,kso2,henry)
                iflgso2 = 0
                iflgo3 = 0
                if (ig.eq.kso2) then
                  iflgso2 = 1
                  henry = henso2
                endif
                if (ig.eq.ko3) iflgo3 = 1
                call vd_gas(i,j,ig,m,istress,iwet,iflgso2,iflgo3,z0,&
                          deltaz,psih,ustar,diffrat(ig),henry,henso2,&
                          f0(ig),rscale(ig),temp0,dstress(0),solflux(i,j),&
                          rj(m,isesn),rlu(m,isesn),rac(m,isesn),&
                          rlcs(m,isesn),rlco(m,isesn),rgss(m,isesn),&
                          rgso(m,isesn),vd)
              else
                vd = 0.
              endif
              vdep(i,j) =  vd
!c
 20     continue
 30   continue
!c
      return
      end subroutine
      
      subroutine drydep_dust_salt(myid,month,tsurf,cellat,cellon,height,&
                            press,windu,windv,landuse,tempk,&
                    pbl_hgt,diam,vdep_aer,sx,ex,sy,ey,nzz,IA)
      integer :: sx,ex,sy,ey
      real tsurf(sx-1:ex+1,sy-1:ey+1),cellat(sx-1:ex+1,sy-1:ey+1),&
          cellon(sx-1:ex+1,sy-1:ey+1)
      real height(sx-1:ex+1,sy-1:ey+1,nzz),&     
          press(sx-1:ex+1,sy-1:ey+1,nzz),&
          windu(sx-1:ex+1,sy-1:ey+1,nzz),&
          windv(sx-1:ex+1,sy-1:ey+1,nzz),&
          tempk(sx-1:ex+1,sy-1:ey+1,nzz)
      real vdep_aer(sx-1:ex+1,sy-1:ey+1),pbl_hgt(sx-1:ex+1,sy-1:ey+1)
      logical lstable
      REAL :: iseason(5,12)
      REAL :: diam  
! diam : log-mean sectional aerosol diameter (m)
      REAL :: rhop ! AEROSOL DENISITY in g/m3
      INTEGER ::  landuse(sx-1:ex+1,sy-1:ey+1)
      REAL,DIMENSION(11,5) :: z0lu     
      INTEGER :: IA ! 2 DUST 1 SEASALT
!
      data pi/3.1415927/
      data cwmin /3.8E-5/
!
!c-----Season indices by month and latitude band
!c     Season Indices            Latitude Bands
!c     1 = summer                1 = <20    Tropical
!c     2 = autumn                2 = 20-35  Sub-tropical
!c     3 = winter w/o snow       3 = 35-50  Temperate
!c     4 = winter w/ snow        4 = 50-75  Cool
!c     5 = spring                5 = >75    Polar
!c                    Latitude Band
      data iseason / 1, 3, 3, 3, 3,& ! Jan
                    1, 5, 3, 3, 3,& ! Feb
                    1, 5, 5, 3, 3,& ! Mar
                    1, 5, 5, 5, 3,& ! Apr
                    1, 1, 5, 5, 3,& ! May
                    1, 1, 1, 1, 5,& ! Jun
                    1, 1, 1, 1, 1,& ! Jul
                    1, 1, 1, 1, 2,& ! Aug
                    1, 1, 2, 2, 3,& ! Sep
                    1, 2, 2, 2, 3,& ! Oct
                    1, 2, 2, 3, 3,& ! Nov
                    1, 2, 3, 3, 3/ ! Dec
!c
!c-----Surface roughness (m) as a function of 11 landuse categories
!c     and 5 seasons; based on AERMET model (ref EPA SCRAM website)
!c
      data z0lu&
      /1.0,0.20,0.100,1.3,1.3,1.30,0.0001,0.002,0.20,0.150,0.30,&
       1.0,0.05,0.010,0.8,1.3,1.05,0.0001,0.002,0.20,0.030,0.30,&
       1.0,0.05,0.010,0.8,1.3,1.05,0.0001,0.002,0.20,0.030,0.30,&
       1.0,0.01,0.001,0.5,1.3,0.90,0.0001,0.002,0.05,0.006,0.15,&
       1.0,0.03,0.050,1.0,1.3,1.15,0.0001,0.002,0.20,0.040,0.30/

       IF(IA==1) rhop = 2.2E06 ! SEASALT
       IF(IA==2) rhop = 2.65E06  ! DUST

      do 30 j = sy,ey
        do 20 i = sx,ex
!c
!c-----Determine season
!c
          mbin = month
          if (cellat(i,j).lt.0.) then
            mbin = mod(month+6,12)
            if (mbin.eq.0) mbin = 12
          endif
          latbin = 1
          if (abs(cellat(i,j)).gt.20.) then
          latbin = 2
          elseif (abs(cellat(i,j)).gt.35.) then
          latbin = 3
          elseif (abs(cellat(i,j)).gt.50.) then
          latbin = 4
          elseif (abs(cellat(i,j)).gt.75.) then
          latbin = 5
          endif
          if ((cellat(i,j).gt.50. .and. cellat(i,j).lt.75.) .and.&
             (cellon(i,j).gt.-15. .and. cellon(i,j).lt.15.)) latbin = 3
           isesn = iseason(latbin,mbin)
!c
!c-----Use input snow cover to set season, if specified
!c
!c
!c-----Load local met variables
          if( height(i,j,1).gt.0 ) then  !! by chenhs
            deltaz = height(i,j,1)/2.
          else
            deltaz=0.
          endif
          temp0 = tsurf(i,j) - 273.15
          prss0 = press(i,j,1) -&
                 2.*deltaz*(press(i,j,2) - press(i,j,1))/height(i,j,2)
          ucomp = (windu(i,j,1) + windu(i-1,j,1))/2.
          vcomp = (windv(i,j,1) + windv(i,j-1,1))/2.
          wind = sqrt(ucomp**2 + vcomp**2)
          wind = amax1(0.1,wind)
!c
!c
!c-----Loop over land use; surface roughness for water is dependent on        
!c     wind speed
         vdep_aer(i,j) = 0.
         totland = 1.
                   
        m=landuse(i,j)
        z0 = z0lu(m,isesn)
        if (m.eq.7) z0 = amax1(z0,2.0e-6*wind**2.5)
!c
!c-----Use input surface roughness, if specified
!c
!c
!c-----Get surface layer micrometeorological parameters for this cell
!and
!c     landuse type
!c
         if (prss0.lt.0) then
           write(iout,'(//,a)') 'ERROR in DRYDEP:'
           write(iout,*) 'Invalid pressure value'
           write(iout,*) 'Cell   Height  Deltaz'
           write(iout,*) i,j,height(i,j,1),deltaz,press(i,j,1),press(i,j,2),prss0
         endif

         pbl = pbl_hgt(i,j)
         call micromet(tempk(i,j,1),tsurf(i,j),press(i,j,1),prss0,&
                    deltaz,wind,z0,pbl,ustar,psih,wstar,lstable)
!cc

         call vd_dust(z0,deltaz,psih,ustar,diam,rhop,temp0,vd)      

         vdep_aer(i,j) =  vd
        
   20     continue
   30   continue
    
                            
      return
     end  subroutine

      subroutine vd_dust(z0,deltaz,psih,ustar,diam,rhop,ts,vd)
!c
!c----CAMx v4.42 070603
!c
!c     VD_AER calculates a deposition velocity for a specific aerosol size
!c     bin, grid cell, and land use category.  The parallel resistance approach
!c     of Slinn and Slinn (1980) is used, as implemented in UAM-AERO
!c     (STI, 1996).
!c
!c     Copyright 1996-2007
!c     ENVIRON International Corporation
!c
!c     Modifications:
!c        3/21/03       Modified Ra calculation to use layer 1 midpoint height
!c                      rather than default 10 m.
!c        9/18/03       Removed small bug in final Vd equation
!c
!c     Input arguments:
!c        z0                  surface roughness length (m)
!c        deltaz              Layer 1 midpoint height (m)
!c        psih                similarity stability correction term
!c        ustar               friction velocity (m/s)
!c        diam                log-mean sectional aerosol diameter (m)
!c        rhop                aerosol density (g/m3)
!c        ts                  surface temperature (K)
!c
!c     Output arguments:
!c        vd                  deposition velocity (m/s)
!c
!c     Routines called:
!c        none
!c
!c     Called by:
!c        DRYDEP

      data vk/0.4/, rmin/1.0/, xmfp/6.5e-8/, g/9.8/, vabs/1.81e-2/
      data boltz/1.38e-20/, pi/3.1415927/, vair/1.5e-5/
!c
!c-----Entry point
!c
!c-----Speed correction factor and sedimendation velocity
!c
      power = amin1(7.6,0.55*diam/xmfp)
      scf = 1. + (2.514 + 0.8*exp(-power))*xmfp/diam
      vsed = rhop*g*(diam*diam)*scf/(18.*vabs)
!c
!c-----Brownian diffusivity and Schmidt number
!c
      difbrwn = boltz*ts*scf/(3.*pi*vabs*diam)
      schmidt = vair/difbrwn
!c
!c-----Stokes number
!c
      stokes = vsed*(ustar*ustar)/(vair*g)
!c
!c-----Compute atmospheric resistance, RA
!c
      ra = (alog(deltaz/z0) - psih)/(vk*ustar)
      ra = amax1(ra,rmin)
!c
!c-----Compute the deposition layer resistance, RD
!c
      sc23 = schmidt**(-2./3.)
      power = -3./stokes
      if (power.lt.-37.) then
        xinert = 10.**(-37.)
      else
        xinert = 10.**(power)
      endif
      rd = 1./(ustar*(sc23 + xinert))
      rd = amax1(rd,rmin)
!c-----Final deposition velocity for this cell, land use, and aerosol size
!c
      vd = vsed + 1./(ra + rd + ra*rd*vsed)
!c
      return
      end subroutine
      
     subroutine vd_gas(i,j,ig,ilu,istress,iwet,iso2,io3,z0,deltaz,psih,ustar,&
                       diffrat,henry,henso2,f0,rscale,ts,dstress,&
                       solflux,rj,rlu,rac,rlcs,rlco,rgss,rgso,vd)
!c
!c----CAMx v4.42 070603
!c
!c     VD_GAS calculates a deposition velocity for a specific gas species,
!c     grid cell, and land use category.  The parallel resistance approach of
!c     Wesely and Hicks (1977) is used with the improvements of Wesely (1989).
!c     Surface resistance (rs) to water (landuse 7) is determined following
!c     Sehmel (1980), as implemented in UAM-AERO (STI, 1996).
!c
!c     Copyright 1996-2007
!c     ENVIRON International Corporation
!c
!c     Modifications:
!c        8/18/03       Relate drought stress to PDI index
!c        3/21/03       Modified Ra calculation to use layer 1 midpoint height
!c                      rather than default 10 m; changed drought stress effects
!c                      to stomatal resistance to be consistent with effects in
!c                      GLOBEIS
!c        3/26/03       Added scaling factor to surface resistance (provided
!c                      on chemparam file)
!c
!c     Input arguments:
!c        ilu                 land use index
!c        istress             vegetation drought stress index
!c        iwet                surface wetness index
!c                            0 = dry
!c                            1 = dew wetted
!c                            2 = rain wetted
!c        iso2                SO2 species flag (1=SO2,0=other)
!c        io3                 O3 species flag (1=O3,0=other)
!c        z0                  surface roughness length (m)
!c        deltaz              Layer 1 midpoint height (m)
!c        psih                similarity stability correction term
!c        ustar               friction velocity (m/s)
!c        diffrat             ratio of molecular diffusivity of water to species
!c        henry               Henry's Law constant (M/atm)
!c        henso2              Henry's Law constant for SO2 (M/atm)
!c        f0                  normalized reactivity parameter
!c        rscale              user-defined surface resistance scaling factor
!c        ts                  surface temperature (C)
!c        dstress             adjustment factors for drought stress
!c        solflux             Solar radiation flux (W/m2)
!c        rj                  baseline minimum stomatal resistance (s/m)
!c        rlu                 baseline upper canopy (cuticle) resistance (s/m)
!c        rac                 baseline canopy height/density resistance (s/m)
!c        rlcs                baseline SO2 lower canopy resistance (s/m)
!c        rlco                baseline O3 lower canopy resistance (s/m)
!c        rgss                baseline SO2 ground surface resistance (s/m)
!c        rgso                baseline O3 ground surface resistance (s/m)
!c
!c     Output arguments:
!c        vd                  deposition velocity (m/s)
!c
!c     Routines called:
!c        none
!c
!c     Called by:
!c        DRYDEP
     real    dstress(0:5)
!c
      data vk/0.4/, rmin/1.0/, rmax/1.e5/, d1/2./, d2/0.667/
      data vair/1.5e-5/, diffh2o/2.3e-5/
!c
!c-----Entry point
!c
!c-----Compute atmospheric resistance, RA
!c
      ra = (alog(deltaz/z0) - psih)/(vk*ustar)
      ra = amax1(ra,rmin)
!c
!c-----Compute the deposition layer resistance, RD
!c
      schmidt = vair*diffrat/diffh2o
      rd = d1*schmidt**d2/(vk*ustar)
      rd = amax1(rd,rmin)
!c
!c-----Compute the surface layer resistance over water, RS
!c
      if (ilu.eq.7) then
        rs = 1./(3.9e-5*henry*ustar*(ts + 273.15))
        rs = amax1(rs,rmin)
        goto 100
      endif
!c
!c
!c-----Compute stomatal resistance, RST
!c     Adjust for vegetation drought stress
!c
      rst = rmax
      if (ts.gt.0. .and. ts.lt.40.) then
        rst = diffrat*rj*(1. + (200./(solflux + 0.1))**2)*&
             (400./(ts*(40.-ts)))
        rst = rst * dstress(istress)
      endif
      if (iwet.gt.0) rst = 3.*rst
      rst = amin1(rmax,rst)
!c
!c-----Compute mesophyll resistance, RM
!c
      rm = 1./(henry/3000. + 100.*f0)
      rm = amax1(rmin,rm)
      rm = amin1(rmax,rm)
!c
!c-----Compute upper canopy resistance, RUC
!c     Adjust for surface wetness
!c
      if (iwet.eq.0) then
        ruc = rlu/(henry/henso2 + f0)
        ruc = amin1(rmax,ruc)
      else
        if (iwet.eq.1) then
          rlus = 100.
          if (ilu.eq.1) rlus = 50.
!c
!c  --- original equations from Wesely 89 ---
!c         rluo = 1./(1./3000. + 1./(3.*rlu))
!c
          rluo = 1000. + rlu
        else
!c
!c  --- original equations from Wesely 89 ---
!c         rlus = 1./(1./5000. + 1./(3.*rlu))
!c
          rlus = 2000. + rlu
          if (ilu.eq.1) rlus = 50.
          rluo = 1./(1./1000. + 1./(3.*rlu))
        endif
        if (iso2.eq.1) then
          ruc = rlus
        elseif (io3.eq.1) then
          ruc = rluo
        else
!c
!c  --- original equations from Wesely 89 ---
!         ruc = 1./(1./(3.*rlu) + 1.e-7*henry + f0/rluo)
!c
       ruc = 1./(henry/(henso2*rlus) + f0/rluo)

          ruc = amin1(rmax,ruc)
        endif
      endif
!         if(i==53.and.j==65.and.ig==11) print*,ruc,rmax,henso2
!c
!c-----Compute resistance to buoyant convection, RDC
!c     (assume effect of local terrain slope is non-resolvable; this
!c     factor is set to 1)
!c
      rdc = 100.*(1. + 1000./(solflux + 10.))
      rdc = amin1(rmax,rdc)
!c
!c-----Compute resistance of exposed surfaces in lower canopy, RCL
!c
      rlc = 1./(henry/(henso2*rlcs) + f0/rlco)
      rlc = amin1(rmax,rlc)
!c
!c-----Compute resistance of ground surface, RGS
!c
      rgs = 1./(henry/(henso2*rgss) + f0/rgso)
      rgs = amin1(rmax,rgs)
!c
!c-----Determine total surface resistance over land, RS
!c
      rs = 1./(rst + rm) + 1./ruc + 1./(rdc + rlc) + 1./(rac + rgs)
      rs = amax1(rmin,1./rs)
!c
!c-----Scale surface resistance for acids according to user-definition
!c
 100  continue
      rs = rs*rscale
!c
!c-----Final deposition velocity for this cell, land use, and species
      vd = 1./(ra + rd + rs)
!c
      if(ig==75) then
       if(vd.gt.0.2) then
        print*,vd,ra,rd,rs,deltaz,psih,ustar
       endif
      endif
      return
      end subroutine

      subroutine micromet(temp,temp0,press,press0,deltaz,wind,z0,pbl,&
                         ustar,psih,wstar,lstable)
!c
!c----CAMx v4.42 070603
!c
!c     MICROMET calculates surface layer micro-meteorological flux-gradient
!c     relationships and variables based on Louis (1979)
!c
!c     Copyright 1996-2007
!c     ENVIRON International Corporation
!c
!c     Modifications:
!c        8/15/02    Specified limits for denominators from 1e-10 to 1e-20
!c        3/21/03    Modified Z/L calculation to use input midpoint height
!c                   (instead of default 10 m) and to increase range from
!c                   (-1 to +1) to (-2.5 to +1.5)
!c        8/20/03    added calculation of w* and logical for stability
!c
!c     Input arguments:
!c        temp                Layer 1 temperature (K)
!c        temp0               Surface temperature (K)
!c        press               Layer 1 pressure (mb)
!c        press0              Surface pressure (mb)
!c        deltaz              Layer 1 midpoint height (m)
!c        wind                Layer 1 total wind speed (m/s)
!c        z0                  Surface roughness length (m)
!c        pbl                 PBL depth (m)
!c
!c     Output arguments:
!c        ustar               friction velocity (m/s)
!c        psih                Similarity stability correction term
!c        wstar               convective velocity scale (m/s)
!c        lstable             stability logical (T=stable)
!c
!c     Routines called:
!c        none
!c
!c     Called by:
!c        DRYDEP
!c        PIGGROW
!c
     logical lstable
      data vk/0.4/, g/9.8/, gamma/0.286/
!c
!c-----Entry point
!c
!c-----Calculate potential temperature and richardson number
!c
      theta = temp*(1000./press)**gamma
      theta0 = temp0*(1000./press0)**gamma
      dtheta = theta - theta0
      thetabar = (theta + theta0)/2.
      ri = (g/thetabar)*deltaz*dtheta/(wind**2 + 1.e-20)
!c
!c-----Determine stability functions
!c
      zscale = vk/alog(deltaz/z0)
      if (ri.le.0.) then
        lstable = .false.
        cm    = 69.56*sqrt(deltaz/z0)*zscale**2
        ch    = 49.82*sqrt(deltaz/z0)*zscale**2
        fm    = 1. - 9.4*ri/(1. + cm*sqrt(abs(ri)))
        fh    = 1. - 9.4*ri/(1. + ch*sqrt(abs(ri)))
      else
        lstable = .true.
        fm = 1./((1. + 4.7*ri)**2)
        fh = fm
      endif
!c
!c-----Calculate micromet variables
!c
      ustar2 = fm*(wind*zscale)**2
      ustar2 = amax1(1.e-20,ustar2)
      ustar = sqrt(ustar2)
      thstar = 1.35*zscale**2*wind*dtheta*fh/ustar
      el = ustar2*temp/(vk*g*thstar + 1.e-20)
      elabs = abs(el)
      wstar = 0.
      if (el.lt.0.) wstar = (pbl*ustar**3./(vk*elabs))**(1./3.)

      if (elabs.ge.1.e4) then
        psih = 0.0
      elseif (el.lt.0.) then
        zoverl = amin1(2.5,deltaz/elabs)
        zmel = alog(zoverl)
        psih = exp(0.598 + 0.39*zmel - 0.090*zmel*zmel)
      else
        zoverl = amin1(1.5,deltaz/elabs)
        psih = -5.*zoverl
      endif
!c
      return
      end subroutine

     subroutine henryfnc(ispc,hlaw0,tfact,temp,ph,knh3,khno3,kso2,hlaw)
!c
!c----CAMx v4.42 070603
!c
!c     HENRYFNC calculates temperature and dissociation adjustments to
!c     baseline Henry's Law constants.
!c
!c     Copyright 2006-2007
!c     ENVIRON International Corporation
!c
!c     Modifications:
!c        None
!c
!c     Input arguments:
!c        ispc                Species index
!c        hlaw0               Baseline Henry's Law constant @298K (M/atm)
!c        tfact               temperature factor
!c        temp                ambient temperature (K)
!c        ph                  pH of liquid
!c        knh3                pointer to NH3
!c        khno3               pointer to HNO3
!c        kso2                pointer to SO2
!c
!c     Output arguments:
!c        hlaw                Adjusted Henry's Law constant (M/atm)
!c
!c     Routines called:
!c        None
!c
!c     Called by:
!c        WETDEP
!c        DRYDEP
!c
      implicit none
      integer ispc,knh3,khno3,kso2
      real hlaw0,tfact,temp,ph,hlaw

      real diss1,diss2
!c
      hlaw = hlaw0*exp(tfact*(1./298. - 1./temp))
      if (ispc.eq.knh3) then
        diss1 = 10.**(-189.1/temp - 4.117)
        diss2 = 10.**(-5839.5/temp - 9.7618*alog(temp) + 61.206)
        hlaw = hlaw*(1. + (diss1/diss2)*10.**(-ph))
      elseif (ispc.eq.khno3) then
       diss1 = 15.4
        hlaw = hlaw*(1. + diss1/(10.**(-ph)))
      elseif (ispc.eq.kso2) then
        diss1 = 10.**(853./temp)/54950.
        diss2 = 10.**(621.9/temp)/1.897e+9
        hlaw = hlaw*(1. + diss1/(10.**(-ph)) +&
                         diss1*diss2/(10.**(-2.*ph)))
      endif

      return
      end subroutine
end module drydep
