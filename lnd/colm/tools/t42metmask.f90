program main

   use precision

   implicit none

   integer, parameter :: nlon = 128
   integer, parameter :: nlat =  64

   integer, parameter :: nlon_met = 360
   integer, parameter :: nlat_met = 180

   integer numlon(nlat)
   real(r8) lat(nlat)
   real(r8) lon(nlon)
   real(r8) lats(nlat+1)
   real(r8) lonw(nlon+1)
   real(r8) lonw2d(nlon+1,nlat)
   real(r8) area(nlon,nlat)
   real(r8) latixy(nlon,nlat)
   real(r8) longxy(nlon,nlat)
   real(r8) landfrac(nlon,nlat)
   integer  landmask(nlon,nlat)

   real(r8) lon_met(nlon_met)
   real(r8) lat_met(nlat_met)
   integer  landmask_met(nlon_met,nlat_met)

   real minlon, maxlon
   integer i, j, i1, i2, j1, j2, n

   integer :: lugrid = 10

   open(lugrid,file="./gridat",form="unformatted",convert="big_endian",status="old")
   read(lugrid) numlon(:)
   read(lugrid) lats(:)
   read(lugrid) lonw2d(:,:)
   read(lugrid) area(:,:)
   read(lugrid) latixy(:,:)
   read(lugrid) longxy(:,:)
   read(lugrid) landfrac(:,:)
   read(lugrid) landmask(:,:)
   close(lugrid)

   open(lugrid,file="./landmask_t42",form="unformatted",convert="little_endian",status="unknown")
   write(lugrid) landmask
   close(lugrid)

   do i = 1, nlat_met
      lat_met(i) = -89.5+(i-1)
   end do

   do i = 1, nlon_met
      lon_met(i) = -179.5+(i-1)
   end do

   lonw(:) = lonw2d(:,1)
   minlon = minval(lonw)
   maxlon = maxval(lonw)

 ! print *, 'lats', lats
 ! print *, 'lonw', lonw

   print *, 'min & max longitude of LSM', minlon, maxlon

   if(maxlon.gt.180.) then
      do i = 1, nlon_met
         if(lon_met(i).lt.minlon) lon_met(i) = lon_met(i)+360.
         if(lon_met(i).gt.maxlon) stop 'error in longitude convertion'
      end do
   end if

   landmask_met(:,:) = 0

   do j1 = 1, nlat
   do i1 = 1, nlon
      n = 0
      do j2 = 1, nlat_met
         if(lats(j1).lt.lat_met(j2) .and. lat_met(j2).lt.lats(j1+1)) then  !* land grids in S->N ascending order
            do i2 = 1, nlon_met
               if(lonw(i1).lt.lon_met(i2) .and. lon_met(i2).lt.lonw(i1+1)) then
                  if(landmask(i1,j1).eq.1) landmask_met(i2,j2) = 1
                  n = n+1
               end if
            end do
         end if
      end do
      if(n.eq.0) print *, 'error', lonw(i1),lonw(i1+1),lats(j1),lats(j1+1)
   end do
   end do

   print *, 'total land model grids', sum(landmask)
   print *, 'total metorological data points', sum(landmask_met)

   do i2 = 1, nlon_met
      if(lon_met(i2).gt.178.59375 .and. lon_met(i2).lt.181.40625) then
         do j2 = 1, nlat_met
            if(lat_met(j2).gt.-19.5342726223644 .and. lat_met(j2).lt.-16.7436678748654) then
          !    print *, lon_met(i2), lat_met(j2), landmask_met(i2,j2)
            end if
         end do
      end if
   end do

   do i2 = 1, nlon_met
      if(lon_met(i2).gt.178.59375 .and. lon_met(i2).lt.181.40625) then
         do j2 = 1, nlat_met
            if(lat_met(j2).gt.-16.7436678748654 .and. lat_met(j2).lt.-13.9530604373611) then
          !    print *, lon_met(i2), lat_met(j2), landmask_met(i2,j2)
            end if
         end do
      end if
   end do

!  178.593750000000    181.406250000000       -19.5342726223644       -16.7436678748654
!  178.593750000000    181.406250000000       -16.7436678748654       -13.9530604373611

   open(lugrid,file="./metmask",form="unformatted",status="unknown")
   write(lugrid) ((landmask_met(i,nlat_met-j+1),i=1,nlon_met),j=1,nlat_met)
   close(lugrid)

end program
