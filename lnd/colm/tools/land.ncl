begin

f1 = addfile("~/swgfs/Princeton_30MIN_GLOBALx/tas/tas_30min_1990-1990.nc","r")

lon = f1->lon
lat = f1->lat
land = f1->land

dimsz = dimsizes(lon)
nlon = dimsz(0)

dimsz = dimsizes(lat)
nlat = dimsz(0)

dimsz = dimsizes(land)
nland = dimsz(0)

metmask = new((/nlat,nlon/),"integer")

metmask = 0

do k = 0, nland-1

   j = land(k)/nlon+1
   i = land(k)%nlon

   if(i.eq.0) then
      i = nlon
      j = j-1
   end if

 ; transform to Fortran style.
   i = i-1
   j = j-1

   metmask(j,i) = 1

end do

fbinrecwrite("metmask_nc",-1,metmask)

end
