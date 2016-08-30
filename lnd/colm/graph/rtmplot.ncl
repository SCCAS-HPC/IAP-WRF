
load "$NCARG_ROOT/lib/ncarg/nclex/gsun/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclex/gsun/gsn_csm.ncl"

begin

   fn1 = "/home/jidy/gfs/colm_rtm/global/global-RTM-1998-06"
   fn2 = "/home/jidy/gfs/clm3run/clm3run.clm2.h0.1998-06.nc"

   f2 = addfile(fn2,"r")

   lnd_globxy1 = fbinrecread(fn1,0,(/360,720/),"float")
   ocn_globxy1 = fbinrecread(fn1,1,(/360,720/),"float")

   lnd_globxy2 = f2->QCHANR(0,:,:)
   ocn_globxy2 = f2->QCHOCNR(0,:,:)

   lnd_globxy1@long_name = "RTM river flow"
   lnd_globxy1@units = "1000 m3/s"
   lnd_globxy1@_FillValue = -9999.

   ocn_globxy1@long_name = "RTM river discharge into ocean"
   ocn_globxy1@units = "1000 m3/s"
   ocn_globxy1@_FillValue = -9999.

   lnd_globxy1 = lnd_globxy1*1.0E-3
   ocn_globxy1 = ocn_globxy1*1.0E-3
   lnd_globxy2 = lnd_globxy2*1.0E-3
   ocn_globxy2 = ocn_globxy2*1.0E-3

   lnd_globxy2@units = "1000 m3/s"
   ocn_globxy2@units = "1000 m3/s"

   rtmlat = new((/360/),"float")
   rtmlon = new((/720/),"float")

   do i = 0, 359
      rtmlat(i) = (i+1)*0.5 - 90.
   end do

   do i = 0, 719
      rtmlon(i) = (i+1)*0.5 - 180.
   end do

   rtmlat@units = "degrees_N"
   rtmlon@units = "degrees_E"

   lnd_globxy1!0 = "lat"
   lnd_globxy1!1 = "lon"
   lnd_globxy1&lat = rtmlat
   lnd_globxy1&lon = rtmlon

   ocn_globxy1!0 = "lat"
   ocn_globxy1!1 = "lon"
   ocn_globxy1&lat = rtmlat
   ocn_globxy1&lon = rtmlon

   wks = gsn_open_wks("pdf","rtm")
 
   res = True

   res@cnMonoLineColor = False

   res@cnFillOn = True

   res@cnFillMode = "RasterFill"
   res@cnLinesOn = False
   res@cnLevelSelectionMode = "ExplicitLevels"
   res@cnLevels = (/0.1, 0.3, 0.6, 1, 2, 4, 6, 8, 10, 13, 16, 20, 30, 40/)

   gsn_define_colormap(wks,"spread_15lev")

   plot = gsn_csm_contour_map_ce(wks,lnd_globxy1,res)

   plot = gsn_csm_contour_map_ce(wks,lnd_globxy2,res)

   delete(res@cnLevels)

  ;res@cnLevels = (/0.001, 0.002, 0.003, 0.004, 0.1, 1, 4,  8, 12, 16, 19/)
   res@cnLevels = (/0.01, 0.02, 0.03, 0.04, 0.07, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 1/)

   gsn_define_colormap(wks,"spread_15lev")

   plot = gsn_csm_contour_map_ce(wks,ocn_globxy1,res)

   plot = gsn_csm_contour_map_ce(wks,ocn_globxy2,res)

end
