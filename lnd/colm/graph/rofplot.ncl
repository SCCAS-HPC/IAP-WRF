
load "$NCARG_ROOT/lib/ncarg/nclex/gsun/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclex/gsun/gsn_csm.ncl"

load "clmshare.ncl"

begin

   fn1 = "/home/jidy/gfs/colm_rtm/global/global-1998-06"
   fn2 = "/home/jidy/gfs/clm3run/clm3run.clm2.h0.1998-06.nc"

   f2 = addfile(fn2,"r")

   qover = f2->QOVER(0,:,:)
   qdrai = f2->QDRAI(0,:,:)

   qover = qover+qdrai

   qover = qover*86400

   wks = gsn_open_wks("pdf","rof")
 
   res = True

   res@cnMonoLineColor = False

   res@cnFillOn = True

  ;res@cnFillMode = "RasterFill"
  ;res@cnLinesOn = False
  ;res@cnLevelSelectionMode = "ExplicitLevels"
  ;res@cnLevels = (/0.1, 0.3, 0.6, 1, 2, 4, 6, 8, 10, 13, 16, 20, 30, 40/)

  ;res@cnLevels = (/0.001, 0.002, 0.003, 0.004, 0.1, 1, 4,  8, 12, 16, 19/)
  ;res@cnLevels = (/0.01, 0.02, 0.03, 0.04, 0.07, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 1/)

   gsn_define_colormap(wks,"spread_15lev")

   plot = gsn_csm_contour_map_ce(wks,qover,res)

end
