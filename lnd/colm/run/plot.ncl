
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"

begin

 wks = gsn_open_wks("x11","mask")

 res = True

 res@cnFillOn = True
 var1 = new((/18,36/),"float")
 var2 = new((/18,36/),"float")

;do i = 0, 50  ;tg
;do i = 0, 85  ;tm
;do i = 0, 13  ;sabg
 do k = 0, 91
   var1 = fbinrecread("../output1/global-1990-01-01-01-00", k, (/18,36/),"float")

   var1@_FillValue = -9999.

   plot = gsn_contour(wks, var1, res)
 end do

end
