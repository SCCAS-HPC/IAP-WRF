
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"

begin

 mask1 = new((/180,360/),"integer")
 mask1 = fbinrecread("../mksrfdata/landmask", 0, (/180,360/),"integer")

 mask2 = new((/180,360/),"integer")
 mask2 = fbinrecread("./gmask",0,(/180,360/),"integer")

 mask3 = new((/180,360/),"integer") 
 mask3 = fbinrecread("../mkinidatacn/landmask2", 0, (/180,360/),"integer")

 wks = gsn_open_wks("x11","mask")

 res = True

 res@cnFillOn = True

 plot = gsn_contour(wks, mask2, res)

end
