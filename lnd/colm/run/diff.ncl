
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"

begin

;wks = gsn_open_wks("x11","mask")

 res = True

 res@cnFillOn = True
 var1 = new((/18,36/),"float")
 var2 = new((/18,36/),"float")

;do i = 0, 50  ;tg
;do i = 0, 85  ;tm
;do i = 0, 13  ; sabg
 do k = 0, 91
   var1 = fbinrecread("../output1/global-1990-01-01-01-00", k, (/18,36/),"float")
   var2 = fbinrecread("../output2/global-1990-01-01-01-00", k, (/18,36/),"float")

;  plot = gsn_contour(wks, var2-var1, res)

   do j = 0, 17
   do i = 0, 35
      if(abs(var2(j,i)-var1(j,i)) .gt. 1.0E-8) then
         print((/k,j,i,var1(j,i),var2(j,i)/))
      end if
   end do
   end do
 end do

end
