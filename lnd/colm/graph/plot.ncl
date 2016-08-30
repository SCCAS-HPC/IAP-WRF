
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"

begin

;f1 = addfile("~/swgfs/colm_esm2/global-1492-12.nc","r")
;f1 = addfile("~/swgfs/colm_esm2/global-1500-12.nc","r")
;f1 = addfile("~/swgfs/colm_esm2/global-1520-12.nc","r")
;f1 = addfile("~/swgfs/colm_esm2/global-1540-12.nc","r")
;f1 = addfile("~/swgfs/colm_esm2/global-1550-12.nc","r")
;f1 = addfile("~/swgfs/colm_esm2/global-1580-12.nc","r")
;f1 = addfile("~/swgfs/colm_esm2/global-1600-12.nc","r")

 f1 = addfile("~/swgfs/CO2C/exe/lnd/CoLM-T42-0001-01-01-01-30.nc","r")


tg = f1->tg

lai = f1->lai

fsena = f1->fsena

lfevpa = f1->lfevpa

assim = f1->assim

respc = f1->respc

fmicr = f1->fmicr

solar = f1->solar

;laimx = f1->laimx

wks=gsn_open_wks("pdf","colmc")

res = True

res@cnFillOn = True

res@tiMainString = "tg"

plot=gsn_csm_contour_map_ce(wks,tg(0,:,:),res)

res@tiMainString = "solar"

plot=gsn_csm_contour_map_ce(wks,solar(0,:,:),res)

res@tiMainString = "lai"

plot=gsn_csm_contour_map_ce(wks,lai(0,:,:),res)

res@tiMainString = "fsena"

plot=gsn_csm_contour_map_ce(wks,fsena(0,:,:),res)

res@tiMainString = "lfevap"

plot=gsn_csm_contour_map_ce(wks,lfevap(0,:,:),res)

res@tiMainString = "assim"

;res@cnLevelSelectionMode = "ExplicitLevels"
;res@cnLevels = (/1.0E-6, 1.E-5, 0.1, 1.0, 2.0, 5.0, 10., 20., 50., 100., 1000./)

plot=gsn_csm_contour_map_ce(wks,assim(0,:,:),res)

;delete(res@cnLevels)

res@tiMainString = "respc"

plot=gsn_csm_contour_map_ce(wks,respc(0,:,:),res)

res@tiMainString = "fmicr"

plot=gsn_csm_contour_map_ce(wks,fmicr(0,:,:),res)

do i = 0, 13
   res@tiMainString = "laimx(" + sprinti("%i",i) + ")"
;  plot=gsn_csm_contour_map_ce(wks,laimx(0,i,:,:),res)
end do

end
