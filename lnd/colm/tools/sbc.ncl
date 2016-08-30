load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

begin

filename = "~/swgfs/t42_offline_restart5/T42-0066-12.nc"

fhist = addfile(filename,"r")

system("rm -f T42-sbc.bin")

setfileoption("bin","WriteByteOrder","BigEndian")

fbinrecwrite("T42-sbc.bin",-1,flt2dble(fhist->albvdir(0,:,:)))
fbinrecwrite("T42-sbc.bin",-1,flt2dble(fhist->albvdif(0,:,:)))
fbinrecwrite("T42-sbc.bin",-1,flt2dble(fhist->albndir(0,:,:)))
fbinrecwrite("T42-sbc.bin",-1,flt2dble(fhist->albndif(0,:,:)))
fbinrecwrite("T42-sbc.bin",-1,flt2dble(fhist->scv(0,:,:)))
fbinrecwrite("T42-sbc.bin",-1,flt2dble(fhist->trad(0,:,:)))

end
