'open usgs.ctl'
'open bats.ctl'
'open sib2.ctl'
'open igbp.ctl'
'open pft1.ctl'
'open pft2.ctl'

'set x 1'
'set y 1'
'set z 1'
'set t 43 162'

'set grid off'
'set grads off'

'set xlint 12'
'set ylint 0.2'
'set vrange 2.4 5.4'

'set cmark  0'
'set ccolor 2'
'd tloop(aave(assim.1,lon=-179.5,lon=179.5,lat=-59.5,lat=89.5)*1.0e6)'

'set cmark  0'
'set ccolor 3'
'd tloop(aave(assim.2,lon=-179.5,lon=179.5,lat=-59.5,lat=89.5)*1.0e6)'

'set cmark  0'
'set ccolor 4'
'd tloop(aave(assim.3,lon=-179.5,lon=179.5,lat=-59.5,lat=89.5)*1.0e6)'

'set cmark  0'
'set ccolor 5'
'd tloop(aave(assim.4,lon=-179.5,lon=179.5,lat=-59.5,lat=89.5)*1.0e6)'

'set cmark  0'
'set ccolor 9'
'd tloop(aave(assim.5,lon=-179.5,lon=179.5,lat=-59.5,lat=89.5)*1.0e6)'

'set cmark  0'
'set ccolor 10'
'd tloop(aave(assim.6,lon=-179.5,lon=179.5,lat=-59.5,lat=89.5)*1.0e6)'

'draw ylab  CO2 Assimilation Rate(1.0e-6mol/(m^2*s))\'
'draw title \Global Land Monthly Mean\'
'cbar_line  -p -c 2 3 4 5 9 10 -l 1 1 1 1 1 1 -t "USGS" "BATS" "SiB2" "IGBP" "PFT" "PFT-DGVM"'

'printim assim.gif gif x1024 y768 white'

'close 6'
'close 5'
'close 4'
'close 3'
'close 2'
'close 1'
