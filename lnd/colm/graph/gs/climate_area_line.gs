'open ../ctl/CN.reg.ctl'

'set grid off'
'set grads off'

*'set xlint 100'
'set ylint 10'
'set vrange 0 100'

'set x 1'
'set y 1'
'set z 1'
'set t 1 500'

'set cmark  0'
'set ccolor 3'
'define tropical=asum((bdttro+bettro+c4)*soil*garea,lon=-179,lon=179,lat=-59,lat=89)'
'd tloop(tropical*1.E-6)'
'undefine tropical'

'set cmark  0'
'set ccolor 4'
'define temperate=asum((bdttem+nettem+bettem+c3)*soil*garea,lon=-179,lon=179,lat=-59,lat=89)'
'd tloop(temperate*1.E-6)'
'undefine temperate'

'set cmark  0'
'set ccolor 7'
'define boreal=asum((netbor+bdtbor+c3arc)*soil*garea,lon=-179,lon=179,lat=-59,lat=89)'
'd tloop(boreal*1.E-6)'
'undefine boreal'

'set t 1 500'
'set cmark  0'
'set ccolor 2'
'define baresoil=asum(bare*soil*garea,lon=-179,lon=179,lat=-59,lat=89))'
'd tloop(baresoil*1.E-6)'
'undefine baresoil'

'draw ylab  total area [1.E6 km2]\'
'draw title \Global climate type cover\'
'cbar_line  -p -c 3 4 7 2 -l 1 1 1 1 -t "tropical" "temperate" "boreal" "bare soil" '

'printim climate_area_line.gif gif x1024 y768 white'

'c'

'set ylint 10'
'set vrange 0 50'

'set cmark  0'
'set ccolor 2'
'define tbdttro=asum(bdttro*soil*garea,lon=-179,lon=179,lat=-59,lat=89)'
'd tloop(tbdttro*1.E-6)'
'undefine tbdttro'

'set cmark  0'
'set ccolor 3'
'define tbettro=asum(bettro*soil*garea,lon=-179,lon=179,lat=-59,lat=89)'
'd tloop(tbettro*1.E-6)'
'undefine tbettro'

'set cmark  0'
'set ccolor 4'
'define tc4=asum(c4*soil*garea,lon=-179,lon=179,lat=-59,lat=89)'
'd tloop(tc4*1.E-6)'
'undefine tc4'

'draw ylab  total area [1.E6 km2]\'
'draw title \Global tropical PFTs cover\'
'cbar_line  -p -c 2 3 4 -l 1 1 1 -t "BDTtro" "BETtro" "C4"'

'printim tropical_area_line.gif gif x1024 y768 white'

'c'

'set ylint 10'
'set vrange 0 50'

'set cmark  0'
'set ccolor 2'
'define tbdttem=asum(bdttem*soil*garea,lon=-179,lon=179,lat=-59,lat=89)'
'd tloop(tbdttem*1.E-6)'
'undefine tbdttem'

'set cmark  0'
'set ccolor 3'
'define tnettem=asum(nettem*soil*garea,lon=-179,lon=179,lat=-59,lat=89)'
'd tloop(tnettem*1.E-6)'
'undefine tnettem'

'set cmark  0'
'set ccolor 4'
'define tbettem=asum(bettem*soil*garea,lon=-179,lon=179,lat=-59,lat=89)'
'd tloop(tbettem*1.E-6)'
'undefine tbettem'

'set cmark  0'
'set ccolor 7'
'define tc3=asum(c3*soil*garea,lon=-179,lon=179,lat=-59,lat=89)'
'd tloop(tc3*1.E-6)'
'undefine tc3'

'draw ylab  total area [1.E6 km2]\'
'draw title \Global temperate PFTs cover\'
'cbar_line  -p -c 2 3 4 7 -l 1 1 1 1 -t "BDTtem" "NETtem" "BETtem" "C3" '

'printim temperate_area_line.gif gif x1024 y768 white'

'c'

'set ylint 10'
'set vrange 0 50'

'set cmark  0'
'set ccolor 2'
'define bnetbor=asum(netbor*soil*garea,lon=-179,lon=179,lat=-59,lat=89)'
'd tloop(bnetbor*1.E-6)'
'undefine bnetbor'

'set cmark  0'
'set ccolor 3'
'define bbdtbor=asum(bdtbor*soil*garea,lon=-179,lon=179,lat=-59,lat=89)'
'd tloop(bbdtbor*1.E-6)'
'undefine bbdtbor'

'set cmark  0'
'set ccolor 4'
'define bc3arc=asum(c3arc*soil*garea,lon=-179,lon=179,lat=-59,lat=89)'
'd tloop(bc3arc*1.E-6)'
'undefine bc3arc'

'draw ylab  total area [1.E6 km2]\'
'draw title \Global boreal PFTs cover\'
'cbar_line  -p -c 2 3 4 -l 1 1 1 -t "NETbor" "BDTbor" "C3-Arctic"'

'printim boreal_area_line.gif gif x1024 y768 white'

'close 1'
