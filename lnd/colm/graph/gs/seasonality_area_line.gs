'open ../ctl/CN.reg.ctl'

'set grid off'
'set grads off'

*'set xlint 100'
'set ylint 20'
'set vrange 0 100'

'set x 1'
'set y 1'
'set z 1'
'set t 1 500'

'set cmark  0'
'set ccolor 3'
'define evergreen=asum((nettem+netbor+bettem+bettro)*soil*garea,lon=-179,lon=179,lat=-59,lat=89)'
'd tloop(evergreen*1.E-6)'
'undefine evergreen'

'set cmark  0'
'set ccolor 4'
'define deciduous=asum((bdttro+bdttem+bdtbor)*soil*garea,lon=-179,lon=179,lat=-59,lat=89)'
'd tloop(deciduous*1.E-6)'
'undefine deciduous'

'set cmark  0'
'set ccolor 7'
'define grass=asum((c3+c3arc+c4)*soil*garea,lon=-179,lon=179,lat=-59,lat=89)'
'd tloop(grass*1.E-6)'
'undefine grass'

'set t 1 500'
'set cmark  0'
'set ccolor 2'
'define baresoil=asum(bare*soil*garea,lon=-179,lon=179,lat=-59,lat=89))'
'd tloop(baresoil*1.E-6)'
'undefine baresoil'

'draw ylab  total area [1.E6 km2]\'
'draw title \Global seasonality cover\'
'cbar_line  -p -c 3 4 7 2 -l 1 1 1 1 -t "evergreen" "deciduous" "grass" "bare soil" '

'printim seasonality_area_line.gif gif x1024 y768 white'

'c'

*plot evergreen trees seperately

'set vrange 0 24'

'set cmark  0'
'set ccolor 3'
'define ebettro=asum(bettro*soil*garea,lon=-179,lon=179,lat=-59,lat=89)'
'd tloop(ebettro*1.E-6)'
'undefine ebettro'

'set cmark  0'
'set ccolor 4'
'define ebettem=asum(bettem*soil*garea,lon=-179,lon=179,lat=-59,lat=89)'
'd tloop(ebettem*1.E-6)'
'undefine ebettem'

'set cmark  0'
'set ccolor 7'
'define enettem=asum(nettem*soil*garea,lon=-179,lon=179,lat=-59,lat=89)'
'd tloop(enettem*1.E-6)'
'undefine enettem'

'set cmark  0'
'set ccolor 2'
'define enetbor=asum(netbor*soil*garea,lon=-179,lon=179,lat=-59,lat=89)'
'd tloop(enetbor*1.E-6)'
'undefine enetbor'

'draw ylab  total area [1.E6 km2]\'
'draw title \Global evergreen trees cover\'
'cbar_line  -p -c 3 4 7 2 -l 1 1 1 1 -t "BETtro" "BETtem" "NETtem" "NETbor" '

'printim evergreen_area_line.gif gif x1024 y768 white'

'c'

*plot deciduous trees seperately

'set vrange 0 14'

'set cmark  0'
'set ccolor 3'
'define dbdttro=asum(bdttro*soil*garea,lon=-179,lon=179,lat=-59,lat=89)'
'd tloop(dbdttro*1.E-6)'
'undefine dbdttro'

'set cmark  0'
'set ccolor 4'
'define dbdttem=asum(bdttem*soil*garea,lon=-179,lon=179,lat=-59,lat=89)'
'd tloop(dbdttem*1.E-6)'
'undefine dbdttem'

'set cmark  0'
'set ccolor 7'
'define dbdtbor=asum(bdtbor*soil*garea,lon=-179,lon=179,lat=-59,lat=89)'
'd tloop(dbdtbor*1.E-6)'
'undefine dbdtbor'

'draw ylab  total area [1.E6 km2]\'
'draw title \Global deciduous trees cover\'
'cbar_line  -p -c 3 4 7  -l  1 1 1 -t "BDTtro" "BDTtem" "BDTbor" '

'printim deciduous_area_line.gif gif x1024 y768 white'

'c'

*plot grasses separately

'set vrange 0 50'

'set cmark  0'
'set ccolor 3'
'define grassc3=asum(c3*soil*garea,lon=-179,lon=179,lat=-59,lat=89)'
'd tloop(grassc3*1.E-6)'
'undefine grassc3'

'set cmark  0'
'set ccolor 2'
'define grassc4=asum(c4*soil*garea,lon=-179,lon=179,lat=-59,lat=89)'
'd tloop(grassc4*1.E-6)'
'undefine grassc4'

'set cmark  0'
'set ccolor 7'
'define grassc3arc=asum(c3arc*soil*garea,lon=-179,lon=179,lat=-59,lat=89)'
'd tloop(grassc3arc*1.E-6)'
'undefine grassc3arc'

'draw ylab  total area [1.E6 km2]\'
'draw title \Global grass cover\'
'cbar_line  -p -c 3 2 7 -l 1 1 1  -t "C3" "C4" "C3-Artic" '

'printim grass_area_line.gif gif x1024 y768 white'

'close 1'
