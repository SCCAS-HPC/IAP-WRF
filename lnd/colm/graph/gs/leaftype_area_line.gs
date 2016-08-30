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
'define broadleaf=asum((bdttro+bdttem+bdtbor+bettem+bettro)*soil*garea,lon=-179,lon=179,lat=-59,lat=89)'
'd tloop(broadleaf*1.E-6)'
'undefine broadleaf'

'set cmark  0'
'set ccolor 4'
'define needleleaf=asum((nettem+netbor)*soil*garea,lon=-179,lon=179,lat=-59,lat=89)'
'd tloop(needleleaf*1.E-6)'
'undefine needleleaf'

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
'draw title \Global leaf type cover\'
'cbar_line  -p -c 3 4 7 2 -l 1 1 1 1 -t "broadleaf" "needleleaf" "grass" "bare soil" '

'printim leaftype_area_line.gif gif x1024 y768 white'

'c'

'set vrange 0 25'
'set ylint 5'

'set cmark  0'
'set ccolor 2'
'define blbdttro=asum(bdttro*soil*garea,lon=-179,lon=179,lat=-59,lat=89)'
'd tloop(blbdttro*1.E-6)'
'undefine blbdttro'

'set cmark  0'
'set ccolor 3'
'define blbdttem=asum(bdttem*soil*garea,lon=-179,lon=179,lat=-59,lat=89)'
'd tloop(blbdttem*1.E-6)'
'undefine blbdttem'

'set cmark  0'
'set ccolor 4'
'define blbdtbor=asum(bdtbor*soil*garea,lon=-179,lon=179,lat=-59,lat=89)'
'd tloop(blbdtbor*1.E-6)'
'undefine blbdtbor'

'set cmark  0'
'set ccolor 5'
'define blbettem=asum(bettem*soil*garea,lon=-179,lon=179,lat=-59,lat=89)'
'd tloop(blbettem*1.E-6)'
'undefine blbettem'

'set cmark  0'
'set ccolor 7'
'define blbettro=asum(bettro*soil*garea,lon=-179,lon=179,lat=-59,lat=89)'
'd tloop(blbettro*1.E-6)'
'undefine blbettro'

'draw ylab  total area [1.E6 km2]\'
'draw title \Global broadleaf trees cover\'
'cbar_line  -p -c 2 3 4 5 7 -l 1 1 1 1 1 -t "BDTtro" "BDTtem" "BDTbor" "BETtem" "BETtro" '

'printim broadleaf_area_line.gif gif x1024 y768 white'

'c'

'set vrange 0 15'
'set ylint 3'

'set cmark  0'
'set ccolor 2'
'define nlnettem=asum(nettem*soil*garea,lon=-179,lon=179,lat=-59,lat=89)'
'd tloop(nlnettem*1.E-6)'
'undefine nlnettem'

'set cmark  0'
'set ccolor 4'
'define nlnetbor=asum(netbor*soil*garea,lon=-179,lon=179,lat=-59,lat=89)'
'd tloop(nlnetbor*1.E-6)'
'undefine nlnetbor'

'draw ylab  total area [1.E6 km2]\'
'draw title \Global needleleaf cover\'
'cbar_line  -p -c 2 4 -l 1 1 -t "NETtem" "NETbor" '

'printim needleleaf_line.gif gif x1024 y768 white'

'close 1'
