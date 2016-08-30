'open ../ctl/CN.amazon.ctl'
'open ../ctl/CN.cbs.ctl'
'open ../ctl/CN.dhs.ctl'
'open ../ctl/CN.hb.ctl'

'set t 1 2000'
'set grads off'
'set grid off'
'set vrange 0 1'
*'set xlab off'

'enable print 4plots.CN.gm'

*1. plot Tropical forest (0.5S,69.5W)
'set cthick 6'
'set ccolor 2'
*'set cstyle 1'
'set cmark  0'
'd bettro.1'

'set ccolor 4'
*'set cstyle 2'
'set cmark  0'
'd bdttro.1'

'set ccolor 3'
*'set cstyle 3'
'set cmark  0'
'd C4.1'

'draw title Tropical forest (1.0S,69.0W)'
'draw xlab model year'
'draw ylab Grid FPC'
'run cbar_line -m 0 0 0  -c 2 4 3  -l 1 1 1  -t "BETtro" "BDTtro" "C4" -p' 
*'printim fpcgrid-tropical-100.gif gif x1280 y768 white'
'print'
'c'

*2. plot Boreal mixed forest (41.5N,127.5E) CBS
'set vrange 0 1'
'set cthick 6'
'set ccolor 2'
*'set cstyle 1'
'set cmark  0'
'd bdtbor.2'

'set ccolor 3'
*'set cstyle 3'
'set cmark  0'
'd netbor.2'

'set ccolor 8'
*'set cstyle 3'
'set cmark  0'
'd c3arc.2'

'draw title Boreal mixed forest (41.5N,127.5E)'
'draw xlab model year'
'draw ylab Grid FPC'
'run cbar_line -m 0 0 0  -c 2 3 8  -l 1 1 1  -t "BDTbor" "NETbor" "C3arc" -p'
*'printim fpcgrid-temperatemix-100.gif gif x1280 y768 white'
'print'
'c'

*3. Temperate evergreen forest (23.5N,112.5E) DHS
'set vrange 0 1'
'set cthick 6'
'set ccolor 2'
*'set cstyle 1'
'set cmark  0'
'd bettem.3'

'set ccolor 4'
*'set cstyle 2'
'set cmark  0'
'd bdttem.3'

'set ccolor 3'
'set cmark  0'
'd nettem.3'

'set ccolor 8'
*'set cstyle 3'
'set cmark  0'
'd c3.3'

'draw title Temperate forest (23.5N,112.5E)'
'draw xlab model year'
'draw ylab Grid FPC'
'run cbar_line -m 0 0 0 0 -c 2 4 3 8  -l 1 1 1 1  -t "BETtem" "BDTtem" "NETtem" "C3" -p'
'print'
'c'

* plot Glassland (43.5N,116.5E) HB
'set vrange 0 1'
'set cthick 6'
'set ccolor 2'
*'set cstyle 1'
'set cmark  0'
'd bare.4'

'set ccolor 4'
*'set cstyle 2'
'set cmark  0'
'd c3.4'

'draw title Glassland (43.5N,116.5E)'
'draw xlab model year'
'draw ylab Grid FPC'
'run cbar_line -m 0 0 -c 2 4  -l 1 1 -t "Bare soil" "C3" -p'
*'printim fpcgrid-grassland-100.gif gif x1280 y768 white'
'print'
'c'

'close 4'
'close 3'
'close 2'
'close 1'
'disable print'
