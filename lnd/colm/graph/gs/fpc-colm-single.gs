'open ../ctl/CN.reg.ctl'
'set t 1 100'
'set grads off'
'set grid off'
'set vrange 0 1'
*'set xlab off'

*1. plot Tropical forest (0.25S,69.75W)
'set x 56'
'set y 30'
'set cthick 6'
'set ccolor 2'
*'set cstyle 1'
'set cmark  0'
'd bettro'

'set ccolor 4'
*'set cstyle 2'
'set cmark  0'
'd bdttro'

'set ccolor 3'
*'set cstyle 3'
'set cmark  0'
'd C4'

'draw title Tropical forest (1.0S,69.0W)'
'draw xlab model year'
'draw ylab Grid FPC'
'run cbar_line -m 0 0 0  -c 2 4 3  -l 1 1 1  -t "BETtro" "BDTtro" "C4" -p' 
*'printim fpcgrid-tropical-100.gif gif x1280 y768 white'
'c'

*2. plot Temperate mixed forest (51.0N,11.0E)
'set vrange 0 1'
'set x 96 '
'set y 56'
'set cthick 6'
'set ccolor 2'
*'set cstyle 1'
'set cmark  0'
'd bdttem'

'set ccolor 4'
*'set cstyle 2'
'set cmark  0'
'd bdtbor'

'set ccolor 3'
*'set cstyle 3'
'set cmark  0'
'd netbor'

'set ccolor 8'
*'set cstyle 3'
'set cmark  0'
'd C3'

'draw title Temperate mixed forest (51.0N,11.0E)'
'draw xlab model year'
'draw ylab Grid FPC'
'run cbar_line -m 0 0 0  -c 2 4 3 8  -l 1 1 1 1  -t "BDTtem" "BDTbor" "NETbor" "C3" -p'
*'printim fpcgrid-temperatemix-100.gif gif x1280 y768 white'
'c'

*3. plot Boreal forest (62.25N,15.75E)
'set vrange 0 1'
'set x 98'
'set y 62'
'set cthick 6'
'set ccolor 2'
*'set cstyle 1'
'set cmark  0'
'd bdtbor'

'set ccolor 4'
*'set cstyle 2'
'set cmark  0'
'd netbor'

'set ccolor 8'
*'set cstyle 3'
'set cmark  0'
'd C3'

'draw title Boreal forest (63.0N,15.0E)'
'draw xlab model year'
'draw ylab Grid FPC'
'run cbar_line -m 0 0 0  -c 2 4 8  -l 1 1 1  -t "BDTbor" "NETbor" "C3" -p'
*'printim fpcgrid-boreal-100.gif gif x1280 y768 white'
'c'

* plot Savanna (14.75S,20.25E)
'set vrange 0 1'
'set x 101'
'set y 23'
'set cthick 6'
'set ccolor 2'
*'set cstyle 1'
'set cmark  0'
'd bettro'

'set ccolor 4'
*'set cstyle 2'
'set cmark  0'
'd bdttro'

'set ccolor 3'
*'set cstyle 2'
'set cmark  0'
'd bdttem'

'set ccolor 5'
*'set cstyle 2'
'set cmark  0'
*'d nettem'
'd C3'

'set ccolor 8'
*'set cstyle 3'
'set cmark  0'
'd C4'

'draw title Savanna (15.0S,21.0E)'
'set vrange 0 1'
'draw xlab model year'
'draw ylab Grid FPC'
'run cbar_line -m 0 0 0 0 0  -c 2 4 3 5 8  -l 1 1 1 1 1  -t "BETtro" "BDTtro" "bdttem" "C3" "C4" -p'
'printim fpcgrid-savanna-100.gif gif x1280 y768 white'
'c'

* plot Glassland (32.25S,120.25E)
'set vrange 0 1'
'set x 151 '
'set y 14'
'set cthick 6'
'set ccolor 2'
*'set cstyle 1'
'set cmark  0'
'd nettem'

'set ccolor 4'
*'set cstyle 2'
'set cmark  0'
'd bdttem'

'set ccolor 4'
*'set cstyle 2'
'set cmark  0'
'd bettem'

'set ccolor 8'
*'set cstyle 3'
'set cmark  0'
'd C3'

'draw title Glassland (33.0S,121.0E)'
'draw xlab model year'
'draw ylab Grid FPC'
'run cbar_line -m 0 0 0  -c 2 4 7  -l 1 1 8  -t "BETtro" "BDTtro" "C4" -p'
*'printim fpcgrid-grassland-100.gif gif x1280 y768 white'
'c'
