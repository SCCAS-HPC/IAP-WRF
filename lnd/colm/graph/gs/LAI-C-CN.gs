'open flx.AMAZON.C.ctl'
'open flx.AMAZON.CN.ctl'
'open flx.DHS.C.ctl'
'open flx.DHS.CN.ctl'
'open flx.CBS.C.ctl'
'open flx.CBS.CN.ctl'
'open flx.HB.C.ctl'
'open flx.HB.CN.ctl'
'set t 12001 14400'
'set z 1'
'set vpage 1 10 4 8.5 '
'set grads off'
'set grid off'
'set cthick 6'
'set ccolor 2'
'set vrange 0 10 '
'set cmark  0'
'd lai.1'
'set ccolor 3'
'set cmark  0'
'd lai.2'
'cbar_line -p -c 2 3 -l 1 1 -t "C-DGVM" "CN-DGVM"'
'draw ylab LAI'
'draw title \tropical forest\'

'set vpage 1 10 0 4.5 '
'set grads off'
'set grid off'
'set cthick 6'
'set ccolor 2'
'set vrange 0 10 '
'set cmark  0'
'd lai.3'
'set ccolor 3'
'set cmark  0'
'd lai.4'
'cbar_line -p -c 2 3 -l 1 1 -t "C-DGVM" "CN-DGVM"'
'draw ylab LAI'
'draw title \temperate forest\'

'printim lai.gif gif x1024 y768 white'

'close 8'
'close 7'
'close 6'
'close 5'
'close 4'
'close 3'
'close 2'
'close 1'


