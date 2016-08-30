'open ../ctl/CN.reg.ctl'

'set grid off'
'set grads off'

*'set xlint 100'
'set ylint 20'
'set vrange -20 100'

'set x 1'
'set y 1'
'set z 1'
'set t 1 500'

'set cmark  0'
'set ccolor 2'
'define anpp=asum((anpp)*soil*garea,lon=-179,lon=179,lat=-59,lat=89)'
'd tloop(anpp*1.E-9)'
'undefine anpp'

'set cmark  0'
'set ccolor 3'
'define amrh=asum((amrh)*soil*garea,lon=-179,lon=179,lat=-89,lat=89)'
'd tloop(amrh*1.E-9)'
'undefine amrh'

'set cmark  0'
'set ccolor 7'
'define nep=asum((anpp-amrh)*soil*garea,lon=-179,lon=179,lat=-59,lat=89)'
'd tloop(nep*1.E-9)'
'undefine nep'

'draw ylab  carbon flux [GT/year]\'
'draw title \Global Carbon Production\'
'cbar_line  -p -c 2 3 7 -l  1 1 1 -t "NPP" "Rh" "NEP" '

*'printim NPP_Rh_NEP_line.gif gif x1024 y768 white'

'close 1'
