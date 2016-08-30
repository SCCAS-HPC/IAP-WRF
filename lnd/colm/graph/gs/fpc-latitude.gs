'open ../ctl/CN.reg.ctl'
*'open co2-355.ctl'
*'open co2-375.ctl'

'set grid off'
'set grads off'
*'set parea 1 9 0.3 8.3'

imagepath='global500/'

'set y 1 75'
'set x 1 '
'set t 500'

t=500

   'set gxout line'
   'set vrange 0 1'
   'set ylint 0.2'
   'set cthick 10'
'set cmark  0'
'set ccolor 4'
   'd ave(nettem+netbor+bettro+bettem+bdttro+bdttem+bdtbor,x=1,x=180)'
'set cmark  0'
'set ccolor 3'
   'd ave(c3+c4+c3arc,x=1,x=180)'
'set cmark  0'
'set ccolor 8'
   'd ave(bare,x=1,x=180)'
   'draw title FPC latitude distribution  spinup year: '%t%''
   'cbar_line  -p -c 4 3 8 -l  1 1 1 -t "tree" "grass" "bare soil" '

   'printim '%imagepath%'FPC-latitude'%t%'.gif gif x1280 y768 white'

'close 1'
