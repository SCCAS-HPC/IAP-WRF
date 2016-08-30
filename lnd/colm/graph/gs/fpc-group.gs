'open ../ctl/CN.reg.ctl'
*'open co2-355.ctl'
*'open co2-375.ctl'

'set grid off'
'set grads off'
*'set parea 1 9 0.3 8.3'

imagepath='global500/'

*'set t 465 512'
*'set t 1425 1472'

t=500
*t=1425

while(t<=500)

   'set t 't

   'set gxout shaded'
   'set clevs 0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0'
   'set ccols 0 4 11 5 13 3 10 7 12 8 2 6'
*   'set cint 2'
   'd (nettem+netbor+bettro+bettem+bdttro+bdttem+bdtbor)'
   'draw title Tree Cover \ spinup year: '%t%''
   'cbarn'
   'printim '%imagepath%'tree'%t%'.gif gif x1024 y768 white'
   'c'

   'set gxout shaded'
*   'set cint 2'
   'set clevs 0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0'
   'set ccols 0 4 11 5 13 3 10 7 12 8 2 6'
   'd (c3+c4+c3arc)'
   'draw title Grass Cover \ spinup year: '%t%''
   'cbarn'
   'printim '%imagepath%'grass'%t%'.gif gif x1024 y768 white'
   'c'

   'set gxout shaded'
*   'set cint 2'
   'set clevs 0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0'
   'set ccols 0 4 11 5 13 3 10 7 12 8 2 6'
   'd bare'
   'draw title Bare Cover \ spinup year: '%t%''
   'cbarn'
   'printim '%imagepath%'bare'%t%'.gif gif x1024 y768 white'
   'c'

    t=t+100; say t
endwhile

'close 1'
