'open ../ctl/CN.reg.ctl'
*'open co2-355.ctl'
*'open co2-375.ctl'

'set grid off'
'set grads off'
*'set parea 1 9 0.3 8.3'

*imagepath='global/'
'set t 1'
'set x 1 180'
'set y 1 75'
'set z 1'
'define landc=bare'

t=100

x=100
while(x<=120)
   'set x 'x
y=45
while(y<=55)
   'set y 'y

   'define tree=(nettem+netbor+bettro+bettem+bdttro+bdttem+bdtbor)'
   'define grass=(c3+c4+c3arc)'
   'define soil=bare'
   say 'tree= 'tree''

   if(tree>grass & tree>soil)
   landc=1.
   endif
   if(grass>tree & grass>soil)
   landc=2.
   endif
   if(soil>tree & soil>grass)
   landc=0.
   endif
*   say landc
   y=y+1 
    'undefine tree'
    'undefine grass'
    'undefine soil'
endwhile

   x=x+1
endwhile

*   'set gxout fgrid'
*   'set fgvals 0. 9  1. 4  2. 2'
   say landc
   'set gxout shaded'
   'draw title Global Land Cover spinup year: '%t%''
   'd landc'
   'cbarn'
*   'printim landcover'%t%'.gif gif x1024 y768 white'
   'undefine landc'

'close 1'
