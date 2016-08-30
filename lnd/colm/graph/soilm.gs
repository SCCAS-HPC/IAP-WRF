'open flx.ctl'

'set t 1 204'

'set z 1'

'define tmp=sumg(wliq+wice,z=1,z=5)'

'set x 1'
'set y 1'

'define tmp2=aave(tmp,lon=-180,lon=180,lat=-90,lat=90)'

'set vrange 70 90'

'd tmp2'

'close 1'
