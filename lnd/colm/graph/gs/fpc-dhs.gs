'reinit';
'open ../ctl/fpc.ctl';
'set t 1 200';
'set grads off';
'set vrange 0 1';
* 'set xlab off';
'set z 1';
'set ccolor 1';
'set cstyle 1';
'set cmark  0';
'd p';
*'set ylpos 0 r';
'set z 5';
'set ccolor 2';
'set cstyle 1';
'set cmark  0';
'd p';
'set z 7';
'set ccolor 7';
'set cstyle 1';
'set cmark  0';
'd p';
'set z 13';
'set ccolor 3';
'set cstyle 1';
'set cmark  0';
'd p';
'set z 17';
'set ccolor 6';
'set cstyle 1';
'd p';
'run cbar_line -m 0 0 0 0 0 -c 1 2 7 3 6 -l 1 1 1 1 1 -t "NETtem" "BETtem" "BDTtem" "C3" "bare soil" -p'; 
'draw title DingHuShan spinup from bare soil';
'draw xlab model year';
'draw ylab fpcgrid';
'enable print fpcgrid-dhs-0424.png';
'print';
'disable print';
'disable fwrite';
;