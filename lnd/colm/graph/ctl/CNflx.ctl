DSET      /p2/zhq/CoLM_PRC/site15/OUT/site15-CN-%y4
DSET      /swgfs/zhangq/CoLM_PRC/global/OUT/global-CN-%y4

TITLE    CLM SIMULATED SURFACE FLUXES
UNDEF    -9999.
OPTIONS  yrev SEQUENTIAL template

XDEF     1 linear 1 1 
YDEF     1 linear 1 1  
XDEF     180 linear -179 2 
YDEF      90 linear  -89 2  
ZDEF     14 linear 1 1 
TDEF      3000  linear 00:30Z01JAN0001 1yr
VARS      38
soil      0  99  01: wind stress: E-W [kg/m/s2] 
urban     0  99  02: wind stress: N-S [kg/m/s2]
wetld     0  99  03: sensible heat from canopy height to atmosphere [W/m2]
ice       0  99  04: latent heat flux from canopy height to atmosphere [W/m2]
lake      0  99  05: evapotranspiration from canopy height to atmosphere [mm/s]
fpc      14  99  06: sensible heat from leaves [W/m2]
bare      0  99  07: respiration (plant+soil) [mol m-2 s-1]
afirec    0  99  08: sunlit leaf temperature [K]
afiref    0  99  09: shaded leaf temperature [K]
avegc     0  99  10: depth of water on foliage [mm]
aestabc   0  99  11: fraction of veg cover, excluding snow-covered veg [-]
anpp      0  99  12: leaf greenness
amrh      0  99  13: leaf area index
alitcag   0  99  14: stem area index
alitcbg   0  99  15: averaged albedo [vis, dir]
asoicf    0  99  16: averaged albedo [vis, dir]
asoics    0  99  17: averaged albedo [vis, dir]
npp_ind  14  99  18: evaporation+transpiration from leaves [mm/s]
lm       14  99  19: evaporation+transpiration from leaves [mm/s]
sm       14  99  20: evaporation+transpiration from leaves [mm/s]
hm       14  99  21: evaporation+transpiration from leaves [mm/s]
rm       14  99  22: evaporation+transpiration from leaves [mm/s]
ca       14  99  23: evaporation+transpiration from leaves [mm/s]
htop     14  99  24: evaporation+transpiration from leaves [mm/s]
nind     14  99  25: evaporation+transpiration from leaves [mm/s]
laimx    14  99  26: evaporation+transpiration from leaves [mm/s]
cnleaf   14  99  27: evaporation+transpiration from leaves [mm/s]
cnsap    14  99  28: evaporation+transpiration from leaves [mm/s]
cnroot   14  99  29: evaporation+transpiration from leaves [mm/s]
an_up     0  99  30: averaged albedo [vis, dir]
stress    0  99  31: averaged albedo [vis, dir]
avegn     0  99  32: averaged albedo [vis, dir]
alitnag   0  99  33: averaged albedo [vis, dir]
alitnbg   0  99  34: averaged albedo [vis, dir]
asoin     0  99  35: averaged albedo [vis, dir]
no3       0  99  36: averaged albedo [vis, dir]
nh4       0  99  37: averaged albedo [vis, dir]
garea     0  99  38: grid area [km^2]
ENDVARS
