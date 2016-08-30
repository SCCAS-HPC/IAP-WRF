DSET    /p2/zhq/CoLM_PRC/site15/OUT/site15-%y4
DSET    /p2/zhq/CoLM_PRC/site15/OUT1/site15-%y4-%m2-%d2
DSET    /p2/zhq/CoLM_PRC/site15/OUT/site15-%y4-%m2-%d2
DSET    /home/jidy/swgfs/colm_esm2/global-%y4-%m2
DSET    /home/jidy/swgfs/colm_esm/global-%y4-%m2

TITLE    CLM SIMULATED SURFACE FLUXES
UNDEF    -9999.
OPTIONS  SEQUENTIAL template big_endian

XDEF     36 linear -175 10 
YDEF     18 linear  -85 10 
ZDEF      10 levels
0.0071006 0.0279250 0.062258 0.118865 0.212193                       
0.3660658 0.6197585 1.038027 1.727635 2.864607
TDEF      99999  linear 00:30Z01JAN1490 1mo
VARS      65
taux      0  99  01: wind stress: E-W [kg/m/s2] 
tauy      0  99  02: wind stress: N-S [kg/m/s2]
fsena     0  99  03: sensible heat from canopy height to atmosphere [W/m2]
lfevpa    0  99  04: latent heat flux from canopy height to atmosphere [W/m2]
fevpa     0  99  05: evapotranspiration from canopy height to atmosphere [mm/s]
fsenl     0  99  06: sensible heat from leaves [W/m2]
fevpl     0  99  07: evaporation+transpiration from leaves [mm/s]
etr       0  99  08: transpiration rate [mm/s]
fseng     0  99  09: sensible heat flux from ground [W/m2]
fevpg     0  99  10: evaporation heat flux from ground [mm/s]
fgrnd     0  99  11: ground heat flux [W/m2]
sabvsun   0  99  12: solar absorbed by sunlit canopy [W/m2]
sabvsha   0  99  13: solar absorbed by shaded [W/m2]
sabg      0  99  14: solar absorbed by ground  [W/m2]
olrg      0  99  15: outgoing long-wave radiation from ground+canopy [W/m2]
rnet      0  99  16: net radiation [W/m2]
zerr      0  99  17: the error of energy balance [W/m2]
assim     0  99  18: canopy assimilation rate [mol m-2 s-1]
respc     0  99  19: respiration (plant+soil) [mol m-2 s-1]
tlsun     0  99  20: sunlit leaf temperature [K]
tlsha     0  99  21: shaded leaf temperature [K]
ldew      0  99  22: depth of water on foliage [mm]
sigf      0  99  23: fraction of veg cover, excluding snow-covered veg [-]
green     0  99  24: leaf greenness
lai       0  99  25: leaf area index
sai       0  99  26: stem area index
albvdir   0  99  27: averaged albedo [vis, dir]
albvdif   0  99  28: averaged albedo [vis, dif]
albndir   0  99  29: averaged albedo [nir, dir]
albndif   0  99  30: averaged albedo [nir, dif]
emis      0  99  31: averaged bulk surface emissivity
z0ma      0  99  32: effective roughness [m]
trad      0  99  33: radiative temperature of surface [K]
ustar     0  99  34: u* in similarity theory [m/s]
tstar     0  99  35: t* in similarity theory [kg/kg]
qstar     0  99  36: q* in similarity theory [kg/kg]
zol       0  99  37: dimensionless height (z/L) used in Monin-Obukhov theory
rib       0  99  38: bulk Richardson number in surface layer
fm        0  99  39: integral of profile function for momentum
fh        0  99  40: integral of profile function for heat
fq        0  99  41: integral of profile function for moisture
tref      0  99  42: 2 m height air temperature [kelvin]
qref      0  99  43: 2 m height air specific humidity [kg/kg]
u10m      0  99  44: 10m u-velocity [m/s]
v10m      0  99  45: 10m v-velocity [m/s]
f10m      0  99  46: integral of profile function for momentum at 10m [-]
xerr      0  99  47: the error of water banace [mm/s]
rsur      0  99  48: surface runoff [mm/s]
rnof      0  99  49: total runoff [mm/s]
tss      10  99  50: soil temperature [K]
wliq     10  99  60: liquid water in soil layers [kg/m2]
wice     10  99  70: ice lens in soil layers [kg/m2]
tg        0  99  80: ground surface temperature [K]
scv       0  99  81: snow cover, water equivalent [mm]
snowdp    0  99  82: snow depth [meter]
fsno      0  99  83: fraction of snow cover on ground
us        0  99  84: wind in eastward direction [m/s]
vs        0  99  85: wind in northward direction [m/s]
tm        0  99  86: temperature at reference height [kelvin]
qm        0  99  87: specific humidity at reference height [kg/kg]
prc       0  99  88: convective precipitation [mm/s]
prl       0  99  89: large scale precipitation [mm/s]
pbot      0  99  90: atmospheric pressure at the surface [pa]
frl       0  99  91: atmospheric infrared (longwave) radiation [W/m2]
solar     0  99  92: downward solar radiation at surface [W/m2]
ENDVARS
