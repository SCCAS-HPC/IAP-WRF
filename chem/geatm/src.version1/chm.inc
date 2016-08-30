      parameter(ngas_com = 32, 
     +          ngas_urb = 21, 
     +          ngas_bio = 10, 
     +          ngas_mar = 11)

      parameter(ngas_max = ngas_com + ngas_urb + ngas_bio + ngas_mar)

      parameter(nmax = ngas_max)

c------------------------------------------------------------------------
c Global Species Indices
c
      common /globalGas/
     + kh2so4,      khno3,       khcl,        knh3,        kno,
     + kno2,        kno3,        kn2o5,       khono,       khno4,
     + ko3,         ko1d,        ko3p,        koh,         kho2,
     + kh2o2,       kco,         kso2,        kch4,        kc2h6,           
     + kch3o2,      kethp,       khcho,       kch3oh,      kanol,
     + kch3ooh,     kethooh,     kald2,       khcooh,      krcooh,
     + kc2o3,       kpan,        kpar,        kaone,       kmgly,
     + keth,        kolet,       kolei,       ktol,        kxyl,
     + kcres,       kto2,        kcro,        kopen,       konit,
     + krooh,       kro2,        kano2,       knap,        kxo2,
     + kxpar,       kisop,       kisoprd,     kisopp,      kisopn,
     + kisopo2,     kterp,       ksv1,        ksv2,        ksv3,
     + ksv4,        ksv5,        ksv6,
     + kdms,        kmsa,        kdmso,       kdmso2,      kch3so2h,
     + kch3sch2oo,  kch3so2,     kch3so3,     kch3so2ch2oo,kch3so2oo,
     + ksulfhox

c-------------------------------------------------------------------------
      character*40 species(nmax)
      common/speciesNames/ species

      common /globalConcentration/ cnn(nmax)

      common /massbalance/ tNi, tSi, tCli, DN, DS, DCl

      common/emissions/ emission(nmax), emit(nmax)

      common /constants/ pi, avogad, deg2rad

      common /timenposition/
     + tbeg_dd,   tbeg_hh,   tbeg_mm,   tbeg_ss,
     + trun_dd,   trun_hh,   trun_mm,   trun_ss,
     + tbeg_sec,  trun_sec,  tcur_sec,  tcur_min,
     + tcur_hrs,  tmid_sec,  tsav_sec,  told_sec, 
     + dt_sec,    dt_min,    rlon,      rlat,
     + zalt_m,    cos_sza  

      common /timecounters/
     + ncstep,     it,        iprint

      common /bulkAirProperties/
     + cair_mlc,  h2o,       o2,      h2,       ppb,
     + te,        pr_atm,    RH,      speed_molec

      common /photolysis/
     + idaytime,  msolar, mphoto

c------------------------------------------------------------------------
