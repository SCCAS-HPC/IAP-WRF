      parameter(nperox = 10)	! total number of alkylperoxy radicals
      parameter(nphoto = 20)	! total number of photolyzing species

      parameter(nrxn_het = 7,
     +          nrxn_com = 74,
     +          nrxn_urb = 44,
     +          nrxn_bio = 20,
     +          nrxn_mar = 35)

      parameter(nreg1 = ngas_com,
     +          nreg2 = ngas_com + ngas_urb,
     +          nreg3 = ngas_com + ngas_urb + ngas_bio,
     +          nreg4 = ngas_com + ngas_mar,
     +          nreg5 = ngas_com + ngas_urb + ngas_mar,
     +          nreg6 = ngas_com + ngas_urb + ngas_bio + ngas_mar)

c------------------------------------------------------------------------

      common /regime/ iregime

      common/gasrateconst/rk_com(nrxn_com),rk_urb(nrxn_urb),
     +                    rk_bio(nrxn_bio),rk_mar(nrxn_mar),
     +                    rk_het(nrxn_het),
     +                    rk_param(nperox),rk_photo(nphoto),
     +                    Aperox(nperox,nperox),Bperox(nperox,nperox),
     +                    uptake_coeff(nrxn_het), hetrxn(nrxn_het)
c
      common/gasratescom/ r_com(nrxn_com),r_urb(nrxn_urb),
     +                    r_bio(nrxn_bio),r_mar(nrxn_mar),
     +                    r_het(nrxn_het)
c
      common/gaspd/ p_com(ngas_max),d_com(ngas_max),
     +              p_urb(ngas_max),d_urb(ngas_max),
     +              p_bio(ngas_max),d_bio(ngas_max),
     +              p_mar(ngas_max),d_mar(ngas_max),
     +              p_het(ngas_max),d_het(ngas_max)

      common/het/ Npcasp(15)
c

c------------------------------------------------------------------------

      common /localGas/
     + ih2so4,      ihno3,       ihcl,        inh3,        ino,
     + ino2,        ino3,        in2o5,       ihono,       ihno4,
     + io3,         io1d,        io3p,        ioh,         iho2,
     + ih2o2,       ico,         iso2,        ich4,        ic2h6,           
     + ich3o2,      iethp,       ihcho,       ich3oh,      ianol,
     + ich3ooh,     iethooh,     iald2,       ihcooh,      ircooh,
     + ic2o3,       ipan,        ipar,        iaone,       imgly,
     + ieth,        iolet,       iolei,       itol,        ixyl,
     + icres,       ito2,        icro,        iopen,       ionit,
     + irooh,       iro2,        iano2,       inap,        ixo2,
     + ixpar,       iisop,       iisoprd,     iisopp,      iisopn,
     + iisopo2,     iterp,       isv1,        isv2,        isv3,
     + isv4,        isv5,        isv6,
     + idms,        imsa,        idmso,       idmso2,      ich3so2h,
     + ich3sch2oo,  ich3so2,     ich3so3,     ich3so2ch2oo,ich3so2oo,
     + isulfhox

      common /gasPeroxyRadicals/
     + jch3o2,      jethp,       jro2,        jc2o3,       jano2,
     + jnap,        jisopp,      jisopn,      jisopo2,     jxo2

      common /photolyzingSpecies/
     + jphoto_no2,    jphoto_no3,   jphoto_hono,   jphoto_hno3, 
     + jphoto_hno4,   jphoto_n2o5,  jphoto_o3a,    jphoto_o3b,    
     + jphoto_h2o2,   jphoto_hchoa, jphoto_hchob,  jphoto_ch3ooh,
     + jphoto_ethooh, jphoto_ald2,  jphoto_aone,   jphoto_mgly,
     + jphoto_open,   jphoto_rooh,  jphoto_onit,   jphoto_isoprd
 
