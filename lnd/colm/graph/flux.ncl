
load "$NCARG_ROOT/lib/ncarg/nclex/gsun/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

procedure rdflux (case:string, date:logical, fsena:float, lfevpa:float, rnet:float, tg:float,\
                  lai:float, sai:float, albvdir:float, albvdif:float, albndir:float, albndif:float,\
                  albdir:float, albdif:float, mois:float, wliq:float, wice:float, rsur:float, rnof:float,\
                  tm:float, qm:float, prc:float, prl:float, solar:float, frl:float, us:float, vs:float, pbot:float)

local fname, year, month, day, second, deltim, ntime, itime, k

begin

   year   = date@year
   month  = date@month
   day    = date@day
   second = date@second
   deltim = date@deltim
   ntime  = date@ntime

   itime = 0

   do while (itime .lt. ntime)

      if (date@monthly) then
         print(sprinti("%0.4d ",year) + sprinti("%0.2d",month))

         fname = "/disk2/jidy/pclm3-lcc/output/" + case + "/GLOBAL-" + \
                 sprinti("%0.4d",year) + "-" + sprinti("%0.2d",month) + "-86400"
      end if

      if (date@daily) then
         print(sprinti("%0.4d ",year) + sprinti("%0.3d",day))

         fname = "/disk2/jidy/pclm3-lcc/output/" + case + "/GLOBAL-" + \
                 sprinti("%0.4d",year) + "-" + sprinti("%0.3d",day) + "-86400"
      end if

      if (date@hourly) then
         print(sprinti("%0.4d ",year) + sprinti("%0.3d ",day) + sprinti("%0.5d",second))

         fname = "/disk2/jidy/pclm3-lcc/output/" + case + "/GLOBAL-" + \
                 sprinti("%0.4d",year) + "-" + sprinti("%0.3d",day) + "-" + sprinti("%0.5d",second)
      end if

      fsena(itime)  = fbinrecread(fname, 2,1,"float")
      lfevpa(itime) = fbinrecread(fname, 3,1,"float")
      rnet(itime)   = fbinrecread(fname,15,1,"float")
      rsur(itime)   = fbinrecread(fname,18,1,"float")
      rnof(itime)   = fbinrecread(fname,19,1,"float")

      do k = 0, 9
         wliq(itime,k) = fbinrecread(fname,32+k,1,"float")
      end do

      do k = 0, 9
         wice(itime,k) = fbinrecread(fname,42+k,1,"float")
      end do

      mois(itime) = sum(wliq(itime,0:9))+sum(wice(itime,0:9))

      tg(itime)  = fbinrecread(fname,52,1,"float")
      lai(itime) = fbinrecread(fname,61,1,"float")
      sai(itime) = fbinrecread(fname,62,1,"float")

      albvdir(itime) = fbinrecread(fname,63,1,"float")
      albvdif(itime) = fbinrecread(fname,65,1,"float")
      albndir(itime) = fbinrecread(fname,64,1,"float")
      albndif(itime) = fbinrecread(fname,66,1,"float")

      albdir(itime) = albvdir(itime)+albndir(itime)
      albdif(itime) = albvdif(itime)+albndif(itime)

      us(itime)    = fbinrecread(fname,83,1,"float")
      vs(itime)    = fbinrecread(fname,84,1,"float")
      tm(itime)    = fbinrecread(fname,85,1,"float")
      qm(itime)    = fbinrecread(fname,86,1,"float")
      prc(itime)   = fbinrecread(fname,87,1,"float")
      prl(itime)   = fbinrecread(fname,88,1,"float")
      pbot(itime)  = fbinrecread(fname,89,1,"float")
      frl(itime)   = fbinrecread(fname,90,1,"float")
      solar(itime) = fbinrecread(fname,91,1,"float")

      itime = itime+1

      if (date@monthly) then
         month = month+1

         if (month .gt. 12) then
            month = 1
            year = year+1
         end if
      end if

      if (date@daily) then
         day = day+1

         if(day .gt. 365) then   ;not complete !!!
            day = 1
            year = year+1
         end if
      end if

      if (date@hourly) then
         second = second+deltim

         if(second .gt. 86400) then  ;not complete !!!
            second = second-86400
            day = day+1
         end if
      end if

   end do

end  ; end of procedure

procedure plotproc(wks:graphic, xi:float, data:float, title:string) 

local xyres, gsres, txres, pxy, xx, yy

begin

   xyres = True
   gsres = True
   txres = True

   xyres@gsnFrame               = False     ; Don't advance the frame.

   xyres@tmXTBorderOn           = False
   xyres@tmXTOn                 = False
   xyres@tmYRBorderOn           = False     ; Don't draw right axis.
   xyres@tmYROn                 = False     ; Don't draw right axis tick marks.
   xyres@tmBorderThicknessF     = 1.0

   xyres@tmXBLabelFont          = 21        ; Change font and size of 
   xyres@tmXBLabelFontHeightF   = 0.015     ; X axis labels.
   xyres@tmXBMajorLengthF       = 0.015     ; Default is 0.02.
   xyres@tmXBMajorThicknessF    = 1.0       ; Default is 2.
   xyres@tmXBMinorLengthF       = 0.0075    ; Default is 0.01.
   xyres@tmXBMinorPerMajor      = 4         ; # of minor tick marks per major.
   xyres@tmXBMode               = "Manual"  ; Set tick mark mode.
   xyres@tmXBTickStartF         = data@xmin
   xyres@tmXBTickEndF           = data@xmax
   xyres@tmXBTickSpacingF       = data@xstep
   xyres@tmXMajorGridThicknessF = 1.0       ; Default is 2.0

   xyres@tmYLLabelFont          = 21        ; See explanations for X axes
   xyres@tmYLLabelFontHeightF   = 0.015     ; resources.
   xyres@tmYLMajorLengthF       = 0.015
   xyres@tmYLMajorThicknessF    = 1.0
   xyres@tmYLMinorPerMajor      = 0
   xyres@tmYLMode               = "Manual"

   xyres@trXMinF                = data@xmin ; Specify data limits for X and Y axes.
   xyres@trXMaxF                = data@xmax

   xyres@vpHeightF              = 0.25      ; Define height, width, and location of plot.
   xyres@vpWidthF               = 0.80
   xyres@vpXF                   = 0.13
 
   xyres@xyMarkLineMode         = "Lines"   ; Set line mode to "Lines".

   txres@txFont                 = 21        ; Change the default font.
   txres@txFontHeightF          = 0.015     ; Set the font height.

   gsres@gsLineColor            = "red"     ; Set polyline color to red.

   xx = (/data@xmin,data@xmax/)             ; Create data for a polyline for marking
   yy = (/0.0,0.0/)                         ; the Y = 0.0 line in each graph.

 ;----------- first plot, same frame --------------------------------
 
   xyres@vpYF             = 0.90

   xyres@trYMinF          = data@ymin1
   xyres@trYMaxF          = data@ymax1

   xyres@tmYLTickStartF   = data@ymin1
   xyres@tmYLTickEndF     = data@ymax1
   xyres@tmYLTickSpacingF = data@ystep1
 
   pxy = gsn_xy(wks,xi,data(0,:),xyres)     ; Plot the interpolated values.

   tx = data@xmax/2.0
   ty = (data@ymax1-data@ymin1)*0.9+data@ymin1
 
   gsn_text(wks,pxy,data@vname1,tx,ty,txres) ; Label the plot.
 
   gsn_polyline(wks,pxy,xx,yy,gsres)        ; Draw polyline at Y=0.
 
 ;----------- second plot, same frame --------------------------------
 
   xyres@vpYF             =  0.60           ; Set Y location of plot.
 
   xyres@trYMinF          =  data@ymin2     ; Set minimum Y axis value.
   xyres@trYMaxF          =  data@ymax2     ; Set maximum Y axis value.
 
   xyres@tmYLTickStartF   =  data@ymin2     ; Define tick mark spacing along Y axis
   xyres@tmYLTickEndF     =  data@ymax2
   xyres@tmYLTickSpacingF =  data@ystep2
 
   pxy = gsn_xy(wks,xi,data(1,:),xyres)     ; Plot the integrals.
 
   tx = data@xmax/2.0
   ty = (data@ymax2-data@ymin2)*0.9+data@ymin2

   gsn_text(wks,pxy,data@vname2,tx,ty,txres) ; Label plot.
 
   gsn_polyline(wks,pxy,xx,yy,gsres)        ; Draw polyline at Y=0.

 ;----------- third plot, same frame --------------------------------

   xyres@vpYF             =  0.30           ; Set Y location of plot.

   xyres@trYMinF          =  data@ymin3
   xyres@trYMaxF          =  data@ymax3

   xyres@tmYLTickStartF   =  data@ymin3
   xyres@tmYLTickEndF     =  data@ymax3
   xyres@tmYLTickSpacingF =  data@ystep3

   pxy = gsn_xy(wks,xi,data(2,:),xyres)

   tx = data@xmax/2.0
   ty = (data@ymax3-data@ymin3)*0.9+data@ymin3

   gsn_text(wks,pxy,data@vname3,tx,ty,txres)

   gsn_polyline(wks,pxy,xx,yy,gsres)

 ;----------- plot elements, same frame -----------------------------
 
   txres@txFontHeightF = 0.02  ; Change the font height.
 
   gsn_text_ndc(wks,title,.5,.95,txres) 
 
   frame(wks)        ; Advance the frame.

end  ; end of procedure


begin

   site = "site2"
   lccs = "usgs"
   delt = "hourly"
   case = "nodynveg-baresoil-" + delt
   satt =  site + "_" + delt

  ;ntime   = 600
   ntime   = 1440
   nl_soil = 10

   date = True

   date@monthly = False
   date@daily   = False
   date@hourly  = True

  ;date@year   = 1951
   date@year   = 1998
   date@month  = 1
   date@day    = 182
   date@second = 1800
   date@deltim = 1800
   date@ntime  = ntime

   time = new((/ntime/),float)

   fsena   = new((/ntime/),float)
   lfevpa  = new((/ntime/),float)
   rnet    = new((/ntime/),float)
   rsur    = new((/ntime/),float)
   rnof    = new((/ntime/),float)
   tg      = new((/ntime/),float)
   lai     = new((/ntime/),float)
   sai     = new((/ntime/),float)
   albvdir = new((/ntime/),float)
   albvdif = new((/ntime/),float)
   albndir = new((/ntime/),float)
   albndif = new((/ntime/),float)
   albdir  = new((/ntime/),float)
   albdif  = new((/ntime/),float)
   mois    = new((/ntime/),float)
   wliq    = new((/ntime,nl_soil/),float)
   wice    = new((/ntime,nl_soil/),float)
   tm      = new((/ntime/),float)
   qm      = new((/ntime/),float)
   frl     = new((/ntime/),float)
   solar   = new((/ntime/),float)
   pbot    = new((/ntime/),float)
   us      = new((/ntime/),float)
   vs      = new((/ntime/),float)
   prc     = new((/ntime/),float)
   prl     = new((/ntime/),float)

   data    = new((/3,ntime/),float)


   wks_type = "pdf"
   wks_type@wkOrientation = "landscape"
  ;wks_type@wkOrientation = "portrait"

   wks = gsn_open_wks(wks_type,"flux-"+site+"-"+case)


   rdflux(site+"-"+lccs+"-"+case,\
          date,fsena,lfevpa,rnet,tg,lai,sai,albvdir,albvdif,albndir,albndif,albdir,\
          albdif,mois,wliq,wice,rsur,rnof,tm,qm,prc,prl,solar,frl,us,vs,pbot)

   do k = 1, ntime
      time(k-1) = k
   end do

   tg_min = min(tg)
   tg_max = max(tg)
   print ("tg:"+sprintf("%f",tg_min)+","+sprintf("%f",tg_max))

   mois_min = min(mois)
   mois_max = max(mois)
   print ("mois:"+sprintf("%f",mois_min)+","+sprintf("%f",mois_max))

   rsur = rsur*3600  ; unit from mm/s to mm/h

   rsur_min = min(rsur)
   rsur_max = max(rsur)
   print ("rsur:"+sprintf("%f",rsur_min)+","+sprintf("%f",rsur_max))

   rnof = rnof*3600  ; unit from mm/s to mm/h

   rnof_min = min(rnof)
   rnof_max = max(rnof)
   print ("rnof:"+sprintf("%f",rnof_min)+","+sprintf("%f",rnof_max))

   fsena_min = min(fsena)
   fsena_max = max(fsena)
   print ("fsena:"+sprintf("%f",fsena_min)+","+sprintf("%f",fsena_max))

   lfevpa_min = min(lfevpa)
   lfevpa_max = max(lfevpa)
   print ("lfevpa:"+sprintf("%f",lfevpa_min)+","+sprintf("%f",lfevpa_max))

   albdir_min = min(albdir)
   albdir_max = max(albdir)
   print ("albdir:"+sprintf("%f",albdir_min)+","+sprintf("%f",albdir_max))

   albdif_min = min(albdif)
   albdif_max = max(albdif)
   print ("albdif:"+sprintf("%f",albdif_min)+","+sprintf("%f",albdif_max))

   rnet_min = min(rnet)
   rnet_max = max(rnet)
   print ("rnet:"+sprintf("%f",rnet_min)+","+sprintf("%f",rnet_max))

   lai_min = min(lai)
   lai_max = max(lai)
   print ("lai:"+sprintf("%f",lai_min)+","+sprintf("%f",lai_max))

   sai_min = min(sai)
   sai_max = max(sai)
   print ("sai:"+sprintf("%f",sai_min)+","+sprintf("%f",sai_max))

   tm_min = min(tm)
   tm_max = max(tm)
   print ("tm:"+sprintf("%f",tm_min)+","+sprintf("%f",tm_max))

   qm_min = min(qm)
   qm_max = max(qm)
   print ("qm:"+sprintf("%f",qm_min)+","+sprintf("%f",qm_max))

   solar_min = min(solar)
   solar_max = max(solar)
   print ("solar:"+sprintf("%f",solar_min)+","+sprintf("%f",solar_max))

   frl_min = min(frl)
   frl_max = max(frl)
   print ("frl:"+sprintf("%f",frl_min)+","+sprintf("%f",frl_max))

   prc = prc*3600  ; unit from kg/m2/s to mm/h
   prl = prl*3600 

   prcp_min = min(prc+prl)
   prcp_max = max(prc+prl)
   print ("prcp:"+sprintf("%f",prcp_min)+","+sprintf("%f",prcp_max))

   us_min = min(us)
   us_max = max(us)
   print ("us:"+sprintf("%f",us_min)+","+sprintf("%f",us_max))

   vs_min = min(vs)
   vs_max = max(vs)
   print ("vs:"+sprintf("%f",vs_min)+","+sprintf("%f",vs_max))

   pbot_min = min(pbot/100.)
   pbot_max = max(pbot/100.)
   print ("pbot:"+sprintf("%f",pbot_min)+","+sprintf("%f",pbot_max))

   data@xmin  =  0.
   data@xmax  =  int2flt(ntime-1)
   data@xstep =  100.

   data(0,:) = tg
   data(1,:) = tm
   data(2,:) = fsena

   data@vname1 = "Soil surface temperature [K]"
   data@vname2 = "Air temperature [K]"
   data@vname3 = "Sensible heat flux [W/m2]"

   data@ymin1  =  280.
   data@ymax1  =  330.
   data@ystep1 =   10.
   data@ymin2  =  280.
   data@ymax2  =  330.
   data@ystep2 =   10.
   data@ymin3  = -250.
   data@ymax3  =  250.
   data@ystep3 =   50.

   plotproc(wks, time, data, "CLM Fluxes") 

   data(0,:) = lfevpa
   data(1,:) = mois
   data(2,:) = prc+prl

   data@vname1 = "Latent heat flux [W/m2]"
   data@vname2 = "Soil moisture [kg/m2]"
   data@vname3 = "Precipitation [mm/h]"

   data@ymin1  = -200.
   data@ymax1  =  600.
   data@ystep1 =  100.
   data@ymin2  =  800.
   data@ymax2  = 1200.
   data@ystep2 =  200.
   data@ymin3  =    0.
   data@ymax3  =   10.
   data@ystep3 =    2.

   plotproc(wks, time, data, "CLM Fluxes") 

   data(0,:) = rsur
   data(1,:) = rnof
   data(2,:) = qm

   data@vname1 = "Surface runoff [mm/h]"
   data@vname2 = "Total runoff [mm/h]"
   data@vname3 = "Specific humidity [kg/kg]"

   data@ymin1  =    0.
   data@ymax1  =    0.2
   data@ystep1 =    0.02
   data@ymin2  =    0.
   data@ymax2  =    0.2
   data@ystep2 =    0.02
   data@ymin3  =    0.
   data@ymax3  =    0.05
   data@ystep3 =    0.01

   plotproc(wks, time, data, "CLM Fluxes") 

   data(0,:) = albdif
   data(1,:) = albdir
   data(2,:) = solar

   data@vname1 = "Diffuse albedo [-]"
   data@vname2 = "Direct albedo [-]"
   data@vname3 = "Solar radiation [W/m2]"

   data@ymin1  =    0.
   data@ymax1  =    1.
   data@ystep1 =    0.2
   data@ymin2  =    0.
   data@ymax2  =    1.
   data@ystep2 =    0.2
   data@ymin3  =    0.
   data@ymax3  = 1000.
   data@ystep3 =  200.

   plotproc(wks, time, data, "CLM Fluxes") 

   data(0,:) = lai
   data(1,:) = sai
   data(2,:) = frl

   data@vname1 = "Leaf area index [-]"
   data@vname2 = "Stem area index [-]"
   data@vname3 = "Longwave radiation [W/m2]"

   data@ymin1  =    0.
   data@ymax1  =    4.
   data@ystep1 =    1.
   data@ymin2  =    0.
   data@ymax2  =    4.
   data@ystep2 =    1.
   data@ymin3  =  250.
   data@ymax3  =  500.
   data@ystep3 =   50.

   plotproc(wks, time, data, "CLM Fluxes") 

   data(0,:) = us
   data(1,:) = vs
   data(2,:) = pbot/100.

   data@vname1 = "Wind u [m/s]"
   data@vname2 = "Wind v [m/s]"
   data@vname3 = "Near surface air pressure [hPa]"

   data@ymin1  =  -10.
   data@ymax1  =   10.
   data@ystep1 =    5.
   data@ymin2  =  -10.
   data@ymax2  =   10.
   data@ystep2 =    5.
   data@ymin3  =  800.
   data@ymax3  = 1100.
   data@ystep3 =  100.

   plotproc(wks, time, data, "CLM Fluxes") 

end
