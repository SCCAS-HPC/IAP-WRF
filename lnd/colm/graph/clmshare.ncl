
procedure looktable(vname:string, p1:integer, p2:integer)

local flux_name, flux_p1, flux_p2, nvar, k

begin

   nvar = 65 + 1   ; my fault ... 

   flux_name = new((/nvar/),string)
   flux_p1 = new((/nvar/),integer)
   flux_p2 = new((/nvar/),integer)


   flux_name( 1) = "taux"      ;0  99  01: wind stress: E-W [kg/m/s2] 
   flux_p1( 1) = 1
   flux_p2( 1) = 1
   flux_name( 2) = "tauy"      ;0  99  02: wind stress: N-S [kg/m/s2]
   flux_p1( 2) = 2
   flux_p2( 2) = 2
   flux_name( 3) = "fsena"     ;0  99  03: sensible heat from canopy height to atmosphere [W/m2]
   flux_p1( 3) = 3
   flux_p2( 3) = 3
   flux_name( 4) = "lfevpa"    ;0  99  04: latent heat flux from canopy height to atmosphere [W/m2]
   flux_p1( 4) = 4
   flux_p2( 4) = 4
   flux_name( 5) = "fevpa"     ;0  99  05: evapotranspiration from canopy to atmosphere [mm/s]
   flux_p1( 5) = 5
   flux_p2( 5) = 5
   flux_name( 6) = "fsenl"     ;0  99  06: sensible heat from leaves [W/m2]
   flux_p1( 6) = 6
   flux_p2( 6) = 6
   flux_name( 7) = "fevpl"     ;0  99  07: evaporation+transpiration from leaves [mm/s]
   flux_p1( 7) = 7
   flux_p2( 7) = 7
   flux_name( 8) = "etr"       ;0  99  08: transpiration rate [mm/s]
   flux_p1( 8) = 8
   flux_p2( 8) = 8
   flux_name( 9) = "fseng"     ;0  99  09: sensible heat flux from ground [W/m2]
   flux_p1( 9) = 9
   flux_p2( 9) = 9
   flux_name(10) = "fevpg"     ;0  99  10: evaporation heat flux from ground [mm/s]
   flux_p1(10) = 10
   flux_p2(10) = 10
   flux_name(11) = "fgrnd"     ;0  99  11: ground heat flux [W/m2]
   flux_p1(11) = 11
   flux_p2(11) = 11
   flux_name(12) = "sabvsun"   ;0  99  12: solar absorbed by sunlit canopy [W/m2]
   flux_p1(12) = 12
   flux_p2(12) = 12
   flux_name(13) = "sabvsha"   ;0  99  13: solar absorbed by shaded [W/m2]
   flux_p1(13) = 13
   flux_p2(13) = 13
   flux_name(14) = "sabg"      ;0  99  14: solar absorbed by ground  [W/m2]
   flux_p1(14) = 14
   flux_p2(14) = 14
   flux_name(15) = "olrg"      ;0  99  15: outgoing long-wave radiation from ground+canopy [W/m2]
   flux_p1(15) = 15
   flux_p2(15) = 15
   flux_name(16) = "rnet"      ;0  99  16: net radiation [W/m2]
   flux_p1(16) = 16
   flux_p2(16) = 16
   flux_name(17) = "xerr"      ;0  99  17: the error of water banace [mm/s]
   flux_p1(17) = 17
   flux_p2(17) = 17
   flux_name(18) = "zerr"      ;0  99  18: the error of energy balance [W/m2]
   flux_p1(18) = 18
   flux_p2(18) = 18
   flux_name(19) = "rsur"      ;0  99  19: surface runoff [mm/s]
   flux_p1(19) = 19
   flux_p2(19) = 19
   flux_name(20) = "rnof"      ;0  99  20: total runoff [mm/s]
   flux_p1(20) = 20
   flux_p2(20) = 20
   flux_name(21) = "assim"     ;0  99  21: canopy assimilation rate [mol m-2 s-1]
   flux_p1(21) = 21
   flux_p2(21) = 21
   flux_name(22) = "respc"     ;0  99  22: respiration (plant+soil) [mol m-2 s-1]
   flux_p1(22) = 22
   flux_p2(22) = 22
   flux_name(23) = "tss"      ;10  99  23: soil temperature [K]
   flux_p1(23) = 23
   flux_p2(23) = 32
   flux_name(24) = "wliq"     ;10  99  24: liquid water in soil layers [kg/m2]
   flux_p1(24) = 33
   flux_p2(24) = 42
   flux_name(25) = "wice"     ;10  99  25: ice lens in soil layers [kg/m2]
   flux_p1(25) = 43
   flux_p2(25) = 52
   flux_name(26) = "tg"        ;0  99  26: ground surface temperature [K]
   flux_p1(26) = 53
   flux_p2(26) = 53
   flux_name(27) = "tlsun"     ;0  99  27: sunlit leaf temperature [K]
   flux_p1(27) = 54
   flux_p2(27) = 54
   flux_name(28) = "tlsha"     ;0  99  28: shaded leaf temperature [K]
   flux_p1(28) = 55
   flux_p2(28) = 55
   flux_name(29) = "ldew"      ;0  99  29: depth of water on foliage [mm]
   flux_p1(29) = 56
   flux_p2(29) = 56
   flux_name(30) = "scv"       ;0  99  30: snow cover, water equivalent [mm]
   flux_p1(30) = 57
   flux_p2(30) = 57
   flux_name(31) = "snowdp"    ;0  99  31: snow depth [meter]
   flux_p1(31) = 58
   flux_p2(31) = 58
   flux_name(32) = "fsno"      ;0  99  32: fraction of snow cover on ground
   flux_p1(32) = 59
   flux_p2(32) = 59
   flux_name(33) = "sigf"      ;0  99  33: fraction of veg cover, excluding snow-covered veg [-]
   flux_p1(33) = 60
   flux_p2(33) = 60
   flux_name(34) = "green"     ;0  99  34: leaf greenness
   flux_p1(34) = 61
   flux_p2(34) = 61
   flux_name(35) = "lai"       ;0  99  35: leaf area index
   flux_p1(35) = 62 
   flux_p2(35) = 62
   flux_name(36) = "sai"       ;0  99  36: stem area index
   flux_p1(36) = 63
   flux_p2(36) = 63
   flux_name(37) = "albvdir"   ;0  99  37: averaged albedo [vis, dir]
   flux_p1(37) = 64
   flux_p2(37) = 64
   flux_name(38) = "albvdif"   ;0  99  38: averaged albedo [vis, dif]
   flux_p1(38) = 65
   flux_p2(38) = 65
   flux_name(39) = "albndir"   ;0  99  39: averaged albedo [nir, dir]
   flux_p1(39) = 66
   flux_p2(39) = 66
   flux_name(40) = "albndif"   ;0  99  40: averaged albedo [nir, dif]
   flux_p1(40) = 67
   flux_p2(40) = 67
   flux_name(41) = "emis"      ;0  99  41: averaged bulk surface emissivity
   flux_p1(41) = 68
   flux_p2(41) = 68
   flux_name(42) = "z0ma"      ;0  99  42: effective roughness [m]
   flux_p1(42) = 69
   flux_p2(42) = 69
   flux_name(43) = "trad"      ;0  99  43: radiative temperature of surface [K]
   flux_p1(43) = 70
   flux_p2(43) = 70
   flux_name(44) = "ustar"     ;0  99  44: u* in similarity theory [m/s]
   flux_p1(44) = 71
   flux_p2(44) = 71
   flux_name(45) = "tstar"     ;0  99  45: t* in similarity theory [kg/kg]
   flux_p1(45) = 72
   flux_p2(45) = 72
   flux_name(46) = "qstar"     ;0  99  46: q* in similarity theory [kg/kg]
   flux_p1(46) = 73
   flux_p2(46) = 73
   flux_name(47) = "zol"       ;0  99  47: dimensionless height (z/L) used in Monin-Obukhov theory
   flux_p1(47) = 74
   flux_p2(47) = 74
   flux_name(48) = "rib"       ;0  99  48: bulk Richardson number in surface layer
   flux_p1(48) = 75
   flux_p2(48) = 75
   flux_name(49) = "fm"        ;0  99  49: integral of profile function for momentum
   flux_p1(49) = 76
   flux_p2(49) = 76
   flux_name(50) = "fh"        ;0  99  50: integral of profile function for heat
   flux_p1(50) = 77
   flux_p2(50) = 77
   flux_name(51) = "fq"        ;0  99  51: integral of profile function for moisture
   flux_p1(51) = 78
   flux_p2(51) = 78
   flux_name(52) = "tref"      ;0  99  52: 2 m height air temperature [kelvin]
   flux_p1(52) = 79
   flux_p2(52) = 79
   flux_name(53) = "qref"      ;0  99  53: 2 m height air specific humidity [kg/kg]
   flux_p1(53) = 80
   flux_p2(53) = 80
   flux_name(54) = "u10m"      ;0  99  54: 10m u-velocity [m/s]
   flux_p1(54) = 81
   flux_p2(54) = 81
   flux_name(55) = "v10m"      ;0  99  55: 10m v-velocity [m/s]
   flux_p1(55) = 82
   flux_p2(55) = 82
   flux_name(56) = "f10m"      ;0  99  56: integral of profile function for momentum at 10m [-]
   flux_p1(56) = 83
   flux_p2(56) = 83
   flux_name(57) = "us"        ;0  99  57: wind in eastward direction [m/s]
   flux_p1(57) = 84
   flux_p2(57) = 84
   flux_name(58) = "vs"        ;0  99  58: wind in northward direction [m/s]
   flux_p1(58) = 85
   flux_p2(58) = 85
   flux_name(59) = "tm"        ;0  99  59: temperature at reference height [kelvin]
   flux_p1(59) = 86
   flux_p2(59) = 86
   flux_name(60) = "qm"        ;0  99  60: specific humidity at reference height [kg/kg]
   flux_p1(60) = 87
   flux_p2(60) = 87
   flux_name(61) = "prc"       ;0  99  61: convective precipitation [mm/s]
   flux_p1(61) = 88
   flux_p2(61) = 88
   flux_name(62) = "prl"       ;0  99  62: large scale precipitation [mm/s]
   flux_p1(62) = 89
   flux_p2(62) = 89
   flux_name(63) = "pbot"      ;0  99  63: atmospheric pressure at the surface [pa]
   flux_p1(63) = 90
   flux_p2(63) = 90
   flux_name(64) = "frl"       ;0  99  64: atmospheric infrared (longwave) radiation [W/m2]
   flux_p1(64) = 91
   flux_p2(64) = 91
   flux_name(65) = "solar"     ;0  99  65: downward solar radiation at surface [W/m2]
   flux_p1(65) = 92
   flux_p2(65) = 92

   do k = 1, nvar
      if(vname .eq. flux_name(k)) then
         p1 = flux_p1(k)-1
         p2 = flux_p2(k)-1
         break
      end if
   end do

   delete(flux_name)
   delete(flux_p1)
   delete(flux_p2)

end ;procedure <looktable>

procedure rdflux (path:string, date:logical, vname:string, vbuf:float)

local fname, year, month, day, hour, second, maxday, deltim, ntime, itime, k, lscalar, lvector, lmatrix

begin

   if(date@delt .eq. "monthly") then
      year = date@year
      month = date@month
   end if

   if(date@delt .eq. "daily") then
      year = date@year
      day = date@day
   end if

   if(date@delt .eq. "6hourly") then
      year = date@year
      day = date@day
      hour = date@hour
   end if

   if(date@delt .eq. "hourly") then
      year = date@year
      day = date@day
      second = date@second
      deltim = date@deltim
   end if

   ntime = date@ntime

   p1 = 1
   p2 = 1

   looktable(vname,p1,p2)

   rank = dimsizes(dimsizes(vbuf))

   if (p1 .eq. p2) then
      rank = rank - 1
   else
      rank = rank - 2
   end if

   lscalar = False
   lvector = False
   lmatrix = False

   if (rank .eq. 0) then 
      lscalar = True
   end if
   if (rank .eq. 1) then
      lvector = True
   end if
   if (rank .eq. 2) then
      lmatrix = True
   end if

   itime = 0

   do while (itime .lt. ntime)

      if(date@delt .eq. "monthly") then
       ; print(sprinti("%0.4d ",year) + sprinti("%0.2d",month))

         fname = path + "/T42-" + sprinti("%0.4d",year) + "-" + sprinti("%0.2d",month) + "-84600"
      end if

      if(date@delt .eq. "daily") then
       ; print(sprinti("%0.4d ",year) + sprinti("%0.3d",day))

         fname = path + "/GLOBAL-" + sprinti("%0.4d",year) + "-" + sprinti("%0.3d",day) + "-86400"
      end if

      if(date@delt .eq. "6hourly") then
       ; print(sprinti("%0.4d ",year) + sprinti("%0.4d ",monthday(year,day)) + sprinti("%0.5d",hour)

         fname = path + "/clmflx." + sprinti("%0.4d",year) + sprinti("%0.4d",monthday(year,day)) + sprinti("%0.5d",hour)
      end if

      if(date@delt .eq. "hourly") then
       ; print(sprinti("%0.4d ",year) + sprinti("%0.3d ",day) + sprinti("%0.5d",second))

         fname = path + "/GLOBAL-" + sprinti("%0.4d",year) + "-" + sprinti("%0.3d",day) + "-" + sprinti("%0.5d",second)
      end if

      if(p1.eq.p2) then
         if(lscalar) then
            vbuf(itime) = fbinrecread(fname,p1,vbuf@dims,"float")
         end if
         if(lvector) then
            vbuf(itime,:) = fbinrecread(fname,p1,vbuf@dims,"float")
         end if
         if(lmatrix) then
            vbuf(itime,:,:) = fbinrecread(fname,p1,vbuf@dims,"float")
         end if
      end if

      if(p1.ne.p2) then
         do k = p1, p2
            if(lscalar) then
               vbuf(itime,k-p1) = fbinrecread(fname,k,vbuf@dims,"float")
            end if
            if(lvector) then
               vbuf(itime,k-p1,:) = fbinrecread(fname,k,vbuf@dims,"float")
            end if
            if(lmatrix) then
               vbuf(itime,k-p1,:,:) = fbinrecread(fname,k,vbuf@dims,"float")
            end if
         end do
      end if

      itime = itime+1

      if(date@delt .eq. "monthly") then
         month = month+1

         if(month .gt. 12) then
            month = 1
            year = year+1
         end if
      end if

      if(date@delt .eq. "daily") then
         day = day+1

         if(isleapyear(year)) then
            maxday = 366
         else
            maxday = 365
         end if

         if(day .gt. maxday) then
            day = 1
            year = year+1
         end if
      end if

      if(date@delt .eq. "6hourly") then
         hour = hour+6

         if(hour .ge. 24) then
            hour = hour-24
            day = day+1
         end if

         if(isleapyear(year)) then
            maxday = 366
         else
            maxday = 365
         end if

         if(day .gt. maxday) then
            day = 1
            year = year+1
         end if
      end if

      if(date@delt .eq. "hourly") then
         second = second+deltim

         if(second .gt. 86400) then
            second = second-86400
            day = day+1
         end if

         if(isleapyear(year)) then
            maxday = 366
         else
            maxday = 365
         end if

         if(day .gt. maxday) then
            day = 1
            year = year+1
         end if
      end if

   end do

end ;procedure <rdflux>


procedure plot_bar(wks:graphic, title:string, \
                   data1:numeric, data2:numeric, data3:numeric, data4:numeric)

local sres, plot1, lbres, labels

begin

  sres = True

  sres@vpWidthF = 0.3
  sres@vpHeightF = 0.3
  sres@trYMinF =  0
  sres@trYMaxF =  100
  sres@gsnDraw = True
  sres@gsnFrame = False
  sres@gsnXYBarChart = True
  sres@tmBorderThicknessF     = 1.0

 ;sres@tmXBMode   = "Explicit"         ; explicit labels
 ;sres@tmXBValues = (/1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25/)
 ;sres@tmXBLabels = (/"Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"/)

;;First frame, first plot

  sres@tmXBLabelFontHeightF = 0.015
  sres@tmXTLabelFontHeightF = 0.015
  sres@tmYLLabelFontHeightF = 0.015
  sres@tiMainFontHeightF = 0.015
  sres@tiMainFont = "helvetica"
  sres@gsnRightString = ""
  sres@tiYAxisString = "Percentage"

  sres@vpXF = 0.15
  sres@vpYF = 0.8
  sres@trXMinF = 0.6
  sres@trXMaxF = data1@size + 0.6
  sres@gsnXYBarChartBarWidth = 0.8*(data1@size/data1@size)                     ; change bar widths
  sres@tiMainString = data1@title

  sres@gsnXYBarChartColors = (/"red"/)
  plot1 = gsn_csm_xy(wks,fspan(1,data1@size,data1@size),data1,sres)

;;First frame, second plot

  sres@vpXF = 0.6
  sres@vpYF = 0.8
  sres@trXMinF = 0.6
  sres@trXMaxF = data2@size + 0.6
  sres@gsnXYBarChartBarWidth = 0.8*(data2@size/data1@size)                      ; change bar widths
  sres@tiMainString = data2@title

  sres@gsnXYBarChartColors = (/"red"/)	
  plot1 = gsn_csm_xy(wks,fspan(1,data2@size,data2@size),data2,sres)

;;First frame, third plot

  sres@vpXF = 0.15
  sres@vpYF = 0.4
  sres@trXMinF = 0.6
  sres@trXMaxF = data3@size + 0.6
  sres@gsnXYBarChartBarWidth = 0.8*(data3@size/data1@size)                      ; change bar widths
  sres@tiMainString = data3@title

  sres@gsnXYBarChartColors = (/"red"/)	
  plot1 = gsn_csm_xy(wks,fspan(1,data3@size,data3@size),data3,sres)

;;First frame, fourth plot

  sres@vpXF = 0.6
  sres@vpYF = 0.4
  sres@trXMinF = 0.6
  sres@trXMaxF = data4@size + 0.6
  sres@gsnXYBarChartBarWidth = 0.8*(data4@size/data3@size)                      ; change bar widths
  sres@tiMainString = data4@title

  sres@gsnXYBarChartColors = (/"red"/)	
  plot1 = gsn_csm_xy(wks,fspan(1,data4@size,data4@size),data4,sres)

  txres = True
  txres@txFontHeightF = 0.02  ; Change the font height.

  gsn_text_ndc(wks,title,0.5,0.9,txres)

  frame(wks)

end ;procedure <plot_bar>


procedure plot_bar2(wks:graphic, title:string, \
                    data1:numeric, data2:numeric, data3:numeric, data4:numeric)

local sres, plot1, plot2, lbres, labels, bar_width

begin

  sres = True

  sres@vpWidthF = 0.3
  sres@vpHeightF = 0.3
  sres@trYMinF =  0
  sres@trYMaxF =  100
  sres@gsnDraw = True
  sres@gsnFrame = False
  sres@gsnXYBarChart = True
  sres@tmBorderThicknessF     = 1.0

 ;sres@tmXBMode   = "Explicit"         ; explicit labels
 ;sres@tmXBValues = (/1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25/)
 ;sres@tmXBLabels = (/"Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"/)

;;First frame, first plot

  sres@tmXBLabelFontHeightF = 0.015
  sres@tmXTLabelFontHeightF = 0.015
  sres@tmYLLabelFontHeightF = 0.015
  sres@tiMainFontHeightF = 0.015
  sres@tiMainFont = "helvetica"
  sres@gsnRightString = ""
  sres@tiYAxisString = "Percentage"

  bar_width = 0.4*(data1@size/data1@size)                      ; change bar widths

  sres@vpXF = 0.15
  sres@vpYF = 0.85
  sres@trXMinF = 0.6
  sres@trXMaxF = data1@size + 0.6
  sres@gsnXYBarChartBarWidth = bar_width 
  sres@tiMainString = data1@title

  sres@gsnXYBarChartColors = (/"red"/)
  plot1 = gsn_csm_xy(wks,fspan(1-bar_width/2,data1@size-bar_width/2,data1@size),data1(0,:),sres)
  sres@gsnXYBarChartColors = (/"blue"/)
  plot2 = gsn_csm_xy(wks,fspan(1+bar_width/2,data1@size+bar_width/2,data1@size),data1(1,:),sres)

;;First frame, second plot

  bar_width = 0.4*(data2@size/data1@size)                      ; change bar widths

  sres@vpXF = 0.6
  sres@vpYF = 0.85
  sres@trXMinF = 0.6
  sres@trXMaxF = data2@size + 0.6
  sres@gsnXYBarChartBarWidth = bar_width
  sres@tiMainString = data2@title

  sres@gsnXYBarChartColors = (/"red"/)
  plot1 = gsn_csm_xy(wks,fspan(1-bar_width/2,data2@size-bar_width/2,data2@size),data2(0,:),sres)
  sres@gsnXYBarChartColors = (/"blue"/)
  plot2 = gsn_csm_xy(wks,fspan(1+bar_width/2,data2@size+bar_width/2,data2@size),data2(1,:),sres)

;;First frame, third plot

  bar_width = 0.4*(data3@size/data1@size)                      ; change bar widths

  sres@vpXF = 0.15
  sres@vpYF = 0.45
  sres@trXMinF = 0.6
  sres@trXMaxF = data3@size + 0.6
  sres@gsnXYBarChartBarWidth = bar_width
  sres@tiMainString = data3@title

  sres@gsnXYBarChartColors = (/"red"/)
  plot1 = gsn_csm_xy(wks,fspan(1-bar_width/2,data3@size-bar_width/2,data3@size),data3(0,:),sres)
  sres@gsnXYBarChartColors = (/"blue"/)
  plot2 = gsn_csm_xy(wks,fspan(1+bar_width/2,data3@size+bar_width/2,data3@size),data3(1,:),sres)

;;First frame, fourth plot

  bar_width = 0.4*(data4@size/data1@size)                      ; change bar widths

  sres@vpXF = 0.6
  sres@vpYF = 0.45
  sres@trXMinF = 0.6
  sres@trXMaxF = data4@size + 0.6
  sres@gsnXYBarChartBarWidth = bar_width
  sres@tiMainString = data4@title

  sres@gsnXYBarChartColors = (/"red"/)
  plot1 = gsn_csm_xy(wks,fspan(1-bar_width/2,data4@size-bar_width/2,data4@size),data4(0,:),sres)
  sres@gsnXYBarChartColors = (/"blue"/)
  plot2 = gsn_csm_xy(wks,fspan(1+bar_width/2,data4@size+bar_width/2,data4@size),data4(1,:),sres)

;;LabelBar

  lbres                    = True          ; labelbar only resources
  lbres@vpWidthF           = 0.3           ; labelbar width
  lbres@vpHeightF          = 0.1           ; labelbar height
  lbres@lbBoxMajorExtentF  = 0.1           ; puts space between color boxes
  lbres@lbFillColors       = (/"blue","red"/)
  lbres@lbMonoFillPattern  = True          ; Solid fill pattern
  lbres@lbLabelFontHeightF = 0.02          ; font height. default is small
  lbres@lbLabelJust        = "CenterLeft"  ; left justify labels
  lbres@lbPerimOn          = False
  lbres@lgPerimColor 	   = "white"
  labels = (/"clay","sand"/)
  gsn_labelbar_ndc(wks,2,labels,0.4,0.1,lbres)	; draw right labelbar column

;;Title

  txres = True
  txres@txFontHeightF = 0.015  ; Change the font height.
  gsn_text_ndc(wks,title,0.5,0.95,txres)

  frame(wks)

end ;procedure <plot_bar2>
