module getclouddepth_c
contains
           SUBROUTINE GETCLOUDDEPTH ( myid,T, PRESS, cwc, pwr, CWC_C, PWR_C,KBOTC, KTOPC, NZZ,&
                       RR,VOLRAT,RAINN,RAINC,LPREC)
!c-----Scan column for layers containing precipitation bottom/top
           INTEGER :: MYID,KBOTC , KTOPC ! THE BOTTOM AND TOP LAYER OF containing precipitation
           REAL    :: cwc(NZZ) ! CLOUD WATER CONTENT  (kg/kg)
           REAL    :: pwr(NZZ) !  rain water content (kg/kg)
           REAL    :: T(NZZ), PRESS(NZZ) ! IN K and HPA
           REAL    :: convfac ! conversion factor: umol/m3 = ppm * convfac   
           LOGICAL :: LPREC
           REAL    :: RR(NZZ),VOLRAT(NZZ)
           REAL    :: CWC_C(NZZ),PWR_C(NZZ) ! g/m3
           REAL    :: CWMIN,rhoh2o
           REAL    :: RAINN,RAINC ! WRF OUTPUT NON-CONVECTION and CONVECTION  RAIN (mm/hr)
           DATA CWMIN/ 0.05 /
           data rhoh2o /1.e6/     ! water density (g/m3) 
            
           LPREC = .TRUE.
           ktopc = 0
           kbotc = 0
           do k = 1, nzz
!--   cpnvert kg/kg to g/m3
           convfac = 44.9 * (273./t(k))*(press(k)/1013.)         
           cwc_c (k) = cwc (k) * 29. * convfac
           pwr_c (k) = pwr (k)* 29. * convfac           
           enddo
 
           do k = 1, nzz
            IF(cwc_c(k).ge.cwmin.or.pwr_c(k).ge.cwmin ) then
             kbotc = k
             goto 25
            ENDIF
           enddo
           
           LPREC = .FALSE.
           GOTO 20 
  25       CONTINUE
           IF (kbotc.eq.nzz) THEN
           LPREC = .FALSE.
           goto 20       
           ENDIF

           NCNT = 1
           DO K = KBOTC+1, NZZ
            IF(cwc_c(k).LT.cwmin.AND.pwr_c(k).LT.cwmin ) then
             KTOPC = K - 1
             GOTO 26
            ENDIF
            NCNT = NCNT + 1
           ENDDO     
           KTOPC = NZZ
  26       CONTINUE
           IF (KBOTC.gt.1 .and. ncnt.eq.1) THEN
            LPREC = .FALSE. 
            goto 20    
           ENDIF
      
  20       CONTINUE

!           IF(LPREC.AND.KBOTC.EQ.1.AND.KTOPC.GE.5) PRINT*,KBOTC,KTOPC
           RR  = 0.0
           VOLRAT = 0.0
           TOTALRR = 0.0

          IF (LPREC) THEN

          DO K = KBOTC, KTOPC
            VOLRAT(K) = ( PWR_C(K)+CWC_C(K) ) / rhoh2o ! drop volume/air volume
            rr(k) = (volrat(k)/1.0e-7)**1.27         ! rainfall rate (mm/hr)
            TOTALRR = TOTALRR + RR(K) 
          ENDDO  

           rate = 10.*ABS( RAINN + RAINC ) / TOTALRR
          IF(rate.EQ.0.0) THEN
            LPREC = .FALSE.
            GOTO 30
          ENDIF

          DO K = KBOTC, KTOPC
            rr(k) = rr(k) * rate
            cwc_c(k) = cwc_c(k) * rate
            pwr_c(k) = pwr_c(k) * rate
            VOLRAT (k) = rr(k)**(1./1.27)*1.E-7
          ENDDO
                  
          GOTO 30 !! by chenhs 
          ENDIF ! LPREC

          IF( 10.*(RAINN + RAINC ).GT.0.1 ) THEN ! >0.1mm/h, by chenhs
             KBOTC = 1
             KTOPC = 11
             LPREC = .TRUE.
           DO K = KBOTC, KTOPC
             rr(k) = 10.0 * (RAINN + RAINC)/11.
             cwc_c(k) = 0.5 * (rr(k)**(1./1.27)*rhoh2o*1.0e-7)
             pwr_c(k) = 0.5 * (rr(k)**(1./1.27)*rhoh2o*1.0e-7)
             VOLRAT (k) = rr(k)**(1./1.27)*1.E-7         
           ENDDO
          ENDIF           
 
  30       CONTINUE         
            
           RETURN
           END SUBROUTINE
end module getclouddepth_c       
