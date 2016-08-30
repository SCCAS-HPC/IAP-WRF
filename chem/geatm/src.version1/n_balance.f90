module nbalance
contains
       SUBROUTINE N_Balance(MYID,C,C_BALANCE,DELTX,DELTY,DELTZ,SX,EX,SY,EY)
         INTEGER MYID,SX,EX,SY,EY
         REAL DT,C_BALANCE
         REAL,DIMENSION(SX-1:EX+1,SY-1:EY+1) :: C,DELTX,DELTY,DELTZ
   
         DO J=SY,EY
         DO I=SX,EX
            IF(C(I,J).GT.0.0) THEN
            C_BALANCE=C_BALANCE+C(I,J)*DELTX(I,J)*DELTY(I,J)*DELTZ(I,J)*1.0E-12 !! Kg
            ENDIF
         ENDDO
         ENDDO

         RETURN
         END SUBROUTINE
        
       SUBROUTINE N_Balance_sum(MYID,C,C_BALANCE,DELTX,DELTY,DELTZ,SX,EX,SY,EY,DT)
         INTEGER MYID,SX,EX,SY,EY
         REAL DT,C_BALANCE
         REAL,DIMENSION(SX-1:EX+1,SY-1:EY+1) :: C,DELTX,DELTY,DELTZ

         DO J=SY,EY
         DO I=SX,EX
            C_BALANCE=C_BALANCE+C(I,J)*DELTX(I,J)*DELTY(I,J)
         ENDDO
         ENDDO

         RETURN
         END SUBROUTINE
         
        SUBROUTINE N_Balance_dep(MYID,C,C_BALANCE,DELTX,DELTY,DELTZ,SX,EX,SY,EY,DT)
         INTEGER MYID,SX,EX,SY,EY
         REAL DT,C_BALANCE
         REAL,DIMENSION(SX-1:EX+1,SY-1:EY+1) :: C,DELTX,DELTY,DELTZ

         DO J=SY,EY
         DO I=SX,EX
            C_BALANCE=C_BALANCE+C(I,J)*DELTX(I,J)*DELTY(I,J)
         ENDDO
         ENDDO

         RETURN
         END SUBROUTINE

       SUBROUTINE make_balance(MYID,C,total_bf,total_af,DELTX,DELTY,DELTZ,SX,EX,SY,EY)
         INTEGER MYID,SX,EX,SY,EY
         REAL total_bf,total_af,grid_amount,vol,ratio
         REAL,DIMENSION(SX-1:EX+1,SY-1:EY+1) :: C,DELTX,DELTY,DELTZ
        
         if(total_af.gt.0.0) then 
           ratio=(total_bf-total_af)/total_af
         else
           ratio=0.0
         endif

         if(abs(ratio).gt.1.0e-7) then  !! this reference ratio can be changed 
         DO J=SY,EY
         DO I=SX,EX
            vol=DELTX(I,J)*DELTY(I,J)*DELTZ(I,J)
            grid_amount=C(I,J)*vol

            C(I,J)=grid_amount*(1+ratio)/vol
            if(C(I,J).LT.0.0) then
               C(I,J)=0.0
            endif
         ENDDO
         ENDDO
         endif
          
         RETURN
         END  SUBROUTINE
end module nbalance         
