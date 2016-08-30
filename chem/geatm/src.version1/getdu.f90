module getdu_model
contains
          SUBROUTINE GETDU(MYID, PPBV, CUL, PRESS, TEMPK)
          INTEGER :: MYID
          REAL    :: PPBV, CUL ! PPBV GAS MIXING RATIO 
                               ! CUL : DENSITIVITY IN MOLESC/M3
          REAL    :: PRESS !HPA
          REAL    :: TEMPK !K
          
          CUL = PRESS*273./1013./TEMPK*1.E+5/22.4*6.02E+12*PPBV

          RETURN
          END subroutine

end module getdu_model
