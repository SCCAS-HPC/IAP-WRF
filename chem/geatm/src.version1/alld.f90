module alld4main
contains
       SUBROUTINE ALLD(MYID, DSO4, DNO3,DSO4SSA, DNO3SSA, DSO4DUST,DNO3DUST,&
                        FSO4_SSA, FNO3_SSA, FSO4_DUST,FNO3_DUST)
! allocate the so4 and no3 from heterogeneous on dust and ssa 4 bins

         REAL :: DSO4SSA,DNO3SSA,DSO4DUST,DNO3DUST,FSO4_DUST,FNO3_DUST
         REAL :: FSO4_SSA,FNO3_SSA
         REAL :: DSO4,DNO3
         
         DSO4SSA  = DSO4SSA  + DSO4 * FSO4_SSA * 4.3 ! ppbv-ug/m3
         DNO3SSA  = DNO3SSA  + DNO3 * FNO3_SSA * 2.8
         DSO4DUST = DSO4DUST + DSO4 * FSO4_DUST* 4.3
         DNO3DUST = DNO3DUST + DNO3 * FNO3_DUST* 2.8
 
         DSO4SSA  = AMAX1 ( DSO4SSA, 1.E-20 )
         DNO3SSA  = AMAX1 ( DNO3SSA, 1.E-20 )
         DSO4DUST = AMAX1 ( DSO4DUST,1.E-20 )
         DNO3DUST = AMAX1 ( DNO3DUST,1.E-20 )


         RETURN
         END SUBROUTINE
end module alld4main        
