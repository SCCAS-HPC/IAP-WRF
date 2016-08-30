!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! MODULE MISC_DEFINITIONS_MODULE
!
! This module defines various non-meteorological constants that are used 
!   by other modules for readability.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module misc_definitions_module

   integer, parameter :: MAX_FILENAME_LEN = 1024

   real, parameter :: NAN=1.E20

   real, parameter :: NOT_MASKED   = -2.,  &
                      MASKED_BOTH  = -1.,  &
                      MASKED_WATER =  0.,  &
                      MASKED_LAND  =  1.

   integer, parameter :: OUTSIDE_DOMAIN=1E8, NOT_PROCESSED=1E9, INVALID=1E9

   integer, parameter :: SIXTEEN_POINT=1, FOUR_POINT=2, N_NEIGHBOR=3, &
                         AVERAGE4=4, AVERAGE16=5, W_AVERAGE4=6, W_AVERAGE16=7, &
                         SEARCH=8

   integer, parameter :: BOTTOM_TOP=1, TOP_BOTTOM=2

   integer, parameter :: CONTINUOUS=0, CATEGORICAL=1, SP_CONTINUOUS=2

   integer, parameter :: ONETWOONE=1, SMTHDESMTH=2, SMTHDESMTH_SPECIAL=3

   integer, parameter :: BIG_ENDIAN=0, LITTLE_ENDIAN=1

end module misc_definitions_module
