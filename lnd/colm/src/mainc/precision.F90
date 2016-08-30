module precision
!-------------------------------------------------------------------------------
! Purpose:
!	Define the precision to use for floating point and integer operations
!	throughout the model.
!-------------------------------------------------------------------------------

   integer, parameter :: r8 = selected_real_kind(12) ! 8 byte real
   integer, parameter :: r4 = selected_real_kind( 6) ! 4 byte real
   integer, parameter :: i8 = selected_int_kind (13) ! 8 byte integer
   integer, parameter :: i4 = selected_int_kind ( 6) ! 4 byte integer

   integer, parameter :: SHR_KIND_RN = kind(1.0)     ! native real
   integer, parameter :: SHR_KIND_IN = kind(1)       ! native integer
   integer, parameter :: SHR_KIND_CL = 256           ! long char

end module precision
