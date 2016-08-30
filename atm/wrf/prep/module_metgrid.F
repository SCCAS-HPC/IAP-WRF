!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Program metgrid
!
! First version: Michael Duda -- January 2006 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  module module_metgrid

   use interp_option_module
   use module_debug
   use process_domain_module
   use module_domain
   use module_configure
   use module_camdomain  
      
   contains

   subroutine metgrid(mgrid, grid, n)
   
   type(metgrid_info) :: mgrid
   type(domain) :: grid
   integer, intent(in) :: n
   
   logical :: extra_col, extra_row
   integer :: ids, ide, jds, jde, kds, kde,    &
              ims, ime, jms, jme, kms, kme,    &    
              ips, ipe, jps, jpe, kps, kpe 
   
!   fort_test=200+mytask  !by Yuzhu Wang
   ! Having determined which processor we are, which grid type we are, and where 
   !   our patch is located in the domain, we can determine if U or V staggered 
   !   fields will have one more row or column than the M staggered fields
   	   
     call get_ijk_from_grid ( grid,    &
                              ids, ide, jds, jde, kds, kde,    &
                              ims, ime, jms, jme, kms, kme,    &    
	                      ips, ipe, jps, jpe, kps, kpe    )    
		       	       
       if (ipe .eq. ide) then ! One more column for U points
            extra_col = .true.
       else
            extra_col = .false.
       end if
  
       if (jpe .eq. jde ) then ! One more row for V points
            extra_row = .true.
       else
            extra_row = .false.
       end if

   ! Get info about which interpolators should be used with each field
       call  read_interp_table()

   ! Now begin the processing work, looping over all domains to be processed    
       call process_domain(mgrid, grid, n, extra_row, extra_col)

   ! Clean up and quit.
       call interp_option_destroy()

!       print *, ' *** Successful completion of metgrid *** '

      end subroutine metgrid

end module module_metgrid
