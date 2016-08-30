!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Program: geogrid
!
! Written by Michael G. Duda
! adapted by Juanxiong He for WRF/CAM coupling 	
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module module_geogrid	

   use gridinfo_module
   use llxy_module
   use module_llxy
   use list_module
   use module_debug
   use process_tile_module
   use source_data_module
   use parallel_module
         
   contains

   subroutine init_geogrid   
     
      ! Have the gridinfo module retrieve description of the grid setup
      call get_grid_params()

      ! Get information about the source data to be processed
      call get_datalist()

      ! Tell the llxy module that it can now compute parameters necessary to do 
      !   transformations for any nest 
      call compute_nest_locations()

   end subroutine init_geogrid

   subroutine geogrid(domain_id, ips, ipe, jps ,jpe, kps, kpe, &
               ids, ide, jds, jde, kds, kde, &
               ims, ime, jms, jme, kms, kme, &
               ht_gc, tmn_gc, snoalb, OC12D, var2d, &
               oa1, oa2, oa3, oa4, ol1, ol2, ol3, ol4, &
               xland, lu_index, toposlpx, toposlpy,&
               sct_dom_gc, scb_dom_gc, slopecat, &
               xlat, xlong, xlat_gc, xlong_gc, clat, clong, &
               f, e, sina, cosa, msftx, msfty, &    
               albedo12m, greenfrac, soilctop, soilcbot, &
               landusef, msfux, msfuy, xlat_u, xlong_u, &
               msfvx, msfvy, xlat_v, xlong_v)
  
   IMPLICIT NONE
   integer, intent(in) ::  domain_id   
   integer, intent(in)  :: ips, ipe, jps ,jpe, kps, kpe, &
               ids, ide, jds, jde, kds, kde, &
               ims, ime, jms, jme, kms, kme   
   real, dimension(ims:ime,jms:jme),intent(inout):: ht_gc, tmn_gc, snoalb, OC12D, var2d, &
                           oa1, oa2, oa3, oa4, ol1, ol2, ol3, ol4, &
                           xland, lu_index, toposlpx, toposlpy,&
                           sct_dom_gc, scb_dom_gc, slopecat
   real, dimension(ims:ime,jms:jme),intent(inout):: xlat, xlong, xlat_gc, xlong_gc, clat, clong, &
                            f, e, sina, cosa, msftx, msfty                                                                                   
   real, dimension(ims:ime,1:12,jms:jme),intent(inout):: albedo12m,  greenfrac
   real, dimension(ims:ime,1:16,jms:jme),intent(inout):: soilctop, soilcbot
   real, dimension(ims:ime,1:24,jms:jme),intent(inout):: landusef
   real, dimension(ims:ime,jms:jme),intent(inout):: msfux, msfuy, xlat_u, xlong_u
   real, dimension(ims:ime,jms:jme),intent(inout):: msfvx, msfvy, xlat_v, xlong_v   
   
   ! Local variables  
   logical :: extra_col, extra_row
   integer :: start_patch_i, end_patch_i, start_patch_j, end_patch_j, end_patch_stag_i, end_patch_stag_j, &
               start_mem_i, end_mem_i, start_mem_j, end_mem_j, end_mem_stag_i, end_mem_stag_j               
   
   integer:: i,j,k
  
            ! Get information about the source data we will use for this nest
            call get_source_params(geog_data_res(domain_id))
  
            ! Set transformations in llxy module to be with respect to current nest
            call select_domain(domain_id)
 
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
            
      ! The following pertains primarily to the C grid
      ! Determine whether only (n-1)th rows/columns should be computed for variables
      !   on staggered grid. In a distributed memory situation, not every tile should
      !   have only (n-1)th rows/columns computed, or we end up with (n-k) 
      !   rows/columns when there are k patches in the y/x direction
            if (extra_col) then
               start_patch_i    = ips    ! The seemingly pointless renaming of start
               end_patch_i      = ipe - 1  !   naming convention with modified end_patch variables, 
               end_patch_stag_i = ipe      !   variables is so that we can maintain consistent
                                                   !   which are marked as intent(in)
               start_mem_i    = start_patch_i    - HALO_WIDTH
               end_mem_i      = end_patch_i      + HALO_WIDTH
               end_mem_stag_i = end_patch_stag_i + HALO_WIDTH
            else                                     
               start_patch_i    = ips
               end_patch_i      = ipe
               end_patch_stag_i = ipe

               start_mem_i    = start_patch_i - HALO_WIDTH
               end_mem_i      = end_patch_i  + HALO_WIDTH
               end_mem_stag_i = end_patch_stag_i  + HALO_WIDTH
             end if
    
           if (extra_row) then
              start_patch_j    = jps
              end_patch_j      = jpe - 1
              end_patch_stag_j = jpe

              start_mem_j    = start_patch_j    - HALO_WIDTH
              end_mem_j      = end_patch_j      + HALO_WIDTH
              end_mem_stag_j = end_patch_stag_j + HALO_WIDTH
           else
              start_patch_j    = jps
              end_patch_j      = jpe
              end_patch_stag_j = jpe

              start_mem_j    = start_patch_j    - HALO_WIDTH
              end_mem_j      = end_patch_j      + HALO_WIDTH
              end_mem_stag_j = end_patch_stag_j + HALO_WIDTH
           end if
              
            ! Process fields for a tile of the current nest
            call process_tile(domain_id, start_patch_i, end_patch_i, start_patch_j, end_patch_j, end_patch_stag_i, end_patch_stag_j, &
                            start_mem_i, end_mem_i, start_mem_j, end_mem_j, end_mem_stag_i, end_mem_stag_j, &
                            ims, ime, jms, jme, extra_col, extra_row, &
                            ht_gc, tmn_gc, snoalb, OC12D, var2d, &
                            oa1, oa2, oa3, oa4, ol1, ol2, ol3, ol4, &
                            xland, lu_index, toposlpx, toposlpy,&
                            sct_dom_gc, scb_dom_gc, slopecat, &
                            xlat, xlong, xlat_gc, xlong_gc, clat, clong, &
                            f, e, sina, cosa, msftx, msfty, &    
                            albedo12m, greenfrac, soilctop, soilcbot, &
                            landusef, msfux, msfuy, xlat_u, xlong_u, &
                            msfvx, msfvy, xlat_v, xlong_v)
                                       
   end subroutine geogrid
   
end module module_geogrid
