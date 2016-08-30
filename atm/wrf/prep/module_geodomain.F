module module_geodomain

   type geogrid_info
    real,dimension(:,:),allocatable:: ht_gc, tmn_gc, snoalb, OC12D, var2d, &
                           oa1, oa2, oa3, oa4, ol1, ol2, ol3, ol4, &
                           xland, lu_index, toposlpx, toposlpy,&
                           sct_dom_gc, scb_dom_gc, slopecat
    real,dimension(:,:),allocatable:: xlat, xlong, xlat_gc, xlong_gc, clat, clong, &
                            f, e, sina, cosa, msftx, msfty          
    real, dimension(:,:,:),allocatable:: albedo12m,  greenfrac
    real, dimension(:,:,:),allocatable:: soilctop, soilcbot
    real, dimension(:,:,:),allocatable:: landusef
    real, dimension(:,:),allocatable:: msfux, msfuy, xlat_u, xlong_u
    real, dimension(:,:),allocatable:: msfvx, msfvy, xlat_v, xlong_v
   end type geogrid_info
   
   type(geogrid_info):: ggrid_info
   
   contains
!===============================================================================   
   subroutine initial_ggrid(ggrid_info,ids, ide, jds, jde, kds, kde,    &
                               ims, ime, jms, jme, kms, kme,    &    
		       	       ips, ipe, jps, jpe, kps, kpe )
      integer,intent(in)::ids, ide, jds, jde, kds, kde,    &
                               ims, ime, jms, jme, kms, kme,    &    
		       	       ips, ipe, jps, jpe, kps, kpe 
      type(geogrid_info)::ggrid_info
   
      if(.not.allocated(ggrid_info%ht_gc))         allocate( ggrid_info%ht_gc(ims:ime, jms:jme) )
      if(.not.allocated(ggrid_info%tmn_gc))        allocate( ggrid_info%tmn_gc(ims:ime, jms:jme) )
      if(.not.allocated(ggrid_info%snoalb))        allocate( ggrid_info%snoalb(ims:ime, jms:jme) )
      if(.not.allocated(ggrid_info%OC12D))         allocate( ggrid_info%OC12D(ims:ime, jms:jme) )
      if(.not.allocated(ggrid_info%var2d))         allocate( ggrid_info%var2d(ims:ime, jms:jme) )
      if(.not.allocated(ggrid_info%oa1))           allocate( ggrid_info%oa1(ims:ime, jms:jme) )
      if(.not.allocated(ggrid_info%oa2))           allocate( ggrid_info%oa2(ims:ime, jms:jme) )
      if(.not.allocated(ggrid_info%oa3))           allocate( ggrid_info%oa3(ims:ime, jms:jme) )
      if(.not.allocated(ggrid_info%oa4))           allocate( ggrid_info%oa4(ims:ime, jms:jme) )
      if(.not.allocated(ggrid_info%ol1))           allocate( ggrid_info%ol1(ims:ime, jms:jme) )
      if(.not.allocated(ggrid_info%ol2))           allocate( ggrid_info%ol2(ims:ime, jms:jme) )
      if(.not.allocated(ggrid_info%ol3))           allocate( ggrid_info%ol3(ims:ime, jms:jme) )
      if(.not.allocated(ggrid_info%ol4))           allocate( ggrid_info%ol4(ims:ime, jms:jme) )
      if(.not.allocated(ggrid_info%xland))         allocate( ggrid_info%xland(ims:ime, jms:jme) )
      if(.not.allocated(ggrid_info%lu_index))      allocate( ggrid_info%lu_index(ims:ime, jms:jme) )
      if(.not.allocated(ggrid_info%toposlpx))      allocate( ggrid_info%toposlpx(ims:ime, jms:jme) )
      if(.not.allocated(ggrid_info%toposlpy))      allocate( ggrid_info%toposlpy(ims:ime, jms:jme) )
      if(.not.allocated(ggrid_info%sct_dom_gc))    allocate( ggrid_info%sct_dom_gc(ims:ime, jms:jme) )
      if(.not.allocated(ggrid_info%scb_dom_gc))    allocate( ggrid_info%scb_dom_gc(ims:ime, jms:jme) )
      if(.not.allocated(ggrid_info%slopecat))      allocate( ggrid_info%slopecat(ims:ime, jms:jme) )
      if(.not.allocated(ggrid_info%xlat))          allocate( ggrid_info%xlat(ims:ime, jms:jme) )
      if(.not.allocated(ggrid_info%xlong))         allocate( ggrid_info%xlong(ims:ime, jms:jme) )
      if(.not.allocated(ggrid_info%xlat_gc))       allocate( ggrid_info%xlat_gc(ims:ime, jms:jme) )
      if(.not.allocated(ggrid_info%xlong_gc))      allocate( ggrid_info%xlong_gc(ims:ime, jms:jme) )
      if(.not.allocated(ggrid_info%clat))          allocate( ggrid_info%clat(ims:ime, jms:jme) )
      if(.not.allocated(ggrid_info%clong))         allocate( ggrid_info%clong(ims:ime, jms:jme) )
      if(.not.allocated(ggrid_info%f))             allocate( ggrid_info%f(ims:ime, jms:jme) )
      if(.not.allocated(ggrid_info%e))             allocate( ggrid_info%e(ims:ime, jms:jme) )
      if(.not.allocated(ggrid_info%sina))          allocate( ggrid_info%sina(ims:ime, jms:jme) )
      if(.not.allocated(ggrid_info%cosa))          allocate( ggrid_info%cosa(ims:ime, jms:jme) )
      if(.not.allocated(ggrid_info%msftx))         allocate( ggrid_info%msftx(ims:ime, jms:jme) )
      if(.not.allocated(ggrid_info%msfty))         allocate( ggrid_info%msfty(ims:ime, jms:jme) )

      if(.not.allocated(ggrid_info%albedo12m))     allocate( ggrid_info%albedo12m(ims:ime, 1:12, jms:jme) )   
      if(.not.allocated(ggrid_info%greenfrac))     allocate( ggrid_info%greenfrac(ims:ime, 1:12, jms:jme) )
      if(.not.allocated(ggrid_info%soilctop))      allocate( ggrid_info%soilctop(ims:ime, 1:16, jms:jme) )
      if(.not.allocated(ggrid_info%soilcbot))      allocate( ggrid_info%soilcbot(ims:ime, 1:16, jms:jme) )
      if(.not.allocated(ggrid_info%landusef))      allocate( ggrid_info%landusef(ims:ime, 1:24, jms:jme) )
   
      if(.not.allocated( ggrid_info%msfux))        allocate( ggrid_info%msfux(ims:ime, jms:jme) )
      if(.not.allocated(ggrid_info%msfuy))         allocate( ggrid_info%msfuy(ims:ime, jms:jme) )
      if(.not.allocated(ggrid_info%xlat_u))        allocate( ggrid_info%xlat_u(ims:ime, jms:jme) )
      if(.not.allocated(ggrid_info%xlong_u))       allocate( ggrid_info%xlong_u(ims:ime, jms:jme) )
      if(.not.allocated(ggrid_info%msfvx))         allocate( ggrid_info%msfvx(ims:ime, jms:jme) )
      if(.not.allocated(ggrid_info%msfvy))         allocate( ggrid_info%msfvy(ims:ime, jms:jme) )
      if(.not.allocated(ggrid_info%xlat_v))        allocate( ggrid_info%xlat_v(ims:ime, jms:jme) )
      if(.not.allocated(ggrid_info%xlong_v))       allocate( ggrid_info%xlong_v(ims:ime, jms:jme) )  
           
  end subroutine initial_ggrid

!===============================================================================  
  subroutine final_ggrid(ggrid_info)
      type(geogrid_info)::ggrid_info
   
      if(allocated(ggrid_info%ht_gc))    deallocate( ggrid_info%ht_gc )
      if(allocated(ggrid_info%tmn_gc))   deallocate( ggrid_info%tmn_gc )
      if(allocated(ggrid_info%snoalb))   deallocate( ggrid_info%snoalb )
      if(allocated(ggrid_info%OC12D))    deallocate( ggrid_info%OC12D )
      if(allocated(ggrid_info%var2d))    deallocate( ggrid_info%var2d )
      if(allocated(ggrid_info%oa1))      deallocate( ggrid_info%oa1 )
      if(allocated(ggrid_info%oa2))      deallocate( ggrid_info%oa2 )
      if(allocated(ggrid_info%oa3))      deallocate( ggrid_info%oa3 )
      if(allocated(ggrid_info%oa4))      deallocate( ggrid_info%oa4 )
      if(allocated(ggrid_info%ol1))      deallocate( ggrid_info%ol1 )
      if(allocated(ggrid_info%ol2))      deallocate( ggrid_info%ol2 )
      if(allocated(ggrid_info%ol3))      deallocate( ggrid_info%ol3 )
      if(allocated(ggrid_info%ol4))      deallocate( ggrid_info%ol4 )
      if(allocated(ggrid_info%xland))     deallocate( ggrid_info%xland )
      if(allocated(ggrid_info%lu_index))  deallocate( ggrid_info%lu_index )
      if(allocated(ggrid_info%toposlpx))  deallocate( ggrid_info%toposlpx )
      if(allocated(ggrid_info%toposlpy))    deallocate( ggrid_info%toposlpy )
      if(allocated(ggrid_info%sct_dom_gc))  deallocate( ggrid_info%sct_dom_gc )
      if(allocated(ggrid_info%scb_dom_gc))  deallocate( ggrid_info%scb_dom_gc )
      if(allocated(ggrid_info%slopecat))    deallocate( ggrid_info%slopecat )
      if(allocated(ggrid_info%xlat))        deallocate( ggrid_info%xlat )
      if(allocated(ggrid_info%xlong))       deallocate( ggrid_info%xlong )
      if(allocated(ggrid_info%xlat_gc))     deallocate( ggrid_info%xlat_gc )
      if(allocated(ggrid_info%xlong_gc))    deallocate( ggrid_info%xlong_gc )
      if(allocated(ggrid_info%clat))        deallocate( ggrid_info%clat )
      if(allocated(ggrid_info%clong))       deallocate( ggrid_info%clong )
      if(allocated(ggrid_info%f))           deallocate( ggrid_info%f )
      if(allocated(ggrid_info%e))           deallocate( ggrid_info%e )
      if(allocated(ggrid_info%sina))        deallocate( ggrid_info%sina )
      if(allocated(ggrid_info%cosa))        deallocate( ggrid_info%cosa )
      if(allocated(ggrid_info%msftx))       deallocate( ggrid_info%msftx )
      if(allocated(ggrid_info%msfty))       deallocate( ggrid_info%msfty )

      if(allocated(ggrid_info%albedo12m))     deallocate( ggrid_info%albedo12m )   
      if(allocated(ggrid_info%greenfrac))     deallocate( ggrid_info%greenfrac )
      if(allocated(ggrid_info%soilctop))      deallocate( ggrid_info%soilctop )
      if(allocated(ggrid_info%soilcbot))      deallocate( ggrid_info%soilcbot )
      if(allocated(ggrid_info%landusef))      deallocate( ggrid_info%landusef )
   
      if(allocated( ggrid_info%msfux))        deallocate( ggrid_info%msfux)
      if(allocated(ggrid_info%msfuy))         deallocate( ggrid_info%msfuy )
      if(allocated(ggrid_info%xlat_u))        deallocate( ggrid_info%xlat_u )
      if(allocated(ggrid_info%xlong_u))       deallocate( ggrid_info%xlong_u )
      if(allocated(ggrid_info%msfvx))         deallocate( ggrid_info%msfvx)
      if(allocated(ggrid_info%msfvy))         deallocate( ggrid_info%msfvy )
      if(allocated(ggrid_info%xlat_v))        deallocate( ggrid_info%xlat_v)
      if(allocated(ggrid_info%xlong_v))       deallocate( ggrid_info%xlong_v)  
           
  end subroutine final_ggrid   
   
end module module_geodomain
