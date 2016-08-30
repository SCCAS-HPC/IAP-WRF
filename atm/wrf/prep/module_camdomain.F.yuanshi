module module_camdomain
  
  use parallel_module 
  use module_dm
  USE module_comm_dm, ONLY : halo_em_phys_a_sub

  use module_domain
  
   type metgrid_info
      integer :: num_metgrid_levels, num_metgrid_soil_levels,iproj
      integer :: ids,  ide,  jds,  jde,  kds,  kde, &
                 ims,  ime,  jms,  jme,  kms,  kme, &
                 imsx, imex, jmsx, jmex, kmsx, kmex, &
                 imsy, imey, jmsy, jmey, kmsy, kmey, &
                 ips,  ipe,  jps,  jpe,  kps,  kpe, &
                 ipsx, ipex, jpsx, jpex, kpsx, kpex, &
                 ipsy, ipey, jpsy, jpey, kpsy, kpey 
      real, pointer, dimension(:,:):: xlat,xlon
      real, pointer, dimension(:,:):: tsk,sst,snow,xice,pslv,psfc,&
                                      xland,ht,sstold,xiceold
      real, pointer, dimension(:,:,:):: u3d,v3d,w3d,q3d,z3d,t3d,p3d,rh3d
      real, pointer, dimension(:,:,:):: soilt,soilm,soildepth,soilthick      
   end type metgrid_info

   type camgrid_info
      integer :: num_metgrid_levels, num_metgrid_soil_levels
      integer :: ids,  ide,  jds,  jde,  kds,  kde, &
                 ims,  ime,  jms,  jme,  kms,  kme, &
                 imsx, imex, jmsx, jmex, kmsx, kmex, &
                 imsy, imey, jmsy, jmey, kmsy, kmey, &
                 ips,  ipe,  jps,  jpe,  kps,  kpe, &
                 ipsx, ipex, jpsx, jpex, kpsx, kpex, &
                 ipsy, ipey, jpsy, jpey, kpsy, kpey
      real(8), pointer, dimension(:,:):: tsk,sst,snow,xice,pslv,psfc,xland,ht
      real(8), pointer, dimension(:,:,:):: u3d,v3d,w3d,q3d,z3d,t3d,p3d
      real(8), pointer, dimension(:,:,:):: uphyten,vphyten,tphyten,qphyten
      real(8), pointer, dimension(:,:,:):: dudt,dvdt,dtdt,dqdt
      real(8), pointer, dimension(:,:,:):: u1,u2,v1,v2,t1,t2,q1,q2
   end type camgrid_info

   type windgrid_info
      integer :: ids,  ide,  jds,  jde,  kds,  kde, &
                 ims,  ime,  jms,  jme,  kms,  kme, &
                 imsx, imex, jmsx, jmex, kmsx, kmex, &
                 imsy, imey, jmsy, jmey, kmsy, kmey, &
                 ips,  ipe,  jps,  jpe,  kps,  kpe, &
                 ipsx, ipex, jpsx, jpex, kpsx, kpex, &
                 ipsy, ipey, jpsy, jpey, kpsy, kpey
      real, pointer, dimension(:,:,:):: u3d,v3d,w3d,t3d,p3d 
      real, pointer, dimension(:,:):: xlat,xlong,psfc
   end type windgrid_info

   type(metgrid_info):: mgrid

   type(camgrid_info):: camgrid, camgrid_ac
   
   type(windgrid_info):: windgrid
   
   contains

!===============================================================================
    subroutine set_metgrid(mgrid)
      use module_configure
      type(metgrid_info) :: mgrid

      mgrid%num_metgrid_soil_levels = head_grid%num_metgrid_soil_levels
      mgrid%num_metgrid_levels =  head_grid%num_metgrid_levels
      mgrid%ids = 1
      mgrid%ide = model_config_rec%cam_nx
      mgrid%jds = 1
      mgrid%jde = model_config_rec%cam_ny
      mgrid%iproj =  model_config_rec%cam_proj

    end subroutine set_metgrid

!===============================================================================

    subroutine set_camgrid(camgrid)
      use module_configure

      type(camgrid_info) :: camgrid

      camgrid%num_metgrid_soil_levels = head_grid%num_metgrid_soil_levels
      camgrid%num_metgrid_levels =  head_grid%num_metgrid_levels
      camgrid%ids = 1
      camgrid%ide = model_config_rec%cam_nx
      camgrid%jds = 1
      camgrid%jde = model_config_rec%cam_ny

    end subroutine set_camgrid

!===============================================================================
    subroutine initial_metgrid(mgrid, ids, ide, jds, jde, &
                               num_metgrid_levels, num_metgrid_soil_levels )
      integer,intent(in)::ids, ide, jds, jde, num_metgrid_levels, num_metgrid_soil_levels 
      type(metgrid_info)::mgrid
   
      if(.not.associated(mgrid%u3d))         allocate( mgrid%u3d(ids:ide, 1:num_metgrid_levels, jds:jde) )
      if(.not.associated(mgrid%v3d))        allocate( mgrid%v3d(ids:ide, 1:num_metgrid_levels, jds:jde)  )
      if(.not.associated(mgrid%w3d))        allocate( mgrid%w3d(ids:ide, 1:num_metgrid_levels, jds:jde)  )
      if(.not.associated(mgrid%q3d))         allocate( mgrid%q3d(ids:ide, 1:num_metgrid_levels, jds:jde)  )
      if(.not.associated(mgrid%t3d))         allocate( mgrid%t3d(ids:ide, 1:num_metgrid_levels, jds:jde)  )
      if(.not.associated(mgrid%z3d))           allocate( mgrid%z3d(ids:ide, 1:num_metgrid_levels, jds:jde)  )
      if(.not.associated(mgrid%p3d))           allocate( mgrid%p3d(ids:ide, 1:num_metgrid_levels, jds:jde)  )
      if(.not.associated(mgrid%rh3d))           allocate( mgrid%rh3d(ids:ide, 1:num_metgrid_levels, jds:jde)  )
      if(.not.associated(mgrid%soilt))           allocate( mgrid%soilt(ids:ide, 1:num_metgrid_soil_levels, jds:jde)  )
      if(.not.associated(mgrid%soilm))          allocate( mgrid%soilm(ids:ide, 1:num_metgrid_soil_levels, jds:jde)  )
      if(.not.associated(mgrid%soilthick))     allocate( mgrid%soilthick(ids:ide, 1:num_metgrid_soil_levels, jds:jde)  )
      if(.not.associated(mgrid%soildepth))     allocate( mgrid%soildepth(ids:ide, 1:num_metgrid_soil_levels, jds:jde)  )
      if(.not.associated(mgrid%tsk))           allocate( mgrid%tsk(ids:ide, jds:jde)  )
      if(.not.associated(mgrid%sst))           allocate( mgrid%sst(ids:ide, jds:jde)  )
      if(.not.associated(mgrid%sstold))           allocate( mgrid%sstold(ids:ide, jds:jde)  )
      if(.not.associated(mgrid%snow))           allocate( mgrid%snow(ids:ide, jds:jde)  )
      if(.not.associated(mgrid%xice))           allocate( mgrid%xice(ids:ide, jds:jde)  )
      if(.not.associated(mgrid%xiceold))        allocate( mgrid%xiceold(ids:ide, jds:jde)  )
      if(.not.associated(mgrid%pslv))           allocate( mgrid%pslv(ids:ide, jds:jde)  )
      if(.not.associated(mgrid%psfc))           allocate( mgrid%psfc(ids:ide, jds:jde)  )
      if(.not.associated(mgrid%xlat))           allocate( mgrid%xlat(ids:ide, jds:jde)  )
      if(.not.associated(mgrid%xlon))           allocate( mgrid%xlon(ids:ide, jds:jde)  )        
      if(.not.associated(mgrid%xland))           allocate( mgrid%xland(ids:ide, jds:jde)  )        
      if(.not.associated(mgrid%ht))           allocate( mgrid%ht(ids:ide, jds:jde)  )        

      mgrid%u3d = 0
      mgrid%v3d = 0
      mgrid%w3d = 0
      mgrid%q3d = 0
      mgrid%t3d = 0
      mgrid%z3d = 0
      mgrid%p3d = 0
      mgrid%rh3d = 0
      mgrid%soilt = 0
      mgrid%soilm = 0
      mgrid%soilthick = 0
      mgrid%soildepth = 0
      mgrid%tsk = 0
      mgrid%sst = 0
      mgrid%sstold = 0
      mgrid%snow = 0
      mgrid%xice = 0
      mgrid%xiceold = 0
      mgrid%pslv = 0
      mgrid%psfc = 0
      mgrid%xlat = 0
      mgrid%xlon = 0
      mgrid%xland = 0    
      mgrid%ht = 0    
           
  end subroutine initial_metgrid

!===============================================================================

    subroutine initial_camgrid(camgrid, ids, ide, jds, jde, &
                               num_metgrid_levels, num_metgrid_soil_levels )
      use module_configure

      integer,intent(in)::ids, ide, jds, jde, num_metgrid_levels, num_metgrid_soil_levels
      type(camgrid_info)::camgrid
      logical :: twoway_coupling
      integer :: twoway_nudging

       twoway_coupling= model_config_rec%twoway_coupling
       twoway_nudging= model_config_rec%twoway_nudging

      if(.not.associated(camgrid%u3d))         allocate( camgrid%u3d(ids:ide, 1:num_metgrid_levels, jds:jde) )
      if(.not.associated(camgrid%v3d))        allocate( camgrid%v3d(ids:ide, 1:num_metgrid_levels, jds:jde)  )
      if(.not.associated(camgrid%w3d))        allocate( camgrid%w3d(ids:ide, 1:num_metgrid_levels, jds:jde)  )
      if(.not.associated(camgrid%q3d))         allocate( camgrid%q3d(ids:ide, 1:num_metgrid_levels, jds:jde)  )
      if(.not.associated(camgrid%t3d))         allocate( camgrid%t3d(ids:ide, 1:num_metgrid_levels, jds:jde)  )
      if(.not.associated(camgrid%z3d))           allocate( camgrid%z3d(ids:ide, 1:num_metgrid_levels, jds:jde)  )
      if(.not.associated(camgrid%p3d))           allocate( camgrid%p3d(ids:ide, 1:num_metgrid_levels, jds:jde)  )
      if (twoway_nudging.gt.0) then
      if(.not.associated(camgrid%dudt))           allocate( camgrid%dudt(ids:ide, 1:num_metgrid_levels, jds:jde)  )
      if(.not.associated(camgrid%dvdt))           allocate( camgrid%dvdt(ids:ide, 1:num_metgrid_levels, jds:jde)  )
      if(.not.associated(camgrid%dtdt))           allocate( camgrid%dtdt(ids:ide, 1:num_metgrid_levels, jds:jde)  )
      if(.not.associated(camgrid%dqdt))           allocate( camgrid%dqdt(ids:ide, 1:num_metgrid_levels, jds:jde)  )
      if(.not.associated(camgrid%uphyten))           allocate( camgrid%uphyten(ids:ide, 1:num_metgrid_levels, jds:jde)  )
      if(.not.associated(camgrid%vphyten))           allocate( camgrid%vphyten(ids:ide, 1:num_metgrid_levels, jds:jde)  )
      if(.not.associated(camgrid%tphyten))           allocate( camgrid%tphyten(ids:ide, 1:num_metgrid_levels, jds:jde)  )
      if(.not.associated(camgrid%qphyten))           allocate( camgrid%qphyten(ids:ide, 1:num_metgrid_levels, jds:jde)  )
      if(.not.associated(camgrid%u1))           allocate( camgrid%u1(ids:ide, 1:num_metgrid_levels, jds:jde)  )
      if(.not.associated(camgrid%u2))           allocate( camgrid%u2(ids:ide, 1:num_metgrid_levels, jds:jde)  )
      if(.not.associated(camgrid%v1))           allocate( camgrid%v1(ids:ide, 1:num_metgrid_levels, jds:jde)  )
      if(.not.associated(camgrid%v2))           allocate( camgrid%v2(ids:ide, 1:num_metgrid_levels, jds:jde)  )
      if(.not.associated(camgrid%t1))           allocate( camgrid%t1(ids:ide, 1:num_metgrid_levels, jds:jde)  )
      if(.not.associated(camgrid%t2))           allocate( camgrid%t2(ids:ide, 1:num_metgrid_levels, jds:jde)  )
      if(.not.associated(camgrid%q1))           allocate( camgrid%q1(ids:ide, 1:num_metgrid_levels, jds:jde)  )
      if(.not.associated(camgrid%q2))           allocate( camgrid%q2(ids:ide, 1:num_metgrid_levels, jds:jde)  )
      end if

      if(.not.associated(camgrid%tsk))           allocate( camgrid%tsk(ids:ide, jds:jde)  )
      if(.not.associated(camgrid%sst))           allocate( camgrid%sst(ids:ide, jds:jde)  )
      if(.not.associated(camgrid%snow))           allocate( camgrid%snow(ids:ide, jds:jde)  )
      if(.not.associated(camgrid%xice))           allocate( camgrid%xice(ids:ide, jds:jde)  )
      if(.not.associated(camgrid%pslv))           allocate( camgrid%pslv(ids:ide, jds:jde)  )
      if(.not.associated(camgrid%psfc))           allocate( camgrid%psfc(ids:ide, jds:jde)  )
      if(.not.associated(camgrid%xland))           allocate( camgrid%xland(ids:ide, jds:jde)  )
      if(.not.associated(camgrid%ht))           allocate( camgrid%ht(ids:ide, jds:jde)  )
 
      camgrid%u3d = 0.0_8
      camgrid%v3d = 0.0_8
      camgrid%w3d = 0.0_8
      camgrid%q3d = 0.0_8
      camgrid%t3d = 0.0_8
      camgrid%z3d = 0.0_8
      camgrid%p3d = 0.0_8

      if (twoway_nudging.gt.0) then
      camgrid%dudt = 0.0_8
      camgrid%dvdt = 0.0_8
      camgrid%dtdt = 0.0_8
      camgrid%dqdt = 0.0_8
      camgrid%uphyten = 0.0_8
      camgrid%vphyten = 0.0_8
      camgrid%tphyten = 0.0_8
      camgrid%qphyten = 0.0_8
      camgrid%u1 = 0.0_8
      camgrid%u2 = 0.0_8
      camgrid%v1 = 0.0_8
      camgrid%v2 = 0.0_8
      camgrid%t1 = 0.0_8
      camgrid%t2 = 0.0_8
      camgrid%q1 = 0.0_8
      camgrid%q2 = 0.0_8
      end if

      camgrid%tsk = 0.0_8
      camgrid%sst = 0.0_8
      camgrid%snow = 0.0_8
      camgrid%xice = 0.0_8
      camgrid%pslv = 0.0_8
      camgrid%psfc = 0.0_8
      camgrid%xland = 0.0_8
      camgrid%ht = 0.0_8

  end subroutine initial_camgrid

!===============================================================================
  subroutine initial_windgrid(windgrid, ids, ide, jds, jde, kps,kpe)
      use module_configure

      integer,intent(in)::ids, ide, jds, jde, kps, kpe
      type(windgrid_info)::windgrid

      if(.not.associated(windgrid%xlat))             allocate( windgrid%xlat(ids:ide, jds:jde) )
      if(.not.associated(windgrid%xlong))            allocate( windgrid%xlong(ids:ide, jds:jde) )
      if(.not.associated(windgrid%psfc))            allocate( windgrid%psfc(ids:ide, jds:jde) )
      if(.not.associated(windgrid%u3d))           allocate( windgrid%u3d(ids:(ide+1), kps:kpe, jds:jde)  )
      if(.not.associated(windgrid%v3d))           allocate( windgrid%v3d(ids:ide, kps:kpe, jds:(jde+1))  )
      if(.not.associated(windgrid%w3d))           allocate( windgrid%w3d(ids:ide, kps:(kpe+1), jds:jde)  )
      if(.not.associated(windgrid%t3d))           allocate( windgrid%t3d(ids:ide, kps:kpe, jds:jde)  )
      if(.not.associated(windgrid%p3d))           allocate( windgrid%p3d(ids:ide, kps:kpe, jds:jde)  )
      windgrid%xlat = 0
      windgrid%xlong = 0
      windgrid%psfc = 0
      windgrid%u3d = 0
      windgrid%v3d = 0
      windgrid%w3d = 0
      windgrid%t3d = 0
      windgrid%p3d = 0

  end subroutine initial_windgrid

!===============================================================================  
  subroutine final_metgrid(mgrid)
      use module_configure
      type(metgrid_info)::mgrid
   
      if(associated(mgrid%u3d))    deallocate( mgrid%u3d )
      if(associated(mgrid%v3d))   deallocate( mgrid%v3d )
      if(associated(mgrid%w3d))   deallocate( mgrid%w3d )
      if(associated(mgrid%q3d))    deallocate( mgrid%q3d )
      if(associated(mgrid%t3d))    deallocate( mgrid%t3d )
      if(associated(mgrid%z3d))      deallocate( mgrid%z3d )
      if(associated(mgrid%p3d))      deallocate( mgrid%p3d )
      if(associated(mgrid%rh3d))      deallocate( mgrid%rh3d )
      if(associated(mgrid%tsk))      deallocate( mgrid%tsk )
      if(associated(mgrid%sst))      deallocate( mgrid%sst )
      if(associated(mgrid%sstold))      deallocate( mgrid%sstold )
      if(associated(mgrid%pslv))      deallocate( mgrid%pslv )
      if(associated(mgrid%psfc))      deallocate( mgrid%psfc )
      if(associated(mgrid%xice))      deallocate( mgrid%xice )
      if(associated(mgrid%xiceold))      deallocate( mgrid%xiceold )
      if(associated(mgrid%snow))      deallocate( mgrid%snow )
      if(associated(mgrid%xlat))     deallocate( mgrid%xlat )
      if(associated(mgrid%xlon))     deallocate( mgrid%xlon )
      if(associated(mgrid%xland))     deallocate( mgrid%xland )
      if(associated(mgrid%soilt))    deallocate( mgrid%soilt)
      if(associated(mgrid%soilm))    deallocate( mgrid%soilm)
      if(associated(mgrid%soilthick)) deallocate( mgrid%soilthick )
      if(associated(mgrid%soildepth))  deallocate( mgrid%soildepth )    
      if(associated(mgrid%ht))  deallocate( mgrid%ht )    
           
  end subroutine final_metgrid

!===============================================================================  
  subroutine final_camgrid(camgrid)
      use module_configure
      type(camgrid_info)::camgrid
      logical :: twoway_coupling
      integer :: twoway_nudging

       twoway_coupling= model_config_rec%twoway_coupling
       twoway_nudging= model_config_rec%twoway_nudging

      if(associated(camgrid%u3d))    deallocate( camgrid%u3d )
      if(associated(camgrid%v3d))   deallocate( camgrid%v3d )
      if(associated(camgrid%w3d))   deallocate( camgrid%w3d )
      if(associated(camgrid%q3d))    deallocate( camgrid%q3d )
      if(associated(camgrid%t3d))    deallocate( camgrid%t3d )
      if(associated(camgrid%z3d))      deallocate( camgrid%z3d )
      if(associated(camgrid%p3d))      deallocate( camgrid%p3d )
      if (twoway_nudging.gt.0) then
      if(associated(camgrid%dudt))      deallocate( camgrid%dudt )
      if(associated(camgrid%dvdt))      deallocate( camgrid%dvdt )
      if(associated(camgrid%dtdt))      deallocate( camgrid%dtdt )
      if(associated(camgrid%dqdt))      deallocate( camgrid%dqdt ) 
      if(associated(camgrid%uphyten))      deallocate( camgrid%uphyten )
      if(associated(camgrid%vphyten))      deallocate( camgrid%vphyten )
      if(associated(camgrid%tphyten))      deallocate( camgrid%tphyten )
      if(associated(camgrid%qphyten))      deallocate( camgrid%qphyten )
      if(associated(camgrid%u1))      deallocate( camgrid%u1 )
      if(associated(camgrid%u2))      deallocate( camgrid%u2 )
      if(associated(camgrid%v1))      deallocate( camgrid%v1 )
      if(associated(camgrid%v2))      deallocate( camgrid%v2 )
      if(associated(camgrid%t1))      deallocate( camgrid%t1 )
      if(associated(camgrid%t2))      deallocate( camgrid%t2 )
      if(associated(camgrid%q1))      deallocate( camgrid%q1 )
      if(associated(camgrid%q2))      deallocate( camgrid%q2 )
      endif
      if(associated(camgrid%tsk))      deallocate( camgrid%tsk )
      if(associated(camgrid%sst))      deallocate( camgrid%sst )
      if(associated(camgrid%pslv))      deallocate( camgrid%pslv )
      if(associated(camgrid%psfc))      deallocate( camgrid%psfc )
      if(associated(camgrid%xice))      deallocate( camgrid%xice )
      if(associated(camgrid%snow))      deallocate( camgrid%snow )
      if(associated(camgrid%xland))     deallocate( camgrid%xland )
      if(associated(camgrid%ht))  deallocate( camgrid%ht )    
           
  end subroutine final_camgrid

!===============================================================================
  subroutine final_windgrid(windgrid)
      use module_configure
      type(windgrid_info)::windgrid

      if(associated(windgrid%xlat))     deallocate( windgrid%xlat )
      if(associated(windgrid%xlong))     deallocate( windgrid%xlong )
      if(associated(windgrid%psfc))     deallocate( windgrid%psfc )
      if(associated(windgrid%u3d))     deallocate( windgrid%u3d )
      if(associated(windgrid%v3d))     deallocate( windgrid%v3d )
      if(associated(windgrid%w3d))     deallocate( windgrid%w3d )
      if(associated(windgrid%t3d))     deallocate( windgrid%t3d )
      if(associated(windgrid%p3d))     deallocate( windgrid%p3d )

  end subroutine final_windgrid
  
!===============================================================================
  subroutine metgrid_decompose(grid, ids,  ide,  jds,  jde,  kds,  kde, &
                               ims,  ime,  jms,  jme,  kms,  kme, &
                               ips,  ipe,  jps,  jpe,  kps,  kpe )
    IMPLICIT NONE
    TYPE(METGRID_INFO)      ::  grid
    INTEGER, INTENT(INOUT)  ::  ids, ide, jds, jde, kds, kde
    INTEGER, INTENT(INOUT)  ::  ims, ime, jms, jme, kms, kme
    INTEGER, INTENT(INOUT)  ::  ips, ipe, jps, jpe, kps, kpe
    
    INTEGER    ::  imsx, imex, jmsx, jmex, kmsx, kmex, &
                    imsy, imey, jmsy, jmey, kmsy, kmey, &    
                    ipsx, ipex, jpsx, jpex, kpsx, kpex, &
                    ipsy, ipey, jpsy, jpey, kpsy, kpey                   


    ids = grid%ids 
    ide = grid%ide 
    jds = grid%jds 
    jde = grid%jde 
    kds = 1 
    kde = 1 

    call compute_memory_dims_rsl_lite( 1, 0,           &
                   0 , 0,  0,              &
                   ids,  ide,  jds,  jde,  kds,  kde, &
                   ims,  ime,  jms,  jme,  kms,  kme, &
                   imsx, imex, jmsx, jmex, kmsx, kmex, &
                   imsy, imey, jmsy, jmey, kmsy, kmey, &
                   ips,  ipe,  jps,  jpe,  kps,  kpe, &
                   ipsx, ipex, jpsx, jpex, kpsx, kpex, &
                   ipsy, ipey, jpsy, jpey, kpsy, kpey )

  end subroutine metgrid_decompose

!===============================================================================
	
   subroutine metgrid_alltogather(grid, ids, ide, jds, jde, kds, kde, &
                              num_metgrid_levels, num_metgrid_soil_levels, &
                              ips, ipe, jps, jpe, kps, kpe)
   implicit none
#if defined(DM_PARALLEL)
     include "mpif.h"
#endif
   
   type(metgrid_info)::grid
   integer,intent(in)::ids, ide, jds, jde, kds, kde,  &   
                       ips, ipe, jps, jpe, kps, kpe,  &
                       num_metgrid_levels, num_metgrid_soil_levels 
                               
   real, dimension(:,:,:), allocatable :: patch_array
   real, dimension(:,:,:), allocatable :: domain_array
   integer :: i, j, k, ierror
   
      allocate(patch_array(ips:ipe,1:num_metgrid_levels,jps:jpe))
      allocate(domain_array(ids:ide,1:num_metgrid_levels,jds:jde))
      
      ! p3d
      do j=jps, jpe
       do k=1,num_metgrid_levels
        do i=ips, ipe
          patch_array(i,k,j)=grid%p3d(i,k,j)
        enddo
       enddo
      enddo
      call gather_whole_field_r(patch_array, ips, ipe, 1, num_metgrid_levels, jps, jpe, &
                                domain_array, ids, ide, 1, num_metgrid_levels, jds, jde)
      call mpi_bcast(domain_array, (ide-ids+1)*(jde-jds+1)*num_metgrid_levels, MPI_REAL, 0, local_communicator, ierror)
      grid%p3d=domain_array

      ! z3d
      do j=jps, jpe
       do k=1,num_metgrid_levels
        do i=ips, ipe
          patch_array(i,k,j)=grid%z3d(i,k,j)
        enddo
       enddo
      enddo 
      call gather_whole_field_r(patch_array, ips, ipe, 1, num_metgrid_levels, jps, jpe, &
                                domain_array, ids, ide, 1, num_metgrid_levels, jds, jde)                        
      call mpi_bcast(domain_array, (ide-ids+1)*(jde-jds+1)*num_metgrid_levels, MPI_REAL, 0, local_communicator, ierror)     
      grid%z3d=domain_array
      
      ! t3d
      do j=jps, jpe
       do k=1,num_metgrid_levels
        do i=ips, ipe
          patch_array(i,k,j)=grid%t3d(i,k,j)
        enddo
       enddo
      enddo 
      call gather_whole_field_r(patch_array, ips, ipe, 1, num_metgrid_levels, jps, jpe, &
                                domain_array, ids, ide, 1, num_metgrid_levels, jds, jde)
      call mpi_bcast(domain_array, (ide-ids+1)*(jde-jds+1)*num_metgrid_levels, MPI_REAL, 0, local_communicator, ierror)   
      grid%t3d=domain_array
      
      ! q3d
      do j=jps, jpe
       do k=1,num_metgrid_levels
        do i=ips, ipe
          patch_array(i,k,j)=grid%q3d(i,k,j)
        enddo
       enddo
      enddo 
      call gather_whole_field_r(patch_array, ips, ipe, 1, num_metgrid_levels, jps, jpe, &
                                domain_array, ids, ide, 1, num_metgrid_levels, jds, jde)
      call mpi_bcast(domain_array, (ide-ids+1)*(jde-jds+1)*num_metgrid_levels, MPI_REAL, 0, local_communicator, ierror)   
      grid%q3d=domain_array
      
     ! u3d
      do j=jps, jpe
       do k=1,num_metgrid_levels
        do i=ips, ipe
          patch_array(i,k,j)=grid%u3d(i,k,j)
        enddo
       enddo
      enddo 
      call gather_whole_field_r(patch_array, ips, ipe, 1, num_metgrid_levels, jps, jpe, &
                                domain_array, ids, ide, 1, num_metgrid_levels, jds, jde)
      call mpi_bcast(domain_array, (ide-ids+1)*(jde-jds+1)*num_metgrid_levels, MPI_REAL, 0, local_communicator, ierror)   
      grid%u3d=domain_array

     ! v3d
      do j=jps, jpe
       do k=1,num_metgrid_levels
        do i=ips, ipe
          patch_array(i,k,j)=grid%v3d(i,k,j)
        enddo
       enddo
      enddo 
      call gather_whole_field_r(patch_array, ips, ipe, 1, num_metgrid_levels, jps, jpe, &
                                domain_array, ids, ide, 1, num_metgrid_levels, jds, jde)
      call mpi_bcast(domain_array, (ide-ids+1)*(jde-jds+1)*num_metgrid_levels, MPI_REAL, 0, local_communicator, ierror)                                                 
      grid%v3d=domain_array
      
     ! rh3d
      do j=jps, jpe
       do k=1,num_metgrid_levels
        do i=ips, ipe
          patch_array(i,k,j)=grid%rh3d(i,k,j)
        enddo
       enddo
      enddo 
      call gather_whole_field_r(patch_array, ips, ipe, 1, num_metgrid_levels, jps, jpe, &
                                domain_array, ids, ide, 1, num_metgrid_levels, jds, jde)
      call mpi_bcast(domain_array, (ide-ids+1)*(jde-jds+1)*num_metgrid_levels, MPI_REAL, 0, local_communicator, ierror)
      grid%rh3d=domain_array

      deallocate(patch_array)
      deallocate(domain_array)
      
      allocate(patch_array(ips:ipe,1:1,jps:jpe))
      allocate(domain_array(ids:ide,1:1,jds:jde))

      ! sst
      do j=jps, jpe
        do i=ips, ipe
          patch_array(i,1,j)=grid%sst(i,j)
        enddo
      enddo 
      call gather_whole_field_r(patch_array, ips, ipe, 1, 1, jps, jpe, &
                                domain_array, ids, ide, 1, 1, jds, jde)
      call mpi_bcast(domain_array, (ide-ids+1)*(jde-jds+1), MPI_REAL, 0, local_communicator, ierror)
      do j=jds, jde
        do i=ids, ide                           
         grid%sst(i,j)=domain_array(i,1,j)
        enddo
      enddo
      
      ! xice
      do j=jps, jpe
        do i=ips, ipe
          patch_array(i,1,j)=grid%xice(i,j)
        enddo
      enddo 
      call gather_whole_field_r(patch_array, ips, ipe, 1, 1, jps, jpe, &
                                domain_array, ids, ide, 1, 1, jds, jde)
      call mpi_bcast(domain_array, (ide-ids+1)*(jde-jds+1), MPI_REAL, 0, local_communicator, ierror)
      do j=jds, jde
        do i=ids, ide                           
         grid%xice(i,j)=domain_array(i,1,j)
        enddo
      enddo
      
      ! psfc
      do j=jps, jpe
        do i=ips, ipe
          patch_array(i,1,j)=grid%psfc(i,j)
        enddo
      enddo 
      call gather_whole_field_r(patch_array, ips, ipe, 1, 1, jps, jpe, &
                                domain_array, ids, ide, 1, 1, jds, jde)
      call mpi_bcast(domain_array, (ide-ids+1)*(jde-jds+1), MPI_REAL, 0, local_communicator, ierror)
      do j=jds, jde
        do i=ids, ide                           
         grid%psfc(i,j)=domain_array(i,1,j)
        enddo
      enddo
      
      ! pslv
      do j=jps, jpe
        do i=ips, ipe
          patch_array(i,1,j)=grid%pslv(i,j)
        enddo
      enddo 
      call gather_whole_field_r(patch_array, ips, ipe, 1, 1, jps, jpe, &
                                domain_array, ids, ide, 1, 1, jds, jde)
      call mpi_bcast(domain_array, (ide-ids+1)*(jde-jds+1), MPI_REAL, 0, local_communicator, ierror)
      do j=jds, jde
        do i=ids, ide                         
         grid%pslv(i,j)=domain_array(i,1,j)
        enddo
      enddo
      
      deallocate(patch_array)
      deallocate(domain_array)

   end subroutine metgrid_alltogather  

!===============================================================================

   subroutine metgrid_oncetogather(grid, ids, ide, jds, jde, kds, kde, &
                                num_metgrid_levels, num_metgrid_soil_levels, &
                                ips, ipe, jps, jpe, kps, kpe)
   implicit none
#if defined(DM_PARALLEL)
     include "mpif.h"
#endif
   
   type(metgrid_info)::grid
   integer,intent(in)::ids, ide, jds, jde, kds, kde,  &   
                       ips, ipe, jps, jpe, kps, kpe,  &
                       num_metgrid_levels, num_metgrid_soil_levels 
                               
   real, dimension(:,:,:), allocatable :: patch_array
   real, dimension(:,:,:), allocatable :: domain_array
   integer :: i, j, k, ierror
   
      allocate(patch_array(ips:ipe,1:num_metgrid_soil_levels,jps:jpe))
      allocate(domain_array(ids:ide,1:num_metgrid_soil_levels,jds:jde))

      ! soilt
      do j=jps, jpe
       do k=1,num_metgrid_soil_levels
        do i=ips, ipe
          patch_array(i,k,j)=grid%soilt(i,k,j)
        enddo
       enddo
      enddo 
      call gather_whole_field_r(patch_array, ips, ipe, 1, num_metgrid_soil_levels, jps, jpe, &
                                 domain_array, ids, ide, 1, num_metgrid_soil_levels, jds, jde)
      call mpi_bcast(domain_array, (ide-ids+1)*(jde-jds+1)*num_metgrid_soil_levels, MPI_REAL, 0, local_communicator, ierror)                                    
      grid%soilt=domain_array
      
      ! soilm
      do j=jps, jpe
       do k=1,num_metgrid_soil_levels
        do i=ips, ipe
          patch_array(i,k,j)=grid%soilm(i,k,j)
        enddo
       enddo
      enddo 
      call gather_whole_field_r(patch_array, ips, ipe, 1, num_metgrid_soil_levels, jps, jpe, &
                                 domain_array, ids, ide, 1, num_metgrid_soil_levels, jds, jde)
      call mpi_bcast(domain_array, (ide-ids+1)*(jde-jds+1)*num_metgrid_soil_levels, MPI_REAL, 0, local_communicator, ierror)
      grid%soilm=domain_array      
      
     ! soilthick
      do j=jps, jpe
       do k=1,num_metgrid_soil_levels
        do i=ips, ipe
          patch_array(i,k,j)=grid%soilthick(i,k,j)
        enddo
       enddo
      enddo 
      call gather_whole_field_r(patch_array, ips, ipe, 1, num_metgrid_soil_levels, jps, jpe, &
                                 domain_array, ids, ide, 1, num_metgrid_soil_levels, jds, jde)
      call mpi_bcast(domain_array, (ide-ids+1)*(jde-jds+1)*num_metgrid_soil_levels, MPI_REAL, 0, local_communicator, ierror)
      grid%soilthick=domain_array      

      ! soildepth
      do j=jps, jpe
       do k=1,num_metgrid_soil_levels
        do i=ips, ipe
          patch_array(i,k,j)=grid%soildepth(i,k,j)
        enddo
       enddo
      enddo 
      call gather_whole_field_r(patch_array, ips, ipe, 1, num_metgrid_soil_levels, jps, jpe, & 
                                domain_array, ids, ide, 1, num_metgrid_soil_levels, jds, jde)
      call mpi_bcast(domain_array, (ide-ids+1)*(jde-jds+1)*num_metgrid_soil_levels, MPI_REAL, 0, local_communicator, ierror)
      grid%soildepth=domain_array
      
      deallocate(patch_array)
      deallocate(domain_array)
      
      allocate(patch_array(ips:ipe,1:1,jps:jpe))
      allocate(domain_array(ids:ide,1:1,jds:jde))
      
      ! tsk
      do j=jps, jpe
        do i=ips, ipe
          patch_array(i,1,j)=grid%tsk(i,j)
        enddo
      enddo 
      call gather_whole_field_r(patch_array, ips, ipe, 1, 1, jps, jpe, &
                                 domain_array, ids, ide, 1, 1, jds, jde)
      call mpi_bcast(domain_array, (ide-ids+1)*(jde-jds+1), MPI_REAL, 0, local_communicator, ierror)
      do j=jds, jde
        do i=ids, ide                           
         grid%tsk(i,j)=domain_array(i,1,j)
        enddo
      enddo
      
      ! snow
      do j=jps, jpe
        do i=ips, ipe
          patch_array(i,1,j)=grid%snow(i,j)
        enddo
      enddo 
      call gather_whole_field_r(patch_array, ips, ipe, 1, 1, jps, jpe, &
                                 domain_array, ids, ide, 1, 1, jds, jde)
      call mpi_bcast(domain_array, (ide-ids+1)*(jde-jds+1), MPI_REAL, 0, local_communicator, ierror)
      do j=jds, jde
        do i=ids, ide                           
         grid%snow(i,j)=domain_array(i,1,j)
        enddo
      enddo
                
      ! xlat
      do j=jps, jpe
        do i=ips, ipe
          patch_array(i,1,j)=grid%xlat(i,j)
        enddo
      enddo 
      call gather_whole_field_r(patch_array, ips, ipe, 1, 1, jps, jpe, &
                                 domain_array, ids, ide, 1, 1, jds, jde)
      call mpi_bcast(domain_array, (ide-ids+1)*(jde-jds+1), MPI_REAL, 0, local_communicator, ierror)
      do j=jds, jde
        do i=ids, ide   
         grid%xlat(i,j)=domain_array(i,1,j)
        enddo
      enddo
      
      ! xlon
      do j=jps, jpe
        do i=ips, ipe
          patch_array(i,1,j)=grid%xlon(i,j)
        enddo
      enddo 
      call gather_whole_field_r(patch_array, ips, ipe, 1, 1, jps, jpe, &
                                 domain_array, ids, ide, 1, 1, jds, jde)
      call mpi_bcast(domain_array, (ide-ids+1)*(jde-jds+1), MPI_REAL, 0, local_communicator, ierror)
      do j=jds, jde
        do i=ids, ide                           
         grid%xlon(i,j)=domain_array(i,1,j)
        enddo
      enddo

     ! xland
      do j=jps, jpe
        do i=ips, ipe
          patch_array(i,1,j)=grid%xland(i,j)
        enddo
      enddo
      call gather_whole_field_r(patch_array, ips, ipe, 1, 1, jps, jpe, &
                                 domain_array, ids, ide, 1, 1, jds, jde)
      call mpi_bcast(domain_array, (ide-ids+1)*(jde-jds+1), MPI_REAL, 0, local_communicator, ierror)
      do j=jds, jde
        do i=ids, ide
         grid%xland(i,j)=domain_array(i,1,j)
        enddo
      enddo

      ! ht
      do j=jps, jpe
        do i=ips, ipe
          patch_array(i,1,j)=grid%ht(i,j)
        enddo
      enddo
      call gather_whole_field_r(patch_array, ips, ipe, 1, 1, jps, jpe, &
                                 domain_array, ids, ide, 1, 1, jds, jde)
      call mpi_bcast(domain_array, (ide-ids+1)*(jde-jds+1), MPI_REAL, 0, local_communicator, ierror)
      do j=jds, jde
        do i=ids, ide
         grid%ht(i,j)=domain_array(i,1,j)
        enddo
      enddo

      deallocate(patch_array)
      deallocate(domain_array)

   end subroutine metgrid_oncetogather  

!===============================================================================

   subroutine windgrid_map_alltogather(windgrid, ids, ide, jds, jde, kds, kde, &
                              num_cam_levels, &
                              ips, ipe, jps, jpe, kps, kpe)
   implicit none
#if defined(DM_PARALLEL)
     include "mpif.h"
#endif

   type(windgrid_info)::windgrid
   integer,intent(in)::ids, ide, jds, jde, kds, kde,  &
                       ips, ipe, jps, jpe, kps, kpe,  &
                       num_cam_levels

   real, dimension(:,:,:), allocatable :: patch_array
   real, dimension(:,:,:), allocatable :: domain_array
   integer :: i, j, k, ierror
   integer :: mpi_ierr

      allocate(patch_array(ips:ipe,1:1,jps:jpe))
      allocate(domain_array(ids:ide,1:1,jds:jde))

      !  latitude
      do j=jps, jpe
        do i=ips, ipe
          patch_array(i,1,j)=windgrid%xlat(i,j)
        enddo
      enddo 
      call gather_whole_field_r(patch_array, ips, ipe, 1, 1, jps, jpe, &
                                 domain_array, ids, ide, 1, 1, jds, jde)
      call mpi_bcast(domain_array, (ide-ids+1)*(jde-jds+1), MPI_REAL, 0, local_communicator, ierror)
      do j=jds, jde
        do i=ids, ide                            
         windgrid%xlat(i,j)=domain_array(i,1,j)
        enddo
      enddo 
      
      ! longitude
      do j=jps, jpe
        do i=ips, ipe
          patch_array(i,1,j)=windgrid%xlong(i,j)
        enddo
      enddo
      call gather_whole_field_r(patch_array, ips, ipe, 1, 1, jps, jpe, &
                                 domain_array, ids, ide, 1, 1, jds, jde)
      call mpi_bcast(domain_array, (ide-ids+1)*(jde-jds+1), MPI_REAL, 0, local_communicator, ierror)
      do j=jds, jde
        do i=ids, ide        
         windgrid%xlong(i,j)=domain_array(i,1,j)
        enddo
      enddo

      deallocate(patch_array)
      deallocate(domain_array)

   end subroutine windgrid_map_alltogather

!===============================================================================
   subroutine windgrid_p_alltogather(windgrid, ids, ide, jds, jde, kds, kde, &
                              num_cam_levels, &
                              ips, ipe, jps, jpe, kps, kpe)
   implicit none
#if defined(DM_PARALLEL)
     include "mpif.h"
#endif

   type(windgrid_info)::windgrid
   integer,intent(in)::ids, ide, jds, jde, kds, kde,  &   
                       ips, ipe, jps, jpe, kps, kpe,  &
                       num_cam_levels
              
   real, dimension(:,:,:), allocatable :: patch_array
   real, dimension(:,:,:), allocatable :: domain_array
   integer :: i, j, k, ierror

      
      allocate(patch_array(ips:ipe,1:1,jps:jpe))
      allocate(domain_array(ids:ide,1:1,jds:jde))

      !  psfc
      do j=jps, jpe
        do i=ips, ipe
          patch_array(i,1,j)=windgrid%psfc(i,j)
        enddo
      enddo 
      call gather_whole_field_r(patch_array, ips, ipe, 1, 1, jps, jpe, &
                                 domain_array, ids, ide, 1, 1, jds, jde)
      call mpi_bcast(domain_array, (ide-ids+1)*(jde-jds+1), MPI_REAL, 0,local_communicator, ierror)
      do j=jds, jde
        do i=ids, ide                           
         windgrid%psfc(i,j)=domain_array(i,1,j)
        enddo
      enddo
      
      deallocate(patch_array)
      deallocate(domain_array)

   end subroutine windgrid_p_alltogather   
!============================================================================== 

 subroutine windgrid_alltogather(windgrid_va, ids, ide, jds, jde, kds, kde, &
                                 ips, ipe, jps, jpe, kps, kpe)
   implicit none
#if defined(DM_PARALLEL)
     include "mpif.h"
#endif

   integer,intent(in)::ids, ide, jds, jde, kds, kde,  &
                       ips, ipe, jps, jpe, kps, kpe
   real, dimension(ids:ide,kps:kpe,jds:jde)::windgrid_va

   real, dimension(:,:,:), allocatable :: patch_array
   real, dimension(:,:,:), allocatable :: domain_array
   integer :: i, j, k, ierror

      allocate(patch_array(ips:ipe,1:(kpe-kps+1),jps:jpe))
      allocate(domain_array(ids:ide,1:(kpe-kps+1),jds:jde))

     ! u3d
      do j=jps, jpe
       do k= 1, kpe-kps+1
        do i=ips, ipe
          patch_array(i,k,j)=windgrid_va(i,k,j)
        enddo
       enddo
      enddo 
      call gather_whole_field_r(patch_array, ips, ipe, 1, kpe-kps+1, jps, jpe, &
                                 domain_array, ids, ide, 1, kpe-kps+1, jds, jde)

      call mpi_bcast(domain_array, (ide-ids+1)*(jde-jds+1)*(kpe-kps+1),MPI_REAL, 0, local_communicator, ierror)
      windgrid_va=domain_array

      deallocate(patch_array)
      deallocate(domain_array)

   end subroutine windgrid_alltogather

!===============================================================================
   subroutine camgrid_gather_state(grid, ids, ide, jds, jde, kds, kde, &
                              num_metgrid_levels, num_metgrid_soil_levels, &
                              ips, ipe, jps, jpe, kps, kpe)
   implicit none
#if defined(DM_PARALLEL)
     include "mpif.h"
#endif
   
   type(camgrid_info)::grid
   integer,intent(in)::ids, ide, jds, jde, kds, kde,  &   
                       ips, ipe, jps, jpe, kps, kpe,  &
                       num_metgrid_levels, num_metgrid_soil_levels 
                               
   real(8), dimension(:,:,:), allocatable :: patch_array
   real(8), dimension(:,:,:), allocatable :: domain_array
   integer :: i, j, k, mpi_ierr
   
      allocate(patch_array(ips:ipe,1:num_metgrid_levels,jps:jpe))
      allocate(domain_array(ids:ide,1:num_metgrid_levels,jds:jde))

      ! p3d
      do j=jps, jpe
       do k=1,num_metgrid_levels
        do i=ips, ipe
          patch_array(i,k,j)=grid%p3d(i,k,j)
        enddo
       enddo
      enddo
      call gather_whole_field_db(patch_array, ips, ipe, 1,num_metgrid_levels,jps, jpe, &
                                 domain_array, ids, ide, 1,num_metgrid_levels,jds, jde)
      if(mytask.eq.0)   grid%p3d=domain_array
      call mpi_barrier(local_communicator, mpi_ierr)

      ! u3d
      do j=jps, jpe
       do k=1,num_metgrid_levels
        do i=ips, ipe
          patch_array(i,k,j)=grid%u3d(i,k,j)
        enddo
       enddo
      enddo
      call gather_whole_field_db(patch_array, ips, ipe, 1, num_metgrid_levels,jps, jpe, &
                                 domain_array, ids, ide, 1, num_metgrid_levels,jds, jde)
      if(mytask.eq.0)   grid%u3d=domain_array
      call mpi_barrier(local_communicator, mpi_ierr)

      ! v3d
      do j=jps, jpe
       do k=1,num_metgrid_levels
        do i=ips, ipe
          patch_array(i,k,j)=grid%v3d(i,k,j)
        enddo
       enddo
      enddo
      call gather_whole_field_db(patch_array, ips, ipe, 1, num_metgrid_levels,jps, jpe, &
                                 domain_array, ids, ide, 1, num_metgrid_levels,jds, jde)
      if(mytask.eq.0)   grid%v3d=domain_array
      call mpi_barrier(local_communicator, mpi_ierr)

      ! t3d
      do j=jps, jpe
       do k=1,num_metgrid_levels
        do i=ips, ipe
          patch_array(i,k,j)=grid%t3d(i,k,j)
        enddo
       enddo
      enddo
      call gather_whole_field_db(patch_array, ips, ipe, 1, num_metgrid_levels,jps, jpe, &
                                 domain_array, ids, ide, 1, num_metgrid_levels,jds, jde)
      if(mytask.eq.0)   grid%t3d=domain_array
      call mpi_barrier(local_communicator, mpi_ierr)

      ! q3d
      do j=jps, jpe
       do k=1,num_metgrid_levels
        do i=ips, ipe
          patch_array(i,k,j)=grid%q3d(i,k,j)
        enddo
       enddo
      enddo
      call gather_whole_field_db(patch_array, ips, ipe, 1, num_metgrid_levels,jps, jpe, &
                                 domain_array, ids, ide, 1, num_metgrid_levels,jds, jde)
      if(mytask.eq.0)   grid%q3d=domain_array
      call mpi_barrier(local_communicator, mpi_ierr)

      ! w3d
      do j=jps, jpe
       do k=1,num_metgrid_levels
        do i=ips, ipe
          patch_array(i,k,j)=grid%w3d(i,k,j)
        enddo
       enddo
      enddo
      call gather_whole_field_db(patch_array, ips, ipe, 1,num_metgrid_levels,jps, jpe, &
                                 domain_array, ids, ide, 1,num_metgrid_levels,jds, jde)
      if(mytask.eq.0)   grid%w3d=domain_array
      call mpi_barrier(local_communicator, mpi_ierr)

      ! z3d
      do j=jps, jpe
       do k=1,num_metgrid_levels
        do i=ips, ipe
          patch_array(i,k,j)=grid%z3d(i,k,j)
        enddo
       enddo
      enddo
      call gather_whole_field_db(patch_array, ips, ipe, 1, num_metgrid_levels,jps, jpe, &
                                 domain_array, ids, ide, 1, num_metgrid_levels,jds, jde)
      if(mytask.eq.0)   grid%z3d=domain_array
      call mpi_barrier(local_communicator, mpi_ierr)

      deallocate(patch_array)
      deallocate(domain_array)

   end subroutine camgrid_gather_state

!===============================================================================
   subroutine camgrid_gather_tend(grid, ids, ide, jds, jde, kds, kde, &
                              num_metgrid_levels, num_metgrid_soil_levels, &
                              ips, ipe, jps, jpe, kps, kpe)
   implicit none
#if defined(DM_PARALLEL)
     include "mpif.h"
#endif

   type(camgrid_info)::grid
   integer,intent(in)::ids, ide, jds, jde, kds, kde,  &
                       ips, ipe, jps, jpe, kps, kpe,  &
                       num_metgrid_levels, num_metgrid_soil_levels

   real(8), dimension(:,:,:), allocatable :: patch_array
   real(8), dimension(:,:,:), allocatable :: domain_array
   integer :: i, j, k, mpi_ierr

      allocate(patch_array(ips:ipe,1:num_metgrid_levels,jps:jpe))
      allocate(domain_array(ids:ide,1:num_metgrid_levels,jds:jde))

      ! u1
      do j=jps, jpe
       do k=1,num_metgrid_levels
        do i=ips, ipe
          patch_array(i,k,j)=grid%u1(i,k,j)
        enddo
       enddo
      enddo 
      call gather_whole_field_db(patch_array, ips, ipe, 1, num_metgrid_levels, jps, jpe, &
                                 domain_array, ids, ide, 1, num_metgrid_levels, jds, jde)
      if(mytask.eq.0)   grid%u1=domain_array                           
      call mpi_barrier(local_communicator, mpi_ierr)
 
      ! u2
      do j=jps, jpe
       do k=1,num_metgrid_levels
        do i=ips, ipe
          patch_array(i,k,j)=grid%u2(i,k,j)
        enddo
       enddo
      enddo 
      call gather_whole_field_db(patch_array, ips, ipe, 1, num_metgrid_levels, jps, jpe, &
                                 domain_array, ids, ide, 1, num_metgrid_levels, jds, jde)                                
      if(mytask.eq.0)   grid%u2=domain_array
      call mpi_barrier(local_communicator, mpi_ierr)
      
      ! v1
      do j=jps, jpe
       do k=1,num_metgrid_levels
        do i=ips, ipe
          patch_array(i,k,j)=grid%v1(i,k,j)
        enddo
       enddo
      enddo 
      call gather_whole_field_db(patch_array, ips, ipe, 1, num_metgrid_levels, jps, jpe, &
                                 domain_array, ids, ide, 1, num_metgrid_levels, jds, jde) 
      if(mytask.eq.0)   grid%v1=domain_array
      call mpi_barrier(local_communicator, mpi_ierr)
      
      ! v2
      do j=jps, jpe
       do k=1,num_metgrid_levels
        do i=ips, ipe
          patch_array(i,k,j)=grid%v2(i,k,j)
        enddo
       enddo
      enddo 
      call gather_whole_field_db(patch_array, ips, ipe, 1, num_metgrid_levels, jps, jpe, &
                                 domain_array, ids, ide, 1, num_metgrid_levels, jds, jde)                             
      if(mytask.eq.0)   grid%v2=domain_array
      call mpi_barrier(local_communicator, mpi_ierr)
      
     ! t1
      do j=jps, jpe
       do k=1,num_metgrid_levels
        do i=ips, ipe
          patch_array(i,k,j)=grid%t1(i,k,j)
        enddo
       enddo
      enddo 
      call gather_whole_field_db(patch_array, ips, ipe, 1, num_metgrid_levels, jps, jpe, &
                                 domain_array, ids, ide, 1, num_metgrid_levels, jds, jde)
      if(mytask.eq.0)   grid%t1=domain_array
      call mpi_barrier(local_communicator, mpi_ierr)

     ! t2
      do j=jps, jpe
       do k=1,num_metgrid_levels
        do i=ips, ipe
          patch_array(i,k,j)=grid%t2(i,k,j)
        enddo
       enddo
      enddo 
      call gather_whole_field_db(patch_array, ips, ipe, 1, num_metgrid_levels, jps, jpe, &
                                 domain_array, ids, ide, 1, num_metgrid_levels, jds, jde)                                     
      if(mytask.eq.0)   grid%t2=domain_array
      call mpi_barrier(local_communicator, mpi_ierr)
      
     ! q1
      do j=jps, jpe
       do k=1,num_metgrid_levels
        do i=ips, ipe
          patch_array(i,k,j)=grid%q1(i,k,j)
        enddo
       enddo
      enddo 
      call gather_whole_field_db(patch_array, ips, ipe, 1, num_metgrid_levels, jps, jpe, &
                                 domain_array, ids, ide, 1, num_metgrid_levels, jds, jde)                                 
      if(mytask.eq.0)   grid%q1=domain_array
      call mpi_barrier(local_communicator, mpi_ierr)
      
     ! q2
      do j=jps, jpe
       do k=1,num_metgrid_levels
        do i=ips, ipe
          patch_array(i,k,j)=grid%q2(i,k,j)
        enddo
       enddo
      enddo 
      call gather_whole_field_db(patch_array, ips, ipe, 1, num_metgrid_levels, jps, jpe, &
                                 domain_array, ids, ide, 1, num_metgrid_levels, jds, jde)
      if(mytask.eq.0)   grid%q2=domain_array
      call mpi_barrier(local_communicator, mpi_ierr)

     !dudt
      do j=jps, jpe
       do k=1,num_metgrid_levels
        do i=ips, ipe
          patch_array(i,k,j)=grid%dudt(i,k,j)
        enddo
       enddo
      enddo
      call gather_whole_field_db(patch_array, ips, ipe, 1, num_metgrid_levels, jps, jpe, &
                                 domain_array, ids, ide, 1, num_metgrid_levels, jds, jde)
      if(mytask.eq.0)   grid%dudt=domain_array
      call mpi_barrier(local_communicator, mpi_ierr)

     !dvdt
      do j=jps, jpe
       do k=1,num_metgrid_levels
        do i=ips, ipe
          patch_array(i,k,j)=grid%dvdt(i,k,j)
        enddo
       enddo
      enddo
      call gather_whole_field_db(patch_array, ips, ipe, 1, num_metgrid_levels, jps, jpe, &
                                 domain_array, ids, ide, 1, num_metgrid_levels, jds, jde)
      if(mytask.eq.0)   grid%dvdt=domain_array
      call mpi_barrier(local_communicator, mpi_ierr)

     !dtdt
      do j=jps, jpe
       do k=1,num_metgrid_levels
        do i=ips, ipe
          patch_array(i,k,j)=grid%dtdt(i,k,j)
        enddo
       enddo
      enddo
      call gather_whole_field_db(patch_array, ips, ipe, 1, num_metgrid_levels, jps, jpe, &
                                 domain_array, ids, ide, 1, num_metgrid_levels, jds, jde)
      if(mytask.eq.0)   grid%dtdt=domain_array
      call mpi_barrier(local_communicator, mpi_ierr)

     !dqdt
      do j=jps, jpe
       do k=1,num_metgrid_levels
        do i=ips, ipe
          patch_array(i,k,j)=grid%dqdt(i,k,j)
        enddo
       enddo
      enddo
      call gather_whole_field_db(patch_array, ips, ipe, 1, num_metgrid_levels, jps, jpe, &
                                 domain_array, ids, ide, 1, num_metgrid_levels, jds, jde)
      if(mytask.eq.0)   grid%dqdt=domain_array
      call mpi_barrier(local_communicator, mpi_ierr)


      deallocate(patch_array)
      deallocate(domain_array)      
      
   end subroutine camgrid_gather_tend  

END MODULE module_camdomain
