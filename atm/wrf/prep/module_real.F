!  Create an initial data set for the WRF model based on real data.  This
!  program is specifically set up for the Eulerian, mass-based coordinate.
module module_real

   
   USE module_domain
   USE module_configure
   USE module_initialize_real
   USE module_optional_input   
   USE module_bdy_prep

   contains
   
   SUBROUTINE init_real(grid)
   
      IMPLICIT NONE

      TYPE(domain)                 :: grid

      INTEGER :: ids , ide , jds , jde , kds , kde
      INTEGER :: ims , ime , jms , jme , kms , kme
      INTEGER :: ips , ipe , jps , jpe , kps , kpe

         !  Various sizes that we need to be concerned about.

         ids = grid%sd31
         ide = grid%ed31
         kds = grid%sd32
         kde = grid%ed32
         jds = grid%sd33
         jde = grid%ed33

         ims = grid%sm31
         ime = grid%em31
         kms = grid%sm32
         kme = grid%em32
         jms = grid%sm33
         jme = grid%em33

         ips = grid%sp31
         ipe = grid%ep31
         kps = grid%sp32
         kpe = grid%ep32
         jps = grid%sp33
         jpe = grid%ep33
   
         !  This is the space needed to save the current 3d data for use in computing
         !  the lateral boundary tendencies.
         IF ( ALLOCATED ( ubdy3dtemp1 ) ) DEALLOCATE ( ubdy3dtemp1 )
         IF ( ALLOCATED ( vbdy3dtemp1 ) ) DEALLOCATE ( vbdy3dtemp1 )
         IF ( ALLOCATED ( tbdy3dtemp1 ) ) DEALLOCATE ( tbdy3dtemp1 )
         IF ( ALLOCATED ( pbdy3dtemp1 ) ) DEALLOCATE ( pbdy3dtemp1 )
         IF ( ALLOCATED ( qbdy3dtemp1 ) ) DEALLOCATE ( qbdy3dtemp1 )
         IF ( ALLOCATED ( mbdy2dtemp1 ) ) DEALLOCATE ( mbdy2dtemp1 )
         IF ( ALLOCATED ( ubdy3dtemp2 ) ) DEALLOCATE ( ubdy3dtemp2 )
         IF ( ALLOCATED ( vbdy3dtemp2 ) ) DEALLOCATE ( vbdy3dtemp2 )
         IF ( ALLOCATED ( tbdy3dtemp2 ) ) DEALLOCATE ( tbdy3dtemp2 )
         IF ( ALLOCATED ( pbdy3dtemp2 ) ) DEALLOCATE ( pbdy3dtemp2 )
         IF ( ALLOCATED ( qbdy3dtemp2 ) ) DEALLOCATE ( qbdy3dtemp2 )
         IF ( ALLOCATED ( mbdy2dtemp2 ) ) DEALLOCATE ( mbdy2dtemp2 )
         ALLOCATE ( ubdy3dtemp1(ims:ime,kms:kme,jms:jme) )
         ALLOCATE ( vbdy3dtemp1(ims:ime,kms:kme,jms:jme) )
         ALLOCATE ( tbdy3dtemp1(ims:ime,kms:kme,jms:jme) )
         ALLOCATE ( pbdy3dtemp1(ims:ime,kms:kme,jms:jme) )
         ALLOCATE ( qbdy3dtemp1(ims:ime,kms:kme,jms:jme) )
         ALLOCATE ( mbdy2dtemp1(ims:ime,1:1,    jms:jme) )
         ALLOCATE ( ubdy3dtemp2(ims:ime,kms:kme,jms:jme) )
         ALLOCATE ( vbdy3dtemp2(ims:ime,kms:kme,jms:jme) )
         ALLOCATE ( tbdy3dtemp2(ims:ime,kms:kme,jms:jme) )
         ALLOCATE ( pbdy3dtemp2(ims:ime,kms:kme,jms:jme) )
         ALLOCATE ( qbdy3dtemp2(ims:ime,kms:kme,jms:jme) )
         ALLOCATE ( mbdy2dtemp2(ims:ime,1:1,    jms:jme) )

         if(allocated(mu0)) deallocate(mu0)
         if(allocated(al)) deallocate(al)
         if(allocated(alb)) deallocate(alb)
         if(allocated(alt)) deallocate(alt)
         if(allocated(pdhs)) deallocate(pdhs)
         if(allocated(p3d)) deallocate(p3d)
         if(allocated(php)) deallocate(php)
         if(allocated(phb)) deallocate(phb)
         if(allocated(ph0)) deallocate(ph0)
         if(allocated(p_hyd)) deallocate(p_hyd)
         allocate(mu0(ims:ime,jms:jme))
         allocate(al(ims:ime,kms:kme,jms:jme))
         allocate(alb(ims:ime,kms:kme,jms:jme))
         allocate(alt(ims:ime,kms:kme,jms:jme))
         allocate(pdhs(ims:ime,kms:kme,jms:jme))
         allocate(p3d(ims:ime,kms:kme,jms:jme))
         allocate(php(ims:ime,kms:kme,jms:jme))
         allocate(phb(ims:ime,kms:kme,jms:jme))
         allocate(ph0(ims:ime,kms:kme,jms:jme))
         allocate(p_hyd(ims:ime,kms:kme,jms:jme))
         allocate(t_init(ims:ime,kms:kme,jms:jme))
         allocate(utemp(ims:ime,kms:kme,jms:jme))
         allocate(vtemp(ims:ime,kms:kme,jms:jme))

   end subroutine init_real      

   SUBROUTINE final_real
   
      IMPLICIT NONE   
         IF ( ALLOCATED ( ubdy3dtemp1 ) ) DEALLOCATE ( ubdy3dtemp1 )
         IF ( ALLOCATED ( vbdy3dtemp1 ) ) DEALLOCATE ( vbdy3dtemp1 )
         IF ( ALLOCATED ( tbdy3dtemp1 ) ) DEALLOCATE ( tbdy3dtemp1 )
         IF ( ALLOCATED ( pbdy3dtemp1 ) ) DEALLOCATE ( pbdy3dtemp1 )
         IF ( ALLOCATED ( qbdy3dtemp1 ) ) DEALLOCATE ( qbdy3dtemp1 )
         IF ( ALLOCATED ( mbdy2dtemp1 ) ) DEALLOCATE ( mbdy2dtemp1 )
         IF ( ALLOCATED ( ubdy3dtemp2 ) ) DEALLOCATE ( ubdy3dtemp2 )
         IF ( ALLOCATED ( vbdy3dtemp2 ) ) DEALLOCATE ( vbdy3dtemp2 )
         IF ( ALLOCATED ( tbdy3dtemp2 ) ) DEALLOCATE ( tbdy3dtemp2 )
         IF ( ALLOCATED ( pbdy3dtemp2 ) ) DEALLOCATE ( pbdy3dtemp2 )
         IF ( ALLOCATED ( qbdy3dtemp2 ) ) DEALLOCATE ( qbdy3dtemp2 )
         IF ( ALLOCATED ( mbdy2dtemp2 ) ) DEALLOCATE ( mbdy2dtemp2 )
         if(allocated(mu0)) deallocate(mu0)
         if(allocated(al)) deallocate(al)
         if(allocated(alb)) deallocate(alb)
         if(allocated(alt)) deallocate(alt)
         if(allocated(pdhs)) deallocate(pdhs)
         if(allocated(p3d)) deallocate(p3d)
         if(allocated(php)) deallocate(php)
         if(allocated(phb)) deallocate(phb)
         if(allocated(ph0)) deallocate(ph0)
         if(allocated(p_hyd)) deallocate(p_hyd)
         if(allocated(t_init)) deallocate(t_init)
         if(allocated(utemp)) deallocate(utemp)
         if(allocated(vtemp)) deallocate(vtemp)
 
   end subroutine final_real   

   subroutine restart_real(grid,loop)
   
   IMPLICIT NONE
   
   TYPE(domain)                 :: grid
   INTEGER , INTENT(IN)         :: loop

   ! prepare optional setting
   flag_qc=0  
   flag_qr=0  
   flag_qi=0  
   flag_qs=0  
   flag_qg=0  
   flag_qni=0  
   flag_qv=0  
   flag_metgrid=1  
   flag_mf_xy=1  
   grid%flag_metgrid=flag_metgrid
   grid%flag_mf_xy=flag_mf_xy
   flag_sst=1
   flag_tsk=1
   flag_tavgsfc=0
   flag_snowh=0
   flag_snow=0
   flag_psfc=1
   flag_soilhgt=1
   flag_toposoil=0
   flag_slp=1

   ! set loop 
   internal_time_loop=loop
  
   ! vertical interpolation
   call bdy_prep(grid)

   ! generate lateral boundary 
   call restart_bdy(grid, loop)

   end subroutine restart_real

   
   SUBROUTINE real_interp(grid, config_flags, loop)

   IMPLICIT NONE


   TYPE(domain) :: grid
   TYPE (grid_config_rec_type)   :: config_flags
   INTEGER ::  loop
   
   CALL       wrf_debug ( 100 , 'calling med_sidata_input' )
   CALL med_sidata_input ( grid , config_flags, loop )
   CALL       wrf_debug ( 100 , 'backfrom med_sidata_input' )

   END subroutine real_interp

   SUBROUTINE med_sidata_input ( grid , config_flags, loop )

   IMPLICIT NONE

  ! Arguments
   TYPE(domain)                :: grid
   TYPE (grid_config_rec_type) :: config_flags
   INTEGER ::  loop

   INTEGER :: idum1, idum2 
   integer::idsi

   CALL set_scalar_indices_from_config ( grid%id , idum1, idum2 )   
   CALL model_to_grid_config_rec ( grid%id , model_config_rec , config_flags )
   internal_time_loop = loop
   ! since bdy_prep doesn't need st_input, sm_input, sw_input, and any surface and soil
   ! information, and flags were also set at loop = 1, so it doesn't need otional_input when loop>1
   if(loop.eq.1) then
   already_been_here = .FALSE.
   CALL wrf_debug ( 100 , 'med_sidata_input: calling init_module_optional_input' )
   CALL init_module_optional_input ( grid , config_flags )
   CALL wrf_debug ( 100 , 'med_sidata_input: calling optional_input' )
   idsi=0 ! can be any value since it's unneccesary in optional_input. Juanxiong he
   CALL optional_input ( grid , idsi , config_flags )   
   endif

   CALL  wrf_debug ( 100 , 'med_sidata_input: calling init_domain' )
   if( loop .eq. 1) then
    CALL init_domain ( grid )
    IF ( ALLOCATED ( st_input ) ) DEALLOCATE ( st_input )
    IF ( ALLOCATED ( sm_input ) ) DEALLOCATE ( sm_input )
    IF ( ALLOCATED ( sw_input ) ) DEALLOCATE ( sw_input )
   else
    call bdy_prep(grid)
   endif 

   if(grid%id.eq.1) then   ! only head grid needs lateral boundary
    CALL assemble_output ( grid , loop  )               
   endif
   
   END SUBROUTINE med_sidata_input

end module module_real
