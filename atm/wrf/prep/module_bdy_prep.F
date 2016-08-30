!REAL:MODEL_LAYER:INITIALIZATION
	
!  This MODULE holds the routines which are used to perform various initializations
!  for the individual domains, specifically for the Eulerian, mass-based coordinate.

!-----------------------------------------------------------------------

MODULE module_bdy_prep
   USE module_bc
   USE module_configure
   USE module_domain
   USE module_model_constants
   USE module_state_description
   USE module_timing
   USE module_date_time
#ifdef DM_PARALLEL
   USE module_dm
   USE module_comm_dm, ONLY : HALO_EM_VINTERP_UV_TEMP_1_sub, HALO_EM_INIT_TEMP_1_sub,&
                              HALO_EM_INIT_TEMP_2_sub,HALO_EM_INIT_TEMP_3_sub,&
                              HALO_EM_INIT_TEMP_4_sub,HALO_EM_INIT_TEMP_5_sub
#endif
   USE module_initialize_real
   USE parallel_module

   REAL , DIMENSION(:,:,:) , ALLOCATABLE :: ubdy3dtemp1 , vbdy3dtemp1 , tbdy3dtemp1 , pbdy3dtemp1 , qbdy3dtemp1
   REAL , DIMENSION(:,:,:) , ALLOCATABLE :: mbdy2dtemp1
   REAL , DIMENSION(:,:,:) , ALLOCATABLE :: ubdy3dtemp2 , vbdy3dtemp2 , tbdy3dtemp2 , pbdy3dtemp2 , qbdy3dtemp2
   REAL , DIMENSION(:,:,:) , ALLOCATABLE :: mbdy2dtemp2
   real , dimension(:,:), allocatable :: mu0
   real , dimension(:,:,:), allocatable :: al, alb, alt, &
                                           pdhs, p3d, phb, php, ph0, p_hyd, t_init,&
                                           utemp, vtemp 
                                           
CONTAINS

!-------------------------------------------------------------------

   SUBROUTINE bdy_prep ( grid )

      IMPLICIT NONE

      !  Input space and data.  No gridded meteorological data has been stored, though.

      TYPE (domain)          :: grid

      !  Local data.

      INTEGER :: idum1, idum2

      CALL bdy_prep_othertime( grid )

   END SUBROUTINE bdy_prep

!-------------------------------------------------------------------
   Subroutine prep_surface_othertime(grid)

      IMPLICIT NONE

      !  Input space and data.  No gridded meteorological data has been stored,
      !  though.
      TYPE (domain)          :: grid

      !  Local domain indices and counters.
      INTEGER :: ids, ide, jds, jde, kds, kde, &
                 ims, ime, jms, jme, kms, kme, &
                 its, ite, jts, jte, kts, kte, &
                 ips, ipe, jps, jpe, kps, kpe, &
                 i, j, k

      INTEGER :: imsx, imex, jmsx, jmex, kmsx, kmex,    &
                 ipsx, ipex, jpsx, jpex, kpsx, kpex,    &
                 imsy, imey, jmsy, jmey, kmsy, kmey,    &
                 ipsy, ipey, jpsy, jpey, kpsy, kpey

      CHARACTER (LEN=19) :: current_timestr

      !  Dimension information stored in grid data structure.
      CALL get_ijk_from_grid (  grid ,                   &
                                ids, ide, jds, jde, kds, kde,    &
                                ims, ime, jms, jme, kms, kme,    &
                                ips, ipe, jps, jpe, kps, kpe,    &
                                imsx, imex, jmsx, jmex, kmsx, kmex,    &
                                ipsx, ipex, jpsx, jpex, kpsx, kpex,    &
                                imsy, imey, jmsy, jmey, kmsy, kmey,    &
                                ipsy, ipey, jpsy, jpey, kpsy, kpey )
      its = ips ; ite = ipe ; jts = jps ; jte = jpe ; kts = kps ; kte = kpe

         !  Get the monthly values interpolated to the current date for the
         !  traditional monthly
         !  fields of green-ness fraction and background albedo.
         call domain_clock_get( grid, current_timestr = current_timestr )
         current_date = trim(current_timestr)//'.0000'
         CALL monthly_interp_to_date ( grid%greenfrac , current_date , grid%vegfra , &
                                       ids , ide , jds , jde , kds , kde , &
                                       ims , ime , jms , jme , kms , kme , &
                                       its , ite , jts , jte , kts , kte )

         CALL monthly_interp_to_date ( grid%albedo12m , current_date , grid%albbck , &
                                       ids , ide , jds , jde , kds , kde , &
                                       ims , ime , jms , jme , kms , kme , &
                                       its , ite , jts , jte , kts , kte )

         !  Get the min/max of each i,j for the monthly green-ness fraction.
         CALL monthly_min_max ( grid%greenfrac , grid%shdmin , grid%shdmax , &
                                ids , ide , jds , jde , kds , kde , &
                                ims , ime , jms , jme , kms , kme , &
                                its , ite , jts , jte , kts , kte )

         !  The model expects the green-ness values in percent, not fraction.
         DO j = jts, MIN(jte,jde-1)
           DO i = its, MIN(ite,ide-1)
              grid%vegfra(i,j) = grid%vegfra(i,j) * 100.
              grid%shdmax(i,j) = grid%shdmax(i,j) * 100.
              grid%shdmin(i,j) = grid%shdmin(i,j) * 100.
           END DO
         END DO

         !  The model expects the albedo fields as a fraction, not a percent.
         !  Set the
         !  water values to 8%.
         DO j = jts, MIN(jte,jde-1)
           DO i = its, MIN(ite,ide-1)
              grid%albbck(i,j) = grid%albbck(i,j) / 100.
              grid%snoalb(i,j) = grid%snoalb(i,j) / 100.
              IF ( grid%landmask(i,j) .LT. 0.5 ) THEN
                 grid%albbck(i,j) = 0.08
                 grid%snoalb(i,j) = 0.08
              END IF
           END DO
         END DO

   end subroutine prep_surface_othertime
!-------------------------------------------------------------------

   SUBROUTINE bdy_prep_othertime ( grid )

      USE module_optional_input
      IMPLICIT NONE

      !  Input space and data.  No gridded meteorological data has been stored, though.

      TYPE (domain)          :: grid


      TYPE (grid_config_rec_type)              :: config_flags

      !  Local domain indices and counters.

      INTEGER :: num_veg_cat , num_soil_top_cat , num_soil_bot_cat
      INTEGER :: loop , num_seaice_changes

      INTEGER :: ids, ide, jds, jde, kds, kde, &
                 ims, ime, jms, jme, kms, kme, &
                 its, ite, jts, jte, kts, kte, &
                 ips, ipe, jps, jpe, kps, kpe, &
                 i, j, k

      INTEGER :: imsx, imex, jmsx, jmex, kmsx, kmex,    &
                 ipsx, ipex, jpsx, jpex, kpsx, kpex,    &
                 imsy, imey, jmsy, jmey, kmsy, kmey,    &
                 ipsy, ipey, jpsy, jpey, kpsy, kpey

      INTEGER :: ns

      !  Local data

      INTEGER :: error
      INTEGER :: im, num_3d_m, num_3d_s
      REAL    :: p_surf, p_level
      REAL    :: cof1, cof2
      REAL    :: qvf , qvf1 , qvf2 , pd_surf
      REAL    :: p00 , t00 , a , tiso
      REAL    :: hold_znw
      LOGICAL :: were_bad

      LOGICAL :: stretch_grid, dry_sounding, debug
      INTEGER IICOUNT

      REAL :: p_top_requested , temp
      INTEGER :: num_metgrid_levels
      REAL , DIMENSION(max_eta) :: eta_levels
      REAL :: max_dz

      REAL    :: zap_close_levels
      INTEGER :: force_sfc_in_vinterp
      INTEGER :: interp_type , lagrange_order , extrap_type , t_extrap_type
      LOGICAL :: lowest_lev_from_sfc , use_levels_below_ground , use_surface
      LOGICAL :: we_have_tavgsfc , we_have_tsk

      INTEGER :: lev500 , loop_count
      REAL    :: zl , zu , pl , pu , z500 , dz500 , tvsfc , dpmu

      LOGICAL , PARAMETER :: want_full_levels = .TRUE.
      LOGICAL , PARAMETER :: want_half_levels = .FALSE.

      CHARACTER (LEN=80) :: a_message
      REAL :: max_mf

!-- Carsel and Parrish [1988]
      REAL , DIMENSION(100) :: lqmi

      REAL :: pexp1,pexp2
      
      CHARACTER (LEN=19) :: current_timestr

      !  Dimension information stored in grid data structure.

      CALL get_ijk_from_grid (  grid ,                   &
                                ids, ide, jds, jde, kds, kde,    &
                                ims, ime, jms, jme, kms, kme,    &
                                ips, ipe, jps, jpe, kps, kpe,    &
                                imsx, imex, jmsx, jmex, kmsx, kmex,    &
                                ipsx, ipex, jpsx, jpex, kpsx, kpex,    &
                                imsy, imey, jmsy, jmey, kmsy, kmey,    &
                                ipsy, ipey, jpsy, jpey, kpsy, kpey )
      its = ips ; ite = ipe ; jts = jps ; jte = jpe ; kts = kps ; kte = kpe

      mu0 = 0.0
      al = 0.0
      alb = 0.0
      alt = 0.0
      pdhs = 0.0
      p3d =0.0
      phb = 0.0
      php = 0.0
      ph0 = 0.0
      p_hyd = 0.0
      t_init = 0.0        
      utemp=0.0
      vtemp=0.0
      
      CALL model_to_grid_config_rec ( grid%id , model_config_rec , config_flags )

      IF ( ( grid%id .EQ. 1 ) .AND. ( .NOT. config_flags%map_proj .EQ. PROJ_CASSINI ) ) THEN
         max_mf = grid%msft(its,jts)
         DO j=jts,MIN(jde-1,jte)
            DO i=its,MIN(ide-1,ite)
               max_mf = MAX ( max_mf , grid%msft(i,j) )
            END DO
         END DO
#if ( defined(DM_PARALLEL)  &&   ! defined(STUBMPI) )
         max_mf = wrf_dm_max_real ( max_mf )
#endif
         WRITE ( a_message , FMT='(A,F5.2,A)' ) 'Max map factor in domain 1 = ',max_mf,'. Scale the dt in the model accordingly.'
         CALL wrf_message ( a_message ) 
      END IF

      !  Check to see if the boundary conditions are set properly in the namelist file.
      !  This checks for sufficiency and redundancy.

      CALL boundary_condition_check( config_flags, bdyzone, error, grid%id )

      !  Pull in the info in the namelist to compare it to the input data.

      grid%real_data_init_type = model_config_rec%real_data_init_type

      !  To define the base state, we call a USER MODIFIED routine to set the three
      !  necessary constants:  p00 (sea level pressure, Pa), t00 (sea level temperature, K),
      !  and A (temperature difference, from 1000 mb to 300 mb, K).
   
      CALL const_module_initialize ( p00 , t00 , a , tiso ) 

      !  Save these constants to write out in model output file

      grid%t00  = t00
      grid%p00  = p00
      grid%tlp  = a
      grid%tiso = tiso
      
      !  Is there any vertical interpolation to do?  The "old" data comes in on the correct
      !  vertical locations already.

      IF ( flag_metgrid .EQ. 1 ) THEN  !   <----- START OF VERTICAL INTERPOLATION PART ---->

        !  If we have any input low-res surface pressure, we store it.

         IF ( flag_psfc .EQ. 1 ) THEN
            DO j = jts, MIN(jte,jde-1)
              DO i = its, MIN(ite,ide-1)
                 grid%p_gc(i,1,j) = grid%psfc_gc(i,j)
              END DO
            END DO
         END IF
        
         !  The number of vertical levels in the input data.  There is no staggering for
         !  different variables.

         num_metgrid_levels = grid%num_metgrid_levels

         !  Check for and semi-fix missing surface fields.

         IF ( grid%p_gc(its,num_metgrid_levels,jts) .LT. grid%p_gc(its,2,jts) ) THEN
            k = 2
         ELSE
            k = num_metgrid_levels
         END IF

         IF ( grid%t_gc(its,1,jts) .EQ. -1.E30 ) THEN
            DO j = jts, MIN(jte,jde-1)
               DO i = its, MIN(ite,ide-1)
                  grid%t_gc(i,1,j) = grid%t_gc(i,k,j)
               END DO
            END DO
            config_flags%use_surface = .FALSE.
            grid%use_surface = .FALSE.
            WRITE ( a_message , * ) 'Missing surface temp, replaced with closest level, use_surface set to false.'
            CALL wrf_message ( a_message ) 
         END IF

         IF ( grid%rh_gc(its,1,jts) .EQ. -1.E30 ) THEN
            DO j = jts, MIN(jte,jde-1)
               DO i = its, MIN(ite,ide-1)
                  grid%rh_gc(i,1,j) = grid%rh_gc(i,k,j)
               END DO
            END DO
            config_flags%use_surface = .FALSE.
            grid%use_surface = .FALSE.
            WRITE ( a_message , * ) 'Missing surface RH, replaced with closest level, use_surface set to false.'
            CALL wrf_message ( a_message ) 
         END IF

         IF ( grid%u_gc(its,1,jts) .EQ. -1.E30 ) THEN
            DO j = jts, MIN(jte,jde-1)
               DO i = its, ite
                  grid%u_gc(i,1,j) = grid%u_gc(i,k,j)
               END DO
            END DO
            config_flags%use_surface = .FALSE.
            grid%use_surface = .FALSE.
            WRITE ( a_message , * ) 'Missing surface u wind, replaced with closest level, use_surface set to false.'
            CALL wrf_message ( a_message ) 
         END IF

         IF ( grid%v_gc(its,1,jts) .EQ. -1.E30 ) THEN
            DO j = jts, jte
               DO i = its, MIN(ite,ide-1)
                  grid%v_gc(i,1,j) = grid%v_gc(i,k,j)
               END DO
            END DO
            config_flags%use_surface = .FALSE.
            grid%use_surface = .FALSE.
            WRITE ( a_message , * ) 'Missing surface v wind, replaced with closest level, use_surface set to false.'
            CALL wrf_message ( a_message ) 
         END IF

         !  Compute the mixing ratio from the input relative humidity.

         IF ( ( flag_qv .NE. 1 ) .AND. ( flag_sh .NE. 1 ) ) THEN
            IF ( grid%p_gc(its,num_metgrid_levels,jts) .LT. grid%p_gc(its,2,jts) ) THEN
               k = 2
            ELSE
               k = num_metgrid_levels
            END IF

            CALL rh_to_mxrat (grid%rh_gc, grid%t_gc, grid%p_gc, grid%qv_gc ,         &
                              config_flags%rh2qv_wrt_liquid ,                        &
                              config_flags%qv_max_p_safe ,                           &
                              config_flags%qv_max_flag , config_flags%qv_max_value , &
                              config_flags%qv_min_p_safe ,                           &
                              config_flags%qv_min_flag , config_flags%qv_min_value , &
                              ids , ide , jds , jde , 1   , num_metgrid_levels ,     &
                              ims , ime , jms , jme , 1   , num_metgrid_levels ,     &
                              its , ite , jts , jte , 1   , num_metgrid_levels )

         ELSE IF ( flag_sh .EQ. 1 ) THEN
            IF ( grid%p_gc(its,num_metgrid_levels,jts) .LT. grid%p_gc(its,2,jts) ) THEN
               k = 2
            ELSE
               k = num_metgrid_levels
            END IF
            IF ( grid%sh_gc(its,kts,jts) .LT. 1.e-6 ) THEN
               DO j = jts, MIN(jte,jde-1)
                  DO i = its, MIN(ite,ide-1)
                     grid%sh_gc(i,1,j) = grid%sh_gc(i,k,j)
                  END DO
               END DO
            END IF

         END IF

         !  Some data sets do not provide a 3d geopotential height field.  

         IF ( grid%ght_gc(its,grid%num_metgrid_levels/2,jts) .LT. 1 ) THEN
            DO j = jts, MIN(jte,jde-1)
               DO k = kts+1 , grid%num_metgrid_levels
                  DO i = its, MIN(ite,ide-1)
                     grid%ght_gc(i,k,j) = grid%ght_gc(i,k-1,j) - &
                        R_d / g * 0.5 * ( grid%t_gc(i,k  ,j) * ( 1 + 0.608 * grid%qv_gc(i,k  ,j) ) +   &
                                          grid%t_gc(i,k-1,j) * ( 1 + 0.608 * grid%qv_gc(i,k-1,j) ) ) * &
                        LOG ( grid%p_gc(i,k,j) / grid%p_gc(i,k-1,j) )
                  END DO
               END DO
            END DO
         END IF


         !  The requested ptop for real data cases.

         p_top_requested = grid%p_top_requested

         grid%p_top = p_top_requested ! for the restart needs p_top_save
         
         p_top_save = grid%p_top

            CALL find_p_top ( grid%p_gc , grid%p_top , &
                              ids , ide , jds , jde , 1   , num_metgrid_levels , &
                              ims , ime , jms , jme , 1   , num_metgrid_levels , &
                              its , ite , jts , jte , 1   , num_metgrid_levels )

#if ( defined( DM_PARALLEL ) && ( ! defined( STUBMPI ) ) )
            grid%p_top = wrf_dm_max_real ( grid%p_top )
#endif
            IF ( grid%p_top .GT. p_top_save ) THEN
               print *,'grid%p_top from last time period = ',p_top_save
               print *,'grid%p_top from this time period = ',grid%p_top
               CALL wrf_error_fatal ( 'grid%p_top > previous value' )
            END IF
            grid%p_top = p_top_save

         !  Get the monthly values interpolated to the current date for the
         !  traditional monthly
         !  fields of green-ness fraction and background albedo.
         call domain_clock_get( grid, current_timestr = current_timestr )
         current_date = trim(current_timestr)//'.0000' 
         CALL monthly_interp_to_date ( grid%greenfrac , current_date , grid%vegfra , &
                                       ids , ide , jds , jde , kds , kde , &
                                       ims , ime , jms , jme , kms , kme , &
                                       its , ite , jts , jte , kts , kte )

         CALL monthly_interp_to_date ( grid%albedo12m , current_date , grid%albbck , &
                                       ids , ide , jds , jde , kds , kde , &
                                       ims , ime , jms , jme , kms , kme , &
                                       its , ite , jts , jte , kts , kte )

         !  Get the min/max of each i,j for the monthly green-ness fraction.
         CALL monthly_min_max ( grid%greenfrac , grid%shdmin , grid%shdmax , &
                                ids , ide , jds , jde , kds , kde , &
                                ims , ime , jms , jme , kms , kme , &
                                its , ite , jts , jte , kts , kte )

         !  The model expects the green-ness values in percent, not fraction.
         DO j = jts, MIN(jte,jde-1)
           DO i = its, MIN(ite,ide-1)
              grid%vegfra(i,j) = grid%vegfra(i,j) * 100.
              grid%shdmax(i,j) = grid%shdmax(i,j) * 100.
              grid%shdmin(i,j) = grid%shdmin(i,j) * 100.
           END DO
         END DO

         !  The model expects the albedo fields as a fraction, not a percent.
         !  Set the
         !  ater values to 8%.
         DO j = jts, MIN(jte,jde-1)
           DO i = its, MIN(ite,ide-1)
              grid%albbck(i,j) = grid%albbck(i,j) / 100.
              grid%snoalb(i,j) = grid%snoalb(i,j) / 100.
              IF ( grid%landmask(i,j) .LT. 0.5 ) THEN
                 grid%albbck(i,j) = 0.08
                 grid%snoalb(i,j) = 0.08
              END IF
           END DO
         END DO
  
         !  If the pressure levels in the middle of the atmosphere are upside down, then
         !  this is hybrid data.  Computing the new surface pressure should use sfcprs2.

         IF ( grid%p_gc(its,num_metgrid_levels/2,jts) .LT. grid%p_gc(its,num_metgrid_levels/2+1,jts) ) THEN
            config_flags%sfcp_to_sfcp = .TRUE.
         END IF

         !  Two ways to get the surface pressure.  1) If we have the low-res input surface
         !  pressure and the low-res topography, then we can do a simple hydrostatic
         !  relation.  2) Otherwise we compute the surface pressure from the sea-level
         !  pressure.
         !  Note that on output, grid%psfc is now hi-res.  The low-res surface pressure and
         !  elevation are grid%psfc_gc and grid%ht_gc (same as grid%ght_gc(k=1)).

         IF      ( ( flag_psfc    .EQ. 1 ) .AND. &
                   ( flag_soilhgt .EQ. 1 ) .AND. &
                   ( flag_slp     .EQ. 1 ) .AND. &
                   ( .NOT. config_flags%sfcp_to_sfcp ) ) THEN
            CALL sfcprs3(grid%ght_gc, grid%p_gc, grid%ht, &
                         grid%pslv_gc, grid%psfc_gc, &
                         ids , ide , jds , jde , 1   , num_metgrid_levels , &
                         ims , ime , jms , jme , 1   , num_metgrid_levels , &
                         its , ite , jts , jte , 1   , num_metgrid_levels )
         ELSE IF ( ( flag_psfc    .EQ. 1 ) .AND. &
                   ( flag_soilhgt .EQ. 1 ) .AND. &
                   ( config_flags%sfcp_to_sfcp ) ) THEN
            CALL sfcprs2(grid%t_gc, grid%qv_gc, grid%ght_gc, grid%psfc_gc, grid%ht, &
                         grid%tavgsfc, grid%p_gc, grid%psfc_gc, we_have_tavgsfc, &
                         ids , ide , jds , jde , 1   , num_metgrid_levels , &
                         ims , ime , jms , jme , 1   , num_metgrid_levels , &
                         its , ite , jts , jte , 1   , num_metgrid_levels )
         ELSE IF ( flag_slp     .EQ. 1 ) THEN
            CALL sfcprs (grid%t_gc, grid%qv_gc, grid%ght_gc, grid%pslv_gc, grid%ht, &
                         grid%tavgsfc, grid%p_gc, grid%psfc_gc, we_have_tavgsfc, &
                         ids , ide , jds , jde , 1   , num_metgrid_levels , &
                         ims , ime , jms , jme , 1   , num_metgrid_levels , &
                         its , ite , jts , jte , 1   , num_metgrid_levels )
         ELSE
            WRITE(a_message,FMT='(3(A,I2),A,L1)') 'ERROR in psfc: flag_psfc = ',flag_psfc, &
                                               ', flag_soilhgt = ',flag_soilhgt , &
                                               ', flag_slp = ',flag_slp , & 
                                               ', sfcp_to_sfcp = ',config_flags%sfcp_to_sfcp 
            CALL wrf_message ( a_message ) 
            CALL wrf_error_fatal ( 'not enough info for a p sfc computation' )
         END IF

         !  If we have no input surface pressure, we'd better stick something in there.

         IF ( flag_psfc .NE. 1 ) THEN
            DO j = jts, MIN(jte,jde-1)
              DO i = its, MIN(ite,ide-1)
                 grid%p_gc(i,1,j) = grid%psfc_gc(i,j)
              END DO
            END DO
         END IF

         !  Integrate the mixing ratio to get the vapor pressure.
         CALL integ_moist ( grid%qv_gc , grid%p_gc , grid%pd_gc , grid%t_gc , grid%ght_gc , grid%intq_gc , &
                            ids , ide , jds , jde , 1   , num_metgrid_levels , &
                            ims , ime , jms , jme , 1   , num_metgrid_levels , &
                            its , ite , jts , jte , 1   , num_metgrid_levels )

         !  Compute the difference between the dry, total surface pressure (input) and the
         !  dry top pressure (constant).

         CALL p_dts ( mu0 , grid%intq_gc , grid%psfc_gc , grid%p_top , &
                      ids , ide , jds , jde , 1   , num_metgrid_levels , &
                      ims , ime , jms , jme , 1   , num_metgrid_levels , &
                      its , ite , jts , jte , 1   , num_metgrid_levels )

         !  Compute the dry, hydrostatic surface pressure.

         CALL p_dhs ( pdhs , grid%ht , p00 , t00 , a , &
                      ids , ide , jds , jde , kds , kde , &
                      ims , ime , jms , jme , kms , kme , &
                      its , ite , jts , jte , kts , kte )

         !  Compute the eta levels if not defined already.

         IF ( grid%znw(1) .NE. 1.0 ) THEN

            eta_levels(1:kde) = model_config_rec%eta_levels(1:kde)
            max_dz            = model_config_rec%max_dz

            CALL compute_eta ( grid%znw , &
                               eta_levels , max_eta , max_dz , &
                               grid%p_top , g , p00 , cvpm , a , r_d , cp , t00 , p1000mb , t0 , tiso , &
                               ids , ide , jds , jde , kds , kde , &
                               ims , ime , jms , jme , kms , kme , &
                               its , ite , jts , jte , kts , kte )
         END IF

         !  The input field is temperature, we want potential temp.

         CALL t_to_theta ( grid%t_gc , grid%p_gc , p00 , &
                           ids , ide , jds , jde , 1   , num_metgrid_levels , &
                           ims , ime , jms , jme , 1   , num_metgrid_levels , &
                           its , ite , jts , jte , 1   , num_metgrid_levels )

         IF ( flag_slp .EQ. 1 ) THEN

            !  On the eta surfaces, compute the dry pressure = mu eta, stored in
            !  grid%pbtemp, since it is a pressure, and we don't need another kms:kme 3d
            !  array floating around.  The grid%pbtemp array is re-computed as the base pressure
            !  later after the vertical interpolations are complete.

            CALL p_dry ( mu0 , grid%znw , grid%p_top , grid%pbtemp , want_full_levels , &
                         ids , ide , jds , jde , kds , kde , &
                         ims , ime , jms , jme , kms , kme , &
                         its , ite , jts , jte , kts , kte )

            !  All of the vertical interpolations are done in dry-pressure space.  The
            !  input data has had the moisture removed (grid%pd_gc).  The target levels (grid%pbtemp)
            !  had the vapor pressure removed from the surface pressure, then they were
            !  scaled by the eta levels.

         interp_type = 2
         lagrange_order = grid%lagrange_order
         lowest_lev_from_sfc = .FALSE.
         use_levels_below_ground = .TRUE.
         use_surface = .TRUE.
         zap_close_levels = grid%zap_close_levels
         force_sfc_in_vinterp = 0
         t_extrap_type = grid%t_extrap_type
         extrap_type = 1

            !  For the height field, the lowest level pressure is the slp (approximately "dry").  The
            !  lowest level of the input height field (to be associated with slp) then is an array
            !  of zeros.

            DO j = jts, MIN(jte,jde-1)
               DO i = its, MIN(ite,ide-1)
                  grid%psfc_gc(i,j) = grid%pd_gc(i,1,j)
                  grid%pd_gc(i,1,j) = grid%pslv_gc(i,j) - ( grid%p_gc(i,1,j) - grid%pd_gc(i,1,j) )
                  grid%ht_gc(i,j) = grid%ght_gc(i,1,j)
                  grid%ght_gc(i,1,j) = 0.
               END DO
            END DO

            CALL vert_interp ( grid%ght_gc , grid%pd_gc , ph0 , grid%pbtemp , &
                               num_metgrid_levels , 'Z' , &
                               interp_type , lagrange_order , extrap_type , &
                               lowest_lev_from_sfc , use_levels_below_ground , use_surface , &
                               zap_close_levels , force_sfc_in_vinterp , &
                               ids , ide , jds , jde , kds , kde , &
                               ims , ime , jms , jme , kms , kme , &
                               its , ite , jts , jte , kts , kte )

            !  Put things back to normal.

            DO j = jts, MIN(jte,jde-1)
               DO i = its, MIN(ite,ide-1)
                  grid%pd_gc(i,1,j) = grid%psfc_gc(i,j)
                  grid%ght_gc(i,1,j) = grid%ht_gc(i,j)
               END DO
            END DO

         END IF

         !  Now the rest of the variables on half-levels to inteprolate.

         CALL p_dry ( mu0 , grid%znw , grid%p_top , grid%pbtemp , want_half_levels , &
                      ids , ide , jds , jde , kds , kde , &
                      ims , ime , jms , jme , kms , kme , &
                      its , ite , jts , jte , kts , kte )

         interp_type = grid%interp_type
         lagrange_order = grid%lagrange_order
         lowest_lev_from_sfc = grid%lowest_lev_from_sfc
         use_levels_below_ground = grid%use_levels_below_ground
         use_surface = grid%use_surface
         zap_close_levels = grid%zap_close_levels
         force_sfc_in_vinterp = grid%force_sfc_in_vinterp
         t_extrap_type = grid%t_extrap_type
         extrap_type = grid%extrap_type

         !  Interpolate RH, diagnose Qv later when have temp and pressure.  Temporarily
         !  store this in the utemp space, for later diagnosis into Qv and stored into moist.

         CALL vert_interp ( grid%rh_gc , grid%pd_gc , utemp , grid%pbtemp , &
                            num_metgrid_levels , 'Q' , &
                            interp_type , lagrange_order , extrap_type , &
                            lowest_lev_from_sfc , use_levels_below_ground , use_surface , &
                            zap_close_levels , force_sfc_in_vinterp , &
                            ids , ide , jds , jde , kds , kde , &
                            ims , ime , jms , jme , kms , kme , &
                            its , ite , jts , jte , kts , kte )

         CALL vert_interp ( grid%t_gc , grid%pd_gc , grid%ttemp      , grid%pbtemp , &
                            num_metgrid_levels , 'T' , &
                            interp_type , lagrange_order , t_extrap_type , &
                            lowest_lev_from_sfc , use_levels_below_ground , use_surface , &
                            zap_close_levels , force_sfc_in_vinterp , &
                            ids , ide , jds , jde , kds , kde , &
                            ims , ime , jms , jme , kms , kme , &
                            its , ite , jts , jte , kts , kte )
     
         !  Do not have full pressure on eta levels, get a first guess at Qv by using
         !  dry pressure.  The use of utemp (rh) and vtemp (temperature) is temporary.
         !  We fix the approximation to Qv after the total pressure is available on
         !  eta surfaces.

         vtemp = grid%ttemp
         CALL theta_to_t ( vtemp , grid%pbtemp , p00 , &
                           ids , ide , jds , jde , kds , kde , &
                           ims , ime , jms , jme , kms , kme , &
                           its , ite , jts , jte , kts , kte )

         CALL rh_to_mxrat (utemp, vtemp, grid%pbtemp, grid%moisttemp,       &
                           config_flags%rh2qv_wrt_liquid ,                        &
                           config_flags%qv_max_p_safe ,                           &
                           config_flags%qv_max_flag , config_flags%qv_max_value , &
                           config_flags%qv_min_p_safe ,                           &
                           config_flags%qv_min_flag , config_flags%qv_min_value , &
                           ids , ide , jds , jde , kds , kde ,                    &
                           ims , ime , jms , jme , kms , kme ,                    &
                           its , ite , jts , jte , kts , kte-1 )

#ifdef DM_PARALLEL
         ips = its ; ipe = ite ; jps = jts ; jpe = jte ; kps = kts ; kpe = kte

         !  For the U and V vertical interpolation, we need the pressure defined
         !  at both the locations for the horizontal momentum, which we get by
         !  averaging two pressure values (i and i-1 for U, j and j-1 for V).  The
         !  pressure field on input (grid%pd_gc) and the pressure of the new coordinate
         !  (grid%pbtemp) are both communicated with an 8 stencil.

#   include "HALO_EM_VINTERP_UV_TEMP_1.inc"
#endif

         CALL vert_interp ( grid%u_gc , grid%pd_gc , grid%utemp          , grid%pbtemp , &
                            num_metgrid_levels , 'U' , &
                            interp_type , lagrange_order , extrap_type , &
                            lowest_lev_from_sfc , use_levels_below_ground , use_surface , &
                            zap_close_levels , force_sfc_in_vinterp , &
                            ids , ide , jds , jde , kds , kde , &
                            ims , ime , jms , jme , kms , kme , &
                            its , ite , jts , jte , kts , kte )

         CALL vert_interp ( grid%v_gc , grid%pd_gc , grid%vtemp           , grid%pbtemp , &
                            num_metgrid_levels , 'V' , &
                            interp_type , lagrange_order , extrap_type , &
                            lowest_lev_from_sfc , use_levels_below_ground , use_surface , &
                            zap_close_levels , force_sfc_in_vinterp , &
                            ids , ide , jds , jde , kds , kde , &
                            ims , ime , jms , jme , kms , kme , &
                            its , ite , jts , jte , kts , kte )

      END IF     !   <----- END OF VERTICAL INTERPOLATION PART ---->

      !  Some of the many weird geopotential initializations that we'll see today: ph0 is total,
      !  and grid%phtemp is a perturbation from the base state geopotential.  We set the base geopotential
      !  at the lowest level to terrain elevation * gravity.

      DO j=jts,jte
         DO i=its,ite
            ph0(i,1,j) = grid%ht(i,j) * g
            grid%phtemp(i,1,j) = 0.
         END DO
      END DO

      !  Base state potential temperature and inverse density (alpha = 1/rho) from
      !  the half eta levels and the base-profile surface pressure.  Compute 1/rho
      !  from equation of state.  The potential temperature is a perturbation from t0.

      DO j = jts, MIN(jte,jde-1)
         DO i = its, MIN(ite,ide-1)

            !  Base state pressure is a function of eta level and terrain, only, plus
            !  the hand full of constants: p00 (sea level pressure, Pa), t00 (sea level
            !  temperature, K), and A (temperature difference, from 1000 mb to 300 mb, K).

            p_surf = p00 * EXP ( -t00/a + ( (t00/a)**2 - 2.*g*grid%ht(i,j)/a/r_d ) **0.5 )


            DO k = 1, kte-1
               php(i,k,j) = grid%znw(k)*(p_surf - grid%p_top) + grid%p_top ! temporary, full lev base pressure
               grid%pbtemp(i,k,j) = grid%znu(k)*(p_surf - grid%p_top) + grid%p_top
               temp = MAX ( tiso, t00 + A*LOG(grid%pbtemp(i,k,j)/p00) )
               t_init(i,k,j) = temp*(p00/grid%pbtemp(i,k,j))**(r_d/cp) - t0
               alb(i,k,j) = (r_d/p1000mb)*(t_init(i,k,j)+t0)*(grid%pbtemp(i,k,j)/p1000mb)**cvpm
            END DO

            !  Base state mu is defined as base state surface pressure minus grid%p_top

            grid%mubtemp(i,j) = p_surf - grid%p_top

            !  Dry surface pressure is defined as the following (this mu is from the input file
            !  computed from the dry pressure).  Here the dry pressure is just reconstituted.

            pd_surf = mu0(i,j) + grid%p_top

            !  Integrate base geopotential, starting at terrain elevation.  This assures that
            !  the base state is in exact hydrostatic balance with respect to the model equations.
            !  This field is on full levels.

            phb(i,1,j) = grid%ht(i,j) * g
            DO k  = 2,kte
               phb(i,k,j) = phb(i,k-1,j) - grid%dnw(k-1)*grid%mubtemp(i,j)*alb(i,k-1,j)
            END DO
         END DO
      END DO

      !  Fill in the outer rows and columns to allow us to be sloppy.

      IF ( ite .EQ. ide ) THEN
      i = ide
      DO j = jts, MIN(jde-1,jte)
         grid%mubtemp(i,j) = grid%mubtemp(i-1,j)
         grid%mutemp(i,j) = grid%mutemp(i-1,j)
         DO k = 1, kte-1
            grid%pbtemp(i,k,j) = grid%pbtemp(i-1,k,j)
            t_init(i,k,j) = t_init(i-1,k,j)
            alb(i,k,j) = alb(i-1,k,j)
         END DO
         DO k = 1, kte
            phb(i,k,j) = phb(i-1,k,j)
         END DO
      END DO
      END IF

      IF ( jte .EQ. jde ) THEN
      j = jde
      DO i = its, ite
         grid%mubtemp(i,j) = grid%mubtemp(i,j-1)
         grid%mutemp(i,j) = grid%mutemp(i,j-1)
         DO k = 1, kte-1
            grid%pbtemp(i,k,j) = grid%pbtemp(i,k,j-1)
            t_init(i,k,j) = t_init(i,k,j-1)
            alb(i,k,j) = alb(i,k,j-1)
         END DO
         DO k = 1, kte
            phb(i,k,j) = phb(i,k,j-1)
         END DO
      END DO
      END IF

      !  Compute the perturbation dry pressure (grid%mubtemp + grid%mutemp + ptop = dry grid%psfc).

      DO j = jts, min(jde-1,jte)
         DO i = its, min(ide-1,ite)
            grid%mutemp(i,j) = mu0(i,j) - grid%mubtemp(i,j)
         END DO
      END DO

      !  Fill in the outer rows and columns to allow us to be sloppy.

      IF ( ite .EQ. ide ) THEN
      i = ide
      DO j = jts, MIN(jde-1,jte)
         grid%mutemp(i,j) = grid%mutemp(i-1,j)
      END DO
      END IF

      IF ( jte .EQ. jde ) THEN
      j = jde
      DO i = its, ite
         grid%mutemp(i,j) = grid%mutemp(i,j-1)
      END DO
      END IF

      lev500 = 0
      DO j = jts, min(jde-1,jte)
         DO i = its, min(ide-1,ite)

            !  Assign the potential temperature (perturbation from t0) and qv on all the mass
            !  point locations.

            DO k =  1 , kde-1
               grid%ttemp(i,k,j)          = grid%ttemp(i,k,j) - t0
            END DO

            dpmu = 10001.
            loop_count = 0

            DO WHILE ( ( ABS(dpmu) .GT. 10. ) .AND. &
                       ( loop_count .LT. 5 ) )

               loop_count = loop_count + 1

               !  Integrate the hydrostatic equation (from the RHS of the bigstep vertical momentum
               !  equation) down from the top to get the pressure perturbation.  First get the pressure
               !  perturbation, moisture, and inverse density (total and perturbation) at the top-most level.

               k = kte-1

               qvf1 = 0.5*(grid%moisttemp(i,k,j)+grid%moisttemp(i,k,j))
               qvf2 = 1./(1.+qvf1)
               qvf1 = qvf1*qvf2

               p3d(i,k,j) = - 0.5*(grid%mutemp(i,j)+qvf1*grid%mubtemp(i,j))/grid%rdnw(k)/qvf2
               qvf = 1. + rvovrd*grid%moisttemp(i,k,j)
               alt(i,k,j) = (r_d/p1000mb)*(grid%ttemp(i,k,j)+t0)*qvf&
                                 *(((p3d(i,k,j)+grid%pbtemp(i,k,j))/p1000mb)**cvpm)
               al(i,k,j) = alt(i,k,j) - alb(i,k,j)
               p_hyd(i,k,j) = p3d(i,k,j) + grid%pbtemp(i,k,j)

               !  Now, integrate down the column to compute the pressure perturbation, and diagnose the two
               !  inverse density fields (total and perturbation).

               DO k=kte-2,1,-1
                  qvf1 = 0.5*(grid%moisttemp(i,k,j)+grid%moisttemp(i,k+1,j))
                  qvf2 = 1./(1.+qvf1)
                  qvf1 = qvf1*qvf2
                  p3d(i,k,j) = p3d(i,k+1,j) - (grid%mutemp(i,j) + qvf1*grid%mubtemp(i,j))/qvf2/grid%rdn(k+1)
                  qvf = 1. + rvovrd*grid%moisttemp(i,k,j)
                  alt(i,k,j) = (r_d/p1000mb)*(grid%ttemp(i,k,j)+t0)*qvf* &
                              (((p3d(i,k,j)+grid%pbtemp(i,k,j))/p1000mb)**cvpm)
                  al(i,k,j) = alt(i,k,j) - alb(i,k,j)
                  p_hyd(i,k,j) = p3d(i,k,j) + grid%pbtemp(i,k,j)
               END DO

#if 1
               !  This is the hydrostatic equation used in the model after the small timesteps.  In
               !  the model, grid%al (inverse density) is computed from the geopotential.

               DO k  = 2,kte
                  grid%phtemp(i,k,j) = grid%phtemp(i,k-1,j) - &
                                grid%dnw(k-1) * ( (grid%mubtemp(i,j)+grid%mutemp(i,j))*al(i,k-1,j) &
                              + grid%mutemp(i,j)*alb(i,k-1,j) )
                  ph0(i,k,j) = grid%phtemp(i,k,j) + phb(i,k,j)
               END DO
#else
               !  Get the perturbation geopotential from the 3d height array from WPS.

               DO k  = 2,kte
                  grid%phtemp(i,k,j) = ph0(i,k,j)*g - phb(i,k,j)
               END DO
#endif

               !  Adjust the column pressure so that the computed 500 mb height is close to the
               !  input value (of course, not when we are doing hybrid input).

               IF ( ( flag_metgrid .EQ. 1 ) .AND. ( i .EQ. its ) .AND. ( j .EQ. jts ) ) THEN
                  DO k = 1 , num_metgrid_levels
                     IF ( ABS ( grid%p_gc(i,k,j) - 50000. ) .LT. 1. ) THEN
                        lev500 = k
                        EXIT
                     END IF
                  END DO
               END IF

               !  We only do the adjustment of height if we have the input data on pressure
               !  surfaces, and folks have asked to do this option.

               IF ( ( flag_metgrid .EQ. 1 ) .AND. &
                    ( config_flags%adjust_heights ) .AND. &
                    ( lev500 .NE. 0 ) ) THEN

                  DO k = 2 , kte-1

                     !  Get the pressures on the full eta levels (grid%php is defined above as
                     !  the full-lev base pressure, an easy array to use for 3d space).

                     pl = php(i,k  ,j) + &
                          ( p3d(i,k-1  ,j) * ( grid%znw(k    ) - grid%znu(k  ) ) + &
                            p3d(i,k    ,j) * ( grid%znu(k-1  ) - grid%znw(k  ) ) ) / &
                          ( grid%znu(k-1  ) - grid%znu(k  ) )
                     pu = php(i,k+1,j) + &
                          ( p3d(i,k-1+1,j) * ( grid%znw(k  +1) - grid%znu(k+1) ) + &
                            p3d(i,k  +1,j) * ( grid%znu(k-1+1) - grid%znw(k+1) ) ) / &
                          ( grid%znu(k-1+1) - grid%znu(k+1) )

                     !  If these pressure levels trap 500 mb, use them to interpolate
                     !  to the 500 mb level of the computed height.

                     IF ( ( pl .GE. 50000. ) .AND. ( pu .LT. 50000. ) ) THEN
                        zl = ( grid%phtemp(i,k  ,j) + phb(i,k  ,j) ) / g
                        zu = ( grid%phtemp(i,k+1,j) + phb(i,k+1,j) ) / g

                        z500 = ( zl * ( LOG(50000.) - LOG(pu    ) ) + &
                                 zu * ( LOG(pl    ) - LOG(50000.) ) ) / &
                               ( LOG(pl) - LOG(pu) )

                        !  Compute the difference of the 500 mb heights (computed minus input), and
                        !  then the change in grid%mutemp.  The grid%php is still full-levels, base pressure.

                        dz500 = z500 - grid%ght_gc(i,lev500,j)
                        tvsfc = ((grid%ttemp(i,1,j)+t0)*((p3d(i,1,j)+php(i,1,j))/p1000mb)**(r_d/cp)) * &
                                (1.+0.6*grid%moisttemp(i,1,j))
                        dpmu = ( php(i,1,j) + p3d(i,1,j) ) * EXP ( g * dz500 / ( r_d * tvsfc ) )
                        dpmu = dpmu - ( php(i,1,j) + p3d(i,1,j) )
                        grid%mutemp(i,j) = grid%mutemp(i,j) - dpmu
                        EXIT
                     END IF

                  END DO
               ELSE
                  dpmu = 0.
               END IF

            END DO

         END DO
      END DO
     
      !  Recompute Qv with full pressure on eta levels.  We use the temporary
      !  space of vtemp for total theta, which gets diagnosed into temperature.
      !  The utemp array is RH on eta levels from the initial interpolation.  The
      !  previous value of moist was computed with dry pressure, not total
      !  pressure.

      vtemp = grid%ttemp + 300
      CALL theta_to_t ( vtemp , grid%pbtemp + p3d , p00 , &
                        ids , ide , jds , jde , kds , kde , &
                        ims , ime , jms , jme , kms , kme , &
                        its , ite , jts , jte , kts , kte )

      CALL rh_to_mxrat (utemp, vtemp, grid%pbtemp + p3d, grid%moisttemp, &
                        config_flags%rh2qv_wrt_liquid ,                           &
                        config_flags%qv_max_p_safe ,                              &
                        config_flags%qv_max_flag , config_flags%qv_max_value ,    &
                        config_flags%qv_min_p_safe ,                              &
                        config_flags%qv_min_flag , config_flags%qv_min_value ,    &
                        ids , ide , jds , jde , kds , kde ,                       &
                        ims , ime , jms , jme , kms , kme ,                       &
                        its , ite , jts , jte , kts , kte-1 )

!  Set flag to denote that we are saving original values of HT, MUB, and
!  PHB for 2-way nesting and cycling.

      grid%save_topo_from_real=1

      ips = its ; ipe = ite ; jps = jts ; jpe = jte ; kps = kts ; kpe = kte
#ifdef DM_PARALLEL
#   include "HALO_EM_INIT_TEMP_1.inc"
#   include "HALO_EM_INIT_TEMP_2.inc"
#   include "HALO_EM_INIT_TEMP_3.inc"
#   include "HALO_EM_INIT_TEMP_4.inc"
#   include "HALO_EM_INIT_TEMP_5.inc"
#endif
      
      RETURN

   END SUBROUTINE bdy_prep_othertime

!-------------------------------------------------------------------
	
   SUBROUTINE assemble_output ( grid , loop )

   USE module_big_step_utilities_em
   USE module_bc
   IMPLICIT NONE

   TYPE(domain)                 :: grid
   INTEGER , INTENT(IN)         :: loop 

   INTEGER :: ids , ide , jds , jde , kds , kde
   INTEGER :: ims , ime , jms , jme , kms , kme
   INTEGER :: ips , ipe , jps , jpe , kps , kpe
   INTEGER :: ijds , ijde , spec_bdy_width
   INTEGER :: i , j , k , idts

   INTEGER :: id1 , interval_seconds , ierr, rc, sst_update, grid_fdda
   INTEGER , SAVE :: id, id2,  id4
   
   LOGICAL :: RESTART

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

   ijds = MIN ( ids , jds )
   ijde = MAX ( ide , jde )

   !  Boundary width, scalar value.

   spec_bdy_width = model_config_rec%spec_bdy_width
   interval_seconds = model_config_rec%interval_seconds
   sst_update = model_config_rec%sst_update
   grid_fdda = model_config_rec%grid_fdda(grid%id)

   IF ( loop .EQ. 1 ) THEN


      IF (  model_config_rec%polar(grid%id) ) THEN

         !  No need to couple data since no lateral BCs required.

      ELSE

         !  We need to save the 3d data to compute a difference during the next loop.  Couple the
         !  3d fields with total mu (grid%mubtemp + grid%mutemp) and the stagger-specific map scale factor.
   
         !  u, theta, h, scalars coupled with my; v coupled with mx
         CALL couple ( grid%mu_2 , grid%mub , ubdy3dtemp1 , grid%u_2                 , 'u' , grid%msfuy , &
                       ids, ide, jds, jde, kds, kde, ims, ime, jms, jme, kms, kme, ips, ipe, jps, jpe, kps, kpe )
         CALL couple ( grid%mu_2 , grid%mub , vbdy3dtemp1 , grid%v_2                 , 'v' , grid%msfvx , &
                       ids, ide, jds, jde, kds, kde, ims, ime, jms, jme, kms, kme, ips, ipe, jps, jpe, kps, kpe )
         CALL couple ( grid%mu_2 , grid%mub , tbdy3dtemp1 , grid%t_2                 , 't' , grid%msfty , &
                       ids, ide, jds, jde, kds, kde, ims, ime, jms, jme, kms, kme, ips, ipe, jps, jpe, kps, kpe )
         CALL couple ( grid%mu_2 , grid%mub , pbdy3dtemp1 , grid%ph_2                , 'h' , grid%msfty , &
                       ids, ide, jds, jde, kds, kde, ims, ime, jms, jme, kms, kme, ips, ipe, jps, jpe, kps, kpe )
         CALL couple ( grid%mu_2 , grid%mub , qbdy3dtemp1 , grid%moist(:,:,:,P_QV)      , 't' , grid%msfty , &
                       ids, ide, jds, jde, kds, kde, ims, ime, jms, jme, kms, kme, ips, ipe, jps, jpe, kps, kpe )
   
         DO j = jps , MIN(jde-1,jpe)
            DO i = ips , MIN(ide-1,ipe)
               mbdy2dtemp1(i,1,j) = grid%mu_2(i,j)
            END DO
         END DO

      END IF

      IF(grid_fdda .GE. 1)THEN
         DO j = jps , jpe
            DO k = kps , kpe
               DO i = ips , ipe
                  grid%fdda3d(i,k,j,p_u_ndg_old) = grid%u_2(i,k,j)
                  grid%fdda3d(i,k,j,p_v_ndg_old) = grid%v_2(i,k,j)
                  grid%fdda3d(i,k,j,p_t_ndg_old) = grid%t_2(i,k,j)
                  grid%fdda3d(i,k,j,p_q_ndg_old) = grid%moist(i,k,j,P_QV)
                  grid%fdda3d(i,k,j,p_ph_ndg_old) = grid%ph_2(i,k,j)
               END DO
            END DO
         END DO

         DO j = jps , jpe
            DO i = ips , ipe
               grid%fdda2d(i,1,j,p_mu_ndg_old) = grid%mu_2(i,j)
            END DO
         END DO
      END IF

      IF (  model_config_rec%polar(grid%id) ) THEN

         !  No need to build boundary arrays, since no lateral BCs are being generated.

      ELSE
   
         !  There are 2 components to the lateral boundaries.  First, there is the starting
         !  point of this time period - just the outer few rows and columns.
   
         CALL stuff_bdy     ( ubdy3dtemp1 , grid%u_bxs, grid%u_bxe, grid%u_bys, grid%u_bye, &
                                                              'U' , spec_bdy_width      , &
                                                                    ids , ide , jds , jde , kds , kde , &
                                                                    ims , ime , jms , jme , kms , kme , &
                                                                    ips , ipe , jps , jpe , kps , kpe )
         CALL stuff_bdy     ( vbdy3dtemp1 , grid%v_bxs, grid%v_bxe, grid%v_bys, grid%v_bye, &
                                                              'V' , spec_bdy_width      , &
                                                                    ids , ide , jds , jde , kds , kde , &
                                                                    ims , ime , jms , jme , kms , kme , &
                                                                    ips , ipe , jps , jpe , kps , kpe )
         CALL stuff_bdy     ( tbdy3dtemp1 , grid%t_bxs, grid%t_bxe, grid%t_bys, grid%t_bye, &
                                                              'T' , spec_bdy_width      , &
                                                                    ids , ide , jds , jde , kds , kde , &
                                                                    ims , ime , jms , jme , kms , kme , &
                                                                    ips , ipe , jps , jpe , kps , kpe )
         CALL stuff_bdy     ( pbdy3dtemp1 , grid%ph_bxs, grid%ph_bxe, grid%ph_bys, grid%ph_bye, &
                                                              'W' , spec_bdy_width      , &
                                                                    ids , ide , jds , jde , kds , kde , &
                                                                    ims , ime , jms , jme , kms , kme , &
                                                                    ips , ipe , jps , jpe , kps , kpe )
         CALL stuff_bdy     ( qbdy3dtemp1 , grid%moist_bxs(:,:,:,P_QV), grid%moist_bxe(:,:,:,P_QV),     &
                                            grid%moist_bys(:,:,:,P_QV), grid%moist_bye(:,:,:,P_QV),     &
                                                              'T' , spec_bdy_width      ,               &
                                                                    ids , ide , jds , jde , kds , kde , &
                                                                    ims , ime , jms , jme , kms , kme , &
                                                                    ips , ipe , jps , jpe , kps , kpe )
         CALL stuff_bdy     ( mbdy2dtemp1 , grid%mu_bxs, grid%mu_bxe, grid%mu_bys, grid%mu_bye, &
                                                              'M' , spec_bdy_width      , &
                                                                    ids , ide , jds , jde , 1 , 1 , &
                                                                    ims , ime , jms , jme , 1 , 1 , &
                                                                    ips , ipe , jps , jpe , 1 , 1 )
      END IF

               grid%dtbc = 0

   ELSE IF ( loop .GT. 1 ) THEN
                  
      IF ( model_config_rec%polar(grid%id) ) THEN

         !  No need to couple fields, since no lateral BCs are required.

      ELSE
   
         !  Couple this time period's data with total mu, and save it in the *bdy3dtemp2 arrays.
   
         !  u, theta, h, scalars coupled with my; v coupled with mx
         CALL couple ( grid%mutemp , grid%mubtemp , ubdy3dtemp2 , grid%utemp                 , 'u' , grid%msfuy , &
                       ids, ide, jds, jde, kds, kde, ims, ime, jms, jme, kms, kme, ips, ipe, jps, jpe, kps, kpe )
         CALL couple ( grid%mutemp , grid%mubtemp , vbdy3dtemp2 , grid%vtemp                 , 'v' , grid%msfvx , &
                       ids, ide, jds, jde, kds, kde, ims, ime, jms, jme, kms, kme, ips, ipe, jps, jpe, kps, kpe )
         CALL couple ( grid%mutemp , grid%mubtemp , tbdy3dtemp2 , grid%ttemp                 , 't' , grid%msfty , &
                       ids, ide, jds, jde, kds, kde, ims, ime, jms, jme, kms, kme, ips, ipe, jps, jpe, kps, kpe )
         CALL couple ( grid%mutemp , grid%mubtemp , pbdy3dtemp2 , grid%phtemp                , 'h' , grid%msfty , &
                       ids, ide, jds, jde, kds, kde, ims, ime, jms, jme, kms, kme, ips, ipe, jps, jpe, kps, kpe )
         CALL couple ( grid%mutemp , grid%mubtemp , qbdy3dtemp2 , grid%moisttemp, 't' , grid%msfty , &
                       ids, ide, jds, jde, kds, kde, ims, ime, jms, jme, kms, kme, ips, ipe, jps, jpe, kps, kpe )
   
         DO j = jps , jpe
            DO i = ips , ipe
               mbdy2dtemp2(i,1,j) = grid%mutemp(i,j)
            END DO
         END DO

      END IF

      IF(grid_fdda .GE. 1)THEN
         DO j = jps , jpe
            DO k = kps , kpe
               DO i = ips , ipe
                  grid%fdda3d(i,k,j,p_u_ndg_new) = grid%utemp(i,k,j)
                  grid%fdda3d(i,k,j,p_v_ndg_new) = grid%vtemp(i,k,j)
                  grid%fdda3d(i,k,j,p_t_ndg_new) = grid%ttemp(i,k,j)
                  grid%fdda3d(i,k,j,p_q_ndg_new) = grid%moisttemp(i,k,j)
                  grid%fdda3d(i,k,j,p_ph_ndg_new) = grid%phtemp(i,k,j)
               END DO
            END DO
         END DO

         DO j = jps , jpe
            DO i = ips , ipe
               grid%fdda2d(i,1,j,p_mu_ndg_new) = grid%mutemp(i,j)
            END DO
         END DO
      END IF

      IF ( model_config_rec%polar(grid%id) ) THEN

         !  No need to build boundary arrays, since no lateral BCs are being generated.

      ELSE

         !  During all of the loops after the first loop, we first compute the boundary
         !  tendencies with the current data values (*bdy3dtemp2 arrays) and the previously 
         !  saved information stored in the *bdy3dtemp1 arrays.
   
         CALL stuff_bdytend ( ubdy3dtemp2 , ubdy3dtemp1 , REAL(interval_seconds) ,                 &
                                                               grid%u_btxs, grid%u_btxe,     &
                                                               grid%u_btys, grid%u_btye,     &
                                                               'U' , &
                                                               spec_bdy_width      , &
                                                               ids , ide , jds , jde , kds , kde , &
                                                               ims , ime , jms , jme , kms , kme , &
                                                               ips , ipe , jps , jpe , kps , kpe )
         CALL stuff_bdytend ( vbdy3dtemp2 , vbdy3dtemp1 , REAL(interval_seconds) ,                 &
                                                               grid%v_btxs, grid%v_btxe,     &
                                                               grid%v_btys, grid%v_btye,     &
                                                               'V' , &
                                                               spec_bdy_width      , &
                                                               ids , ide , jds , jde , kds , kde , &
                                                               ims , ime , jms , jme , kms , kme , &
                                                               ips , ipe , jps , jpe , kps , kpe )
         CALL stuff_bdytend ( tbdy3dtemp2 , tbdy3dtemp1 , REAL(interval_seconds) ,                 &
                                                               grid%t_btxs, grid%t_btxe,     &
                                                               grid%t_btys, grid%t_btye,     &
                                                               'T' , &
                                                               spec_bdy_width      , &
                                                               ids , ide , jds , jde , kds , kde , &
                                                               ims , ime , jms , jme , kms , kme , &
                                                               ips , ipe , jps , jpe , kps , kpe )
         CALL stuff_bdytend ( pbdy3dtemp2 , pbdy3dtemp1 , REAL(interval_seconds) ,                 &
                                                               grid%ph_btxs, grid%ph_btxe,   &
                                                               grid%ph_btys, grid%ph_btye,   &
                                                               'W' , &
                                                               spec_bdy_width      , &
                                                               ids , ide , jds , jde , kds , kde , &
                                                               ims , ime , jms , jme , kms , kme , &
                                                               ips , ipe , jps , jpe , kps , kpe )
         CALL stuff_bdytend ( qbdy3dtemp2 , qbdy3dtemp1 , REAL(interval_seconds) ,                 &
                                                               grid%moist_btxs(:,:,:,P_QV), grid%moist_btxe(:,:,:,P_QV), &
                                                               grid%moist_btys(:,:,:,P_QV), grid%moist_btye(:,:,:,P_QV), &
                                                               'T' , &
                                                               spec_bdy_width      , &
                                                               ids , ide , jds , jde , kds , kde , &
                                                               ims , ime , jms , jme , kms , kme , &
                                                               ips , ipe , jps , jpe , kps , kpe )
         CALL stuff_bdytend ( mbdy2dtemp2 , mbdy2dtemp1 , REAL(interval_seconds) ,                 &
                                                               grid%mu_btxs, grid%mu_btxe,   &
                                                               grid%mu_btys, grid%mu_btye,   &
                                                               'M' , &
                                                               spec_bdy_width      , &
                                                               ids , ide , jds , jde , 1 , 1 , &
                                                               ims , ime , jms , jme , 1 , 1 , &
                                                               ips , ipe , jps , jpe , 1 , 1 )
      END IF
    
      
      !  Is this or is this not the last time time?  We can remove some unnecessary
      !  stores if it is not.

         IF ( model_config_rec%polar(grid%id) ) THEN

            !  No need to build boundary arrays, since no lateral BCs are being generated.

         ELSE

            !  There are 2 components to the lateral boundaries.  First, there is the starting
            !  point of this time period - just the outer few rows and columns.
   
            CALL stuff_bdy     ( ubdy3dtemp1 , grid%u_bxs, grid%u_bxe, grid%u_bys, grid%u_bye, &
                                                                 'U' , spec_bdy_width      , &
                                                                       ids , ide , jds , jde , kds , kde , &
                                                                       ims , ime , jms , jme , kms , kme , &
                                                                       ips , ipe , jps , jpe , kps , kpe )
            CALL stuff_bdy     ( vbdy3dtemp1 , grid%v_bxs, grid%v_bxe, grid%v_bys, grid%v_bye, &
                                                                 'V' , spec_bdy_width      , &
                                                                       ids , ide , jds , jde , kds , kde , &
                                                                       ims , ime , jms , jme , kms , kme , &
                                                                       ips , ipe , jps , jpe , kps , kpe )
            CALL stuff_bdy     ( tbdy3dtemp1 , grid%t_bxs, grid%t_bxe, grid%t_bys, grid%t_bye, &
                                                                 'T' , spec_bdy_width      , &
                                                                       ids , ide , jds , jde , kds , kde , &
                                                                       ims , ime , jms , jme , kms , kme , &
                                                                       ips , ipe , jps , jpe , kps , kpe )
            CALL stuff_bdy     ( pbdy3dtemp1 , grid%ph_bxs, grid%ph_bxe, grid%ph_bys, grid%ph_bye, &
                                                                 'W' , spec_bdy_width      , &
                                                                       ids , ide , jds , jde , kds , kde , &
                                                                       ims , ime , jms , jme , kms , kme , &
                                                                       ips , ipe , jps , jpe , kps , kpe )
            CALL stuff_bdy     ( qbdy3dtemp1 , grid%moist_bxs(:,:,:,P_QV), grid%moist_bxe(:,:,:,P_QV),     &
                                               grid%moist_bys(:,:,:,P_QV), grid%moist_bye(:,:,:,P_QV),     &
                                                                 'T' , spec_bdy_width      ,               &
                                                                       ids , ide , jds , jde , kds , kde , &
                                                                       ims , ime , jms , jme , kms , kme , &
                                                                       ips , ipe , jps , jpe , kps , kpe )
            CALL stuff_bdy     ( mbdy2dtemp1 , grid%mu_bxs, grid%mu_bxe, grid%mu_bys, grid%mu_bye, &
                                                                 'M' , spec_bdy_width      , &
                                                                       ids , ide , jds , jde , 1 , 1 , &
                                                                       ims , ime , jms , jme , 1 , 1 , &
                                                                       ips , ipe , jps , jpe , 1 , 1 )
   
         END IF

         IF ( model_config_rec%polar(grid%id) ) THEN

            !  No need to swap old for new for the boundary data, it is not required.

         ELSE

            !  We need to save the 3d data to compute a difference during the next loop.  Couple the
            !  3d fields with total mu (grid%mubtemp + grid%mutemp) and the stagger-specific map scale factor.
            !  We load up the boundary data again for use in the next loop.

            DO j = jps , jpe
               DO k = kps , kpe
                  DO i = ips , ipe
                     ubdy3dtemp1(i,k,j) = ubdy3dtemp2(i,k,j)
                     vbdy3dtemp1(i,k,j) = vbdy3dtemp2(i,k,j)
                     tbdy3dtemp1(i,k,j) = tbdy3dtemp2(i,k,j)
                     pbdy3dtemp1(i,k,j) = pbdy3dtemp2(i,k,j)
                     qbdy3dtemp1(i,k,j) = qbdy3dtemp2(i,k,j)
                  END DO
               END DO
            END DO

            DO j = jps , jpe
               DO i = ips , ipe
                  mbdy2dtemp1(i,1,j) = mbdy2dtemp2(i,1,j)
               END DO
            END DO

         END IF

         IF(grid_fdda .GE. 1)THEN
            DO j = jps , jpe
               DO k = kps , kpe
                  DO i = ips , ipe
                     grid%fdda3d(i,k,j,p_u_ndg_old) = grid%fdda3d(i,k,j,p_u_ndg_new)
                     grid%fdda3d(i,k,j,p_v_ndg_old) = grid%fdda3d(i,k,j,p_v_ndg_new)
                     grid%fdda3d(i,k,j,p_t_ndg_old) = grid%fdda3d(i,k,j,p_t_ndg_new)
                     grid%fdda3d(i,k,j,p_q_ndg_old) = grid%fdda3d(i,k,j,p_q_ndg_new)
                     grid%fdda3d(i,k,j,p_ph_ndg_old) = grid%fdda3d(i,k,j,p_ph_ndg_new)
                  END DO
               END DO
            END DO

            DO j = jps , jpe
               DO i = ips , ipe
                  grid%fdda2d(i,1,j,p_mu_ndg_old) = grid%fdda2d(i,1,j,p_mu_ndg_new)
               END DO
            END DO
         END IF
   
         grid%dtbc = 0
                 
   END IF

      call exchange_bdy_halo(grid,spec_bdy_width,loop)  

END SUBROUTINE assemble_output

   SUBROUTINE restart_bdy(grid,loop)

   USE module_big_step_utilities_em
   USE module_bc
   IMPLICIT NONE

   TYPE(domain)                 :: grid
   INTEGER , INTENT(IN)         :: loop
   
   INTEGER :: ids , ide , jds , jde , kds , kde
   INTEGER :: ims , ime , jms , jme , kms , kme
   INTEGER :: ips , ipe , jps , jpe , kps , kpe
   INTEGER :: ijds , ijde , spec_bdy_width
   INTEGER :: i , j , k , idts

   INTEGER :: id1 , interval_seconds , ierr, rc, sst_update, grid_fdda
   INTEGER , SAVE :: id, id2,  id4 

   !  Boundary width, scalar value.

   spec_bdy_width = model_config_rec%spec_bdy_width
   interval_seconds = model_config_rec%interval_seconds
   sst_update = model_config_rec%sst_update
   grid_fdda = model_config_rec%grid_fdda(grid%id)

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

   ijds = MIN ( ids , jds )
   ijde = MAX ( ide , jde )

   !  Boundary width, scalar value.

   spec_bdy_width = model_config_rec%spec_bdy_width
   interval_seconds = model_config_rec%interval_seconds
   sst_update = model_config_rec%sst_update
   grid_fdda = model_config_rec%grid_fdda(grid%id)
   
   IF (  model_config_rec%polar(grid%id) ) THEN

         !  No need to couple data since no lateral BCs required.

      ELSE

         !  We need to save the 3d data to compute a difference during the next loop.  Couple the
         !  3d fields with total mu (grid%mubtemp + grid%mu_2) and the stagger-specific map scale factor.
   
         !  u, theta, h, scalars coupled with my; v coupled with mx
         CALL couple ( grid%mutemp , grid%mubtemp , ubdy3dtemp1 , grid%utemp                 , 'u' , grid%msfuy , &
                       ids, ide, jds, jde, kds, kde, ims, ime, jms, jme, kms, kme, ips, ipe, jps, jpe, kps, kpe )
         CALL couple ( grid%mutemp , grid%mubtemp , vbdy3dtemp1 , grid%vtemp                 , 'v' , grid%msfvx , &
                       ids, ide, jds, jde, kds, kde, ims, ime, jms, jme, kms, kme, ips, ipe, jps, jpe, kps, kpe )
         CALL couple ( grid%mutemp , grid%mubtemp , tbdy3dtemp1 , grid%ttemp                 , 't' , grid%msfty , &
                       ids, ide, jds, jde, kds, kde, ims, ime, jms, jme, kms, kme, ips, ipe, jps, jpe, kps, kpe )
         CALL couple ( grid%mutemp , grid%mubtemp , pbdy3dtemp1 , grid%phtemp                , 'h' , grid%msfty , &
                       ids, ide, jds, jde, kds, kde, ims, ime, jms, jme, kms, kme, ips, ipe, jps, jpe, kps, kpe )
         CALL couple ( grid%mutemp , grid%mubtemp , qbdy3dtemp1 , grid%moisttemp             , 't' , grid%msfty , &
                       ids, ide, jds, jde, kds, kde, ims, ime, jms, jme, kms, kme, ips, ipe, jps, jpe, kps, kpe )
   
         DO j = jps , MIN(jde-1,jpe)
            DO i = ips , MIN(ide-1,ipe)
               mbdy2dtemp1(i,1,j) = grid%mutemp(i,j)
            END DO
         END DO

      END IF

      IF(grid_fdda .GE. 1)THEN
         DO j = jps , jpe
            DO k = kps , kpe
               DO i = ips , ipe
                  grid%fdda3d(i,k,j,p_u_ndg_old) = grid%utemp(i,k,j)
                  grid%fdda3d(i,k,j,p_v_ndg_old) = grid%vtemp(i,k,j)
                  grid%fdda3d(i,k,j,p_t_ndg_old) = grid%ttemp(i,k,j)
                  grid%fdda3d(i,k,j,p_q_ndg_old) = grid%moisttemp(i,k,j)
                  grid%fdda3d(i,k,j,p_ph_ndg_old) = grid%phtemp(i,k,j)
               END DO
            END DO
         END DO

         DO j = jps , jpe
            DO i = ips , ipe
               grid%fdda2d(i,1,j,p_mu_ndg_old) = grid%mu_2(i,j)
            END DO
         END DO
      END IF

      IF (  model_config_rec%polar(grid%id) ) THEN

         !  No need to build boundary arrays, since no lateral BCs are being generated.

      ELSE
   
         !  There are 2 components to the lateral boundaries.  First, there is the starting
         !  point of this time period - just the outer few rows and columns.
   
         CALL stuff_bdy     ( ubdy3dtemp1 , grid%u_bxs, grid%u_bxe, grid%u_bys, grid%u_bye, &
                                                              'U' , spec_bdy_width      , &
                                                                    ids , ide , jds , jde , kds , kde , &
                                                                    ims , ime , jms , jme , kms , kme , &
                                                                    ips , ipe , jps , jpe , kps , kpe )
         CALL stuff_bdy     ( vbdy3dtemp1 , grid%v_bxs, grid%v_bxe, grid%v_bys, grid%v_bye, &
                                                              'V' , spec_bdy_width      , &
                                                                    ids , ide , jds , jde , kds , kde , &
                                                                    ims , ime , jms , jme , kms , kme , &
                                                                    ips , ipe , jps , jpe , kps , kpe )
         CALL stuff_bdy     ( tbdy3dtemp1 , grid%t_bxs, grid%t_bxe, grid%t_bys, grid%t_bye, &
                                                              'T' , spec_bdy_width      , &
                                                                    ids , ide , jds , jde , kds , kde , &
                                                                    ims , ime , jms , jme , kms , kme , &
                                                                    ips , ipe , jps , jpe , kps , kpe )
         CALL stuff_bdy     ( pbdy3dtemp1 , grid%ph_bxs, grid%ph_bxe, grid%ph_bys, grid%ph_bye, &
                                                              'W' , spec_bdy_width      , &
                                                                    ids , ide , jds , jde , kds , kde , &
                                                                    ims , ime , jms , jme , kms , kme , &
                                                                    ips , ipe , jps , jpe , kps , kpe )
         CALL stuff_bdy     ( qbdy3dtemp1 , grid%moist_bxs(:,:,:,P_QV), grid%moist_bxe(:,:,:,P_QV),     &
                                            grid%moist_bys(:,:,:,P_QV), grid%moist_bye(:,:,:,P_QV),     &
                                                              'T' , spec_bdy_width      ,               &
                                                                    ids , ide , jds , jde , kds , kde , &
                                                                    ims , ime , jms , jme , kms , kme , &
                                                                    ips , ipe , jps , jpe , kps , kpe )
         CALL stuff_bdy     ( mbdy2dtemp1 , grid%mu_bxs, grid%mu_bxe, grid%mu_bys, grid%mu_bye, &
                                                              'M' , spec_bdy_width      , &
                                                                    ids , ide , jds , jde , 1 , 1 , &
                                                                    ims , ime , jms , jme , 1 , 1 , &
                                                                    ips , ipe , jps , jpe , 1 , 1 )
      END IF
  
      call exchange_bdy_halo(grid,spec_bdy_width,loop)  
      
      grid%dtbc = 0
            
   end subroutine restart_bdy

      subroutine exchange_bdy_halo(grid,spec_bdy_width,loop)
   
       TYPE (domain)          :: grid
       integer, intent(in)   :: spec_bdy_width,loop

      !  Local domain indices and counters.
      
      INTEGER :: ids, ide, jds, jde, kds, kde, &
                 ims, ime, jms, jme, kms, kme, &
                 its, ite, jts, jte, kts, kte, &
                 ips, ipe, jps, jpe, kps, kpe, &
                 i, j, k
                 
      if(loop.gt.1) then

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

        ! u
        call exchange_ybdy_halo_r(grid%u_bxs, spec_bdy_width, &
                              ims, ime, jms, jme, kms, kme, &
                              ips, ipe, jps, jpe, kps, kpe, &
                              ids, ide, jds, jde, kds, kde)

        call exchange_ybdy_halo_r(grid%u_bxe, spec_bdy_width, &
                              ims, ime, jms, jme, kms, kme, &
                              ips, ipe, jps, jpe, kps, kpe, &
                              ids, ide, jds, jde, kds, kde)

        call exchange_xbdy_halo_r(grid%u_bys, spec_bdy_width, &
                              ims, ime, jms, jme, kms, kme, &
                              ips, ipe, jps, jpe, kps, kpe, &
                              ids, ide, jds, jde, kds, kde)
                              
        call exchange_xbdy_halo_r(grid%u_bye, spec_bdy_width, &
                              ims, ime, jms, jme, kms, kme, &
                              ips, ipe, jps, jpe, kps, kpe, &
                              ids, ide, jds, jde, kds, kde)

        call exchange_ybdy_halo_r(grid%u_btxs, spec_bdy_width, &
                              ims, ime, jms, jme, kms, kme, &
                              ips, ipe, jps, jpe, kps, kpe, &
                              ids, ide, jds, jde, kds, kde)
                              
        call exchange_ybdy_halo_r(grid%u_btxe, spec_bdy_width, &
                              ims, ime, jms, jme, kms, kme, &
                              ips, ipe, jps, jpe, kps, kpe, &
                              ids, ide, jds, jde, kds, kde)

        call exchange_xbdy_halo_r(grid%u_btys, spec_bdy_width, &
                              ims, ime, jms, jme, kms, kme, &
                              ips, ipe, jps, jpe, kps, kpe, &
                              ids, ide, jds, jde, kds, kde)
                              
        call exchange_xbdy_halo_r(grid%u_btye, spec_bdy_width, &
                              ims, ime, jms, jme, kms, kme, &
                              ips, ipe, jps, jpe, kps, kpe, &
                              ids, ide, jds, jde, kds, kde)
        ! v
        call exchange_ybdy_halo_r(grid%v_bxs, spec_bdy_width, &
                              ims, ime, jms, jme, kms, kme, &
                              ips, ipe, jps, jpe, kps, kpe, &
                              ids, ide, jds, jde, kds, kde)
                              
        call exchange_ybdy_halo_r(grid%v_bxe, spec_bdy_width, &
                              ims, ime, jms, jme, kms, kme, &
                              ips, ipe, jps, jpe, kps, kpe, &
                              ids, ide, jds, jde, kds, kde)

        call exchange_xbdy_halo_r(grid%v_bys, spec_bdy_width, &
                              ims, ime, jms, jme, kms, kme, &
                              ips, ipe, jps, jpe, kps, kpe, &
                              ids, ide, jds, jde, kds, kde)
                              
        call exchange_xbdy_halo_r(grid%v_bye, spec_bdy_width, &
                              ims, ime, jms, jme, kms, kme, &
                              ips, ipe, jps, jpe, kps, kpe, &
                              ids, ide, jds, jde, kds, kde)

        call exchange_ybdy_halo_r(grid%v_btxs, spec_bdy_width, &
                              ims, ime, jms, jme, kms, kme, &
                              ips, ipe, jps, jpe, kps, kpe, &
                              ids, ide, jds, jde, kds, kde)
                              
        call exchange_ybdy_halo_r(grid%v_btxe, spec_bdy_width, &
                              ims, ime, jms, jme, kms, kme, &
                              ips, ipe, jps, jpe, kps, kpe, &
                              ids, ide, jds, jde, kds, kde)

        call exchange_xbdy_halo_r(grid%v_btys, spec_bdy_width, &
                              ims, ime, jms, jme, kms, kme, &
                              ips, ipe, jps, jpe, kps, kpe, &
                              ids, ide, jds, jde, kds, kde)
                              
        call exchange_xbdy_halo_r(grid%v_btye, spec_bdy_width, &
                              ims, ime, jms, jme, kms, kme, &
                              ips, ipe, jps, jpe, kps, kpe, &
                              ids, ide, jds, jde, kds, kde)

        ! ph                       
        call exchange_ybdy_halo_r(grid%ph_bxs, spec_bdy_width, &
                              ims, ime, jms, jme, kms, kme, &
                              ips, ipe, jps, jpe, kps, kpe, &
                              ids, ide, jds, jde, kds, kde)
                              
        call exchange_ybdy_halo_r(grid%ph_bxe, spec_bdy_width, &
                              ims, ime, jms, jme, kms, kme, &
                              ips, ipe, jps, jpe, kps, kpe, &
                              ids, ide, jds, jde, kds, kde)

        call exchange_xbdy_halo_r(grid%ph_bys, spec_bdy_width, &
                              ims, ime, jms, jme, kms, kme, &
                              ips, ipe, jps, jpe, kps, kpe, &
                              ids, ide, jds, jde, kds, kde)
                              
        call exchange_xbdy_halo_r(grid%ph_bye, spec_bdy_width, &
                              ims, ime, jms, jme, kms, kme, &
                              ips, ipe, jps, jpe, kps, kpe, &
                              ids, ide, jds, jde, kds, kde)

        call exchange_ybdy_halo_r(grid%ph_btxs, spec_bdy_width, &
                              ims, ime, jms, jme, kms, kme, &
                              ips, ipe, jps, jpe, kps, kpe, &
                              ids, ide, jds, jde, kds, kde)
                              
        call exchange_ybdy_halo_r(grid%ph_btxe, spec_bdy_width, &
                              ims, ime, jms, jme, kms, kme, &
                              ips, ipe, jps, jpe, kps, kpe, &
                              ids, ide, jds, jde, kds, kde)

        call exchange_xbdy_halo_r(grid%ph_btys, spec_bdy_width, &
                              ims, ime, jms, jme, kms, kme, &
                              ips, ipe, jps, jpe, kps, kpe, &
                              ids, ide, jds, jde, kds, kde)
                              
        call exchange_xbdy_halo_r(grid%ph_btye, spec_bdy_width, &
                              ims, ime, jms, jme, kms, kme, &
                              ips, ipe, jps, jpe, kps, kpe, &
                              ids, ide, jds, jde, kds, kde)

       ! t
        call exchange_ybdy_halo_r(grid%t_bxs, spec_bdy_width, &
                              ims, ime, jms, jme, kms, kme, &
                              ips, ipe, jps, jpe, kps, kpe, &
                              ids, ide, jds, jde, kds, kde)
                              
        call exchange_ybdy_halo_r(grid%t_bxe, spec_bdy_width, &
                              ims, ime, jms, jme, kms, kme, &
                              ips, ipe, jps, jpe, kps, kpe, &
                              ids, ide, jds, jde, kds, kde)

        call exchange_xbdy_halo_r(grid%t_bys, spec_bdy_width, &
                              ims, ime, jms, jme, kms, kme, &
                              ips, ipe, jps, jpe, kps, kpe, &
                              ids, ide, jds, jde, kds, kde)
                              
        call exchange_xbdy_halo_r(grid%t_bye, spec_bdy_width, &
                              ims, ime, jms, jme, kms, kme, &
                              ips, ipe, jps, jpe, kps, kpe, &
                              ids, ide, jds, jde, kds, kde)

        call exchange_ybdy_halo_r(grid%t_btxs, spec_bdy_width, &
                              ims, ime, jms, jme, kms, kme, &
                              ips, ipe, jps, jpe, kps, kpe, &
                              ids, ide, jds, jde, kds, kde)
                              
        call exchange_ybdy_halo_r(grid%t_btxe, spec_bdy_width, &
                              ims, ime, jms, jme, kms, kme, &
                              ips, ipe, jps, jpe, kps, kpe, &
                              ids, ide, jds, jde, kds, kde)

        call exchange_xbdy_halo_r(grid%t_btys, spec_bdy_width, &
                              ims, ime, jms, jme, kms, kme, &
                              ips, ipe, jps, jpe, kps, kpe, &
                              ids, ide, jds, jde, kds, kde)
                              
        call exchange_xbdy_halo_r(grid%t_btye, spec_bdy_width, &
                              ims, ime, jms, jme, kms, kme, &
                              ips, ipe, jps, jpe, kps, kpe, &
                              ids, ide, jds, jde, kds, kde)
                              
        ! moisture
        call exchange_ybdy_halo_r(grid%moist_bxs(:,:,:,P_QV), spec_bdy_width, &
                              ims, ime, jms, jme, kms, kme, &
                              ips, ipe, jps, jpe, kps, kpe, &
                              ids, ide, jds, jde, kds, kde)
                              
        call exchange_ybdy_halo_r(grid%moist_bxe(:,:,:,P_QV), spec_bdy_width, &
                              ims, ime, jms, jme, kms, kme, &
                              ips, ipe, jps, jpe, kps, kpe, &
                              ids, ide, jds, jde, kds, kde)

        call exchange_xbdy_halo_r(grid%moist_bys(:,:,:,P_QV), spec_bdy_width, &
                              ims, ime, jms, jme, kms, kme, &
                              ips, ipe, jps, jpe, kps, kpe, &
                              ids, ide, jds, jde, kds, kde)
                              
        call exchange_xbdy_halo_r(grid%moist_bye(:,:,:,P_QV), spec_bdy_width, &
                              ims, ime, jms, jme, kms, kme, &
                              ips, ipe, jps, jpe, kps, kpe, &
                              ids, ide, jds, jde, kds, kde)

        call exchange_ybdy_halo_r(grid%moist_btxs(:,:,:,P_QV), spec_bdy_width, &
                              ims, ime, jms, jme, kms, kme, &
                              ips, ipe, jps, jpe, kps, kpe, &
                              ids, ide, jds, jde, kds, kde)
                              
        call exchange_ybdy_halo_r(grid%moist_btxe(:,:,:,P_QV), spec_bdy_width, &
                              ims, ime, jms, jme, kms, kme, &
                              ips, ipe, jps, jpe, kps, kpe, &
                              ids, ide, jds, jde, kds, kde)

        call exchange_xbdy_halo_r(grid%moist_btys(:,:,:,P_QV), spec_bdy_width, &
                              ims, ime, jms, jme, kms, kme, &
                              ips, ipe, jps, jpe, kps, kpe, &
                              ids, ide, jds, jde, kds, kde)
                              
        call exchange_xbdy_halo_r(grid%moist_btye(:,:,:,P_QV), spec_bdy_width, &
                              ims, ime, jms, jme, kms, kme, &
                              ips, ipe, jps, jpe, kps, kpe, &
                              ids, ide, jds, jde, kds, kde)  
                              
        !mu
        call exchange_ybdy_halo_r(grid%mu_bxs, spec_bdy_width, &
                              ims, ime, jms, jme, 1, 1, &
                              ips, ipe, jps, jpe, 1, 1, &
                              ids, ide, jds, jde, 1, 1)
                              
        call exchange_ybdy_halo_r(grid%mu_bxe, spec_bdy_width, &
                              ims, ime, jms, jme, 1, 1, &
                              ips, ipe, jps, jpe, 1, 1, &
                              ids, ide, jds, jde, 1, 1)

        call exchange_xbdy_halo_r(grid%mu_bys, spec_bdy_width, &
                              ims, ime, jms, jme, 1, 1, &
                              ips, ipe, jps, jpe, 1, 1, &
                              ids, ide, jds, jde, 1, 1)
                              
        call exchange_xbdy_halo_r(grid%mu_bye, spec_bdy_width, &
                              ims, ime, jms, jme, 1, 1, &
                              ips, ipe, jps, jpe, 1, 1, &
                              ids, ide, jds, jde, 1, 1)

        call exchange_ybdy_halo_r(grid%mu_btxs, spec_bdy_width, &
                              ims, ime, jms, jme, 1, 1, &
                              ips, ipe, jps, jpe, 1, 1, &
                              ids, ide, jds, jde, 1, 1)
                              
        call exchange_ybdy_halo_r(grid%mu_btxe, spec_bdy_width, &
                              ims, ime, jms, jme, 1, 1, &
                              ips, ipe, jps, jpe, 1, 1, &
                              ids, ide, jds, jde, 1, 1)

        call exchange_xbdy_halo_r(grid%mu_btys, spec_bdy_width, &
                              ims, ime, jms, jme, 1, 1, &
                              ips, ipe, jps, jpe, 1, 1, &
                              ids, ide, jds, jde, 1, 1)
                              
        call exchange_xbdy_halo_r(grid%mu_btye, spec_bdy_width, &
                              ims, ime, jms, jme, 1, 1, &
                              ips, ipe, jps, jpe, 1, 1, &
                              ids, ide, jds, jde, 1, 1)
   
       endif
                              
      end subroutine exchange_bdy_halo
   
END MODULE module_bdy_prep
