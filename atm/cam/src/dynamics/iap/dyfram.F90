SUBROUTINE DYFRAM( NENG,NFST,NSEQ, NSLT, adv_state, detam, etamid, cwava,   &
                  t3, u3, v3, q3, ps, omga, phis, n3, npt )    
!------------------------------------------------------------------------------------------------
! Purpose: Main routine of the dynamical frame calculation
! Original version: DYFRAM.f (IAP 9L)
! Reconstructed & modified : ZhangHe
! Completed : 2005.9.19
! Update : 2006.12.25, ZhangHe, removed dumb parameter MLF & DTSM
!          2007.4.25, ZhangHe,
!          2007.12.15, ZhangHe, changed 'trans_IAP' to 'trans_IAP_tend'
!                               'trans_antiIAP'  
!          ZhangHe, 2007.12.21
!          ZhangHe, 2008.04.18, forecast q by SLT from CAM3.1
!          ZhangHe, 2008.04.26, added dumb parameter 't3, u3, v3 ... '
!          ZhangHe, 2008.04.30, moved calling of sub. DIAGPP, DIAGBB & DIAGHI to sub. trans_af_mass
!          ZhangHe, 2008.06.03, for parallel version
!          juanxiong he, 2010.08, renamed trans_grid_dp to copy_grid_dp	
!          ZhangHe, 2010.8.10, added NSLT 
!          ZhangHe, 2010.8.12, modified UVW0  
!          ZhangHe, 2011.05.04, do not call scanslt_run if adiabatic run
!          Jiang Jinrong, October 2012, for 2D parellel
! Reviewed: ZhangHe, 2011-11-18
!           ZhangHe, 2012-10-23
!           ZhangHe, 2013-02-05
! Modified: ZhangHe, 2013-03-29, use dynamic arrays 
!------------------------------------------------------------------------------------------------
   use shr_kind_mod, only: r8 => shr_kind_r8
   use IAP_grid,   only: NX, NY, NL, NZ, period
   use Dyn_const,  only: SIG
   use flexib,     only: IBCFFT, DTadv, DTIMEQ 
   use IAP_prog,   only: WPA, WS, P, U, V, T, UVW0    
   use tendency,   only: DPsa, NADDSS, SU, SV, ST
   use Trans_coef, only: PHYSDM, trans_IAP_tend, trans_antiIAP  !zhh 2007.12.15
   use stdatm,     only: DIAGBB
   use Engy_const, only: ENERGY 
   use pmgrid,     only: plon, plat, plev, beglat, endlat,   &   
                         beglatdyn, endlatdyn,beglev,endlev,endlevp   !jjr 2012.09
   use spmd_utils, only: masterproc, iam   
   use constituents,  only: pcnst                    ! juanxiong he, 201008
   use scanslt,     only: scanslt_run, advection_state
   use perf_mod, only : t_startf, t_stopf ! juanxiong he
   use cam_control_mod, only: adiabatic, ideal_phys    !zhh, 2011-11-18
!
   IMPLICIT NONE
!------------------------------------Arguments--------------------------------------------------
   integer,  intent(in) :: NENG     ! switch to compute energy budget
   integer,  intent(in) :: NFST     ! switch to do smoothing 
   integer,  intent(in) :: NSEQ     ! dynamical computation times of every hour
   integer,  intent(in) :: NSLT     ! slt computation times of every DTadv
   integer,  intent(in) :: n3       ! current time level
   integer,  intent(in) :: npt      ! number of time levels in the dycore
!
   real(r8), intent(inout) :: t3(plon,beglev:endlev,beglat:endlat,npt)    ! temperature
   real(r8), intent(inout) :: u3(plon, beglev:endlev,beglat:endlat,npt)    ! u-wind component
   real(r8), intent(inout) :: v3(plon,beglev:endlev,beglat:endlat,npt)    ! v-wind component
   real(r8), intent(inout) :: q3(plon,beglev:endlev,pcnst,beglat:endlat,npt)   ! specific humidity
   real(r8), intent(inout) :: omga(plon,beglev:endlev,beglat:endlat)      ! p-surface vertical velocity
   real(r8), intent(inout) :: ps(plon,beglat:endlat,npt)         ! surface pressure
   real(r8), intent(inout) :: phis(plon,beglat:endlat)           ! surface geopotential
!
   type(advection_state), intent(inout) :: adv_state        ! Advection state data
   real(r8), intent(in) :: etamid(plev)     ! vertical coords at midpoints 
   real(r8), intent(inout) :: cwava(plat)   ! weight applied to global integrals
   real(r8), intent(inout) :: detam(plev)   ! intervals between vert full levs.
!----------------------------------Local workspace----------------------------------------------
   real(r8), allocatable :: Un(:,:,:,:)
   real(r8), allocatable :: Vn(:,:,:,:)
   real(r8), allocatable :: Wn(:,:,:,:)
   integer  :: I, J, K, NCYC, NI, jd        ! loop index
!------------------------------------------------------------------------------------------------
! allocate arrays
   allocate ( Un(NX,beglev:endlev,beglatdyn:endlatdyn,0:NSLT) )
   allocate ( Vn(NX,beglev:endlev,beglatdyn:endlatdyn,0:NSLT) )
   allocate ( Wn(NX,beglev:endlev,beglatdyn:endlatdyn,0:NSLT) )
!
!  PREPARE FOR THE CYCLE OF THE RUNNING
!
! -------------------- BACK FROM PHYSICS ROUTINE TO DYNAMICS ROUTINE --------------------------
   IF ( NADDSS == 0 ) THEN
      call PHYSDM(IBCFFT, U, V, T)    !zhh  2008.4.21
   ELSE 
      CALL trans_IAP_tend(IBCFFT, SU, SV, ST)    !zhh 2008.4.21
   END IF
!----------------------------------------------------------------------------------
!
!     START THE DYNAMICAL INTEGRATING CYCLE
   DO NCYC = 1, NSEQ   
! ============================ zhh 2010.8.12 =============================
      DO J = beglatdyn, endlatdyn
         DO K = beglev ,endlev
            DO I = 1 ,NX
               UVW0(I,K,J,1)  = U (I,K,J)
               UVW0(I,K,J,2)  = V (I,K,J)
               UVW0(I,K,J,3)  = WS(I,K,J)
            ENDDO
         ENDDO
      ENDDO
! ============================ zhh 2010.8.12 =============================
!     PREDICT DRY-ADIABATIC SYSTEM
!-------------------------- solve the dry dynamical framework -----------------------------
      call t_startf('NLITTI')
      CALL NLITTI        
      call t_stopf('NLITTI')
!-----------------------------------------------------------------------------------------
      call copy_grid_dp( 2, t3, u3, v3, q3, omga, ps, phis, n3, npt )
!
! SLT scan from south to north for forecasting moisture

      do NI = 0, NSLT    !zhh 2010.08.12
         do J = beglatdyn, endlatdyn
            do K = beglev ,endlev
               do I = 1 ,NX
                  Un(I,K,J,NI) = UVW0(I,K,J,1) + dble(NI)/dble(NSLT) * (U(I,K,J)-UVW0(I,K,J,1))
                  Vn(I,K,J,NI) = UVW0(I,K,J,2) + dble(NI)/dble(NSLT) * (V(I,K,J)-UVW0(I,K,J,2))
                  Wn(I,K,J,NI) = UVW0(I,K,J,3) + dble(NI)/dble(NSLT) * (WS(I,K,J)-UVW0(I,K,J,3))
               end do
            end do
         end do
      end do
!
      call t_startf('scanslt_run')
      do NI = 1, NSLT    !zhh 2010.08.10
         do J = beglatdyn, endlatdyn
            do K =  beglev,endlev
               do I = 1 ,NX
                  UVW0(I,K,J,1) = Un(I,K,J,NI-1) + Un(I,K,J,NI) - U(I,K,J)
                  UVW0(I,K,J,2) = Vn(I,K,J,NI-1) + Vn(I,K,J,NI) - V(I,K,J)
                  UVW0(I,K,J,3) = Wn(I,K,J,NI-1) + Wn(I,K,J,NI) - WS(I,K,J)
               end do
            end do
         end do
!
         if ((.not. ideal_phys) .and. (.not. adiabatic)) then
            call scanslt_run( adv_state, DTIMEQ, detam, etamid, cwava ) !zhh 2010.08.10
         end if 
      end do
      call t_stopf('scanslt_run')
!
   ENDDO
!     CLOSE THE DYNAMICAL INTEGRATING CYCLE
!
!     GET THE TIME AVERAGED P-SURFACE VERTICAL VELOCITY
   DO J = beglatdyn, endlatdyn
      DO K = beglev ,endlev+1
         DO I = 1 ,NX
            WPA(I,K,J) = WS(I,K,J)*P(I,J) + DPsa(I,J)*SIG(K)  
         ENDDO
      ENDDO
   ENDDO
   DO J = beglatdyn ,endlatdyn
      DO K = beglev ,endlev+1
         call period( WPA(1,K,J) )
      ENDDO
   ENDDO
!
!---------------------------------------------------------------------------------
   CALL DIAGPP
!---------------------------------------------------------------------------------
!---------------------SET MSA  TEMPERATURE & CHARACTERISTIC PHASE-SPEED-----------
   CALL DIAGBB( 2 )
!---------------------------------------------------------------------------------
!--------------- COMPUTE GEOPOTENTIAL HEIGHT DEPARTURE ---------------------------
   CALL DIAGHI( 0 )
!---------------------------------------------------------------------------------
!----------------------  CALCULATE ENERGY IF NECESSARY ---------------------------
   IF ( NENG.EQ.1 ) CALL ENERGY   
!---------------------------------------------------------------------------------
! deallocate arrays
   deallocate(Un)
   deallocate(Vn)
   deallocate(Wn)
!
   RETURN      ! end of call dyfram     
END
