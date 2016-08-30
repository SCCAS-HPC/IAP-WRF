module sm9h
!------------------------------------------------------------------------------------------------
! Purpose: Set the consts for 9-POINT HORIZONTAL AREAL SMOOTHING
! Original version : SM9HAS.f (IAP 9L)
! Reconstructed to module : ZhangHe
! Completed : 2005.9.13
!             2008.5, Wujianping, parallel version
!             2008.6.9, ZhangHe, available for both serial & parallel
!             2012 October, Jiang Jinrong, for 2D parallel
!             2013-01-28, Zhang He, new mp_send3d, dyn_state is added
!------------------------------------------------------------------------------------------------
   use shr_kind_mod, only: r8 => shr_kind_r8
   use IAP_grid, only: NX, NY, IB, IE, JB, JE, period
   use pmgrid,   only: beglatdyn, endlatdyn, loc_JB, loc_JE,npr_y,&
                       beglatdynex,endlatdynex
   use spmd_utils, only: iam
   
   implicit none

   save
   public
 
   real(r8) :: DAP(NY,2)     ! U or V grid area weight  
   real(r8) :: RA00(NY,2)    ! weight factor of central value
   real(r8) :: RA01(NY,2)    ! weight factor of NORTH/SOUTH/EAST/WEST value
   real(r8) :: RA02(NY,2)    ! weight factor of NE/NW/SE/SW value

!================================================================================================
CONTAINS
!================================================================================================

!================================================================================================
   SUBROUTINE SM9HAS(F,IV)
!------------------------------------------------------------------------------------------------
!  Purpose: PERFORM 9-POINT HORIZONTAL AREAL SMOOTHING
!           IV = 1/2  FOR THE FIELD ON THE P/V GRID
!------------------------------------------------------------------------------------------------
      use dynamics_vars,   only: T_FVDYCORE_STATE         !zhh 2013-01-28
      use dyn_internal_state,   only : get_dyn_state      !zhh 2013-01-28
#if ( defined SPMD )
      use mod_comm, only: mp_send3d, mp_recv3d
#endif
      implicit none
!----------------------------------------Arguments-----------------------------------------------
      integer,  intent(in)    :: IV          ! index which denotes U or V grid
      real(r8), intent(inout) :: F(NX,NY)    ! input variable need smooth
!-------------------------------------Local workspace-------------------------------------------
      real(r8) :: W(NX,NY)
      real(r8) :: DJ, RD0J, RD1J, RD2J     
      type (T_FVDYCORE_STATE), pointer :: dyn_state   !zhh 2013-01-28
      integer  :: JF, JV, I, J, J1, JJ, I1, II, dest, src
      integer  :: commyz     !zhh 2013-01-28
!------------------------------------------------------------------------------------------------
!
      dyn_state => get_dyn_state()
      commyz = dyn_state%grid%commyz

!     WEIGHTED BY THE AREA SUROUNDING THE CELL
      JF           = NY - IV
      JV           = JF + 1
      DO J = beglatdyn ,min(JV,endlatdyn)
         DJ        = DAP(J,IV)
         DO I = IB, IE
            W(I,J) = F(I,J) * DJ
         END DO
         call period( W(1,J) )
      END DO
!
#if (defined SPMD)
! Latitude distribution for plat=16 on 3 processors.             !wjp 2007.05
!  procId:      iam=0     |      iam=1     |      iam=2          !wjp 2007.05
!  CAM:    01 02 03 04 05 | 06 07 08 09 10 | 11 12 13 14 15 16   !wjp 2007.05
!  IAP:    16 15 14 13 12 | 11 10 09 08 07 | 06 05 04 03 02 01   !wjp 2007.05
!                left     |     current    |      right          !wjp 2007.05
!!      call mpi_move_left(W(1,endlatdyn), W(1,beglatdyn-1),NX)
!!      call mpi_move_right( W(1,beglatdyn), W(1,endlatdyn+1),NX)
      src = iam+1
      dest  = iam-1
      if ( mod(iam,npr_y) == 0 ) dest = -1
      if ( mod(iam+1,npr_y) == 0 ) src = -1
      call mp_send3d( commyz, dest, src, NX, NY,1,                  &
                      1, NX, 1, NY,1,1,       &
                      1, NX, endlatdyn, endlatdyn, 1, 1, w )
      call mp_recv3d( commyz, src, NX,  NY,1,                       &
                      1, NX, 1, NY,1,1,       &
                      1, NX, beglatdynex, beglatdynex,1,1, w )
!
      src = iam-1
      dest  = iam+1
      if ( mod(iam,npr_y) == 0 ) src = -1
      if ( mod(iam+1,npr_y) == 0 ) dest = -1
      call mp_send3d( commyz, dest, src, NX,  NY,1,                 &
                      1, NX, 1,NY ,1,1,       &
                      1, NX, beglatdyn, beglatdyn,1,1,  w )
      call mp_recv3d( commyz, src, NX,  NY,1,                       &
                      1, NX, 1, NY,1,1,       &
                      1, NX, endlatdynex, endlatdynex,1,1, w)

#endif
!
      DO J = max(JB,beglatdyn),min(JF,endlatdyn)
         JJ        = J + 1
         J1        = J - 1
         RD0J      = RA00(J,IV)
         RD1J      = RA01(J,IV)
         RD2J      = RA02(J,IV)
         DO I = IB,IE
            II     = I + 1
            I1     = I - 1
            F(I,J) = ( W(II,J ) + W(I1,J ) + W(I ,J1) + W(I ,JJ) )*RD1J        &
                   + ( W(II,JJ) + W(II,J1) + W(I1,JJ) + W(I1,J1) )*RD2J + W(I,J)*RD0J         
         END DO
		 call period( F(1,J) )
      END DO
      RETURN
   END SUBROUTINE 

!================================================================================================
   SUBROUTINE STSM9C(IDWEAK)
!------------------------------------------------------------------------------------------------
!  Popurse: Set const used in sub. SM9HAS         
!------------------------------------------------------------------------------------------------
      use Dyn_const, only: DXYP, DXYV 
      implicit none
!----------------------------------------Arguments-----------------------------------------------
      integer, intent(in) :: IDWEAK   ! index of weight factor
!-------------------------------------Local workspace-------------------------------------------
      real(r8) :: FTRUNC              ! coefficient for reducing trunction error
      real(r8) :: DJ                  ! 
      integer  :: I, J                ! loop index
!------------------------------------------------------------------------------------------------

      FTRUNC      = 1.0E-10
      DO J = beglatdyn, endlatdyn
         DAP(J,1) = DXYP(J) * FTRUNC
         DAP(J,2) = DXYV(J) * FTRUNC
      END DO
!
      IF ( IDWEAK.EQ.+1 ) THEN  
!       ++++  WITH  THE CENTRAL VALUE WEIGHTED BY 1 / 2
!                   NORTH/SOUTH/EAST/WEST      BY 3 /32
!                   NE/NW/SE/SW                BY 1 /32
         DO I = 1, 2
            DO J = loc_JB, loc_JE
               DJ        = DAP(J,I)
               RA00(J,I) = 1.0D0 /( 2.00D0*DJ )
               RA01(J,I) = 3.0D0 /( 32.0D0*DJ )
               RA02(J,I) = 1.0D0 /( 32.0D0*DJ )
            END DO
		 END DO
      ELSE
!       ++++  WITH  THE CENTRAL VALUE WEIGHTED BY 1 / 4
!                   NORTH/SOUTH/EAST/WEST      BY 1 / 8
!                   NE/NW/SE/SW                BY 1 /16
         DO I = 1, 2
            DO J = loc_JB, loc_JE
               DJ        = DAP(J,I)
               RA00(J,I) = 1.0D0 /( 4.00D0*DJ )
               RA01(J,I) = 1.0D0 /( 8.00D0*DJ )
               RA02(J,I) = 1.0D0 /( 16.0D0*DJ )
            END DO
		 END DO
      ENDIF

      RETURN
   END SUBROUTINE

end module
