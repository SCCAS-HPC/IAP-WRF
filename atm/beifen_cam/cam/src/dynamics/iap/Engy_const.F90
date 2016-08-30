module Engy_const
!------------------------------------------------------------------------------------------------
! Purpose:  set constants related to energy   
! Author :  ZhangHe
! Completed : 2005.8.30
! Update: 2007.5.14, ZhangHe, update the name of using moduls &
!                             specify which variables to use only
!         2008.5, Wujianping, parellel version
!         2008.6.9, ZhangHe, available for both serial & parellel
!         2010.08, Juanxiong He, available for both serial & parallel
!         2011.5.4, ZhangHe, write GMKENG, GMAPEA, GMAPES and GMASS
!         OCT 2012, Jiang Jinrong, 2D parellel
! Reviewed: ZhangHe, 2011-11-18
!           ZhangHe, 2012-10-25
!           ZhangHe, 2013-01-24
! Modified: Zhang He, 2013-03-29, removed Q and added qtmp in sub. GMPRFL
!------------------------------------------------------------------------------------------------
   use shr_kind_mod, only: r8 => shr_kind_r8
   use IAP_grid,  only: NX, NY, NL, IB, IE, JB, JE, EX, IM
   use mathconst, only: ZERO, HALF, FOUR
   use pmgrid,    only: beglatdyn, endlatdyn, loc_JB, loc_JE, beglev, endlev, plon, plat, plev, &
                        beglatdynex, endlatdynex
   use Dyn_const,    only: PMTOP
#if (defined SPMD)
   use mpishorthand, only: mpir8, mpicom
   use spmd_utils,   only: iam, masterproc, npes
#endif

   implicit none
#if (defined SPMD)
   include 'mpif.h'
#endif
   save
   public

   real(r8) :: PESBGM     ! mean model value of surface pressure
   real(r8) :: PESBGA     ! global atmosphere climatic mass
   real(r8) :: PESBAM     ! climatic atmosphere mass unit area
   real(r8) :: GBAREA     ! global area
   real(r8) :: DAP(NY)    ! area weight of P&U grid
   real(r8) :: DAV(NY)    ! area weight of  V  grid
   real(r8) :: DF(NY,2)   ! DF(J,1) = DAP(J);  DF(J,2) = DAV(J)
   real(r8) :: DSGH(NL)   ! DSGH(K) = 0.5 * DSIG(K)
   real(r8) :: GMKENG     ! global mean kinetic energy per mass
   real(r8) :: GMAPEA     ! global mean available potential energy per mass
   real(r8) :: GMAPES     ! global mean available surface potential energy per mass
   real(r8) :: GMTAE      ! global mean total available energy per mass
   real(r8) :: GMPENA     ! global mean generalized internal energy per mass
   real(r8) :: GMPENS     ! global mean surface potential energy per mass
   real(r8) :: GMPSFC     ! global mean surface pressure 
   real(r8) :: GMASS      ! atmosphere mass unit area
   integer, parameter :: NGENG = 6  ! kinds of energy
   real(r8) :: ENGREF(NGENG)     ! reference energy which the increase rates are calculated
   real(r8) :: ENGCUR(NGENG)     ! current energy
   real(r8) :: ENGRAT(NGENG)     ! the increase rates of energy
   private NF, NLF
   integer, parameter :: NF  = 4  ! four prognostic variables: U, V, T, Q
   integer, parameter :: NLF = NL * NF
      
!================================================================================================
CONTAINS
!================================================================================================

!================================================================================================
   SUBROUTINE ENERGY
!-----------------------------------------------------------------------------------------------
!  COMPUTE GLOBAL DYNAMIC ENERGY BUGET
!-----------------------------------------------------------------------------------------------
      use stdatm,    only: H0B, P00, PSB
      use physconst, only: RAD, CP
      use Dyn_const, only: DSIG
      use IAP_prog,  only: Psa, P, GHS, UT, VT, TT, T
      use cam_control_mod, only: adiabatic    !zhh 2011-11-18
      use time_manager, only: get_step_size, &
                              is_first_step, is_first_restart_step, &
                              get_curr_date

      IMPLICIT NONE
!----------------------------------Local workspace----------------------------------------------      
      real(r8) :: WA(NX,NY)   ! available surface potential energy
      real(r8) :: WB(NX,NY)   ! surface potential energy
      real(r8) :: WU(NL)      ! WU = 0.5 * U^2 * DSIG(K) * P0
      real(r8) :: WV(NL)      ! WV = 0.5 * V^2 * DSIG(K) * P0
      real(r8) :: WT(NL)      ! WT = 0.5 * (Phi)^2 * DSIG(K) * P0
      real(r8) :: WE(NL)      ! WE = Pes * Cp*T * DSIG(K)
      real(r8) :: TKENG       ! total kinetic energy
      real(r8) :: TAPEA
      real(r8) :: TAPES
      real(r8) :: TPENA
      real(r8) :: TPENS
      real(r8) :: TMASS
      real(r8) :: DPS0
      real(r8) :: WAJ, WBJ
      real(r8) :: EUJ         ! kinetic energy (U component)
      real(r8) :: EVJ         ! kinetic energy (V component)
      real(r8) :: EAJ
      real(r8) :: EASJ        ! available surface potential energy
      real(r8) :: EEJ         ! generalized internal energy
      real(r8) :: ESJ         ! surface potential energy
      real(r8) :: TMJ
      real(r8) :: UT0, VT0, TT0, T0
      real(r8) :: DSGK
      real(r8) :: DSHK
      real(r8) :: EAK         ! available potential energy
      real(r8) :: EASI        ! available surface potential energy
      real(r8) :: DAPJ
      INTEGER  :: I,K,J       ! loop index
      real(r8) :: t1(6), t2(6)
      integer  :: ierr
      integer  :: yr, mon, day, ncsec, dtime, ymd, n

!-----------------------------------------------------------------------------------------------
!
      TKENG      = ZERO
      TAPEA      = ZERO
      TAPES      = ZERO
      TPENA      = ZERO
      TPENS      = ZERO
      TMASS      = ZERO
!!      DO J = loc_JB,loc_JE
      DO J = beglatdyn, endlatdyn
         if ( j >= JB .and. j <= JE ) then 
            DO I = IB,IE
               DPS0    = Psa (I,J) 
               WA(I,J) = H0B(I,J) * DPS0*DPS0 * HALF  ! WA = Ees, available surface potential energy
               WB(I,J) = P  (I,J) * GHS(I,J)          ! WB : surface potential energy
!       GHS : surface geopotential 
            END DO
!
         else   ! at north and south pole         
            DPS0 = Psa(IB,J)
            WAJ  = H0B(IB,J) * DPS0*DPS0 * HALF
            WBJ  = P  (IB,J) * GHS(IB,J)
            DO I   = IB,IE
               WA(I,J)    = WAJ
               WB(I,J)    = WBJ
            END DO
         end if
      END DO
!      
      DO J = beglatdyn, endlatdyn
         EUJ        = ZERO
         EVJ        = ZERO
         EAJ        = ZERO
         EASJ       = ZERO
         EEJ        = ZERO
         ESJ        = ZERO
         TMJ        = ZERO
         DO I = IB,IE
            DO K = beglev ,endlev           !jjr
               UT0     = UT(I,K,J)
               VT0     = VT(I,K,J)
               TT0     = TT(I,K,J)
               T0      = T (I,K,J)
               DSGK    = DSIG(K)
               DSHK    = DSGH(K)                     ! DSGH(K)=0.5 * DSIG(K)
               WU(K)   = DSHK  * ( UT0*UT0 ) * P00   ! WU = 0.5 * U^2 * DSIG(K) * P0
               WV(K)   = DSHK  * ( VT0*VT0 ) * P00   ! WV = 0.5 * V^2 * DSIG(K) * P0
               WT(K)   = DSHK  * ( TT0*TT0 ) * P00   ! WT = 0.5 * (Phi)^2 * DSIG(K) * P0
               WE(K)   = DSGK  * P(I,J) * CP * T0    ! WE = Pes * Cp*T * DSIG(K)
            END DO
!
            EAK        = ZERO
            DO K   =  beglev,endlev  
               EAK     = EAK   + WT(K)  ! available potential energy   \
               EUJ     = EUJ   + WU(K)  ! kinetic energy (U component)  |__ sum for K
               EVJ     = EVJ   + WV(K)  ! kinetic energy (V component)  |
               EEJ     = EEJ   + WE(K)  ! generalized internal energy  /
            END DO
!
            EASI       = WA(I,J)         !available surface potential energy \
            TMJ        = TMJ   +  P(I,J) !                                    |
            ESJ        = ESJ   + WB(I,J) !surface potential energy            |--> sum for I
            EASJ       = EASJ  + EASI    !available surface potential energy  |
            EAJ        = EAJ   + EAK     !available potential energy         /
         END DO  ! I = IB,IE
!
         DAPJ       = DAP(J)             ! DAP: area weight of P & U grid 
         TKENG      = TKENG + (DAPJ*EUJ + DAV(J)*EVJ)         !total kinetic energy, Ek \
         TAPEA      = TAPEA +  DAPJ*EAJ           !total available potential energy, Eep |
         TAPES      = TAPES +  DAPJ*EASJ  !total available surface potential energy, Ees |_ sum for J
         TPENA      = TPENA +  DAPJ*EEJ          !total generalized internal energy, Cp*T| 
         TPENS      = TPENS +  DAPJ*ESJ             !total surface potential energy, Es  | 
         TMASS      = TMASS +  DAPJ*TMJ                        ! total surface pressure / 
!        TMASS/g:  total atmospheric mass
! the actual total energy is shoud be divided by g, for ex. TKENG = TKENG / g
      END DO   ! J = 1 ,NY
!
      t1(1) = TKENG
      t1(2) = TAPEA
      t1(3) = TAPES
      t1(4) = TPENA
      t1(5) = TPENS
      t1(6) = TMASS
#if (defined SPMD)
      call mpi_reduce(t1,t2,6,mpir8,mpi_sum,0,mpicom,ierr)
      if (iam.eq.0) then
         GMKENG = t2(1) / PESBGA   ! mean kinetic energy per mass
         GMAPEA = t2(2) / PESBGA   ! mean available potential energy per mass
         GMAPES = t2(3) / PESBGA   ! mean available surface potential energy per mass
         GMPENA = t2(4) / PESBGA   ! mean generalized internal energy per mass
         GMPENS = t2(5) / PESBGA   ! mean surface potential energy per mass
         GMPSFC = t2(6) / GBAREA + PMTOP   ! mean surface pressure
         GMASS  = t2(6) / GBAREA
         GMTAE  = GMKENG + GMAPEA + GMAPES ! global mean total available energy per mass
! ========================= zhh 2011.05.04 ===============================
         if (adiabatic) then
! open the file
            if (is_first_step() .or. is_first_restart_step()) then
               open (unit=31,file='energy.txt') 
               dtime = get_step_size()
               n = 1
               print*, '================================================'
               print*, 'dtime =', dtime
               print*, '================================================'
            end if
! output energy
            call get_curr_date(yr, mon, day, ncsec)
            ymd = yr*10000 + mon*100 + day         
            print*, 'ncsec =', ncsec
            print*, 'GMKENG =', GMKENG
            if (ncsec==dtime) then
               write (31,*)  n, ymd, real(GMKENG), real(GMAPEA), real(GMAPES), real(GMTAE), real(GMASS) 
               n = n+1
               print*, 'GMKENG =', GMKENG
               print*, 'GMAPEA =', GMAPEA
               print*, 'GMAPES =', GMAPES
            end if
         end if
! ========================= zhh 2011.05.04 ===============================
      endif                         
#else
      GMKENG = t1(1) / PESBGA   ! mean kinetic energy per mass
      GMAPEA = t1(2) / PESBGA   ! mean available potential energy per mass
      GMAPES = t1(3) / PESBGA   ! mean available surface potential energy per mass
      GMPENA = t1(4) / PESBGA   ! mean generalized internal energy per mass
      GMPENS = t1(5) / PESBGA   ! mean surface potential energy per mass
      GMPSFC = t1(6) / GBAREA + PMTOP   ! mean surface pressure
      GMASS  = t1(6) / GBAREA
      GMTAE  = GMKENG + GMAPEA + GMAPES ! global mean total available energy per mass
#endif
      RETURN
   END SUBROUTINE

!================================================================================================
   SUBROUTINE GMPRFL(GMUVTQ)
!-----------------------------------------------------------------------------------------------
!  COMPUTE GLOBAL DYNAMIC PROFILES
!-----------------------------------------------------------------------------------------------
      use IAP_prog,  only: U, V, T
      use prognostics, only: q3, n3

      IMPLICIT NONE
!------------------------------------Arguments--------------------------------------------------
      real(r8), intent(out) :: GMUVTQ(NL,NF)  ! GLOBAL PROFILES of U, V, T, Q
!----------------------------------Local workspace----------------------------------------------
      real(r8) :: WW( NX, NLF, NY ) 
      real(r8) :: WM( NLF)
      real(r8) :: WD( NL, NF )        ! area weight of U, V, T, Q
      integer  :: L, K, KL, K1, K2, K3, K4, I, J, jp, ierr
#if (defined SPMD)
      real(r8) :: temp(NL,NF)         ! temporary storage for GMUVTQ.
#endif
      real(r8), allocatable :: qtmp(:,:,:)    ! specific humidity 
!-----------------------------------------------------------------------------------------------

      allocate(qtmp(NX,beglev:endlev,beglatdynex:endlatdynex))
      do j = beglatdyn, endlatdyn
         jp = plat + 1 - j
         do k = beglev,endlev
            do i = 1, plon
               qtmp(I+EX,K,J) = q3(I,K,1,Jp,n3)
            end do
         end do
      end do

      DO L = 1 ,NF   ! NF=4, for U, V, T, Q
         DO K = 1 ,NL
            GMUVTQ(K,L)= ZERO
         END DO
      END DO
      DO J = beglatdyn, endlatdyn
         DO K = 1 ,NL
            WD(K,1)    = DAP(J)  ! U
            WD(K,2)    = DAV(J)  ! V
            WD(K,3)    = DAP(J)  ! T
            WD(K,4)    = DAP(J)  ! Q
         END DO
         DO K   = 1 ,NLF  !NLF=NL*NF
            WM(K)      = ZERO
         END DO
         DO I = IB,IE
!------------------------------------- define  WW(I,K,J) --------------------------------
            DO K = 1, NL
               K1 = K
               WW(I,K,J) = U(I,K1,J)
            END DO
            DO K  = NL+1, 2*NL
               K2 = K - NL
               WW(I,K,J) = V(I,K2,J)
            END DO
            DO K = 2*NL+1, 3*NL
               K3 = K - 2*NL
               WW(I,K,J) = T(I,K3,J)
            END DO
            DO K = 3*NL+1, 4*NL
               K4 = K - 3*NL
               WW(I,K,J) = qtmp(I,K4,J)
            END DO
!------------------------------------------------------------------------------------------
            DO K = 1 ,NLF
               WM(K)   = WM(K) + WW(I,K,J) 
            END DO
         END DO
         DO L = 1 ,NF
            DO K = 1 ,NL
               KL      = (L-1)*NL+K        ! equal KL = 1,NF*NL,1
               GMUVTQ(K,L) = GMUVTQ(K,L) + WM(KL)*WD(K,L)
            END DO
         END DO
      END DO
#if (defined SPMD)
      call mpi_allreduce(GMUVTQ, temp, NLF, mpir8, mpi_sum, mpicom, ierr)
      DO L = 1 ,NF
         DO K = 1 ,NL
            GMUVTQ(K,L) = temp(K,L) / GBAREA  !GBAREA: total global area
!           GMUVTQ(K,1) = U(K) , U PROFILE
!           GMUVTQ(K,2) = V(K) , V PROFILE
!           GMUVTQ(K,3) = T(K) , T PROFILE
!           GMUVTQ(K,4) = Q(K) , Q PROFILE
         END DO
      END DO
#else
      DO L = 1 ,NF
         DO K = 1 ,NL
            GMUVTQ(K,L) = GMUVTQ(K,L) / GBAREA  !GBAREA: total global area
         END DO
      END DO
#endif

      deallocate(qtmp)

      RETURN
   END SUBROUTINE

!================================================================================================
   SUBROUTINE GAREAM(F,IX,IP,GMF)
!-----------------------------------------------------------------------------------------------
!  COMPUTE THE GLOBAL AREA MEAN OF A VARIABLE F
!-----------------------------------------------------------------------------------------------
      IMPLICIT NONE
!------------------------------------Arguments--------------------------------------------------
      integer,  intent(in)  :: IX          ! number of F's first dimension
      integer,  intent(in)  :: IP          ! index of P or V grid
      real(r8), intent(in)  :: F(IX,NY)    ! input variable
      real(r8), intent(out) :: GMF         ! GLOBAL AREA MEAN OF F
!----------------------------------Local workspace----------------------------------------------
      real(r8) :: GMJ, temp
      integer  :: I, J, ierr
!-----------------------------------------------------------------------------------------------

      GMF        = ZERO
      DO J = beglatdyn, endlatdyn
         GMJ     = ZERO
         DO I = 1 ,IM
            GMJ  = GMJ + F(I,J)
         END DO
         GMF     = GMF + GMJ*DF(J,IP)  !DF(J,1) = DAP(J);  DF(J,2) = DAV(J)
      END DO
#if (defined SPMD)
      call mpi_allreduce(GMF, temp, 1, mpir8, mpi_sum, mpicom, ierr)
      GMF = temp / GBAREA
#else
      GMF = GMF / GBAREA
#endif
      RETURN
   END SUBROUTINE

!================================================================================================
   SUBROUTINE GPAREA(FFF,IX,IP,GMF)
!-----------------------------------------------------------------------------------------------
!  COMPUTE THE GLOBAL AREA MASS MEAN OF A VARIABLE FFF
!-----------------------------------------------------------------------------------------------
      use Dyn_const, only: DSIG
      use IAP_prog,  only: P

      IMPLICIT NONE
!------------------------------------Arguments--------------------------------------------------
      integer,  intent(in)  :: IX              ! number of FFF's first dimension
      integer,  intent(in)  :: IP              ! index of P or V grid
      real(r8), intent(in)  :: FFF(IX,NL,NY)   ! input variable
      real(r8), intent(out) :: GMF             ! GLOBAL AREA MASS MEAN OF FFF
!----------------------------------Local workspace----------------------------------------------
      real(r8) :: GMJ, temp
      real(r8) :: W( NX, NY )
      integer  :: I, J, K
      integer  :: ierr
!-----------------------------------------------------------------------------------------------

      GMF        = ZERO
      DO K = 1 ,NL
         DO J = beglatdyn, endlatdyn
            DO I = 1 ,IM
               W(I,J)  = FFF(I,K,J) * DSIG(K) * P(I,J)
            END DO
         END DO
         DO J = beglatdyn, endlatdyn
            GMJ        = ZERO
            DO I = 1 ,IM
               GMJ     = GMJ + W(I,J)
            END DO
            GMF        = GMF + GMJ*DF(J,IP)
         END DO
      END DO
#if (defined SPMD)
      call mpi_allreduce(GMF, temp, 1, mpir8, mpi_sum, mpicom, ierr)
      GMF = temp / PESBGA
#else
      GMF = GMF / PESBGA
#endif
      RETURN
   END SUBROUTINE

!================================================================================================
   SUBROUTINE STENGC
!-----------------------------------------------------------------------------------------------
!  SET CONSTANTS USED BY SUB.ENERGY
!-----------------------------------------------------------------------------------------------
      use stdatm,    only: PSB
      use physconst, only: RAD
      use Dyn_const, only: DSIG, DLON, DLAT, PMTOP, FIM

      IMPLICIT NONE
!----------------------------------Local workspace----------------------------------------------
      real(r8) :: FENGC      ! coefficient for reducing trunction error
      real(r8) :: DS00
      real(r8) :: HALFDT
      real(r8) :: DAJ
      real(r8) :: YU
      real(r8) :: YV
      real(r8) :: DS11
      real(r8) :: DSJJ
      real(r8) :: SUMPES     ! global climatic value of surface pressure
      real(r8) :: GPOINT     ! total global points
      integer  :: I, J, K    ! loop index
#if (defined SPMD)
      integer numeles,rdispls(npes),rcounts(npes),ibeg,ilen,ierr
      real*8 ts(3,NY), tr(3,NY)
#endif
!-----------------------------------------------------------------------------------------------

      FENGC      = 1.0E-10      ! FOR REDUCING TRUNCATION ERROR
      DS00       = RAD*DLON * RAD*DLAT * FENGC
      HALFDT     = DLAT * HALF
      DO J   = loc_JB,loc_JE
         YU      = DLAT * DBLE( J-1 )
         YV      = YU   + HALFDT
         DAP(J)  = DS00 * SIN( YU )  ! DAP(J)=a^2 * d(sita) * d(lamda)*sin[sita(J)]
         DAV(J)  = DS00 * SIN( YV )  ! DAV(J)=a^2 * d(sita) * d(lamda)*sin[sita(J+1/2)]
      END DO

! for DAP & DVP at polar 
      DS11       = DS00 * SIN( HALFDT )
      YU         = DLAT * DBLE( JE )  
      DSJJ       = DS00 * SIN( YU - HALFDT ) 
      DAP(1 )    = DS11 / FOUR  
      DAP(NY)    = DSJJ / FOUR       ! DAP(NY) = DAP(1 )
      DAV(1 )    = DS11
      DAV(NY)    = ZERO
!  define DF 
      DO J =1,NY 
         DF(J,1) = DAP(J)
         DF(J,2) = DAV(J)
      END DO
! 
      GBAREA     = ZERO
      PESBGM     = ZERO
      PESBGA     = ZERO
      DO J = 1,NY 
         SUMPES     = ZERO
         DO I = IB,IE
            SUMPES  = SUMPES + ( PSB(I,J) - PMTOP )  ! global climatic value of surface pressure
         END DO
         DAJ     = DAP(J)
!jjr#if (defined SPMD)
!jjr         ts(1,J) = DAJ
!jjr         ts(2,J) = DAJ*SUMPES
!jjr         ts(3,J) =     SUMPES
!jjr#else
         GBAREA  = GBAREA + DAJ         ! global area , GBAREA = 4 * PI * a^2
         PESBGA  = PESBGA + DAJ*SUMPES  ! global atmosphere climatic mass
         PESBGM  = PESBGM +     SUMPES  ! sum of global model value of surface pressure
!jjr#endif
      END DO
!jjr#if (defined SPMD)
!      numeles = 3*(endlatdyn-beglatdyn+1)
!      call mpi_gather(numeles, 1, mpi_integer, rcounts, 1, mpi_integer, 0, mpicom, ierr)
!      rdispls(1) = 0
!      do i=1,npes-1
!         rdispls(i+1)=rdispls(i)+rcounts(i)
!      enddo
!      call mpi_gatherv(ts(1,beglatdyn),numeles,mpir8,tr,rcounts,rdispls,mpir8,0,mpicom,ierr)
!      if (iam.eq.0) then
!         do j=1,NY
!            GBAREA=GBAREA+tr(1,J)
!            PESBGA=PESBGA+tr(2,J)
!            PESBGM=PESBGM+tr(3,J)
!         enddo
!      endif
!      call mpi_bcast(GBAREA,1,mpir8,0,mpicom,ierr)
!      call mpi_bcast(PESBGA,1,mpir8,0,mpicom,ierr)
!      call mpi_bcast(PESBGM,1,mpir8,0,mpicom,ierr)
!#endif
      GPOINT     = FIM    * DBLE( NY )  ! FIM = DBLE( NLON )
      GBAREA     = GBAREA * FIM         ! equal sum from i=1 to i=IM
      PESBAM     = PESBGA / GBAREA      ! atmosphere climatic mass unit area
      PESBGM     = PESBGM / GPOINT      ! mean model value of surface pressure
      DO K = 1 ,NL
         DSGH(K) = HALF   * DSIG(K)
      END DO

      RETURN
   END SUBROUTINE

end module
