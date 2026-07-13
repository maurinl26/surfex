!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!   ############################################################################
SUBROUTINE SNOW_LOAD_MEB(IO, PK, PEK, DEK, PTSTEP, PSR, PWRVNMAX, PCHEATV, PMELTVN, &
                         PVELC, PSUBVCOR, PCHIP, PWRMAX)
!   ############################################################################
!
!!****  *SNOW_LOAD_MEB*
!!
!!    PURPOSE
!!    -------
!
!     Calculate temporal evolution of canopy-intercepted snow
!     
!!**  METHOD
!!    ------
!
!
!!    EXTERNAL
!!    --------
!!
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!
!!      
!!    REFERENCE
!!    ---------
!!
!!      
!!    AUTHOR
!!    ------
!!
!!      P. Samuelsson           * SMHI *
!!      A. Boone                * CNRM-GAME, Meteo-France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    02/2011
!                   09/2022 (A. Bouchet, A. Boone) new snow unloading option added
!-------------------------------------------------------------------------------
!
!*       0.     DECLARATIONS
!               ------------
!
USE MODD_ISBA_OPTIONS_n,   ONLY : ISBA_OPTIONS_t
USE MODD_ISBA_n, ONLY : ISBA_PE_t, ISBA_P_t
USE MODD_DIAG_EVAP_ISBA_n, ONLY : DIAG_EVAP_ISBA_t
USE MODD_CSTS,     ONLY : XTT, XLMTT
USE MODE_MEB,              ONLY : SNOW_INTERCEPT_EFF
!
USE MODD_SNOW_PAR, ONLY : XRHOSMAX_ES
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*      0.1    Declaration of Arguments
!
TYPE(ISBA_OPTIONS_t), INTENT(INOUT) :: IO
TYPE(ISBA_P_t), INTENT(INOUT) :: PK
TYPE(ISBA_PE_t), INTENT(INOUT) :: PEK
TYPE(DIAG_EVAP_ISBA_t), INTENT(INOUT) :: DEK
!
REAL,               INTENT(IN)    :: PTSTEP
!
REAL, DIMENSION(:), INTENT(IN)    :: PSR, PCHEATV, PVELC, PMELTVN, PWRVNMAX, PCHIP, &
                                     PWRMAX
!
REAL, DIMENSION(:), INTENT(OUT)   :: PSUBVCOR
!
!
!*      0.2    declarations of local variables
!
REAL, DIMENSION(SIZE(PSR))        :: ZSRINT, ZWRVN, ZSUB, ZKVN, ZUNLOAD_POT
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!*      0.3    declarations of local parameters
!
! Snow unloading parameters (Roesch el al., Clim. Dyn., 2001)
!
REAL, PARAMETER                   :: ZUNLOAD_T     = 1.87E+5  ! K s ! initial : 1.87E+5
REAL, PARAMETER                   :: ZUNLOAD_TT    = 270.15   ! K   ! initial : 270.15
REAL, PARAMETER                   :: ZUNLOAD_V     = 1.56E+5  ! m   ! initial : 1.56E+5
!
! Snow unloading parameters (Lundquist, 2021)
!
REAL, PARAMETER                   :: ZUNLOAD_MT    = .25          ! --
REAL, PARAMETER                   :: ZUNLOAD_MV    = .25          ! --
REAL, PARAMETER                   :: ZUNLOAD_MTD   = 0.           ! --
REAL, PARAMETER                   :: ZUNLOAD_CT    = 1.87E+5      ! s-1
REAL, PARAMETER                   :: ZUNLOAD_CV    = 1.56E+5      ! s-1
REAL, PARAMETER                   :: ZUNLOAD_CTD   = 1.2861E-6    ! s-1
!
!-------------------------------------------------
! 0) Initialization
!
IF (LHOOK) CALL DR_HOOK('SNOW_LOAD_MEB',0,ZHOOK_HANDLE)
!
ZSRINT(:)      = 0.0
ZWRVN(:)       = 0.0
ZSUB(:)        = 0.0
ZKVN(:)        = 0.0
DEK%XUNLOAD(:) = 0.0
!
!
! 1) Potential Unloading (falling off branches, etc...): loss
! Note, the temperature effect is assumed to vanish for cold temperatures.
! This is the initial value which may be adjusted after WRVN has been updated.

IF(IO%LMEB_INT_ULUN)THEN
  ZUNLOAD_POT(:) = PTSTEP * PEK%XWRVN(:)*( 	                               &
                      ZUNLOAD_MT  * (MAX(0.0, PEK%XTV(:) - XTT) / ZUNLOAD_CT)  &
                    + ZUNLOAD_MV  * (MAX(0.0, PVELC(:)) / ZUNLOAD_CV)	       &
                    + ZUNLOAD_MTD * ZUNLOAD_CTD ) 
ELSE
  ZUNLOAD_POT(:) = PTSTEP * PEK%XWRVN(:)*(                           &
                       (MAX(0.0, PEK%XTV(:)-ZUNLOAD_TT) / ZUNLOAD_T) &
                     + (PVELC(:) / ZUNLOAD_V))
ENDIF

! 2) Consider the case when maximum interception is zero...
! this only occurs when vegetation canopy is *totally* buried. The follwing line
! results in non-zero snow loading (total removal of intercepted snow) 
! only during the timestep when vegetation has just been buried:
!
WHERE(PWRVNMAX(:) == 0.0)
  !
  DEK%XSR_GN(:) = PSR(:) + PEK%XWRVN(:)/PTSTEP    ! kg m-2 s-1
  PEK%XWRVN(:)    = 0.0

  ! for a totally buried canopy, the following are zero:
  DEK%XMELT_CV(:) = 0.0
  DEK%XFRZ_CV(:)  = 0.0
  PSUBVCOR(:)    = 0.0
  !
  ELSEWHERE
  !
  ! 2) Case for snow beneath or only partially covering the vegetation canopy:
  ! The following are computed as steps to ensure mass conservation.
  !
  ! Interception: gain
  ZKVN(:)        = SNOW_INTERCEPT_EFF(PCHIP,PVELC,PWRVNMAX)
  ZSRINT(:)      = MAX(0.0,PWRVNMAX(:)-PEK%XWRVN(:))*(1.0-EXP(-ZKVN(:)*PSR(:)*PTSTEP)) ! kg m-2
  ZSRINT(:)      = MIN(PSR(:)*PTSTEP, ZSRINT(:))  ! kg m-2 
  ZWRVN(:)       = PEK%XWRVN(:) + ZSRINT(:)           ! kg m-2 
  DEK%XSR_GN(:)  = MAX(0.0, PSR(:) - ZSRINT(:)/PTSTEP) ! kg m-2 s-1
  !
  ! Sublimation: gain or loss
  ! NOTE for the rare case that sublimation exceeds snow mass (possible as traces of snow disappear)
  ! compute a mass correction to be removed from soil (to conserve mass): PSUBVCOR
  ZSUB(:)        = DEK%XLES_CV(:)*(PTSTEP/PK%XLSTT(:))    ! kg m-2
  PSUBVCOR(:)    = MAX(0.0, ZSUB(:) - ZWRVN(:))/PTSTEP  ! kg m-2 s-1
  ZWRVN(:)       = MAX(0.0, ZWRVN(:) - ZSUB(:))         ! kg m-2
  !
  ! Phase change: loss (melt of snow mass)
  DEK%XMELT_CV(:) = PTSTEP*MAX(0.0, PMELTVN(:))         ! kg m-2  
  DEK%XMELT_CV(:) = MIN(DEK%XMELT_CV(:), ZWRVN(:))
  ZWRVN(:)        = ZWRVN(:)    - DEK%XMELT_CV(:)
  PEK%XWR(:)      = PEK%XWR(:)  + DEK%XMELT_CV(:)      ! NOTE: liq reservoir can exceed maximum holding
  !                                                    !       capacity here, but this is accounted for
  !                                                    !       in main prognostic PWRV routine.
  ! Phase change: gain (freeze of intercepted water) 
  ! Note, to get a better estimate of water available for freezing, remove Er in 
  ! estimation of water for freezing:
  ! Also, update liquid water stored on the canopy here:
  DEK%XFRZ_CV(:) = PTSTEP*MAX(0.0, -PMELTVN(:))        ! kg m-2  
  DEK%XFRZ_CV(:)  = MIN(DEK%XFRZ_CV(:), MAX(0.0,PEK%XWR(:) - &
                     DEK%XLER_CV(:)*(PTSTEP/PK%XLVTT(:))))
  ZWRVN(:)       = ZWRVN(:)   + DEK%XFRZ_CV(:)
  PEK%XWR(:)     = PEK%XWR(:) - DEK%XFRZ_CV(:)
  ! Unloading (possibly limited by existing intercepted snow):
  DEK%XUNLOAD(:)  = MIN(ZWRVN(:), ZUNLOAD_POT(:))         ! kg m-2 
  ZWRVN(:)       = ZWRVN(:) - DEK%XUNLOAD(:)                                           ! kg m-2 
  DEK%XSR_GN(:)   = DEK%XSR_GN(:) + DEK%XUNLOAD(:)/PTSTEP ! kg m-2 s-1
  !
  ! Diagnostic updates:
  ! final phase change (units)
  DEK%XMELT_CV(:) = DEK%XMELT_CV(:)/PTSTEP ! kg m-2 s-1
  DEK%XUNLOAD(:) = DEK%XUNLOAD(:)/PTSTEP ! kg m-2 s-1
  DEK%XFRZ_CV(:)  = DEK%XFRZ_CV(:) /PTSTEP ! kg m-2 s-1
  !
  ! Prognostic Updates:
  PEK%XWRVN(:)       = ZWRVN(:)
  PEK%XTV(:)         = PEK%XTV(:) + (DEK%XFRZ_CV(:) - DEK%XMELT_CV(:))*(XLMTT*PTSTEP)/PCHEATV(:) ! K
  !
END WHERE
!
IF (LHOOK) CALL DR_HOOK('SNOW_LOAD_MEB',1,ZHOOK_HANDLE)
!
END SUBROUTINE SNOW_LOAD_MEB
