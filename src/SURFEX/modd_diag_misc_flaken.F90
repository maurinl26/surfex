!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!     ############################
      MODULE MODD_DIAG_MISC_FLAKE_n
!     ############################
!
!!****  *MODD_DIAG_MISC_FLAKE - declaration of diagnostic variables for FLAKE scheme
!!
!!    PURPOSE
!!    -------
!
!!
!!**  IMPLICIT ARGUMENTS
!!    ------------------
!!      None 
!!
!!    REFERENCE
!!    ---------
!!
!!    AUTHOR
!!    ------
!!      P. Le Moigne   *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original     10/2004
!!      P. Le Moigne 05/2023 : temperature profile in sediments
!!      P. Le Moigne 08/2023 : heat storage flux
!!      P. Le Moigne 11/2023 : effective rain P-E
!------------------------------------------------------------------------------------
!
!*       0.   DECLARATIONS
!             ------------
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE

TYPE DIAG_MISC_FLAKE_t
!------------------------------------------------------------------------------
!
  LOGICAL :: LWATER_PROFILE   ! flag for water diagnostics
  LOGICAL :: LSEDIM_PROFILE   ! flag for sediment diagnostics
  LOGICAL :: LFLKFLUX         ! flag for heat and radiative diagnostics
  LOGICAL :: LFLKWATER        ! flag for water budget P-E
!
!* miscellaneous variables
!
  REAL, POINTER, DIMENSION(:)   :: XZWAT_PROFILE ! depth of output levels (m) in namelist
  REAL, POINTER, DIMENSION(:)   :: XZWPROF       ! depth of output levels (m)
  REAL, POINTER, DIMENSION(:,:) :: XTWPROF       ! Water temperature in output levels (K)
  
  REAL, POINTER, DIMENSION(:)   :: XZSED_PROFILE ! depth of output levels (m) in namelist
  REAL, POINTER, DIMENSION(:)   :: XZSPROF       ! depth of output levels in the sediment layer (m)
  REAL, POINTER, DIMENSION(:,:) :: XTSPROF       ! Temperature in the sediment layer in output levels (K)
!
  REAL, POINTER, DIMENSION(:)   :: XQBOT_FLK     ! Heat flux through the water-bottom sediment interface [W m^{-2}]
  REAL, POINTER, DIMENSION(:)   :: XQW_FLK       ! Heat flux through the ice-water or air-water interface [W m^{-2}]
  REAL, POINTER, DIMENSION(:)   :: XQSNOW_FLK    ! Heat flux through the air-snow interface [W m^{-2}]
  REAL, POINTER, DIMENSION(:)   :: XQICE_FLK     ! Heat flux through the snow-ice or air-ice interface [W m^{-2}]
  REAL, POINTER, DIMENSION(:)   :: XIATM_FLK     ! Radiation flux at the lower boundary of the atmosphere [W m^{-2}]
  REAL, POINTER, DIMENSION(:)   :: XISNOW_FLK    ! Radiation flux through the air-snow interface [W m^{-2}]
  REAL, POINTER, DIMENSION(:)   :: XIICE_FLK     ! Radiation flux through the snow-ice or air-ice interface [W m^{-2}]
  REAL, POINTER, DIMENSION(:)   :: XIW_FLK       ! Radiation flux through the ice-water or air-water interface [W m^{-2}]
  REAL, POINTER, DIMENSION(:)   :: XIH_FLK       ! Radiation flux through the mixed-layer-thermocline interface [W m^{-2}]
  REAL, POINTER, DIMENSION(:)   :: XIBOT_FLK     ! Radiation flux through the water-bottom sediment interface [W m^{-2}]
  REAL, POINTER, DIMENSION(:)   :: XIINTM0H_FLK  ! Mean radiation flux over the mixed layer [W m^{-1}]
  REAL, POINTER, DIMENSION(:)   :: XIINTMHD_FLK  ! Mean radiation flux over the thermocline [W m^{-1}]
  REAL, POINTER, DIMENSION(:)   :: XQSTAR_FLK    ! A generalized heat flux scale [W m^{-2}]
  REAL, POINTER, DIMENSION(:)   :: XQSTO_FLK     ! Heat storage flux [W m^{-2}]
  !
  REAL, POINTER, DIMENSION(:)   :: XQSTOC_FLK    ! Accumulated heat storage flux [J m^{-2}]
  REAL, POINTER, DIMENSION(:)   :: XIBOTC_FLK    ! Accumulated radiation flux through the water-bottom sediment interface [J m^{-2}]
  REAL, POINTER, DIMENSION(:)   :: XQBOTC_FLK    ! Accumulated heat flux through the water-bottom sediment interface [J m^{-2}]
  REAL, POINTER, DIMENSION(:)   :: XIWC_FLK      ! Accumulated radiation flux through the ice-water or air-water interface [J m^{-2}]
  REAL, POINTER, DIMENSION(:)   :: XQWC_FLK      ! Accumulated heat flux through the ice-water or air-water interface [J m^{-2}]
!
  REAL, POINTER, DIMENSION(:)   :: XPME          ! Effective rain P-E (kg/m2/s)
  REAL, POINTER, DIMENSION(:)   :: XPMEC         ! Accumulated effective rain P-E (kg/m2)
  REAL, POINTER, DIMENSION(:)   :: XPRE          ! Total precipitation (kg/m2/s)
  REAL, POINTER, DIMENSION(:)   :: XPREC         ! Accumulated precipitation P-E (kg/m2)
!------------------------------------------------------------------------------
!
END TYPE DIAG_MISC_FLAKE_t

CONTAINS
!
SUBROUTINE DIAG_MISC_FLAKE_INIT(DMF)
TYPE(DIAG_MISC_FLAKE_t), INTENT(INOUT) :: DMF
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK("MODD_DIAG_MISC_FLAKE_N:DIAG_MISC_FLAKE_INIT",0,ZHOOK_HANDLE)
!
NULLIFY(DMF%XZWAT_PROFILE)
NULLIFY(DMF%XZWPROF)
NULLIFY(DMF%XTWPROF)
NULLIFY(DMF%XZSED_PROFILE)
NULLIFY(DMF%XZSPROF)
NULLIFY(DMF%XTSPROF)
NULLIFY(DMF%XQBOT_FLK)
NULLIFY(DMF%XQW_FLK)  
NULLIFY(DMF%XQSNOW_FLK)  
NULLIFY(DMF%XQICE_FLK)   
NULLIFY(DMF%XIATM_FLK)   
NULLIFY(DMF%XISNOW_FLK)  
NULLIFY(DMF%XIICE_FLK)   
NULLIFY(DMF%XIW_FLK)     
NULLIFY(DMF%XIH_FLK)     
NULLIFY(DMF%XIBOT_FLK)     
NULLIFY(DMF%XIINTM0H_FLK)
NULLIFY(DMF%XIINTMHD_FLK)
NULLIFY(DMF%XQSTAR_FLK)  
NULLIFY(DMF%XQSTO_FLK)  
NULLIFY(DMF%XPME)  
NULLIFY(DMF%XPMEC)  
NULLIFY(DMF%XPRE)  
NULLIFY(DMF%XPREC)  
NULLIFY(DMF%XQSTOC_FLK)
NULLIFY(DMF%XIBOTC_FLK)
NULLIFY(DMF%XQBOTC_FLK)
NULLIFY(DMF%XIWC_FLK)
NULLIFY(DMF%XQWC_FLK)
!
DMF%LWATER_PROFILE=.FALSE.
DMF%LSEDIM_PROFILE=.FALSE.
DMF%LFLKFLUX=.FALSE.
DMF%LFLKWATER=.FALSE.
!
IF (LHOOK) CALL DR_HOOK("MODD_DIAG_MISC_FLAKE_N:DIAG_MISC_FLAKE_INIT",1,ZHOOK_HANDLE)
END SUBROUTINE DIAG_MISC_FLAKE_INIT
!
END MODULE MODD_DIAG_MISC_FLAKE_n
