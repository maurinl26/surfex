!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!     #########
      SUBROUTINE DEFAULT_PREP_ISBA
!     ###########################
!
!!****  *DEFAULT_PREP_ISBA* - routine to set default values for the configuration for ISBA fields preparation
!!
!!    PURPOSE
!!    -------
!!
!!**  METHOD
!!    ------
!!
!!    EXTERNAL
!!    --------
!!
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!    REFERENCE
!!    ---------
!!
!!
!!    AUTHOR
!!    ------
!!      V. Masson   *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original     01/2004 
!!      P Le Moigne  03/2007  Initialization using ASCLLV format
!!      P Samuelsson 10/2014  MEB
!!      P Le Moigne  08/2023  Initialization of soil temperature and moisture profiles
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_PREP_ISBA,  ONLY : CFILE_ISBA, CTYPE, CFILEPGD_ISBA, CTYPEPGD,     &
                            CFILE_HUG, CTYPE_HUG, CFILE_TG, CTYPE_TG,       &
                            XUNIF_HUG_SOIL, XUNIF_HUGI_SOIL, XUNIF_TG_SOIL, &
                            XWR_DEF, LEXTRAP_TG, LEXTRAP_WG, LEXTRAP_WGI,   &
                            LEXTRAP_SN, XWRV_DEF, XWRVN_DEF, XQC_DEF,       &
                            CFILE_HUG_SOIL, CFILE_TG_SOIL
!
USE MODN_PREP_ISBA
!
USE MODD_SURF_PAR,   ONLY : XUNDEF
USE MODD_SNOW_PAR,   ONLY : XANSMIN, XRHOSMAX
USE MODD_CSTS,       ONLY : XTT
!
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*       0.1   Declarations of arguments
!              -------------------------
!
!
!
!*       0.2   Declarations of local variables
!              -------------------------------
!
!-------------------------------------------------------------------------------
!

REAL(KIND=JPRB) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('DEFAULT_PREP_ISBA',0,ZHOOK_HANDLE)
CFILE_ISBA = '                          '
CTYPE      = 'GRIB  '
CFILEPGD_ISBA = '                          '
CTYPEPGD      = '      '
CFILE_HUG  = '                          '
CTYPE_HUG  = '      '
CFILE_TG   = '                          '
CTYPE_TG   = '      '
!
CFILE_HUG_SOIL(:,:)  = '                          '
CFILE_TG_SOIL (:,:)  = '                          '
!
XUNIF_HUG_SOIL(:)  = XUNDEF
XUNIF_HUGI_SOIL(:) = XUNDEF
XUNIF_TG_SOIL(:)   = XUNDEF
!
XWR_DEF   = 0.
XWRV_DEF  = 0.
XWRVN_DEF = 0.
XQC_DEF   = 0.
!
LISBA_CANOPY = .FALSE.
LEXTRAP_TG   = .FALSE.
LEXTRAP_WG   = .FALSE.
LEXTRAP_WGI  = .FALSE.
LEXTRAP_SN   = .FALSE.
 
IF (LHOOK) CALL DR_HOOK('DEFAULT_PREP_ISBA',1,ZHOOK_HANDLE)

!-------------------------------------------------------------------------------
!
END SUBROUTINE DEFAULT_PREP_ISBA
