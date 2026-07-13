!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!     #########
SUBROUTINE PREP_ISBA_UNIF(DTCO, U, IO, KLUOUT,HSURF,PFIELD)
!     #################################################################################
!
!!****  *PREP_ISBA_UNIF* - prepares ISBA field from prescribed values
!!
!!    PURPOSE
!!    -------
!
!!**  METHOD
!!    ------
!!
!!    REFERENCE
!!    ---------
!!      
!!
!!    AUTHOR
!!    ------
!!     V. Masson 
!!
!!    MODIFICATIONS
!!    -------------
!!      Original     01/2004
!!      P Samuelsson 02/2012 MEB
!!      P Le Moigne  08/2023 Initialization of soil temperature and moisture profiles
!!------------------------------------------------------------------
!
!
USE MODD_DATA_COVER_n,   ONLY : DATA_COVER_t
USE MODD_SURF_ATM_n,     ONLY : SURF_ATM_t
USE MODD_ISBA_OPTIONS_n, ONLY : ISBA_OPTIONS_t
!
USE MODD_PREP,           ONLY : CINTERP_TYPE
USE MODD_SURF_PAR,       ONLY : XUNDEF
USE MODD_PREP_ISBA,      ONLY : XUNIF_HUG_SOIL, XUNIF_HUGI_SOIL, XUNIF_TG_SOIL,    &
                                XWR_DEF, XWRV_DEF, XWRVN_DEF, XQC_DEF
!
USE MODI_ABOR1_SFX
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
TYPE(DATA_COVER_t), INTENT(INOUT) :: DTCO
TYPE(SURF_ATM_t), INTENT(INOUT)   :: U
TYPE(ISBA_OPTIONS_t), INTENT(INOUT) :: IO
!
INTEGER,            INTENT(IN)  :: KLUOUT    ! output listing logical unit
 CHARACTER(LEN=7),   INTENT(IN) :: HSURF     ! type of field
REAL, POINTER, DIMENSION(:,:,:) :: PFIELD    ! field to interpolate horizontally
!
!*      0.2    declarations of local variables
!
INTEGER :: IL, JLAYER
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('PREP_ISBA_UNIF',0,ZHOOK_HANDLE)
!
 CALL GET_TYPE_DIM_n(DTCO, U, 'NATURE',IL)
!
SELECT CASE(HSURF)
!
!*      3.0    Orography
!
  CASE('ZS     ')
    ALLOCATE(PFIELD(1,1,1))
    PFIELD = 0.
   
!
!*      3.1    Profile of soil relative humidity
!
  CASE('WG     ') 
    IF(ANY(XUNIF_HUG_SOIL(1:IO%NGROUND_LAYER)==XUNDEF))THEN
      CALL ABOR1_SFX('PREP_ISBA_UNIF: No values for '//TRIM(HSURF)//' check your namelist and NAM_PREP_ISBA !')
    ENDIF          
    ALLOCATE(PFIELD(IL,IO%NGROUND_LAYER,1))
    DO JLAYER=1,IO%NGROUND_LAYER
      PFIELD(:,JLAYER,1) = XUNIF_HUG_SOIL(JLAYER)
    ENDDO

!*      3.2    Profile of soil humidity for ice

  CASE('WGI    ')
    ALLOCATE(PFIELD(IL,IO%NGROUND_LAYER,1))
    DO JLAYER=1,IO%NGROUND_LAYER
      PFIELD(:,JLAYER,1) = XUNIF_HUGI_SOIL(JLAYER)
    ENDDO

!*      3.3    Profile of temperatures

  CASE('TG     ')
    IF(ANY(XUNIF_TG_SOIL(1:IO%NGROUND_LAYER)==XUNDEF))THEN
      CALL ABOR1_SFX('PREP_ISBA_UNIF: No values for '//TRIM(HSURF)//' check your namelist and NAM_PREP_ISBA !')
    ENDIF
    ALLOCATE(PFIELD(IL,IO%NGROUND_LAYER,1))
    DO JLAYER=1,IO%NGROUND_LAYER
      PFIELD(:,JLAYER,1) = XUNIF_TG_SOIL(JLAYER)
    ENDDO

!*      3.4    Other quantities

  CASE('WR     ')
    ALLOCATE(PFIELD(IL,1,1))
    PFIELD = XWR_DEF

  CASE('WRL    ')
    ALLOCATE(PFIELD(IL,1,1))
    PFIELD = XWRV_DEF

  CASE('WRLI    ')
    ALLOCATE(PFIELD(IL,1,1))
    PFIELD = XWRV_DEF

  CASE('WRVN   ')
    ALLOCATE(PFIELD(IL,1,1))
    PFIELD = XWRVN_DEF

  CASE('TV     ')
    ALLOCATE(PFIELD(IL,1,1))
    PFIELD = XUNIF_TG_SOIL(1)

  CASE('TL     ')
    ALLOCATE(PFIELD(IL,1,1))
    PFIELD = XUNIF_TG_SOIL(1)

  CASE('TC     ')
    ALLOCATE(PFIELD(IL,1,1))
    PFIELD = XUNIF_TG_SOIL(1)

  CASE('QC     ')
    ALLOCATE(PFIELD(IL,1,1))
    PFIELD = XQC_DEF

  CASE('LAI    ')
    ALLOCATE(PFIELD(IL,1,1))
    PFIELD = XUNDEF

  CASE('ICE_STO')
    ALLOCATE(PFIELD(IL,1,1))
    PFIELD = 0.0
!
  CASE DEFAULT
    CALL ABOR1_SFX('PREP_ISBA_UNIF: '//TRIM(HSURF)//" initialization not implemented !")
!
END SELECT
!
!*      4.     Interpolation method
!              --------------------
!
CINTERP_TYPE='UNIF  '
!
IF (LHOOK) CALL DR_HOOK('PREP_ISBA_UNIF',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------------
END SUBROUTINE PREP_ISBA_UNIF
