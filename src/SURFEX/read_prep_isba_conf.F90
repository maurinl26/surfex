!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!     #########
      SUBROUTINE READ_PREP_ISBA_CONF(HPROGRAM,HVAR,HFILE,HFILETYPE,HFILEPGD,     &
                                     HFILEPGDTYPE,HATMFILE,HATMFILETYPE,HPGDFILE,&
                                     HPGDFILETYPE,KLUOUT,OUNIF,KPATCH,KGROUND_LAYER)
!     #######################################################
!
!!****  *READ_PREP_ISBA_CONF* - routine to read the configuration for ISBA
!!                              fields preparation
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
!!      P Samuelsson 02/2012  MEB
!!      P Le Moigne  08/2023  Initialization of soil temperature and moisture profiles
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODI_READ_PREP_SURF_ATM_CONF
USE MODI_READ_SURF
!
USE MODN_PREP_ISBA
USE MODD_PREP_ISBA,  ONLY : CFILE_ISBA, CTYPE, CFILEPGD_ISBA, CTYPEPGD,        &
                            CFILE_HUG, CTYPE_HUG, CFILE_TG, CTYPE_TG,          &
                            XUNIF_HUG_SOIL, XUNIF_HUGI_SOIL, XUNIF_TG_SOIL,    &
                            XWSNOW, XTSNOW, XRSNOW, XASNOW,                    &
                            CFILE_HUG_SOIL, CFILE_TG_SOIL
!
USE MODD_SURF_PAR,   ONLY : XUNDEF
!
USE MODD_SURF_PAR,        ONLY : XUNDEF
USE MODD_DATA_COVER_PAR,  ONLY : NVEGTYPE
USE MODD_ISBA_OPTIONS_n, ONLY : ISBA_OPTIONS_t
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
USE MODI_ABOR1_SFX
!
IMPLICIT NONE
!
!*       0.1   Declarations of arguments
!              -------------------------
!
CHARACTER(LEN=6),  INTENT(IN)  :: HPROGRAM    ! program calling ISBA
CHARACTER(LEN=7),  INTENT(IN)  :: HVAR        ! variable treated
CHARACTER(LEN=28), INTENT(OUT) :: HFILE       ! file name
CHARACTER(LEN=6),  INTENT(OUT) :: HFILETYPE   ! file type
CHARACTER(LEN=28), INTENT(OUT) :: HFILEPGD    ! file name
CHARACTER(LEN=6),  INTENT(OUT) :: HFILEPGDTYPE! file type
CHARACTER(LEN=28), INTENT(IN)  :: HATMFILE    ! atmospheric file name
CHARACTER(LEN=6),  INTENT(IN)  :: HATMFILETYPE! atmospheric file type
CHARACTER(LEN=28), INTENT(IN)  :: HPGDFILE    ! atmospheric file name
CHARACTER(LEN=6),  INTENT(IN)  :: HPGDFILETYPE! atmospheric file type
INTEGER,           INTENT(IN)  :: KLUOUT      ! logical unit of output listing
LOGICAL,           INTENT(OUT) :: OUNIF       ! flag for prescribed uniform field

INTEGER,  INTENT(IN), OPTIONAL :: KPATCH      ! number of patches
INTEGER,  INTENT(IN), OPTIONAL :: KGROUND_LAYER      ! number of soil layers
!
!*       0.2   Declarations of local variables
!              -------------------------------
!
INTEGER :: JLAYER
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
!
!
IF (LHOOK) CALL DR_HOOK('READ_PREP_ISBA_CONF',0,ZHOOK_HANDLE)
HFILE    = '                         '
HFILETYPE    = '      '
!
HFILEPGD = '                         '
HFILEPGDTYPE = '      '
!
OUNIF    = .FALSE.
!
!-------------------------------------------------------------------------------
!
!* choice of input file
!  --------------------
!
SELECT CASE (HVAR)
  CASE ('WG     ','WGI    ')
    IF (LEN_TRIM(CFILE_HUG)>0 .AND. LEN_TRIM(CTYPE_HUG)>0 ) THEN
      HFILE     = CFILE_HUG
      HFILETYPE = CTYPE_HUG
    END IF
  CASE ('TG     ','TV     ','TC     ','TL     ')
    IF (LEN_TRIM(CFILE_TG)>0 .AND. LEN_TRIM(CTYPE_TG)>0 ) THEN
      HFILE     = CFILE_TG
      HFILETYPE = CTYPE_TG
    END IF
END SELECT
!
IF (LEN_TRIM(HFILE)==0 .AND. LEN_TRIM(CFILE_ISBA)>0 .AND. LEN_TRIM(CTYPE)>0) THEN
  HFILE     = CFILE_ISBA
  HFILETYPE = CTYPE
END IF
!
IF (LEN_TRIM(HFILEPGD)==0 .AND. LEN_TRIM(CFILEPGD_ISBA)>0 .AND. LEN_TRIM(CTYPEPGD)>0) THEN
  HFILEPGD     = CFILEPGD_ISBA
  HFILEPGDTYPE = CTYPEPGD
END IF
!
!! If no file name in the scheme namelist,
!! try to find a name in NAM_SURF_ATM
!
IF (LEN_TRIM(HFILE)==0) THEN
!
 CALL READ_PREP_SURF_ATM_CONF(HPROGRAM,HFILE,HFILETYPE,HFILEPGD,HFILEPGDTYPE,&
                             HATMFILE,HATMFILETYPE,HPGDFILE,HPGDFILETYPE,KLUOUT)
!
END IF
!
!! If no file name in the scheme namelist,
!! nor in NAM_SURF_ATM, look if ascii input files are present
!
SELECT CASE (HVAR)
  CASE ('WG     ','WGI    ')
    IF ( LEN_TRIM(CTYPE_HUG )>0 .AND. ALL(LEN_TRIM(CFILE_HUG_SOIL(1:KPATCH,1:KGROUND_LAYER))>0) ) THEN  
       HFILETYPE = CTYPE_HUG 
    END IF
    IF (HVAR=='WGI    ' .AND. HFILETYPE=='ASCLLV') THEN
       OUNIF = .TRUE.
       DO JLAYER=1,KGROUND_LAYER
         IF (XUNIF_HUGI_SOIL(JLAYER)==XUNDEF) XUNIF_HUGI_SOIL(JLAYER) = 0.
         IF (LHOOK) CALL DR_HOOK('READ_PREP_ISBA_CONF',1,ZHOOK_HANDLE)
       ENDDO
       RETURN
    ENDIF
  CASE ('TG     ','TV     ','TC     ','TL     ')
    IF ( LEN_TRIM(CTYPE_TG )>0 .AND. ALL(LEN_TRIM(CFILE_TG_SOIL(1:KPATCH,1:KGROUND_LAYER))>0) ) THEN  
       HFILETYPE = CTYPE_TG 
    END IF
END SELECT
!
!-------------------------------------------------------------------------------
!
!* Is an uniform field prescribed?
!  ------------------------------
!
SELECT CASE (HVAR)
  CASE ('WG     ')
    OUNIF = ANY(XUNIF_HUG_SOIL(:)/=XUNDEF)
    DO JLAYER=1,KGROUND_LAYER
      IF (OUNIF .AND. (XUNIF_HUG_SOIL(JLAYER)==XUNDEF)) THEN
         WRITE(KLUOUT,*)'ONE OF XUNIF_HUG_SOIL LAYER IS GIVEN'
         CALL ABOR1_SFX('READ_PREP_ISBA_CONF: XUNIF_HUG_SOIL MUST BE SET FOR ALL LAYERS')
      END IF
    ENDDO
    !
  CASE ('WGI    ')
    OUNIF = ANY(XUNIF_HUGI_SOIL(:)/=XUNDEF)
    DO JLAYER=1,KGROUND_LAYER
      IF (OUNIF .AND. (XUNIF_HUGI_SOIL(JLAYER)==XUNDEF)) THEN
         WRITE(KLUOUT,*)'ONE OF XUNIF_HUGI_SOIL LAYER IS GIVEN'
         CALL ABOR1_SFX('READ_PREP_ISBA_CONF: XUNIF_HUGI_SOIL MUST BE SET FOR ALL LAYERS')
      END IF
    ENDDO
   !
  CASE ('TG     ','TV     ','TC     ')
    OUNIF = ANY(XUNIF_TG_SOIL(:)/=XUNDEF)
    DO JLAYER=1,KGROUND_LAYER
      IF (OUNIF .AND. (XUNIF_TG_SOIL(JLAYER)==XUNDEF)) THEN
         WRITE(KLUOUT,*)'ONE OF XUNIF_TG_SOIL LAYER IS GIVEN'
         CALL ABOR1_SFX('READ_PREP_ISBA_CONF: XUNIF_TG_SOIL MUST BE SET FOR ALL LAYERS')
      END IF
    ENDDO
    !
END SELECT
!
!-------------------------------------------------------------------------------
!
!* no file given ? nor specific value in namelist? One takes the default value.
!
IF (HFILETYPE=='      ' .AND. .NOT. OUNIF) THEN
  IF (HVAR(1:2)/='ZS') WRITE(KLUOUT,*) 'NO FILE FOR FIELD ',HVAR, &
                                        ': UNIFORM DEFAULT FIELD IS PRESCRIBED'
  IF (HVAR(1:3)=='WGI') THEN
    XUNIF_HUGI_SOIL(:) = 0.
  ENDIF                                     
  OUNIF = .TRUE.
END IF
IF (LHOOK) CALL DR_HOOK('READ_PREP_ISBA_CONF',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE READ_PREP_ISBA_CONF
