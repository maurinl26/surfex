!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!     #############################################################
      SUBROUTINE WRITE_HEADER_MNH
!     #############################################################
!
!!****  * - routine to header-type fields in a lfi file to emulate a MesoNH file
!!
!!    PURPOSE
!!    -------
!
!
!!**  METHOD
!!    ------
!!
!!    EXTERNAL
!!    --------
!!
!!     
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!
!!    REFERENCE
!!    ---------
!!
!!
!!    AUTHOR
!!    ------
!!
!!      V. Masson      *METEO-FRANCE*
!!
!!    MODIFICATIONS
!!    -------------
!!
!!      original         21/05/08
!!      Y. Seity          09/2018: add 'FA' case                                  
!!      A.Druel           08/2019: Permit to change the size of caracters (write / read) with constant
!!
!----------------------------------------------------------------------------
!
!*      0.    DECLARATIONS
!             ------------
!
#ifdef SFX_LFI
USE MODI_FMWRIT
#endif
!
USE MODD_IO_SURF_LFI,        ONLY : CFILEOUT_LFI, CLUOUT_LFI, LMNH_COMPATIBLE, LCARTESIAN

#ifdef SFX_FA
USE MODE_WRITE_SURF_FA
#endif
USE MODN_IO_OFFLINE
USE MODD_IO_SURF_FA
!
USE MODD_DATA_COVER_PAR,     ONLY : NCAR_FILES
USE MODD_SURF_PAR, ONLY : LEN_HREC
!
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*      0.1   Declarations of arguments
!
!*      0.2   Declarations of local variables
!
 CHARACTER(LEN=100)        :: YCOMMENT=' '
INTEGER                    :: IRESP
INTEGER                    :: INB ! number of articles in the file
 CHARACTER(LEN=28)         :: YNAME
 CHARACTER(LEN=10)         :: YBIBUSER =' '
 CHARACTER(LEN=LEN_HREC)         :: YREC
 CHARACTER(LEN=NCAR_FILES) :: YFIELD
INTEGER :: J
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!----------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('WRITE_HEADER_MNH',0,ZHOOK_HANDLE)
  SELECT CASE (CSURF_FILETYPE)
    CASE ('LFI   ')
#ifdef SFX_LFI
      CALL FMOPEN(CFILEOUT_LFI,'UNKNOWN',CLUOUT_LFI,0,1,1,INB,IRESP)
!
      CALL FMWRITN0(CFILEOUT_LFI,'MASDEV',CLUOUT_LFI,1,47,4,100,YCOMMENT,IRESP)
      CALL FMWRITN0(CFILEOUT_LFI,'BUGFIX',CLUOUT_LFI,1,0,4,100,YCOMMENT,IRESP)
      CALL FMWRITC0(CFILEOUT_LFI,'BIBUSER',CLUOUT_LFI,1,YBIBUSER,4,100,YCOMMENT,IRESP)
     YNAME=CFILEOUT_LFI
      CALL FMWRITC0(CFILEOUT_LFI,'MY_NAME',CLUOUT_LFI,1,YNAME,4,100,YCOMMENT,IRESP)
     YNAME=' '
      CALL FMWRITC0(CFILEOUT_LFI,'DAD_NAME',CLUOUT_LFI,1,YNAME,4,100,YCOMMENT,IRESP)
      CALL FMWRITC0(CFILEOUT_LFI,'PROGRAM',CLUOUT_LFI,1,'SURFEX',4,100,YCOMMENT,IRESP)
      CALL FMWRITN0(CFILEOUT_LFI,'KMAX',CLUOUT_LFI,1,0,4,100,YCOMMENT,IRESP)
      CALL FMWRITC0(CFILEOUT_LFI,'STORAGE_TYPE',CLUOUT_LFI,1,'SU    ',4,100,YCOMMENT,IRESP)
      CALL FMWRITL0(CFILEOUT_LFI,'CARTESIAN       ',CLUOUT_LFI,1,LCARTESIAN,4,100,YCOMMENT,IRESP)
      CALL FMWRITL0(CFILEOUT_LFI,'THINSHELL       ',CLUOUT_LFI,1,.TRUE.,4,100,YCOMMENT,IRESP)
!
      CALL FMCLOS(CFILEOUT_LFI,'KEEP',CLUOUT_LFI,IRESP)
#endif
#ifdef SFX_FA
    CASE ('FA    ')
      YREC =       'MASDEV';                     ; CALL WRITE_SURF0_FA (YREC, 47,         IRESP, YCOMMENT)
      YREC =       'BUGFIX';                     ; CALL WRITE_SURF0_FA (YREC, 0,          IRESP, YCOMMENT)
      YREC =      'BIBUSER'; YFIELD = YBIBUSER   ; CALL WRITE_SURF0_FA (YREC, YFIELD,     IRESP, YCOMMENT)
      YFIELD = CFILEOUT_FA;
      J = LEN (TRIM (YFIELD))
      IF (J > 3) THEN
        IF (YFIELD (J-2:J) == '.fa') YFIELD (J-2:J) = '   '
      ENDIF
      YREC =      'MY_NAME';                       CALL WRITE_SURF0_FA (YREC, YFIELD,     IRESP, YCOMMENT)
      YREC =     'DAD_NAME'; YFIELD = ' '        ; CALL WRITE_SURF0_FA (YREC, YFIELD,     IRESP, YCOMMENT)
      YREC =      'PROGRAM'; YFIELD = 'SURFEX'   ; CALL WRITE_SURF0_FA (YREC, YFIELD,     IRESP, YCOMMENT)
      YREC =         'KMAX';                     ; CALL WRITE_SURF0_FA (YREC, 0,          IRESP, YCOMMENT)
      YREC = 'STORAGE_TYPE'; YFIELD = 'SU    '   ; CALL WRITE_SURF0_FA (YREC, YFIELD,     IRESP, YCOMMENT)
      YREC =    'CARTESIAN';                     ; CALL WRITE_SURF0_FA (YREC, LCARTESIAN, IRESP, YCOMMENT)
      YREC =    'THINSHELL';                     ; CALL WRITE_SURF0_FA (YREC, .TRUE.,     IRESP, YCOMMENT)
#endif
  END SELECT

IF (LHOOK) CALL DR_HOOK('WRITE_HEADER_MNH',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
END SUBROUTINE WRITE_HEADER_MNH

