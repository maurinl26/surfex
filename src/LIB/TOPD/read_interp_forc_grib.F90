!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!     #######################
      SUBROUTINE READ_INTERP_FORC_GRIB(YSC,HFILE,PFIELD_OUT,KNUM_GRIB,KLTYPE_GRIB,KLEV1,OINTERP)
!     #######################
!
!!****  *READ_INTERP_FORC_GRIB*  
!!
!!    PURPOSE
!!    -------
!!     This routine aims at reading forcing variables from grib files
!!     and interpolates values on SURFEX domain 
!!**  METHOD
!!    ------
!
!!    EXTERNAL
!!    --------
!!
!!    none
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------ 
!!
!!    
!!    
!!
!!      
!!    REFERENCE
!!    ---------
!!
!!    
!!      
!!    AUTHOR
!!    ------
!!
!!      B. Vincendon    * Meteo-France *
!!
!!    MODIFICATIONS
!!    -------------
!!
!!      Original   11/2007
!-------------------------------------------------------------------------------
!
!*       0.     DECLARATIONS
!               ------------
USE MODD_SURFEX_n,   ONLY : SURFEX_t
USE MODD_PREP,       ONLY : CINGRID_TYPE, CINTERP_TYPE, LINTERP
USE MODD_GRID_GRIB,  ONLY : NNI, CINMODEL
USE MODD_TYPE_DATE_SURF
USE MODD_SURF_PAR,   ONLY : XUNDEF
!
USE MODI_GET_LUOUT
USE MODI_OPEN_FILE
USE MODI_CLOSE_FILE
USE MODE_READ_GRIB
USE MODI_HOR_INTERPOL
USE MODI_PREP_GRIB_GRID
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
TYPE(SURFEX_t), INTENT(INOUT)     :: YSC
!
CHARACTER(LEN=*),   INTENT(IN)    :: HFILE       ! Grib file to be read
REAL, DIMENSION(:), INTENT(OUT)   :: PFIELD_OUT  ! Forcing parameter read on file
INTEGER,            INTENT(IN)    :: KNUM_GRIB
INTEGER,            INTENT(INOUT) :: KLTYPE_GRIB
INTEGER,            INTENT(INOUT) :: KLEV1
LOGICAL,            INTENT(IN)    :: OINTERP
!!!
!*      0.2    declarations of local variables
!
INTEGER                   :: IUNIT       ! Unit of the files
INTEGER                   :: ILUOUT      ! Unit of the files
INTEGER                   :: IRESP    ! Return of reading
TYPE (DATE_TIME)          :: TZTIME_GRIB    ! current date and time
CHARACTER(LEN=6)          :: YINMODEL       ! model from which GRIB file originates
REAL, DIMENSION(:),   ALLOCATABLE         :: ZFLD_LECT  ! Field read
REAL, DIMENSION(:,:), ALLOCATABLE         :: ZFLD_INTERPOL !Field interpolated
REAL, DIMENSION(:),   POINTER :: ZPNT_LECT
REAL, DIMENSION(:,:), POINTER :: ZPNT_INTERPOL
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('READ_INTERP_FORC_GRIB',0,ZHOOK_HANDLE)
!
!
!*       1.0    preparing file openning
!               ----------------------
CALL PREP_GRIB_GRID(HFILE,IUNIT,YINMODEL,CINGRID_TYPE,CINTERP_TYPE,TZTIME_GRIB)
CINTERP_TYPE='HORIBL'
!
!
ALLOCATE(ZPNT_LECT(NNI))
ALLOCATE(ZFLD_LECT(NNI))
ZFLD_LECT(:)=XUNDEF
!
!*       2.0    reading file
!               ----------------------
CALL READ_GRIB(HFILE,CINMODEL,IUNIT,KNUM_GRIB,IRESP,PFIELD=ZPNT_LECT,&
               KLTYPE=KLTYPE_GRIB,KLEV1=KLEV1)

IF (IRESP /=0) THEN
  WRITE(*,*) 'WARNING'
  WRITE(*,*) '-------'
  WRITE(*,*) 'error when reading article ', KNUM_GRIB,'IRESP=',IRESP
  WRITE(*,*) 'default value may be used, who knows???'
  WRITE(*,*) ' '
  ZPNT_LECT(:)=XUNDEF
ENDIF
ZFLD_LECT(:)=ZPNT_LECT(:)
!
!
!*       3.0    interpolation on good grid
!               ----------------------
IF (OINTERP) THEN
!
  IF (ALL(ZFLD_LECT(:)==XUNDEF)) THEN
    PFIELD_OUT(:)=XUNDEF
  ELSE
    PFIELD_OUT(:)=XUNDEF
    ALLOCATE(ZPNT_INTERPOL(SIZE(ZFLD_LECT,1),1))
    ZPNT_INTERPOL(:,1)=ZFLD_LECT(:)
    DEALLOCATE(ZFLD_LECT)
    !
    ALLOCATE(ZFLD_INTERPOL(YSC%U%NDIM_FULL,1))
    ZFLD_INTERPOL(:,:)=XUNDEF 
    !
    CALL GET_LUOUT('ASCII ',ILUOUT)
    CALL HOR_INTERPOL(YSC%DTCO,YSC%U,YSC%GCP,ILUOUT,ZPNT_INTERPOL,ZFLD_INTERPOL)
    PFIELD_OUT(:)=ZFLD_INTERPOL(:,1)
  ENDIF
  !
ELSE
  PFIELD_OUT(:)=XUNDEF
  PFIELD_OUT(:)=ZFLD_LECT(:)
ENDIF
!
DEALLOCATE(ZPNT_LECT)
!
IF (LHOOK) CALL DR_HOOK('READ_INTERP_FORC_GRIB',1,ZHOOK_HANDLE)
!
END SUBROUTINE READ_INTERP_FORC_GRIB
