!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!     #######################
      SUBROUTINE READ_INTERP_FORC_SAF(YSC,HFILE,PFIELD_OUT,KDATE_DEB,KNUM_GRIB,KLTYPE_GRIB,KLEV1,KNB_HOUR_SAF)
!     #######################
!
!!****  *READ_INTERP_FORC_SAFRAN*  
!!
!!    PURPOSE
!!    -------
!!     This routine aims at reading forcing variables from grib files
!!     and interpolates values on SURFEX domain
!
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
USE MODD_SURFEX_n, ONLY : SURFEX_t
USE MODD_PREP,       ONLY : CINGRID_TYPE, CINTERP_TYPE
USE MODD_GRID_GRIB,  ONLY : NNI
USE MODD_TYPE_DATE_SURF
USE MODD_SURF_PAR,   ONLY : XUNDEF
!
USE MODI_GET_LUOUT
USE MODI_OPEN_FILE
USE MODI_CLOSE_FILE
USE MODI_READ_GRIB_SAF
USE MODI_HOR_INTERPOL
USE MODI_PREP_GRIB_GRID
!
USE YOMHOOK   ,ONLY : LHOOK, DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
TYPE(SURFEX_t), INTENT(INOUT)                   :: YSC
!
CHARACTER(LEN=*),                 INTENT(IN)    :: HFILE       ! Grib file to be read
REAL, DIMENSION(:,:),ALLOCATABLE, INTENT(OUT)   :: PFIELD_OUT  ! Forcing parameter read on file
INTEGER,                          INTENT(IN)    :: KDATE_DEB
INTEGER,                          INTENT(IN)    :: KNUM_GRIB
INTEGER,                          INTENT(INOUT) :: KLTYPE_GRIB
INTEGER,                          INTENT(INOUT) :: KLEV1
INTEGER,                          INTENT(IN)    :: KNB_HOUR_SAF
!!
!*      0.2    declarations of local variables
!
!

INTEGER                   :: JWRK           ! loop control
INTEGER                   :: IUNIT          ! Unit of the files
INTEGER                   :: ILUOUT         ! Unit of the files
INTEGER                   :: IRESP          ! Return of reading
TYPE (DATE_TIME)          :: TZTIME_GRIB    ! current date and time
CHARACTER(LEN=6)          :: YINMODEL       ! model from which GRIB file originates
REAL, DIMENSION(:),   ALLOCATABLE     :: ZFLD_LECT     ! Field read
REAL, DIMENSION(:,:), ALLOCATABLE     :: ZFLD_INTERPOL !Field interpolated
REAL, DIMENSION(:),   POINTER         :: ZPNT_LECT
REAL,   POINTER,   DIMENSION(:,:)     :: ZPNT_INTERPOL
REAL,   POINTER,   DIMENSION(:,:)     :: ZPNT_LECT_2D !Field not yet interpolated
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('READ_INTERP_FORC_SAF',0,ZHOOK_HANDLE)
!
!
!*       1.0    preparing file openning
!               ----------------------
CALL PREP_GRIB_GRID(HFILE,IUNIT,YINMODEL,CINGRID_TYPE,CINTERP_TYPE,TZTIME_GRIB)
CINTERP_TYPE='HORIBL'
!
ALLOCATE (PFIELD_OUT(KNB_HOUR_SAF,YSC%U%NDIM_FULL))
PFIELD_OUT(:,:)=XUNDEF
!
ALLOCATE(ZPNT_LECT(NNI))
ALLOCATE(ZPNT_LECT_2D(KNB_HOUR_SAF,NNI))
!
!*       2.0    reading file
!               ----------------------
CALL READ_GRIB_SAF(HFILE,IUNIT,KNUM_GRIB,KDATE_DEB,IRESP,KNI=NNI,PDOUT=ZPNT_LECT_2D,&
               KLTYPE=KLTYPE_GRIB,KLEV1=KLEV1,KNB_HOUR_SAF=KNB_HOUR_SAF)

IF (IRESP /=0) THEN
  WRITE(*,*) 'WARNING'
  WRITE(*,*) '-------'
  WRITE(*,*) 'error when reading article ', KNUM_GRIB,'IRESP=',IRESP
  WRITE(*,*) 'default value may be used, who knows???'
  WRITE(*,*) ' '
ENDIF
!
WHERE(ZPNT_LECT_2D<10E-10)
        ZPNT_LECT_2D=XUNDEF
ENDWHERE
!
!*       3.0    interpolation on good grid
!               ----------------------
!
ALLOCATE(ZFLD_LECT(NNI))
ALLOCATE(ZPNT_INTERPOL(SIZE(ZFLD_LECT,1),1))
ALLOCATE(ZFLD_INTERPOL(YSC%U%NDIM_FULL,1))
ZFLD_LECT(:)=XUNDEF
!
DO JWRK =1,KNB_HOUR_SAF
  ZFLD_LECT(:)=ZPNT_LECT_2D(JWRK,:)
  ZPNT_INTERPOL(:,1)=ZFLD_LECT(:)
  ZFLD_INTERPOL(:,:)=XUNDEF 
!
  CALL GET_LUOUT('ASCII ',ILUOUT)
  CALL HOR_INTERPOL(YSC%DTCO,YSC%U,YSC%GCP,ILUOUT,ZPNT_INTERPOL,ZFLD_INTERPOL)
  PFIELD_OUT(JWRK,:)=ZFLD_INTERPOL(:,1)
ENDDO
!
DEALLOCATE(ZFLD_LECT)
DEALLOCATE(ZPNT_INTERPOL)
DEALLOCATE(ZFLD_INTERPOL)
DEALLOCATE(ZPNT_LECT)
DEALLOCATE(ZPNT_LECT_2D)
!
IF (LHOOK) CALL DR_HOOK('READ_INTERP_FORC_SAF',1,ZHOOK_HANDLE)
!
END SUBROUTINE READ_INTERP_FORC_SAF
