!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!     #######################
      SUBROUTINE READ_GRIB_LARGE_RAIN(HFILE,HGRIB_TYPE,PLARGE_RAIN)
!     #######################
!
!!****  *READ_GRIB_LARGE_RAIN*  
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
USE MODD_PREP,       ONLY : CINGRID_TYPE, CINTERP_TYPE
USE MODD_GRID_GRIB,  ONLY : NNI, CINMODEL
USE MODD_TYPE_DATE_SURF
USE MODD_SURF_PAR,   ONLY : XUNDEF
!
USE MODI_GET_LUOUT
USE MODI_OPEN_FILE
USE MODI_CLOSE_FILE
USE MODI_HOR_INTERPOL
USE MODI_PREP_GRIB_GRID
USE MODI_GET_GRIB_CODE
!
USE MODE_READ_GRIB
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
CHARACTER(LEN=*),   INTENT(IN)  :: HFILE        ! Grib file to be read
CHARACTER(LEN=*),   INTENT(IN)  :: HGRIB_TYPE   ! Grib type of file to be read
REAL, DIMENSION(:), INTENT(OUT) :: PLARGE_RAIN  ! Rainfall on large grid
!
!*      0.2    declarations of local variables
!
TYPE (DATE_TIME)                :: TZTIME_GRIB    ! current date and time
INTEGER                   :: INUM_GRIB   ! number of parameter in grib file
INTEGER                   :: ITYPE_GRIB  ! number of level in grib file
INTEGER                   :: ILEV1       ! level in grib file
INTEGER                   :: IUNIT       ! Unit of the files
INTEGER                   :: IRESP    ! Return of reading
INTEGER, DIMENSION(10)    :: ITAB_NUM_GRIB        ! Code of the parameter to get
INTEGER, DIMENSION(10)    :: ITAB_TYPE_GRIB
INTEGER, DIMENSION(10)    :: ITAB_LEV1
CHARACTER(LEN=5), DIMENSION(10)  :: YREC_NAME       ! Grib type of file to be read
CHARACTER(LEN=6)                 :: YINMODEL       ! model from which GRIB file originates
REAL, DIMENSION(:),  ALLOCATABLE :: ZFLD_LECT  ! Field read
REAL, DIMENSION(:),  POINTER     :: ZPNT_LECT
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('READ_GRIB_LARGE_RAIN',0,ZHOOK_HANDLE)
!
!
!*       1.0    preparing file openning
!               ----------------------
CALL GET_GRIB_CODE(HGRIB_TYPE,YREC_NAME,ITAB_NUM_GRIB,ITAB_TYPE_GRIB,ITAB_LEV1)
INUM_GRIB  = ITAB_NUM_GRIB(8)
ITYPE_GRIB = ITAB_TYPE_GRIB(8)
ILEV1      = ITAB_LEV1(8)
write(*,*) '********************************************'
write(*,*) 'Treatment of Large rain : INUM_GRIB=',INUM_GRIB 
write(*,*) '********************************************'

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
CALL READ_GRIB(HFILE,CINMODEL,IUNIT,INUM_GRIB,IRESP,PFIELD=ZPNT_LECT,&
               KLTYPE=ITYPE_GRIB,KLEV1=ILEV1)

IF (IRESP /=0) THEN
  WRITE(*,*) 'WARNING'
  WRITE(*,*) '-------'
  WRITE(*,*) 'error when reading article ', INUM_GRIB,'IRESP=',IRESP
  WRITE(*,*) 'default value may be used, who knows???'
  WRITE(*,*) ' '
ENDIF
ZFLD_LECT(:)  = ZPNT_LECT(:)
PLARGE_RAIN(:)= ZFLD_LECT(:)
!
DEALLOCATE(ZPNT_LECT)
DEALLOCATE(ZFLD_LECT)
!
IF (LHOOK) CALL DR_HOOK('READ_GRIB_LARGE_RAIN',1,ZHOOK_HANDLE)
!
END SUBROUTINE READ_GRIB_LARGE_RAIN
