!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!     #######################
      SUBROUTINE READ_FORC_SAF(YSC,KDATE_DEB,PT2M,PQ2M,PRADSW,PRADLW,PFF10M,PRR1,KNB_HOUR_SAF)
!     #######################
!
!!****  *READ_FORC_SAFRAN*  
!!
!!    PURPOSE
!!    -------
!     This routine aims at reading forcing variables from grib files
!!**  METHOD
!!    ------
!!    Usage of table grid 2 :
!!    001 : P         - Pression                                  - Pa       - niveau 1/0
!!    011 : T         - Temperature                               - K        - niveau 105/2
!!    033 : U         - Premiere composante (zonale) du vent      - m s**-1  - niveau 105/10
!!    034 : V         - Sec. composante (meridienne) du vent      - m s**-1  - niveau 105/10
!!    051 : Q         - Humidite specifique                       - kg kg**-1- niveau 105/2
!!    061 : PRECIP    - Precip. totales (toutes formes)           - kg m**-2 - niveau 1/1
!!    111 : FLSOLAIRE - Bilan ray. courtes long. d'ondes (au sol) - W m**-2  - niveau 105/2
!!    112 : FLTHERM   - Bilan ray. grandes long. d'ondes (au sol) - W m**-2  - niveau 105/2
!!    (codes of parameters read thanks to routine GET_GRIB_CODE)
!!
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
USE MODD_SURF_PAR,    ONLY : XUNDEF
!
USE MODI_READ_INTERP_FORC_SAF
USE MODI_GET_GRIB_CODE
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
TYPE(SURFEX_t), INTENT(INOUT) :: YSC
!
INTEGER,                           INTENT(IN)  :: KDATE_DEB
REAL, DIMENSION(:,:),ALLOCATABLE,  INTENT(OUT) :: PT2M  ! forcing variables read on file
REAL, DIMENSION(:,:),ALLOCATABLE,  INTENT(OUT) :: PQ2M
REAL, DIMENSION(:,:),ALLOCATABLE,  INTENT(OUT) :: PRADSW
REAL, DIMENSION(:,:),ALLOCATABLE,  INTENT(OUT) :: PRADLW
REAL, DIMENSION(:,:),ALLOCATABLE,  INTENT(OUT) :: PFF10M
REAL, DIMENSION(:,:),ALLOCATABLE,  INTENT(OUT) :: PRR1
INTEGER,                           INTENT(IN)  :: KNB_HOUR_SAF
!
!
!*      0.2    declarations of local variables
!
!
INTEGER                   :: INUM_GRIB   ! number of parameter in grib file
INTEGER                   :: ITYPE_GRIB  ! number of level in grib file
INTEGER                   :: ILEV1       ! level in grib file
CHARACTER(LEN=28)         :: YFILE
CHARACTER(LEN=5), DIMENSION(10) :: YREC_NAME       ! Grib type of file to be read
INTEGER,          DIMENSION(10) :: ITAB_NUM_GRIB        ! Code of the parameter to get
INTEGER,          DIMENSION(10) :: ITAB_TYPE_GRIB
INTEGER,          DIMENSION(10) ::   ITAB_LEV1
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('READ_FORC_SAF',0,ZHOOK_HANDLE)
!
!
ALLOCATE (PT2M(KNB_HOUR_SAF,YSC%U%NDIM_FULL))
ALLOCATE (PQ2M(KNB_HOUR_SAF,YSC%U%NDIM_FULL))
ALLOCATE (PRADSW(KNB_HOUR_SAF,YSC%U%NDIM_FULL))
ALLOCATE (PRADLW(KNB_HOUR_SAF,YSC%U%NDIM_FULL))
ALLOCATE (PFF10M(KNB_HOUR_SAF,YSC%U%NDIM_FULL))
ALLOCATE (PRR1(KNB_HOUR_SAF,YSC%U%NDIM_FULL))
!
!*       1.0    default values
!               ----------------------
write(*,*) 'Dans read_forc_saf, date :',KDATE_DEB

PT2M(:,:)   = XUNDEF
PQ2M(:,:)   = XUNDEF
PRADSW(:,:) = XUNDEF
PRADLW(:,:) = XUNDEF
PFF10M(:,:) = XUNDEF
PRR1(:,:)   = XUNDEF
CALL GET_GRIB_CODE('SAFRAN',YREC_NAME,ITAB_NUM_GRIB,ITAB_TYPE_GRIB,ITAB_LEV1)
YFILE='datafile'

!*       2.0    reading values in grib file
!               ----------------------
!
!
! *** Treatment of 2m-temperature ***
INUM_GRIB  = ITAB_NUM_GRIB(1)
ITYPE_GRIB = ITAB_TYPE_GRIB(1)
ILEV1      = ITAB_LEV1(1)
write(*,*) '********************************************'
write(*,*) 'Treatment of 2m-temperature : INUM_GRIB=',INUM_GRIB
write(*,*) '********************************************'
CALL READ_INTERP_FORC_SAF(YSC,YFILE,PT2M,KDATE_DEB,INUM_GRIB,ITYPE_GRIB,ILEV1,KNB_HOUR_SAF)

! *** Treatment of 2m-specific humidity***
!
INUM_GRIB  = ITAB_NUM_GRIB(2)
ITYPE_GRIB = ITAB_TYPE_GRIB(2)
ILEV1      = ITAB_LEV1(2)
write(*,*) '********************************************'
write(*,*) 'Treatment of 2m-specific humidity : INUM_GRIB=',INUM_GRIB
write(*,*) '********************************************'
CALL READ_INTERP_FORC_SAF(YSC,YFILE,PQ2M,KDATE_DEB,INUM_GRIB,ITYPE_GRIB,ILEV1,KNB_HOUR_SAF)
!
! *** Treatment of SW radiation***
!
INUM_GRIB  = ITAB_NUM_GRIB(3)
ITYPE_GRIB = ITAB_TYPE_GRIB(3)
ILEV1      = ITAB_LEV1(3)
write(*,*) '********************************************'
write(*,*) 'Treatment of SW radiation : INUM_GRIB=',INUM_GRIB
write(*,*) '********************************************'
CALL READ_INTERP_FORC_SAF(YSC,YFILE,PRADSW,KDATE_DEB,INUM_GRIB,ITYPE_GRIB,ILEV1,KNB_HOUR_SAF)
!
! *** Treatment of LW radiation***
!
INUM_GRIB  = ITAB_NUM_GRIB(4)
ITYPE_GRIB = ITAB_TYPE_GRIB(4)
ILEV1      = ITAB_LEV1(4)
write(*,*) '********************************************'
write(*,*) 'Treatment of LW radiation : INUM_GRIB=',INUM_GRIB
write(*,*) '********************************************'
CALL READ_INTERP_FORC_SAF(YSC,YFILE,PRADLW,KDATE_DEB,INUM_GRIB,ITYPE_GRIB,ILEV1,KNB_HOUR_SAF)
!
! *** Treatment of wind velocity***
!
INUM_GRIB  = ITAB_NUM_GRIB(5)
ITYPE_GRIB = ITAB_TYPE_GRIB(5)
ILEV1      = ITAB_LEV1(5)
write(*,*) '********************************************'
write(*,*) 'Treatment of wind  velocity : INUM_GRIB=',INUM_GRIB
write(*,*) '********************************************'
CALL READ_INTERP_FORC_SAF(YSC,YFILE,PFF10M,KDATE_DEB,INUM_GRIB,ITYPE_GRIB,ILEV1,KNB_HOUR_SAF)
!
! *** Treatment of hourly rainfall***
!
INUM_GRIB  = ITAB_NUM_GRIB(8)
ITYPE_GRIB = ITAB_TYPE_GRIB(8)
ILEV1      = ITAB_LEV1(8)
write(*,*) '********************************************'
write(*,*) 'Treatment of hourly precipitation : INUM_GRIB=',INUM_GRIB
write(*,*) '********************************************'
CALL READ_INTERP_FORC_SAF(YSC,YFILE,PRR1,KDATE_DEB,INUM_GRIB,ITYPE_GRIB,ILEV1,KNB_HOUR_SAF)
!
WHERE (PT2M>=XUNDEF)
  PT2M=293.
ENDWHERE
WHERE (PQ2M>=XUNDEF)
  PQ2M=0.005
ENDWHERE
WHERE (PRADSW>=XUNDEF)
  PRADSW=10.
ENDWHERE
WHERE (PRADLW>=XUNDEF)
  PRADLW=0.
ENDWHERE
WHERE (PFF10M>=XUNDEF)
  PFF10M=1.
ENDWHERE
WHERE (PRR1>=XUNDEF)
  PRR1=0.
ENDWHERE
!
IF (LHOOK) CALL DR_HOOK('READ_FORC_SAF',1,ZHOOK_HANDLE)
!
END SUBROUTINE READ_FORC_SAF
