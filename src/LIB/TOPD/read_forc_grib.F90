!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!     #######################
      SUBROUTINE READ_FORC_GRIB(YSC,HFILE,CGRIB_TYPE,PT2M,PQ2M,PRADSW, &
                      PRADLW,PU10M,PV10M,PPS,PRR1,PRS1,OINTERP)
!     #######################
!
!!****  *READ_FORC_GRIB*  
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
USE MODD_SURF_PAR, ONLY   :XUNDEF
!
USE MODI_READ_INTERP_FORC_GRIB
USE MODI_GET_GRIB_CODE
!
USE MODE_THERMOS
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
CHARACTER(LEN=*),   INTENT(IN)  :: HFILE      ! Grib file to be read
CHARACTER(LEN=*),   INTENT(IN)  :: CGRIB_TYPE ! Grib type of file to be read
REAL, DIMENSION(:), INTENT(OUT) :: PT2M       ! forcing variables read on file
REAL, DIMENSION(:), INTENT(OUT) :: PQ2M
REAL, DIMENSION(:), INTENT(OUT) :: PRADSW
REAL, DIMENSION(:), INTENT(OUT) :: PRADLW
REAL, DIMENSION(:), INTENT(OUT) :: PU10M
REAL, DIMENSION(:), INTENT(OUT) :: PV10M
REAL, DIMENSION(:), INTENT(OUT) :: PPS
REAL, DIMENSION(:), INTENT(OUT) :: PRR1
REAL, DIMENSION(:), INTENT(OUT) :: PRS1

LOGICAL, INTENT(IN) :: OINTERP
!
!
!*      0.2    declarations of local variables
!
TYPE(SURFEX_t), INTENT(INOUT) :: YSC
!
INTEGER                   :: INUM_GRIB   ! number of parameter in grib file
INTEGER                   :: ITYPE_GRIB  ! number of level in grib file
INTEGER                   :: ILEV1       ! level in grib file
CHARACTER(LEN=5), DIMENSION(10) :: YREC_NAME       ! Grib type of file to be read
INTEGER,          DIMENSION(10) :: ITAB_NUM_GRIB        ! Code of the parameter to get
INTEGER,          DIMENSION(10) :: ITAB_TYPE_GRIB
INTEGER,          DIMENSION(10) :: ITAB_LEV1
REAL,      DIMENSION(SIZE(PRR1)):: ZH2M,ZQSATS,ZRG1,ZRS1
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('READ_FORC_GRIB',0,ZHOOK_HANDLE)
!
!*       1.0    default values
!               ----------------------
PT2M(:)   = XUNDEF
PQ2M(:)   = XUNDEF
PRADSW(:) = XUNDEF
PRADLW(:) = XUNDEF
PU10M(:)  = XUNDEF
PV10M(:)  = XUNDEF
PPS(:)    = XUNDEF
PRR1(:)   = XUNDEF
PRS1(:)   = XUNDEF
CALL GET_GRIB_CODE(CGRIB_TYPE,YREC_NAME,ITAB_NUM_GRIB,ITAB_TYPE_GRIB,ITAB_LEV1)
!
!*       2.0    reading values in grib file
!               ----------------------
! *** Treatment of 2m-temperature ***
INUM_GRIB  = ITAB_NUM_GRIB(1)
ITYPE_GRIB = ITAB_TYPE_GRIB(1)
ILEV1      = ITAB_LEV1(1)
write(*,*) '********************************************'
write(*,*) 'Treatment of 2m-temperature : INUM_GRIB=',INUM_GRIB
write(*,*) '********************************************'
CALL READ_INTERP_FORC_GRIB(YSC,HFILE,PT2M,INUM_GRIB,ITYPE_GRIB,ILEV1,OINTERP)

! *** Treatment of 2m-specific humidity***
!
INUM_GRIB  = ITAB_NUM_GRIB(2)
ITYPE_GRIB = ITAB_TYPE_GRIB(2)
ILEV1      = ITAB_LEV1(2)
write(*,*) '********************************************'
write(*,*) 'Treatment of 2m-specific humidity :INUM_GRIB=',INUM_GRIB
write(*,*) '********************************************'
CALL READ_INTERP_FORC_GRIB(YSC,HFILE,PQ2M,INUM_GRIB,ITYPE_GRIB,ILEV1,OINTERP)
!
! *** Treatment of SW radiation***
!
INUM_GRIB  = ITAB_NUM_GRIB(3)
ITYPE_GRIB = ITAB_TYPE_GRIB(3)
ILEV1      = ITAB_LEV1(3)
write(*,*) '********************************************'
write(*,*) 'Treatment of SW radiation : INUM_GRIB=',INUM_GRIB
write(*,*) '********************************************'
CALL READ_INTERP_FORC_GRIB(YSC,HFILE,PRADSW,INUM_GRIB,ITYPE_GRIB,ILEV1,OINTERP)
!
! *** Treatment of LW radiation***
!
INUM_GRIB  = ITAB_NUM_GRIB(4)
ITYPE_GRIB = ITAB_TYPE_GRIB(4)
ILEV1      = ITAB_LEV1(4)
write(*,*) '********************************************'
write(*,*) 'Treatment of LW radiation : INUM_GRIB=',INUM_GRIB
write(*,*) '********************************************'
CALL READ_INTERP_FORC_GRIB(YSC,HFILE,PRADLW,INUM_GRIB,ITYPE_GRIB,ILEV1,OINTERP)
!
! *** Treatment of u-wind component***
!
INUM_GRIB  = ITAB_NUM_GRIB(5)
ITYPE_GRIB = ITAB_TYPE_GRIB(5)
ILEV1      = ITAB_LEV1(5)
write(*,*) '********************************************'
write(*,*) 'Treatment of u-wind component : INUM_GRIB=',INUM_GRIB
write(*,*) '********************************************'
CALL READ_INTERP_FORC_GRIB(YSC,HFILE,PU10M,INUM_GRIB,ITYPE_GRIB,ILEV1,OINTERP)
!
! *** Treatment of v-wind component***
!
INUM_GRIB  = ITAB_NUM_GRIB(6)
ITYPE_GRIB = ITAB_TYPE_GRIB(6)
ILEV1      = ITAB_LEV1(6)
write(*,*) '********************************************'
write(*,*) 'Treatment of v-wind component : INUM_GRIB=',INUM_GRIB
write(*,*) '********************************************'
CALL READ_INTERP_FORC_GRIB(YSC,HFILE,PV10M,INUM_GRIB,ITYPE_GRIB,ILEV1,OINTERP)
!
! *** Treatment of hourly rainfall***
!
INUM_GRIB  = ITAB_NUM_GRIB(8)
ITYPE_GRIB = ITAB_TYPE_GRIB(8)
ILEV1      = ITAB_LEV1(8)
write(*,*) '********************************************'
write(*,*) 'Treatment of hourly precipitation : INUM_GRIB=',INUM_GRIB
write(*,*) '********************************************'
CALL READ_INTERP_FORC_GRIB(YSC,HFILE,PRR1,INUM_GRIB,ITYPE_GRIB,ILEV1,OINTERP)
!
!*** Treatment of pressure***
!
INUM_GRIB  = ITAB_NUM_GRIB(7)
ITYPE_GRIB = ITAB_TYPE_GRIB(7)
ILEV1      = ITAB_LEV1(7)
write(*,*) '********************************************'
write(*,*) 'Treatment of surface pressure : INUM_GRIB=',INUM_GRIB
write(*,*) '********************************************'
CALL READ_INTERP_FORC_GRIB(YSC,HFILE,PPS,INUM_GRIB,ITYPE_GRIB,ILEV1,OINTERP)
!
! *** Treatment of hourly solid precipitation (snow)***
!
INUM_GRIB  = ITAB_NUM_GRIB(9)
ITYPE_GRIB = ITAB_TYPE_GRIB(9)
ILEV1      = ITAB_LEV1(9)
write(*,*) '********************************************'
write(*,*) 'Treatment of hourly solid precip. (snow): INUM_GRIB=',INUM_GRIB
write(*,*) '********************************************'
CALL READ_INTERP_FORC_GRIB(YSC,HFILE,ZRS1,INUM_GRIB,ITYPE_GRIB,ILEV1,OINTERP)
!
! *** Treatment of hourly solid precipitation (snow)***
!
INUM_GRIB  = ITAB_NUM_GRIB(10)
ITYPE_GRIB = ITAB_TYPE_GRIB(10)
ILEV1      = ITAB_LEV1(10)
write(*,*) '********************************************'
write(*,*) 'Treatment of hourly solid precip. (graupel): INUM_GRIB=',INUM_GRIB
write(*,*) '********************************************'
CALL READ_INTERP_FORC_GRIB(YSC,HFILE,ZRG1,INUM_GRIB,ITYPE_GRIB,ILEV1,OINTERP)
!
! *** Treatment of hourly solid precipitation (total)***
!
PRS1(:) = ZRG1(:)+ZRS1(:)
!
! *** Treatment of 2m-relative humidity ***
! *** if the specific humidity is missing***
IF (MINVAL(PQ2M)==XUNDEF) THEN
!
  IF (CGRIB_TYPE/='PEAROP')THEN
    INUM_GRIB=52
    ITYPE_GRIB=ITAB_TYPE_GRIB(2)
    ILEV1=ITAB_LEV1(2)
  ELSE
    INUM_GRIB=254
    ITYPE_GRIB=1
    ILEV1=0
  ENDIF
  write(*,*) '********************************************'
  write(*,*) 'Treatment of 2m-relative humidity :INUM_GRIB=',INUM_GRIB
  write(*,*) '********************************************'
  CALL READ_INTERP_FORC_GRIB(YSC,HFILE,ZH2M,INUM_GRIB,ITYPE_GRIB,ILEV1,OINTERP)
  ZQSATS(:) = QSAT(PT2M(:),PPS(:))
  PQ2M(:)   = ZH2M(:) *ZQSATS(:)/100.
!
ENDIF
!
IF (LHOOK) CALL DR_HOOK('READ_FORC_GRIB',1,ZHOOK_HANDLE)
!
END SUBROUTINE READ_FORC_GRIB
