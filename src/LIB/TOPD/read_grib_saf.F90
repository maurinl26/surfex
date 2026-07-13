!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!     ##########################################################################
    SUBROUTINE READ_GRIB_SAF(HGRIB,KLUOUT,KPARAM,KDATE_DEB,KRET, &
                         KNI,PDOUT,KLTYPE,KLEV1,KLEV2,KNB_HOUR_SAF)
!     ##########################################################################
!!****  *READ_GRIB_SAF* - Read a grib field & performs interpolation (optional)
!!
!!    PURPOSE
!!    -------
!!
!!    Searchs & reads a field in a grib file. Returns the field or an interpolated
!!    one.
!!
!!    METHOD
!!    ------
!!
!!   The field to read is defined by :
!!    . its number (required),
!!    . its level type (may be set to -1 to accept any type),
!!    . its level id 1 (may be set to -1 to accept any value),
!!    . its level id 2 (may be set to -1 to accept any value).
!!   If '-1' values are used, the routine returns the founded values.
!!
!!
!!   EXTERNAL
!!   --------
!!   
!!   subroutine PBOPEN        : open a grib file
!!   subroutine PBGRIB        : read datas from a grib file
!!   subroutine PBCLOSE       : close a gribfile
!!   subroutine GRIBEX        : decodes grib data
!!
!!   IMPLICIT ARGUMENTS
!!   ------------------
!!
!!   REFERENCE
!!   ---------
!!
!!   This routine is based on books describing the Grib file format :
!!     'Encoding and Decoding Grib data', John D.Chambers(ECMWF), October 1995
!!     'Accessing GRIB and BUFR data', John D.Chambers(ECMWF), May 1994
!!     'A guide to Grib' Edition 1, John D.Stackpole(NOAA), March 1994
!!
!!   AUTHOR
!!   ------
!!
!!   V.Bousquet
!!
!!   MODIFICATIONS
!!   -------------
!!
!!   Original       07/01/1999
!!
!!   Modification 06/2003 (V. Masson) simplification for externalization
!!   Modification 11/2005 (I. Mallet) increase IPACK to read CEP files 
!!                                    from cy29r1 (24bits coding)
!!   Modification 03/2006 (I. Mallet) increase IPACK to read CEP files 
!!-----------------------------------------------------------------------------------
!
! 0. DECLARATIONS
! ---------------
!
USE MODD_SURF_PAR,  ONLY : XUNDEF
USE MODD_TYPE_DATE_SURF
USE MODD_GRID_GRIB, ONLY : NIDX
!
USE MODI_ADD_FORECAST_TO_DATE_SURF
USE MODE_READ_GRIB
!
USE YOMHOOK   ,ONLY : LHOOK, DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
! 0.1. Declaration of arguments
! -----------------------------
 CHARACTER(LEN=*),                   INTENT(IN)    :: HGRIB  ! Grib file name
 INTEGER,                            INTENT(IN)    :: KLUOUT ! logical unit of output listing
 INTEGER,                            INTENT(IN)    :: KPARAM ! Parameter to read
 INTEGER,                            INTENT(IN)    :: KDATE_DEB
 INTEGER,                            INTENT(OUT)   :: KRET   ! Result
 INTEGER, OPTIONAL,                  INTENT(IN)    :: KNI    ! Number of input points
 REAL,    OPTIONAL, DIMENSION(:,:),    POINTER     :: PDOUT  ! Output datas
 INTEGER, OPTIONAL,                  INTENT(INOUT) :: KLTYPE ! type of level (Grib code table 3)
 INTEGER, OPTIONAL,                  INTENT(INOUT) :: KLEV1  ! level definition
 INTEGER, OPTIONAL,                  INTENT(INOUT) :: KLEV2  ! level definition
 INTEGER,                            INTENT(IN)    :: KNB_HOUR_SAF ! number of steps to get
! 
! 0.2. Declaration of local variables
! -----------------------------------
! Variable involved in the task of reading the grib file
INTEGER                            :: IPACK     ! Size of the different buffer arrays
PARAMETER (IPACK=2800000)                       ! used to store the Grib informations
REAL,    DIMENSION(:), POINTER     :: ZPT_FIELD
INTEGER                            :: IRET      ! Return code from subroutines
INTEGER                            :: IUNIT     ! Unit number attached to the file
INTEGER                            :: IPARAM
INTEGER                            :: ICNT ! date in a interger format
INTEGER(KIND=kindOfInt)            :: IGRIB
!
! Local variables
INTEGER                            :: IFOUND    ! Number of correct parameters
INTEGER                            :: IYEAR_DEB,IMONTH_DEB,IDAY_DEB,IHOUR_DEB
INTEGER, DIMENSION(:), ALLOCATABLE :: IYEAR,IMONTH,IDAY,IHOUR
INTEGER                            :: IY,IM,ID
REAL                               :: ZTIME
INTEGER                            :: JWRK2,JSTP
INTEGER, DIMENSION(:), ALLOCATABLE :: IREAD
INTEGER                            :: ILTYPE, ILEV1, ILEV2
INTEGER                            :: IUNITTIME,IP1
INTEGER                            :: ITIME,ISIZE
TYPE (DATE_TIME)                   :: TPTIME_GRIB    ! current date and time
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-----------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('READ_GRIB_SAF',0,ZHOOK_HANDLE)
!
!* 0. Stores dates
!----------------
!
IYEAR_DEB=FLOOR(KDATE_DEB/1000000.)
IMONTH_DEB=FLOOR(KDATE_DEB/10000.)-IYEAR_DEB*100
IDAY_DEB=FLOOR(KDATE_DEB/100.)-IYEAR_DEB*10000-IMONTH_DEB*100
IHOUR_DEB=KDATE_DEB-IYEAR_DEB*1000000-IMONTH_DEB*10000-IDAY_DEB*100
!
IY=IYEAR_DEB
IM=IMONTH_DEB
ID=IDAY_DEB
ZTIME=IHOUR_DEB*3600.
!
ALLOCATE(IYEAR(KNB_HOUR_SAF))
ALLOCATE(IMONTH(KNB_HOUR_SAF))
ALLOCATE(IDAY(KNB_HOUR_SAF))
ALLOCATE(IHOUR(KNB_HOUR_SAF))

DO JSTP=1,KNB_HOUR_SAF
  IYEAR(JSTP)=IY
  IMONTH(JSTP)=IM
  IDAY(JSTP) =ID
  IHOUR(JSTP)=ZTIME/3600.
  ZTIME= ZTIME+3600.
  IF (86400.-ZTIME < 1.E-6) THEN
    CALL ADD_FORECAST_TO_DATE_SURF(IY,IM,ID,ZTIME)
  ENDIF
ENDDO!JSTP
!
ALLOCATE (IREAD(KNB_HOUR_SAF))
IREAD(:)=0
!
IF (PRESENT(KNI) .AND. PRESENT(PDOUT)) THEN
  ALLOCATE (PDOUT(KNB_HOUR_SAF,KNI))
  PDOUT(:,:)=XUNDEF
ENDIF
!
!* 1. Search for the wanted field
!--------------------------------
!
write(*,*)  'READ_GRIB_SAF PARAM: ',KPARAM
ICNT=0
!
! 1.1 Open the file
!
CALL GRIB_OPEN_FILE(IUNIT,HGRIB,'R',IRET)
IF (IRET /= 0) THEN
  CALL ABOR1_SFX('READ_GRIB_SAF: Error opening the grib file '//HGRIB)
END IF
!
CALL GRIB_MULTI_SUPPORT_ON()
CALL GRIB_NEW_FROM_FILE(IUNIT,IGRIB,IRET)
IF (IRET /= 0) THEN
  CALL ABOR1_SFX('READ_GRIB_SAF: Error in reading the grib file')
END IF

DO WHILE (IRET /= GRIB_END_OF_FILE ) 
  ILTYPE=-2
  IF (PRESENT(KLTYPE)) ILTYPE=KLTYPE
  ILEV1=-2
  IF (PRESENT(KLEV1)) ILEV1=KLEV1
  ILEV2=-2
  IF (PRESENT(KLEV2)) ILEV2=KLEV2
  IFOUND = 0
  KRET=0
  !
  IPARAM=KPARAM
  CALL GET_GRIB_MESSAGE_FROMFILE(IUNIT,KLUOUT,IPARAM,ILTYPE,ILEV1,ILEV2,IGRIB,IFOUND)
  IF (IFOUND==4) THEN
    IF (.NOT.ASSOCIATED(ZPT_FIELD)) THEN
      CALL GRIB_GET_SIZE(IGRIB,'values',ISIZE,KRET)
      IF (KRET.NE.0) CALL ABOR1_SFX(" READ_GRIB: Problem getting size of values")
      ALLOCATE(ZPT_FIELD(ISIZE))
    ENDIF
    WRITE (KLUOUT,'(A)') ' | Reading date'
    !
    CALL GRIB_GET(IGRIB,'year',TPTIME_GRIB%TDATE%YEAR,IRET)
    CALL GRIB_GET(IGRIB,'month',TPTIME_GRIB%TDATE%MONTH,IRET)
    CALL GRIB_GET(IGRIB,'day',TPTIME_GRIB%TDATE%DAY,IRET)
    CALL GRIB_GET(IGRIB,'time',ITIME,IRET)
    TPTIME_GRIB%TIME=INT(ITIME/100)*3600+(ITIME-INT(ITIME/100)*100)*60
    !
    !  
    CALL GRIB_GET(IGRIB,'P1',IP1,IRET)
    IF ( IP1>0 ) THEN
      CALL GRIB_GET(IGRIB,'unitOfTimeRange',IUNITTIME,IRET)      
      SELECT CASE (IUNITTIME)       ! Time unit indicator
        CASE (1)                    !hour
          TPTIME_GRIB%TIME   = TPTIME_GRIB%TIME + IP1*3600.
        CASE (0)                    !minute
          TPTIME_GRIB%TIME   = TPTIME_GRIB%TIME + IP1*60.
      END SELECT
    ENDIF
    !
    write(*,*) 'READ_GRIB_SAF DATE : ',TPTIME_GRIB%TDATE%YEAR,TPTIME_GRIB%TDATE%MONTH,TPTIME_GRIB%TDATE%DAY,ITIME/100
    DO JSTP=1,KNB_HOUR_SAF
      IF ((TPTIME_GRIB%TDATE%YEAR ==IYEAR(JSTP) ) .AND.&
        (TPTIME_GRIB%TDATE%MONTH==IMONTH(JSTP)) .AND.&
        (TPTIME_GRIB%TDATE%DAY  ==IDAY(JSTP)  ) .AND.&
        (ITIME/100==IHOUR(JSTP) ) ) THEN ! same date
      
        CALL GRIB_GET(IGRIB,'values',ZPT_FIELD,KRET)

        IF (KRET.NE.0) CALL ABOR1_SFX(" READ_GRIB: Problem getting values")

        write(*,*) 'min Max zpt_field :',MINVAL(ZPT_FIELD),MAXVAL(ZPT_FIELD)

        PDOUT(JSTP,1:KNI)=ZPT_FIELD(1:KNI)

        write(*,*) 'Field filled for time :',JSTP

        IREAD(JSTP)=1
      ENDIF
      IF (COUNT(IREAD==1)==KNB_HOUR_SAF) GOTO 60
    ENDDO
  ENDIF
  IRET=KRET
  CALL GRIB_NEW_FROM_FILE(IUNIT,IGRIB,IRET) 
END DO
60 CONTINUE

!-----------------------------------------------------------------------------------
!
!* 2. Closes file
!----------------
!
CALL GRIB_RELEASE(IGRIB,KRET)
IF (KRET.NE.0) CALL ABOR1_SFX("READ_GRIB_SAF: Problem releasing memory")
CALL GRIB_CLOSE_FILE(IUNIT)
!
DO JSTP=1,KNB_HOUR_SAF
  DO JWRK2=1,KNI
    IF (PDOUT(JSTP,JWRK2)==XUNDEF) THEN
      PDOUT(JSTP,JWRK2)=PDOUT(JSTP-1,JWRK2)
    ENDIF
  ENDDO
ENDDO
!
DEALLOCATE(IYEAR)
DEALLOCATE(IMONTH)
DEALLOCATE(IDAY)
DEALLOCATE(IHOUR)
DEALLOCATE (IREAD)
!-----------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('READ_GRIB_SAF',1,ZHOOK_HANDLE)
!
CONTAINS
!     ####################
      SUBROUTINE GET_GRIB_MESSAGE_FROMFILE(KUNIT,KLUOUT,KPARAM,KLTYPE,KLEV1,KLEV2,KGRIB,KFOUND)
!     ####################
!
USE MODD_GRID_GRIB, ONLY : NIDX
!
IMPLICIT NONE
!
INTEGER, INTENT(IN) :: KUNIT!unit of input grib file
INTEGER, INTENT(IN) :: KLUOUT!unit of output  file
INTEGER, INTENT(INOUT)    :: KPARAM ! Parameter to read
INTEGER, INTENT(INOUT)  :: KLTYPE
INTEGER, INTENT(INOUT)  :: KLEV1
INTEGER, INTENT(INOUT)  :: KLEV2
INTEGER(KIND=kindOfInt), INTENT(INOUT) :: KGRIB
INTEGER, INTENT(OUT) :: KFOUND
!
INTEGER :: ILTYPE,IPARAM
INTEGER :: ILEV1
INTEGER :: ILEV2
INTEGER(KIND=kindOfInt) :: IRET
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
IF (LHOOK) CALL DR_HOOK(' GET_GRIB_MESSAGE',0,ZHOOK_HANDLE)
!
  IRET = 0
  KFOUND=0
  !
  CALL GRIB_GET(KGRIB,'indicatorOfParameter',IPARAM,IRET)
  CALL TEST_IRET(KLUOUT,IPARAM,KPARAM,IRET)
  IF (IRET.EQ.0)  KFOUND = KFOUND + 1
  IF (KLTYPE/=-2) THEN
    CALL GRIB_GET(KGRIB,'indicatorOfTypeOfLevel',ILTYPE,IRET)
    CALL TEST_IRET(KLUOUT,ILTYPE,KLTYPE,IRET)
  ENDIF
  !
  IF (IRET.EQ.0) THEN
    !
    KFOUND = KFOUND + 1
    !
    IF (KLEV1/=-2) THEN
      CALL GRIB_GET(KGRIB,'topLevel',ILEV1,IRET)
      CALL TEST_IRET(KLUOUT,ILEV1,KLEV1,IRET)
    ENDIF
    !
    IF (IRET.EQ.0) THEN
      !
      KFOUND = KFOUND + 1
      !
      IF (KLEV2/=-2) THEN
        CALL GRIB_GET(KGRIB,'bottomLevel',ILEV2,IRET)
        CALL TEST_IRET(KLUOUT,ILEV2,KLEV2,IRET)
      ENDIF
      !
      IF (IRET.EQ.0) KFOUND = KFOUND + 1
      !
    ENDIF
    !
  ENDIF
  !
!
IF (LHOOK) CALL DR_HOOK(' GET_GRIB_MESSAGE',1,ZHOOK_HANDLE)
!
END SUBROUTINE GET_GRIB_MESSAGE_FROMFILE
!
!-----------------------------------------------------------------------------------
!       ##############
        SUBROUTINE TEST_IRET(KLUOUT,VAL1,VAL0,KRET)
!       ##############
!
IMPLICIT NONE
!
INTEGER, INTENT(IN) :: KLUOUT ! logical unit of output listing
INTEGER, INTENT(IN) :: VAL1
INTEGER, INTENT(INOUT) :: VAL0
INTEGER(KIND=kindOfInt), INTENT(INOUT) :: KRET   ! number of the message researched
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
IF (LHOOK) CALL DR_HOOK(' TEST_IRET',0,ZHOOK_HANDLE)
!
IF (KRET > 0) THEN
  WRITE (KLUOUT,'(A)')' | Error encountered in the Grib file, skipping field'
ELSE IF (KRET == -6) THEN
  WRITE (KLUOUT,'(A)')' | ECMWF pseudo-Grib data encountered, skipping field'
ELSEIF (VAL1 /= VAL0) THEN
  IF (VAL0 == -1) THEN
    VAL0 = VAL1
  ELSE
    KRET=1
  ENDIF
ENDIF
!
IF (LHOOK) CALL DR_HOOK(' TEST_IRET',1,ZHOOK_HANDLE)
!
END SUBROUTINE TEST_IRET
!
!-----------------------------------------------------------
!
END SUBROUTINE READ_GRIB_SAF
