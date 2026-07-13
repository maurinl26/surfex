!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!     #########
      SUBROUTINE READ_ISBA_NUDGING_n (IO, S, K, NK, NP, NPGLO, HPROGRAM, KI, ODIM)
!     ##################################
!
!!****  *READ_ISBA_NUDGING_n* - routine to initialise ISBA Snow and soil moisture nudging
!!                         
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
!!      B. Decharme   *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    12/2021
!!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_ISBA_OPTIONS_n, ONLY : ISBA_OPTIONS_t
USE MODD_ISBA_n,         ONLY : ISBA_S_t, ISBA_K_t, ISBA_NK_t, ISBA_NP_t, ISBA_P_t, ISBA_PE_t
!                                
USE MODD_SURF_PAR,       ONLY : XUNDEF, LEN_HREC
!
USE MODI_MAKE_CHOICE_ARRAY
USE MODI_PACK_SAME_RANK
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
TYPE(ISBA_OPTIONS_t), INTENT(INOUT) :: IO
TYPE(ISBA_S_t),       INTENT(INOUT) :: S
TYPE(ISBA_K_t),       INTENT(INOUT) :: K
TYPE(ISBA_NK_t),      INTENT(INOUT) :: NK
TYPE(ISBA_NP_t),      INTENT(INOUT) :: NP, NPGLO
!
CHARACTER(LEN=6),     INTENT(IN)    :: HPROGRAM ! calling program
INTEGER,              INTENT(IN)    :: KI       ! number of points
LOGICAL,              INTENT(IN)    :: ODIM
!
!*       0.2   Declarations of local variables
!              -------------------------------
!
TYPE(ISBA_K_t), POINTER  :: KK
TYPE(ISBA_P_t), POINTER  :: PK, PKGLO
!
INTEGER, PARAMETER ::   INMONTH=3
!
CHARACTER(LEN=LEN_HREC) :: YRECFM         ! Name of the article to be read
CHARACTER(LEN=LEN_HREC) :: YWORK          ! Work variable
!
REAL, DIMENSION(:,:)  ,ALLOCATABLE :: ZWORK2D  ! 2D array to write data in file
!
CHARACTER(LEN=4) :: YLVL
CHARACTER(LEN=2) :: YDAY, YMTH
!
INTEGER :: INDAYS   ! Number of days in the months
INTEGER :: IMONTH   ! Current month
INTEGER :: IYEAR    ! Current year
INTEGER :: INDAYTOT
!
INTEGER :: JI, JP, JL, JDAY, JMTH   ! loop counter
!
LOGICAL :: GYEARBI
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('READ_ISBA_NUDGING_N',0,ZHOOK_HANDLE)
!-------------------------------------------------------------------------------
!
!*       1.     Physical dimension
!               -----------------
!
!Compute the number of days in the month (every february is 29 days long)
!
IMONTH=S%TTIME%TDATE%MONTH
SELECT CASE (IMONTH)
  CASE(4,6,9,11)
    INDAYS=30
  CASE(1,3,5,7:8,10,12)
    INDAYS=31
  CASE(2)
    INDAYS=28
END SELECT
!
IYEAR=S%TTIME%TDATE%YEAR
GYEARBI=(((MOD(IYEAR,4)==0).AND.(MOD(IYEAR,100)/=0)).OR.(MOD(IYEAR,400)==0))
!
IF(GYEARBI.AND.IMONTH==2)THEN
  INDAYTOT=INDAYS+1
ELSE
  INDAYTOT=INDAYS
ENDIF
!
!-------------------------------------------------------------------------------
!
CALL ALLOC_READ_ISBA_NUDGING(NP)
IF(IO%LLULU)THEN
  CALL ALLOC_READ_ISBA_NUDGING(NPGLO)        
ENDIF
!
!-------------------------------------------------------------------------------
ALLOCATE(ZWORK2D(KI,IO%NPATCH))
!-------------------------------------------------------------------------------
!
!
!*       2.     Nudging mask:
!               -------------
!
IF(IO%LNUDG_SWE_MASK.OR.IO%LNUDG_WG_MASK) THEN
  DO JP = 1, IO%NPATCH
     KK => NK%AL(JP)  
     PK => NP%AL(JP)
     ALLOCATE(KK%XNUDG_MASK(PK%NSIZE_P))
     CALL PACK_SAME_RANK(PK%NR_P,K%XNUDG_MASK,KK%XNUDG_MASK)
  ENDDO
ENDIF
!
!-------------------------------------------------------------------------------
!
!
!*       3.     Snow mass nudging:
!               -----------------
!
IF(IO%LNUDG_SWE)THEN
  !
  !  The nudging is applied separately on each patch
  !
  DO JDAY=1,INDAYS
     WRITE(YDAY,'(I2)')JDAY
     YRECFM='N_SWE_DD'//ADJUSTL(YDAY(:LEN_TRIM(YDAY)))
     CALL MAKE_CHOICE_ARRAY(HPROGRAM, IO%NPATCH, ODIM, YRECFM, ZWORK2D)
     DO JP = 1,IO%NPATCH
        PK => NP%AL(JP)
        CALL PACK_SAME_RANK(PK%NR_P,ZWORK2D(:,JP),PK%XNUDG_SWE(:,JDAY))
        IF(IO%LLULU) NPGLO%AL(JP)%XNUDG_SWE(:,JDAY)=ZWORK2D(:,JP)
     ENDDO
  ENDDO
  !
  IF(INDAYTOT/=INDAYS)THEN
    DO JP = 1,IO%NPATCH
       PK => NP%AL(JP)
       PK%XNUDG_SWE(:,INDAYTOT)=PK%XNUDG_SWE(:,INDAYS) ! February 28 value for the 29th
       IF(IO%LLULU) NPGLO%AL(JP)%XNUDG_SWE(:,INDAYTOT)=NPGLO%AL(JP)%XNUDG_SWE(:,INDAYS)
    ENDDO
  ENDIF
  !
  !  Make sure undefined values are set to XUNDEF
  !
  DO JP = 1,IO%NPATCH
     PK => NP%AL(JP)
     DO JDAY=1,INDAYTOT
        DO JI=1,KI
           PK%XNUDG_SWE(JI,JDAY) = MAX(PK%XNUDG_SWE(JI,JDAY),0.0)
           IF(PK%XNUDG_SWE(JI,JDAY)>=1.E+15)THEN
             PK%XNUDG_SWE(JI,JDAY) = XUNDEF
           ENDIF
        ENDDO
     ENDDO
  ENDDO
  !
  IF(IO%LLULU)THEN
    DO JP = 1,IO%NPATCH
       PKGLO => NPGLO%AL(JP)
       DO JDAY=1,INDAYTOT
          DO JI=1,KI
             PKGLO%XNUDG_SWE(JI,JDAY) = MAX(PKGLO%XNUDG_SWE(JI,JDAY),0.0)
             IF(PKGLO%XNUDG_SWE(JI,JDAY)>=1.E+15)THEN
               PKGLO%XNUDG_SWE(JI,JDAY) = XUNDEF
             ENDIF
          ENDDO
       ENDDO
    ENDDO
  ENDIF
  !
ENDIF
!
!-------------------------------------------------------------------------------
!
!*       4.     Soil Moisture Nudging:
!               ---------------------
!
!The nudging is applied separately on each patch
!
IF(IO%CNUDG_WG=='DAY')THEN   
  !
  !  Nudging values are daily
  !
  DO JDAY=1,INDAYS
     WRITE(YDAY,'(I2)')JDAY
     YWORK='_DD'//ADJUSTL(YDAY(:LEN_TRIM(YDAY)))
     DO JL=1,IO%NGROUND_LAYER
        WRITE(YLVL,'(I4)') JL
        YRECFM='N_WG'//ADJUSTL(YLVL(:LEN_TRIM(YLVL)))
        YRECFM=ADJUSTR(YRECFM(:LEN_TRIM(YRECFM)))//ADJUSTL(YWORK(:LEN_TRIM(YWORK)))
        CALL MAKE_CHOICE_ARRAY(HPROGRAM, IO%NPATCH, ODIM, YRECFM, ZWORK2D)
        DO JP = 1,IO%NPATCH
           PK => NP%AL(JP)
           CALL PACK_SAME_RANK(PK%NR_P,ZWORK2D(:,JP),PK%XNUDG_WGTOT(:,JL,JDAY))
           IF(IO%LLULU) NPGLO%AL(JP)%XNUDG_WGTOT(:,JL,JDAY)=ZWORK2D(:,JP)
        ENDDO
     ENDDO
  ENDDO
  !
  IF(INDAYTOT/=INDAYS)THEN
    DO JP = 1,IO%NPATCH
       PK => NP%AL(JP)
       PK%XNUDG_WGTOT(:,:,INDAYS+1)=PK%XNUDG_WGTOT(:,:,INDAYS) ! February 28 value for the 29th
       IF(IO%LLULU) NPGLO%AL(JP)%XNUDG_WGTOT(:,:,INDAYTOT)=NPGLO%AL(JP)%XNUDG_WGTOT(:,:,INDAYS)
    ENDDO
  ENDIF
  !
ELSE
  !
  !  Nudging values are montly
  !
  DO JMTH=1,3
     WRITE(YMTH,'(I2)')JMTH-1
     YWORK='_MTH'//ADJUSTL(YMTH(:LEN_TRIM(YMTH)))
     DO JL=1,IO%NGROUND_LAYER
        WRITE(YLVL,'(I4)') JL
        YRECFM='N_WG'//ADJUSTL(YLVL(:LEN_TRIM(YLVL)))
        YRECFM=ADJUSTR(YRECFM(:LEN_TRIM(YRECFM)))//ADJUSTL(YWORK(:LEN_TRIM(YWORK)))
        CALL MAKE_CHOICE_ARRAY(HPROGRAM, IO%NPATCH, ODIM, YRECFM, ZWORK2D)
        DO JP = 1,IO%NPATCH
           PK => NP%AL(JP)
           CALL PACK_SAME_RANK(PK%NR_P,ZWORK2D(:,JP),PK%XNUDG_WGTOT(:,JL,JMTH))
           IF(IO%LLULU) NPGLO%AL(JP)%XNUDG_WGTOT(:,JL,JMTH)=ZWORK2D(:,JP)
        ENDDO
     ENDDO
  ENDDO
!
ENDIF
!
!Make sure undefined values are set to XUNDEF
!
IF (IO%CNUDG_WG/='DEF') THEN   
   DO JP = 1,IO%NPATCH
      PK => NP%AL(JP)
      WHERE((PK%XNUDG_WGTOT(:,:,:)>1.00).OR.(PK%XNUDG_WGTOT(:,:,:)<0.00))
             PK%XNUDG_WGTOT(:,:,:)=XUNDEF
      ENDWHERE
   ENDDO
ENDIF
!
!-------------------------------------------------------------------------------
DEALLOCATE(ZWORK2D)
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('READ_ISBA_NUDGING_N',1,ZHOOK_HANDLE)
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
CONTAINS
!-------------------------------------------------------------------------------
!
SUBROUTINE ALLOC_READ_ISBA_NUDGING(NA)
!
IMPLICIT NONE
!
TYPE(ISBA_NP_t),      INTENT(INOUT) :: NA
!
TYPE(ISBA_P_t),  POINTER :: PA
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
IF (LHOOK) CALL DR_HOOK('READ_ISBA_NUDGING_N:ALLOC_READ_ISBA_NUDGING',0,ZHOOK_HANDLE)
!
IF(IO%LNUDG_SWE)THEN
  DO JP = 1,IO%NPATCH
     PA => NA%AL(JP)
     ALLOCATE(PA%XNUDG_SWE(PA%NSIZE_P,INDAYTOT))
     PA%XNUDG_SWE(:,:) = XUNDEF
  ENDDO
ENDIF
!
!
IF(IO%CNUDG_WG=='DAY')THEN
  !
  !Nudging values are daily
  !
  DO JP = 1,IO%NPATCH
     PA => NA%AL(JP)
     ALLOCATE(PA%XNUDG_WGTOT(PA%NSIZE_P,IO%NGROUND_LAYER,INDAYS))
     PA%XNUDG_WGTOT(:,:,:)=XUNDEF
  ENDDO
ELSE
  !
  !Nudging values are monthly
  !
  DO JP = 1,IO%NPATCH
     PA => NA%AL(JP)
     ALLOCATE(PA%XNUDG_WGTOT(PA%NSIZE_P,IO%NGROUND_LAYER,INMONTH))
     PA%XNUDG_WGTOT(:,:,:)=XUNDEF
  ENDDO
  !
ENDIF
!
IF (LHOOK) CALL DR_HOOK('READ_ISBA_NUDGING_N:ALLOC_READ_ISBA_NUDGING',1,ZHOOK_HANDLE)
!
END SUBROUTINE ALLOC_READ_ISBA_NUDGING
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE READ_ISBA_NUDGING_n
