!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!     #########
      SUBROUTINE READ_ISBA_LANDUSE_n (IO, S, HPROGRAM, KI, KVERSION, ODIM   )
!     ##################################
!
!!****  *READ_ISBA_LANDUSE_n* - routine to initialise ISBA variables for Land-use Land Cover change
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
!!      Original    12/2023 Split from previous read_isban.F90 routine
!!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_ISBA_OPTIONS_n, ONLY : ISBA_OPTIONS_t
USE MODD_ISBA_n,         ONLY : ISBA_S_t
!
USE MODD_SURF_PAR,       ONLY : XUNDEF, NUNDEF, LEN_HREC
!
USE MODI_READ_SURF
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
!
CHARACTER(LEN=6),     INTENT(IN)    :: HPROGRAM ! calling program
INTEGER,              INTENT(IN)    :: KI       ! number of points
INTEGER,              INTENT(IN)    :: KVERSION
LOGICAL,              INTENT(IN)    :: ODIM
!
!*       0.2   Declarations of local variables
!              -------------------------------
!
CHARACTER(LEN=LEN_HREC) :: YRECFM         ! Name of the article to be read
CHARACTER(LEN=4)        :: YLVL
!
INTEGER :: IRESP             ! Error code after redding
INTEGER :: IWORK             ! Work integer
!
INTEGER :: JP, JL, JNL, JNS, JNB  ! loop counter on layers
INTEGER :: JI, JNCANT
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('READ_ISBA_LANDUSE_N',0,ZHOOK_HANDLE)
!-------------------------------------------------------------------------------
!
!*       1.     Physical dimension
!               -----------------
!
CALL ALLOC_READ_LULCC
!
!-------------------------------------------------------------------------------
!
!*       2. LULCC variables
!        ------------------
!
!
IF(KVERSION>=9)THEN
  !
  YRECFM='DTLUL'
  CALL READ_SURF(HPROGRAM,'DTLUL',S%TLULCC,IRESP)
  !
  ! * Water conservation
  !
  IF(IO%CISBA=='DIF')THEN
      YRECFM = 'WCONSRV'
      CALL READ_SURF(HPROGRAM,YRECFM,S%XWCONSRV(:),IRESP)
  ENDIF
  !
  ! * Biomass and carbon specific treatment
  !
  IF(IO%CPHOTO=='NCB')THEN
    !
    ! * Carbon conservation
    !   
    YRECFM = 'CCONSRV'
    CALL READ_SURF(HPROGRAM,YRECFM,S%XCCONSRV(:),IRESP)
    !
    ! * Land-use Land Cover change flux
    !
    YRECFM='FLUATM'
    CALL MAKE_CHOICE_ARRAY(HPROGRAM, IO%NPATCH, ODIM, YRECFM, S%XFLUATM)
    YRECFM='FLURES'
    CALL MAKE_CHOICE_ARRAY(HPROGRAM, IO%NPATCH, ODIM, YRECFM, S%XFLURES)
    !
  ENDIF
  !
ENDIF
!
!-------------------------------------------------------------------------------
!
!
!*       3. Carbon Managing variables
!        ----------------------------
!
IF(IO%LLULCC_MANAGE.AND.KVERSION>=9)THEN
  !
  YRECFM='FLUANT'
  CALL MAKE_CHOICE_ARRAY(HPROGRAM, IO%NPATCH, ODIM, YRECFM, S%XFLUANT)
  YRECFM='FANTATM'
  CALL MAKE_CHOICE_ARRAY(HPROGRAM, IO%NPATCH, ODIM, YRECFM, S%XFANTATM)    
  !
  DO JNCANT=1,IO%NNDECADAL
     WRITE(YLVL,'(I4)')  JNCANT
     YRECFM='CANTD'//ADJUSTL(YLVL(:LEN_TRIM(YLVL)))
     CALL MAKE_CHOICE_ARRAY(HPROGRAM, IO%NPATCH, ODIM, YRECFM, S%XCSTOCK_DECADAL(:,JNCANT,:))
     YRECFM='CEXPD'//ADJUSTL(YLVL(:LEN_TRIM(YLVL)))
     CALL MAKE_CHOICE_ARRAY(HPROGRAM, IO%NPATCH, ODIM, YRECFM, S%XEXPORT_DECADAL(:,JNCANT,:))
  ENDDO
  !
  DO JNCANT=1,IO%NNCENTURY
     WRITE(YLVL,'(I4)')  JNCANT
     YRECFM='CANTC'//ADJUSTL(YLVL(:LEN_TRIM(YLVL)))
     CALL MAKE_CHOICE_ARRAY(HPROGRAM, IO%NPATCH, ODIM, YRECFM, S%XCSTOCK_CENTURY(:,JNCANT,:))
     YRECFM='CEXPC'//ADJUSTL(YLVL(:LEN_TRIM(YLVL)))
     CALL MAKE_CHOICE_ARRAY(HPROGRAM, IO%NPATCH, ODIM, YRECFM, S%XEXPORT_CENTURY(:,JNCANT,:))
  ENDDO
  !
ENDIF
!
!-------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('READ_ISBA_LANDUSE_N',1,ZHOOK_HANDLE)
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
CONTAINS
!-------------------------------------------------------------------------------
!
SUBROUTINE ALLOC_READ_LULCC
!
IMPLICIT NONE
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
IF (LHOOK) CALL DR_HOOK('READ_ISBA_LANDUSE_N:ALLOC_READ_LULCC',0,ZHOOK_HANDLE)
!
S%TLULCC%TDATE%DAY=0
S%TLULCC%TDATE%MONTH=0
S%TLULCC%TDATE%YEAR=0
!
! * Water conservation
!
IF(IO%CISBA=='DIF')THEN
  ALLOCATE(S%XWCONSRV(KI))
  S%XWCONSRV(:) = 0.0
ENDIF
!
! * Biomass and carbon specific treatment
!
IF(IO%CPHOTO=='NCB')THEN
  !
  ! * Carbon conservation
  !   
  ALLOCATE(S%XCCONSRV(KI))
  S%XCCONSRV(:) = 0.0
  !
  ! * Land-use Land Cover change flux
  !
  ALLOCATE(S%XFLUATM (KI,IO%NPATCH))
  ALLOCATE(S%XFLURES (KI,IO%NPATCH))
  S%XFLUATM (:,:) = 0.0
  S%XFLURES (:,:) = 0.0
  !
ENDIF
!
! * Land-use Land Cover change carbon Managing variables 
!   (global because can not be packed)
!
IF(IO%LLULCC_MANAGE)THEN
  !
  ALLOCATE(S%XFLUANT (KI,IO%NPATCH))
  ALLOCATE(S%XFANTATM(KI,IO%NPATCH))
  S%XFLUANT (:,:) = 0.0
  S%XFANTATM(:,:) = 0.0
  !
  ALLOCATE(S%XEXPORT_DECADAL(KI,IO%NNDECADAL,IO%NPATCH))
  ALLOCATE(S%XCSTOCK_DECADAL(KI,IO%NNDECADAL,IO%NPATCH))
  ALLOCATE(S%XEXPORT_CENTURY(KI,IO%NNCENTURY,IO%NPATCH))
  ALLOCATE(S%XCSTOCK_CENTURY(KI,IO%NNCENTURY,IO%NPATCH))
  S%XEXPORT_DECADAL(:,:,:) = 0.0
  S%XCSTOCK_DECADAL(:,:,:) = 0.0
  S%XEXPORT_CENTURY(:,:,:) = 0.0
  S%XCSTOCK_CENTURY(:,:,:) = 0.0
  !
ENDIF
!
IF (LHOOK) CALL DR_HOOK('READ_ISBA_LANDUSE_N:ALLOC_READ_LULCC',1,ZHOOK_HANDLE)
!
END SUBROUTINE ALLOC_READ_LULCC
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE READ_ISBA_LANDUSE_n
