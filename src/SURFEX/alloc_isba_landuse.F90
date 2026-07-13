!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!###############################################################################
SUBROUTINE ALLOC_ISBA_LANDUSE (U, IO, DTI, S, SOLD, NPGLO, NPEGLO, HPROGRAM, KI)  
!###############################################################################
!
!!****  *ALLOC_ISBA_LANDUSE* - routine to allocate required variables for land use for ISBA field
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
!!      Original    12/2023
!!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_SURF_ATM_n,     ONLY : SURF_ATM_t
USE MODD_ISBA_OPTIONS_n, ONLY : ISBA_OPTIONS_t
USE MODD_ISBA_n,         ONLY : ISBA_S_t, ISBA_P_t, ISBA_PE_t, ISBA_NP_t, ISBA_NPE_t, &
                                ISBA_S_INIT, ISBA_NP_INIT, ISBA_NPE_INIT
USE MODD_DATA_ISBA_n,    ONLY : DATA_ISBA_t
!
USE MODD_SURF_PAR,       ONLY : XUNDEF, XSURF_EPSILON, LEN_HREC
!
USE MODD_DATA_COVER_PAR, ONLY : NVEGTYPE
!
USE MODD_ASSIM,          ONLY : LASSIM, CASSIM_ISBA
!
USE MODI_ABOR1_SFX
!
USE MODI_SURF_PATCH
USE MODI_READ_SURF
!
USE YOMHOOK,  ONLY : LHOOK,   DR_HOOK
USE PARKIND1, ONLY : JPRB
!
IMPLICIT NONE
!
!*       0.1   Declarations of arguments
!              -------------------------
!
TYPE(SURF_ATM_t),      INTENT(INOUT) :: U
TYPE(ISBA_OPTIONS_t),  INTENT(INOUT) :: IO
TYPE(DATA_ISBA_t),     INTENT(INOUT) :: DTI
TYPE(ISBA_S_t),        INTENT(INOUT) :: S, SOLD
TYPE(ISBA_NP_t),       INTENT(INOUT) :: NPGLO
TYPE(ISBA_NPE_t),      INTENT(INOUT) :: NPEGLO
!
CHARACTER(LEN=6),      INTENT(IN)    :: HPROGRAM          ! program calling surf. schemes
INTEGER,               INTENT(IN)    :: KI
!
!
!*       0.2   Declarations of local arguments on Patch grid
!
TYPE(ISBA_P_t),  POINTER :: PKGLO
TYPE(ISBA_PE_t), POINTER :: PEKGLO
!
!*       0.3   Declarations of local arguments on complete ISBA grid
!
REAL, DIMENSION(KI,NVEGTYPE) :: ZVEGTYPE ! work array (ISBA grid)
!
CHARACTER(LEN=LEN_HREC)      :: YRECFM         ! Name of the article to be read
CHARACTER(LEN=4)             :: YLVL
!
REAL    :: ZEPSILON
INTEGER :: IRESP        ! Error code after redding
INTEGER :: JI, JP, JVEG ! loop counter
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('ALLOC_ISBA_LANDUSE',0,ZHOOK_HANDLE)
!-------------------------------------------------------------------------------
!
!*       1.     Control checks
!               --------------
!
IF(TRIM(CASSIM_ISBA)/='OI'.OR.LASSIM) THEN
  CALL ABOR1_SFX('ALLOC_ISBA_LANDUSE: ABORT because Assimilation procedure not implemented under Land-use change case')
ENDIF
!
IF(U%LECOSG) THEN
  CALL ABOR1_SFX('ALLOC_ISBA_LANDUSE: ABORT because ECOSG not yet implemented under Land-use change case')
ENDIF
!
IF(HPROGRAM=='ASCII'.OR.HPROGRAM=='TEXTE')THEN
  ZEPSILON=1.0E-8
ELSE
  ZEPSILON=XSURF_EPSILON
ENDIF
!
!-------------------------------------------------------------------------------
!
!*       2.     Old vegtype and patch grids 
!               ---------------------------
!  
!* Work pointer init
!
CALL ISBA_S_INIT(SOLD)  
!  
!* Fraction of each vegetation type for each grid mesh previous year (-)
!
ZVEGTYPE(:,:)=XUNDEF
!
DO JVEG=1,NVEGTYPE
   WRITE(YLVL,'(I4)') JVEG
   YRECFM='VEGTYPE'//ADJUSTL(YLVL(:LEN_TRIM(YLVL)))
   CALL READ_SURF(HPROGRAM,YRECFM,ZVEGTYPE(:,JVEG),IRESP)
ENDDO
!
ALLOCATE(SOLD%XPATCH        (KI,         IO%NPATCH))
ALLOCATE(SOLD%XVEGTYPE_PATCH(KI,NVEGTYPE,IO%NPATCH))
!
CALL SURF_PATCH(IO%NPATCH,ZVEGTYPE,DTI%NPAR_VEG_IRR_USE,SOLD%XPATCH,SOLD%XVEGTYPE_PATCH)
!
!-------------------------------------------------------------------------------
!
!*       3.     Work pointers to read ISBA variable on the global grid 
!               ------------------------------------------------------
!
!* Find if PATCH distribution has changed due to LULCC
!
IO%LLULU=ANY(ABS(S%XPATCH(:,:)-SOLD%XPATCH(:,:))>ZEPSILON)
!
IF(IO%LLULU)THEN
  !
  CALL CHECK_DATE
  !
  CALL ISBA_NP_INIT(NPGLO,IO%NPATCH)
  CALL ISBA_NPE_INIT(NPEGLO,IO%NPATCH)
  !
  DO JP = 1,IO%NPATCH
     !
     PKGLO => NPGLO%AL(JP)
     !
     !dimension of the patch set to global dimenssion
     PKGLO%NSIZE_P = KI   
     !
     !mask of the patch in tile nature
     ALLOCATE(PKGLO%NR_P(KI))
     DO JI = 1,KI
        PKGLO%NR_P(JI)=JI
     ENDDO
     !    
     !patch in tile nature
     ALLOCATE(PKGLO%XPATCH(KI))
     PKGLO%XPATCH(:)=SOLD%XPATCH(:,JP)
     !
  ENDDO  
  !
ENDIF
!
!-------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('ALLOC_ISBA_LANDUSE',1,ZHOOK_HANDLE)
!-------------------------------------------------------------------------------
!
CONTAINS
!
!-------------------------------------------------------------------------------
!
SUBROUTINE CHECK_DATE
!
USE MODI_GET_LUOUT
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
INTEGER :: ILUOUT
!
IF (LHOOK) CALL DR_HOOK('ALLOC_ISBA_LANDUSE:CHECK_DATE',0,ZHOOK_HANDLE)
!
CALL GET_LUOUT(HPROGRAM,ILUOUT)
!
S%TLULCC%TDATE%DAY=S%TTIME%TDATE%DAY
S%TLULCC%TDATE%MONTH=S%TTIME%TDATE%MONTH
!
IF(S%TLULCC%TDATE%YEAR==0)THEN
  S%TLULCC%TDATE%YEAR=S%TTIME%TDATE%YEAR
ELSE
  S%TLULCC%TDATE%YEAR=S%TLULCC%TDATE%YEAR+1
ENDIF
!
IF(S%TLULCC%TDATE%YEAR/=S%TTIME%TDATE%YEAR)THEN
  WRITE(ILUOUT,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'  
  WRITE(ILUOUT,*)'!!  LAND USE CHECK_DATE: ABORT LAND USE DATE INCONSISTENCY' 
  WRITE(ILUOUT,*)'!!  LAND USE SHOULD BE ONCE BY YEAR' 
  WRITE(ILUOUT,*)'!!  LAND USE YEAR SHOULD BE :',S%TLULCC%TDATE%YEAR 
  WRITE(ILUOUT,*)'!!  ISBA CURRENT YEAR IS    :',S%TTIME%TDATE%YEAR
  WRITE(ILUOUT,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'  
  CALL ABOR1_SFX('ALLOC_ISBA_LANDUSE:CHECK_DATE: ABORT LAND USE DATE INCONSISTENCY')
ENDIF
!
IF (LHOOK) CALL DR_HOOK('ALLOC_ISBA_LANDUSE:CHECK_DATE',1,ZHOOK_HANDLE)
!
END SUBROUTINE CHECK_DATE
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE ALLOC_ISBA_LANDUSE
