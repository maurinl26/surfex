!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!     #########
      SUBROUTINE READ_ISBA_MEB_n (IO, NP, NPE, NPGLO, NPEGLO, &
                                  HPROGRAM, KI, ODIM          )
!     ##################################
!
!!****  *READ_ISBA_MEB_n* - routine to initialise ISBA physicals variables
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
USE MODD_ISBA_n,         ONLY : ISBA_NP_t, ISBA_NPE_t, ISBA_P_t, ISBA_PE_t
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
TYPE(ISBA_NP_t),      INTENT(INOUT) :: NP, NPGLO
TYPE(ISBA_NPE_t),     INTENT(INOUT) :: NPE, NPEGLO
!
CHARACTER(LEN=6),     INTENT(IN)    :: HPROGRAM ! calling program
INTEGER,              INTENT(IN)    :: KI       ! number of points
LOGICAL,              INTENT(IN)    :: ODIM
!
!*       0.2   Declarations of local variables
!              -------------------------------
!
TYPE(ISBA_P_t), POINTER  :: PK
TYPE(ISBA_PE_t), POINTER :: PEK
!
CHARACTER(LEN=LEN_HREC) :: YRECFM         ! Name of the article to be read
!
REAL, DIMENSION(:,:)  ,ALLOCATABLE :: ZWORK2D  ! 2D array to write data in file
!
INTEGER :: JP      ! loop counter
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('READ_ISBA_MEB_N',0,ZHOOK_HANDLE)
!-------------------------------------------------------------------------------
!
!*       1.     Physical dimension
!               -----------------
!
CALL ALLOC_READ_ISBA_MEB(NP,NPE)
IF(IO%LLULU)THEN
  CALL ALLOC_READ_ISBA_MEB(NPGLO,NPEGLO)        
ENDIF
!
!-------------------------------------------------------------------------------
ALLOCATE(ZWORK2D(KI,IO%NPATCH))
!-------------------------------------------------------------------------------
!
!
!*       2.  MEB Prognostic or Semi-prognostic variables
!            -------------------------------------------
!
!* water intercepted on litter
!
YRECFM = 'WRL'
CALL MAKE_CHOICE_ARRAY(HPROGRAM, IO%NPATCH, ODIM, YRECFM, ZWORK2D)
!
DO JP = 1,IO%NPATCH
   PK => NP%AL(JP)
   PEK => NPE%AL(JP)
   CALL PACK_SAME_RANK(PK%NR_P,ZWORK2D(:,JP),PEK%XWRL(:))
   IF(IO%LLULU) NPEGLO%AL(JP)%XWRL(:)=ZWORK2D(:,JP)
ENDDO
!
YRECFM = 'WRLI'
CALL MAKE_CHOICE_ARRAY(HPROGRAM, IO%NPATCH, ODIM, YRECFM, ZWORK2D)
!
DO JP = 1,IO%NPATCH
   PK => NP%AL(JP)
   PEK => NPE%AL(JP)
   CALL PACK_SAME_RANK(PK%NR_P,ZWORK2D(:,JP),PEK%XWRLI(:))
   IF(IO%LLULU) NPEGLO%AL(JP)%XWRLI(:)=ZWORK2D(:,JP)
ENDDO 
!
!* snow intercepted on vegetation canopy leaves
!
YRECFM = 'WRVN'
CALL MAKE_CHOICE_ARRAY(HPROGRAM, IO%NPATCH, ODIM, YRECFM, ZWORK2D)
!
DO JP = 1,IO%NPATCH
   PK => NP%AL(JP)
   PEK => NPE%AL(JP)
   CALL PACK_SAME_RANK(PK%NR_P,ZWORK2D(:,JP),PEK%XWRVN(:))
   IF(IO%LLULU) NPEGLO%AL(JP)%XWRVN(:)=ZWORK2D(:,JP)
ENDDO   
!
!* vegetation canopy temperature
!
YRECFM = 'TV'
CALL MAKE_CHOICE_ARRAY(HPROGRAM, IO%NPATCH, ODIM, YRECFM, ZWORK2D)
!
DO JP = 1,IO%NPATCH
   PK => NP%AL(JP)
   PEK => NPE%AL(JP)
   CALL PACK_SAME_RANK(PK%NR_P,ZWORK2D(:,JP),PEK%XTV(:))
   IF(IO%LLULU) NPEGLO%AL(JP)%XTV(:)=ZWORK2D(:,JP)
ENDDO    
!
!* litter temperature
!
YRECFM = 'TL'
CALL MAKE_CHOICE_ARRAY(HPROGRAM, IO%NPATCH, ODIM, YRECFM, ZWORK2D)
!
DO JP = 1,IO%NPATCH
   PK => NP%AL(JP)
   PEK => NPE%AL(JP)
   CALL PACK_SAME_RANK(PK%NR_P,ZWORK2D(:,JP),PEK%XTL(:))
   IF(IO%LLULU) NPEGLO%AL(JP)%XTL(:)=ZWORK2D(:,JP)
ENDDO  
!
!* vegetation canopy air temperature
!
YRECFM = 'TC'
CALL MAKE_CHOICE_ARRAY(HPROGRAM, IO%NPATCH, ODIM, YRECFM, ZWORK2D)
!
DO JP = 1,IO%NPATCH
   PK => NP%AL(JP)
   PEK => NPE%AL(JP)
   CALL PACK_SAME_RANK(PK%NR_P,ZWORK2D(:,JP),PEK%XTC(:))
   IF(IO%LLULU) NPEGLO%AL(JP)%XTC(:)=ZWORK2D(:,JP)
ENDDO    
!
!-------------------------------------------------------------------------------
DEALLOCATE(ZWORK2D)
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('READ_ISBA_MEB_N',1,ZHOOK_HANDLE)
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
CONTAINS
!-------------------------------------------------------------------------------
!
SUBROUTINE ALLOC_READ_ISBA_MEB(NA,NAE)
!
IMPLICIT NONE
!
TYPE(ISBA_NP_t),      INTENT(INOUT) :: NA
TYPE(ISBA_NPE_t),     INTENT(INOUT) :: NAE
!
TYPE(ISBA_P_t),  POINTER :: PA
TYPE(ISBA_PE_t), POINTER :: PEA
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
IF (LHOOK) CALL DR_HOOK('READ_ISBA_MEB_N:ALLOC_READ_ISBA_MEB',0,ZHOOK_HANDLE)
!
DO JP = 1,IO%NPATCH
   !
   PEA => NAE%AL(JP)
   PA => NA%AL(JP)
   !
   ALLOCATE(PEA%XWRL (PA%NSIZE_P))
   ALLOCATE(PEA%XWRLI(PA%NSIZE_P))
   ALLOCATE(PEA%XWRVN(PA%NSIZE_P))
   ALLOCATE(PEA%XTV  (PA%NSIZE_P))
   ALLOCATE(PEA%XTL  (PA%NSIZE_P))
   ALLOCATE(PEA%XTC  (PA%NSIZE_P))
   !
   PEA%XWRL (:) = XUNDEF
   PEA%XWRLI(:) = XUNDEF
   PEA%XWRVN(:) = XUNDEF
   PEA%XTV  (:) = XUNDEF
   PEA%XTL  (:) = XUNDEF
   PEA%XTC  (:) = XUNDEF
   !
ENDDO
!
IF (LHOOK) CALL DR_HOOK('READ_ISBA_MEB_N:ALLOC_READ_ISBA_MEB',1,ZHOOK_HANDLE)
!
END SUBROUTINE ALLOC_READ_ISBA_MEB
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE READ_ISBA_MEB_n
