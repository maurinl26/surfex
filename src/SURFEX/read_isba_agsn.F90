!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!     #########
      SUBROUTINE READ_ISBA_AGS_n (IO, NP, NPE, NPGLO, NPEGLO,          &
                                  HPROGRAM, KI, KVERSION, KBUGFIX, ODIM)
!     ##################################
!
!!****  *READ_ISBA_AGS_n* - routine to initialise ISBA variables for interactive vegetation
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
USE MODD_CO2V_PAR,       ONLY : XANFMINIT
!
USE MODD_SURF_PAR,       ONLY : XUNDEF, NUNDEF, LEN_HREC
!
USE MODE_READ_SURF_LAYERS
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
INTEGER,              INTENT(IN)    :: KVERSION
INTEGER,              INTENT(IN)    :: KBUGFIX
LOGICAL,              INTENT(IN)    :: ODIM
!
!*       0.2   Declarations of local variables
!              -------------------------------
!
TYPE(ISBA_P_t), POINTER  :: PK, PKGLO
TYPE(ISBA_PE_t), POINTER :: PEK, PEKGLO
!
!
CHARACTER(LEN=LEN_HREC) :: YRECFM         ! Name of the article to be read
CHARACTER(LEN=4)        :: YLVL
!
REAL, DIMENSION(:,:,:),ALLOCATABLE :: ZWORK3D  ! 3D array to write data in file
REAL, DIMENSION(:,:)  ,ALLOCATABLE :: ZWORK2D  ! 2D array to write data in file
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
IF (LHOOK) CALL DR_HOOK('READ_ISBA_AGS_N',0,ZHOOK_HANDLE)
!-------------------------------------------------------------------------------
!
!*       1.     Physical dimension
!               -----------------
!
CALL ALLOC_READ_ISBA_AGS
!
!-------------------------------------------------------------------------------
ALLOCATE(ZWORK2D(KI,IO%NPATCH))
!-------------------------------------------------------------------------------
!
!*       2.  Assimilation
!           ---------
!
!* Assimilation
!
YRECFM = 'AN'
CALL MAKE_CHOICE_ARRAY(HPROGRAM, IO%NPATCH, ODIM, YRECFM, ZWORK2D)
DO JP = 1,IO%NPATCH
   CALL PACK_SAME_RANK(NP%AL(JP)%NR_P,ZWORK2D(:,JP),NPE%AL(JP)%XAN(:))
ENDDO
!
YRECFM = 'ANDAY'
CALL MAKE_CHOICE_ARRAY(HPROGRAM, IO%NPATCH, ODIM, YRECFM, ZWORK2D)
DO JP = 1,IO%NPATCH
   CALL PACK_SAME_RANK(NP%AL(JP)%NR_P,ZWORK2D(:,JP),NPE%AL(JP)%XANDAY(:))
ENDDO  
!
YRECFM = 'ANFM'
CALL MAKE_CHOICE_ARRAY(HPROGRAM, IO%NPATCH, ODIM, YRECFM, ZWORK2D)
DO JP = 1,IO%NPATCH
   CALL PACK_SAME_RANK(NP%AL(JP)%NR_P,ZWORK2D(:,JP),NPE%AL(JP)%XANFM(:))
ENDDO  
!
!-------------------------------------------------------------------------------
!
!*       3.   LAI, Biomass and Respiration
!           -------------------------
!
IF (IO%CPHOTO=='NIT'.OR.IO%CPHOTO=='NCB') THEN
  !
  !* LAI
  !
  YRECFM = 'LAI'
  CALL MAKE_CHOICE_ARRAY(HPROGRAM, IO%NPATCH, ODIM, YRECFM, ZWORK2D)  
  DO JP = 1,IO%NPATCH
     PK  => NP%AL(JP)
     PEK => NPE%AL(JP)
     CALL PACK_SAME_RANK(PK%NR_P,ZWORK2D(:,JP),PEK%XLAI(:))
     IF(IO%LLULU) NPEGLO%AL(JP)%XLAI(:)=ZWORK2D(:,JP)
  ENDDO
  !
  !* Biomass
  !
  IF (KVERSION>7 .OR. KVERSION==7 .AND. KBUGFIX>=3) THEN
    YRECFM='BIOMA'
  ELSE
    YRECFM='BIOMASS'
  ENDIF
  !
  ALLOCATE(ZWORK3D(KI,IO%NNBIOMASS,IO%NPATCH))    
  !  
  CALL READ_SURF_LAYERS(HPROGRAM,YRECFM,ODIM,ZWORK3D,IRESP)
  !
  DO JNB=1,IO%NNBIOMASS
     DO JP = 1,IO%NPATCH
        PK => NP%AL(JP)
        PEK => NPE%AL(JP)
        CALL PACK_SAME_RANK(PK%NR_P,ZWORK3D(:,JNB,JP),PEK%XBIOMASS(:,JNB))
        IF(IO%LLULU) NPEGLO%AL(JP)%XBIOMASS(:,JNB)=ZWORK3D(:,JNB,JP)
     ENDDO 
  ENDDO 
  !
  DEALLOCATE(ZWORK3D)
  !
  !* Respiration
  !
  IWORK=0
  IF(IO%CPHOTO=='NCB'.OR.KVERSION<8)THEN
    IWORK=2
  ENDIF
  !
  DO JNB=2,IO%NNBIOMASS-IWORK
     !
     WRITE(YLVL,'(I1)') JNB
     !
     IF (KVERSION>7 .OR. (KVERSION==7 .AND. KBUGFIX>=3)) THEN
        YRECFM='RESPI'//ADJUSTL(YLVL(:LEN_TRIM(YLVL)))
     ELSE
        YRECFM='RESP_BIOM'//ADJUSTL(YLVL(:LEN_TRIM(YLVL)))
     ENDIF    
     !
     CALL MAKE_CHOICE_ARRAY(HPROGRAM, IO%NPATCH, ODIM, YRECFM, ZWORK2D)
     !
     DO JP = 1,IO%NPATCH
        PK => NP%AL(JP)
        PEK => NPE%AL(JP)
        CALL PACK_SAME_RANK(PK%NR_P,ZWORK2D(:,JP),PEK%XRESP_BIOMASS(:,JNB))
     ENDDO
     !
  ENDDO
  !
ENDIF
!
!-------------------------------------------------------------------------------
DEALLOCATE(ZWORK2D)
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('READ_ISBA_AGS_N',1,ZHOOK_HANDLE)
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
CONTAINS
!-------------------------------------------------------------------------------
!
SUBROUTINE ALLOC_READ_ISBA_AGS
!
IMPLICIT NONE
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
IF (LHOOK) CALL DR_HOOK('READ_ISBA_AGS_N:ALLOC_READ_ISBA_AGS',0,ZHOOK_HANDLE)
!
!* Assimilation
!        
DO JP = 1,IO%NPATCH
   !
   PK => NP%AL(JP)
   PEK => NPE%AL(JP)
   ! 
   ALLOCATE(PEK%XANFM  (PK%NSIZE_P))
   ALLOCATE(PEK%XAN    (PK%NSIZE_P))
   ALLOCATE(PEK%XANDAY (PK%NSIZE_P))
   PEK%XANFM (:) = XANFMINIT  
   PEK%XAN   (:) = 0.
   PEK%XANDAY(:) = 0.
   !
ENDDO
!
!* Alloc biomass and respiration (and global read LAI)
!
IF(IO%CPHOTO=='NIT'.OR.IO%CPHOTO=='NCB')THEN
  !
  DO JP = 1,IO%NPATCH
     !
     PK => NP%AL(JP)
     PEK => NPE%AL(JP)
     !
     ALLOCATE(PEK%XBIOMASS     (PK%NSIZE_P,IO%NNBIOMASS))
     ALLOCATE(PEK%XRESP_BIOMASS(PK%NSIZE_P,IO%NNBIOMASS))
     !
     PEK%XBIOMASS     (:,:) = 0.
     PEK%XRESP_BIOMASS(:,:) = 0.    
     !
  ENDDO
  !
  IF(IO%LLULU)THEN
    DO JP = 1,IO%NPATCH
       PKGLO => NPGLO%AL(JP)
       PEKGLO => NPEGLO%AL(JP)
       ALLOCATE(PEKGLO%XLAI    (PKGLO%NSIZE_P)             )
       ALLOCATE(PEKGLO%XBIOMASS(PKGLO%NSIZE_P,IO%NNBIOMASS))
       PEKGLO%XLAI    (:  ) = 0. 
       PEKGLO%XBIOMASS(:,:) = 0. 
    ENDDO
  ENDIF
  !
ENDIF
!
IF (LHOOK) CALL DR_HOOK('READ_ISBA_AGS_N:ALLOC_READ_ISBA_AGS',1,ZHOOK_HANDLE)
!
END SUBROUTINE ALLOC_READ_ISBA_AGS
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE READ_ISBA_AGS_n
