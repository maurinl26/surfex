!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!     #########
      SUBROUTINE READ_ISBA_CC_n (IO, NP, NPE, NPGLO, NPEGLO, &
                                 HPROGRAM, KI, KVERSION, ODIM)
!     ##################################
!
!!****  *READ_ISBA_CC_n* - routine to initialise ISBA variables for Carbon Cycle
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
USE MODD_SURF_PAR,       ONLY : XUNDEF, NUNDEF, LEN_HREC
!
USE MODE_READ_SURF_LAYERS
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
TYPE(ISBA_NP_t),      INTENT(INOUT) :: NP, NPGLO
TYPE(ISBA_NPE_t),     INTENT(INOUT) :: NPE, NPEGLO
!
CHARACTER(LEN=6),     INTENT(IN)    :: HPROGRAM ! calling program
INTEGER,              INTENT(IN)    :: KI       ! number of points
INTEGER,              INTENT(IN)    :: KVERSION
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
IF (LHOOK) CALL DR_HOOK('READ_ISBA_CC_N',0,ZHOOK_HANDLE)
!-------------------------------------------------------------------------------
!
!*       1.     Physical dimension
!               -----------------
!
CALL ALLOC_READ_ISBA_CC(NP,NPE)
IF(IO%LLULU)THEN
  CALL ALLOC_READ_ISBA_CC(NPGLO,NPEGLO)        
ENDIF
!
!-------------------------------------------------------------------------------
ALLOCATE(ZWORK2D(KI,IO%NPATCH))
!-------------------------------------------------------------------------------
!
!
!*       2. Bulk Soil carbon
!        -------------------
!
!
IF(IO%CRESPSL=='CNT')THEN
  !
  DO JNL=1,IO%NNLITTER
    DO JNS=1,IO%NNLITTLEVS
      WRITE(YLVL,'(I1,A1,I1)') JNL,'_',JNS
      YRECFM='LITTER'//ADJUSTL(YLVL(:LEN_TRIM(YLVL)))
      CALL MAKE_CHOICE_ARRAY(HPROGRAM, IO%NPATCH, ODIM, YRECFM, ZWORK2D)
      DO JP = 1,IO%NPATCH
        CALL PACK_SAME_RANK(NP%AL(JP)%NR_P,ZWORK2D(:,JP),NPE%AL(JP)%XLITTER(:,JNL,JNS))
        IF(IO%LLULU) NPEGLO%AL(JP)%XLITTER(:,JNL,JNS)=ZWORK2D(:,JP)
      ENDDO       
    END DO
  END DO
  !
  DO JNS=1,IO%NNSOILCARB
    WRITE(YLVL,'(I4)') JNS
    YRECFM='SOILCARB'//ADJUSTL(YLVL(:LEN_TRIM(YLVL)))
    CALL MAKE_CHOICE_ARRAY(HPROGRAM, IO%NPATCH, ODIM, YRECFM, ZWORK2D)
    DO JP = 1,IO%NPATCH
      CALL PACK_SAME_RANK(NP%AL(JP)%NR_P,ZWORK2D(:,JP),NPE%AL(JP)%XSOILCARB(:,JNS))
      IF(IO%LLULU) NPEGLO%AL(JP)%XSOILCARB(:,JNS)=ZWORK2D(:,JP)
    ENDDO      
  END DO
  !
  DO JNS=1,IO%NNLITTLEVS
    WRITE(YLVL,'(I4)') JNS
    YRECFM='LIGN_STR'//ADJUSTL(YLVL(:LEN_TRIM(YLVL)))
    CALL MAKE_CHOICE_ARRAY(HPROGRAM, IO%NPATCH, ODIM, YRECFM, ZWORK2D)
    DO JP = 1,IO%NPATCH
      CALL PACK_SAME_RANK(NP%AL(JP)%NR_P,ZWORK2D(:,JP),NPE%AL(JP)%XLIGNIN_STRUC(:,JNS))
      IF(IO%LLULU) NPEGLO%AL(JP)%XLIGNIN_STRUC(:,JNS)=ZWORK2D(:,JP)
    ENDDO     
  END DO
  !
ENDIF
!
!-------------------------------------------------------------------------------
!
!
!*       3. Multi-layer Soil carbon
!        --------------------------
!
!
IF(IO%CRESPSL=='DIF'.AND.KVERSION>=9)THEN
  !
  YRECFM='SFLIGN'
  CALL MAKE_CHOICE_ARRAY(HPROGRAM, IO%NPATCH, ODIM, YRECFM, ZWORK2D)
  DO JP = 1,IO%NPATCH
     CALL PACK_SAME_RANK(NP%AL(JP)%NR_P,ZWORK2D(:,JP),NPE%AL(JP)%XSURFACE_LIGNIN_STRUC(:))
     IF(IO%LLULU) NPEGLO%AL(JP)%XSURFACE_LIGNIN_STRUC(:)=ZWORK2D(:,JP)
  ENDDO
  !
  DO JNL=1,IO%NNLITTER
     WRITE(YLVL,'(I1)') JNL
     YRECFM='SFLIT'//ADJUSTL(YLVL(:LEN_TRIM(YLVL)))
     CALL MAKE_CHOICE_ARRAY(HPROGRAM, IO%NPATCH, ODIM, YRECFM, ZWORK2D)
     DO JP = 1,IO%NPATCH
        CALL PACK_SAME_RANK(NP%AL(JP)%NR_P,ZWORK2D(:,JP),NPE%AL(JP)%XSURFACE_LITTER(:,JNL))
        IF(IO%LLULU) NPEGLO%AL(JP)%XSURFACE_LITTER(:,JNL)=ZWORK2D(:,JP)
     ENDDO
  END DO
  !
  ALLOCATE(ZWORK3D(KI,IO%NGROUND_LAYER,IO%NPATCH))
  !
  CALL READ_SURF_LAYERS(HPROGRAM,'DFLIGN',ODIM,ZWORK3D,IRESP)
  DO JL=1,IO%NGROUND_LAYER
     DO JP = 1,IO%NPATCH
        CALL PACK_SAME_RANK(NP%AL(JP)%NR_P,ZWORK3D(:,JL,JP),NPE%AL(JP)%XSOILDIF_LIGNIN_STRUC(:,JL))
        IF(IO%LLULU) NPEGLO%AL(JP)%XSOILDIF_LIGNIN_STRUC(:,JL)=ZWORK3D(:,JL,JP)
     ENDDO
  ENDDO
  !
  CALL READ_SURF_LAYERS(HPROGRAM,'DFLIT1L',ODIM,ZWORK3D,IRESP)
  DO JL=1,IO%NGROUND_LAYER
     DO JP = 1,IO%NPATCH
        CALL PACK_SAME_RANK(NP%AL(JP)%NR_P,ZWORK3D(:,JL,JP),NPE%AL(JP)%XSOILDIF_LITTER(:,JL,1))
        IF(IO%LLULU) NPEGLO%AL(JP)%XSOILDIF_LITTER(:,JL,1)=ZWORK3D(:,JL,JP)
     ENDDO
  ENDDO
  CALL READ_SURF_LAYERS(HPROGRAM,'DFLIT2L',ODIM,ZWORK3D,IRESP)
  DO JL=1,IO%NGROUND_LAYER
     DO JP = 1,IO%NPATCH
        CALL PACK_SAME_RANK(NP%AL(JP)%NR_P,ZWORK3D(:,JL,JP),NPE%AL(JP)%XSOILDIF_LITTER(:,JL,2))
        IF(IO%LLULU) NPEGLO%AL(JP)%XSOILDIF_LITTER(:,JL,2)=ZWORK3D(:,JL,JP)
     ENDDO
  ENDDO
  !
  CALL READ_SURF_LAYERS(HPROGRAM,'DFSOC1L',ODIM,ZWORK3D,IRESP)
  DO JL=1,IO%NGROUND_LAYER
     DO JP = 1,IO%NPATCH
        CALL PACK_SAME_RANK(NP%AL(JP)%NR_P,ZWORK3D(:,JL,JP),NPE%AL(JP)%XSOILDIF_CARB(:,JL,1))
        IF(IO%LLULU) NPEGLO%AL(JP)%XSOILDIF_CARB(:,JL,1)=ZWORK3D(:,JL,JP)
     ENDDO
  ENDDO
  CALL READ_SURF_LAYERS(HPROGRAM,'DFSOC2L',ODIM,ZWORK3D,IRESP)
  DO JL=1,IO%NGROUND_LAYER
     DO JP = 1,IO%NPATCH
        CALL PACK_SAME_RANK(NP%AL(JP)%NR_P,ZWORK3D(:,JL,JP),NPE%AL(JP)%XSOILDIF_CARB(:,JL,2))
        IF(IO%LLULU) NPEGLO%AL(JP)%XSOILDIF_CARB(:,JL,2)=ZWORK3D(:,JL,JP)
     ENDDO
  ENDDO
  CALL READ_SURF_LAYERS(HPROGRAM,'DFSOC3L',ODIM,ZWORK3D,IRESP)
  DO JL=1,IO%NGROUND_LAYER
     DO JP = 1,IO%NPATCH
        CALL PACK_SAME_RANK(NP%AL(JP)%NR_P,ZWORK3D(:,JL,JP),NPE%AL(JP)%XSOILDIF_CARB(:,JL,3))
        IF(IO%LLULU) NPEGLO%AL(JP)%XSOILDIF_CARB(:,JL,3)=ZWORK3D(:,JL,JP)
     ENDDO
  ENDDO
  !
  DEALLOCATE(ZWORK3D)
  !
ENDIF
!
!-------------------------------------------------------------------------------
!
!
!*       4. Multi-layer Soil gas
!        -----------------------
!
!
IF(IO%LSOILGAS.AND.KVERSION>=9)THEN
  !
  ALLOCATE(ZWORK3D(KI,IO%NGROUND_LAYER,IO%NPATCH))
  !
  CALL READ_SURF_LAYERS(HPROGRAM,'GASO2L',ODIM,ZWORK3D,IRESP)
  DO JL=1,IO%NGROUND_LAYER
     DO JP = 1,IO%NPATCH
        CALL PACK_SAME_RANK(NP%AL(JP)%NR_P,ZWORK3D(:,JL,JP),NPE%AL(JP)%XSGASO2(:,JL))
        IF(IO%LLULU) NPEGLO%AL(JP)%XSGASO2(:,JL)=ZWORK3D(:,JL,JP)
     ENDDO
  ENDDO
  !
  CALL READ_SURF_LAYERS(HPROGRAM,'GASCO2L',ODIM,ZWORK3D,IRESP)
  DO JL=1,IO%NGROUND_LAYER
     DO JP = 1,IO%NPATCH
        CALL PACK_SAME_RANK(NP%AL(JP)%NR_P,ZWORK3D(:,JL,JP),NPE%AL(JP)%XSGASCO2(:,JL))
        IF(IO%LLULU) NPEGLO%AL(JP)%XSGASCO2(:,JL)=ZWORK3D(:,JL,JP)
     ENDDO
  ENDDO
  !
  CALL READ_SURF_LAYERS(HPROGRAM,'GASCH4L',ODIM,ZWORK3D,IRESP)
  DO JL=1,IO%NGROUND_LAYER
     DO JP = 1,IO%NPATCH
        CALL PACK_SAME_RANK(NP%AL(JP)%NR_P,ZWORK3D(:,JL,JP),NPE%AL(JP)%XSGASCH4(:,JL))
        IF(IO%LLULU) NPEGLO%AL(JP)%XSGASCH4(:,JL)=ZWORK3D(:,JL,JP)
     ENDDO
  ENDDO
  !
  DEALLOCATE(ZWORK3D)
  !
ENDIF
!
!-------------------------------------------------------------------------------
!
!
!*       5. Fire scheme
!        --------------
!
!
IF(IO%LFIRE.AND.KVERSION>=9)THEN
  !
  YRECFM='FIREIND'
  CALL MAKE_CHOICE_ARRAY(HPROGRAM, IO%NPATCH, ODIM, YRECFM, ZWORK2D)
  DO JP = 1,IO%NPATCH
     CALL PACK_SAME_RANK(NP%AL(JP)%NR_P,ZWORK2D(:,JP),NPE%AL(JP)%XFIREIND(:))
     IF(IO%LLULU) NPEGLO%AL(JP)%XFIREIND(:)=ZWORK2D(:,JP)
  ENDDO
  !
  YRECFM='MOISTLITFIRE'
  CALL MAKE_CHOICE_ARRAY(HPROGRAM, IO%NPATCH, ODIM, YRECFM, ZWORK2D)
  DO JP = 1,IO%NPATCH
     CALL PACK_SAME_RANK(NP%AL(JP)%NR_P,ZWORK2D(:,JP),NPE%AL(JP)%XMOISTLIT_FIRE(:))
     IF(IO%LLULU) NPEGLO%AL(JP)%XMOISTLIT_FIRE(:)=ZWORK2D(:,JP)
  ENDDO
  !
  YRECFM='TEMPLITFIRE'
  CALL MAKE_CHOICE_ARRAY(HPROGRAM, IO%NPATCH, ODIM, YRECFM, ZWORK2D)
  DO JP = 1,IO%NPATCH
     CALL PACK_SAME_RANK(NP%AL(JP)%NR_P,ZWORK2D(:,JP),NPE%AL(JP)%XTEMPLIT_FIRE(:))
     IF(IO%LLULU) NPEGLO%AL(JP)%XTEMPLIT_FIRE(:)=ZWORK2D(:,JP)
  ENDDO
  !
ENDIF
!
!-------------------------------------------------------------------------------
DEALLOCATE(ZWORK2D)
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('READ_ISBA_CC_N',1,ZHOOK_HANDLE)
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
CONTAINS
!-------------------------------------------------------------------------------
!
SUBROUTINE ALLOC_READ_ISBA_CC(NA,NAE)
!
IMPLICIT NONE
!
TYPE(ISBA_NP_t),  INTENT(INOUT) :: NA
TYPE(ISBA_NPE_t), INTENT(INOUT) :: NAE
!
TYPE(ISBA_P_t),  POINTER :: PA
TYPE(ISBA_PE_t), POINTER :: PEA
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
IF (LHOOK) CALL DR_HOOK('READ_ISBA_CC_N:ALLOC_READ_ISBA_CC',0,ZHOOK_HANDLE)
!
! * Soil carbon
!
IF(IO%CRESPSL=='CNT')THEN
  !
  DO JP = 1,IO%NPATCH
     !
     PEA => NAE%AL(JP)
     PA => NA%AL(JP)
     !
     ALLOCATE(PEA%XLITTER      (PA%NSIZE_P,IO%NNLITTER,IO%NNLITTLEVS))
     ALLOCATE(PEA%XSOILCARB    (PA%NSIZE_P,IO%NNSOILCARB))
     ALLOCATE(PEA%XLIGNIN_STRUC(PA%NSIZE_P,IO%NNLITTLEVS))
     !
     PEA%XLITTER      (:,:,:) = 0.
     PEA%XSOILCARB    (:,:) = 0. 
     PEA%XLIGNIN_STRUC(:,:) = 0.  
     !
  ENDDO
  !
ELSEIF (IO%CRESPSL=='DIF') THEN
  !
  DO JP = 1,IO%NPATCH
     !
     PEA => NAE%AL(JP)
     PA => NA%AL(JP)
     !
     ALLOCATE(PEA%XSURFACE_LIGNIN_STRUC(PA%NSIZE_P))
     ALLOCATE(PEA%XSURFACE_LITTER      (PA%NSIZE_P,IO%NNLITTER))
     PEA%XSURFACE_LIGNIN_STRUC(:  ) = 0.0
     PEA%XSURFACE_LITTER      (:,:) = 0.0
     !
     ALLOCATE(PEA%XSOILDIF_LIGNIN_STRUC(PA%NSIZE_P,IO%NGROUND_LAYER))
     ALLOCATE(PEA%XSOILDIF_LITTER      (PA%NSIZE_P,IO%NGROUND_LAYER,IO%NNLITTER))
     ALLOCATE(PEA%XSOILDIF_CARB        (PA%NSIZE_P,IO%NGROUND_LAYER,IO%NNSOILCARB))
     PEA%XSOILDIF_LIGNIN_STRUC(:,:  ) = 0.0
     PEA%XSOILDIF_LITTER      (:,:,:) = 0.0
     PEA%XSOILDIF_CARB        (:,:,:) = 0.0
     !   
  ENDDO
  !
ENDIF
!
! * Soil gas
!
IF(IO%CRESPSL=='DIF'.AND.IO%LSOILGAS)THEN
  !
  DO JP = 1,IO%NPATCH
     !
     PEA => NAE%AL(JP)
     PA => NA%AL(JP)
     !
     ALLOCATE(PEA%XSGASO2 (PA%NSIZE_P,IO%NGROUND_LAYER))
     ALLOCATE(PEA%XSGASCO2(PA%NSIZE_P,IO%NGROUND_LAYER))
     ALLOCATE(PEA%XSGASCH4(PA%NSIZE_P,IO%NGROUND_LAYER))
     !
     PEA%XSGASO2 (:,:) = 0.0
     PEA%XSGASCO2(:,:) = 0.0
     PEA%XSGASCH4(:,:) = 0.0
     !
  ENDDO
  !
ENDIF 
!
! * Fire scheme
!
IF(IO%LFIRE)THEN
  !
  DO JP = 1,IO%NPATCH
     !
     PEA => NAE%AL(JP)
     PA => NA%AL(JP)
     !
     ALLOCATE(PEA%XFIREIND      (PA%NSIZE_P))
     ALLOCATE(PEA%XMOISTLIT_FIRE(PA%NSIZE_P))
     ALLOCATE(PEA%XTEMPLIT_FIRE (PA%NSIZE_P))
     !
     PEA%XFIREIND      (:) = 0.0
     PEA%XMOISTLIT_FIRE(:) = PEA%XWG(:,1)+PEA%XWGI(:,1)
     PEA%XTEMPLIT_FIRE (:) = PEA%XTG(:,1)
     !
  ENDDO
  !
ENDIF
!
IF (LHOOK) CALL DR_HOOK('READ_ISBA_CC_N:ALLOC_READ_ISBA_CC',1,ZHOOK_HANDLE)
!
END SUBROUTINE ALLOC_READ_ISBA_CC
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE READ_ISBA_CC_n
