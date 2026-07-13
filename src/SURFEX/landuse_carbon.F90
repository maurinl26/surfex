!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!     #########
    SUBROUTINE LANDUSE_CARBON(IO, S, NP, NPE, NPGLO, NPEGLO, TLU, KI)
!   #################################################################
!!****  *LAND USE CARBON*
!!
!!    PURPOSE
!!    -------
!
!     Update and conserv carbon stocks after land-use change for ISBA-CC
!               
!!**  METHOD
!!    ------
!!
!!    EXTERNAL
!!    --------
!!    none
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      
!!    none
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
!!      Original    11/2020    
!!
!-------------------------------------------------------------------------------
!
!*       0.     DECLARATIONS
!               ------------
!
USE MODD_ISBA_OPTIONS_n, ONLY : ISBA_OPTIONS_t
USE MODD_ISBA_n,         ONLY : ISBA_S_t, ISBA_P_t, ISBA_PE_t, ISBA_K_t, ISBA_NP_t, ISBA_NPE_t
USE MODD_INIT_LANDUSE,   ONLY : LULCC_t
!
USE MODD_CO2V_PAR, ONLY : XGTOKG
!
USE MODI_BIOMASS_TO_SURFACE_LITTER
USE MODI_BIOMASS_TO_SOIL_LITTER
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
TYPE(ISBA_OPTIONS_t),  INTENT(INOUT) :: IO
TYPE(ISBA_S_t),        INTENT(INOUT) :: S
TYPE(ISBA_NP_t),       INTENT(INOUT) :: NP, NPGLO
TYPE(ISBA_NPE_t),      INTENT(INOUT) :: NPE, NPEGLO
TYPE(LULCC_t),         INTENT(INOUT) :: TLU
!
INTEGER,               INTENT(IN)    :: KI
!
!*      0.2    declarations of local parameter
!
TYPE(ISBA_P_t),  POINTER :: PK, PKGLO
TYPE(ISBA_PE_t), POINTER :: PEK, PEKGLO
!
REAL, DIMENSION(KI)                :: ZOLDVEG_GRID
REAL, DIMENSION(KI)                :: ZNEWVEG_GRID
!
REAL, DIMENSION(KI)                :: ZCSURF_LIGNIN_NEW    ! current year surface lignin C stock
REAL, DIMENSION(KI)                :: ZCSURF_LIGNIN_OLD    ! previous year surface lignin C stock
REAL, DIMENSION(KI)                :: ZCSOIL_LIGNIN_NEW    ! current year lignin C stock
REAL, DIMENSION(KI)                :: ZCSOIL_LIGNIN_OLD    ! current year lignin C stock
!
REAL, DIMENSION(KI,IO%NNLITTER)   :: ZCSURF_LITTER_NEW    ! current year surface litter C stock
REAL, DIMENSION(KI,IO%NNLITTER)   :: ZCSURF_LITTER_OLD    ! previous year surface litter C stock
REAL, DIMENSION(KI,IO%NNLITTER)   :: ZCSOIL_LITTER_NEW    ! current year litter C stock
REAL, DIMENSION(KI,IO%NNLITTER)   :: ZCSOIL_LITTER_OLD    ! previous year litter C stock
!
REAL, DIMENSION(KI,IO%NNSOILCARB) :: ZCSOIL_RESERV_NEW  ! current  year soil C stock
REAL, DIMENSION(KI,IO%NNSOILCARB) :: ZCSOIL_RESERV_OLD  ! previous year soil C stock
!
REAL, DIMENSION(KI,1,IO%NNLITTER) :: ZSOIL_LITTER
REAL, DIMENSION(KI,1)             :: ZSOIL_LIGNIN_STRUC
!
REAL, DIMENSION(KI,IO%NNLITTER,IO%NPATCH) :: ZHARVEST
!
INTEGER :: INP, INL, INC, JGLO ! loop counter on levels
INTEGER :: JI, JP, JNL, JNC     ! loop counter on levels
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!-----------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('LANDUSE_CARBON',0,ZHOOK_HANDLE)
!-----------------------------------------------------------------
!
!*      1.     Preliminaries
!              -------------
!
INP=IO%NPATCH
INL=IO%NNLITTER
INC=IO%NNSOILCARB
!
ZHARVEST(:,:,:)=0.0
!
ZCSURF_LIGNIN_NEW(:)   = 0.0
ZCSURF_LIGNIN_OLD(:)   = 0.0
ZCSOIL_LIGNIN_NEW(:)   = 0.0
ZCSOIL_LIGNIN_OLD(:)   = 0.0
!
ZCSURF_LITTER_OLD(:,:) = 0.0
ZCSOIL_LITTER_OLD(:,:) = 0.0
ZCSOIL_RESERV_OLD(:,:) = 0.0
!
ZCSURF_LITTER_NEW(:,:) = 0.0
ZCSOIL_LITTER_NEW(:,:) = 0.0
ZCSOIL_RESERV_NEW(:,:) = 0.0
!
!* No carbon where new patch = 0.0
!
DO JP=1,INP
   PK => NP%AL(JP)
   PEK => NPE%AL(JP)
   DO JI=1,PK%NSIZE_P 
      IF(PK%XPATCH(JI)==0.0)THEN
        PEK%XLIGNIN_STRUC(JI,:  )=0.0
        PEK%XLITTER      (JI,:,:)=0.0
      ENDIF
   ENDDO
ENDDO
!
DO JP=1,INP  
   PK => NP%AL(JP)
   PEK => NPE%AL(JP)
   DO JNC=1,INC
      DO JI=1,PK%NSIZE_P
         IF(PK%XPATCH(JI)==0.0)THEN
           PEK%XSOILCARB(JI,JNC)=0.0
         ENDIF
      ENDDO
   ENDDO
ENDDO
!
!* Old carbon grid
!
ZOLDVEG_GRID(:) = 0.0
DO JP=1,INP  
   PKGLO => NPGLO%AL(JP)
   PEKGLO => NPEGLO%AL(JP)
   DO JGLO=1,KI
      IF(PEKGLO%XSOILCARB(JGLO,1)>0.0)THEN
        ZOLDVEG_GRID(JGLO)=ZOLDVEG_GRID(JGLO)+PKGLO%XPATCH(JGLO)
      ENDIF       
   ENDDO
ENDDO
!
!* New carbon grid
!
ZNEWVEG_GRID(:) = 0.0
DO JP=1,INP  
   PK => NP%AL(JP)
   PEK => NPE%AL(JP)
   DO JI=1,PK%NSIZE_P
      JGLO = PK%NR_P(JI)
      IF(PEK%XSOILCARB(JI,1)>0.0)THEN
        ZNEWVEG_GRID(JGLO)=ZNEWVEG_GRID(JGLO)+PK%XPATCH(JI)
      ENDIF       
   ENDDO
ENDDO
!
!* To compute carbon budget
!
DO JP=1,INP
   PKGLO => NPGLO%AL(JP)
   PEKGLO => NPEGLO%AL(JP)
   DO JNL=1,INL
      DO JGLO=1,KI 
         TLU%XLITTER_GRID_OLD(JGLO) = TLU%XLITTER_GRID_OLD(JGLO) + (PEKGLO%XLITTER(JGLO,JNL,1)+PEKGLO%XLITTER(JGLO,JNL,2)) &
                                                                 * PKGLO%XPATCH(JGLO)
      ENDDO
   ENDDO
ENDDO
!
DO JP=1,INP
   PKGLO => NPGLO%AL(JP)
   PEKGLO => NPEGLO%AL(JP)
   DO JNC=1,INC
      DO JGLO=1,KI
         TLU%XCSOIL_GRID_OLD(JGLO) = TLU%XCSOIL_GRID_OLD(JGLO) + PEKGLO%XSOILCARB(JGLO,JNC)*PKGLO%XPATCH(JGLO)
      ENDDO
   ENDDO
ENDDO
!
!-----------------------------------------------------------------
!
!*      2.     Patch disappear = surface litter is harvested
!              ---------------------------------------------
!
DO JP=1,INP  
   PKGLO => NPGLO%AL(JP)
   PEKGLO => NPEGLO%AL(JP)
   DO JNL=1,INL
      DO JGLO=1,KI 
         IF(PKGLO%XPATCH(JGLO)>0.0.AND.S%XPATCH(JGLO,JP)==0.0)THEN 
           ZHARVEST(JGLO,JNL,JP)=ZHARVEST(JGLO,JNL,JP)+PEKGLO%XLITTER(JGLO,JNL,1)*PKGLO%XPATCH(JGLO)
         ENDIF
      ENDDO
   ENDDO
ENDDO
!
!-----------------------------------------------------------------
!
!*      3.     Compute previous and current year total grid-cell biomass and carbon stocks
!              ---------------------------------------------------------------------------
!
!Carbon stock in gC/m2
!
DO JP=1,INP
   PKGLO => NPGLO%AL(JP)
   PEKGLO => NPEGLO%AL(JP)
   DO JGLO=1,KI 
      ZCSURF_LIGNIN_OLD(JGLO) = ZCSURF_LIGNIN_OLD(JGLO) + PEKGLO%XLIGNIN_STRUC(JGLO,1)*PKGLO%XPATCH(JGLO)
      ZCSOIL_LIGNIN_OLD(JGLO) = ZCSOIL_LIGNIN_OLD(JGLO) + PEKGLO%XLIGNIN_STRUC(JGLO,2)*PKGLO%XPATCH(JGLO)  
   ENDDO
ENDDO 
!
DO JP=1,INP
   PK => NP%AL(JP)
   PEK => NPE%AL(JP)
   DO JI=1,PK%NSIZE_P 
      JGLO = PK%NR_P(JI)
      ZCSURF_LIGNIN_NEW(JGLO) = ZCSURF_LIGNIN_NEW(JGLO) + PEK%XLIGNIN_STRUC(JI,1)*PK%XPATCH(JI)
      ZCSOIL_LIGNIN_NEW(JGLO) = ZCSOIL_LIGNIN_NEW(JGLO) + PEK%XLIGNIN_STRUC(JI,2)*PK%XPATCH(JI)         
   ENDDO
ENDDO
!
DO JP=1,INP
   PKGLO => NPGLO%AL(JP)
   PEKGLO => NPEGLO%AL(JP)
   DO JNL=1,INL
      DO JI=1,KI
         ZCSURF_LITTER_OLD(JI,JNL) = ZCSURF_LITTER_OLD(JI,JNL) + PEKGLO%XLITTER(JI,JNL,1)*PKGLO%XPATCH(JI) - ZHARVEST(JI,JNL,JP)
         ZCSOIL_LITTER_OLD(JI,JNL) = ZCSOIL_LITTER_OLD(JI,JNL) + PEKGLO%XLITTER(JI,JNL,2)*PKGLO%XPATCH(JI)
      ENDDO
   ENDDO
ENDDO
! 
DO JP=1,INP
   PK => NP%AL(JP)
   PEK => NPE%AL(JP)
   DO JNL=1,INL
      DO JI=1,PK%NSIZE_P 
         JGLO = PK%NR_P(JI)
         ZCSURF_LITTER_NEW(JGLO,JNL) = ZCSURF_LITTER_NEW(JGLO,JNL) + PEK%XLITTER(JI,JNL,1)*PK%XPATCH(JI)
         ZCSOIL_LITTER_NEW(JGLO,JNL) = ZCSOIL_LITTER_NEW(JGLO,JNL) + PEK%XLITTER(JI,JNL,2)*PK%XPATCH(JI)           
      ENDDO
   ENDDO
ENDDO
!
DO JP=1,INP
   PKGLO => NPGLO%AL(JP)
   PEKGLO => NPEGLO%AL(JP)
   DO JNC=1,INC
      DO JGLO=1,KI
         ZCSOIL_RESERV_OLD(JGLO,JNC) = ZCSOIL_RESERV_OLD(JGLO,JNC) + PEKGLO%XSOILCARB(JGLO,JNC)*PKGLO%XPATCH(JGLO)
      ENDDO
   ENDDO
ENDDO
!
DO JP=1,INP
   PK => NP%AL(JP)
   PEK => NPE%AL(JP)
   DO JNC=1,INC
      DO JI=1,PK%NSIZE_P
         JGLO = PK%NR_P(JI)
         ZCSOIL_RESERV_NEW(JGLO,JNC) = ZCSOIL_RESERV_NEW(JGLO,JNC) + PEK%XSOILCARB(JI,JNC)*PK%XPATCH(JI)
      ENDDO
   ENDDO
ENDDO
!
!-----------------------------------------------------------------
!
!*      4.      Ensures grid-cell conservation for carbon and litter stock
!              -----------------------------------------------------------
!
DO JP=1,INP
   PK => NP%AL(JP)
   PEK => NPE%AL(JP)
   DO JI=1,PK%NSIZE_P
      JGLO = PK%NR_P(JI)
      IF(ZCSURF_LIGNIN_NEW(JGLO)/=ZCSURF_LIGNIN_OLD(JGLO).AND.ZCSURF_LIGNIN_NEW(JGLO)/=0.0)THEN
        PEK%XLIGNIN_STRUC(JI,1) = PEK%XLIGNIN_STRUC(JI,1) * ZCSURF_LIGNIN_OLD(JGLO)/ZCSURF_LIGNIN_NEW(JGLO)
      ENDIF
      IF(ZCSOIL_LIGNIN_NEW(JGLO)/=ZCSOIL_LIGNIN_OLD(JGLO).AND.ZCSOIL_LIGNIN_NEW(JGLO)/=0.0)THEN
        PEK%XLIGNIN_STRUC(JI,2) = PEK%XLIGNIN_STRUC(JI,2) * ZCSOIL_LIGNIN_OLD(JGLO)/ZCSOIL_LIGNIN_NEW(JGLO)
      ENDIF
   ENDDO
ENDDO
!
DO JP=1,INP  
   PK => NP%AL(JP)
   PEK => NPE%AL(JP)
   DO JNL=1,INL
      DO JI=1,PK%NSIZE_P
         JGLO = PK%NR_P(JI)
         IF(ZCSURF_LITTER_NEW(JGLO,JNL)/=ZCSURF_LITTER_OLD(JGLO,JNL).AND.ZCSURF_LITTER_NEW(JGLO,JNL)/=0.0 ) THEN
           PEK%XLITTER(JI,JNL,1) = PEK%XLITTER(JI,JNL,1) * ZCSURF_LITTER_OLD(JGLO,JNL)/ZCSURF_LITTER_NEW(JGLO,JNL)
         ENDIF
         IF(ZCSOIL_LITTER_NEW(JGLO,JNL)/=ZCSOIL_LITTER_OLD(JGLO,JNL).AND.ZCSOIL_LITTER_NEW(JGLO,JNL)/=0.0 ) THEN
           PEK%XLITTER(JI,JNL,2) = PEK%XLITTER(JI,JNL,2) * ZCSOIL_LITTER_OLD(JGLO,JNL)/ZCSOIL_LITTER_NEW(JGLO,JNL)
         ENDIF
      ENDDO
   ENDDO
ENDDO
!
DO JP=1,INP 
   PK => NP%AL(JP)
   PEK => NPE%AL(JP)
   DO JNC=1,INC
      DO JI=1,PK%NSIZE_P
         JGLO = PK%NR_P(JI)
         IF(ZCSOIL_RESERV_NEW(JGLO,JNC)/=ZCSOIL_RESERV_OLD(JGLO,JNC).AND.ZCSOIL_RESERV_NEW(JGLO,JNC)/= 0.0)THEN
           PEK%XSOILCARB(JI,JNC) = PEK%XSOILCARB(JI,JNC) * ZCSOIL_RESERV_OLD(JGLO,JNC)/ZCSOIL_RESERV_NEW(JGLO,JNC)
         ENDIF
      ENDDO
   ENDDO
ENDDO
!
! If all vegetated patches disappeared -> try to conserv
!
DO JP=1,INP
   PKGLO => NPGLO%AL(JP)
   PEKGLO => NPEGLO%AL(JP)
   DO JNL=1,INL
      DO JGLO=1,KI
         IF(ZOLDVEG_GRID(JGLO)>0.0.AND.ZNEWVEG_GRID(JGLO)==0.0)THEN ! all vegetated patches disappeared = try to conserv
           S%XCCONSRV(JGLO)=S%XCCONSRV(JGLO)+PEKGLO%XLITTER(JGLO,JNL,2)*PKGLO%XPATCH(JGLO)*XGTOKG
         ENDIF
      ENDDO
   ENDDO
ENDDO
!
DO JP=1,INP
   PKGLO => NPGLO%AL(JP)
   PEKGLO => NPEGLO%AL(JP)
   DO JNC=1,INC
      DO JGLO=1,KI
         IF(ZOLDVEG_GRID(JGLO)>0.0.AND.ZNEWVEG_GRID(JGLO)==0.0)THEN ! all vegetated patches disappeared = try to conserv
           S%XCCONSRV(JGLO)=S%XCCONSRV(JGLO)+PEKGLO%XSOILCARB(JGLO,JNC)*PKGLO%XPATCH(JGLO)*XGTOKG
         ENDIF
      ENDDO
   ENDDO
ENDDO
!-----------------------------------------------------------------
!
!*      5.     Dead roots added to soil litter (gC/m2)
!              ---------------------------------------
!
DO JP=1,INP  
   !
   PK => NP%AL(JP)
   PEK => NPE%AL(JP)
   !
   CALL BIOMASS_TO_SURFACE_LITTER(PK%XTURNOVER(:,:),PEK%XLITTER(:,:,1),PEK%XLIGNIN_STRUC(:,1))
   !
   ZSOIL_LITTER      (:,:,:) = 0.0
   ZSOIL_LIGNIN_STRUC(:,:  ) = 0.0
   !
   ZSOIL_LITTER      (1:PK%NSIZE_P,1,:) = PEK%XLITTER      (:,:,2)
   ZSOIL_LIGNIN_STRUC(1:PK%NSIZE_P,1  ) = PEK%XLIGNIN_STRUC(:,  2)
   !
   CALL BIOMASS_TO_SOIL_LITTER(PK%XTURNOVER(:,:),ZSOIL_LITTER(1:PK%NSIZE_P,:,:),ZSOIL_LIGNIN_STRUC(1:PK%NSIZE_P,:))
   !
   PEK%XLITTER      (:,:,2) = ZSOIL_LITTER      (1:PK%NSIZE_P,1,:)
   PEK%XLIGNIN_STRUC(:,  2) = ZSOIL_LIGNIN_STRUC(1:PK%NSIZE_P,1  )   
   !
ENDDO
!
!-----------------------------------------------------------------
!
!*      6.     Compute current year total grid-cell carbon stocks in gC/m2
!              ---------------------------------------------------------------
!
DO JP=1,INP 
   DO JNL=1,INL
      DO JGLO=1,KI
         TLU%XLULCC_HARVEST     (JGLO,JP)=TLU%XLULCC_HARVEST     (JGLO,JP)+ZHARVEST(JGLO,JNL,JP)*XGTOKG
         TLU%XLULCC_HARVEST_GRID(JGLO   )=TLU%XLULCC_HARVEST_GRID(JGLO   )+ZHARVEST(JGLO,JNL,JP)*XGTOKG
      ENDDO
   ENDDO
ENDDO
!
DO JP=1,INP
   PK => NP%AL(JP)
   PEK => NPE%AL(JP)
   DO JNL=1,INL
      DO JI=1,PK%NSIZE_P 
         JGLO = PK%NR_P(JI)
         TLU%XLITTER_GRID_NEW(JGLO) = TLU%XLITTER_GRID_NEW(JGLO) + (PEK%XLITTER(JI,JNL,1)+PEK%XLITTER(JI,JNL,2))*PK%XPATCH(JI)
      ENDDO
   ENDDO
ENDDO
DO JP=1,INP
   PK => NP%AL(JP)
   PEK => NPE%AL(JP)
   DO JNC=1,INC
      DO JI=1,PK%NSIZE_P
           JGLO = PK%NR_P(JI)
           TLU%XCSOIL_GRID_NEW(JGLO) = TLU%XCSOIL_GRID_NEW(JGLO) + PEK%XSOILCARB(JI,JNC)*PK%XPATCH(JI)
      ENDDO
   ENDDO
ENDDO
!
IF (LHOOK) CALL DR_HOOK('LANDUSE_CARBON',1,ZHOOK_HANDLE)
!-----------------------------------------------------------------
!
END SUBROUTINE LANDUSE_CARBON
