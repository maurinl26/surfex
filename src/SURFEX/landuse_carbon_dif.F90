!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!     #########
    SUBROUTINE LANDUSE_CARBON_DIF(IO, S, NP, NPE, NPGLO, NPEGLO, TLU, KI)
!   #####################################################################
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
USE MODD_SURF_PAR,       ONLY : XUNDEF,NUNDEF
!
USE MODD_CO2V_PAR,       ONLY : XGTOKG
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
!
!*      0.2    declarations of local parameter
!
TYPE(ISBA_P_t),  POINTER :: PK, PKGLO
TYPE(ISBA_PE_t), POINTER :: PEK, PEKGLO
!
REAL, PARAMETER                                    :: ZLIMIT = 1.E-12
!
REAL, DIMENSION(KI)                                :: ZOLDVEG_GRID
REAL, DIMENSION(KI)                                :: ZNEWVEG_GRID
!
REAL, DIMENSION(KI,IO%NPATCH)                      :: ZP_CSOIL    ! current year litter C stock per patch
!
REAL, DIMENSION(KI,IO%NGROUND_LAYER,IO%NNSOILCARB) :: ZCSOIL_NEW  ! current  year soil C stock
REAL, DIMENSION(KI,IO%NGROUND_LAYER,IO%NNSOILCARB) :: ZCSOIL_OLD  ! previous year soil C stock
!
REAL, DIMENSION(KI,IO%NGROUND_LAYER)               :: ZCSOIL_GRID_OLD, ZCSOIL_GRID_NEW
!
INTEGER,DIMENSION(KI)                              :: INWG_LAYER_MAX
INTEGER,DIMENSION(KI,IO%NPATCH)                    :: INWG_LAYER, INWG_LAYER_VEG
!
REAL    :: ZSCARB
!
INTEGER :: INL, INP, INC, JLMAX  ! levels
INTEGER :: JI, JL, JP, JNC, JGLO ! loop counter on levels
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!-----------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('LANDUSE_CARBON_DIF',0,ZHOOK_HANDLE)
!-----------------------------------------------------------------
!
!*      1.     Preliminaries
!              -------------
!
INL  =IO%NGROUND_LAYER
INP  =IO%NPATCH
INC  =IO%NNSOILCARB
!
ZCSOIL_OLD(:,:,:) = 0.0
ZCSOIL_NEW(:,:,:) = 0.0
!
ZCSOIL_GRID_OLD (:,:) = 0.0
ZCSOIL_GRID_NEW (:,:) = 0.0
!
!* No carbon where new patch = 0.0
!
!
DO JP=1,INP
   PK => NP%AL(JP)
   PEK => NPE%AL(JP)
   DO JNC=1,INC
      DO JL=1,INL
         DO JI=1,PK%NSIZE_P
            IF(PK%XPATCH(JI)==0.0)THEN
               PEK%XSOILDIF_CARB(JI,JL,JNC)=XUNDEF
            ENDIF
         ENDDO
      ENDDO
   ENDDO
ENDDO
!
!-----------------------------------------------------------------
!
!*      2.     If new carbon layers are very low, set to 0 and conserv
!              -------------------------------------------------------
!
!
DO JP=1,INP
   PK => NP%AL(JP)
   PEK => NPE%AL(JP)
   DO JL=1,INL
      DO JI=1,PK%NSIZE_P
         JGLO=PK%NR_P(JI)
         ZSCARB=SUM(PEK%XSOILDIF_CARB(JI,JL,:))
         IF(PK%XPATCH(JI)>0.0.AND.JL<=PK%NWG_LAYER(JI).AND.ZSCARB>0.0.AND.ZSCARB<ZLIMIT)THEN
           S%XCCONSRV(JGLO)=S%XCCONSRV(JGLO)+ZSCARB*PK%XPATCH(JI)*XGTOKG
           PEK%XSOILDIF_CARB(JI,JL,:)=0.0
         ENDIF
      ENDDO
   ENDDO
ENDDO
!
!-----------------------------------------------------------------
!
!*      3.     Soil carbon grids
!              -----------------
!
ZP_CSOIL(:,:) = 0.0
DO JP=1,INP
   PK => NP%AL(JP)
   PEK => NPE%AL(JP)
   DO JNC=1,INC
      DO JL=1,INL
         DO JI=1,PK%NSIZE_P
            JGLO=PK%NR_P(JI)
            IF(JL<=PK%NWG_LAYER(JI))THEN
               ZP_CSOIL(JGLO,JP)=ZP_CSOIL(JGLO,JP)+PEK%XSOILDIF_CARB(JI,JL,JNC)*PK%XPATCH(JI)
            ENDIF        
         ENDDO
      ENDDO
   ENDDO
ENDDO
!
!Old carbon grid
!
ZOLDVEG_GRID(:) = 0.0
DO JP=1,INP
   PKGLO => NPGLO%AL(JP)
   PEKGLO => NPEGLO%AL(JP)
   DO JGLO=1,KI
      IF(PEKGLO%XLAI(JGLO)/=XUNDEF)THEN
        ZOLDVEG_GRID(JGLO)=ZOLDVEG_GRID(JGLO)+PKGLO%XPATCH(JGLO)
      ENDIF 
   ENDDO
ENDDO
!
!New carbon grid
!
INWG_LAYER(:,:) = NUNDEF
ZNEWVEG_GRID(:) = 0.0
DO JP=1,INP
   PK => NP%AL(JP)
   PEK => NPE%AL(JP)
   DO JI=1,PK%NSIZE_P
      JGLO=PK%NR_P(JI)
      INWG_LAYER(JGLO,JP)=PK%NWG_LAYER(JI)      
      IF(PEK%XLAI(JI)/=XUNDEF.AND.ZP_CSOIL(JGLO,JP)>0.0)THEN
        ZNEWVEG_GRID(JGLO)=ZNEWVEG_GRID(JGLO)+PK%XPATCH(JI)
      ENDIF       
   ENDDO
ENDDO
!
!New grid carbon layers depth (accounting the case rootfrac layer depth < nwg_layer)
!
INWG_LAYER_VEG(:,:) = NUNDEF
DO JP=1,INP
   PK => NP%AL(JP)
   PEK => NPE%AL(JP)
   DO JL=1,INL
      DO JI=1,PK%NSIZE_P
         JGLO=PK%NR_P(JI)
         IF(PEK%XSOILDIF_CARB(JI,JL,1)>0.0.AND.PEK%XSOILDIF_CARB(JI,JL,1)<XUNDEF.AND.PEK%XWG(JI,JL)/=XUNDEF)THEN
           INWG_LAYER_VEG(JGLO,JP)=JL
         ENDIF       
      ENDDO
   ENDDO
ENDDO
!
DO JGLO=1,KI
   IF(ALL(INWG_LAYER_VEG(JGLO,:)==NUNDEF))THEN
     INWG_LAYER_MAX(JGLO)=MAXVAL(INWG_LAYER(JGLO,:),INWG_LAYER(JGLO,:)/=NUNDEF)
   ELSE
     INWG_LAYER_MAX(JGLO)=MAXVAL(INWG_LAYER_VEG(JGLO,:),INWG_LAYER_VEG(JGLO,:)/=NUNDEF)
   ENDIF
ENDDO
!
!* To compute carbon budget
!
DO JP=1,INP
   PKGLO => NPGLO%AL(JP)
   PEKGLO => NPEGLO%AL(JP)
   DO JNC=1,INC
      DO JL=1,INL
         DO JGLO=1,KI     
            IF(PKGLO%XPATCH(JGLO)>0.0.AND.JL<=PKGLO%NWG_LAYER(JGLO))THEN
              ZCSOIL_GRID_OLD(JGLO,JL) = ZCSOIL_GRID_OLD(JGLO,JL) + PEKGLO%XSOILDIF_CARB(JGLO,JL,JNC)*PKGLO%XPATCH(JGLO)
            ENDIF
         ENDDO
      ENDDO
   ENDDO
ENDDO
!
!-----------------------------------------------------------------
!
!*      4.     Compute previous and current year total grid-cell carbon stocks
!              ---------------------------------------------------------------
!
!
!Soil carbon stock in gC/m2
!
DO JP=1,INP
   PKGLO => NPGLO%AL(JP)
   PEKGLO => NPEGLO%AL(JP)
   DO JNC=1,INC
      DO JL=1,INL
         DO JGLO=1,KI
            IF(PKGLO%XPATCH(JGLO)>0.0.AND.JL<=PKGLO%NWG_LAYER(JGLO))THEN  
              JLMAX=MIN(JL,INWG_LAYER_MAX(JGLO))
              ZCSOIL_OLD(JGLO,JLMAX,JNC) = ZCSOIL_OLD(JGLO,JLMAX,JNC) + PEKGLO%XSOILDIF_CARB(JGLO,JL,JNC)*PKGLO%XPATCH(JGLO)
            ENDIF 
         ENDDO
      ENDDO
   ENDDO
ENDDO
!
DO JP=1,INP
   PK => NP%AL(JP)
   PEK => NPE%AL(JP)
   DO JNC=1,INC
      DO JL=1,INL
         DO JI=1,PK%NSIZE_P
            JGLO=PK%NR_P(JI)
            IF(PK%XPATCH(JI)>0.0.AND.JL<=PK%NWG_LAYER(JI))THEN      
              ZCSOIL_NEW(JGLO,JL,JNC) = ZCSOIL_NEW(JGLO,JL,JNC) + PEK%XSOILDIF_CARB(JI,JL,JNC)*PK%XPATCH(JI)        
            ENDIF
         ENDDO
      ENDDO
   ENDDO
ENDDO
!
!-----------------------------------------------------------------
!
!*      5.      Ensures grid-cell conservation for carbon stock
!              ------------------------------------------------
!
!
DO JP=1,INP
   PK => NP%AL(JP)
   PEK => NPE%AL(JP)
   DO JNC=1,INC
      DO JL=1,INL
         DO JI=1,PK%NSIZE_P
            JGLO=PK%NR_P(JI)
            IF(PK%XPATCH(JI)>0.0.AND.JL<=PK%NWG_LAYER(JI).AND.ZCSOIL_NEW(JGLO,JL,JNC)/=ZCSOIL_OLD(JGLO,JL,JNC) &
                                                         .AND.ZCSOIL_NEW(JGLO,JL,JNC)/=0.0)THEN
              PEK%XSOILDIF_CARB(JI,JL,JNC) = PEK%XSOILDIF_CARB(JI,JL,JNC) * ZCSOIL_OLD(JGLO,JL,JNC) / ZCSOIL_NEW(JGLO,JL,JNC)
            ENDIF
         ENDDO
      ENDDO
   ENDDO
ENDDO
!
!-----------------------------------------------------------------
!
!*      6.     If all vegetated patches disappeared -> try to conserv
!              ------------------------------------------------------
!
DO JP=1,INP
   PKGLO => NPGLO%AL(JP)
   PEKGLO => NPEGLO%AL(JP)
   DO JNC=1,INC
      DO JL=1,INL
         DO JGLO=1,KI
            IF(PKGLO%XPATCH(JGLO)>0.0.AND.JL<=PKGLO%NWG_LAYER(JGLO).AND.ZOLDVEG_GRID(JGLO)>0.0.AND.ZNEWVEG_GRID(JGLO)==0.0)THEN
              S%XCCONSRV(JGLO)=S%XCCONSRV(JGLO)+PEKGLO%XSOILDIF_CARB(JGLO,JL,JNC)*PKGLO%XPATCH(JGLO)*XGTOKG
            ENDIF
         ENDDO
      ENDDO
   ENDDO
ENDDO
!
!-----------------------------------------------------------------
!
!*      7.     Compute current year total grid-cell carbon stocks in gC/m2
!              ---------------------------------------------------------------
!
!
DO JP=1,INP
   PK => NP%AL(JP)
   PEK => NPE%AL(JP)
   DO JNC=1,INC
      DO JL=1,INL
         DO JI=1,PK%NSIZE_P
            JGLO=PK%NR_P(JI)
            IF(PK%XPATCH(JI)>0.0.AND.JL<=PK%NWG_LAYER(JI))THEN
              ZCSOIL_GRID_NEW(JGLO,JL) = ZCSOIL_GRID_NEW(JGLO,JL) + PEK%XSOILDIF_CARB(JI,JL,JNC)*PK%XPATCH(JI)
            ENDIF
         ENDDO
      ENDDO
   ENDDO
ENDDO
!
DO JL=1,INL
   DO JGLO=1,KI
      TLU%XCSOIL_GRID_OLD (JGLO) = TLU%XCSOIL_GRID_OLD (JGLO) + ZCSOIL_GRID_OLD (JGLO,JL)
      TLU%XCSOIL_GRID_NEW (JGLO) = TLU%XCSOIL_GRID_NEW (JGLO) + ZCSOIL_GRID_NEW (JGLO,JL)
   ENDDO
ENDDO
!
IF (LHOOK) CALL DR_HOOK('LANDUSE_CARBON_DIF',1,ZHOOK_HANDLE)
!-----------------------------------------------------------------
!
END SUBROUTINE LANDUSE_CARBON_DIF
