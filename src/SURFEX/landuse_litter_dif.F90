!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!     #########
    SUBROUTINE LANDUSE_LITTER_DIF(IG, IO, S, NP, NPE, NPGLO, NPEGLO, TLU, KI, KLUOUT)
!   #####################################################################
!!****  *LAND USE LITTER*
!!
!!    PURPOSE
!!    -------
!
!     Update and conserv litter stocks after land-use change for ISBA-CC
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
!!      Original    04/2024    
!!
!-------------------------------------------------------------------------------
!
!*       0.     DECLARATIONS
!               ------------
!
USE MODD_SFX_GRID_n,     ONLY : GRID_t
USE MODD_ISBA_OPTIONS_n, ONLY : ISBA_OPTIONS_t
USE MODD_ISBA_n,         ONLY : ISBA_S_t, ISBA_P_t, ISBA_PE_t, ISBA_K_t, ISBA_NP_t, ISBA_NPE_t
USE MODD_INIT_LANDUSE,   ONLY : LULCC_t
!
USE MODD_SURF_PAR,       ONLY : XUNDEF,NUNDEF
!
USE MODD_CO2V_PAR,       ONLY : XGTOKG
!
USE MODI_BIOMASS_TO_SURFACE_LITTER
USE MODI_BIOMASS_TO_SOIL_LITTER
!
USE MODI_ABOR1_SFX
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
TYPE(GRID_t),          INTENT(INOUT) :: IG
TYPE(ISBA_OPTIONS_t),  INTENT(INOUT) :: IO
TYPE(ISBA_S_t),        INTENT(INOUT) :: S
TYPE(ISBA_NP_t),       INTENT(INOUT) :: NP, NPGLO
TYPE(ISBA_NPE_t),      INTENT(INOUT) :: NPE, NPEGLO
TYPE(LULCC_t),         INTENT(INOUT) :: TLU
!
INTEGER,               INTENT(IN)    :: KLUOUT
INTEGER,               INTENT(IN)    :: KI
!
!
!*      0.2    declarations of local parameter
!
TYPE(ISBA_P_t),  POINTER :: PK, PKGLO
TYPE(ISBA_PE_t), POINTER :: PEK, PEKGLO
!
REAL, PARAMETER                                    :: ZLIMIT = 1.0E-12
!
REAL, DIMENSION(KI)                                :: ZSURF_LIGNIN_NEW  ! current year surface lignin C stock
REAL, DIMENSION(KI)                                :: ZSURF_LIGNIN_OLD  ! previous year surface lignin C stock
!
REAL, DIMENSION(KI,IO%NNLITTER)                    :: ZSURF_LITTER_NEW  ! current year surface litter C stock
REAL, DIMENSION(KI,IO%NNLITTER)                    :: ZSURF_LITTER_OLD  ! previous year surface litter C stock
!
REAL, DIMENSION(KI,IO%NPATCH)                      :: ZP_SOIL_LITTER    ! current year litter C stock per patch
!
REAL, DIMENSION(KI,IO%NGROUND_LAYER)               :: ZSOIL_LIGNIN_NEW  ! current year lignin C stock
REAL, DIMENSION(KI,IO%NGROUND_LAYER)               :: ZSOIL_LIGNIN_OLD  ! current year lignin C stock
!
REAL, DIMENSION(KI,IO%NGROUND_LAYER,IO%NNLITTER)   :: ZSOIL_LITTER_NEW  ! current year litter C stock
REAL, DIMENSION(KI,IO%NGROUND_LAYER,IO%NNLITTER)   :: ZSOIL_LITTER_OLD  ! previous year litter C stock
!
REAL, DIMENSION(KI,IO%NNLITTER,IO%NPATCH)          :: ZHARVEST
!
REAL, DIMENSION(KI,IO%NGROUND_LAYER)               :: ZLITTER_GRID_OLD, ZLITTER_GRID_NEW
!
REAL, DIMENSION(KI)                                :: ZOLDVEG_GRID, ZBUDGET, ZHARVEST_GRID
REAL, DIMENSION(KI)                                :: ZNEWVEG_GRID, ZCONSERV, ZTURNOVER
!
INTEGER,DIMENSION(KI)                              :: INWG_LAYER_MAX
INTEGER,DIMENSION(KI,IO%NPATCH)                    :: INWG_LAYER, INWG_LAYER_VEG
!
LOGICAL :: LSTOP
!
REAL    :: ZLITTER
!
INTEGER :: INL, INP, INLIT, JLMAX ! levels
INTEGER :: JI, JL, JP, JNL, JGLO ! loop counter on levels
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!-----------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('LANDUSE_LITTER_DIF',0,ZHOOK_HANDLE)
!-----------------------------------------------------------------
!
!*      1.     Preliminaries
!              -------------
!
INL  =IO%NGROUND_LAYER
INP  =IO%NPATCH
INLIT=IO%NNLITTER
!
ZSURF_LIGNIN_OLD(:)   = 0.0
ZSOIL_LIGNIN_OLD(:,:) = 0.0
!
ZSURF_LIGNIN_NEW(:)   = 0.0
ZSOIL_LIGNIN_NEW(:,:) = 0.0
!
ZSURF_LITTER_OLD(:,:) = 0.0
ZSOIL_LITTER_OLD(:,:,:) = 0.0
!
ZSURF_LITTER_NEW(:,:) = 0.0
ZSOIL_LITTER_NEW(:,:,:) = 0.0
!
ZHARVEST  (:,:,:)=0.0
!
ZLITTER_GRID_OLD(:,:) = 0.0
ZLITTER_GRID_NEW(:,:) = 0.0
!
ZCONSERV(:) = 0.0
!
!* No litter where new patch = 0.0
!
DO JP=1,INP
   PK => NP%AL(JP)
   PEK => NPE%AL(JP)
   DO JI=1,PK%NSIZE_P 
      IF(PK%XPATCH(JI)==0.0)THEN
        PEK%XSURFACE_LIGNIN_STRUC(JI)=XUNDEF              
      ENDIF      
   ENDDO
ENDDO
!
DO JP=1,INP
   PK => NP%AL(JP)
   PEK => NPE%AL(JP)
   DO JL=1,INL
      DO JI=1,PK%NSIZE_P
         IF(PK%XPATCH(JI)==0.0)THEN
           PEK%XSOILDIF_LIGNIN_STRUC(JI,JL)=XUNDEF
         ENDIF
      ENDDO
   ENDDO
ENDDO  
!  
DO JP=1,INP
   PK => NP%AL(JP)
   PEK => NPE%AL(JP)
   DO JNL=1,INLIT
      DO JI=1,PK%NSIZE_P 
         IF(PK%XPATCH(JI)==0.0)THEN
            PEK%XSURFACE_LITTER(JI,JNL)=XUNDEF
         ENDIF
      ENDDO
   ENDDO
ENDDO  
!
DO JP=1,INP
   PK => NP%AL(JP)
   PEK => NPE%AL(JP)
   DO JNL=1,INLIT
      DO JL=1,INL
         DO JI=1,PK%NSIZE_P
            IF(PK%XPATCH(JI)==0.0)THEN
               PEK%XSOILDIF_LITTER(JI,JL,JNL)=XUNDEF
            ENDIF
         ENDDO
      ENDDO
   ENDDO
ENDDO
!
!
!-----------------------------------------------------------------
!
!*      2.     If new litter layers are very low, set to 0 and conserv
!              -------------------------------------------------------
!
!
DO JP=1,INP
   PK => NP%AL(JP)
   PEK => NPE%AL(JP)
   DO JL=1,INL
      DO JI=1,PK%NSIZE_P
         JGLO=PK%NR_P(JI)
         ZLITTER=SUM(PEK%XSOILDIF_LITTER(JI,JL,:))
         IF(PK%XPATCH(JI)>0.0.AND.JL<=PK%NWG_LAYER(JI).AND.ZLITTER>0.0.AND.ZLITTER<ZLIMIT)THEN
           ZCONSERV(JGLO)=ZCONSERV(JGLO)+ZLITTER*PK%XPATCH(JI)
           PEK%XSOILDIF_LITTER(JI,JL,:)=0.0
         ENDIF
      ENDDO
   ENDDO
ENDDO
!
!-----------------------------------------------------------------
!
!*      3.     Soil litter grids
!              -----------------
!
ZP_SOIL_LITTER(:,:) = 0.0
DO JP=1,INP
   PK => NP%AL(JP)
   PEK => NPE%AL(JP)
   DO JNL=1,INLIT
      DO JL=1,INL
         DO JI=1,PK%NSIZE_P
            JGLO=PK%NR_P(JI)
            IF(JL<=PK%NWG_LAYER(JI))THEN
               ZP_SOIL_LITTER(JGLO,JP)=ZP_SOIL_LITTER(JGLO,JP)+PEK%XSOILDIF_LITTER(JI,JL,JNL)*PK%XPATCH(JI)
            ENDIF        
         ENDDO
      ENDDO
   ENDDO
ENDDO
!
!Old litter grid
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
!New litter grid
!
INWG_LAYER(:,:) = NUNDEF
ZNEWVEG_GRID(:) = 0.0
DO JP=1,INP
   PK => NP%AL(JP)
   PEK => NPE%AL(JP)
   DO JI=1,PK%NSIZE_P
      JGLO=PK%NR_P(JI)
      INWG_LAYER(JGLO,JP)=PK%NWG_LAYER(JI)   
      IF(PEK%XLAI(JI)/=XUNDEF.AND.ZP_SOIL_LITTER(JGLO,JP)>0.0)THEN
        ZNEWVEG_GRID(JGLO)=ZNEWVEG_GRID(JGLO)+PK%XPATCH(JI)
      ENDIF       
   ENDDO
ENDDO
!
!New grid litter layers depth (accounting the case rootfrac layer depth < nwg_layer)
!
INWG_LAYER_VEG(:,:) = NUNDEF
DO JP=1,INP
   PK => NP%AL(JP)
   PEK => NPE%AL(JP)
   DO JL=1,INL
      DO JI=1,PK%NSIZE_P
         JGLO=PK%NR_P(JI)
         IF(PEK%XSOILDIF_LITTER(JI,JL,1)>0.0.AND.PEK%XSOILDIF_LITTER(JI,JL,1)<XUNDEF.AND.PEK%XWG(JI,JL)/=XUNDEF)THEN
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
!* To compute litter budget
!
DO JP=1,INP
   PKGLO => NPGLO%AL(JP)
   PEKGLO => NPEGLO%AL(JP)
   DO JNL=1,INLIT
      DO JGLO=1,KI
         TLU%XLITTER_GRID_OLD(JGLO) = TLU%XLITTER_GRID_OLD(JGLO) + PEKGLO%XSURFACE_LITTER(JGLO,JNL)*PKGLO%XPATCH(JGLO)
      ENDDO
   ENDDO
ENDDO
!
DO JP=1,INP
   PKGLO => NPGLO%AL(JP)
   PEKGLO => NPEGLO%AL(JP)
   DO JNL=1,INLIT
      DO JL=1,INL
         DO JGLO=1,KI     
            IF(PKGLO%XPATCH(JGLO)>0.0.AND.JL<=PKGLO%NWG_LAYER(JGLO))THEN      
              ZLITTER_GRID_OLD(JGLO,JL) = ZLITTER_GRID_OLD(JGLO,JL) + PEKGLO%XSOILDIF_LITTER(JGLO,JL,JNL)*PKGLO%XPATCH(JGLO)
            ENDIF
         ENDDO
      ENDDO
   ENDDO
ENDDO
!
!-----------------------------------------------------------------
!
!*      4.     Patch disappear = surface litter is harvested
!              ---------------------------------------------
!
DO JP=1,INP
   PKGLO => NPGLO%AL(JP)
   PEKGLO => NPEGLO%AL(JP)
   DO JNL=1,INLIT
      DO JGLO=1,KI      
         IF(PKGLO%XPATCH(JGLO)>0.0.AND.S%XPATCH(JGLO,JP)==0.0)THEN
            ZHARVEST(JGLO,JNL,JP)=ZHARVEST(JGLO,JNL,JP)+PEKGLO%XSURFACE_LITTER(JGLO,JNL)*PKGLO%XPATCH(JGLO)
         ENDIF         
      ENDDO
   ENDDO
ENDDO  
!
!-----------------------------------------------------------------
!
!*      5.     Compute previous and current year total grid-cell biomass and litter stocks
!              ---------------------------------------------------------------------------
!
!Surface lignin stock in gC/m2
!
DO JP=1,INP
   PKGLO => NPGLO%AL(JP)
   PEKGLO => NPEGLO%AL(JP)
   DO JGLO=1,KI
      ZSURF_LIGNIN_OLD(JGLO) = ZSURF_LIGNIN_OLD(JGLO) + PEKGLO%XSURFACE_LIGNIN_STRUC(JGLO)*PKGLO%XPATCH(JGLO)
   ENDDO
ENDDO
!
DO JP=1,INP
   PK => NP%AL(JP)
   PEK => NPE%AL(JP)
   DO JI=1,PK%NSIZE_P 
      JGLO=PK%NR_P(JI)
      ZSURF_LIGNIN_NEW(JGLO) = ZSURF_LIGNIN_NEW(JGLO) + PEK%XSURFACE_LIGNIN_STRUC(JI)*PK%XPATCH(JI)
   ENDDO
ENDDO
!
!Soil lignin stock in gC/m2
!
DO JP=1,INP
   PKGLO => NPGLO%AL(JP)
   PEKGLO => NPEGLO%AL(JP)
   DO JL=1,INL
      DO JGLO=1,KI
         IF(PKGLO%XPATCH(JGLO)>0.0.AND.JL<=PKGLO%NWG_LAYER(JGLO))THEN      
           JLMAX=MIN(JL,INWG_LAYER_MAX(JGLO))
           ZSOIL_LIGNIN_OLD(JGLO,JLMAX) = ZSOIL_LIGNIN_OLD(JGLO,JLMAX) + PEKGLO%XSOILDIF_LIGNIN_STRUC(JGLO,JL)*PKGLO%XPATCH(JGLO)  
         ENDIF
      ENDDO
   ENDDO
ENDDO  
!
DO JP=1,INP
   PK => NP%AL(JP)
   PEK => NPE%AL(JP)
   DO JL=1,INL
      DO JI=1,PK%NSIZE_P   
         JGLO=PK%NR_P(JI)
         IF(PK%XPATCH(JI)>0.0.AND.JL<=PK%NWG_LAYER(JI))THEN      
           ZSOIL_LIGNIN_NEW(JGLO,JL) = ZSOIL_LIGNIN_NEW(JGLO,JL) + PEK%XSOILDIF_LIGNIN_STRUC(JI,JL)*PK%XPATCH(JI)
         ENDIF
      ENDDO
   ENDDO
ENDDO  
!
!Surface litter stock in gC/m2
!
DO JP=1,INP
   PKGLO => NPGLO%AL(JP)
   PEKGLO => NPEGLO%AL(JP)
   DO JNL=1,INLIT
      DO JGLO=1,KI 
         IF(PKGLO%XPATCH(JGLO)>0.0)THEN
            ZSURF_LITTER_OLD(JGLO,JNL) = ZSURF_LITTER_OLD(JGLO,JNL) + PEKGLO%XSURFACE_LITTER(JGLO,JNL)*PKGLO%XPATCH(JGLO) &
                                                                    - ZHARVEST(JGLO,JNL,JP)
         ENDIF
      ENDDO
   ENDDO
ENDDO
!
DO JP=1,INP
   PK => NP%AL(JP)
   PEK => NPE%AL(JP)
   DO JNL=1,INLIT
      DO JI=1,PK%NSIZE_P 
         JGLO=PK%NR_P(JI)
         IF(PK%XPATCH(JI)>0.0)THEN
            ZSURF_LITTER_NEW(JGLO,JNL) = ZSURF_LITTER_NEW(JGLO,JNL) + PEK%XSURFACE_LITTER(JI,JNL)*PK%XPATCH(JI)
         ENDIF
      ENDDO
   ENDDO
ENDDO
!
!Soil litter stock in gC/m2
!
DO JP=1,INP
   PKGLO => NPGLO%AL(JP)
   PEKGLO => NPEGLO%AL(JP)
   DO JNL=1,INLIT
      DO JL=1,INL
         DO JGLO=1,KI        
            IF(PKGLO%XPATCH(JGLO)>0.0.AND.JL<=PKGLO%NWG_LAYER(JGLO))THEN
              JLMAX=MIN(JL,INWG_LAYER_MAX(JGLO))
              ZSOIL_LITTER_OLD(JGLO,JLMAX,JNL) = ZSOIL_LITTER_OLD(JGLO,JLMAX,JNL) + PEKGLO%XSOILDIF_LITTER(JGLO,JL,JNL) &
                                                                                  * PKGLO%XPATCH(JGLO)
            ENDIF              
         ENDDO
      ENDDO
   ENDDO
ENDDO
!
DO JP=1,INP
   PK => NP%AL(JP)
   PEK => NPE%AL(JP)
   DO JNL=1,INLIT
      DO JL=1,INL
         DO JI=1,PK%NSIZE_P        
            JGLO=PK%NR_P(JI)
            IF(PK%XPATCH(JI)>0.0.AND.JL<=PK%NWG_LAYER(JI))THEN      
              ZSOIL_LITTER_NEW(JGLO,JL,JNL) = ZSOIL_LITTER_NEW(JGLO,JL,JNL) + PEK%XSOILDIF_LITTER(JI,JL,JNL)*PK%XPATCH(JI)     
            ENDIF             
         ENDDO
      ENDDO
   ENDDO
ENDDO
!
!-----------------------------------------------------------------
!
!*      6.      Ensures grid-cell conservation for litter stock
!              ------------------------------------------------
!
!
DO JP=1,INP
   PK => NP%AL(JP)
   PEK => NPE%AL(JP)
   DO JI=1,PK%NSIZE_P
      JGLO=PK%NR_P(JI)
      IF(PK%XPATCH(JI)>0.0.AND.ZSURF_LIGNIN_NEW(JGLO)/=ZSURF_LIGNIN_OLD(JGLO).AND.ZSURF_LIGNIN_NEW(JGLO)/=0.0)THEN
        PEK%XSURFACE_LIGNIN_STRUC(JI) = PEK%XSURFACE_LIGNIN_STRUC(JI) * ZSURF_LIGNIN_OLD(JGLO)/ZSURF_LIGNIN_NEW(JGLO)
      ENDIF
   ENDDO
ENDDO
!
DO JP=1,INP
   PK => NP%AL(JP)
   PEK => NPE%AL(JP)
   DO JL=1,INL
      DO JI=1,PK%NSIZE_P
         JGLO=PK%NR_P(JI)
         IF(JL<=PK%NWG_LAYER(JI).AND.ZSOIL_LIGNIN_NEW(JGLO,JL)/=ZSOIL_LIGNIN_OLD(JGLO,JL).AND.ZSOIL_LIGNIN_NEW(JGLO,JL)/=0.0)THEN
           PEK%XSOILDIF_LIGNIN_STRUC(JI,JL) = PEK%XSOILDIF_LIGNIN_STRUC(JI,JL) * ZSOIL_LIGNIN_OLD(JGLO,JL) &
                                                                               / ZSOIL_LIGNIN_NEW(JGLO,JL)
         ENDIF
      ENDDO
   ENDDO
ENDDO
!
DO JP=1,INP  
   PK => NP%AL(JP)
   PEK => NPE%AL(JP)
   DO JNL=1,INLIT
      DO JI=1,PK%NSIZE_P
        JGLO=PK%NR_P(JI)
        IF(PK%XPATCH(JI)>0.0.AND.ZSURF_LITTER_NEW(JGLO,JNL)/=ZSURF_LITTER_OLD(JGLO,JNL) &
                            .AND.ZSURF_LITTER_NEW(JGLO,JNL)/=0.0                        )THEN
           PEK%XSURFACE_LITTER(JI,JNL) = PEK%XSURFACE_LITTER(JI,JNL) * ZSURF_LITTER_OLD(JGLO,JNL)/ZSURF_LITTER_NEW(JGLO,JNL)
        ENDIF
      ENDDO
   ENDDO
ENDDO
!
DO JP=1,INP
   PK => NP%AL(JP)
   PEK => NPE%AL(JP)
   DO JNL=1,INLIT
      DO JL=1,INL
         DO JI=1,PK%NSIZE_P
            JGLO=PK%NR_P(JI)
            IF(PK%XPATCH(JI)>0.0.AND.JL<=PK%NWG_LAYER(JI).AND.ZSOIL_LITTER_NEW(JGLO,JL,JNL)/=ZSOIL_LITTER_OLD(JGLO,JL,JNL) &
                                                         .AND.ZSOIL_LITTER_NEW(JGLO,JL,JNL)/=0.0                           )THEN
              PEK%XSOILDIF_LITTER(JI,JL,JNL) = PEK%XSOILDIF_LITTER(JI,JL,JNL) * ZSOIL_LITTER_OLD(JGLO,JL,JNL) &
                                                                              / ZSOIL_LITTER_NEW(JGLO,JL,JNL)
            ENDIF
         ENDDO
      ENDDO
   ENDDO
ENDDO
!
!-----------------------------------------------------------------
!
!*      7.     If all vegetated patches disappeared -> try to conserv
!              ------------------------------------------------------
!
DO JP=1,INP
   PKGLO => NPGLO%AL(JP)
   PEKGLO => NPEGLO%AL(JP)
   DO JNL=1,INLIT
      DO JL=1,INL
         DO JGLO=1,KI
            IF(PKGLO%XPATCH(JGLO)>0.0.AND.JL<=PKGLO%NWG_LAYER(JGLO).AND.ZOLDVEG_GRID(JGLO)>0.0.AND.ZNEWVEG_GRID(JGLO)==0.0)THEN
              ZCONSERV(JGLO)=ZCONSERV(JGLO)+PEKGLO%XSOILDIF_LITTER(JGLO,JL,JNL)*PKGLO%XPATCH(JGLO)
            ENDIF
         ENDDO
      ENDDO
   ENDDO
ENDDO
!
!-----------------------------------------------------------------
!
!*      8.     Dead roots added to soil litter (gC/m2)
!              ---------------------------------------
!
DO JP=1,INP
   PK => NP%AL(JP)
   PEK => NPE%AL(JP)
   CALL BIOMASS_TO_SURFACE_LITTER(PK%XTURNOVER,PEK%XSURFACE_LITTER,PEK%XSURFACE_LIGNIN_STRUC)
   CALL BIOMASS_TO_SOIL_LITTER   (PK%XTURNOVER,PEK%XSOILDIF_LITTER,PEK%XSOILDIF_LIGNIN_STRUC,&
                                  KWG_LAYER=PK%NWG_LAYER,PROOTFRAC=PK%XROOTFRAC)  
ENDDO
!
!-----------------------------------------------------------------
!
!*      9.     Compute current year surface litter stocks in gC/m2
!              ---------------------------------------------------
!
!
ZHARVEST_GRID(:) = 0.0
DO JP=1,INP  
   DO JNL=1,INLIT
      DO JGLO=1,KI
         ZHARVEST_GRID     (JGLO   )=ZHARVEST_GRID     (JGLO   )+ZHARVEST(JGLO,JNL,JP)
         TLU%XLULCC_HARVEST(JGLO,JP)=TLU%XLULCC_HARVEST(JGLO,JP)+ZHARVEST(JGLO,JNL,JP)*XGTOKG
      ENDDO
   ENDDO
ENDDO
!  
TLU%XLULCC_HARVEST_GRID(:)=TLU%XLULCC_HARVEST_GRID(:)+ZHARVEST_GRID(:)*XGTOKG
!
DO JP=1,INP
   PK => NP%AL(JP)
   PEK => NPE%AL(JP)
   DO JNL=1,INLIT
      DO JI=1,PK%NSIZE_P 
         JGLO=PK%NR_P(JI)
         IF(PK%XPATCH(JI)>0.0)THEN
           TLU%XLITTER_GRID_NEW(JGLO) = TLU%XLITTER_GRID_NEW(JGLO) + PEK%XSURFACE_LITTER(JI,JNL)*PK%XPATCH(JI)
         ENDIF
      ENDDO
   ENDDO
ENDDO
!
ZTURNOVER(:) = 0.0
DO JP=1,INP  
   PK => NP%AL(JP)
   DO JI=1,PK%NSIZE_P
      JGLO=PK%NR_P(JI)
      ZTURNOVER(JGLO) = ZTURNOVER(JGLO) + (PK%XTURNOVER(JI,1)+PK%XTURNOVER(JI,2))*PK%XPATCH(JI)
   ENDDO
ENDDO
!
ZBUDGET(:) = TLU%XLITTER_GRID_NEW(:)-TLU%XLITTER_GRID_OLD(:)+ZHARVEST_GRID(:)-ZTURNOVER(:)
!
LSTOP=.FALSE.
DO JGLO=1,KI
   IF(ABS(ZBUDGET(JGLO))>1.0E-10)THEN
     WRITE(KLUOUT,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
     WRITE(KLUOUT,*)'LANDUSE_LITTER_DIF: NO SURFACE CONSERVATION IN AT LEAST ONE GRID CELL'
     WRITE(KLUOUT,*)'LON = ',IG%XLON(JGLO),' LAT =',IG%XLAT(JGLO),'BUDGET =',ZBUDGET(JGLO),'g/m2'
     WRITE(KLUOUT,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
     LSTOP=.TRUE.
   ENDIF           
ENDDO
IF(LSTOP) CALL ABOR1_SFX('LANDUSE_LITTER_DIF: INCONSISTENCY IN SURFACE LITTER BUDGET')
!
!
!-----------------------------------------------------------------
!
!*      10.     Compute current year total grid-cell litter stocks in gC/m2
!              ---------------------------------------------------------------
!
S%XCCONSRV(:) = S%XCCONSRV(:) + ZCONSERV(:)*XGTOKG
!
DO JP=1,INP
   PK => NP%AL(JP)
   PEK => NPE%AL(JP)
   DO JNL=1,INLIT
      DO JL=1,INL
         DO JI=1,PK%NSIZE_P
            JGLO=PK%NR_P(JI)
            IF(PK%XPATCH(JI)>0.0.AND.JL<=PK%NWG_LAYER(JI))THEN      
              ZLITTER_GRID_NEW(JGLO,JL) = ZLITTER_GRID_NEW(JGLO,JL) + PEK%XSOILDIF_LITTER(JI,JL,JNL)*PK%XPATCH(JI)
            ENDIF
         ENDDO
      ENDDO
   ENDDO
ENDDO
!
ZBUDGET(:) = 0.0
DO JL=1,INL
   DO JGLO=1,KI
      !
      TLU%XLITTER_GRID_OLD(JGLO) = TLU%XLITTER_GRID_OLD(JGLO) + ZLITTER_GRID_OLD(JGLO,JL)
      TLU%XLITTER_GRID_NEW(JGLO) = TLU%XLITTER_GRID_NEW(JGLO) + ZLITTER_GRID_NEW(JGLO,JL)
      !
      ZBUDGET(JGLO) = ZBUDGET(JGLO) + (ZLITTER_GRID_NEW(JGLO,JL)-ZLITTER_GRID_OLD(JGLO,JL))
      !
   ENDDO
ENDDO
!
ZTURNOVER(:) = 0.0
DO JP=1,INP  
   PK => NP%AL(JP)
   DO JI=1,PK%NSIZE_P
      JGLO=PK%NR_P(JI)
      ZTURNOVER(JGLO) = ZTURNOVER(JGLO) + (PK%XTURNOVER(JI,4)+PK%XTURNOVER(JI,6))*PK%XPATCH(JI)
   ENDDO
ENDDO
!
ZBUDGET(:) = ZBUDGET(:) + ZCONSERV(:) - ZTURNOVER(:)
!
LSTOP=.FALSE.
DO JGLO=1,KI
   IF(ABS(ZBUDGET(JGLO))>1.0E-10)THEN
     WRITE(KLUOUT,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
     WRITE(KLUOUT,*)'LANDUSE_LITTER_DIF: NO SOIL LITTER CONSERVATION IN AT LEAST ONE GRID CELL'
     WRITE(KLUOUT,*)'LON = ',IG%XLON(JGLO),' LAT =',IG%XLAT(JGLO),'BUDGET =',ZBUDGET(JGLO),'g/m2'
     WRITE(KLUOUT,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
     LSTOP=.TRUE.
   ENDIF           
ENDDO
IF(LSTOP) CALL ABOR1_SFX('LANDUSE_LITTER_DIF: INCONSISTENCY IN SOIL LITTER BUDGET')
!
IF (LHOOK) CALL DR_HOOK('LANDUSE_LITTER_DIF',1,ZHOOK_HANDLE)
!-----------------------------------------------------------------
!
END SUBROUTINE LANDUSE_LITTER_DIF
