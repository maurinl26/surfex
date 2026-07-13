!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!#############################################################
SUBROUTINE LANDUSE_NEWPATCH (IO, K, NK, NP, NPE, SOLD, NPGLO, NPEGLO, KI)  
!#############################################################
!
!!****  *LANDUSE_NEWPATCH* - The algorithm ATTRIBUTE_CLOSEST_VEGTYPE is simple:
!!                           It uses attribution rules for each kind of vegetation
!!                           to find the closest patch matching when a new patch 
!!                           is created in the grid cell
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
USE MODD_SFX_GRID_n,     ONLY : GRID_t
USE MODD_ISBA_OPTIONS_n, ONLY : ISBA_OPTIONS_t
USE MODD_ISBA_n,         ONLY : ISBA_S_t, ISBA_P_t, ISBA_PE_t, ISBA_K_t, ISBA_NK_t, &
                                ISBA_NP_t, ISBA_NPE_t, ISBA_S_INIT, ISBA_NP_INIT, ISBA_NPE_INIT
!
USE MODD_INIT_LANDUSE
!
USE MODD_SURF_PAR,       ONLY : XUNDEF, NUNDEF
!
USE MODD_CSTS,           ONLY : XDAY
!
USE MODD_DATA_COVER_PAR, ONLY : NVEGTYPE
!
USE MODD_CO2V_PAR,       ONLY : XANFMINIT,XCA_NIT,XCC_NIT
!
USE MODI_ATTRIBUTE_CLOSEST_VEGTYPE
!
USE YOMHOOK,  ONLY : LHOOK,   DR_HOOK
USE PARKIND1, ONLY : JPRB
!
IMPLICIT NONE
!
!*       0.1   Declarations of arguments
!              -------------------------
!
TYPE(ISBA_OPTIONS_t),  INTENT(INOUT) :: IO
TYPE(ISBA_S_t),        INTENT(INOUT) :: SOLD
TYPE(ISBA_K_t),        INTENT(INOUT) :: K
TYPE(ISBA_NK_t),       INTENT(INOUT) :: NK
TYPE(ISBA_NP_t),       INTENT(INOUT) :: NP, NPGLO
TYPE(ISBA_NPE_t),      INTENT(INOUT) :: NPE, NPEGLO
!
INTEGER,               INTENT(IN)    :: KI
!
!
!*       0.2   Declarations of local arguments on Patch grid
!
TYPE(LULCC_t) :: TLU
!
TYPE(ISBA_K_t),  POINTER :: KK
TYPE(ISBA_P_t),  POINTER :: PK, PKGLO
TYPE(ISBA_PE_t), POINTER :: PEK, PEKGLO, PNEAR
!
!*       0.3   Declarations of local arguments on complete ISBA grid
!
INTEGER :: JI, JL, JNL, JNLS, JNC, JP, JP_NEAR ! loop counter
INTEGER :: INL, INS, INB, INLIT, INLITS, INC   ! dimension
INTEGER :: IGLO, ISIZE_LMEB_PATCH              ! Work integer
REAL    :: ZCC_CA, ZINVCA                      ! Work real
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('LANDUSE_NEWPATCH',0,ZHOOK_HANDLE)
!-------------------------------------------------------------------------------
!        
!
!*       0. Global initialisation
!        ------------------------
!
INL   =IO%NGROUND_LAYER
INS   =NPE%AL(1)%TSNOW%NLAYER
INB   =IO%NNBIOMASS
INLIT =IO%NNLITTER
INLITS=IO%NNLITTLEVS
INC   =IO%NNSOILCARB
!
ISIZE_LMEB_PATCH=COUNT(IO%LMEB_PATCH(:))
!
IF(IO%CPHOTO=='NIT'.OR.IO%CPHOTO=='NCB')THEN   
  ZCC_CA=(XCC_NIT/EXP(XCA_NIT*LOG(10.)))  
  ZINVCA=(1.0/(1.0-XCA_NIT))
ENDIF
!
DO JP=1,IO%NPATCH
   !
   KK => NK%AL(JP)
   PK => NP%AL(JP)
   PEK => NPE%AL(JP)
   PKGLO => NPGLO%AL(JP)
   PEKGLO => NPEGLO%AL(JP)
   !
   DO JI=1,PK%NSIZE_P
      !
      !* Ensures that solely created patches are treated
      !     
      IGLO = PK%NR_P(JI)
      !
      IF((PK%XPATCH(JI)>0.).AND.(PKGLO%XPATCH(IGLO)==0.0))THEN
        !
        !* New patch appears
        !
        CALL ATTRIBUTE_CLOSEST_VEGTYPE(IO%NPATCH,NVEGTYPE,SOLD%XPATCH(IGLO,:),JP,JP_NEAR)
        !
        PNEAR => NPEGLO%AL(JP_NEAR)
        !
        !* Soil temperature
        !
        PEK%XTG(JI,:) = PNEAR%XTG(IGLO,:)
        !
        !* soil liquid and ice water contents
        !
        PEK%XWG (JI,:) = PNEAR%XWG (IGLO,:)
        PEK%XWGI(JI,:) = PNEAR%XWGI(IGLO,:)
        !
        IF(IO%CISBA/='DIF')THEN
          PEK%XWGI(JI,3) = 0.0 !Only 2 soil ice layers in ISBA-FR
        ENDIF
        !
        DO JL=1,INL
           IF(PEK%XWG(JI,JL)/=XUNDEF.AND.(PEK%XWG(JI,JL)+PEK%XWGI(JI,JL))>KK%XWSAT(JI,JL))THEN
              PEK%XWGI(JI,JL) = KK%XWSAT(JI,JL) - PEK%XWG(JI,JL)
           ENDIF
        ENDDO
        !
        !* water intercepted on leaves
        !
        PEK%XWR(JI) = 0.0
        !
        !* Canopy air specific humidity
        !
        PEK%XQC(JI) = PNEAR%XQC(IGLO)
        ! 
        !* glacier ice storage (semi-pro)
        !
        IF(IO%LGLACIER)THEN
          PEK%XICE_STO(JI) = 0.0
        ENDIF
        !
        !* snow Albedo
        !
        PEK%TSNOW%ALB(JI) = PNEAR%TSNOW%ALB(IGLO)
        !
        !* snow water equivalent and density
        !
        PEK%TSNOW%WSNOW(JI,:) = PNEAR%TSNOW%WSNOW(IGLO,:)
        PEK%TSNOW%RHO  (JI,:) = PNEAR%TSNOW%RHO  (IGLO,:)
        !
        !* Heat content and age
        !
        IF(PEK%TSNOW%SCHEME=='3-L'.OR.PEK%TSNOW%SCHEME=='CRO')THEN
          DO JL = 1,INS
             PEK%TSNOW%HEAT(JI,JL) = PNEAR%TSNOW%HEAT(IGLO,JL)
             PEK%TSNOW%AGE (JI,JL) = PNEAR%TSNOW%AGE (IGLO,JL)
          ENDDO
        ENDIF
        !
        !* Optical Diameter, Sphericity and History
        !
        IF(PEK%TSNOW%SCHEME=='CRO')THEN
          DO JL = 1,INS
             PEK%TSNOW%DIAMOPT(JI,JL) = PNEAR%TSNOW%DIAMOPT(IGLO,JL)
             PEK%TSNOW%SPHERI (JI,JL) = PNEAR%TSNOW%SPHERI (IGLO,JL)
             PEK%TSNOW%HIST   (JI,JL) = PNEAR%TSNOW%HIST   (IGLO,JL)
          ENDDO
        ENDIF
        !
        !* aerodynamical resistance
        !
        PEK%XRESA(JI) = PNEAR%XRESA(IGLO)
        !
        !* Leaf Area Index
        !
        IF(IO%CPHOTO/='NON'.AND.IO%CPHOTO/='AGS'.AND.IO%CPHOTO/='AST')THEN
           PEK%XLAI(JI) = PEK%XLAIMIN(JI)
        ENDIF
        !
        !* Assimilation and evapotranspiration
        !
        IF (IO%CPHOTO/='NON') THEN
           PEK%XAN   (JI) = 0.0 
           PEK%XANDAY(JI) = 0.0 
           PEK%XANFM (JI) = XANFMINIT 
        ENDIF
        !
        !* biomass (similar computation as in prep are done using LAIMIN)
        !
        IF(IO%CPHOTO=='NIT'.OR.IO%CPHOTO=='NCB')THEN
          PEK%XRESP_BIOMASS(JI,:) = 0.0 
          PEK%XBIOMASS     (JI,:) = 0.0
          PEK%XBIOMASS     (JI,1) = PEK%XLAIMIN(JI) * PK%XBSLAI_NITRO(JI)                                 ! Parameter initialized
          PEK%XBIOMASS     (JI,2) = MAX(0.,EXP(ZINVCA*LOG(PEK%XBIOMASS(JI,1)/ZCC_CA))-PEK%XBIOMASS(JI,1)) ! Optimization : X**n = EXP(n*LOG(X))
        ENDIF
        !
        !* MEB Prognostic or Semi-prognostic variables
        !
        IF(ISIZE_LMEB_PATCH>0)THEN
          !
          !* liquid water retained on litter
          PEK%XWRL(JI) = PNEAR%XWRL(IGLO)
          !
          !* ice retained on litter
          PEK%XWRLI(JI) = PNEAR%XWRLI(IGLO)
          !
          !* snow retained on the foliage (as for WR)
          PEK%XWRVN(JI) = 0.0
          !
          !* canopy vegetation temperature
          PEK%XTV(JI) = PNEAR%XTV(IGLO)
          !
          !* litter temperature
          PEK%XTL(JI) = PNEAR%XTL(IGLO)
          !
          !* canopy air temperature
          PEK%XTC(JI) = PNEAR%XTC(IGLO)
          !
        ENDIF
        !
        !* Fire scheme
        !
        IF(IO%LFIRE)THEN
          PEK%XFIREIND      (JI)= PNEAR%XFIREIND(IGLO)
          PEK%XMOISTLIT_FIRE(JI)= PNEAR%XMOISTLIT_FIRE(IGLO)
          PEK%XTEMPLIT_FIRE (JI)= PNEAR%XTEMPLIT_FIRE(IGLO)
        ENDIF
        !
        !* litter and soil carbon (previous mean litter and soil carbon are attributed to the new patch)
        !
        IF(IO%CRESPSL=='CNT')THEN
          !       
          DO JNL=1,INLIT
             DO JNLS=1,INLITS
                PEK%XLITTER(JI,JNL,JNLS) = PNEAR%XLITTER(IGLO,JNL,JNLS)
             ENDDO
          ENDDO
          !
          DO JNLS=1,INLITS
             PEK%XLIGNIN_STRUC(JI,JNLS) = PNEAR%XLIGNIN_STRUC(IGLO,JNLS)
          ENDDO
          !
          DO JNC=1,INC
             PEK%XSOILCARB(JI,JNC) = PNEAR%XSOILCARB(IGLO,JNC) 
          ENDDO
          !
        ELSEIF(IO%CRESPSL=='DIF') THEN
          !
          DO JNL=1,INLIT
             PEK%XSURFACE_LITTER(JI,JNL)=PNEAR%XSURFACE_LITTER(IGLO,JNL)
          ENDDO
          !
          DO JNL=1,INLIT
             DO JL=1,INL
                IF(JL<=PK%NWG_LAYER(JI).AND.PK%NWG_LAYER(JI)/=NUNDEF.AND.PNEAR%XSOILDIF_LITTER(IGLO,JL,JNL)==XUNDEF)THEN
                  PEK%XSOILDIF_LITTER(JI,JL,JNL) = 0.0
                ELSE
                  PEK%XSOILDIF_LITTER(JI,JL,JNL) = PNEAR%XSOILDIF_LITTER(IGLO,JL,JNL)
                ENDIF
             ENDDO
          ENDDO
          !
          PEK%XSURFACE_LIGNIN_STRUC(JI) = PNEAR%XSURFACE_LIGNIN_STRUC(IGLO)
          !
          DO JL=1,INL
             IF(JL<=PK%NWG_LAYER(JI).AND.PK%NWG_LAYER(JI)/=NUNDEF.AND.PNEAR%XSOILDIF_LIGNIN_STRUC(IGLO,JL)==XUNDEF)THEN
               PEK%XSOILDIF_LIGNIN_STRUC(JI,JL) = 0.0
              ELSE
               PEK%XSOILDIF_LIGNIN_STRUC(JI,JL) = PNEAR%XSOILDIF_LIGNIN_STRUC(IGLO,JL)
              ENDIF
          ENDDO
          !
          DO JNC=1,INC
             DO JL=1,INL
                IF(JL<=PK%NWG_LAYER(JI).AND.PK%NWG_LAYER(JI)/=NUNDEF.AND.PNEAR%XSOILDIF_CARB(IGLO,JL,JNC)==XUNDEF)THEN
                  PEK%XSOILDIF_CARB(JI,JL,JNC) = 0.0
                ELSE
                  PEK%XSOILDIF_CARB(JI,JL,JNC) = PNEAR%XSOILDIF_CARB(IGLO,JL,JNC)
                ENDIF
             ENDDO
          ENDDO
          !
        ENDIF
        !
        !* Soil gas scheme
        !
        IF(IO%LSOILGAS)THEN
          PEK%XSGASO2 (JI,:) = PNEAR%XSGASO2 (IGLO,:)
          PEK%XSGASCO2(JI,:) = PNEAR%XSGASCO2(IGLO,:)
          PEK%XSGASCH4(JI,:) = PNEAR%XSGASCH4(IGLO,:)
        ENDIF       
        !
      ENDIF ! end of emerging or luluccf-driven changing pfts distribution
      !
   ENDDO  ! end of grid-cell loop
   !
ENDDO ! end of patch loop
!
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('LANDUSE_NEWPATCH',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
END SUBROUTINE LANDUSE_NEWPATCH
