!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!     #########
    SUBROUTINE LANDUSE_BIOMASS(IG, IO, S, NP, NPE, NPGLO, NPEGLO, TLU, KI, KLUOUT)
!   ##############################################################################
!!****  *LAND USE BIOMASS*
!!
!!    PURPOSE
!!    -------
!
!     Update and conserv biomass after land-use change
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
USE MODD_SFX_GRID_n,     ONLY : GRID_t
USE MODD_ISBA_OPTIONS_n, ONLY : ISBA_OPTIONS_t
USE MODD_ISBA_n,         ONLY : ISBA_S_t, ISBA_P_t, ISBA_PE_t, ISBA_K_t, ISBA_NP_t, ISBA_NPE_t
USE MODD_INIT_LANDUSE,   ONLY : LULCC_t
!
USE MODD_CO2V_PAR, ONLY : XPCCO2, XKGTOG, XGTOKG
!
USE MODD_SURF_PAR,ONLY : XUNDEF                 
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
!*      0.2    declarations of local parameter
!
TYPE(ISBA_P_t),  POINTER :: PK, PKGLO
TYPE(ISBA_PE_t), POINTER :: PEK, PEKGLO 
!
REAL, PARAMETER                            :: ZLIMIT = 1.0E-12
!
REAL, DIMENSION(KI,IO%NNBIOMASS,IO%NPATCH) :: ZTURNOVER    ! gC/m2
!
REAL, DIMENSION(KI,IO%NNBIOMASS)   :: ZBIOMASS_RESERV_OLD  ! Kg/m2
REAL, DIMENSION(KI,IO%NNBIOMASS)   :: ZBIOMASS_RESERV_NEW  ! Kg/m2
REAL, DIMENSION(KI,IO%NNBIOMASS)   :: ZTURNOVER_GRID       ! gC/m2

REAL, DIMENSION(KI)                :: ZBUDGET, ZCONSERV
REAL, DIMENSION(KI)                :: ZTURNOVER_TOT
!
REAL, DIMENSION(KI)                :: ZOLDVEG_GRID
REAL, DIMENSION(KI)                :: ZNEWVEG_GRID
!
LOGICAL :: LSTOP
!
INTEGER :: INP, INB, JGLO ! loop counter on levels
INTEGER :: JI, JP, JNB    ! loop counter on levels
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!-----------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('LANDUSE_BIOMASS',0,ZHOOK_HANDLE)
!-----------------------------------------------------------------
!
!*      1.     Preliminaries
!              -------------
!
INB=IO%NNBIOMASS
INP=IO%NPATCH
!
ZCONSERV     (:) = 0.0
ZBUDGET      (:) = 0.0
ZTURNOVER_TOT(:) = 0.0
ZNEWVEG_GRID (:) = 0.0
ZOLDVEG_GRID (:) = 0.0
!
ZBIOMASS_RESERV_OLD(:,:) = 0.0
ZBIOMASS_RESERV_NEW(:,:) = 0.0
ZTURNOVER_GRID     (:,:) = 0.0
!
ZTURNOVER(:,:,:) = 0.0
!
!Biomass in kg/m2
!
DO JP=1,INP
   PK => NP%AL(JP)
   PEK => NPE%AL(JP)
   DO JNB=1,INB
      DO JI=1,PK%NSIZE_P
         IF(PK%XPATCH(JI)==0.OR.PEK%XBIOMASS(JI,JNB)<ZLIMIT)THEN
            PEK%XBIOMASS     (JI,JNB)=0.0
            PEK%XRESP_BIOMASS(JI,JNB)=0.0
         ENDIF
      ENDDO
   ENDDO
ENDDO
!
DO JP=1,INP
   PKGLO => NPGLO%AL(JP)
   PEKGLO => NPEGLO%AL(JP)
   DO JNB=1,INB
      DO JGLO=1,KI
         IF(PKGLO%XPATCH(JGLO)==0.OR.PEKGLO%XBIOMASS(JGLO,JNB)<ZLIMIT)THEN
            PEKGLO%XBIOMASS(JGLO,JNB)=0.0
         ENDIF      
         TLU%XBIOM_GRID_OLD(JGLO) = TLU%XBIOM_GRID_OLD(JGLO) + PEKGLO%XBIOMASS(JGLO,JNB)*PKGLO%XPATCH(JGLO)
      ENDDO
   ENDDO
ENDDO
!
!-----------------------------------------------------------------
!
!*      2.     NIT case : Simple global grid-cell conservation for biomass
!              -----------------------------------------------------------
!
IF(IO%CPHOTO=='NIT')THEN
  !
  DO JP=1,INP
     PKGLO => NPGLO%AL(JP)
     PEKGLO => NPEGLO%AL(JP)
     DO JNB=1,INB
        DO JGLO=1,KI
           ZBIOMASS_RESERV_OLD(JGLO,JNB) = ZBIOMASS_RESERV_OLD(JGLO,JNB) + PEKGLO%XBIOMASS(JGLO,JNB)*PKGLO%XPATCH(JGLO)
        ENDDO
     ENDDO
  ENDDO
  !
  DO JP=1,INP
     PK => NP%AL(JP)
     PEK => NPE%AL(JP)
     DO JNB=1,INB
        DO JI=1,PK%NSIZE_P
           JGLO = PK%NR_P(JI)
           ZBIOMASS_RESERV_NEW(JGLO,JNB) = ZBIOMASS_RESERV_NEW(JGLO,JNB) + PEK%XBIOMASS(JI,JNB)*PK%XPATCH(JI)
        ENDDO
     ENDDO
  ENDDO
  !
  DO JP=1,INP 
     PK => NP%AL(JP)
     PEK => NPE%AL(JP)
     DO JNB=1,INB
        DO JI=1,PK%NSIZE_P
           JGLO = PK%NR_P(JI)
           IF(ZBIOMASS_RESERV_NEW(JGLO,JNB)/=ZBIOMASS_RESERV_OLD(JGLO,JNB).AND.ZBIOMASS_RESERV_NEW(JGLO,JNB)/= 0.0)THEN
             PEK%XBIOMASS(JI,JNB) = PEK%XBIOMASS(JI,JNB) * ZBIOMASS_RESERV_OLD(JGLO,JNB)/ZBIOMASS_RESERV_NEW(JGLO,JNB)
          ENDIF
        ENDDO
     ENDDO
  ENDDO
  !
ENDIF
!
!-----------------------------------------------------------------
!
!*      3.     ISBA-CC case : Update current year biomass
!              ------------------------------------------
!   
IF(IO%CPHOTO=='NCB')THEN  
  !
  !* The patch has grown = we "sow" a seed of zero biomass
  !  -----------------------------------------------------
  !
  !  Biomass=[ Biomass*Patch_old + 0.0 * (Patch_new-Patch_old) ] / (Patch_old+Patch_new-Patch_old) = Biomass*Patch_old/Patch_new                     
  !
  DO JP=1,INP  
     PK => NP%AL(JP)
     PEK => NPE%AL(JP)
     PKGLO => NPGLO%AL(JP)
     PEKGLO => NPEGLO%AL(JP)
     DO JNB=1,INB
        DO JI=1,PK%NSIZE_P
           JGLO = PK%NR_P(JI)
           IF(PK%XPATCH(JI)>PKGLO%XPATCH(JGLO).AND.PKGLO%XPATCH(JGLO)>0.0)THEN
             PEK%XBIOMASS(JI,JNB) = PEKGLO%XBIOMASS(JGLO,JNB) * PKGLO%XPATCH(JGLO) / PK%XPATCH(JI)
           ENDIF
        ENDDO
     ENDDO
  ENDDO
  !
  !* The patch has diminished or disappear
  !  -------------------------------------
  !
  !(1) we conserve biomass density                      
  !    Bio_new =[ Bio_old*Patch_old - Bio_old * (Patch_old-Patch_new) ] / Patch_new = Bio_old 
  !
  DO JP=1,INP 
     !
     PKGLO => NPGLO%AL(JP)
     PEKGLO => NPEGLO%AL(JP)
     !
     DO JGLO=1,KI
        !
        IF(S%XPATCH(JGLO,JP)<PKGLO%XPATCH(JGLO).AND.PKGLO%XPATCH(JGLO)>0.0)THEN
          !          
          !(2) we harvest the leaves and stems based on the old grid (kgC/m2 including patch)
          !    Carbon due to biomass harvested after land use change including patch fraction (kgC/m2)
          TLU%XLULCC_HARVEST(JGLO,JP) = (PEKGLO%XBIOMASS(JGLO,1)+PEKGLO%XBIOMASS(JGLO,2)+ &
                                         PEKGLO%XBIOMASS(JGLO,3)+PEKGLO%XBIOMASS(JGLO,5)) &
                                      * XPCCO2 * (PKGLO%XPATCH(JGLO)-S%XPATCH(JGLO,JP))
          !     
          !(3) we leave roots to die (gC/m2 including patch)
          ZTURNOVER(JGLO,4,JP) = PEKGLO%XBIOMASS(JGLO,4) * (PKGLO%XPATCH(JGLO)-S%XPATCH(JGLO,JP)) * (XPCCO2*XKGTOG)
          ZTURNOVER(JGLO,6,JP) = PEKGLO%XBIOMASS(JGLO,6) * (PKGLO%XPATCH(JGLO)-S%XPATCH(JGLO,JP)) * (XPCCO2*XKGTOG)
          !
        ENDIF
        !
     ENDDO
     !
  ENDDO
  !
  ! Biomass carbon transferred to litter pools due to lulucf processes (kgC m-2)
  !
  S%XFLURES(:,:) = (ZTURNOVER(:,4,:)+ZTURNOVER(:,6,:)) * XGTOKG
  !
  ! We add roots to turnover in the new grid (gC/m2)
  !
  ZTURNOVER_GRID(:,:) = 0.0
  DO JP=1,INP  
     DO JGLO=1,KI
        ZTURNOVER_GRID(JGLO,4)=ZTURNOVER_GRID(JGLO,4)+ZTURNOVER(JGLO,4,JP)
        ZTURNOVER_GRID(JGLO,6)=ZTURNOVER_GRID(JGLO,6)+ZTURNOVER(JGLO,6,JP)
     ENDDO
  ENDDO
  !
  ZNEWVEG_GRID(:) = 0.0
  DO JP=1,INP  
     PK => NP%AL(JP)
     PEK => NPE%AL(JP)
     DO JI=1,PK%NSIZE_P       
        IF(PEK%XLAI(JI)/=XUNDEF)THEN
          JGLO = PK%NR_P(JI)
          ZNEWVEG_GRID(JGLO)=ZNEWVEG_GRID(JGLO)+PK%XPATCH(JI)
        ENDIF 
     ENDDO
  ENDDO     
  !
  DO JP=1,INP  
     PK => NP%AL(JP)
     PEK => NPE%AL(JP)
     DO JI=1,PK%NSIZE_P
        JGLO = PK%NR_P(JI)
        IF(ZNEWVEG_GRID(JGLO)>0.0.AND.PEK%XLAI(JI)/=XUNDEF)THEN
          PK%XTURNOVER(JI,4)=ZTURNOVER_GRID(JGLO,4)/ZNEWVEG_GRID(JGLO)
          PK%XTURNOVER(JI,6)=ZTURNOVER_GRID(JGLO,6)/ZNEWVEG_GRID(JGLO)
        ENDIF
     ENDDO
  ENDDO 
  !
  ! all vegetated patches disappeared -> try to conserv
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
  DO JP=1,INP
     DO JGLO=1,KI
        IF(ZOLDVEG_GRID(JGLO)>0.0.AND.ZNEWVEG_GRID(JGLO)==0.0)THEN 
          ZCONSERV(JGLO) = ZCONSERV(JGLO) + (ZTURNOVER(JGLO,4,JP)+ZTURNOVER(JGLO,6,JP))* XGTOKG
        ENDIF
     ENDDO
  ENDDO 
  !
  ! * The patch is new : 
  !  -------------------
  !
  DO JP=1,INP  
     PK => NP%AL(JP)
     PEK => NPE%AL(JP)
     PKGLO => NPGLO%AL(JP)
     DO JI=1,PK%NSIZE_P
        JGLO = PK%NR_P(JI)
        !We use a negative turnover to conserv global carbon stock
        IF(PK%XPATCH(JI)>0.0.AND.PKGLO%XPATCH(JGLO)==0.0)THEN
           PK%XTURNOVER(JI,1) = - PEK%XBIOMASS(JI,1) * (XPCCO2*XKGTOG)
           PK%XTURNOVER(JI,2) = - PEK%XBIOMASS(JI,2) * (XPCCO2*XKGTOG)
        ENDIF
     ENDDO
  ENDDO 
  !
  ! Biomass carbon transferred to litter pools due to lulucf processes (kgC m-2)
  !
  DO JP=1,INP  
     PK => NP%AL(JP)
     DO JI=1,PK%NSIZE_P
        JGLO = PK%NR_P(JI)
        S%XFLURES(JGLO,JP) = S%XFLURES(JGLO,JP)+(PK%XTURNOVER(JI,1)+PK%XTURNOVER(JI,2))*PK%XPATCH(JI)*XGTOKG
     ENDDO
  ENDDO 
  !
  DO JP=1,INP  
     DO JGLO=1,KI
        TLU%XLULCC_HARVEST_GRID(JGLO) = TLU%XLULCC_HARVEST_GRID(JGLO)+TLU%XLULCC_HARVEST(JGLO,JP)
     ENDDO
  ENDDO 
  !
  ZTURNOVER_TOT(:)=0.0
  DO JP=1,INP  
     PK => NP%AL(JP)
     DO JI=1,PK%NSIZE_P
        JGLO = PK%NR_P(JI)
        ZTURNOVER_TOT(JGLO)=ZTURNOVER_TOT(JGLO)+(PK%XTURNOVER(JI,1)+PK%XTURNOVER(JI,2) &
                                               + PK%XTURNOVER(JI,4)+PK%XTURNOVER(JI,6))&
                                               * PK%XPATCH(JI)*(XGTOKG/XPCCO2)
     ENDDO
  ENDDO
!
ENDIF
!
!-----------------------------------------------------------------
!
!*      4.     Compute current year total grid-cell biomass
!              --------------------------------------------
!
! * Biomass in kg/m2
!
DO JP=1,INP
   PK => NP%AL(JP)
   PEK => NPE%AL(JP)
   DO JNB=1,INB
      DO JI=1,PK%NSIZE_P
         JGLO = PK%NR_P(JI)
         TLU%XBIOM_GRID_NEW(JGLO) = TLU%XBIOM_GRID_NEW(JGLO) + PEK%XBIOMASS(JI,JNB)*PK%XPATCH(JI)
      ENDDO
   ENDDO
ENDDO
!
!-----------------------------------------------------------------
!
!*      7.     Test total biomass conservation
!              -------------------------------
!
!* Biomass conservation (kg/m2)
!
S%XCCONSRV(:) = ZCONSERV(:)
!
ZBUDGET(:) = TLU%XBIOM_GRID_NEW(:)-TLU%XBIOM_GRID_OLD(:)
!
IF(IO%CPHOTO=='NCB')THEN  
  ZBUDGET(:) = ZBUDGET(:) + ZTURNOVER_TOT(:) + (ZCONSERV(:)+TLU%XLULCC_HARVEST_GRID(:))/XPCCO2
ENDIF
!
LSTOP=.FALSE.
DO JGLO=1,KI
   IF(ABS(ZBUDGET(JGLO))>1.E-12)THEN
     WRITE(KLUOUT,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
     WRITE(KLUOUT,*)'LANDUSE_BIOMASS: NO CONSERVATION IN AT LEAST ONE GRID CELL'
     WRITE(KLUOUT,*)'LON = ',IG%XLON(JGLO),' LAT =',IG%XLAT(JGLO),'BUDGET =',ZBUDGET(JGLO),'kg/m2'
     WRITE(KLUOUT,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
     LSTOP=.TRUE.
   ENDIF           
ENDDO
IF(LSTOP) CALL ABOR1_SFX('LANDUSE_BIOMASS: INCONSISTENCY IN BIOMASS BUDGET')
!
IF (LHOOK) CALL DR_HOOK('LANDUSE_BIOMASS',1,ZHOOK_HANDLE)
!-----------------------------------------------------------------
!
END SUBROUTINE LANDUSE_BIOMASS
