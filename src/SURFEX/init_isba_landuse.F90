!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!#############################################################
SUBROUTINE INIT_ISBA_LANDUSE (IG, IO, S, K, NK, NP, NPE, SOLD, NPGLO, NPEGLO, HPROGRAM, KI)  
!#############################################################
!
!!****  *INIT_ISBA_LANDUSE* - routine to initialize land use for ISBA field
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
!!      Original    07/2011
!!      Completelly reframed 08/2016 R. Séférian
!!      R. Seferian 10/2016 correct error in landuse computation fields
!!      R. Seferian 11/2016 : add cmip6 diagnostics
!!      J. Colin    12/2017 : add computations in case the water or snow is
!!                            nudged seperately on each patch
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
USE MODD_SURF_PAR,       ONLY : XUNDEF, NUNDEF, LEN_HREC
!
USE MODD_CSTS,           ONLY : XDAY
!
USE MODD_DATA_COVER_PAR, ONLY : NVEGTYPE
!
USE MODD_CO2V_PAR,       ONLY : XPCCO2,XKGTOG
!
USE MODD_SFX_OASIS,       ONLY :LCPL_RIVCARB
!
USE MODI_ABOR1_SFX
!
USE MODI_READ_SURF
USE MODE_READ_SURF_LAYERS
!
USE MODI_GET_LUOUT
!
USE MODI_LANDUSE_NEWPATCH
USE MODI_LANDUSE_HYDRO 
USE MODI_LANDUSE_BIOMASS
USE MODI_LANDUSE_CARBON
USE MODI_LANDUSE_LITTER_DIF
USE MODI_LANDUSE_CARBON_DIF
USE MODI_LANDUSE_SOILGAS
USE MODI_LANDUSE_CARBON_MANAGING
USE MODI_LANDUSE_NUDGING 
!
USE YOMHOOK,  ONLY : LHOOK,   DR_HOOK
USE PARKIND1, ONLY : JPRB
!
IMPLICIT NONE
!
!*       0.1   Declarations of arguments
!              -------------------------
!
TYPE(GRID_t),          INTENT(INOUT) :: IG
TYPE(ISBA_OPTIONS_t),  INTENT(INOUT) :: IO
TYPE(ISBA_S_t),        INTENT(INOUT) :: S, SOLD
TYPE(ISBA_K_t),        INTENT(INOUT) :: K
TYPE(ISBA_NK_t),       INTENT(INOUT) :: NK
TYPE(ISBA_NP_t),       INTENT(INOUT) :: NP, NPGLO
TYPE(ISBA_NPE_t),      INTENT(INOUT) :: NPE, NPEGLO
!
CHARACTER(LEN=6),      INTENT(IN)    :: HPROGRAM          ! program calling surf. schemes
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
REAL, DIMENSION(KI,IO%NGROUND_LAYER,IO%NPATCH) :: ZWORK3D        ! previous year soil layer depth
REAL, DIMENSION(KI,                 IO%NPATCH) :: ZWORK          ! work array (ISBA grid)
REAL, DIMENSION(KI                           ) :: ZBUDGET        ! Carbon budget
!
CHARACTER(LEN=LEN_HREC)                        :: YRECFM         ! Name of the article to be read
CHARACTER(LEN=4)                               :: YLVL
INTEGER                                        :: IRESP          ! Error code after redding
!
LOGICAL :: GDIM, GSTOP, GLULU_LAYER
INTEGER :: ILUOUT                                  ! unit of output listing file
INTEGER :: JI, JJ, JL, JNL, JNLS, JNC, JP, JP_NEAR ! loop counter
INTEGER :: INP, INL, INS, INB, INLIT, INLITS, INC  ! dimension
INTEGER :: JT, INTIME                              ! loop on time (nudging) and size
INTEGER :: JGLO, ISIZE_LMEB_PATCH, IWORK           ! Work integer
REAL    :: ZNDAYS                                  ! Work real
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('INIT_ISBA_LANDUSE',0,ZHOOK_HANDLE)
!-------------------------------------------------------------------------------
!        
!
!*       0. Global initialisation
!        ------------------------
!
CALL GET_LUOUT(HPROGRAM,ILUOUT)
!
INP   =IO%NPATCH
INL   =IO%NGROUND_LAYER
INS   =NPE%AL(1)%TSNOW%NLAYER
INB   =IO%NNBIOMASS
INLIT =IO%NNLITTER
INLITS=IO%NNLITTLEVS
INC   =IO%NNSOILCARB
!
ISIZE_LMEB_PATCH=COUNT(IO%LMEB_PATCH(:))
!
!* number of day by year
!
ZNDAYS=365.
IF(((MOD(S%TTIME%TDATE%YEAR,4)==0).AND.(MOD(S%TTIME%TDATE%YEAR,100)/=0)).OR.(MOD(S%TTIME%TDATE%YEAR,400)==0))THEN
  ZNDAYS=366.
ENDIF
!
!* initialize total co2 land use flux to atm
!
ALLOCATE(S%XFLUCLEARTOATM(KI))
ALLOCATE(S%XFHARVESTTOATM(KI))
S%XFLUCLEARTOATM(:)=0.0
S%XFHARVESTTOATM(:)=0.0
!
!* find vanishing layer
!
GLULU_LAYER=IO%LLULU
DO JP=1,INP
   PK => NP%AL(JP)
   PEK => NPE%AL(JP)
   DO JI=1,PK%NSIZE_P
      IWORK=PK%NWG_LAYER(JI)
      IF(PK%NWG_LAYER(JI)/=NUNDEF)THEN
        IF(.NOT.IO%LLULU.AND.PEK%XWG(JI,IWORK)==XUNDEF)GLULU_LAYER=.TRUE.
      ENDIF
   ENDDO
ENDDO
!
!-------------------------------------------------------------------------------
!
IF(IO%LLULU)THEN
  !
  WRITE(ILUOUT,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
  WRITE(ILUOUT,*)'!!!                                                        !!!'
  WRITE(ILUOUT,*)'!!!                    WARNING    WARNING                  !!!'
  WRITE(ILUOUT,*)'!!!                                                        !!!'
  WRITE(ILUOUT,*)'!!!              Patches Distribution has changed          !!!'
  WRITE(ILUOUT,*)'!!!  Land-use Land Cover Change computation are performed  !!!'
  WRITE(ILUOUT,*)'!!!                                                        !!!'
  WRITE(ILUOUT,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
  !
  WRITE(*,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
  WRITE(*,*)'!!!                                                        !!!'
  WRITE(*,*)'!!!                    WARNING    WARNING                  !!!'
  WRITE(*,*)'!!!                                                        !!!'
  WRITE(*,*)'!!!              Patches Distribution has changed          !!!'
  WRITE(*,*)'!!!  Land-use Land Cover Change computation are performed  !!!'
  WRITE(*,*)'!!!                                                        !!!'
  WRITE(*,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'     
  !
  !
  !*       1. Previous vertical grid
  !        -------------------------
  !
  !* Initialise or allocate some parameters
  !
  CALL INIT_LULCC(TLU)  
  !
  CALL ALLOC_LULCC(TLU,KI,IO%NPATCH)  
  !
  DO JP = 1,INP
     PKGLO => NPGLO%AL(JP)
     ALLOCATE(PKGLO%XDG (KI,INL))
     ALLOCATE(PKGLO%XDZG(KI,INL))
     ALLOCATE(PKGLO%NWG_LAYER(KI))
     PKGLO%XDG    (:,:) = XUNDEF
     PKGLO%XDZG   (:,:) = XUNDEF
     PKGLO%NWG_LAYER(:) = NUNDEF
  ENDDO
  !
  !*  Read previous DG 
  !
  CALL READ_SURF(HPROGRAM,'SPLIT_PATCH',GDIM,IRESP)
  CALL READ_SURF_LAYERS(HPROGRAM,'DG',GDIM,ZWORK3D(:,:,:),IRESP)
  !
  DO JP = 1,INP
     PKGLO => NPGLO%AL(JP)
     PKGLO%XDG (:,:) = ZWORK3D(:,:,JP)
     PKGLO%XDZG(:,1) = PKGLO%XDG(:,1)
  ENDDO
  !
  DO JL=2,INL
     DO JP = 1,INP
        PKGLO => NPGLO%AL(JP)
        PKGLO%XDZG(:,JL) = PKGLO%XDG(:,JL)-PKGLO%XDG(:,JL-1)
     ENDDO
  ENDDO
  !
  !* Previous NWG_LAYER
  !
  DO JP = 1,INP
     PKGLO => NPGLO%AL(JP)
     PEKGLO => NPEGLO%AL(JP)
     DO JL=1,INL
        DO JGLO=1,KI
           IF(PKGLO%XPATCH(JGLO)>0.0.AND.PEKGLO%XWG(JGLO,JL)/=XUNDEF)THEN
              PKGLO%NWG_LAYER(JGLO) = JL
           ENDIF
        ENDDO
     ENDDO
  ENDDO
  !
  !* Vegetation interactive case
  !
  IF(IO%CPHOTO=='NIT'.OR.IO%CPHOTO=='NCB')THEN
    DO JP = 1,INP
       PK => NP%AL(JP)
       PK%XTURNOVER(:,:) = 0.0
    ENDDO
  ENDIF
  !
  !-----------------------------------------------------------------------------
  !
  !*       2. Treat case when new PFTs is created due to LULCC
  !        ---------------------------------------------------
  !
  CALL LANDUSE_NEWPATCH (IO, K, NK, NP, NPE, SOLD, NPGLO, NPEGLO, KI)
  !
  !-----------------------------------------------------------------------------
  !
ENDIF ! End of Land-use Land Cover Change case
!
!
!*       3. Conserv mass and energy
!        --------------------------
!
!
!*ISBA-DF case : Update water content and conserv
!
IF(IO%CISBA=='DIF'.AND.GLULU_LAYER)THEN
  !
  !initialize annual resiual conservation flux
  S%XWCONSRV(:) = 0.0
  !
  !compute conservation
  CALL LANDUSE_HYDRO(IO, S, NK, NP, NPE, NPGLO, NPEGLO, KI)
  !
  !annual resiual water conservation flux : kg/m2 by year to kg/m2/s
  S%XWCONSRV(:) = S%XWCONSRV(:)/(ZNDAYS*XDAY)
  !
ENDIF
!
!*ISBA-Ags case : Update biomass
!not done
!
!*ISBA-CC case : Update biomass and carbon stocks and conserv
!
IF(IO%LLULU)THEN
  !
  IF(IO%CPHOTO=='NCB')THEN
    !
    !initialize annual fluxes
    S%XFLURES (:,:) = 0.0
    S%XFLUATM (:,:) = 0.0
    S%XCCONSRV(:  ) = 0.0
    !
    !compute conservation
    CALL LANDUSE_BIOMASS(IG, IO, S, NP, NPE, NPGLO, NPEGLO, TLU, KI, ILUOUT)
    !
    IF(IO%CRESPSL=='CNT')THEN
      CALL LANDUSE_CARBON(IO, S, NP, NPE, NPGLO, NPEGLO, TLU, KI)
    ELSEIF(IO%CRESPSL=='DIF')THEN
      CALL LANDUSE_LITTER_DIF(IG, IO, S, NP, NPE, NPGLO, NPEGLO, TLU, KI, ILUOUT)   
      CALL LANDUSE_CARBON_DIF(IO, S, NP, NPE, NPGLO, NPEGLO, TLU, KI)   
    ENDIF
    !
    IF(IO%LSOILGAS)THEN
      CALL LANDUSE_SOILGAS(IO, S, NP, NPE, NPGLO, NPEGLO, KI)
    ENDIF
    !
    ZBUDGET(:) = (TLU%XBIOM_GRID_NEW     (:)-TLU%XBIOM_GRID_OLD  (:))*XPCCO2*XKGTOG &
               + (TLU%XLITTER_GRID_NEW   (:)-TLU%XLITTER_GRID_OLD(:))               &
               + (TLU%XCSOIL_GRID_NEW    (:)-TLU%XCSOIL_GRID_OLD (:))               &
               + (TLU%XLULCC_HARVEST_GRID(:)+S%XCCONSRV          (:))*XKGTOG
    !
    GSTOP=.FALSE.
    DO JGLO=1,KI
       IF(ABS(ZBUDGET(JGLO))>1.0E-10)THEN
         WRITE(ILUOUT,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
         WRITE(ILUOUT,*)'INIT_ISBA_LANDUSE: NO CARBON CONSERVATION IN AT LEAST ONE GRID CELL'
         WRITE(ILUOUT,*)'LON = ',IG%XLON(JGLO),' LAT =',IG%XLAT(JGLO),'RESIDUE =',ZBUDGET(JGLO),'gC/m2'
         WRITE(ILUOUT,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
         GSTOP=.TRUE.
       ENDIF  
    ENDDO
    IF(GSTOP) CALL ABOR1_SFX('INIT_ISBA_LANDUSE: INCONSISTENCY IN CARBON BUDGET, SEE LISTING_OFF')
    !
    !*ISBA-CC case : annual resiual carbon conservation flux : kgC/m2 by year to kgC/m2/s
    !
    S%XCCONSRV(:) = S%XCCONSRV(:)/(ZNDAYS*XDAY)  
    !
    !*Managing land-use scheme for carbon cycle
    !
    IF(IO%LLULCC_MANAGE)THEN
      CALL LANDUSE_CARBON_MANAGING(IO, S, SOLD, TLU, KI)  ! Flux of carbon due to LULCC
    ELSE
      S%XFLUATM (:,:) = TLU%XLULCC_HARVEST(:,:)
    ENDIF
    !
    !*kgC/m2 by year to kgC/m2/s
    !
    S%XFLURES (:,:) = S%XFLURES (:,:)/(ZNDAYS*XDAY)
    S%XFLUATM (:,:) = S%XFLUATM (:,:)/(ZNDAYS*XDAY)
    !
    IF(IO%LLULCC_MANAGE)THEN
      S%XFLUANT (:,:) = S%XFLUANT (:,:)/(ZNDAYS*XDAY)
      S%XFANTATM(:,:) = S%XFANTATM(:,:)/(ZNDAYS*XDAY)
    ENDIF  
    !
  ENDIF
  !
  !-----------------------------------------------------------------------------
  !
  !*       4. Re-initialisation for normal ISBA computation
  !        ------------------------------------------------
  !
  CALL INIT_LULCC(TLU)   
  !
  IF(IO%CPHOTO=='NIT'.OR.IO%CPHOTO=='NCB')THEN           
    DO JP=1,INP
       PK => NP%AL(JP)
       PK%XTURNOVER(:,:) = 0.0
    ENDDO
  ENDIF
  !
  !           
ENDIF 
!
!
!-------------------------------------------------------------------------------
!        
!
!*       6. Coupling surface-atmosphere for CO2 flux
!        -------------------------------------------
!
IF(IO%CPHOTO=='NCB')THEN
  DO JP=1,INP
     DO JGLO=1,KI
        S%XFLUCLEARTOATM(JGLO)=S%XFLUCLEARTOATM(JGLO)+S%XFLUATM(JGLO,JP)
     ENDDO
  ENDDO
ENDIF
!
IF(IO%LLULCC_MANAGE)THEN
  DO JP=1,INP
     DO JGLO=1,KI
        S%XFLUCLEARTOATM(JGLO)=S%XFLUCLEARTOATM(JGLO)+S%XFANTATM(JGLO,JP)
     ENDDO
  ENDDO
ENDIF
!
! Carbon conservation in case where soil lixivation not done 
!
IF(IO%CPHOTO=='NCB'.AND..NOT.(LCPL_RIVCARB.AND.IO%LCLEACH))THEN
!IF(.NOT.(LCPL_RIVCARB.AND.IO%LCLEACH))THEN
  S%XFLUCLEARTOATM(:)=S%XFLUCLEARTOATM(:)+S%XCCONSRV(:)
ENDIF
!
!
!-------------------------------------------------------------------------------
!        
!
!*       7. Nudging case
!        ---------------
!
IF(IO%LNUDG_SWE.OR.IO%CNUDG_WG/='DEF')THEN
  CALL LANDUSE_NUDGING(IO, S, NK, NP, NPE, SOLD, NPGLO, HPROGRAM, KI)
ENDIF
!
!-------------------------------------------------------------------------------
!        
!
!*       8. Reset LULCC specific pointer
!        -------------------------------
!
CALL ISBA_S_INIT(SOLD)
IF(IO%LLULU)THEN
  CALL ISBA_NP_INIT(NPGLO,0)
  CALL ISBA_NPE_INIT(NPEGLO,0)
ENDIF
!
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('INIT_ISBA_LANDUSE',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
END SUBROUTINE INIT_ISBA_LANDUSE
