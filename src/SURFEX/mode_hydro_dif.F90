!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!     ######spl
      MODULE MODE_HYDRO_DIF 
!     ################
!
!!****  *MODE_HYDRO_DIF * - pedo-transfert functions
!!
!!    PURPOSE
!!    -------
!
!!**  METHOD
!!    ------
!!    
!!
!!    EXTERNAL
!!    --------
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
!!      B. Decharme       * Meteo France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original        11/2010
!!      A. Boone        02/2025 Added function to compute the molecular diffusivity of water,
!!                              and baresoil evap soil resistance options (Marti-Lopez et al., 2025)
!-----------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!
!
INTERFACE VAPCONDCF
  MODULE PROCEDURE VAPCONDCF       
END INTERFACE
!
INTERFACE INFMAX_FUNC
  MODULE PROCEDURE INFMAX_FUNC
END INTERFACE
!
INTERFACE TRIDIAG_DIF
  MODULE PROCEDURE TRIDIAG_DIF
END INTERFACE
!
INTERFACE MDIFFUSIVITY_WATER
  MODULE PROCEDURE MDIFFUSIVITY_WATER_0D
  MODULE PROCEDURE MDIFFUSIVITY_WATER_1D
END INTERFACE 
!
INTERFACE SOIL_RES_EG_S92
  MODULE PROCEDURE SOIL_RES_EG_S92
END INTERFACE SOIL_RES_EG_S92
!
INTERFACE SOIL_RES_EG_DSL
  MODULE PROCEDURE SOIL_RES_EG_DSL
END INTERFACE SOIL_RES_EG_DSL
!
!-------------------------------------------------------------------------------
CONTAINS
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
!vapor conductivity (m s-1)
!-------------------------------------------------------------------------------
!
FUNCTION VAPCONDCF(PTG,PPS,PWG,PWGI,PPSIA,PWSAT,PWFC,PQSAT,PQSATI,KWG_LAYER,KNL) RESULT(PVAPCOND)
!
! Uses method of Braud et al. (1993) for
!
USE MODD_CSTS,       ONLY : XMV, XMD, XTT, XP00, XG, XRV, XRHOLW
USE MODD_ISBA_PAR,   ONLY : XWGMIN
USE MODD_SURF_PAR,   ONLY : XUNDEF
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!
!*      0.1    declarations of arguments
!
REAL, DIMENSION(:  ), INTENT(IN)         :: PPS
REAL, DIMENSION(:,:), INTENT(IN)         :: PWG,PWGI,PPSIA,PWSAT, &
                                            PWFC,PTG,PQSAT,PQSATI
INTEGER, DIMENSION(:), INTENT(IN)        :: KWG_LAYER       !Moisture layer
INTEGER,               INTENT(IN)        :: KNL             ! number of vertical levels
!
REAL, DIMENSION(SIZE(PWG,1),SIZE(PWG,2)) :: PVAPCOND
!
!
!*      0.2    declarations of local variables
!
REAL    :: ZDVA, ZFVA, ZCHI, ZHUM, ZWORK,  &
           ZPV, ZESAT, ZESATI, ZWG, ZVC
!
INTEGER :: INI, JJ, JL, IDEPTH
!
! Parameters:
!
REAL, PARAMETER                     :: ZTORTY = 0.66         ! (-)
REAL, PARAMETER                     :: ZNV    = 1.88         ! (-)
REAL, PARAMETER                     :: ZCV    = 2.17e-5      ! (m2/s)
REAL, PARAMETER                     :: ZWK    = 0.05         ! (m3 m-3)
REAL, PARAMETER                     :: ZLIM   = TINY(1.0)
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('MODE_HYDRO_DIF:VAPCONDCF',0,ZHOOK_HANDLE)
!
INI = SIZE(PWG,1)
!
PVAPCOND(:,:) = 0.0
!
! Only perform this computation if the soil is sufficiently
! dry (as otherwise the hydraulic conductivity dominates
! the diffusion coefficient). Arbitrarily base threshold on field
! capacity water content:
!
DO JL=1,KNL
   DO JJ=1,INI
!
      IDEPTH = KWG_LAYER(JJ)
      ZWG    = PWG (JJ,JL) + PWGI(JJ,JL)
!      
      IF(JL<=IDEPTH .AND. ZWG < PWFC(JJ,JL) .AND. ZWG > XWGMIN)THEN
!
!        Vapor pressure over liquid and solid water surfaces (Pa), respectively:
!
         ZESAT  = PQSAT(JJ,JL)* PPS(JJ)/((XMV/XMD)+PQSAT(JJ,JL) *(1.-(XMV/XMD)))
!
         ZESATI = PQSATI(JJ,JL)*PPS(JJ)/((XMV/XMD)+PQSATI(JJ,JL)*(1.-(XMV/XMD)))
!
!        molecular diffusivity of water vapor (m2 s-1):
!
         ZDVA   = MDIFFUSIVITY_WATER(PTG(JJ,JL),PPS(JJ))
!
!        function of pore space: 
!
         ZFVA   = (PWSAT(JJ,JL) - ZWG)*(1.+(ZWG/(PWSAT(JJ,JL)-ZWK)))
         ZFVA   = MIN(ZFVA,PWSAT(JJ,JL))
!
!        relative humidity of air in soil pores:
!
         ZHUM   = MAX(ZLIM,EXP(PPSIA(JJ,JL)*XG/(XRV*PTG(JJ,JL))))
!
!        fraction of frozen water:
!
         ZCHI   = PWGI(JJ,JL)/ZWG
!
!        vapor pressure within pore space (Pa):
!
         ZPV    = ZHUM*(ZCHI*ZESAT + (1.-ZCHI)*ZESATI)
!
!        vapor conductivity (kg m-2 s-1)
!
         ZVC    = ZTORTY*PPS(JJ)*ZDVA*ZFVA*XG*ZPV/                  &
                  ((PPS(JJ)-ZPV)*(XRV*XRV*PTG(JJ,JL)*PTG(JJ,JL)))
!
!        vapor conductivity (m s-1)
!
         PVAPCOND(JJ,JL) = ZVC/XRHOLW
!
      ENDIF
!
   ENDDO
ENDDO
!
IF (LHOOK) CALL DR_HOOK('MODE_HYDRO_DIF:VAPCONDCF',1,ZHOOK_HANDLE)
!
END FUNCTION VAPCONDCF
!
!-------------------------------------------------------------------------------
!Green-Ampt approximation (derived form) for maximum infiltration
!-------------------------------------------------------------------------------
!
FUNCTION INFMAX_FUNC(PWG,PWSAT,PFRZ,PCONDSAT,PMPOTSAT,PBCOEF,PDZG,PDG,KLAYER_HORT)
USE YOMHOOK      ,ONLY : LHOOK,   DR_HOOK
USE MODD_SGH_PAR, ONLY : XHORT_DEPTH
USE PARKIND1     ,ONLY : JPRB
IMPLICIT NONE
REAL, DIMENSION(:,:), INTENT(IN) :: PWG
REAL, DIMENSION(:,:), INTENT(IN) :: PWSAT           
REAL, DIMENSION(:,:), INTENT(IN) :: PFRZ
REAL, DIMENSION(:,:), INTENT(IN) :: PCONDSAT            
REAL, DIMENSION(:,:), INTENT(IN) :: PMPOTSAT    
REAL, DIMENSION(:,:), INTENT(IN) :: PBCOEF 
REAL, DIMENSION(:,:), INTENT(IN) :: PDZG
REAL, DIMENSION(:,:), INTENT(IN) :: PDG
INTEGER,              INTENT(IN) :: KLAYER_HORT   
!
REAL, DIMENSION(SIZE(PWG,1)) :: ZGREEN_AMPT, ZDEPTH
REAL                         :: ZS, ZCOEF
INTEGER                      :: JJ,JL,INI
!
REAL, DIMENSION(SIZE(PWG,1)) :: INFMAX_FUNC
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
IF (LHOOK) CALL DR_HOOK('MODE_HYDRO_DIF:INFMAX_FUNC',0,ZHOOK_HANDLE)
!
INI   =SIZE(PWG,1)
!
ZGREEN_AMPT(:) = 0.0
ZDEPTH     (:) = 0.0
!
DO JL=1,KLAYER_HORT
   DO JJ=1,INI  
      IF(ZDEPTH(JJ)<XHORT_DEPTH)THEN
         ZS              = MIN(1.0,PWG(JJ,JL)/PWSAT(JJ,JL))
         ZCOEF           = PBCOEF(JJ,JL)*PMPOTSAT(JJ,JL)*(ZS-1.0)/PDZG(JJ,JL)        
         ZGREEN_AMPT(JJ) = ZGREEN_AMPT(JJ)+PDZG(JJ,JL)*PFRZ(JJ,JL)*PCONDSAT(JJ,JL)*(ZCOEF+1.0)
         ZDEPTH     (JJ) = PDG(JJ,JL)
      ENDIF      
   ENDDO
ENDDO
!
INFMAX_FUNC(:) = ZGREEN_AMPT(:)/ZDEPTH(:)
!
IF (LHOOK) CALL DR_HOOK('MODE_HYDRO_DIF:INFMAX_FUNC',1,ZHOOK_HANDLE)
END FUNCTION INFMAX_FUNC
!
!-------------------------------------------------------------------------------
!Solve tridiagonal matrix (for method see tridiag_ground.F90)
!-------------------------------------------------------------------------------
!
SUBROUTINE TRIDIAG_DIF(PAMTRX,PBMTRX,PCMTRX,PFRC,KWG_LAYER,KNL,PSOL)
USE MODD_SURF_PAR, ONLY : XUNDEF
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
IMPLICIT NONE
REAL,    DIMENSION(:,:), INTENT(IN)  :: PAMTRX    ! lower diag. elements of A matrix
REAL,    DIMENSION(:,:), INTENT(IN)  :: PBMTRX    ! main  diag. elements of A matrix
REAL,    DIMENSION(:,:), INTENT(IN)  :: PCMTRX    ! upper diag. elements of A matrix
REAL,    DIMENSION(:,:), INTENT(IN)  :: PFRC      ! Forcing term
INTEGER, DIMENSION(:),   INTENT(IN)  :: KWG_LAYER !Moisture layer
INTEGER,                 INTENT(IN)  :: KNL       ! number of vertical levels
REAL,    DIMENSION(:,:), INTENT(OUT) :: PSOL      ! solution of A.SOL = FRC
!
REAL, DIMENSION(SIZE(PFRC,1),SIZE(PFRC,2)) :: ZWORK! work array
REAL, DIMENSION(SIZE(PFRC,1))              :: ZDET ! work array
INTEGER                                    :: JL   ! vertical loop control
INTEGER :: JJ, INI, IDEPTH
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
IF (LHOOK) CALL DR_HOOK('MODE_HYDRO_DIF:TRIDIAG_DIF',0,ZHOOK_HANDLE)
!
INI=SIZE(PFRC,1)
!
PSOL (:,:)=XUNDEF
ZWORK(:,:)=XUNDEF
!
!first level
ZDET(:)   = PBMTRX(:,1)
PSOL(:,1) = PFRC(:,1) / ZDET(:)
!
!other levels
DO JL=2,KNL
   DO JJ=1,INI
      IDEPTH=KWG_LAYER(JJ)
      IF(JL<=IDEPTH)THEN
        ZWORK(JJ,JL) = PCMTRX(JJ,JL-1)/ZDET(JJ)
        ZDET (JJ)    = PBMTRX(JJ,JL) - PAMTRX(JJ,JL)*ZWORK(JJ,JL)
        PSOL (JJ,JL) = (PFRC (JJ,JL) - PAMTRX(JJ,JL)*PSOL(JJ,JL-1))/ZDET(JJ)  
      ENDIF
   ENDDO 
ENDDO        
!
!levels going down
DO JL=KNL-1,1,-1
   DO JJ=1,INI
      IDEPTH=KWG_LAYER(JJ)
      IF(JL<IDEPTH)THEN
         PSOL(JJ,JL) = PSOL(JJ,JL)-ZWORK(JJ,JL+1)*PSOL(JJ,JL+1)
      ENDIF
   ENDDO
ENDDO
!
IF (LHOOK) CALL DR_HOOK('MODE_HYDRO_DIF:TRIDIAG_DIF',1,ZHOOK_HANDLE)
!
END SUBROUTINE TRIDIAG_DIF
!
!-------------------------------------------------------------------------------
! Molecular diffusivity of water vapor (m2 s-1)
!-------------------------------------------------------------------------------
!
FUNCTION MDIFFUSIVITY_WATER_1D(PT,PPS) RESULT(ZDVA)
!
USE MODD_CSTS,    ONLY : XTT, XP00  
!  
USE YOMHOOK      ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1     ,ONLY : JPRB
!
IMPLICIT NONE

REAL, DIMENSION(:), INTENT(IN) :: PT
REAL, DIMENSION(:), INTENT(IN) :: PPS
!
REAL, DIMENSION(SIZE(PT))      :: ZWORK, ZDVA
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
REAL, PARAMETER                :: ZNV    = 1.88         ! (-)
REAL, PARAMETER                :: ZCV    = 2.17e-5      ! (m2/s)
!
IF (LHOOK) CALL DR_HOOK('MODE_HYDRO_DIF:MDIFFUSIVITY_WATER_1D',0,ZHOOK_HANDLE)
!
ZWORK(:)  = ZNV*LOG(PT(:)/XTT)
ZDVA(:)   = ZCV*(XP00/PPS(:))*EXP(ZWORK(:))
!
IF (LHOOK) CALL DR_HOOK('MODE_HYDRO_DIF:MDIFFUSIVITY_WATER_1D',1,ZHOOK_HANDLE)

END FUNCTION MDIFFUSIVITY_WATER_1D
!-------------------------------------------------------------------------------
FUNCTION MDIFFUSIVITY_WATER_0D(PT,PPS) RESULT(ZDVA)
!
USE MODD_CSTS,    ONLY : XTT, XP00  
!  
USE YOMHOOK      ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1     ,ONLY : JPRB
!
IMPLICIT NONE

REAL, INTENT(IN) :: PT
REAL, INTENT(IN) :: PPS
!
REAL             :: ZWORK, ZDVA
!
REAL(KIND=JPRB)  :: ZHOOK_HANDLE
!
REAL, PARAMETER                :: ZNV    = 1.88         ! (-)
REAL, PARAMETER                :: ZCV    = 2.17e-5      ! (m2/s)
!
IF (LHOOK) CALL DR_HOOK('MODE_HYDRO_DIF:MDIFFUSIVITY_WATER_0D',0,ZHOOK_HANDLE)
!
ZWORK  = ZNV*LOG(PT/XTT)
ZDVA   = ZCV*(XP00/PPS)*EXP(ZWORK)
!
IF (LHOOK) CALL DR_HOOK('MODE_HYDRO_DIF:MDIFFUSIVITY_WATER_0D',1,ZHOOK_HANDLE)

END FUNCTION MDIFFUSIVITY_WATER_0D
!
!-------------------------------------------------------------------------------
!Soil resistance for baresoil evaporation: Sellers et al (1992)
!-------------------------------------------------------------------------------
!
FUNCTION SOIL_RES_EG_S92(PWG,PWSAT) RESULT(ZRSOIL)
USE MODD_SURF_PAR, ONLY : XUNDEF
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
IMPLICIT NONE
REAL,    DIMENSION(:), INTENT(IN)  :: PWG    ! surface soil moisture content (m3/m3)
REAL,    DIMENSION(:), INTENT(IN)  :: PWSAT  ! saturation volumetric water content (m3/m3)
REAL,    DIMENSION(SIZE(PWSAT))    :: ZRSOIL ! s/m
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
! Parameterization based on Sellers et al., 1992
!
REAL, PARAMETER            :: ZRG_COEF1    = 8.206  ! Ground/litter resistance coefficient 
REAL, PARAMETER            :: ZRG_COEF2    = 4.255  ! Ground/litter resistance coefficient 
!
IF (LHOOK) CALL DR_HOOK('MODE_HYDRO_DIF:SOIL_RES_EG_S92',0,ZHOOK_HANDLE)
!
ZRSOIL(:) = EXP(ZRG_COEF1 - ZRG_COEF2 * PWG(:) / PWSAT(:))
!
IF (LHOOK) CALL DR_HOOK('MODE_HYDRO_DIF:SOIL_RES_EG_S92',1,ZHOOK_HANDLE)
!
END FUNCTION SOIL_RES_EG_S92
!-------------------------------------------------------------------------------
!Soil resistance for baresoil evaporation: DSL (Swenson and Lawrence, 2014)
!-------------------------------------------------------------------------------
!
FUNCTION SOIL_RES_EG_DSL(PTG,PPS,PWG,PWGI,PWSAT,PBCOEF,PMPOTSAT) RESULT(ZRSOIL)
USE MODD_SURF_PAR, ONLY : XUNDEF
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
IMPLICIT NONE
REAL,    DIMENSION(:), INTENT(IN)  :: PTG      ! soil temperature (K)
REAL,    DIMENSION(:), INTENT(IN)  :: PPS      ! surface pressure (Pa)
REAL,    DIMENSION(:), INTENT(IN)  :: PWG      ! surface soil moisture content (m3/m3)
REAL,    DIMENSION(:), INTENT(IN)  :: PWGI     ! surface soil frozen moisture content (m3/m3)
REAL,    DIMENSION(:), INTENT(IN)  :: PWSAT    ! saturation volumetric water content (m3/m3)
REAL,    DIMENSION(:), INTENT(IN)  :: PMPOTSAT ! saturation matric potential (m)
REAL,    DIMENSION(:), INTENT(IN)  :: PBCOEF   ! slope of soil water rentention curve (-)
REAL,    DIMENSION(SIZE(PWSAT))    :: ZRSOIL   ! s/m
!
REAL,    DIMENSION(SIZE(PWSAT))    :: ZWSAT, ZW_AIR, ZDSL, ZDVA, ZPHI_AIR, ZW_DSL0, ZB_INV
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE

!*      0.3    declarations of local parameters

REAL, PARAMETER            :: ZPSI_AIR = -10000.0      ! DSL: matric potential in m of air (m)
REAL, PARAMETER            :: ZKP_DSL  =      0.8      ! DSL: (-) 
REAL, PARAMETER            :: ZDELTAZ  =      0.015    ! DSL: thickness of layer for resistance computation
!
IF (LHOOK) CALL DR_HOOK('MODE_HYDRO_DIF:SOIL_RES_EG_DSL',0,ZHOOK_HANDLE)
!
! init:
!
ZB_INV(:)   = 0.
ZW_AIR(:)   = 0.
ZDSL(:)     = 0.
ZDVA(:)     = 0.
ZPHI_AIR(:) = 0.
ZRSOIL(:)   = 0. 

! Ice-modified satuation volumetric water content (VWC) (m3/m3)

ZWSAT(:)   = PWSAT(:)-PWGI(:)

ZW_DSL0(:) = ZKP_DSL*PWSAT(:) !threshold VWC for the formation of a DSL

WHERE( PWG(:) < ZW_DSL0(:)) ! Only where sufficiently dry (below saturation)

   ZB_INV(:)   = 1/PBCOEF(:)
   ZW_AIR(:)   = PWSAT(:)*(PMPOTSAT(:)/ZPSI_AIR)**ZB_INV(:)
   ZDSL(:)     = ZDELTAZ*(ZW_DSL0(:)-PWG(:))/(ZW_DSL0(:)-ZW_AIR(:))

! molecular diffusivity of water vapor (m2 s-1)
! from the article of Swenson and Lawrence (2014)

   ZDVA(:)     = MDIFFUSIVITY_WATER(PTG(:),PPS(:))

   ZPHI_AIR(:) = PWSAT(:)-ZW_AIR(:)  ! matric potential depression (relative to air)
   
   ZRSOIL(:)   = ZDSL(:)/(ZDVA(:)*                                              &
                 ZPHI_AIR(:)*ZPHI_AIR(:)*(ZPHI_AIR(:)/PWSAT(:))**(3*ZB_INV(:))) ! tau

END WHERE

IF (LHOOK) CALL DR_HOOK('MODE_HYDRO_DIF:SOIL_RES_EG_DSL',1,ZHOOK_HANDLE)

END FUNCTION SOIL_RES_EG_DSL
!-------------------------------------------------------------------------------
!
END MODULE MODE_HYDRO_DIF
