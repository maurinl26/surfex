!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!     #########
    SUBROUTINE LANDUSE_CARBON_MANAGING(IO, S, SOLD, TLU, KI)       
!   ###############################################################
!!
!!    PURPOSE
!!    -------
!
!     Performs land use land cover change managing computation at yearly time step
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
!!    R. Séférian 08/2015
!!    B. Decharme 08/2021
!!
!-------------------------------------------------------------------------------
!
!*       0.     DECLARATIONS
!               ------------
!
USE MODD_ISBA_n,         ONLY : ISBA_S_t
USE MODD_ISBA_OPTIONS_n, ONLY : ISBA_OPTIONS_t
USE MODD_INIT_LANDUSE,   ONLY : LULCC_t
!
USE MODD_DATA_COVER_PAR, ONLY : NVEGTYPE
USE MODD_LANDUSE_PAR
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
TYPE(ISBA_OPTIONS_t),  INTENT(INOUT) :: IO
TYPE(ISBA_S_t),        INTENT(INOUT) :: S, SOLD
TYPE(LULCC_t),         INTENT(INOUT) :: TLU
!
INTEGER,               INTENT(IN)    :: KI
!
!*      0.2    declarations of local parameter
!
REAL, DIMENSION(KI,IO%NPATCH)                             :: ZFLUATM_ANNUAL          
REAL, DIMENSION(KI,IO%NPATCH)                             :: ZFLUATM_DECADAL          
REAL, DIMENSION(KI,IO%NPATCH)                             :: ZFLUATM_CENTURY         
REAL, DIMENSION(KI,IO%NPATCH)                             :: ZFLUANT          
REAL, DIMENSION(KI,IO%NPATCH)                             :: ZEXPORT_COEF_ANNUAL
REAL, DIMENSION(KI,IO%NPATCH)                             :: ZEXPORT_COEF_DECADAL
REAL, DIMENSION(KI,IO%NPATCH)                             :: ZEXPORT_COEF_CENTURY
!
REAL, DIMENSION(KI,IO%NNDECADAL+1,IO%NPATCH)   :: ZCSTOCK_DECADAL 
REAL, DIMENSION(KI,IO%NNCENTURY+1,IO%NPATCH)   :: ZCSTOCK_CENTURY
!
INTEGER :: JI, JIND, JL, JP, JYEAR, JVEG, IDECA2, ICENT2
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!-----------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('LANDUSE_CARBON_MANAGING',0,ZHOOK_HANDLE)
!-----------------------------------------------------------------
!
!*      1.     Preliminaries
!              -------------
!
IDECA2=IO%NNDECADAL+1
ICENT2=IO%NNCENTURY+1
!
! 1.1 Initialize matrix 
!
ZFLUATM_ANNUAL (:,:) = 0.
ZFLUATM_DECADAL(:,:) = 0.
ZFLUATM_CENTURY(:,:) = 0.
ZFLUANT        (:,:) = 0.
!
! 1.2 Initialize stock matrix
! Current year (year=1) anthropogenic stock of carbon is set to zero
! other years (2...N) conserve the previous stock of anth carbon
!
ZCSTOCK_DECADAL(:,:,:) = 0.
ZCSTOCK_DECADAL(:,2:IDECA2,:) = S%XCSTOCK_DECADAL(:,:,:)
!
ZCSTOCK_CENTURY(:,:,:) = 0.
ZCSTOCK_CENTURY(:,2:ICENT2,:) = S%XCSTOCK_CENTURY(:,:,:)
!
! 1.3 Characterized biomass export properties for each patch
!
ZEXPORT_COEF_ANNUAL (:,:)   = 0.
ZEXPORT_COEF_DECADAL(:,:)   = 0.
ZEXPORT_COEF_CENTURY(:,:)   = 0.
!
DO JP=1,IO%NPATCH
   DO JVEG=1,NVEGTYPE
      DO JI=1,KI
         ZEXPORT_COEF_ANNUAL(JI,JP)  = ZEXPORT_COEF_ANNUAL(JI,JP)  &
                                     + XEXPORT_COEF_ANNUAL(JVEG) * SOLD%XVEGTYPE_PATCH(JI,JVEG,JP)
         ZEXPORT_COEF_DECADAL(JI,JP) = ZEXPORT_COEF_DECADAL(JI,JP) &
                                     + XEXPORT_COEF_DECADAL(JVEG)* SOLD%XVEGTYPE_PATCH(JI,JVEG,JP)
         ZEXPORT_COEF_CENTURY(JI,JP) = ZEXPORT_COEF_CENTURY(JI,JP) &
                                     + XEXPORT_COEF_CENTURY(JVEG)* SOLD%XVEGTYPE_PATCH(JI,JVEG,JP)
      ENDDO
   ENDDO
ENDDO
!
!-----------------------------------------------------------------
!
!*      3.     LULCC Carbon Fluxes
!              -------------------
!
!* Annual export of carbon
ZFLUATM_ANNUAL(:,:) = ZEXPORT_COEF_ANNUAL(:,:) * TLU%XLULCC_HARVEST(:,:)
!
!* Decadal export of carbon
ZCSTOCK_DECADAL(:,1,:) = ZEXPORT_COEF_DECADAL(:,:) * TLU%XLULCC_HARVEST(:,:)
!
!* Carbon flux to the decadal anthropogenic C pool
ZFLUANT(:,:) = ZFLUANT(:,:) + ZEXPORT_COEF_DECADAL(:,:) * TLU%XLULCC_HARVEST(:,:)
!
!temporal decay
DO JP=1,IO%NPATCH
   DO JI=1,KI
      DO JIND=1,IO%NNDECADAL-1
         JYEAR = IO%NNDECADAL - JIND + 1
         ZFLUATM_DECADAL   (JI,      JP) = ZFLUATM_DECADAL  (JI,        JP) + S%XEXPORT_DECADAL(JI,JYEAR  ,JP)
         ZCSTOCK_DECADAL   (JI,JYEAR,JP) = ZCSTOCK_DECADAL  (JI,JYEAR-1,JP) - S%XEXPORT_DECADAL(JI,JYEAR-1,JP)
         S%XEXPORT_DECADAL (JI,JYEAR,JP) = S%XEXPORT_DECADAL(JI,JYEAR-1,JP)
         IF(ZCSTOCK_DECADAL(JI,JYEAR,JP) < 0.)THEN
            ZCSTOCK_DECADAL(JI,JYEAR,JP) = 0.
         ENDIF
      ENDDO
   ENDDO
ENDDO
!
!update of the stock  
ZFLUATM_DECADAL    (:,:) = ZFLUATM_DECADAL(:,:) + S%XEXPORT_DECADAL(:,2,:)  
S%XEXPORT_DECADAL(:,2,:) = ZCSTOCK_DECADAL(:,1,:)/REAL(IO%NNDECADAL)
ZCSTOCK_DECADAL  (:,2,:) = ZCSTOCK_DECADAL(:,1,:)
!
!
!* Centennial export of carbon
ZCSTOCK_CENTURY(:,1,:) = ZCSTOCK_CENTURY(:,1,:) + ZEXPORT_COEF_CENTURY(:,:) * TLU%XLULCC_HARVEST(:,:)
!
!* Carbon flux to the centennial anthropogenic C pool
ZFLUANT(:,:) = ZFLUANT(:,:) + ZEXPORT_COEF_CENTURY(:,:) * TLU%XLULCC_HARVEST(:,:)
!
!  temporal decay
DO JP=1,IO%NPATCH
   DO JI=1,KI
      DO JIND=1,IO%NNCENTURY-1
         JYEAR = IO%NNCENTURY - JIND + 1
         ZFLUATM_CENTURY   (JI,      JP) = ZFLUATM_CENTURY  (JI,        JP) + S%XEXPORT_CENTURY(JI,JYEAR  ,JP)
         ZCSTOCK_CENTURY   (JI,JYEAR,JP) = ZCSTOCK_CENTURY  (JI,JYEAR-1,JP) - S%XEXPORT_CENTURY(JI,JYEAR-1,JP)
         S%XEXPORT_CENTURY (JI,JYEAR,JP) = S%XEXPORT_CENTURY(JI,JYEAR-1,JP)
         IF(ZCSTOCK_CENTURY(JI,JYEAR,JP) < 0.)THEN
            ZCSTOCK_CENTURY(JI,JYEAR,JP) = 0.
         ENDIF
      ENDDO
   ENDDO
ENDDO
!
!update of the stock  
ZFLUATM_CENTURY    (:,:) = ZFLUATM_CENTURY(:,:) + S%XEXPORT_CENTURY(:,2,:)  
S%XEXPORT_CENTURY(:,2,:) = ZCSTOCK_CENTURY(:,1,:)/REAL(IO%NNCENTURY)
ZCSTOCK_CENTURY  (:,2,:) = ZCSTOCK_CENTURY(:,1,:)
!
!reset current year reservoir
ZCSTOCK_DECADAL(:,1,:) = 0.
ZCSTOCK_CENTURY(:,1,:) = 0.
!
!* Finalize
! 
! Carbon flux due to vegetation clearance
S%XFLUATM(:,:)  =  ZFLUATM_ANNUAL(:,:) 
!
! Carbon flux from the anthropogenic carbon pool
S%XFANTATM(:,:) = (ZFLUATM_DECADAL(:,:)+ ZFLUATM_CENTURY(:,:)) 
!
! Carbon flux from natural carbon stock to anthropogenic carbon pool
S%XFLUANT(:,:)  =  ZFLUANT(:,:)
!
! Conserve history of stock for the following years
S%XCSTOCK_DECADAL(:,:,:) = ZCSTOCK_DECADAL(:,2:IDECA2,:) 
S%XCSTOCK_CENTURY(:,:,:) = ZCSTOCK_CENTURY(:,2:ICENT2,:) 
!
IF (LHOOK) CALL DR_HOOK('LANDUSE_CARBON_MANAGING',1,ZHOOK_HANDLE)
!-----------------------------------------------------------------
!
END SUBROUTINE LANDUSE_CARBON_MANAGING
