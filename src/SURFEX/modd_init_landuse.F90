!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!     #####################
      MODULE MODD_INIT_LANDUSE
!     #####################
!
!!****  *MODD_INIT_LANDUSE* - declaration of landuse types used for landuse initialisation
!!
!!    PURPOSE
!!    -------
!       The purpose of this declarative module is to define landuse types used for landuse initialisation
!
!!
!!**  IMPLICIT ARGUMENTS
!!    ------------------
!!      NONE 
!!
!!    REFERENCE
!!    --------- 
!!       
!!    AUTHOR
!!    ------
!!      B. Decharme   *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    20/12/2021                      
!-------------------------------------------------------------------------------
!
!*       0.   DECLARATIONS
!             ------------
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
TYPE LULCC_t
!
!* Carbon due to biomass harvested after land use change by patch including patch fraction (kgC/m2)
!
  REAL, DIMENSION(:,:), POINTER :: XLULCC_HARVEST
!
! Carbon due to biomass harvested after land use change including patch fraction (kgC/m2)
!
  REAL, DIMENSION(:),   POINTER :: XLULCC_HARVEST_GRID
!
!* To compute carbon budget after LULCC
!
  REAL, DIMENSION(:),   POINTER :: XBIOM_GRID_OLD 
  REAL, DIMENSION(:),   POINTER :: XBIOM_GRID_NEW 
  REAL, DIMENSION(:),   POINTER :: XLITTER_GRID_OLD
  REAL, DIMENSION(:),   POINTER :: XLITTER_GRID_NEW
  REAL, DIMENSION(:),   POINTER :: XCSOIL_GRID_OLD
  REAL, DIMENSION(:),   POINTER :: XCSOIL_GRID_NEW
!
END TYPE LULCC_t
!
!-------------------------------------------------------------------------------
!
CONTAINS
!
!-------------------------------------------------------------------------------
!
SUBROUTINE INIT_LULCC(YLULCC)
!
TYPE(LULCC_t), INTENT(INOUT) :: YLULCC
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
IF (LHOOK) CALL DR_HOOK("MODD_INIT_LANDUSE:INIT_LULCC_P",0,ZHOOK_HANDLE)
!
NULLIFY(YLULCC%XLULCC_HARVEST)
NULLIFY(YLULCC%XLULCC_HARVEST_GRID)
NULLIFY(YLULCC%XBIOM_GRID_OLD)
NULLIFY(YLULCC%XBIOM_GRID_NEW)
NULLIFY(YLULCC%XLITTER_GRID_OLD)
NULLIFY(YLULCC%XLITTER_GRID_NEW)
NULLIFY(YLULCC%XCSOIL_GRID_OLD)
NULLIFY(YLULCC%XCSOIL_GRID_NEW)
!
IF (LHOOK) CALL DR_HOOK("MODD_INIT_LANDUSE:INIT_LULCC_P",1,ZHOOK_HANDLE)
!
END SUBROUTINE INIT_LULCC
!
!-------------------------------------------------------------------------------
!
SUBROUTINE ALLOC_LULCC(YLULCC,KI,KP)
!
TYPE(LULCC_t), INTENT(INOUT) :: YLULCC
!
INTEGER, INTENT(IN) :: KI
INTEGER, INTENT(IN) :: KP
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
IF (LHOOK) CALL DR_HOOK("MODD_INIT_LANDUSE:ALLOC_LULCC_P",0,ZHOOK_HANDLE)
!
ALLOCATE(YLULCC%XLULCC_HARVEST(KI,KP))
YLULCC%XLULCC_HARVEST(:,:) = 0.0
!
ALLOCATE(YLULCC%XLULCC_HARVEST_GRID(KI))
YLULCC%XLULCC_HARVEST_GRID(:) = 0.0
!
ALLOCATE(YLULCC%XBIOM_GRID_OLD  (KI))
ALLOCATE(YLULCC%XBIOM_GRID_NEW  (KI))
ALLOCATE(YLULCC%XLITTER_GRID_OLD(KI))
ALLOCATE(YLULCC%XLITTER_GRID_NEW(KI))
ALLOCATE(YLULCC%XCSOIL_GRID_OLD (KI))
ALLOCATE(YLULCC%XCSOIL_GRID_NEW (KI))
YLULCC%XBIOM_GRID_OLD  (:) = 0.0
YLULCC%XBIOM_GRID_NEW  (:) = 0.0
YLULCC%XLITTER_GRID_OLD(:) = 0.0
YLULCC%XLITTER_GRID_NEW(:) = 0.0
YLULCC%XCSOIL_GRID_OLD (:) = 0.0
YLULCC%XCSOIL_GRID_NEW (:) = 0.0
!
IF (LHOOK) CALL DR_HOOK("MODD_INIT_LANDUSE:ALLOC_LULCC_P",1,ZHOOK_HANDLE)
!
END SUBROUTINE ALLOC_LULCC
!
!-------------------------------------------------------------------------------
!
END MODULE MODD_INIT_LANDUSE

