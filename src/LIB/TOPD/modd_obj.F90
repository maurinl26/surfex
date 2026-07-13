!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!     ###########################
      MODULE MODD_OBJ
!     ###########################
!
!!****  *MODD_OBJ
!!
!!    PURPOSE
!!    -------
!
!!
!!**  IMPLICIT ARGUMENTS
!!    ------------------
!!      None 
!!
!!    REFERENCE
!!    ---------
!!
!!    AUTHOR
!!    ------
!!    B. Vincendon
!!
!!    MODIFICATIONS
!!    -------------
!!      Original       11/2007
!
!*       0.   DECLARATIONS
!             ------------
!
!
IMPLICIT NONE
 
!-------------------------------------------------------------------------------
! For each grid-point, each time step, index of the objet 
INTEGER, DIMENSION(:,:), ALLOCATABLE :: NN_REG_RAD ! for radar fields
INTEGER, DIMENSION(:,:), ALLOCATABLE :: NN_REG_ARO ! for arome fields
! For each object, X and Y coordinates of the barycenter
INTEGER, DIMENSION(:,:), ALLOCATABLE :: NX_BARY_RAD,NY_BARY_RAD ! for radar fields
INTEGER, DIMENSION(:,:), ALLOCATABLE :: NX_BARY_ARO,NY_BARY_ARO ! for arome fields
INTEGER, DIMENSION(:),   ALLOCATABLE :: NX_BARY_MAX_RAD,NY_BARY_MAX_RAD!barycentres des objets radar plus étendus
INTEGER, DIMENSION(:),   ALLOCATABLE :: NX_BARY_MAX_ARO,NY_BARY_MAX_ARO!barycentres des objets arome plus étendus
INTEGER, DIMENSION(:),   ALLOCATABLE :: NREGMAX_ARO,NREGMAX_RAD ! Numero of the biggest region
INTEGER, DIMENSION(:),   ALLOCATABLE :: NMAX_ARO,NMAX_RAD ! number of pixels in the biggest region
INTEGER, DIMENSION(:),   ALLOCATABLE :: NTOT_REG_ARO,NTOT_REG_RAD ! Total Number of regions
REAL,    DIMENSION(:),   ALLOCATABLE :: XAVGMAX_ARO,XAVGMAX_RAD ! average rainfall in the biggest region
LOGICAL, DIMENSION(:),   ALLOCATABLE :: LRADOK !true if mesh with the radar "portée"
!***
INTEGER :: NTEST
!
!***
LOGICAL, DIMENSION(:,:),  ALLOCATABLE :: LINAREA ! LOGICAL to delineate areas (TRUE if in interresting aeras)
INTEGER, DIMENSION(:,:),  ALLOCATABLE :: NNUMREG  ! Numero of  the region for each grid-point
!-------------------------------------------------------------------------------------
!
END MODULE MODD_OBJ

