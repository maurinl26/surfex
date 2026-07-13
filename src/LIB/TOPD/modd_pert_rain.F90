!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!     ###########################
      MODULE MODD_PERT_RAIN
!     ###########################
!
!!****  *MODD_PERT_RAIN
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
INTEGER                               :: NNB_STEPS_MODIF ! Number of timesteps to modify
INTEGER                               :: NNB_MEMBERS ! Final number of members
INTEGER                               :: NNB_MEMBERS_LOC ! Number of members with modified location
INTEGER                               :: NNB_MAX_MEMBERS ! Maximal number of members
REAL,    DIMENSION(:,:,:),ALLOCATABLE :: XRAIN_NEW
REAL,    DIMENSION(:,:,:),ALLOCATABLE :: XRAIN_NEWSEL
REAL,    DIMENSION(:),    ALLOCATABLE :: XPROB_DECAL
INTEGER , DIMENSION(:),   ALLOCATABLE ::NNB_MEMBERS_SEL ! Number of selected members
INTEGER , DIMENSION(:),   ALLOCATABLE ::NMEMBER_SEL ! Members where convective objects were changed
!-------------------------------------------------------------------------------------
!
END MODULE MODD_PERT_RAIN

