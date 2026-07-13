!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!     ######################################
      SUBROUTINE PGD_GRIDTYPE_INIT(HGRID,KGRID_PAR,PGRID_PAR)
!     ######################################
!!
!!    PURPOSE
!!    -------
!!
!!
!!    METHOD
!!    ------
!!   
!!    EXTERNAL
!!    --------
!!
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!
!!    REFERENCE
!!    ---------
!!
!!    AUTHOR
!!    ------
!!
!!    V. Masson                   Meteo-France
!!
!!    MODIFICATION
!!    ------------
!!
!!    Original     13/10/03
!----------------------------------------------------------------------------
!
!*    0.     DECLARATION
!            -----------
!
USE MODD_PGD_GRID,   ONLY : CGRID, XGRID_PAR, NGRID_PAR
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*    0.1    Declaration of dummy arguments
!            ------------------------------
!
CHARACTER(LEN=10),           INTENT(IN)  :: HGRID       ! type of grid
INTEGER                                  :: KGRID_PAR   ! size of PGRID_PAR
REAL,    DIMENSION(:),       POINTER     :: PGRID_PAR   ! parameters defining this grid
!
!
!*    0.2    Declaration of local variables
!            ------------------------------
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('PGD_GRIDTYPE_INIT',0,ZHOOK_HANDLE)
!
ALLOCATE(XGRID_PAR(KGRID_PAR))
XGRID_PAR=PGRID_PAR
CGRID=HGRID
NGRID_PAR=KGRID_PAR
!
IF (LHOOK) CALL DR_HOOK('PGD_GRIDTYPE_INIT',1,ZHOOK_HANDLE)
!-------------------------------------------------------------------------------
!
END SUBROUTINE PGD_GRIDTYPE_INIT
