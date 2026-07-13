!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!     ########################
      SUBROUTINE GAUSS_1REAL(PAVG,PSTDDEV,PRAIN,PDENSITY_FCT)
!     ########################
!
!!*****    * GAUSS *
!
!!    PURPOSE
!!    -------
!   This routine gives for a set of values of rain, the cooresponding density
!   values given a value for the shape and scale parameters.
!     
!         
!     
!!**  METHOD
!!    ------
!
!!    EXTERNAL
!!    --------
!!
!!    none
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------ 
!!
!!    
!!    
!!
!!      
!!    REFERENCE
!!    ---------
!!
!!    
!!      
!!    AUTHOR
!!    ------
!!
!!     Vincendon
!!
!!    MODIFICATIONS
!!    -------------
!!
!!      Original   12/2008
!----------------------------------------------------------------------
!*       0.      DECLARATIONS
!                ------------
!
USE MODD_CSTS, ONLY : XPI
!
USE YOMHOOK   ,ONLY : LHOOK, DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
REAL,INTENT(IN)   :: PAVG ! average
REAL,INTENT(IN)   :: PSTDDEV! standard deviation
REAL, INTENT(IN)  :: PRAIN
REAL, INTENT(OUT) :: PDENSITY_FCT
!*      0.2    declarations of local variables
! 
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!-----------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('GAUSS_1REAL',0,ZHOOK_HANDLE)
!
!*       0.     Initialization:
! 
PDENSITY_FCT=0.

IF (PSTDDEV/=0.) PDENSITY_FCT=EXP(-0.5*((PRAIN-PAVG)/PSTDDEV)**2.)/(PSTDDEV*(2.*XPI)**(0.5))
!
IF (LHOOK) CALL DR_HOOK('GAUSS_1REAL',1,ZHOOK_HANDLE)
!
END SUBROUTINE GAUSS_1REAL
