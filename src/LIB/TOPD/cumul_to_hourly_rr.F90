!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!     #######################
      SUBROUTINE CUMUL_TO_HOURLY_RR(PRRC,PRR1,KNB_STEPS)
!     #######################
!
!!****  *CUMUL_TO_HOURLY_RR*  
!!
!!    PURPOSE
!!    -------
!     This routine aims at computing hourly rainfall from accumulated ones and
!     to convert mm in kg/m2/s.
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
!!      B. Vincendon    * Meteo-France *
!!
!!    MODIFICATIONS
!!    -------------
!!
!!      Original   11/2007
!-------------------------------------------------------------------------------
!
!*       0.     DECLARATIONS
!               ------------
USE MODD_SURF_PAR,   ONLY : XUNDEF
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
REAL, DIMENSION(:,:) ,INTENT(IN)    :: PRRC     ! Input cumulated field
REAL, DIMENSION(:,:) ,INTENT(INOUT) :: PRR1     ! Output hourly field
INTEGER              ,INTENT(IN)    ::KNB_STEPS ! Total number of time steps
!
!*      0.2    declarations of local variables
!
INTEGER :: JWRK, JWRK2
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('CUMUL_TO_HOURLY_RR',0,ZHOOK_HANDLE)
!
WHERE (PRRC(1,:)/=XUNDEF)
  PRR1(1,:)=PRRC(1,:)/3600.
ENDWHERE
!
DO JWRK=2,KNB_STEPS
  WHERE (PRRC(JWRK,:)/=XUNDEF)
    PRR1(JWRK,:)=(PRRC(JWRK,:)-PRRC(JWRK-1,:))/3600.
  ENDWHERE
ENDDO
!
IF (LHOOK) CALL DR_HOOK('CUMUL_TO_HOURLY_RR',1,ZHOOK_HANDLE)
!
END SUBROUTINE CUMUL_TO_HOURLY_RR
