!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!     #######################
      SUBROUTINE FIND_POINTS(PRR1,PTHRESHOLD,KI,KJ,OINAREA)
!     #######################
!
!!****  *FIND_POINTS*  
!!
!!    PURPOSE
!!    -------
!!    This routine aims at delineating areas corresponding to a thresold exeedance 
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
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
REAL, DIMENSION(:)                 , INTENT(IN)  :: PRR1     ! Imput hourly field
REAL                               , INTENT(IN)  :: PTHRESHOLD ! Threshold from which points are in an area
INTEGER                            , INTENT(IN)  :: KI
INTEGER                            , INTENT(IN)  :: KJ
LOGICAL, DIMENSION(:,:),ALLOCATABLE, INTENT(OUT) :: OINAREA ! LOGICAL to delineate areas (TRUE if in interresting aeras)
!
!
!*      0.2    declarations of local variables
!
!
INTEGER                              :: JWRK, JWRK1, JWRK2
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('FIND_POINTS',0,ZHOOK_HANDLE)
!
IF (.NOT.ALLOCATED(OINAREA)) ALLOCATE (OINAREA(KI,KJ))
OINAREA(:,:)=.FALSE.
DO JWRK2=1,KJ
  DO JWRK1=1,KI
    JWRK=(JWRK2-1)*KI+JWRK1
    IF(PRR1(JWRK)>=PTHRESHOLD) OINAREA(JWRK1,JWRK2)=.TRUE.
  ENDDO
ENDDO
!
IF (LHOOK) CALL DR_HOOK('FIND_POINTS',1,ZHOOK_HANDLE)
!
END SUBROUTINE FIND_POINTS
