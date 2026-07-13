!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!-----------------------------------------------------------------
!     ########################
      SUBROUTINE SAT_AREA_FRAC(PDEF,PAS,PAS_MESHBV,GTOPD)
!     ########################
!
!!*****    * SAT_AREA_FRAC *
!
!!    PURPOSE
!!    -------
!    
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
!!      K. Chancibault  * LTHE / Meteo-France *
!!
!!    MODIFICATIONS
!!    -------------
!!
!!      Original   27/11/2006
!!                 03/2014 (B. Vincendon) computation based of pixels counts instead of areas 
!!                 07/2017 (B. Vincendon) adding a variable of saturation by mesh and catchment
!
!----------------------------------------------------------------------
!*       0.      DECLARATIONS
!                ------------
!
USE MODD_TOPODYN,       ONLY : NNCAT, NNMC
USE MODD_COUPLING_TOPD, ONLY : NMASKT, NNPIX, NNBV_IN_MESH
USE MODD_SURF_PAR,        ONLY : XUNDEF, NUNDEF
!
USE YOMHOOK   ,ONLY : LHOOK, DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
REAL,    DIMENSION(:,:),INTENT(IN)     :: PDEF    ! deficit
REAL,    DIMENSION(:)  , INTENT(OUT)   :: PAS     !contributive area fraction in Isba meshes
REAL,    DIMENSION(:,:), INTENT(OUT)   :: PAS_MESHBV     !contributive area fraction in Isba meshes by catchment
LOGICAL, DIMENSION(:)  , INTENT(INOUT) :: GTOPD     ! 
!
!*      0.2    declarations of local variables
INTEGER               :: JJ, JI,JMESH
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-----------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('SAT_AREA_FRAC',0,ZHOOK_HANDLE)
!
! full grid only
!*       0.     Initialization:
! 
PAS(:)=0.0
PAS_MESHBV(:,:)=0.0
!
DO JJ=1,NNCAT
  DO JI=1,NNMC(JJ)
    IF (PDEF(JJ,JI).LE.0.001 .AND. NMASKT(JJ,JI)/=NUNDEF .AND. NMASKT(JJ,JI)/=0) THEN 
      PAS(NMASKT(JJ,JI)) = PAS(NMASKT(JJ,JI)) +1. 
      PAS_MESHBV(NMASKT(JJ,JI),JJ) = PAS_MESHBV(NMASKT(JJ,JI),JJ) +1. 
    ENDIF
  ENDDO
ENDDO
!
! Calculation of the saturated area ratio in each Isba mesh
DO JMESH=1,SIZE(NNBV_IN_MESH,1)
  !
  IF(NNPIX(JMESH)/=0.AND.NNPIX(JMESH)/=XUNDEF.AND.PAS(JMESH)/=XUNDEF)THEN
    PAS(JMESH) = PAS(JMESH) /REAL(NNPIX(JMESH))
  ELSE
    PAS(JMESH) = XUNDEF
  ENDIF
  !
  DO JJ=1,NNCAT
    IF(NNPIX(JMESH)/=0.AND.NNPIX(JMESH)/=XUNDEF)THEN
      PAS_MESHBV(JMESH,JJ)= PAS_MESHBV(JMESH,JJ)/REAL(NNPIX(JMESH)) !full grid
    ELSE
      PAS_MESHBV(JMESH,JJ)= XUNDEF
    ENDIF
  ENDDO
  !
ENDDO
!
IF (LHOOK) CALL DR_HOOK('SAT_AREA_FRAC',1,ZHOOK_HANDLE)
!
END SUBROUTINE SAT_AREA_FRAC
