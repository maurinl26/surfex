!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!-----------------------------------------------------------------
!     #####################
      SUBROUTINE CONTROL_WATER_BUDGET_TOPD (IO,S, U, PWGM,PWG,PDG,PMESH_SIZE,PAVG_MESH_SIZE,PWSAT,PWOVSAT_IBV)
!     #####################
!
!!****  *CONTROL_WATER_BUDGET_TOPD*  
!!
!!    PURPOSE
!!    -------
!     To control water budget after topodyn_lat lateral distribution
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
!!      B. Vincendon *  Meteo-France *
!!
!!    MODIFICATIONS
!!    -------------
!!
!!      Original   : Out of COUPL_TOPD in february 2014
!!      Modif 07/2017 : water over saturation saved for runoff computation
!-------------------------------------------------------------------------------
!
!*       0.     DECLARATIONS
!               ------------
!
!
USE MODD_ISBA_OPTIONS_n, ONLY : ISBA_OPTIONS_t
USE MODD_ISBA_n, ONLY : ISBA_S_t
USE MODD_SURF_ATM_n, ONLY : SURF_ATM_t
!
USE MODD_SURF_PAR,         ONLY : XUNDEF, NUNDEF
USE MODD_COUPLING_TOPD,    ONLY : XTOTBV_IN_MESH, NNBV_IN_MESH
USE MODD_TOPODYN,          ONLY : NNCAT
USE MODD_ISBA_PAR,         ONLY : XWGMIN
USE MODI_AVG_PATCH_WG
!
USE MODI_PACK_SAME_RANK
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
!
!
TYPE(ISBA_OPTIONS_t), INTENT(INOUT) :: IO
TYPE(ISBA_S_t), INTENT(INOUT) :: S
TYPE(SURF_ATM_t), INTENT(INOUT) :: U
!
REAL, DIMENSION(:,:), INTENT(IN)    :: PWGM
REAL, DIMENSION(:,:), INTENT(INOUT) :: PWG
REAL, DIMENSION(:,:), INTENT(IN)    :: PDG
REAL, DIMENSION(:),   INTENT(IN)    :: PMESH_SIZE
REAL,                 INTENT(IN)    :: PAVG_MESH_SIZE
REAL, DIMENSION(:),   INTENT(IN)    :: PWSAT
REAL, DIMENSION(:,:), INTENT(INOUT) :: PWOVSAT_IBV
!
!
!*      0.2    declarations of local variables
!
!
REAL, DIMENSION(SIZE(PWG,1),3)     :: ZWG_3L, ZWGI_3L, ZDG_3L          
REAL                               :: ZSTOCK_WGM, ZSTOCK_WG
REAL                               :: ZAVG_DGALL, ZCONTROL_WATER_BUDGET_TOPD
REAL                               :: ZTMP, ZTMP2
INTEGER                            :: JMESH, JP, JJ, JCAT
REAL,    DIMENSION(U%NDIM_NATURE) :: ZSUMPATCH
REAL,    DIMENSION(U%NDIM_NATURE) :: ZWG_CORR, ZAVG_WGM, ZAVG_WG, ZAVG_DG
REAL,    DIMENSION(U%NDIM_NATURE) :: ZTOTBV_IN_MESH
INTEGER, DIMENSION(U%NDIM_NATURE,NNCAT) :: INBV_IN_MESH
LOGICAL, DIMENSION(U%NDIM_NATURE) :: LMODIF
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('CONTROL_WATER_BUDGET_TOPD',0,ZHOOK_HANDLE)
!
ZTMP = 0.0
IF(IO%NPATCH/=1) THEN
 !
  ZSUMPATCH(:) = 0.0
  DO JP=1,IO%NPATCH
    DO JJ=1,U%NDIM_NATURE
      ZSUMPATCH(JJ) = ZSUMPATCH(JJ) + S%XPATCH(JJ,JP)
    ENDDO
  ENDDO
  ZAVG_WGM(:)  = 0.
  ZAVG_WG(:)  = 0.
  ZAVG_DG(:)  = 0.
  !
  DO JP=1,IO%NPATCH
    DO JJ=1,U%NDIM_NATURE     
      IF(ZSUMPATCH(JJ) > 0..AND.PWGM(JJ,JP)/=XUNDEF.AND.PWG(JJ,JP)/=XUNDEF.AND.PDG (JJ,JP)/=XUNDEF)THEN
        !
        ZAVG_WGM(JJ)  = ZAVG_WGM(JJ)  + S%XPATCH(JJ,JP) * PWGM(JJ,JP)  * PDG (JJ,JP) 
        ZAVG_WG(JJ)  = ZAVG_WG(JJ)  + S%XPATCH(JJ,JP) * PWG(JJ,JP)  * PDG (JJ,JP) 
        ZAVG_DG(JJ) = ZAVG_DG(JJ) + S%XPATCH(JJ,JP) * PDG (JJ,JP)
        !          
      ENDIF
    ENDDO
  ENDDO     
  !     
  WHERE (ZAVG_DG(:)>0.0.AND.ZSUMPATCH(:)>0.)
    ZAVG_WGM(:) = ZAVG_WGM(:) / ZAVG_DG(:)
    ZAVG_WG(:)  = ZAVG_WG(:)  / ZAVG_DG(:)
  ENDWHERE
  !
ELSE
  ZAVG_WGM(:)= PWGM(:,1)
  ZAVG_WG(:) = PWG(:,1) 
  ZAVG_DG(:) = PDG(:,1)
  ZSUMPATCH(:) = 1.0
ENDIF
!
ZSTOCK_WGM = SUM(ZAVG_WGM(:)*ZAVG_DG(:)*PMESH_SIZE(:),&
            MASK=(ZAVG_WGM(:)/=XUNDEF.AND.&
                  ZAVG_DG(:)/=XUNDEF.AND.&
                  PMESH_SIZE(:)/=XUNDEF.AND.&
                  ZSUMPATCH(:)>0.))    ! water stocked in the ground (m3)
!
ZSTOCK_WG  = SUM(ZAVG_WG(:)*ZAVG_DG(:)*PMESH_SIZE(:),&
            MASK=(ZAVG_WG(:)/=XUNDEF.AND.&
                  ZAVG_DG(:)/=XUNDEF.AND.&
                  PMESH_SIZE(:)/=XUNDEF.AND.&
                  ZSUMPATCH(:)>0.))    ! water stocked in the ground (m3)
!
IF ( COUNT(ZAVG_DG(:)/=XUNDEF.AND.ZSUMPATCH(:)>0.)/=0. )&
ZAVG_DGALL = SUM(ZAVG_DG(:),MASK=(ZAVG_DG(:)/=XUNDEF.AND.ZSUMPATCH(:)>0.))&
          / COUNT(ZAVG_DG(:)/=XUNDEF.AND.ZSUMPATCH(:)>0.)

IF (ZAVG_DGALL/=0.) THEN
  ZCONTROL_WATER_BUDGET_TOPD = ( ZSTOCK_WG - ZSTOCK_WGM )/ ZAVG_DGALL / PAVG_MESH_SIZE
77 CONTINUE
  !
  IF (ZCONTROL_WATER_BUDGET_TOPD==0.0) GOTO 66
  !
  ZTMP  = COUNT( ZAVG_WG(:)/=ZAVG_WGM(:).AND.ZAVG_WG(:)/=XUNDEF.AND.ZAVG_WGM(:)/=XUNDEF.AND.ZSUMPATCH(:)>0. )
  !
  LMODIF(:)=.FALSE.
  CALL PACK_SAME_RANK(U%NR_NATURE,XTOTBV_IN_MESH,ZTOTBV_IN_MESH)
  DO JCAT=1,NNCAT
    CALL PACK_SAME_RANK(U%NR_NATURE,NNBV_IN_MESH(:,JCAT),INBV_IN_MESH(:,JCAT))
  ENDDO
  IF (ZTMP/=0.) THEN
     !
     WHERE (ZTOTBV_IN_MESH(:)/=0.0.AND.ZAVG_WGM(:)/=XUNDEF.AND.ZAVG_WG(:)/=XUNDEF.AND.&
       ZAVG_WG(:)/=ZAVG_WGM(:) .AND. ZAVG_WG(:)>XWGMIN+(ZCONTROL_WATER_BUDGET_TOPD/ZTMP).AND.&
       ZAVG_WG(:)<=PWSAT(:)+(ZCONTROL_WATER_BUDGET_TOPD/ZTMP).AND.ZSUMPATCH(:)>0.)
       LMODIF(:)=.TRUE.
     ENDWHERE
     !
     WHERE (LMODIF)
       ZAVG_WG(:) = ZAVG_WG(:) - (ZCONTROL_WATER_BUDGET_TOPD/ZTMP)
     ENDWHERE
     ! 
  ENDIF
  !
  ZSTOCK_WG  = SUM(ZAVG_WG(:)*ZAVG_DG(:)*PMESH_SIZE(:),&
            MASK=(ZAVG_WG(:)/=XUNDEF.AND.&
                  ZAVG_DG(:)/=XUNDEF.AND.&
                  PMESH_SIZE(:)/=XUNDEF.AND.&
                  ZSUMPATCH(:)>0.))    ! water stocked in the ground (m3)

  !
  IF (ZAVG_DGALL/=0. .AND. PAVG_MESH_SIZE/=0.) THEN
    ZCONTROL_WATER_BUDGET_TOPD = ( ZSTOCK_WG - ZSTOCK_WGM )/ ZAVG_DGALL / PAVG_MESH_SIZE
  ENDIF
  !
  IF ((ABS(ZCONTROL_WATER_BUDGET_TOPD)>100.).AND.(COUNT(LMODIF)/=0)) GOTO 77
  !
ENDIF!ZAVG_DGALL/=0.
!
66 CONTINUE
! Adding excess in runoff
PWOVSAT_IBV(:,:)=0.
IF (1==1) THEN
  !
  DO JP=1,IO%NPATCH
    !
    DO JJ=1,U%NDIM_NATURE     
      !
      IF ((PWG(JJ,JP)/=XUNDEF).AND.(S%XPATCH(JJ,JP)>0.)&
         .AND.(S%XPATCH(JJ,JP)/=XUNDEF).AND.(ZTOTBV_IN_MESH(JJ)/=0.0))THEN
        !
        PWG(JJ,JP)=MAX(ZAVG_WG(JJ),XWGMIN)
        DO JCAT=1,NNCAT
          IF(INBV_IN_MESH(JJ,JCAT)/=0.AND.ZTMP/=0..AND.&
             ZAVG_WG(JJ)/=XUNDEF.AND.PWSAT(JJ)/=XUNDEF.AND.&
             ZCONTROL_WATER_BUDGET_TOPD/=XUNDEF)&
            PWOVSAT_IBV(JJ,JCAT)=PWOVSAT_IBV(JJ,JCAT)+ MAX(0.0,ZAVG_WG(JJ)-PWSAT(JJ)) +&
                          (ZCONTROL_WATER_BUDGET_TOPD/ZTMP)
        ENDDO
        !
      ENDIF
      !
      PWG(JJ,JP)=MIN(PWG(JJ,JP),PWSAT(JJ))
      !
    ENDDO
  ENDDO
ENDIF
!
IF (LHOOK) CALL DR_HOOK('CONTROL_WATER_BUDGET_TOPD',1,ZHOOK_HANDLE)
!
END SUBROUTINE CONTROL_WATER_BUDGET_TOPD

