!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!-------------------------------------------------------------------------------
!     ####################
      SUBROUTINE ISBA_TO_TOPDSAT(PKAPPA,PKAPPAC,KI,PRO_I,PRO_T)
!     ####################
!
!!****  *ISBA_TO_TOPDSAT*  
!!
!!    PURPOSE
!!    -------
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
!!      Original   23/11/2005
!!      Modif B Vincendon 07/2017 : more control for UNDEF variables
!-------------------------------------------------------------------------------
!
!*       0.     DECLARATIONS
!               ------------
!
USE MODD_SURF_PAR,  ONLY : XUNDEF,NUNDEF
!
USE MODD_TOPODYN, ONLY : NNCAT, NNMC, NMESHT
USE MODD_COUPLING_TOPD,ONLY: NMASKI, NMASKT, NNBV_IN_MESH
!
USE YOMHOOK   ,ONLY : LHOOK, DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
INTEGER, INTENT(IN)               :: KI      ! Number of Isba meshes
REAL, DIMENSION(:,:), INTENT(IN)  :: PKAPPA  ! Hydrological indexes on the catchments 
                                             ! at the previous time step
REAL, DIMENSION(:), INTENT(IN)    :: PKAPPAC ! Hydrological index at saturation at the 
                                             ! previous time step
REAL, DIMENSION(:,:), INTENT(IN)  :: PRO_I   ! Runoff on Isba grid for each catchment
REAL, DIMENSION(:,:), INTENT(OUT) :: PRO_T   ! Runoff on TOPODYN grid
!
!
!*      0.2    declarations of local variables
!
INTEGER                          :: JCAT, JPIX, JMESH_ISBA,JJ ! Loop indexes
INTEGER, DIMENSION(NNCAT,KI)     :: INSAT       ! number of saturated pixels in an ISBA mesh
INTEGER, DIMENSION(NNCAT,KI)     :: INDRY       ! Number of non-saturated pixels in an ISBA mesh
REAL, DIMENSION(NNCAT,NMESHT)    :: ZROSAT      ! 
REAL, DIMENSION(NNCAT,NMESHT)    :: ZRODRY      ! 
!
REAL::ZSMALL,ZTMP,ZTMP2
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('ISBA_TO_TOPDSAT',0,ZHOOK_HANDLE)
!
!*       0.     Initialization :
!               --------------
!
INSAT(:,:)=0
INDRY(:,:)=0
ZROSAT(:,:)=0.0
ZRODRY(:,:)=0.0
!Full grid
! Only Isba meshes over studied catchments are scanned
DO JMESH_ISBA = 1,KI
  !
  DO JCAT = 1,NNCAT
    !
    JJ=1
    JPIX=NMASKI(JMESH_ISBA,JCAT,JJ)
    !
    DO WHILE (JPIX/=NUNDEF .AND.(JJ<=SIZE(NMASKI,3)))
      !
      IF (PKAPPA(JCAT,JPIX)/=XUNDEF .AND. PKAPPAC(JCAT)/=XUNDEF.AND.&
          NMASKT(JCAT,JPIX)/=NUNDEF .AND. NMASKT(JCAT,JPIX)/=0.AND.&
          PRO_I(NMASKT(JCAT,JPIX),JCAT)/=XUNDEF) THEN
        ! Calculation of the saturated and dry catchment pixels in each Isba mesh
        IF (PKAPPA(JCAT,JPIX).GE.PKAPPAC(JCAT)) THEN
          INSAT(JCAT,NMASKT(JCAT,JPIX)) = INSAT(JCAT,NMASKT(JCAT,JPIX)) + 1
          ZROSAT(JCAT,JPIX) = PRO_I(NMASKT(JCAT,JPIX),JCAT)
        ELSE
          INDRY(JCAT,NMASKT(JCAT,JPIX)) = INDRY(JCAT,NMASKT(JCAT,JPIX)) + 1 
          ZRODRY(JCAT,JPIX) = PRO_I(NMASKT(JCAT,JPIX),JCAT)
        ENDIF
      ENDIF
      !
      JJ=JJ+1
      IF (JJ<=SIZE(NMASKI,3)) JPIX=NMASKI(JMESH_ISBA,JCAT,JJ)
      !
    ENDDO
    !
  ENDDO
  !
ENDDO
!
!
DO JCAT = 1,NNCAT
  !
  DO JPIX = 1,NNMC(JCAT)
    !
    IF ( NMASKT(JCAT,JPIX)/=NUNDEF.AND.NMASKT(JCAT,JPIX)/=0) THEN
      ! calculation of the runoff and deep drainage to rout in each Isba mesh, for each catchment
      IF (INSAT(JCAT,NMASKT(JCAT,JPIX))>0 .AND. INSAT(JCAT,NMASKT(JCAT,JPIX))/=XUNDEF.AND.&
          ZROSAT(JCAT,JPIX)/=XUNDEF ) THEN
        PRO_T(JCAT,JPIX) = ZROSAT(JCAT,JPIX) / INSAT(JCAT,NMASKT(JCAT,JPIX))

        ! if no runoff : calculation of the deep drainage to rout in each Isba mesh for each catchment
      ELSEIF (INDRY(JCAT,NMASKT(JCAT,JPIX))>0 .AND. INDRY(JCAT,NMASKT(JCAT,JPIX))/=XUNDEF.AND.&
           ZRODRY(JCAT,JPIX)/=XUNDEF) THEN
        PRO_T(JCAT,JPIX) = ZRODRY(JCAT,JPIX) / INDRY(JCAT,NMASKT(JCAT,JPIX))
      ELSE
        PRO_T(JCAT,JPIX) = 0.
      ENDIF
    ENDIF
    !
  ENDDO
  !
  ! budget control 
  ZTMP=0.
  ZTMP2=0.
  !
  DO JPIX = 1,NNMC(JCAT)
    !
    IF (PRO_T(JCAT,JPIX)/=XUNDEF) ZTMP = ZTMP + PRO_T(JCAT,JPIX)
    IF (NMASKT(JCAT,JPIX)/=NUNDEF.AND.NMASKT(JCAT,JPIX)/=0 .AND. &
        PRO_I(NMASKT(JCAT,JPIX),JCAT)/=XUNDEF .AND. &
        NNBV_IN_MESH(NMASKT(JCAT,JPIX),JCAT)/=0 .AND.&
        NNBV_IN_MESH(NMASKT(JCAT,JPIX),JCAT)/=NUNDEF) &
      ZTMP2 = ZTMP2 + (PRO_I(NMASKT(JCAT,JPIX),JCAT)/NNBV_IN_MESH(NMASKT(JCAT,JPIX),JCAT))
    !
  ENDDO!JPIX
  !
  ZSMALL=ABS(ZTMP2*0.001)
  !
  IF( ABS(ZTMP-ZTMP2) > ZSMALL ) THEN
    WHERE ( PRO_T(JCAT,:)/=XUNDEF )
      PRO_T(JCAT,:) = PRO_T(JCAT,:)- ((ZTMP-ZTMP2)/NNMC(JCAT))
    ENDWHERE
  ENDIF
  !
ENDDO
!
IF (LHOOK) CALL DR_HOOK('ISBA_TO_TOPDSAT',1,ZHOOK_HANDLE)
!
END SUBROUTINE ISBA_TO_TOPDSAT
