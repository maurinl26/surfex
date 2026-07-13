!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!-----------------------------------------------------------------
!     ####################
      SUBROUTINE TOPD_TO_ISBA (K, PEK, UG, U, KI,KSTEP,GTOPD)
!     ####################
!
!!****  *TOPD_TO_ISBA*  
!!
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
!!      Original   09/10/2003
!!                 03/2014 (B. Vincendon) correction for meshes covered by several watersheds
!!                 03/2015 (E. Artinyan) YSTEP jusqu'a 99999 steps
!!                 07/2022 (B. Decharme) MOD(KSTEP,NFREQ_MAPS_WG) crash if NFREQ_MAPS_WG==0
!-------------------------------------------------------------------------------
!
!*       0.     DECLARATIONS
!               ------------
!
USE MODD_ISBA_n, ONLY : ISBA_K_t, ISBA_PE_t
USE MODD_SURF_ATM_GRID_n, ONLY : SURF_ATM_GRID_t
USE MODD_SURF_ATM_n, ONLY : SURF_ATM_t
!
USE MODI_UNPACK_SAME_RANK
USE MODI_PACK_SAME_RANK
!
USE MODI_WRITE_FILE_ISBAMAP
USE MODI_OPEN_FILE
USE MODI_CLOSE_FILE
!
USE MODD_TOPD_PAR, ONLY : NUNIT
USE MODD_TOPODYN,         ONLY : NNCAT, NNMC, NNB_TOPD_STEP
USE MODD_COUPLING_TOPD,   ONLY : XWG_FULL, XDTOPT, XWTOPT, XWOVSATI_P,&
                                 NMASKT, XTOTBV_IN_MESH,NFREQ_MAPS_WG, XBV_IN_MESH
!
USE MODD_SURF_PAR,        ONLY : XUNDEF,NUNDEF
USE MODD_ISBA_PAR,        ONLY : XWGMIN
!
USE YOMHOOK   ,ONLY : LHOOK, DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
!
!
TYPE(ISBA_K_t), INTENT(INOUT) :: K
TYPE(ISBA_PE_t), POINTER :: PEK
TYPE(SURF_ATM_GRID_t), INTENT(INOUT) :: UG
TYPE(SURF_ATM_t), INTENT(INOUT) :: U
!
INTEGER, INTENT(IN)                 :: KI      ! Grid dimensions
INTEGER, INTENT(IN)                 :: KSTEP   ! Topodyn current time step
LOGICAL, DIMENSION(:), INTENT(INOUT)   :: GTOPD     ! 
!
!*      0.2    declarations of local variables
!
!
INTEGER                            :: JJ, JI , JMESH, JCAT         ! loop control 
REAL, DIMENSION(U%NDIM_FULL)       :: ZW              ! TOPODYN water content on ISBA grid (mm)
REAL, DIMENSION(U%NDIM_FULL)       :: ZWSAT_FULL      ! Water content at saturation on the layer 2 
                                                      ! on the full grid
REAL, DIMENSION(U%NDIM_FULL)       :: ZWG_OLD

REAL, DIMENSION(U%NDIM_FULL)       :: ZWOVSATI_F
!
REAL, DIMENSION(U%NDIM_FULL,NNCAT) :: ZCOUNT, ZW_CAT
!
CHARACTER(LEN=5)                   :: YSTEP
!
LOGICAL                            :: LWORK
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('TOPD_TO_ISBA',0,ZHOOK_HANDLE)
!
!*       0.     Initialization:
!               ---------------
!
ZW(:)= 0.0
ZW_CAT(:,:)= 0.0
ZCOUNT(:,:)=0.0
!
ZWOVSATI_F(:)=0.
XWOVSATI_P(:)=0.
CALL UNPACK_SAME_RANK(U%NR_NATURE,K%XWSAT(:,2),ZWSAT_FULL)
CALL UNPACK_SAME_RANK(U%NR_NATURE,PEK%XWG(:,2),ZWG_OLD)
!
!*       1.     TOPODYN-LAT => ISBA
!               -------------------
!*       1.1    mobilizable water
!               -----------------
!
DO JJ=1,NNCAT
  DO JI=1,NNMC(JJ)
    IF ( (XDTOPT(JJ,JI) /= XUNDEF).AND. (NMASKT(JJ,JI) /= NUNDEF) .AND.XWTOPT(JJ,JI)/=XUNDEF)THEN
      ZW_CAT(NMASKT(JJ,JI),JJ) = ZW_CAT(NMASKT(JJ,JI),JJ) + XWTOPT(JJ,JI)
      ZCOUNT(NMASKT(JJ,JI),JJ) = ZCOUNT(NMASKT(JJ,JI),JJ) + 1.0
    ENDIF
  ENDDO
ENDDO
!
!
DO JMESH=1,U%NDIM_FULL
  !
  IF (XTOTBV_IN_MESH(JMESH)>0.0 .AND. XTOTBV_IN_MESH(JMESH)/=XUNDEF ) THEN
  ! at least 1 catchment over mesh
    DO JCAT=1,NNCAT
      !
      IF (XBV_IN_MESH(JMESH,JCAT)>0.) THEN ! Catchment JCAT is on Jmesh
        IF(ZW_CAT(JMESH,JCAT)/=XUNDEF .AND.ZCOUNT(JMESH,JCAT)/=0.)THEN
          ZW(JMESH) = ZW(JMESH) +ZW_CAT(JMESH,JCAT) /  ZCOUNT(JMESH,JCAT)*&
                 MIN(1.0,(XBV_IN_MESH(JMESH,JCAT)/UG%G%XMESH_SIZE(JMESH)))
        ENDIF
      ENDIF
      !
    ENDDO !JCAT
    !
    ZW(JMESH) = ZW(JMESH) + ZWG_OLD(JMESH) *&
                 MIN(1.0,(UG%G%XMESH_SIZE(JMESH)-XTOTBV_IN_MESH(JMESH))/UG%G%XMESH_SIZE(JMESH))
  ENDIF ! (XTOTBV_IN_MESH(JMESH)>0.0 .AND. XTOTBV_IN_MESH(JMESH)/=XUNDEF )
  !
  IF (ZW(JMESH)/=0.) THEN
    XWG_FULL(JMESH) = MAX(ZW(JMESH),XWGMIN)
  ENDIF
  !
  IF ( XWG_FULL(JMESH) > ZWSAT_FULL(JMESH) .AND.&
      XWG_FULL(JMESH)/=XUNDEF .AND. ZWSAT_FULL(JMESH)/=XUNDEF )THEN
     !ludo calcul sat avant wg
    ZWOVSATI_F(JMESH) = ZWOVSATI_F(JMESH)+MAX(0.,XWG_FULL(JMESH) - ZWSAT_FULL(JMESH))
    IF (ZWOVSATI_F(JMESH) <  0.0)ZWOVSATI_F(JMESH)=0.
    XWG_FULL(JMESH) = ZWSAT_FULL(JMESH)
  ENDIF
  !
ENDDO!JMESH
!
CALL PACK_SAME_RANK(U%NR_NATURE,ZWOVSATI_F,XWOVSATI_P)
!
LWORK=.FALSE.
IF (NFREQ_MAPS_WG/=0) THEN
   LWORK=MOD(KSTEP,NFREQ_MAPS_WG)==0
ENDIF
!
!IF ( (NFREQ_MAPS_WG/=0 .AND. MOD(KSTEP,NFREQ_MAPS_WG)==0) .OR. (KSTEP==NNB_TOPD_STEP)) THEN
IF ( LWORK .OR. (KSTEP==NNB_TOPD_STEP)) THEN
  ! writing of YSTEP to be able to write maps
  IF (KSTEP<10) THEN
    WRITE(YSTEP,'(I1)') KSTEP
  ELSEIF (KSTEP < 100) THEN
    WRITE(YSTEP,'(I2)') KSTEP
  ELSEIF (KSTEP < 1000) THEN
    WRITE(YSTEP,'(I3)') KSTEP
  ELSEIF (KSTEP < 10000) THEN
    WRITE(YSTEP,'(I4)') KSTEP
  ELSE
    WRITE(YSTEP,'(I5)') KSTEP
  ENDIF
  !  
  CALL OPEN_FILE('ASCII ',NUNIT,HFILE='carte_w'//YSTEP,HFORM='FORMATTED',HACTION='WRITE')
  CALL WRITE_FILE_ISBAMAP(UG,NUNIT,XWG_FULL,U%NDIM_FULL)
  CALL CLOSE_FILE('ASCII ',NUNIT)
  !
ENDIF
!
IF (LHOOK) CALL DR_HOOK('TOPD_TO_ISBA',1,ZHOOK_HANDLE)
!
END SUBROUTINE TOPD_TO_ISBA
