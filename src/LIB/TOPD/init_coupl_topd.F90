!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!-----------------------------------------------------------------
!     #######################
      SUBROUTINE INIT_COUPL_TOPD (DEC, IO, S, K, NP, NPE, UG, U, HPROGRAM )
!     #######################
!
!!****  *INIT_COUPL_TOPD*  
!!
!!    PURPOSE
!!    -------
!!     This routine aims at initialising the variables 
!     needed for coupling with Topmodel.
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
!!      Original   16/10/2003
!!      Modif BV : supression of variables specific to Topmodel
!!      20/12/2007 - mll : Adaptation between a lonlat grid system for ISBA
!!                         and lambert II projection for topmodel
!!      11/2011: Modif BV : Creation of masks between ISBA and TOPODYN
!                transfered in PGD step (routine init_pgd_topd)
!!      03/2014: Modif BV : New organisation for first time step (displacement
!!                          in coupl_topd)
!!      07/2017 : Modif BV : Adding new variables for satuared area and saturation excess 
!-------------------------------------------------------------------------------
!
!*       0.     DECLARATIONS
!               ------------
!
! Modules
!
USE MODD_DIAG_EVAP_ISBA_n, ONLY : DIAG_EVAP_ISBA_t
USE MODD_ISBA_OPTIONS_n, ONLY : ISBA_OPTIONS_t
USE MODD_ISBA_n, ONLY : ISBA_S_t, ISBA_K_t, ISBA_NP_t, ISBA_NPE_t

!
USE MODD_SURF_ATM_GRID_n, ONLY : SURF_ATM_GRID_t
USE MODD_SURF_ATM_n, ONLY : SURF_ATM_t
!
!
USE MODD_COUPLING_TOPD, ONLY : XWSTOPI, XWFCTOPI, XDTOPI, &
                               XAS_NATURE, XAIBV_F, XATOP, XATOP_NATURE,XAS_IBV_P,&
                               XWWTOPI, XWTOPT, XAVG_RUNOFFCM, XAVG_DRAINCM,&
                               XDTOPT, XKA_PRE, XKAC_PRE, NMASKI, XDMAXFC, &
                               XWG_FULL, XWSTOPT, XWFCTOPT, XWWTOPT, NMASKT, & 
                               NNBV_IN_MESH, XBV_IN_MESH, XTOTBV_IN_MESH,&
                               XRUNOFF_IBV_P, NNPIX,&
                               XFRAC_D2, XFRAC_D3, XWGI_FULL,&
                               XRUN_TOROUT, XDR_TOROUT,&
                               LSTOCK_TOPD,NNB_STP_RESTART,NMASKT_PATCH,&
                               LPERT_PARAM, XAS_IBV_P,XWOVSATI_P
USE MODD_DUMMY_EXP_PROFILE,ONLY :XF_PARAM, XC_DEPTH_RATIO
USE MODD_TOPODYN,       ONLY : NNCAT, XMPARA, NMESHT, XDXT,&
                                 NNMC, XRTOP_D2, NNB_TOPD_STEP,  XDMAXT
!
USE MODD_SURF_PAR,         ONLY : XUNDEF, NUNDEF

!
USE MODI_OPEN_FILE
USE MODI_WRITE_FILE_ISBAMAP
USE MODI_CLOSE_FILE
! Interfaces
USE MODI_GET_LUOUT
USE MODI_READ_FILE_MASKTOPD
USE MODI_PACK_SAME_RANK
USE MODI_UNPACK_SAME_RANK
USE MODI_ISBA_TO_TOPD
USE MODI_RESTART_COUPL_TOPD
USE MODI_AVG_PATCH_WG
USE MODI_DG_DFTO3L
!
USE MODE_SOIL
USE MODE_SOIL_PERT
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
!
TYPE(DIAG_EVAP_ISBA_t), INTENT(INOUT) :: DEC
TYPE(ISBA_OPTIONS_t), INTENT(INOUT) :: IO
TYPE(ISBA_S_t), INTENT(INOUT) :: S
TYPE(ISBA_K_t), INTENT(INOUT) :: K
TYPE(ISBA_NP_t), INTENT(INOUT) :: NP
TYPE(ISBA_NPE_t), INTENT(INOUT) :: NPE
!
TYPE(SURF_ATM_GRID_t), INTENT(INOUT) :: UG
TYPE(SURF_ATM_t), INTENT(INOUT) :: U
!
 CHARACTER(LEN=*), INTENT(IN) :: HPROGRAM   
!
!*      0.2    declarations of local variables
!
REAL, DIMENSION(:), ALLOCATABLE   :: ZSAND_FULL, ZCLAY_FULL, ZDG_FULL ! Isba variables on the full domain
REAL, DIMENSION(:), ALLOCATABLE   :: ZFRAC    ! fraction of SurfEx mesh that covers one or several catchments
REAL, DIMENSION(:),ALLOCATABLE    :: ZSANDTOPI, ZCLAYTOPI!, ZWWILTTOPI !sand and clay fractions on TOPMODEL layers
!
!ludo
REAL, DIMENSION(:), ALLOCATABLE   :: ZKSAT_FULL  !ksat surf 
REAL, DIMENSION(:), ALLOCATABLE   :: ZDG2_FULL, ZDG3_FULL, ZWG2_FULL, ZWG3_FULL, ZRTOP_D2
!                                          
REAL                              :: ZAVG_CLAY,ZAVG_SAND
INTEGER                   :: JJ,JI            ! loop control 
INTEGER                   :: JCAT,JMESH      ! loop control 
INTEGER                   :: ILUOUT           ! Logical unit for output filr
!
REAL, DIMENSION(U%NDIM_NATURE)  :: ZRUNOFFD_NAT   
REAL, DIMENSION(U%NDIM_NATURE,3)  :: ZWG_3L, ZWGI_3L, ZDG_3L   
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('INIT_COUPL_TOPD',0,ZHOOK_HANDLE)
!
 CALL GET_LUOUT(HPROGRAM,ILUOUT)
!
WRITE(ILUOUT,*) 'INITIALISATION INIT_COUPL_TOPD'
!
ALLOCATE(NMASKT(NNCAT,NMESHT))
NMASKT(:,:) = NUNDEF
ALLOCATE(XWOVSATI_P(U%NDIM_NATURE))
XWOVSATI_P(:)=0.
!
!*       1    Initialization:
!               ---------------
ALLOCATE(NMASKT_PATCH(U%NDIM_NATURE))
!
IF (IO%CISBA=='DIF') THEN
  CALL DG_DFTO3L(IO, NP, ZDG_3L)
ELSEIF (IO%CISBA=='3-L') THEN
  CALL AVG_PATCH_WG(IO, NP, NPE, ZWG_3L, ZWGI_3L, ZDG_3L)
ENDIF
!
! la surface saturee, à l'initialisation est nulle, donc on initialise les 
! lambdas de telle sorte qu'aucun pixel ne soit sature
!
ALLOCATE(XKA_PRE (NNCAT,NMESHT))
ALLOCATE(XKAC_PRE(NNCAT))
XKA_PRE(:,:) = 0.0
XKAC_PRE(:)  = MAXVAL(XKA_PRE) + 1.
!
!Cumulated runoff and drainage initialisation
!
IF(.NOT.ALLOCATED(XAVG_RUNOFFCM)) ALLOCATE(XAVG_RUNOFFCM(U%NDIM_NATURE))
XAVG_RUNOFFCM(:) = DEC%XRUNOFF(:)
!
IF(.NOT.ALLOCATED(XAVG_DRAINCM )) ALLOCATE(XAVG_DRAINCM (U%NDIM_NATURE))
XAVG_DRAINCM (:) = DEC%XDRAIN(:)
!
!
! Reading masks
 CALL READ_FILE_MASKTOPD(U%NDIM_FULL)
!
WRITE(ILUOUT,*) ' INIT_COUPL_TOPD 1. ok'
!*      2.1     Fraction of SurfEx mesh with TOPMODEL
!               -------------------------------------
!
ALLOCATE(NNBV_IN_MESH  (U%NDIM_FULL,NNCAT))
ALLOCATE(XBV_IN_MESH   (U%NDIM_FULL,NNCAT))
ALLOCATE(XTOTBV_IN_MESH(U%NDIM_FULL))
XTOTBV_IN_MESH(:) = 0.0
ALLOCATE(XAIBV_F(U%NDIM_FULL,NNCAT)) ! fraction covered by each catchment FULL GRID
XAIBV_F(:,:)=0.0
!
!
DO JJ=1,U%NDIM_FULL !full grid
  !
  XBV_IN_MESH(JJ,:) = 0.0
  !
  DO JI=1,NNCAT
    NNBV_IN_MESH(JJ,JI) = COUNT( NMASKI(JJ,JI,:)/=NUNDEF )
    XBV_IN_MESH (JJ,JI) = REAL(NNBV_IN_MESH(JJ,JI)) * XDXT(JI)**2
    XTOTBV_IN_MESH (JJ) = XTOTBV_IN_MESH(JJ) + XBV_IN_MESH(JJ,JI)
    !
    IF (XTOTBV_IN_MESH(JJ)> UG%G%XMESH_SIZE(JJ).AND.UG%G%XMESH_SIZE(JJ)/=XUNDEF) THEN
      XBV_IN_MESH(JJ,JI) = XBV_IN_MESH(JJ,JI) * UG%G%XMESH_SIZE(JJ)/XTOTBV_IN_MESH(JJ)
      XTOTBV_IN_MESH (JJ) = UG%G%XMESH_SIZE(JJ)
    ENDIF
    !
    IF(UG%G%XMESH_SIZE(JJ)/=XUNDEF)&
    XAIBV_F(JJ,JI)= XBV_IN_MESH(JJ,JI)/UG%G%XMESH_SIZE(JJ)
  ENDDO
  !
ENDDO
!
WRITE(ILUOUT,*) ' INIT_COUPL_TOPD 2.1 ok'
!
!*      2.2     Fraction of SurfEx mesh with each catchment
!               -------------------------------------------
!
ALLOCATE(XRUNOFF_IBV_P(U%NDIM_NATURE,NNCAT))
XRUNOFF_IBV_P(:,:) = 0.0
!
ALLOCATE(ZFRAC(U%NDIM_FULL))  ! fraction not covered by catchments
ZFRAC(:) = ( UG%G%XMESH_SIZE(:)-XTOTBV_IN_MESH(:) ) / UG%G%XMESH_SIZE(:)
ZFRAC(:) = MIN(MAX(ZFRAC(:),0.),1.)
!
ALLOCATE(XATOP(U%NDIM_FULL)) ! fraction covered by catchments FULL GRID
 XATOP=1.-ZFRAC
ALLOCATE(XATOP_NATURE(U%NDIM_NATURE)) ! fraction covered by catchments NATURE GRID
 CALL PACK_SAME_RANK(U%NR_NATURE,(1.-ZFRAC),XATOP_NATURE)
!
!
IF (HPROGRAM=='POST  ') GOTO 10
WRITE(ILUOUT,*) ' INIT_COUPL_TOPD 2.2 ok'
!
!*      3.0     Wsat, Wfc and depth for TOPODYN on ISBA grid
!               --------------------------------------------
!*      3.1     clay, sand fraction, depth hydraulic conductivity at saturation of the layer for TOPODYN
!               ---------------------------------------------------------
!
ALLOCATE(ZSAND_FULL(U%NDIM_FULL))
ALLOCATE(ZCLAY_FULL(U%NDIM_FULL))
CALL UNPACK_SAME_RANK(U%NR_NATURE,K%XSAND(:,2),ZSAND_FULL)
CALL UNPACK_SAME_RANK(U%NR_NATURE,K%XCLAY(:,2),ZCLAY_FULL)
!
ZAVG_SAND=SUM(ZSAND_FULL(:),MASK=ZSAND_FULL(:)/=XUNDEF)/COUNT(ZSAND_FULL(:)/=XUNDEF)
ZAVG_CLAY=SUM(ZCLAY_FULL(:),MASK=ZCLAY_FULL(:)/=XUNDEF)/COUNT(ZCLAY_FULL(:)/=XUNDEF)
!
DO JJ=1,U%NDIM_FULL
  IF(ZSAND_FULL(JJ)==XUNDEF) ZSAND_FULL(JJ)=ZAVG_SAND
  IF(ZCLAY_FULL(JJ)==XUNDEF) ZCLAY_FULL(JJ)=ZAVG_CLAY
ENDDO
!
!ludo prof variable pour tr lat (OK car sol homogene verticalement, faux sinon)
ALLOCATE(ZDG2_FULL(U%NDIM_FULL))
ALLOCATE(ZDG3_FULL(U%NDIM_FULL))
CALL UNPACK_SAME_RANK(U%NR_NATURE,ZDG_3L(:,2),ZDG2_FULL)
CALL UNPACK_SAME_RANK(U%NR_NATURE,ZDG_3L(:,3),ZDG3_FULL)
!
!
ALLOCATE(ZRTOP_D2(U%NDIM_FULL))
ZRTOP_D2(:) = 0.
!
DO JMESH=1,U%NDIM_FULL
  IF ( ZDG2_FULL(JMESH)/=XUNDEF .AND. ZFRAC(JMESH)<1. ) THEN
    DO JCAT=1,NNCAT
     !moyenne ponderee pour cas ou plusieurs BV sur maille
       ZRTOP_D2(JMESH) = ZRTOP_D2(JMESH) + XRTOP_D2(JCAT)*MIN(XBV_IN_MESH(JMESH,JCAT)/XTOTBV_IN_MESH(JMESH),1.)    
    END DO
  ENDIF   
ENDDO
!ZTOP_D2 * D2 < D3 : the depth concerned by lateral transfers is lower than D2
WHERE( ZDG2_FULL/=XUNDEF .AND. ZRTOP_D2*ZDG2_FULL>ZDG3_FULL ) ZRTOP_D2(:) = ZDG3_FULL(:)/ZDG2_FULL(:)
!
DEALLOCATE(ZFRAC)
!
ALLOCATE(XFRAC_D2 (U%NDIM_FULL))
ALLOCATE(XFRAC_D3 (U%NDIM_FULL))
XFRAC_D2(:)=1.
XFRAC_D3(:)=0.
!
IF (IO%CISBA=='3-L') THEN
  !
  WHERE( ZDG2_FULL/=XUNDEF  ) ! if the depth is < D2
    XFRAC_D2(:) = MIN(1.,ZRTOP_D2(:))
  END WHERE
  !
  WHERE( ZDG2_FULL/=XUNDEF .AND. ZRTOP_D2*ZDG2_FULL>ZDG2_FULL  ) ! if the depth is > D2
    XFRAC_D3(:) = (ZRTOP_D2(:)*ZDG2_FULL(:)-ZDG2_FULL(:)) / (ZDG3_FULL(:)-ZDG2_FULL(:))
    XFRAC_D3(:) = MAX(0.,XFRAC_D3(:))
  END WHERE
  !
ENDIF
 !
ALLOCATE(ZDG_FULL(U%NDIM_FULL))
 CALL UNPACK_SAME_RANK(NP%AL(1)%NR_P,NP%AL(1)%XRUNOFFD(:),ZRUNOFFD_NAT)
 CALL UNPACK_SAME_RANK(U%NR_NATURE,ZRUNOFFD_NAT,ZDG_FULL)
!
ALLOCATE(ZSANDTOPI(U%NDIM_FULL))
ALLOCATE(ZCLAYTOPI(U%NDIM_FULL))
ZSANDTOPI(:)=0.0
ZCLAYTOPI(:)=0.0
ALLOCATE(XDTOPI(U%NDIM_FULL))
XDTOPI(:)=0.0
!
WHERE ( ZDG_FULL/=XUNDEF .AND. ZDG_FULL/=0. )
  XDTOPI = ZDG_FULL
  ZSANDTOPI = ZSANDTOPI + ZSAND_FULL * ZDG_FULL
  ZCLAYTOPI = ZCLAYTOPI + ZCLAY_FULL * ZDG_FULL
  ZSANDTOPI = ZSANDTOPI / XDTOPI
  ZCLAYTOPI = ZCLAYTOPI / XDTOPI
ELSEWHERE
  ZSANDTOPI = XUNDEF
  ZCLAYTOPI = XUNDEF
  XDTOPI = XUNDEF
END WHERE
!
DEALLOCATE(ZSAND_FULL)
DEALLOCATE(ZCLAY_FULL)
WRITE(ILUOUT,*) ' INIT_COUPL_TOPD 3.1 ok'
!
!*      4.1     depth of the Isba layer on TOP-LAT grid
!               ---------------------------------------
!
ALLOCATE(XDTOPT(NNCAT,NMESHT))
XDTOPT(:,:) = 0.0
CALL ISBA_TO_TOPD(XDTOPI,XDTOPT)
WRITE(ILUOUT,*) ' INIT_COUPL_TOPD 4.1 ok'
!!
!
!*      3.2     Wsat and Wfc on TOPODYN layer
!               -----------------------------
!
ALLOCATE(XWSTOPI   (U%NDIM_FULL))
ALLOCATE(XWFCTOPI  (U%NDIM_FULL))
XWSTOPI (:) = 0.0
XWFCTOPI(:) = 0.0
!
IF (LPERT_PARAM) THEN
  !
  XWSTOPI    = WSAT_FUNC_PERT_1D (ZCLAYTOPI,ZSANDTOPI,IO%CPEDOTF)
  IF (IO%CISBA=='2-L' .OR. IO%CISBA=='3-L') THEN
    !  field capacity at hydraulic conductivity = 0.1mm/day
    XWFCTOPI   = WFC_FUNC_PERT_1D(ZCLAYTOPI,ZSANDTOPI,IO%CPEDOTF)
  ELSE IF (IO%CISBA=='DIF') THEN
    !  field capacity at water potential = 0.33bar        
    XWFCTOPI   = W33_FUNC_PERT_1D(ZCLAYTOPI,ZSANDTOPI,IO%CPEDOTF)
  END IF
  !
ELSE
  !
  XWSTOPI    = WSAT_FUNC_1D (ZCLAYTOPI,ZSANDTOPI,IO%CPEDOTF)
  IF (IO%CISBA=='2-L' .OR. IO%CISBA=='3-L') THEN
    !  field capacity at hydraulic conductivity = 0.1mm/day
    XWFCTOPI   = WFC_FUNC_1D  (ZCLAYTOPI,ZSANDTOPI,IO%CPEDOTF)
  ELSE IF (IO%CISBA=='DIF') THEN
    !  field capacity at water potential = 0.33bar         
    XWFCTOPI   = WFC_FUNC_1D  (ZCLAYTOPI,ZSANDTOPI,IO%CPEDOTF)
  END IF
  !
END IF
!
!modif ludo test ksat exp
WRITE(ILUOUT,*) 'CKSAT==',IO%CKSAT

ALLOCATE(ZKSAT_FULL(U%NDIM_FULL))
ZKSAT_FULL(:) = 0.0
ALLOCATE(XWWTOPI(U%NDIM_FULL))
XWWTOPI(:) = 0.0

IF (LPERT_PARAM) THEN
  XWWTOPI(:) = WWILT_FUNC_PERT(ZCLAYTOPI,ZSANDTOPI,IO%CPEDOTF)
ELSE
  XWWTOPI(:) = WWILT_FUNC(ZCLAYTOPI,ZSANDTOPI,IO%CPEDOTF)
END IF
!
DEALLOCATE(ZSANDTOPI)
DEALLOCATE(ZCLAYTOPI)
DEALLOCATE(ZRTOP_D2)
!
WRITE(ILUOUT,*) ' INIT_COUPL_TOPD 3.2 ok'
!*      4.3     Ko on TOP-LAT grid
!               ------------------
!
ALLOCATE(XWWTOPT(NNCAT,NMESHT))
 CALL ISBA_TO_TOPD(XWWTOPI,XWWTOPT)
WHERE (XWWTOPT == XUNDEF) XWWTOPT = 0.0
!
ALLOCATE(XWG_FULL(U%NDIM_FULL))
ALLOCATE(XWGI_FULL(U%NDIM_FULL))
ALLOCATE(XWTOPT(NNCAT,NMESHT))
XWTOPT(:,:) = 0.0
ALLOCATE(XWSTOPT (NNCAT,NMESHT))
ALLOCATE(XWFCTOPT(NNCAT,NMESHT))
ALLOCATE(XDMAXFC(NNCAT,NMESHT))
XDMAXFC(:,:) = XUNDEF
ALLOCATE(XDMAXT(NNCAT,NMESHT))
XDMAXT(:,:)=XUNDEF
!
ALLOCATE(ZWG2_FULL(U%NDIM_FULL))
ALLOCATE(ZWG3_FULL(U%NDIM_FULL))
!
IF (IO%CISBA=='3-L')THEN
  !
  CALL UNPACK_SAME_RANK(U%NR_NATURE,ZWG_3L(:,2),ZWG2_FULL)
  CALL UNPACK_SAME_RANK(U%NR_NATURE,ZWG_3L(:,3),ZWG3_FULL)
  !
  DO JMESH=1,U%NDIM_FULL
    IF ( ZDG_FULL(JMESH)/=XUNDEF .AND. ZDG_FULL(JMESH)/=0. )THEN
      XWG_FULL(JMESH) = XFRAC_D2(JMESH)*(ZDG2_FULL(JMESH)/ZDG_FULL(JMESH))*ZWG2_FULL(JMESH)&
                  + XFRAC_D3(JMESH)*((ZDG3_FULL(JMESH)-ZDG2_FULL(JMESH))/ZDG_FULL(JMESH))*ZWG3_FULL(JMESH)
    ELSE
      XWG_FULL(JMESH) = XUNDEF
    ENDIF
  ENDDO
  !
ENDIF
!
XWGI_FULL = 0.
WRITE(ILUOUT,*) ' INIT_COUPL_TOPD 4.3 ok'
!
!
!*      4.4     Initialisation of the previous time step water storage on topodyn-lat grid
!               --------------------------------------------------------------------------
!*      4.5     M parameter on TOPODYN grid
!               ------------------------
!*      4.5.1   Mean depth soil on catchment
!
ALLOCATE(XMPARA (NNCAT))
XMPARA  (:) = 0.0
IF (.NOT.ALLOCATED(XF_PARAM)) ALLOCATE(XF_PARAM(SIZE(S%XF_PARAM)))
WRITE(ILUOUT,*) ' INIT_COUPL_TOPD 4.5 ok'
!
! 
!*      5.0      Initial saturated area computation
!               -----------------------------------------------------------
!
ALLOCATE(XAS_NATURE(U%NDIM_NATURE))
XAS_NATURE(:) = 0.0
ALLOCATE(XAS_IBV_P(U%NDIM_NATURE,NNCAT))
XAS_IBV_P(:,:) = 0.0
WRITE(ILUOUT,*) ' INIT_COUPL_TOPD 5.0 ok'
!
!*      6.0     Stock management in case of restart
!               -----------------------------------------------------------
!
10 CONTINUE
!
!stock
ALLOCATE(XRUN_TOROUT(NNCAT,NNB_TOPD_STEP+NNB_STP_RESTART))
ALLOCATE(XDR_TOROUT (NNCAT,NNB_TOPD_STEP+NNB_STP_RESTART))
XRUN_TOROUT(:,:) = 0.
XDR_TOROUT (:,:) = 0.
!
IF (HPROGRAM=='POST  ') GOTO 20
!
IF (LSTOCK_TOPD) CALL RESTART_COUPL_TOPD(UG, U%NR_NATURE, HPROGRAM,U%NDIM_FULL)
WRITE(ILUOUT,*) ' INIT_COUPL_TOPD 6.0 ok'
!
!*      7.0     deallocate
!               ----------
!
DEALLOCATE(ZDG2_FULL)
DEALLOCATE(ZDG3_FULL)
DEALLOCATE(ZWG2_FULL)
DEALLOCATE(ZWG3_FULL)
!
20 CONTINUE
!
IF (LHOOK) CALL DR_HOOK('INIT_COUPL_TOPD',1,ZHOOK_HANDLE)
!
END SUBROUTINE INIT_COUPL_TOPD
