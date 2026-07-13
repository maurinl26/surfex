!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!-----------------------------------------------------------------
!     #####################
      SUBROUTINE COUPL_TOPD (DEC, DC, DMI, PMESH_SIZE, IO, S, K, NK, NP, NPE, &
                             UG, U, HPROGRAM, HSTEP, KI, KSTEP)
!     #####################
!
!!****  *COUPL_TOPD*  
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
!!      Original   15/10/2003
!!      09/2007 : New organisation of exfiltration, computation of saturated
!!                area, routing.
!!                Soil ice content taken into account
!!      09/2013 : Modifications to be able to run with ISBA-DF and more than 1
!!                patch
!!      03/2014: Modif BV : New organisation for first time step (displacement
!!                          from init_coupl_topd)
!!      07/2015: Modif BV : modification of recharge computation
!!      07/2017: Modif BV : change name of variables packed and on full grid
!!               + computation of runoff by mesh and catchment to avoid problems on catchments interfaces
!!      07/2022 (B. Decharme) MOD(KSTEP,NFREQ_MAPS_ASAT) crash if NFREQ_MAPS_ASAT==0
!-------------------------------------------------------------------------------
!
!*       0.     DECLARATIONS
!               ------------
!
!
USE MODD_DIAG_n, ONLY : DIAG_t
USE MODD_DIAG_EVAP_ISBA_n, ONLY : DIAG_EVAP_ISBA_t
USE MODD_DIAG_MISC_ISBA_n, ONLY : DIAG_MISC_ISBA_t
!
USE MODD_ISBA_OPTIONS_n, ONLY : ISBA_OPTIONS_t
USE MODD_ISBA_n, ONLY : ISBA_S_t, ISBA_K_t, ISBA_NK_t, ISBA_NP_t, ISBA_NPE_t, ISBA_PE_t
!
USE MODD_SURF_ATM_GRID_n, ONLY : SURF_ATM_GRID_t
USE MODD_SURF_ATM_n, ONLY : SURF_ATM_t
!
USE MODD_TOPD_PAR, ONLY : NUNIT
USE MODD_TOPODYN,        ONLY   : NNCAT, NMESHT, NNMC, XMPARA, XDMAXT, XQTOT,XQB_RUN,XQB_DR
USE MODD_COUPLING_TOPD,  ONLY   : XWG_FULL, XDTOPI, XKAC_PRE, XDTOPT, XWTOPT, XWSTOPT, XWWTOPT,&
                                  XAS_NATURE,&
                                  XKA_PRE, NMASKT, XWOVSATI_P, XAS_IBV_P,&
                                  XRUNOFF_IBV_P, XAIBV_F,XATOP, XWFCTOPI, NNPIX,&
                                  XFRAC_D2, XFRAC_D3, XWSTOPI, XDMAXFC, XWFCTOPT, XWGI_FULL,&
                                  NFREQ_MAPS_ASAT, XAVG_RUNOFFCM,&
                                  XAVG_DRAINCM,XRAINFALLCM,XAVG_HORTCM,&
                                  LBUDGET_TOPD, LPERT_PARAM, LPERT_INIT,&
                                  NNBV_IN_MESH,XTOTBV_IN_MESH                                  !
USE MODD_CSTS,             ONLY : XRHOLW, XRHOLI
USE MODD_SURF_PAR,         ONLY : XUNDEF, NUNDEF
USE MODD_ISBA_PAR,         ONLY : XWGMIN

USE MODD_DUMMY_EXP_PROFILE,ONLY :XF_PARAM, XC_DEPTH_RATIO
!
USE MODI_GET_LUOUT
USE MODI_UNPACK_SAME_RANK
USE MODI_PACK_SAME_RANK
USE MODI_ISBA_TO_TOPD
USE MODI_RECHARGE_SURF_TOPD
USE MODI_TOPODYN_LAT
USE MODI_SAT_AREA_FRAC
USE MODI_TOPD_TO_ISBA
USE MODI_DIAG_ISBA_TO_ROUT
USE MODI_ISBA_TO_TOPDSAT
USE MODI_ROUTING
USE MODI_OPEN_FILE
USE MODI_WRITE_FILE_ISBAMAP
USE MODI_CLOSE_FILE
USE MODI_DG_DFTO3L
USE MODI_AVG_PATCH_WG
USE MODI_DISPATCH_WG
USE MODI_TOPD_TO_DF
USE MODI_INIT_BUDGET_COUPL_ROUT
USE MODI_CONTROL_WATER_BUDGET_TOPD
!
USE MODE_RANDOM_PERT
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
TYPE(DIAG_t), INTENT(INOUT) :: DC
TYPE(DIAG_MISC_ISBA_t), INTENT(INOUT) :: DMI
REAL, DIMENSION(:), INTENT(IN) :: PMESH_SIZE
TYPE(ISBA_OPTIONS_t), INTENT(INOUT) :: IO
TYPE(ISBA_S_t), INTENT(INOUT) :: S
TYPE(ISBA_K_t), INTENT(INOUT) :: K
TYPE(ISBA_NK_t), INTENT(INOUT) :: NK
TYPE(ISBA_NP_t), INTENT(INOUT) :: NP
TYPE(ISBA_NPE_t), INTENt(INOUT) :: NPE
TYPE(SURF_ATM_GRID_t), INTENT(INOUT) :: UG
TYPE(SURF_ATM_t), INTENT(INOUT) :: U
!
CHARACTER(LEN=6), INTENT(IN) :: HPROGRAM ! program calling surf. schemes
CHARACTER(LEN=*), INTENT(IN) :: HSTEP  ! atmospheric loop index
INTEGER, INTENT(IN)          :: KI    ! Grid dimensions
INTEGER, INTENT(IN)          :: KSTEP ! current time step 
!
!*      0.2    declarations of local variables
!
TYPE(ISBA_PE_t), POINTER :: PEK
!On pixels
REAL, DIMENSION(NNCAT,NMESHT) :: ZRT             ! recharge on TOP-LAT grid (m)
REAL, DIMENSION(NNCAT,NMESHT) :: ZDEFT           ! local deficits on TOPODYN grid (m)
REAL, DIMENSION(NNCAT,NMESHT) :: ZRI_WGIT        ! water changing of phase on TOPMODEL grid
REAL, DIMENSION(NNCAT,NMESHT) :: ZRUNOFF_TOPD    ! Runoff on the Topodyn grid (m3/s)
REAL, DIMENSION(NNCAT,NMESHT) :: ZDRAIN_TOPD     ! Drainage from Isba on Topodyn grid (m3/s)
REAL, DIMENSION(NNCAT,NMESHT) :: ZKAPPA          ! topographic index
REAL, DIMENSION(NNCAT)        :: ZKAPPAC         ! critical topographic index
!On full grid
REAL, DIMENSION(U%NDIM_FULL)           :: ZRI             ! recharge on ISBA grid (m)
REAL, DIMENSION(U%NDIM_FULL)           :: ZRI_WGI         ! water changing of phase on ISBA grid
REAL, DIMENSION(U%NDIM_FULL)           :: ZWM,ZWIM        ! Water content on SurfEx grid after the previous topodyn time step
REAL, DIMENSION(U%NDIM_FULL)           :: Z_WSTOPI, Z_WFCTOPI
REAL, DIMENSION(U%NDIM_FULL)           :: ZRUNOFFC_FULL   ! Cumulated runoff from isba on the full domain (kg/m2)
REAL, DIMENSION(U%NDIM_FULL)           :: ZRUNOFFC_FULLM  ! Cumulated runoff from isba on the full domain (kg/m2) at t-dt
REAL, DIMENSION(U%NDIM_FULL)           :: ZRUNOFF_ISBA_F    ! Runoff from Isba (kg/m2)
REAL, DIMENSION(U%NDIM_FULL)           :: ZDRAINC_FULL    ! Cumulated drainage from Isba on the full domain (kg/m2)
REAL, DIMENSION(U%NDIM_FULL)           :: ZDRAINC_FULLM   ! Cumulated drainage from Isba on the full domain (kg/m2) at t-dt
REAL, DIMENSION(U%NDIM_FULL)           :: ZDRAIN_ISBA     ! Drainage from Isba (m3/s)
REAL, DIMENSION(U%NDIM_FULL)           :: ZDG_FULL
REAL, DIMENSION(U%NDIM_FULL)           :: ZWG2_FULL, ZWG3_FULL, ZDG2_FULL, ZDG3_FULL
REAL, DIMENSION(U%NDIM_FULL)           :: ZWGI_FULL !(m3/m3)
REAL, DIMENSION(U%NDIM_FULL)           :: ZWOVSATI_F !(m3/m3)
REAL, DIMENSION(U%NDIM_FULL)           :: ZRUNOFFD_F !(kg/m2)
REAL, DIMENSION(U%NSIZE_NATURE) :: ZRUNOFFD_NAT!(kg/m2)
REAL, DIMENSION(U%NDIM_FULL)           :: ZASI_F             ! Saturated area fraction for each Isba meshes on full grid
REAL, DIMENSION(U%NDIM_FULL,NNCAT)     :: ZAS_IBV_F      ! Saturated area fraction for each Isba meshes on catchment on full grid
REAL, DIMENSION(U%NDIM_FULL,NNCAT)     :: ZRUNOFF_IBV_F    ! Runoff from Isba (kg/m2)
REAL, DIMENSION(NNCAT)        :: Z_DW1,Z_DW2     ! Wsat-Wfc to actualise M in fonction of WI
REAL                          :: ZAVG_MESH_SIZE, ZWSATMAX
LOGICAL, DIMENSION(NNCAT)     :: GTOPD           ! logical variable = true if topodyn_lat runs
INTEGER                       :: JJ, JI, JL, JP    ! loop control 
INTEGER                       :: ILUOUT          ! unit number of listing file
INTEGER                       :: IACT_GROUND_LAYER, IDEPTH,IMASK, ISUM
!
REAL, DIMENSION(U%NDIM_FULL)            :: ZF_PARAM_FULL
REAL, DIMENSION(NNCAT,NMESHT)  :: ZF_PARAMT
!On isba grid (packed)
REAL,    DIMENSION(U%NSIZE_NATURE,3)     :: ZWG_3L,ZWGI_3L,ZDG_3L          
REAL,    DIMENSION(U%NSIZE_NATURE)       :: ZMESH_SIZE, ZWSAT
REAL, DIMENSION(U%NSIZE_NATURE,IO%NGROUND_LAYER,IO%NPATCH) :: ZWG_TMP
REAL, DIMENSION(U%NSIZE_NATURE,IO%NPATCH) :: ZWG, ZDG
REAL,    DIMENSION(U%NSIZE_NATURE,NNCAT) :: ZWOVSAT_IBV_P    ! Temporary variable for runoff from Isba (kg/m2)
INTEGER, DIMENSION(U%NSIZE_NATURE)       :: INPIXI_P
INTEGER, DIMENSION(U%NSIZE_NATURE,NNCAT) :: INBV_IN_MESH
! Taking several patches into account
REAL, DIMENSION(U%NSIZE_NATURE)   :: ZSUMFRD2, ZSUMFRD3
REAL, DIMENSION(U%NSIZE_NATURE,3) :: ZWG_CTL
! Perturbating the initial soil moisture field
REAL, DIMENSION(U%NDIM_FULL) :: ZRAND_MAP
REAL, DIMENSION(7)  :: ZRANDOM
!
LOGICAL                      :: LWORK
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('COUPL_TOPD',0,ZHOOK_HANDLE)
!
CALL GET_LUOUT(HPROGRAM,ILUOUT)
!
!
!*       0.     Initialization:
!               ---------------
!Nature grid
ZWSATMAX=MAXVAL(XWSTOPI,MASK=XWSTOPI/=XUNDEF)
!
ZWG_TMP(:,:,:) = 0.
DO JP = 1,IO%NPATCH
  !
  PEK => NPE%AL(JP)
  !
  DO JJ = 1,NP%AL(JP)%NSIZE_P
    IMASK = NP%AL(JP)%NR_P(JJ)
    NPE%AL(JP)%XWG(JJ,:) = MAX(NPE%AL(JP)%XWG(JJ,:),XWGMIN)
    ZWG_TMP(IMASK,:,JP)  = NPE%AL(JP)%XWG(JJ,:)
  ENDDO
ENDDO
!
DO JJ=1,NNCAT
 CALL PACK_SAME_RANK(U%NR_NATURE,NNBV_IN_MESH(:,JJ),INBV_IN_MESH(:,JJ))
ENDDO
!               ---------------
IF (IO%CISBA=='DIF') THEN
  CALL DG_DFTO3L(IO, NP, ZDG_3L)
  ZWG_3L(:,2)  = DMI%XFRD2_TWG (:)
  ZWG_3L(:,3)  = DMI%XFRD3_TWG (:)
  ZWGI_3L(:,2) = DMI%XFRD2_TWGI(:)
  ZWGI_3L(:,3) = DMI%XFRD3_TWGI(:)
ELSEIF (IO%CISBA=='3-L') THEN
  CALL AVG_PATCH_WG(IO, NP, NPE, ZWG_3L, ZWGI_3L, ZDG_3L)
ENDIF
!
!
!*       1.     ISBA => TOPODYN
!               ---------------
!*       1.1    Computation of the useful depth and water for lateral transfers
!               -----------------------------------
!
CALL UNPACK_SAME_RANK(U%NR_NATURE,ZDG_3L(:,2),ZDG2_FULL)
CALL UNPACK_SAME_RANK(U%NR_NATURE,ZDG_3L(:,3),ZDG3_FULL)
WHERE ( ZDG2_FULL/=XUNDEF )
  ZDG_FULL = XFRAC_D2*ZDG2_FULL + XFRAC_D3*(ZDG3_FULL-ZDG2_FULL)
ELSEWHERE
  ZDG_FULL = XUNDEF
END WHERE
!
CALL UNPACK_SAME_RANK(U%NR_NATURE,ZWG_3L(:,2),ZWG2_FULL)
CALL UNPACK_SAME_RANK(U%NR_NATURE,ZWG_3L(:,3),ZWG3_FULL)
!
IF (KSTEP==1) THEN 
 ZWM (1:U%NDIM_FULL) = ZWG2_FULL(:)
 ZWIM(1:U%NDIM_FULL) = 0.0
 ALLOCATE(XRAINFALLCM(U%NSIZE_NATURE))
 ALLOCATE(XAVG_HORTCM(U%NSIZE_NATURE))
 XRAINFALLCM(:)=0.0
 XAVG_HORTCM(:)=0.0
ELSE
 ZWM (1:U%NDIM_FULL) = XWG_FULL(1:U%NDIM_FULL)
 ZWIM(1:U%NDIM_FULL) = XWGI_FULL(1:U%NDIM_FULL)
ENDIF
!
IF (KSTEP==1) THEN 
CALL UNPACK_SAME_RANK(U%NR_NATURE,ZWG_3L(:,2),XWG_FULL)
write(*,*)'XWG_FULL TSTP1',SUM(XWG_FULL,MASK=XWG_FULL/=XUNDEF)
!
! Full grid
!
  IF (LBUDGET_TOPD) CALL INIT_BUDGET_COUPL_ROUT(DEC, DC, DMI, PMESH_SIZE, IO, NP, NPE, U, U%NDIM_FULL)
 CALL ISBA_TO_TOPD(XWG_FULL,XWTOPT)
        !
ELSEIF (KSTEP==48) THEN 
        !
 ZRAND_MAP(:)=1.
 IF (LPERT_INIT) THEN
  CALL CREATE_RANDOM_MAP(UG,U,ZRAND_MAP)
  WRITE(*,*) 'BE CAREFUL PERTURBATION OF INITIAL SWI IS ACTIVATED' 
  WRITE(*,*) 'THE RANDOM COEFS ARE FROM ', MINVAL(ZRAND_MAP),'FROM ', MAXVAL(ZRAND_MAP) 
 ENDIF
  WHERE ( ZDG_FULL/=XUNDEF .AND. ZDG_FULL/=0. )
   XWG_FULL = (XFRAC_D2*(ZDG2_FULL/ZDG_FULL)*ZWG2_FULL + XFRAC_D3*((ZDG3_FULL-ZDG2_FULL)/ZDG_FULL)*ZWG3_FULL)*ZRAND_MAP
  ELSEWHERE
   XWG_FULL = XUNDEF
  END WHERE
  CALL ISBA_TO_TOPD(XWG_FULL,XWTOPT)
ELSE
 WHERE ( ZDG_FULL/=XUNDEF .AND. ZDG_FULL/=0. .AND.&
         ZWG2_FULL/=XUNDEF .AND. ZDG2_FULL/=XUNDEF.AND. XFRAC_D2/=XUNDEF)
  XWG_FULL = (XFRAC_D2*(ZDG2_FULL/ZDG_FULL)*ZWG2_FULL + XFRAC_D3*((ZDG3_FULL-ZDG2_FULL)/ZDG_FULL)*ZWG3_FULL)
 ELSEWHERE
  XWG_FULL = XUNDEF
 END WHERE
  CALL ISBA_TO_TOPD(XWG_FULL,XWTOPT)
ENDIF
!
!ludo prise en compte glace (pas de glace dans 3e couche)
CALL UNPACK_SAME_RANK(U%NR_NATURE,ZWGI_3L(:,2),ZWGI_FULL)
WHERE ( ZWGI_FULL/=XUNDEF .AND. XFRAC_D2>0 .AND. ZDG_FULL/=0. )
  XWGI_FULL = XFRAC_D2*(ZDG2_FULL/ZDG_FULL)*ZWGI_FULL
ELSEWHERE
  XWGI_FULL = XUNDEF
END WHERE
!
WHERE ( (XDTOPI/=XUNDEF).AND.(XWGI_FULL/=XUNDEF).AND.(ZWIM/=XUNDEF))
  ZRI_WGI = ( (XWGI_FULL - ZWIM)  ) * XDTOPI!old code
ELSEWHERE
  ZRI_WGI = 0.0
END WHERE
!
WHERE ( XDTOPI==XUNDEF ) 
  ZRI_WGI = 0.0
END WHERE
!
CALL ISBA_TO_TOPD(ZRI_WGI,ZRI_WGIT)
!
!!!!!!!!!!!!!!!!!
!Determination of Wsat, Wfc, Dmax
!!!!!!!!!!!!!!!
!test reservoir top=eau+glace -> pas de modif Wsat et Wfc
  Z_WSTOPI  = XWSTOPI 
  Z_WFCTOPI = XWFCTOPI
WHERE ( XWGI_FULL/=0. .AND.XWGI_FULL/=XUNDEF .AND. XWSTOPI/=0. )
  Z_WSTOPI  = XWSTOPI - XWGI_FULL
  Z_WFCTOPI = XWFCTOPI  * Z_WSTOPI / XWSTOPI
END WHERE
!ludo calcul en fct teneur glace
!
 CALL ISBA_TO_TOPD(Z_WSTOPI,XWSTOPT)
 CALL ISBA_TO_TOPD(Z_WFCTOPI,XWFCTOPT)
!
!ludo test empeche erreur num chgt phase
WHERE ( ABS(XWSTOPT-XWTOPT) < 0.0000000001 .AND. XWTOPT/=XUNDEF) XWSTOPT = XWTOPT
!
WHERE ( XWTOPT>XWSTOPT .AND. XWTOPT/=XUNDEF) XWTOPT = XWSTOPT
!
WHERE ( XWFCTOPT/= XUNDEF .AND. XWSTOPT/=XUNDEF .AND. XDTOPT/=XUNDEF)&
                XDMAXFC = (XWSTOPT - XWFCTOPT) * XDTOPT ! (m)
XDMAXT=XDMAXFC

!
!actualisation M
IF( IO%CKSAT=='EXP' .OR. IO%CKSAT=='SGH' ) THEN
  !ludo test
  XF_PARAM(:) = S%XF_PARAM(:)
  CALL UNPACK_SAME_RANK(U%NR_NATURE,XF_PARAM(:),ZF_PARAM_FULL)
  CALL ISBA_TO_TOPD(ZF_PARAM_FULL,ZF_PARAMT)
  !
  !passage de f a M (M=Wsat-Wfc/f)
  !ludo test ksat exp
  WHERE( ZF_PARAMT/=XUNDEF .AND. ZF_PARAMT/=0. ) ZF_PARAMT = (XWSTOPT-XWFCTOPT)/ZF_PARAMT
  !
  DO JJ=1,NNCAT
    XMPARA(JJ) = SUM(ZF_PARAMT(JJ,:),MASK=ZF_PARAMT(JJ,:)/=XUNDEF) / NNMC(JJ)
  ENDDO
  !
ELSE
  !
 DO JJ=1,NNCAT
  ZRANDOM(7)=1.
  IF (LPERT_PARAM)CALL READ_RANDOM_NUMBER(ZRANDOM)
  XMPARA(JJ) = (SUM( XDMAXFC(JJ,:),MASK=XDMAXFC(JJ,:)/=XUNDEF )/NNMC(JJ)/4.)*ZRANDOM(7)
 ENDDO
  !
ENDIF
!
!!!!!!!!!!!!!!!
!*       1.2    Water recharge 
!               ---------------
! Topodyn uses :
! - a water recharge = water added since last time step to compute hydrological similarity indexes
! - the total water content to compute a deficit
!
! This recharge is computed without regarding the changing of phase of water
! and the lateral transfers are performed regarding wsat et Wfc of last time step
!
!Full grid
!
WHERE ( (XDTOPI/=XUNDEF).AND.(XWG_FULL/=XUNDEF).AND.(ZWM/=XUNDEF).AND.(ZRI_WGI/=XUNDEF))
  ZRI = ( (XWG_FULL - ZWM)  ) * XDTOPI+ ZRI_WGI
ELSEWHERE
  ZRI = 0.0
ENDWHERE
!
! The water recharge on ISBA grid is computed on TOPMODEL grid
CALL RECHARGE_SURF_TOPD(ZRI,ZRT,U%NDIM_FULL)
!
!*       2.     Lateral distribution
!               --------------------
!*       2.1    Computation of local deficits on TOPODYN grid
!               ----------------------------------------
!
CALL TOPODYN_LAT(ZRT(:,:),ZDEFT(:,:),ZKAPPA(:,:),ZKAPPAC(:),GTOPD)
!
!*       2.2    Computation of contributive area on ISBA grid
!               ----------------------------------------
!
ZASI_F(:)=0.
ZAS_IBV_F(:,:)=0.
XAS_NATURE(:)=0.
XAS_IBV_P(:,:)=0.
CALL SAT_AREA_FRAC(ZDEFT,ZASI_F,ZAS_IBV_F,GTOPD)!work on full grid
!
!from full to nature grid
CALL PACK_SAME_RANK(U%NR_NATURE,NNPIX,INPIXI_P)
CALL PACK_SAME_RANK(U%NR_NATURE,ZASI_F,XAS_NATURE)

DO JJ=1,NNCAT
 CALL PACK_SAME_RANK(U%NR_NATURE,ZAS_IBV_F(:,JJ),XAS_IBV_P(:,JJ))
ENDDO
!
!*       2.3    Runoff from contributive area on ISBA grid
!               ----------------------------------------
!
XRUNOFF_IBV_P(:,:)=0.
DO JJ=1,NNCAT
 DO JI=1,U%NSIZE_NATURE
  IF (XAS_IBV_P(JI,JJ)/=XUNDEF    .AND. INPIXI_P(JI)/=0. .AND.&
      INBV_IN_MESH(JI,JJ)/=XUNDEF .AND. INPIXI_P(JI)/=XUNDEF .AND.&
      DEC%XRAINFALL(JI)/=XUNDEF .AND. XRAINFALLCM(JI)/=XUNDEF .AND.&
      DEC%XHORT(JI)/=XUNDEF .AND. XAVG_HORTCM(JI)/=XUNDEF) THEN
          ! kg/m2
      XRUNOFF_IBV_P(JI,JJ) = MAX(DEC%XRAINFALL(JI)-XRAINFALLCM(JI),0.0) * MAX(XAS_IBV_P(JI,JJ),0.0)+&
                            (DEC%XHORT(JI)-XAVG_HORTCM(JI))*INBV_IN_MESH(JI,JJ)/INPIXI_P(JI) 
  ELSE
 XRUNOFF_IBV_P(JI,JJ) = 0.0
  ENDIF
 ENDDO
ENDDO
XRAINFALLCM(1:U%NSIZE_NATURE) = DEC%XRAINFALL(1:U%NSIZE_NATURE)
XAVG_HORTCM(1:U%NSIZE_NATURE) = DEC%XHORT(1:U%NSIZE_NATURE)
!
!*       3.    Deficit (m) -> water storage (m3/m3) and changing of phase
!               ------------------------------------
!
! Full grid
DO JJ=1,NNCAT
  WHERE ( XDTOPT(JJ,:)/=XUNDEF .AND. XDTOPT(JJ,:)/=0. .AND. ZDEFT(JJ,:)/=XUNDEF)
    XWTOPT(JJ,:) = XWSTOPT(JJ,:) - ( ZDEFT(JJ,:) / XDTOPT(JJ,:) )      
   !changing phase
    XWTOPT(JJ,:) = XWTOPT(JJ,:) - ZRI_WGIT(JJ,:)
  END WHERE
ENDDO
!'
!
!*       3.    Deficit (m) -> water storage (m3/m3) and changing of phase
!               ------------------------------------
!
! Full grid
DO JJ=1,NNCAT
  WHERE ( XDTOPT(JJ,:)/=XUNDEF .AND. XDTOPT(JJ,:)/=0. .AND. ZDEFT(JJ,:)/=XUNDEF)
    XWTOPT(JJ,:) = XWSTOPT(JJ,:) - ( ZDEFT(JJ,:) / XDTOPT(JJ,:) )      
   !changing phase
    XWTOPT(JJ,:) = XWTOPT(JJ,:) - ZRI_WGIT(JJ,:)
  END WHERE
ENDDO
!*       4.     TOPODYN => ISBA
!               ---------------
!*       4.1    Calculation of water storage on ISBA grid
!               -----------------------------------------
!
CALL TOPD_TO_ISBA(K, PEK,UG, U, U%NDIM_FULL,KSTEP,GTOPD)!=modif of XWG_FULL from XWTOPT
! Nature grid
CALL PACK_SAME_RANK(U%NR_NATURE, (1-XFRAC_D2)*ZWG2_FULL + XFRAC_D2*XWG_FULL, ZWG_3L(:,2))
CALL PACK_SAME_RANK(U%NR_NATURE, (1-XFRAC_D3)*ZWG3_FULL + XFRAC_D3*XWG_FULL, ZWG_3L(:,3))
!
!*       4.2    Budget correction
!  -----------------------------------------
!
 CALL PACK_SAME_RANK(U%NR_NATURE,UG%G%XMESH_SIZE,ZMESH_SIZE)
 ZAVG_MESH_SIZE = SUM(ZMESH_SIZE(:),MASK=ZMESH_SIZE(:)/=XUNDEF) / COUNT(ZMESH_SIZE(:)/=XUNDEF)
!
IF (IO%CISBA=='DIF') THEN
 CALL TOPD_TO_DF(IO, NK, NP, NPE, ZWG_3L)
ELSEIF (IO%CISBA=='3-L') THEN
 CALL DISPATCH_WG(S, NP, NPE, ZWG_3L, ZWGI_3L, ZDG_3L)
ENDIF
!
!
IACT_GROUND_LAYER=3

IF (IO%CISBA=='DIF') THEN

DO JL=2,IO%NGROUND_LAYER
  ISUM = 0
  DO JP = 1,IO%NPATCH
    IF (ALL(NPE%AL(JP)%XWG(:,JL)==XUNDEF)) ISUM = ISUM + 1
  ENDDO
  IF (ISUM==IO%NPATCH) THEN
         IACT_GROUND_LAYER=JL-1
         !WRITE(ILUOUT,*) 'IACT_GROUND_LAYER=',IACT_GROUND_LAYER
         EXIT
 ENDIF
ENDDO
!
ENDIF
!
ZWOVSAT_IBV_P(:,:)=0.
CALL PACK_SAME_RANK(U%NR_NATURE,Z_WSTOPI,ZWSAT)
!!!!! Budget ctrl on layer 2
ZWG(:,:) = 0.
ZDG(:,:) = 0.
DO JP = 1,IO%NPATCH
  DO JJ = 1,NP%AL(JP)%NSIZE_P
    IMASK = NP%AL(JP)%NR_P(JJ)
    ZWG(IMASK,JP) = NPE%AL(JP)%XWG(JJ,2)
    ZDG(IMASK,JP) = NP%AL(JP)%XDG(JJ,2)
  ENDDO
ENDDO
CALL CONTROL_WATER_BUDGET_TOPD(IO, S, U, ZWG_TMP(:,2,:), ZWG, ZDG,&
           ZMESH_SIZE,ZAVG_MESH_SIZE,ZWSAT(:),ZWOVSAT_IBV_P(:,:))
DO JP = 1,IO%NPATCH
  DO JJ = 1,NP%AL(JP)%NSIZE_P
    IMASK = NP%AL(JP)%NR_P(JJ)
    NPE%AL(JP)%XWG(JJ,2) = ZWG(IMASK,JP)
  ENDDO
ENDDO   
!
!!!!! Budget ctrl on layer JL 
DO JL = 3,IACT_GROUND_LAYER
  ZWG(:,:) = 0.
  ZDG(:,:) = 0.
  DO JP = 1,IO%NPATCH
    DO JJ = 1,NP%AL(JP)%NSIZE_P
      IMASK = NP%AL(JP)%NR_P(JJ)
      ZWG(IMASK,JP) = NPE%AL(JP)%XWG(JJ,JL)
      ZDG(IMASK,JP) = NP%AL(JP)%XDG(JJ,JL)-NP%AL(JP)%XDG(JJ,JL-1)
    ENDDO
  ENDDO
 CALL CONTROL_WATER_BUDGET_TOPD(IO, S, U, ZWG_TMP(:,JL,:), ZWG, ZDG, &
           ZMESH_SIZE,ZAVG_MESH_SIZE,ZWSAT(:),ZWOVSAT_IBV_P(:,:))
  DO JP = 1,IO%NPATCH
    DO JJ = 1,NP%AL(JP)%NSIZE_P
      IMASK = NP%AL(JP)%NR_P(JJ)
      NPE%AL(JP)%XWG(JJ,JL) = ZWG(IMASK,JP)
    ENDDO
  ENDDO          

ENDDO
!
DO JP = 1,IO%NPATCH
  WHERE(NPE%AL(JP)%XWG(:,:)>ZWSATMAX.AND.NPE%AL(JP)%XWG(:,:)/=XUNDEF)
   NPE%AL(JP)%XWG(:,:)=ZWSATMAX
  ENDWHERE
  WHERE(NPE%AL(JP)%XWG(:,:)<XWGMIN)
   NPE%AL(JP)%XWG(:,:)=XWGMIN
  ENDWHERE
ENDDO
!
!
!*      5.0    Total discharge
!              ---------------
!
!*      5.1    Total water for runoff on TOPODYN grid
!              ---------------------------------------
!
!
! Full grid
!
ZWOVSATI_F(:)=0.
XWOVSATI_P(:)=XWOVSATI_P(:)
CALL UNPACK_SAME_RANK(U%NR_NATURE,XWOVSATI_P(:),ZWOVSATI_F(:))!from isba grid to full grid
 CALL UNPACK_SAME_RANK(NP%AL(1)%NR_P,NP%AL(1)%XRUNOFFD(:),ZRUNOFFD_NAT)
 CALL UNPACK_SAME_RANK(U%NR_NATURE,ZRUNOFFD_NAT,ZRUNOFFD_F)!from isba grid to full grid

ZRUNOFF_IBV_F(:,:)=0.
DO JJ=1,NNCAT
  XRUNOFF_IBV_P(:,JJ)=XRUNOFF_IBV_P(:,JJ)+ZWOVSAT_IBV_P(:,JJ)*ZRUNOFFD_NAT(:)*XRHOLW !isba grid !kg/m2
  CALL UNPACK_SAME_RANK(U%NR_NATURE,XRUNOFF_IBV_P(:,JJ),ZRUNOFF_IBV_F(:,JJ))! from isba grid to full grid
!
 DO JI=1,U%NDIM_FULL! full_grid
   IF(ZWOVSATI_F(JI)/=XUNDEF.AND.XAIBV_F(JI,JJ)/=XUNDEF.AND.ZWOVSATI_F(JI)>=0.)THEN
   ZRUNOFF_IBV_F(JI,JJ)=ZRUNOFF_IBV_F(JI,JJ)+ZWOVSATI_F(JI)* XAIBV_F(JI,JJ)*ZRUNOFFD_F(JI)*XRHOLW!kg/m2

   ENDIF
 ENDDO
ENDDO
!
!
ZRUNOFF_TOPD(:,:) = 0.0
!
!
DO JJ=1,NNCAT
 WHERE (ZRUNOFF_IBV_F(:,JJ)/=XUNDEF.AND.UG%G%XMESH_SIZE(:)/=XUNDEF)
  ZRUNOFF_IBV_F(:,JJ)=ZRUNOFF_IBV_F(:,JJ)*UG%G%XMESH_SIZE(:)/XRHOLW/3600.!! from kg/m2 to m3/s
 ENDWHERE
ENDDO
CALL ISBA_TO_TOPDSAT(XKA_PRE,XKAC_PRE,U%NDIM_FULL,ZRUNOFF_IBV_F,ZRUNOFF_TOPD)
!
!
!
!*      5.2    Total water for drainage on TOPODYN grid
!              ----------------------------------------
!In XAVG_DRAINC, the paches have been average
CALL UNPACK_SAME_RANK(U%NR_NATURE,DEC%XDRAIN,ZDRAINC_FULL)
CALL UNPACK_SAME_RANK(U%NR_NATURE,XAVG_DRAINCM,ZDRAINC_FULLM)
!
CALL DIAG_ISBA_TO_ROUT(UG%G%XMESH_SIZE,ZDRAINC_FULL,ZDRAINC_FULLM,ZDRAIN_ISBA)
!
WHERE (ZDRAIN_ISBA==XUNDEF) ZDRAIN_ISBA=0.
!
XAVG_DRAINCM(:)  = DEC%XDRAIN(:)
!
ZDRAIN_TOPD(:,:) = 0.0
ZDRAIN_ISBA=ZDRAIN_ISBA*XATOP
!
CALL ISBA_TO_TOPD(ZDRAIN_ISBA,ZDRAIN_TOPD)
!
DO JJ=1,NNCAT
  DO JI=1,NNMC(JJ)
    IF (NMASKT(JJ,JI)/=NUNDEF) &
      ZDRAIN_TOPD(JJ,JI) = ZDRAIN_TOPD(JJ,JI) / NNPIX(NMASKT(JJ,JI))
  ENDDO
ENDDO   
!
!*      6    Routing (runoff + drainage + exfiltration)
!
CALL ROUTING(ZRUNOFF_TOPD,ZDRAIN_TOPD,KSTEP)
!
XKA_PRE(:,:) = ZKAPPA(:,:)
XKAC_PRE(:) = ZKAPPAC(:)
!!
!*      7.0    Computing Alert levels
!              ----------------------------

!*      8.0    Writing results in map files
!              ----------------------------
!
LWORK=.FALSE.
IF (NFREQ_MAPS_ASAT/=0) THEN
   LWORK=MOD(KSTEP,NFREQ_MAPS_ASAT)==0
ENDIF
!
IF (LWORK) THEN
  CALL OPEN_FILE('ASCII ',NUNIT,HFILE='carte_surfcont'//HSTEP,HFORM='FORMATTED',HACTION='WRITE')
  CALL WRITE_FILE_ISBAMAP(UG, &
                          NUNIT,ZASI_F,U%NDIM_FULL)
  CALL CLOSE_FILE('ASCII ',NUNIT)
ENDIF
!
IF (LHOOK) CALL DR_HOOK('COUPL_TOPD',1,ZHOOK_HANDLE)
!
!
END SUBROUTINE COUPL_TOPD
