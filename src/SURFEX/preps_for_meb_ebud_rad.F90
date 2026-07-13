!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!SFX_LIC for details. version 1.
!   ############################################################################
SUBROUTINE PREPS_FOR_MEB_EBUD_RAD(IO, PPS,                                     &
     PLAICV,PSNOWRHO,PSNOWSWE,PSNOWHEAT,PSNOWLIQ,                              &
     PSNOWTEMP,PSNOWDZ,PSCOND,PHEATCAPS,PEMISNOW,PSIGMA_F,PCHIP,               &
     PTSTEP,PSR,PTA,PVMOD,PSNOWAGE,PSNOWDIAMOPT,PSNOWSPHERI,PSNOWHIST,         &
     PPERMSNOWFRAC,HSNOW_ISBA, PUREF, PZ0EFF          )
!   ############################################################################
!
!!****  *PREPS_FOR_MEB_EBUD_RAD*
!!
!!    PURPOSE
!!    -------
!
!     Get preliminary estimates of certain parameters needed for energy budget
!     solution of snowpack, and some other misc inputs needed by radiation
!     routines for MEB.
!
!!**  METHOD
!!    ------
!
!
!!    EXTERNAL
!!    --------
!!
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
!!    AUTHOR
!!    ------
!!
!!    A. Boone                * CNRM-GAME, Meteo-France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    02/2011
!-------------------------------------------------------------------------------
!
!*       0.     DECLARATIONS
!               ------------
!
USE MODD_ISBA_OPTIONS_n,   ONLY : ISBA_OPTIONS_t
!
USE MODD_SNOW_PAR,            ONLY : XRHOSMAX_ES, XRHOSMIN_ES, XEMISSN, XSNOWDMIN, &
                                     XSNOWTHRMCOND1
!
USE MODD_CSTS,                ONLY : XTT, XLMTT, XRHOLW, XCI
!
USE MODD_SURF_PAR,            ONLY : XUNDEF
!
USE MODD_SNOW_METAMO,         ONLY : XSNOWDZMIN, XUEPSI
!
USE MODE_SNOW3L,              ONLY : SNOW3LTHRM, SNOW3LSCAP, SNOW3LFALL,         &
                                     SNOW3LTRANSF, SNOW3LGRID, SNOW3LCOMPACTN,   &
                                     SNOWCROTHRM
USE MODE_SNOWCRO,             ONLY:  SNOWNLFALL_UPGRID, SNOWNLGRIDFRESH_1D
!
USE MODE_MEB,                 ONLY : MEB_SHIELD_FACTOR
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!
!*      0.1    Declaration of Arguments
!
REAL                                :: PTSTEP   ! time step (s)
TYPE(ISBA_OPTIONS_t), INTENT(INOUT) :: IO
REAL, DIMENSION(:),   INTENT(IN)    :: PLAICV
REAL, DIMENSION(:),   INTENT(IN)    :: PPS
REAL, DIMENSION(:),   INTENT(IN)    :: PSR
REAL, DIMENSION(:),   INTENT(IN)    :: PTA
REAL, DIMENSION(:),   INTENT(IN)    :: PVMOD
REAL, DIMENSION(:),   INTENT(IN)    :: PPERMSNOWFRAC
REAL, DIMENSION(:,:), INTENT(IN)    :: PSNOWHEAT
REAL, DIMENSION(:),   INTENT(IN)    :: PUREF, PZ0EFF
CHARACTER(LEN=*),     INTENT(IN)    :: HSNOW_ISBA
REAL, DIMENSION(:,:), INTENT(INOUT) :: PSNOWSWE, PSNOWAGE, PSNOWRHO
REAL, DIMENSION(:,:), INTENT(INOUT) :: PSNOWDIAMOPT, PSNOWSPHERI, PSNOWHIST ! crocus
!
REAL, DIMENSION(:),   INTENT(OUT)   :: PSIGMA_F, PCHIP
REAL, DIMENSION(:),   INTENT(OUT)   :: PEMISNOW
REAL, DIMENSION(:,:), INTENT(OUT)   :: PSNOWDZ, PSCOND, PHEATCAPS, PSNOWTEMP, PSNOWLIQ
!
!
!*      0.2    declarations of local variables
!
INTEGER                                            :: JI, JK, JJ, INLVLS, ISIZE_SNOW, INI
INTEGER, DIMENSION(SIZE(PTA))                      :: NMASK      ! indices correspondance between arrays
REAL, DIMENSION(SIZE(PLAICV,1))                    :: ZPSNA
REAL, DIMENSION(SIZE(PTA))                         :: ZSNOW, ZSNOWFALL
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!----------------------------------------------------
! 0) Initialization
!
IF (LHOOK) CALL DR_HOOK('PREPS_FOR_MEB_EBUD_RAD',0,ZHOOK_HANDLE)
!
INI             = SIZE(PSNOWRHO,1)
INLVLS          = SIZE(PSNOWRHO,2)
!
! Initialize some output variables:
! where snow depth below threshold, set to non-snow values
!
PSNOWTEMP(:,:)      = XTT
PSNOWLIQ (:,:)      = 0.0
PSCOND   (:,:)      = XSNOWTHRMCOND1
PHEATCAPS(:,:)      = XRHOSMIN_ES*XCI
!
! Test variables to check for existance of snow:
!
WHERE(PSNOWRHO(:,:)==XUNDEF)
   PSNOWDZ(:,:)     = 0.
   PSNOWRHO(:,:) = XRHOSMIN_ES
! The previous initialization of density is necessary
! at the first time step with snow on the ground in case of
! divergence of input between MEB and the snow scheme. The problem
! should be reduced after this version which implements here Crocus
! discretization routines, but it might still occur again if combined
! with blowing snow or machine made snow.

ELSEWHERE
   PSNOWDZ(:,:)     = PSNOWSWE(:,:)/PSNOWRHO(:,:)
END WHERE
!
ZSNOWFALL(:)     = PSR(:)*PTSTEP/XRHOSMAX_ES
!
ZSNOW(:)         = 0.0
DO JK=1,INLVLS
   DO JI=1,INI
      ZSNOW(JI)  = ZSNOW(JI) + PSNOWDZ(JI,JK)
   ENDDO
ENDDO
!
! Here, as in snow3l (ISBA-ES), we account for several processes
! on the snowpack before surface energy budget computations
! (i.e. snowfall on albedo, density, thickness, and compaction etc...)

! ===============================================================
! === Packing: Only call snow model routines when there is snow on the surface
!              exceeding a minimum threshold OR if the equivalent
!              snow depth falling during the current time step exceeds 
!              this limit.
!
! counts the number of points where the computations will be made
!
!
ISIZE_SNOW = 0
NMASK(:)   = 0
!
DO JJ=1,INI
   IF (ZSNOW(JJ) >= XSNOWDMIN .OR. ZSNOWFALL(JJ) >= XSNOWDMIN) THEN
      ISIZE_SNOW = ISIZE_SNOW + 1
      NMASK(ISIZE_SNOW) = JJ
   ENDIF
ENDDO
!
IF (ISIZE_SNOW>0) THEN
   CALL CALL_SNOW_ROUTINES(ISIZE_SNOW,INLVLS,NMASK)
ENDIF
!
! ===============================================================
!
!
! View factor: (1 - shielding factor)
!
ZPSNA(:)          = 0.
PCHIP(:)          = MEB_SHIELD_FACTOR(PLAICV,ZPSNA)
PSIGMA_F(:)       = 1.0 - PCHIP(:)
!
! snow emissivity
!
PEMISNOW(:)       = XEMISSN
!
IF (LHOOK) CALL DR_HOOK('PREPS_FOR_MEB_EBUD_RAD',1,ZHOOK_HANDLE)
!
!
CONTAINS
!================================================================
SUBROUTINE CALL_SNOW_ROUTINES(KSIZE1,KSIZE2,KMASK)
!
! Make some snow computations only over regions with snow cover or snow falling
!
IMPLICIT NONE
!
INTEGER,               INTENT(IN) :: KSIZE1
INTEGER,               INTENT(IN) :: KSIZE2
INTEGER, DIMENSION(:), INTENT(IN) :: KMASK
!
REAL, DIMENSION(KSIZE1,KSIZE2) :: ZP_SNOWSWE
REAL, DIMENSION(KSIZE1,KSIZE2) :: ZP_SNOWRHO
REAL, DIMENSION(KSIZE1,KSIZE2) :: ZP_SNOWHEAT
REAL, DIMENSION(KSIZE1,KSIZE2) :: ZP_SNOWTEMP
REAL, DIMENSION(KSIZE1,KSIZE2) :: ZP_SNOWLIQ
REAL, DIMENSION(KSIZE1,KSIZE2) :: ZP_SNOWDZ
REAL, DIMENSION(KSIZE1,KSIZE2) :: ZP_SCOND
REAL, DIMENSION(KSIZE1,KSIZE2) :: ZP_SNOWAGE
REAL, DIMENSION(KSIZE1,KSIZE2) :: ZP_SNOWDIAMOPT, ZP_SNOWSPHERI, ZP_SNOWHIST
REAL, DIMENSION(KSIZE1,KSIZE2) :: ZP_SNOWDZO, ZP_SNOWDZN
REAL, DIMENSION(KSIZE1,KSIZE2) :: ZP_HEATCAPS
REAL, DIMENSION(KSIZE1)        :: ZP_SNOW
REAL, DIMENSION(KSIZE1)        :: ZP_SNOWHMASS
REAL, DIMENSION(KSIZE1)        :: ZP_PERMSNOWFRAC
REAL, DIMENSION(KSIZE1)        :: ZP_PS
REAL, DIMENSION(KSIZE1)        :: ZP_SR
REAL, DIMENSION(KSIZE1)        :: ZP_TA
REAL, DIMENSION(KSIZE1)        :: ZP_VMOD
REAL, DIMENSION(KSIZE1)        :: ZP_WORK
REAL, DIMENSION(KSIZE1)        :: ZP_UNLOAD
REAL, DIMENSION(KSIZE1)        :: ZP_SNOWMASSNEW
!
! For Crocus discretization routines
LOGICAL, DIMENSION(KSIZE1)     :: GP_SNOWFALL
INTEGER, DIMENSION(KSIZE1)     :: INLVLS_USE, INLVLS_USE_OLD
INTEGER :: IMAX_USE
! Fresh snow characteristics
REAL, DIMENSION(KSIZE1)        :: ZP_SNOWRHOF, ZP_SNOWDZF
REAL, DIMENSION(KSIZE1)        :: ZP_SNOWDIAMOPTF, ZP_SNOWSPHERIF, ZP_SNOWHISTF
REAL, DIMENSION(KSIZE1)        :: ZP_SNOWAGEF
REAL, DIMENSION(KSIZE1,1)      :: ZP_WETCOEF
REAL, DIMENSION(KSIZE1,1)      :: ZP_SNOWIMPURF
REAL, DIMENSION(KSIZE1,KSIZE2,1):: ZP_SNOWIMPUR
REAL, DIMENSION(KSIZE1)        :: ZP_SNOWMAK
REAL, DIMENSION(KSIZE1)        :: ZP_SNOWALB
REAL, DIMENSION(KSIZE1, 4)     :: ZP_BLOWSNW
!
LOGICAL, DIMENSION(KSIZE1)     :: GP_MODIF_GRID
!
REAL, DIMENSION(KSIZE1)        :: ZP_UREF, ZP_Z0EFF
!
LOGICAL :: GSUCCESS ! Flag to test the success of Crocus regridding
!
INTEGER         :: JWRK, JJ, JI
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!----------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('PREPS_FOR_MEB_EBUD_RAD:CALL_SNOW_ROUTINES',0,ZHOOK_HANDLE)
!
! pack the variables
!
DO JWRK=1,KSIZE2
   DO JJ=1,KSIZE1
      JI = KMASK(JJ)
      ZP_SNOWSWE (JJ,JWRK) = PSNOWSWE (JI,JWRK)
      ZP_SNOWRHO (JJ,JWRK) = PSNOWRHO (JI,JWRK)
      ZP_SNOWHEAT(JJ,JWRK) = PSNOWHEAT(JI,JWRK)
      ZP_SNOWAGE (JJ,JWRK) = PSNOWAGE (JI,JWRK)
      ZP_SNOWDZ  (JJ,JWRK) = PSNOWDZ  (JI,JWRK)
   ENDDO
ENDDO
!
IF (HSNOW_ISBA =="CRO") THEN
  DO JWRK=1,KSIZE2
    DO JJ=1,KSIZE1
      JI = KMASK(JJ)
      ZP_SNOWDIAMOPT (JJ,JWRK) = PSNOWDIAMOPT (JI,JWRK)
      ZP_SNOWSPHERI  (JJ,JWRK) = PSNOWSPHERI  (JI,JWRK)
      ZP_SNOWHIST    (JJ,JWRK) = PSNOWHIST    (JI,JWRK)
    ENDDO
  ENDDO
ELSE
  ! For MEB 3-L compatibility
  DO JWRK=1,KSIZE2
    DO JJ=1,KSIZE1
      ZP_SNOWDIAMOPT (JJ,JWRK) = XUNDEF
      ZP_SNOWSPHERI  (JJ,JWRK) = XUNDEF
      ZP_SNOWHIST    (JJ,JWRK) = XUNDEF
    ENDDO
  ENDDO
ENDIF
!
DO JJ=1,KSIZE1
   JI = KMASK(JJ)
   ZP_SNOW        (JJ) = ZSNOW(JI)
   ZP_PS          (JJ) = PPS  (JI)
   ZP_SR          (JJ) = PSR  (JI)
   ZP_TA          (JJ) = PTA  (JI)
   ZP_VMOD        (JJ) = PVMOD(JI)
   ZP_PERMSNOWFRAC(JJ) = PPERMSNOWFRAC(JI)
   ZP_UNLOAD      (JJ) = 0.0 ! Currently set to zero
   ZP_UREF        (JJ) = PUREF (JI)
   ZP_Z0EFF       (JJ) = PZ0EFF(JI)
ENDDO
!
!---------------------------------------------------------------
!
! Local working:
!
!
IF (HSNOW_ISBA /="CRO") THEN
  !
   WHERE(ZP_SNOWHEAT(:,:) == XUNDEF) ! when snow first falls, possibly
                                     ! remove FLAG values and initialize:
      ZP_SNOWHEAT(:,:)   = 0.
      ZP_SNOWAGE(:,:)    = 0.
      ZP_SNOWDZ(:,:)     = 0.
      ZP_SNOWSWE(:,:)    = 0.
      ZP_SNOWRHO(:,:)    = XRHOSMIN_ES
   ELSEWHERE
      ZP_SNOWHEAT(:,:)   = ZP_SNOWHEAT(:,:)*ZP_SNOWDZ(:,:) ! J/m3 to J/m2
   END WHERE
  !
   CALL SNOW3LFALL(PTSTEP, ZP_SR, ZP_TA, ZP_VMOD, ZP_SNOW, ZP_SNOWRHO, ZP_SNOWDZ,  &
                  ZP_SNOWHEAT, ZP_SNOWHMASS, ZP_WORK, ZP_SNOWAGE, ZP_PERMSNOWFRAC, &
                  ZP_UNLOAD, ZP_SNOWMASSNEW)
  !
  CALL SNOW3LGRID(ZP_SNOWDZN,ZP_SNOW,PSNOWDZ_OLD=ZP_SNOWDZ)
  !
  CALL SNOW3LTRANSF(ZP_SNOW,ZP_SNOWDZ,ZP_SNOWDZN,ZP_SNOWRHO,ZP_SNOWHEAT,ZP_SNOWAGE)
  !
  ! NOTE, in ISBA-ES, the number of active snow layers is fixed:
  !
  INLVLS_USE(:) = KSIZE2
  !
ELSE
  ZP_SNOWHEAT(:,:)   = ZP_SNOWHEAT(:,:)*ZP_SNOWDZ(:,:) ! J/m3 to J/m2
  !
  ! don't care about that for now
  ZP_WETCOEF = 0.
  ZP_BLOWSNW(:,:) = 0.
  ZP_SNOWMAK = 0.
  ZP_SNOWIMPUR (:,:,:) = 0.
  ! number of active layers
  INLVLS_USE(:) = 0
  DO JWRK = 1,KSIZE2
    DO JJ = 1,KSIZE1
      IF ( ZP_SNOWSWE(JJ,JWRK)>0. ) THEN
        INLVLS_USE(JJ) = JWRK
      ENDIF
    ENDDO  !  end loop snow layers
  ENDDO    ! end loop grid points
  !
  !
  IMAX_USE = MAXVAL(INLVLS_USE)
  INLVLS_USE_OLD = INLVLS_USE
  ZP_SNOWDZO(:,:) = ZP_SNOWDZ(:,:)
  !
  CALL SNOWNLFALL_UPGRID(PTSTEP,ZP_SR,ZP_TA,ZP_VMOD,ZP_SNOW, ZP_SNOWRHO,                   &
                         ZP_SNOWDZ, ZP_SNOWHEAT,ZP_SNOWHMASS,ZP_SNOWALB,ZP_PERMSNOWFRAC,   &
                         ZP_SNOWDIAMOPT,ZP_SNOWSPHERI,ZP_SNOWHIST,ZP_SNOWAGE,GP_SNOWFALL,  &
                         ZP_SNOWDZN, ZP_SNOWRHOF, ZP_SNOWDZF,                              &
                         ZP_SNOWDIAMOPTF, ZP_SNOWSPHERIF, ZP_SNOWHISTF,                    &
                         ZP_SNOWAGEF,ZP_WETCOEF, ZP_SNOWIMPURF,GP_MODIF_GRID,INLVLS_USE,   &
                         IO%CSNOWDRIFT,IO%CSNOWFPAPPUS,ZP_Z0EFF,ZP_UREF,                   &
                         ZP_BLOWSNW, IO%CSNOWFALL,                                         &
                         ZP_SNOWMAK, .FALSE., .FALSE.,IMAX_USE)
  !
  ! Update grid/discretization
  ! Reset grid to conform to Crocus specifications:
  !
  DO JJ=1,KSIZE1
    !
    IF ( GP_MODIF_GRID(JJ) ) THEN
      CALL SNOWNLGRIDFRESH_1D(JJ,ZP_SNOW(JJ),ZP_SNOWDZ(JJ,:),ZP_SNOWDZN(JJ,:),ZP_SNOWRHO(JJ,:),                       &
                              ZP_SNOWHEAT(JJ,:),ZP_SNOWDIAMOPT(JJ,:),ZP_SNOWSPHERI(JJ,:),                             &
                              ZP_SNOWHIST(JJ,:),ZP_SNOWAGE(JJ,:),ZP_SNOWIMPUR(JJ,:,:),GP_SNOWFALL(JJ),ZP_SNOWRHOF(JJ),&
                              ZP_SNOWDZF(JJ),ZP_SNOWHMASS(JJ),ZP_SNOWDIAMOPTF(JJ),ZP_SNOWSPHERIF(JJ),                 &
                              ZP_SNOWHISTF(JJ),ZP_SNOWAGEF(JJ),ZP_SNOWIMPURF(JJ,:),INLVLS_USE(JJ), INLVLS_USE_OLD(JJ),&
                              GSUCCESS)
      ! In case of troubles, this could help
      !  IF (.NOT. GSUCCESS) THEN
      !    PRINT*, 'regridding problem in preps_for_meb_ebud_rad'
      !    PRINT*, 'JJ=', JJ
      !    PRINT*, 'ZP_SNOWDZO=',ZP_SNOWDZO(JJ,:)
      !    PRINT*, 'ZP_SNOWDZN=',ZP_SNOWDZN(JJ,:)
      !  END IF                 
    ENDIF
    ! To avoid tests below for inactive layers and avoid division by 0, initialize density for inactive layers
    ZP_SNOWRHO(JJ, INLVLS_USE(JJ)+1:KSIZE2) = XRHOSMIN_ES
    !
  ENDDO
END IF
!
! Snow heat capacity (J m-3 K-1)
ZP_HEATCAPS(:,:)   = SNOW3LSCAP(ZP_SNOWRHO)
!
! initialisations:
ZP_SNOWTEMP(:,:) = XTT
ZP_SNOWLIQ(:,:) = 0.0
ZP_SNOWSWE(:,:) = 0
!
DO JJ=1,KSIZE1
  ! Snow temperature (K)
  ZP_SNOWTEMP(JJ,1:INLVLS_USE(JJ)) = XTT + &
                                    ( ((ZP_SNOWHEAT(JJ,1:INLVLS_USE(JJ))/MAX(1.E-10,ZP_SNOWDZ(JJ,1:INLVLS_USE(JJ)))) &
                                    + XLMTT*ZP_SNOWRHO(JJ,1:INLVLS_USE(JJ)))/ZP_HEATCAPS(JJ,1:INLVLS_USE(JJ)) )
  !
  ZP_SNOWLIQ(JJ,1:INLVLS_USE(JJ)) = MAX(0.0,ZP_SNOWTEMP(JJ,1:INLVLS_USE(JJ))-XTT)*ZP_HEATCAPS(JJ,1:INLVLS_USE(JJ)) * &
                                   ZP_SNOWDZ(JJ,1:INLVLS_USE(JJ))/(XLMTT*XRHOLW) 
  !
  ! SWE:
  ZP_SNOWSWE(JJ,1:INLVLS_USE(JJ)) = ZP_SNOWDZ(JJ,1:INLVLS_USE(JJ))*ZP_SNOWRHO(JJ,1:INLVLS_USE(JJ))
END DO
!
ZP_SNOWTEMP(:,:)   = MIN(XTT,ZP_SNOWTEMP(:,:))
!
IF(HSNOW_ISBA=="CRO") THEN
  !
  ! à voir ce qu'on garde ici
  CALL SNOWCROTHRM(ZP_SNOWRHO, ZP_SCOND, ZP_SNOWTEMP, ZP_PS, ZP_SNOWLIQ, IO%CSNOWCOND)
ELSE
  !
  CALL SNOW3LCOMPACTN(PTSTEP, XSNOWDZMIN, ZP_SNOWRHO, ZP_SNOWDZ, ZP_SNOWTEMP, ZP_SNOW, ZP_SNOWLIQ)
  !
  ! Snow thermal conductivity:
  !
  CALL SNOW3LTHRM(ZP_SNOWRHO, ZP_SCOND, ZP_SNOWTEMP, ZP_PS)
  !
ENDIF
!
!----------------------------------------------------------------
!
! Unpack:
!
DO JWRK=1,KSIZE2
   DO JJ=1,KSIZE1
      JI = KMASK(JJ)
      PSNOWSWE (JI,JWRK) = ZP_SNOWSWE (JJ,JWRK)
      PSNOWRHO (JI,JWRK) = ZP_SNOWRHO (JJ,JWRK)
      PSNOWAGE (JI,JWRK) = ZP_SNOWAGE (JJ,JWRK)
      PSNOWDZ  (JI,JWRK) = ZP_SNOWDZ  (JJ,JWRK)
      PSNOWTEMP(JI,JWRK) = ZP_SNOWTEMP(JJ,JWRK)
      PSNOWLIQ (JI,JWRK) = ZP_SNOWLIQ (JJ,JWRK)
      PSCOND   (JI,JWRK) = ZP_SCOND   (JJ,JWRK)
      PHEATCAPS(JI,JWRK) = ZP_HEATCAPS(JJ,JWRK)
   ENDDO
ENDDO
!
IF (LHOOK) CALL DR_HOOK('PREPS_FOR_MEB_EBUD_RAD:CALL_SNOW_ROUTINES',1,ZHOOK_HANDLE)
!
END SUBROUTINE CALL_SNOW_ROUTINES
!================================================================  
!
END SUBROUTINE PREPS_FOR_MEB_EBUD_RAD
