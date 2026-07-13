!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!SFX_LIC for details. version 1.
!     #########
SUBROUTINE SNOW3L_ISBA(IO, G, PK, PEK, DK, DEK, DMK, OMEB, HIMPLICIT_WIND,        &
                       TPTIME, PTSTEP, PTIMEC, PVEGTYPE, PTG, PCT, PSOILHCAPZ,    &
                       PSOILCONDZ, PPS, PTAR, PTAC, PSW_RAD, PQA, PVMOD, PVDIR,   &
                       PLW_RAD, PRR, PSR, PRHOA, PUREF,PEXNS, PEXNA,              &
                       PDIRCOSZW, PSLOPEDIR, PZREF, PALB, PD_G, PDZG,             &
                       PPEW_A_COEF, PPEW_B_COEF, PPET_A_COEF, PPEQ_A_COEF,        &
                       PPET_B_COEF, PPEQ_B_COEF, PTHRUFAL, PGRNDFLUX, PFLSN_COR,  &
                       PEVAPCOR, PLES3L, PLEL3L, PEVAP, PSNOWSFCH, PRI,PZENITH,   &
                       PAZIM, PQS, NPAR_VEG_IRR_USE, KTAB_SYT,                    &
                       P_DIR_SW, P_SCA_SW, PIMPWET, PIMPDRY,                      &
                       PBLOWSNW_FLUX, PBLOWSNW_CONC                               )
!     ######################################################################################
!
!!****  *SNOW3L_ISBA*
!!
!!    PURPOSE
!!    -------
!
!     3-Layer snow scheme option (Boone and Etchevers 1999)
!     This routine is NOT called as snow depth goes below
!     a critical threshold which is vanishingly small.
!     This routine acts as an interface between SNOW3L and ISBA.
!
!!**  METHOD
!!    ------
!
!     Direct calculation
!
!!    EXTERNAL
!!    --------
!
!     None
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!
!!    REFERENCE
!!    ---------
!!
!!    Boone and Etchevers (1999)
!!    Belair (1995)
!!    Noilhan and Planton (1989)
!!    Noilhan and Mahfouf (1996)
!!
!!    AUTHOR
!!    ------
!!	A. Boone           * Meteo-France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original        7/99  Boone
!!      Packing added   4/00  Masson & Boone
!!      z0h and snow    2/06  LeMoigne
!!
!!      Modified by B. Decharme   (03/2009): Consistency with Arpege permanent
!!                                           snow/ice treatment
!!      Modified by A. Boone      (04/2010): Implicit coupling with atmosphere permitted.
!!
!!      Modified by B. Decharme   (04/2010): check suspicious low temperature for ES and CROCUS
!!      Modified by B. Decharme   (08/2013): Qsat as argument (needed for coupling with atm)
!!      Modified by A. Boone      (10/2014): MEB: pass in fluxes when using MEB
!!      Modified by M. Lafaysse   (08/2015): MEB-Crocus coupling
!!      Modified by B. Decharme   (03/2016): No snowdrift under forest
!!      Modified by P. Hagenmuller(09/2017): Mepra outputs
!!      Modified by A. Druel      (02/2019): Streamlines the code and adapt it to be compatible with new irrigation
!!      Modified by B. Decharme   (07/2019): add many diag for water and energy balance computation 
!!      Modified by B. Decharme   (07/2023): CHORT='CM6' -> Tuning for CNRM-CM/ESM to reduce the too earlier 
!!                                                          spring time snowmelt due to the non representation 
!!                                                          of vegetation-snow interaction. Obsolete with MEB. 
!!      Modified by L. Viallon G. (03/2024): Add SNOWHIST diagnostic for Crocus
!!
!-------------------------------------------------------------------------------
!
USE MODD_ISBA_OPTIONS_n,   ONLY : ISBA_OPTIONS_t
USE MODD_SFX_GRID_n,       ONLY : GRID_t
USE MODD_ISBA_n,           ONLY : ISBA_PE_t, ISBA_P_t
USE MODD_DIAG_n,           ONLY : DIAG_t
USE MODD_DIAG_EVAP_ISBA_n, ONLY : DIAG_EVAP_ISBA_t
USE MODD_DIAG_MISC_ISBA_n, ONLY : DIAG_MISC_ISBA_t
USE MODD_CSTS,             ONLY : XTT, XPI, XDAY, XLMTT, XLSTT, XRHOLW, XCI
USE MODD_SNOW_PAR,         ONLY : XRHOSMAX_ES, XSNOWDMIN, XRHOSMIN_ES, XEMISSN, XZ0SN
USE MODD_PREP_SNOW,        ONLY : NIMPUR, LSNOW_FRAC_TOT
USE MODD_SURF_PAR,         ONLY : XUNDEF
USE MODD_TYPE_DATE_SURF,   ONLY : DATE_TIME
USE MODD_DATA_COVER_PAR,   ONLY : NVT_SNOW, NVEGTYPE,             &
                                  NVT_TEBD, NVT_TRBE, NVT_BONE,   &
                                  NVT_TRBD, NVT_TEBE, NVT_TENE,   &
                                  NVT_BOBD, NVT_BOND, NVT_SHRB
USE MODD_AGRI,             ONLY : NVEG_IRR
USE MODD_CONST_TARTES,     ONLY : XPWAVEIND_MODIS, NPNBANDS_MODIS
!
USE MODD_BLOWSNW_SURF
USE MODI_SNOW3L
USE MODI_SNOWCRO
USE MODI_SNOWCRO_DIAG
USE MODI_SNOW_SYTRON
USE MODI_SNOW_MAKING
USE SNOWPAPPUS_ENGINE
!
#ifdef SFX_OL
USE MODN_IO_OFFLINE,       ONLY : XTSTEP_OUTPUT
#endif
!
USE MODI_ABOR1_SFX
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
USE MODD_SLOPE_EFFECT, ONLY: NIX, NIY
!
#ifdef SFX_MPI
USE MODD_SURFEX_MPI, ONLY: NSIZE_TASK, NRANK
! Warning NSIZE_TASK is an array with an index starting from 0 not 1 !!
#endif
!
IMPLICIT NONE
!
!
!
!*      0.1    declarations of arguments
!
TYPE(ISBA_OPTIONS_t),   INTENT(INOUT) :: IO
TYPE(GRID_t),           INTENT(INOUT) :: G
TYPE(ISBA_P_t),         INTENT(INOUT) :: PK
TYPE(ISBA_PE_t),        INTENT(INOUT) :: PEK
TYPE(DIAG_t),           INTENT(INOUT) :: DK
TYPE(DIAG_EVAP_ISBA_t), INTENT(INOUT) :: DEK
TYPE(DIAG_MISC_ISBA_t), INTENT(INOUT) :: DMK
!
LOGICAL, INTENT(IN)                 :: OMEB       ! True = coupled to MEB. This means surface fluxes ae IMPOSED
!                                                 ! as an upper boundary condition to the explicit snow schemes.
!                                                 ! If = False, then energy
!                                                 ! budget and fluxes are computed herein.
!
 CHARACTER(LEN=*),     INTENT(IN)    :: HIMPLICIT_WIND   ! wind implicitation option
!                                                       ! 'OLD' = direct
!                                                       ! 'NEW' = Taylor serie, order 1
!
TYPE(DATE_TIME), INTENT(IN)         :: TPTIME     ! current date and time
REAL, INTENT(IN)                    :: PTSTEP
!                                      PTSTEP    = time step of the integration
REAL,                 INTENT(IN)    :: PTIMEC    ! cumulated time since beginning
!
REAL, DIMENSION(:,:), INTENT(IN)    :: PVEGTYPE ! fraction of each vegetation
!
!
REAL, DIMENSION(:,:), INTENT(INOUT) :: PTG
!                                      PTG       = Soil temperature profile (K)
!
REAL, DIMENSION(:,:), INTENT(IN)    :: PSOILHCAPZ, PD_G, PDZG
REAL, DIMENSION(:),   INTENT(IN)    :: PCT, PSOILCONDZ
!                                      PD_G      = Depth to bottom of each soil layer (m)
!                                      PDZG      = Soil layer thicknesses (m)
!                                      DMK%XCG       = area-averaged soil heat capacity [(K m2)/J]
!                                      PCT       = area-averaged surface heat capacity [(K m2)/J]
!                                      PSOILCONDZ= soil thermal conductivity (W m-1 K-1)
!                                      PSOILHCAPZ= soil heat capacity (J m-3 K-1)
!
REAL, DIMENSION(:), INTENT(IN)      :: PPS, PTAR, PTAC, PSW_RAD, PQA,                       &
                                       PVMOD, PVDIR, PLW_RAD, PSR, PRR
!                                      PSW_RAD = incoming solar radiation (W/m2)
!                                      PLW_RAD = atmospheric infrared radiation (W/m2)
!                                      PRR     = rain rate [kg/(m2 s)]
!                                      PSR     = total snow rate (SWE) [kg/(m2 s)] (including PUNLOAD)
!                                      PUNLOAD = snow unloading rate (SWE) [kg/(m2 s)]
!                                      PTAC    = atmospheric temperature at level za (K)
!                                                NOTE, when MEB used, it corresponds to the canopy air T
!                                      PTAR    = reference air T above the surface/canopy (K)
!                                      PVMOD   = modulus of the wind parallel to the orography (m/s)
!                                      PVDIR   = wind direction (rad)
!                                      PPS     = surface pressure
!                                      PQA     = atmospheric specific humidity
!                                                at level za
REAL, DIMENSION(:,:), INTENT(IN)    :: P_DIR_SW, P_SCA_SW
!                                      P_DIR_SW, P_SCA_SW = direct and diffuse spectral solar irradiance (W/m2/um)
REAL, DIMENSION(:,:), INTENT(IN)    :: PIMPWET,PIMPDRY !(g m-2 s-1)
!
REAL, DIMENSION(:), INTENT(IN)      :: PZREF, PUREF, PEXNS, PEXNA, PDIRCOSZW, PSLOPEDIR, PRHOA, PALB
!                                      PZREF     = reference height of the first
!                                                  atmospheric level
!                                      PUREF     = reference height of the wind
!                                      PRHOA     = air density
!                                      PEXNS     = Exner function at surface
!                                      PEXNA     = Exner function at lowest atmos level
!                                      PDIRCOSZW = Cosinus of the angle between the
!                                                  normal to the surface and the vertical
!                                      PSLOPEDIR = Slope direction
!                                      PALB      = soil/vegetation albedo

!
REAL, DIMENSION(:), INTENT(IN)      :: PPEW_A_COEF, PPEW_B_COEF,                   &
                                       PPET_A_COEF, PPEQ_A_COEF, PPET_B_COEF,      &
                                       PPEQ_B_COEF
!                                      PPEW_A_COEF = wind coefficient
!                                      PPEW_B_COEF = wind coefficient
!                                      PPET_A_COEF = A-air temperature coefficient
!                                      PPET_B_COEF = B-air temperature coefficient
!                                      PPEQ_A_COEF = A-air specific humidity coefficient
!                                      PPEQ_B_COEF = B-air specific humidity coefficient
!
INTEGER , DIMENSION(:), INTENT(IN)  ::  KTAB_SYT    ! Array of index defining
                                                     ! opposite points for Sytron
REAL, DIMENSION(:), INTENT(INOUT)   :: PLES3L, PLEL3L, PEVAP, PGRNDFLUX
!                                      PLEL3L        = evaporation heat flux from snow (W/m2)
!                                      PLES3L        = sublimation (W/m2)
!                                      PEVAP         = total evaporative flux from snow (kg/m2/s)
!                                      PGRNDFLUX     = soil/snow interface heat flux (W/m2)
!
REAL, DIMENSION(:), INTENT(INOUT)   :: PRI
!                                      PRI        = Richardson number (-)
!
REAL, DIMENSION(:), INTENT(OUT)     :: PTHRUFAL, PFLSN_COR, PEVAPCOR
!                                      PTHRUFAL  = rate that liquid water leaves snow pack:
!                                                  paritioned into soil infiltration/runoff
!                                                  by ISBA [kg/(m2 s)]
!                                      PFLSN_COR = soil/snow correction heat flux (W/m2) 
!                                      PEVAPCOR  = evaporation/sublimation correction term:
!                                                  extract any evaporation exceeding the
!                                                  actual snow cover (as snow vanishes)
!                                                  and apply it as a surface soil water
!                                                  sink. [kg/(m2 s)]
!
REAL, DIMENSION(:), INTENT(OUT)     :: PSNOWSFCH
!
REAL, DIMENSION(:), INTENT(OUT)     :: PQS
!                                      PQS = surface humidity (kg/kg)
!
! ajout_EB pour prendre en compte angle zenithal du soleil dans LRAD
! puis plus tard dans LALB
REAL, DIMENSION(:), INTENT(IN)      :: PZENITH   ! solar zenith angle
REAL, DIMENSION(:), INTENT(IN)      :: PAZIM     ! azimuthal angle      (radian from North, clockwise)
!
INTEGER,DIMENSION(:), INTENT(IN)    :: NPAR_VEG_IRR_USE ! vegtype with irrigation
!
REAL, DIMENSION(:,:), INTENT(INOUT) :: PBLOWSNW_FLUX
!                                      PBLOWSNW_FLUX  = Blowing snow particles flux:
!                                       1: Number (#/m2/s) 2: Mass (kg/m2/s)
!                                       IN : contains sedimentation flux
!                                       OUT : contains emitted turbulent flux towards the atmosphere
REAL, DIMENSION(:,:), INTENT(IN)    :: PBLOWSNW_CONC
!                                      PBLOWSNW_CONC = Blowing snow particles concentration:
!                                       1: Number (#/m3) 2: Mass (kg/m3)
!
!*      0.2    declarations of local variables
!
LOGICAL, DIMENSION(SIZE(PPS))      ::  LCOLD ! flag for excessively cold and thin layers
!
INTEGER                             :: JWRK, JJ,JIMP ! Loop control
!
INTEGER                             :: JVEG, JK ! loop on vegtypes
!
INTEGER                             :: INLVLS   ! maximum number of snow layers
INTEGER                             :: INLVLG   ! number of ground layers
INTEGER                             :: IBLOWSNW ! number of blowing snow variables
INTEGER                             :: ILOCNIY  ! Processor local NIY value (number of lines for each proc)
!
REAL, DIMENSION(SIZE(PTG,1),SIZE(PTG,2)) :: ZTG0 ! Initial soil temperature profile
!
REAL, DIMENSION(SIZE(PTAR))         :: ZRRSNOW, ZSOILCOND, ZSNOW, ZSNOWFALL,  &
                                       ZSNOWABLAT_DELTA, ZSNOWSWE_1D, ZSNOWD, &
                                       ZSNOWH, ZSNOWH1, ZGRNDFLUXN, ZPSN,     &
                                       ZSOILCOR, ZSNOWSWE_OUT, ZTHRUFAL,      &
                                       ZSNOW_MASS_BUDGET, ZWGHT, ZWORK, ZC2,  &
                                       ZSNOW_ENERGY_BUDGET
!                                      ZSOILCOND    = soil thermal conductivity [W/(m K)]
!                                      ZRRSNOW      = rain rate over snow [kg/(m2 s)]
!                                      ZSNOW        = snow depth (m)
!                                      ZSNOWFALL    = minimum equivalent snow depth
!                                                     for snow falling during the
!                                                     current time step (m)
!                                      ZSNOWABLAT_DELTA = FLAG =1 if snow ablates completely
!                                                     during current time step, else=0
!                                      ZSNOWSWE_1D  = TOTAL snowpack SWE (kg m-2)
!                                      ZSNOWD       = snow depth
!                                      ZSNOWH       = snow total heat content (J m-2)
!                                      ZSNOWH1      = snow surface layer heat content (J m-2)
!                                      ZGRNDFLUXN   = corrected snow-ground flux (if snow fully ablated during timestep)
!                                      ZPSN         = snow fraction working array
!                                      ZSOILCOR = for vanishingy thin snow cover,
!                                                 allow any excess evaporation
!                                                 to be extracted from the soil
!                                                 to maintain an accurate water
!                                                 balance [kg/(m2 s)]
!                                      ZSNOW_MASS_BUDGET = snow water equivalent budget (kg/m2/s)
!                                      ZWGHT        = MEB surface layer weight for distributing energy
!                                                     between litter and ground layers for the case
!                                                     of total ablation during a timestep (-).
!                                      ZWORK        = local working variable (*)
!                                      ZC2          = sub-surface heat capacity [(K m2)/J]
!                                      ZSNOW_ENERGY_BUDGET = snow energy budget (W/m2)

REAL, DIMENSION(SIZE(PVMOD),4)      :: ZBLOWSNW   ! Properties of deposited blowing snow
!                                                   1 : Deposition flux (kg/m2/s)
!                                                   2 : Density of deposited snow (kg/m3)
!                                                   3 : SGRA1 of deposited snow diam opt
!                                                   4 : SGRA2 of deposited snow Sphericité
REAL, DIMENSION(SIZE(PVMOD),3)      :: ZQ_OUT
!
REAL, DIMENSION(SIZE(PVMOD),11)     :: ZPAPPUS_DEBUG   ! 1 ZFRIC 2 ZFRIC_T 3 PZ0 4 ZHSALT_SUSP 5 ZHSALT1 6 ZHSALT2 7 ZVFALL
!
REAL, DIMENSION(SIZE(PTAR))         :: ZBLOWSNW_ACC, ZBLOWSNW_DEPFLUX
!                                      ZBLOWSNW_ACC  = minimum equivalent snow depth
!                                                      for deposition of blown snow particles
!                                                      during the current time step (m)
!                                      ZBLOWSNW_DEPFLUX = deposition flux of blowing snow (kg/m2/s)
REAL, DIMENSION(SIZE(DK%XZ0))       :: ZZ0
REAL                                :: ZSNOWDMIN_S, ZCDT_MIN, ZSNOWSWE_MIN ! variables for checking cold/thin snowpack
REAL, DIMENSION(SIZE(PPS))          :: ZSNOWHEAT_S, ZSNOWSWE_S, ZSNOWDZ_S, ZSNOWAGE_S, ZSNOWAGE_X, ZSNOWRHO_X, &
                                       ZSNOWDZ_X, ZSNOWTEMP_X, ZMASS_COR, ZHEAT_COR
!                                      These are bulk snowpack variables used to correct for cold-thin snowpacks
!                                      ZSNOWHEAT_S = total snow heat content (J m-2)
!                                      ZSNOWSWE_S  = total snow mass (kg m-2)
!                                      ZSNOWDZ_S   = total snow depth (m)
!                                      ZSNOWAGE_S  = average snowpack age (s)
!                                      ZSNOWAGE_X  = working snowpack age (s)
!                                      ZSNOWRHO_X  = working snow density (kg m-3)
!                                      ZSNOWDZ_X   = working snow depth (m)
!                                      ZSNOWTEMP_X = working snow T (K)
!                                      ZMASS_COR   = possible snow mass correction (kg m-2)
!                                      ZHEAT_COR   = possible snow heat correction (J m-2)
!
!*      0.3    declarations of packed  variables
!
INTEGER                            :: ISIZE_SNOW ! number of points where computations are done
INTEGER, DIMENSION(SIZE(PTAR))     :: NMASK      ! indices correspondance between arrays
!
LOGICAL, DIMENSION(SIZE(PTAR))     :: GREMOVE_SNOW
!
REAL, DIMENSION(SIZE(PTAR))        :: ZSWNET_N, ZSWNET_NS, ZLWNET_N
!
REAL, DIMENSION(SIZE(PPS))         :: ZSNOWMAK
!
!*      0.4    declarations of local parameters
!
INTEGER, PARAMETER                 ::  JNSOIL_COR     = 3  ! spread any vanishing snow over
                                                           ! DIFF soil layers. Must be >= 1
                                                           ! currently just sfc for FR soil
REAL,    PARAMETER                 :: ZCHECK_TEMP = 150.0  ! K
!                                      Limit to check suspicious low temperature (K)
REAL,    PARAMETER                 :: ZCHECK_TEMP_DT = 5. ! safely above ZCHECK_TEMP (K)

LOGICAL                            :: GREG_GRID ! Boolean to set if grid is regular
REAL                               :: ZCELLSIZE ! Value of cell size for regular grid
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
! - - ---------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('SNOW3L_ISBA',0,ZHOOK_HANDLE)
!
!*       0.     Initialize variables:
!               ---------------------
!
IF (SIZE(DMK%XSNOWDEND)>0) THEN
  DMK%XSNOWDEND (:,:) = XUNDEF
  DMK%XSNOWSPHER(:,:) = XUNDEF
  DMK%XSNOWSIZE (:,:) = XUNDEF
  DMK%XSNOWSSA  (:,:) = XUNDEF
  DMK%XSNOWHIST (:,:) = XUNDEF
  DMK%XSNOWRAM  (:,:) = XUNDEF
  DMK%XSNOWSHEAR(:,:) = XUNDEF
  DMK%XSNOWTYPEMEPRA(:,:) = XUNDEF
  DMK%XACC_RAT      (:,:) = XUNDEF
  DMK%XNAT_RAT      (:,:) = XUNDEF
ENDIF
!
IF (ASSOCIATED(G%XMESH_SIZE) .AND. (NIX .NE. 0)) THEN ! check if regular rectangular grid
  GREG_GRID = .TRUE.
  ZCELLSIZE = G%XMESH_SIZE(1)
ELSE
  GREG_GRID=.FALSE. ! not grid initialisation
  ZCELLSIZE = 0
END IF
!
IF (PEK%TSNOW%SCHEME=='CRO') THEN
  DMK%XIMPUR_CONC (:,:,:) = XUNDEF
  DMK%XSPEC_ALB(:,:) = XUNDEF
  DMK%XDIFF_RATIO(:,:) = XUNDEF
  DMK%XSPEC_TOT(:,:) = XUNDEF
  IF (SIZE(DMK%XSPECMOD)>1) THEN
    DMK%XSPECMOD(:,:) = XUNDEF
  ENDIF
ENDIF
!
DMK%XSNOWHMASS(:)  = 0.0
DMK%XSRSFC(:)      = PSR(:)         ! these are snow and rain rates passed to ISBA,
DMK%XRRSFC(:)      = PRR(:)         ! so initialize here if SNOW3L not used:
!
PFLSN_COR(:)   = 0.0
PTHRUFAL(:)    = 0.0
PEVAPCOR(:)    = 0.0
PQS(:)         = XUNDEF
!
DEK%XSNDRIFT  (:) = 0.0
DEK%XMELTSTOT (:) = 0.0
DEK%XSNREFREEZ(:) = 0.0
!
IF (PEK%TSNOW%SCHEME=='CRO' .AND. IO%LSNOWSYTRON) THEN
   DEK%XSYTMASS(:) = 0.0
ENDIF
!
ZSNOW(:)       = 0.0
ZSNOWD(:)      = 0.0
ZGRNDFLUXN(:)  = 0.0
ZSNOWH(:)      = 0.0
ZSNOWH1(:)     = 0.0
ZSNOWSWE_1D(:) = 0.0
ZSNOWSWE_OUT(:)= 0.0
ZSOILCOND(:)   = 0.0
ZRRSNOW(:)     = 0.0
ZSNOWFALL(:)   = 0.0
ZSNOWABLAT_DELTA(:) = 0.0
ZWGHT(:)       = 0.0
ZWORK(:)       = 0.0
ZC2(:)         = PCT(:)
!
DMK%XSNOWLIQ(:,:) = 0.0
DMK%XSNOWDZ (:,:) = 0.0
ZTG0        (:,:) = PTG(:,:)
!
ZBLOWSNW(:,:)  = 0.0
ZBLOWSNW_ACC(:)  = 0.0
ZBLOWSNW_DEPFLUX(:) = 0.0
ZSNOWMAK(:) = 0.0
!
INLVLS   = SIZE(PEK%TSNOW%WSNOW(:,:),2)
INLVLG   = MIN(SIZE(PD_G(:,:),2),SIZE(PTG(:,:),2))
IBLOWSNW = SIZE(ZBLOWSNW(:,:),2)
!
!
IF(.NOT.OMEB)THEN
  !
  ! If MEB activated, these values are input, else initialize here:
  !
  PGRNDFLUX(:)      = 0.0
  PLES3L(:)         = 0.0
  PLEL3L(:)         = 0.0
  PEVAP(:)          = 0.0
  PRI(:)            = XUNDEF
  PEK%TSNOW%EMIS(:) = XEMISSN
  DMK%XRNSNOW(:)    = 0.0
  DMK%XHSNOW(:)     = 0.0
  DMK%XGFLUXSNOW(:) = 0.0
  DMK%XHPSNOW(:)    = 0.0
  DMK%XUSTARSNOW(:) = 0.0
  DMK%XCDSNOW(:)    = 0.0
  DMK%XCHSNOW(:)    = 0.0
  !
ENDIF
!
! Use ISBA-SNOW3L or NOT: NOTE that if explicit soil diffusion method in use,
! then *must* use explicit snow model:
!
IF (PEK%TSNOW%SCHEME=='3-L' .OR. IO%CISBA == 'DIF' .OR. PEK%TSNOW%SCHEME == 'CRO') THEN
!
  IF(.NOT.OMEB)THEN
    !
    ! If MEB activated, these values are input, else initialize here:
    !
    ZSWNET_N (:) = 0.0
    ZSWNET_NS(:) = 0.0
    ZLWNET_N (:) = 0.0
  ELSE
    ZSWNET_N (:) = DEK%XSWNET_N (:)
    ZSWNET_NS(:) = DEK%XSWNET_NS(:)
    ZLWNET_N (:) = DEK%XLWNET_N (:)
  END IF
  !
  ! - Snow and rain falling onto the 3-L grid space:
  !
  DMK%XSRSFC(:) = 0.0
  !
  DO JJ=1,SIZE(PSR)
    ZRRSNOW(JJ)        = PEK%XPSN(JJ)*PRR(JJ)
    DMK%XRRSFC(JJ)     = PRR(JJ) - ZRRSNOW(JJ)
    ZSNOWFALL(JJ)      = PSR(JJ)*PTSTEP/XRHOSMAX_ES  ! maximum possible snowfall depth (m) does not take account redistribution scheme (check ZBLOWSNW_ACC)
  ENDDO

  IF (IO%CISBA == 'DIF') THEN
    ZSOILCOND(:) = PSOILCONDZ(:)
  ELSE
    !
    ! - Soil thermal conductivity
    !   is implicit in Force-Restore soil method, so it
    !   must be backed-out of surface thermal coefficients
    !   (Etchevers and Martin 1997):
    !
    ZSOILCOND(:) = 4.*XPI/( DMK%XCG(:)*DMK%XCG(:)*XDAY/(PD_G(:,1)*PCT(:)) )
    !
  ENDIF
  !
  ! ===============================================================
  ! HERE STARTS THE SNOW REDISTRIBUTION: 3 Schemes are possible
  ! - Sytron (LSNOWSYTRON in namelist)
  ! - Meso-NH coupling (PBLOWSNW_FLUX from where ?)  <- TO BE CHECKED
  ! - SnowPappus (LSNOWPAPPUS in namelist)
  !
  ! General idea of this part of the code:
  ! - snow redistribution change the value of ZBLOWSNW_ACC
  ! - snow machine made can change ZSNOWFALL
  ! - packing: call the model where there is snow (due to snowdz or snowfall or zblowsnw_acc)
  !
  !
  ! ===============================================================
  !        Snow redistribution scheme Sytron
  !
  IF (PEK%TSNOW%SCHEME=='CRO' .AND. IO%LSNOWSYTRON) THEN

    CALL SNOW_SYTRON(PTSTEP,PPS,PTAC,PQA,PVMOD,PVDIR,PSLOPEDIR,PDIRCOSZW,  &
                     PEK%TSNOW%HEAT,PEK%TSNOW%WSNOW,PEK%TSNOW%RHO,         &
                     PEK%TSNOW%DIAMOPT,PEK%TSNOW%SPHERI,PEK%TSNOW%HIST,    &
                     PEK%TSNOW%AGE,KTAB_SYT, ZBLOWSNW,DEK%XSYTMASS,IO%CSNOWMOB)

    !
    ! Calculate maximum snow depth (m) of deposited blown snow particles
    !
    WHERE(ZBLOWSNW(:,1)> 0.)
      ZBLOWSNW_ACC(:)=ZBLOWSNW(:,1)*PTSTEP/ZBLOWSNW(:,2)
    END WHERE
  ENDIF
  !
  ! ===============================================================
  !        Snow redistribution when coupled to Meso-NH
  !
  !
  ! TO BE CHECKED: two times "Calculate maximum deposited..." in that case
  !
  ! Calculate maximum deposited snow depth (m) of blown snow particles
  !
  IF(SIZE(PBLOWSNW_FLUX,2) /= 0) THEN
    DO JJ=1,SIZE(PSR)
      ZBLOWSNW_ACC(JJ)=(PBLOWSNW_FLUX(JJ,2)+PBLOWSNW_FLUX(JJ,3))*PTSTEP/XRHO_DEP
      IF (PBLOWSNW_FLUX(JJ,2)+PBLOWSNW_FLUX(JJ,3)>0.) THEN
        ZBLOWSNW_DEPFLUX(JJ) = (PBLOWSNW_FLUX(JJ,2)+PBLOWSNW_FLUX(JJ,3))
      ENDIF
    ENDDO
  ENDIF
  !
  IF (PEK%TSNOW%SCHEME =='CRO' .AND. SIZE(PBLOWSNW_FLUX,2) /= 0) THEN
    CALL SNOWPACK_EVOL(IO%CSNOWRES,PBLOWSNW_FLUX,PEK%TSNOW%HEAT,          &
                       PEK%TSNOW%WSNOW,PEK%TSNOW%RHO,                     &
                       PEK%TSNOW%DIAMOPT,PEK%TSNOW%SPHERI,PEK%TSNOW%HIST, &
                       PEK%TSNOW%AGE, PTSTEP,PRHOA,PTAC,                  &
                       PBLOWSNW_CONC,PVMOD, PQA,PPS,                      &
                       PUREF,PEXNS,PDIRCOSZW,                             &
                       PZREF,DK%XZ0EFF,DK%XZ0H,ZBLOWSNW,                  &
                       PTG(:,1))
    !
    ! Calculate maximum snow depth (m) of deposited blown snow particles
    !
    WHERE(ZBLOWSNW(:,1)> 0.) !traitement similaire a sytron masque sur la deposition
      ZBLOWSNW_ACC(:)=ZBLOWSNW(:,1)*PTSTEP/ZBLOWSNW(:,2)
    END WHERE
  ENDIF
  !
  ! ===============================================================
  !        Snow redistribution scheme Pappus
  ! ===============================================================
  !
  IF (PEK%TSNOW%SCHEME=='CRO' .AND. IO%LSNOWPAPPUS) THEN
    IF (IO%LPAPPUDEBUG) THEN
      WRITE(*,*)"ffffffffffffffffffffffffffffffffffffffffff"
      WRITE(*,*) "date", TPTIME%TDATE%YEAR, TPTIME%TDATE%MONTH, TPTIME%TDATE%DAY
      WRITE(*,*) "heure", TPTIME%TIME
      WRITE(*,*)"GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG"
    ENDIF
    ! WRITE(*,*) G%XMESH_SIZE ! not allocated
    ! WRITE(*,*)"ffffffffffffffffffffffffffffffffffffffffff"
    ! WRITE(*,*)"wsnow :", PEK%TSNOW%WSNOW(1,:)
    ! WRITE(*,*)"HEAT :", PEK%TSNOW%HEAT(1,:)
    ! WRITE(*,*)"RHO :", PEK%TSNOW%RHO(1,:)
    ! WRITE(*,*)"DEP_SUP :", PEK%TSNOW%DEP_SUP(1)
    ! WRITE(*,*)"DEP_TOT :", PEK%TSNOW%DEP_TOT(1)
    ! WRITE(*,*) (associated(G%XMESH_SIZE))
    ! WRITE(*,*) "nlayer", PEK%TSNOW%NLAYER
    ! WRITE(*,*) "depth",PEK%TSNOW%DEP_SUP(:)
    ! WRITE(*,*) "depth tot", PEK%TSNOW%DEP_TOT(1)
    !=============================================================================================================================>
    ZZ0(:)=XZ0SN ! Pour resoudre le bug des z0 initialisé dans snow3L A LA PLACE DE DK%XZ0
    ! initialisation diag cumul
    DEK%XQDEP_TOT(:)  = 0.0
    DEK%XQ_OUT_SUBL(:)= 0.0
    DEK%XQT_TOT(:)    = 0.0
    DEK%XSNOWDEBTC(:) = 0.0
    DEK%XBLOWSNWFLUX_1M(:)=0.0
    DEK%XBLOWSNWFLUXINT(:)=0.0
    !=============================================================================================================================>
    CALL SNOWPAPPUS(IO%CSNOWPAPPUSERODEPO, IO%CSALTPAPPUS, IO%CPAPPUSSUBLI, IO%CLIMVFALL, IO%LPAPPULIMTFLUX,     &
                    IO%LPAPPUDEBUG, IO%CSNOWMOB, PTSTEP, ZCELLSIZE, ZZ0, PUREF, PVMOD, PVDIR, PRHOA, PPS,  &
                    PQA,PTAC, PDIRCOSZW, PEK, DEK%XSNOWDEBTC, PTHRUFAL, ZQ_OUT, DEK%XQT_TOT, DEK%XBLOWSNWFLUX_1M,&
                    DEK%XBLOWSNWFLUXINT, ZBLOWSNW, ZPAPPUS_DEBUG, GREG_GRID)
    ! XMESH_SIZE(1) SUPPOSE QUE XMESH_SIZE CONSTANT SUR TOUT (:)
    ! NE PAS OUBLIER DE METTRE A JOURS TOUTES LES VALEUR LOCALE QUI ON ETE INITIALISÉ AVEC PEK AVANT SNOWPAPPUS
    ! Initisalisation de diags :
    DMK%XQ_OUT_SALT(:)=ZQ_OUT(:,1)
    DMK%XQ_OUT_SUSP(:)=ZQ_OUT(:,2)
    DEK%XQ_OUT_SUBL(:)=ZQ_OUT(:,3)
    DMK%XVFRIC_PAPPUS=ZPAPPUS_DEBUG(:,1)
    DMK%XVFRIC_T_PAPPUS= ZPAPPUS_DEBUG(:,2)
    DMK%XPZ0_PAPPUS= ZPAPPUS_DEBUG(:,3)
    DMK%XVFALL_PAPPUS= ZPAPPUS_DEBUG(:,7)
    DEK%XQDEP_TOT(:)=ZPAPPUS_DEBUG(:,4)
    !
    WHERE(ZBLOWSNW(:,1)> 0.)
      ZBLOWSNW_ACC(:)=ZBLOWSNW(:,1)*PTSTEP/ZBLOWSNW(:,2)
      ZBLOWSNW_DEPFLUX(:) = ZBLOWSNW(:,1)
    END WHERE
  ENDIF
  !
  ! ===============================================================
  ! END OF SNOW REDISTRIBUTION
  ! ===============================================================
  ! AFTER Snow redistribution : Calculate preliminary snow depth (m)
  ! TAKE CARE this MUST be done AFTER snow transport modules
  !
  ZSNOW(:)      = 0.
  ZSNOWH(:)     = 0.
  ZSNOWSWE_1D(:)= 0.
  ZSNOWH1(:)    = PEK%TSNOW%HEAT(:,1)*PEK%TSNOW%WSNOW(:,1)/PEK%TSNOW%RHO(:,1) ! sfc layer only
  !
  DO JWRK=1,SIZE(PEK%TSNOW%WSNOW(:,:),2)
    DO JJ=1,SIZE(PEK%TSNOW%WSNOW(:,:),1)
      ZSNOWSWE_1D(JJ) = ZSNOWSWE_1D(JJ) + PEK%TSNOW%WSNOW(JJ,JWRK)
      ZSNOW(JJ)       = ZSNOW(JJ)       + PEK%TSNOW%WSNOW(JJ,JWRK)/PEK%TSNOW%RHO(JJ,JWRK)
      ZSNOWH(JJ)      = ZSNOWH(JJ)      + PEK%TSNOW%HEAT (JJ,JWRK)*PEK%TSNOW%WSNOW(JJ,JWRK)/PEK%TSNOW%RHO(JJ,JWRK)
    END DO
  ENDDO
  !
  ! ===============================================================
  ! MPI DEBUG OUTPUT
  ! ===============================================================
  ! OUTPUT AVALAIBLE ONLY FOR CROCUS AND IGN RECTANGULAR GRID
  ! OUTPUT CHANGED IF MPI COMPILATION
  !
  IF (PEK%TSNOW%SCHEME=='CRO') THEN
    IF (GREG_GRID) THEN ! if rectangular grid
#ifdef SFX_MPI
      ILOCNIY = NSIZE_TASK(NRANK)/NIX
      DO JJ = 1, SIZE(PDIRCOSZW)
        DMK%XJJ(JJ)=JJ
        DMK%XNRANK(JJ)=NRANK
        DMK%XILOCNIY(JJ)=ILOCNIY
        DMK%XSIZE_TASK(JJ)=NSIZE_TASK(NRANK)
      END DO
#else
      ILOCNIY = NIY
      DO JJ = 1, SIZE(PDIRCOSZW)
        DMK%XJJ(JJ)=JJ
        DMK%XNRANK(JJ)=XUNDEF
        DMK%XILOCNIY(JJ)=ILOCNIY
        DMK%XSIZE_TASK(JJ)=XUNDEF
      END DO
#endif
    END IF
  END IF
  !
  ! ===============================================================
  ! Snow machine made can change the value of snowfall
  ! ===============================================================
  !
  IF (PEK%TSNOW%SCHEME=='CRO' .AND. IO%LSNOWMAK_BOOL) THEN
    CALL SNOW_MAKING(IO, TPTIME, PTSTEP, PEK, PSR, PTAC, PVMOD, PPS, PQA, PRHOA, ZSNOW, ZSNOWFALL, ZSNOWMAK, DMK%XWBT)
  ENDIF
  !
  ! ===============================================================
  ! PACKING: Only call snow model when:
  !          - there is snow on the surface exceeding a minimum threshold
  !          - OR if the equivalent snow depth falling during the current time
  !            step exceeds this limit
  !          - OR if the blowing snow from snow redistribution exceeds this limit
  !
  ! counts the number of points where the computations will be made
  !
  !
  ISIZE_SNOW = 0
  NMASK(:) = 0
  !
  IF (PEK%TSNOW%SCHEME =='CRO' .AND. SIZE(PBLOWSNW_FLUX,2) /= 0) THEN
    DO JJ=1,SIZE(ZSNOW)
      IF (ZSNOW(JJ) >= XSNOWDMIN .OR. ZSNOWFALL(JJ) >= XSNOWDMIN .OR. ZBLOWSNW_ACC(JJ) >= XSNOWDMIN) THEN
        ISIZE_SNOW = ISIZE_SNOW + 1
        NMASK(ISIZE_SNOW) = JJ
      ENDIF
    ENDDO
  ELSE
    DO JJ=1,SIZE(ZSNOW)
      IF (ZSNOW(JJ) >= XSNOWDMIN .OR. ZSNOWFALL(JJ) >= XSNOWDMIN) THEN
        ISIZE_SNOW = ISIZE_SNOW + 1
        NMASK(ISIZE_SNOW) = JJ
      ENDIF
    ENDDO
  ENDIF
  !
  IF (ISIZE_SNOW>0) CALL CALL_MODEL(ISIZE_SNOW,INLVLS,INLVLG,IBLOWSNW,NMASK)
  !
  ! ===============================================================
  !
  ! Remove trace amounts of snow and reinitialize snow prognostic variables
  ! if snow cover is ablated.
  ! If MEB used, soil T already computed, therefore correct heating/cooling
  ! effect of updated snow-soil flux
  !
  ZSNOWD(:) = 0.
  ZSNOWSWE_OUT(:) = 0.
  DO JWRK=1,SIZE(PEK%TSNOW%WSNOW(:,:),2)
    DO JJ=1,SIZE(PEK%TSNOW%WSNOW(:,:),1)
      ZSNOWD      (JJ) = ZSNOWD      (JJ) + PEK%TSNOW%WSNOW(JJ,JWRK)/PEK%TSNOW%RHO(JJ,JWRK)
      ZSNOWSWE_OUT(JJ) = ZSNOWSWE_OUT(JJ) + PEK%TSNOW%WSNOW(JJ,JWRK)
    ENDDO
  END DO
  !
  GREMOVE_SNOW(:)=(ZSNOWD(:)<XSNOWDMIN*1.1) !eqv ZSNOWD< 1.1E-6
  !
  IF(OMEB)THEN
    ZPSN(:) = 1.0
    IF(IO%CISBA == 'DIF')THEN
      ZWGHT(:) = PSOILHCAPZ(:,2)*PDZG(:,2)/(PSOILHCAPZ(:,1)*PDZG(:,1) + PSOILHCAPZ(:,2)*PDZG(:,2))
      ZC2(:)   = 1.0/(PSOILHCAPZ(:,2)*PDZG(:,2))
    ELSE
      ZWGHT(:) = (PD_G(:,2)-PD_G(:,1))/PD_G(:,2)
    ENDIF
  ELSE
    !  To Conserve mass in ISBA without MEB,
    !  EVAP must be weignted by the snow fraction
    !  in the calulation of THRUFAL
    ZPSN(:) = PEK%XPSN(:)
  ENDIF
  !
  ZSNOWABLAT_DELTA(:) = 0.0
  ZTHRUFAL        (:) = PTHRUFAL(:)
  !
  WHERE(GREMOVE_SNOW(:))
    !
    ZSNOWSWE_OUT(:)     = 0.0
    PLES3L(:)           = MIN(PLES3L(:), XLSTT*(ZSNOWSWE_1D(:)/PTSTEP + PSR(:)+ZBLOWSNW_DEPFLUX(:)))
    PLEL3L(:)           = 0.0
    PEVAP(:)            = PLES3L(:)/PK%XLSTT(:)
    PTHRUFAL(:)         = MAX(0.0, ZSNOWSWE_1D(:)/PTSTEP + PSR(:) + ZBLOWSNW_DEPFLUX(:) - PEVAP(:)*ZPSN(:) + ZRRSNOW(:)) ! kg m-2 s-1 manque snowmaking
    ZTHRUFAL(:)         = MAX(0.0, ZSNOWSWE_1D(:)/PTSTEP + PSR(:) + ZBLOWSNW_DEPFLUX(:) - PEVAP(:)         + ZRRSNOW(:)) ! kg m-2 s-1 manque snowmaking
    DEK%XMELTSTOT(:)    = PTHRUFAL(:)
    DEK%XSNREFREEZ(:)   = 0.0
    !
    DMK%XSRSFC(:)       = 0.0
    DMK%XRRSFC(:)       = DMK%XRRSFC(:)
    !
    ZSNOWABLAT_DELTA(:) = 1.0
    !
    PEK%TSNOW%ALB(:)    = XUNDEF
    !
    PEVAPCOR(:)         = 0.0
    ZSOILCOR(:)         = 0.0
    !
    DMK%XGFLUXSNOW(:)   = DMK%XRNSNOW(:) - DMK%XHSNOW(:) - PLES3L(:) - PLEL3L(:)
    DMK%XSNOWHMASS(:)   = -(PSR(:)+ZBLOWSNW_DEPFLUX(:))*(XLMTT*PTSTEP) ! manque snowmaking pour le bilan
    !
    DEK%XRESTOREN(:)     = 0.0
    DEK%XDELHEATN(:)     = -ZSNOWH(:) /PTSTEP
    DEK%XDELHEATN_SFC(:) = -ZSNOWH1(:)/PTSTEP
    DEK%XDELPHASEN(:)    = 0.0
    DEK%XDELPHASEN_SFC(:)= 0.0 
    !
    PSNOWSFCH(:)        = DEK%XDELHEATN_SFC(:) - (ZSWNET_NS(:) + ZLWNET_N(:) - DMK%XHSNOW(:) &
                        - PLES3L(:) - PLEL3L(:)) + DEK%XRESTOREN(:) - DMK%XSNOWHMASS(:)/PTSTEP 
    !
    ZGRNDFLUXN(:)       = (ZSNOWH(:)+DMK%XSNOWHMASS(:))/PTSTEP + DMK%XGFLUXSNOW(:)
    ZWORK(:)            = ZPSN(:)*(ZGRNDFLUXN(:) - PGRNDFLUX(:) - PFLSN_COR(:))
    PTG(:,1)            = ZTG0(:,1) + PTSTEP*ZWORK(:)*(1.-ZWGHT(:))*PCT(:)
    PTG(:,2)            = ZTG0(:,2) + PTSTEP*ZWORK(:)*    ZWGHT(:) *ZC2(:)
    PGRNDFLUX(:)        = ZGRNDFLUXN(:)
    DEK%XDELHEATG_SFC(:)= DEK%XDELHEATG_SFC(:) + (PTG(:,1)-ZTG0(:,1))/(PTSTEP*PCT(:)) - ZWORK(:)*(1.-ZWGHT(:))
    DEK%XDELHEATG    (:)= DEK%XDELHEATG    (:) + (PTG(:,2)-ZTG0(:,2))/(PTSTEP*ZC2(:)) - ZWGHT(:) * ZWORK(:) + DEK%XDELHEATG_SFC(:)  
    !
  END WHERE
  !
  IF (PEK%TSNOW%SCHEME=='CRO') THEN 
    WHERE(GREMOVE_SNOW(:))
      PEK%TSNOW%DEP_SUP(:) = 0
      PEK%TSNOW%DEP_TOT(:) = 0
      PEK%TSNOW%DEP_HUM(:) = 0
      PEK%TSNOW%NAT_LEV(:) = 6
      PEK%TSNOW%AVA_TYP(:) = 6
      PEK%TSNOW%PRO_SUP_TYP(:) = 6
    END WHERE
  ENDIF
  !
  !
  DO JWRK=1,INLVLS
    DO JJ=1,SIZE(PEK%TSNOW%WSNOW(:,:),1)
      PEK%TSNOW%WSNOW(JJ,JWRK) = (1.0-ZSNOWABLAT_DELTA(JJ))*PEK%TSNOW%WSNOW(JJ,JWRK)
      PEK%TSNOW%HEAT (JJ,JWRK) = (1.0-ZSNOWABLAT_DELTA(JJ))*PEK%TSNOW%HEAT (JJ,JWRK)
      PEK%TSNOW%RHO  (JJ,JWRK) = (1.0-ZSNOWABLAT_DELTA(JJ))*PEK%TSNOW%RHO  (JJ,JWRK) + ZSNOWABLAT_DELTA(JJ)*XRHOSMIN_ES
      PEK%TSNOW%AGE  (JJ,JWRK) = (1.0-ZSNOWABLAT_DELTA(JJ))*PEK%TSNOW%AGE  (JJ,JWRK)
      DMK%XSNOWTEMP  (JJ,JWRK) = (1.0-ZSNOWABLAT_DELTA(JJ))*DMK%XSNOWTEMP  (JJ,JWRK) + ZSNOWABLAT_DELTA(JJ)*XTT
      DMK%XSNOWLIQ   (JJ,JWRK) = (1.0-ZSNOWABLAT_DELTA(JJ))*DMK%XSNOWLIQ   (JJ,JWRK)
      DMK%XSNOWDZ    (JJ,JWRK) = (1.0-ZSNOWABLAT_DELTA(JJ))*DMK%XSNOWDZ    (JJ,JWRK)
    ENDDO
  ENDDO
  DO JIMP=1,NIMPUR
    DO JWRK=1,INLVLS
      DO JJ=1,SIZE(PEK%TSNOW%WSNOW,1)
        PEK%TSNOW%IMPUR (JJ,JWRK,JIMP)=(1.0-ZSNOWABLAT_DELTA(JJ))*PEK%TSNOW%IMPUR(JJ,JWRK,JIMP) !F.T
      ENDDO
    ENDDO
  ENDDO
  !
  IF (PEK%TSNOW%SCHEME=='CRO') THEN
    DO JWRK=1,INLVLS
      DO JJ=1,SIZE(PEK%TSNOW%DIAMOPT(:,:),1)
        PEK%TSNOW%DIAMOPT(JJ,JWRK) = (1.0-ZSNOWABLAT_DELTA(JJ))*PEK%TSNOW%DIAMOPT(JJ,JWRK)
        PEK%TSNOW%SPHERI (JJ,JWRK) = (1.0-ZSNOWABLAT_DELTA(JJ))*PEK%TSNOW%SPHERI (JJ,JWRK)
        PEK%TSNOW%HIST   (JJ,JWRK) = (1.0-ZSNOWABLAT_DELTA(JJ))*PEK%TSNOW%HIST   (JJ,JWRK)
      ENDDO
    ENDDO
  ENDIF
  !
  !  ===============================================================
  !
  !  Compute snow mass budget and energy budget
  !
  ZSNOW_MASS_BUDGET(:) = (ZSNOWSWE_1D(:)-ZSNOWSWE_OUT(:))/PTSTEP + PSR     (:)+ZRRSNOW (:) &
                                                                 - PEVAP   (:)-ZTHRUFAL(:) &
                                                                 + PEVAPCOR(:)+ZSOILCOR(:)
  !
  ZSNOW_ENERGY_BUDGET(:) = DEK%XDELHEATN(:)-DMK%XSNOWHMASS(:)/PTSTEP-DMK%XGFLUXSNOW(:)+PGRNDFLUX(:)+PFLSN_COR(:)
  !
  !  ===============================================================
  !
  !  To Conserve mass in ISBA, the latent heat flux part of
  !  the EVAPCOR term must be weignted by the snow fraction
  !
  PEVAPCOR (:) = PEVAPCOR(:)*ZPSN(:) + ZSOILCOR(:)
  !
  ! ===============================================================
  !
! If approaching COLD limit, remix properties and regrid
! NOTE happens when snow is very thin, so to conserve energy highly accurately,
! put any excess energy (should be quite small if present) into ground.

! Identify excessive cooling:
  
  LCOLD(:) = .false.
  DO JWRK=1,INLVLS
     DO JJ=1,SIZE(PEK%TSNOW%WSNOW,1)
        IF (PEK%TSNOW%WSNOW(JJ,JWRK)>0.0) THEN
           IF (DMK%XSNOWTEMP(JJ,JWRK)<ZCHECK_TEMP)LCOLD(JJ) = .true.
        ENDIF
     ENDDO
  ENDDO

! sum properties for cold points:
  
  ZSNOWHEAT_S(:) = 0.
  ZSNOWSWE_S(:)  = 0.
  ZSNOWDZ_S(:)   = 0.
  ZSNOWAGE_S(:)  = 0.
  ZSNOWAGE_X(:)  = 0.
  ZSNOWDZ_X(:)   = 0.
  DO JWRK=1,INLVLS
     DO JJ=1,SIZE(PEK%TSNOW%WSNOW(:,:),1)
        IF(LCOLD(JJ))THEN
           ZSNOWAGE_X(JJ)  = MAX(ZSNOWAGE_X(JJ),PEK%TSNOW%AGE (JJ,JWRK)) ! max for each point
           ZSNOWDZ_X(JJ)   = PEK%TSNOW%WSNOW(JJ,JWRK)/ &
                             MAX(XRHOSMIN_ES,MIN(XRHOSMAX_ES,PEK%TSNOW%RHO  (JJ,JWRK)))
! sums:
           ZSNOWHEAT_S(JJ) = ZSNOWHEAT_S(JJ) + PEK%TSNOW%HEAT (JJ,JWRK)
           ZSNOWSWE_S(JJ)  = ZSNOWSWE_S(JJ)  + PEK%TSNOW%WSNOW(JJ,JWRK)
           ZSNOWDZ_S(JJ)   = ZSNOWDZ_S(JJ)   + ZSNOWDZ_X(JJ)
           ZSNOWAGE_S(JJ)  = ZSNOWAGE_S(JJ)  + PEK%TSNOW%AGE (JJ,JWRK)*ZSNOWDZ_X(JJ)
        ENDIF
     ENDDO
  ENDDO

! bulk properties:

  ZSNOWDMIN_S    = XSNOWDMIN*INLVLS         ! min total depth
  ZSNOWSWE_MIN   = ZSNOWDMIN_S*XRHOSMIN_ES  ! min total SWE

  ZCDT_MIN       = XCI*(ZCHECK_TEMP + ZCHECK_TEMP_DT -XTT) 
  
  ZSNOWRHO_X(:)  = 0.
  ZSNOWTEMP_X(:) = 0.
                         
  DO JJ=1,SIZE(PEK%TSNOW%WSNOW(:,:),1)
     IF(LCOLD(JJ))THEN
           
        ZSNOWRHO_X(JJ) = MAX(XRHOSMIN_ES,MIN(XRHOSMAX_ES, ZSNOWSWE_S(JJ)/  &
                         MAX(ZSNOWDMIN_S,ZSNOWDZ_S(JJ)) ))  ! limited
        ZSNOWDZ_S(JJ)  = ZSNOWSWE_S(JJ)/ZSNOWRHO_X(JJ)      ! can be small (true SWE sum)
        ZSNOWAGE_X(JJ) = MIN(ZSNOWAGE_X(JJ), MAX( 0., ZSNOWAGE_S(JJ)/      &
                         MAX(ZSNOWDMIN_S,ZSNOWDZ_S(JJ)) ) ) ! limited 

        ZSNOWTEMP_X(JJ)= XTT + MIN(0., MAX(                                            & 
                         ((ZSNOWHEAT_S(JJ)/MAX(ZSNOWSWE_S(JJ),ZSNOWSWE_MIN))           & 
                         + XLMTT), ZCDT_MIN ) )/XCI         ! limited                               
     ENDIF
  ENDDO

! dispatch avg bulk properites to each layer asumming const grid thickness (thin snowpack):
  
  DO JWRK=1,INLVLS
     DO JJ=1,SIZE(PEK%TSNOW%WSNOW(:,:),1)
        IF(LCOLD(JJ))THEN  
           PEK%TSNOW%WSNOW(JJ,JWRK)  = ZSNOWSWE_S(JJ) /INLVLS ! conserved
           PEK%TSNOW%HEAT (JJ,JWRK)  = ZSNOWHEAT_S(JJ)/INLVLS ! conserved
           PEK%TSNOW%RHO  (JJ,JWRK)  = ZSNOWRHO_X(JJ)         
           PEK%TSNOW%AGE(JJ,JWRK)    = ZSNOWAGE_X(JJ)
           DMK%XSNOWDZ  (JJ,JWRK)    = ZSNOWDZ_S(JJ)  /INLVLS
           DMK%XSNOWTEMP(JJ,JWRK)    = ZSNOWTEMP_X(JJ)
           DMK%XSNOWLIQ (JJ,JWRK)    = 0.
        ENDIF
     ENDDO
  ENDDO

! Remove snow if excessively cold & thin: no numerical shocks to soil
! will occur as long as snow is *thin*. If not thin, then let model stop
! because cold problem lies elsewhere.  
! NOTE: strictly speaking the min SWE occurs for XRHOSMIN_ES, or ould use
!       XRHOSMAX_ES to be less restrictive, or an avg value like 300 kg/m3...
!       But *more* restrictive (lower density threshold) means smaller flux
!       and mass correction added to soil!!

! first, save any mass or heat to possibly be removed:
  
  WHERE(LCOLD(:) .AND. ZSNOWSWE_S(:) <= ZSNOWSWE_MIN)
     ZMASS_COR(:) = ZSNOWSWE_S(:) /XRHOLW ! m
     ZHEAT_COR(:) = ZSNOWHEAT_S(:)        ! J m-2
  ELSEWHERE
     ZMASS_COR(:) = 0.
     ZHEAT_COR(:) = 0.
  ENDWHERE
  !
  !  Update snow mass budget and energy budget diags
  !  EVAPCOR is a soil mass sink, thus here it is a source (neg sign)
  !
  ZSNOW_MASS_BUDGET(:)   = ZSNOW_MASS_BUDGET(:) - ZMASS_COR(:)/PTSTEP
  PEVAPCOR(:)            = PEVAPCOR(:) - ZMASS_COR(:)/PTSTEP
  !
  ! PFLSN_COR is a diag when using MEB & CEB to close e-budget.
  ! When using CEB, it is also used to heat/cool the soil...but when using MEB,
  ! do that within this routine (below)
  !
  ZSNOW_ENERGY_BUDGET(:) = ZSNOW_ENERGY_BUDGET(:) - ZHEAT_COR(:)/PTSTEP
  PFLSN_COR(:)           = PFLSN_COR(:) - ZHEAT_COR(:)/PTSTEP
  !
  ! Add what is lost by the snow to the ground if MEB on:  
  ! Also, divide by the number of layers over which the
  ! correction is to be spread (NSOIL_COR >= 1)

  IF(OMEB)THEN
     IF(IO%CISBA == 'DIF')THEN
        DO JWRK=1,JNSOIL_COR
           DO JJ=1,SIZE(PEK%XTG,1)
              PEK%XTG(JJ,JWRK) = PEK%XTG(JJ,JWRK) + (ZHEAT_COR(JJ)/            &
                                (PSOILHCAPZ(JJ,JWRK)*PDZG(JJ,JWRK)))/JNSOIL_COR
           ENDDO
        ENDDO
     ELSE
        PEK%XTG(:,1) = PEK%XTG(:,1) +  ZHEAT_COR(:)/PCT(:)
     ENDIF
  ENDIF
  !
  ! remove snow and reset snow properties:  
  
  DO JWRK=1,INLVLS
     DO JJ=1,SIZE(PEK%TSNOW%WSNOW(:,:),1)
        IF(LCOLD(JJ) .AND. ZSNOWSWE_S(JJ) <= ZSNOWSWE_MIN)THEN  
           PEK%TSNOW%WSNOW(JJ,JWRK)  = 0.
           PEK%TSNOW%HEAT (JJ,JWRK)  = 0.
           PEK%TSNOW%RHO  (JJ,JWRK)  = XRHOSMIN_ES
           PEK%TSNOW%AGE(JJ,JWRK)    = 0.
           DMK%XSNOWDZ  (JJ,JWRK)    = 0.
           DMK%XSNOWTEMP(JJ,JWRK)    = XTT
           DMK%XSNOWLIQ (JJ,JWRK)    = 0.
        ENDIF
     ENDDO
  ENDDO
  !
  ! check suspicious low temperature
  !
  DO JWRK=1,INLVLS
    !
    DO JJ=1,SIZE(PEK%TSNOW%WSNOW,1)
      !
      IF (PEK%TSNOW%WSNOW(JJ,JWRK)>0.0) THEN
        !
        IF (DMK%XSNOWTEMP(JJ,JWRK)<ZCHECK_TEMP) THEN
          WRITE(*,*) 'Suspicious low temperature :',DMK%XSNOWTEMP(JJ,JWRK)
          WRITE(*,*) 'At point and location      :',JJ,'LAT=',G%XLAT(JJ),'LON=',G%XLON(JJ)
          WRITE(*,*) 'At snow level / total layer:',JWRK,'/',INLVLS
          WRITE(*,*) 'SNOW MASS BUDGET (kg/m2/s) :',ZSNOW_MASS_BUDGET(JJ)
          WRITE(*,*) 'SNOW ENERGY BUDGET (W/m2)  :',ZSNOW_ENERGY_BUDGET(JJ)
          WRITE(*,*) 'SWE BY LAYER      (kg/m2)  :',PEK%TSNOW%WSNOW (JJ,1:INLVLS)
          WRITE(*,*) 'DEKTH BY LAYER      (m)    :',DMK%XSNOWDZ  (JJ,1:INLVLS)
          WRITE(*,*) 'DENSITY BY LAYER   (kg/m3) :',PEK%TSNOW%RHO(JJ,1:INLVLS)
          WRITE(*,*) 'TEMPERATURE BY LAYER (K)   :',DMK%XSNOWTEMP(JJ,1:INLVLS)
          CALL ABOR1_SFX('SNOW3L_ISBA: Suspicious low temperature')
        ENDIF
        !
      ELSE
        !
        !Prognostic variables forced to XUNDEF for correct outputs
        DMK%XSNOWDZ(JJ,JWRK)=XUNDEF
        ! Careful : to compute average surface temperature in ISBA_SNOW_AGR
        ! PSNOWTEMP(JJ,1) is required when PPSN(JJ)>0 even if PSNOWSWE(JJ,1)==0
        ! (vanishing snowpack)
        IF (.NOT.(PEK%XPSN(JJ)>0.0.AND.JWRK==1)) DMK%XSNOWTEMP(JJ,JWRK) = XUNDEF
        DMK%XSNOWLIQ  (JJ,JWRK) = XUNDEF
        PEK%TSNOW%HEAT(JJ,JWRK) = XUNDEF
        PEK%TSNOW%RHO (JJ,JWRK) = XUNDEF
        PEK%TSNOW%AGE (JJ,JWRK) = XUNDEF
        IF (PEK%TSNOW%SCHEME=='CRO') THEN
          PEK%TSNOW%DIAMOPT(JJ,JWRK) = XUNDEF
          PEK%TSNOW%SPHERI(JJ,JWRK) = XUNDEF
          PEK%TSNOW%HIST (JJ,JWRK) = XUNDEF
        END IF
      ENDIF
    ENDDO
  ENDDO
  !
  DEK%XSWNET_N(:)  = ZSWNET_N(:) 
  DEK%XSWNET_NS(:) = ZSWNET_NS(:)
  DEK%XLWNET_N(:)  = ZLWNET_N(:)
  !
  ! ===============================================================
  !
ENDIF
!
IF (LHOOK) CALL DR_HOOK('SNOW3L_ISBA',1,ZHOOK_HANDLE)
!
 CONTAINS
!
!================================================================
SUBROUTINE CALL_MODEL(KSIZE1,KSIZE2,KSIZE3,KSIZE4,KMASK)
!
IMPLICIT NONE
!
INTEGER, INTENT(IN) :: KSIZE1
INTEGER, INTENT(IN) :: KSIZE2
INTEGER, INTENT(IN) :: KSIZE3
INTEGER, INTENT(IN) :: KSIZE4
INTEGER, DIMENSION(:), INTENT(IN) :: KMASK
!
REAL, DIMENSION(KSIZE1,KSIZE2) :: ZP_SNOWSWE
REAL, DIMENSION(KSIZE1,KSIZE2) :: ZP_SNOWDZ
REAL, DIMENSION(KSIZE1,KSIZE2) :: ZP_SNOWRHO
REAL, DIMENSION(KSIZE1,KSIZE2) :: ZP_SNOWHEAT
REAL, DIMENSION(KSIZE1,KSIZE2) :: ZP_SNOWTEMP
REAL, DIMENSION(KSIZE1,KSIZE2) :: ZP_SNOWLIQ
REAL, DIMENSION(KSIZE1,KSIZE2) :: ZP_SNOWDIAMOPT
REAL, DIMENSION(KSIZE1,KSIZE2) :: ZP_SNOWSPHERI
REAL, DIMENSION(KSIZE1,KSIZE2) :: ZP_SNOWHIST
REAL, DIMENSION(KSIZE1,KSIZE2) :: ZP_SNOWAGE
REAL, DIMENSION(KSIZE1,KSIZE2,NIMPUR) :: ZP_SNOWIMPUR
REAL, DIMENSION(KSIZE1,NIMPUR) :: ZP_IMPWET
REAL, DIMENSION(KSIZE1,NIMPUR) :: ZP_IMPDRY
REAL, DIMENSION(KSIZE1,KSIZE4) :: ZP_BLOWSNW
REAL, DIMENSION(KSIZE1)        :: ZP_VFRIC_T
REAL, DIMENSION(KSIZE1)        :: ZP_SNOWALB
REAL, DIMENSION(KSIZE1)        :: ZP_SWNETSNOW
REAL, DIMENSION(KSIZE1)        :: ZP_SWNETSNOWS
REAL, DIMENSION(KSIZE1)        :: ZP_LWNETSNOW
REAL, DIMENSION(KSIZE1)        :: ZP_PS
REAL, DIMENSION(KSIZE1)        :: ZP_SRSNOW
REAL, DIMENSION(KSIZE1)        :: ZP_UNLOAD
REAL, DIMENSION(KSIZE1)        :: ZP_RRSNOW
REAL, DIMENSION(KSIZE1)        :: ZP_PSN3L
REAL, DIMENSION(KSIZE1)        :: ZP_TAR
REAL, DIMENSION(KSIZE1)        :: ZP_TAC
REAL, DIMENSION(KSIZE1)        :: ZP_CT
REAL, DIMENSION(KSIZE1,KSIZE3) :: ZP_TG
REAL, DIMENSION(KSIZE1,KSIZE3) :: ZP_D_G
REAL, DIMENSION(KSIZE1,KSIZE3) :: ZP_DZG
REAL, DIMENSION(KSIZE1,KSIZE3) :: ZP_SOILHCAPZ
REAL, DIMENSION(KSIZE1)        :: ZP_SOILD
REAL, DIMENSION(KSIZE1)        :: ZP_DELHEATG
REAL, DIMENSION(KSIZE1)        :: ZP_DELHEATG_SFC
REAL, DIMENSION(KSIZE1)        :: ZP_SW_RAD
REAL, DIMENSION(KSIZE1)        :: ZP_QA
REAL, DIMENSION(KSIZE1)        :: ZP_LVTT
REAL, DIMENSION(KSIZE1)        :: ZP_LSTT
REAL, DIMENSION(KSIZE1)        :: ZP_VMOD
REAL, DIMENSION(KSIZE1)        :: ZP_LW_RAD
REAL, DIMENSION(KSIZE1)        :: ZP_RHOA
REAL, DIMENSION(KSIZE1)        :: ZP_UREF
REAL, DIMENSION(KSIZE1)        :: ZP_EXNS
REAL, DIMENSION(KSIZE1)        :: ZP_EXNA
REAL, DIMENSION(KSIZE1)        :: ZP_DIRCOSZW
REAL, DIMENSION(KSIZE1)        :: ZP_SLOPEDIR
REAL, DIMENSION(KSIZE1)        :: ZP_ZREF
REAL, DIMENSION(KSIZE1)        :: ZP_Z0NAT
REAL, DIMENSION(KSIZE1)        :: ZP_Z0HNAT
REAL, DIMENSION(KSIZE1)        :: ZP_Z0EFF
REAL, DIMENSION(KSIZE1)        :: ZP_ALB
REAL, DIMENSION(KSIZE1)        :: ZP_SOILCOND
REAL, DIMENSION(KSIZE1)        :: ZP_THRUFAL
REAL, DIMENSION(KSIZE1)        :: ZP_MELTSTOT
REAL, DIMENSION(KSIZE1)        :: ZP_SNREFREEZ
REAL, DIMENSION(KSIZE1)        :: ZP_GRNDFLUX
REAL, DIMENSION(KSIZE1)        :: ZP_FLSN_COR
REAL, DIMENSION(KSIZE1)        :: ZP_RESTOREN
REAL, DIMENSION(KSIZE1)        :: ZP_EVAPCOR
REAL, DIMENSION(KSIZE1)        :: ZP_SOILCOR
REAL, DIMENSION(KSIZE1)        :: ZP_GFLXCOR
REAL, DIMENSION(KSIZE1)        :: ZP_RNSNOW
REAL, DIMENSION(KSIZE1)        :: ZP_HSNOW
REAL, DIMENSION(KSIZE1)        :: ZP_GFLUXSNOW
REAL, DIMENSION(KSIZE1)        :: ZP_DELHEATN
REAL, DIMENSION(KSIZE1)        :: ZP_DELHEATN_SFC
REAL, DIMENSION(KSIZE1)        :: ZP_DELPHASEN
REAL, DIMENSION(KSIZE1)        :: ZP_DELPHASEN_SFC
REAL, DIMENSION(KSIZE1)        :: ZP_SNOWSFCH
REAL, DIMENSION(KSIZE1)        :: ZP_HPSNOW
REAL, DIMENSION(KSIZE1)        :: ZP_LES3L
REAL, DIMENSION(KSIZE1)        :: ZP_LEL3L
REAL, DIMENSION(KSIZE1)        :: ZP_EVAP
REAL, DIMENSION(KSIZE1)        :: ZP_SNDRIFT
REAL, DIMENSION(KSIZE1)        :: ZP_RI
REAL, DIMENSION(KSIZE1)        :: ZP_QS
REAL, DIMENSION(KSIZE1)        :: ZP_EMISNOW
REAL, DIMENSION(KSIZE1)        :: ZP_CDSNOW
REAL, DIMENSION(KSIZE1)        :: ZP_USTARSNOW
REAL, DIMENSION(KSIZE1)        :: ZP_CHSNOW
REAL, DIMENSION(KSIZE1)        :: ZP_SNOWHMASS
REAL, DIMENSION(KSIZE1)        :: ZP_VEGTYPE
REAL, DIMENSION(KSIZE1)        :: ZP_FOREST
REAL, DIMENSION(KSIZE1)        :: ZP_PEW_A_COEF
REAL, DIMENSION(KSIZE1)        :: ZP_PEW_B_COEF
REAL, DIMENSION(KSIZE1)        :: ZP_PET_A_COEF
REAL, DIMENSION(KSIZE1)        :: ZP_PET_B_COEF
REAL, DIMENSION(KSIZE1)        :: ZP_PEQ_A_COEF
REAL, DIMENSION(KSIZE1)        :: ZP_PEQ_B_COEF
REAL, DIMENSION(KSIZE1)        :: ZP_ZENITH
REAL, DIMENSION(KSIZE1)        :: ZP_AZIM
REAL, DIMENSION(KSIZE1)        :: ZP_LAT,ZP_LON
REAL, DIMENSION(KSIZE1)        :: ZP_PSN_INV
REAL, DIMENSION(KSIZE1)        :: ZP_PSN
REAL, DIMENSION(KSIZE1)        :: ZP_PSN_GFLXCOR
REAL, DIMENSION(KSIZE1)        :: ZP_WORK
REAL, DIMENSION(KSIZE1)        :: ZP_TG_CPL
REAL, DIMENSION(KSIZE1)        :: ZP_VEG
!
REAL, DIMENSION(KSIZE1,KSIZE2) :: ZP_SNOWDEND
REAL, DIMENSION(KSIZE1,KSIZE2) :: ZP_SNOWSPHER
REAL, DIMENSION(KSIZE1,KSIZE2) :: ZP_SNOWSIZE
REAL, DIMENSION(KSIZE1,KSIZE2) :: ZP_SNOWSSA
REAL, DIMENSION(KSIZE1,KSIZE2) :: ZP_SNOWTYPEMEPRA
REAL, DIMENSION(KSIZE1,KSIZE2) :: ZP_SNOWRAM
REAL, DIMENSION(KSIZE1,KSIZE2) :: ZP_SNOWSHEAR
REAL, DIMENSION(KSIZE1,KSIZE2) :: ZP_ACC_RAT
REAL, DIMENSION(KSIZE1,KSIZE2) :: ZP_NAT_RAT
!
REAL, DIMENSION(KSIZE1) :: ZP_SNDPT_12H
REAL, DIMENSION(KSIZE1) :: ZP_SNDPT_1DY
REAL, DIMENSION(KSIZE1) :: ZP_SNDPT_3DY
REAL, DIMENSION(KSIZE1) :: ZP_SNDPT_5DY
REAL, DIMENSION(KSIZE1) :: ZP_SNDPT_7DY
REAL, DIMENSION(KSIZE1) :: ZP_SNSWE_1DY
REAL, DIMENSION(KSIZE1) :: ZP_SNSWE_3DY
REAL, DIMENSION(KSIZE1) :: ZP_SNSWE_5DY
REAL, DIMENSION(KSIZE1) :: ZP_SNSWE_7DY
REAL, DIMENSION(KSIZE1) :: ZP_SNRAM_SONDE
REAL, DIMENSION(KSIZE1) :: ZP_SN_WETTHCKN
REAL, DIMENSION(KSIZE1) :: ZP_SN_REFRZNTHCKN
REAL, DIMENSION(KSIZE1) :: ZP_DEP_HIG
REAL, DIMENSION(KSIZE1) :: ZP_DEP_MOD
REAL, DIMENSION(KSIZE1) :: ZP_DEP_SUP
REAL, DIMENSION(KSIZE1) :: ZP_DEP_TOT
REAL, DIMENSION(KSIZE1) :: ZP_DEP_HUM
REAL, DIMENSION(KSIZE1) :: ZP_ACC_LEV
REAL, DIMENSION(KSIZE1) :: ZP_NAT_LEV
REAL, DIMENSION(KSIZE1) :: ZP_PRO_SUP_TYP
REAL, DIMENSION(KSIZE1) :: ZP_PRO_INF_TYP
REAL, DIMENSION(KSIZE1) :: ZP_AVA_TYP
REAL, DIMENSION(KSIZE1) :: ZP_SNOWMAK
REAL, DIMENSION(KSIZE1,KSIZE2,NIMPUR) :: ZP_SNOWIMP_CONC !F.T
REAL, DIMENSION(KSIZE1,SIZE(P_DIR_SW,2)) :: ZP_DIR_SW !F.T
REAL, DIMENSION(KSIZE1,SIZE(P_DIR_SW,2)) :: ZP_SCA_SW !F.T
REAL, DIMENSION(KSIZE1,SIZE(P_DIR_SW,2)) :: ZP_SPEC_ALB !F.T
REAL, DIMENSION(KSIZE1,SIZE(P_DIR_SW,2)) :: ZP_DIFF_RATIO !F.T
REAL, DIMENSION(KSIZE1,SIZE(P_DIR_SW,2)) :: ZP_SPEC_TOT !F.T
!
REAL, PARAMETER :: ZDEPTHABS = 0.60 ! m
!
INTEGER :: JWRK, JJ, JI
!
LOGICAL :: GCOMPUTECRODIAG ! flag to compute Crocus-MEPRA diagnostics
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
IF (LHOOK) CALL DR_HOOK('SNOW3L_ISBA:CALL_MODEL',0,ZHOOK_HANDLE)
!
! Initialize:
!
GCOMPUTECRODIAG = .FALSE.
ZP_PSN_GFLXCOR(:)  = 0.
ZP_WORK(:)         = 0.
ZP_SOILD(:)        = 0.
ZP_UNLOAD(:)       = 0.
!
!
! pack the variables
!
DO JWRK=1,KSIZE2
   DO JJ=1,KSIZE1
      JI = KMASK(JJ)
      ZP_SNOWSWE (JJ,JWRK) = PEK%TSNOW%WSNOW(JI,JWRK)
      ZP_SNOWRHO (JJ,JWRK) = PEK%TSNOW%RHO  (JI,JWRK)
      ZP_SNOWHEAT(JJ,JWRK) = PEK%TSNOW%HEAT (JI,JWRK)
      ZP_SNOWAGE (JJ,JWRK) = PEK%TSNOW%AGE  (JI,JWRK)
      ZP_SNOWTEMP(JJ,JWRK) = DMK%XSNOWTEMP  (JI,JWRK)
      ZP_SNOWLIQ (JJ,JWRK) = DMK%XSNOWLIQ   (JI,JWRK)
      ZP_SNOWDZ  (JJ,JWRK) = DMK%XSNOWDZ    (JI,JWRK)
   ENDDO
ENDDO
!
IF (PEK%TSNOW%SCHEME=='CRO') THEN
  IF (IO%LSNOWPAPPUS) THEN
    DO JJ=1,KSIZE1
      JI = KMASK(JJ)
      ZP_VFRIC_T    (JJ) = DMK%XVFRIC_T_PAPPUS   (JI)
    ENDDO
  ENDIF
  DO JWRK=1,KSIZE2
    DO JJ=1,KSIZE1
      JI = KMASK(JJ)
      ZP_SNOWDIAMOPT(JJ,JWRK) = PEK%TSNOW%DIAMOPT(JI,JWRK)
      ZP_SNOWSPHERI (JJ,JWRK) = PEK%TSNOW%SPHERI (JI,JWRK)
      ZP_SNOWHIST   (JJ,JWRK) = PEK%TSNOW%HIST   (JI,JWRK)
    ENDDO
  ENDDO

  DO JJ=1,KSIZE1
    JI = KMASK(JJ)
    ZP_DEP_SUP    (JJ) = PEK%TSNOW%DEP_SUP    (JI)
    ZP_DEP_TOT    (JJ) = PEK%TSNOW%DEP_TOT    (JI)
    ZP_DEP_HUM    (JJ) = PEK%TSNOW%DEP_HUM    (JI)
    ZP_NAT_LEV    (JJ) = PEK%TSNOW%NAT_LEV    (JI)
    ZP_PRO_SUP_TYP(JJ) = PEK%TSNOW%PRO_SUP_TYP(JI)
    ZP_AVA_TYP    (JJ) = PEK%TSNOW%AVA_TYP    (JI)
  ENDDO

  DO JIMP=1,NIMPUR
    DO JWRK=1,KSIZE2
      DO JJ=1,KSIZE1
        JI = KMASK(JJ)
        ZP_SNOWIMPUR(JJ,JWRK,JIMP) =PEK%TSNOW%IMPUR(JI,JWRK,JIMP)
      ENDDO
    ENDDO
  ENDDO

  DO JIMP=1,NIMPUR
    DO JJ=1,KSIZE1
      JI = KMASK(JJ)
      ZP_IMPWET(JJ,JIMP)=PIMPWET(JI,JIMP)
      ZP_IMPDRY(JJ,JIMP)=PIMPDRY(JI,JIMP)
    ENDDO
  ENDDO !end BC merge oubli

  DO JWRK=1,KSIZE4
    DO JJ=1,KSIZE1
      JI = KMASK(JJ)
      ZP_BLOWSNW(JJ,JWRK) = ZBLOWSNW(JI,JWRK)
    ENDDO
  ENDDO

ELSE
  DO JWRK=1,KSIZE2
    DO JJ=1,KSIZE1
      ZP_SNOWDIAMOPT(JJ,JWRK) = XUNDEF
      ZP_SNOWSPHERI (JJ,JWRK) = XUNDEF
      ZP_SNOWHIST   (JJ,JWRK) = XUNDEF
    ENDDO
  ENDDO

  DO JIMP=1,NIMPUR
    DO JJ=1,KSIZE1
      ZP_IMPWET(JJ,JIMP)=XUNDEF
      ZP_IMPDRY(JJ,JIMP)=XUNDEF
      DO JWRK=1,KSIZE2
        ZP_SNOWIMPUR(JJ,JWRK,JIMP) = XUNDEF
      ENDDO
    ENDDO
  ENDDO
  !
  DO JWRK=1,KSIZE4
    DO JJ=1,KSIZE1
      ZP_BLOWSNW(JJ,JWRK) = XUNDEF
    ENDDO
  ENDDO
  !
ENDIF
! 
DO JWRK=1,KSIZE3
   DO JJ=1,KSIZE1
      JI                    = KMASK           (JJ)
      ZP_TG       (JJ,JWRK) = PTG        (JI,JWRK)
      ZP_D_G      (JJ,JWRK) = PD_G       (JI,JWRK)
      ZP_SOILHCAPZ(JJ,JWRK) = PSOILHCAPZ (JI,JWRK)
   ENDDO
ENDDO
!
IF (IO%CISBA=='DIF'.OR.OMEB) THEN
  DO JWRK=1,KSIZE3
    DO JJ=1,KSIZE1
      JI                    = KMASK           (JJ)
      ZP_DZG      (JJ,JWRK) = PDZG       (JI,JWRK)
    ENDDO
  ENDDO
ENDIF
!
IF(IO%LMEB_INT_USFC)THEN
  DO JJ=1,KSIZE1
    JI = KMASK(JJ)
    ZP_UNLOAD(JJ) = DEK%XUNLOAD(JI)
  ENDDO
ENDIF
!
DO JJ=1,KSIZE1
   !
   JI = KMASK(JJ)
   !
   ZP_LVTT        (JJ) = PK%XLVTT         (JI)
   ZP_LSTT        (JJ) = PK%XLSTT         (JI)
   ZP_SNOWALB     (JJ) = PEK%TSNOW%ALB    (JI)
   ZP_PSN3L       (JJ) = PEK%XPSN         (JI)
   ZP_Z0NAT       (JJ) = DK%XZ0           (JI)
   ZP_Z0HNAT      (JJ) = DK%XZ0H          (JI)
   ZP_Z0EFF       (JJ) = DK%XZ0EFF        (JI)
   ZP_PS          (JJ) = PPS              (JI)
   ZP_SRSNOW      (JJ) = PSR              (JI)
   ZP_CT          (JJ) = PCT              (JI)
   ZP_TAR         (JJ) = PTAR             (JI)
   ZP_TAC         (JJ) = PTAC             (JI)
   ZP_DELHEATG    (JJ) = DEK%XDELHEATG    (JI)
   ZP_DELHEATG_SFC(JJ) = DEK%XDELHEATG_SFC(JI)
   ZP_SW_RAD      (JJ) = PSW_RAD          (JI)
   ZP_QA          (JJ) = PQA              (JI)
   ZP_VMOD        (JJ) = PVMOD            (JI)
   ZP_LW_RAD      (JJ) = PLW_RAD          (JI)
   ZP_RHOA        (JJ) = PRHOA            (JI)
   ZP_UREF        (JJ) = PUREF            (JI)
   ZP_EXNS        (JJ) = PEXNS            (JI)
   ZP_EXNA        (JJ) = PEXNA            (JI)
   ZP_DIRCOSZW    (JJ) = PDIRCOSZW        (JI)
   ZP_SLOPEDIR    (JJ) = PSLOPEDIR        (JI)
   ZP_ZREF        (JJ) = PZREF            (JI)
   ZP_ALB         (JJ) = PALB             (JI)
   ZP_RRSNOW      (JJ) = ZRRSNOW          (JI)
   ZP_SOILCOND    (JJ) = ZSOILCOND        (JI)
   !  
   ZP_PEW_A_COEF(JJ) = PPEW_A_COEF(JI)
   ZP_PEW_B_COEF(JJ) = PPEW_B_COEF(JI)
   ZP_PET_A_COEF(JJ) = PPET_A_COEF(JI)
   ZP_PEQ_A_COEF(JJ) = PPEQ_A_COEF(JI)
   ZP_PET_B_COEF(JJ) = PPET_B_COEF(JI)
   ZP_PEQ_B_COEF(JJ) = PPEQ_B_COEF(JI)
   !
   ZP_LAT(JJ) = G%XLAT(JI)
   ZP_LON(JJ) = G%XLON(JI)

   ZP_ZENITH    (JJ) = PZENITH    (JI)
   ZP_AZIM      (JJ) = PAZIM      (JI)
   !
   ZP_GRNDFLUX    (JJ) = PGRNDFLUX        (JI)
   ZP_RNSNOW      (JJ) = DMK%XRNSNOW      (JI)
   ZP_HSNOW       (JJ) = DMK%XHSNOW       (JI)
   ZP_DELHEATN    (JJ) = DEK%XDELHEATN    (JI)
   ZP_DELHEATN_SFC(JJ) = DEK%XDELHEATN_SFC(JI)
   ZP_SNOWSFCH    (JJ) = PSNOWSFCH        (JI)
   ZP_HPSNOW      (JJ) = DMK%XHPSNOW      (JI)
   ZP_LES3L       (JJ) = PLES3L           (JI)
   ZP_LEL3L       (JJ) = PLEL3L           (JI)
   ZP_EVAP        (JJ) = PEVAP            (JI)
   ZP_EMISNOW     (JJ) = PEK%TSNOW%EMIS   (JI)
   ZP_SWNETSNOW   (JJ) = ZSWNET_N         (JI)
   ZP_SWNETSNOWS  (JJ) = ZSWNET_NS        (JI)
   ZP_LWNETSNOW   (JJ) = ZLWNET_N         (JI)
   !
   ZP_GFLUXSNOW    (JJ) = DMK%XGFLUXSNOW    (JI)
   ZP_SNOWHMASS    (JJ) = DMK%XSNOWHMASS    (JI)
   ZP_RESTOREN     (JJ) = DEK%XRESTOREN     (JI)
   ZP_DELPHASEN    (JJ) = DEK%XDELPHASEN    (JI)
   ZP_DELPHASEN_SFC(JJ) = DEK%XDELPHASEN_SFC(JI)
   !
   ZP_SNOWMAK(JJ)  = ZSNOWMAK(JI)
   !
ENDDO
!
DO JWRK=1,SIZE(P_DIR_SW,2)
  DO JJ=1,KSIZE1
    JI = KMASK(JJ)
    ZP_DIR_SW(JJ,JWRK)=P_DIR_SW(JI,JWRK)
    ZP_SCA_SW(JJ,JWRK)=P_SCA_SW(JI,JWRK)
  ENDDO
ENDDO
!
DO JJ=1,KSIZE1
   JI = KMASK(JJ)
   ZP_VEGTYPE (JJ) = 0.
   ZP_FOREST  (JJ) = 0.
   DO JVEG = 1, NVEGTYPE+NVEG_IRR
     JK = JVEG
     IF (JVEG > NVEGTYPE) JK = NPAR_VEG_IRR_USE( JVEG - NVEGTYPE )
     IF ( JK == NVT_SNOW ) ZP_VEGTYPE (JJ) = ZP_VEGTYPE (JJ) + PVEGTYPE (JI,JVEG)
     IF ( JK == NVT_TEBD .OR. JK == NVT_TRBE .OR. JK == NVT_BONE .OR. JK == NVT_TRBD .OR. JK == NVT_TEBE .OR. &
          JK == NVT_TENE .OR. JK == NVT_BOBD .OR. JK == NVT_BOND .OR. JK == NVT_SHRB )                        &
       ZP_FOREST  (JJ) = ZP_FOREST  (JJ) + PVEGTYPE (JI,JVEG)
   ENDDO
   !
   ZP_VEG(JJ) = PEK%XVEG(JI)
   !
ENDDO
!
!
! ===============================================================
! conversion of snow heat from J/m3 into J/m2
WHERE(ZP_SNOWSWE(:,:)>0.) &
  ZP_SNOWHEAT(:,:) = ZP_SNOWHEAT(:,:) / ZP_SNOWRHO (:,:) * ZP_SNOWSWE (:,:)
! ===============================================================
!
ZP_PSN_INV(:)       = 0.0
ZP_PSN(:)           = ZP_PSN3L(:)
!
ZP_TG_CPL(:)        = ZP_TG(:,1)
!
IF(OMEB)THEN
!
!   MEB (case of imposed surface fluxes)
!   - Prepare inputs for explicit snow scheme(s):
!     If using MEB, these are INPUTs ONLY:
!     divide fluxes by snow fraction to make "snow-relative"
!
   ZP_PSN(:)          = MAX(1.E-4, ZP_PSN3L(:))
   ZP_PSN_INV(:)      = 1.0/ZP_PSN(:)
   ZP_SRSNOW(:)       = ZP_SRSNOW(:) - ZP_UNLOAD(:) ! PSR = only atmosphere snowfall (not unloading)
!
   ZP_RNSNOW(:)       = ZP_RNSNOW(:)       *ZP_PSN_INV(:)
   ZP_SWNETSNOW(:)    = ZP_SWNETSNOW(:)    *ZP_PSN_INV(:)
   ZP_SWNETSNOWS(:)   = ZP_SWNETSNOWS(:)   *ZP_PSN_INV(:)
   ZP_LWNETSNOW(:)    = ZP_LWNETSNOW(:)    *ZP_PSN_INV(:)
   ZP_HSNOW(:)        = ZP_HSNOW(:)        *ZP_PSN_INV(:)
   ZP_GFLUXSNOW(:)    = ZP_GFLUXSNOW(:)    *ZP_PSN_INV(:)
   ZP_RESTOREN(:)     = ZP_RESTOREN(:)     *ZP_PSN_INV(:)
   ZP_SNOWHMASS(:)    = ZP_SNOWHMASS(:)    *ZP_PSN_INV(:)
   ZP_LES3L(:)        = ZP_LES3L(:)        *ZP_PSN_INV(:)
   ZP_LEL3L(:)        = ZP_LEL3L(:)        *ZP_PSN_INV(:)
   ZP_GRNDFLUX(:)     = ZP_GRNDFLUX(:)     *ZP_PSN_INV(:)
   ZP_EVAP(:)         = ZP_EVAP(:)         *ZP_PSN_INV(:)
   ZP_HPSNOW(:)       = ZP_HPSNOW(:)       *ZP_PSN_INV(:)
   ZP_DELHEATN(:)     = ZP_DELHEATN(:)     *ZP_PSN_INV(:)
   ZP_DELHEATN_SFC(:) = ZP_DELHEATN_SFC(:) *ZP_PSN_INV(:)
   ZP_DELPHASEN(:)    = ZP_DELPHASEN(:)    *ZP_PSN_INV(:)
   ZP_DELPHASEN_SFC(:)= ZP_DELPHASEN_SFC(:)*ZP_PSN_INV(:)
   ZP_SNOWSFCH(:)     = ZP_SNOWSFCH(:)     *ZP_PSN_INV(:)
!
   ZP_SRSNOW(:)       = ZP_SRSNOW(:)       *ZP_PSN_INV(:)
   ZP_UNLOAD(:)       = ZP_UNLOAD(:)       *ZP_PSN_INV(:)
   ZP_RRSNOW(:)       = ZP_RRSNOW(:)       *ZP_PSN_INV(:)
!
   DO JJ=1,KSIZE2
      DO JI=1,KSIZE1
         ZP_SNOWSWE(JI,JJ)  = ZP_SNOWSWE(JI,JJ) *ZP_PSN_INV(JI)
         ZP_SNOWHEAT(JI,JJ) = ZP_SNOWHEAT(JI,JJ)*ZP_PSN_INV(JI)
         ZP_SNOWDZ(JI,JJ)   = ZP_SNOWDZ(JI,JJ)  *ZP_PSN_INV(JI)
      ENDDO
   ENDDO
!
ELSEIF(IO%CHORT=='CM6'.AND.IO%CISBA=='DIF'.AND.(.NOT.LSNOW_FRAC_TOT))THEN
   !
   ZP_SOILD(:) = 0.0
   ZP_WORK (:) = 0.0
   DO JJ=1,KSIZE3
      DO JI=1,KSIZE1
         IF(ZP_SOILD(JI) < ZDEPTHABS)THEN
              ZP_WORK (JI) = ZP_WORK(JI) + ZP_DZG(JI,JJ)*ZP_TG(JI,JJ) ! K
              ZP_SOILD(JI) = ZP_D_G(JI,JJ)
         ENDIF
      ENDDO
   ENDDO
   !
   ZP_TG_CPL(:)=(1.0-ZP_VEG(:))*ZP_TG(:,1)+ZP_VEG(:)*ZP_WORK(:)/ZP_SOILD(:)
   !
ENDIF
!
! Call snow schemes :
!  
IF (PEK%TSNOW%SCHEME=='CRO') THEN 
   !
   ! ------------------------
   ! Main call to Crocus
      CALL SNOWCRO(IO%CSNOWRES, TPTIME, OMEB, HIMPLICIT_WIND,    &
                ZP_PEW_A_COEF, ZP_PEW_B_COEF, ZP_PET_A_COEF, ZP_PEQ_A_COEF,   &
                ZP_PET_B_COEF, ZP_PEQ_B_COEF, ZP_SNOWSWE, ZP_SNOWRHO,         &
                ZP_SNOWHEAT, ZP_SNOWALB, ZP_SNOWDIAMOPT, ZP_SNOWSPHERI,       &
                ZP_SNOWHIST, ZP_SNOWAGE,ZP_SNOWIMPUR, PTSTEP, ZP_PS,ZP_SRSNOW,&
                ZP_UNLOAD, ZP_RRSNOW, ZP_PSN3L, ZP_TAR, ZP_TG(:,1), ZP_SW_RAD,&
                ZP_QA,ZP_VMOD, ZP_LW_RAD, ZP_RHOA, ZP_UREF, ZP_EXNS, ZP_EXNA, &
                ZP_DIRCOSZW, ZP_SLOPEDIR, ZP_ZREF, ZP_Z0NAT, ZP_Z0EFF,        &
                ZP_Z0HNAT, ZP_ALB, ZP_SOILCOND, ZP_D_G(:,1), ZP_SNOWLIQ,      &
                ZP_SNOWTEMP, ZP_SNOWDZ, ZP_THRUFAL, ZP_GRNDFLUX, ZP_EVAPCOR,  &
                ZP_GFLXCOR, ZP_SWNETSNOW, ZP_SWNETSNOWS, ZP_LWNETSNOW,        &
                ZP_RNSNOW, ZP_HSNOW, ZP_GFLUXSNOW, ZP_HPSNOW, ZP_LES3L,       &
                ZP_LEL3L, ZP_EVAP, ZP_SNDRIFT, ZP_RI, ZP_EMISNOW, ZP_CDSNOW,  &
                ZP_USTARSNOW, ZP_CHSNOW, ZP_SNOWHMASS, ZP_QS, ZP_VEGTYPE,     &
                ZP_ZENITH, ZP_AZIM, ZP_LAT, ZP_LON, ZP_BLOWSNW, IO%CSNOWDRIFT,&
                IO%CSNOWFPAPPUS, IO%LSNOWDRIFT_SUBLIM, IO%LSNOW_ABS_ZENITH,   &
                IO%CSNOWMETAMO,IO%CSNOWRAD,IO%LATMORAD,ZP_DIR_SW, ZP_SCA_SW,  &
                ZP_SPEC_ALB, ZP_DIFF_RATIO,ZP_SPEC_TOT, ZP_RESTOREN,ZP_IMPWET,&
                ZP_IMPDRY, IO%CSNOWFALL, IO%CSNOWCOND, IO%CSNOWHOLD,          &
                IO%CSNOWCOMP, IO%CSNOWZREF,ZP_SNOWMAK, IO%LSNOWCOMPACT_BOOL,  &
                IO%LSNOWMAK_BOOL,IO%LSNOWTILLER,IO%LSELF_PROD,                &
                IO%LSNOWMAK_PROP, ZP_VFRIC_T)
!
  ZP_FLSN_COR(:) = 0.0
  ZP_SOILCOR (:) = 0.0
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Bertrand : Diag should be coded in crocus but I don't have time
  ZP_MELTSTOT (:) = 0.0
  ZP_SNREFREEZ(:) = 0.0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
#ifndef SFX_OL
  ! En couplé il faudra voir si on veut virer les diagnostics, les calculer tout le temps, ou trouver une autre solution
  GCOMPUTECRODIAG = (SIZE(DMK%XSNOWDEND)>0)
#else
  ! This condition should be applied relatively to the time since beginning of simulation, not to the absolute time.
  ! It is temporarily removed to allow for example daily outputs with runs starting at 6h
  ! However a correct condition should be implemented in a near future to reduce computation time
  GCOMPUTECRODIAG = (SIZE(DMK%XSNOWDEND)>0).AND.(MOD(PTIMEC,XTSTEP_OUTPUT) == 0. )
#endif
  !
  !Ajout test sur pas de temps de sortie
  IF (GCOMPUTECRODIAG) THEN
    CALL SNOWCRO_DIAG(IO%CSNOWHOLD, IO%CSNOWMETAMO, ZP_SNOWDZ, ZP_SNOWSWE, ZP_SNOWRHO, ZP_SNOWDIAMOPT, ZP_SNOWSPHERI,       &
                      ZP_SNOWAGE, ZP_SNOWHIST, ZP_SNOWTEMP, ZP_SNOWLIQ, ZP_DIRCOSZW, ZP_SNOWIMPUR, ZP_SNOWDEND,             &
                      ZP_SNOWSPHER,ZP_SNOWSIZE, ZP_SNOWSSA, ZP_SNOWTYPEMEPRA, ZP_SNOWRAM, ZP_SNOWSHEAR,                     &
                      ZP_ACC_RAT, ZP_NAT_RAT,                                                                               &
                      ZP_SNDPT_12H, ZP_SNDPT_1DY, ZP_SNDPT_3DY, ZP_SNDPT_5DY, ZP_SNDPT_7DY,                                 &
                      ZP_SNSWE_1DY, ZP_SNSWE_3DY, ZP_SNSWE_5DY, ZP_SNSWE_7DY,                                               &
                      ZP_SNRAM_SONDE, ZP_SN_WETTHCKN, ZP_SN_REFRZNTHCKN,ZP_SNOWIMP_CONC,                                    &
                      ZP_DEP_HIG, ZP_DEP_MOD, ZP_DEP_SUP, ZP_DEP_TOT, ZP_DEP_HUM,                                           &
                      ZP_ACC_LEV, ZP_NAT_LEV, ZP_PRO_SUP_TYP, ZP_PRO_INF_TYP, ZP_AVA_TYP)
  ENDIF
  !
ELSE
!
  CALL SNOW3L(IO%CSNOWRES, OMEB, HIMPLICIT_WIND,                           &
             ZP_PEW_A_COEF, ZP_PEW_B_COEF,                                 &
             ZP_PET_A_COEF, ZP_PEQ_A_COEF,ZP_PET_B_COEF, ZP_PEQ_B_COEF,    &
             ZP_SNOWSWE, ZP_SNOWRHO, ZP_SNOWHEAT, ZP_SNOWALB,              &
             ZP_SNOWAGE, PTSTEP, ZP_PS, ZP_SRSNOW, ZP_UNLOAD,              &
             ZP_RRSNOW, ZP_PSN3L, ZP_TAR, ZP_TAC, ZP_TG_CPL(:),            &
             ZP_SW_RAD, ZP_QA, ZP_VMOD, ZP_LW_RAD, ZP_RHOA, ZP_UREF,       &
             ZP_EXNS, ZP_EXNA, ZP_DIRCOSZW, ZP_ZREF, ZP_Z0NAT, ZP_Z0EFF,   &
             ZP_Z0HNAT, ZP_ALB, ZP_SOILCOND, ZP_D_G(:,1),                  &
             ZP_LVTT, ZP_LSTT, ZP_SNOWLIQ,                                 &
             ZP_SNOWTEMP, ZP_SNOWDZ, ZP_THRUFAL, ZP_MELTSTOT, ZP_SNREFREEZ,&
             ZP_GRNDFLUX, ZP_EVAPCOR, ZP_SOILCOR, ZP_GFLXCOR, ZP_SNOWSFCH, &
             ZP_DELHEATN, ZP_DELHEATN_SFC, ZP_DELPHASEN, ZP_DELPHASEN_SFC, &
             ZP_SWNETSNOW, ZP_SWNETSNOWS, ZP_LWNETSNOW, ZP_RESTOREN,       &
             ZP_RNSNOW, ZP_HSNOW, ZP_GFLUXSNOW, ZP_HPSNOW, ZP_LES3L,       &
             ZP_LEL3L, ZP_EVAP, ZP_SNDRIFT, ZP_RI,                         &
             ZP_EMISNOW, ZP_CDSNOW, ZP_USTARSNOW,                          &
             ZP_CHSNOW, ZP_SNOWHMASS, ZP_QS, ZP_VEGTYPE,  ZP_FOREST,       &
             ZP_ZENITH, IO%CSNOWDRIFT, IO%LSNOWDRIFT_SUBLIM                )
ENDIF
!
  IF(OMEB)THEN
!
! - reverse transform: back to surface-relative
!
     ZP_RNSNOW(:)       = ZP_RNSNOW(:)       /ZP_PSN_INV(:)
     ZP_SWNETSNOW(:)    = ZP_SWNETSNOW(:)    /ZP_PSN_INV(:)
     ZP_SWNETSNOWS(:)   = ZP_SWNETSNOWS(:)   /ZP_PSN_INV(:)
     ZP_LWNETSNOW(:)    = ZP_LWNETSNOW(:)    /ZP_PSN_INV(:)
     ZP_HSNOW(:)        = ZP_HSNOW(:)        /ZP_PSN_INV(:)
     ZP_LES3L(:)        = ZP_LES3L(:)        /ZP_PSN_INV(:)
     ZP_LEL3L(:)        = ZP_LEL3L(:)        /ZP_PSN_INV(:)
     ZP_GRNDFLUX(:)     = ZP_GRNDFLUX(:)     /ZP_PSN_INV(:)
     ZP_EVAP(:)         = ZP_EVAP(:)         /ZP_PSN_INV(:)
     ZP_HPSNOW(:)       = ZP_HPSNOW(:)       /ZP_PSN_INV(:)
     ZP_GFLUXSNOW(:)    = ZP_GFLUXSNOW(:)    /ZP_PSN_INV(:)
     ZP_DELHEATN(:)     = ZP_DELHEATN(:)     /ZP_PSN_INV(:)
     ZP_DELHEATN_SFC(:) = ZP_DELHEATN_SFC(:) /ZP_PSN_INV(:)
     ZP_DELPHASEN(:)    = ZP_DELPHASEN(:)    /ZP_PSN_INV(:)
     ZP_DELPHASEN_SFC(:)= ZP_DELPHASEN_SFC(:)/ZP_PSN_INV(:)
     ZP_SNOWSFCH(:)     = ZP_SNOWSFCH(:)     /ZP_PSN_INV(:)
     ZP_RESTOREN(:)     = ZP_RESTOREN(:)     /ZP_PSN_INV(:)
!
     ZP_SRSNOW(:)       = ZP_SRSNOW(:)       /ZP_PSN_INV(:)
     ZP_RRSNOW(:)       = ZP_RRSNOW(:)       /ZP_PSN_INV(:)
     ZP_SNDRIFT(:)      = ZP_SNDRIFT(:)      /ZP_PSN_INV(:)
!
     DO JJ=1,KSIZE2
        DO JI=1,KSIZE1
           ZP_SNOWSWE (JI,JJ) = ZP_SNOWSWE (JI,JJ) /ZP_PSN_INV(JI)
           ZP_SNOWHEAT(JI,JJ) = ZP_SNOWHEAT(JI,JJ)/ZP_PSN_INV(JI)
           ZP_SNOWDZ  (JI,JJ) = ZP_SNOWDZ  (JI,JJ)  /ZP_PSN_INV(JI)
        ENDDO
     ENDDO
!     
     ZP_SNOWHMASS(:)  = ZP_SNOWHMASS(:)/ZP_PSN_INV(:)
     ZP_THRUFAL(:)    = ZP_THRUFAL(:)  /ZP_PSN_INV(:)
     ZP_MELTSTOT(:)   = ZP_MELTSTOT(:) /ZP_PSN_INV(:)
     ZP_SNREFREEZ(:)  = ZP_SNREFREEZ(:)/ZP_PSN_INV(:)
!
!    Final Adjustments:
!    ------------------
!    Add cooling/heating flux correction to underlying soil.
!    This term is usually active for vanishingly thin snowpacks..
!    it is put outside of the snow scheme owing to it's dependence on
!    snow fraction. It is related to a possible correction to the ground-snow
!    heat flux when it is imposed (using MEB).
!    Also, it is added as a heat sink/source here since
!    fluxes have already be computed and should not be adjusted at this point:
!    applying it to the soil has the same impact as soil freeze-thaw, in the
!    sense it is computed after the fluxes have been updated.
!    (and update heat storage diagnostic in a consistent manner)
!
!    Energy is thickness weighted, thus thicker layers receive more energy and energy
!    is evenly distributed to depth ZDEPTHABS. An
!    alternate method is to weight near surface layers more and diminish weights
!    (thus eenrgy received by each layer) with depth. Both methods conserve energy as
!    long as vertical weights are normalized.

!    i) Determine soil depth for energy absorption:

     ZP_SOILD(:) = ZP_DZG(:,1)
     DO JJ=2,KSIZE3
        DO JI=1,KSIZE1
           IF(ZP_DZG(JI,JJ) <= ZDEPTHABS)THEN
              ZP_SOILD(JI) = ZP_DZG(JI,JJ)
           ENDIF
        ENDDO
     ENDDO

!    ii) Distribute (possible) energy to absorb vertically over some layer (defined above):

     ZP_PSN_GFLXCOR(:)  = ZP_PSN(:)*ZP_GFLXCOR(:)                                ! (W/m2)
     ZP_WORK(:)         = ZP_PSN_GFLXCOR(:)*PTSTEP/ZP_SOILD(:)

     ZP_TG(:,1)         = ZP_TG(:,1)         + ZP_WORK(:)*ZP_CT(:)*ZP_D_G(:,1)   ! (K)
     DO JJ=2,KSIZE3
        DO JI=1,KSIZE1
           IF (ZP_SOILD(JI) <= ZDEPTHABS) THEN
              ZP_TG(JI,JJ) = ZP_TG(JI,JJ)    + ZP_WORK(JI)/ZP_SOILHCAPZ(JI,JJ)   ! K
           ENDIF
        ENDDO
     ENDDO

     ZP_DELHEATG(:)     = ZP_DELHEATG(:)     + ZP_PSN_GFLXCOR(:)                 ! (W/m2)
     ZP_DELHEATG_SFC(:) = ZP_DELHEATG_SFC(:) + ZP_PSN_GFLXCOR(:)                 ! (W/m2)
!
     ZP_FLSN_COR(:)     = 0.0
!
  ELSE
!
!    To conserve energy in ISBA, the correction flux must be distributed at least
!    over the first 60cm depth. This method prevent numerical oscillations
!    especially when explicit snow vanishes. Final Adjustments are done in ISBA_CEB
!
     ZP_DELHEATG    (:) = 0.0
     ZP_DELHEATG_SFC(:) = 0.0
     ZP_FLSN_COR    (:) = ZP_GFLXCOR(:) ! (W/m2)
!
  ENDIF
!
!
!===============================================================
!conversion of snow heat from J/m2 into J/m3
WHERE(ZP_SNOWSWE (:,:)>0.)
      ZP_SNOWHEAT(:,:)=ZP_SNOWHEAT(:,:)*ZP_SNOWRHO(:,:)/ZP_SNOWSWE(:,:)
ENDWHERE
!===============================================================
!
! === Packing:
!
! unpack variables
!
DO JWRK=1,KSIZE2
  DO JJ=1,KSIZE1
    JI = KMASK(JJ)
    PEK%TSNOW%WSNOW(JI,JWRK) = ZP_SNOWSWE  (JJ,JWRK)
    PEK%TSNOW%RHO  (JI,JWRK) = ZP_SNOWRHO  (JJ,JWRK)
    PEK%TSNOW%HEAT (JI,JWRK) = ZP_SNOWHEAT (JJ,JWRK)
    PEK%TSNOW%AGE  (JI,JWRK) = ZP_SNOWAGE  (JJ,JWRK)
    DMK%XSNOWTEMP(JI,JWRK)   = ZP_SNOWTEMP (JJ,JWRK)
    DMK%XSNOWLIQ (JI,JWRK)   = ZP_SNOWLIQ  (JJ,JWRK)
    DMK%XSNOWDZ  (JI,JWRK)   = ZP_SNOWDZ   (JJ,JWRK)
  ENDDO
ENDDO
!
IF (GCOMPUTECRODIAG) THEN
  PEK%TSNOW%DEP_SUP    (:) = 0
  PEK%TSNOW%DEP_TOT    (:) = 0
  PEK%TSNOW%DEP_HUM    (:) = 0
  PEK%TSNOW%NAT_LEV    (:) = 6
  PEK%TSNOW%PRO_SUP_TYP(:) = 6
  PEK%TSNOW%AVA_TYP    (:) = 6

  DO JJ=1,KSIZE1
    JI = KMASK(JJ)
    PEK%TSNOW%DEP_SUP    (JI) = ZP_DEP_SUP    (JJ)
    PEK%TSNOW%DEP_TOT    (JI) = ZP_DEP_TOT    (JJ)
    PEK%TSNOW%DEP_HUM    (JI) = ZP_DEP_HUM    (JJ)
    PEK%TSNOW%NAT_LEV    (JI) = ZP_NAT_LEV    (JJ)
    PEK%TSNOW%PRO_SUP_TYP(JI) = ZP_PRO_SUP_TYP(JJ)
    PEK%TSNOW%AVA_TYP    (JI) = ZP_AVA_TYP    (JJ)
  ENDDO
ENDIF
!
IF (PEK%TSNOW%SCHEME=='CRO') THEN
!
  DO JWRK=1,KSIZE2
    DO JJ=1,KSIZE1
      JI = KMASK(JJ)
      PEK%TSNOW%DIAMOPT(JI,JWRK) = ZP_SNOWDIAMOPT(JJ,JWRK)
      PEK%TSNOW%SPHERI (JI,JWRK) = ZP_SNOWSPHERI (JJ,JWRK)
      PEK%TSNOW%HIST   (JI,JWRK) = ZP_SNOWHIST   (JJ,JWRK)
    ENDDO
  ENDDO
  !
  IF(SIZE(PBLOWSNW_FLUX,2) == 4)THEN
    DO JWRK=1,KSIZE4
      DO JJ=1,KSIZE1
        JI = KMASK(JJ)
        PBLOWSNW_FLUX(JI,JWRK) = ZP_BLOWSNW(JJ,JWRK)
      ENDDO
    ENDDO
  ENDIF
  !
  DO JIMP=1,NIMPUR
    DO JWRK=1,KSIZE2
      DO JJ=1,KSIZE1
        JI = KMASK(JJ)
        PEK%TSNOW%IMPUR(JI,JWRK,JIMP) = ZP_SNOWIMPUR(JJ,JWRK,JIMP)
      ENDDO
    ENDDO
  ENDDO

!
IF (PEK%TSNOW%SCHEME=='CRO' .AND. GCOMPUTECRODIAG) THEN
   DO JIMP=1,NIMPUR
     DO JWRK=1,KSIZE2
       DO JJ=1,KSIZE1
         JI = KMASK(JJ)
         DMK%XIMPUR_CONC (JI,JWRK,JIMP) = ZP_SNOWIMP_CONC (JJ,JWRK,JIMP)
       ENDDO
     ENDDO
   ENDDO
   DO JWRK=1,SIZE(P_DIR_SW,2)
      DO JJ=1,KSIZE1
         JI = KMASK(JJ)
         DMK%XDIFF_RATIO(JI,JWRK) = ZP_DIFF_RATIO(JJ,JWRK)
         DMK%XSPEC_ALB  (JI,JWRK) = ZP_SPEC_ALB  (JJ,JWRK) 
         DMK%XSPEC_TOT  (JI,JWRK) = ZP_SPEC_TOT  (JJ,JWRK)
      ENDDO
   ENDDO
   IF ((SIZE(DMK%XSPECMOD)>1) .AND. (SIZE(DMK%XSPEC_ALB,2)>1))THEN
     DO JWRK=1, NPNBANDS_MODIS
       DO JJ=1, KSIZE1
         JI = KMASK(JJ)
         DMK%XSPECMOD(JI, JWRK) = DMK%XSPEC_ALB(JI, XPWAVEIND_MODIS(JWRK))
       ENDDO
     ENDDO
   ENDIF
ENDIF
!
  IF (GCOMPUTECRODIAG)THEN
  ! This is equivalent to test the value of DGMI%LPROSNOW which does not enter in ISBA
    DO JWRK = 1,KSIZE2
      DO JJ=1,KSIZE1
        JI = KMASK(JJ)
        DMK%XSNOWDEND     (JI,JWRK) = ZP_SNOWDEND     (JJ,JWRK)
        DMK%XSNOWSPHER    (JI,JWRK) = ZP_SNOWSPHER    (JJ,JWRK)
        DMK%XSNOWSIZE     (JI,JWRK) = ZP_SNOWSIZE     (JJ,JWRK)
        DMK%XSNOWSSA      (JI,JWRK) = ZP_SNOWSSA      (JJ,JWRK)
        DMK%XSNOWHIST     (JI,JWRK) = ZP_SNOWHIST     (JJ,JWRK)
        DMK%XSNOWTYPEMEPRA(JI,JWRK) = ZP_SNOWTYPEMEPRA(JJ,JWRK)
        DMK%XSNOWRAM      (JI,JWRK) = ZP_SNOWRAM      (JJ,JWRK)
        DMK%XSNOWSHEAR    (JI,JWRK) = ZP_SNOWSHEAR    (JJ,JWRK)
        DMK%XACC_RAT      (JI,JWRK) = ZP_ACC_RAT      (JJ,JWRK)
        DMK%XNAT_RAT      (JI,JWRK) = ZP_NAT_RAT      (JJ,JWRK)
      ENDDO
    ENDDO
  ENDIF
ENDIF
!
DO JWRK=1,KSIZE3
   DO JJ=1,KSIZE1
      JI              = KMASK          (JJ)
      PTG    (JI,JWRK)= ZP_TG        (JJ,JWRK)
   ENDDO
ENDDO
!
DO JJ=1,KSIZE1
  !
  JI                 = KMASK          (JJ)
  !
  PEK%TSNOW%ALB (JI) = ZP_SNOWALB     (JJ)
  PEK%TSNOW%EMIS(JI) = ZP_EMISNOW     (JJ)
  !
  DMK%XCDSNOW   (JI) = ZP_CDSNOW      (JJ)
  DMK%XUSTARSNOW(JI) = ZP_USTARSNOW   (JJ)
  DMK%XCHSNOW   (JI) = ZP_CHSNOW      (JJ)
  DMK%XSNOWHMASS(JI) = ZP_SNOWHMASS   (JJ)
  DMK%XRNSNOW   (JI) = ZP_RNSNOW      (JJ)
  DMK%XHSNOW    (JI) = ZP_HSNOW       (JJ)
  DMK%XHPSNOW   (JI) = ZP_HPSNOW      (JJ)
  DMK%XGFLUXSNOW(JI) = ZP_GFLUXSNOW   (JJ)
  !
  DEK%XSNDRIFT      (JI)   = ZP_SNDRIFT      (JJ)
  !
  DEK%XDELHEATG     (JI)   = ZP_DELHEATG     (JJ)
  DEK%XDELHEATG_SFC (JI)   = ZP_DELHEATG_SFC (JJ)
  DEK%XMELTSTOT     (JI)   = ZP_MELTSTOT     (JJ)
  DEK%XSNREFREEZ    (JI)   = ZP_SNREFREEZ    (JJ)
  DEK%XDELHEATN     (JI)   = ZP_DELHEATN     (JJ)
  DEK%XDELHEATN_SFC (JI)   = ZP_DELHEATN_SFC (JJ)
  DEK%XDELPHASEN    (JI)   = ZP_DELPHASEN    (JJ)
  DEK%XDELPHASEN_SFC(JI)   = ZP_DELPHASEN_SFC(JJ)
  DEK%XRESTOREN     (JI)   = ZP_RESTOREN     (JJ)
  !
  PTHRUFAL     (JI)   = ZP_THRUFAL     (JJ)
  PEVAPCOR     (JI)   = ZP_EVAPCOR     (JJ)
  PRI          (JI)   = ZP_RI          (JJ)
  PQS          (JI)   = ZP_QS          (JJ)
  PGRNDFLUX    (JI)   = ZP_GRNDFLUX    (JJ)
  PFLSN_COR    (JI)   = ZP_FLSN_COR    (JJ)
  PSNOWSFCH    (JI)   = ZP_SNOWSFCH    (JJ)
  PLES3L       (JI)   = ZP_LES3L       (JJ)
  PLEL3L       (JI)   = ZP_LEL3L       (JJ)
  PEVAP        (JI)   = ZP_EVAP        (JJ)
  ZSOILCOR     (JI)   = ZP_SOILCOR     (JJ)
  !
  ZSWNET_N      (JI) = ZP_SWNETSNOW   (JJ)
  ZSWNET_NS     (JI) = ZP_SWNETSNOWS  (JJ)
  ZLWNET_N      (JI) = ZP_LWNETSNOW   (JJ)
  !
ENDDO
!
IF (GCOMPUTECRODIAG)THEN
  ! This is equivalent to test the value of DGMI%LPROSNOW which does not enter in ISBATHEN
  DMK%XSNDPT_12H     (:) = XUNDEF
  DMK%XSNDPT_1DY     (:) = XUNDEF
  DMK%XSNDPT_3DY     (:) = XUNDEF
  DMK%XSNDPT_5DY     (:) = XUNDEF
  DMK%XSNDPT_7DY     (:) = XUNDEF
  DMK%XSNSWE_1DY     (:) = XUNDEF
  DMK%XSNSWE_3DY     (:) = XUNDEF
  DMK%XSNSWE_5DY     (:) = XUNDEF
  DMK%XSNSWE_7DY     (:) = XUNDEF
  DMK%XSNRAM_SONDE   (:) = XUNDEF
  DMK%XSN_WETTHCKN   (:) = XUNDEF
  DMK%XSN_REFRZNTHCKN(:) = XUNDEF
  DMK%XDEP_HIG       (:) = XUNDEF
  DMK%XDEP_MOD       (:) = XUNDEF
  DMK%XACC_LEV       (:) = 4
  DMK%XPRO_INF_TYP   (:) = 6
  DO JJ=1,KSIZE1
    JI = KMASK(JJ)
    DMK%XSNDPT_12H     (JI) = ZP_SNDPT_12H     (JJ)
    DMK%XSNDPT_1DY     (JI) = ZP_SNDPT_1DY     (JJ)
    DMK%XSNDPT_3DY     (JI) = ZP_SNDPT_3DY     (JJ)
    DMK%XSNDPT_5DY     (JI) = ZP_SNDPT_5DY     (JJ)
    DMK%XSNDPT_7DY     (JI) = ZP_SNDPT_7DY     (JJ)
    DMK%XSNSWE_1DY     (JI) = ZP_SNSWE_1DY     (JJ)
    DMK%XSNSWE_3DY     (JI) = ZP_SNSWE_3DY     (JJ)
    DMK%XSNSWE_5DY     (JI) = ZP_SNSWE_5DY     (JJ)
    DMK%XSNSWE_7DY     (JI) = ZP_SNSWE_7DY     (JJ)
    DMK%XSNRAM_SONDE   (JI) = ZP_SNRAM_SONDE   (JJ)
    DMK%XSN_WETTHCKN   (JI) = ZP_SN_WETTHCKN   (JJ)
    DMK%XSN_REFRZNTHCKN(JI) = ZP_SN_REFRZNTHCKN(JJ)
    DMK%XDEP_HIG       (JI) = ZP_DEP_HIG       (JJ)
    DMK%XDEP_MOD       (JI) = ZP_DEP_MOD       (JJ)
    DMK%XACC_LEV       (JI) = ZP_ACC_LEV       (JJ)
    DMK%XPRO_INF_TYP   (JI) = ZP_PRO_INF_TYP   (JJ)
  ENDDO
ENDIF
!
IF (LHOOK) CALL DR_HOOK('SNOW3L_ISBA:CALL_MODEL',1,ZHOOK_HANDLE)
!
END SUBROUTINE CALL_MODEL
!
END SUBROUTINE SNOW3L_ISBA
