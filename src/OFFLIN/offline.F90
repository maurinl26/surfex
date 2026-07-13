 !SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!SFX_LIC for details. version 1.
! -------------------------------------------------
PROGRAM OFFLINE
!
! -------------------------------------------------
! Driver structure
! ----------------
! 1. Initializations
! 2. Temporal loops
!   2.a Read forcing
!   2.b Interpolate forcing in time
!   2.c Run surface
!   2.d Write prognostics and diagnostics variables
!
! modifications
! 09/2012 G. Pigeon: coherence between radiation and zenith angle because of
!                    trouble with radiation received by wall in TEB
! 03/2014 E. Martin change indices names in OMP module according to GMAP changes
! 05/2014 B. Decharme delete trip; I5 format to print DAY
! 04/2013 P. Le Moigne Add XDELTA_OROG to fix the maximum difference allowed between
!                      forcing and surface file orographies if LSET_FORC_ZS=.F
! 12/2013 S.Senesi     Add call to Gelato diag files init and close
! 02/2016: replace DOUBLE PRECISION by REAL to handle problem for promotion of real with GMKPACK or IBM SP
! 06/2016 S.Senesi     Use XIOS for diags output
! 08/2016 R. Seferian  New implentation of landuse change
! 03/2020 B. Decharme  New forcing interpolation
! 05/2023 P. Le Moigne temperature profile in sediments
! 11/2023 P. Le Moigne effective rain P-E 
! -------------------------------------------------
!
USE MODD_OFF_SURFEX_n
!
USE MODD_TYPE_DATE_SURF, ONLY : DATE
!
USE MODD_WRITE_SURF_ATM, ONLY : LFIRST_WRITE, NCPT_WRITE
!
USE MODD_FORC_ATM,  ONLY: CSV         ,&! name of all scalar variables
                            XDIR_ALB    ,&! direct albedo for each band
                            XSCA_ALB    ,&! diffuse albedo for each band
                            XEMIS       ,&! emissivity
                            XTSRAD      ,&! radiative temperature
                            XTSUN       ,&! solar time                    (s from midnight)
                            XZS         ,&! orography                             (m)
                            XZREF       ,&! height of T,q forcing                 (m)
                            XUREF       ,&! height of wind forcing                (m)
                            XTA         ,&! air temperature forcing               (K)
                            XQA         ,&! air humidity forcing                  (kg/kg)
                            XSV         ,&! scalar variables
                            XU          ,&! zonal wind                            (m/s)
                            XV          ,&! meridian wind                         (m/s)
                            XDIR_SW     ,&! direct  solar radiation (on horizontal surf.)
                            XSCA_SW     ,&! diffuse solar radiation (on horizontal surf.)
                            XSW_BANDS   ,&! mean wavelength of each shortwave band (m)
                            XZENITH     ,&! zenithal angle       (radian from the vertical)
                            XZENITH2    ,&! zenithal angle       (radian from the vertical)
                            XAZIM       ,&! azimuthal angle      (radian from North, clockwise)
                            XLW         ,&! longwave radiation (on horizontal surf.)
                            XPS         ,&! pressure at atmospheric model surface (Pa)
                            XPA         ,&! pressure at forcing level             (Pa)
                            XRHOA       ,&! density at forcing level              (kg/m3)
                            XCO2        ,&! CO2 concentration in the air          (kg/m3)
                            XIMPWET     ,&! wet deposit coef for each type of impurity (g/m2/s)
                            XIMPDRY     ,&! dry deposit coef for each type of impurity (g/m2/s)
                            XO3         ,&! Ozone
                            XAE         ,&! Aerosol optical depth
                            XSNOW       ,&! snow precipitation                    (kg/m2/s)
                            XRAIN       ,&! liquid precipitation                  (kg/m2/s)
                            XSFTH       ,&! flux of heat                          (W/m2)
                            XSFTQ       ,&! flux of water vapor                   (kg/m2/s)
                            XSFU        ,&! zonal momentum flux                   (m/s)
                            XSFV        ,&! meridian momentum flux                (m/s)
                            XSFCO2      ,&! flux of CO2                           (kg/m2/s)
                            XSFTS       ,&! flux of scalar var.                   (kg/m2/s)
                            XPEW_A_COEF ,&! implicit coefficients
                            XPEW_B_COEF ,&! needed if HCOUPLING='I'
                            XPET_A_COEF ,&
                            XPEQ_A_COEF ,&
                            XPET_B_COEF ,&
                            XPEQ_B_COEF ,&
                            XTSURF      ,&! effective temperature                  (K)
                            XZ0         ,&! surface roughness length for momentum  (m)
                            XZ0H        ,&! surface roughness length for heat      (m)
                            XQSURF      ,&! specific humidity at surface           (kg/kg)
                            XZWS          ! significant wave height                (m)
!
USE MODD_SURF_PAR,   ONLY : NUNDEF
!
USE MODD_SURF_CONF,  ONLY : CPROGNAME, CSOFTWARE
USE MODD_CSTS,       ONLY : XPI, XDAY, XRV, XRD, XG
USE MODD_SYTRON_PAR, ONLY : NTAB_SYT
USE MODD_IO_SURF_ASC,ONLY : CFILEIN,CFILEIN_SAVE,CFILEOUT,CFILEPGD
USE MODD_IO_SURF_FA, ONLY : CFILEIN_FA, CFILEIN_FA_SAVE,       &
                            CFILEOUT_FA, NUNIT_FA, CDNOMC,     &
                            IVERBFA, LFANOCOMPACT, CFILEPGD_FA
USE MODD_IO_SURF_LFI,ONLY : CFILEIN_LFI, CFILEIN_LFI_SAVE, CLUOUT_LFI, CFILEOUT_LFI, &
                            LMNH_COMPATIBLE, CFILEPGD_LFI
USE MODD_IO_SURF_NC, ONLY : CFILEIN_NC, CFILEIN_NC_SAVE, CFILEOUT_NC, CLUOUT_NC, &
                            CFILEPGD_NC, LDEF_nc=>LDEF, LRESET_DIAG_nc=>LRESET_DIAG
USE MODI_ENDDEF_IO_SURF_NC_n
USE MODI_CLOSE_IO_SURF_NC_n
USE MODD_IO_SURF_OL, ONLY : XSTART, XCOUNT, XSTRIDE, LPARTW,   &
                            XSTARTW, XCOUNTW, LTIME_WRITTEN,   &
                            NSTEP_OUTPUT, NEND_ATM, XTIMEC,    &
                            LDEF_ol=>LDEF, LRESET_DIAG_ol=>LRESET_DIAG
USE MODD_WRITE_BIN,  ONLY : NWRITE
!
USE MODD_SURFEX_MPI, ONLY : NCOMM, NPROC, NRANK, NPIO, WLOG_MPI, PREP_LOG_MPI,   &
                            NINDEX, NSIZE_TASK, XTIME_NPIO_READ, &
                            XTIME_COMM_READ, XTIME_WRITE, XTIME_CALC, IDX_W, END_LOG_MPI
!
USE MODD_SURFEX_OMP, ONLY :  NBLOCK, NBLOCKTOT
!
USE MODD_COUPLING_TOPD, ONLY : NNB_TOPD, NNB_STP_RESTART, LBUDGET_TOPD, LTOPD_STEP, &
                               LCOUPL_TOPD, NTOPD_STEP, NYEAR, NMONTH, NDAY, NH, NM, &
                               LSUBCAT
USE MODD_TOPODYN, ONLY : XTOPD_STEP, NNB_TOPD_STEP, XQTOT, XQB_RUN, XQB_DR
!
USE MODD_SLOPE_EFFECT, ONLY: XZS_THREAD,XZS_XY_THREAD,XSLOPANG_THREAD,&
                             XSLOPAZI_THREAD,XSURF_TRIANGLE_THREAD
!
USE MODD_SFX_OASIS, ONLY : LOASIS, XRUNTIME
!
USE MODD_XIOS, ONLY : LXIOS, TXIOS_CONTEXT, LXIOS_DEF_CLOSED, LADD_DIM=>LALLOW_ADD_DIM, NTIMESTEP
!
USE MODD_SURF_ATM_TURB_n, ONLY : SURF_ATM_TURB_t
!
USE MODE_POS_SURF
!
USE MODE_CRODEBUG
!
USE MODE_DATES_NETCDF
!
USE MODN_IO_OFFLINE
USE MODD_ASSIM, ONLY : LCROCO
!
USE MODI_GET_LUOUT
USE MODI_OPEN_NAMELIST
USE MODI_TEST_NAM_VAR_SURF
USE MODI_CLOSE_NAMELIST
USE MODI_READ_ALL_NAMELISTS
USE MODI_OPEN_CLOSE_BIN_ASC_FORC
USE MODI_OPEN_FILEIN_OL
USE MODI_OL_READ_ATM_CONF
USE MODI_ABOR1_SFX
USE MODI_OL_ALLOC_ATM
USE MODI_COMPARE_OROGRAPHY
USE MODI_SUNPOS
USE MODI_INIT_INDEX_MPI
USE MODI_OL_READ_ATM
USE MODI_IO_BUFF_CLEAN
USE MODI_INIT_SURF_ATM_n
USE MODI_OL_TIME_INTERP_ATM
USE MODI_OL_PRECIP_FRC_PDF
USE MODI_COUPLING_SURF_ATM_n
USE MODI_ADD_FORECAST_TO_DATE_SURF
USE MODI_WRITE_SURF_ATM_n
USE MODI_WRITE_HEADER_MNH
USE MODI_FLAG_UPDATE
USE MODI_FLAG_DIAG_UPDATE
USE MODI_DIAG_SURF_ATM_n
USE MODI_WRITE_DIAG_SURF_ATM_n
USE MODI_GET_SURF_VAR_n
USE MODI_GATHER_AND_WRITE_MPI
USE MODI_CLOSE_FILEIN_OL
USE MODI_CLOSE_FILEOUT_OL
USE MODI_INIT_OUTPUT_OL_n
USE MODI_SFX_XIOS_INIT_OUTPUT_OL
USE MODI_INIT_OUTPUT_NC_n
!
USE MODI_WRITE_HEADER_FA
USE MODI_ABOR1_SFX
!
USE MODI_WRITE_DISCHARGE_FILE
USE MODI_WRITE_DISCHARGE_FILE_SUB
USE MODI_WRITE_BUDGET_COUPL_ROUT
USE MODI_PREP_RESTART_COUPL_TOPD
!
USE MODI_INIT_SLOPE_PARAM
USE MODI_SLOPE_RADIATIVE_EFFECT
!
USE MODI_SFX_OASIS_READ_NAM
USE MODI_SFX_OASIS_INIT
USE MODI_SFX_OASIS_DEF_OL
USE MODI_SFX_OASIS_RECV_OL
USE MODI_SFX_OASIS_SEND_OL
USE MODI_SFX_OASIS_END
!RJ: missing modi
USE MODI_LOCAL_SLOPE_PARAM
!
#ifdef WXIOS
USE XIOS, ONLY : XIOS_CONTEXT_FINALIZE, XIOS_CLOSE_CONTEXT_DEFINITION, XIOS_UPDATE_CALENDAR, &
                 XIOS_INITIALIZE
#endif
USE MODI_SFX_XIOS_READNAM_OL
USE MODI_SFX_XIOS_SETUP_OL
!
USE MODE_GLT_DIA_LU
!
USE MODI_INIT_SYTRON_TABLE
!
! spectral repartition of irradiance
USE MODD_CONST_ATM, ONLY : JPNBANDS_ATM
USE MODE_ATMO_TARTES, ONLY : RADIANCE
#ifdef SFX_MPI
#ifdef SFX_MPL
USE MPL_DATA_MODULE, ONLY : LMPLUSERCOMM, MPLUSERCOMM
#endif
#endif
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
#ifdef AIX64
!$ USE OMP_LIB
#endif
!
IMPLICIT NONE
!
#ifdef SFX_MPI
INCLUDE 'mpif.h'
#endif
!
#ifndef AIX64
!$ INCLUDE 'omp_lib.h'
#endif
!
!*      0.    declarations of local variables
!
 CHARACTER(LEN=3), PARAMETER       :: YINIT     = 'ALL'
!
 CHARACTER(LEN=28)                 :: YLUOUT    = 'LISTING_OFFLINE             '
!
INTEGER                           :: IYEAR               ! current year (UTC)
INTEGER                           :: IMONTH              ! current month (UTC)
INTEGER                           :: IDAY                ! current day (UTC)
!

! ###########################################################################
INTEGER                           :: IYEAR_NEXT          ! current year (UTC)
INTEGER                           :: IMONTH_NEXT         ! current month (UTC)
INTEGER                           :: IDAY_NEXT           ! current day (UTC)
!
INTEGER                           :: IYEAR_PREV          ! current year (UTC)
INTEGER                           :: IMONTH_PREV         ! current month (UTC)
INTEGER                           :: IDAY_PREV           ! current day (UTC)
! ###########################################################################


INTEGER                           :: IYEAR2              ! current year at end of timestep(UTC)
INTEGER                           :: IMONTH2             ! current month at end of timestep(UTC)
INTEGER                           :: IDAY2               ! current day at end of timestep(UTC)
REAL                              :: ZTIME               ! current time since start of the day (s)
REAL                              :: ZTIME2              ! current time since start of the day at end of timestep (s)
REAL                              :: ZTIMEC              ! current duration since start of the run (s)



! ###########################################################################
REAL                              :: ZTIME_NEXT          ! current time since start of the day (s)
REAL                              :: ZTIME_PREV          ! current time since start of the day (s)
! ###########################################################################



!
INTEGER                           :: IYEAR_OUT           ! output year name
INTEGER                           :: IMONTH_OUT          ! output month name
INTEGER                           :: IDAY_OUT            ! output day name
REAL                              :: ZTIME_OUT           ! output time since start of the run (s)
!
INTEGER, DIMENSION(11)  :: IDATEF
!
 CHARACTER(LEN=28), PARAMETER      :: YATMFILE     = '                            '
 CHARACTER(LEN=6),  PARAMETER      :: YATMFILETYPE = '      '
 CHARACTER(LEN=2),  PARAMETER      :: YTEST        = 'OK'          ! must be equal to 'OK'
!
REAL, DIMENSION(:), POINTER       :: ZLAT                ! latitude                         (rad)
REAL, DIMENSION(:), POINTER       :: ZLON                ! longitude                        (rad)
REAL, DIMENSION(:), POINTER       :: ZZS_FORC            ! orography                        (m)
REAL, DIMENSION(:), POINTER       :: ZZREF               ! Forcing level for T
REAL, DIMENSION(:), POINTER       :: ZUREF               ! Forcing level for U
!
REAL                              :: ZTSTEP              ! atmospheric time-step            (s)
!
INTEGER                           :: INI                 ! grid dimension
INTEGER                           :: JLOOP               ! loop counter
INTEGER                           :: IBANDS              ! Number of radiative bands
INTEGER                           :: INB_STEP_ATM        ! Number of atmospheric time-steps
INTEGER                           :: INB_ATM             ! Number of Isba time-steps
                                                         ! within a forcing time-step
INTEGER                           :: ID_FORC             ! indice of forcing in the file
INTEGER                           :: INB_LINES           ! nb of lines to read in the forcing file
INTEGER                           :: IDMAX               ! nb of lines to read in the forcing file at last
INTEGER                           :: JFORC_STEP          ! atmospheric loop index
INTEGER                           :: IFORC_STEP          ! atmospheric count index
INTEGER                           :: JSURF_STEP          ! isba loop index
INTEGER                           :: ICOUNT              ! day counter
LOGICAL                           :: LLAST_TIMESTEP      ! .True. for Last timestep
INTEGER                           :: ITIMESTARTINDEX
REAL                              :: ZDURATION           ! duration of run                       (s)
!
REAL, DIMENSION(:,:), ALLOCATABLE :: ZTA                 ! air temperature forcing               (K)
REAL, DIMENSION(:,:), ALLOCATABLE :: ZQA                 ! air humidity forcing                  (kg/m3)
REAL, DIMENSION(:,:), ALLOCATABLE :: ZWIND               ! wind speed                            (m/s)
REAL, DIMENSION(:,:), ALLOCATABLE :: ZSCA_SW             ! diffuse solar radiation (on horizontal surf.)
REAL, DIMENSION(:,:), ALLOCATABLE :: ZDIR_SW             ! direct  solar radiation (on horizontal surf.)
REAL, DIMENSION(:,:), ALLOCATABLE :: ZLW                 ! longwave radiation (on horizontal surf.)
REAL, DIMENSION(:,:), ALLOCATABLE :: ZSNOW               ! snow precipitation                    (kg/m2/s)
REAL, DIMENSION(:,:), ALLOCATABLE :: ZRAIN               ! liquid precipitation                  (kg/m2/s)
REAL, DIMENSION(:,:), ALLOCATABLE :: ZPS                 ! pressure at forcing level             (Pa)
REAL, DIMENSION(:,:), ALLOCATABLE :: ZCO2                ! CO2 concentration in the air          (kg/m3)
REAL, DIMENSION(:,:), ALLOCATABLE :: ZO3                 ! Ozone
REAL, DIMENSION(:,:), ALLOCATABLE :: ZAE                 ! Aerosol optical depth
REAL, DIMENSION(:,:,:), ALLOCATABLE :: ZIMPWET           ! wet deposit coefficient for each impurity type (g/m²/s)
REAL, DIMENSION(:,:,:), ALLOCATABLE :: ZIMPDRY           ! dry deposit coefficient for each impurity type (g/m²/s)
REAL, DIMENSION(:,:), ALLOCATABLE :: ZDIR                ! wind direction
INTEGER                           :: ILUOUT              ! ascii output unit number
INTEGER                           :: ILUNAM              ! namelist unit number
INTEGER                           :: IRET                ! error return code
INTEGER                           :: INB
INTEGER                           :: INW, JNW
CHARACTER(LEN=20)                 :: YTAG
LOGICAL                           :: GFOUND              ! return logical when reading namelist
LOGICAL                           :: GSHADOWS
REAL, DIMENSION(:),   ALLOCATABLE :: ZSW                 ! total solar radiation (on horizontal surf.)
REAL, DIMENSION(:),   ALLOCATABLE :: ZCOEF               ! coefficient for solar radiation interpolation near sunset/sunrise
REAL, DIMENSION(:),   ALLOCATABLE :: ZPDISTRIB           ! PDF for precipitation time interpolation
!
TYPE(DATE_TIME) :: TDATE_END
! Flag diag :
!
INTEGER                           :: I2M, IBEQ, IDSTEQ
LOGICAL                           :: GFRAC, GDIAG_GRID, GSURF_BUDGET, GRAD_BUDGET, GCOEF,    &
                                     GSURF_VARS, GDIAG_OCEAN, GDIAG_SEAICE, GWATER_PROFILE,  &
                                     GSEDIM_PROFILE, GSURF_EVAP_BUDGET, GFLOOD, GPGD_ISBA,   &
                                     GCH_NO_FLUX_ISBA, GSURF_MISC_BUDGET_ISBA, GPGD_TEB,     &
                                     GSURF_MISC_BUDGET_TEB, GLUTILES_BUDGET, GDIAG_RESTART,  &
                                     GFLKFLUX,GFLKWATER, GCROCO, GINTERPOL_TS
!
! Inquiry mode arrays:
!
REAL, DIMENSION(:), ALLOCATABLE   :: ZSEA, ZWATER, ZNATURE, ZTOWN
REAL, DIMENSION(:), ALLOCATABLE   :: ZSEA_FULL, ZWATER_FULL, ZNATURE_FULL, ZTOWN_FULL
REAL, DIMENSION(:), ALLOCATABLE   :: ZT2M, ZQ2M
REAL, DIMENSION(:), ALLOCATABLE   :: ZZ0, ZZ0H, ZQS
REAL, DIMENSION(:), ALLOCATABLE   :: ZQS_SEA, ZQS_WATER, ZQS_NATURE, ZQS_TOWN
REAL, DIMENSION(:), ALLOCATABLE   :: ZPSNG, ZPSNV
REAL, DIMENSION(:), ALLOCATABLE   :: ZZ0EFF
REAL, DIMENSION(:), ALLOCATABLE   :: ZZS
REAL, DIMENSION(:), ALLOCATABLE   :: ZZ0_FULL, ZZ0EFF_FULL, ZZS_FULL
REAL, DIMENSION(:), ALLOCATABLE   :: ZSUMZEN




! ###########################################################################
! Local variables for radiation at previous and next forcing
REAL, DIMENSION(:), ALLOCATABLE :: ZTSUN_PREV_FORC
REAL, DIMENSION(:), ALLOCATABLE :: ZZENITH_PREV_FORC
REAL, DIMENSION(:), ALLOCATABLE :: ZAZIM_PREV_FORC
!
REAL, DIMENSION(:), ALLOCATABLE :: ZTSUN_NEXT_FORC
REAL, DIMENSION(:), ALLOCATABLE :: ZZENITH_NEXT_FORC
REAL, DIMENSION(:), ALLOCATABLE :: ZAZIM_NEXT_FORC
! ###########################################################################



INTEGER :: ISERIES, ISIZE
!
! MPI variables
!
CHARACTER(LEN=100) :: YNAME
CHARACTER(LEN=10)  :: YRANK
INTEGER :: ILEVEL, INFOMPI, INKPROMA, JBLOCK,JIMP
INTEGER, DIMENSION(:), ALLOCATABLE :: ISIZE_OMP
DOUBLE PRECISION :: XTIME0, XTIME1, XTIME
!
! SFX - OASIS coupling variables
!
LOGICAL :: GSAVHOOK
INTEGER :: IBLOCKTOT, IBLOCK
!
TYPE(SURF_ATM_TURB_t) :: AT         ! atmospheric turbulence parameters
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
! --------------------------------------------------------------------------------------
!
!*     0.1.   MPI, OASIS, XIOS and dr_hook initializations
!
CSOFTWARE='OFFLINE'
!
INFOMPI=1
!
#ifdef WXIOS
CALL SFX_XIOS_READNAM_OL(CNAMELIST)
#else
LXIOS=.FALSE.
#endif
!
#ifdef CPLOASIS
!Must be call before DRHOOK !
CALL SFX_OASIS_INIT(CNAMELIST,NCOMM)
#else
LOASIS   = .FALSE.
XRUNTIME = 0.0
#ifdef WXIOS
IF (LXIOS) THEN
  CALL XIOS_INITIALIZE('surfex',return_comm=NCOMM)
ENDIF
#endif
#endif
!
!disabling DR_HOOK should not be necessary if a newer version of DR_HOOK is used
IF (LOASIS.OR.LXIOS) LHOOK=.FALSE.
!
#ifdef SFX_MPI
#ifdef SFX_MPL
IF (LOASIS.OR.LXIOS) THEN
  LMPLUSERCOMM = .TRUE.
  MPLUSERCOMM = NCOMM
ENDIF
#endif
#endif
!
#ifdef SFX_MPI
IF(.NOT.LOASIS.AND..NOT.LXIOS)THEN
 CALL MPI_INIT_THREAD(MPI_THREAD_MULTIPLE,ILEVEL,INFOMPI)
 IF (INFOMPI /= MPI_SUCCESS) THEN
    CALL ABOR1_SFX('OFFLINE: ERROR WHEN INITIALIZING MPI')
 ENDIF
 NCOMM=MPI_COMM_WORLD
ENDIF
CALL MPI_COMM_SIZE(NCOMM,NPROC,INFOMPI)
CALL MPI_COMM_RANK(NCOMM,NRANK,INFOMPI)
#endif
!
IF (LHOOK) CALL DR_HOOK('OFFLINE',0,ZHOOK_HANDLE)
!
!RJ: init modd_surefx_omp
!$OMP PARALLEL
!$ NBLOCKTOT = OMP_GET_NUM_THREADS()
!$ NBLOCK = OMP_GET_THREAD_NUM()
!$OMP END PARALLEL
!
IBLOCKTOT = 1
IBLOCK = 0
!
CALL PREP_LOG_MPI
!
CALL WLOG_MPI(' ')
!
CALL WLOG_MPI('NBLOCKTOT ',KLOG=NBLOCKTOT)
!
#ifdef SFX_MPI
XTIME0 = MPI_WTIME()
#endif
!
!
!*      0.3.   Open ascii file for writing
!
WRITE(YRANK,FMT='(I10)') NRANK
YNAME=TRIM(YLUOUT)//ADJUSTL(YRANK)
!
CLUOUT_LFI =  ADJUSTL(ADJUSTR(YNAME)//'.txt')
CLUOUT_NC  =  ADJUSTL(ADJUSTR(YNAME)//'.txt')
!
CALL GET_LUOUT('ASCII ',ILUOUT)
OPEN(UNIT=ILUOUT,FILE=ADJUSTL(ADJUSTR(YNAME)//'.txt'),FORM='FORMATTED',ACTION='WRITE')
!
!
IF ( NRANK==NPIO ) THEN
  !
!RJ: be verbose just for openmp
  IF(NBLOCKTOT==1) THEN
!$  WRITE(*,*) "CAUTION: DID YOU THINK TO SET OMP_NUM_THREADS=1?"
!$  WRITE(*,*) "PLEASE VERIFY OMP_NUM_THREADS IS INITIALIZED : TYPE ECHO $OMP_NUM_THREADS IN A TERMINAL"
  !
!$  WRITE(ILUOUT,*) "CAUTION: DID YOU THINK TO SET OMP_NUM_THREADS=1?"
!$  WRITE(ILUOUT,*) "PLEASE VERIFY OMP_NUM_THREADS IS INITIALIZED : TYPE ECHO $OMP_NUM_THREADS IN A TERMINAL"
  ENDIF
  !
ENDIF
!
!*      0.4.   Reads namelists
!
 CALL OPEN_NAMELIST('ASCII ',ILUNAM,CNAMELIST)
!
 CALL POSNAM(ILUNAM,'NAM_IO_OFFLINE',GFOUND,ILUOUT)
IF (GFOUND) READ (UNIT=ILUNAM,NML=NAM_IO_OFFLINE)
 CALL CLOSE_NAMELIST('ASCII ',ILUNAM)
!
IF (NPROC==1) THEN
  XIO_FRAC=1.
ELSE
  XIO_FRAC = MAX(MIN(XIO_FRAC,1.),0.)
ENDIF
!
CALL TEST_NAM_VAR_SURF(ILUOUT,'CSURF_FILETYPE',CSURF_FILETYPE,'ASCII ','LFI   ','FA    ','NC    ')
#ifdef WXIOS
 CALL TEST_NAM_VAR_SURF(ILUOUT,'CTIMESERIES_FILETYPE',CTIMESERIES_FILETYPE,'NETCDF','TEXTE ','BINARY',&
                                                                            'ASCII ','LFI   ','FA    ',&
                                                                            'NONE  ','OFFLIN','NC    '&
                                                                            ,'XIOS  ')
#else
 CALL TEST_NAM_VAR_SURF(ILUOUT,'CTIMESERIES_FILETYPE',CTIMESERIES_FILETYPE,'NETCDF','TEXTE ','BINARY',&
                                                                            'ASCII ','LFI   ','FA    ',&
                                                                            'NONE  ','OFFLIN','NC    ')  
#endif
CALL TEST_NAM_VAR_SURF(ILUOUT,'CFORCING_FILETYPE',CFORCING_FILETYPE,'NETCDF','ASCII ','BINARY')
CALL TEST_NAM_VAR_SURF(ILUOUT,'CTIME_INTERP_PRCP',CTIME_INTERP_PRCP,'DEF','OLD','PDF')
!
IF (NSCAL>59) CALL ABOR1_SFX("OFFLINE: NSCAL MUST BE LOWER THAN OR EQUAL TO 59")
!
!
IF (CTIMESERIES_FILETYPE=='NETCDF') CTIMESERIES_FILETYPE='OFFLIN'
!
IF ((TRIM(CTIMESERIES_FILETYPE) /= 'XIOS') .AND. LADD_DIM) THEN
      CALL ABOR1_SFX('CANNOT YET SET LALLOW_ADD_DIM TO .TRUE. WITHOUT SETTING CTIMESERIES_FILETYPE to XIOS ')
ENDIF
!
!
CFILEPGD = ADJUSTL(ADJUSTR(CPGDFILE)//'.txt')
CFILEIN  = ADJUSTL(ADJUSTR(CPREPFILE)//'.txt')
CFILEIN_SAVE = CFILEIN
!
CFILEPGD_LFI = CPGDFILE
CFILEIN_LFI  = CPREPFILE
CFILEIN_LFI_SAVE = CFILEIN_LFI
!
CFILEPGD_FA = ADJUSTL(ADJUSTR(CPGDFILE)//'.fa')
CFILEIN_FA  = ADJUSTL(ADJUSTR(CPREPFILE)//'.fa')
CFILEIN_FA_SAVE  = CFILEIN_FA
!
CFILEPGD_NC = ADJUSTL(ADJUSTR(CPGDFILE)//'.nc')
CFILEIN_NC  = ADJUSTL(ADJUSTR(CPREPFILE)//'.nc')
CFILEIN_NC_SAVE  = CFILEIN_NC
!
!     Allocations of Surfex Types
CALL SURFEX_ALLOC_LIST(1)
YSC => YSURF_LIST(1)
!
!     Reading all namelist (also assimilation)
CALL READ_ALL_NAMELISTS(YSC, CSURF_FILETYPE,'ALL',.FALSE.)
!
!
!*      0.5.   Reads SFX - OASIS coupling namelists
!
CALL SFX_OASIS_READ_NAM(CSURF_FILETYPE,XTSTEP_SURF)
!
!*      0.6   Assume FA filetype consistency
!
CPROGNAME = CSURF_FILETYPE
!
! --------------------------------------------------------------------------------------
!
!*      1.    Initializations
!
!       netcdf file handling
!
IF (NRANK==NPIO) THEN
  !
  XSTART            = NUNDEF
  XSTRIDE           = NUNDEF
  XCOUNT            = NUNDEF
  XSTARTW           = 0
  XCOUNTW           = 1
  LPARTW            = .TRUE.
  !
ENDIF
!
!New interpolation method
IF(LNEW_TIME_INTERP_ATM)THEN
  NEND_ATM=2
ELSE
  NEND_ATM=1
ENDIF
!
#ifdef SFX_MPI
XTIME = (MPI_WTIME() - XTIME0)
#endif
 CALL WLOG_MPI('READ NAMELISTS ',PLOG=XTIME)
#ifdef SFX_MPI
XTIME0 = MPI_WTIME()
#endif
!
!
!       splitting of the grid
!
GSHADOWS = LSHADOWS_SLOPE .OR. LSHADOWS_OTHER
CALL INIT_INDEX_MPI(YSC%DTCO, YSC%U, YSC%UG, YSC%GCP, CSURF_FILETYPE, 'OFF', YALG_MPI, XIO_FRAC, GSHADOWS)
!
 CALL WLOG_MPI(' ')
 CALL WLOG_MPI('TIME_NPIO_READ init_index ',PLOG=XTIME_NPIO_READ)
 CALL WLOG_MPI('TIME_COMM_READ init_index ',PLOG=XTIME_COMM_READ)
XTIME_NPIO_READ = 0.
XTIME_COMM_READ = 0.
!
#ifdef SFX_MPI
XTIME = (MPI_WTIME() - XTIME0)
#endif
 CALL WLOG_MPI(' ')
 CALL WLOG_MPI('INIT_INDEX_MPI ',PLOG=XTIME)
 CALL WLOG_MPI(' ')
#ifdef SFX_MPI
XTIME0 = MPI_WTIME()
#endif
!
!       forcing file handling
!
IF (CFORCING_FILETYPE=='ASCII ' .OR. CFORCING_FILETYPE=='BINARY') CALL OPEN_CLOSE_BIN_ASC_FORC('CONF ',CFORCING_FILETYPE,'R')
IF (CFORCING_FILETYPE=='NETCDF') CALL OPEN_FILEIN_OL
!
!       configuration of run
!
CALL OL_READ_ATM_CONF(YSC%DTCO, YSC%U, YSC%UG%G%CGRID, CSURF_FILETYPE, CFORCING_FILETYPE,  &
                      LDELAYEDSTART_NC, NDATESTOP, ZDURATION, ZTSTEP, INI,  &
                      IYEAR, IMONTH, IDAY, ZTIME, ZLAT, ZLON, ZZS_FORC,     &
                      ZZREF, ZUREF, ITIMESTARTINDEX     )
!
TDATE_END%TDATE%YEAR  = IYEAR
TDATE_END%TDATE%MONTH = IMONTH
TDATE_END%TDATE%DAY   = IDAY
TDATE_END%TIME        = ZTIME/3600
CALL ADDHOURS(TDATE_END,INT(ZDURATION/3600 ))
!
 CALL WLOG_MPI(' ')
 CALL WLOG_MPI('TIME_NPIO_READ forc conf ',PLOG=XTIME_NPIO_READ)
 CALL WLOG_MPI('TIME_COMM_READ forc conf ',PLOG=XTIME_COMM_READ)
XTIME_NPIO_READ = 0.
XTIME_COMM_READ = 0.
!
#ifdef SFX_MPI
XTIME = (MPI_WTIME() - XTIME0)
#endif
 CALL WLOG_MPI('OL_READ_ATM_CONF ',PLOG=XTIME)
 CALL WLOG_MPI(' ')
#ifdef SFX_MPI
XTIME0 = MPI_WTIME()
#endif
!
!*     time steps coherence check
!
IF ( (MOD(XTSTEP_OUTPUT,ZTSTEP)*MOD(ZTSTEP,XTSTEP_OUTPUT) /= 0) .OR. (MOD(ZTSTEP,XTSTEP_SURF) /= 0) ) THEN
   WRITE(ILUOUT,*)' FORCING  AND OUTPUT/SURFACE TIME STEP SHOULD BE MULTIPLE', &
     NINT(ZTSTEP),NINT(XTSTEP_OUTPUT),NINT(XTSTEP_SURF)
   CALL ABOR1_SFX('OFFLINE: FORCING  AND OUTPUT/SURFACE TIME STEP SHOULD BE MULTIPLE')
ENDIF
!
IF ( ZTIME /= 0. .AND. MOD(ZTIME,XTSTEP_SURF) /= 0  ) THEN
   WRITE(ILUOUT,*)' INITIAL AND SURFACE TIME STEP SHOULD BE MULTIPLE', &
   NINT(ZTIME),NINT(XTSTEP_SURF)
   CALL ABOR1_SFX('OFFLINE: INITIAL AND SURFACE TIME STEP SHOULD BE MULTIPLE')
ENDIF
!
IF(LOASIS.AND.ZDURATION/=XRUNTIME)THEN
   WRITE(ILUOUT,*)'Total simulated time given by Forcing field and OASIS namcouple are different'
   WRITE(ILUOUT,*)'From Forcing (s) : ',ZDURATION, 'From OASIS   (s) : ',XRUNTIME
   CALL ABOR1_SFX('OFFLINE: TOTAL SIMULATED TIME DIFFERENT BETWEEN FORCING AND OASIS')
ENDIF
!
INB_STEP_ATM  = INT(ZDURATION / ZTSTEP)
INB_ATM       = INT(ZTSTEP / XTSTEP_SURF)
NSTEP_OUTPUT  = INT(ZDURATION / XTSTEP_OUTPUT)
!
XTOPD_STEP = 0
NNB_TOPD_STEP = 0
NTOPD_STEP = 0
IF ( LCOUPL_TOPD ) THEN
  !
  XTOPD_STEP = FLOAT(NNB_TOPD)* XTSTEP_SURF
  NNB_TOPD_STEP = INT( ZDURATION / XTOPD_STEP )
  !
  IF ( NNB_STP_RESTART==0 .AND. .NOT.LRESTART ) NNB_STP_RESTART = -1
  !
  NTOPD_STEP = 1
  !
ENDIF
!
!       allocation of variables
!
IBANDS = 1
! special case for snowcro tartes!!!!!!

IF (LSPECSNOW) THEN
  IBANDS=JPNBANDS_ATM
ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
CALL OL_ALLOC_ATM(INI,IBANDS,NSCAL)
!XSW_BANDS=1.2E-6
!
XZS   = ZZS_FORC
XZREF = ZZREF
XUREF = ZUREF
!
!       compare orography
!
CALL COMPARE_OROGRAPHY(YSC%DTCO, YSC%U, CSURF_FILETYPE, LSET_FORC_ZS, XDELTA_OROG)
!
!       miscellaneous initialization
!
ICOUNT = 0
ZTIMEC = 0.
!
CALL SUNPOS(IYEAR, IMONTH, IDAY, ZTIME, ZLON, ZLAT, XTSUN, XZENITH, XAZIM)
!
!number of lines read in forcing files
INB_LINES=1
IF (NB_READ_FORC.EQ.1) THEN
  INB_LINES=INB_STEP_ATM
ELSEIF (NB_READ_FORC.NE.0) THEN
  !to be sure the number of readings will be NB_READ_FORC as a maximum
  INB_LINES=CEILING(1.*(INB_STEP_ATM+1)/NB_READ_FORC)
ENDIF
!number of lines to be read effectively
IDMAX=INB_LINES+NEND_ATM
!effective number of readings of the forcing files
NB_READ_FORC=CEILING(1.*(INB_STEP_ATM+1)/INB_LINES)
!
!Compute idealized PDF of precipitation over the forçing timestep
IF(LNEW_TIME_INTERP_ATM.AND.CTIME_INTERP_PRCP=='PDF')THEN
  ALLOCATE(ZPDISTRIB(INB_ATM))
  ZPDISTRIB(:) = 0.0
  CALL OL_PRECIP_FRC_PDF(INB_ATM,LPRINT,ZTSTEP,ZPDISTRIB)
ENDIF
!
!open Gelato specific diagnostic files (if requested by Gelato wizzard user)
!
#if ! defined in_arpege
CALL OPNDIA()
#endif
!
!       allocate local atmospheric variables
!
IF (.NOT.ALLOCATED(ZTA))    ALLOCATE(ZTA    (INI,IDMAX))
IF (.NOT.ALLOCATED(ZQA))    ALLOCATE(ZQA    (INI,IDMAX))
IF (.NOT.ALLOCATED(ZWIND))  ALLOCATE(ZWIND  (INI,IDMAX))
IF (.NOT.ALLOCATED(ZDIR_SW))ALLOCATE(ZDIR_SW(INI,IDMAX))
IF (.NOT.ALLOCATED(ZSCA_SW))ALLOCATE(ZSCA_SW(INI,IDMAX))
IF (.NOT.ALLOCATED(ZLW))    ALLOCATE(ZLW    (INI,IDMAX))
IF (.NOT.ALLOCATED(ZSNOW))  ALLOCATE(ZSNOW  (INI,IDMAX))
IF (.NOT.ALLOCATED(ZRAIN))  ALLOCATE(ZRAIN  (INI,IDMAX))
IF (.NOT.ALLOCATED(ZPS))    ALLOCATE(ZPS    (INI,IDMAX))
IF (.NOT.ALLOCATED(ZCO2))   ALLOCATE(ZCO2   (INI,IDMAX))
IF (.NOT.ALLOCATED(ZO3))    ALLOCATE(ZO3    (INI,IDMAX))
IF (.NOT.ALLOCATED(ZAE))    ALLOCATE(ZAE    (INI,IDMAX))
IF (.NOT.ALLOCATED(ZIMPWET))ALLOCATE(ZIMPWET(INI,NIMPUROF,IDMAX))
IF (.NOT.ALLOCATED(ZIMPDRY))ALLOCATE(ZIMPDRY(INI,NIMPUROF,IDMAX))
IF (.NOT.ALLOCATED(ZDIR))   ALLOCATE(ZDIR   (INI,IDMAX))
IF (.NOT.ALLOCATED(ZCOEF))  ALLOCATE(ZCOEF  (INI))
IF (.NOT.ALLOCATED(ZSUMZEN))ALLOCATE(ZSUMZEN(INI))
IF (.NOT.ALLOCATED(ZSW))    ALLOCATE(ZSW    (INI))
!
! ###########################################################################
IF (.NOT.ALLOCATED(ZTSUN_PREV_FORC))   ALLOCATE(ZTSUN_PREV_FORC(INI))
IF (.NOT.ALLOCATED(ZZENITH_PREV_FORC)) ALLOCATE(ZZENITH_PREV_FORC(INI))
IF (.NOT.ALLOCATED(ZAZIM_PREV_FORC))   ALLOCATE(ZAZIM_PREV_FORC(INI))
IF (.NOT.ALLOCATED(ZTSUN_NEXT_FORC))   ALLOCATE(ZTSUN_NEXT_FORC(INI))
IF (.NOT.ALLOCATED(ZZENITH_NEXT_FORC)) ALLOCATE(ZZENITH_NEXT_FORC(INI))
IF (.NOT.ALLOCATED(ZAZIM_NEXT_FORC))   ALLOCATE(ZAZIM_NEXT_FORC(INI))
! ###########################################################################
!
ZIMPWET   (:,:,:) = NUNDEF
ZIMPDRY   (:,:,:) = NUNDEF
ZO3        (:,:)  = NUNDEF
ZAE        (:,:)  = NUNDEF
!      computes initial air co2 concentration and  density
!
#ifdef SFX_MPI
XTIME = (MPI_WTIME() - XTIME0)
#endif
 CALL WLOG_MPI('COMPARE_OROGRAPHY SUNPOS ',PLOG=XTIME)
#ifdef SFX_MPI
XTIME0 = MPI_WTIME()
#endif
!
!* opens forcing files (if ASCII or BINARY)
!
IF (CFORCING_FILETYPE=='ASCII ' .OR. CFORCING_FILETYPE=='BINARY') &
        CALL OPEN_CLOSE_BIN_ASC_FORC('OPEN ',CFORCING_FILETYPE,'R')
!
CALL OL_READ_ATM(CSURF_FILETYPE, CFORCING_FILETYPE, ITIMESTARTINDEX,&
                 ZTA,ZQA,ZWIND,ZDIR_SW,ZSCA_SW,ZLW,ZSNOW,ZRAIN,ZPS, &
                 ZCO2,ZIMPWET,ZIMPDRY,ZO3,ZAE,ZDIR                  )
!
 CALL WLOG_MPI(' ')
 CALL WLOG_MPI('TIME_NPIO_READ forc ',PLOG=XTIME_NPIO_READ)
 CALL WLOG_MPI('TIME_COMM_READ forc ',PLOG=XTIME_COMM_READ)
XTIME_NPIO_READ = 0.
XTIME_COMM_READ = 0.
!
#ifdef SFX_MPI
XTIME = (MPI_WTIME() - XTIME0)
#endif
 CALL WLOG_MPI(' ')
 CALL WLOG_MPI('OL_READ_ATM0 ',PLOG=XTIME)
 CALL WLOG_MPI(' ')
#ifdef SFX_MPI
XTIME0 = MPI_WTIME()
#endif
!
XCO2 (:) = ZCO2(:,1)
XRHOA(:) = ZPS(:,1) / (XRD * ZTA(:,1) * ( 1.+((XRV/XRD)-1.)*ZQA(:,1) ) + XG * XZREF )
!Set the value of impur deposit coef
IF (LFORCATMOTARTES) THEN
  XO3(:)  = ZO3(:,1)
  XAE(:)  = ZAE(:,1)
ENDIF
IF (LFORCIMP) THEN
  DO JIMP=1,NIMPUROF
    XIMPWET(:,JIMP)=ZIMPWET(:,JIMP,1)
    XIMPDRY(:,JIMP)=ZIMPDRY(:,JIMP,1)
  ENDDO
ENDIF
!
!       surface Initialisation
!
#ifdef SFX_MPI
XTIME = (MPI_WTIME() - XTIME0)
#endif
 CALL WLOG_MPI('CO2 RHOA ',PLOG=XTIME)
!
 CALL IO_BUFF_CLEAN
!
 !CALL SURFEX_DEALLO_LIST
 !CALL SURFEX_ALLOC_LIST(IBLOCKTOT)
!
#ifdef SFX_MPI
XTIME0 = MPI_WTIME()
#endif
!
CALL GOTO_MODEL(1)
!
CALL INIT_SURF_ATM_n(YSC, CSURF_FILETYPE, YINIT, INI, NSCAL, IBANDS, CSV,   &
                     XCO2(:), XRHOA(:), XZENITH(:),XAZIM(:),XSW_BANDS,      &
                     XDIR_ALB(:,:), XSCA_ALB(:,:), XEMIS(:), XTSRAD(:),     &
                     XTSURF(:), IYEAR, IMONTH, IDAY, ZTIME, TDATE_END%TDATE,&
                     AT,YATMFILE, YATMFILETYPE, YTEST                       )
!
! initialization routines to compute shadows
IF (GSHADOWS) THEN
  IF (IBLOCK==0) THEN
    CALL INIT_SLOPE_PARAM(YSC%UG%G, YSC%UG%XGRID_FULL_PAR,ZZS_FORC,INI,ZLAT)
  END IF
  CALL LOCAL_SLOPE_PARAM(1,INI)
END IF
!
! initialization routines to define sytron grid
IF (NBLOCK==0) THEN
  IF(YSC%IM%O%LSNOWSYTRON) THEN
    CALL INIT_SYTRON_TABLE(YSC%USS,ZZS_FORC,INI,ZLAT,ZLON)
  ELSE
    ALLOCATE(NTAB_SYT(INI))
    NTAB_SYT(:)=-999
  ENDIF
ENDIF

#ifdef SFX_MPI
XTIME = (MPI_WTIME() - XTIME0)
#endif
CALL WLOG_MPI(' ')
CALL WLOG_MPI('INIT_SURF_ATM ',PLOG=XTIME)
CALL WLOG_MPI(' ')
!
CALL WLOG_MPI('TIME_NPIO_READ init ',PLOG=XTIME_NPIO_READ)
CALL WLOG_MPI('TIME_COMM_READ init ',PLOG=XTIME_COMM_READ)
CALL WLOG_MPI(' ')
!
XTIME_NPIO_READ = 0.
XTIME_COMM_READ = 0.
!

#ifdef SFX_MPI
XTIME0 = MPI_WTIME()
#endif
!
CALL INIT_CRODEBUG(YSC%IM%NPE%AL(1)%TSNOW%SCHEME)
!
! * SURFEX - OASIS  grid, partitions and local field definitions
!
IF(LOASIS)THEN
  CALL SFX_OASIS_DEF_OL(YSC%IM%O, YSC%U, YSC%UG, CSURF_FILETYPE,YALG_MPI)
ENDIF
!
! --------------------------------------------------------------------------------------
!
IF (LXIOS) THEN 
  XTSTEP_OUTPUT = XTSTEP_SURF
ENDIF
!
CALL SFX_XIOS_SETUP_OL(YSC,ILUOUT,IYEAR,IMONTH,IDAY,ZTIME,XTSTEP_OUTPUT,XSW_BANDS)
!
NWRITE = 0
!
#ifdef SFX_MPI
XTIME = (MPI_WTIME() - XTIME0)
#endif
 CALL WLOG_MPI('INIT FINISHED ',PLOG=XTIME)
#ifdef SFX_MPI
XTIME0 = MPI_WTIME()
#endif
!*      2.    Temporal loops
!
XTIME_CALC(:) = 0.
XTIME_WRITE(:) = 0.
!
LFIRST_WRITE = .TRUE.
LDEF_ol = .TRUE.
IF (CTIMESERIES_FILETYPE=="OFFLIN") CALL INIT_OUTPUT_OL_n (YSC)
!
NCPT_WRITE = 0
!
DO JFORC_STEP=1,INB_STEP_ATM
  !
#ifdef SFX_MPI
  XTIME1 = MPI_WTIME()
#endif
  ! read Forcing
  !
  !indice of forcing line in forcing arrays
  ID_FORC=JFORC_STEP-INT(JFORC_STEP/INB_LINES)*INB_LINES
  IF (ID_FORC==0) ID_FORC=INB_LINES
  !new forcings to read
  IF (ID_FORC==1 .AND. JFORC_STEP.NE.1) THEN
    !if last part of forcing, the last point has to be adjusted on the end of
    !files
    IF (JFORC_STEP/INB_LINES==NB_READ_FORC-1) THEN 
      IDMAX=INB_STEP_ATM-JFORC_STEP+1+NEND_ATM
      !for ascii and binary forcing files
      ZTA    (:,IDMAX)=ZTA    (:,SIZE(ZTA,2))
      ZQA    (:,IDMAX)=ZQA    (:,SIZE(ZTA,2))
      ZWIND  (:,IDMAX)=ZWIND  (:,SIZE(ZTA,2))
      ZDIR_SW(:,IDMAX)=ZDIR_SW(:,SIZE(ZTA,2))
      ZSCA_SW(:,IDMAX)=ZSCA_SW(:,SIZE(ZTA,2))
      ZLW    (:,IDMAX)=ZLW    (:,SIZE(ZTA,2))
      ZSNOW  (:,IDMAX)=ZSNOW  (:,SIZE(ZTA,2))
      ZRAIN  (:,IDMAX)=ZRAIN  (:,SIZE(ZTA,2))
      ZPS    (:,IDMAX)=ZPS    (:,SIZE(ZTA,2))
      ZCO2   (:,IDMAX)=ZCO2   (:,SIZE(ZTA,2))
      ZDIR   (:,IDMAX)=ZDIR   (:,SIZE(ZTA,2))
      IF (LFORCATMOTARTES) THEN  
        ZO3(:,IDMAX)=ZO3(:,SIZE(ZTA,2))
        ZAE(:,IDMAX)=ZAE(:,SIZE(ZTA,2))
      ENDIF
      IF (LFORCIMP) THEN
        DO JIMP=1,NIMPUROF
          ZIMPWET(:,JIMP,IDMAX)=ZIMPWET(:,JIMP,SIZE(ZTA,2))
          ZIMPDRY(:,JIMP,IDMAX)=ZIMPDRY(:,JIMP,SIZE(ZTA,2))
        ENDDO
      ENDIF
    ENDIF
    IFORC_STEP=ITIMESTARTINDEX+JFORC_STEP-1
    CALL OL_READ_ATM(CSURF_FILETYPE, CFORCING_FILETYPE, IFORC_STEP,       &
                     ZTA(:,1:IDMAX),ZQA(:,1:IDMAX),ZWIND(:,1:IDMAX),      &
                     ZDIR_SW(:,1:IDMAX),ZSCA_SW(:,1:IDMAX),ZLW(:,1:IDMAX),&
                     ZSNOW(:,1:IDMAX),ZRAIN(:,1:IDMAX),ZPS(:,1:IDMAX),    &
                     ZCO2(:,1:IDMAX),ZIMPWET(:,:,1:IDMAX),                &
                     ZIMPDRY(:,:,1:IDMAX),ZO3(:,1:IDMAX),ZAE(:,1:IDMAX),  &
                     ZDIR(:,1:IDMAX) )
  ENDIF

#ifdef SFX_MPI
  XTIME_CALC(1) = XTIME_CALC(1) + (MPI_WTIME() - XTIME1)
  XTIME1 = MPI_WTIME()
#endif
  !
  !COMPUTE SUM ZENITH angle between 2 timestepA
  ZSUMZEN(:)=0.0
  DO JSURF_STEP = 1,INB_ATM
    IDAY2  = IDAY
    ZTIME2 = ZTIME + (JSURF_STEP-1.)*XTSTEP_SURF
    IF (ZTIME2>86400.) THEN
      ZTIME2 = ZTIME2-86400
      IDAY2  = IDAY+1
    ENDIF
    CALL SUNPOS(IYEAR, IMONTH, IDAY2, ZTIME+(JSURF_STEP-1.)*XTSTEP_SURF, &
                ZLON, ZLAT, XTSUN, XZENITH, XAZIM)
    !
    ZSUMZEN(:)= ZSUMZEN(:) + MAX(COS(XZENITH(:)+0.1),0.)/(INB_ATM*1.0)
    !
  ENDDO
  WHERE ( ZSUMZEN<0.01 ) ZSUMZEN = 0.0
  !

! ###########################################################################
  !COMPUTE SUM ZENITH angle between 2 timestep
  !ZSUMZEN(:)=0.0
  !DO JSURF_STEP = 1,INB_ATM
  !  IDAY2  = IDAY
  !  ZTIME2 = ZTIME + (JSURF_STEP-1.)*XTSTEP_SURF
  !  IF (ZTIME2>86400.) THEN
  !    ZTIME2 = ZTIME2-86400
  !    IDAY2  = IDAY+1
  !  ENDIF
  !
  !
  ! Calculation of zenith angle for previous and next forcing time step
  !
  !
  ! Previous forcing time step
  ! 
  IYEAR_PREV  = IYEAR
  IMONTH_PREV = IMONTH
  IDAY_PREV   = IDAY
  ZTIME_PREV  = ZTIME
  !
  CALL SUNPOS(IYEAR_PREV, IMONTH_PREV, IDAY_PREV, ZTIME_PREV, ZLON, ZLAT, ZTSUN_PREV_FORC, ZZENITH_PREV_FORC, ZAZIM_PREV_FORC)
  !
  !
  ! Next forcing time step
  !
  IYEAR_NEXT  = IYEAR
  IMONTH_NEXT = IMONTH
  IDAY_NEXT   = IDAY
  ZTIME_NEXT  = ZTIME + INB_ATM*XTSTEP_SURF
  !
  IF (ZTIME_NEXT .GT. 86400.) THEN
     ZTIME_NEXT = ZTIME_NEXT - 86400.0
     IDAY_NEXT  = IDAY_NEXT + 1
  ENDIF
  !
  !
  CALL SUNPOS(IYEAR_NEXT, IMONTH_NEXT, IDAY_NEXT, ZTIME_NEXT, ZLON, ZLAT, ZTSUN_NEXT_FORC, ZZENITH_NEXT_FORC, ZAZIM_NEXT_FORC)
  !
  ! ###########################################################################

  DO JSURF_STEP=1,INB_ATM
    !
    !
    LLAST_TIMESTEP = (JFORC_STEP==INB_STEP_ATM) .AND. (JSURF_STEP==INB_ATM)
    !
    !
    ! time interpolation of the forcing
    !
#ifdef SFX_MPI
    XTIME1 = MPI_WTIME()
#endif
    !
    CALL SUNPOS(IYEAR, IMONTH, IDAY, ZTIME, ZLON, ZLAT, XTSUN, XZENITH, XAZIM)
    IYEAR2 = IYEAR
    IMONTH2= IMONTH
    IDAY2  = IDAY
    ZTIME2 = ZTIME+XTSTEP_SURF
    CALL ADD_FORECAST_TO_DATE_SURF(IYEAR2, IMONTH2, IDAY2, ZTIME2)
    CALL SUNPOS(IYEAR2, IMONTH2, IDAY2, ZTIME2, ZLON, ZLAT, XTSUN, XZENITH2, XAZIM)
    !
#ifdef SFX_MPI
    XTIME_CALC(2) = XTIME_CALC(2) + (MPI_WTIME() - XTIME1)
    XTIME1 = MPI_WTIME()
#endif

! ###########################################################################
    !interpolation between beginning and end of current forcing time step
    CALL OL_TIME_INTERP_ATM(JSURF_STEP,INB_ATM,ZPDISTRIB,           &
                            ZTA    (:,ID_FORC:ID_FORC+NEND_ATM),    &
                            ZQA    (:,ID_FORC:ID_FORC+NEND_ATM),    &
                            ZWIND  (:,ID_FORC:ID_FORC+NEND_ATM),    &
                            ZDIR_SW(:,ID_FORC:ID_FORC+NEND_ATM),    &
                            ZSCA_SW(:,ID_FORC:ID_FORC+NEND_ATM),    &
                            ZLW    (:,ID_FORC:ID_FORC+NEND_ATM),    &
                            ZSNOW  (:,ID_FORC:ID_FORC+NEND_ATM),    &
                            ZRAIN  (:,ID_FORC:ID_FORC+NEND_ATM),    &
                            ZPS    (:,ID_FORC:ID_FORC+NEND_ATM),    &
                            ZCO2   (:,ID_FORC:ID_FORC+NEND_ATM),    &
                            ZDIR   (:,ID_FORC:ID_FORC+NEND_ATM),    &
                            ZO3    (:,ID_FORC:ID_FORC+NEND_ATM),    &
                            ZAE    (:,ID_FORC:ID_FORC+NEND_ATM),    &
                            ZIMPWET(:,:,ID_FORC:ID_FORC+NEND_ATM),  &
                            ZIMPDRY(:,:,ID_FORC:ID_FORC+NEND_ATM),  &
                            ZSUMZEN,ZZENITH_PREV_FORC,XZENITH2,ZZENITH_NEXT_FORC )
! ###########################################################################

#ifdef SFX_MPI
    XTIME_CALC(3) = XTIME_CALC(3) + (MPI_WTIME() - XTIME1)
    XTIME1 = MPI_WTIME()
#endif
    !
    IF(LADAPT_SW)THEN
      !
      ! coherence between solar zenithal angle and radiation
      ! when solar beam close to horizontal -> reduction of direct radiation to
      ! the benefit of scattered radiation
      ! when pi/2 - 0.1 < ZENITH < pi/2 - 0.05 => weight of direct to scattered radiation decreases linearly with zenith
      ! when pi/2 - 0.05 < ZENITH => all the direct radiation is converted to scattered radiation
      ! coherence between solar zenithal angle and radiation
      !
      ZCOEF(:) = (XPI/2. - XZENITH(:) - 0.05) / 0.05
      ZCOEF(:) = MAX(MIN(ZCOEF,1.),0.)
      DO JLOOP=1,SIZE(XDIR_SW,2)
        XSCA_SW(:,JLOOP) = XSCA_SW(:,JLOOP) + XDIR_SW(:,JLOOP) * (1 - ZCOEF)
        XDIR_SW(:,JLOOP) = XDIR_SW(:,JLOOP) * ZCOEF(:)
      ENDDO
      !
    ELSE
      !
      ZSW(:) = 0.
      DO JLOOP=1,SIZE(XDIR_SW,2)
        ZSW(:) = ZSW(:) + XDIR_SW(:,JLOOP) + XSCA_SW(:,JLOOP)
      END DO
      WHERE (ZSW(:)>0.)
        XZENITH  = MIN (XZENITH ,XPI/2.-0.01)
        XZENITH2 = MIN (XZENITH2,XPI/2.-0.01)
      ELSEWHERE
        XZENITH  = MAX (XZENITH ,XPI/2.)
        XZENITH2 = MAX (XZENITH2,XPI/2.)
      END WHERE
      !
    ENDIF
    !
    ! updates time
    ZTIMEC= ZTIMEC+XTSTEP_SURF
    IF (LCOUPL_TOPD) LTOPD_STEP = ( MOD((((JFORC_STEP-1)*INB_ATM)+JSURF_STEP),NNB_TOPD) == 0 )
    !
    ! run Surface
    !
#ifdef SFX_MPI
    XTIME_CALC(4) = XTIME_CALC(4) + (MPI_WTIME() - XTIME1)
#endif
    !
    CALL IO_BUFF_CLEAN
    !
    IF(LOASIS)THEN
     ! Receive fields to other models proc by proc
     CALL SFX_OASIS_RECV_OL(YSC%FM%F, YSC%IM, YSC%SM%S, YSC%U, YSC%WM%W, &
                            YSC%TM, YSC%GDM, YSC%GRM, &
                            CSURF_FILETYPE, INI, IBANDS, ZTIMEC, XTSTEP_SURF, XZENITH, &
                            XSW_BANDS, XTSRAD, XDIR_ALB, XSCA_ALB, XEMIS, XTSURF   )
    ENDIF
    !
#ifdef SFX_MPI
    XTIME1 = MPI_WTIME()
#endif
    !
! IF SNOWCRORAD=TARTES then spectral repartition of direct and diffuse radiation
    ! for any spectral calculation regarding snow (tartes and/or atmotartes)
    IF (LSPECSNOW) THEN
      XSCA_SW(:,2:IBANDS)=0.
      XDIR_SW(:,2:IBANDS)=0.
    ENDIF
!	write(*,*) IMONTH, IDAY, ZTIME/(3600.)
    IF (YSC%IM%O%LATMORAD) THEN
      CALL RADIANCE (INI,IYEAR, IMONTH, IDAY, &
                     IBANDS,LFORCATMOTARTES,XSCA_SW,XDIR_SW)
    ENDIF
    !
    IF(GSHADOWS) THEN
      CALL SLOPE_RADIATIVE_EFFECT(XTSTEP_SURF, XZENITH, XAZIM, XPS, XTA, XRAIN, XDIR_SW, XLW, &
                                  XZS_THREAD, XZS_XY_THREAD, XSLOPANG_THREAD, XSLOPAZI_THREAD,&
                                  XSURF_TRIANGLE_THREAD)
    END IF
    !
    CALL COUPLING_SURF_ATM_n(YSC, CSURF_FILETYPE, 'E', ZTIMEC, XTSTEP_SURF,  &
                            IYEAR, IMONTH, IDAY, ZTIME, INI, NSCAL, IBANDS,  &
                            XTSUN, XZENITH, XZENITH2, XAZIM, XZREF, XUREF,   &
                            XZS ,XU, XV, XQA, XTA, XRHOA, XSV, XCO2, XIMPWET,&
                            XIMPDRY, CSV, XRAIN, XSNOW, XLW, XDIR_SW,        &
                            XSCA_SW, XSW_BANDS, XPS, XPA, XSFTQ, XSFTH,      &
                            XSFTS, XSFCO2, XSFU, XSFV, XTSRAD, XDIR_ALB,     &
                            XSCA_ALB, XEMIS, XTSURF, XZ0, XZ0H, XQSURF,      &
                            XPEW_A_COEF, XPEW_B_COEF,XPET_A_COEF,XPEQ_A_COEF,&
                            XPET_B_COEF,XPEQ_B_COEF, XZWS, YTEST      )
    !
#ifdef SFX_MPI
    XTIME_CALC(5) = XTIME_CALC(5) + (MPI_WTIME() - XTIME1)
#endif
    !
#ifdef SFX_MPI
    XTIME1 = MPI_WTIME()
#endif
    !
    IF(LOASIS)THEN
     ! Send fields to other models proc by proc
     CALL SFX_OASIS_SEND_OL(YSC%FM%F, YSC%IM, YSC%SM%S, YSC%U, YSC%WM%W, YSC%SM%SD%D, CSURF_FILETYPE,INI,ZTIMEC,XTSTEP_SURF)
    ENDIF
    !
    ZTIME = ZTIME + XTSTEP_SURF
    CALL ADD_FORECAST_TO_DATE_SURF(IYEAR, IMONTH, IDAY, ZTIME)
#ifdef SFX_MPI
    XTIME_CALC(6) = XTIME_CALC(6) + (MPI_WTIME() - XTIME1)
    !
    XTIME1 =  MPI_WTIME()
#endif
    ! ecrit Surface
    !
    IF ( LCOUPL_TOPD .AND. LTOPD_STEP ) THEN
      !
      IF (.NOT.ALLOCATED(NYEAR))  ALLOCATE(NYEAR(NNB_TOPD_STEP))
      IF (.NOT.ALLOCATED(NMONTH)) ALLOCATE(NMONTH(NNB_TOPD_STEP))
      IF (.NOT.ALLOCATED(NDAY))   ALLOCATE(NDAY(NNB_TOPD_STEP))
      IF (.NOT.ALLOCATED(NH))     ALLOCATE(NH(NNB_TOPD_STEP))
      IF (.NOT.ALLOCATED(NM))     ALLOCATE(NM(NNB_TOPD_STEP))
      !
      NYEAR (NTOPD_STEP) = IYEAR
      NMONTH(NTOPD_STEP) = IMONTH
      NDAY  (NTOPD_STEP) = IDAY
      NH    (NTOPD_STEP) = INT(ZTIME/3600.)
      NM    (NTOPD_STEP) = INT((ZTIME-NH(NTOPD_STEP)*3600.)/60.)
      !
      IF ( NM(NTOPD_STEP)==60 ) THEN
        !
        NM(NTOPD_STEP) = 0
        NH(NTOPD_STEP) = NH(NTOPD_STEP)+1
        !
      ENDIF
      !
      IF ( NH(NTOPD_STEP)==24 ) THEN
        !
        NH  (NTOPD_STEP) = 0
        NDAY(NTOPD_STEP) = NDAY(NTOPD_STEP)+1
        !
        !!AJOUT BEC
        SELECT CASE (NMONTH(NTOPD_STEP))
          CASE(4,6,9,11)
            IF ( NDAY(NTOPD_STEP)==31 ) THEN
              NMONTH(NTOPD_STEP) = NMONTH(NTOPD_STEP)+1
              NDAY  (NTOPD_STEP) = 1
            ENDIF
          CASE(1,3,5,7:8,10)
            IF ( NDAY(NTOPD_STEP)==32 ) THEN
              NMONTH(NTOPD_STEP) = NMONTH(NTOPD_STEP)+1
              NDAY  (NTOPD_STEP) = 1
            ENDIF
          CASE(12)
            IF ( NDAY(NTOPD_STEP)==32 ) THEN
              NYEAR (NTOPD_STEP) = NYEAR(NTOPD_STEP)+1
              NMONTH(NTOPD_STEP) = 1
              NDAY  (NTOPD_STEP) = 1
            ENDIF
          CASE(2)
            IF( MOD(NYEAR(NTOPD_STEP),4)==0 .AND. MOD(NYEAR(NTOPD_STEP),100)/=0 .OR. MOD(NYEAR(NTOPD_STEP),400)==0 ) THEN
              IF (NDAY(NTOPD_STEP)==30) THEN
                NMONTH(NTOPD_STEP) = NMONTH(NTOPD_STEP)+1
                NDAY  (NTOPD_STEP) = 1
              ENDIF
            ELSE
              IF (NDAY(NTOPD_STEP)==29) THEN
                NMONTH(NTOPD_STEP) = NMONTH(NTOPD_STEP)+1
                NDAY  (NTOPD_STEP) = 1
              ENDIF
            ENDIF
        END SELECT
        !
      ENDIF
      !
      ! * 2. Stocking date of each time step
      !
      NTOPD_STEP = NTOPD_STEP + 1
      !
    ENDIF
    !
    ZTIME_OUT  = ZTIME
    IDAY_OUT   = IDAY
    IMONTH_OUT = IMONTH
    IYEAR_OUT  = IYEAR
    !
    IF(ZTIME==0.0)THEN
      ZTIME_OUT = 86400.
      IDAY_OUT   = IDAY-1
      IF(IDAY_OUT==0)THEN
        IMONTH_OUT = IMONTH - 1
        IF(IMONTH_OUT==0)THEN
          IMONTH_OUT=12
          IYEAR_OUT = IYEAR - 1
        ENDIF
        SELECT CASE(IMONTH_OUT)
               CASE(4,6,9,11)
                 IDAY_OUT=30
               CASE(1,3,5,7:8,10,12)
                 IDAY_OUT=31
               CASE(2)
                 IF(((MOD(IYEAR_OUT,4)==0).AND.(MOD(IYEAR_OUT,100)/=0)).OR.(MOD(IYEAR_OUT,400)==0))THEN
                   IDAY_OUT=29
                 ELSE
                   IDAY_OUT=28
                 ENDIF
        END SELECT
      ENDIF
    ENDIF
    !
    IF ((MOD(ZTIMEC,XTSTEP_OUTPUT) == 0. .AND. CTIMESERIES_FILETYPE/='NONE  '.AND. .NOT. (LFIX_OUTPUT)) .OR. &
    (MOD(ZTIME,XTSTEP_OUTPUT) == 0. .AND. CTIMESERIES_FILETYPE/='NONE  '.AND. (LFIX_OUTPUT))) THEN !B. Cluzet added a case for fixed output at 06 12 18 and 00 hours or more generally at 00 + n*XTIMESTEP
      !
      IF (LFIX_OUTPUT) THEN
        XTIMEC = ZTIMEC
      ENDIF
      !
      IF (NRANK==NPIO) THEN
        !
        !* name of the file
        IF (CTIMESERIES_FILETYPE=="ASCII " .OR. &
            CTIMESERIES_FILETYPE=="LFI   " .OR. &
            CTIMESERIES_FILETYPE=="FA    " .OR. &
            CTIMESERIES_FILETYPE=="NC    "    ) THEN  
          !
          IF(LOUT_TIMENAME)THEN
            ! if true, change the name of output file at the end of a day
            ! (ex: 19860502_00h00 -> 19860501_24h00)
            WRITE(YTAG,FMT='(I4.4,I2.2,I2.2,A1,I2.2,A1,I2.2)') IYEAR_OUT,IMONTH_OUT,IDAY_OUT,&
                 '_',INT(ZTIME_OUT/3600.),'h',NINT(ZTIME_OUT)/60-60*INT(ZTIME_OUT/3600.)
            !
          ELSE
            ! if false, default
            WRITE(YTAG,FMT='(I4.4,I2.2,I2.2,A1,I2.2,A1,I2.2)') IYEAR,IMONTH,IDAY,&
                 '_',INT(ZTIME/3600.),'h',NINT(ZTIME)/60-60*INT(ZTIME/3600.)                  
          ENDIF
          !
          CFILEOUT    = ADJUSTL(ADJUSTR(CSURFFILE)//'.'//YTAG(1:LEN_TRIM(YTAG))//'.txt')
          CFILEOUT_LFI= ADJUSTL(ADJUSTR(CSURFFILE)//'.'//YTAG(1:LEN_TRIM(YTAG)))
          CFILEOUT_FA = ADJUSTL(ADJUSTR(CSURFFILE)//'.'//YTAG(1:LEN_TRIM(YTAG))//'.fa')
          CFILEOUT_NC = ADJUSTL(ADJUSTR(CSURFFILE)//'.'//YTAG(1:LEN_TRIM(YTAG))//'.nc')
          !
          IF (CTIMESERIES_FILETYPE=='FA    ') THEN
#ifdef SFX_FA
            LFANOCOMPACT = LDIAG_FA_NOCOMPACT
            IDATEF(1)= IYEAR
            IDATEF(2)= IMONTH
            IDATEF(3)= IDAY
            !ZTIME instead of ZTIME_OUT (FA XRD do not like 24h)
            IDATEF(4)= FLOOR(ZTIME/3600.)
            IDATEF(5)= FLOOR(ZTIME/60.) - IDATEF(4) * 60
            IDATEF(6)= NINT(ZTIME) - IDATEF(4) * 3600 - IDATEF(5) * 60
            IDATEF(7:11) = 0
            NUNIT_FA = 19
            IF (CSURF_FILETYPE/='FA    ') THEN
              CALL WRITE_HEADER_FA(YSC%GCP, YSC%UG%G%CGRID, YSC%UG%XGRID_FULL_PAR, CSURF_FILETYPE,'ALL')
            ELSE
              CALL FAITOU(IRET,NUNIT_FA,.TRUE.,CFILEOUT_FA,'UNKNOWN',.TRUE.,.FALSE.,IVERBFA,0,INB,CDNOMC)
            ENDIF
            CALL FANDAR(IRET,NUNIT_FA,IDATEF)
#endif
          END IF
          !
        END IF
        !
        XSTARTW = XSTARTW + 1
        NWRITE  = NWRITE  + 1
        LTIME_WRITTEN=.FALSE.
        !
      ENDIF
      !
#ifdef SFX_MPI
      XTIME_WRITE(1) = XTIME_WRITE(1) + (MPI_WTIME() - XTIME1)
#endif
      !
      INW = 1
      !
      IF ( LXIOS .AND. .NOT. LXIOS_DEF_CLOSED ) INW = 2
      !
      IF (LXIOS) THEN
        NTIMESTEP=INT(ZTIMEC/XTSTEP_OUTPUT)
      ENDIF
      !
      LDEF_nc = .FALSE.
      IF (CTIMESERIES_FILETYPE=="NC    ") THEN
        LDEF_nc = .TRUE.
        INW = 2
        CALL INIT_OUTPUT_NC_n (YSC%TM%BDD, YSC%CHE, YSC%CHN, YSC%CHU, YSC%SM%DTS, &
                               YSC%TM%DTT, YSC%DTZ, YSC%IM, YSC%UG, YSC%U, YSC%DUO%CSELECT)
      ENDIF
      !
      LDEF_ol = .FALSE.
      IF (CTIMESERIES_FILETYPE=="OFFLIN".AND.LFIRST_WRITE) THEN
        LDEF_ol = .TRUE.
        INW = 2
      ENDIF
      !
      IDX_W = 0
      !
      DO JNW = 1,INW
        !
        CALL IO_BUFF_CLEAN
        !
#ifdef SFX_MPI
        XTIME1 =  MPI_WTIME()
#endif
        !
        IF (JNW.EQ.INW) THEN
           LRESET_DIAG_ol = .TRUE.
           LRESET_DIAG_nc = .TRUE.
        ELSE
           LRESET_DIAG_ol = .FALSE.
           LRESET_DIAG_nc = .FALSE.
        ENDIF
        !
        CALL WRITE_SURF_ATM_n(YSC, CTIMESERIES_FILETYPE,'ALL')
#ifdef SFX_MPI
        XTIME_WRITE(2) = XTIME_WRITE(2) + (MPI_WTIME() - XTIME1)
        XTIME1 =  MPI_WTIME()
#endif
        CALL DIAG_SURF_ATM_n(YSC, CTIMESERIES_FILETYPE)
#ifdef SFX_MPI
        XTIME_WRITE(3) = XTIME_WRITE(3) + (MPI_WTIME() - XTIME1)
        XTIME1 =  MPI_WTIME()
#endif
        !
#ifdef WXIOS
        IF (LXIOS_DEF_CLOSED) THEN
          CALL XIOS_UPDATE_CALENDAR(NTIMESTEP)
        ENDIF
#endif
        !
        CALL WRITE_DIAG_SURF_ATM_n(YSC, CTIMESERIES_FILETYPE,'ALL')
        !
        IF (LXIOS) THEN
#ifdef WXIOS
          IF (.NOT. LXIOS_DEF_CLOSED) THEN
            CALL XIOS_CLOSE_CONTEXT_DEFINITION()
            LXIOS_DEF_CLOSED=.TRUE.
          ENDIF
#endif
        ENDIF
        !
#ifdef SFX_MPI
        XTIME_WRITE(4) = XTIME_WRITE(4) + (MPI_WTIME() - XTIME1)
#endif
        !
        IF (CTIMESERIES_FILETYPE=="NC    ") THEN
            IF (JNW==1) THEN
                CALL ENDDEF_IO_SURF_NC_n()
            ELSE
                CALL CLOSE_IO_SURF_NC_n()
            ENDIF
        ENDIF
        LDEF_nc = .FALSE.
        LDEF_ol = .FALSE.
        !
        NCPT_WRITE = 0
        !
        LFIRST_WRITE = .FALSE.
        !
        IF (LLAST_TIMESTEP .AND. (.NOT. LRESTART)) THEN
           LRESET_DIAG_ol = .TRUE.
           LRESET_DIAG_nc = .TRUE.
        ENDIF
      ENDDO
      !
      IF (LCOUPL_TOPD .AND. NTOPD_STEP > NNB_TOPD_STEP) THEN
        !
        ! Writing of file resulting of coupling with TOPMODEL or routing ****
        CALL WRITE_DISCHARGE_FILE(CSURF_FILETYPE,'q_total.txt','FORMATTED',&
                                  NYEAR,NMONTH,NDAY,NH,NM,XQTOT)
        CALL WRITE_DISCHARGE_FILE(CSURF_FILETYPE,'q_runoff.txt','FORMATTED',&
                                  NYEAR,NMONTH,NDAY,NH,NM,XQB_RUN)
        CALL WRITE_DISCHARGE_FILE(CSURF_FILETYPE,'q_drainage.txt','FORMATTED',&
                                  NYEAR,NMONTH,NDAY,NH,NM,XQB_DR)
        ! Writing of budget files
        IF (LBUDGET_TOPD) CALL WRITE_BUDGET_COUPL_ROUT
        ! Writing results on sub-catchments
        IF (LSUBCAT) CALL WRITE_DISCHARGE_FILE_SUB(CSURF_FILETYPE,&
                                'q_total_sub.txt','FORMATTED',&
                                NYEAR,NMONTH,NDAY,NH,NM)
        !
      ENDIF
      !
#ifdef SFX_MPI
      XTIME1 =  MPI_WTIME()
#endif
      !
      IF (NRANK==NPIO) THEN
        IF (CTIMESERIES_FILETYPE=='FA    ') THEN
#ifdef SFX_FA
          CALL FAIRME(IRET,NUNIT_FA,'UNKNOWN')
#endif
        END IF
        !* add informations in the file
        IF (CTIMESERIES_FILETYPE=='LFI   ' .AND. LMNH_COMPATIBLE) CALL WRITE_HEADER_MNH
      ENDIF
#ifdef SFX_MPI
      XTIME_WRITE(5) = XTIME_WRITE(5) + (MPI_WTIME() - XTIME1)
#endif
      !
      !
       ELSEIF (MOD(ZTIMEC,XTSTEP_OUTPUT) == 0. .AND. CTIMESERIES_FILETYPE=='NONE  '&
       .AND.LCOUPL_TOPD .AND. NTOPD_STEP > NNB_TOPD_STEP) THEN
         !
         ! Writing of file resulting of coupling with TOPMODEL or routing ****
         CALL WRITE_DISCHARGE_FILE(CSURF_FILETYPE,'q_total.txt','FORMATTED',&
                                   NYEAR,NMONTH,NDAY,NH,NM,XQTOT)
         CALL WRITE_DISCHARGE_FILE(CSURF_FILETYPE,'q_runoff.txt','FORMATTED',&
                                   NYEAR,NMONTH,NDAY,NH,NM,XQB_RUN)
         CALL WRITE_DISCHARGE_FILE(CSURF_FILETYPE,'q_drainage.txt','FORMATTED',&
                                   NYEAR,NMONTH,NDAY,NH,NM,XQB_DR)
         ! Writing of budget files
         IF (LBUDGET_TOPD) CALL WRITE_BUDGET_COUPL_ROUT
         !
         ! Writing results on sub-catchments
         IF (LSUBCAT) CALL WRITE_DISCHARGE_FILE_SUB(CSURF_FILETYPE,&
                                'q_total_sub.txt','FORMATTED',&
                                NYEAR,NMONTH,NDAY,NH,NM)
       !
    ENDIF
    !
  END DO
  !
  IF (NRANK==NPIO) THEN
    IF (LPRINT) THEN
      IF (MOD(ZTIMEC,XDAY) == 0.) THEN
        ICOUNT = ICOUNT + 1
        CALL WLOG_MPI('SFX DAY :',KLOG=ICOUNT,KLOG2=INT(ZDURATION/XDAY))
        WRITE(YTAG,FMT='(I4.4,A1,I2.2,A1,I2.2,A4,I2.2,A1,I2.2)') IYEAR_OUT,'-',IMONTH_OUT,'-',IDAY_OUT,&
                       ' at ',INT(ZTIME_OUT/3600.),'h',NINT(ZTIME_OUT)/60-60*INT(ZTIME_OUT/3600.)
        WRITE(*,'(A10,I5,A2,I5,A17,A20)')'SFX  DAY :',ICOUNT,' /',INT(ZDURATION/XDAY), &
                                         ' ; DATE ENDING : ',ADJUSTL(ADJUSTR(YTAG))
      ENDIF
    ENDIF
  ENDIF
  !
END DO

IF ( (TRIM(CTIMESERIES_FILETYPE)=="XIOS") ) THEN
    CALL SFX_XIOS_INIT_OUTPUT_OL(YSC)
END IF
!
 CALL WLOG_MPI(' ')
 CALL WLOG_MPI('OL_READ_ATM ',PLOG=XTIME_CALC(1))
 CALL WLOG_MPI('SUNPOS ',PLOG=XTIME_CALC(2))
 CALL WLOG_MPI('OL_TIME_INTERP_ATM ',PLOG=XTIME_CALC(3))
 CALL WLOG_MPI('')
 CALL WLOG_MPI('ZENITH ',PLOG=XTIME_CALC(4))
 CALL WLOG_MPI('')
 CALL WLOG_MPI('COUPLING_SURF_ATM ',PLOG=XTIME_CALC(5))
 CALL WLOG_MPI('')
 CALL WLOG_MPI('ADD_FORECAST_TO_DATE_SURF ',PLOG=XTIME_CALC(6))
 CALL WLOG_MPI('DEF_DATE ',PLOG=XTIME_WRITE(1))
 CALL WLOG_MPI('')
 CALL WLOG_MPI('WRITE_SURF_ATM ',PLOG=XTIME_WRITE(2))
 CALL WLOG_MPI('DIAG_SURF_ATM ',PLOG=XTIME_WRITE(3))
 CALL WLOG_MPI('WRITE_DIAG_SURF_ATM ',PLOG=XTIME_WRITE(4))
 CALL WLOG_MPI('')
 CALL WLOG_MPI('CLOSE FILES ',PLOG=XTIME_WRITE(5))
 CALL WLOG_MPI('')
!
IF (CFORCING_FILETYPE=='ASCII ' .OR. CFORCING_FILETYPE=='BINARY') &
        CALL OPEN_CLOSE_BIN_ASC_FORC('CLOSE',CFORCING_FILETYPE,'R')
!
IF (CFORCING_FILETYPE=='NETCDF') CALL CLOSE_FILEIN_OL
IF (CTIMESERIES_FILETYPE=='OFFLIN') CALL CLOSE_FILEOUT_OL
!
! --------------------------------------------------------------------------------------
!
!*    3.     write restart file
!            ------------------
!
IF ( LRESTART ) THEN
  !
  LFIRST_WRITE = .TRUE.
  !
  IF(CSURF_FILETYPE=='FA'.OR.YSC%DUO%LDIAG_MIP)THEN
    GDIAG_RESTART = .FALSE.
  ELSE
    GDIAG_RESTART = .TRUE.
  ENDIF
  !
  IF (NRANK==NPIO) THEN
    !* name of the file
    CFILEOUT    = ADJUSTL(ADJUSTR(CSURFFILE)//'.txt')
    CFILEOUT_LFI= CSURFFILE
    CFILEOUT_FA = ADJUSTL(ADJUSTR(CSURFFILE)//'.fa')
    CFILEOUT_NC = ADJUSTL(ADJUSTR(CSURFFILE)//'.nc')

    !* opens the file
    IF (CSURF_FILETYPE=='FA    ') THEN
#ifdef SFX_FA
      LFANOCOMPACT = .TRUE.
      IDATEF(1)= IYEAR
      IDATEF(2)= IMONTH
      IDATEF(3)= IDAY
      IDATEF(4)= FLOOR(ZTIME/3600.)
      IDATEF(5)= FLOOR(ZTIME/60.) - IDATEF(4) * 60
      IDATEF(6)= NINT(ZTIME) - IDATEF(4) * 3600 - IDATEF(5) * 60
      IDATEF(7:11) = 0
      NUNIT_FA = 19
      CALL FAITOU(IRET,NUNIT_FA,.TRUE.,CFILEOUT_FA,'UNKNOWN',.TRUE.,.FALSE.,IVERBFA,0,INB,CDNOMC)
      CALL FANDAR(IRET,NUNIT_FA,IDATEF)
#endif
    END IF
    !
  ENDIF
  !
  INW = 1
  IF (CSURF_FILETYPE=="NC    ") INW = 2
  !
  LDEF_nc = .TRUE.
  LDEF_ol = .TRUE.
  !
  IF (ASSOCIATED(YSC%DUO%CSELECT)) DEALLOCATE(YSC%DUO%CSELECT)
  ALLOCATE(YSC%DUO%CSELECT(0))
  !
  IF (CSURF_FILETYPE=="NC    ") THEN
    CALL INIT_OUTPUT_NC_n (YSC%TM%BDD, YSC%CHE, YSC%CHN, YSC%CHU, YSC%SM%DTS, YSC%TM%DTT, &
                           YSC%DTZ, YSC%IM, YSC%UG, YSC%U, YSC%DUO%CSELECT)
  ENDIF
  !
  DO JNW = 1,INW
    !
    CALL IO_BUFF_CLEAN
    !
    CALL FLAG_UPDATE(YSC%IM%ID%O, YSC%DUO,.FALSE.,.TRUE.,.FALSE.,.FALSE.,.FALSE.)
    !
    IF (LRESTART_2M) THEN
      I2M       = 1
      GPGD_ISBA = .TRUE.
    ELSE
      I2M       = 0
      GPGD_ISBA = .FALSE.
    ENDIF
    GFRAC                  = .TRUE.
    GDIAG_GRID             = .TRUE.
    GSURF_BUDGET           = .FALSE.
    GLUTILES_BUDGET        = .FALSE.
    GRAD_BUDGET            = .FALSE.
    GCOEF                  = .FALSE.
    GSURF_VARS             = .FALSE.
    IBEQ                   = 0
    IDSTEQ                 = 0
    GDIAG_OCEAN            = .FALSE.
    GDIAG_SEAICE           = .FALSE.
    GWATER_PROFILE         = .FALSE.
    GSEDIM_PROFILE         = .FALSE.
    GSURF_EVAP_BUDGET      = .FALSE.
    GFLOOD                 = .FALSE.
    GCH_NO_FLUX_ISBA       = .FALSE.
    GSURF_MISC_BUDGET_ISBA = .FALSE.
    GPGD_TEB               = .FALSE.
    GSURF_MISC_BUDGET_TEB  = .FALSE.  
    GFLKFLUX               = .FALSE.
    GFLKWATER              = .FALSE.
    !
    CALL FLAG_DIAG_UPDATE(YSC%FM, YSC%IM, YSC%SM, YSC%TM, YSC%WM, YSC%DUO, YSC%U, YSC%SV, &
                          GFRAC, GDIAG_GRID, I2M, GSURF_BUDGET, GRAD_BUDGET, GCOEF,       &
                          GSURF_VARS, IBEQ, IDSTEQ, GDIAG_OCEAN, GDIAG_SEAICE,            &
                          GWATER_PROFILE, GSEDIM_PROFILE, GSURF_EVAP_BUDGET, GFLOOD,      &
                          GPGD_ISBA, GCH_NO_FLUX_ISBA, GSURF_MISC_BUDGET_ISBA, GPGD_TEB,  &
                          GSURF_MISC_BUDGET_TEB, GLUTILES_BUDGET, GFLKFLUX, GFLKWATER, LCROCO)
    !
    YSC%DUO%LSNOWDIMNC = .FALSE.
    !
    !* writes into the file
    CALL WRITE_SURF_ATM_n(YSC, CSURF_FILETYPE,'ALL')
    IF(GDIAG_RESTART.OR.LRESTART_2M.OR.LCROCO) THEN
       CALL WRITE_DIAG_SURF_ATM_n(YSC, CSURF_FILETYPE,'ALL')
    ENDIF
    !
    IF (CSURF_FILETYPE=="NC    ") THEN
      IF (JNW==1) THEN
        CALL ENDDEF_IO_SURF_NC_n()
      ELSE
        CALL CLOSE_IO_SURF_NC_n()
      ENDIF
    ENDIF
    IF (LLAST_TIMESTEP) THEN 
       LRESET_DIAG_ol = .TRUE.
       LRESET_DIAG_nc = .TRUE.
    ENDIF 
    !
    LDEF_nc = .FALSE.
    LDEF_ol = .FALSE.
    !
    NCPT_WRITE = 0
    !
    LFIRST_WRITE = .FALSE.
    !
  ENDDO
  !
  !* closes the file
  IF (NRANK==0 ) THEN
    IF (CSURF_FILETYPE=='FA    ') THEN
#ifdef SFX_FA
      CALL FAIRME(IRET,NUNIT_FA,'UNKNOWN')
#endif
    END IF
    !* add informations in the file
    IF (CSURF_FILETYPE=='LFI   ' .AND. LMNH_COMPATIBLE) CALL WRITE_HEADER_MNH
    !
  ENDIF
  !
  IF (LCOUPL_TOPD .AND. NTOPD_STEP > NNB_TOPD_STEP) &
          CALL PREP_RESTART_COUPL_TOPD(YSC%UG, YSC%U, CSURF_FILETYPE,INI)
  !
END IF
!
! --------------------------------------------------------------------------------------
!
!*    4.     inquiry mode
!            ------------
!
IF ( LINQUIRE ) THEN
  !
  ALLOCATE( ZSEA       ( INI ) )
  ALLOCATE( ZWATER     ( INI ) )
  ALLOCATE( ZNATURE    ( INI ) )
  ALLOCATE( ZTOWN      ( INI ) )
  ALLOCATE( ZT2M       ( INI ) )
  ALLOCATE( ZQ2M       ( INI ) )
  ALLOCATE( ZZ0        ( INI ) )
  ALLOCATE( ZZ0H       ( INI ) )
  ALLOCATE( ZQS_SEA    ( INI ) )
  ALLOCATE( ZQS_WATER  ( INI ) )
  ALLOCATE( ZQS_NATURE ( INI ) )
  ALLOCATE( ZQS_TOWN   ( INI ) )
  ALLOCATE( ZQS        ( INI ) )
  ALLOCATE( ZPSNG      ( INI ) )
  ALLOCATE( ZPSNV      ( INI ) )
  ALLOCATE( ZZ0EFF     ( INI ) )
  ALLOCATE( ZZS        ( INI ) )
  !
  ISERIES = 0
  CALL GET_SURF_VAR_n(YSC%FM, YSC%IM, YSC%SM, YSC%TM, YSC%GDM, YSC%WM, YSC%DUO, YSC%DU,  YSC%UG, YSC%U, YSC%USS, &
                      CSURF_FILETYPE,INI,ISERIES,PSEA=ZSEA,PWATER=ZWATER,PNATURE=ZNATURE,PTOWN=ZTOWN, &
                      PT2M=ZT2M,PQ2M=ZQ2M,PQS=ZQS,PZ0=ZZ0,PZ0H=ZZ0H,PZ0EFF=ZZ0EFF,PQS_SEA=ZQS_SEA,  &
                      PQS_WATER=ZQS_WATER,PQS_NATURE=ZQS_NATURE,PQS_TOWN=ZQS_TOWN,                  &
                      PPSNG=ZPSNG,PPSNV=ZPSNV,PZS=ZZS                                         )
  !
  ISIZE = SIZE(NINDEX)
  IF (NRANK==NPIO) THEN
    ALLOCATE(ZSEA_FULL   (ISIZE))
    ALLOCATE(ZWATER_FULL (ISIZE))
    ALLOCATE(ZNATURE_FULL(ISIZE))
    ALLOCATE(ZTOWN_FULL  (ISIZE))
    ALLOCATE(ZZ0_FULL    (ISIZE))
    ALLOCATE(ZZ0EFF_FULL (ISIZE))
    ALLOCATE(ZZS_FULL    (ISIZE))
  ELSE
    ALLOCATE(ZSEA_FULL   (0))
    ALLOCATE(ZWATER_FULL (0))
    ALLOCATE(ZNATURE_FULL(0))
    ALLOCATE(ZTOWN_FULL  (0))
    ALLOCATE(ZZ0_FULL    (0))
    ALLOCATE(ZZ0EFF_FULL (0))
    ALLOCATE(ZZS_FULL    (0))
  ENDIF
  CALL GATHER_AND_WRITE_MPI(ZSEA,ZSEA_FULL)
  CALL GATHER_AND_WRITE_MPI(ZWATER,ZWATER_FULL)
  CALL GATHER_AND_WRITE_MPI(ZNATURE,ZNATURE_FULL)
  CALL GATHER_AND_WRITE_MPI(ZTOWN,ZTOWN_FULL)
  CALL GATHER_AND_WRITE_MPI(ZZ0,ZZ0_FULL)
  CALL GATHER_AND_WRITE_MPI(ZZ0EFF,ZZ0EFF_FULL)
  CALL GATHER_AND_WRITE_MPI(ZZS,ZZS_FULL)

  IF (NRANK==NPIO) THEN
    WRITE(ILUOUT,'(A32,I4,A3,I4)') ' GRID BOXES CONTAINING SEA    : ',COUNT( ZSEA_FULL    (:) > 0. ),' / ',ISIZE
    WRITE(ILUOUT,'(A32,I4,A3,I4)') ' GRID BOXES CONTAINING WATER  : ',COUNT( ZWATER_FULL  (:) > 0. ),' / ',ISIZE
    WRITE(ILUOUT,'(A32,I4,A3,I4)') ' GRID BOXES CONTAINING NATURE : ',COUNT( ZNATURE_FULL (:) > 0. ),' / ',ISIZE
    WRITE(ILUOUT,'(A32,I4,A3,I4)') ' GRID BOXES CONTAINING TOWN   : ',COUNT( ZTOWN_FULL   (:) > 0. ),' / ',ISIZE
    WRITE(ILUOUT,*)'ZZ0    = ',ZZ0_FULL
    WRITE(ILUOUT,*)'ZZ0EFF = ',ZZ0EFF_FULL
    WRITE(ILUOUT,*)'ZZS = ',ZZS_FULL
    WRITE(ILUOUT,*)'MINVAL(ZZS) = ',MINVAL(ZZS_FULL),' MAXVAL(ZZS) = ',MAXVAL(ZZS_FULL)
  ENDIF
  !
  DEALLOCATE( ZSEA       )
  DEALLOCATE( ZWATER     )
  DEALLOCATE( ZNATURE    )
  DEALLOCATE( ZTOWN      )
  DEALLOCATE( ZT2M       )
  DEALLOCATE( ZQ2M       )
  DEALLOCATE( ZZ0        )
  DEALLOCATE( ZZ0H       )
  DEALLOCATE( ZQS_SEA    )
  DEALLOCATE( ZQS_WATER  )
  DEALLOCATE( ZQS_NATURE )
  DEALLOCATE( ZQS_TOWN   )
  DEALLOCATE( ZQS        )
  DEALLOCATE( ZPSNG      )
  DEALLOCATE( ZPSNV      )
  DEALLOCATE( ZZ0EFF     )
  DEALLOCATE( ZZS        )
  DEALLOCATE( NTAB_SYT   )
  !
  IF (NRANK==NPIO) THEN
    DEALLOCATE(ZSEA_FULL   )
    DEALLOCATE(ZWATER_FULL )
    DEALLOCATE(ZNATURE_FULL)
    DEALLOCATE(ZTOWN_FULL  )
    DEALLOCATE(ZZ0_FULL    )
    DEALLOCATE(ZZ0EFF_FULL )
    DEALLOCATE(ZZS_FULL    )
  ENDIF
  !
ENDIF
!
! --------------------------------------------------------------------------------------
!
!    4'    Close Gelato specific diagnostic
#if ! defined in_arpege
CALL CLSDIA()
#endif
!
!
!*    5.     Close parallelized I/O
!            ----------------------
!
IF (NRANK==NPIO) THEN
  !
  WRITE(YTAG,FMT='(I4.4,A1,I2.2,A1,I2.2,A4,I2.2,A1,I2.2)') IYEAR_OUT,'-',IMONTH_OUT,'-',IDAY_OUT,&
                  ' at ',INT(ZTIME_OUT/3600.),'h',NINT(ZTIME_OUT)/60-60*INT(ZTIME_OUT/3600.)     
  !
  WRITE(ILUOUT,*) ' '
  WRITE(ILUOUT,*) '    --------------------------'
  WRITE(ILUOUT,*) '    | OFFLINE ENDS CORRECTLY |'
  WRITE(ILUOUT,*) '    --------------------------'
  WRITE(ILUOUT,'(A14,A20)')'DATE ENDING : ',ADJUSTL(ADJUSTR(YTAG))
  WRITE(ILUOUT,*) ' '
  CLOSE(ILUOUT)
  WRITE(*,*) ' '
  WRITE(*,*) '    --------------------------'
  WRITE(*,*) '    | OFFLINE ENDS CORRECTLY |'
  WRITE(*,*) '    --------------------------'
  WRITE(*,'(A14,A20)')'DATE ENDING : ',ADJUSTL(ADJUSTR(YTAG))
  WRITE(*,*) ' '
ENDIF
!
 CALL SURFEX_DEALLO_LIST
!
IF (ALLOCATED(NINDEX)) DEALLOCATE(NINDEX)
IF (ALLOCATED(NSIZE_TASK)) DEALLOCATE(NSIZE_TASK)
!
 CALL END_LOG_MPI
!
IF (LHOOK) CALL DR_HOOK('OFFLINE',1,ZHOOK_HANDLE)
!
! * MPI and OASIS must be finalized after the last DR_HOOK call
!
IF (LXIOS) THEN
#ifdef WXIOS
  CALL XIOS_CONTEXT_FINALIZE()
#endif
ENDIF
!
CALL SFX_OASIS_END
!
#ifdef SFX_MPI
IF(.NOT. LOASIS .AND. .NOT. LXIOS) THEN
  CALL MPI_FINALIZE(INFOMPI)
ENDIF
#endif
!
! --------------------------------------------------------------------------------------
!
END PROGRAM OFFLINE
