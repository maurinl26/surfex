!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!------------------------
PROGRAM CREATE_FORC
!
!------------------------
!!
!!    PURPOSE
!!    -------
!!   This program prepares the forcing files for offline 
!!   SURFEX/TOPMODEL run : it works with complete GRIB files
!!
!!    METHOD
!!    ------
!!   
!!
!!    EXTERNAL
!!    --------
!!
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!
!!    REFERENCE
!!    ---------
!!    From routines availables in SURFEX OFFLINE:
!!     - prep.f90
!!     - prep_input_experiment.f90
!!
!!    AUTHOR
!!    ------
!!
!!    B. Vincendon                 Meteo-France
!!
!!    MODIFICATION
!!    ------------
!!
!!    Original     11/2014
!!
!----------------------------------------------------------------------------
!
! Modules to read grib files and to interpolate
USE MODE_POS_SURF
USE MODD_IO_SURF_ASC,     ONLY : NUNIT, CFILEIN, CFILEOUT, NMASK, NLUOUT, NFULL, CMASK
USE MODD_SURF_CONF
USE MODD_SURF_PAR
USE MODD_SURF_ATM_n, ONLY : SURF_ATM_t
USE MODD_SURF_ATM_GRID_n, ONLY : SURF_ATM_GRID_t
USE MODD_DATA_COVER_n, ONLY : DATA_COVER_t
!
USE MODD_OFF_SURFEX_n
USE MODD_TYPE_DATE_SURF
USE MODN_PREP_ISBA
USE MODN_IO_OFFLINE
!
USE MODI_OPEN_NAMELIST
USE MODI_CLOSE_NAMELIST
USE MODI_READ_FORC_SAF
USE MODI_READ_SURF
USE MODI_PREP_OUTPUT_GRID
USE MODI_CUMUL_TO_HOURLY_RR
USE MODI_ADD_FORECAST_TO_DATE_SURF
USE MODI_INIT_FORC_GRIB
!
USE MODI_TEST_NAM_VAR_SURF 
! Modules from PREP_INPUT_EXPERIMENT
USE MODD_CSTS
USE MODD_TYPE_DATE_SURF
USE MODD_OL_FILEID,      ONLY : XNETCDF_FILENAME_OUT
!
USE MODI_CREATE_FILE
USE MODI_SUNPOS
USE MODI_WRITE_SURF
USE MODI_WRITE_NETCDF
USE MODI_DEF_VAR_NETCDF
USE MODI_READ_INTERP_QPE_RADAR
USE MODI_READ_INTERP_QPE_PLUVIOS
USE MODI_OPEN_FILE
USE MODI_CLOSE_FILE
!
USE MODI_READ_FORC_GRIB
USE MODI_GET_LUOUT
USE MODI_SELECT_MEMBER
USE MODI_WRITE_TIME
USE MODE_POS_SURF,  ONLY : POSNAM
!
USE MODD_FORC_OFFLINE_NC
! To introduce a perturbation into the rainfall field
USE MODI_READ_GRIB_LARGE_RAIN
USE MODI_HOR_INTERPOL
USE MODD_OBJ   ,         ONLY : LINAREA, NNUMREG, LRADOK
USE MODD_GRID_AROME,     ONLY : NX, NY
USE MODD_GRID_GRIB,      ONLY : NNI
USE MODD_PREP,           ONLY : CINGRID_TYPE
! 
USE MODD_PERT_RAIN
USE MODI_CHANGE_LOC_OBJECTS
USE MODI_CHANGE_AMPLI_RAIN
USE MODI_GET_SIZES_PARALLEL
USE MODD_SURFEX_OMP,     ONLY : NBLOCKTOT
!
USE YOMHOOK   ,ONLY : LHOOK, DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!------------------------------------------------------------------------------
!
IMPLICIT NONE
!
include 'netcdf.inc'
!
!*    0.     Declaration of local variables
!            ------------------------------
!
INTEGER           :: ILUOUT ! unit number of the listing file
INTEGER           :: ILUNAM ! unit number of the namelist file
INTEGER           :: IYEAR, IMONTH, IDAY !date
REAL              :: ZTIME               !time
LOGICAL           :: GFOUND
!
CHARACTER(LEN=28) :: YLUOUT    ='LISTING_CREATE_FORC         '  ! name of listing
CHARACTER(LEN=28) :: YNAMELIST ='OPTIONS.nam                 '  ! name of the entire namelist
!
INTEGER           :: JSTP ! loop index for time steps
INTEGER           :: JWRK, JSEQ ! loop index 
INTEGER           :: JWRK1, JWRK2, JMB ! loop indexes
INTEGER           :: JK, JJ, JN ! variable to change order of radar data!*      
!
REAL              :: XLATCEN  ! latitude  of center point
REAL              :: XLONCEN  ! longitude of center point
INTEGER           :: NIMAX    ! number of points in I direction
INTEGER           :: NJMAX    ! number of points in J direction
REAL              :: XDX      ! increment in X direction (in meters)
REAL              :: XDY      ! increment in Y direction (in meters)
INTEGER           :: INI      ! number of grid points
!
NAMELIST/NAM_CONF_PROJ_GRID/NIMAX,NJMAX,XLATCEN,XLONCEN,XDX,XDY
!
!
!     0.1 Declaration of namelist NAM_FORC_OFFLINE_NC and related variables
!
NAMELIST/NAM_FORC_OFFLINE_NC/NNB_FORC_STP, NNB_FORC_SEQUENCES,&
                             NSTP_BEG, NSTP_END, CTYPE_SEQUENCES,&
                             CGRIB_BASE_NAME, CRAD_BASE_NAME,&
                             CGRIB_TYPE,&
                             XTA_IDEA, XQA_IDEA, XDIRSW_IDEA, XSCASW_IDEA,&
                             XLW_IDEA, XPS_IDEA, XRAIN_IDEA, XSNOW_IDEA,&
                             XWINDSPEED_IDEA, XWINDDIR_IDEA,&
                             NNB_MEMBERS_ENS, NNB_STEPS_MODIF
!                             
INTEGER                         :: INB_FORC_STP       ! Number of forcing time steps
!                                                     where model's gribs are used
INTEGER                         :: INB_FORC_SEQUENCES ! Number (up to 5) of sequences
!                                                     with different sources in
!                                                     the forcing serie
INTEGER,           DIMENSION(:), ALLOCATABLE :: IBEG   ! Time step of a sequence beginning
INTEGER,           DIMENSION(:), ALLOCATABLE :: IEND   ! Time step of a sequence end
CHARACTER(LEN=6),  DIMENSION(:), ALLOCATABLE :: YTYPE_SEQUENCES! Type of forcing sequences
!                                                     can be 'MODEL ','NORAIN','RADAR ',&
!                                                            'SAFRAN','IDEA  'and 'PERTRR'
CHARACTER(LEN=18), DIMENSION(15) :: YGRIB_BASE_NAME    ! Base name of grib file
CHARACTER(LEN=18), DIMENSION(15) :: YRAD_BASE_NAME     ! Base name of radar file
CHARACTER(LEN=6)                 :: YGRIB_TYPE         ! Model source
REAL                             :: ZSTEP_ISBA         ! ISBA time step
REAL                             :: ZSTEP_FORCING      ! forcing time step
!
INTEGER                          :: INB_MEMBERS
!
INTEGER                          :: INB_HOUR_SAF ! number of hours for safran data
INTEGER                          :: INB_STEP_MOD ! number of total number of model data
INTEGER                          :: INB_STEP_RAD ! number of total number of radar data
CHARACTER(LEN=28)                :: YFILE_GRB ! name of grib file to be read
CHARACTER(LEN=29)                :: YFILE_RAD ! Name of the radar data file

CHARACTER(LEN=4)                 :: YSTEP
INTEGER                          :: IUNIT, ITMP
TYPE (DATE_TIME)                 :: TZTIME_GRIB    ! current date and time
CHARACTER(LEN=100)               :: YCOMMENT ! comment
!
!     0.2 Work 1D variables
!
! 1D tables (NB_points)
REAL, DIMENSION(:), ALLOCATABLE :: ZWRK_T2M 
REAL, DIMENSION(:), ALLOCATABLE :: ZWRK_Q2M
REAL, DIMENSION(:), ALLOCATABLE :: ZWRK_RADSW
REAL, DIMENSION(:), ALLOCATABLE :: ZWRK_RADLW
REAL, DIMENSION(:), ALLOCATABLE :: ZWRK_U10M
REAL, DIMENSION(:), ALLOCATABLE :: ZWRK_V10M
REAL, DIMENSION(:), ALLOCATABLE :: ZWRK_FF10M
REAL, DIMENSION(:), ALLOCATABLE :: ZWRK_PS
REAL, DIMENSION(:), ALLOCATABLE :: ZWRK_RR1
REAL, DIMENSION(:), ALLOCATABLE :: ZWRK_RS1
REAL, DIMENSION(:), ALLOCATABLE :: ZRAIN_TMP
!
!
!     0.3 Reading input files
!
CHARACTER(LEN=6)                 :: YPROG = 'OFFLIN'
!
TYPE (DATE_TIME)                 :: TDTCUR
!
INTEGER                          :: IFILE_ID
INTEGER                          :: IVAR_ID
INTEGER                          :: IRES
INTEGER                          :: IRET
INTEGER,            DIMENSION(6) :: IDDIM
INTEGER,            DIMENSION(2) :: KDIMS
CHARACTER(LEN=100) ,DIMENSION(6) :: YNAME_DIM
CHARACTER(LEN=100) ,DIMENSION(2) :: YATT_TITLE,YATT,YATB
!
REAL,                DIMENSION(:), ALLOCATABLE :: ZTM
!
REAL                             :: ZSTEP_OUTPUT        ! output time step (s)
REAL                             :: ZTMP, ZDEN       ! trash variable
INTEGER                          :: IDATE_DEB, ISEC
INTEGER, DIMENSION(3)            :: ITIME
LOGICAL :: GCUMUL    ! true if file 
!                      contains cumulated values of rain, radiations
!                      then a file is indexed 00 to permit hourly values
!                      computation
! In case the user wants to write at the screen the cumulated values
LOGICAL,                          PARAMETER :: GWRITE_CUMULATED_RAIN=.FALSE.
LOGICAL :: LCPL_ESM!
! In case SAFRAN data are read from an text file
LOGICAL,                          PARAMETER :: GSAF_TXT=.FALSE.
CHARACTER(LEN=18),                PARAMETER :: YFILE_FRC_SAFRAN='forcing_sept2002.txt'
!
!     0.4 Writing NETCDF forcing file
!
! In case the user wants to write rainfall into a llv text file
LOGICAL,                          PARAMETER  :: GWRITE_LLV=.FALSE.
!
CHARACTER(LEN=10), PARAMETER                 :: YFILE_FORCING_OUT = 'FORCING.nc'
CHARACTER(LEN=14), DIMENSION(:), ALLOCATABLE :: YFILE_FORCING_NEW !Netcdf forcing files
CHARACTER(LEN=14)                            :: YFILE_ASC         ! Ascii output file
!
!  Variables useful in FORCING file
! 2D tables (NB_points,NB_time_steps)
REAL, DIMENSION(:,:),  ALLOCATABLE :: ZU10M     ! U component of wind
REAL, DIMENSION(:,:),  ALLOCATABLE :: ZV10M     ! V component of wind
REAL, DIMENSION(:,:),  ALLOCATABLE :: ZRR1      ! Accumulated precipitations
REAL, DIMENSION(:,:),  ALLOCATABLE :: ZRS1      ! Accumulated snowfall (snow+graupel)
REAL, DIMENSION(:,:),  ALLOCATABLE :: ZCO2      ! CO2 concentration (kg/m3) 
REAL, DIMENSION(:,:),  ALLOCATABLE :: ZDIR_SW   ! Solar direct   radiation (W/m2)
REAL, DIMENSION(:,:),  ALLOCATABLE :: ZSCA_SW   ! Solar diffused radiation (W/m2)
REAL, DIMENSION(:,:),  ALLOCATABLE :: ZLW       ! Longwave radiation (W/m2)
REAL, DIMENSION(:,:),  ALLOCATABLE :: ZDIR_SW_1H! Solar direct   radiation (W/m2 anytime) 
REAL, DIMENSION(:,:),  ALLOCATABLE :: ZLW_1H    ! Longwave radiation (W/m2 anytime )
REAL, DIMENSION(:,:),  ALLOCATABLE :: ZWINDSPEED! Wind speed (m/s)
REAL, DIMENSION(:,:),  ALLOCATABLE :: ZWINDDIR  ! Wind dir. (deg. from N, clockwise)
REAL, DIMENSION(:,:),  ALLOCATABLE :: ZRAIN     ! rain rate (kg/m2/s)
REAL, DIMENSION(:,:),  ALLOCATABLE :: ZSNOW     ! snow rate (kg/m2/s)
REAL, DIMENSION(:,:),  ALLOCATABLE :: ZTA       ! temperature (K)
REAL, DIMENSION(:,:),  ALLOCATABLE :: ZQA       ! humidity (kg/m3)
REAL, DIMENSION(:,:),  ALLOCATABLE :: ZPS       ! pressure (Pa)
LOGICAL                            :: GINTERP = .TRUE.
! 1D tables (NB_points)
REAL, DIMENSION(:),    ALLOCATABLE :: ZZREF     ! height of temperature forcing (m)
REAL, DIMENSION(:),    ALLOCATABLE :: ZUREF     ! height of wind forcing (m)
REAL, DIMENSION(:),    ALLOCATABLE :: ZZS       ! orography (m)
REAL, DIMENSION(:),    ALLOCATABLE :: ZLON      ! longitude (degrees)
REAL, DIMENSION(:),    ALLOCATABLE :: ZLAT      ! latitude  (degrees)
REAL, DIMENSION(:),    ALLOCATABLE :: ZRR1_0,ZRS1_0,ZDIR_SW_0,ZLW_0
!
!Safran readed variables
LOGICAL                            :: GSAFRAN = .TRUE.
CHARACTER(LEN=14)                  :: YDATE
REAL, DIMENSION(:,:),  ALLOCATABLE :: ZSAF_T2M,ZSAF_Q2M,ZSAF_RADSW,ZSAF_RADLW
REAL, DIMENSION(:,:),  ALLOCATABLE :: ZSAF_FF10M,ZSAF_PS,ZSAF_RR1,ZRAD_RR1
!
! Variables for ideal cases
REAL,  DIMENSION(15)  :: ZTA_IDEA, ZQA_IDEA, ZDIRSW_IDEA, ZSCASW_IDEA
REAL,  DIMENSION(15)  :: ZLW_IDEA, ZPS_IDEA, ZRAIN_IDEA, ZSNOW_IDEA
REAL,  DIMENSION(15)  :: ZWINDSPEED_IDEA, ZWINDDIR_IDEA
!
! Variables to generate ensemble of rainfall introducing perturbation                        
LOGICAL                                         :: GPERTRAIN !true if one sequence has to be perturbated
CHARACTER(LEN=4)                                :: YMEMBER
REAL,             DIMENSION(:,:,:), ALLOCATABLE :: ZRAIN_SEL
REAL,   POINTER,  DIMENSION(:,:)                :: ZPNT_INTERPOL
REAL,             DIMENSION(:,:),   ALLOCATABLE :: ZFLD_INTERPOL !Field interpolated
REAL,             DIMENSION(:,:),   ALLOCATABLE :: ZRAIN_IN, ZRAIN_OUT
!
!
REAL,             DIMENSION(:,:),   ALLOCATABLE :: ZLARGE_RAIN     ! rain on large domain
REAL,             DIMENSION(:),     ALLOCATABLE :: ZLARGE_RAIN_0   ! at 0
REAL,             DIMENSION(:,:),   ALLOCATABLE :: ZLARGE_RR1      ! hourly rain on large domain
REAL,             DIMENSION(:,:),   ALLOCATABLE :: ZLARGE_RAIN_NEW ! new hourly rain on large domain
REAL,             DIMENSION(:),     ALLOCATABLE :: ZCUMRAIN 
REAL                                            :: ZCUM 
INTEGER,          DIMENSION(:),     ALLOCATABLE :: ISIZE_OMP
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
IF (LHOOK) CALL DR_HOOK('CREATE_FORC',0,ZHOOK_HANDLE)
!------------------------------------------------------------------------------
!
!*    1.      Set default names and parallelized I/O
!             --------------------------------------
!
!     Allocations of Surfex Types
CALL SURFEX_ALLOC_LIST(1)
YSC => YSURF_LIST(1)
!
!     1.1     initializations
!             ---------------
!
CFILEIN = ADJUSTL(ADJUSTR(CPGDFILE)//'.txt') ! output of PGD program
IYEAR   = NUNDEF
IMONTH  = NUNDEF
IDAY    = NUNDEF
ZTIME   = XUNDEF
!
LCPL_ESM = .FALSE.
!
!     1.2     output listing
!             --------------
CALL GET_LUOUT('ASCII ',ILUOUT)
OPEN (UNIT=ILUOUT,FILE=ADJUSTL(ADJUSTR(YLUOUT)//'.txt'),FORM='FORMATTED',ACTION='WRITE')
!
!     1.3     output file name read in namelist
!             ---------------------------------
! Default values
INB_FORC_STP       = 1 
INB_FORC_SEQUENCES = 1 
!
ZSTEP_ISBA      = 3600.
ZSTEP_FORCING   = 3600. 
CGRIB_BASE_NAME = 'GRIBIN'
CRAD_BASE_NAME  = 'RADAR'
CGRIB_TYPE      = 'AROME'
!
XTA_IDEA(:)        = 291.00
XQA_IDEA(:)        = 0.013
XDIRSW_IDEA(:)     = 200.00
XSCASW_IDEA(:)     = 0.
XLW_IDEA(:)        = 400.0
XPS_IDEA(:)        = 1015.0
XRAIN_IDEA(:)      = 0.
XSNOW_IDEA(:)      = 0.
XWINDSPEED_IDEA(:) = 1.0
XWINDDIR_IDEA(:)   = 0.
! 
INB_MEMBERS        = 50
NNB_STEPS_MODIF    = 0
!
CALL OPEN_NAMELIST('ASCII ',ILUNAM,CNAMELIST)
!
!Reading times steps duration
CALL POSNAM(ILUNAM,'NAM_IO_OFFLINE',GFOUND,ILUOUT)
IF (GFOUND) READ(UNIT=ILUNAM,NML=NAM_IO_OFFLINE)
!
ZSTEP_ISBA    = XTSTEP_SURF
ZSTEP_FORCING = XTSTEP_OUTPUT
!
!Reading space dimensions
CALL POSNAM(ILUNAM,'NAM_CONF_PROJ_GRID',GFOUND,ILUOUT)
IF (GFOUND) READ(UNIT=ILUNAM,NML=NAM_CONF_PROJ_GRID)
!
!Reading date
CALL POSNAM(ILUNAM,'NAM_PREP_ISBA',GFOUND,ILUOUT)
IF (GFOUND) READ (UNIT=ILUNAM,NML=NAM_PREP_ISBA)
!
TDTCUR%TDATE%YEAR  = NYEAR
TDTCUR%TDATE%MONTH = NMONTH
TDTCUR%TDATE%DAY   = NDAY
TDTCUR%TIME        = XTIME
IYEAR=NYEAR
IMONTH=NMONTH
IDAY=NDAY
ZTIME=XTIME
!
!Reading informations to force offline coupling system
!
CALL POSNAM(ILUNAM,'NAM_FORC_OFFLINE_NC',GFOUND,ILUOUT)
IF (GFOUND) READ (UNIT=ILUNAM,NML=NAM_FORC_OFFLINE_NC)
!
INB_FORC_STP       = NNB_FORC_STP
INB_FORC_SEQUENCES = NNB_FORC_SEQUENCES
INB_MEMBERS        = NNB_MEMBERS_ENS
write(*,*) 'NNB_STEPS_MODIF:',NNB_STEPS_MODIF
!
ALLOCATE(YTYPE_SEQUENCES(INB_FORC_SEQUENCES))
ALLOCATE(IBEG(INB_FORC_SEQUENCES))
ALLOCATE(IEND(INB_FORC_SEQUENCES))
!
YTYPE_SEQUENCES(:) = 'SAFRAN'
IBEG(:)            = 1 
IEND(:)            = INB_FORC_STP 
YTYPE_SEQUENCES(1:INB_FORC_SEQUENCES) = CTYPE_SEQUENCES(1:INB_FORC_SEQUENCES)
!
DO JSEQ=1,INB_FORC_SEQUENCES
  CALL TEST_NAM_VAR_SURF(ILUOUT,'CTYPE_SEQUENCES',CTYPE_SEQUENCES(JSEQ),&
                      'SAFRAN','RADAR ','MODEL ','NORAIN','PERTRR','MODRAD','IDEA  ')
  IBEG(JSEQ) = NSTP_BEG(JSEQ)
  IEND(JSEQ) = NSTP_END(JSEQ) 
ENDDO
!
YGRIB_BASE_NAME(1:INB_FORC_SEQUENCES) = CGRIB_BASE_NAME(1:INB_FORC_SEQUENCES)
YRAD_BASE_NAME(1:INB_FORC_SEQUENCES)  = CRAD_BASE_NAME(1:INB_FORC_SEQUENCES)
!
!
ZTA_IDEA(1:INB_FORC_SEQUENCES)       = XTA_IDEA(1:INB_FORC_SEQUENCES)
ZQA_IDEA(1:INB_FORC_SEQUENCES)       = XQA_IDEA(1:INB_FORC_SEQUENCES)
ZDIRSW_IDEA(1:INB_FORC_SEQUENCES)    = XDIRSW_IDEA(1:INB_FORC_SEQUENCES)
ZSCASW_IDEA(1:INB_FORC_SEQUENCES)    = XSCASW_IDEA(1:INB_FORC_SEQUENCES)
ZLW_IDEA(1:INB_FORC_SEQUENCES)       = XLW_IDEA(1:INB_FORC_SEQUENCES)
ZPS_IDEA(1:INB_FORC_SEQUENCES)       = XPS_IDEA(1:INB_FORC_SEQUENCES)
ZRAIN_IDEA(1:INB_FORC_SEQUENCES)     = XRAIN_IDEA(1:INB_FORC_SEQUENCES)
ZSNOW_IDEA(1:INB_FORC_SEQUENCES)     = XSNOW_IDEA(1:INB_FORC_SEQUENCES)
ZWINDSPEED_IDEA(1:INB_FORC_SEQUENCES)= XWINDSPEED_IDEA(1:INB_FORC_SEQUENCES)
ZWINDDIR_IDEA(1:INB_FORC_SEQUENCES)  = XWINDDIR_IDEA(1:INB_FORC_SEQUENCES)
!
CALL CLOSE_NAMELIST('ASCII ',ILUNAM)
!
DO JSEQ=1,INB_FORC_SEQUENCES
  IF (IEND(JSEQ)>INB_FORC_STP) THEN
    WRITE(*,*) 'Sequence ',JSEQ,'inconstistancy INB_FORC_STP vs IEND : ',INB_FORC_STP,IEND(JSEQ)
    stop
  ENDIF
  IF (IBEG(JSEQ)>IEND(JSEQ)) THEN
    WRITE(*,*) 'Sequence ',JSEQ,'inconstistancy IBEG vs IEND : ',IBEG(JSEQ),IEND(JSEQ)
    stop
  ENDIF
ENDDO
!
WRITE(ILUOUT,*) 'Creating FORCING for SURFEX V8'
write(*,*)      'Creating FORCING for SURFEX V8'
!
WRITE(ILUOUT,*) 'Number of sequences :',INB_FORC_SEQUENCES
write(*,*)      'Number of sequences :',INB_FORC_SEQUENCES
!
DO JSEQ=1,INB_FORC_SEQUENCES
  WRITE(ILUOUT,*)'SEQUENCE :',JSEQ,'TYPE : ',YTYPE_SEQUENCES(JSEQ)
  write(*,*)     'SEQUENCE :',JSEQ,'TYPE : ',YTYPE_SEQUENCES(JSEQ)
  WRITE(ILUOUT,*)'Beginning at step:',IBEG(JSEQ),'/ End at step:',IEND(JSEQ)
  write(*,*)     'Beginning at step:',IBEG(JSEQ),'/ End at step:',IEND(JSEQ)
ENDDO
!----------------------------------------------------------------------------
!
!*    2.      Preparation of surface physiographic fields
!             -------------------------------------------
CALL INIT_FORC_GRIB(YSC,'ASCII ')
write(*,*) 'init forc grib'
!
INI=YSC%U%NDIM_FULL
ALLOCATE(ZZS(INI))
ZZS=YSC%U%XZS
!
CALL PREP_OUTPUT_GRID(YSC%UG%G, YSC%UG%G, YSC%U%NSIZE_FULL, ILUOUT)
!
write(*,*) 'preparation ok'
!----------------------------------------------------------------------------
!
!*    3.      Reading variables and interpolation
!             -------------------------------------------
! Allocation of variables and initialisation

ALLOCATE (ZU10M      (INB_FORC_STP,INI))
ALLOCATE (ZV10M      (INB_FORC_STP,INI))
ALLOCATE (ZTA        (INB_FORC_STP,INI))
ALLOCATE (ZQA        (INB_FORC_STP,INI))
ALLOCATE (ZPS        (INB_FORC_STP,INI))
ALLOCATE (ZWINDSPEED (INB_FORC_STP,INI))
ALLOCATE (ZWINDDIR   (INB_FORC_STP,INI))
ALLOCATE (ZSCA_SW    (INB_FORC_STP,INI))
ALLOCATE (ZRAIN      (INB_FORC_STP,INI))
ALLOCATE (ZSNOW      (INB_FORC_STP,INI))
ALLOCATE (ZCO2       (INB_FORC_STP,INI))
ALLOCATE (ZDIR_SW_1H (INB_FORC_STP,INI))
ALLOCATE (ZLW_1H     (INB_FORC_STP,INI))
ALLOCATE (ZRR1       (INB_FORC_STP,INI))
ALLOCATE (ZRS1       (INB_FORC_STP,INI))
ALLOCATE (ZDIR_SW    (INB_FORC_STP,INI))
ALLOCATE (ZLW        (INB_FORC_STP,INI))
!
ZU10M(:,:)      = XUNDEF
ZV10M(:,:)      = XUNDEF
ZRR1(:,:)       = 0.
ZRS1(:,:)       = 0.
ZTA(:,:)        = XUNDEF
ZQA(:,:)        = XUNDEF
ZDIR_SW(:,:)    = XUNDEF
ZSCA_SW(:,:)    = 0.
ZLW(:,:)        = XUNDEF
ZWINDSPEED(:,:) = XUNDEF
ZWINDDIR(:,:)   = XUNDEF
ZPS(:,:)        = XUNDEF
ZRAIN(:,:)      = XUNDEF
ZSNOW(:,:)      = XUNDEF
ZCO2(:,:)       = 0.000620   ! (kg/m3, equivalent to 350 ppm) 
!
write(*,*) 'allocations ok',SIZE(ZTA,1),SIZE(ZTA,2)
!    3.1     Reading variables from simulations(grib)
!             -------------------------------------------
!
DO JSEQ=1,INB_FORC_SEQUENCES
  WRITE(ILUOUT,*) 'Sequence',JSEQ,'=',YTYPE_SEQUENCES(JSEQ)
  write(*,*)      'Sequence',JSEQ,'=',YTYPE_SEQUENCES(JSEQ)
  !
  ! Determining  what is the beginning date
  IF (JSEQ==1) THEN
    ZTIME=ZTIME+IBEG(JSEQ)*ZSTEP_FORCING
  ELSE
    ZTIME=ZTIME+(IBEG(JSEQ)-IBEG(JSEQ-1))*ZSTEP_FORCING
  ENDIF
  !
  IF (86400.-ZTIME < 1.E-6) THEN
    CALL ADD_FORECAST_TO_DATE_SURF(IYEAR,IMONTH,IDAY,ZTIME)
  ENDIF
  !
  IDATE_DEB=(IYEAR*1000000)+(IMONTH*10000)+(IDAY*100)+FLOOR(ZTIME/3600.)
  write(*,*) 'IDATE_DEB =',IDATE_DEB
  !
  IF (YTYPE_SEQUENCES(JSEQ)=='IDEA  ') THEN  
    !Idealized cases
    ZTA(IBEG(JSEQ):IEND(JSEQ),:)        = ZTA_IDEA(JSEQ)
    ZQA(IBEG(JSEQ):IEND(JSEQ),:)        = ZQA_IDEA(JSEQ)
    ZDIR_SW_1H(IBEG(JSEQ):IEND(JSEQ),:) = ZDIRSW_IDEA(JSEQ)
    ZSCA_SW(IBEG(JSEQ):IEND(JSEQ),:)    = ZSCASW_IDEA(JSEQ)
    ZLW_1H(IBEG(JSEQ):IEND(JSEQ),:)     = ZLW_IDEA(JSEQ)
    ZPS(IBEG(JSEQ):IEND(JSEQ),:)        = ZPS_IDEA(JSEQ)
    ZRAIN(IBEG(JSEQ):IEND(JSEQ),:)      = ZRAIN_IDEA(JSEQ)
    ZSNOW(IBEG(JSEQ):IEND(JSEQ),:)      = ZSNOW_IDEA(JSEQ)
    ZWINDSPEED(IBEG(JSEQ):IEND(JSEQ),:) = ZWINDSPEED_IDEA(JSEQ)
    ZWINDDIR(IBEG(JSEQ):IEND(JSEQ),:)   = ZWINDDIR_IDEA(JSEQ)
  ELSE IF((YTYPE_SEQUENCES(JSEQ)=='MODEL ').OR.& 
        (YTYPE_SEQUENCES(JSEQ)=='PERTRR').OR.&
        (YTYPE_SEQUENCES(JSEQ)=='MODRAD')) THEN
  !Cases with grib files to read
    YGRIB_TYPE=CGRIB_TYPE(JSEQ)
    INB_STEP_MOD=IEND(JSEQ)-IBEG(JSEQ)+1
    ALLOCATE(ZWRK_T2M   (INI))  
    ALLOCATE(ZWRK_Q2M   (INI))  
    ALLOCATE(ZWRK_RADSW (INI))  
    ALLOCATE(ZWRK_RADLW (INI))  
    ALLOCATE(ZWRK_U10M  (INI))  
    ALLOCATE(ZWRK_V10M  (INI))  
    ALLOCATE(ZWRK_PS    (INI))  
    ALLOCATE(ZWRK_RR1   (INI)) 
    ALLOCATE(ZWRK_RS1   (INI)) 
    !
    DO JSTP=IBEG(JSEQ),IEND(JSEQ)
      JWRK1=JSTP-(IBEG(JSEQ)-1)
      IF (JWRK1<10) THEN
        WRITE(YSTEP,'(I1)') JWRK1
        YSTEP='0'//YSTEP
      ELSEIF (JWRK1 < 100) THEN
        WRITE(YSTEP,'(I2)') JWRK1
      ELSE
        WRITE(*,*) 'You cannot use more than 99 hours of model hourly forecast'
        STOP
      ENDIF
      !
      YFILE_GRB=TRIM(YGRIB_BASE_NAME(JSEQ))//'_'//TRIM(YSTEP)//'.grb'
      write(*,*) YFILE_GRB
      CALL READ_FORC_GRIB(YSC,YFILE_GRB,YGRIB_TYPE,ZWRK_T2M(:),ZWRK_Q2M(:),&
                       ZWRK_RADSW(:),ZWRK_RADLW(:),&
                       ZWRK_U10M(:),ZWRK_V10M(:),&
                       ZWRK_PS(:),ZWRK_RR1(:),ZWRK_RS1(:),GINTERP)
      write(*,*) 'RAIN from grib : ',&
          MINVAL(ZWRK_RR1),MAXVAL(ZWRK_RR1,MASK=ZWRK_RR1/=XUNDEF)
      !
      ZTA(JSTP,:)     = ZWRK_T2M(:)
      ZQA(JSTP,:)     = ZWRK_Q2M(:)
      ZDIR_SW(JSTP,:) = ZWRK_RADSW(:)
      ZLW(JSTP,:)     = ZWRK_RADLW(:)
      ZU10M(JSTP,:)   = ZWRK_U10M(:)
      ZV10M(JSTP,:)   = ZWRK_V10M(:)
      ZPS(JSTP,:)     = ZWRK_PS(:)
      ZRR1(JSTP,:)    = ZWRK_RR1(:)
      ZRS1(JSTP,:)    = ZWRK_RS1(:)
      !
      IF (YTYPE_SEQUENCES(JSEQ)=='PERTRR') THEN
        IF (.NOT.ALLOCATED(ZLARGE_RAIN)) ALLOCATE (ZLARGE_RAIN(INB_STEP_MOD,NNI))
        CALL READ_GRIB_LARGE_RAIN(YFILE_GRB,YGRIB_TYPE,ZLARGE_RAIN(JWRK1,:))
      ENDIF
    ENDDO!JSTP

! Gestion des cas où on ne débute pas au premier pas de temps AROME
! La première pluie lue est cumulée depuis le début de la simulation
! Donc on est obligés de récupérer le fichier AROME précédent
! pour obtenir une pluie horaire
    YFILE_GRB=TRIM(YGRIB_BASE_NAME(JSEQ))//'_00.grb'
    INQUIRE(FILE=YFILE_GRB, EXIST=GCUMUL)
    IF (GCUMUL) THEN
      ALLOCATE (ZRR1_0(INI))
      ALLOCATE (ZRS1_0(INI))
      ALLOCATE(ZDIR_SW_0(INI))
      ALLOCATE(ZLW_0(INI))
      YFILE_GRB=TRIM(YGRIB_BASE_NAME(JSEQ))//'_00.grb'
      write(*,*) YFILE_GRB
      ZRR1_0(:)=0.0
      ZRS1_0(:)=0.0
      ZDIR_SW_0(:)=0.0
      ZLW_0(:)=0.0
      CALL READ_FORC_GRIB(YSC,YFILE_GRB,YGRIB_TYPE,ZWRK_T2M(:),ZWRK_Q2M(:),&
                          ZWRK_RADSW(:),ZWRK_RADLW(:),&
                          ZWRK_U10M(:),ZWRK_V10M(:),&
                          ZWRK_PS(:),ZWRK_RR1(:),ZWRK_RS1(:),GINTERP)
      ZRR1_0(:)=ZWRK_RR1(:)
      ZRS1_0(:)=ZWRK_RS1(:)
      ZDIR_SW_0(:)=ZWRK_RADSW(:)
      ZLW_0(:)=ZWRK_RADLW(:)
      write(*,*) 'RAIN from grib : ',&
          MINVAL(ZWRK_RR1),MAXVAL(ZWRK_RR1,MASK=ZWRK_RR1/=XUNDEF)
      !
      IF (YTYPE_SEQUENCES(JSEQ)=='PERTRR') THEN
        IF (.NOT.ALLOCATED(ZLARGE_RAIN_0)) ALLOCATE(ZLARGE_RAIN_0(NNI))
        ZLARGE_RAIN_0(:)=0.
        CALL READ_GRIB_LARGE_RAIN(YFILE_GRB,YGRIB_TYPE,ZLARGE_RAIN_0(:))
      ENDIF
    ENDIF!GCUMUL
    !
    DEALLOCATE(ZWRK_T2M)  
    DEALLOCATE(ZWRK_Q2M)  
    DEALLOCATE(ZWRK_RADSW)  
    DEALLOCATE(ZWRK_RADLW)  
    DEALLOCATE(ZWRK_U10M)  
    DEALLOCATE(ZWRK_V10M)  
    DEALLOCATE(ZWRK_PS)  
    DEALLOCATE(ZWRK_RR1)
    DEALLOCATE(ZWRK_RS1)
    !
    ! Values of model rainfall are cumulated from the beginning of the
    ! simulation so we decumulate them.
    ! Moreover they are in mm so we pass them in kg.m-2
    CALL CUMUL_TO_HOURLY_RR(ZRR1(IBEG(JSEQ):IEND(JSEQ),:),&
                           ZRAIN(IBEG(JSEQ):IEND(JSEQ),:),&
                           INB_STEP_MOD)
    CALL CUMUL_TO_HOURLY_RR(ZRS1(IBEG(JSEQ):IEND(JSEQ),:),&
                           ZSNOW(IBEG(JSEQ):IEND(JSEQ),:),&
                           INB_STEP_MOD)
    IF (GCUMUL) THEN 
      ZRAIN(IBEG(JSEQ),:)=ZRAIN(IBEG(JSEQ),:)-(ZRR1_0(:)/3600.)
      ZSNOW(IBEG(JSEQ),:)=ZSNOW(IBEG(JSEQ),:)-(ZRS1_0(:)/3600.)
    ENDIF
    IF (CGRIB_TYPE(JSEQ)=='AROMAN') THEN
      ZRAIN(IBEG(JSEQ):IEND(JSEQ),:)=ZRR1(IBEG(JSEQ):IEND(JSEQ),:)
      ZSNOW(IBEG(JSEQ):IEND(JSEQ),:)=ZRS1(IBEG(JSEQ):IEND(JSEQ),:)
    ENDIF
    IF (CGRIB_TYPE(JSEQ)=='ARP3H ') THEN
      ZRAIN(IBEG(JSEQ):IEND(JSEQ),:)=ZRAIN(IBEG(JSEQ):IEND(JSEQ),:)/7200.
      ZSNOW(IBEG(JSEQ):IEND(JSEQ),:)=ZSNOW(IBEG(JSEQ):IEND(JSEQ),:)/7200.
    ENDIF
    IF (YTYPE_SEQUENCES(JSEQ)=='PERTRR') THEN 
      IF (.NOT.ALLOCATED(ZLARGE_RR1)) ALLOCATE (ZLARGE_RR1(INB_STEP_MOD,NNI)) 
      CALL CUMUL_TO_HOURLY_RR(ZLARGE_RAIN(:,:),ZLARGE_RR1(:,:),INB_STEP_MOD)
      DEALLOCATE(ZLARGE_RAIN)
      IF (GCUMUL) THEN
        write(*,*) MINVAL(ZLARGE_RAIN_0),MAXVAL(ZLARGE_RAIN_0)
        WHERE (ZLARGE_RAIN_0(:)/=XUNDEF .AND. ZLARGE_RR1(1,:)/=XUNDEF)
          ZLARGE_RR1(1,:)=ZLARGE_RR1(1,:)-(ZLARGE_RAIN_0(:)/3600.)
        ENDWHERE
      ENDIF
      IF (ALLOCATED(ZLARGE_RAIN_0)) DEALLOCATE(ZLARGE_RAIN_0)
    ENDIF !IF (YTYPE_SEQUENCES(JSEQ)=='PERTRR')

    ! To compute cumulated rain
    IF (GWRITE_CUMULATED_RAIN) THEN
      IF (.NOT.ALLOCATED(ZCUMRAIN))ALLOCATE (ZCUMRAIN(INI))
      ZCUMRAIN(:)=0.
      ZCUM=0.
      DO JWRK=1,INI
        DO JSTP=IBEG(JSEQ),IEND(JSEQ)
          ZCUMRAIN(JWRK)=ZCUMRAIN(JWRK)+ZRAIN(JSTP,JWRK)
        ENDDO
      ENDDO
      ZCUM = SUM(ZCUMRAIN(:),MASK=ZCUMRAIN(:)/=XUNDEF)/&
          COUNT(ZCUMRAIN(:)/=XUNDEF)
      write(*,*) 'Cumul sequence ',JSEQ,ZCUM*3600.,&
         MINVAL(ZRAIN(JSTP,:)*3600.),MAXVAL(ZRAIN(IBEG(JSEQ):IEND(JSEQ),:)*3600.,MASK=ZRAIN(IBEG(JSEQ):IEND(JSEQ),:)/=XUNDEF)
    ENDIF !(GWRITE_CUMULATED_RAIN)
    !
    IF (YTYPE_SEQUENCES(JSEQ)=='PERTRR') THEN
      WHERE (ZLARGE_RR1(:,:)<0.0000001)ZLARGE_RR1(:,:)=0.
    ENDIF
    !
    ZWINDSPEED(IBEG(JSEQ):IEND(JSEQ),:) = SQRT(ZV10M(IBEG(JSEQ):IEND(JSEQ),:)**2+&
                                            ZU10M(IBEG(JSEQ):IEND(JSEQ),:)**2)
    WHERE (ZV10M(IBEG(JSEQ):IEND(JSEQ),:)/=0.)                                            
      ZWINDDIR  (IBEG(JSEQ):IEND(JSEQ),:) = ATAN(ZU10M(IBEG(JSEQ):IEND(JSEQ),:)/ZV10M(IBEG(JSEQ):IEND(JSEQ),:))&
                                    *180./XPI
    ELSEWHERE ( ZU10M(IBEG(JSEQ):IEND(JSEQ),:)>0.)
      ZWINDDIR  (IBEG(JSEQ):IEND(JSEQ),:) = 90.
    ELSEWHERE 
      ZWINDDIR  (IBEG(JSEQ):IEND(JSEQ),:) = 270.
    ENDWHERE
    !
    IF ((CGRIB_TYPE(JSEQ)=='AROME ').OR.(CGRIB_TYPE(JSEQ)=='PEAROM').OR.(CGRIB_TYPE(JSEQ)=='PEAROP')&
    .OR.(CGRIB_TYPE(JSEQ)=='ARPEGE').OR.(CGRIB_TYPE(JSEQ)=='MESONH').OR.(CGRIB_TYPE(JSEQ)=='ARP3H ')) THEN
         ! Values of AROME radiations are cumulated from the beginning of the
         ! simulation so we decumulate them
         ! Moreover they are in J.m-2 so wa pass them in W.m-2
      ZDIR_SW_1H(:,:) = 0.
      ZLW_1H(:,:)     = 0.
      CALL CUMUL_TO_HOURLY_RR(ZDIR_SW(IBEG(JSEQ):IEND(JSEQ),:),&
                             ZDIR_SW_1H(IBEG(JSEQ):IEND(JSEQ),:),&
                             INB_STEP_MOD)
      CALL CUMUL_TO_HOURLY_RR(ZLW(IBEG(JSEQ):IEND(JSEQ),:),&
                             ZLW_1H(IBEG(JSEQ):IEND(JSEQ),:),&
                             INB_STEP_MOD)
      IF (GCUMUL) THEN
        ZDIR_SW_1H(IBEG(JSEQ),:) = ZDIR_SW_1H(IBEG(JSEQ),:)-(ZDIR_SW_0(:)/3600.)
        ZLW_1H(IBEG(JSEQ),:)     = ZLW_1H(IBEG(JSEQ),:)-(ZLW_0(:)/3600.)
      ENDIF  
      IF (CGRIB_TYPE(JSEQ)=='ARP3H ') THEN   
        ZDIR_SW_1H(IBEG(JSEQ):IEND(JSEQ),:) = ZDIR_SW_1H(IBEG(JSEQ):IEND(JSEQ),:)/7200.
        ZLW_1H(IBEG(JSEQ):IEND(JSEQ),:)     = ZLW_1H(IBEG(JSEQ):IEND(JSEQ),:)/7200. 
      ENDIF
    ELSE
      ZDIR_SW_1H(IBEG(JSEQ):IEND(JSEQ),:) = ZDIR_SW(IBEG(JSEQ):IEND(JSEQ),:)/3600.
      ZLW_1H(IBEG(JSEQ):IEND(JSEQ),:)     = ZLW(IBEG(JSEQ):IEND(JSEQ),:)/3600.  
    ENDIF !(GRIB_TYPE == ...) 
  ELSE ! Generally, safran data are needed (else default values are used)
    !
    INQUIRE(FILE='datafile', EXIST=GSAFRAN)
    write(*,*) 'GSAFRAN',GSAFRAN
    INB_HOUR_SAF = FLOOR((IEND(JSEQ)-IBEG(JSEQ)+1)*ZSTEP_FORCING/3600.) 
    write(*,*) 'INB_HOUR_SAF =',INB_HOUR_SAF
    !
    ZTA(IBEG(JSEQ):IEND(JSEQ),:)        = 291.
    ZQA(IBEG(JSEQ):IEND(JSEQ),:)        = 0.013
    ZDIR_SW_1H(IBEG(JSEQ):IEND(JSEQ),:) = 200.
    ZSCA_SW(IBEG(JSEQ):IEND(JSEQ),:)    = 0.
    ZLW_1H(IBEG(JSEQ):IEND(JSEQ),:)     = 400.
    ZPS(IBEG(JSEQ):IEND(JSEQ),:)        = 101500.
    ZRAIN(IBEG(JSEQ):IEND(JSEQ),:)      = 0.
    ZSNOW(IBEG(JSEQ):IEND(JSEQ),:)      = 0.
    ZWINDSPEED(IBEG(JSEQ):IEND(JSEQ),:) = 1. 
    ZWINDDIR(IBEG(JSEQ):IEND(JSEQ),:)   = 0.
    ZCO2(IBEG(JSEQ):IEND(JSEQ),:)       = 0.000620 

    IF (GSAFRAN) THEN
      write(*,*) 'Sequence safran'
!*    3.1      Determining number of hours to read 
!             -------------------------------------------
      ALLOCATE (ZSAF_T2M  (INB_HOUR_SAF,INI))
      ALLOCATE (ZSAF_Q2M  (INB_HOUR_SAF,INI))
      ALLOCATE (ZSAF_RADSW(INB_HOUR_SAF,INI))
      ALLOCATE (ZSAF_RADLW(INB_HOUR_SAF,INI))
      ALLOCATE (ZSAF_FF10M(INB_HOUR_SAF,INI))
      ALLOCATE (ZSAF_RR1  (INB_HOUR_SAF,INI))
      ZSAF_RR1(:,:)=0.

      IF (GSAF_TXT) THEN
        CALL OPEN_FILE('ASCII ',IUNIT,YFILE_FRC_SAFRAN,'FORMATTED',HACTION='READ ')
        DO JWRK1=1,INI
          DO JSTP=1,INB_HOUR_SAF
            READ(IUNIT,*) ZSAF_RADSW(JSTP,JWRK1),ZSAF_RADLW(JSTP,JWRK1),ZSAF_T2M(JSTP,JWRK1),&
                      ZSAF_FF10M(JSTP,JWRK1),ZCO2(JSTP,JWRK1),&
                      ZPS(JSTP,JWRK1),ZSAF_Q2M(JSTP,JWRK1)
          ENDDO
        ENDDO
        CALL CLOSE_FILE('ASCII ',IUNIT) 
      ELSE
        CALL READ_FORC_SAF(YSC,IDATE_DEB,ZSAF_T2M,ZSAF_Q2M,ZSAF_RADSW,ZSAF_RADLW,&
                         ZSAF_FF10M,ZSAF_RR1,INB_HOUR_SAF)
      ENDIF!(GSAF_TXT)
      !
      IF (ZSTEP_FORCING>=3600.) THEN! cases where time step is not an hour
        DO JSTP=IBEG(JSEQ),IEND(JSEQ)
          JWRK=JSTP-IBEG(JSEQ)+1
          ZTA(JSTP,:)       = ZSAF_T2M  (JWRK*INB_HOUR_SAF/(IEND(JSEQ)-IBEG(JSEQ)+1),:)
          ZQA(JSTP,:)       = ZSAF_Q2M  (JWRK*INB_HOUR_SAF/(IEND(JSEQ)-IBEG(JSEQ)+1),:)
          ZDIR_SW_1H(JSTP,:)= ZSAF_RADSW(JWRK*INB_HOUR_SAF/(IEND(JSEQ)-IBEG(JSEQ)+1),:)
          ZLW_1H(JSTP,:)    = ZSAF_RADLW(JWRK*INB_HOUR_SAF/(IEND(JSEQ)-IBEG(JSEQ)+1),:)
          ZWINDSPEED(JSTP,:)= ZSAF_FF10M(JWRK*INB_HOUR_SAF/(IEND(JSEQ)-IBEG(JSEQ)+1),:)
          ZRAIN(JSTP,:)     = ZSAF_RR1  (JWRK*INB_HOUR_SAF/(IEND(JSEQ)-IBEG(JSEQ)+1),:)
          !
          IF ((INB_HOUR_SAF/(IEND(JSEQ)-IBEG(JSEQ)+1))>1) THEN
            DO  JWRK1=1,(INB_HOUR_SAF/(IEND(JSEQ)-IBEG(JSEQ)+1))-1
              ZRAIN(JSTP,:)=ZRAIN(JSTP,:)+&
                   ZSAF_RR1((JWRK*INB_HOUR_SAF/(IEND(JSEQ)-IBEG(JSEQ)+1))-JWRK1,:)
                   ! Cumulation on time steps
            ENDDO!JWRK1
          ENDIF
        ENDDO!JSTP
      ELSE !
        JWRK1=0
        ITMP=FLOOR(3600./ZSTEP_FORCING)-1!3
        DO JSTP=IBEG(JSEQ),IEND(JSEQ)
          JWRK=JSTP-IBEG(JSEQ)+1
          DO JWRK2=JSTP+JWRK1,JSTP+JWRK1+ITMP
            ZTA(JWRK2,:)       = ZSAF_T2M(JWRK,:)
            ZQA(JWRK2,:)       = ZSAF_Q2M(JWRK,:)
            ZDIR_SW_1H(JWRK2,:)= ZSAF_RADSW(JWRK,:)
            ZLW_1H(JWRK2,:)    = ZSAF_RADLW(JWRK,:)
            ZWINDSPEED(JWRK2,:)= ZSAF_FF10M(JWRK,:)
           ZRAIN(JWRK2,:)     = ZSAF_RR1(JWRK,:)
          ENDDO!JWRK2
          JWRK1=JWRK1+ITMP
        ENDDO!JSTP
      ENDIF!(ZSTEP_FORCING>=3600.)
      DEALLOCATE(ZSAF_T2M)  
      DEALLOCATE(ZSAF_Q2M)  
      DEALLOCATE(ZSAF_RADSW)  
      DEALLOCATE(ZSAF_RADLW)  
      DEALLOCATE(ZSAF_FF10M)  
      DEALLOCATE(ZSAF_RR1)  
    ENDIF ! GSAFRAN
  !
  ENDIF ! YTYPE_SEQUENCES??
  !
  !
  ! rain can also be taken to 0 or to radar data
  ! 
  IF (YTYPE_SEQUENCES(JSEQ)=='NORAIN') ZRAIN(IBEG(JSEQ):IEND(JSEQ),:)=0.
  !
  IF ((YTYPE_SEQUENCES(JSEQ)=='RADAR ') .OR.&
     (YTYPE_SEQUENCES(JSEQ)=='MODRAD').OR.&
     (YTYPE_SEQUENCES(JSEQ)=='WRFLLV')) THEN
    write(*,*) 'Sequence ',YTYPE_SEQUENCES(JSEQ)
    !
    ALLOCATE(ZRAIN_TMP(INI))
    ZRAIN_TMP(:)=0.
    INB_STEP_RAD=IEND(JSEQ)-IBEG(JSEQ)+1
    !
    ALLOCATE (ZRAD_RR1(INB_STEP_RAD,INI))
    ZRAD_RR1(:,:)=0.
    DO JSTP=IBEG(JSEQ),IEND(JSEQ)
      JWRK1=JSTP-(IBEG(JSEQ)-1)
      IF (JWRK1<10) THEN
        WRITE(YSTEP,'(I1)') JWRK1
        YSTEP='0'//YSTEP
      ELSEIF (JWRK1 < 100) THEN
        WRITE(YSTEP,'(I2)') JWRK1
      ELSEIF (JWRK1 < 1000) THEN
        WRITE(YSTEP,'(I3)') JWRK1
      ELSE
        WRITE(YSTEP,'(I4)') JWRK1
      ENDIF
      !
      YFILE_RAD=TRIM(YRAD_BASE_NAME(JSEQ))//'_'//TRIM(YSTEP)//'.txt'
      !
      write(*,*) YFILE_RAD,'=fichier radar traite' !
      !
      IF ((INDEX(YRAD_BASE_NAME(JSEQ),'PLUVIOS')==0).AND.&
        (INDEX(YRAD_BASE_NAME(JSEQ),'WRFLLV')==0)) THEN
        write(*,*) 'SIZE(XGRID_PAR):',SIZE(YSC%UG%G%XGRID_PAR)
        CALL READ_INTERP_QPE_RADAR(YSC,YFILE_RAD,ZRAIN_TMP)
      ELSE
        CALL READ_INTERP_QPE_PLUVIOS(YSC,YFILE_RAD,ZRAIN_TMP)
      ENDIF
      ! If the rain provided by radar is undefined or equal to 0
      ! we keep the safran data
      WHERE ((ZRAIN_TMP<=XUNDEF).AND.(ZRAIN_TMP/=0.)) ZRAD_RR1(JWRK1,:)=ZRAIN_TMP(:)/ZSTEP_FORCING
      WHERE (ZRAIN_TMP>=XUNDEF)                       ZRAD_RR1(JWRK1,:)=XUNDEF
      IF (((INB_HOUR_SAF/INB_STEP_RAD)>1).AND.(YTYPE_SEQUENCES(JSEQ)/='MODRAD')) THEN
        ZRAIN(JSTP,:)= ZRAD_RR1(JWRK1*INB_HOUR_SAF/INB_STEP_RAD,:)
        DO  JWRK1=1,(INB_HOUR_SAF/INB_STEP_RAD)-1
          ZRAIN(JSTP,:)=ZRAIN(JSTP,:)+&
                   ZRAD_RR1((JWRK1*INB_HOUR_SAF/INB_STEP_RAD)-JWRK1,:)
                     ! Cumulation on time steps
        ENDDO!JWRK1
      ELSE
        WHERE (ZRAD_RR1(JWRK1,:)/=XUNDEF) ZRAIN(JSTP,:)=ZRAD_RR1(JWRK1,:)
      ENDIF
    ENDDO!JSTP
    DEALLOCATE(ZRAD_RR1)
    DEALLOCATE(ZRAIN_TMP)
    !
    IF (YTYPE_SEQUENCES(JSEQ)=='WRFLLV') THEN
      DO JWRK=IEND(JSEQ),IBEG(JSEQ)+1,-1
        ZRAIN(JWRK,:)=ZRAIN(JWRK,:)-ZRAIN(JWRK-1,:)
      ENDDO
    ENDIF
  ENDIF !(YTYPE_SEQUENCES(JSEQ)=='RADAR ' .OR. 'MODRAD' .OR. 'WRFLLV')
ENDDO !JSEQ
!
! some corrections
WHERE (ZTA>=XUNDEF)       ZTA        = 291.00
WHERE (ZQA>=XUNDEF)       ZQA        = 0.013
WHERE (ZDIR_SW_1H>=XUNDEF) ZDIR_SW_1H = 200.00
WHERE (ZSCA_SW>=XUNDEF)   ZSCA_SW    = 0.
WHERE (ZLW_1H>=XUNDEF)     ZLW_1H     = 400.0
WHERE (ZPS>=XUNDEF)       ZPS        = 101500.0
WHERE (ZWINDSPEED>=XUNDEF)ZWINDSPEED = 1.000
WHERE (ZWINDDIR>=XUNDEF)  ZWINDDIR   = 0.000
WHERE (ZRAIN>=XUNDEF)     ZRAIN(:,:) = 0.
WHERE (ZSNOW>=XUNDEF)      ZSNOW(:,:) = 0.
WHERE (ZDIR_SW_1H>10000) ZDIR_SW_1H = 200.
WHERE (ZDIR_SW_1H<0.001) ZDIR_SW_1H = 0.
WHERE (ZLW_1H<0.001)     ZLW_1H     = 0.
WHERE (ZLW_1H>10000)     ZLW_1H     = 400.
WHERE (ZPS<80000.0)      ZPS        = 101500.0
WHERE (ZRAIN<0.00000001) ZRAIN      = 0.
WHERE (ZSNOW<0.0000001)  ZSNOW      = 0.
!
write(*,*) ' FIN LECTURE FORCAGE'
!
GPERTRAIN=.FALSE.
!
DO JSEQ=1,INB_FORC_SEQUENCES
  IF (YTYPE_SEQUENCES(JSEQ)=='PERTRR') GPERTRAIN=.TRUE.
ENDDO
write(*,*) 'GPERTRAIN',GPERTRAIN
!
IF (GPERTRAIN) THEN
!----------------------------------------------------------------------------
!
!*    6.      Modification of AROME (NB_STEPS_MODIF hours of rain)
!             -------------------------------------------
!----------------------------------------------------------------------------
!*         6.1 Finding rain field structure for NB_STEPS_MODIF hours of AROME
!          and modifing rainfall location and amplitude
!              -------------------------------------------
  ALLOCATE(YFILE_FORCING_NEW(INB_MEMBERS))
  NNB_MAX_MEMBERS=500
  ALLOCATE (ZLARGE_RAIN_NEW(NNI,NNB_MAX_MEMBERS))
  ZLARGE_RAIN_NEW(:,:)=XUNDEF
  ALLOCATE (ZFLD_INTERPOL(INI,1))
  ALLOCATE (ZPNT_INTERPOL(NNI,1))
  ALLOCATE (XRAIN_NEW(NNB_STEPS_MODIF,INI,NNB_MAX_MEMBERS))
  ALLOCATE (XPROB_DECAL(NNB_MAX_MEMBERS))
  write(*,*) 'NNB_MAX_MEMBERS ',NNB_MAX_MEMBERS!,NX,NY
  XPROB_DECAL(:)=0.
!
!!!!!!!!!!!!!!!!!!!!!!!location!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  DO JSTP=1,NNB_STEPS_MODIF 
    ZLARGE_RAIN_NEW(:,:)=0.
    CALL  CHANGE_LOC_OBJECTS(ZLARGE_RR1(JSTP,:),NX,NY,JSTP,&
             ZLARGE_RAIN_NEW(:,:))!O2
    DO JMB=1,NNB_MEMBERS_LOC
      ZPNT_INTERPOL(:,1)=ZLARGE_RAIN_NEW(:,JMB)
      CALL HOR_INTERPOL(YSC%DTCO,YSC%U,YSC%GCP,ILUOUT,ZPNT_INTERPOL,ZFLD_INTERPOL)
      XRAIN_NEW(JSTP,:,JMB)=ZFLD_INTERPOL(:,1)
    ENDDO !jmb=membres
  ENDDO !jstep
  write(*,*) 'avt select'
  !
  CALL SELECT_MEMBER(YSC%U)
  !
  ALLOCATE (ZRAIN_SEL(NNB_STEPS_MODIF,INI,INB_MEMBERS))
  DO JMB=1,INB_MEMBERS
     ZRAIN_SEL(:,:,JMB)=XRAIN_NEW(:,:,NMEMBER_SEL(JMB))
  ENDDO
  write(*,*) 'zrainsel ok'
  !
!!!!!!!!!!!!!!!!!!!!!!!amplitude!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ALLOCATE (XRAIN_NEWSEL(NNB_STEPS_MODIF,INI,INB_MEMBERS))
  ALLOCATE(ZRAIN_IN(INI,NNB_STEPS_MODIF))
  ALLOCATE(ZRAIN_OUT(INI,NNB_STEPS_MODIF))
  !
  DO JMB=1,INB_MEMBERS
    DO JSTP=1,NNB_STEPS_MODIF 
      ZRAIN_IN(:,JSTP)=ZRAIN_SEL(JSTP,:,JMB)
    ENDDO
    CALL CHANGE_AMPLI_RAIN(ZRAIN_IN(:,:),NIMAX,NJMAX,ZRAIN_OUT(:,:))
    ZRAIN_OUT(:,:)=ZRAIN_IN(:,:)
    DO JSTP=1,NNB_STEPS_MODIF 
      XRAIN_NEWSEL(JSTP,:,JMB)=ZRAIN_OUT(:,JSTP)
    ENDDO
    write(*,*) 'member',JMB,', min1:',MINVAL(XRAIN_NEWSEL(1:NNB_STEPS_MODIF,:,JMB))*3600.,&
                        ', max1:',MAXVAL(XRAIN_NEWSEL(1:NNB_STEPS_MODIF,:,JMB))*3600.,&
             ', min2:',MINVAL(ZRAIN(NNB_STEPS_MODIF:INB_FORC_STP,:))*3600.,&
             ', max2:',MAXVAL(ZRAIN(NNB_STEPS_MODIF:INB_FORC_STP,:))*3600.
  ENDDO
  DEALLOCATE(ZRAIN_IN)
  DEALLOCATE(ZRAIN_OUT)
ELSE !GPERRAIN
  INB_MEMBERS=1 !(to enter the JMB loop)
ENDIF !GPERTRAIN
!----------------------------------------------------------------------------
!
!*    4.      Writing OUTPUT file(s)
!----------------------------------------------------------------------------
XNETCDF_FILENAME_OUT= &
                                       (/'ISBA_VEG_EVOLUTION.OUT.nc  ',&
                                         'ISBA_VEG_EVOLUTION_P.OUT.nc', &
                                         'ISBA_VEG_EVOLUTION_A.OUT.nc', &
                                         'ISBA_PROGNOSTIC.OUT.nc     ',&
                                         'ISBA_DIAGNOSTICS.OUT.nc    ',&
                                         'ISBA_DIAG_CUMUL.OUT.nc     ',&
                                         'SEAFLUX_PROGNOSTIC.OUT.nc  ',&
                                         'SEAFLUX_DIAGNOSTICS.OUT.nc ',&
                                         'SEAFLUX_DIAG_CUMUL.OUT.nc  ',&
                                         'WATFLUX_PROGNOSTIC.OUT.nc  ',&
                                         'WATFLUX_DIAGNOSTICS.OUT.nc ',&
                                         'WATFLUX_DIAG_CUMUL.OUT.nc  ',&
                                         'FLAKE_PROGNOSTIC.OUT.nc    ',&
                                         'FLAKE_DIAGNOSTICS.OUT.nc   ',&
                                         'FLAKE_DIAG_CUMUL.OUT.nc    ',&
                                         'TEB_PROGNOSTIC.OUT.nc      ',&
                                         'TEB_DIAGNOSTICS.OUT.nc     ',&
                                         'TEB_CANOPY.OUT.nc          ',&
                                         'PARAMS.nc                  ',&
                                         'FORCING.nc                 ',&
                                         'SURF_ATM_DIAGNOSTICS.OUT.nc',&
                                         'SURF_ATM_DIAG_CUMUL.OUT.nc ',&
                                         'SURF_ATM_DIAGNOSTIC1.OUT.nc',&
                                         'SURF_ATM_DIAGNOSTIC2.OUT.nc',&
                                         'SURF_ATM_DIAGNOSTIC3.OUT.nc'/)  

ZSTEP_OUTPUT=ZSTEP_FORCING
ALLOCATE(ZZREF(INI))
ALLOCATE(ZUREF(INI))
ZZREF(:)=2.
ZUREF(:)=10.
!----------------------------------------------------------------------------
!      
!*       4.2     Writing of FORCING.nc file
!               --------------------------
!
!----------------------------------------------------------------------------
!      
!      
!      4.2.01    grid definition
!            ---------------
!
ALLOCATE(ZLAT(INI))
ALLOCATE(ZLON(INI))
ZLAT(:)=YSC%UG%G%XLAT(:)
ZLON(:)=YSC%UG%G%XLON(:)
!
!
!----------------------------------------------------------------------------
!        4.2.1    define dimensions
!               -----------------
!
KDIMS(1) = INB_FORC_STP  ! time dimension
KDIMS(2) = INI    ! space dimension
!
!----------------------------------------------------------------------------
!      
!        4.2.2    define dimension names
!               ----------------------
!
YNAME_DIM(1) = 'time'
YNAME_DIM(2) = 'Number_of_points'
!
ALLOCATE(ZTM(INB_FORC_STP))
ZDEN = 1.
YATB(1)='seconds since '
!
DO JSTP = 1, INB_FORC_STP
  ZTM(JSTP) = ZSTEP_FORCING/ZDEN * (JSTP-1)
ENDDO
ISEC=MAX(0,NINT(TDTCUR%TIME))
ITIME(1)=FLOOR(ISEC/3600.)
ITIME(2)=FLOOR((ISEC-ITIME(1)*3600)/60.)
ITIME(3)=ISEC-ITIME(1)*3600-ITIME(2)*60 
!
CALL WRITE_TIME(TDTCUR%TDATE%YEAR,1,"-",YATB(1))
CALL WRITE_TIME(TDTCUR%TDATE%MONTH,0,"-",YATB(1))
CALL WRITE_TIME(TDTCUR%TDATE%DAY,0,"",YATB(1))
CALL WRITE_TIME(ITIME(1),1,":",YATB(1))
CALL WRITE_TIME(ITIME(2),0,":",YATB(1))
CALL WRITE_TIME(ITIME(3),0,"",YATB(1))
!
! CAUTION : Here treatment differs if PERTURBAtion of rainfall is required
!----------------------------------------------------------------------------
!      
!        4.2.3    create file
!               -----------
IF (GPERTRAIN) &
  write(*,*) 'member 0 , min1:',MINVAL(ZRAIN(1:NNB_STEPS_MODIF,:))*3600.,', max1:',MAXVAL(ZRAIN(1:NNB_STEPS_MODIF,:))*3600.,&
             ', min2:',MINVAL(ZRAIN(NNB_STEPS_MODIF:INB_FORC_STP,:)),', max2:',MAXVAL(ZRAIN(NNB_STEPS_MODIF:INB_FORC_STP,:))

!
WRITE(YMEMBER,'(I1)') 0
DO JMB=1,INB_MEMBERS
  IF (.NOT.GPERTRAIN) THEN !(if not gpertrain nor gpecosm, INB_MEMBERS=1)
    CALL CREATE_FILE(YFILE_FORCING_OUT,KDIMS,YNAME_DIM,IFILE_ID,IDDIM)
  ELSE
    IF (JMB<10) THEN
      WRITE(YMEMBER,'(I1)') JMB
      YMEMBER='P0'//YMEMBER
    ELSEIF (JMB < 100) THEN
      WRITE(YMEMBER,'(I2)') JMB
      YMEMBER='P'//YMEMBER
    ENDIF
    YFILE_FORCING_NEW(JMB)=YFILE_FORCING_OUT//'.'//YMEMBER
    CALL CREATE_FILE(YFILE_FORCING_NEW(JMB),KDIMS,YNAME_DIM,IFILE_ID,IDDIM)
    write(*,*) JMB,YFILE_FORCING_NEW(JMB)!,KDIMS,YNAME_DIM,IFILE_ID,IDDIM
    !
    DO JSEQ=1,INB_FORC_SEQUENCES
      IF (YTYPE_SEQUENCES(JSEQ)=='PERTRR') THEN
        write(*,*) "Modification of RAIN by PERTRR method for sequence JSEQ =",JSEQ
        write(*,*) "First time step of modification:",IBEG(JSEQ)
        write(*,*) "Last time step of modification:",MIN(IBEG(JSEQ)+NNB_STEPS_MODIF-1,IEND(JSEQ))
        ZRAIN(IBEG(JSEQ):MIN(IBEG(JSEQ)+NNB_STEPS_MODIF-1,IEND(JSEQ)),:) = &
                                XRAIN_NEWSEL(1:NNB_STEPS_MODIF,:,JMB)
      ENDIF
    ENDDO !JSEQ
  ENDIF!(.NOT.GPERTRAIN)
!----------------------------------------------------------------------------
  IF(GWRITE_LLV) THEN 
    YFILE_ASC='llv.'//YMEMBER
    CALL OPEN_FILE('ASCII ',IUNIT,YFILE_ASC,'FORMATTED',HACTION='WRITE')
    DO JWRK1=1,INI
      DO JSTP=1,INB_FORC_STP
        WRITE(IUNIT,*) ZLAT(JWRK1),ZLON(JWRK1),ZRAIN(1:INB_FORC_STP,JWRK1)
      ENDDO
    ENDDO
    CALL CLOSE_FILE('ASCII ',IUNIT) 
    GOTO 111
  ENDIF
  !      
  !        4.2.4    Write into forcing file
  !               -----------------------
  !
  YATT_TITLE(1) = 'units'
  !
  CALL WRITE_NETCDF(IFILE_ID,'FRC_TIME_STP','Forcing_Time_Step',ZSTEP_FORCING)
  !
  YATT(1)='seconds since '
  CALL WRITE_NETCDF(IFILE_ID,'time','Time_since_beginning',ZTM,IDDIM(1),YATT_TITLE(1:1),YATB(1:1))
  !
  CALL WRITE_NETCDF(IFILE_ID,'LON','Longitude', ZLON, IDDIM(2))
  CALL WRITE_NETCDF(IFILE_ID,'LAT','Latitude',  ZLAT, IDDIM(2))
  !
  CALL WRITE_NETCDF(IFILE_ID,'ZS','Surface_Orography',ZZS, IDDIM(2))
  !
  YATT_TITLE(1)='units'
  YATT(1) = 'm'
  CALL WRITE_NETCDF(IFILE_ID,'ZREF','Reference_Height',ZZREF,IDDIM(2),YATT_TITLE(1:1),YATT(1:1))
  CALL WRITE_NETCDF(IFILE_ID,'UREF','Reference_Height_for_Wind',ZUREF,IDDIM(2),YATT_TITLE(1:1),YATT(1:1))
  !
  ! 2D VARIABLES WITH 2 COMMENTS
  !
  YATT_TITLE(1) = 'measurement_height' 
  YATT      (1) = '2m'
  YATT_TITLE(2) = 'units'
  !
  YATT      (2) = 'K'
  CALL WRITE_NETCDF(IFILE_ID,'Tair','Near_Surface_Air_Temperature',TRANSPOSE(ZTA(1:INB_FORC_STP,:)),&
        IDDIM(2),IDDIM(1),YATT_TITLE(1:2),YATT(1:2))
  !
  YATT      (2) = 'Kg/Kg'
  CALL WRITE_NETCDF(IFILE_ID,'Qair','Near_Surface_Specific_Humidity',TRANSPOSE(ZQA(1:INB_FORC_STP,:)),&
        IDDIM(2),IDDIM(1),YATT_TITLE(1:2),YATT(1:2))
  !
  YATT_TITLE(1) = 'measurement height' 
  YATT      (1) = '10m'
  YATT_TITLE(2) = 'unit'
  YATT      (2) = 'm/s'
  CALL WRITE_NETCDF(IFILE_ID,'Wind','Wind_Speed',TRANSPOSE(ZWINDSPEED(1:INB_FORC_STP,:)),&
        IDDIM(2),IDDIM(1),YATT_TITLE(1:2),YATT(1:2))
  !
  !
  ! 2D VARIABLES WITH 1 COMMENT
  !
  YATT_TITLE(1) = 'units' 
  !
  YATT(1) = 'W/m2'
  CALL WRITE_NETCDF(IFILE_ID,'DIR_SWdown','Surface_Indicent_Direct_Shortwave_Radiation' ,TRANSPOSE(ZDIR_SW_1H(1:INB_FORC_STP,:)),&
        IDDIM(2),IDDIM(1),YATT_TITLE(1:2),YATT(1:2))
  !
  YATT(1) = 'W/m2'
  CALL WRITE_NETCDF(IFILE_ID,'SCA_SWdown','Surface_Incident_Diffuse_Shortwave_Radiation' ,TRANSPOSE(ZSCA_SW(1:INB_FORC_STP,:)),&
        IDDIM(2),IDDIM(1),YATT_TITLE(1:2),YATT(1:2))
  !
  YATT(1) = 'W/m2'
  CALL WRITE_NETCDF(IFILE_ID,'LWdown','Surface_Incident_Longwave_Radiation' ,TRANSPOSE(ZLW_1H(1:INB_FORC_STP,:)),&
        IDDIM(2),IDDIM(1),YATT_TITLE(1:2),YATT(1:2))
  !
  YATT(1) = 'Pa'
  CALL WRITE_NETCDF(IFILE_ID,'PSurf','Surface_Pressure',TRANSPOSE(ZPS(1:INB_FORC_STP,:)),&
        IDDIM(2),IDDIM(1),YATT_TITLE(1:2),YATT(1:2))
  !
  YATT(1) = 'Kg/m2/s'
  CALL WRITE_NETCDF(IFILE_ID,'Rainf','Rainfall_Rate',TRANSPOSE(ZRAIN(1:INB_FORC_STP,:)),&
        IDDIM(2),IDDIM(1),YATT_TITLE(1:2),YATT(1:2))
  !
  YATT(1) = 'Kg/m2/s'
  CALL WRITE_NETCDF(IFILE_ID,'Snowf','Snowfall_Rate',TRANSPOSE(ZSNOW(1:INB_FORC_STP,:)),&
        IDDIM(2),IDDIM(1),YATT_TITLE(1:2),YATT(1:2))
  !
  YATT(1) = 'Kg/m3'
  CALL WRITE_NETCDF(IFILE_ID,'CO2air','Near_Surface_CO2_Concentration',TRANSPOSE(ZCO2(1:INB_FORC_STP,:)),&
        IDDIM(2),IDDIM(1),YATT_TITLE(1:2),YATT(1:2))
  !
  YATT(1) = 'deg'
  CALL WRITE_NETCDF(IFILE_ID,'Wind_DIR','Wind_Direction',TRANSPOSE(ZWINDDIR(1:INB_FORC_STP,:)),&
        IDDIM(2),IDDIM(1),YATT_TITLE(1:2),YATT(1:2))
  !
  !              4.2.5 closing file
  !                  ------------
  !
  IRET=NF_CLOSE(IFILE_ID)
  IF(.NOT.GPERTRAIN) THEN 
    write(*,*) 'Sortie boucle'
  GOTO 111
  ENDIF
  write(*,*) JMB
ENDDO!INB_MEMBERS
111 CONTINUE

write(*,*) ' FIN ECRITURE FORCING.nc'
WRITE(ILUOUT,*) ' '
WRITE(ILUOUT,*) '    ------------------------------------'
WRITE(ILUOUT,*) '    | CREATE_FORCING ENDS CORRECTLY     |'
WRITE(ILUOUT,*) '    ------------------------------------'
!
CLOSE(ILUOUT)

WRITE(*,*) '    ------------------------------------'
WRITE(*,*) '    | CREATE_FORCING ENDS CORRECTLY     |'
WRITE(*,*) '    ------------------------------------'
!
!
IF (LHOOK) CALL DR_HOOK('CREATE_FORC',1,ZHOOK_HANDLE)
!
END PROGRAM CREATE_FORC
