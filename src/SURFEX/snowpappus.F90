!SFX_LIC Copyright 2001-2019 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!SFX_LIC for details. version 1.
MODULE SNOWPAPPUS_ENGINE
  !
  !  ▄▀▀▀▀▄  ▄▀▀▄ ▀▄  ▄▀▀▀▀▄   ▄▀▀▄    ▄▀▀▄  ▄▀▀▄▀▀▀▄  ▄▀▀█▄   ▄▀▀▄▀▀▀▄  ▄▀▀▄▀▀▀▄  ▄▀▀▄ ▄▀▀▄  ▄▀▀▀▀▄
  ! █ █   ▐ █  █ █ █ █      █ █   █    ▐  █ █   █   █ ▐ ▄▀ ▀▄ █   █   █ █   █   █ █   █    █ █ █   ▐
  !    ▀▄   ▐  █  ▀█ █      █ ▐  █        █ ▐  █▀▀▀▀    █▄▄▄█ ▐  █▀▀▀▀  ▐  █▀▀▀▀  ▐  █    █     ▀▄
  ! ▀▄   █    █   █  ▀▄    ▄▀   █   ▄    █     █       ▄▀   █    █         █        █    █   ▀▄   █
  !  █▀▀▀   ▄▀   █     ▀▀▀▀      ▀▄▀ ▀▄ ▄▀   ▄▀       █   ▄▀   ▄▀        ▄▀          ▀▄▄▄▄▀   █▀▀▀
  !  ▐      █    ▐                     ▀    █         ▐   ▐   █         █                     ▐
  !         ▐                               ▐                 ▐         ▐
  !
  IMPLICIT NONE
CONTAINS

  !SNOWPAPPUS SUBROUTINE
  !SUSPENSION_PAPPUS SUBROUTINE
  !VFALL_SUSPENSION SUBROUTINE
  !SALTATION_PAPPUS SUBRTOUINE
  !WINDSPEEDTHRESOLD SUBROUTINE
  !SNOWAPPETIZER SUBROUTINE
  !SNOWEATER SUBROUTINE

  SUBROUTINE SNOWPAPPUS(HSNOWPAPPUSERODEPO,HSALTPAPPUS,HPAPPUSSUBLI,HLIMVFALL,OPAPPULIMTFLUX,OPAPPUDEBUG,HSNOWMOB,&
                        PTSTEP,PMESH_SIZE,PZ0,PUREF, PVMOD,PVDIR,PRHOA,PPS,PQA,PTA,PDIRCOSZW, PEK, PSNOWDEBTC,&
                        PTHRUFAL,XQ_OUT,PQT_TOT,PBLOWSNWFLUX_1M,PBLOWSNWFLUXINT,PBLOWSNW,PAPPUS_DEBUG,OREG_GRID)
    !
    ! ROUTINE SNOW PAPPUS - ROUTINE DE TRANSPORT DE NEIGE
    !      HADDJERI Ange and BARON Matthieu
    !
    ! Snowpappus subroutine is cut in differents blocks
    !  1 -   Initialisation Block
    !          Varibale initialisation
    !  2 -   Computation Block
    !          Loop on point and compute different fluxes
    !          flux limitation sub-Block
    !  3 -   Comunication Block
    !          MPI communication
    !  4 -   Divergence Block
    !          Compute divergence on domain or subdomaine (if MPI)
    !  5 -   Erosion Block
    !          Modification of snowpack
    !  6 -   Ending Block
    !          prepare output variable
    !
    ! lister les options possibles
    !
    !
    !           `::`
    !            /
    !           `     __/\__
    !                 \_\/_/
    !                 /_/\_\
    !                   \/
    !                  /
    !                '     `::`
    !     .;:;             /
    !      ::;            '
    !   _ ';:;;'
    !   >'. ||  _
    !   `> \||.'<
    !     `>|/ <`
    !      `||/`
    !


    USE MODD_SURF_PAR,   ONLY : XUNDEF ! NETCDF FILL VALUE
    USE MODD_ISBA_n,      ONLY : ISBA_PE_t
    !  USE MODD_CSTS,        ONLY : XLMTT, XTT, XCI,XRHOLW
    USE MODD_SLOPE_EFFECT, ONLY: NIX,NIY ! DANS OFFLINE
    USE MODE_SNOW3L
#ifdef SFX_MPI
    USE MODD_SURFEX_MPI, ONLY: NSIZE_TASK,NRANK,NPROC,NCOMM !POSSIBLE BUG SI NON DEFINI
#endif
    USE YOMHOOK,       ONLY : LHOOK,   DR_HOOK
    USE PARKIND1,      ONLY : JPRB
    USE OMP_LIB
    USE MODD_SNOW_METAMO, ONLY : XSNOWDZMIN,XUEPSI ! MINIMUM DZ D'UNE COUCHE
    USE MODI_ABOR1_SFX
    USE MODD_SNOW_PAR, ONLY :  XRHODEPPAPPUS, XDIAMDEPPAPPUS, XSPHDEPPAPPUS, XLFETCHPAPPUS
    USE MODD_CSTS,     ONLY : XTT,XCI,XRHOLW,XLMTT
    !  USE MODD_TYPE_SNOW

    IMPLICIT NONE

#ifdef SFX_MPI
    INCLUDE "mpif.h"
#endif
    !
    CHARACTER(3),INTENT(IN)          :: HSNOWPAPPUSERODEPO ! deposition rate calculation option$
    LOGICAL,INTENT(IN)               :: OPAPPULIMTFLUX ! true = activate snow flux limitation
    LOGICAL,INTENT(IN)               :: OPAPPUDEBUG ! TRUE = active debug mode warnings and debug output in snowpappus
    CHARACTER(4),INTENT(IN)          :: HSNOWMOB ! mobility index calculation option
    CHARACTER(4),INTENT(IN)          :: HPAPPUSSUBLI ! NAMELIST SUBLIMATION CONTROLING OPTION
    CHARACTER(3),INTENT(IN)          :: HSALTPAPPUS ! namelist option controling saltation
    CHARACTER(4),INTENT(IN)          :: HLIMVFALL ! gives the option to decide what is old or new snow for fall speed calculation
    REAL, INTENT(IN)                 :: PTSTEP ! time step of calculation (s)
    REAL, INTENT(IN)                 :: PMESH_SIZE ! mesh (pixel) size (m)
    REAL, DIMENSION(:),INTENT(IN)    :: PZ0 !roughness length for momentum (m) CHANGEMENT PAR LA RUGOSITÉ DE LA NEIGE
    REAL, DIMENSION(:),INTENT(IN)    :: PUREF !reference height of the wind (m)
    REAL, DIMENSION(:),INTENT(IN)    :: PVMOD !modulus of the wind parallel to the orography (m/s)
    REAL, DIMENSION(:),INTENT(IN)    :: PRHOA !air density(kg/m3)
    REAL, DIMENSION(:),INTENT(IN)    :: PPS !AIR PRESURE
    REAL, DIMENSION(:),INTENT(IN)    :: PQA !AIR SPECIFIC HUMIDITY
    REAL, DIMENSION(:), INTENT(IN)   :: PTA    ! atmospheric temperature at level za (K)
    REAL, DIMENSION(:),INTENT(IN)    :: PVDIR !wind direction (rad)
    REAL, DIMENSION(:), INTENT(IN)   :: PDIRCOSZW !cosine of the slope angle (i.e. projection of the normal vector on the vertical axis)
    LOGICAL, INTENT(IN)              :: OREG_GRID ! bool for the activation of the 2d part of snowpapp necessiting regular grid so the G%XMESH_SIZE pointer defined and NIX!=0
    !
    TYPE(ISBA_PE_t), INTENT(INOUT)             :: PEK ! PEK variable type tsnow contain snow layers properties see MODD_TYPE_SNOW.F90
    REAL, DIMENSION(SIZE(PVMOD)), INTENT(INOUT):: PSNOWDEBTC !cumulated amount of snow which was not removed on the point because it became snowfree (kg/m2)
    REAL, DIMENSION(SIZE(PVMOD)), INTENT(INOUT):: PTHRUFAL ! PTHRUFAL  = rate that liquid water leaves snow pack: paritioned into soil infiltration/runoff by ISBA [kg/(m2 s)]
    !
    REAL, DIMENSION(SIZE(PVMOD),3), INTENT(OUT):: XQ_OUT ! ued for output and is passend to diagnostics contain ZQSALT in (:,1) and ZQSUSP in (:,2) (kg/ms) and zqsubl in (:)
    REAL, DIMENSION(SIZE(PVMOD)), INTENT(OUT)  :: PQT_TOT ! sum of saltating and suspensive transport rate (ZQSALT + ZQSUSP) (kg/ms)
    REAL, DIMENSION(SIZE(PVMOD)), INTENT(OUT)  :: PBLOWSNWFLUX_1M! mass flux 1 meter above the snow surface (kg/m2/s) (for comp with
    !CLB SPC data from datapaper ( Vionnet 2016 ))
    REAL, DIMENSION(SIZE(PVMOD)), INTENT(OUT)  :: PBLOWSNWFLUXINT! mass flux 0.2-1.2  meter above the snow surface (kg/m2/s) (for comp with
    !CLB SPC data from datapaper ( Vionnet 2016 ))
    REAL, DIMENSION(SIZE(PVMOD),4), INTENT(OUT):: PBLOWSNW   ! Properties of deposited blowing snow
    !    1 : Deposition flux (kg/m2/s)
    !    2 : Density of deposited snow (kg/m3)
    !    3 : SGRA1 of deposited snow
    !    4 : SGRA2 of deposited snow
    REAL,DIMENSION(SIZE(PVMOD),11),INTENT(OUT)  :: PAPPUS_DEBUG !Output fordebugging => (JJ,n) avec n= 1 ZFRIC 2 ZFRIC_T 3 PZ0 4 ZHSALT_SUSP 5 ZHSALT1 6 ZHSALT2  7 ZVFALL 8 JJ 9 Nrank 10 ILOCNIY 11 NSIZE_TASK
    !
    REAL, DIMENSION(SIZE(PVMOD))                :: ZQDEP_TOT ! total blowing snow deposition rate ( kg/m2/s )
    REAL, DIMENSION(SIZE(PVMOD))                :: ZQDEP_SALT ! divergence result for saltation (kg/ms)
    REAL, DIMENSION(SIZE(PVMOD))                :: ZQDEP_SUSP ! divergence result for suspension (kg/ms)
    REAL, DIMENSION(3,SIZE(PVMOD))              :: ZQMPI  !TABLEAU DES FLUX DE TRANSPORT USED FOR MPI (S,JJ) JJ=POINT NUMBER S=1 Q SALTATION S=2 Q SUSPENSION S=3 PVDIR wind direction
    ! the shape is choosen like this (3,jj) for JJ point to be continuous in memory and thus give faster mpi comm
    REAL, DIMENSION(SIZE(PEK%TSNOW%WSNOW,1),SIZE(PEK%TSNOW%WSNOW,2))  :: ZSNOWSWE ! Snow Water Equivalent of the snow layers (kg/m2)
    REAL, DIMENSION(SIZE(PEK%TSNOW%WSNOW,1),SIZE(PEK%TSNOW%WSNOW,2))  :: ZSNOWHEAT !SNOW LAYER HEAT CONTENT (J/m3) !
    REAL, DIMENSION(SIZE(PEK%TSNOW%WSNOW,1),SIZE(PEK%TSNOW%WSNOW,2))  :: ZSNOWRHO ! Snow layer(s) averaged density (kg/m3)
    REAL, DIMENSION(SIZE(PEK%TSNOW%WSNOW,1),SIZE(PEK%TSNOW%WSNOW,2))  :: ZSNOWHIST ! Snow layer(s) grain historical parameter
    REAL, DIMENSION(SIZE(PEK%TSNOW%WSNOW,1),SIZE(PEK%TSNOW%WSNOW,2))  :: ZSNOWAGE ! Snow layer(s) age
    REAL, DIMENSION(SIZE(PEK%TSNOW%WSNOW,1),SIZE(PEK%TSNOW%WSNOW,2))  :: ZSNOWSPHERI !snow sphericity
    REAL, DIMENSION(SIZE(PEK%TSNOW%WSNOW,1),SIZE(PEK%TSNOW%WSNOW,2))  :: ZSNOWDIAMOPT !snow optical diameter (m)
    REAL, DIMENSION(SIZE(PEK%TSNOW%WSNOW,1),SIZE(PEK%TSNOW%WSNOW,2))  :: ZSNOWDZ ! DEPTH OF SNOW LAYER (m)
    REAL, DIMENSION(SIZE(PEK%TSNOW%WSNOW,1),SIZE(PEK%TSNOW%WSNOW,2))  :: ZSNOWDRYMASS ! dry snow mass (kg/m2)
    REAL, DIMENSION(SIZE(PEK%TSNOW%WSNOW,1),SIZE(PEK%TSNOW%WSNOW,2))  :: ZSNOWLIQ ! snow liquid mass (kg/m2)
    REAL, DIMENSION(SIZE(PEK%TSNOW%WSNOW,1),SIZE(PEK%TSNOW%WSNOW,2))  :: ZSCAP ! snow heat capacity (J/Km^3)
    REAL, DIMENSION(SIZE(PEK%TSNOW%WSNOW,1),SIZE(PEK%TSNOW%WSNOW,2))  :: ZSNOWTEMP ! snow temperature (K)
    ! local copies of variables for PEK%type_snow, geting states to snow layers
    REAL, DIMENSION(SIZE(PVMOD))     :: ZQSUSP !transport rate in suspension (kg/m/s)
    REAL, DIMENSION(SIZE(PVMOD))     :: ZQSUBRATE !transport rate in suspension (kg/m2/s) (defined >0)
    REAL, DIMENSION(SIZE(PVMOD))     :: ZVFRIC_T !threshold friction velocity (hauteur du sol) for snow transport (m/s)
    REAL                             :: ZVFRIC! wind friction velocity (m/s)
    REAL                             :: ZVFALL !ass averaged terminal fall velocity of snow particles at the bottom of the suspension layer (m/s)
    REAL                             :: ZCSALT !snow concentration in the saltation layer (kg/m3)
    REAL                             :: ZQSALT !TRANSPORT RATE OF SALTATING SNOW (KG/(M.S))
    REAL                             :: ZHSALT1 !hsalt selon pomeroy90 (m)
    REAL                             :: ZHSALT2 !hsalt selon pomeroy92 (m)
    REAL                             :: ZHMIN_SUSP !minimum height where flux is computed with suspension routine (m)
    REAL                             :: ZHSALT_SUSP !heigth of the bottom of suspension layer (m) A DESACTIVER
    REAL                             :: ZLFETCH ! length of upwind fetch (m)
    REAL                             :: ZFETCH_COEFF ! effect of fetch on saltation transport
    REAL,PARAMETER                   :: ZFMAX=500.! fetch distance necessary to reach equilibrium saltation concentration (see Liston 1998)
    REAL                             :: ZSNOWDINNER !SNOW MASS TO BE REMOVED BY SNOWEATER MODULE kg/m2
    INTEGER,DIMENSION(SIZE(PVMOD))   :: INLVLS_USE ! number of active snow layers
    REAL,DIMENSION(SIZE(PVMOD))      :: ZWSNOW_T ! total SWE (kg/m2)
    REAL,DIMENSION(SIZE(PVMOD))      :: ZSNOW ! total snow height (m)
    INTEGER                          :: INB_PTS !number of points => was necessary to pass lists as function arguments (I don't understand why )
    INTEGER                          :: JJ,JST ! main iterator, iteratorS for PTEMPQ_OUT COMPLETION
    INTEGER                          :: ILOCNIY ! local value of NIX for own cpu
    INTEGER                          :: IMAX_USE ! max number of layer over the domain
    REAL, DIMENSION(3,NIX)           :: Z_MPI_NORTH! NORTHERN boundary values for subdomain (/!\ always check the domain and array "orientation")
    REAL, DIMENSION(3,NIX)           :: Z_MPI_SOUTH! SOUTHERN boundary values for subdomain (/!\ always check the domain and array "orientation")
    REAL, DIMENSION(3,NIY)           :: Z_MPI_WEST! WESTERN boundary values for subdomain (/!\ always check the domain and array "orientation")
    REAL, DIMENSION(3,NIY)           :: Z_MPI_EAST !EASTERN boundary values for subdomain (/!\ always check the domain and array "orientation")
    REAL(KIND=JPRB)                  :: ZHOOK_HANDLE,ZHOOK_HANDLE2
    LOGICAL, DIMENSION(SIZE(PVMOD))  :: GTRANSPORTABLE !marks if snow can be transported on the point (take account of icy crust)
    LOGICAL                          :: GVERIF ! bool for expected result coherence verification and error raising
    REAL                             :: ZSNOWHEATJM3TOJM2 ! used to convert snowheat from Jm-3 to Jm-2
#ifdef SFX_MPI
    !MPI LOCAL VARIABLES
    INTEGER                             :: IERR !IERR
    INTEGER, DIMENSION(MPI_STATUS_SIZE) :: ISTATUS_RECV_UP !ISTATUS_SEND_UP
    INTEGER, DIMENSION(MPI_STATUS_SIZE) :: ISTATUS_RECV_DOWN!ISTATUS_SEND_DOWN,
#endif
#ifdef SFX_MPI
    REAL, ALLOCATABLE, DIMENSION(:,:,:)  :: ZDOMAINE_2D !DOMAINE REDRESSÉ EN 2D AVEC LES BORD MPI ,ILOCNIY=NSIZE_TASK(NRANK)/NIX
    REAL, ALLOCATABLE, DIMENSION(:,:,:)  :: ZDIV_2D !sous domaine mpi transposé en 2D avec les bords mpi
#endif
#ifndef SFX_MPI
    REAL, DIMENSION(NIX+2, NIY+2,3) :: ZDOMAINE_2D !DOMAINE REDRESSÉ EN 2D AVEC LES BORD MPI ,ILOCNIY=NSIZE_TASK(NRANK)/NIX
    REAL, DIMENSION(NIX+2,NIY+2,2)  :: ZDIV_2D !sous domaine mpi transposé en 2D avec les bords du domaine
#endif

    IF (LHOOK) CALL DR_HOOK('SNOWPAPPUS',0,ZHOOK_HANDLE2)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !                                            INITIALISATION BLOCK                                                   !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !initialization of variables locales

#ifdef SFX_MPI
    IF (OREG_GRID) THEN
      ALLOCATE(ZDOMAINE_2D(NIX+2,(NSIZE_TASK(NRANK)/NIX)+2,3))
      ALLOCATE(ZDIV_2D(NIX+2,(NSIZE_TASK(NRANK)/NIX)+2,3))
    ELSE
      ALLOCATE(ZDOMAINE_2D(0,0,0)) ! should no be used
      ALLOCATE(ZDIV_2D(0,0,0)) ! should not be used
    END IF
#endif

    ZSNOWSWE(:,:)=PEK%TSNOW%WSNOW(:,:)
    ZSNOWRHO(:,:)=PEK%TSNOW%RHO(:,:)
    ZSNOWHIST(:,:)=PEK%TSNOW%HIST(:,:)
    ZSNOWAGE(:,:)=PEK%TSNOW%AGE(:,:)
    ZSNOWSPHERI(:,:)=PEK%TSNOW%SPHERI(:,:)
    ZSNOWDIAMOPT(:,:)=PEK%TSNOW%DIAMOPT(:,:)
    ZSNOWHEAT(:,:)=PEK%TSNOW%HEAT(:,:)
    ZSNOWDZ(:,:)=0.
    ZSNOWDRYMASS(:,:)=0.
    ZSNOWLIQ(:,:)=0.
    ZSCAP(:,:)=0.
    ZSNOWTEMP(:,:)=0.
    INB_PTS=SIZE(PVMOD)
    ZWSNOW_T(:)=SUM(ZSNOWSWE,DIM=2)
    ZQDEP_SALT=0.
    ZQDEP_SUSP=0.
    ZLFETCH=XLFETCHPAPPUS
    ZQSALT=0.
    ZQSUSP=0.
    ZQMPI=0.
    XQ_OUT(:,:)=0.
    PQT_TOT(:)=0.
    PBLOWSNW(:,:)=0.
    PBLOWSNWFLUX_1M(:)=0.
    PBLOWSNWFLUXINT(:)=0.
    PSNOWDEBTC=0.
    ZVFRIC=XUNDEF ! fill value in netcdf
    PAPPUS_DEBUG=XUNDEF
    ZVFRIC_T(:)=XUNDEF ! fill value in netcdf
    ZHSALT1=XUNDEF ! fill value in netcdf
    ZHSALT2=XUNDEF ! fill value in netcdf
    ZVFALL=XUNDEF ! fill value in netcdf
    ZQDEP_TOT=0.
    GVERIF=.TRUE.
    ZDOMAINE_2D=0.
    ZDIV_2D=0.
    ZSNOWDINNER=0.
    ZQSUBRATE=0.
    ZSNOWHEATJM3TOJM2=0
    !INITIALISATION DE LA VALEUR LOCALE DU NOMBRE DE LIGNE POUR SOUS DOMAINE DU RuN

    IF (OREG_GRID) THEN ! Cas 2D IGN + quelque soit la compil (MPI/NoMPI)
#ifndef SFX_MPI
      !! cas 2D et no mpi
      ILOCNIY=NIY
#endif
#ifdef SFX_MPI
      ! cas 2D et cas mpi
      ! variable speciale mpi
      !INITIALISATION DE LA VALEUR LOCALE DU NOMBRE DE LIGNE POUR SOUS DOMAINE DU RANK
      ILOCNIY=NSIZE_TASK(NRANK)/NIX
      !initialisation des bord locaux
      IF (MOD(NSIZE_TASK(NRANK),NIX).NE.0)THEN
        WRITE(*,*)"DOMAINE NON CARRÉ :( BIG ERROR IN SUBDOMAIN DIVISION => MESSAGE FROM SNOWPAPPUS"
        STOP
      END IF
      IF (SIZE(PVMOD).NE.NIX*ILOCNIY) THEN
        WRITE(*,*)"BIG ERROR IN SUBDOMAIN DIVISION :( SIZE(PVMOD).NE.NIX*ILOCNIY => MESSAGE FROM SNOWPAPPUS"
        WRITE(*,*)"NIX=",NIX
        WRITE(*,*)"LOCAL NIY=",ILOCNIY
        WRITE(*,*)"SIZE(PVMOD)=",SIZE(PVMOD)
        WRITE(*,*)"NSIZE_TASK(NRANK)=",NSIZE_TASK(NRANK)
      END IF !ajouter une option de definition pou rle cas non mpi et non grille 2D
#endif
    ELSE  ! cas 1D quelque soit la compilation ou grille non IGN
      ILOCNIY=1
    END IF

    Z_MPI_NORTH(:,:)=0.
    Z_MPI_SOUTH(:,:)=0.
    Z_MPI_WEST(:,:)=0.
    Z_MPI_EAST(:,:)=0.

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !                                      END OF INITIALISATION BLOCK                                                  !
    !                                       RECOMPUTE SNOW PARAMETRES                                                   !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    INLVLS_USE(:) = 0
    ZSNOW(:) = 0
    DO JST = 1,SIZE(ZSNOWSWE(:,:),2)
      DO JJ = 1,SIZE(PVMOD)
        IF ( ZSNOWSWE(JJ,JST)>0. ) THEN
          ZSNOWDZ(JJ,JST) = ZSNOWSWE(JJ,JST) / ZSNOWRHO(JJ,JST) ! compute snow height
          ZSNOW(JJ) = ZSNOW(JJ) + ZSNOWDZ(JJ,JST)
          INLVLS_USE(JJ) = JST
        ELSE
          ZSNOWDZ(JJ,JST) = 0.
        ENDIF
      ENDDO  !  end loop snow layers
    ENDDO    ! end loop grid points
    IMAX_USE = MAXVAL(INLVLS_USE)
    ! active layers
    DO JST = 1, IMAX_USE
      DO JJ = 1,SIZE(PVMOD)
        IF (JST <= INLVLS_USE(JJ)) THEN
          !
          ZSCAP    (JJ,JST) = ZSNOWRHO(JJ,JST) * XCI
          !
          ZSNOWTEMP(JJ,JST) = XTT + &
          ( ( ZSNOWHEAT(JJ,JST)/ZSNOWDZ(JJ,JST) + XLMTT*ZSNOWRHO(JJ,JST) )/ZSCAP(JJ,JST) )
          !
          ZSNOWLIQ (JJ,JST) = MAX( 0.0, ZSNOWTEMP(JJ,JST)-XTT ) * ZSCAP(JJ,JST) * &
          ZSNOWDZ(JJ,JST) / (XLMTT*XRHOLW)
          !
          ZSNOWTEMP(JJ,JST) = MIN( XTT, ZSNOWTEMP(JJ,JST) )
          !
          ZSNOWDRYMASS (JJ, JST) = ZSNOWSWE(JJ,JST) - ZSNOWLIQ(JJ,JST)*XRHOLW
        ENDIF
      ENDDO  !  end loop active snow layers
    ENDDO

    !Determines points covered by snow
    GTRANSPORTABLE(:)=(ZSNOW(:)>0. .AND. ZSNOWHIST(:,1)<2) ! empêche le transport en cas de croûte de regel ou d'absence de neige
    !ne marche pas bien

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !                                           END OF RECOMPUTE SNOW PARAMETRES                                        !
    !                                               BEGIN COMPUTATION BLOCK                                             !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    IF (LHOOK) CALL DR_HOOK('SNOWPAPPUS_COMPT',0,ZHOOK_HANDLE)

    ! 1) Compute threshold friction velocity PVFRIC_T

    CALL WINDSPEEDTHRESOLD(ZSNOWDIAMOPT(:,1),ZSNOWSPHERI(:,1),PEK%TSNOW%AGE(:,1),HSNOWMOB,ZSNOWRHO(:,1),PZ0,&
    GTRANSPORTABLE,INB_PTS,ZVFRIC_T)

    !variable output
    PAPPUS_DEBUG(:,2)=ZVFRIC_T(:)
    PAPPUS_DEBUG(:,3)=PZ0(:)

    ! 2') computes the effect of fetch on saltation transport
    ZFETCH_COEFF=1-EXP(-3*ZLFETCH/ZFMAX)!Liston 1998, from Takeuchi 1980 measurements

    ! 2)  Loop on GRID points SIZE(PVMOD)
    DO JJ=1, SIZE(PVMOD) ! Loop on GRID points SIZE(PVMOD) CAR PVMOD EST DE LA TAILLE NSIZE_P
      !initialisation pour chaque pas de temps necessaire pour les sorties diag
      ZQSALT=0.
      ZHSALT1=XUNDEF ! fill value in netcdf
      ZHSALT2=XUNDEF ! fill value in netcdf
      ZVFALL=XUNDEF ! fill value in netcdf
      ZHSALT_SUSP=XUNDEF ! fill value in netcdf

      !compute friction velocity
      ZVFRIC=0.41*MAX(0.,PVMOD(JJ))/LOG(MAX(1.,PUREF(JJ)/PZ0(JJ)))

      IF (GTRANSPORTABLE(JJ)) THEN  !checks if there is snow on the point
        ! WRITE(*,*) "== POINT N° ", JJ, "=="
        !WRITE(*,*)ZVFRIC,PVMOD(JJ),PUREF(JJ)

        IF (ZVFRIC>=ZVFRIC_T(JJ)) THEN ! check if transport occur
          ! 3) compute saltation transport flux and boundary condition for suspension
          CALL SALTATION_P90_P92(ZVFRIC,ZVFRIC_T(JJ),PRHOA(JJ),ZQSALT,ZCSALT,ZHSALT1,ZHSALT2)
          IF (HSALTPAPPUS=='S04') THEN
            ! in this case, the computed flux up to 15 cm is replaced with Sorensen 04 flux
            CALL SALTATION_S04(ZVFRIC,ZVFRIC_T(JJ),PRHOA(JJ),ZQSALT)
            ZHMIN_SUSP=0.15
          ELSEIF(HSALTPAPPUS=='P90') THEN
            ZHMIN_SUSP=ZHSALT2
          ENDIF
          !Applies the effect of the fetch
          ZQSALT=ZFETCH_COEFF*ZQSALT
          ZCSALT=ZFETCH_COEFF*ZCSALT
          ! 4) compute suspension transport flux
          CALL VFALL_SUSPENSION(ZVFRIC,ZSNOWSPHERI(JJ,1),ZSNOWDIAMOPT(JJ,1),ZSNOWAGE(JJ,1),HLIMVFALL,ZVFALL)!terminal fall speed parametrization
          CALL SUSPENSION_PAPPUS(ZVFRIC,ZVFRIC_T(JJ),ZHSALT2,ZHMIN_SUSP,ZVFALL,ZLFETCH,ZCSALT,PZ0(JJ),ZQSUSP(JJ),&
          PBLOWSNWFLUX_1M(JJ),PBLOWSNWFLUXINT(JJ))

          ! 5) Total transport fluxes
          PQT_TOT(JJ)=ZQSUSP(JJ)+ZQSALT

          ! 6) compute transport sublimation flux
          IF (HPAPPUSSUBLI=='SBSM') THEN
            ZSNOWHEATJM3TOJM2 = ZSNOWHEAT(JJ,1) / ZSNOWRHO (JJ,1) * ZSNOWSWE (JJ,1)  ! convert heat content from J/m3 into J/m2
            CALL SBSM (PVMOD(JJ),PUREF(JJ),PZ0(JJ),PTA(JJ),PPS(JJ),PQA(JJ),INLVLS_USE(JJ),ZSNOWRHO(JJ,1),ZSNOWHEAT(JJ,1), &
            ZSNOWDZ(JJ,1),ZQSUBRATE(JJ))
          ELSEIF ( HPAPPUSSUBLI=='BJ03' ) THEN
            CALL BINTANJA_98_3M(PVMOD(JJ),PUREF(JJ),PZ0(JJ),PTA(JJ),ZQSUBRATE(JJ))
          ELSEIF ( HPAPPUSSUBLI=='BJ10' ) THEN
            CALL BINTANJA_98_10M(PVMOD(JJ),PUREF(JJ),PZ0(JJ),PTA(JJ),ZQSUBRATE(JJ))
          ELSEIF ( HPAPPUSSUBLI=='GR06' ) THEN
            CALL GORDON06 (ZVFRIC_T(JJ),PVMOD(JJ),PUREF(JJ),PZ0(JJ),PTA(JJ),PPS(JJ),PQA(JJ),PRHOA(JJ),ZQSUBRATE(JJ))
          END IF

          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          !    VERIFICATION BLOCK     !
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          IF (GVERIF) THEN
            IF (ZQSUSP(JJ)<0) THEN
              CALL ABOR1_SFX('QSUSP<0 Forbiden value')
            END IF
            IF (ZQSALT<0) THEN
              CALL ABOR1_SFX('QSALT<0 Forbiden value')
            END IF
          END IF
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          !  END VERIFICATION BLOCK   !
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        ENDIF
      ENDIF! end if snow on point

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !                  (STILL IN DO LOOP)                        !
      !                BEGIN LIMITATION BLOCK                      !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      ! ANOTHER LIMITATION OF OUR PAPPUS SCHEME IS THAT WE ARE NOT LIMITING THE EROSION FLUX AT THE SNOW QUANTITY ON POINT
      !THEREFORE MASS CONSERVATION IS NOT RESPECTED AND WE ARE CREATING SNOW IF WE ERODE MORE THAN THE SNOW ON POINT
      ! A WAY OT SOLVE THIS PROBLEM IS TO COMPUTE DEPOSITION FLOOR IT TO THE ON POINT SNOW MASS AND COMMUNICATE CHANGES IN RE-PROCESSED
      !FLUXES TO NEIBOURS. THIS MEAN TWO TIMES MORE COMMUNICATIONS
      ! ANOTHER APPROXIMATION IS TO FLOOR QT_TOT FLUXES TO THE MAXIMUM MASS OF SNOW ON POINTS * TIMESTEP
      !
      !  WITH X BEING THE SNOW FLUX COMMING IN THE PIXEL
      !  WITH dZ BEING THE SNOW MASS VARIATION ON PIXEL
      !  WITH Y BEING THE SNOW FLUX COMMING OUT THE PIXEL
      ! X=>  |  Z  |  Y=>            THE CONSERVATION EQUATION IS     Y = X+dZ*TIMESTEP
      ! BUT WITH THE UPWIND DIVERGENCE SCHEME WE HAVE BIG HYPOTHESIS THAT TRANSPORT IS LIMITED TO 1 PIXEL TRAVEL IN ONE TIME STEP
      ! THEREOFRE WE CAN SAY THAT SNOW INCOMMING CANNOT GO OUT THE PIXEL IN THE SAME TIMESTEP
      ! THIS HYPOTHESIS IS EQUIVALIENT TO X=0 FOR THE CURRENT TIME STEP AND THE CONSERVATION EQUATION IS Y=dZ*TIMESTEP
      !SO     MAX(Y)<=MAX(dZ)

      IF ( OPAPPULIMTFLUX ) THEN
        ! apply two condition to limit transport flux to swe present on pixel
        ! mean deposited snow cannot be transported again in same time step
        ! should make PSNOWDEBTC obsolete if working properly
        ! limiting errosion :
        ! Qsub must be <or = to sum(swe)/ptstep
        ! Qt must be < or =  to (sum(swe)/ptstep - qsubl)*(pmesh/cos)
        !
        IF (ZQSUBRATE(JJ)>ZWSNOW_T(JJ)/PTSTEP) THEN
          IF (OPAPPUDEBUG) THEN
            WRITE(*,*) "Warning SnowPappus : Qsubrate > SWE/Ptstep, check the code. Limiting sublimation to Qsubrate == SWE/Ptstep"
          ENDIF
          ZQSUBRATE(JJ)=ZWSNOW_T(JJ)/PTSTEP
        ENDIF
        IF (PQT_TOT(JJ)/=0) THEN
          ZQSALT=ZQSALT*MIN(((ZWSNOW_T(JJ)/PTSTEP)-ZQSUBRATE(JJ))*PMESH_SIZE/PDIRCOSZW(JJ),PQT_TOT(JJ))/PQT_TOT(JJ) ! LIMITATION OF FLUX TO MAX OF TRANSPORTABLE SNOW mass ON POINT
          ZQSUSP(JJ)=ZQSUSP(JJ)*MIN(((ZWSNOW_T(JJ)/PTSTEP)-ZQSUBRATE(JJ))*PMESH_SIZE/PDIRCOSZW(JJ),PQT_TOT(JJ))/PQT_TOT(JJ)
        END IF
        PQT_TOT(JJ)=ZQSUSP(JJ)+ZQSALT
      END IF

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !                  (STILL IN DO LOOP)                        !
      !                 END LIMITATION BLOCK                       !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !  output computation variables
      XQ_OUT(JJ,1)=ZQSALT
      XQ_OUT(JJ,2)=ZQSUSP(JJ)
      XQ_OUT(JJ,3)=ZQSUBRATE(JJ)
      ZQMPI(1,JJ)=ZQSALT ! NORME DU VECTEUR
      ZQMPI(2,JJ)=ZQSUSP(JJ)!! voir les valeur du tableau de suspension zonale
      ZQMPI(3,JJ)=PVDIR(JJ) !ANGLE DU VECTEUR
      PAPPUS_DEBUG(JJ,1)=ZVFRIC
      PAPPUS_DEBUG(JJ,5)=ZHSALT1
      PAPPUS_DEBUG(JJ,6)=ZHSALT2
      PAPPUS_DEBUG(JJ,7)=ZVFALL
      !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !END DO LOOP on points
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    END DO

    IF (LHOOK) CALL DR_HOOK('SNOWPAPPUS_COMPT',1,ZHOOK_HANDLE)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !                                                                                                                   !
    !                                              END COMPUTATION BLOCK                                                !
    !                                           BEGIN 2D PART OF SNOWAPPPUS (IF)                                        !
    !                                             BEGIN COMUNICATION BLOCK                                              !
    !                                                                                                                   !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !CONSEQUENCE DU DECOUPAGE AVEC OFFLIN/INITINDEXMPI => SET_NB_POINTS_LIN()
    !remplissage des rank MPI debut du domaine commence par N-1 VERS 0

    !             ^      -----------------------!                             Pour chaque bande de chaque processeur :
    !             |                0                                          ex: Rank = 0   ------------------------------!
    !             |       ----------------------!         North                               NIX+ILOCNIY .. JJ ..  NIX*NIY
    !             |                1                        |
    !             |      -----------------------!   West --   -- East                         1 2 3 4 5 6 7 ..JJ.. NIX-1 NIX
    !          Yaxis              ...                       |                                ------------------------------!
    !             |      -----------------------!         South                  Rank = 1     NIX+ILOCNIY .. JJ ..  NIX*NIY
    !             |              NPROC-2
    !             |      -----------------------!                                              1 2 3 4 5 6 7 ..JJ.. NIX-1 NIX
    !             |             NPROC-1                                                      ------------------------------!
    !             |     -----------------------!
    !             0------------Xaxis--------------->

    ! pour visualiser les decoupes voir les outputs suivants : 'XJJ_MPIOUT','XNRANK_MPIOUT','XILOCNIY_MPIOUT','XNSIZETASK_MPIOUT'
    ! LE DOMAINE DE SIMULATION EST DECOUPÉ EN SOUS DOMAINES ASSIGNÉ A UN PROCESSEUR
    !
    ! PRINCIPE DE BASE: RECUPERER LES VALEUR DES BORDS DES SOUS DOMAINES ADJACENT SI NECESSAIRE
    ! LE SOUS DOMAINE PARTAGE SA DERNIERE LIGNE AVEC LE DOMAINE DU DESSOUS (NRANK-1) ET ATTENDS LA PREMIERE LIGNE DU DOMAINE DU DESSOUS (SAUF SI NRANK=0)
    ! LE SOUS DOMAINE PARTAGE SA PREMIERE LIGNE AVEC LE DOMAINE DU DESSUS (NRANK +1) ET ATTENDS LE DERNIERE LIGNE DU DOMAINE DU DESSUS (SAUF SI NRANK=NPROC-1)
    ! CE FONCTIONNEMENT ENTRAIRE UN PARTAGE ET LA RECEPTION DE 4 * TAILLE DE LA LIGNE * NB BYTE D'UN REAL A CHAQUE PAS DE TEMPS
    ! POUR EVITER CELA LES LOGICALS G_ATLEASTONETRANSPORT G_ATLEASTONETRANSPORT_UP G_ATLEASTONETRANSPORT_DOWN QUI REPRESENTENT RESPECTIVEMENT SI IL Y A UN EVENEMENT
    ! DE TRANSPORT DANS LE BORDS HAUT ET LE BORD BAS A CE PAS DE TEMPS POUR LE SOUS DOMAINE ET LES SOUS DOMAINES ADJACENTS SI .TRUE;
    ! l'OBJECTIF EST DE PARTAGER LE LOGICAL G_ATLEASTONETRANSPORT ET RECUPERER CELUI DES VOISINS A CHAQUE PAS DE TEMPS
    ! SI LES VOISINS OU NOUS MEME AVONS DU TRANSPORT AU PAS DE TEMPS EN COURS ALORS ON RENTRE DANS LA BOUCLE MPI ET ON PARTAGE NOS BORDS
    !
    !DEADLOCK PROTECTION DANS LE CAS OU LE TRANSPORT EST TRES LOCALISÉ ET QU'ON RENTRE DANS LE MPI ALORS QUE LES VOISINS NON CAR ILS N'ONT PAS DE NEIGE
    !MOYEN DE FAIRE AVEC DE MPI_PROBE SANS DOUTE PLUS EFFICACE MAIS BON ...
    !COMMUNICATION DE SI IL Y A EU DU TRANSPORT A CE PAS DE TEMPS ET RECUPERATION DU STATUS DES AUTRES


    !      2D part of snowpappus
    IF (OREG_GRID) THEN ! only on 2D grids
    !
    !
#ifdef SFX_MPI
      IF (LHOOK) CALL DR_HOOK('SNOWPAPPUS_MPI',0,ZHOOK_HANDLE)
      IF (NPROC>1)THEN
        IF (NRANK==0)THEN
          CALL MPI_RECV(Z_MPI_SOUTH,SIZE(Z_MPI_SOUTH)*KIND(Z_MPI_SOUTH)/4,MPI_REAL,NRANK+1,0,NCOMM,ISTATUS_RECV_UP,IERR)
          !
          CALL MPI_SEND(ZQMPI(:,1:NIX),SIZE(ZQMPI(:,1:NIX))*KIND(ZQMPI(:,1:NIX))/4,MPI_REAL,NRANK+1,0,NCOMM,IERR)
          !
        ELSE IF (NRANK==NPROC-1) THEN
          CALL MPI_SEND(ZQMPI(:,(NIX*(ILOCNIY-1)+1):NIX*ILOCNIY),SIZE(ZQMPI(:,(NIX*(ILOCNIY-1)+1):NIX*ILOCNIY))*&
          KIND(ZQMPI(:,(NIX*(ILOCNIY-1)+1):NIX*ILOCNIY))/4,MPI_REAL,NRANK-1,0,NCOMM,IERR)
          !
          CALL MPI_RECV(Z_MPI_NORTH,SIZE(Z_MPI_NORTH)*KIND(Z_MPI_NORTH)/4,MPI_REAL,NRANK-1,0,NCOMM,ISTATUS_RECV_DOWN,IERR)
          !
        ELSE
          CALL MPI_SEND(ZQMPI(:,(NIX*(ILOCNIY-1)+1):NIX*ILOCNIY),SIZE(ZQMPI(:,(NIX*(ILOCNIY-1)+1):NIX*ILOCNIY))*&
          KIND(ZQMPI(:,(NIX*(ILOCNIY-1)+1):NIX*ILOCNIY))/4,MPI_REAL,NRANK-1,0,NCOMM,IERR)
          !
          CALL MPI_RECV(Z_MPI_SOUTH,SIZE(Z_MPI_SOUTH)*KIND(Z_MPI_SOUTH)/4,MPI_REAL,NRANK+1,0,NCOMM,ISTATUS_RECV_UP,IERR)
          !
          CALL MPI_RECV(Z_MPI_NORTH,SIZE(Z_MPI_NORTH)*KIND(Z_MPI_NORTH)/4,MPI_REAL,NRANK-1,0,NCOMM,ISTATUS_RECV_DOWN,IERR)
          !
          CALL MPI_SEND(ZQMPI(:,1:NIX),SIZE(ZQMPI(:,1:NIX))*KIND(ZQMPI(:,1:NIX))/4,MPI_REAL,NRANK+1,0,NCOMM,IERR)
          !
        END IF
      END IF
      IF (LHOOK) CALL DR_HOOK('SNOWPAPPUS_MPI',1,ZHOOK_HANDLE)
#endif
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !                                             END COMUNICATION BLOCK                                                !
      !                                             BEGIN DIVERGENCE BLOCK                                                !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      !             "UPWIND SCHEME" (1D example)
      !
      !     \          |     * ====|          |
      !     \     * ===|     :     |          |
      !     \     :    |     :     |    *==== |
      !     \_____:____|_____:_____|____:_____|
      !        pixel A   pixel B      pixel C
      !
      ! The problem is that fluxes are computed at the pixel center not at the faces/side
      ! Here the aim a numerical scheme is to compute transport fluxes at pixels sides
      ! the use of an "upwind scheme" mean that the computed flux at the center IS the flux that goes thru the pixel face
      ! this is made to prevent from interpolation problems using other scheme

      !  ┬ ┬┬ ┬┌─┐┌─┐┌┬┐┬ ┬┌─┐┌─┐┬┌─┐
      !  ├─┤└┬┘├─┘│ │ │ ├─┤├┤ └─┐│└─┐
      !  ┴ ┴ ┴ ┴  └─┘ ┴ ┴ ┴└─┘└─┘┴└─┘ of divergence
      !
      ! divergence computation are done using an idealised flat surface and then reprojected accorded to each point z slope
      ! (in reality slope isn't isotropic, for a flat 250 side pixel projected on real topography, the pixel side are bigger than 250m)

      ! the 'upwind' scheme implemented here is valid only for folowing hypothesis:
      ! IMPORTANT :: The maximum snow displacement in one timestep is one pixel !

      ! link to improve scheme:
      !https://ocw.mit.edu/courses/earth-atmospheric-and-planetary-sciences/12-950-atmospheric-and-oceanic-modeling-spring-2004/lecture-notes/lec10.pdf
      !https://en.wikipedia.org/wiki/False_diffusion
      !https://pdf.sciencedirectassets.com/271589/1-s2.0-S0307904X00X0145X/1-s2.0-0307904X9390048L/main.pdf?X-Amz-Security-Token=IQoJb3JpZ2luX2VjEEkaCXVzLWVhc3QtMSJHMEUCIQCPgXnzSYqeTcjaKsp6FZAvB17CiQ9m41XYYz0BWg0EZQIgYrqQk54P9THrckvhgXqh%2BxZT5TGDP%2FHpB8V%2FCpXplrgqgwQI0v%2F%2F%2F%2F%2F%2F%2F%2F%2F%2FARAEGgwwNTkwMDM1NDY4NjUiDJeNcfoCynFVGYS5birXA0PD6hv2Bjh5ZdVZa4NwUze5vo7KM95WVtg5mxDhYT95FyqjVPNcbFrzfpe9czEuw22yELtNacG%2BozFzufdWwasnpwQLD5VINbGMeaEZCnueiYTCscoRE5AxCMoOCyVVFLhnoWqnKj%2FwBA6tUqOrp%2F2MtbaLkdjrdvLtBDLqaX1hy3E7g8ir0sBOt6X1gQju8Hyw1af1TL2LrVkatXYsUtttvFHqTFrzQcqwVSUL6sgtr0ytBuiH4BDP0Lz%2FM5RokbT8oJv3lOz16rwVDlQPABD6256MS%2FgWEzO2QZoSTFyTieYJ79Es3zEqfRteiovU%2FtxCmVcRDC8rhAiyi9lCDqeEAEOG9OQMsvADuNF3ZRJ38HsYQZxamEo5%2FgjMDZu0oONLTq1DKbCUgOC8k2jPj9H%2FXakgBgrKAterA36p%2Fz5XwieUjTW43N7LnyGdq10BPezqnrZ8ofBLmanMsQbz%2FRMXcxLhdvKGPxRT6u4RkvfPPnfrq%2F6fbWvt0OFdAKoVE0Y7vfJuLKgph0VgU3bduGxq2GdR9U4bLCkLH0rZOyoarhlQ2Z7BDG5CU6XhxB0SW0TmC5%2FkCw4xLhBiZ5vU64BtkrWv9YI864Qs23KW3pkIPJvKDgKafTCmjIuSBjqlAXY0dDu5Rbo0kb3VhLXEyHVQsTlCD7YSq8EoZnI5OdGgyZTysHY%2FcZSHXBwgOy%2FWK8c%2FVWZ8OlZXmA5FeWoKa0EFmLkwE0HikgNKYXexjAwb5GWb7RpB%2BF5nmtg%2Fgn%2Bo4UBfyf2jHbwbRqZ2yBS2A9k81qr771B0MUn5EqbuKxzHi9swE6KJ0QXh7mNaXemRSNPNMPvqIzyHMa432E2Kpv5p7nyjUw%3D%3D&X-Amz-Algorithm=AWS4-HMAC-SHA256&X-Amz-Date=20220329T092934Z&X-Amz-SignedHeaders=host&X-Amz-Expires=300&X-Amz-Credential=ASIAQ3PHCVTYRNBCFCLL%2F20220329%2Fus-east-1%2Fs3%2Faws4_request&X-Amz-Signature=d7636f5d02d5ae5020bc0868a7e832388618b3ac57d22de678eb4b191515b945&hash=e8744df9126369ca6be3f51e5d25e2b292c75fb9754d4f962f62dec93507251a&host=68042c943591013ac2b2430a89b270f6af2c76d8dfd086a07176afe7c76c2c61&pii=0307904X9390048L&tid=spdf-86ee6440-2043-49bf-bbb2-d6b7283cddd1&sid=bdae24e639b278497d891794f998452facfagxrqb&type=client&ua=52055752565a00505d&rr=6f3794f69cde3ba3
      !https://en.wikipedia.org/wiki/Iterative_Stencil_Loops
      !https://math.stackexchange.com/questions/2916234/how-to-obtain-the-9-point-laplacian-formula
      !http://acoustics.ae.illinois.edu/pdfs/lele-1992.pdf
      !https://scicomp.stackexchange.com/questions/2114/implicit-finite-difference-schemes-for-advection-equation
      !https://www.uni-muenster.de/imperia/md/content/physik_tp/lectures/ws2016-2017/num_methods_i/advection.pdf
      !https://docs.lib.purdue.edu/cgi/viewcontent.cgi?article=1928&context=cstech
      !
      IF (LHOOK) CALL DR_HOOK('SNOWPAPPUS_DIV',0,ZHOOK_HANDLE)    !
      !#ifdef SFX_MPI
      IF (SIZE(PVMOD).EQ.NIX*ILOCNIY) THEN ! Reshape 1D point list to the 2D spatial domain
        !table core
        ZDOMAINE_2D(2:NIX+1,2:ILOCNIY+1,1) = RESHAPE(ZQMPI(1,:),(/NIX,ILOCNIY/)) ! QSALT
        ZDOMAINE_2D(2:NIX+1,2:ILOCNIY+1,2) = RESHAPE(ZQMPI(2,:),(/NIX,ILOCNIY/)) ! QSUSP
        ZDOMAINE_2D(2:NIX+1,2:ILOCNIY+1,3) = RESHAPE(ZQMPI(3,:),(/NIX,ILOCNIY/)) ! RAW PVDIR
        !southern side (MPI)
        ZDOMAINE_2D(2:NIX+1,1,1) = Z_MPI_SOUTH(1,:) ! MPI SOUTH QSALT
        ZDOMAINE_2D(2:NIX+1,1,2) = Z_MPI_SOUTH(2,:) ! MPI SOUTH QSUSP
        ZDOMAINE_2D(2:NIX+1,1,3) = Z_MPI_SOUTH(3,:) ! MPI SOUTH RAW PVDIR
        !northern side (MPI)
        ZDOMAINE_2D(2:NIX+1,ILOCNIY+2,1) = Z_MPI_NORTH(1,:) ! MPI NORTH QSALT
        ZDOMAINE_2D(2:NIX+1,ILOCNIY+2,2) = Z_MPI_NORTH(2,:) ! MPI NORTH QSUSP
        ZDOMAINE_2D(2:NIX+1,ILOCNIY+2,3) = Z_MPI_NORTH(3,:) ! MPI NORTH RAW PVDIR
        ! no east and west mpi boundary using SET_NB_POINTS_LIN() domain cutting

        !COMPUTE DIVERGENCE FOR SALTATION
        ZDIV_2D(2:NIX+1,2:ILOCNIY+1,1)= -MAX(-ZDOMAINE_2D(2:NIX+1,2:ILOCNIY+1,1)*COS(ZDOMAINE_2D(2:NIX+1,2:ILOCNIY+1,3)),0.)& ! OUT OF NORTH SIDE
        -MAX(ZDOMAINE_2D(2:NIX+1,2:ILOCNIY+1,1)*COS(ZDOMAINE_2D(2:NIX+1,2:ILOCNIY+1,3)),0.)& !OUT OF SOUTH SIDE
        -MAX(ZDOMAINE_2D(2:NIX+1,2:ILOCNIY+1,1)*SIN(ZDOMAINE_2D(2:NIX+1,2:ILOCNIY+1,3)),0.)& !OUT OF WEST SIDE
        -MAX(-ZDOMAINE_2D(2:NIX+1,2:ILOCNIY+1,1)*SIN(ZDOMAINE_2D(2:NIX+1,2:ILOCNIY+1,3)),0.)& !OUT OS EAST SIDE
        +MAX(ZDOMAINE_2D(2:NIX+1,3:ILOCNIY+2,1)*COS(ZDOMAINE_2D(2:NIX+1,3:ILOCNIY+2,3)),0.)& ! NORTH
        +MAX(ZDOMAINE_2D(3:NIX+2,2:ILOCNIY+1,1)*SIN(ZDOMAINE_2D(3:NIX+2,2:ILOCNIY+1,3)),0.)& !EAST
        +MAX(-ZDOMAINE_2D(2:NIX+1,1:ILOCNIY,1)*COS(ZDOMAINE_2D(2:NIX+1,1:ILOCNIY,3)),0.)& !SOUTH
        +MAX(-ZDOMAINE_2D(1:NIX,2:ILOCNIY+1,1)*SIN(ZDOMAINE_2D(1:NIX,2:ILOCNIY+1,3)),0.) !WEST

        !COMPUTE DIVERGENCE FOR SUSPENSION
        ZDIV_2D(2:NIX+1,2:ILOCNIY+1,2)= -MAX(-ZDOMAINE_2D(2:NIX+1,2:ILOCNIY+1,2)*COS(ZDOMAINE_2D(2:NIX+1,2:ILOCNIY+1,3)),0.)& ! OUT OF NORTH SIDE
        -MAX(ZDOMAINE_2D(2:NIX+1,2:ILOCNIY+1,2)*COS(ZDOMAINE_2D(2:NIX+1,2:ILOCNIY+1,3)),0.)& !OUT OF NORTH SIDE
        -MAX(ZDOMAINE_2D(2:NIX+1,2:ILOCNIY+1,2)*SIN(ZDOMAINE_2D(2:NIX+1,2:ILOCNIY+1,3)),0.)& !OUT OF WEST SIDE
        -MAX(-ZDOMAINE_2D(2:NIX+1,2:ILOCNIY+1,2)*SIN(ZDOMAINE_2D(2:NIX+1,2:ILOCNIY+1,3)),0.)& !OUT OF EAST SIDE
        +MAX(ZDOMAINE_2D(2:NIX+1,3:ILOCNIY+2,2)*COS(ZDOMAINE_2D(2:NIX+1,3:ILOCNIY+2,3)),0.)& ! NORTH
        +MAX(ZDOMAINE_2D(3:NIX+2,2:ILOCNIY+1,2)*SIN(ZDOMAINE_2D(3:NIX+2,2:ILOCNIY+1,3)),0.)& !EAST
        +MAX(-ZDOMAINE_2D(2:NIX+1,1:ILOCNIY,2)*COS(ZDOMAINE_2D(2:NIX+1,1:ILOCNIY,3)),0.)& !SOUTH
        +MAX(-ZDOMAINE_2D(1:NIX,2:ILOCNIY+1,2)*SIN(ZDOMAINE_2D(1:NIX,2:ILOCNIY+1,3)),0.) !WEST

        !RESHAPING THE COMPUTED DIVERGENCE FIELDS IN 1D
        ZQDEP_SALT(:)=RESHAPE(ZDIV_2D(2:NIX+1,2:ILOCNIY+1,1),(/NIX*ILOCNIY/))
        ZQDEP_SUSP(:)=RESHAPE(ZDIV_2D(2:NIX+1,2:ILOCNIY+1,2),(/NIX*ILOCNIY/))

      ELSE IF (SIZE(PVMOD).NE.NIX*ILOCNIY) THEN ! erreur dans les declarations
        WRITE(*,*)"BIG ERROR IN DIVERGENCE => MESSAGE FROM SNOWPAPPUS.F90"
        STOP
      END IF
      !END IF

      IF (LHOOK) CALL DR_HOOK('SNOWPAPPUS_DIV',1,ZHOOK_HANDLE)
      !#endif

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !                                           END DIVERGENCE BLOCK                                                    !
      !                                           BEGIN EROSION BLOCK                                                     !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      IF (LHOOK) CALL DR_HOOK('SNOWPAPPUS_ERODE',0,ZHOOK_HANDLE)
      !
      !0.computation of deposition rate as a function of the selected option
      !
      IF (HSNOWPAPPUSERODEPO == 'ERO') THEN
        ! represents a "pure erosion case" assuming computed blowing snow flux is built over a distance of 250 meters
        ! rq : this option is mainly for numerical test purpose and does not necessarily represents a realistic case
        ZQDEP_TOT(:) = - PQT_TOT(:)/250. - ZQSUBRATE(:)
        !
      ELSEIF (HSNOWPAPPUSERODEPO == 'DEP') THEN
        ! represents a systematic deposition of the same strength as the erosion in the above option
        ! rq : this option is mainly for numerical test purpose and does not necessarily represents a realistic case
        ZQDEP_TOT(:) = PQT_TOT(:)/250. - ZQSUBRATE(:)

      ELSEIF (HSNOWPAPPUSERODEPO == 'DIV') THEN
        ! Deposition rate is taken as the result of the divergence block
        ZQDEP_TOT(:) =(ZQDEP_SALT(:)+ZQDEP_SUSP(:))*PDIRCOSZW(:)/(PMESH_SIZE) - ZQSUBRATE(:)
        !
        !/!\ The given deposition rate is in surfacic mass, the surface being the surface occupied by the surface of the
        ! sloping surface
        !/!\ Here a big simplification is done : it is considered that the transport rate QT_TOT which are computed
        ! are the ones applying to a projection of the 250x 250 m pixel over an horizontal plane. This allows to close
        ! the mass balance. Then the computed deposition rate is projected  along the normal vector of the slope.
        ! This is not physical and should ideally be changed for complex environments.
        !
      ELSEIF (HSNOWPAPPUSERODEPO == 'NON') THEN
        ZQDEP_TOT(:) = 0.

      ENDIF
      !
      !Limiting the ero/dep flux to the minimum flux alowed by SNOWCRO:SNOWNLFALL_UPGRID
      IF (OPAPPULIMTFLUX) THEN
        WHERE (ZQDEP_TOT(:)<XUEPSI .AND. ZQDEP_TOT(:)>-XUEPSI)
          ZQDEP_TOT(:) = 0.
        ENDWHERE
      ENDIF

      !Sortie de ZQDEP_TOT en diag
      PAPPUS_DEBUG(:,4)=ZQDEP_TOT(:)
      !
      !1.ablation/erosion of snow
      !
      DO JJ=1, SIZE(PVMOD)! Loop on spatial points
        !
        !1.1 case when there is snow net ablation
        !
        ZSNOWDINNER=0.

        IF (ZQDEP_TOT(JJ)<0. .AND. ZWSNOW_T(JJ)>0.) THEN

          ZSNOWDINNER = -ZQDEP_TOT(JJ)*PTSTEP !conversion kg/m2/s => kg/m2

          DO WHILE (ZSNOWDINNER>0. .AND. INLVLS_USE(JJ)>0.)

            IF (ZSNOWDINNER<ZSNOWDRYMASS(JJ,1)) THEN ! if less than one dry mass layer is to be removed

              CALL SNOWAPPETIZER(ZSNOWDINNER,ZSNOWRHO(JJ,1),ZSNOWHEAT(JJ,1),ZSNOWDRYMASS(JJ,1),ZSNOWLIQ(JJ,1),ZSNOWSWE(JJ,1),&
              ZSNOWDZ(JJ,1),ZSCAP(JJ,1),ZSNOWTEMP(JJ,1))! supprime une partie de couche
              ZSNOWDINNER=0. ! pb ca devrait etre

            ELSE IF (ZSNOWDINNER >= ZSNOWDRYMASS(JJ,1)) THEN ! if a complete layer has to be removed

              ZSNOWDINNER=ZSNOWDINNER-ZSNOWDRYMASS(JJ,1) !
              CALL SNOWEATER(ZSNOWHEAT(JJ,:), ZSNOWDZ(JJ,:),ZSNOWSWE(JJ,:),ZSNOWRHO(JJ,:),ZSNOWSPHERI(JJ, :),&
              ZSNOWDIAMOPT(JJ,:),ZSNOWHIST(JJ,:),ZSNOWAGE(JJ,:),ZSNOWDRYMASS(JJ,:), ZSNOWLIQ(JJ,:),&
              ZSCAP(JJ,:), ZSNOWTEMP(JJ,:),PTHRUFAL(JJ), PTSTEP, INLVLS_USE(JJ))

            ENDIF
          ENDDO

          ! alerts and records when there is less snow than what must be removes ( breaks the mass conservation in the grid )
          IF (ZSNOWDINNER>0.) THEN
            IF (OPAPPUDEBUG) THEN
              write(*,*)'Warning SnowPappus :', ZSNOWDINNER, 'kg/m2 of snow were not removed on point', JJ
            ENDIF
            PSNOWDEBTC(JJ)=ZSNOWDINNER/PTSTEP !conversion kg/m2 => kg/m2/s
          ENDIF
          !
          !1.2 case when there is net snow deposition
          !
        ELSEIF (ZQDEP_TOT(JJ)>0) THEN !si déposition sup a XUEPSI (condition limite sur chute de neige dans SNOCRO : SNOWNLFALL_UPGRID )

          PBLOWSNW(JJ,1) = ZQDEP_TOT(JJ)     ! deposition flux (kg/m2/s)
          PBLOWSNW(JJ,2) = XRHODEPPAPPUS     ! density of deposited snow (kg/m3)
          PBLOWSNW(JJ,3) = XDIAMDEPPAPPUS     !optical diameter of deposited snow (m)
          PBLOWSNW(JJ,4) = XSPHDEPPAPPUS     ! sphericity of deposited snow

        END IF
        ! following is a technical condition for closing the mass bilan when
        ! prescribing ZQDEP_TOT<0 when no snow is on the ground
        IF (ZQDEP_TOT(JJ)<0 .AND. ZWSNOW_T(JJ)==0. .AND. ZSNOWDINNER==0. ) THEN
          PSNOWDEBTC(JJ)=-ZQDEP_TOT(JJ) ! already in kg/m2/s
        ENDIF

        IF (OPAPPUDEBUG) THEN
          !the aim is to check mass convervation in snowpappus.
          !we compare swe before and after computation and look for error > 1E-10
          ! error < 1E-10 are considered rounding errors
          IF (ZQDEP_TOT(JJ)<=0) THEN
            IF (abs(sum(PEK%TSNOW%WSNOW(JJ,:))+(ZQDEP_TOT(JJ)*900)+PSNOWDEBTC(JJ)*900-sum(ZSNOWSWE(JJ,:)))<1E-10) THEN
              !              WRITE(*,*) 'bILAN ero OK'
            ELSE
              write(*,*)'Bilan ero FAUX :', sum(PEK%TSNOW%WSNOW(JJ,:))+(ZQDEP_TOT(JJ)*900)+PSNOWDEBTC(JJ)*900&
              -sum(ZSNOWSWE(JJ,:))
              write(*,*)'detail', sum(ZSNOWSWE(JJ,:)),ZQDEP_TOT(JJ)*900,PSNOWDEBTC(JJ)*900,sum(PEK%TSNOW%WSNOW(JJ,:))
            ENDIF
          ELSE
            IF(abs(sum(ZSNOWSWE(JJ,:))-sum(PEK%TSNOW%WSNOW(JJ,:)))<1E-10) THEN
              !             WRITE(*,*) 'bILAN depo OK'
            ELSE
              write(*,*)'Bilan depo FAUX:',sum(ZSNOWSWE(JJ,:))-sum(PEK%TSNOW%WSNOW(JJ,:))
            ENDIF
          ENDIF
        ENDIF
      END DO
      ! end pappus errode

      IF (LHOOK) CALL DR_HOOK('SNOWPAPPUS_ERODE',1,ZHOOK_HANDLE)

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !                                       END EROSION BLOCK                                                           !
      !                                                                                                                   !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      ! writes prognstic variables in PEK
      PEK%TSNOW%WSNOW(:,:)=ZSNOWSWE(:,:)
      PEK%TSNOW%HEAT(:,:)=ZSNOWHEAT(:,:)
      PEK%TSNOW%RHO(:,:)=ZSNOWRHO(:,:)
      PEK%TSNOW%HIST(:,:)=ZSNOWHIST(:,:)
      PEK%TSNOW%AGE(:,:)=ZSNOWAGE(:,:)
      PEK%TSNOW%SPHERI(:,:)=ZSNOWSPHERI(:,:)
      PEK%TSNOW%DIAMOPT(:,:)=ZSNOWDIAMOPT(:,:)

    END IF
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !                                       END 2D PART OF SNOWPAPPUS (IF)                                              !
      !                                             ENDING BLOCK                                                          !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    IF (LHOOK) CALL DR_HOOK('SNOWPAPPUS',1,ZHOOK_HANDLE2)

  END SUBROUTINE SNOWPAPPUS

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !   _____                     _____
  !  / ____|                   |  __ \
  ! | (___  _ __   _____      _| |__) |_ _ _ __  _ __  _   _ ___
  !  \___ \| '_ \ / _ \ \ /\ / /  ___/ _` | '_ \| '_ \| | | / __|
  !  ____) | | | | (_) \ V  V /| |  | (_| | |_) | |_) | |_| \__ \
  ! |_____/|_| |_|\___/ \_/\_/ |_|   \__,_| .__/| .__/ \__,_|___/
  !                                       | |   | |
  !                                       |_|   |_|
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  SUBROUTINE SUSPENSION_PAPPUS(PVFRIC,PVFRIC_T,PHSALT_SUSP,PHMIN_SUSP,PVFALL,PLFETCH,PCSALT,PZ0,PQSUSP,PBLOWSNWFLUX_1M,&
                               PBLOWSNWFLUXINT)
    ! Subroutine used un SNOWPAPPUS.F90 used to calculates the suspension flux
    !
    !
    !
    REAL,INTENT(IN)  ::PVFRIC!wind friction speed (m/s)
    REAL,INTENT(IN)  ::PVFRIC_T! thresold wind friction speed to trigger saltation phenomenon (m/s)
    REAL,INTENT(IN)  ::PHSALT_SUSP ! HAUTEUR DE SALTATION 2  (m)
    REAL,INTENT(IN)  ::PHMIN_SUSP ! height of the lower limit of flux taken into account in suspension flux calculation (m)
    REAL,INTENT(IN)  ::PVFALL ! fall velocity (m/s)
    REAL,INTENT(IN)  ::PLFETCH
    REAL,INTENT(IN)  ::PCSALT ! concentration in the saltation layer (kg/m3)
    REAL,INTENT(IN)  ::PZ0 ! roughness length for momentum
    REAL,INTENT(OUT) ::PQSUSP ! transport rate or linar flux of suspensive snow (kg/ms)
    REAL,INTENT(OUT) ::PBLOWSNWFLUX_1M ! snow flux 1 m above snow surface (kg/m2/s)
    REAL,INTENT(OUT) ::PBLOWSNWFLUXINT ! snow flux integrated between 0.2 and 1.2 m above the surface (kg/m(2)/s)
    REAL,PARAMETER    :: ZK=0.41! von Karman's constant
    REAL,PARAMETER   :: ZHMAX_INT=1.2
    REAL,PARAMETER   :: ZHMIN_INT=0.2
    REAL    :: ZEXPO!exponent of the power law describing the concentration profile
    REAL    ::ZHMAX_SUSP,ZHSALT_SUSP !maximum reached by suspended particles

    IF (PVFRIC>=PVFRIC_T) THEN ! checks if transport occurs
      ZEXPO=PVFALL/(ZK*PVFRIC)
      !maximum heigth of the suspension layer  (simplification of Pomeroy 1993 formula)
      ZHMAX_SUSP=PHSALT_SUSP+PLFETCH*ZK**2/SQRT(LOG(5./PZ0)*LOG(PHSALT_SUSP/PZ0))
     !calculation of the suspension mass flux, in a monodisperse hypothesis

      ZHSALT_SUSP=MAX(PZ0,PHSALT_SUSP)! /!\ Blindage pour eviter le cas ou hsalt < z0 qui entraire des flux negatifs

      IF (ZEXPO /= 1) THEN
        ! DERIVATION DE LA SUSPENTION PAR ANGE TENTATIVE DE CORRECTION POUR DES SORTIES QSUSP NEGATIVE
        ! PQSUSP=PCSALT*(1/PHSALT_SUSP)**(-ZEXPO)*(PVFRIC/ZK)&
        ! *( ( ( LOG(ZHMAX_SUSP/PZ0)*ZHMAX_SUSP**(1-ZEXPO)/ (1-ZEXPO)) - (LOG(PHSALT_SUSP/PZ0)*PHSALT_SUSP**(1-ZEXPO)/(1-ZEXPO)))&
        ! -( ( ZHMAX_SUSP**(1-ZEXPO)/((ZEXPO-1)**2) ) - ( PHSALT_SUSP**(1-ZEXPO)/((ZEXPO-1)**2) ) ) )

        ! LE DENOMINATEUR EST NEGATIF (ZK*(1-ZEXPO)) QUAND ZEXPO >1 DONC QSUSP NEGATIF AUSSI => ERREUR
        PQSUSP = PCSALT*ZHSALT_SUSP*PVFRIC*(((ZHMAX_SUSP/ZHSALT_SUSP)**(-ZEXPO+1)*(LOG(ZHMAX_SUSP/PZ0)-1/(1-ZEXPO))) &
        -((PHMIN_SUSP/ZHSALT_SUSP)**(-ZEXPO+1)*(LOG(PHMIN_SUSP/PZ0)-1/(1-ZEXPO))))/(ZK*(1-ZEXPO))
        PBLOWSNWFLUXINT = PCSALT*ZHSALT_SUSP*PVFRIC*(((ZHMAX_INT/ZHSALT_SUSP)**(-ZEXPO+1)*(LOG(ZHMAX_INT/PZ0)-1/(1-ZEXPO))) &
        -((ZHMIN_INT/ZHSALT_SUSP)**(-ZEXPO+1)*(LOG(ZHMIN_INT/PZ0)-1/(1-ZEXPO))))/(ZK*(1-ZEXPO))
      ELSE !case when the above expression is not defined
        PQSUSP = PCSALT*ZHSALT_SUSP*PVFRIC*(LOG(ZHMAX_SUSP/PZ0)**2-LOG(PHMIN_SUSP/PZ0)**2)/(2*ZK)
        PBLOWSNWFLUXINT = PCSALT*ZHSALT_SUSP*PVFRIC*(LOG(ZHMAX_INT/PZ0)**2-LOG(ZHMIN_INT/PZ0)**2)/(2*ZK)
      ENDIF
      PBLOWSNWFLUX_1M=PCSALT*(1./PHSALT_SUSP)**(-ZEXPO)*PVFRIC*LOG(1./PZ0)/ZK
    ENDIF

  END SUBROUTINE SUSPENSION_PAPPUS

  SUBROUTINE VFALL_SUSPENSION(PVFRIC,PSPHERI,PDIAMOPT,PSNOWAGE,HLIMVFALL,PVFALL)
    USE MODE_SNOW3L,ONLY: CHECK_DENDRITIC, GETDENDRICITY
    USE MODD_SNOW_METAMO, ONLY: XVDIAM6, XUEPSI  !parameters necessary to determine if snow is dendritic
    uSE MODD_SNOW_PAR, ONLY: XAGELIMPAPPUS2, XDEMAXVFALL
    !parametrization of the mass-averaged terminal fall speed of suspended particles ( multiplied by schmidt number )
    REAL,INTENT(IN)    ::PVFRIC ! friction velocity (m/s)
    REAL,INTENT(IN)    ::PSPHERI !snow sphericity
    REAL,INTENT(IN)    ::PDIAMOPT !snow optical diameter (m)
    REAL,INTENT(IN)    ::PSNOWAGE ! snow age (day)
    CHARACTER(4),INTENT(IN)::HLIMVFALL
    REAL,INTENT(OUT)   ::PVFALL ! suspeded particles terminal fall speed (m/s)
    !local variables
    LOGICAL            ::GCOND ! bool for deciding if it is old or new snow, ! may be used uninitialized
    LOGICAL            ::ZDENDRITRUE ! bool beeing true if snow is dendritic
    CHARACTER(3)       ::HVFALL='MIX'
    REAL               ::ZDENDRICITY ! snow dendricity ( old metamorphism )
    REAL               ::ZF ! factor
    REAL               ::ZVFALLY ! terminal fall speed of young snow
    REAL               ::ZVFALLO ! terminal fall speed of old snow
    ZVFALLY=MIN(0.38*PVFRIC+0.12,0.8)
    ZVFALLO=0.8
    CALL CHECK_DENDRITIC(PDIAMOPT,PSPHERI,ZDENDRITRUE)
    IF (HLIMVFALL =='DEND') THEN
      GCOND=ZDENDRITRUE
    ELSEIF (HLIMVFALL == 'PREC'.OR. HLIMVFALL =='MIXT') THEN
      GCOND=(PSNOWAGE<XAGELIMPAPPUS2 .AND. ZDENDRITRUE)
    ENDIF
    IF (HVFALL=='MIX')THEN
      IF (GCOND) THEN
        PVFALL=ZVFALLY
      ELSE
        IF (ZDENDRITRUE .AND. HLIMVFALL=='MIXT') THEN
          CALL GETDENDRICITY(PDIAMOPT,PSPHERI,ZDENDRICITY)
          ZF=MIN(1.,ZDENDRICITY/XDEMAXVFALL)
          PVFALL = ZVFALLO*(1-ZF)+(ZF*ZVFALLY)
        ELSE
          PVFALL= ZVFALLO
        ENDIF
      !(A faire : mettre une option namelist pour pouvoir jouer sur cette paramétrisation)
      ENDIF
      !
    ELSEIF (HVFALL=='OLD') THEN
      PVFALL=0.8
    !
    ELSEIF (HVFALL=='N98') THEN
      PVFALL=0.38*PVFRIC+0.12
    !
    ENDIF
  END SUBROUTINE VFALL_SUSPENSION


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE SALTATION_P90_P92(PVFRIC,PVFRIC_T,PRHOA,PQSALT,PCSALT,PHSALT1,PHSALT2)
    ! Subroutine used in SNOWPAPPUS.F90 to compute the saltations parameters
    ! don't forget to check the parameters
    ! using equations from POMEROY 90 and 92
    !
    REAL,INTENT(IN)  ::PVFRIC !wind friction speed (M/s)
    REAL,INTENT(IN)  ::PVFRIC_T ! thresold wind friction speed to trigger saltation phenomenon (M/s)
    REAL,INTENT(IN)  ::PRHOA ! air density (kg/m3)
    REAL,INTENT(OUT) ::PHSALT1 ! Saltation height 1 as defined in POMEROY 90 paper (M)
    REAL,INTENT(OUT) ::PHSALT2 ! Saltation height 2 as defined in POMEROY 92 paper (M)
    REAL,INTENT(OUT) ::PQSALT ! Transport rate or Linear saltating flux (kg/ms)
    REAL,INTENT(OUT) ::PCSALT ! Concentration in the saltation layer (kg/m3)
    REAL,PARAMETER :: PVFNEROD=0. ! NON ERRODIBLE FRICTION VELOCITY AS DEFINED IN POMEROY90 U*n (M/S) DEFAULT =0
    REAL,PARAMETER :: PPSALTEFFIC=.68 != (C . E )* U* = SALTATION EFFICIENCY COEFFICIENT AS DEFINED IN POMEROY90 (DEFAULT SET AS 0.68 /U*) (M/S), neeed to divide by PVFRIC in the expression
    REAL,PARAMETER :: ZG=9.80665 ! ZG = STANDARD GRAVITY ACCELERATION (M.S-2)

    PHSALT1=(1.6*PVFRIC*PVFRIC)/(2.*ZG) ! hsalt 1 from Pomeroy 1990
    PHSALT2=0.08436*(PVFRIC**1.27) ! hasalt2 from Pomeroy 1992 => taken as the top of the saltation layer here

    PQSALT = (PPSALTEFFIC/PVFRIC) * (PRHOA/ZG)*((PVFRIC_T*PVFRIC**2) - (PVFRIC_T*PVFNEROD**2) -(PVFRIC_T**3) )*(PHSALT2/PHSALT1)
    ! here, we have Qsalt from Pomeroy 1990, multiplied by the ratio hsalt2/hsalt1, see Baron, Haddjeri et al. 2023 for more details

    PCSALT = PQSALT/(PHSALT2*2.8*PVFRIC_T)!Concentration in the saltation layer, as in Pomeroy & Gray 1990
  END SUBROUTINE SALTATION_P90_P92

  SUBROUTINE SALTATION_S04(PVFRIC,PVFRIC_T,PRHOA,PQSALT)
    !computes saltation transport according to Sorensen 2004 parameterization refined by Vionnet et al. 2012 ( PhD Thesis )
    REAL,INTENT(IN)  ::PVFRIC !wind friction speed (m/s)
    REAL,INTENT(IN)  ::PVFRIC_T ! thresold wind friction speed to trigger saltation phenomenon (m/s)
    REAL,INTENT(IN)  ::PRHOA ! air density (kg/m3)
    REAL,INTENT(OUT) ::PQSALT !vertically integrated total saltation flux (kg/m/s)
    REAL             ::ZV
    REAL,PARAMETER   :: ZG=9.80665 ! ZG = STANDARD GRAVITY ACCELERATION (M.S-2)
    REAL,PARAMETER   :: A=2.6 ! A, B and C are constants for the Sorensen 2004 formula
    REAL,PARAMETER   :: B=2.5
    REAL,PARAMETER   :: C=2
    ZV=PVFRIC/PVFRIC_T
    PQSALT = (PRHOA*PVFRIC**3/ZG)*(1-1/ZV**2)*(A+B/ZV**2+C/ZV)
  END SUBROUTINE SALTATION_S04


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  SUBROUTINE WINDSPEEDTHRESOLD (PDIAMOPT,PSPHERI,PSNOWAGE,HSNOWMOB,PSNOWRHO,PZ0,GTRANSPORTABLE,KNB_PTS,PVFRIC_T)
    !Subroutine calculating threshold wind velocity for transport ( if snow is considered "transportable" )
    USE MODE_SNOW3L, ONLY : CHECK_DENDRITIC,GETDENDRICITY,GETGRAINSIZE_B21
    USE MODD_SNOW_PAR, ONLY : XAGELIMPAPPUS, XWINDTHRFRESH

    IMPLICIT NONE
    INTEGER, INTENT(IN)            ::KNB_PTS
    REAL, DIMENSION(KNB_PTS), INTENT(IN)  :: PDIAMOPT ! optical diameter of the surface snow layer (m)
    REAL, DIMENSION(KNB_PTS), INTENT(IN)  :: PSPHERI ! sphericity of the surface snow layer (m)
    REAL, DIMENSION(KNB_PTS), INTENT(IN)  :: PSNOWRHO ! density of the surface snow layer (m)
    REAL, DIMENSION(KNB_PTS), INTENT(IN)  :: PZ0 ! surface roughness length
    REAL, DIMENSION(KNB_PTS), INTENT(IN)  :: PSNOWAGE ! age of the surface snow layer ( day )
    CHARACTER(4),INTENT(IN):: HSNOWMOB ! option for wind speed threshold calculation
    LOGICAL, DIMENSION(KNB_PTS), INTENT(IN):: GTRANSPORTABLE !Logical telling if the snow is transportable on the point
    REAL, DIMENSION(KNB_PTS), INTENT(OUT) :: PVFRIC_T ! threshold friction velocity

    !LOCAL VARIABLE
    REAL              :: ZDRIFTINDEX ! driftability index ( see definition in Vionnet et al. 2012 )
    REAL, PARAMETER   :: ZRHOSNOWMIN = 50. ! minimum density (kg/m3) for HSNOWMOB='VI12' option, see Vionnet et al. 2012
    REAL              :: ZDENDRICITY ! snow surface dendricity ( old version of microstructure representation, see Vionnet 2012)
    REAL              :: ZGRAINSIZE ! snow surface grain size ( """"" )
    REAL              :: PVFRIC_T5 ! wind speed threshold at 5 meter heigth
    REAL              :: XWINDTHR = 9. ! constant threshold wind velocity used when HSNOWMOB='CONS'
    LOGICAL           :: GDENDRITRUE ! logical testing if the considered layer is dendritic
    INTEGER           :: JJ ! loop indexer (on points, spatial)
    IF (HSNOWMOB=='CONS') THEN
      DO JJ=1,SIZE(PVFRIC_T)
        IF (GTRANSPORTABLE(JJ)) THEN
          ! case 1 : young snow => fixed threshold wind speed
          IF (PSNOWAGE(JJ)<XAGELIMPAPPUS) THEN
            PVFRIC_T5 = XWINDTHRFRESH
          ! case 2 : not young snow => threshold wind speed depends on HSNOWMOB option
          ELSE
            PVFRIC_T5 = XWINDTHR
          ENDIF
        ENDIF
        PVFRIC_T(JJ)=0.4*PVFRIC_T5/LOG(5/PZ0(JJ)) ! calcul de vfric_t au sol
      ENDDO
    ELSEIF (HSNOWMOB=='COGM') THEN
      DO JJ=1,SIZE(PVFRIC_T)
        IF (GTRANSPORTABLE(JJ)) THEN
          IF (PSNOWAGE(JJ)<XAGELIMPAPPUS) THEN
            PVFRIC_T5 = XWINDTHRFRESH
          ELSE
            CALL CHECK_DENDRITIC(PDIAMOPT(JJ),PSPHERI(JJ),GDENDRITRUE)
            IF (GDENDRITRUE)THEN
              ! calculation of driftability index from snow surface microstructure
              CALL GETDENDRICITY(PDIAMOPT(JJ),PSPHERI(JJ),ZDENDRICITY)
              ZDRIFTINDEX=.75*ZDENDRICITY-.5*PSPHERI(JJ)+.5
              ZDRIFTINDEX=MAX(-0.9999,ZDRIFTINDEX)!limits the drift index to avoid log(0)=>limits 5m threshold to ~ 120 m/s
              ! calculation of 5m threshold wind speed from driftability index
              PVFRIC_T5=-11.7*LOG((1+ZDRIFTINDEX)/2.868)
            ELSE
              PVFRIC_T5 = XWINDTHR
            ENDIF
          ENDIF
          ! caculation of threshold wind friction velocity from 5m wind speed
          PVFRIC_T(JJ)=0.4*PVFRIC_T5/LOG(5/PZ0(JJ)) ! calcul de vfric_t au sol
        ENDIF
      ENDDO
    ELSEIF (HSNOWMOB=='LI07') THEN
      DO JJ=1,SIZE(PVFRIC_T)
        IF (GTRANSPORTABLE(JJ)) THEN
          IF (PSNOWAGE(JJ)<XAGELIMPAPPUS) THEN
            PVFRIC_T5 = XWINDTHRFRESH
            PVFRIC_T(JJ)=0.4*PVFRIC_T5/LOG(5/PZ0(JJ))
          ELSEIF (PSNOWRHO(JJ)<=300) THEN
            PVFRIC_T(JJ)=0.1*EXP(0.003*PSNOWRHO(JJ))
          ELSE
            PVFRIC_T(JJ)=0.005*EXP(0.013*PSNOWRHO(JJ))
          ENDIF
        ENDIF
      ENDDO
    ELSEIF (HSNOWMOB=='VI12') THEN
      DO JJ=1,SIZE(PVFRIC_T)
        IF (GTRANSPORTABLE(JJ)) THEN
          IF (PSNOWAGE(JJ)<XAGELIMPAPPUS) THEN
            PVFRIC_T5 = XWINDTHRFRESH
          ELSE
            ! calculation of driftability index from snow surface microstructure
            CALL CHECK_DENDRITIC(PDIAMOPT(JJ),PSPHERI(JJ),GDENDRITRUE)
            IF (GDENDRITRUE)THEN
              CALL GETDENDRICITY(PDIAMOPT(JJ),PSPHERI(JJ),ZDENDRICITY)
              ZDRIFTINDEX=.34*(.75*ZDENDRICITY-.5*PSPHERI(JJ)+.5)+.66*(1.25-.0042*(PSNOWRHO(JJ)-ZRHOSNOWMIN))
            ELSE
              CALL GETGRAINSIZE_B21(PDIAMOPT(JJ),PSPHERI(JJ),ZGRAINSIZE)
              ZDRIFTINDEX=.34*(-.583*ZGRAINSIZE-.833*PSPHERI(JJ)+.833)+.66*(1.25-.0042*(PSNOWRHO(JJ)-ZRHOSNOWMIN))
              ! error in vionnet 2012 not -.583
            ENDIF
            ZDRIFTINDEX=MAX(-0.9999,ZDRIFTINDEX)!limits the drift index to avoid log(0)=>limits 5m threshold to ~ 120 m/s

            PVFRIC_T5=-11.7*LOG((1+ZDRIFTINDEX)/2.868) ! calcul de VFRIC_T
          ENDIF
          ! caculation of threshold wind friction velocity from 5m wind speed
          PVFRIC_T(JJ)=0.4*PVFRIC_T5/LOG(5/PZ0(JJ)) ! calcul de vfric_t au sol
        ENDIF
      ENDDO
    ELSEIF (HSNOWMOB=='GM98') THEN
      DO JJ=1,SIZE(PVFRIC_T)
        IF (GTRANSPORTABLE(JJ)) THEN
          IF (PSNOWAGE(JJ)<XAGELIMPAPPUS) THEN
            PVFRIC_T5 = XWINDTHRFRESH
          ELSE
            CALL CHECK_DENDRITIC(PDIAMOPT(JJ),PSPHERI(JJ),GDENDRITRUE)
            ! calculation of ZDRIFTINDEX en fonction du type de neige
            IF (GDENDRITRUE)THEN
              CALL GETDENDRICITY(PDIAMOPT(JJ),PSPHERI(JJ),ZDENDRICITY)
              ZDRIFTINDEX=.75*ZDENDRICITY-.5*PSPHERI(JJ)+.5
            ELSE
              CALL GETGRAINSIZE_B21(PDIAMOPT(JJ),PSPHERI(JJ),ZGRAINSIZE)
              ZDRIFTINDEX=-.583*ZGRAINSIZE-.833*PSPHERI(JJ)+.833
              ! error in vionnet 2012 not -.583
            ENDIF
            ZDRIFTINDEX=MAX(-0.9999,ZDRIFTINDEX)!limiting the drift index in order to avoid log(0), limits wind speed threshold to 120 m/s

            PVFRIC_T5=-11.7*LOG((1+ZDRIFTINDEX)/2.868)
          ENDIF
          PVFRIC_T(JJ)=0.4*PVFRIC_T5/LOG(5/PZ0(JJ))
        ENDIF
      ENDDO
    ENDIF

  END SUBROUTINE WINDSPEEDTHRESOLD

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    SUBROUTINE SNOWAPPETIZER(PEROD_SURFACICMASS,PSNOWRHO,PSNOWHEAT,PSNOWDRYMASS,PSNOWLIQ,PSNOWSWE,PSNOWDZ,PSCAP,PSNOWTEMP)
      !Routine to remove dry part of a snow layer
      ! (no complete layer removal is made here)
      !
      !USE MODE_SNOW3L,   ONLY : SNOW3LSCAP
      !USE MODD_SNOW_PAR, ONLY : XSNOWDMIN
      USE MODD_CSTS,     ONLY : XLMTT, XTT, XCI,XRHOLW
      !
      IMPLICIT NONE
      !
      REAL, INTENT(IN) :: PEROD_SURFACICMASS ! snow SURFACIC mass to eat KG/M^2
      REAL, INTENT(INOUT) :: PSNOWRHO ! snow layer density
      REAL, INTENT(INOUT) :: PSNOWHEAT ! WARNING : Heat content of the snow layer (J/m3)
      REAL, INTENT(INOUT) :: PSNOWDRYMASS ! SNOW DRY MASS KG/M2
      REAL, INTENT(INOUT) :: PSNOWLIQ ! snow liquid mass (kg/m2)
      REAL, INTENT(OUT) :: PSNOWSWE ! Snow water equivalent of the snow layer ( kg/m2 )
      REAL, INTENT(OUT) :: PSNOWDZ ! snow layer thikness (m)
      REAL, INTENT(OUT) :: PSCAP ! snow heat capacity (J/K)
      REAL, INTENT(OUT) :: PSNOWTEMP
      !
      !1. computation snow layer temperature ,energy and width to be removed
      !
      ! HYPOTHESIS WARNING:
      !  ERODING ONLY DRY SNOW !!!
      !  consequences accumulate water in current layer
      !
      !COMPUTE NEW SNOWDRY MASS BY REMOVING PEROD_SURFACICMASS
      PSNOWDRYMASS = PSNOWDRYMASS - PEROD_SURFACICMASS
      PSNOWDZ = PSNOWDRYMASS/PSNOWRHO
      PSNOWRHO = (PSNOWDRYMASS + PSNOWLIQ*XRHOLW)/PSNOWDZ
      PSNOWSWE = PSNOWRHO * PSNOWDZ
      !COMPUTE NEW PARAMETERS
      PSCAP = PSNOWRHO * XCI
      PSNOWTEMP = XTT + (( PSNOWHEAT/PSNOWDZ + XLMTT*PSNOWRHO )/PSCAP )
      !PSNOWLIQ = MAX( 0.0, PSNOWTEMP-XTT ) * ZSCAP * PSNOWDZ / (XLMTT*XRHOLW)
      !PSNOWDRYMASS = PSNOWSWE - PSNOWLIQ*XRHOLW ! NEW DRY MASS
      PSNOWHEAT = PSNOWDZ * ( PSCAP*(PSNOWTEMP-XTT) - XLMTT*PSNOWRHO ) + XLMTT * XRHOLW * PSNOWLIQ
      PSNOWTEMP = MIN( XTT, PSNOWTEMP )

    END SUBROUTINE SNOWAPPETIZER


    SUBROUTINE SNOWEATER(PSNOWHEAT,PSNOWDZ,PSNOWSWE,PSNOWRHO,PSNOWSPHERI,PSNOWDIAMOPT,PSNOWHIST,&
      PSNOWAGE,PSNOWDRYMASS, PSNOWLIQ, PSCAP, PSNOWTEMP,PTHRUFAL, PTSTEP, KNLVLS_USE)
      !Routine for removal of complete snow layer
      !Snow prognostic variables depend on the layer, their values are organized as in this example, with a N-Layer snowpack :
      ! PSNOWDZ(1) = surface layer heigth
      ! ...
      ! PSNOWDZ (I) = Ith layer heigth under surface layer
      ! ...
      ! PSNOWDZ(N) = heigth of first snow layer above ground
      ! PSNOWDZ(N+1) = 0. ( or other initialization value )
      ! ...
      ! PSNOWDZ(KNLVLS_USE) = 0.
      ! note that, while the size of the table is fixed, the number of layers in the snowpack is not ( and inferior to its size)
      USE MODD_CSTS,     ONLY : XTT, XLMTT, XRHOLW, XCI,XRHOLI
      USE MODD_SNOW_PAR, ONLY : XSNOWDMIN
      IMPLICIT NONE
      INTEGER,INTENT(INOUT)       :: KNLVLS_USE ! number of active snow layers
      INTEGER :: JLAYER !loop index on snow layers
      REAL, DIMENSION(KNLVLS_USE), INTENT(INOUT) :: PSNOWHEAT ! heat content of snow layers (J/m-3)
      REAL, DIMENSION(KNLVLS_USE), INTENT(INOUT) :: PSNOWDZ ! heigth of snow layers (m)
      REAL, DIMENSION(KNLVLS_USE), INTENT(INOUT) :: PSNOWSWE ! snow water equivalent of snow layers
      REAL, DIMENSION(KNLVLS_USE), INTENT(INOUT) :: PSNOWRHO ! density of snow layers
      REAL, DIMENSION(KNLVLS_USE), INTENT(INOUT) :: PSNOWSPHERI ! sphericity of snow layers
      REAL, DIMENSION(KNLVLS_USE), INTENT(INOUT) :: PSNOWDIAMOPT ! optical diameter of snow layers
      REAL, DIMENSION(KNLVLS_USE), INTENT(INOUT) :: PSNOWHIST ! history of snow layer (see more in snowcro.f90 => SNOWCROMETAMO)
      REAL, DIMENSION(KNLVLS_USE), INTENT(INOUT) :: PSNOWAGE ! age of snow layers
      REAL, DIMENSION(KNLVLS_USE),INTENT(INOUT)  :: PSNOWDRYMASS ! SNOW DRY MASS KG/M2
      REAL, DIMENSION(KNLVLS_USE),INTENT(INOUT)  :: PSNOWLIQ ! snow liquid mass (kg/m2)
      REAL, DIMENSION(KNLVLS_USE),INTENT(INOUT)  :: PSCAP ! snow heat capacity (J/K)
      REAL, DIMENSION(KNLVLS_USE),INTENT(INOUT)  :: PSNOWTEMP
      REAL, INTENT(INOUT)                        :: PTHRUFAL ! rate that liquid water leaves snow pack: paritioned into soil infiltration/runoff by ISBA [kg/(m2 s)]
      REAL, INTENT(IN)                           :: PTSTEP ! time step in (s)
      REAL                                       :: ZPHASE !
      REAL                                       :: ZSNOWDZ
      REAL                                       :: ZSNOWLIQ, ZSOLIDMASS_BEFORE_REFR, ZSOLIDMASS_AFTER_REFR
      !
      !1. removal of surface layer => shift of the index of all layers
      !1a. only one layer => remove dry mass liquid mass go to PTHRUFAL
      !1b. 2 or more layer move liquid mass down + remove layer

      IF (KNLVLS_USE==1) THEN
      !
      !1a remove last snow layer put liquid mass to ptstep

        PTHRUFAL = PTHRUFAL + PSNOWLIQ(1)*XRHOLW / PTSTEP ! check units

        PSNOWSWE(KNLVLS_USE)=0.0
        PSNOWRHO(KNLVLS_USE)=999.
        PSNOWDZ(KNLVLS_USE)=0.
        PSNOWSPHERI(KNLVLS_USE)=0.
        PSNOWDIAMOPT(KNLVLS_USE)=0.
        PSNOWHIST(KNLVLS_USE)=0.
        PSNOWAGE(KNLVLS_USE)=0.
        PSNOWHEAT(KNLVLS_USE)=0.
        PSCAP(KNLVLS_USE)=0.
        PSNOWLIQ(KNLVLS_USE)=0.
        PSNOWDRYMASS(KNLVLS_USE)=0.
        PSNOWTEMP(KNLVLS_USE)=0.
        KNLVLS_USE=KNLVLS_USE-1

      ELSEIF(KNLVLS_USE>=2) THEN
        ! 1b update liquid content layer
        ! add PSNOWLIQ(1) to 2nd layer and recompute parameters
        PSNOWLIQ(2) = PSNOWLIQ(2)+PSNOWLIQ(1)
        PSNOWSWE(2) = PSNOWSWE(2)+PSNOWLIQ(1)*XRHOLW
        !
        !COMPUTE NEW PARAMETERS
        ! extracted from snowcro SNOWCROREFRZ routine
        !
        ! Calculate the maximum possible refreezing
        ZPHASE = MIN( PSCAP(2)* MAX(0.0, XTT - PSNOWTEMP(2)) * PSNOWDZ(2), &
                          PSNOWLIQ(2) * XLMTT * XRHOLW )
        ! Reduce liquid content if freezing occurs:
        ZSNOWLIQ = PSNOWLIQ(2) - ZPHASE/(XLMTT*XRHOLW)
        ! Warm layer and reduce liquid if freezing occurs:
        ZSNOWDZ = MAX(XSNOWDMIN/(SIZE(PSNOWDZ(:))-1), PSNOWDZ(2))
        ! Difference with ISBA-ES: a possible cooling of current refreezing water
        !                          is taken into account to calculate temperature change
        ZSOLIDMASS_BEFORE_REFR =  ( PSNOWRHO(2) * ZSNOWDZ - &
                                ( PSNOWLIQ(2) - PSNOWLIQ(1)  ) * XRHOLW )
        ZSOLIDMASS_AFTER_REFR =  ( PSNOWRHO(2) * ZSNOWDZ - &
                               ( ZSNOWLIQ - PSNOWLIQ(1) ) * XRHOLW )
        !
        PSNOWTEMP(2) = XTT + ( PSNOWTEMP(2)-XTT )*ZSOLIDMASS_BEFORE_REFR/ZSOLIDMASS_AFTER_REFR + &
                            ZPHASE/( XCI*ZSOLIDMASS_AFTER_REFR )
        ! Density is adjusted to conserve the mass
        PSNOWRHO(2) = PSNOWSWE(2)/ZSNOWDZ
        ! keeps snow denisty below ice density
        IF ( PSNOWRHO(2)>XRHOLI ) THEN
          WRITE(*,*) "Warning SnowPappus : SNOWEATER, snow density above ice, correcting density the hard way, can leads to pb..."
          PSNOWDZ (2) = PSNOWDZ(2) * PSNOWRHO(2) / XRHOLI
          PSNOWRHO(2) = XRHOLI
        ENDIF
        PSNOWLIQ(2) = ZSNOWLIQ
        PSCAP(2) = PSNOWRHO(2) * XCI
        PSNOWHEAT(2) = PSNOWDZ(2) * ( PSCAP(2)*(PSNOWTEMP(2)-XTT) - XLMTT*PSNOWRHO(2) ) + XLMTT * XRHOLW * PSNOWLIQ(2)
        PSNOWTEMP(2) = MIN( XTT, PSNOWTEMP(2) )
        !
        !1b move layer up
        DO JLAYER=1,KNLVLS_USE-1
        !
          PSNOWSWE(JLAYER)=PSNOWSWE(JLAYER+1)
          PSNOWHEAT(JLAYER)=PSNOWHEAT(JLAYER+1)
          PSNOWDZ(JLAYER)=PSNOWDZ(JLAYER+1)
          PSNOWRHO(JLAYER)=PSNOWRHO(JLAYER+1)
          PSNOWSPHERI(JLAYER)=PSNOWSPHERI(JLAYER+1)
          PSNOWDIAMOPT(JLAYER)=PSNOWDIAMOPT(JLAYER+1)
          PSNOWHIST(JLAYER)=PSNOWHIST(JLAYER+1)
          PSNOWAGE(JLAYER)=PSNOWAGE(JLAYER+1)
          PSNOWLIQ(JLAYER)=PSNOWLIQ(JLAYER+1)
          PSCAP(JLAYER)=PSCAP(JLAYER+1)
          PSNOWDRYMASS(JLAYER)=PSNOWDRYMASS(JLAYER+1)
          PSNOWTEMP(JLAYER)=PSNOWTEMP(JLAYER+1)
        ENDDO
        !
        !1b. set to zero last layer
        !
        PSNOWSWE(KNLVLS_USE)=0.0
        PSNOWRHO(KNLVLS_USE)=999.
        PSNOWDZ(KNLVLS_USE)=0.
        PSNOWSPHERI(KNLVLS_USE)=0.
        PSNOWDIAMOPT(KNLVLS_USE)=0.
        PSNOWHIST(KNLVLS_USE)=0.
        PSNOWAGE(KNLVLS_USE)=0.
        PSNOWHEAT(KNLVLS_USE)=0.
        PSCAP(KNLVLS_USE)=0.
        PSNOWLIQ(KNLVLS_USE)=0.
        PSNOWDRYMASS(KNLVLS_USE)=0.
        PSNOWTEMP(KNLVLS_USE)=0.
        KNLVLS_USE=KNLVLS_USE-1
      ENDIF

  END SUBROUTINE SNOWEATER



  SUBROUTINE BINTANJA_98_3M (PVMOD,PUREF,PZ0,PTA,PSUBL_RATE)

    !SNOW DRIFT SUBLIMATION RATE ACCORDING TO PARAMETRISATION OF BINTANJA_98 (EQUATION 10)
    ! USE OF COEFICIENTS FROM TABLE 3 AT 3M
    !Bintanja, R. (1998). The contribution of snowdrift sublimation to the surface mass balance of Antarctica. Annals of Glaciology, 27, 251-259.
    !DOI: https://doi.org/10.3189/1998AoG27-1-251-259
    !

    IMPLICIT NONE
    REAL, INTENT(IN):: PVMOD !modulus of the wind parallel to the orography (m/s)
    REAL, INTENT(IN):: PUREF ! height of wind
    REAL, INTENT(IN):: PZ0   ! atmospheric roughness of wind
    REAL, INTENT(IN):: PTA    ! atmospheric temperature at level za (K)
    REAL, PARAMETER :: PPCA0 = -137.517
    REAL, PARAMETER :: PPCA1 = .184875
    REAL, PARAMETER :: PPCA2 = -3.00521E-4
    REAL, PARAMETER :: PPCA3 = 144.087
    REAL, PARAMETER :: PPCA4 = -78.1198
    REAL, PARAMETER :: PPCA5 = 20.5968
    REAL, PARAMETER :: PPCA6 = -2.5627
    REAL, PARAMETER :: PPCA7 = .113710
    REAL ::ZCOEF_SUBL_RATE,ZGAMMA,PVMOD_3M
    REAL,INTENT(OUT):: PSUBL_RATE

    PVMOD_3M=PVMOD*LOG(3./PZ0)/LOG(PUREF/PZ0) ! CHANGE IN WIND SPEED TO 3 M ACCORDING TO LOG WIND PROFILE
    ZGAMMA=LOG(PVMOD_3M)
    ZCOEF_SUBL_RATE=PPCA0+PPCA1*PTA+PPCA2*PTA**2+PPCA3*ZGAMMA+PPCA4*ZGAMMA**2&
                       +PPCA5*ZGAMMA**3+PPCA6*ZGAMMA**4+PPCA7*ZGAMMA**5

    PSUBL_RATE=10**(ZCOEF_SUBL_RATE)
    IF (PSUBL_RATE<0.) THEN
      PSUBL_RATE=0.
    END IF
  END SUBROUTINE

  SUBROUTINE BINTANJA_98_10M (PVMOD,PUREF,PZ0,PTA,PSUBL_RATE)

    !SNOW DRIFT SUBLIMATION RATE ACCORDING TO PARAMETRISATION OF BINTANJA_98 (EQUATION 10)
    ! USE OF COEFICIENTS FROM TABLE 3 AT 10M
    !Bintanja, R. (1998). The contribution of snowdrift sublimation to the surface mass balance of Antarctica. Annals of Glaciology, 27, 251-259.
    !DOI: https://doi.org/10.3189/1998AoG27-1-251-259
    !

    IMPLICIT NONE
    REAL,INTENT(IN) :: PVMOD !modulus of the wind parallel to the orography (m/s)
    REAL, INTENT(IN):: PUREF ! height of wind
    REAL, INTENT(IN):: PZ0   ! atmospheric roughness of wind
    REAL,INTENT(IN) :: PTA    ! atmospheric temperature at level za (K)
    REAL, PARAMETER :: PPCA0 = -50.5902
    REAL, PARAMETER :: PPCA1 = .183630
    REAL, PARAMETER :: PPCA2 = -2.96572E-4
    REAL, PARAMETER :: PPCA3 = 9.38304
    REAL, PARAMETER :: PPCA4 = -3.57458E-3
    REAL, PARAMETER :: PPCA5 = -.249308
    REAL, PARAMETER :: PPCA6 = -.127787
    REAL, PARAMETER :: PPCA7 = 2.99190E-2
    REAL ::ZCOEF_SUBL_RATE,ZGAMMA,PVMOD_10M
    REAL, INTENT(OUT):: PSUBL_RATE !SUBLIMATION RATE

    PVMOD_10M=PVMOD*LOG(10./PZ0)/LOG(PUREF/PZ0) ! CHANGE IN WIND SPEED TO 10 M ACCORDING TO LOG WIND PROFILE
    ZGAMMA=LOG(PVMOD_10M)
    ZCOEF_SUBL_RATE=PPCA0+PPCA1*PTA+PPCA2*PTA**2+PPCA3*ZGAMMA+PPCA4*ZGAMMA**2&
                       +PPCA5*ZGAMMA**3+PPCA6*ZGAMMA**4+PPCA7*ZGAMMA**5

    PSUBL_RATE=10**(ZCOEF_SUBL_RATE)
    IF (PSUBL_RATE<0.) THEN
      PSUBL_RATE=0.
    END IF
  END SUBROUTINE

  SUBROUTINE SBSM (PVMOD,PUREF,PZ0,PTA,PPS,PQA,KNLVLS_USE,PSNOWRHO,PSNOWHEAT, PSNOWDZ,PSUBL_RATE)

    !SNOW DRIFT SUBLIMATION RATE ACCORDING TO PARAMETRISATION OF SBSM ESSERY 1999
    ! https://documentation.help/CRHM_Borland/modules_sbsm.htm
    ! http://www.merrittnet.org/Papers/Essery_et_al_1999.pdf

    USE MODE_SNOW3L
    USE MODE_THERMOS
    USE MODD_SNOW_PAR, ONLY : XSNOWDMIN
    USE MODD_CSTS,     ONLY : XLMTT, XTT
    IMPLICIT NONE
    REAL,  INTENT(IN)    :: PVMOD  !modulus of the wind parallel to the orography (m/s)
    REAL,  INTENT(IN)    :: PUREF  !heigth of the wind (m)
    REAL, INTENT(IN)     :: PZ0   ! atmospheric roughness of wind
    REAL,  INTENT(IN)    :: PTA    ! atmospheric temperature at level za (K)
    REAL,  INTENT(IN)    :: PPS    ! atmospheric PRESURE at level za
    REAL,  INTENT(IN)    :: PQA    ! atmospheric SPECIFIC HUMIDITY at level za

    INTEGER, INTENT(IN) :: KNLVLS_USE
    REAL,  INTENT(IN) :: PSNOWHEAT, PSNOWDZ,PSNOWRHO
    REAL               :: PDIFF,PRSAT,PCOND
    REAL              :: ZSCAP, ZSNOWTEMP,PSCALE,PVP,PRMIX,PUNDERSAT,PPQSAT,ZVMOD
    REAL, PARAMETER :: PPM=18.01 !Molecular weight of water (kg/kmole)
    REAL, PARAMETER :: PPLS=2.838E6 !Latent heat of sublimation (J/kg)
    REAL, PARAMETER :: PPR=8313.0 !Universal gas constant (J/kmole/K)
    REAL,INTENT(OUT):: PSUBL_RATE !SUBLIMATION RATE

    ZVMOD=PVMOD*LOG(10./PZ0)/LOG(PUREF/PZ0)
    ZSCAP      = SNOW3LSCAP(PSNOWRHO)
    ZSNOWTEMP  = XTT + (PSNOWHEAT +XLMTT*PSNOWRHO*PSNOWDZ)/(ZSCAP*MAX(XSNOWDMIN/KNLVLS_USE,PSNOWDZ))
    ZSNOWTEMP  = MIN(XTT, ZSNOWTEMP)
    PPQSAT=PSAT_0D(PTA)
    PRMIX=PQA/(1-PQA)
    PVP=PRMIX*PPS/(0.622+PRMIX)
    PUNDERSAT=MAX((1-(PVP/PPQSAT)),0.) ! LIMITATION DE LA SUBLIMATION EN CAS DE SURSATUration
    PDIFF = 2.06e-5*(PTA/273.0)**1.75
    PRSAT = PPM*611.15*EXP(22.45*(PTA-273.15)/ZSNOWTEMP)/(PPR*(ZSNOWTEMP))
    PCOND = 0.00063*ZSNOWTEMP + 0.0673
    PSCALE= ((PPLS*PPM/(PPR*ZSNOWTEMP)) - 1.0)/(PCOND*PTA) + 1.0/(PPLS*PDIFF*PRSAT)
    PSUBL_RATE =(PUNDERSAT/PSCALE)*137.6*((ZVMOD/25.0)**5)/1000.0
    PSUBL_RATE=-PSUBL_RATE ! defined positively so we inverse
    IF (PSUBL_RATE<0.) THEN
      PSUBL_RATE=0.
    END IF


  END SUBROUTINE

  SUBROUTINE GORDON06 (PVFRIC_T,PVMOD,PUREF,PZ0,PTA,PPS,PQA,PRHOA,PQS)
    ! SNOW DRIFT SUBLIMATION RATE ACCORDING TO PARAMETRISATION OF GORDON
    ! https://doi.org/10.3137/ao.440303
    ! https://www.tandfonline.com/doi/pdf/10.3137/ao.440303?needAccess=true
    USE MODD_CSTS,     ONLY : XTT
    USE MODE_THERMOS, ONLY : QSATI

    REAL, INTENT(IN) :: PVFRIC_T ! wind friction thresold veolity
    REAL, INTENT(IN) :: PVMOD ! wind speed
    REAL, INTENT(IN) :: PUREF ! wind height
    REAL, INTENT(IN) :: PZ0 ! atmospheric roughness
    REAL, INTENT(IN) :: PTA ! air temperature
    REAL, INTENT(IN) :: PPS ! air presure
    REAL, INTENT(IN) :: PQA ! air specific humidity
    REAL, INTENT(IN) :: PRHOA ! air density
    REAL, INTENT(OUT):: PQS ! sublimation rate
    REAL             ::ZVFRIC_T5 ! wind speed thresold at 5m
    REAL             ::ZVMOD_5 ! wind speed at 5m
    REAL             ::ZQSATI
    REAL             ::ZRHI

    ZVFRIC_T5  = LOG(5./PZ0)*PVFRIC_T/0.4
    ZVMOD_5=PVMOD*LOG(5./PZ0)/LOG(PUREF/PZ0)
    ZQSATI = QSATI( PTA,PPS )
    ZRHI = PQA / ZQSATI
    ! computation of sublimation rate according to Gordon's PhD
    PQS = 0.0018 * (XTT/PTA)**4. * ZVFRIC_T5 * PRHOA * ZQSATI * (1.-ZRHI) * (ZVMOD_5/ZVFRIC_T5)**3.6
    IF (PQS<0) THEN
      PQS=0.
    END IF
  END SUBROUTINE

END MODULE SNOWPAPPUS_ENGINE
