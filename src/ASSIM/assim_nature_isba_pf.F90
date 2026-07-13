SUBROUTINE ASSIM_NATURE_ISBA_PF(IO, S, HPROGRAM, KI, HTEST, PLON_IN, PLAT_IN)

! -----------------------------------------------------------------------------
!
! Point-by-point Particle filter data assimilation of snow observations.
! fractions into Crocus Snowpack simulations 
!
! based on assim_nature_isba_enkf (Trygve Aspelien and Alina Barbu)
!
! The control vector is a subvector of PSPEC_ALB corresponding to the chosen observation bands
! - Choice in namelist and PSWE (for Fractional Snow Cover Area assimilation)
!
! The observations can be any element of (PB1,...,PB7, PSCF) - Choice in namelist
!

!
! First version : B. Cluzet, 2018
  
! -----------------------------------------------------------------------------
USE MODD_ISBA_OPTIONS_n, ONLY : ISBA_OPTIONS_t
USE MODD_ISBA_n,         ONLY : ISBA_S_t
!
USE MODD_SURFEX_MPI,     ONLY : NRANK, NPIO, NCOMM, NPROC
USE MODD_ASSIM,          ONLY : NVAR, NPRINTLEV, CVAR,    &
                               XF,COBS,NENS,              &
                               NOBS, NOBSTYPE, NNCO, XYO, &
                               NPFPART, NECHGU,CPF_CROCUS,&
                               LGLOBAL_PF, NEFF_PF,       &
                               XDLOC_PF
!
USE MODE_PF
!
USE MODD_SURFEX_MPI,     ONLY : NRANK, NPIO
USE MODD_SURF_PAR,       ONLY : XUNDEF
!
USE MODI_GATHER_AND_WRITE_MPI
!
#ifdef SFX_ARO
USE YOMMP,               ONLY : MYPROC
#endif
!
USE YOMHOOK,             ONLY : LHOOK,DR_HOOK
USE PARKIND1,            ONLY : JPRB
!
USE MODI_ABOR1_SFX
!
!
! -----------------------------------------------------------
!
IMPLICIT NONE
!
#ifdef SFX_MPI
INCLUDE "mpif.h"
#endif
!
TYPE(ISBA_OPTIONS_t), INTENT(INOUT) :: IO
TYPE(ISBA_S_t), INTENT(INOUT) :: S
!
CHARACTER(LEN=6),   INTENT(IN) :: HPROGRAM     ! program calling surf. schemes
INTEGER,            INTENT(IN) :: KI           ! = isize_full : nb of points
CHARACTER(LEN=2),   INTENT(IN) :: HTEST        ! must be equal to 'OK'
REAL(KIND=JPRB), DIMENSION (KI), INTENT(IN) ::  PLON_IN   ! coordinates of the points
REAL(KIND=JPRB), DIMENSION (KI), INTENT(IN) ::  PLAT_IN
!
!    Declarations of local variables
!

CHARACTER(LEN=7)   :: YMYPROC          
!
REAL, DIMENSION(:,:), ALLOCATABLE :: ZRINV_LOC             ! inverse of covariance matrix of obs. errors (localized cases)
REAL, DIMENSION(:,:), ALLOCATABLE :: ZRINV_G ! inverse of covariance matrix of obs. errors (global case)
REAL, DIMENSION(:,:), ALLOCATABLE :: ZINNOV       ! member innovation
!
REAL, DIMENSION(KI, KI, NOBSTYPE) :: ZBG_CORR   ! Spatial Correlation Matrix of the ensemble for each var (k-localization)
!
REAL, DIMENSION(IO%NPATCH)                :: ZLIKE        ! member likelihood
REAL, DIMENSION(IO%NPATCH, NENS)          :: ZLIKE_ENS    ! vector of members likelihood
REAL, DIMENSION(IO%NPATCH, NENS)          :: ZWEIGHT      ! members weight
REAL, DIMENSION(KI)                       :: ZALPHA_VEC   ! alpha vector for outputs
!
INTEGER :: IMYPROC
!
INTEGER :: JPT, JOBS, JLOC, JJ
INTEGER :: IDIM1
!
INTEGER, DIMENSION(KI, NOBSTYPE) :: IMASK_G ! globalization mask for each var at the given location.
INTEGER, DIMENSION(KI, NOBSTYPE) :: IMASK_LOC  ! localization mask for each var at the given location. (/!\in localized case,  will be affected until NLOC ONLY
INTEGER, DIMENSION(KI, KI, NOBSTYPE) :: IMASK_LOC_WRITE
INTEGER, DIMENSION(NOBSTYPE) :: ICOUNT_LOC ! counter for the number of points in IMASK at the given location(<=NLOC in loc case, <=NDIM_FULL in glob case)
!
!INTEGER, DIMENSION(KI,KI) :: ICOUNT_COMMON ! count the number of common members btw i and j
!INTEGER, DIMENSION(KI,KI, SIZE(XF, 3)) :: IMASK_COMMON ! mask of common defined members between i and j points.
!
INTEGER :: ISTAT
!
CHARACTER(LEN=14) :: YFMT
CHARACTER(LEN=14) :: YFMTG
CHARACTER(LEN=18) :: YFMTA
!
LOGICAL :: GDIAGONAL_R, GCOND_PM, GCOND_NK
!
INTEGER :: IKLOC, ILEN_COND, JK
LOGICAL, DIMENSION(KI) :: GMASK_COND ! LOGICAL Mask of the points where the Neff cond has not been reached yet.
REAL :: ZEFF_LOC, ZALPHA, ZWEIGHT_KEEP
!
REAL(KIND=JPRB)                            :: ZHOOK_HANDLE
!
!
IF (LHOOK) CALL DR_HOOK('ASSIM_NATURE_ISBA_PF',0,ZHOOK_HANDLE)
!
!############################# BEGINNING ###############################
!
IF (HTEST/='OK') THEN
  CALL ABOR1_SFX('ASSIM_NATURE_ISBA_PF: FATAL ERROR DURING ARGUMENT TRANSFER')
END IF
!
IF ( NPRINTLEV>0 .AND. NRANK==NPIO) THEN
  WRITE(*,*)
  WRITE(*,*) '   ---------------------------------'
  WRITE(*,*) '   |   ENTERING  PARTICLE-FILTER   |'
  WRITE(*,*) '   ---------------------------------'
  WRITE(*,*)
ENDIF
!
#ifdef SFX_ARO
IF ( MYPROC > 0 ) THEN 
  IMYPROC = MYPROC
ELSE
  IMYPROC = 1
ENDIF
#else
IMYPROC = NRANK+1
#endif
!
WRITE(YMYPROC(1:7),'(I7.7)') IMYPROC
!
IF ( NPRINTLEV > 0 .AND. NRANK==NPIO ) WRITE(*,*) 'number of patches =',IO%NPATCH
!
!############################# ANALYSIS ###############################
IF ( NPRINTLEV > 0 ) THEN
  IF (NRANK==NPIO) THEN
    WRITE(*,*) 'PERFORMING ANALYSIS'
  ENDIF
ENDIF
!!-----------------------LOCALIZED CASES-----------------------------------
IF (.NOT.LGLOBAL_PF) THEN 
  PRINT*,  '    --------------------------------'
  PRINT*,  '    |   PERFORMING LOCALIZED PF    |'
  PRINT*,  '    ================================'
  PRINT*,  '    | algo type   :', CPF_CROCUS
  IF (CPF_CROCUS == 'RLOCAL') THEN
    PRINT*,  '    | angle (deg) :', XDLOC_PF
  ENDIF
  PRINT*,  '    | target NEFF :', NEFF_PF
  PRINT*,  '    --------------------------------'
  !
  ! (KLOCAL): compute the ensemble bg correlation matrix
  IF (CPF_CROCUS == 'KLOCAL') THEN
    ZBG_CORR(:,:,:) = XUNDEF
    CALL BG_CORR_MATRIX(ZBG_CORR)  ! to define properly
  ENDIF
  ! BC 22/05/19 : iterative process on IKLOC for the localization, with Neff (Namelist) as a target
  ! BC 06/20 count the max number of valid obs per variable (for the KLOCAL)
  ! if no obs available (IKLOC==0), force IKLOC to 1
  ! to ensure that all the particles are properly replicated.
  CALL COUNT_MAX_OBS(IKLOC)
  !
  IF (IKLOC == 0) IKLOC=1
  !
  GMASK_COND(:) = .TRUE. ! if True, klocalize/inflate
  ILEN_COND = COUNT(GMASK_COND)  ! counter for the points still left to assimilate.
  !
  IMASK_LOC_WRITE(:,:,:)=0
  DO WHILE((ILEN_COND >= 1) .AND. (IKLOC >= 1))
    DO JPT=1,KI ! loop on grid points
      IF (GMASK_COND(JPT)) THEN
        ! init
        IMASK_LOC(:,:) = 0 ! localization mask for each var at the given location.
        ICOUNT_LOC(:) = 0  ! counter for the number of loc. obs for each var in IMASK at the given location.
        IDIM1=0            ! IDIM1 is equal to the total number of assimilated obs at the given location.
        !
        ! perform localization (select the obs).
        IF (CPF_CROCUS == 'KLOCAL') THEN
          CALL KLOCALIZE(IMASK_LOC, ICOUNT_LOC, IKLOC, JPT, ZBG_CORR, PLON_IN, PLAT_IN)
        ELSE
          CALL RLOCALIZE(IMASK_LOC, ICOUNT_LOC, JPT, PLON_IN, PLAT_IN)
        ENDIF
        !
        DO  JOBS=1, SIZE(ICOUNT_LOC)
          IDIM1 = IDIM1 + ICOUNT_LOC(JOBS)
        END DO
        !
        ALLOCATE(ZINNOV(IDIM1, NENS))
        !
        CALL INNOVATION(ZINNOV, IMASK_LOC, ICOUNT_LOC)
        !
        ALLOCATE(ZRINV_LOC(IDIM1, IDIM1))
        ZRINV_LOC(:,:) = 0
        !
        CALL SET_OBS_ERROR(ZRINV_LOC, IMASK_LOC, ICOUNT_LOC, GDIAGONAL_R)
        !
        CALL LIKELIHOOD(ZINNOV,ZRINV_LOC, ZLIKE_ENS)  ! BC 10/05/19 reference to JPT seems now useless.
        !
        ! perform machine-precision inflation in any case
        ! because where likelihoods are under the machine precision,
        ! LIKELIHOOD returns a default value. 
        ZALPHA = 1.
        CALL PM_INFLATE(ZALPHA, ZRINV_LOC, IMASK_LOC, ICOUNT_LOC, GDIAGONAL_R, ZINNOV, ZLIKE_ENS)
        CALL NORMALIZE(ZLIKE_ENS, ZWEIGHT)
        !
        CALL EFFWEIGHTS(ZWEIGHT, ZEFF_LOC)
        ! the inflation should be performed in the rlocal case or when the kloc degenerates with IKLOC = 1
        IF   ((  (CPF_CROCUS == 'RLOCAL') .OR. ((CPF_CROCUS == 'KLOCAL') .AND.(IKLOC==1))  ) .AND. (ZEFF_LOC < NEFF_PF)) THEN
          CALL INFLATE(ZALPHA, ZWEIGHT, ZRINV_LOC, IMASK_LOC, ICOUNT_LOC, GDIAGONAL_R, ZINNOV, ZLIKE_ENS, ZEFF_LOC)
        ENDIF
        ! overwrite the ZALPHA_VEC value (for output)
        ZALPHA_VEC(JPT) = ZALPHA
        DEALLOCATE(ZINNOV)
        DEALLOCATE(ZRINV_LOC)
        !
        ! actualize the mask and "check" IKLOC stopping condition :
        ! either we locally >= Neff OR we reached minimal IKLOC value (1)
        ! BC june 2020: be careful, the first term of ths condition must
        ! be TRUE if ZEFF_LOC is strictly equal to NEFF_PF (happens when both =1.)
        ! and the (absulute) epsilon must be superior
        ! to the epsilon in the convergence condition of the inflate ago.
        IF (((ZEFF_LOC - FLOAT(NEFF_PF)) > -1E-3) .OR. (CPF_CROCUS == 'KLOCAL' .AND. IKLOC == 1)) THEN
          GMASK_COND(JPT) = .FALSE.
          CALL RESAMPLING(ZWEIGHT, JPT)
          IMASK_LOC_WRITE(JPT, :,:) = IMASK_LOC(:,:)
        ENDIF
        !
      ENDIF
    ENDDO
    !
    ! stopping condition setting
    IF (CPF_CROCUS =='KLOCAL') IKLOC = IKLOC -1
    IF (CPF_CROCUS =='RLOCAL') IKLOC=0  ! kill the while in the RLOCAL case.
    !
    ILEN_COND = COUNT(GMASK_COND)
    !
  ENDDO
  !
  !*****write local selection********
  OPEN (unit=241,file='IMASK',status='unknown',IOSTAT=ISTAT)
  WRITE(YFMT, '(A1,I5,A8)')'(', KI,'(I5,A1))'
  DO JPT=1, KI
    DO JOBS=1, NOBSTYPE
      WRITE(241,YFMT) (IMASK_LOC_WRITE(JPT, JLOC, JOBS),',', JLOC=1, KI)
    END DO
  END DO
  CLOSE(241)
  !**********************************
  !*****write inflation factors******
  OPEN (unit=242,file='ALPHA',status='unknown',IOSTAT=ISTAT)
  WRITE(YFMTA, '(A1,I5,A12)') '(', KI,'(ES11.4,A1))'
  WRITE(242,YFMTA) (ZALPHA_VEC(JPT),',', JPT=1, KI)
  CLOSE(242)
  !**********************************
  !
!!-----------------------GLOBAL CASE------------------------------
ELSE
  ! for the sake of clarity, global case is outside the loop 
  ! while it could have been put in a IF JPT==1 case
  ! indeed, it's calling the same routines, except the Localize/globalize case, 
  ! which depends on point in the local case, while is done only once in the global case.
  PRINT*,  '    -------------------------------------'
  PRINT*,  '    |        PERFORMING GLOBAL PF       |'
  PRINT*,  '    ====================================='
  PRINT*,  '    | target NEFF :', NEFF_PF
  PRINT*,  '    -------------------------------------'

  IMASK_G(:,:) = 0
  ICOUNT_LOC(:)= 0
  CALL GLOBALIZE(IMASK_G, ICOUNT_LOC)
  !
  ! allocate ZRINV_G
  ! ZRINV_loc has a variable shape <= (NPTS*NOBSTYPE,NPTS*NOBSTYPE) depending on the definited/comparable points for each var
  IDIM1=0
  DO  JOBS=1, SIZE(ICOUNT_LOC)
    IDIM1 = IDIM1 + ICOUNT_LOC(JOBS)
  END DO
  !
  ALLOCATE(ZINNOV(IDIM1, NENS))
  !
  CALL INNOVATION(ZINNOV, IMASK_G, ICOUNT_LOC)
  !
  ALLOCATE(ZRINV_G(IDIM1, IDIM1))
  ZRINV_G(:,:) = 0 ! and not XUNDEF
  !
  ! compute weights using an inflation method
  ZALPHA=1 ! (ZALPHA <=1) R will be multiplied by 1/ZALPHA -> R^-1 mult by ZALPHA
  ! to ensure that the inflation method is proprly initiated,
  ! check the machine precision condition of the likelihoods computations
  !
  CALL PM_INFLATE(ZALPHA, ZRINV_G, IMASK_G, ICOUNT_LOC, GDIAGONAL_R, ZINNOV, ZLIKE_ENS)
  !
  ! check the efficient number of particles if >NEFF_PF, do not try to inflate.
  CALL NORMALIZE(ZLIKE_ENS, ZWEIGHT)
  CALL EFFWEIGHTS(ZWEIGHT, ZEFF_LOC)
  IF (ZEFF_LOC < NEFF_PF) THEN
    CALL INFLATE (ZALPHA, ZWEIGHT, ZRINV_G, IMASK_G, ICOUNT_LOC, GDIAGONAL_R, ZINNOV, ZLIKE_ENS, ZEFF_LOC)
  ENDIF
  !
  DEALLOCATE(ZINNOV)
  !
  CALL RESAMPLING(ZWEIGHT)
  !
  DEALLOCATE(ZRINV_G)
  !
  !*****write inflation factors******
  OPEN (unit=242,file='ALPHA',status='unknown',IOSTAT=ISTAT)
  WRITE(242,'(ES11.4)') ZALPHA
  CLOSE(242)
  !**********************************
  
ENDIF  

!2- write the weights
IF (LGLOBAL_PF) THEN
  ALLOCATE(S%NPART(1,SIZE(NPFPART, 3))) ! (1,1,nmbrs) in local

  S%NPART(:,:) = 0  ! ok for initialization ?
  S%NPART = NPFPART(:,1,:)
ELSE
  ALLOCATE(S%NPART(SIZE(NPFPART, 1), SIZE(NPFPART, 3))) ! (1,1,nmbrs) in global

  S%NPART(:,:) = 0  ! ok for initialization ?
  S%NPART = NPFPART(:,1,:)
ENDIF
OPEN (unit=118,file='PART',status='unknown',IOSTAT=ISTAT)
IF (LGLOBAL_PF) THEN
  WRITE(YFMTG, '(A1,I5,A8)') '(', SIZE(S%NPART,2),'(I5,A1))'
  WRITE(118,YFMTG) (S%NPART(1,JPT),',', JPT=1, SIZE(S%NPART,2))
ELSE
  WRITE(YFMTG, '(A1,I5,A8)') '(', SIZE(S%NPART,2),'(I5,A1))'
  DO JPT=1, SIZE(NPFPART, 1)
   WRITE(118,YFMTG) (S%NPART(JPT,JJ),',', JJ=1, SIZE(S%NPART,2))
  END DO
ENDIF
 CLOSE(118)
IF (LHOOK) CALL DR_HOOK('ASSIM_NATURE_ISBA_PF',1,ZHOOK_HANDLE)
!
END SUBROUTINE ASSIM_NATURE_ISBA_PF


