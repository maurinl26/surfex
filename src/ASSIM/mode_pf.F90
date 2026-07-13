MODULE MODE_PF

CONTAINS
!********************************************************************
!********************************************************************
!********************************************************************
SUBROUTINE COUNT_MAX_OBS(KKLOC)
!********************************************************************
! -Purpose
!    Compute the number of available observations for each variable in XYO
!    and return the maximal value.
!
! -Author
!     B. Cluzet 10/06/20
!********************************************************************
USE MODD_ASSIM, ONLY: XYO
USE MODD_SURF_PAR,   ONLY : XUNDEF
!
IMPLICIT NONE
!
INTEGER, INTENT(OUT) :: KKLOC     ! ids of Points 1 and 2 in the geometry
!
INTEGER, DIMENSION(SIZE(XYO,2))  :: KCOUNT
INTEGER :: JPT, JOBS
KKLOC = 0
DO JOBS=1, SIZE(XYO,2)
  KCOUNT(JOBS) = 0
  DO JPT=1, SIZE(XYO,1)
    IF (XYO(JPT, JOBS) .NE. XUNDEF) THEN
      KCOUNT(JOBS) = KCOUNT(JOBS) + 1
    ENDIF
  END DO
END DO
KKLOC = MAXVAL(KCOUNT)
!
END SUBROUTINE COUNT_MAX_OBS
!********************************************************************
!********************************************************************
!********************************************************************
SUBROUTINE BG_CORR_MATRIX(PBG_CORR)
!********************************************************************
! -Purpose
!    Compute Spatial Correlation matrix of the ensemble
!    Separately for each var
! -Author
!     B. Cluzet 09/05/19
!
!********************************************************************
USE MODD_SURFEX_MPI,    ONLY : NRANK, NPIO
USE MODD_ASSIM, ONLY : NENS, XF, NOBSTYPE, XDLOC_PF   ! XF is (NPTS, NPATCH,NENS, NOBSTYPE)
USE MODD_SURF_PAR,   ONLY : XUNDEF

IMPLICIT NONE

REAL, DIMENSION(:,:,:), INTENT(INOUT) :: PBG_CORR ! (NPTS, NPTS, NOBSTYPE) Spatial  (not cross-variables) Covariance matrix of the ensemble for each variable 

REAL, DIMENSION(SIZE(PBG_CORR,1)) :: ZBG_MEAN ! (NPTS)

REAL, DIMENSION(SIZE(PBG_CORR,1), SIZE(XF,3)) :: ZBG_ANOM ! (NPTS, NENS) normalized anomalies
REAL, DIMENSION(SIZE(PBG_CORR,1), SIZE(PBG_CORR,2)) :: ZBG_COV ! (NPTS, NPTS) : covariance matrix.
REAL, DIMENSION(SIZE(XF,3)) :: ZLOC_COV
REAL, DIMENSION(:), ALLOCATABLE :: ZDEVI_SQ, ZDEVJ_SQ
REAL :: ZSUM, ZSTDI, ZSTDJ

INTEGER, DIMENSION(SIZE(PBG_CORR,1),SIZE(PBG_CORR,2)) :: ICOUNT_COMMON
INTEGER, DIMENSION(SIZE(PBG_CORR,1),SIZE(PBG_CORR,2), SIZE(XF,3)) :: IMASK_COMMON ! mask of common defined members between i and j points.
                                                                                       ! BC 09/05/19: expensive but hard to avoid
INTEGER, DIMENSION(SIZE(PBG_CORR,1)) :: ICOUNT_VALID
INTEGER, DIMENSION(SIZE(PBG_CORR,1), SIZE(XF,3)) :: IMASK_VALID
INTEGER  :: JOBS, JPT, JPTI, JPTJ, JC, JENS, JENSV
!
INTEGER ::ISTAT
!
CHARACTER(LEN=18) :: YFMT
!********************************************************************
! loop on vars
PRINT*, '----------------------------------------'
PRINT*, '   computing correlation matrix         '

DO JOBS=1, NOBSTYPE
  PRINT*, '    var :', JOBS
  !0. ------------ count the defined members and the common members ------------------
  IF (JOBS == 1) THEN  ! all def/not def a the same time
    ICOUNT_VALID(:)=0
    IMASK_VALID(:,:)=0
    ICOUNT_COMMON(:,:)=0
    IMASK_COMMON(:,:,:)=0
    !
    DO JPTI=1,SIZE(PBG_CORR,1)
      ! count valid
      DO JENS=1, SIZE(XF,3)
        IF (XF(JPTI,1,JENS, JOBS) .NE. XUNDEF ) THEN ! f... the patches for now
          ICOUNT_VALID(JPTI) = ICOUNT_VALID(JPTI) + 1
          IMASK_VALID(JPTI, ICOUNT_VALID(JPTI)) = JENS
        ENDIF
      END DO
      !
      ! count common
      DO JPTJ=1,SIZE(PBG_CORR,2)
        DO JENS=1, SIZE(XF,3)
          IF ((XF(JPTI,1,JENS, JOBS) .NE. XUNDEF ) .AND. (XF(JPTJ,1,JENS, JOBS) .NE. XUNDEF )) THEN
            ICOUNT_COMMON(JPTI, JPTJ) = ICOUNT_COMMON(JPTI, JPTJ) + 1
            IMASK_COMMON(JPTI, JPTJ, ICOUNT_COMMON(JPTI, JPTJ)) = JENS
          ENDIF
        END DO
      END DO
    END DO
  ENDIF
  !
  !1. -------------- compute the mean-------------------------
  ZBG_MEAN(:) = XUNDEF
  DO JPT=1, SIZE(PBG_CORR,1)
    IF (ICOUNT_VALID(JPT) .NE. 0) THEN
      ZSUM = 0
      !
      DO JENSV=1, ICOUNT_VALID(JPT)
        ZSUM = ZSUM + XF(JPT,1,IMASK_VALID(JPT, JENSV), JOBS)
      END DO
      ZBG_MEAN(JPT) = ZSUM / ICOUNT_VALID(JPT)
    ENDIF
  END DO
  !
  !2. ------------- compute anoms X---------
  ZBG_ANOM(:,:) = XUNDEF
  DO JENS=1, SIZE(XF, 3)
    DO JPT=1, SIZE(PBG_CORR,1)
      IF (ICOUNT_VALID(JPT) .NE. 0) THEN
        ZBG_ANOM(JPT,JENS) = (XF(JPT,1,JENS,JOBS) - ZBG_MEAN(JPT))! / ICOUNT_VALID
      END IF
    END DO
  END DO
  !
  !3. ------------- compute correlation --------------
  !3.1 compute the spread of each sample
  !3.2 compute the covariance of the two samples
  ! 3.2.1 check the number of common points
  ! 3.2.2 if (check) compute the spread of each sample (redundancies...)
  ! 3.2.3 if (check) compute the covariance
  !3.3. normalize by the spreads
  DO JPTI=1, SIZE(PBG_CORR,1)
    !
    DO JPTJ=1, SIZE(PBG_CORR,1)
      ZSTDI = XUNDEF
      ZSTDJ = XUNDEF
      PBG_CORR(JPTI, JPTJ, JOBS) = XUNDEF
      ! 3.2.1
      IF (ICOUNT_COMMON(JPTI, JPTJ) .GT. NINT(0.1*SIZE(XF,3))) THEN  ! if count is lower than 10 percent of the pop, reject.
        ! 3.2.2
        ZLOC_COV(:) = 0.
        ALLOCATE(ZDEVI_SQ(ICOUNT_COMMON(JPTI, JPTJ)))
        ALLOCATE(ZDEVJ_SQ(ICOUNT_COMMON(JPTI, JPTJ)))
        DO JC=1, ICOUNT_COMMON(JPTI, JPTJ)
          ZDEVI_SQ(JC) = ZBG_ANOM(JPTI, IMASK_COMMON(JPTI, JPTJ, JC)) * &
                         ZBG_ANOM(JPTI, IMASK_COMMON(JPTI, JPTJ, JC))
          ZDEVJ_SQ(JC) = ZBG_ANOM(JPTJ, IMASK_COMMON(JPTI, JPTJ, JC)) * &
                         ZBG_ANOM(JPTJ, IMASK_COMMON(JPTI, JPTJ, JC))
          ! hand-made X.X^T computation -_-
          ZLOC_COV(JC) = (ZBG_ANOM(JPTI, IMASK_COMMON(JPTI, JPTJ, JC)) * ZBG_ANOM(JPTJ, IMASK_COMMON(JPTI, JPTJ, JC)))
        END DO
        ZSTDI = SQRT(SUM(ZDEVI_SQ))
        ZSTDJ = SQRT(SUM(ZDEVJ_SQ))
        DEALLOCATE(ZDEVI_SQ)
        DEALLOCATE(ZDEVJ_SQ)
        !
        ! corr = cov(1,2) / (std1*std2)
        IF ((ZSTDI < 1E-12).NEQV.(ZSTDJ < 1E-12)) THEN !  .NEQV. is exclusive or
          PBG_CORR(JPTI, JPTJ, JOBS) = XUNDEF
        ELSEIF ((ZSTDI < 1E-12) .AND. (ZSTDJ < 1E-12)) THEN
          PBG_CORR(JPTI, JPTJ, JOBS) = 1.
        ELSE
          PBG_CORR(JPTI, JPTJ, JOBS) = SUM(ZLOC_COV)/ (ZSTDI * ZSTDJ)
        ENDIF
      ENDIF
    END DO
  END DO
END DO
!! (dirty) writing of correlation matrix :
WRITE(YFMT, '(A1,I5,A12)') '(', SIZE(PBG_CORR,2),'(ES11.4,A1))'
OPEN (unit=117,file='BG_CORR',status='unknown',IOSTAT=ISTAT)
DO JPTI=1, SIZE(PBG_CORR, 1)
  DO JOBS=1, NOBSTYPE
    WRITE(117,YFMT) (PBG_CORR(JPTI, JPTJ,JOBS),',', JPTJ=1,SIZE(PBG_CORR,2))
  END DO
END DO
 CLOSE(117)
!! (dirty) writing of counts matrix :
!OPEN (unit=119,file='COUNT_COMMON',status='unknown',IOSTAT=ISTAT)
!DO JPTI=1, SIZE(PBG_CORR, 1)
!   WRITE(119,*) ICOUNT_COMMON(JPTI,:)
!END DO
! CLOSE(119) 

!OPEN (unit=120,file='ANOM',status='unknown',IOSTAT=ISTAT)
!DO JPTI=1, SIZE(ZBG_ANOM, 1)
!   WRITE(120,*) ZBG_ANOM(JPTI,:)
!END DO
! CLOSE(120) 

!
PRINT*, '--END COMPUTING CORRELATION--'
!
END SUBROUTINE BG_CORR_MATRIX

!********************************************************************
!********************************************************************
!********************************************************************
SUBROUTINE KLOCALIZE(KMASK, KCOUNT_LOC, KKLOC, KPT, PBG_CORR, PLON_IN, PLAT_IN)
!********************************************************************
! -Purpose
!    Compute Mask of k-localisation on a given point.
!    Taking account of the availability of obs and ensemble.
! -Author
!     B. Cluzet 09/05/19
!
!********************************************************************
USE MODD_SURFEX_MPI,    ONLY : NRANK, NPIO
USE MODD_ASSIM, ONLY : NENS, XYO, XF, CPF_CROCUS, LLOO_PF
USE MODD_SURF_PAR,   ONLY : XUNDEF
USE PARKIND1,            ONLY : JPRB
!
IMPLICIT NONE
INTEGER, DIMENSION(:,:), INTENT(INOUT) :: KMASK   ! Mask of the local selection of points to assimilate
INTEGER, DIMENSION(:), INTENT(INOUT) :: KCOUNT_LOC   ! counter for KMASK

INTEGER, INTENT(IN) :: KKLOC         !  Current localization size
INTEGER, INTENT(IN) :: KPT         !  Index of the current local point "around" which the selection is performed
REAL, DIMENSION(:,:,:), INTENT(IN) :: PBG_CORR ! Spatial  (not cross-variables) Covariance matrix of the ensemble for each variable (k-localization)
REAL(KIND=JPRB), DIMENSION (:), INTENT(IN) ::  PLON_IN  ! coordinate vectors of the geometry
REAL(KIND=JPRB), DIMENSION (:), INTENT(IN) ::  PLAT_IN
REAL, DIMENSION(:), ALLOCATABLE :: ZCOMP_CORR ! selection of the corrs of comparable points (absolute value).
LOGICAL, DIMENSION(:), ALLOCATABLE :: GPOP_SEL

INTEGER, DIMENSION(SIZE(PBG_CORR, 1)):: IMASK_COMP  ! mask where ensemble and obs are comparable 
                                                    ! e.g. cross-cov computable AND obs exists.
INTEGER ::ICOUNT_COMP                               ! len of this mask. 

INTEGER :: JOBS, JPTJ, JCOMP, IMAX

REAL :: ZDEG  ! angle between pt1 and pt2 (degrees)

LOGICAL ::GCOND_RLOC
!
!********************************************************************
DO JOBS=1, SIZE(PBG_CORR,3)
  !
  ! find the points where you can compare ensemble and obs.
  ! BC update 25/11/19 : discard correlations inferior to 0.3
  ICOUNT_COMP = 0
  IMASK_COMP(:) = 0
  DO JPTJ=1, SIZE(PBG_CORR,2) ! work by lines 
    !
    ! 0st cond : in the localisation radius
    ! 1st cond eq >10% of members in common btw KPT and JPTJ
    ! 2nd that obs exists
    ! 3rd that correlation is significant
    CALL IN_RADIUS(GCOND_RLOC, KPT, JPTJ, PLON_IN, PLAT_IN, ZDEG)
    IF ((GCOND_RLOC) .AND. (PBG_CORR(KPT,JPTJ,JOBS) .NE. XUNDEF) .AND. &
    (XYO(JPTJ, JOBS) .NE. XUNDEF) .AND. (ABS(PBG_CORR(KPT,JPTJ,JOBS)) > 0.3)) THEN
      ICOUNT_COMP = ICOUNT_COMP + 1
      IMASK_COMP(ICOUNT_COMP) = JPTJ
    ENDIF
  END DO
  !
  ! prepare the selection
  ALLOCATE(ZCOMP_CORR(ICOUNT_COMP))
  ALLOCATE(GPOP_SEL(ICOUNT_COMP))
  DO JCOMP =1, ICOUNT_COMP
    ZCOMP_CORR(JCOMP) = ABS(PBG_CORR(KPT, IMASK_COMP(JCOMP), JOBS))  ! absolute value.
    GPOP_SEL(JCOMP) = .TRUE.
  END DO 
  !
  ! in that set of points, select up to KKLOC biggest correlations.
  KCOUNT_LOC(JOBS) = 0
  DO WHILE((KCOUNT_LOC(JOBS) < KKLOC) .AND. (KCOUNT_LOC(JOBS) < ICOUNT_COMP))
    IMAX = MAXLOC(ZCOMP_CORR, DIM=1, MASK=GPOP_SEL) ! Second arg is a mask where to search.
    KCOUNT_LOC(JOBS) = KCOUNT_LOC(JOBS) + 1
    KMASK(KCOUNT_LOC(JOBS), JOBS) = IMASK_COMP(IMAX)
    GPOP_SEL(IMAX) = .FALSE. ! set false to the mask where the value has been selected
  END DO
  !
  DEALLOCATE(ZCOMP_CORR)
  DEALLOCATE(GPOP_SEL)
END DO
END SUBROUTINE KLOCALIZE

!********************************************************************
!********************************************************************
!********************************************************************
SUBROUTINE RLOCALIZE(KMASK, KCOUNT_LOC, KPT, PLON_IN, PLAT_IN)
!********************************************************************
! -Purpose
!    Compute Mask of r-localisation on a given point.
!    Taking account of the availability of obs and ensemble.
!
! -Author
!     B. Cluzet 09/05/19
!
!********************************************************************
USE MODD_SURFEX_MPI,    ONLY : NRANK, NPIO
USE MODD_ASSIM, ONLY : NENS, XYO, XF, CPF_CROCUS, XDLOC_PF, LLOO_PF
USE MODD_SURF_PAR,   ONLY : XUNDEF
USE PARKIND1,            ONLY : JPRB
USE MODI_ABOR1_SFX
!
IMPLICIT NONE
!
INTEGER, DIMENSION(:,:), INTENT(INOUT) :: KMASK   ! (NOBS x NOBSTYPE) Mask of the local selection of points to assimilate
INTEGER, DIMENSION(:), INTENT(INOUT) :: KCOUNT_LOC   ! (NOBSTYPE) counter for KMASK
INTEGER, INTENT(IN) :: KPT         !  Index of the current local point "around" which the selection is performed
REAL(KIND=JPRB), DIMENSION (:), INTENT(IN) ::  PLON_IN  ! coordinate vectors of the geometry
REAL(KIND=JPRB), DIMENSION (:), INTENT(IN) ::  PLAT_IN

REAL :: ZDEG  ! angle between pt1 and pt2 (degrees)
INTEGER :: JOBS, JPTOBS, JENS

LOGICAL :: G1VALID ! True if at lest 1 member is defined at this point
LOGICAL :: GCOND_RLOC ! True if the considered point is within a radius of XDLOC_PF
!
!********************************************************************
! purely localised case (Cluzet et al., (submitted), Deschamps-Berger et al., ...)
! characterised by a default value of 1E-8 on XDLOC_PF.
IF (XDLOC_PF < 1E-7) THEN ! 1e-7 degrees approx. 1cm on the earth's surface
  IF (LLOO_PF) THEN
    CALL ABOR1_SFX('Leave-one-out with purely localised does not make much sense...')
  ELSE
    DO JOBS=1, SIZE(KMASK, 2)
      CALL TEST_DEF_ENS_OBS(KCOUNT_LOC, KMASK, KPT, JOBS)
    END DO
  ENDIF
ELSE  !localised within a radius of XDLOC_PF
  DO JOBS=1, SIZE(KMASK, 2)  ! loop on obs. types
    KCOUNT_LOC(JOBS)=0
    DO JPTOBS=1, SIZE(KMASK, 1)  ! loop on all obs locations to check their distance with the considered point
      CALL IN_RADIUS(GCOND_RLOC, KPT, JPTOBS, PLON_IN, PLAT_IN, ZDEG)
      IF (GCOND_RLOC) THEN
        CALL TEST_DEF_ENS_OBS(KCOUNT_LOC, KMASK, JPTOBS, JOBS)
      ENDIF
    END DO
  END DO
ENDIF
!
END SUBROUTINE RLOCALIZE

!********************************************************************
!********************************************************************
!********************************************************************
SUBROUTINE GLOBALIZE(KMASK, KCOUNT_LOC)
!********************************************************************
! -Purpose
!    Compute Global mask (excluding not snowy/scf criteria)
!    Taking account of the availability of obs and ensemble.
! -Author
!     B. Cluzet 10/05/19
!
!********************************************************************
USE MODD_ASSIM, ONLY : XYO, XF
USE MODD_SURF_PAR,   ONLY : XUNDEF

IMPLICIT NONE
INTEGER, DIMENSION(:,:), INTENT(INOUT) :: KMASK   ! Mask of the local selection of points to assimilate
INTEGER, DIMENSION(:), INTENT(INOUT) :: KCOUNT_LOC   ! counter for KMASK
INTEGER :: JPT, JOBS, JENS

LOGICAL :: G1VALID ! for a given point in the loop, True if at least 1 member is defined.
!
!********************************************************************
!
DO JOBS=1,SIZE(KMASK, 2)
  DO JPT=1,SIZE(KMASK,1)
    CALL TEST_DEF_ENS_OBS(KCOUNT_LOC, KMASK, JPT, JOBS)
  END DO
END DO
!
END SUBROUTINE GLOBALIZE

!********************************************************************
!********************************************************************
!********************************************************************
SUBROUTINE TEST_DEF_ENS_OBS(KCOUNT_LOC, KMASK, KPT, KOBS)
!********************************************************************
! -Purpose
!    Small useful piece of code. At location KPT,
!    check the definition of at least 1 member, 
!    and the existence of an observation of type KOBS
!    Update KCOUNT_LOC and KMASK accordingly
!
! -Author
!     B. Cluzet 10/06/20
!
!********************************************************************

USE MODD_ASSIM,    ONLY: XF, XYO
USE MODD_SURF_PAR, ONLY: XUNDEF
!
IMPLICIT NONE
!
INTEGER, DIMENSION(:,:), INTENT(INOUT) :: KMASK   ! (NOBS x NOBSTYPE) Mask of the local selection of points to assimilate
INTEGER, DIMENSION(:), INTENT(INOUT) :: KCOUNT_LOC   ! (NOBSTYPE) counter for KMASK
INTEGER, INTENT(IN)  :: KPT
INTEGER, INTENT(IN)  :: KOBS
!
INTEGER :: JENS
LOGICAL :: G1VALID
!
! test if at least 1 defined member
G1VALID= .FALSE.
JENS = 0
DO WHILE( (JENS < SIZE(XF,3)) .AND. .NOT. G1VALID)
  JENS=JENS+1
  IF(XF(KPT,1,JENS, KOBS) .NE. XUNDEF) G1VALID = .TRUE.
  !
END DO
!
!account for obs. availability and at least 1 mb.:
IF (G1VALID) THEN
  IF (XYO(KPT, KOBS) .NE. XUNDEF) THEN
    KCOUNT_LOC(KOBS)= KCOUNT_LOC(KOBS) + 1
    KMASK(KCOUNT_LOC(KOBS),KOBS) = KPT
  ENDIF
ENDIF
!
END SUBROUTINE TEST_DEF_ENS_OBS

!********************************************************************
!********************************************************************
!********************************************************************
SUBROUTINE INNOVATION(PINNOV, KMASK_LOC, KCOUNT_LOC)
!********************************************************************
! -Purpose
!    Compute Innovation vector (loc and glob)
!    @TODO : IMPLEMENT PROPERLY for glob.
! -Author
!     B. Cluzet 10/05/19
!
!********************************************************************
USE MODD_ASSIM, ONLY : XF, XYO, NOBSTYPE, NENS
USE MODD_SURF_PAR,   ONLY : XUNDEF
IMPLICIT NONE

INTEGER, DIMENSION(:,:), INTENT(IN) :: KMASK_LOC  ! localization mask
INTEGER, DIMENSION(:), INTENT(IN) :: KCOUNT_LOC   ! counter for KMASK
REAL, DIMENSION(:,:), INTENT(INOUT) :: PINNOV

INTEGER :: JLOC, JOBS, JENS, JPT, JVAR
!
!********************************************************************
JVAR = 0
DO JOBS=1, NOBSTYPE
  IF (KCOUNT_LOC(JOBS)  .NE. 0) THEN  ! else loop from 1 to 0 -_-
    DO JPT=1, KCOUNT_LOC(JOBS)
      JVAR = JVAR + 1
      JLOC = KMASK_LOC(JPT, JOBS)
      DO JENS =1,SIZE(PINNOV,2)
        ! be careful, if a member is not defined, then innov will be XUNDEF.
        !
        ! => enforce 0.2 (soil albedo) when the member is not defined (makes some sense for the reflectances.)
        !@ TODO : BE CAREFUL!! there should be a test if we are assimilating refl., this doesn't make any sense otherwise
        IF (XF(JLOC,1, JENS, JOBS) .NE. XUNDEF) THEN
          PINNOV(JVAR, JENS) = XYO(JLOC, JOBS) - XF(JLOC,1, JENS, JOBS) ! !f... the patches.
        ELSE
          PINNOV(JVAR, JENS) = XYO(JLOC, JOBS) - 0.2
        ENDIF
      END DO
    END DO
  ENDIF
END DO
!
END SUBROUTINE INNOVATION
!********************************************************************
!********************************************************************
!********************************************************************
SUBROUTINE SET_OBS_ERROR(PRINV, KMASK, KCOUNT_LOC, LDIAGONAL_R, PALPHA)
!********************************************************************
! -Purpose
!    Compute Observation error matrix in global and local cases
! -Author
!     B. Cluzet 30/04/19
!
!********************************************************************
USE MODD_SURFEX_MPI,    ONLY : NRANK, NPIO
USE MODD_ASSIM, ONLY : NENS, XYO, XYO_SCF, NOBSTYPE, LGLOBAL_PF, COBS, XERROBS, &
                       XERROBS_FACTOR, NRETM, NRETS, NRETD,NRETR,NCOUNTM
IMPLICIT NONE

REAL, DIMENSION(:,:), INTENT(INOUT)    :: PRINV ! Observation error matrix and its inverse.
INTEGER, DIMENSION(:,:), INTENT(IN) :: KMASK
INTEGER, DIMENSION(:), INTENT(IN) :: KCOUNT_LOC   ! counter for KMASK
LOGICAL, INTENT(OUT) :: LDIAGONAL_R
REAL, INTENT(INOUT), OPTIONAL :: PALPHA  ! optional : inflate the R matrix by 1/ZALPHA (following Larue et al., HESS, 2018)
REAL, DIMENSION(SIZE(PRINV,1), SIZE(PRINV,2)) :: ZR

REAL :: ZALPHA
INTEGER :: JPT, JLOC, JJ, JOBS, JVAR
LOGICAL :: GASSIMALL
!
!********************************************************************
!in this routine, the assumption is made that R is diagonal.
LDIAGONAL_R = .TRUE.
! default value for ZALPHA is 1
IF (.NOT.(PRESENT(PALPHA))) THEN
  ZALPHA = 1.
ELSE
  ZALPHA = PALPHA
ENDIF
!
! ---------- LOCAL AND GLOBAL CASES ASAME ! -----------------
JVAR=0
DO JOBS = 1,NOBSTYPE
  !
  IF (KCOUNT_LOC(JOBS)  .NE. 0) THEN  ! else loop from 1 to 0 -_-
    DO JPT = 1, KCOUNT_LOC(JOBS)
      !!-----------------*****SET OBSERVATION ERROR ********------------------
      !
      ! B. Cluzet on each point setting localized observation error matrix
      ! according to values in XERROBS and multiplying by XERROBS_FACTOR (namelist)
      ! be careful SCF error is relative so far.
      !
      ! But first if SCF assimilation is asked, we need to check SCF obs value at point JPT:
      !  - if SCF > 0.8 : assimilate reflectance and SCF IF ONLY there is snow in the model
      !  - else assimilate SCF only.

      ! Two different cases : 1) assimilation of MODIS and/or SCF, 2) assimilation of snow heights
      !-----------------------------------------------------------------------

      ! get current localized pt
      JLOC = KMASK(JPT, JOBS)
      ! set current index in PRINV
      JVAR = JVAR + 1
      ! init
      !
      !1-- assimilation of MODIS and/or SCF
      IF ((NRETM==1) .OR. (NRETS==1)) THEN
        IF (NRETS==1) THEN
          IF (XYO(JLOC,NCOUNTM+1)>0.8) THEN
            GASSIMALL = .TRUE.  ! assimilate reflectance and SCF......
          ELSE
            GASSIMALL = .FALSE. ! assimilate  scf only.
          ENDIF
        ELSE
          GASSIMALL = .TRUE.    !...... or reflectance only, if IRETS==0
        ENDIF
        !
        !1a-- assimilation of MODIS and SCF
        IF (GASSIMALL) THEN ! assimilate reflectance and SCF
          !
          IF (COBS(JOBS) == "SCF" ) THEN
            ! relative error for fractional snow cover measurement only
            ! this one has to be updated at each point
            ! need to prevent errors =0. for no snow in obs :cut at the value for scf=0.1
            IF (XYO(JLOC, JOBS)>0.1) THEN
              ZR(JVAR, JVAR) = XERROBS(JOBS)*XYO(JLOC,JOBS)*XERROBS_FACTOR(JOBS) ! BC bug here
            ELSE
              ZR(JVAR, JVAR) = XERROBS(JOBS)*0.1*XERROBS_FACTOR(JOBS)
            ENDIF

          ELSE ! usual cases (modis bands, reflectance)
            ! where there is no snow in the obs, increase the errors by a factor of 10.
            IF (XYO_SCF(JLOC) ==0) THEN
              ZR(JVAR, JVAR) = XERROBS(JOBS)*XERROBS_FACTOR(JOBS) * 10.
            ELSE
              ZR(JVAR, JVAR) = XERROBS(JOBS)*XERROBS_FACTOR(JOBS)
            ENDIF
          ENDIF
        ENDIF
        !
      !2-- assimilation of snow depth
      ELSEIF (NRETD==1) THEN
        ZR(JVAR,JVAR) = XERROBS(JOBS)*XERROBS_FACTOR(JOBS)
      ENDIF
      PRINV(JVAR, JVAR) = 1./((1./ZALPHA) * ZR(JVAR, JVAR))
    END DO
  ENDIF
END DO
!
END SUBROUTINE SET_OBS_ERROR
!********************************************************************
!********************************************************************
!********************************************************************
SUBROUTINE LIKELIHOOD(PINNOV,PRINV, PLIKE_ENS)
!********************************************************************
! -Purpose
!    Compute likelihoods
! -Author
!     B. Cluzet 30/04/19
!
!********************************************************************

USE MODD_ASSIM, ONLY : NENS, XYO, XF, LGLOBAL_PF
IMPLICIT NONE

REAL, DIMENSION(:,:), INTENT(IN)    :: PINNOV, PRINV ! Observation error matrix and its inverse.
REAL, DIMENSION(:,:), INTENT(OUT)     :: PLIKE_ENS
REAL, DIMENSION(SIZE(PINNOV,1)) :: ZLOC1
REAL, DIMENSION(SIZE(PINNOV,1), SIZE(PINNOV,2)) :: ZLIKE_G, ZQ
REAL :: ZLOC2
!REAL, DIMENSION(1,1,1): ZLOC3
REAL :: ZLOC3
INTEGER :: JENS, JVAR1, JVAR2, JVAR
!
!********************************************************************
!**** Likelihoods (weights) computation**
!      WITH NO ASSUMPTION ON R
! /!\ be careful : there is a risk of manipulating really tiny values when:
!  |  IKLOC AND NVAR are high (high pb dimension) AND prescribed errors are too low.
!  |  as a consequence, the member will be ruled out from the sample (likelihood=0)
!  |
!********************************************************************
DO JENS = 1,NENS
  !DO JJ=1, IM%O%NPATCH ! f... the patches
  !1. ZLOC1 = MATMUL(PINNOV, PRINV) -> TODO : check efficiency
  DO JVAR1=1, SIZE(PRINV,1)
    ZLOC1(JVAR1)=0
    DO JVAR2=1, SIZE(PRINV,2)
      ZLOC1(JVAR1) = ZLOC1(JVAR1) + PRINV(JVAR1, JVAR2) * PINNOV(JVAR2, JENS)
    END DO
  END DO
  !
  !2. ZLOC2 = MATMUL(ZLOC1, PINNOV.T)
  ZLOC2 = 0
  DO JVAR1 = 1, SIZE(PINNOV,1)
    ZLOC2 = ZLOC2 + ZLOC1(JVAR1)* PINNOV(JVAR1, JENS)
  END DO
  !3. LIKELIHOOD = e-1/2*PINNOV.PRINV.PINNOV.T
  ! /!\ sometimes reaching the limits of double accuracy.
  IF (ZLOC2 > 1300.) THEN !e(-745) = 5e-324; e(-746) = 0. (machine precision limit, take some margins for the weighting)
   ! PRINT*, 'WARNINGLIKELIHOOD : dealing with too low likelihoods : cutting.'
   ! PRINT*, 'Check _FillValue. Increase errors/reduce the problem dimension.'
   ! PRINT*, 'ZLOC2', ZLOC2
   ! PRINT*, 'JENS', JENS
   ! PRINT*,'INNOV', PINNOV(:,JENS)
    ZLOC2 = 1300.
    !CALL ABOR1_SFX('LIKELIHOOD : dealing with too low likelihoods. Check _FillValue. Increase errors/reduce the problem dimension.')
  ENDIF
  PLIKE_ENS(1,JENS)=EXP(-0.5*ZLOC2)

  !END DO ! f... the patches
ENDDO

END SUBROUTINE LIKELIHOOD
!********************************************************************
!********************************************************************
!********************************************************************
SUBROUTINE NORMALIZE(PLIKE_ENS, PWEIGHT)
!********************************************************************
! -Purpose
!    Compute PF weights (local or global) by normalization of lilkelihoods
! -Author
!     B. Cluzet 24/04/18
! -Modifications
!     B. Cluzet, 2019 : dev mutualization between local/global
!********************************************************************

USE MODD_ASSIM, ONLY : NENS, LGLOBAL_PF
IMPLICIT NONE


REAL, DIMENSION(:,:), INTENT(IN)    :: PLIKE_ENS    ! vector of members likelihood
REAL, DIMENSION(:,:), INTENT(INOUT) :: PWEIGHT      ! members weight
!
REAL, DIMENSION(SIZE(PLIKE_ENS, 1)) :: ZLIKESUM
INTEGER :: JP, JENS
!
!*********************************************************************
!
ZLIKESUM = SUM(PLIKE_ENS,2) ! sum of all likelihoods for normalization of weights
!


!**** Weighting (likelihood normalisation)**************
DO JP=1,SIZE(PLIKE_ENS, 1)
  IF(ZLIKESUM(JP) == 0.) THEN
    PWEIGHT(JP,:)=1./NENS
  ELSE
    DO JENS = 1,NENS
      ! here, the denominator is protected from 0 by the abort in likelihood
      PWEIGHT(JP,JENS) = PLIKE_ENS(JP,JENS)/ZLIKESUM(JP)
    END DO
  ENDIF
END DO


END SUBROUTINE NORMALIZE
!********************************************************************
!********************************************************************
!********************************************************************
SUBROUTINE INFLATE_FANNY(PWEIGHT,PALPHA, PLIKE_ENS, PLIKE_KEEP)
!********************************************************************
! -Status : to check (not finalized)
! -Purpose
!    Compute the R-inflation factor PALPHA 
!    to select NEFF_PF particles with certainty.
!    Based on Larue et al., HESS, 2018 (p.5730)
! -Author
!     B. Cluzet 20/09/19
!********************************************************************

USE MODD_ASSIM, ONLY : NEFF_PF, NENS
IMPLICIT NONE
REAL, DIMENSION(:,:), INTENT(INOUT) ::PWEIGHT
REAL, INTENT(INOUT)              :: PALPHA       ! inflation factor
REAL, DIMENSION(:,:), INTENT(IN) :: PLIKE_ENS    ! vector of members likelihood
REAL, INTENT(IN)                 :: PLIKE_KEEP   ! Nkeepth biggest likelihood  
!
REAL, DIMENSION(SIZE(PLIKE_ENS,2)) :: ZLIKE_EXP  ! like to the exponent alpha
REAL, DIMENSION(1,SIZE(PLIKE_ENS,2)) :: ZLIKE_EXP_ARG  ! like to the exponent alpha
REAL :: ZTOL, ZALPHA_0
!
LOGICAL :: GCOND
!
INTEGER :: JENS
!
!
!*********************************************************************
ZTOL = 0.002
ZALPHA_0=PALPHA + 2*ZTOL

DO WHILE (ABS(ZALPHA_0 - PALPHA) .GT. ZTOL)
  ZALPHA_0 = PALPHA
  DO JENS = 1,SIZE(PLIKE_ENS,2)
    ZLIKE_EXP(JENS) = EXP(ZALPHA_0 * LOG(PLIKE_ENS(1, JENS)))  ! au diable les patches ! BC FPE ici
  END DO
  PALPHA = (LOG(1. / FLOAT(NENS)) + LOG(SUM(ZLIKE_EXP)))/ LOG(PLIKE_KEEP)  ! (A.5)
  PRINT*, PALPHA
END DO  
!
DO JENS = 1,SIZE(PLIKE_ENS,2)
    ZLIKE_EXP_ARG(1, JENS) = EXP(PALPHA * LOG(PLIKE_ENS(1, JENS)))  ! au diable les patches ! BC FPE ici
END DO
!
CALL NORMALIZE(ZLIKE_EXP_ARG, PWEIGHT)

END SUBROUTINE INFLATE_FANNY
!********************************************************************
!********************************************************************
!********************************************************************
SUBROUTINE PM_INFLATE(PALPHA, PRINV, KMASK, KCOUNT_LOC, LDIAGONAL_R, PINNOV, PLIKE_ENS)
!********************************************************************
! 
! -Purpose
! ensure that the computation of the likelihoods do not go beyond machine-precision.
!     B. Cluzet 24/09/19
!********************************************************************
USE MODD_ASSIM, ONLY : NEFF_PF 
!
USE MODE_QUICK_SORT, ONLY:  QUICK_SORT

IMPLICIT NONE
REAL, INTENT(INOUT)              :: PALPHA       ! inflation factor
REAL, DIMENSION(:,:), INTENT(INOUT)    :: PRINV ! Observation error matrix and its inverse.
INTEGER, DIMENSION(:,:), INTENT(IN) :: KMASK
INTEGER, DIMENSION(:), INTENT(INOUT) :: KCOUNT_LOC   ! counter for KMASK
LOGICAL, INTENT(INOUT) :: LDIAGONAL_R
REAL, DIMENSION(:,:), INTENT(IN)    :: PINNOV
REAL, DIMENSION(:,:), INTENT(INOUT) :: PLIKE_ENS    ! vector of members likelihood
!
REAL, DIMENSION(SIZE(PLIKE_ENS, 2))    :: ZSORTED    ! sorted likelihoods or weights, au diable les patches
REAL :: ZLIKE_KEEP, ZALPHA
!
INTEGER, DIMENSION(SIZE(PLIKE_ENS, 2)) :: IORDER
INTEGER :: JITER
LOGICAL :: GCOND_PM
!
GCOND_PM = .TRUE.
ZALPHA = PALPHA
JITER = 0
DO WHILE(GCOND_PM)
  JITER = JITER + 1
  CALL SET_OBS_ERROR(PRINV, KMASK, KCOUNT_LOC, LDIAGONAL_R, PALPHA = ZALPHA)
  !
  CALL LIKELIHOOD(PINNOV,PRINV, PLIKE_ENS)
  !
  ZSORTED = PLIKE_ENS(1,:)
  CALL QUICK_SORT(ZSORTED, IORDER)

  ! at least Nkeepth biggest LIKE should exceed machine precision
  ZLIKE_KEEP = ZSORTED(SIZE(ZSORTED) - NEFF_PF + 1)
  IF ((ZLIKE_KEEP > EXP(-649.)) .OR. (JITER > 1000)) THEN
    GCOND_PM = .FALSE.
  ELSE
    ZALPHA = 0.1 * ZALPHA ! inflate R because too many weights below the machine precision.
  ENDIF
END DO
IF(JITER > 1000) THEN
  CALL ABOR1_SFX('INFLATE : failed to reach machine precision in likelihoods. Increase obs errors.')
ENDIF
!
PALPHA = ZALPHA
!
END SUBROUTINE PM_INFLATE
!********************************************************************
!********************************************************************
!********************************************************************
SUBROUTINE INFLATE(PALPHA, PWEIGHT, PRINV, KMASK, KCOUNT_LOC, LDIAGONAL_R, PINNOV, PLIKE_ENS, PEFF_LOC)
!********************************************************************
! 
! -Purpose
    ! optimization of the inflation factor of the matrix
    ! we want to find (1/ALPHA) inflation factor on R such that Neff = Neff_PF.
    ! Neff is define as 1/sum(ZWEIGHTS^2) (Doucet et al., Sequential Monte-Carlo Methods, Springer, 2001)
    !
    ! for that, we apply the secant-bisect hybrid method to F : x : |-> Neff(x) - Neff_PF in ]0, 1[
    ! inspired on rtsafe algo from numerical recipes in FORTRAN 77, 
! -Author
!     B. Cluzet 24/09/19
!********************************************************************
USE MODD_ASSIM, ONLY : NEFF_PF
IMPLICIT NONE
REAL, DIMENSION(:,:), INTENT(INOUT) ::PWEIGHT
REAL, INTENT(INOUT)              :: PALPHA       ! inflation factor
REAL, DIMENSION(:,:), INTENT(INOUT)    :: PRINV ! Observation error matrix and its inverse.
INTEGER, DIMENSION(:,:), INTENT(IN) :: KMASK
INTEGER, DIMENSION(:), INTENT(INOUT) :: KCOUNT_LOC   ! counter for KMASK
LOGICAL, INTENT(INOUT) :: LDIAGONAL_R
REAL, DIMENSION(:,:), INTENT(IN)    :: PINNOV
REAL, DIMENSION(:,:), INTENT(INOUT) :: PLIKE_ENS    ! vector of members likelihood
REAL, INTENT(INOUT)               :: PEFF_LOC
!
REAL, DIMENSION(SIZE(PWEIGHT,1), SIZE(PWEIGHT,2)) :: ZWEIGHT ! local PWEIGHT
REAL :: ZTOL, ZALPHA_0, ZEFF_LOC, ZEFF_LOC_1, ZEFF_LOC_2,ZALPHA, ZALPHA_1, ZALPHA_2
REAL :: ZL, ZH ! lower and upper boundaries such as f(zl)<0 and f(zh) >0
!
LOGICAL :: GCOND_NK
!
INTEGER :: JITER
!
!
!*********************************************************************
    ! INITIALIZATION
    GCOND_NK= .TRUE. ! get sure you keep at least Nk particles
    ZALPHA_1 = 0.0000001*PALPHA ! second initial value needed for the secant method
    ZEFF_LOC = PEFF_LOC
    ZALPHA = PALPHA
    JITER=0
    ! bisection init
    ZL = ZALPHA ! f is a decreasing function (f(1)<0, f(0)>0)
    ZH = ZALPHA_1
    !
    ! evaluate the function on the two provided values
    ! first value (x0) has just been computed
    DO WHILE ((GCOND_NK) .AND. (JITER <1000))
      JITER=JITER+1
      ! the second value is necessary only for the first call (otherwise it has already bee computed)
      IF(JITER <2) THEN
        ! second value (x1)
        CALL SET_OBS_ERROR(PRINV, KMASK, KCOUNT_LOC, LDIAGONAL_R, PALPHA = ZALPHA_1)
        CALL LIKELIHOOD(PINNOV,PRINV, PLIKE_ENS)
        CALL NORMALIZE(PLIKE_ENS,ZWEIGHT)
        CALL EFFWEIGHTS(ZWEIGHT, ZEFF_LOC_1)
      ENDIF
      !
      ! if f(x1) = f(x0), this is pathological and we should exit.
      !pathological case 1: if the function is flat, bisect.
      IF (ABS(ZEFF_LOC_1 - ZEFF_LOC) < 1E-12) THEN
        ZALPHA_2 = ZL + 0.5*(ZH - ZL)
      ELSE 
        ! secant(x2). NEFF_PF simplifies in the denominator.
        ZALPHA_2 = ZALPHA_1 - (ZEFF_LOC_1 - NEFF_PF) * (ZALPHA_1 - ZALPHA) / (ZEFF_LOC_1 - ZEFF_LOC)
        !
        ! pathological case 2: overshooting might occur : if so, jump to a bisection method.
        IF ((ZALPHA_2 < ZL) .OR. (ZALPHA_2 > ZH)) THEN
          ZALPHA_2 = ZL + 0.5*(ZH - ZL)
        ENDIF
      ENDIF

      !
      ! evaluate the function and continue if necessary
      CALL SET_OBS_ERROR(PRINV, KMASK, KCOUNT_LOC, LDIAGONAL_R, PALPHA = ZALPHA_2)
      CALL LIKELIHOOD(PINNOV,PRINV, PLIKE_ENS)
      CALL NORMALIZE(PLIKE_ENS, ZWEIGHT)
      CALL EFFWEIGHTS(ZWEIGHT, ZEFF_LOC_2)
      ! convergence condition on the function value
      IF (ABS(ZEFF_LOC_2 - NEFF_PF) < 1E-4) THEN
        GCOND_NK = .FALSE.
        ! ZALPHA2 is the solution : return it
        PALPHA = ZALPHA_2
        PWEIGHT = ZWEIGHT
        PEFF_LOC = ZEFF_LOC_2
      ELSE
        ZALPHA = ZALPHA_1
        ZALPHA_1=ZALPHA_2
        
        ZEFF_LOC = ZEFF_LOC_1        
        ZEFF_LOC_1 = ZEFF_LOC_2
        ! bisection : update the brackets
        IF ((ZEFF_LOC_2 - NEFF_PF) < 0.) THEN
          ZL = ZALPHA_2
        ELSE
          ZH = ZALPHA_2
        ENDIF
      ENDIF
    END DO
    IF(JITER == 1000.) THEN
      CALL ABOR1_SFX('INFLATE : failed to converge.')
    ENDIF
END SUBROUTINE INFLATE
!********************************************************************
!********************************************************************
!********************************************************************
SUBROUTINE EFFWEIGHTS(PWEIGHT, PEFF_LOC)
!********************************************************************
! -Purpose
!    Compute efficient weights (local) (Liu and Chen 1995)
! -@TODO:
!    Chech for global case.
! -Author
!     B. Cluzet 22/05/19
!********************************************************************
USE MODD_ASSIM, ONLY : NENS
IMPLICIT NONE

REAL, DIMENSION(:,:), INTENT(IN) :: PWEIGHT
REAL, INTENT(OUT) :: PEFF_LOC

REAL :: ZSQ_SUM
INTEGER :: JENS
!*********************************************************************
PEFF_LOC = 0.
ZSQ_SUM = 0.
!DO JP=1,SIZE(PLIKE_ENS, 1) ! F... the patches
DO JENS=1, NENS
  ZSQ_SUM = ZSQ_SUM + PWEIGHT(1, JENS) * PWEIGHT(1, JENS)
END DO
!
! ZSQ_SUM is always >0. (see weighting)
PEFF_LOC = 1./ ZSQ_SUM
! END DO
END SUBROUTINE EFFWEIGHTS
!********************************************************************
!********************************************************************
!********************************************************************
SUBROUTINE RESAMPLING(PWEIGHT, KPOINT)
!********************************************************************
! -Purpose
!    Resample (local or global)
! -Author
!     B. Cluzet 24/04/18
! -Modifications
!     B. Cluzet 24/04/18 : dev mutualization between local/global
!********************************************************************

USE MODD_ASSIM, ONLY : NENS, LGLOBAL_PF, NPFPART
IMPLICIT NONE

REAL, DIMENSION(:,:), INTENT(IN)    :: PWEIGHT    !(NPATCHES, NENS) members weights, f... the patches
INTEGER, INTENT(IN), OPTIONAL       :: KPOINT     !id of point at which to apply the filter (localised only)

!- Internal vars
!REAL, DIMENSION(:)                :: ZRAND        ! random draw for resampling f... the patches
REAL :: ZRAND
REAL, DIMENSION(SIZE(PWEIGHT, 1), SIZE(PWEIGHT,2)) :: ZWEIGHTCUMUL ! vector of cumulations of weights, f... the patches
INTEGER :: JP, JENS, JID, JRANK
!
!*********************************************************************
!
IF (.NOT. LGLOBAL_PF) THEN
  JID=KPOINT
ELSE
  JID=1
ENDIF
!
!**** Resampling (Kitagawa) (Kitagawa et al., 1996)*****
!
! memberwise cumulation (after normalisation)
DO JENS = 1, NENS
  ZWEIGHTCUMUL(:,JENS) = SUM(PWEIGHT(:,1:JENS),2)
END DO
!
! random draw in ]0, 1/NENS[
CALL RANDOM_NUMBER(ZRAND)
ZRAND = ZRAND/NENS
!
! NENS draws of ZRAND will give JRANK, the ids of particules to replicate at this location
DO JENS = 1, NENS
  !
  !"searchsorted function" : find JRANK, the rank of ZRAND in ZWEIGHT_CUMUL
  DO JP=1, SIZE(PWEIGHT,1)  ! f... the patches
    IF (ZRAND < ZWEIGHTCUMUL(1,1)) THEN
      NPFPART(JID,JP,JENS) = 1
    ENDIF
    DO JRANK = 2, NENS
      IF (ZWEIGHTCUMUL(JP,JRANK-1) < ZRAND .AND. ZRAND < ZWEIGHTCUMUL(JP,JRANK)) THEN
        NPFPART(JID,JP,JENS) = JRANK
      ENDIF
    END DO
  END DO
  !
  ! increment ZRAND
  !DO JP=1, IM%O%NPATCH !  f... the patches
  ZRAND = ZRAND + 1./NENS
  !END DO
END DO

!######## VERY IMPORTANT RESORTING STEP !!!!
CALL RESORT(JID)
!################################
END SUBROUTINE RESAMPLING
!********************************************************************
!********************************************************************
!********************************************************************
SUBROUTINE RESORT(KID)
!********************************************************************
! -Purpose
!    Re-sort the NPFPART to ensure that ith-element of NPFPART equals i
!    for each replicated particle
!    /!\ assuming that NPFPART is SORTED.
!    ex1 :(/1,1,2,3,3,3,8,8,9,9,9,9,9,16,16,16/)
!      =>(/1,2,3,1,3,3,8,8,9,9,9,9,9,16,16,16/)
!    ex2 :(/9,9,9,9,13,14,15,15,16,16,16,16,16,16,16,16/)
!      =>(/9,9,9,15,16,16,16,16,9,16,16,16,13,14,15,16/)
!
! -Author
!     B. Cluzet 23/06/19
!********************************************************************
USE MODD_ASSIM, ONLY : NPFPART, LGLOBAL_PF
IMPLICIT NONE
INTEGER, INTENT(IN) :: KID
INTEGER, DIMENSION(SIZE(NPFPART,1), SIZE(NPFPART,2), SIZE(NPFPART,3)) :: NPFPART_2
INTEGER, DIMENSION(SIZE(NPFPART, 3)) :: IUN  ! unique replicated particles at point i.
INTEGER, DIMENSION(SIZE(NPFPART, 3)) :: IUNIND  ! loc of first ocur. of unique particles in NPFPART

INTEGER, DIMENSION(SIZE(NPFPART, 3)) :: INOT_UNIND  ! locs of not unique particles (complementary of iunind)
INTEGER, DIMENSION(SIZE(NPFPART, 3)) :: IFREE_IND ! free indices (indices of not rep particles)
INTEGER :: JENS, JCT_UN, JCT_NOT_UNIND, JUN, JCT_FREE
INTEGER :: ISTAT

NPFPART_2 = 0
!

  ! in each point, find the unique elements 
  ! /!\ assuming NPFPART is SORTED !!!
  IUN(:) = 0
  IUN(1) = NPFPART(KID,1,1)
  IUNIND(:) =  0
  IUNIND(1) = 1

  INOT_UNIND(:) = 0

  JCT_UN = 2
  JCT_NOT_UNIND = 1
  ! find unique/not unique particles
  DO JENS = 2, SIZE(NPFPART, 3)
    IF (NPFPART(KID,1, JENS) > IUN(JCT_UN-1)) THEN
      IUN(JCT_UN) = NPFPART(KID, 1, JENS)
      IUNIND(JCT_UN) = JENS
      JCT_UN = JCT_UN + 1
    ELSE
      INOT_UNIND(JCT_NOT_UNIND) = JENS
      JCT_NOT_UNIND = JCT_NOT_UNIND + 1
    ENDIF
  END DO
  !
  JCT_UN = JCT_UN-1
  JCT_NOT_UNIND = JCT_NOT_UNIND -1
  !
  ! find free indices
  IFREE_IND(:) = 0
  JCT_FREE = 0
  DO JENS = 1, SIZE(NPFPART,3)
    IF (.NOT. (ANY(NPFPART(KID,1,:) == JENS)))THEN
      JCT_FREE = JCT_FREE + 1
      IFREE_IND(JCT_FREE)=JENS
    END IF
  END DO
  ! 
  !put not unique to free
  DO JUN = 1, JCT_NOT_UNIND
    NPFPART_2(KID,1, IFREE_IND(JUN)) = NPFPART(KID,1,INOT_UNIND(JUN))
  END DO
  !
  !put unique where they belong
  DO JUN=1, JCT_UN
    NPFPART_2(KID,1,IUN(JUN))=IUN(JUN)
  END DO
!
! actualize NPFPART
DO JENS=1, SIZE(NPFPART,3)
  NPFPART(KID,1,JENS) = NPFPART_2(KID,1,JENS)
END DO
END SUBROUTINE RESORT
!********************************************************************
!********************************************************************
!********************************************************************
SUBROUTINE COMPUTE_DIST_DEGREES(KPT1, KPT2, PLON_IN, PLAT_IN, PDEG)
!********************************************************************
! -Purpose
!    Compute angular distance (degrees)
!    between PT1 and PT2 of the LATLON grid
!
! -Author
!     B. Cluzet 08/06/20
!********************************************************************
USE MODD_CSTS, ONLY : XPI
USE PARKIND1,  ONLY : JPRB
!
IMPLICIT NONE
!
INTEGER, INTENT(IN) :: KPT1, KPT2     ! ids of Points 1 and 2 in the geometry
REAL(KIND=JPRB), DIMENSION (:), INTENT(IN) ::  PLON_IN  ! coordinate vectors of the geometry
REAL(KIND=JPRB), DIMENSION (:), INTENT(IN) ::  PLAT_IN
REAL, INTENT(OUT) :: PDEG           ! distance in degrees
!
REAL :: ZLAT1, ZLAT2, ZLON1, ZLON2  ! point coordinates in degrees
ZLAT1 = PLAT_IN(KPT1) * XPI / 180.
ZLAT2 = PLAT_IN(KPT2) * XPI / 180.
ZLON1 = PLON_IN(KPT1) * XPI / 180.
ZLON2 = PLON_IN(KPT2) * XPI / 180.
!
! In the semi-distributed geometry, lat/lon of all points are equal.
IF ((ZLAT1 - ZLAT2 <1E-7) .AND. (ZLON1 - ZLON2 <1E-7)) THEN
  PDEG = 0.
ELSE
  ! simple haversine formula (rounding error is not an issue).
  PDEG = ACOS(COS(ZLAT1) * COS(ZLAT2) * COS(ZLON1 - ZLON2) + SIN(ZLAT1) * SIN(ZLAT2))
  PDEG = PDEG * 180. / XPI
ENDIF
!
END SUBROUTINE COMPUTE_DIST_DEGREES
!********************************************************************
!********************************************************************
!********************************************************************
SUBROUTINE IN_RADIUS(LCOND_RLOC, KPT, KPTOBS, PLON_IN, PLAT_IN, PDEG)
!********************************************************************
! -Purpose
!    Check if the location of obs. point KPTOBS 
!    is within a disk of ZDEG RADIUS around KPT.
!    /!\ this scripts rejects the center of the disk if LLOO_PF = .TRUE.
!
! -Author
!     B. Cluzet 27/08/20
!********************************************************************
USE MODD_ASSIM, ONLY : XDLOC_PF, LLOO_PF
USE PARKIND1,  ONLY : JPRB
!
IMPLICIT NONE
! 
LOGICAL, INTENT(OUT) :: LCOND_RLOC  ! True if inside the disk (except if at the center in the LLOO_PF).
INTEGER, INTENT(IN)  :: KPT, KPTOBS ! id of the disk center and the observed point
REAL(KIND=JPRB), DIMENSION (:), INTENT(IN) ::  PLON_IN  ! coordinate vectors of the geometry
REAL(KIND=JPRB), DIMENSION (:), INTENT(IN) ::  PLAT_IN
REAL, INTENT(OUT) :: PDEG           ! distance in degrees
!
LCOND_RLOC = .FALSE.
IF (KPT .NE. KPTOBS) THEN
  CALL COMPUTE_DIST_DEGREES(KPT, KPTOBS, PLON_IN, PLAT_IN, PDEG)
  IF (PDEG < XDLOC_PF) THEN
    LCOND_RLOC =.TRUE.
  ENDIF
ELSEIF (.NOT. LLOO_PF) THEN ! if LOO_PF we want to reject local obs.
  LCOND_RLOC = .TRUE.
ENDIF

END SUBROUTINE IN_RADIUS
!********************************************************************
!********************************************************************
!********************************************************************
END MODULE MODE_PF
