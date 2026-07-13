!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!     #########
      SUBROUTINE HOR_EXTRAPOL_SURF(KLUOUT,HCOORTYPE,KILEN,PILA1,PILA2,PILO1,PILO2,&
                                   KINLA,KINLO,KP,PFIELD_IN,PLAT,PLON,PFIELD,OINTERP,&
                                   PILATARRAY)
!     ###################################################################
!
!!**** *HOR_EXTRAPOL_SURF* extrapolate a surface field
!!
!!    PURPOSE
!!    -------
!!
!!    METHOD
!!    ------
!!       For each point to interpolate, the nearest valid point value is set.
!!
!!    EXTERNAL
!!    --------
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!    REFERENCE
!!    ---------
!!
!!    AUTHOR
!!    ------
!!
!!    V. Masson          Meteo-France
!!
!!    MODIFICATION
!!    ------------
!!
!!    Original     01/12/98
!!     V. Masson    01/2004 extrapolation in latitude and longitude
!!     M. Jidane    11/2013 add OpenMP directives
!!     Q. Rodier    06/2021 avoid abort for interpolation of ALL(PFIELD)=XUNDEF with ECOSG
!!     A. Napoly    10/2022 add OpenMP directives and optimisations in loops
!----------------------------------------------------------------------------
!
!*    0.     DECLARATION
!            -----------
!
USE MODD_SURFEX_MPI, ONLY : NRANK, NPROC, NPIO, NCOMM, IDX_I
USE MODD_SURF_PAR,   ONLY : XUNDEF
USE MODD_CSTS,       ONLY : XPI
USE MODN_PREP_SURF_ATM, ONLY : NHALO_PREP
!
USE MODI_ABOR1_SFX
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB

IMPLICIT NONE
!
#ifdef SFX_MPI
INCLUDE "mpif.h"
#endif
!
!*    0.1    Declaration of arguments
!            ------------------------
!
INTEGER,               INTENT(IN)     :: KLUOUT,KILEN,KINLA
CHARACTER(LEN=4),      INTENT(IN)     :: HCOORTYPE! type of coordinate
REAL, INTENT(IN) :: PILA1,PILA2,PILO1,PILO2
INTEGER, DIMENSION(:), INTENT(IN) :: KINLO
INTEGER, DIMENSION(:,:), INTENT(IN) :: KP
REAL,   DIMENSION(:,:),  INTENT(IN)     :: PFIELD_IN! input field on grid mesh
REAL,   DIMENSION(:),  INTENT(IN)     :: PLAT     ! latitude of each grid mesh.
REAL,   DIMENSION(:),  INTENT(IN)     :: PLON     ! longitude of each grid mesh.
REAL,   DIMENSION(:,:),  INTENT(INOUT)  :: PFIELD   ! field on grid mesh
LOGICAL,DIMENSION(:),  INTENT(IN)     :: OINTERP  ! .true. where physical value is needed
REAL, DIMENSION(:), INTENT(IN), OPTIONAL :: PILATARRAY
!
!*    0.2    Declaration of local variables
!            ------------------------------
!
integer,parameter :: iint4=kind(0)/4,ireal4=kind(0.)/4
INTEGER, DIMENSION(:), ALLOCATABLE :: IMASK
INTEGER, DIMENSION(:,:), ALLOCATABLE :: ipos
INTEGER  :: INO     ! output array size
INTEGER, DIMENSION(2) :: ITSIZE
INTEGER :: isize(2,0:NPROC-1),ibor(0:nproc-1,2)
INTEGER :: J, ID0, ICOMPT, nxtr,n1,n2
INTEGER :: INFOMPI, IDX, INL
INTEGER  :: JI, JL, JLAT, JLON, JP   ! loop index on points
INTEGER  :: JISC  ! loop index on valid points
INTEGER, DIMENSION(KINLA) :: ioff
INTEGER, DIMENSION(SIZE(KP,2)) :: IKP
#ifdef SFX_MPI
INTEGER, DIMENSION(MPI_STATUS_SIZE) :: ISTATUS
#endif
LOGICAL  :: GLALO ! flag true is second coordinate is a longitude or pseudo-lon.
                  !      false if metric coordinates
REAL, DIMENSION(:), ALLOCATABLE :: ZTLONMIN, ZTLONMAX, ZTLATMIN, ZTLATMAX
REAL, DIMENSION(:,:), ALLOCATABLE :: ZFIELD,ZNDIST
REAL :: ZLAT  ! latitude of point to define
REAL :: ZLON  ! longitude of point to define
REAL :: ZDIST ! current distance to valid point (in lat/lon grid)
REAL :: ZDLONSC! longitude of valid point
REAL :: ZIDLOMAX, ZIDLOMIN, ZIDLAMAX, ZIDLAMIN
REAL, DIMENSION(:,:), ALLOCATABLE :: ZCOOR
REAL, DIMENSION(:), ALLOCATABLE :: ZIDLA,ZCOSLA
REAL, DIMENSION(:), ALLOCATABLE :: ZLA       ! input "latitude"  coordinate
REAL, DIMENSION(:), ALLOCATABLE :: ZLO       ! input "longitude" coordinate
REAL :: ZLONN,ZLONX,ZLATN,ZLATX,ZDLON,ZDLAT
REAL, DIMENSION(KINLA) :: ZINLAT,ZIDLO
REAL(KIND=JPRB) :: ZRAD ! conversion degrees to radians
REAL(KIND=JPRB) :: ZHOOK_HANDLE,ZH

IF (LHOOK) CALL DR_HOOK('HOR_EXTRAPOL_SURF',0,ZHOOK_HANDLE)

IF (LHOOK) CALL DR_HOOK("HOR_EXTRAPOL_SURF_LOLA",0,ZH)

GLALO = HCOORTYPE=='LALO'

IF (PRESENT(PILATARRAY)) THEN
  ALLOCATE(ZIDLA(KINLA))

  ZIDLA(1) = 0.
  DO JLAT=2,KINLA
    ZIDLA(JLAT) = PILATARRAY(JLAT) - PILATARRAY(JLAT-1)
  ENDDO

  ZINLAT(1) = PILA1
  DO JLAT=2,KINLA
    ZINLAT(JLAT) = ZINLAT(JLAT-1)+ZIDLA(JLAT)
  ENDDO

  ZIDLAMAX = MAXVAL(ABS(ZIDLA))
  ZIDLAMIN = MINVAL(ABS(ZIDLA(2:KINLA)))

  DEALLOCATE(ZIDLA)
ELSE
  ZIDLAMIN = (PILA2-PILA1)/(KINLA-1)

  DO JLAT=1,KINLA
    ZINLAT(JLAT) = PILA1+(JLAT-1)*ZIDLAMIN
  ENDDO

  ZIDLAMIN = ABS(ZIDLAMIN)
  ZIDLAMAX = ZIDLAMIN
ENDIF

IF (KINLA >= 1) THEN
  IOFF(1) = 0
  DO JLAT = 1, KINLA - 1
    IOFF(JLAT + 1) = IOFF(JLAT) + KINLO(JLAT)
  ENDDO
ENDIF

ALLOCATE(ZLA(KILEN),ZLO(KILEN))

IF (GLALO) THEN
  ZIDLO(:) = (PILO2-PILO1)/KINLO(1:KINLA)
ELSE
  ZIDLO(:) = (PILO2-PILO1)/(KINLO(1:KINLA)-1)
END IF

DO JLAT=1,KINLA
  ZLA(IOFF(JLAT)+1:IOFF(JLAT)+KINLO(JLAT)) = ZINLAT(JLAT)

  DO JLON=1,KINLO(JLAT)
    ZLO(IOFF(JLAT)+JLON) = PILO1+(JLON-1)*ZIDLO(JLAT)
  END DO
END DO

ZIDLOMAX = MAXVAL(ABS(ZIDLO))
ZIDLOMIN = MINVAL(ABS(ZIDLO))

IF (LHOOK) CALL DR_HOOK("HOR_EXTRAPOL_SURF_LOLA",1,ZH)

!-------------------------------------------------------------------------------
!
!*      4.   Loop on points to define
!            ------------------------
!
IF (LHOOK) CALL DR_HOOK("HOR_EXTRAPOL_SURF_MASK",0,ZH)

INO = SIZE(PFIELD,1)

NXTR = 0
ALLOCATE(IMASK(INO))
IMASK(:)=0
DO JI=1,INO
  IF (.NOT.OINTERP(JI).OR.ALL(PFIELD(JI,:) /= XUNDEF)) CYCLE

  NXTR = NXTR + 1
  IMASK(NXTR) = JI
ENDDO

IF (LHOOK) CALL DR_HOOK("HOR_EXTRAPOL_SURF_MASK",1,ZH)

ZRAD=XPI/180

!1: ZTLONMIN, ZTLONMAX, ZTLATMIN, ZTLATMAX contain for each point to extrapol 
! the limits of the domain where to search for the valid points, according
! to NHALO_PREP

IF (LHOOK) CALL DR_HOOK("HOR_EXTRAPOL_SURF_ZNX",0,ZH)

ITSIZE(1) = NXTR
ALLOCATE(ZTLONMIN(NXTR),ZTLONMAX(NXTR),ZTLATMIN(NXTR),ZTLATMAX(NXTR))

IF (NHALO_PREP > 0) THEN
  ZDLON=ZIDLOMAX*NHALO_PREP
  ZDLAT=ZIDLAMAX*NHALO_PREP

  !$OMP PARALLEL DO PRIVATE(IKP)
  DO JP=1,NXTR
    IKP(:) = KP(IMASK(JP),:)
    ZTLONMIN(JP) = MINVAL(ZLO(IKP))-ZDLON
    ZTLONMAX(JP) = MAXVAL(ZLO(IKP))+ZDLON
    ZTLATMIN(JP) = MINVAL(ZLA(IKP))-ZDLAT
    ZTLATMAX(JP) = MAXVAL(ZLA(IKP))+ZDLAT
  ENDDO
  !$OMP END PARALLEL DO

  ITSIZE(2) = MAXVAL(CEILING((ZTLONMAX-ZTLONMIN+1)/ZIDLOMIN)*&
    CEILING((ZTLATMAX-ZTLATMIN+1)/ZIDLAMIN))
ELSE
  ZLONN=MINVAL(ZLO(:))
  ZLONX=MAXVAL(ZLO(:))
  ZLATN=MINVAL(ZLA(:))
  ZLATX=MAXVAL(ZLA(:))

  !$OMP PARALLEL DO
  DO JP=1,NXTR
    ZTLONMIN(JP) = ZLONN
    ZTLONMAX(JP) = ZLONX
    ZTLATMIN(JP) = ZLATN
    ZTLATMAX(JP) = ZLATX
  ENDDO
  !$OMP END PARALLEL DO

  ITSIZE(2) = CEILING((ZLONX-ZLONN+1)/ZIDLOMIN)*CEILING((ZLATX-ZLATN+1)/ZIDLAMIN)
ENDIF

IF (ITSIZE(2) < 0) ITSIZE(2) = 0

IF (LHOOK) CALL DR_HOOK("HOR_EXTRAPOL_SURF_ZNX",1,ZH)

!NPIO knows the numbers of points to extrapolate for all tasks
IF (NPROC > 1) THEN
#ifdef SFX_MPI
  CALL MPI_GATHER(ITSIZE,2*IINT4,MPI_INTEGER,ISIZE,2*IINT4,MPI_INTEGER,NPIO,NCOMM,INFOMPI)
  IBOR(:,1) = ISIZE(1,:)
  IBOR(:,2) = ISIZE(2,:)
#endif
ELSE
  IBOR(0,:) = ITSIZE(:)
ENDIF

IF (LHOOK) CALL DR_HOOK("HOR_EXTRAPOL_SURF_POS",0,ZH)
IF (NRANK == NPIO) THEN
  ALLOCATE(IPOS(MAXVAL(IBOR(:,2)),MAXVAL(IBOR(:,1))))
  ALLOCATE(ZCOOR(MAXVAL(IBOR(:,1)),2))
ELSE
  ALLOCATE(IPOS(ITSIZE(2),NXTR))
  ALLOCATE(ZCOOR(NXTR,2))
ENDIF

! useless:
!ipos(:,:) = 0
!ZCOOR(:,:) = 0

!$OMP PARALLEL DO PRIVATE(ICOMPT,JLAT,JLON,ZLON,J)
DO JP=1,NXTR
  !coordinates of the point in the grid 
  ZCOOR(JP,1) = PLAT(IMASK(JP))
  ZCOOR(JP,2) = PLON(IMASK(JP))

  ICOMPT = 0

  !loop on the grid
  DO JLAT = 1,KINLA
    IF (ZINLAT(JLAT) < ZTLATMIN(JP).OR.ZINLAT(JLAT) > ZTLATMAX(JP)) CYCLE

    DO JLON = 1,KINLO(JLAT)
      J = IOFF(JLAT)+JLON
      IF (ZTLONMIN(JP) <= ZLO(J).AND.ZLO(J) <= ZTLONMAX(JP)) THEN
        ICOMPT = ICOMPT + 1
        !ipos: indexes of the points needed to interpolate in the complete grid
        IPOS(ICOMPT,JP) = J
      ENDIF
    ENDDO
  ENDDO

  IF (ICOMPT < ITSIZE(2)) IPOS(ICOMPT+1,JP) = 0
ENDDO
!$OMP END PARALLEL DO

DEALLOCATE(ZTLONMIN,ZTLONMAX,ZTLATMIN,ZTLATMAX)
IF (LHOOK) CALL DR_HOOK("HOR_EXTRAPOL_SURF_POS",1,ZH)

INL = SIZE(PFIELD,2)

IF (NRANK == NPIO) THEN
  IF (LHOOK) CALL DR_HOOK("HOR_EXTRAPOL_SURF_MINVAL",0,ZH)
  DO JP=NPIO,NPROC-1+NPIO
    J = JP
    IF (JP > NPROC-1) J = JP-NPROC

    IF (IBOR(J,1) == 0) CYCLE

    N1 = IBOR(J,1)
    N2 = IBOR(J,2)

    ALLOCATE(ZFIELD(INL,N1),ZNDIST(INL,N1))

    ZFIELD(:,:) = XUNDEF

    IF (J /= NPIO) THEN
      ! receive positions and coordinates
      ! note: NPIO values get lost, this is why first jp is NPIO
#ifdef SFX_MPI        
      CALL MPI_RECV(IPOS(1:N2,1:N1),N1*N2*IINT4,MPI_INTEGER,J,IDX_I+1,NCOMM,ISTATUS,INFOMPI)
      CALL MPI_RECV(ZCOOR(1:N1,:),2*N1*IREAL4,MPI_REAL,J,IDX_I+2,NCOMM,ISTATUS,INFOMPI)
#endif
    ENDIF

    IF (GLALO) THEN
      ALLOCATE(ZCOSLA(N1))
      ZCOSLA(:) = COS(ZCOOR(1:N1,1)*ZRAD)

      !$OMP PARALLEL DO PRIVATE(JISC,JL,ID0,ZDLONSC,ZDIST)
      DO JI=1,N1
        ZNDIST(:,JI) = HUGE(0.)

        DO JISC=1,N2
          !index in the whole grid of the point used to interpolate
          ID0 = IPOS(JISC,JI)
          IF (ID0 == 0) EXIT

          IF (ALL(PFIELD_IN(ID0,:) == XUNDEF)) CYCLE

          ZDLONSC = ZLO(ID0)-ZCOOR(JI,2)

          IF (ZDLONSC > 180) THEN
             ZDLONSC = ZDLONSC - 360
          ELSE IF (ZDLONSC < -180) THEN
             ZDLONSC = ZDLONSC + 360
          END IF

          ZDIST = (ZLA(ID0)-ZCOOR(JI,1))**2+(ZDLONSC*ZCOSLA(JI))**2

          DO JL=1,INL
            IF (ZNDIST(JL,JI) >= ZDIST.AND.PFIELD_IN(ID0,JL) /= XUNDEF) THEN
              ZFIELD(JL,JI) = PFIELD_IN(ID0,JL)
              ZNDIST(JL,JI) = ZDIST
            ENDIF
          ENDDO
        END DO   
      ENDDO
      !$OMP END PARALLEL DO

      DEALLOCATE(ZCOSLA)
    ELSE
      !$OMP PARALLEL DO PRIVATE(JISC,JL,ID0,ZDIST)
      DO JI=1,N1
        ZNDIST(:,JI) = HUGE(0.)

        DO JISC=1,N2
          !index in the whole grid of the point used to interpolate
          ID0 = IPOS(JISC,JI)
          IF (ID0 == 0) EXIT

          IF (ALL(PFIELD_IN(ID0,:) == XUNDEF)) CYCLE

          ZDIST = (ZLA(ID0)-ZCOOR(JI,1))**2+(ZLO(ID0)-ZCOOR(JI,2))**2

          DO JL=1,INL
            IF (ZNDIST(JL,JI) >= ZDIST.AND.PFIELD_IN(ID0,JL) /= XUNDEF) THEN
              ZFIELD(JL,JI) = PFIELD_IN(ID0,JL)
              ZNDIST(JL,JI) = ZDIST
            ENDIF
          ENDDO
        END DO   
      ENDDO
      !$OMP END PARALLEL DO
    ENDIF

    IF (J == NPIO) THEN
      ! note: n1 is nxtr (NPIO's one) since j is NPIO
      DO JI=1,N1
        PFIELD(IMASK(JI),:) = ZFIELD(:,JI)
      ENDDO
    ELSE
      !send values found to extrapolate
#ifdef SFX_MPI        
      CALL MPI_SEND(ZFIELD,SIZE(ZFIELD)*IREAL4,MPI_REAL,J,IDX_I+3,NCOMM,INFOMPI)
#endif
    ENDIF

    ! optim: deallocate zndist after send
    DEALLOCATE(ZFIELD,ZNDIST)
  ENDDO
  IF (LHOOK) CALL DR_HOOK("HOR_EXTRAPOL_SURF_MINVAL",1,ZH)
ELSE IF (NXTR > 0) THEN
#ifdef SFX_MPI    
  IDX = IDX_I + 1
  CALL MPI_SEND(IPOS,SIZE(IPOS)*IINT4,MPI_INTEGER,NPIO,IDX,NCOMM,INFOMPI)

  IDX = IDX_I + 2
  CALL MPI_SEND(ZCOOR,2*NXTR*IREAL4,MPI_REAL,NPIO,IDX,NCOMM,INFOMPI)

  ! optim: allocate after sends
  ALLOCATE(ZFIELD(INL,NXTR))

  IDX = IDX_I + 3
  CALL MPI_RECV(ZFIELD,INL*NXTR*IREAL4,MPI_REAL,NPIO,IDX,NCOMM,ISTATUS,INFOMPI)
#else
  ALLOCATE(ZFIELD(INL,NXTR))
#endif

  DO JI=1,NXTR
    PFIELD(IMASK(JI),:) = ZFIELD(:,JI)
  ENDDO

  DEALLOCATE(ZFIELD)
ENDIF

IDX_I = IDX_I + 3

DEALLOCATE(ZCOOR,ZLA,ZLO)
DEALLOCATE(IPOS,IMASK)

DO JL=1,INL
  IF (ANY(PFIELD(:,JL) == XUNDEF.AND.OINTERP(:)).AND..NOT.ALL(PFIELD(:,JL)==XUNDEF)) THEN
    WRITE(*,*) 'LAYER ',JL,': NO EXTRAPOLATION : INCREASE YOUR HALO_PREP IN NAM_PREP_SURF_ATM'
    CALL ABOR1_SFX('NO EXTRAPOLATION : INCREASE YOUR HALO_PREP IN NAM_PREP_SURF_ATM')
  ENDIF
  WHERE (.NOT.OINTERP(:)) PFIELD(:,JL) = XUNDEF
ENDDO

IF (LHOOK) CALL DR_HOOK("HOR_EXTRAPOL_SURF",1,ZHOOK_HANDLE)
END SUBROUTINE HOR_EXTRAPOL_SURF
