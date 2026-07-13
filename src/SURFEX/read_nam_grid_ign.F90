!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!     ################################################################
      SUBROUTINE READ_NAM_GRID_IGN(PGRID_FULL_PAR,KDIM_FULL,HPROGRAM,KGRID_PAR,KL,PGRID_PAR,HDIR)
!     ################################################################
!
!!****  *READ_NAM_GRID_IGN* - routine to read in namelist the horizontal grid
!!
!!    PURPOSE
!!    -------
!!
!!**  METHOD
!!    ------
!!
!!    EXTERNAL
!!    --------
!!
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!    REFERENCE
!!    ---------
!!
!!
!!    AUTHOR
!!    ------
!!      E. Martin   *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    10/2007 
!!      07/2011     add maximum domain dimension for output (B. Decharme)
!!      08/2022     When the grid is prescribed from LL point and nrows / ncols,
!!                  remove the dependency to the fixed allocation of XX and YY
!!                  by assigning directly ZX and ZY (M. Lafaysse)
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_SURFEX_MPI, ONLY : NRANK, NSIZE_TASK
!
USE MODD_SURF_PAR, ONLY : XUNDEF
!
USE MODE_POS_SURF
!
USE MODI_OPEN_NAMELIST
USE MODI_CLOSE_NAMELIST
USE MODI_GET_LUOUT
USE MODI_TEST_NAM_VAR_SURF
!
USE MODE_GRIDTYPE_IGN
USE MODI_GET_XYALL_IGN
!
USE MODI_READ_AND_SEND_MPI
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*       0.1   Declarations of arguments
!              -------------------------
!
REAL, DIMENSION(:), POINTER :: PGRID_FULL_PAR
INTEGER, INTENT(IN) :: KDIM_FULL
!
CHARACTER(LEN=6),           INTENT(IN)    :: HPROGRAM   ! calling program
INTEGER,                    INTENT(INOUT) :: KGRID_PAR  ! size of PGRID_PAR
INTEGER,                    INTENT(OUT)   :: KL         ! number of points
REAL, DIMENSION(KGRID_PAR), INTENT(OUT)   :: PGRID_PAR  ! parameters defining this grid
CHARACTER(LEN=1), INTENT(IN) :: HDIR
!
!*       0.2   Declarations of local variables
!              -------------------------------
!
INTEGER :: ILUOUT ! output listing logical unit
INTEGER :: ILUNAM ! namelist file  logical unit
INTEGER :: ILAMBERT ! Lambert type

REAL, DIMENSION(:),   ALLOCATABLE :: ZX, ZX0       ! X conformal coordinate of grid mesh
REAL, DIMENSION(:),   ALLOCATABLE :: ZY, ZY0       ! Y conformal coordinate of grid mesh
REAL, DIMENSION(:),   ALLOCATABLE :: ZDX, ZDX0     ! X grid mesh size
REAL, DIMENSION(:),   ALLOCATABLE :: ZDY, ZDY0     ! Y grid mesh size
!
!*       0.3   Declarations of namelist
!              ------------------------
!
CHARACTER(LEN=3) :: CLAMBERT  ! Lambert type
INTEGER :: NPOINTS  ! number of points
REAL, DIMENSION(1000000) :: XX  ! X coordinate of grid mesh center (in meters)
REAL, DIMENSION(1000000) :: XY  ! Y coordinate of grid mesh center (in meters)
REAL, DIMENSION(1000000) :: XDX ! X mesh size (in meters)
REAL, DIMENSION(1000000) :: XDY ! Y mesh size (in meters)
!
REAL :: XX_LLCORNER ! X coordinate of left  side of the domain
REAL :: XY_LLCORNER ! Y coordinate of lower side of the domain
REAL :: XCELLSIZE   ! size of the cell (equal in X and Y)
INTEGER :: NCOLS    ! number of columns
INTEGER :: NROWS    ! number of rows
!
REAL, DIMENSION(:), ALLOCATABLE :: ZXALL  ! maximum domain X coordinate of grid mesh
REAL, DIMENSION(:), ALLOCATABLE :: ZYALL  ! maximum domain Y coordinate of grid mesh
INTEGER                         :: IDIMX  ! maximum domain length in X
INTEGER                         :: IDIMY  ! maximum domain length in Y
!
REAL, DIMENSION(:), POINTER :: ZGRID_PAR
!
INTEGER :: JCOLS, JROWS, IINDEX ! loop counters
LOGICAL :: GFOUND
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
NAMELIST/NAM_IGN/CLAMBERT,NPOINTS,XX,XY,XDX,XDY,      &
                 XX_LLCORNER, XY_LLCORNER, XCELLSIZE, &
                 NCOLS, NROWS
!
!------------------------------------------------------------------------------
!
!*       1.    opening of namelist
! 
IF (LHOOK) CALL DR_HOOK('READ_NAM_GRID_IGN',0,ZHOOK_HANDLE)
 CALL GET_LUOUT(HPROGRAM,ILUOUT)
!
IF (HDIR/='H') THEN
  !
  CALL OPEN_NAMELIST(HPROGRAM,ILUNAM)
  !
  XX_LLCORNER = XUNDEF
  XY_LLCORNER = XUNDEF
  XCELLSIZE   = XUNDEF
  NCOLS = 0
  NROWS = 0
  !
  !---------------------------------------------------------------------------
  !
  !*       2.    Reading of projection parameters
  !              --------------------------------
  !
  CALL POSNAM(ILUNAM,'NAM_IGN',GFOUND,ILUOUT)
  IF (GFOUND) READ(UNIT=ILUNAM,NML=NAM_IGN)
  !
  !---------------------------------------------------------------------------
  CALL CLOSE_NAMELIST(HPROGRAM,ILUNAM)
  !---------------------------------------------------------------------------
  !
  !*       3.    Number of points
  !              ---------------------------------
  !
  IF ( XX_LLCORNER/=XUNDEF .AND. XY_LLCORNER/=XUNDEF &
              .AND. NCOLS>0 .AND. NROWS>0 ) THEN
    NPOINTS = NCOLS * NROWS
    ! else npoints must be defined in the namelist
  ENDIF
  !
  KL = NPOINTS
  !
  !---------------------------------------------------------------------------
  !
  !*       4.    Allocate Array of X and Y coordinates and increments
  !              --------------------------------------------------------
  !
  ALLOCATE(ZX(KL))
  ALLOCATE(ZY(KL))  
  ALLOCATE(ZDX(KL))
  ALLOCATE(ZDY(KL))
  !
  IF (XCELLSIZE/=XUNDEF) THEN
    !
    WRITE(ILUOUT,*) 'Initialisation of IGN Coordinates for a regular grid'
    !      
    ZDX(:) = XCELLSIZE
    ZDY(:) = XCELLSIZE
    !
    IF ( XX_LLCORNER/=XUNDEF .AND. XY_LLCORNER/=XUNDEF &
              .AND. NCOLS>0 .AND. NROWS>0 ) THEN
      ! M Lafaysse 08/2022
      ! There is no guarantee that NPOINTS < 1000000 (prescribed size of XX and YY)
      ! so it is necessary to initialize directly ZX and ZY which are allocated to the appropriate size.
      ! XX and XY must be used when the user defined himself all coordinates (e.g. SIM).
      DO JROWS=1,NROWS
        DO JCOLS=1,NCOLS
          !
          IINDEX = JCOLS + (JROWS-1) * NCOLS
          ZX(IINDEX) = XX_LLCORNER + (JCOLS-0.5) * XCELLSIZE
          ZY(IINDEX) = XY_LLCORNER + (JROWS-0.5) * XCELLSIZE
          !
        END DO
      END DO
      !
    ELSE
      !
      ! List of all x and y coordinates prescribed in namelist
      ZX(:) = XX(:KL)
      ZY(:) = XY(:KL)
      !
    ENDIF
  ELSE
    ! List of all x and y coordinates prescribed in namelist
    ZX(:) = XX(:KL)
    ZY(:) = XY(:KL)
    ! List of all dx and dy increments prescribed in namelist
    ZDX(:) = XDX(:KL)
    ZDY(:) = XDY(:KL)
  END IF
  !
  !---------------------------------------------------------------------------
  !
  !*       5.    Lambert type
  !              ------------
  !
  CALL TEST_NAM_VAR_SURF(ILUOUT,'CLAMBERT',CLAMBERT,'L1 ','L2 ','L3 ',&
                         'L4 ','L2E','L93' )  
  !
  SELECT CASE (CLAMBERT)
    CASE ('L1 ')
      ILAMBERT=1
    CASE ('L2 ')
      ILAMBERT=2
    CASE ('L3 ')
      ILAMBERT=3
    CASE ('L4 ')
      ILAMBERT=4
    CASE ('L2E')
      ILAMBERT=5
    CASE ('L93')
      ILAMBERT=6
  END SELECT
  !
  !---------------------------------------------------------------------------
  !
  !*       7.    maximum domain lengths
  !              ----------------------
  !
  ! If regular grid, squeeze get_xy_all_ign which is very slow for large grid
  IF ( XX_LLCORNER/=XUNDEF .AND. XY_LLCORNER/=XUNDEF .AND. NCOLS>0 .AND. NROWS>0 ) THEN
    ALLOCATE(ZXALL(NCOLS))
    ALLOCATE(ZYALL(NROWS))
    DO JCOLS=1,NCOLS
      ZXALL(JCOLS) = XX_LLCORNER + (JCOLS-0.5) * XCELLSIZE
    END DO
    DO JROWS=1,NROWS
      ZYALL(JROWS) = XY_LLCORNER + (JROWS-0.5) * XCELLSIZE
    END DO
    IDIMX = NCOLS
    IDIMY = NROWS
  ELSE
    ! usual case
    ALLOCATE(ZXALL(KL*5))
    ALLOCATE(ZYALL(KL*5))
    CALL GET_XYALL_IGN(ZX,ZY,ZDX,ZDY,ZXALL,ZYALL,IDIMX,IDIMY)
  END IF

  !
  !---------------------------------------------------------------------------
  !
  !*       8.    All this information stored into pointer PGRID_PAR
  !              --------------------------------------------------
  !
  CALL PUT_GRIDTYPE_IGN(ZGRID_PAR,ILAMBERT,ZX,ZY,ZDX,ZDY,        &
                      IDIMX,IDIMY,ZXALL(1:IDIMX),ZYALL(1:IDIMY))
  !
ELSE
  !
  ALLOCATE(ZX0(KDIM_FULL),ZY0(KDIM_FULL),ZDX0(KDIM_FULL),ZDY0(KDIM_FULL))
  !
  CALL GET_GRIDTYPE_IGN(PGRID_FULL_PAR,KLAMBERT=ILAMBERT,&
                        PX=ZX0,PY=ZY0,PDX=ZDX0,PDY=ZDY0)
  !
  KL = NSIZE_TASK(NRANK)
  ALLOCATE(ZX(KL),ZY(KL),ZDX(KL),ZDY(KL))
  ALLOCATE(ZXALL(KL*3),ZYALL(KL*3))
  IDIMX=0
  IDIMY=0
  !
  CALL READ_AND_SEND_MPI(ZX0,ZX)
  CALL READ_AND_SEND_MPI(ZY0,ZY)
  CALL READ_AND_SEND_MPI(ZDX0,ZDX)
  CALL READ_AND_SEND_MPI(ZDY0,ZDY)
  !
  DEALLOCATE(ZX0,ZY0,ZDX0,ZDY0)    
  !
  CALL PUT_GRIDTYPE_IGN(ZGRID_PAR,ILAMBERT,ZX,ZY,ZDX,ZDY,        &
                        IDIMX,IDIMY,ZXALL,ZYALL)
  !
ENDIF              
!
!---------------------------------------------------------------------------
DEALLOCATE(ZX)
DEALLOCATE(ZY)
DEALLOCATE(ZDX)
DEALLOCATE(ZDY)
DEALLOCATE(ZXALL)
DEALLOCATE(ZYALL)
!---------------------------------------------------------------------------
!
!* 1st call : initializes dimension
!
IF (KGRID_PAR==0) THEN
  KGRID_PAR = SIZE(ZGRID_PAR)
!
ELSE
!
!* 2nd call : initializes grid array
!
  PGRID_PAR(:) = 0.
  PGRID_PAR(:) = ZGRID_PAR
END IF
!
DEALLOCATE(ZGRID_PAR)
IF (LHOOK) CALL DR_HOOK('READ_NAM_GRID_IGN',1,ZHOOK_HANDLE)
!
!---------------------------------------------------------------------------
!
END SUBROUTINE READ_NAM_GRID_IGN
