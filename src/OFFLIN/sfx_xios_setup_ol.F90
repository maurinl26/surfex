SUBROUTINE SFX_XIOS_SETUP_OL(YSC, KLUOUT, KYEAR, KMONTH, KDAY, PTIME, PSTEP, PSW_BANDS)
!
!**** *SFX_XIOS_SETUP_OL*  - 
!
!     Purpose.
!     --------
!
!       Call SFX_XIOS_SETUP, providing it with Arpege MPI communicator
!       and passing args about : output logical unit, date, model
!       timestep and, for the whole MPI task : lat/lon of centers and
!       corners, cell index and mask, tile masks
!
!**   Interface.
!     ----------
!       *CALL*  *SFX_XIOS_SETUP_OL*
!
!     Input:
!     -----
!
!     Output:
!     ------
!
!
!     Method:
!     ------
!
!     Externals:
!     ---------
!
!     Reference:
!     ---------
!
!     Author:
!     -------
!      S.Senesi, aug 2015
!
!     Modifications.
!     --------------
!      S.Senesi, aug 2016
!
!     -----------------------------------------------------------
!
USE MODN_IO_OFFLINE, ONLY  : LALLOW_ADD_DIM, LGRID_MODE
!
USE MODD_SURFEX_n,    ONLY : SURFEX_t
USE MODD_XIOS,  ONLY       : LXIOS, LADD_DIM=>LALLOW_ADD_DIM, NBLOCK_XIOS=>NBLOCK, LGRID=>LGRID_MODE
USE MODD_SURFEX_MPI, ONLY  : NRANK, NPROC, NCOMM, NINDEX
USE MODD_IO_SURF_OL, ONLY  : NMASK_IGN
!
USE MODI_GET_SURF_GRID_DIM_n
USE MODI_GET_MESH_CORNER
USE MODI_SFX_XIOS_SETUP
USE MODI_SFX_XIOS_GAUSS_GRID
USE MODI_ABOR1_SFX
!
USE MODE_GRIDTYPE_GAUSS
USE MODI_GET_IGN_MASKALL
!
USE PARKIND1  ,ONLY : JPRB, JPIM
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
!
IMPLICIT NONE
!
TYPE(SURFEX_t),  INTENT(INOUT) :: YSC
INTEGER,            INTENT(IN) :: KLUOUT
INTEGER(KIND=JPIM), INTENT(IN) :: KYEAR      ! Current Year
INTEGER(KIND=JPIM), INTENT(IN) :: KMONTH     ! Current Month
INTEGER(KIND=JPIM), INTENT(IN) :: KDAY       ! Current Day 
REAL              , INTENT(IN) :: PTIME      ! Time in the day
REAL              , INTENT(IN) :: PSTEP      ! Atmospheric time step
REAL,DIMENSION(:),  INTENT(IN) :: PSW_BANDS
!
! Local variables
!
INTEGER                              :: IX, JX
INTEGER, ALLOCATABLE, DIMENSION(:)   :: IXINDEX ! Index of the grid meshes for the 
                                                ! current MPI-task  in global 1D grid (start at 0)
LOGICAL, ALLOCATABLE, DIMENSION(:)   :: LLXMASK ! Cells mask
REAL,    ALLOCATABLE, DIMENSION(:,:) :: ZCORLON
REAL,    ALLOCATABLE, DIMENSION(:,:) :: ZCORLAT
INTEGER, ALLOCATABLE, DIMENSION(:)   :: IMN,IMS,IMW,IMT ! Tile masks, 0-based
REAL, DIMENSION(:), ALLOCATABLE      :: ZX, ZY
INTEGER                              :: NLATI
INTEGER, ALLOCATABLE, DIMENSION(:)   :: NLOPA
!
CHARACTER(LEN=10)   :: YGRID      ! grid type 
INTEGER             :: IDIM1, IDIM2
LOGICAL             :: LLRECT, LLGAUSS
!
REAL(KIND=JPRB)     :: ZHOOK_HANDLE
!
!------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('SFX_XIOS_SETUP_OL',0,ZHOOK_HANDLE)
!
#ifdef WXIOS
IF (LXIOS) THEN
   !
   LADD_DIM=LALLOW_ADD_DIM
   LGRID=LGRID_MODE
   NBLOCK_XIOS=1
   !
   !Build XIOS index from Surfex's NINDEX
   ALLOCATE(IXINDEX(YSC%U%NSIZE_FULL))
   JX=0
   DO IX=1,YSC%U%NDIM_FULL
      IF(NINDEX(IX)==NRANK)THEN
         JX=JX+1
         IF(JX>YSC%U%NSIZE_FULL)THEN
            CALL ABOR1_SFX('sfx_xios_setup_ol : internal error with XIOS index')
         ENDIF
         IXINDEX(JX) = IX-1
      ENDIF
   ENDDO
   !
   ALLOCATE(LLXMASK(YSC%U%NSIZE_FULL))
   LLXMASK(:)=.TRUE.
   !
   ALLOCATE(IMN(SIZE(YSC%U%NR_NATURE))); IMN(:)=YSC%U%NR_NATURE-1
   ALLOCATE(IMS(SIZE(YSC%U%NR_SEA)))   ; IMS(:)=YSC%U%NR_SEA   -1
   ALLOCATE(IMW(SIZE(YSC%U%NR_WATER ))); IMW(:)=YSC%U%NR_WATER -1
   ALLOCATE(IMT(SIZE(YSC%U%NR_TOWN  ))); IMT(:)=YSC%U%NR_TOWN  -1
   !
   ALLOCATE(ZCORLAT(YSC%U%NSIZE_FULL,4))
   ALLOCATE(ZCORLON(YSC%U%NSIZE_FULL,4))
   CALL GET_MESH_CORNER(YSC%UG, KLUOUT,ZCORLAT(:,:),ZCORLON(:,:))
   !
   CALL GET_SURF_GRID_DIM_N(YSC%UG, YGRID, LLRECT, IDIM1, IDIM2)
   !
   LLGAUSS=.FALSE.
   IF (YSC%UG%G%CGRID == "GAUSS     ") THEN
     LLGAUSS=.TRUE.
     CALL GET_GRIDTYPE_GAUSS(YSC%UG%G%XGRID_PAR,KNLATI=NLATI)
     ALLOCATE(NLOPA(NLATI))
     CALL GET_GRIDTYPE_GAUSS(YSC%UG%G%XGRID_PAR,KNLOPA=NLOPA)
     IDIM1=SUM(NLOPA)
     IDIM2=1
   ELSEIF (TRIM(YSC%UG%G%CGRID)=='IGN') THEN
     ALLOCATE(ZX(IDIM1),ZY(IDIM2))
     CALL GET_IGN_MASKALL(YSC%UG,YSC%U%NDIM_FULL,ZX,ZY)
     IXINDEX(:) = NMASK_IGN(IXINDEX(:)+1) -1
     DEALLOCATE(NMASK_IGN,ZX,ZY)
   ENDIF
   !
   CALL SFX_XIOS_SETUP(YSC%IM,YSC%UG,NCOMM,KLUOUT,KYEAR,KMONTH,KDAY,PTIME,PSTEP,     &
                       IDIM1,IDIM2,0,PSW_BANDS,TRANSPOSE(ZCORLAT),TRANSPOSE(ZCORLON),&
                       IXINDEX,LLXMASK,IMN,IMS,IMW,IMT)
   !
   DEALLOCATE(IMN,IMS,IMW,IMT)
   DEALLOCATE(ZCORLAT,ZCORLON)
   DEALLOCATE(IXINDEX)
   DEALLOCATE(LLXMASK)
   IF(LLGAUSS)THEN
     CALL SFX_XIOS_GAUSS_GRID(NLATI,MAXVAL(NLOPA),NRANK,NPROC)
     DEALLOCATE(NLOPA)
   ENDIF
   !
ENDIF
#endif
!
IF (LHOOK) CALL DR_HOOK('SFX_XIOS_SETUP_OL',1,ZHOOK_HANDLE)
!
END SUBROUTINE SFX_XIOS_SETUP_OL
