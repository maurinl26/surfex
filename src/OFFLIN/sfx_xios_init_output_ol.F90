!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
SUBROUTINE SFX_XIOS_INIT_OUTPUT_OL(YSC)
!
!!****  *SFX_XIOS_INIT_OUTPUT_OL* writes coordonates values
!!
!!    PURPOSE
!!    -------
!       Write if asked the values of the coordonates (for IGN grid) with Lat, Lon values also.
!       It also writes the altitude for all the points (ZS).
!!
!!**  IMPLICIT ARGUMENTS
!!    ------------------
!!      None
!!
!!    REFERENCE
!!    ---------
!!
!!    AUTHOR
!!    ------
!!      Mathieu Fructus   *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!
!!=================================================================
!
!*       0.   DECLARATIONS
!             ------------
!
USE MODD_SURFEX_n, ONLY : SURFEX_t
USE MODE_GRIDTYPE_IGN
USE MODI_WRITE_SURF
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
TYPE(SURFEX_t), INTENT(IN) :: YSC
!
! Local variables
REAL, DIMENSION(:), ALLOCATABLE   :: ZLAT, ZLON       ! Latitude and Longitude
REAL, DIMENSION(:), ALLOCATABLE   :: ZX               ! X conformal coordinate of grid mesh (dim IIMAX)
REAL, DIMENSION(:), ALLOCATABLE   :: ZY               ! Y conformal coordinate of grid mesh (dim IJMAX)
REAL(KIND=JPRB)                   :: ZHOOK_HANDLE     ! in order to profile with Dr Hook
CHARACTER(LEN=100)                :: ZCOMMENT         ! comment linked to the netcdf output field
INTEGER                           :: ZRESP            ! error return code from write_surf
INTEGER                           :: IL, ZDIMX, ZDIMY ! index for number of points, dimension in x and y direction
!
IF (LHOOK) CALL DR_HOOK('SFX_XIOS_INIT_OUTPUT_OL',0,ZHOOK_HANDLE)
!
IL = NINT(YSC%UG%G%XGRID_PAR(2))
ZDIMX = NINT(YSC%UG%G%XGRID_PAR(3+4*IL))
ZDIMY = NINT(YSC%UG%G%XGRID_PAR(4+4*IL))
ALLOCATE(ZX(IL))
ALLOCATE(ZY(IL))
ALLOCATE(ZLAT(IL))
ALLOCATE(ZLON(IL))
CALL GET_GRIDTYPE_IGN(YSC%UG%G%XGRID_PAR, PX=ZX, PY=ZY)
CALL LATLON_IGN(6,ZX,ZY,PLAT=ZLAT,PLON=ZLON)
!
ZCOMMENT='elevation'
CALL WRITE_SURF((/'ZS'/),'XIOS  ','ZS',YSC%U%XZS,ZRESP,ZCOMMENT)
!
ZCOMMENT='XX Lambert93'
CALL WRITE_SURF((/'XX'/),'XIOS  ','XX',ZX,ZRESP,ZCOMMENT)
!
ZCOMMENT='YY Lambert93'
CALL WRITE_SURF((/'YY'/),'XIOS  ','YY',ZY,ZRESP,ZCOMMENT)
!
ZCOMMENT='latitude'
CALL WRITE_SURF((/'LAT'/),'XIOS  ','LAT',ZLAT,ZRESP,ZCOMMENT)
!
ZCOMMENT='longitude'
CALL WRITE_SURF((/'LON'/),'XIOS  ','LON',ZLON,ZRESP,ZCOMMENT)
!
DEALLOCATE(ZLAT,ZLON,ZX,ZY)
!
IF (LHOOK) CALL DR_HOOK('SFX_XIOS_INIT_OUTPUT_OL',1,ZHOOK_HANDLE)
!
END SUBROUTINE SFX_XIOS_INIT_OUTPUT_OL
