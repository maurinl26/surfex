!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!     #########
      SUBROUTINE AVERAGE1_OROGRAPHY (UG,USS, &
                                     KLUOUT,KNBLINES,PLAT,PLON,PVALUE,PNODATA)
!     #######################################################
!
!!**** *AVERAGE1_OROGRAPHY* computes the sum of orography, squared orography
!!                              and subgrid orography characteristics
!!
!!    PURPOSE
!!    -------
!!
!!    METHOD
!!    ------
!!   
!!    EXTERNAL
!!    --------
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!
!!    REFERENCE
!!    ---------
!!
!!    AUTHOR
!!    ------
!!
!!    V. Masson         Meteo-France
!!
!!    MODIFICATION
!!    ------------
!!
!!    Original    12/09/95
!!
!----------------------------------------------------------------------------
!
!*    0.     DECLARATION
!            -----------
!
!
USE MODD_SURF_ATM_GRID_n, ONLY : SURF_ATM_GRID_t
USE MODD_SSO_n, ONLY : SSO_t
!
USE MODD_PGDWORK,       ONLY : XALL, XEXT_ALL, NSIZE_ALL, XSSO_ALL, NSSO_ALL, NSSO, &
                               NFSSO_ALL, XFSSO_ALL, NFSSO, LORORAD
!
USE MODI_GET_MESH_INDEX
USE MODD_POINT_OVERLAY, ONLY : NOVMX
!
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*    0.1    Declaration of arguments
!            ------------------------
!
!
TYPE(SURF_ATM_GRID_t), INTENT(INOUT) :: UG
TYPE(SSO_t), INTENT(INOUT) :: USS
!
INTEGER,                 INTENT(IN)    :: KLUOUT
INTEGER,                 INTENT(IN)    :: KNBLINES
REAL, DIMENSION(:),      INTENT(IN)    :: PLAT    ! latitude of the point to add
REAL, DIMENSION(:),      INTENT(IN)    :: PLON    ! longitude of the point to add
REAL, DIMENSION(:),      INTENT(IN)    :: PVALUE  ! value of the point to add
REAL, OPTIONAL, INTENT(IN) :: PNODATA
!
!*    0.2    Declaration of other local variables
!            ------------------------------------
!
INTEGER, DIMENSION(NOVMX,SIZE(PLAT)) :: IINDEX  ! mesh index of all input points
                                         ! 0 indicates the point is out of the domain
INTEGER, DIMENSION(NOVMX,SIZE(PLAT)) :: ISSOX   ! X submesh index in their mesh of all input points
INTEGER, DIMENSION(NOVMX,SIZE(PLAT)) :: ISSOY   ! Y submesh index in their mesh of all input points
INTEGER, DIMENSION(NOVMX,SIZE(PLAT)) :: IFSSOX  ! X fractional-submesh index in their mesh of all input points
INTEGER, DIMENSION(NOVMX,SIZE(PLAT)) :: IFSSOY  ! Y fractional-submesh index in their mesh of all input points
!
INTEGER :: JL, JOV,IP        ! loop index on input arrays
REAL(KIND=JPRB) :: ZHOOK_HANDLE,ZH
!----------------------------------------------------------------------------
!
!
!*    1.     Get position
!            ------------
!     
IF (LHOOK) CALL DR_HOOK('AVERAGE1_OROGRAPHY',0,ZHOOK_HANDLE)
!      
IF (PRESENT(PNODATA)) THEN
  IF (LORORAD) THEN
    CALL GET_MESH_INDEX(UG,KLUOUT,KNBLINES,PLAT,PLON,IINDEX,PVALUE,PNODATA,NSSO,ISSOX,ISSOY, &
                        NFSSO,IFSSOX,IFSSOY)
  ELSE
    CALL GET_MESH_INDEX(UG,KLUOUT,KNBLINES,PLAT,PLON,IINDEX,PVALUE,PNODATA,NSSO,ISSOX,ISSOY)
  ENDIF
ELSE
  IF (LORORAD) THEN
    CALL GET_MESH_INDEX(UG,KLUOUT,KNBLINES,PLAT,PLON,IINDEX,KSSO=NSSO,KISSOX=ISSOX,KISSOY=ISSOY, & 
                        KFSSO=NFSSO,KFISSOX=IFSSOX,KFISSOY=IFSSOY)
  ELSE
    CALL GET_MESH_INDEX(UG,KLUOUT,KNBLINES,PLAT,PLON,IINDEX,KSSO=NSSO,KISSOX=ISSOX,KISSOY=ISSOY)
  ENDIF
ENDIF
!
!*    2.     Loop on all input data points
!            -----------------------------

DO JOV = 1, NOVMX
  DO JL = 1 , SIZE(PLAT)
!
!*    3.     Tests on position
!            -----------------
!     
    IF (IINDEX(JOV,JL)==0) CYCLE
!
!*    4.     Summation
!            ---------
!
    IP = IINDEX(JOV,JL)
    NSIZE_ALL(IP,1) = NSIZE_ALL(IP,1)+1

!*    5.     Orography
!            ---------
!
    XALL(IP,1,1) = XALL(IP,1,1)+PVALUE(JL)
!
!*    6.     Square of Orography
!            -------------------
!
    XALL(IP,2,1) = XALL(IP,2,1)+PVALUE(JL)**2
!
!*    7.     Maximum orography in a subgrid square
!            -------------------------------------
!
    NSSO_ALL(IP,ISSOX(JOV,JL),ISSOY(JOV,JL)) = 1
    XSSO_ALL(IP,ISSOX(JOV,JL),ISSOY(JOV,JL)) = &
         MAX (  XSSO_ALL(IP,ISSOX(JOV,JL),ISSOY(JOV,JL)) , PVALUE(JL) )   
!
!
!*    8.     Maximum orography in the mesh
!            -----------------------------
!
    XEXT_ALL(IP,1) = MAX(XEXT_ALL(IP,1),PVALUE(JL))
!
!
!*    9.     Minimum orography in the mesh
!            -----------------------------
!
    XEXT_ALL(IP,2) = MIN(XEXT_ALL(IP,2),PVALUE(JL))
  END DO
ENDDO
!
!
!*    10.    Orographic radiation parameters
!            -------------------------------
!
IF (LORORAD) THEN
  IF (LHOOK) CALL DR_HOOK('AVERAGE1_OROGRAPHY_FSSO',0,ZH)
  DO JOV = 1, NOVMX
    DO JL = 1 , SIZE(PLAT)
      IF (IINDEX(JOV,JL)==0) CYCLE

      XFSSO_ALL(IINDEX(JOV,JL),IFSSOX(JOV,JL),IFSSOY(JOV,JL)) = &
         XFSSO_ALL(IINDEX(JOV,JL),IFSSOX(JOV,JL),IFSSOY(JOV,JL)) + PVALUE(JL)
      NFSSO_ALL(IINDEX(JOV,JL),IFSSOX(JOV,JL),IFSSOY(JOV,JL)) = &
         NFSSO_ALL(IINDEX(JOV,JL),IFSSOX(JOV,JL),IFSSOY(JOV,JL)) + 1
    END DO
  ENDDO
  IF (LHOOK) CALL DR_HOOK('AVERAGE1_OROGRAPHY_FSSO',1,ZH)
ENDIF
IF (LHOOK) CALL DR_HOOK('AVERAGE1_OROGRAPHY',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE AVERAGE1_OROGRAPHY
