!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!-------------------------------------------------------------------------------
!     ####################
      SUBROUTINE GET_UPSLOPE(PCONN,KFILE,KNMC,KOUTLET,PDX,KNXC,KNYC,PTOPD,KLIST_PIX,KB_SUBPIX)
!     ####################
!
!!****  *GET_UPSLOPE*  
!!    ------------------ 
!!    PURPOSE
!!    ---------
!     This routine aims at findind a sub-catchment mask and surface
!!    
!!    AUTHOR 
!!    ------
!!
!!      B. Vincendon    * Meteo-France *
!!
!!    MODIFICATIONS
!!    -------------
!!
!!      Original   03/2013
!-------------------------------------------------------------------------------
!
!*       0.     DECLARATIONS
!               ------------
!
USE MODD_SURF_PAR  ,ONLY : XUNDEF, NUNDEF
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
REAL   , DIMENSION(:,:), INTENT(IN) :: PCONN
INTEGER, DIMENSION(:)  , INTENT(IN) :: KFILE
INTEGER                , INTENT(IN) :: KNMC
INTEGER                , INTENT(IN) :: KOUTLET
REAL                   , INTENT(IN) :: PDX
INTEGER                , INTENT(IN) :: KNXC
INTEGER                , INTENT(IN) :: KNYC
REAL   , DIMENSION(:)  , INTENT(IN) :: PTOPD
INTEGER, DIMENSION(:)  , INTENT(OUT):: KLIST_PIX ! mask of the sub-catchment
INTEGER                , INTENT(OUT):: KB_SUBPIX ! number of pixels in the sub-catchment
!
!*      0.2    declarations of local variables
!
INTEGER                            :: JI,JN,JITER
INTEGER                            :: ILINE,ILINE_OLD
INTEGER                            :: IREF
LOGICAL                            :: GFOUND,GTAGGED
INTEGER, DIMENSION(:), ALLOCATABLE :: INEW, IOLD
INTEGER                            :: INB_NEW, INB_OLD
REAL                               :: ZSURFACE, ZDX2
REAL, DIMENSION(:,:), ALLOCATABLE  :: ZHISTORY
INTEGER                            :: I_X0,I_Y0,I_Z0
INTEGER                            :: I_X,I_Y,I_Z

!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('GET_UPSLOPE',0,ZHOOK_HANDLE)
!
ALLOCATE(INEW(KNMC))
ALLOCATE(IOLD(KNMC))
ALLOCATE(ZHISTORY(KNMC,2))
!
!       1. Initiliaze
!
ZDX2=PDX*PDX
!
DO JI=1,KNMC
  KLIST_PIX(JI)=0
  ZHISTORY(JI,1)=PCONN(JI,1)
  ZHISTORY(JI,2)=0
ENDDO
!
!       2.  Looking for upslopes pixels
!
GFOUND  = .TRUE.
INB_NEW = 0
INB_OLD = 1
KB_SUBPIX = 1
IOLD(1) = KOUTLET
ILINE = KFILE(KOUTLET)
ZHISTORY(JI,2) = 1
KLIST_PIX(KB_SUBPIX) = KOUTLET
JITER = 0
!
DO WHILE (GFOUND)
  ! Update N° Iteration
  JITER=JITER+1
  ! Update Screen Text Control
  ZSURFACE=REAL(KB_SUBPIX)*ZDX2
  ZSURFACE=ZSURFACE/1000000.
  !
  GFOUND=.FALSE.
  ! Explore 'previous New Pixels' (stored in Old Vector) for upslope neighboorhing
  INB_NEW=0
  !
  DO JI=1,INB_OLD
    ILINE=KFILE(IOLD(JI))
    I_Z0=IOLD(JI)
    I_Y0=(I_Z0/KNXC)+1
    I_X0=I_Z0-(I_Y0-1)*KNXC
    !
    DO I_X=I_X0-1,I_X0+1
      !
      DO I_Y=I_Y0-1,I_Y0+1
        !
        IF(I_X>0.AND.I_Y>0.AND.I_X<=KNXC.AND.I_Y<=KNYC) THEN
          !
          I_Z=I_X+(I_Y-1)*KNXC
          IF (I_Z/=I_Z0.AND.PTOPD(I_Z)/=XUNDEF.AND.PTOPD(I_Z)>=PTOPD(I_Z0))THEN
            !
            IREF=I_Z
            IF (IREF/=0) THEN
              ILINE_OLD=KFILE(IREF)
              IF (ZHISTORY(ILINE_OLD,2)==0.) THEN
                !
                GTAGGED=.FALSE.
                DO JN=1,INB_NEW
                  !
                  IF (INEW(JN)==IREF) THEN
                    GTAGGED=.TRUE.
                  ENDIF
                  !
                ENDDO
                IF (.NOT.GTAGGED) THEN
                  GFOUND=.TRUE.
                  INB_NEW=INB_NEW+1
                  INEW(INB_NEW)=IREF
                ENDIF
                !
              ENDIF
            ENDIF
          ENDIF
        ENDIF
      ENDDO
    ENDDO
  ENDDO  
!!
! Check if these new pixels were not already memorized
! Memorize them if not.
  INB_OLD=0
  DO JI=1,INB_NEW
    IREF=INEW(JI)
    ILINE=KFILE(IREF)
    !
    IF (ZHISTORY(ILINE,2)==0.) THEN
      !
      INB_OLD=INB_OLD+1
      IOLD(INB_OLD)=IREF     
      KB_SUBPIX=MIN(KB_SUBPIX+1,KNMC)
      KLIST_PIX(KB_SUBPIX)=IREF
      ZHISTORY(ILINE,2)=1
      !
    ENDIF
    !
  ENDDO
  !
ENDDO
!
!
IF (LHOOK) CALL DR_HOOK('GET_UPSLOPE',1,ZHOOK_HANDLE)
!
END SUBROUTINE GET_UPSLOPE

