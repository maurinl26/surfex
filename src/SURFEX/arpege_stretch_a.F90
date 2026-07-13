!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!     ##########################################################################
SUBROUTINE ARPEGE_STRETCH_A(KN,PLAP,PLOP,PCOEF,PLAR,PLOR,PLAC,PLOC)
!     ##########################################################################
!!****  *ARPEGE_STRETCH_A* - Projection to Arpege stretched grid
!!
!!   PURPOSE
!!   -------
!!
!!   Projection from standard Lat,Lon grid to Arpege stretched grid
!!
!!   METHOD
!!   ------
!!
!!   The projection is defined in two steps :
!!    1. A rotation to place the stretching pole at the north pole
!!    2. The stretching
!!   This routine is a basic implementation of the informations founded in 
!!     'Note de travail Arpege n#3'
!!     'Transformation de coordonnees'
!!     J.F.Geleyn 1988
!!   This document describes a slightly different transformation in 3 steps. Only the
!!   two first steps are to be taken in account (at the time of writing this paper has
!!   not been updated).
!!
!!   EXTERNAL
!!   --------
!!
!!
!!   IMPLICIT ARGUMENTS
!!   ------------------
!!
!!   REFERENCE
!!   ---------
!!
!!   This routine is based on : 
!!     'Note de travail ARPEGE' number 3
!!     by J.F. GELEYN (may 1988)
!!
!!   AUTHOR
!!   ------
!!
!!   V.Bousquet
!!
!!   MODIFICATIONS
!!   -------------
!!
!!   Original       07/01/1999
!!   V. Masson      01/2004    Externalization of surface
!!
!
! 0. DECLARATIONS
! ---------------
!
  USE MODD_CSTS,ONLY : XPI
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
  IMPLICIT NONE
!
! 0.1. Declaration of arguments
! -----------------------------

  INTEGER,             INTENT(IN)  :: KN            ! Number of points to convert
  REAL,                INTENT(IN)  :: PLAP          ! Latitude of stretching pole
  REAL,                INTENT(IN)  :: PLOP          ! Longitude of stretching pole
  REAL,                INTENT(IN)  :: PCOEF         ! Stretching coefficient
  REAL, DIMENSION(KN), INTENT(IN)  :: PLAR          ! Lat. of points
  REAL, DIMENSION(KN), INTENT(IN)  :: PLOR          ! Lon. of points
  REAL, DIMENSION(KN), INTENT(OUT) :: PLAC          ! Computed pseudo-lat. of points
  REAL, DIMENSION(KN), INTENT(OUT) :: PLOC          ! Computed pseudo-lon. of points
!
  REAL                             :: ZSINSTRETCHLA ! Sine of stretching point lat.
  REAL                             :: ZSINSTRETCHLO ! Sine of stretching point lon.
  REAL                             :: ZCOSSTRETCHLA ! Cosine of stretching point lat.
  REAL                             :: ZCOSSTRETCHLO ! Cosine of stretching point lon.
  REAL                             :: ZSINLA        ! Sine of computed point latitude
  REAL                             :: ZSINLO        ! Sine of computed point longitude
  REAL                             :: ZCOSLA        ! Cosine of computed point latitude
  REAL                             :: ZCOSLO        ! Cosine of computed point longitude
  REAL                             :: ZSINLAS       ! Sine of point's pseudo-latitude
  REAL                             :: ZSINLOS       ! Sine of point's pseudo-longitude
  REAL                             :: ZCOSLOS       ! Cosine of point's pseudo-lon.
  REAL                             :: ZA,ZB,ZD      ! Dummy variables used for 
  REAL                             :: ZX,ZY,ZMU ! computations
!  
  INTEGER                          :: JP        ! Dummy loop counter
  REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
  IF (LHOOK) CALL DR_HOOK('ARPEGE_STRETCH_A',0,ZHOOK_HANDLE)

  ZSINSTRETCHLA = SIN(PLAP*XPI/180.)
  ZCOSSTRETCHLA = COS(PLAP*XPI/180.)
  ZSINSTRETCHLO = SIN(PLOP*XPI/180.)
  ZCOSSTRETCHLO = COS(PLOP*XPI/180.)
  ! L = longitude (0 = Greenwich, + toward east)
  ! l = latitude (90 = N.P., -90 = S.P.)
  ! p stands for stretching pole
  PLAC(:) = PLAR(:) * XPI / 180.
  PLOC(:) = PLOR(:) * XPI / 180.
  ! A = 1 + c.c
  ZA = 1. + PCOEF**2
  ! B = 1 - c.c
  ZB = 1. - PCOEF**2

  DO JP=1, KN
    ZSINLA = SIN(PLAC(JP))
    ZCOSLA = COS(PLAC(JP))
    ZSINLO = SIN(PLOC(JP))
    ZCOSLO = COS(PLOC(JP))
    ! X = cos(Lp-L)
    ZX = ZCOSLO*ZCOSSTRETCHLO + ZSINLO*ZSINSTRETCHLO
    ! Y = sin(Lp-L)
    ZY = ZCOSLO*ZSINSTRETCHLO - ZSINLO*ZCOSSTRETCHLO

    ! D = (1+c.c) + (1-c.c)(sin lp.sin l + cos lp.cos l.cos(Lp-L))
    ZMU = ZSINSTRETCHLA*ZSINLA+ZCOSSTRETCHLA*ZCOSLA*ZX
    ZD = ZA + ZB*ZMU
    !          (1-c.c)+(1+c.c)((sin lp.sin l + cos lp.cos l.cos(Lp-L))
    ! sin lr = -------------------------------------------------------
    !                                  D
    ZSINLAS = (ZB + ZA*ZMU) / ZD
    ! D' = D * cos lr
    ZD = ZD * AMAX1(1e-6,SQRT(1.-ZSINLAS**2))
    !          2.c.(cos lp.sin l - sin lp.cos l.cos(Lp-L))
    ! cos Lr = -------------------------------------------
    !                              D'
    ZCOSLOS = 2.*PCOEF*(ZCOSSTRETCHLA*ZSINLA-ZSINSTRETCHLA*ZCOSLA*ZX) / ZD
    !          2.c.cos l.cos(Lp-L)
    ! sin Lr = -------------------
    !                  D'
    ZSINLOS = 2.*PCOEF*(ZCOSLA*ZY) / ZD

    ! saturations (corrects calculation errors)
    ZSINLAS = MAX(ZSINLAS,-1.)
    ZSINLAS = MIN(ZSINLAS, 1.)
    ZCOSLOS = MAX(ZCOSLOS,-1.)
    ZCOSLOS = MIN(ZCOSLOS, 1.)

    ! back from sine & cosine
    PLAC(JP) = ASIN(ZSINLAS)*180./XPI
    PLOC(JP) = SIGN(ACOS(ZCOSLOS),ZSINLOS)*180./XPI
  ENDDO

  IF (LHOOK) CALL DR_HOOK('ARPEGE_STRETCH_A',1,ZHOOK_HANDLE)
END SUBROUTINE ARPEGE_STRETCH_A
!-------------------------------------------------------------------------------
