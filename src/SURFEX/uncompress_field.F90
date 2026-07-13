!SFX_LIC Copyright 1994-2019 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!SFX_LIC for details. version 1.
!------------------------------------------------------------------------------------------------------------
! Modifications:
!  P. Wautelet 19/09/2019: correct support of 64bit integers (MNH_INT=8)
!  H. Petithomme 06/2023: optimization of do loop and vector assignments
!------------------------------------------------------------------------------------------------------------
SUBROUTINE UNCOMPRESS_FIELD(KLONG,PSEUIL,PFIELD_IN,PFIELD_OUT)

USE PARKIND1,ONLY: JPRB
USE YOMHOOK,ONLY: LHOOK,DR_HOOK

IMPLICIT NONE
 
INTEGER, INTENT(IN) :: KLONG
REAL, INTENT(IN) :: PSEUIL
REAL, DIMENSION(:), INTENT(IN) :: PFIELD_IN
REAL, DIMENSION(:), INTENT(OUT) :: PFIELD_OUT
INTEGER :: ICPT, I, N
REAL(KIND=JPRB) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('UNCOMPRESS_FIELD',0,ZHOOK_HANDLE)
ICPT = 0

PFIELD_OUT(:) = 0.

! boucle sur les colonnes
DO I=1,SIZE(PFIELD_IN,1)

  ! si la valeur est valide
  IF (PFIELD_IN(I)<PSEUIL) THEN

    ! on la met dans lwrite à l'indice icpt
    ICPT = ICPT + 1
    PFIELD_OUT(ICPT) = PFIELD_IN(I)

  ELSE

    !ideb = icpt + 1
    N = MIN(KLONG-ICPT,NINT(PFIELD_IN(I)-PSEUIL))
    PFIELD_OUT(ICPT+1:ICPT+N) = 0
    ICPT = ICPT+N

  ENDIF

  ! si on a dépassé la dernière colonne, on sort de la boucle
  IF (ICPT >= KLONG) EXIT

ENDDO

! test temporary: to remove after
WHERE (MOD(PFIELD_OUT(1:ICPT),100.) == 0) PFIELD_OUT(1:ICPT) = 0

IF (LHOOK) CALL DR_HOOK('UNCOMPRESS_FIELD',1,ZHOOK_HANDLE)
END SUBROUTINE UNCOMPRESS_FIELD
