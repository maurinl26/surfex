!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!     #########
      SUBROUTINE GET_NEAR_MESHES_IGN(KGRID_PAR,KL,PGRID_PAR,KNEAR_NBR,KNEAR)
!     ##############################################################
!
!!**** *GET_NEAR_MESHES_IGN* get the near grid mesh indices
!!
!!    PURPOSE
!!    -------
!!
!!    METHOD
!!    ------
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
!!    Original    03/2004
!!
!
!     Modifié par Renaud Lestrigant (02/2016) : changement complet de l'algo
!     de recherche des plus proches voisins.
!     Récupération d'un code sur Internet et adaptation locale.
!     (http://jblevins.org/mirror/amiller/qsort.f90)
!     
!     Modifié par B. Cluzet (20/19/19): déplacement des sous-routines dans mode_quick_sort.F90 (créé) pour mutualisation
!
!----------------------------------------------------------------------------
!
!*    0.     DECLARATION
!            -----------
!
USE MODD_SURFEX_MPI, ONLY : NINDEX, NRANK, NNUM
USE MODE_GRIDTYPE_IGN
USE MODE_QUICK_SORT, ONLY : QUICK_SORT 
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*    0.1    Declaration of arguments
!            ------------------------
!
INTEGER,                         INTENT(IN)    :: KGRID_PAR ! size of PGRID_PAR
INTEGER,                         INTENT(IN)    :: KL        ! number of points
INTEGER,                         INTENT(IN)    :: KNEAR_NBR ! number of nearest points wanted
REAL,    DIMENSION(KGRID_PAR),   INTENT(IN)    :: PGRID_PAR ! grid parameters
INTEGER, DIMENSION(:,:),POINTER :: KNEAR    ! near mesh indices
!
!*    0.2    Declaration of other local variables
!            ------------------------------------
!
REAL,DIMENSION(KL)  :: ZX
REAL,DIMENSION(KL)  :: ZY
REAL,DIMENSION(KL)  :: ZDX
REAL,DIMENSION(KL)  :: ZDY
REAL,DIMENSION(KL)  :: ZDIS
INTEGER,DIMENSION(KL) :: INDZDIS

REAL :: ZMAXVALDIS

INTEGER :: JP, ID, ISIZE
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!----------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('GET_NEAR_MESHES_IGN_1',0,ZHOOK_HANDLE)
!
 CALL GET_GRIDTYPE_IGN(PGRID_PAR,PX=ZX,PY=ZY,PDX=ZDX,PDY=ZDY)
!
KNEAR(:,:) = 0
!
ISIZE = MIN(KNEAR_NBR,KL)
!
! calcul de la distance de tous les points 2 à 2
!
!
IF (LHOOK) CALL DR_HOOK('GET_NEAR_MESHES_IGN_1',1,ZHOOK_HANDLE)
IF (LHOOK) CALL DR_HOOK('GET_NEAR_MESHES_IGN_2',0,ZHOOK_HANDLE)
!
!$OMP PARALLEL DO PRIVATE(JP,ZDIS,ZMAXVALDIS,INDZDIS,ID)
!
DO JP=1,KL
  !
  IF (NINDEX(JP)==NRANK) THEN 
    !
    ID = NNUM(JP)
    !
    ! distance du point JP à tous les autres points
    ZDIS(:) = SQRT((ZX(:)-ZX(JP))**2 + (ZY(:)-ZY(JP))**2)
    ! distance maximale entre JP et les autres moints
    ZMAXVALDIS = 2. * MAXVAL(ZDIS)
    ZDIS(JP) = ZMAXVALDIS

    CALL QUICK_SORT(ZDIS, INDZDIS)
    KNEAR(ID,1:ISIZE) = INDZDIS(1:ISIZE)
    !
  ENDIF
  !
ENDDO
!$OMP END PARALLEL DO
!
IF (LHOOK) CALL DR_HOOK('GET_NEAR_MESHES_IGN_2',1,ZHOOK_HANDLE)
END SUBROUTINE GET_NEAR_MESHES_IGN
