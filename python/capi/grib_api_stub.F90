!SFX_KARPOS
! grib_api_stub.F90 — Stub GRIB_API pour build OFFLINE/SODA NetCDF-only.
!
! Cette distribution SURFEX suppose grib_api toujours disponible (modd_grid_grib
! et mode_read_grib font `USE GRIB_API`). Les forcings Karpos sont en NetCDF :
! les chemins de lecture GRIB ne sont JAMAIS exercés. On fournit donc un stub
! qui débloque la compilation sans dépendance eccodes.
!
! Stratégie : le module n'exporte que les données (kindOfInt, GRIB_END_OF_INDEX).
! Les procédures grib_* sont laissées EXTERNES (interface implicite) : les appels
! à signatures variées de mode_read_grib compilent grâce à -fallow-argument-mismatch,
! et les corps stub ci-dessous s'arrêtent net si jamais ils sont atteints.
!
MODULE GRIB_API
  IMPLICIT NONE
  PUBLIC
  ! kindOfInt : grib_api le dérive du C int ; le kind entier par défaut convient
  ! (SURFEX déclare INTEGER(KIND=kindOfInt) et passe des entiers par défaut).
  INTEGER, PARAMETER :: kindOfInt = KIND(0)
  ! Sentinelles de fin d'index / fichier (valeurs grib_api ; jamais atteintes ici).
  INTEGER(kindOfInt), PARAMETER :: GRIB_END_OF_INDEX = -43
  INTEGER(kindOfInt), PARAMETER :: GRIB_END_OF_FILE  = -1
END MODULE GRIB_API

! --- Corps stub des procédures grib_api (interface implicite) -----------------
! Appelées uniquement sur un chemin GRIB, absent des runs OFFLINE/SODA NetCDF.

SUBROUTINE GRIB_ABORT_KARPOS(HWHO)
  CHARACTER(LEN=*), INTENT(IN) :: HWHO
  WRITE(0,*) 'FATAL_KARPOS_GRIB: ', TRIM(HWHO), ' appelée (build NetCDF-only).'
  FLUSH(0)
  ERROR STOP 'GRIB non compilé'
END SUBROUTINE GRIB_ABORT_KARPOS

! Jeu complet des externes grib_api référencés par le code SURFEX compilé
! (mode_read_grib, prep_grib_grid, read_grib_saf, read_grib_large_rain, …).
! Fonctions de LIBÉRATION/NETTOYAGE : no-op (appelées inconditionnellement en fin
! de prep_surf_atm même sans GRIB ouvert — inoffensives à ignorer).
SUBROUTINE GRIB_CLOSE_FILE();      END SUBROUTINE
SUBROUTINE GRIB_INDEX_RELEASE();   END SUBROUTINE
SUBROUTINE GRIB_RELEASE();         END SUBROUTINE
! Fonctions de LECTURE : abort (indiquent une vraie tentative de lecture GRIB,
! bypassée par les patches prep_hor_*_field / read_prep_file_date).
SUBROUTINE GRIB_OPEN_FILE();       CALL GRIB_ABORT_KARPOS('grib_open_file');       END SUBROUTINE
SUBROUTINE GRIB_NEW_FROM_FILE();   CALL GRIB_ABORT_KARPOS('grib_new_from_file');   END SUBROUTINE
SUBROUTINE GRIB_MULTI_SUPPORT_ON();END SUBROUTINE
SUBROUTINE GRIB_IS_MISSING();      CALL GRIB_ABORT_KARPOS('grib_is_missing');      END SUBROUTINE
SUBROUTINE GRIB_INDEX_CREATE();    CALL GRIB_ABORT_KARPOS('grib_index_create');    END SUBROUTINE
SUBROUTINE GRIB_INDEX_GET();       CALL GRIB_ABORT_KARPOS('grib_index_get');       END SUBROUTINE
SUBROUTINE GRIB_INDEX_GET_SIZE();  CALL GRIB_ABORT_KARPOS('grib_index_get_size');  END SUBROUTINE
SUBROUTINE GRIB_INDEX_SELECT();    CALL GRIB_ABORT_KARPOS('grib_index_select');    END SUBROUTINE
SUBROUTINE GRIB_NEW_FROM_INDEX();  CALL GRIB_ABORT_KARPOS('grib_new_from_index');  END SUBROUTINE
SUBROUTINE GRIB_GET();             CALL GRIB_ABORT_KARPOS('grib_get');             END SUBROUTINE
SUBROUTINE GRIB_GET_SIZE();        CALL GRIB_ABORT_KARPOS('grib_get_size');        END SUBROUTINE
