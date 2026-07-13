!SFX_KARPOS
! sfx_capi.F90 — Interface C stable (iso_c_binding) autour de SURFEX pour Karpos.
!
! Expose la surface OFFLINE (intégration forward) ET SODA (assimilation de surface),
! plus l'injection d'observations Sencrop et la lecture des champs. C'est cette
! couche que Cython enveloppe (_surfex.pyx) → numpy/xarray côté Python.
!
! Codes retour : 0 = OK, <0 = erreur (voir SFX_* ci-dessous).
!
! ⚠️  SCAFFOLD : la surface (ABI) est définitive et compile ; le câblage profond
!     vers les drivers SURFEX (INIT_SURF_ATM_n, boucle temporelle, SODA) est
!     marqué TODO et renvoie SFX_ENOTIMPL tant qu'il n'est pas branché.
!
MODULE SFX_CAPI
  USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_INT, C_DOUBLE, C_CHAR, C_NULL_CHAR
  IMPLICIT NONE
  PRIVATE

  ! Codes retour exposés côté C
  INTEGER(C_INT), PARAMETER, PUBLIC :: SFX_OK        =  0
  INTEGER(C_INT), PARAMETER, PUBLIC :: SFX_ENOTIMPL  = -1   ! pas encore câblé
  INTEGER(C_INT), PARAMETER, PUBLIC :: SFX_EINVAL    = -2   ! argument invalide
  INTEGER(C_INT), PARAMETER, PUBLIC :: SFX_ENOSTATE  = -3   ! sfx_init non appelé
  INTEGER(C_INT), PARAMETER, PUBLIC :: SFX_ERUN      = -4   ! échec run SURFEX

  LOGICAL, SAVE :: GINIT = .FALSE.

CONTAINS

  ! Convertit une chaîne C (terminée \0) en chaîne Fortran.
  FUNCTION c2f(cstr) RESULT(fstr)
    CHARACTER(KIND=C_CHAR), DIMENSION(*), INTENT(IN) :: cstr
    CHARACTER(LEN=:), ALLOCATABLE :: fstr
    INTEGER :: i
    fstr = ''
    i = 1
    DO
      IF (cstr(i) == C_NULL_CHAR) EXIT
      fstr = fstr // cstr(i)
      i = i + 1
      IF (i > 4096) EXIT
    END DO
  END FUNCTION c2f

  ! --- Cycle de vie -----------------------------------------------------------

  FUNCTION sfx_init() RESULT(ir) BIND(C, NAME='sfx_init')
    INTEGER(C_INT) :: ir
    ! TODO: init sans MPI (SFX_INIT / mise en place environnement offline)
    GINIT = .TRUE.
    ir = SFX_OK
  END FUNCTION sfx_init

  FUNCTION sfx_finalize() RESULT(ir) BIND(C, NAME='sfx_finalize')
    INTEGER(C_INT) :: ir
    GINIT = .FALSE.
    ir = SFX_OK
  END FUNCTION sfx_finalize

  ! --- Programmes maîtres SURFEX ---------------------------------------------
  ! Chaque étape prend le chemin d'un répertoire de namelists/inputs.

  FUNCTION sfx_run_pgd(cdir) RESULT(ir) BIND(C, NAME='sfx_run_pgd')
    CHARACTER(KIND=C_CHAR), DIMENSION(*), INTENT(IN) :: cdir
    INTEGER(C_INT) :: ir
    CHARACTER(LEN=:), ALLOCATABLE :: dir
    IF (.NOT. GINIT) THEN; ir = SFX_ENOSTATE; RETURN; END IF
    dir = c2f(cdir)
    ! TODO: appeler le driver PGD sur `dir` (physiographie : cover, sol, orographie)
    ir = SFX_ENOTIMPL
  END FUNCTION sfx_run_pgd

  FUNCTION sfx_run_prep(cdir) RESULT(ir) BIND(C, NAME='sfx_run_prep')
    CHARACTER(KIND=C_CHAR), DIMENSION(*), INTENT(IN) :: cdir
    INTEGER(C_INT) :: ir
    CHARACTER(LEN=:), ALLOCATABLE :: dir
    IF (.NOT. GINIT) THEN; ir = SFX_ENOSTATE; RETURN; END IF
    dir = c2f(cdir)
    ! TODO: appeler le driver PREP (état initial : T_sol, humidité, neige)
    ir = SFX_ENOTIMPL
  END FUNCTION sfx_run_prep

  FUNCTION sfx_run_offline(cdir) RESULT(ir) BIND(C, NAME='sfx_run_offline')
    CHARACTER(KIND=C_CHAR), DIMENSION(*), INTENT(IN) :: cdir
    INTEGER(C_INT) :: ir
    CHARACTER(LEN=:), ALLOCATABLE :: dir
    IF (.NOT. GINIT) THEN; ir = SFX_ENOSTATE; RETURN; END IF
    dir = c2f(cdir)
    ! TODO: intégration forward OFFLINE (forcing NetCDF → T2m, bilan surface)
    ir = SFX_ENOTIMPL
  END FUNCTION sfx_run_offline

  ! SODA : assimilation de surface. C'est le moteur du "nudging Sencrop"
  ! (OI / SEKF / ENKF / PF selon CASSIM_ISBA de la namelist).
  FUNCTION sfx_run_soda(cdir) RESULT(ir) BIND(C, NAME='sfx_run_soda')
    CHARACTER(KIND=C_CHAR), DIMENSION(*), INTENT(IN) :: cdir
    INTEGER(C_INT) :: ir
    CHARACTER(LEN=:), ALLOCATABLE :: dir
    IF (.NOT. GINIT) THEN; ir = SFX_ENOSTATE; RETURN; END IF
    dir = c2f(cdir)
    ! TODO: appeler le driver SODA (analyse : first-guess + obs → état corrigé)
    ir = SFX_ENOTIMPL
  END FUNCTION sfx_run_soda

  ! --- Assimilation Sencrop ---------------------------------------------------

  ! Injecte n observations de stations (lat/lon/valeur) pour la prochaine analyse
  ! SODA. Typiquement T2m Sencrop au temps d'analyse.
  FUNCTION sfx_set_obs(n, lats, lons, vals) RESULT(ir) BIND(C, NAME='sfx_set_obs')
    INTEGER(C_INT), VALUE, INTENT(IN) :: n
    REAL(C_DOUBLE), DIMENSION(*), INTENT(IN) :: lats, lons, vals
    INTEGER(C_INT) :: ir
    IF (.NOT. GINIT) THEN; ir = SFX_ENOSTATE; RETURN; END IF
    IF (n < 0) THEN; ir = SFX_EINVAL; RETURN; END IF
    ! TODO: stocker les obs → format attendu par SODA (OBS/CANARI) ; peupler
    !       l'opérateur d'observation H aux points-stations Sencrop.
    ir = SFX_ENOTIMPL
  END FUNCTION sfx_set_obs

  ! Lit un champ diagnostique dans un buffer fourni par l'appelant.
  ! name : 'T2M', 'WG1', 'TG1', ... ; nout : taille du buffer ; out : sortie.
  FUNCTION sfx_get_field(cname, out, nout) RESULT(ir) BIND(C, NAME='sfx_get_field')
    CHARACTER(KIND=C_CHAR), DIMENSION(*), INTENT(IN) :: cname
    INTEGER(C_INT), VALUE, INTENT(IN) :: nout
    REAL(C_DOUBLE), DIMENSION(*), INTENT(OUT) :: out
    INTEGER(C_INT) :: ir
    CHARACTER(LEN=:), ALLOCATABLE :: name
    IF (.NOT. GINIT) THEN; ir = SFX_ENOSTATE; RETURN; END IF
    name = c2f(cname)
    ! TODO: copier le champ SURFEX `name` (grille courante) dans out(1:nout)
    ir = SFX_ENOTIMPL
  END FUNCTION sfx_get_field

  ! Taille (nb de points de grille) du champ courant — pour dimensionner out.
  FUNCTION sfx_field_size() RESULT(nsize) BIND(C, NAME='sfx_field_size')
    INTEGER(C_INT) :: nsize
    ! TODO: renvoyer NDIM de la grille SURFEX active
    nsize = 0
  END FUNCTION sfx_field_size

END MODULE SFX_CAPI
