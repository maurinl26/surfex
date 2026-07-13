# cython: language_level=3
"""Liaison Cython bas-niveau vers l'ABI C de SURFEX (voir capi/sfx_capi.F90)."""

import numpy as np
cimport numpy as cnp

cnp.import_array()

# --- Déclarations de l'ABI C (BIND(C) côté Fortran) --------------------------
cdef extern from *:
    """
    extern int sfx_init(void);
    extern int sfx_finalize(void);
    extern int sfx_run_pgd(const char *dir);
    extern int sfx_run_prep(const char *dir);
    extern int sfx_run_offline(const char *dir);
    extern int sfx_run_soda(const char *dir);
    extern int sfx_set_obs(int n, const double *lats, const double *lons, const double *vals);
    extern int sfx_get_field(const char *name, double *out, int nout);
    extern int sfx_field_size(void);
    """
    int sfx_init()
    int sfx_finalize()
    int sfx_run_pgd(const char *dir)
    int sfx_run_prep(const char *dir)
    int sfx_run_offline(const char *dir)
    int sfx_run_soda(const char *dir)
    int sfx_set_obs(int n, const double *lats, const double *lons, const double *vals)
    int sfx_get_field(const char *name, double *out, int nout)
    int sfx_field_size()


class SurfexError(RuntimeError):
    """Erreur renvoyée par une routine SURFEX (code retour < 0)."""


_ERR = {
    -1: "non implémenté (câblage driver SURFEX à venir)",
    -2: "argument invalide",
    -3: "sfx_init() non appelé",
    -4: "échec du run SURFEX",
}


cdef _check(int rc, str what):
    if rc != 0:
        raise SurfexError(f"{what}: {_ERR.get(rc, f'code {rc}')}")


def init():
    _check(sfx_init(), "init")


def finalize():
    _check(sfx_finalize(), "finalize")


def run_pgd(str directory):
    _check(sfx_run_pgd(directory.encode()), "run_pgd")


def run_prep(str directory):
    _check(sfx_run_prep(directory.encode()), "run_prep")


def run_offline(str directory):
    _check(sfx_run_offline(directory.encode()), "run_offline")


def run_soda(str directory):
    """Étape d'assimilation SODA (moteur du nudging Sencrop)."""
    _check(sfx_run_soda(directory.encode()), "run_soda")


def set_obs(cnp.ndarray[cnp.double_t, ndim=1] lats,
            cnp.ndarray[cnp.double_t, ndim=1] lons,
            cnp.ndarray[cnp.double_t, ndim=1] vals):
    """Injecte des observations de stations (Sencrop) pour la prochaine analyse."""
    cdef int n = lats.shape[0]
    if lons.shape[0] != n or vals.shape[0] != n:
        raise ValueError("lats/lons/vals doivent avoir la même longueur")
    lats = np.ascontiguousarray(lats)
    lons = np.ascontiguousarray(lons)
    vals = np.ascontiguousarray(vals)
    _check(sfx_set_obs(n, &lats[0], &lons[0], &vals[0]), "set_obs")


def get_field(str name):
    """Retourne un champ diagnostique SURFEX (ex. 'T2M') en numpy 1-D."""
    cdef int n = sfx_field_size()
    if n <= 0:
        raise SurfexError("get_field: grille non initialisée (field_size=0)")
    cdef cnp.ndarray[cnp.double_t, ndim=1] out = np.empty(n, dtype=np.double)
    _check(sfx_get_field(name.encode(), &out[0], n), "get_field")
    return out
