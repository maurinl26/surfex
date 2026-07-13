"""karpos_surfex — SURFEX OFFLINE + assimilation (SODA) empaqueté pour Karpos.

API fonctionnelle (pilote les exécutables SURFEX construits par ``build.sh``) :

    import karpos_surfex as ks

    ks.run_pgd("domain/")          # physiographie (une fois)
    ks.run_prep("prep/")           # état initial
    ks.run_offline("run/")         # intégration forward → sortie NetCDF
    t2m = ks.get_field("T2M", "run/SURF_ATM_DIAGNOSTICS.OUT.nc")

    # Assimilation Sencrop via SODA :
    ks.set_obs(lats, lons, sencrop_t2m, "assim/OBS_sencrop.nc")
    ks.run_soda("assim/")          # analyse (OI / SEKF selon namelist)

Les exécutables sont localisés via ``$SURFEX_EXE_DIR`` (ou un ``exe/`` ancêtre).
L'assimilation Sencrop passe par le moteur natif SURFEX (SODA), pas un module
ad hoc — cf. ``karpos-downscaling#30`` (EPIC D).

``init``/``finalize`` et le module ``_surfex`` (Cython/Fortran) constituent la
voie in-process (tight-loop) pour l'avenir ; le driver subprocess rend l'API
fonctionnelle dès aujourd'hui (modèle batch nocturne / job-array HPC).
"""

from ._surfex import SurfexError, init, finalize  # noqa: F401
from .driver import (  # noqa: F401
    find_exe,
    get_field,
    run_offline,
    run_pgd,
    run_prep,
    run_soda,
    set_obs,
)

__all__ = [
    "SurfexError",
    "init",
    "finalize",
    "find_exe",
    "run_pgd",
    "run_prep",
    "run_offline",
    "run_soda",
    "set_obs",
    "get_field",
]

__version__ = "0.0.1.dev0"
