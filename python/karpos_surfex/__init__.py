"""karpos_surfex — SURFEX OFFLINE + assimilation (SODA) empaqueté pour Karpos.

Surface haut-niveau au-dessus de l'ABI C (voir ``capi/sfx_capi.F90``) :

    import karpos_surfex as ks

    ks.init()
    ks.run_pgd("domain/")          # physiographie (une fois)
    ks.run_prep("prep/")           # état initial
    ks.run_offline("run/")         # intégration forward (first-guess)

    # Assimilation Sencrop via SODA :
    ks.set_obs(lats, lons, t2m)    # obs stations Sencrop
    ks.run_soda("assim/")          # analyse (OI / SEKF selon namelist)
    t2m_grid = ks.get_field("T2M") # champ analysé → numpy

L'assimilation Sencrop n'est PAS un module séparé : elle passe par le moteur
d'assimilation de surface natif de SURFEX (SODA — OI / SEKF / ENKF / PF). Voir
``karpos-downscaling#30`` (EPIC D) pour le contexte scientifique.
"""

from ._surfex import (  # noqa: F401
    SurfexError,
    init,
    finalize,
    run_pgd,
    run_prep,
    run_offline,
    run_soda,
    set_obs,
    get_field,
)

__all__ = [
    "SurfexError",
    "init",
    "finalize",
    "run_pgd",
    "run_prep",
    "run_offline",
    "run_soda",
    "set_obs",
    "get_field",
]

__version__ = "0.0.1.dev0"
