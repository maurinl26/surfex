"""Orchestration de la chaîne SURFEX OFFLINE : PGD → PREP → OFFLINE (→ SODA).

Assemble un répertoire de run (données ECOCLIMAP liées) et enchaîne les
programmes maîtres via le driver. Pensé pour un usage batch/CI : idempotent,
échoue franchement à la première étape en erreur, retourne un récapitulatif.
"""

from __future__ import annotations

import os
from pathlib import Path

from ._surfex import SurfexError
from . import driver

DEFAULT_STEPS = ("pgd", "prep", "offline")
_ECOCLIMAP_COVERS = (
    "ecoclimapI_covers_param.bin",
    "ecoclimapII_eu_covers_param.bin",
    "ecoclimapII_af_covers_param.bin",
)


def link_ecoclimap(workdir, ecoclimap_dir):
    """Lie les fichiers de cover ECOCLIMAP dans le répertoire de run (requis par PGD)."""
    workdir = Path(workdir)
    src = Path(ecoclimap_dir)
    if not src.is_dir():
        raise SurfexError(f"répertoire ECOCLIMAP inexistant : {src}")
    linked = []
    for cover in _ECOCLIMAP_COVERS:
        s = src / cover
        if not s.is_file():
            continue
        dst = workdir / cover
        if dst.exists() or dst.is_symlink():
            dst.unlink()
        os.symlink(s.resolve(), dst)
        linked.append(cover)
    if not linked:
        raise SurfexError(f"aucun cover ECOCLIMAP (*.bin) trouvé dans {src}")
    return linked


def run_chain(workdir, steps=DEFAULT_STEPS, ecoclimap_dir=None):
    """Enchaîne les étapes SURFEX dans workdir.

    Parameters
    ----------
    workdir : chemin du répertoire de run (namelists + forcing présents).
    steps   : sous-ensemble ordonné de ('pgd', 'prep', 'offline', 'soda').
    ecoclimap_dir : si fourni, lie les covers ECOCLIMAP avant PGD.

    Returns
    -------
    dict : {'workdir', 'steps', 'ecoclimap', 'outputs'} en cas de succès.
    Lève SurfexError (avec l'étape et la sortie SURFEX) à la première erreur.
    """
    workdir = Path(workdir)
    if not workdir.is_dir():
        raise SurfexError(f"répertoire de run inexistant : {workdir}")

    valid = {"pgd", "prep", "offline", "soda"}
    steps = tuple(s.lower() for s in steps)
    bad = [s for s in steps if s not in valid]
    if bad:
        raise SurfexError(f"étapes inconnues : {bad} (attendu ⊂ {sorted(valid)})")

    ecoclimap = []
    if ecoclimap_dir is not None and "pgd" in steps:
        ecoclimap = link_ecoclimap(workdir, ecoclimap_dir)

    for step in steps:
        getattr(driver, f"run_{step}")(workdir)   # lève SurfexError si échec

    outputs = sorted(
        f.name for f in workdir.iterdir()
        if f.suffix in (".nc", ".OUT") or f.name.startswith("LISTING")
    )
    return {
        "workdir": str(workdir),
        "steps": list(steps),
        "ecoclimap": ecoclimap,
        "outputs": outputs,
    }
