"""Driver fonctionnel : pilote les exécutables SURFEX construits + lit les sorties.

Approche opérationnelle (batch nocturne Karpos / job-array HPC) : chaque étape lance
le programme maître réel (PGD/PREP/OFFLINE/SODA) dans un répertoire de travail
contenant ses namelists (`OPTIONS.nam`, …) et forcings, puis on relit les sorties
NetCDF en numpy.

Le shim Fortran in-process (`_surfex`) reste la voie « tight-loop » pour plus tard ;
ce driver rend l'API fonctionnelle dès maintenant, sans ré-embarquer le driver
offline.F90 (1845 lignes, branches MPI/XIOS/OASIS).
"""

from __future__ import annotations

import os
import shutil
import subprocess
from pathlib import Path

from ._surfex import SurfexError

MASTERS = ("PGD", "PREP", "OFFLINE", "SODA")


def _exe_dir() -> Path | None:
    """Répertoire des exécutables : $SURFEX_EXE_DIR, sinon un exe/ ancêtre."""
    d = os.environ.get("SURFEX_EXE_DIR")
    if d:
        return Path(d)
    here = Path(__file__).resolve()
    for parent in here.parents:
        cand = parent / "exe"
        if cand.is_dir() and any(cand.glob("OFFLINE*")):
            return cand
    return None


def find_exe(name: str) -> str:
    """Localise un programme maître (ex. 'OFFLINE'). Gère le suffixe -ARCH-...

    Cherche dans $SURFEX_EXE_DIR, un exe/ ancêtre, puis le PATH.
    """
    name = name.upper()
    d = _exe_dir()
    if d is not None:
        matches = [
            m for m in sorted(d.glob(f"{name}*"))
            if m.is_file() and os.access(m, os.X_OK)
        ]
        if matches:
            return str(matches[0])
    found = shutil.which(name)
    if found:
        return found
    raise SurfexError(
        f"exécutable {name} introuvable — lancer ../build.sh puis "
        f"export SURFEX_EXE_DIR=<repo>/exe (ou mettre les binaires dans le PATH)"
    )


def _run(name: str, workdir: str | os.PathLike) -> int:
    """Exécute un maître SURFEX dans workdir (qui doit contenir ses namelists)."""
    exe = find_exe(name)
    wd = Path(workdir)
    if not wd.is_dir():
        raise SurfexError(f"répertoire de travail inexistant : {wd}")
    proc = subprocess.run(
        [exe], cwd=str(wd), capture_output=True, text=True, check=False
    )
    if proc.returncode != 0:
        tail = "\n".join((proc.stdout + proc.stderr).splitlines()[-20:])
        raise SurfexError(f"{name} a échoué (rc={proc.returncode}) dans {wd} :\n{tail}")
    return proc.returncode


def run_pgd(workdir):
    """Physiographie (cover, sol, orographie). Lit PRE_PGD1.nam / OPTIONS.nam."""
    return _run("PGD", workdir)


def run_prep(workdir):
    """État initial (T_sol, humidité, neige)."""
    return _run("PREP", workdir)


def run_offline(workdir):
    """Intégration forward OFFLINE (forcing → T2m, bilan de surface)."""
    return _run("OFFLINE", workdir)


def run_soda(workdir):
    """Assimilation de surface (SODA — OI/SEKF/ENKF/PF selon CASSIM_ISBA)."""
    return _run("SODA", workdir)


def get_field(name: str, ncfile: str | os.PathLike):
    """Lit une variable d'un fichier de sortie SURFEX NetCDF → numpy (squeezé)."""
    try:
        from netCDF4 import Dataset
    except ImportError as exc:  # pragma: no cover
        raise SurfexError(
            "get_field requiert netCDF4 : pip install 'karpos-surfex[netcdf]'"
        ) from exc
    import numpy as np

    with Dataset(str(ncfile)) as ds:
        if name not in ds.variables:
            dispo = list(ds.variables)[:25]
            raise SurfexError(f"variable '{name}' absente de {ncfile}. Dispo : {dispo}")
        return np.asarray(ds.variables[name][:]).squeeze()


def set_obs(lats, lons, vals, path: str | os.PathLike = "OBS_sencrop.nc"):
    """Écrit des observations de stations (Sencrop) dans un NetCDF pour SODA.

    Produit un fichier (lat, lon, value) que la namelist SODA doit référencer.
    Le format exact d'ingestion est propre au cas (NAM_OBS / CANARI) — voir
    docs et la calibration EPIC D ; ici on matérialise les obs de façon portable.
    """
    try:
        from netCDF4 import Dataset
    except ImportError as exc:  # pragma: no cover
        raise SurfexError(
            "set_obs requiert netCDF4 : pip install 'karpos-surfex[netcdf]'"
        ) from exc
    import numpy as np

    lats = np.ascontiguousarray(lats, dtype="f8")
    lons = np.ascontiguousarray(lons, dtype="f8")
    vals = np.ascontiguousarray(vals, dtype="f8")
    if not (lats.shape == lons.shape == vals.shape) or lats.ndim != 1:
        raise SurfexError("lats/lons/vals doivent être 1-D de même longueur")

    with Dataset(str(path), "w", format="NETCDF4") as ds:
        ds.createDimension("station", lats.size)
        for nm, arr, unit in (
            ("lat", lats, "degrees_north"),
            ("lon", lons, "degrees_east"),
            ("obs", vals, "K"),
        ):
            v = ds.createVariable(nm, "f8", ("station",))
            v.units = unit
            v[:] = arr
        ds.source = "karpos-surfex set_obs (Sencrop)"
    return str(path)
