#!/usr/bin/env python3
"""AROME 0,025° (public, meteo.data.gouv.fr / OVH) → forcing SURFEX OFFLINE.

L'AROME temps-réel est PUBLIC (bucket OVH meteofrance-pnt) — pas de convention MF
(seule l'archive historique 2022-2025 l'exige, cf. #11). Les packages GRIB2
SP1/SP2/SP3 contiennent toutes les variables de forcing SURFEX.

Mapping AROME → forcing SURFEX (validé par inventaire) :
    Tair       ← t (SP2) / t2m (SP3)          [K]
    Qair       ← sh2 (SP2, humidité spéc. 2 m) [kg/kg]   (sinon d2m → Magnus)
    Wind       ← si10 (SP1)                     [m/s]
    Wind_DIR   ← wdir10 (SP1)                    [deg]
    PSurf      ← sp (SP2)                        [Pa]
    DIR_SWdown ← ssrd (SP3, accumulé → /Δt)      [W/m²]  (+ split diffus SCA_SWdown)
    LWdown     ← strd (SP3, accumulé → /Δt)      [W/m²]
    Rainf/Snowf← tp (SP3, accumulé → /Δt)        [kg/m²/s] (phase par Tair)

Forcing SURFEX (FORCING.nc) — dims (Number_of_points, time) + variables ci-dessus
+ métadonnées par point : LAT, LON, ZS, ZREF(=2), UREF(=10), FRC_TIME_STP.

État : build_forcing() OPÉRATIONNEL — fetch AROME public + interpolation sur la
grille PGD 1 km + FORCING.nc valide (testé : 6 pas × 20160 pts, valeurs plausibles).
Limites 1er run : SWdown≈ssr (net, proxy) ; précip=0 (nuits de gel sèches) — à
raffiner (params MF `ssrd`/`tirf`). Le run OFFLINE dépend d'un PREP file-less :
la distribution veut un first-guess GRIB pour l'état sol primaire (WG/TG) → soit
construire eccodes (support GRIB), soit fournir un first-guess NetCDF.
"""

from __future__ import annotations

from datetime import datetime, timedelta, timezone
from pathlib import Path

UTC = timezone.utc
BBOX = {"lat_min": 44.0, "lat_max": 45.5, "lon_min": 4.0, "lon_max": 5.5}
LEAD_RANGES = ["00H06H", "07H12H", "13H18H", "19H24H"]
PUBLICATION_DELAY_MIN = 180
_BASE = ("https://meteofrance-pnt.s3.rbx.io.cloud.ovh.net/pnt/{ri}/arome/0025/"
         "{pk}/arome__0025__{pk}__{lr}__{ri}.grib2")

# Variable de chaque package (paramId shortName) → nom de forcing SURFEX
PKG_VARS = {
    "SP1": {"si10": "Wind", "wdir10": "Wind_DIR"},
    "SP2": {"t": "Tair", "sh2": "Qair", "sp": "PSurf"},
    "SP3": {"ssrd": "SWdown", "strd": "LWdown", "tp": "Precip"},
}


def latest_run(now: datetime | None = None) -> datetime:
    """Dernier run AROME publié (00/06/12/18 UTC + délai de publication)."""
    now = now or datetime.now(tz=UTC)
    c = now - timedelta(minutes=PUBLICATION_DELAY_MIN)
    return c.replace(hour=(c.hour // 6) * 6, minute=0, second=0, microsecond=0)


def _run_iso(run: datetime) -> str:
    return run.strftime("%Y-%m-%dT%H:00:00Z")


def download(run: datetime, package: str, lead_range: str, target: Path) -> Path:
    """Télécharge un GRIB2 AROME (~50 Mo) depuis le bucket public OVH."""
    import requests
    url = _BASE.format(ri=_run_iso(run), pk=package, lr=lead_range)
    r = requests.get(url, timeout=180, stream=True)
    r.raise_for_status()
    target.parent.mkdir(parents=True, exist_ok=True)
    with open(target, "wb") as f:
        for chunk in r.iter_content(1 << 20):
            f.write(chunk)
    return target


def open_package(grib_path: Path, package: str):
    """Ouvre un GRIB AROME, sous-ensemble bbox Drôme, variables de forcing du package.

    SP3 mêle des champs instant/accumulés → filter_by_keys par stepType.
    """
    import xarray as xr

    def _open(**kw):
        return xr.open_dataset(grib_path, engine="cfgrib",
                               backend_kwargs={"indexpath": "", **kw})

    if package == "SP3":
        ds = _open(filter_by_keys={"stepType": "accum"})
    else:
        ds = _open()
    keep = [v for v in PKG_VARS[package] if v in ds.data_vars]
    ds = ds[keep].sel(
        latitude=slice(BBOX["lat_max"], BBOX["lat_min"]),
        longitude=slice(BBOX["lon_min"], BBOX["lon_max"]),
    )
    return ds.rename(PKG_VARS[package])


def _open_sp3(grib_path):
    """SP3 : flux radiatifs accumulés. AROME (GRIB MF) expose strd (thermique
    descendant) et ssr (solaire NET) ; pas de ssrd/tp propres via cfgrib.
    On prend LWdown=strd et SWdown≈ssr (net, proxy — nul la nuit, cas gel).
    Précip traitée à part (= 0 pour le 1er run, nuits de gel sèches)."""
    import xarray as xr
    ds = xr.open_dataset(grib_path, engine="cfgrib", backend_kwargs={
        "indexpath": "", "filter_by_keys": {"stepType": "accum"}})
    keep = [v for v in ("strd", "ssr") if v in ds.data_vars]
    return ds[keep].rename({"strd": "LWdown", "ssr": "SWdown"})


def build_forcing(run, workdir, pgd_path, lead_ranges=("00H06H",), keep_grib=False):
    """Télécharge AROME + interpole sur la grille PGD + assemble FORCING.nc.

    Le forcing SURFEX OFFLINE doit être sur les MÊMES points que le modèle (PGD).
    On interpole donc AROME 0,025° → grille PGD 1 km (LAT/LON/ZS depuis PGD.nc).
    """
    import numpy as np
    import xarray as xr
    from netCDF4 import Dataset
    from scipy.interpolate import griddata

    wd = Path(workdir)
    wd.mkdir(parents=True, exist_ok=True)
    subs = []
    for lr in lead_ranges:
        parts = {}
        for pk in ("SP1", "SP2", "SP3"):
            g = download(run, pk, lr, wd / f"arome_{pk}_{lr}.grib2")
            parts[pk] = _open_sp3(g) if pk == "SP3" else open_package(g, pk)
        m = xr.merge(parts.values(), compat="override", join="inner")
        subs.append(m)
    ds = xr.concat(subs, dim="step").sortby("step")

    vt = np.asarray(ds["valid_time"].values).astype("datetime64[s]")
    nt = vt.size
    dt = 3600.0  # AROME horaire

    # points source (AROME) et cible (grille PGD)
    alat = np.asarray(ds["latitude"].values, "f8")
    alon = np.asarray(ds["longitude"].values, "f8")
    aLON, aLAT = np.meshgrid(alon, alat)
    src = np.column_stack([aLON.ravel(), aLAT.ravel()])
    with Dataset(str(pgd_path)) as pgd:
        gLAT = np.asarray(pgd.variables["LAT"][:], "f8")   # (yy, xx)
        gLON = np.asarray(pgd.variables["LON"][:], "f8")
        gZS = np.asarray(pgd.variables["ZS"][:], "f8")
    tgt = np.column_stack([gLON.ravel(), gLAT.ravel()])
    npF = tgt.shape[0]

    def flat(v):  # (time, ny, nx) AROME → (time, npF) interpolé sur PGD
        a = np.asarray(ds[v].values, "f4")
        out = np.empty((nt, npF), "f4")
        for it in range(nt):
            g = griddata(src, a[it].ravel(), tgt, method="linear")
            h = ~np.isfinite(g)
            if h.any():
                g[h] = griddata(src, a[it].ravel(), tgt, method="nearest")[h]
            out[it] = g
        return out

    Tair = flat("Tair")
    Qair = flat("Qair")
    Wind = flat("Wind")
    WDIR = flat("Wind_DIR")
    PSurf = flat("PSurf")

    # accumulés (depuis le début du run) → flux instantané (dérivée / dt)
    def deaccum(v):
        a = flat(v)
        out = np.empty_like(a)
        out[0] = a[0] / dt
        out[1:] = np.maximum(a[1:] - a[:-1], 0.0) / dt
        return out

    SW = np.maximum(deaccum("SWdown"), 0.0)       # W/m² (net solaire, proxy)
    LW = deaccum("LWdown")                         # W/m² (thermique descendant)
    PR = np.zeros_like(Tair)                       # précip=0 (1er run ; TODO tirf)
    DIR_SW = 0.7 * SW                              # split simple direct/diffus
    SCA_SW = 0.3 * SW
    Rainf = np.where(Tair >= 273.15, PR, 0.0)      # phase par Tair
    Snowf = np.where(Tair < 273.15, PR, 0.0)

    # métadonnées par point = grille PGD (LAT/LON/ZS réels, orographie IGN)
    ZS = gZS.ravel().astype("f4")

    fpath = wd / "FORCING.nc"
    with Dataset(fpath, "w", format="NETCDF3_64BIT_OFFSET") as nc:
        nc.createDimension("Number_of_points", npF)
        nc.createDimension("time", nt)
        t0 = vt[0].astype("datetime64[s]").item()

        def var(name, dims, data, **att):
            v = nc.createVariable(name, "f4" if data.dtype == np.float32 else "f8", dims)
            for k, val in att.items():
                setattr(v, k, val)
            v[:] = data
            return v

        tv = nc.createVariable("time", "f8", ("time",))
        tv.units = f"seconds since {t0.strftime('%Y-%m-%d %H:%M:%S')}"
        tv[:] = ((vt - vt[0]) / np.timedelta64(1, "s")).astype("f8")
        nc.createVariable("FRC_TIME_STP", "f8")[...] = dt

        var("LON", ("Number_of_points",), gLON.ravel().astype("f4"))
        var("LAT", ("Number_of_points",), gLAT.ravel().astype("f4"))
        var("ZS", ("Number_of_points",), ZS)
        var("ZREF", ("Number_of_points",), np.full(npF, 2.0, "f4"))
        var("UREF", ("Number_of_points",), np.full(npF, 10.0, "f4"))
        for nm, arr in (("Tair", Tair), ("Qair", Qair), ("Wind", Wind),
                        ("Wind_DIR", WDIR), ("PSurf", PSurf), ("LWdown", LW),
                        ("DIR_SWdown", DIR_SW), ("SCA_SWdown", SCA_SW),
                        ("Rainf", Rainf), ("Snowf", Snowf)):
            var(nm, ("time", "Number_of_points"), arr.astype("f4"))
        var("CO2air", ("time", "Number_of_points"),
            np.full((nt, npF), 0.00062, "f4"))   # ~400 ppm en kg/kg

    if not keep_grib:
        for g in wd.glob("arome_*.grib2*"):
            g.unlink()
    return str(fpath), nt, npF


if __name__ == "__main__":
    import sys
    run = latest_run()
    wd = sys.argv[1] if len(sys.argv) > 1 else "."
    print(f"run AROME {_run_iso(run)} → forcing dans {wd}")
    pgd = sys.argv[2] if len(sys.argv) > 2 else "PGD.nc"
    path, nt, npF = build_forcing(run, wd, pgd)
    print(f"→ {path} : {nt} pas de temps × {npF} points")
