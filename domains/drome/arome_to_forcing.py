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

⚠️  SCAFFOLD : le fetch + l'inventaire sont opérationnels ; l'assemblage du
FORCING.nc (conversions unités + accumulés→flux + split SW + phase précip) et le
run PREP/OFFLINE sont la prochaine étape (voir write_forcing / TODO).
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


# TODO — assemblage FORCING.nc :
#   * concaténer les LEAD_RANGES (0..24h+), aligner le temps
#   * accumulés (ssrd/strd/tp) → flux instantanés (diff / Δt)
#   * split DIR_SWdown / SCA_SWdown (ex. modèle Erbs sur SWdown)
#   * phase précip Rainf/Snowf selon Tair (seuil ~273.15 K)
#   * aplatir (lat,lon)→Number_of_points ; LAT/LON/ZS(depuis DEM)/ZREF=2/UREF=10
#   * FRC_TIME_STP = 3600 s
#   → écrire FORCING.nc (format SURFEX OFFLINE), puis PREP (init uniforme) + OFFLINE.


if __name__ == "__main__":
    run = latest_run()
    print(f"dernier run AROME : {_run_iso(run)}")
    print("packages/variables :", PKG_VARS)
