"""Ingestion des observations Sencrop → obs SURFEX (set_obs / SODA).

Lit le bulk Sencrop Karpos (S23) et le catalogue stations, joint les coordonnées,
sélectionne les obs à un temps d'analyse, et alimente `set_obs` (→ SODA).

Layout attendu (local ou s3://) :
    <root>/stations_integrated.csv        # bucket_id, latitude, longitude
    <root>/<year>.csv/part-*.csv          # station_id, timestamp, temperature (°C)

Le join se fait sur `timeseries.station_id == catalogue.bucket_id`.
Dépendances (extra `sencrop`) : pandas, s3fs (pour s3://).

Exemple :
    import karpos_surfex as ks
    from karpos_surfex import sencrop
    lats, lons, tK = sencrop.load_observations(
        "s3://karpos-backtest-data/sencrop", "2023-04-05T04:00:00")
    ks.set_obs(lats, lons, tK, "assim/OBS_sencrop.nc")
"""

from __future__ import annotations

from urllib.parse import urlparse

from ._surfex import SurfexError

# bbox du domaine Drôme-Ardèche (cf. domains/drome)
DROME_BBOX = {"lat_min": 44.0, "lat_max": 45.5, "lon_min": 4.0, "lon_max": 5.5}
STATIONS_FILE = "stations_integrated.csv"


def _pd():
    try:
        import pandas as pd
        return pd
    except ImportError as e:  # pragma: no cover
        raise SurfexError(
            "ingestion Sencrop : pip install 'karpos-surfex[sencrop]' (pandas, s3fs)"
        ) from e


def _is_remote(root: str) -> bool:
    return urlparse(str(root)).scheme in ("s3", "gs", "gcs", "az", "abfs")


def _join(root: str, *parts: str) -> str:
    if _is_remote(root):
        return "/".join([str(root).rstrip("/"), *parts])
    from pathlib import Path
    return str(Path(root).joinpath(*parts))


def _year_partition(root: str, year: int) -> str:
    """Unique part-*.csv dans <root>/<year>.csv/ (Spark)."""
    import fsspec

    pattern = _join(root, f"{year}.csv", "part-*.csv")
    fs, _ = fsspec.core.url_to_fs(str(root))
    matches = sorted(fs.glob(pattern))
    if not matches:
        raise SurfexError(f"aucune partition Sencrop pour {year} sous {root}")
    m = matches[0]
    if _is_remote(root) and "://" not in m:
        m = f"{urlparse(str(root)).scheme}://{m}"
    return m


def load_stations_catalog(root: str, bbox: dict | None = DROME_BBOX):
    """Catalogue stations (bucket_id, latitude, longitude), filtré bbox."""
    pd = _pd()
    df = pd.read_csv(_join(root, STATIONS_FILE))
    need = {"bucket_id", "latitude", "longitude"}
    missing = need - set(df.columns)
    if missing:
        raise SurfexError(f"{STATIONS_FILE} : colonnes manquantes {sorted(missing)}")
    if bbox:
        df = df[
            (df.latitude >= bbox["lat_min"]) & (df.latitude <= bbox["lat_max"])
            & (df.longitude >= bbox["lon_min"]) & (df.longitude <= bbox["lon_max"])
        ]
    return df[["bucket_id", "latitude", "longitude"]].dropna()


def load_observations(
    root: str,
    timestamp,
    bbox: dict | None = DROME_BBOX,
    tol_minutes: int = 30,
    station_only: bool = True,
):
    """Obs Sencrop à `timestamp` (± tol) → (lats, lons, temperatures_K) en numpy.

    Args:
        root: racine du bulk Sencrop (local ou s3://).
        timestamp: instant d'analyse (str ISO8601 ou pd.Timestamp, UTC).
        bbox: filtre spatial (défaut : domaine Drôme).
        tol_minutes: tolérance temporelle autour de `timestamp`.
        station_only: ne garder que temperature_source == 'station'.
    """
    import numpy as np
    pd = _pd()

    ts = pd.Timestamp(timestamp)
    if ts.tzinfo is None:
        ts = ts.tz_localize("UTC")
    year = ts.year

    df = pd.read_csv(_year_partition(root, year))
    df["timestamp"] = pd.to_datetime(df["timestamp"], utc=True, errors="coerce")
    if station_only and "temperature_source" in df.columns:
        df = df[df["temperature_source"] == "station"]

    window = pd.Timedelta(minutes=tol_minutes)
    df = df[(df["timestamp"] >= ts - window) & (df["timestamp"] <= ts + window)]
    if df.empty:
        raise SurfexError(f"aucune obs Sencrop à {ts} (± {tol_minutes} min)")

    # obs la plus proche du temps cible par station
    df = df.assign(_dt=(df["timestamp"] - ts).abs())
    df = df.sort_values("_dt").drop_duplicates("station_id", keep="first")

    cat = load_stations_catalog(root, bbox=bbox)
    merged = df.merge(cat, left_on="station_id", right_on="bucket_id", how="inner")
    merged = merged.dropna(subset=["latitude", "longitude", "temperature"])
    if merged.empty:
        raise SurfexError("aucune obs Sencrop après jointure catalogue + bbox")

    lats = merged["latitude"].to_numpy(dtype="f8")
    lons = merged["longitude"].to_numpy(dtype="f8")
    tK = merged["temperature"].to_numpy(dtype="f8") + 273.15   # °C → K
    return lats, lons, tK


def ingest_to_obs(root, timestamp, out_path="OBS_sencrop.nc", **kwargs):
    """Charge les obs Sencrop et écrit un fichier d'obs points (lat/lon/K)."""
    from . import driver
    lats, lons, tK = load_observations(root, timestamp, **kwargs)
    driver.set_obs(lats, lons, tK, out_path)
    return out_path, len(lats)


# --- Opérateur d'observation : stations → grille modèle (pour SODA) ----------

def to_grid(lats, lons, vals, pgd_path, method: str = "linear"):
    """Interpole des obs de stations sur la grille du PGD → champ 2D (yy, xx).

    SODA assimile des champs T2m/HU2m **grillés** (READ_SURF('NC','T2M',…)), pas
    des points. Cet opérateur projette les stations Sencrop sur la grille modèle.
    Interpolation `method` à l'intérieur, complétée au plus proche voisin ailleurs.
    """
    import numpy as np
    from netCDF4 import Dataset
    from scipy.interpolate import griddata

    with Dataset(str(pgd_path)) as ds:
        glat = np.asarray(ds.variables["LAT"][:], dtype="f8")   # (yy, xx)
        glon = np.asarray(ds.variables["LON"][:], dtype="f8")

    pts = np.column_stack([np.asarray(lons, "f8"), np.asarray(lats, "f8")])
    tgt = (glon, glat)
    grid = griddata(pts, np.asarray(vals, "f8"), tgt, method=method)
    holes = ~np.isfinite(grid)
    if holes.any():   # hors enveloppe convexe → plus proche voisin
        grid[holes] = griddata(pts, np.asarray(vals, "f8"), tgt, method="nearest")[holes]
    return grid


def write_obs_file(pgd_path, field, out_path, varname: str = "T2M"):
    """Écrit un fichier d'obs SODA : copie du PGD + champ `varname` (yy, xx).

    Réutilise la structure SURFEX valide du PGD.nc (grille, en-têtes) → lisible
    tel quel par READ_SURF('NC', varname, …).
    """
    import shutil
    import numpy as np
    from netCDF4 import Dataset

    shutil.copyfile(str(pgd_path), str(out_path))
    with Dataset(str(out_path), "a") as ds:
        if varname in ds.variables:
            v = ds.variables[varname]
        else:
            v = ds.createVariable(varname, "f8", ("yy", "xx"))
        v[:] = np.asarray(field, "f8")
        v.long_name = f"{varname} observation (Sencrop, grillé)"
    return str(out_path)


def ingest_to_soda_obs(root, timestamp, pgd_path, out_path="OBSERVATIONS.nc",
                       varname="T2M", method="linear", **kwargs):
    """Pipeline complet Sencrop → champ d'obs grillé SURFEX pour SODA.

    Charge les obs, les grille sur le PGD, écrit le NC d'obs (variable `varname`).
    Namelist SODA : &NAM_OBS COBS_M='T2M', NNCO(1)=1, CFILE_FORMAT_OBS='NC'.
    """
    lats, lons, tK = load_observations(root, timestamp, **kwargs)
    grid = to_grid(lats, lons, tK, pgd_path, method=method)
    write_obs_file(pgd_path, grid, out_path, varname=varname)
    return out_path, len(lats)
