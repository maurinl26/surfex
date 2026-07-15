#!/usr/bin/env python3
"""CERRA (Copernicus CDS, public) → forcing SURFEX OFFLINE.

CERRA `reanalysis-cerra-single-levels` est **public** (pas de convention MF) et
couvre 1984–2021 (+ extensions) → débloque les saisons de gel historiques, là où
l'AROME public n'est que temps-réel et l'AROME archive est verrouillé MF (#11).

On demande `grid=[0.05,0.05]` : CDS regrille la projection Lambert CERRA en
lat/lon régulier → NetCDF simple (pas de GRIB projeté à reprojeter).

Deux flux (comme la doc CERRA) :
  - **analyse** (instantané, 3-horaire) : t2m, r2, si10, wdir10, sp
  - **forecast** (accumulé, leadtimes) : ssrd, strd, tp   → flux instantanés

Mapping → forcing SURFEX :
  Tair=t2m · Qair=q(r2,t2m,sp) · Wind=si10 · Wind_DIR=wdir10 · PSurf=sp
  DIR/SCA_SWdown=ssrd · LWdown=strd · Rainf/Snowf=tp (phase par Tair)

Prérequis : ~/.cdsapirc (clé CDS), cdsapi, xarray, netCDF4, scipy.

    python cerra_to_forcing.py <workdir> PGD.nc 2023-04-04 2023-04-06
"""

from __future__ import annotations

import sys
from pathlib import Path

BBOX = {"lat_min": 44.0, "lat_max": 45.5, "lon_min": 4.0, "lon_max": 5.5}
AREA = [BBOX["lat_max"], BBOX["lon_min"], BBOX["lat_min"], BBOX["lon_max"]]  # N,W,S,E
GRID = [0.05, 0.05]
DATASET = "reanalysis-cerra-single-levels"
ANALYSIS_VARS = ["2m_temperature", "2m_relative_humidity",
                 "10m_wind_speed", "10m_wind_direction", "surface_pressure"]
FORECAST_VARS = ["surface_solar_radiation_downwards",
                 "surface_thermal_radiation_downwards", "total_precipitation"]
ANALYSIS_TIMES = ["00:00", "03:00", "06:00", "09:00", "12:00", "15:00", "18:00", "21:00"]


def _days(d0: str, d1: str):
    from datetime import date, timedelta
    a, b = date.fromisoformat(d0), date.fromisoformat(d1)
    out = []
    while a <= b:
        out.append(a); a += timedelta(days=1)
    return out


def download(workdir, d0, d1):
    """Télécharge analyse + forecast CERRA (bbox Drôme, 0,05°) sur [d0, d1]."""
    import cdsapi
    wd = Path(workdir); wd.mkdir(parents=True, exist_ok=True)
    days = _days(d0, d1)
    years = sorted({str(d.year) for d in days})
    months = sorted({f"{d.month:02d}" for d in days})
    daynums = sorted({f"{d.day:02d}" for d in days})
    c = cdsapi.Client()

    fa = wd / "cerra_analysis.nc"
    if not fa.exists():
        c.retrieve(DATASET, {
            "variable": ANALYSIS_VARS, "level_type": "surface_or_atmosphere",
            "data_type": "reanalysis", "product_type": "analysis",
            "year": years, "month": months, "day": daynums, "time": ANALYSIS_TIMES,
            "area": AREA, "grid": GRID, "data_format": "netcdf",
        }, str(fa))

    ff = wd / "cerra_forecast.nc"
    if not ff.exists():
        c.retrieve(DATASET, {
            "variable": FORECAST_VARS, "level_type": "surface_or_atmosphere",
            "data_type": "reanalysis", "product_type": "forecast",
            "year": years, "month": months, "day": daynums,
            "time": ["00:00", "06:00", "12:00", "18:00"],
            "leadtime_hour": ["1", "2", "3", "4", "5", "6"],
            "area": AREA, "grid": GRID, "data_format": "netcdf",
        }, str(ff))
    return fa, ff


def build_forcing(workdir, pgd_path, d0, d1):
    """Télécharge CERRA + interpole sur la grille PGD → FORCING.nc SURFEX."""
    import numpy as np
    import xarray as xr
    from netCDF4 import Dataset
    from scipy.interpolate import griddata

    wd = Path(workdir)
    fa, ff = download(wd, d0, d1)
    an = xr.open_dataset(fa)
    fc = xr.open_dataset(ff)

    # AXE TEMPS = forecast horaire (leadtimes 1-6 → continu). L'analyse 3-horaire
    # est interpolée temporellement dessus. La déaccumulation gère la remise à zéro
    # à chaque init (00/06/12/18) via le leadtime déduit de l'heure.
    from scipy.interpolate import interp1d

    def _t(ds):
        n = "valid_time" if "valid_time" in ds.coords else "time"
        return np.asarray(ds[n].values).astype("datetime64[s]")

    ft = _t(fc); order = np.argsort(ft); ft = ft[order]
    nt = ft.size; dt = 3600.0
    at = _t(an); aord = np.argsort(at); at = at[aord]

    def _latlon(ds):
        la = ds["latitude"].values; lo = ds["longitude"].values
        lo2, la2 = np.meshgrid(lo, la) if la.ndim == 1 else (lo, la)
        return la2.astype("f8"), lo2.astype("f8")

    sLAT, sLON = _latlon(fc)                       # analyse et forecast : même grille CDS 0,05°
    src = np.column_stack([sLON.ravel(), sLAT.ravel()])
    with Dataset(str(pgd_path)) as pgd:
        gLAT = np.asarray(pgd.variables["LAT"][:], "f8")
        gLON = np.asarray(pgd.variables["LON"][:], "f8")
        gZS = np.asarray(pgd.variables["ZS"][:], "f8")
    tgt = np.column_stack([gLON.ravel(), gLAT.ravel()])
    npF = tgt.shape[0]

    def _sp(field2d):                              # (ny,nx) → (npF,) sur PGD
        g = griddata(src, np.asarray(field2d, "f8").ravel(), tgt, method="linear")
        h = ~np.isfinite(g)
        if h.any():
            g[h] = griddata(src, np.asarray(field2d, "f8").ravel(), tgt, method="nearest")[h]
        return g.astype("f4")

    def analysis(var):                             # instant 3h → temps ft → PGD
        a = np.asarray(an[var].values)[aord]       # (nt_a, ny, nx)
        ta = (at - at[0]) / np.timedelta64(1, "s")
        tf = (ft - at[0]) / np.timedelta64(1, "s")
        ai = interp1d(ta, a, axis=0, bounds_error=False, fill_value=(a[0], a[-1]))(tf)
        return np.stack([_sp(ai[i]) for i in range(nt)]).astype("f4")

    def deaccum(var):                              # accumulé/init → flux horaire → PGD
        a = np.asarray(fc[var].values)[order]      # (nt, ny, nx)
        hod = ((ft - ft.astype("datetime64[D]")) / np.timedelta64(1, "h")).astype(int)
        lt = ((hod - 1) % 6) + 1                    # leadtime 1..6 depuis l'heure
        flux = np.empty_like(a, dtype="f8")
        for i in range(nt):
            flux[i] = a[i] / dt if (lt[i] == 1 or i == 0) else np.maximum(a[i] - a[i - 1], 0.0) / dt
        return np.stack([_sp(flux[i]) for i in range(nt)]).astype("f4")

    Tair = analysis(_pick(an, "t2m", "2m_temperature"))
    R2 = analysis(_pick(an, "r2", "2m_relative_humidity"))
    Wind = analysis(_pick(an, "si10", "10m_wind_speed"))
    WDIR = analysis(_pick(an, "wdir10", "10m_wind_direction"))
    PSurf = analysis(_pick(an, "sp", "surface_pressure"))
    Qair = _rh_to_q(R2, Tair, PSurf)

    SW = deaccum(_pick(fc, "ssrd", "surface_solar_radiation_downwards"))
    LW = deaccum(_pick(fc, "strd", "surface_thermal_radiation_downwards"))
    PR = deaccum(_pick(fc, "tp", "total_precipitation"))
    DIR_SW, SCA_SW = 0.85 * SW, 0.15 * SW
    Rainf = np.where(Tair >= 273.15, PR, 0.0)
    Snowf = np.where(Tair < 273.15, PR, 0.0)

    fpath = wd / "FORCING.nc"
    _write_forcing(fpath, ft, dt, gLON, gLAT, gZS, npF, nt,
                   Tair, Qair, Wind, WDIR, PSurf, LW, DIR_SW, SCA_SW, Rainf, Snowf)
    return str(fpath), nt, npF


def _pick(ds, *names):
    for n in names:
        if n in ds.variables:
            return n
    raise KeyError(f"aucune de {names} dans {list(ds.variables)}")


def _rh_to_q(rh_pct, t_k, p_pa):
    import numpy as np
    es = 611.2 * np.exp(17.67 * (t_k - 273.15) / (t_k - 29.65))   # Pa (Bolton)
    e = np.clip(rh_pct, 0, 100) / 100.0 * es
    return (0.622 * e / (p_pa - 0.378 * e)).astype("f4")



def _write_forcing(fpath, vt, dt, gLON, gLAT, gZS, npF, nt, Tair, Qair, Wind,
                   WDIR, PSurf, LW, DIR_SW, SCA_SW, Rainf, Snowf):
    import numpy as np
    from netCDF4 import Dataset
    with Dataset(str(fpath), "w", format="NETCDF3_64BIT_OFFSET") as nc:
        nc.createDimension("Number_of_points", npF)
        nc.createDimension("time", nt)
        t0 = vt[0].astype("datetime64[s]").item()
        tv = nc.createVariable("time", "f8", ("time",))
        tv.units = f"seconds since {t0.strftime('%Y-%m-%d %H:%M:%S')}"
        tv[:] = ((vt - vt[0]) / np.timedelta64(1, "s")).astype("f8")
        nc.createVariable("FRC_TIME_STP", "f8")[...] = dt

        def v(name, dims, data):
            x = nc.createVariable(name, "f4", dims); x[:] = np.asarray(data, "f4")
        v("LON", ("Number_of_points",), gLON.ravel())
        v("LAT", ("Number_of_points",), gLAT.ravel())
        v("ZS", ("Number_of_points",), gZS.ravel())
        v("ZREF", ("Number_of_points",), np.full(npF, 2.0))
        v("UREF", ("Number_of_points",), np.full(npF, 10.0))
        for nm, arr in (("Tair", Tair), ("Qair", Qair), ("Wind", Wind),
                        ("Wind_DIR", WDIR), ("PSurf", PSurf), ("LWdown", LW),
                        ("DIR_SWdown", DIR_SW), ("SCA_SWdown", SCA_SW),
                        ("Rainf", Rainf), ("Snowf", Snowf)):
            v(nm, ("time", "Number_of_points"), arr)
        v("CO2air", ("time", "Number_of_points"), np.full((nt, npF), 0.00062, "f4"))


if __name__ == "__main__":
    wd, pgd, d0, d1 = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4]
    path, nt, npF = build_forcing(wd, pgd, d0, d1)
    print(f"→ {path} : {nt} pas × {npF} points")
