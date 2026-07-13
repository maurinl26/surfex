#!/usr/bin/env python3
"""Convertit un DEM NetCDF (elevation, lat, lon) en orographie SURFEX 'DIRECT'.

Produit <out>.dir (INTEGER*2 big-endian, nord→sud, ouest→est) + <out>.hdr, format
lu par read_direct_gauss.F90 (YZSFILETYPE='DIRECT').

    python dem_to_surfex.py dem_attributes.nc DROME_ZS
    # → DROME_ZS.dir + DROME_ZS.hdr  (namelist : &NAM_ZS YZS='DROME_ZS', YZSFILETYPE='DIRECT')

Source DEM Karpos : s3://karpos-backtest-data/downscaling/dem/dem_attributes.nc
(IGN BD ALTI, grille Drôme 118×167).
"""

from __future__ import annotations

import sys
import numpy as np
from netCDF4 import Dataset

NODATA = -9999


def convert(nc_path: str, out_base: str, var: str = "elevation") -> None:
    with Dataset(nc_path) as ds:
        z = np.asarray(ds.variables[var][:], dtype="f8")          # (y, x)
        lat = np.asarray(ds.variables["lat"][:], dtype="f8")       # (y,)
        lon = np.asarray(ds.variables["lon"][:], dtype="f8")       # (x,)

    nrows, ncols = z.shape
    assert lat.size == nrows and lon.size == ncols, "dims lat/lon incohérentes"

    # Ordonner nord→sud (ligne 1 = latitude max) et ouest→est.
    if lat[0] < lat[-1]:            # lat ascendante (sud→nord) → on retourne
        z = z[::-1, :]
        lat = lat[::-1]
    if lon[0] > lon[-1]:
        z = z[:, ::-1]
        lon = lon[::-1]

    dlat = abs(lat[1] - lat[0])
    dlon = abs(lon[1] - lon[0])
    north = lat.max() + dlat / 2.0
    south = lat.min() - dlat / 2.0
    west = lon.min() - dlon / 2.0
    east = lon.max() + dlon / 2.0

    # int16 big-endian ; NaN/inf → nodata ; arrondi au mètre.
    zi = np.where(np.isfinite(z), np.rint(z), NODATA)
    zi = np.clip(zi, -32000, 32000).astype(">i2")
    zi[~np.isfinite(z)] = NODATA

    with open(f"{out_base}.dir", "wb") as f:
        f.write(zi.tobytes(order="C"))     # ligne par ligne, nord→sud

    # NB : read_direct_gauss cherche 'recordtype' en 1re ligne ; READHEAD ignore
    # la 1re ligne (commentaire) puis lit les clés → recordtype DOIT être la 1re.
    hdr = (
        f"recordtype: integer 16 bits\n"
        f"nodata: {NODATA}\n"
        f"north: {north:.6f}N\n"
        f"south: {south:.6f}N\n"
        f"east: {east:.6f}E\n"
        f"west: {west:.6f}E\n"
        f"rows: {nrows}\n"
        f"cols: {ncols}\n"
    )
    with open(f"{out_base}.hdr", "w") as f:
        f.write(hdr)

    print(f"→ {out_base}.dir ({zi.nbytes} o) + {out_base}.hdr")
    print(f"   grille {nrows}×{ncols}, N={north:.3f} S={south:.3f} "
          f"W={west:.3f} E={east:.3f}, z∈[{np.nanmin(z):.0f},{np.nanmax(z):.0f}] m")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit("usage: dem_to_surfex.py <dem.nc> <out_base>")
    convert(sys.argv[1], sys.argv[2])
