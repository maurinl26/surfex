# Domaine Karpos — Drôme-Ardèche (SURFEX 1 km)

Grille SURFEX pour la détection de gel arboricole Drôme-Ardèche (EPIC D).

| Paramètre | Valeur |
|---|---|
| Projection | Lambert conforme (`CONF PROJ`), tangente 45°N |
| Centre | 44,75°N / 4,75°E |
| Résolution | **1 km** (`XDX=XDY=1000 m`) |
| Grille | **120 × 168** points (~120 × 168 km) |
| Emprise | ~44–45,5°N × 4–5,5°E |
| Cover | ECOCLIMAP-II (573 classes) |

## État

- ✅ **Grille + orographie réelle validées** : `PGD` produit `PGD.nc` (120×168) avec
  l'orographie **IGN BD ALTI** — ZS ∈ [20,6 ; 1826,9] m (fonds de vallée → Diois/Vercors).
- ⏳ **Cover réel** : encore uniforme (`NAM_COVER XUNIF_COVER(4)`) — attend la
  carte ECOCLIMAP-II.

## Orographie — DEM IGN (opérationnel)

```bash
# DEM sur S3 : s3://karpos-backtest-data/downscaling/dem/dem_attributes.nc (IGN BD ALTI)
python dem_to_surfex.py dem_attributes.nc DROME_ZS   # → DROME_ZS.dir + .hdr
# namelist : &NAM_ZS YZS='DROME_ZS', YZSFILETYPE='DIRECT'
```

`dem_to_surfex.py` convertit le NetCDF (elevation, lat, lon) en raster SURFEX
`DIRECT` (INTEGER*2 big-endian) ; PGD (`average_orography`) l'agrège à 1 km.

## Données restantes

- **Carte de cover ECOCLIMAP-II** (`LECOCLIMAP=T` + `NAM_COVER YCOVER=...`) : raster
  global (~Go), source Météo-France / CNRM (cf. licence #10 pour SG). En attendant,
  cover uniforme prescrit.

## Lancer

```bash
export SURFEX_EXE_DIR=<repo>/exe
python -m karpos_surfex run \
  --workdir domains/drome \
  --ecoclimap <repo>/MY_RUN/ECOCLIMAP \
  --steps pgd
```

Puis `--steps pgd,prep,offline` une fois PREP (init NetCDF/uniforme) et le
forcing NetCDF (CERRA, #4) en place.

## Forcing AROME (#4) & run OFFLINE

- `arome_to_forcing.py` — AROME public (OVH) → `FORCING.nc` sur la grille PGD 1 km.
  **Opérationnel** (testé : 6 pas × 20160 pts, valeurs plausibles). SWdown≈ssr
  (proxy), précip=0 (nuits de gel sèches) — à raffiner.
- `OPTIONS_run.nam` — namelist PGD + PREP (uniforme) + OFFLINE (forcing NC).

### ⛔ Blocage OFFLINE : PREP veut un first-guess GRIB

PGD ✅. PREP initialise WGI/WR/SN_VEG/LAI en uniforme, mais l'**état sol primaire
(WG/TG)** réclame un fichier atmosphérique GRIB (`HFILETYPE`), non contournable via
namelist (testé : XUNIF, ISBA-seul, `CISBA='3-L'`, `CFILE/CFILETYPE` blancs).
Comportement de la distribution — même le cas de référence `cdp9697` gribe.

**Débloquage (prochaine session focalisée)** :
1. Construire **eccodes** (support GRIB réel — supprime le stub) puis fournir un
   first-guess GRIB (analyse ECMWF/AROME sol, ou climatologie), OU
2. Patcher SURFEX pour autoriser un PREP 100 % uniforme (init sans fichier).

Tout l'amont est prêt (grille + orographie + forcing AROME + obs Sencrop).

## Fichiers

- `OPTIONS.nam` — namelist PGD (smoke). `OPTIONS_run.nam` — chaîne complète.
