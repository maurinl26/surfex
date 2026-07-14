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

### ✅ Chaîne complète PGD → PREP → OFFLINE (voie 2 : PREP file-less patché)

**La chaîne tourne de bout en bout** avec forcing AROME réel :
```bash
python domains/drome/arome_to_forcing.py <workdir> PGD.nc   # FORCING.nc AROME
python -m karpos_surfex run --workdir <workdir> --ecoclimap MY_RUN/ECOCLIMAP \
       --steps pgd,prep,offline
# → SURF_ATM_DIAGNOSTICS.OUT.nc : TS, TSRAD (T_skin), T2M, HU2M, RN, H, LE… 120×168
```

**PREP file-less** obtenu par patch SURFEX (repli init uniforme au lieu du GRIB) :
`prep_hor_isba_field`, `prep_hor_snow_field`, `read_prep_file_date`,
`prep_surf_atm` (CLEAR_GRIB_INDEX désactivé) + stub grib (release no-op) +
`build.sh` (relink forcé — sinon la modif du stub seul n'est pas relinkée).

Validé : run AROME 2026-07-13 13Z → TS/TSRAD ∈ [23,9 ; 38,9] °C (après-midi été).
Le `PREP.nc` produit sert de first-guess pour SODA (#7).

⚠️ Le temps du PREP (`XTIME`) doit == le 1er pas du forcing (le join SP3 accumulé
décale d'1 h). Limites forcing : SWdown proxy, précip=0 (à raffiner, #4).

## Assimilation Sencrop (#7) — SODA

État : **mécanisme SODA validé**, blocage sur le forcing de période de gel (data, pas code).

- Ingestion Sencrop → obs grillées → **injectées dans l'état SURFEX** (`sencrop.py`).
- **SODA tourne** (`OPTIONS_soda.nam`, nuit de gel 2023-04-05) : lit les obs T2M/HU2M
  sur la vraie grille Drôme, atteint l'étage d'analyse.
- **Blocage analyse** — les deux méthodes butent sur la **même contrainte** :
  - **OI** (`CASSIM_ISBA='OI'`) : first-guess **FA** opérationnel MF (précip/nébulosité/flux CANARI).
  - **EKF/SEKF** (`CASSIM_ISBA='EKF'`) : `PREP_INIT.nc` + prévisions **perturbées** → il faut
    faire tourner OFFLINE sur la période de gel.
  - Or **pas de forcing pour les saisons de gel** : AROME public = temps-réel seul,
    AROME **archive** = verrou MF (#11), CERRA S3 = t2m-only.

**Donc #7 (comme #5) dépend de l'archive AROME (MF) pour l'historique.** Le socle
technique (SURFEX + SODA + ingestion Sencrop) est complet et validé ; la validation
sur événements de gel réels attend le forcing archive.

## ✅ Baseline CERRA vs Sencrop (#16/#5) — validé sur gel réel

CERRA `reanalysis-cerra-single-levels` est **public** (Copernicus CDS, pas de MF) →
débloque les saisons de gel historiques (`cerra_to_forcing.py`).

**Run 2023-04-04→06 (épisode de gel Drôme), forcing CERRA → SURFEX 1 km :**

| @ 2023-04-05 04h | Sencrop (air 2 m) | T_skin SURFEX (CERRA) |
|---|---|---|
| Plage | −2,8 → 12,5 °C | −0,8 → 3,7 °C |
| Stations en gel (<0°C) | 9 / 37 | **12 / 37** |
| — | biais T_skin−air = **−1,37 °C** · RMSE 2,49 °C · **r = 0,73** | T_skin min domaine −4,2 °C |

Le T_skin descend plus bas que l'air — **le signal de gel radiatif**, exactement ce
que Karpos vise (cf. `product.md` « CERRA-Land T_skin »). Résultat sur un **premier
run non calibré**, rayonnement clear-sky (forecast CERRA `ssrd`/`strd` en secours,
cf. ci-dessous).

**Limites / suites** : `T2M` (diag 2 m air) sort en fill — comparaison faite sur
`TSRAD` (T_skin, la variable opérante) ; rayonnement = clear-sky proxy (Brutsaert LW
+ solaire clair) en attendant le forecast CERRA ; un seul épisode (baseline saisons
#5 + métriques POD/FAR #6 à dérouler) ; pas de calibration.

## Fichiers

- `OPTIONS.nam` — PGD smoke. `OPTIONS_run.nam` — OFFLINE (AROME). `OPTIONS_cerra.nam` — OFFLINE (CERRA).
- `OPTIONS_soda.nam` — chaîne SODA (PGD+PREP+assimilation OI/EKF).
- `arome_to_forcing.py`, `dem_to_surfex.py` — converters forcing/DEM.
