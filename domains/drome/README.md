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

**Run 2023-04-04→06 (épisode de gel Drôme), forcing CERRA (rayonnement réel) → SURFEX 1 km.**
Comparaison @2023-04-05 04h vs 37 stations Sencrop :

| Variable modèle | biais | RMSE | r | gel <0°C (modèle / obs) |
|---|---|---|---|---|
| **T2m air** (air-vs-air) | **−0,54 °C** | 2,14 °C | **0,72** | 12 / 9 |
| **T_skin** (TSRAD) | −2,12 °C | 3,01 °C | 0,72 | **18** / 9 |

- **T2m quasi non biaisé** (−0,54°C) vs le vrai air Sencrop — sur un run **non calibré**.
- **T_skin ~1,6°C plus froid que l'air** → détecte **18/37 stations en gel** (vs 9 en air) :
  le signal de gel radiatif que le T2m rate, **exactement ce que Karpos vise**
  (`product.md` « CERRA-Land T_skin »).

Détails techniques :
- Forcing = **CERRA réel** : analyse 3-horaire (t2m/r2/si10/wdir10/sp) + forecast leadtimes
  **1-6** (ssrd/strd/tp), déaccumulés par init, interpolés sur la grille PGD 1 km.
- `T2M` diag : fix **`N2M=2`** (méthode Paulson `CLS_TQ`) — sinon le 2 m sort en fill.

**Suites** : dérouler la **saison Feb–Avr** (#5) ; **métriques POD/FAR/CSI** (#6) ; boucler
**SODA** forcé CERRA (#7). Split direct/diffus SWdown encore simple (0,85/0,15).

## Fichiers

- `OPTIONS.nam` — PGD smoke. `OPTIONS_run.nam` — OFFLINE (AROME). `OPTIONS_cerra.nam` — OFFLINE (CERRA).
- `OPTIONS_soda.nam` — chaîne SODA (PGD+PREP+assimilation OI/EKF).
- `arome_to_forcing.py`, `dem_to_surfex.py` — converters forcing/DEM.
