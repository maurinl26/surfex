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

- ✅ **Géométrie de grille validée** : `PGD` produit `PGD.nc` (xx=120, yy=168).
  Smoke avec **cover uniforme prescrit** (`NAM_COVER XUNIF_COVER(4)`) et
  **orographie uniforme** (`XUNIF_ZS=300`) — valide la grille sans données externes.
- ⏳ **Physiographie réelle** : nécessite deux jeux de données à acquérir.

## Données à acquérir (physiographie réelle)

Pour un `PGD` physiquement correct, remplacer les valeurs uniformes par :

1. **Carte de cover ECOCLIMAP-II** (`LECOCLIMAP=T` + `NAM_COVER YCOVER=...`).
   Raster global (~Go). Source : Météo-France / CNRM (cf. licence #10 pour SG).
2. **DEM** pour l'orographie (`NAM_ZS YZS=...`) :
   - SRTM 30 m (OpenTopography / USGS) — libre, global.
   - IGN BD ALTI 25 m — libre, France, meilleure qualité locale.
   Agrégé à 1 km par les routines `average_orography` de PGD.

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

## Fichiers

- `OPTIONS.nam` — namelist PGD/PREP/OFFLINE du domaine.
