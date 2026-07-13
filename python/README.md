# karpos-surfex

Package Python (CPU) exposant **SURFEX OFFLINE + assimilation de surface (SODA)**
pour Karpos. Backend de build : **meson-python** ; interop **Cython** via une
couche C-ABI Fortran (`iso_c_binding`).

> Statut : **scaffold**. L'ABI et l'API Python sont définitives et compilent ;
> le câblage profond du shim vers les drivers SURFEX est en cours (les appels
> renvoient `SurfexError: non implémenté` tant qu'ils ne sont pas branchés).

## Pourquoi cette architecture

SURFEX ne se compile pas par simple `cpp` : il génère 334+ modules d'interface
explicites (`MODI_*`) via perl au build. Réécrire ça en Meson serait lourd et
risqué. On adopte donc une **façade hybride** :

```
meson-python (façade, produit le wheel)
  ├─ capi/grib_api_stub.F90   → stub grib_api (build NetCDF-only, sans eccodes)
  ├─ capi/sfx_capi.F90        → ABI C stable (bind(C)) : OFFLINE + SODA + obs + champs
  ├─ karpos_surfex/_surfex.pyx→ Cython → numpy/xarray
  └─ libsurfex.a (optionnel)  → produit par ../build.sh (chaîne SURFEX éprouvée)
```

Le gros Fortran reste compilé par la chaîne SURFEX (`../build.sh`), qui sait déjà
générer les interfaces et gérer les flags par-fichier. Meson compile le shim +
Cython et **lie** `libsurfex.a`. Migration vers un build 100 % Meson-natif
possible plus tard si justifié.

## La surface exposée

| Python | Rôle |
|---|---|
| `init()` / `finalize()` | cycle de vie |
| `run_pgd(dir)` | physiographie (cover, sol, orographie) |
| `run_prep(dir)` | état initial (T_sol, humidité, neige) |
| `run_offline(dir)` | intégration forward (forcing NetCDF → T2m) |
| `run_soda(dir)` | **assimilation** (OI / SEKF / ENKF / PF) |
| `set_obs(lats, lons, vals)` | injecte les obs **Sencrop** pour l'analyse |
| `get_field(name)` | lit un champ (`T2M`, `WG1`, `TG1`…) → numpy |

**Le « nudging Sencrop » = `set_obs()` + `run_soda()`** : l'assimilation passe par
le moteur natif SURFEX (SODA), pas par un module ad hoc. Cf. `karpos-downscaling#30`.

## Build (uv — env géré, aligné sur karpos-engine)

Environnement Python géré par **uv** (Python épinglé dans `.python-version` = 3.11).

```bash
# 1) produire libsurfex.a via la chaîne SURFEX
cd .. && ./build.sh --compiler gfortran

# 2a) build du wheel
cd python
uv build --wheel --config-setting=setup-args=-Dsurfex_lib=$PWD/../exe/libsurfex.a

# 2b) ou installer en editable pour développer
uv venv
uv pip install -e . --config-setting=setup-args=-Dsurfex_lib=$PWD/../exe/libsurfex.a
uv run --with netCDF4 python -c "import karpos_surfex as ks; print(ks.__version__)"
```

Repli pip : `pip install . -Csetup-args=-Dsurfex_lib=…` (identique).

Prérequis système : `gfortran`, `netcdf-fortran` (`nf-config`). uv gère Python +
`meson-python`, `cython`, `numpy` automatiquement.

## Exemple

```python
import karpos_surfex as ks
import numpy as np

ks.init()
ks.run_offline("run/")                     # first-guess
ks.set_obs(lats, lons, sencrop_t2m)        # obs Sencrop
ks.run_soda("assim/")                      # analyse
t2m = ks.get_field("T2M")                  # grille analysée
ks.finalize()
```
