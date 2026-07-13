# open-SURFEX V9.1.0 — surface physique pour Karpos

Distribution open-source de **SURFEX** (SURFace EXternalisée, CNRM / Météo-France), packagée
comme **modèle de surface** du moteur de détection de gel [Karpos](https://github.com/maurinl26/karpos-engine).

L'objectif de ce dépôt n'est pas de redistribuer SURFEX pour lui-même, mais d'en faire une
**brique reproductible** : SURFEX tourne en mode `OFFLINE` (ISBA forcé hors-ligne), forcé par
**AROME 0,025°** (public), sur la grille **1 km Drôme-Ardèche**, avec orographie **IGN BD ALTI**,
et **assimile les capteurs Sencrop in-verger** (via SODA) pour produire une T_skin / T2m
nocturne physiquement cohérente.

> ✅ **Statut : brique CPU debout et validée.** La chaîne **PGD → PREP → OFFLINE tourne de bout
> en bout** sur données réelles et produit le **T_skin (TSRAD) 1 km** — voir [Bilan](#-bilan--état)
> ci-dessous. Correspond à l'EPIC D de Karpos
> ([`maurinl26/karpos-downscaling#30`](https://github.com/maurinl26/karpos-downscaling/issues/30)),
> voie « sans Météo-France ». Le cap *opposable* final reste conditionné à la Phase 0 (partenariat MF).

---

## 🏔️ Bilan / état

Point de départ : distribution SURFEX brute. Arrivée : **chaîne SURFEX 1 km opérationnelle
sur la Drôme, produisant le T_skin à partir de données AROME réelles, packagée en wheel
Python, avec DevOps complet.**

**Ce qui tourne aujourd'hui**
```
Forcing AROME 0,025° public (SP1/2/3, temps-réel, meteo.data.gouv.fr)
   ↓  domains/drome/arome_to_forcing.py   (conversions + interpolation → grille PGD)
FORCING.nc (20 160 pts, grille 1 km Drôme)
   ↓  PGD (orographie IGN BD ALTI réelle) + PREP (init uniforme) + OFFLINE (ISBA)
SURF_ATM_DIAGNOSTICS.OUT.nc  →  TS, TSRAD (T_skin), T2M, HU2M, RN, H, LE…  [120×168]
```
Lancement : `python -m karpos_surfex run --workdir <dir> --ecoclimap MY_RUN/ECOCLIMAP --steps pgd,prep,offline`

**État de la roadmap** ([issues](https://github.com/maurinl26/surfex/issues))

| Statut | Issues |
|---|---|
| ✅ **Fermées** | #1 build · #2 Docker GHCR · #3 domaine Drôme+IGN · #4 forcing AROME→OFFLINE · #13 CI · #14 orchestration · #15 packaging |
| 🟢 **Quasi** | #5 baseline (chaîne ✅, dérouler saisons) · #7 nudging Sencrop (obs grillées + first-guess ✅, reste le run SODA) |
| ⚪ À faire | #6 harness POD/FAR/CSI |
| 🔒 Verrou MF | #8–12 (Phase 0 : méthodo + partenariat Météo-France) |

**Livré**
- `build.sh` — build portable **macOS/Linux × gfortran/flang/ifort**, NetCDF système.
- Package **`karpos-surfex`** (meson-python + Cython + ABI C `iso_c_binding`), wheel CPU, **uv**.
- **Orchestration** `driver.py`/`orchestrate.py` + CLI `python -m karpos_surfex`.
- **Ingestion Sencrop** (S3 → obs → opérateur d'obs grillé pour SODA) — `python/karpos_surfex/sencrop.py`.
- **Domaine Drôme** (`domains/drome/`) : grille 1 km, DEM IGN, forcing AROME, namelists.
- **CI/CD** : GitHub Actions (build/test), **image GHCR** (`ghcr.io/maurinl26/surfex`), release HPC.
- **Débloquage PREP file-less** : patches SURFEX (repli init uniforme au lieu du GRIB) — voir `domains/drome/README.md`.

**Note GPU** : hors-scope. SURFEX/ISBA est un modèle physique CPU (parallélisme par point/tuile
→ OpenMP + job-array HPC, pas GPU). Le GPU concerne le downscaling DL (`karpos-downscaling`), pas
cette brique. Axe de montée en charge éventuel : **MPI multi-nœuds** (support présent, laissé `NOMPI`).

---

## Pourquoi SURFEX pour Karpos

Les baselines actuelles de Karpos butent sur le **gel radiatif de floraison** :

| Méthode | Régime advectif (hiver) | Régime radiatif (flo) |
|---|---|---|
| Lot B — QDM/RBF Sencrop sur CERRA | POD 86 % / FAR 14 % | POD 32 % / FAR 51 % |
| Lot C — U-Net FiLM | régresse vers la climatologie | sous-performe Lot B |

SURFEX est le **bon outil scientifique** pour passer ce cap :

- schéma de surface **physique** (ISBA — bilan d'énergie/eau) ;
- **cold-pools** et drainage froid topographique nativement modélisés ;
- **inversions nocturnes** via bilan radiatif ;
- le **nudging Sencrop** préserve la cohérence physique, là où un RBF spatial mélange les régimes ;
- **opposable** : SURFEX est l'opérationnel Météo-France → argument pour le certificat DEP.

Voir la note produit : [`karpos-engine/docs/product.md`](https://github.com/maurinl26/karpos-engine/blob/main/docs/product.md).

---

## Contenu du dépôt

| Répertoire | Rôle |
|---|---|
| **`build.sh`** | **Build portable** (macOS/Linux × gfortran/flang/ifort) — point d'entrée unique |
| **`python/`** | **Package `karpos-surfex`** (meson-python + Cython + ABI C) : driver, orchestration, ingestion Sencrop |
| **`domains/drome/`** | **Domaine Drôme 1 km** : namelists, DEM IGN, forcing AROME, SODA — voir son `README.md` |
| **`.github/workflows/`** | **CI/CD** : build+test, PyPI, Docker GHCR, release HPC, runner mac mini |
| **`Dockerfile`** | Image reproductible `ghcr.io/maurinl26/surfex` (multi-stage) |
| `src/` | Sources Fortran SURFEX (`SURFEX/`, `OFFLIN/`, `ASSIM/`, `FORC/`, `LIB/`) + `Makefile` |
| `src/Rules.*.mk` | Règles de compilation par plateforme (+ compat gfortran ≥10, profils flang) |
| `bin/` | Outils auxiliaires (conversion LFI/NetCDF, préprocesseurs Perl `splr.pl`…) |
| `MY_RUN/` | Cas-tests OFFLINE + covers ECOCLIMAP-II versionnés (`*.bin`) |
| `exe/` | Cible des exécutables compilés + `libsurfex.a` — vide au checkout |
| `Licence_CeCILL-C_V1-*.txt` | Licence **CeCILL-C** (compatible exploitation commerciale) |

> Version : **SURFEX V9.1.0** (`NVERSION = 9`, `SFX-V9-1-0`). La roadmap Karpos mentionnait v8.1/8.2 ;
> la V9.1.0 convient pour le cap EPIC D.

---

## Build & run

**Build** (portable, NetCDF système détecté via `nf-config`) :
```bash
./build.sh --compiler gfortran     # ou --compiler ifort (HPC) ; macOS/Linux
# → exe/{PGD,PREP,OFFLINE,SODA} + exe/libsurfex.a
```

**Package Python** (uv-managed, cf. `python/README.md`) :
```bash
cd python && uv build --wheel --config-setting=setup-args=-Dsurfex_lib=$PWD/../exe/libsurfex.a
```

**Chaîne complète sur le domaine Drôme** (forcing AROME public → T_skin) :
```bash
python domains/drome/arome_to_forcing.py <workdir> PGD.nc        # FORCING.nc AROME
export SURFEX_EXE_DIR=$PWD/exe
python -m karpos_surfex run --workdir <workdir> \
       --ecoclimap MY_RUN/ECOCLIMAP --steps pgd,prep,offline
# → SURF_ATM_DIAGNOSTICS.OUT.nc : TS, TSRAD (T_skin), T2M…  [120×168]
```

Détails domaine, DEM et assimilation : [`domains/drome/README.md`](domains/drome/README.md).
Déploiement HPC : [`docs/HPC.md`](docs/HPC.md). Image conteneur : `ghcr.io/maurinl26/surfex`.

---

## Chaîne cible Karpos

```
Forcing AROME 0,025° (public) :  rayonnement, vent, humidité, pression
            ↓  (converter + interpolation grille PGD)
SURFEX OFFLINE (ISBA) — orographie IGN BD ALTI à 1 km
  → TS, TSRAD (T_skin), T2M, HU2M, bilan d'énergie cohérent
            ↓
Assimilation Sencrop via SODA (OI / SEKF)  ← src/ASSIM
  obs stations grillées (opérateur d'obs) → état corrigé
            ↓
T_skin / T2m nocturne 1 km, physique-consistant → first-guess opposable
```

**Domaine** : Drôme-Ardèche, bbox 44–45,5°N × 4–5,5°E, **1 km** (120×168, CONF PROJ).

---

## Roadmap (EPIC D — `karpos-downscaling#30`)

Voie **« sans Météo-France »** (démarrable) — cf. [Bilan](#-bilan--état) :

- [x] **Build reproductible** (`build.sh`, CI, Docker GHCR) — #1/#2/#13
- [x] **Packaging CPU** (meson-python + Cython) — #15
- [x] **Domaine Drôme 1 km + orographie IGN** — #3
- [x] **Forcing AROME public → OFFLINE** — #4
- [x] **Chaîne forward PGD→PREP→OFFLINE** (T_skin produit) — #5 (baseline saisons à dérouler)
- [x] **Ingestion + opérateur d'obs Sencrop** — #7 (reste le run SODA)
- [ ] **Harness POD/FAR/CSI** stratifié — #6

Voie **« attente Météo-France »** (verrou Phase 0) — #8–12 : méthodologie signée + partenariat
MF (mentor CNRM/GAME, licences), forcings AROME **archive**, validation hold-out + rapport opposable (DEP).

**Hors-scope** : AROME-Carmel, Méso-NH, couplage atmosphère interactif (SURFEX reste **forced offline**) ; GPU (cf. Bilan).

---

## Licence

SURFEX est distribué sous **CeCILL-C V1** (voir `Licence_CeCILL-C_V1-en.txt` / `-fr.txt`),
compatible avec une exploitation commerciale sous réserve des obligations de la licence.
Vérifier séparément les licences **ECOCLIMAP-SG** et de l'**archive AROME** (Phase 0).

## Références

- Documentation SURFEX : <https://www.umr-cnrm.fr/surfex/>
- Source CNRM : <https://git.umr-cnrm.fr/git/Surfex_Git2.git>
- EPIC D SURFEX + nudging Sencrop : `maurinl26/karpos-downscaling#30`
- Note produit Karpos : `karpos-engine/docs/product.md`
