# open-SURFEX V9.1.0 — surface physique pour Karpos

Distribution open-source de **SURFEX** (SURFace EXternalisée, CNRM / Météo-France), packagée
comme **modèle de surface** du moteur de détection de gel [Karpos](https://github.com/maurinl26/karpos-engine).

L'objectif de ce dépôt n'est pas de redistribuer SURFEX pour lui-même, mais d'en faire une
**brique reproductible** : SURFEX tourne en mode `OFFLINE` (ISBA forcé hors-ligne), forcé par
réanalyse (CERRA 5,5 km, puis AROME 1,3 km), et **assimile les capteurs Sencrop in-verger**
pour produire une T_min nocturne 1 km physiquement cohérente sur la Drôme-Ardèche.

> ⚠️ **Statut : cap scientifique, non lancé.** Ce travail correspond à l'EPIC D de Karpos
> ([`maurinl26/karpos-downscaling#30`](https://github.com/maurinl26/karpos-downscaling/issues/30)),
> bloqué par une **Phase 0 non négociable** (méthodologie signée + partenariat Météo-France).
> Ne rien lancer côté build/run avant clôture de la Phase 0.

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
| `src/` | Sources Fortran SURFEX (`SURFEX/`, `OFFLIN/`, `ASSIM/`, `FORC/`, `LIB/`) + `Makefile` |
| `src/configure` | Génère l'environnement de build (`VERSION_MASTER=SFX-V9-1`) |
| `src/Rules.*.mk` | Règles de compilation par plateforme (gfortran, ifort, pgi…) |
| `conf/` | Fichiers de configuration compilateur (`conf_LXg95`, `profile_surfex.ihm`…) |
| `bin/` | Outils auxiliaires (conversion LFI/NetCDF, préprocesseurs Perl) |
| `MY_RUN/` | Cas-tests OFFLINE (namelists, forcings, ECOCLIMAP) — `hapex`, `cdp9697`, `me93`… |
| `exe/` | Cible des exécutables compilés (`PREP`, `PGD`, `OFFLINE`) — vide au checkout |
| `Licence_CeCILL-C_V1-*.txt` | Licence **CeCILL-C** (compatible exploitation commerciale) |

> Version : **SURFEX V9.1.0** (`NVERSION = 9`, `SFX-V9-1-0`). La roadmap Karpos mentionnait v8.1/8.2 ;
> la V9.1.0 convient pour le cap EPIC D.

---

## Build (aperçu)

Build natif recommandé sous **Linux (Ubuntu 22.04 LTS)** ; Apple Silicon en best-effort.
Nécessite un compilateur Fortran (gfortran) et NetCDF-Fortran.

```bash
cd src
export SFX_LIB=$PWD/../exe        # cible des binaires
./configure                       # génère le profil de build
# éditer le profil / choisir la Rules.*.mk adaptée (ex. Rules.LXgfortran.mk)
source ../conf/profile_surfex.<host>
make -j          # compile PGD, PREP, OFFLINE
```

Les exécutables `PGD`, `PREP`, `OFFLINE` atterrissent dans `exe/`.

> 🎯 **À faire (Phase 1, EPIC D)** : figer une chaîne reproductible → `Dockerfile` `karpos-surfex`
> (image privée GHCR), Ubuntu 22.04 + gfortran + NetCDF, tag de version épinglé.

---

## Cas-test OFFLINE

Les cas de `MY_RUN/` valident la chaîne PGD → PREP → OFFLINE sur des forcings 1D/petit domaine :

```bash
cd MY_RUN/KTEST/<cas>     # ex. hapex
# adapter les namelists (OPTIONS.nam, PRE_*.nam) et FORCING
../../../exe/PGD
../../../exe/PREP
../../../exe/OFFLINE
```

Ils servent de **smoke-test** avant de passer au domaine Drôme-Ardèche réel.

---

## Chaîne cible Karpos

```
Forcings (CERRA 5,5 km → AROME 1,3 km) :
  rayonnement net, vent, humidité, pression, précip
            ↓
SURFEX standalone OFFLINE (ISBA + FLAKE + TEB)
  ECOCLIMAP-SG land cover + DEM (SRTM 30 m / IGN BD ALTI) agrégés à 1 km
  → T_surface, T_air 2 m à 1 km, bilan d'énergie cohérent
            ↓
Nudging Sencrop (IAU / 2D-VAR simplifié)  ← src/ASSIM
  Pousse l'état SURFEX vers les obs station, correction étalée dans le temps
            ↓
T2m_min nocturne 1 km, physique-consistante → API /v1/forecast/surfex/baronnies
```

**Domaine** : Drôme-Ardèche, bbox 44–45,5°N × 4–5,5°E, résolution cible **1 km**.
**Schémas** : ISBA-DIF (physique) vs Force-Restore (rapide) — arbitrage Phase 0 ; TEB (urbain local),
FLAKE (lac de Serre-Ponçon).

---

## Roadmap (EPIC D — `karpos-downscaling#30`)

- [ ] **Phase 0 — Verrou** : méthodologie signée + **partenariat Météo-France** (convention, mentor CNRM/GAME, licences ECOCLIMAP-SG & AROME archive). *Sans Phase 0, on ne lance pas la Phase 1.*
- [ ] **Phase 1 — Setup** : récup code CNRM ([git.umr-cnrm.fr](https://git.umr-cnrm.fr/git/Surfex_Git2.git)), build natif, `Dockerfile` reproductible.
- [ ] **Phase 2 — Forcings + domaine** : ECOCLIMAP-SG + DEM Drôme, PREP_PGD (invariants) + PREP_OFFLINE (init saison gel).
- [ ] **Phase 3 — Baseline forward** : run Feb–Avr 2023 sans nudging ; biais résiduels par régime ; ≥ Lot B.
- [ ] **Phase 4 — Nudging Sencrop** : IAU/2D-VAR, calibration du rappel + portée spatiale ; saisons 2024–2025.
- [ ] **Phase 5 — Validation hold-out** : test pur 2026, POD ≥ 90 % / FAR ≤ 20 % stratifiés, rapport opposable (DEP).
- [ ] **Phase 6 — Production** : orchestration nocturne, hébergement (Scaleway CPU / Mac mini), endpoint API PWA.

**Hors-scope** : AROME-Carmel, Méso-NH, couplage atmosphère interactif (SURFEX reste **forced offline**).

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
