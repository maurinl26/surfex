# Build & déploiement sur HPC classique

Guide pour construire `karpos-surfex` (SURFEX OFFLINE + SODA + wrapper Python)
sur un cluster HPC classique piloté par modules d'environnement (Lmod / Environment
Modules), sans droits admin ni Homebrew.

## 1. Récupérer la release

```bash
tar xzf karpos-surfex-<version>-src.tar.gz
cd karpos-surfex-<version>
```

Le tarball est **auto-suffisant** (sources SURFEX + package Python + `build.sh`).
L'archive XIOS de 85 Mo n'y est pas (inutile : build `OFFLINE` sans serveur d'I/O).

## 2. Charger la toolchain (modules)

`build.sh` détecte le NetCDF via `nf-config` : il suffit qu'un module NetCDF-Fortran
soit chargé. Deux variantes courantes :

**GNU (gfortran)**
```bash
module load gcc/12          # gfortran >= 10 requis
module load netcdf-fortran  # doit fournir nf-config dans le PATH
module load python/3.12
```

**Intel (ifort)** — fréquent sur HPC
```bash
module load intel           # ifort + icc
module load netcdf-fortran  # build Intel
module load python/3.12
```

Vérifier : `which gfortran nf-config` (ou `ifort`) et `nf-config --flibs`.

## 3. Construire la bibliothèque SURFEX

```bash
./build.sh --compiler gfortran      # GNU
./build.sh --compiler ifort         # Intel (Linux uniquement)
./build.sh --compiler gfortran --jobs 16   # paralléliser
```

Produit :
- `exe/libsurfex.a` — bibliothèque (dépendance du package Python) ;
- `exe/{PGD,PREP,OFFLINE,SODA}-*` — exécutables des programmes maîtres.

> Le build lie le NetCDF **du module** (pas de recompilation des tarballs bundlés).
> Un stub `grib_api` évite toute dépendance eccodes (forcings OFFLINE = NetCDF).

## 4. Installer le package Python

```bash
python -m venv ~/venvs/surfex && source ~/venvs/surfex/bin/activate
pip install ./python -Csetup-args=-Dsurfex_lib=$PWD/exe/libsurfex.a
```

### Nœuds sans accès Internet (air-gapped)

`pip install` utilise l'isolation de build → télécharge `meson-python`, `cython`,
`numpy`, `ninja`, `meson`. Sur un nœud sans réseau :

```bash
# 1) sur le nœud de connexion (avec réseau) : pré-installer les outils de build
pip install meson-python cython numpy ninja meson

# 2) puis build sans isolation
pip install ./python --no-build-isolation \
  -Csetup-args=-Dsurfex_lib=$PWD/exe/libsurfex.a
```

## 5. Vérifier

```bash
python -c "import karpos_surfex as ks; ks.init(); ks.finalize(); print('OK', ks.__version__)"
```

## 6. Exécution — parallélisme

Le build est **NOMPI** (SURFEX OFFLINE mono-nœud, forcé hors-ligne). Pour couvrir un
domaine, on parallélise par **job array** (une tuile / un point par tâche) plutôt que
par MPI — cohérent avec la chaîne nocturne Karpos.

```bash
#SBATCH --array=0-99
#SBATCH --cpus-per-task=4        # OpenMP intra-tâche (build compilé -fopenmp)
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
python run_tile.py --tile $SLURM_ARRAY_TASK_ID
```

## Dépannage

| Symptôme | Cause / correction |
|---|---|
| `nf-config: command not found` | module NetCDF-Fortran non chargé |
| `library 'netcdf' not found` au link | lib C NetCDF hors PATH : charger le module `netcdf` (C) en plus de `netcdf-fortran` |
| `ifort: command not found` | `module load intel` ; ifort = Linux uniquement |
| erreurs d'arguments (gfortran ≥10) | déjà géré (`-fallow-argument-mismatch` dans les Rules) |
| build lent | `--jobs N` ; le stub grib évite de compiler eccodes |
