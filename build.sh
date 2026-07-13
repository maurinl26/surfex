#!/usr/bin/env bash
#
# build.sh — Build portable de SURFEX OFFLINE pour Karpos.
#
# Objectif : contourner la détection par hostname du `src/configure` d'origine
# (qui ne reconnaît que les machines CNRM) et fournir un point d'entrée unique,
# reproductible, macOS + Linux, gfortran + flang, lié au NetCDF système.
#
# Stratégie « hybride staged » : ce wrapper pilote le build natif .mk (éprouvé)
# pour produire les binaires OFFLINE/PGD/PREP. Un CMakeLists.txt mince viendra
# ensuite envelopper la lib statique pour le packaging Python (scikit-build).
#
# Usage :
#   ./build.sh [--compiler gfortran|flang] [--debug] [--jobs N] [--clean]
#
set -euo pipefail

# --- Défauts -----------------------------------------------------------------
COMPILER=gfortran
OPTLEVEL=O2
JOBS="$( (command -v nproc >/dev/null && nproc) || sysctl -n hw.ncpu 2>/dev/null || echo 4)"
DO_CLEAN=0

# --- Parsing arguments -------------------------------------------------------
while [ $# -gt 0 ]; do
  case "$1" in
    --compiler) COMPILER="$2"; shift 2 ;;
    --debug)    OPTLEVEL=DEBUG; shift ;;
    --jobs)     JOBS="$2"; shift 2 ;;
    --clean)    DO_CLEAN=1; shift ;;
    -h|--help)  grep '^#' "$0" | sed 's/^#\s\?//'; exit 0 ;;
    *) echo "Argument inconnu : $1" >&2; exit 2 ;;
  esac
done

ROOT="$(cd "$(dirname "$0")" && pwd)"

# --- Détection OS ------------------------------------------------------------
case "$(uname -s)" in
  Darwin) OSPFX=MC ;;
  Linux)  OSPFX=LX ;;
  *) echo "OS non supporté : $(uname -s) (macOS/Linux uniquement)" >&2; exit 1 ;;
esac

# --- Sélection compilateur / ARCH -------------------------------------------
case "$COMPILER" in
  gfortran)
    command -v gfortran >/dev/null || { echo "gfortran introuvable dans le PATH" >&2; exit 1; }
    ;;
  flang)
    command -v flang >/dev/null || command -v flang-new >/dev/null || {
      echo "flang introuvable dans le PATH" >&2; exit 1; }
    echo "⚠️  Profil flang NON VALIDÉ (scaffold) — repli conseillé : --compiler gfortran"
    ;;
  *) echo "Compilateur non supporté : $COMPILER (gfortran|flang)" >&2; exit 1 ;;
esac
export ARCH="${OSPFX}${COMPILER}"   # MCgfortran, LXgfortran, MCflang, LXflang

if [ ! -f "$ROOT/src/Rules.${ARCH}.mk" ]; then
  echo "Profil de règles manquant : src/Rules.${ARCH}.mk" >&2; exit 1
fi

# --- Localisation du NetCDF système -----------------------------------------
# On lie le NetCDF de l'environnement (Homebrew sur macOS, apt/module sur Linux)
# plutôt que de reconstruire les tarballs bundlés (VER_CDF=CDFAUTO, lent).
if command -v nf-config >/dev/null 2>&1; then
  INC_NETCDF="$(nf-config --fflags)"
  LIB_NETCDF="$(nf-config --flibs)"
  # nc-config complète les chemins de la lib C si distincts (cas Homebrew)
  if command -v nc-config >/dev/null 2>&1; then
    INC_NETCDF="$INC_NETCDF $(nc-config --cflags)"
    LIB_NETCDF="$LIB_NETCDF $(nc-config --libs)"
  fi
else
  echo "nf-config introuvable — installe netcdf-fortran (brew install netcdf-fortran" >&2
  echo "  ou apt-get install libnetcdff-dev)." >&2
  exit 1
fi
export INC_NETCDF LIB_NETCDF

# --- Variables de configuration SURFEX --------------------------------------
export VER_MPI=NOMPI       # OFFLINE mono-nœud : pas de MPI
export VER_OMP=OMP         # OpenMP activé
export VER_XIOS=0          # pas de serveur d'I/O XIOS
export VER_CDF=CDFKARPOS   # valeur neutre : INC/LIB fournis en ligne de commande
export VER_GRIBAPI=NONE    # pas de GRIB (forcings OFFLINE en NetCDF)
export VER_ECCODES=NONE    # idem eccodes
export OPTLEVEL
export NEED_TOOLS=NO NEED_NCARG=NO MVWORK=NO

echo "──────────────────────────────────────────────────────────"
echo " SURFEX build — ARCH=$ARCH  OPT=$OPTLEVEL  jobs=$JOBS"
echo " NetCDF inc : $INC_NETCDF"
echo " NetCDF lib : $LIB_NETCDF"
echo "──────────────────────────────────────────────────────────"

cd "$ROOT/src"

# --- Configure (génère conf/profile_surfex-<XYZ>) ---------------------------
./configure

# Source le profil généré le plus récent (exporte SRC_SURFEX, chemins OBJDIR, PATH bin/)
PROFILE="$(ls -t "$ROOT"/conf/profile_surfex-* 2>/dev/null | head -1)"
[ -n "$PROFILE" ] || { echo "Profil non généré par configure" >&2; exit 1; }
echo "→ source $PROFILE"
# Le profil CNRM n'est pas écrit pour le mode strict (SRC_SURFEX non initialisée,
# dernière instruction à rc≠0) : on relâche errexit+nounset le temps du sourcing.
set +eu
# shellcheck disable=SC1090
. "$PROFILE"
set -eu

# --- Compilation (INC/LIB en ligne de commande = surcharge toute règle) -----
MK=(make -j"$JOBS" INC_NETCDF="$INC_NETCDF" LIB_NETCDF="$LIB_NETCDF")
[ "$DO_CLEAN" = 1 ] && "${MK[@]}" clean || true

"${MK[@]}"
"${MK[@]}" installmaster

echo "──────────────────────────────────────────────────────────"
echo " Build terminé. Binaires dans exe/ :"
ls -1 "$ROOT"/exe/ 2>/dev/null | grep -Ev '^empty.txt$' || echo "  (aucun binaire — voir le log)"
