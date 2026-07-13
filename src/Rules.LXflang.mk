##########################################################
#                                                        #
# Compiler Options — Linux / LLVM flang                  #
#                                                        #
# ⚠️  SCAFFOLD NON VALIDÉ — flang n'était pas disponible  #
#     au moment de l'écriture. Flags à confirmer sur une  #
#     machine avec LLVM flang (>= 19 recommandé).         #
#     Piste de repli : gfortran (Rules.LXgfortran.mk).    #
##########################################################
#
# Génération des dépendances via splr.pl (comme les autres profils)
USE_SPLR = YES
#
# flang (LLVM). Ancien nom : flang-new.
# Promotion des réels en double (SURFEX est écrit en real*8 implicite).
OPT_BASE  = -fdefault-real-8 -fdefault-double-8 -fPIC
#
OPT_PERF0 = -O0
OPT_PERF2 = -O2
OPT_CHECK = -fcheck=all
OPT_I8    = -fdefault-integer-8
#
ifeq "$(MNH_INT)" "I8"
OPT_BASE         += $(OPT_I8)
MNH_MPI_RANK_KIND ?=8
else
MNH_MPI_RANK_KIND ?=4
endif
#
OPT       = $(OPT_BASE) $(OPT_PERF2)
OPT0      = $(OPT_BASE) $(OPT_PERF0)
OPT_NOCB  = $(OPT_BASE) $(OPT_PERF2)
#
ifeq "$(OPTLEVEL)" "DEBUG"
OPT       = $(OPT_BASE) $(OPT_PERF0) $(OPT_CHECK)
OPT0      = $(OPT_BASE) $(OPT_PERF0) $(OPT_CHECK)
OPT_NOCB  = $(OPT_BASE) $(OPT_PERF0)
endif
#
ifneq "$(OPTLEVEL)" "DEBUG"
OBJSD += spll_teb_garden.o
$(OBJSD) : OPT = $(OPT_BASE) $(OPT_PERF0)
endif
#
ifeq "$(VER_MPI)" "NOMPI"
F90= flang
CC = clang
else
F90 = mpif90
CC  = mpicc
endif
#
REALFC=flang
#
FCFLAGS_OMP= -fopenmp
CFLAGS_OMP= -fopenmp
ifeq "$(VER_OMP)" "NOOMP"
FCFLAGS_OMP=
CFLAGS_OMP=
endif
#
F90FLAGS      = $(FCFLAGS_OMP) $(OPT)
F77 = $(F90)
F77FLAGS      = $(FCFLAGS_OMP) $(OPT)
FX90 = $(F90)
FX90FLAGS     = $(FCFLAGS_OMP) $(OPT)
#
LDFLAGS   =  $(FCFLAGS_OMP)
#
# preprocessing flags
CPP = cpp -P -traditional -Wcomment
#
FPPFLAGS_SURFEX    =
#
TARGET_GRIBEX=linux
CNAME_GRIBEX=_flang
##########################################################
include Makefile.SURFEX.mk
