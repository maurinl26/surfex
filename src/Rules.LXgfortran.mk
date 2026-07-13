#SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
#SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
#SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
#SFX_LIC for details. version 1.
##########################################################
#                                                        #
# Compiler Options                                       #
#                                                        #
##########################################################
#
#   Gfortran version
GFV=$(shell  gfortran --version | grep -E -m1 -o ' [[:digit:]\.]{2,}( |$$)'  | sed 's/\.//g' )
# use splr.pl script for dependency generation
USE_SPLR = YES
#
#OBJDIR_PATH=/home/escj/azertyuiopqsdfghjklm/wxcvbn/azertyuiopqsdfghjklmwxcvbn
#
# -fallow-argument-mismatch / -fallow-invalid-boz : requis pour compiler SURFEX V9
# avec gfortran >= 10 (appels type-MPI non conformes, désormais erreurs par défaut).
OPT_BASE  = -fdefault-real-8 -fdefault-double-8 -g -fno-second-underscore -fpic  -ffpe-trap=overflow,zero,invalid -fbacktrace -fconvert=swap -fallow-argument-mismatch -fallow-invalid-boz
#
OPT_PERF0 = -O0
OPT_PERF2 = -O2
OPT_CHECK = -fcheck=bounds,do,mem,pointer,recursion -finit-real=nan
OPT_I8    = -fdefault-integer-8
#
#
# Integer 4/8 option
#
#MNH_INT   ?=I4
#RJ LFI_RECL  ?=512
#
ifeq "$(MNH_INT)" "I8"
OPT_BASE         += $(OPT_I8)
#RJ LFI_INT           ?=8
MNH_MPI_RANK_KIND ?=8
else
MNH_MPI_RANK_KIND ?=4
#RJ LFI_INT           ?=4
endif
#
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
F90= gfortran -pg
CC = gcc
else
F90 = mpif90
CC  = mpicc
endif
#
REALFC=gfortran
#
FCFLAGS_OMP= -fopenmp
CFLAGS_OMP= -fopenmp
ifeq "$(VER_OMP)" "NOOMP"
FCFLAGS_OMP= -pg
CFLAGS_OMP=
endif
#
F90FLAGS      = $(FCFLAGS_OMP) $(OPT)
F77 = $(F90)
F77FLAGS      = $(FCFLAGS_OMP) $(OPT)
FX90 = $(F90)
FX90FLAGS     = $(FCFLAGS_OMP) $(OPT)
#
LDFLAGS   =  $(FCFLAGS_OMP) -Wl,-warn-once -ldl -lrt
#
# preprocessing flags
#
CPP = cpp -P -traditional -Wcomment
#
FPPFLAGS_SURFEX    =
#RJ FPPFLAGS_SURCOUCHE = -DMNH_MPI_DOUBLE_PRECISION -DMNH_LINUX -DMNH_MPI_BSEND -DDEV_NULL  -DMNH_MPI_RANK_KIND=$(MNH_MPI_RANK_KIND)
#RJ FPPFLAGS_RAD       =
#RJ FPPFLAGS_NEWLFI    = -DSWAPIO -DLINUX -DLFI_INT=${LFI_INT} -DLFI_RECL=${LFI_RECL}
#RJ FPPFLAGS_MNH       = -DMNH -DAINT=INT -DAMOD=MOD
#
#
# ecCodes or grib_api selection
#SFX_GRIBAPI: if set to no:  use ecCodes
#             if set to yes: use grib_api (deprecated library)
#
SFX_GRIBAPI=no
#
# Gribex flags
#
TARGET_GRIBEX=linux
CNAME_GRIBEX=_gfortran
#
# Force -fallow-argument-mismatch option for gcc >= 10.1
# Necessary because some subroutines may be called with different datatypes
# Known list: MPI_Allgatherv,MPI_Allreduce,MPI_Bcast,MPI_Bsend,MPI_Gather,MPI_Gatherv,MPI_Recv,LEPOLY,EXTRACT_BBUFF,FILL_BBUFF
# + ecCodes + netCDF-fortran < 4.5.3
#
ifeq ($(shell test $(GFV) -ge 1010 ; echo $$?),0)
OPT_BASE += -fallow-argument-mismatch
GRIB_FLAGS += -fallow-argument-mismatch
#NETCDF_SUPPFLAGS += -fallow-argument-mismatch
ECCODES_FFLAGS += -fallow-argument-mismatch
endif
#
##########################################################
#                                                        #
# Source of MESONH PACKAGE  Distribution                 #
#                                                        #
##########################################################
#
include Makefile.SURFEX.mk
#
##########################################################
#                                                        #
# extra VPATH, Compilation flag modification             #
#         systeme module , etc ...                       #
#         external precompiled module librairie          #
#         etc ...                                        #
#                                                        #
##########################################################

#RJ ifneq "$(findstring 8,$(LFI_INT))" ""
#RJ OBJS_I8=spll_NEWLFI_ALL.o
#RJ $(OBJS_I8) : OPT = $(OPT_BASE) $(OPT_PERF2) $(OPT_I8)
#RJ endif
