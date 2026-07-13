!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!     ######################
      MODULE MODD_IO_SURF_OL
!     ######################
!
!!****  *MODD_IO_SURF_OL* Keep in memory the netcdf ID of the output files
!!
!!    PURPOSE
!!    -------
!
!!
!!**  IMPLICIT ARGUMENTS
!!    ------------------
!!      None 
!!
!!    REFERENCE
!!    ---------
!!
!!    AUTHOR
!!    ------
!!      F. Habets   *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      modified 04/04 by P. LeMoigne: add logical for town, sea and water
!!      modified 07/11 by B. Decharme: add mask for IGN grid
!!      modified 03/20 by B. Decharme: add NEND_ATM for forcing interpolation
!
!*       0.   DECLARATIONS
!             ------------
!
!
IMPLICIT NONE
!------------------------------------------------------------------------------
!
!* variables for each patch
!
INTEGER, DIMENSION(:),POINTER :: NMASK
!$OMP THREADPRIVATE(NMASK)
INTEGER                       :: NSTEP_OUTPUT
LOGICAL                       :: LMASK = .FALSE.
LOGICAL                       :: LPARTR,LPARTW
LOGICAL                       :: LTIME_WRITTEN
INTEGER                       :: XTYPE
INTEGER                       :: XSTART,XCOUNT,XSTRIDE
INTEGER                       :: XSTARTW,XCOUNTW
INTEGER                       :: NEND_ATM       
!
LOGICAL :: LDEF 
INTEGER :: NID_NC
LOGICAL :: LRESET_DIAG = .TRUE.  ! if true diagnostics are reinitialised
!
REAL                          :: XTIMEC   ! cumulated time since start of the run
INTEGER, DIMENSION(:),ALLOCATABLE :: NMASK_IGN

!------------------------------------------------------------------------------
!
END MODULE MODD_IO_SURF_OL

