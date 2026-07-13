!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!     #######################################################
      SUBROUTINE CLOSE_IO_SURF_NC_n()
!     #######################################################
!
!!****  *CLOSE_AUX_IO_SURF* - chooses the routine to OPENialize IO
!!
!!    PURPOSE
!!    -------
!!
!!**  METHOD
!!    ------
!!
!!    EXTERNAL
!!    --------
!!
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!    REFERENCE
!!    ---------
!!
!!
!!    AUTHOR
!!    ------
!!      V. Masson   *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODI_HANDLE_ERR
!
USE MODD_SURFEX_MPI, ONLY : NRANK, NPIO
!
USE MODD_IO_SURF_NC, ONLY: NID_NCOUT
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
USE NETCDF
!
IMPLICIT NONE
!
!
!*       0.1   Declarations of arguments
!              -------------------------
!
!*       0.2   Declarations of local variables
!              -------------------------------
!
INTEGER :: IRET ! return code
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('CLOSE_IO_SURF_NC',0,ZHOOK_HANDLE)
!
IF (NRANK == NPIO) THEN
  IRET = NF90_CLOSE(NID_NCOUT)
  IF (IRET.NE.NF90_NOERR) CALL HANDLE_ERR(IRET,'CLOSE_IO_SURF_NC_n')
ENDIF
!
IF (LHOOK) CALL DR_HOOK('CLOSE_IO_SURF_NC',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE CLOSE_IO_SURF_NC_n
