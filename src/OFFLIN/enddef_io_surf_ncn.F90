!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!     #######################################################
      SUBROUTINE ENDDEF_IO_SURF_NC_n()
!     #######################################################
!
!!****  *ENDDEF_IO_SURF_NC_n* - routine to finish definition mode of IO netcdf files
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
!!      M .Lafaysse   *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    03/2020 
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_SURFEX_MPI, ONLY : NRANK, NPIO
!
USE MODD_IO_SURF_NC, ONLY : NID_NCOUT
USE MODI_HANDLE_ERR
!
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
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!*       0.2   Declarations of local variables
!              -------------------------------
!
INTEGER :: IRET
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('ENDDEF_IO_SURF_NC_N',0,ZHOOK_HANDLE)
!
IF (NRANK==NPIO) THEN
    IRET = NF90_ENDDEF(NID_NCOUT)
    IF (IRET.NE.NF90_NOERR) CALL HANDLE_ERR(IRET,'ENDDEF_IO_SURF_NC_n')
ENDIF
!
IF (LHOOK) CALL DR_HOOK('ENDDEF_IO_SURF_NC_N',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE ENDDEF_IO_SURF_NC_n
