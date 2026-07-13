!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!     ###########################
      MODULE MODD_FORC_OFFLINE_NC
!     ###########################
!
!!****  *MODD_FORC_OFFLINE_NC - declaration of forcing variables for ISBA-TOPMODEL
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
!!    B. Vincendon
!!
!!    MODIFICATIONS
!!    -------------
!!      Original       11/2014
!
!*       0.   DECLARATIONS
!             ------------
!
!
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------
INTEGER                         :: NNB_FORC_STP       ! Number of total forcing time steps
INTEGER                         :: NNB_FORC_SEQUENCES ! Number (up to 5) of sequences
!                                                     with different sources in
!                                                     the forcing serie
INTEGER,          DIMENSION(15) :: NSTP_BEG           ! Time step of a sequence beginning
INTEGER,          DIMENSION(15) :: NSTP_END           ! Time step of a sequence end
CHARACTER(LEN=6), DIMENSION(15) :: CTYPE_SEQUENCES    ! Type of forcing sequences
!                                                     can be 'MODEL ','NORAIN','RADAR ','SAFRAN','MODRAD'
CHARACTER(LEN=18),DIMENSION(15) :: CGRIB_BASE_NAME    ! Base name of grib file
CHARACTER(LEN=18),DIMENSION(15) :: CRAD_BASE_NAME     ! Base name of radar file
CHARACTER(LEN=6), DIMENSION(15) :: CGRIB_TYPE         ! Model source
REAL                            :: XSTEP_ISBA         ! ISBA time step
REAL                            :: XSTEP_FORCING      ! forcing time step
INTEGER                         :: NNB_MEMBERS_ENS    ! Number of members of the ensemble 
!***
! Variables for ideal cases
REAL,              DIMENSION(15) :: XTA_IDEA
REAL,              DIMENSION(15) :: XQA_IDEA
REAL,              DIMENSION(15) :: XDIRSW_IDEA
REAL,              DIMENSION(15) :: XSCASW_IDEA
REAL,              DIMENSION(15) :: XLW_IDEA
REAL,              DIMENSION(15) :: XPS_IDEA
REAL,              DIMENSION(15) :: XRAIN_IDEA
REAL,              DIMENSION(15) :: XSNOW_IDEA
REAL,              DIMENSION(15) :: XWINDSPEED_IDEA
REAL,              DIMENSION(15) :: XWINDDIR_IDEA
!
!
!-------------------------------------------------------------------------------------
!
END MODULE MODD_FORC_OFFLINE_NC

