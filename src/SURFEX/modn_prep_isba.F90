!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!     ##################
      MODULE MODN_PREP_ISBA
!     ##################
!
!!****  *MODN_PREP_ISBA* - declaration of namelist NAM_PREP_ISBA
!!
!!    PURPOSE
!!    -------
!       The purpose of this module is to specify  the namelist NAM_PREP_ISBA
!
!!
!!**  IMPLICIT ARGUMENTS
!!    ------------------
!!
!!    REFERENCE
!!    ---------
!!
!!       
!!    AUTHOR
!!    ------
!!      V. Masson    *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original     01/2004                   
!!      P Le Moigne  08/2023  Initialization of soil temperature and moisture profiles
!-------------------------------------------------------------------------------
!
!*       0.   DECLARATIONS
!             ------------
!
USE MODD_PREP_ISBA,  ONLY : CFILE_ISBA, CTYPE, CFILEPGD_ISBA, CTYPEPGD,          &
                            CFILE_HUG, CTYPE_HUG, CFILE_TG, CTYPE_TG,            &
                            XUNIF_HUG_SOIL, XUNIF_HUGI_SOIL, XUNIF_TG_SOIL,      &
                            LEXTRAP_TG,LEXTRAP_WG, LEXTRAP_WGI,LEXTRAP_SN,       &
                            CFILE_HUG_SOIL, CFILE_TG_SOIL

!
IMPLICIT NONE
INTEGER           :: NYEAR        ! YEAR for surface
INTEGER           :: NMONTH       ! MONTH for surface
INTEGER           :: NDAY         ! DAY for surface
REAL              :: XTIME        ! TIME for surface
LOGICAL           :: LISBA_CANOPY !flag to use air layers inside the canopy
!
NAMELIST/NAM_PREP_ISBA/CFILE_ISBA, CTYPE, CFILEPGD_ISBA, CTYPEPGD,         &
                       CFILE_HUG, CTYPE_HUG, CFILE_TG, CTYPE_TG,           &
                       XUNIF_HUG_SOIL, XUNIF_HUGI_SOIL, XUNIF_TG_SOIL,     &
                       NYEAR, NMONTH, NDAY, XTIME, LISBA_CANOPY,LEXTRAP_TG,&
                       LEXTRAP_WG,LEXTRAP_WGI,LEXTRAP_SN,                  &
                       CFILE_HUG_SOIL, CFILE_TG_SOIL
  
!
END MODULE MODN_PREP_ISBA
