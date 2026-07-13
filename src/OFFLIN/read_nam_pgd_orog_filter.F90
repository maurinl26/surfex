!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!     #########
      SUBROUTINE READ_NAM_PGD_OROG_FILTER(HPROGRAM, KOPTFILTER, KZSFILTER, PCOFILTER, PTHFILTER)  
!     ##############################################################
!
!!**** *READ_NAM_PGD_OROG_FILTER* reads namelist for Orography
!!
!!    PURPOSE
!!    -------
!!
!!    METHOD
!!    ------
!!   
!     Remark about RTHFILTER
!      NOPTFILTER      = 0  : filtering is done everywhere
!      NOPTFILTER      = 1  : filtering is done at locations where orography is above a threshold
!
!!    EXTERNAL
!!    --------
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!    REFERENCE
!!    ---------
!!
!!    AUTHOR
!!    ------
!!
!!    B. Decharme        Meteo-France
!!
!!    MODIFICATION
!!    ------------
!!
!!    Original    02/2010
!!    09/2018 : Y. Seity : Add new filtering options
!----------------------------------------------------------------------------
!
!*    0.     DECLARATION
!            -----------
!
USE MODD_SURF_PAR,       ONLY : XUNDEF
!
USE MODI_GET_LUOUT
USE MODI_OPEN_NAMELIST
USE MODI_CLOSE_NAMELIST
!
USE MODE_POS_SURF
!
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*    0.1    Declaration of arguments
!            ------------------------
!                                   
CHARACTER(LEN=6),    INTENT(IN)    :: HPROGRAM    ! Type of program
INTEGER,             INTENT(OUT)   :: KOPTFILTER  ! Filtering option
INTEGER,             INTENT(OUT)   :: KZSFILTER   ! number of orographic spatial filter iterations                     
REAL(KIND=JPRB),     INTENT(OUT)   :: PCOFILTER   ! Filtering coefficient
REAL(KIND=JPRB),     INTENT(OUT)   :: PTHFILTER   ! Filtering threshold
!
!*    0.2    Declaration of local variables
!            ------------------------------
!
INTEGER                           :: ILUOUT    ! output listing logical unit
INTEGER                           :: ILUNAM    ! namelist file logical unit
LOGICAL                           :: GFOUND    ! flag when namelist is present
!
!*    0.3    Declaration of namelists
!            ------------------------
!
INTEGER                  :: NOPTFILTER  ! Filtering option
INTEGER                  :: NZSFILTER   ! number of orographic spatial filter iterations
REAL(KIND=JPRB)          :: RCOFILTER   ! Filtering coefficient
REAL(KIND=JPRB)          :: RTHFILTER   ! Filtering threshold

REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
NAMELIST/NAM_ZS_FILTER/NOPTFILTER,NZSFILTER,RCOFILTER,RTHFILTER
!
!-------------------------------------------------------------------------------
!
!*    1.      Initializations of defaults
!             ---------------------------
!
IF (LHOOK) CALL DR_HOOK('READ_NAM_PGD_OROG_FILTER',0,ZHOOK_HANDLE)

NOPTFILTER     = 0
NZSFILTER      = 1
RCOFILTER      = 1.
RTHFILTER      = 0.
!
 CALL GET_LUOUT(HPROGRAM,ILUOUT)
!
!-------------------------------------------------------------------------------
!
!*    2.      Reading of namelist
!             -------------------
!
 CALL OPEN_NAMELIST(HPROGRAM,ILUNAM)
!
 CALL POSNAM(ILUNAM,'NAM_ZS_FILTER',GFOUND,ILUOUT)
IF (GFOUND) READ(UNIT=ILUNAM,NML=NAM_ZS_FILTER)
!
 CALL CLOSE_NAMELIST(HPROGRAM,ILUNAM)
!
!-------------------------------------------------------------------------------
!
KOPTFILTER = NOPTFILTER
KZSFILTER = NZSFILTER ! number of orographic spatial filter iterations
PCOFILTER = RCOFILTER
PTHFILTER = RTHFILTER

IF (LHOOK) CALL DR_HOOK('READ_NAM_PGD_OROG_FILTER',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE READ_NAM_PGD_OROG_FILTER
