!----------------------------------------------------------------------------!
!     ##############################################################
      SUBROUTINE READ_NAM_PGD_TOPD(HPROGRAM,OCOUPL_TOPD,HCAT,PF_PARAM_BV,PC_DEPTH_RATIO_BV,&
                                   ODUMMY_SUBCAT,OSUBCAT,KSUBCAT,PLX,PLY,&
                                   HSUBCAT,HFILE_SUBCAT,OWRITE_SEVERITY_MAPS)
!     ##############################################################
!
!!**** *READ_NAM_TOPD_n* reads namelist NAM_TOPD
!!
!!    PURPOSE
!!    -------
!!    NAM_TOPD is a namelist used only for Topmodel coupling
!!    It permits to define the different catchments studied.
!!    This routine aims at reading and initialising those names.
!!
!!    METHOD
!!    ------
!!   
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
!!    B. Vincendon        Meteo-France
!!
!!    MODIFICATION
!!    ------------
!!
!!    Original   11/2006
!!
!----------------------------------------------------------------------------
!
!*    0.     DECLARATION
!            -----------
!
USE MODI_GET_LUOUT
USE MODI_OPEN_NAMELIST
USE MODI_CLOSE_NAMELIST
!
USE MODD_TOPD_PAR, ONLY : JPCAT
USE MODD_TOPODYN, ONLY : NNCAT
!
USE MODE_POS_SURF
!
USE YOMHOOK   ,ONLY : LHOOK, DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*    0.1    Declaration of arguments
!            ------------------------
!
 CHARACTER(LEN=6),                    INTENT(IN)  :: HPROGRAM     ! Type of program
 LOGICAL,                             INTENT(OUT) :: OCOUPL_TOPD
 CHARACTER(LEN=15), DIMENSION(JPCAT), INTENT(OUT) :: HCAT         ! Names of catchments         
 REAL,              DIMENSION(JPCAT), INTENT(OUT) :: PF_PARAM_BV
 REAL,              DIMENSION(JPCAT), INTENT(OUT) :: PC_DEPTH_RATIO_BV 
 LOGICAL,                             INTENT(OUT) :: ODUMMY_SUBCAT
 LOGICAL,                             INTENT(OUT) :: OSUBCAT
 INTEGER,           DIMENSION(JPCAT),      INTENT(OUT) :: KSUBCAT
 REAL,              DIMENSION(JPCAT,JPCAT),INTENT(OUT) :: PLX 
 REAL,              DIMENSION(JPCAT,JPCAT),INTENT(OUT) :: PLY
 CHARACTER(LEN=15), DIMENSION(JPCAT,JPCAT), INTENT(OUT) :: HSUBCAT         ! Names of catchments         
 CHARACTER(LEN=15), DIMENSION(JPCAT), INTENT(OUT) :: HFILE_SUBCAT!
 LOGICAL,                             INTENT(OUT) :: OWRITE_SEVERITY_MAPS
!
!*    0.2    Declaration of local variables
!            ------------------------------
!
CHARACTER(LEN=15), DIMENSION(JPCAT) :: CCAT
LOGICAL                           :: LCOUPL_TOPD, LDUMMY_SUBCAT
LOGICAL                           :: LSUBCAT, LWRITE_SEVERITY_MAPS
REAL, DIMENSION(JPCAT)            :: XF_PARAM_BV
REAL, DIMENSION(JPCAT)            :: XC_DEPTH_RATIO_BV
INTEGER, DIMENSION(JPCAT)         :: NSUBCAT
REAL, DIMENSION(JPCAT,JPCAT)      :: XLX,XLY
CHARACTER(LEN=15), DIMENSION(JPCAT,JPCAT) :: CSUBCAT
CHARACTER(LEN=15), DIMENSION(JPCAT) :: CFILE_SUBCAT
INTEGER                           :: ILUOUT    ! output listing logical unit
INTEGER                           :: ILUNAM    ! namelist file logical unit
LOGICAL                           :: GFOUND    ! flag when namelist is present
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!*    0.3    Declaration of namelists
!  
!
NAMELIST/NAM_PGD_TOPD/CCAT, LCOUPL_TOPD,XF_PARAM_BV, XC_DEPTH_RATIO_BV,  &
                      LDUMMY_SUBCAT,LSUBCAT,NSUBCAT,XLX,XLY,CFILE_SUBCAT,&
                      CSUBCAT,LWRITE_SEVERITY_MAPS
!-------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('READ_NAM_PGD_TOPD',0,ZHOOK_HANDLE)
!
!*    1.      Initializations of defaults
!             ---------------------------
!
LCOUPL_TOPD = .FALSE.
CCAT(:) = '   '
CSUBCAT(:,:) = '   '
CFILE_SUBCAT(:) = '   '
XF_PARAM_BV(:) = 2.5
XC_DEPTH_RATIO_BV(:) = 1.
LDUMMY_SUBCAT=.FALSE.
LSUBCAT=.FALSE.
NSUBCAT(:)=0
XLX(:,:)=0.
XLY(:,:)=0.
LWRITE_SEVERITY_MAPS=.FALSE.
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
 CALL POSNAM(ILUNAM,'NAM_PGD_TOPD',GFOUND,ILUOUT)
IF (GFOUND) READ(UNIT=ILUNAM,NML=NAM_PGD_TOPD)
!
 CALL CLOSE_NAMELIST(HPROGRAM,ILUNAM)
!
!         2.   Initialises number of catchments and time step variables
!              -------------------------------------------------------
!
NNCAT=COUNT(CCAT(:)/='   ')
!
!-------------------------------------------------------------------------------
!
!*    3.      Fills output arguments
!             ----------------------
!
OCOUPL_TOPD = LCOUPL_TOPD
HCAT(1:NNCAT) = CCAT(1:NNCAT)
HFILE_SUBCAT(:) = CFILE_SUBCAT(:)
HSUBCAT(:,:) = CSUBCAT(:,:)
PF_PARAM_BV(1:NNCAT) = XF_PARAM_BV(1:NNCAT)
PC_DEPTH_RATIO_BV(1:NNCAT) = XC_DEPTH_RATIO_BV(1:NNCAT)
ODUMMY_SUBCAT=LDUMMY_SUBCAT
OSUBCAT=LSUBCAT
KSUBCAT(:)=NSUBCAT(:)
PLX(:,:)=XLX(:,:)
PLY(:,:)=XLY(:,:)
OWRITE_SEVERITY_MAPS=LWRITE_SEVERITY_MAPS
!
IF (LHOOK) CALL DR_HOOK('READ_NAM_PGD_TOPD',1,ZHOOK_HANDLE)
!-------------------------------------------------------------------------------
!
END SUBROUTINE READ_NAM_PGD_TOPD
