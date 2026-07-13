!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!----------------------------------------------------------------------------!
!     ##############################################################
       SUBROUTINE READ_NAM_TOPD(HPROGRAM,OBUDGET_TOPD,KNB_TOPD,OSTOCK_TOPD,&
                                KNB_STOCK,KNB_RESTART,KFREQ_MAPS_WG, &
                                KFREQ_MAPS_ASAT,KFREQ_MAPS_RUNOFF,PSPEEDR, &
                                PSPEEDG,PSPEEDH,OSPEEDR_VAR,PQINIT,PRTOP_D2,&
                                OPERT_PARAM, OPERT_INIT)
!     ##############################################################
!
!!**** *READ_NAM TOPD* reads namelist NAM_TOPD
!!
!!    PURPOSE
!!    -------
!!
!!    NAM_TOPD is a namelist used to define whether Topmodel coupling
!!    is performed or not and the time step ratio between hydrological 
!!    model and ISBA.
!!    This routine aims at reading and initialising those variables.
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
!!    Original    11/2006
!!    B. Vincendon 02/2014 : adding possibility to choose the speed of water on hillslopes and
!!                           to write runoff maps on watersheds
!!
!----------------------------------------------------------------------------
!
!*    0.     DECLARATION
!            -----------
!
USE MODD_TOPD_PAR, ONLY : JPCAT
USE MODD_TOPODYN, ONLY : NNCAT
!
USE MODE_POS_SURF
!
USE MODI_GET_LUOUT
USE MODI_OPEN_NAMELIST
USE MODI_CLOSE_NAMELIST
!
USE YOMHOOK   ,ONLY : LHOOK, DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*    0.1    Declaration of arguments
!            ------------------------
!
 CHARACTER(LEN=6),                   INTENT(IN)   :: HPROGRAM     ! Type of program
LOGICAL,                            INTENT(OUT)  :: OBUDGET_TOPD ! budget computation
INTEGER,                            INTENT(OUT)  :: KNB_TOPD     ! Ratio between Topmodel time step and ISBA time step
LOGICAL,                            INTENT(OUT)  :: OSTOCK_TOPD  ! T if use of stock from previous simulation
INTEGER,                            INTENT(OUT)  :: KNB_STOCK    ! number of time step to read in previous simulation
INTEGER,                            INTENT(OUT)  :: KNB_RESTART  ! number of time step to write for next simulation
INTEGER,                            INTENT(OUT)  :: KFREQ_MAPS_WG! 
INTEGER,                            INTENT(OUT)  :: KFREQ_MAPS_ASAT! 
INTEGER,                            INTENT(OUT)  :: KFREQ_MAPS_RUNOFF
REAL, DIMENSION(JPCAT),INTENT(OUT)               :: PSPEEDR ! River speed
REAL, DIMENSION(JPCAT),INTENT(OUT)               :: PSPEEDG ! Ground speed
REAL, DIMENSION(JPCAT),INTENT(OUT)               :: PSPEEDH ! Hillslope speed
LOGICAL,                            INTENT(OUT)  :: OSPEEDR_VAR  ! T to modulate river speed according to the discharge
REAL, DIMENSION(JPCAT),INTENT(OUT)               :: PQINIT  ! Initial discharge at catchments outlet
REAL, DIMENSION(JPCAT),INTENT(OUT)               :: PRTOP_D2
LOGICAL,                            INTENT(OUT)  :: OPERT_PARAM ! T to introduce perturbation in ISBA-TOP parameters
LOGICAL,                            INTENT(OUT)  :: OPERT_INIT  ! T to introduce perturbation in ISBA-TOP initial soil moisture
!
!*    0.2    Declaration of local variables
!            ------------------------------
!
LOGICAL                           :: LBUDGET_TOPD
LOGICAL                           :: LSTOCK_TOPD
INTEGER                           :: NNB_TOPD
INTEGER                           :: NFREQ_MAPS_WG
INTEGER                           :: NFREQ_MAPS_ASAT
INTEGER                           :: NFREQ_MAPS_RUNOFF
INTEGER                           :: NNB_STP_STOCK
INTEGER                           :: NNB_STP_RESTART
REAL, DIMENSION(JPCAT)            :: XSPEEDR
REAL, DIMENSION(JPCAT)            :: XSPEEDG
REAL, DIMENSION(JPCAT)            :: XSPEEDH
LOGICAL                           :: LSPEEDR_VAR
REAL, DIMENSION(JPCAT)            :: XQINIT
REAL, DIMENSION(JPCAT)            :: XRTOP_D2
LOGICAL                           :: LPERT_PARAM
LOGICAL                           :: LPERT_INIT
!
INTEGER                           :: ILUOUT    ! output listing logical unit
INTEGER                           :: ILUNAM    ! namelist file logical unit
LOGICAL                           :: GFOUND    ! flag when namelist is present
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!*    0.3    Declaration of namelists
!  
NAMELIST/NAM_TOPD/LBUDGET_TOPD, LSTOCK_TOPD, NNB_TOPD, &
                  NFREQ_MAPS_WG, NFREQ_MAPS_ASAT, NFREQ_MAPS_RUNOFF,&
                  NNB_STP_STOCK, NNB_STP_RESTART, &
                  XSPEEDR, XSPEEDG, XSPEEDH, LSPEEDR_VAR,&
                  XQINIT, XRTOP_D2,&
                  LPERT_PARAM, LPERT_INIT
!-------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('READ_NAM_TOPD',0,ZHOOK_HANDLE)
!
!*    1.      Initializations of defaults
!             ---------------------------
!
LBUDGET_TOPD = .FALSE.
LSTOCK_TOPD = .FALSE.
NNB_TOPD = 1
NFREQ_MAPS_WG = 0
NFREQ_MAPS_ASAT = 0
NFREQ_MAPS_RUNOFF = 0
NNB_STP_STOCK = 1
NNB_STP_RESTART = 1
XSPEEDR(:) = 1.0 ! default value of river speed 
XSPEEDG(:) = 0.1 ! default value of speed in the ground
XSPEEDH(:) = 0.1 ! default value of hillspeed
LSPEEDR_VAR = .FALSE.
XQINIT(:) = 0.
XRTOP_D2(:) = 1.
LPERT_PARAM = .FALSE.
LPERT_INIT  = .FALSE.
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
 CALL POSNAM(ILUNAM,'NAM_TOPD',GFOUND,ILUOUT)
IF (GFOUND) READ(UNIT=ILUNAM,NML=NAM_TOPD)
!
 CALL CLOSE_NAMELIST(HPROGRAM,ILUNAM)
!
!-------------------------------------------------------------------------------
!
!*    3.      Fills output arguments
!             ----------------------
!
OBUDGET_TOPD = LBUDGET_TOPD
OSTOCK_TOPD = LSTOCK_TOPD
KNB_TOPD = NNB_TOPD
KNB_STOCK = NNB_STP_STOCK
KNB_RESTART = NNB_STP_RESTART
KFREQ_MAPS_WG = NFREQ_MAPS_WG
KFREQ_MAPS_ASAT = NFREQ_MAPS_ASAT
KFREQ_MAPS_RUNOFF = NFREQ_MAPS_RUNOFF
PSPEEDR(1:NNCAT) = XSPEEDR(1:NNCAT)
PSPEEDG(1:NNCAT) = XSPEEDG(1:NNCAT)
WHERE(XSPEEDH(1:NNCAT)/=0.3)
  PSPEEDH(1:NNCAT) = XSPEEDH(1:NNCAT)
ELSEWHERE
  PSPEEDH(1:NNCAT) = XSPEEDR(1:NNCAT)/10.
ENDWHERE
OSPEEDR_VAR = LSPEEDR_VAR
PQINIT(1:NNCAT) = XQINIT(1:NNCAT)
PRTOP_D2(1:NNCAT) = XRTOP_D2(1:NNCAT)
OPERT_PARAM = LPERT_PARAM
OPERT_INIT  = LPERT_INIT
!
WRITE(ILUOUT,*) 'NAM_TOPD'
WRITE(ILUOUT,*) 'LBUDGET = ',LBUDGET_TOPD
WRITE(ILUOUT,*) 'NNB_TOP = ',NNB_TOPD
WRITE(ILUOUT,*) 'LSTOCK = ',LSTOCK_TOPD
WRITE(ILUOUT,*) 'NNB_RESTART,NNB_STOCK = ',NNB_STP_RESTART,NNB_STP_STOCK
WRITE(ILUOUT,*) 'NFREQ_MAPS_WG,NFREQ_MAPS_ASAT = ',NFREQ_MAPS_WG,NFREQ_MAPS_ASAT
WRITE(ILUOUT,*) 'NFREQ_MAPS_RUNOFF = ',NFREQ_MAPS_RUNOFF
WRITE(ILUOUT,*) 'XSPEEDR = ',XSPEEDR
WRITE(ILUOUT,*) 'LSPEEDR_VAR = ',LSPEEDR_VAR
WRITE(ILUOUT,*) 'LPERT_PARAM = ',LPERT_PARAM
WRITE(ILUOUT,*) 'LPERT_INIT  = ',LPERT_INIT 
!
IF (LHOOK) CALL DR_HOOK('READ_NAM_TOPD',1,ZHOOK_HANDLE)
!-------------------------------------------------------------------------------
!
END SUBROUTINE READ_NAM_TOPD
