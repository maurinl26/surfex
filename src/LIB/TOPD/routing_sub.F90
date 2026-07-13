!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!-------------------------------------------------------------------------------
!     ####################
      SUBROUTINE ROUTING_SUB(PRO,PDR,KSTEP)
!     ####################
!
!!****  *ROUTING*  
!!
!!    PURPOSE
!!    -------
!     To route the runoff and the exfiltration discharge to the catchment outlet
!         
!     
!!**  METHOD
!!    ------
!
!!    EXTERNAL
!!    --------
!!
!!    none
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------ 
!!
!!
!!      
!!    REFERENCE
!!    ---------
!!
!!    
!!     
!!    AUTHOR
!!    ------
!!       B Vincendon	* Meteo-France *
!!
!!    MODIFICATIONS
!!    -------------
!!
!!      Original   04/2014
!!      Modif B Vincendon 07/2017 : adding possibility of varing river speed according
!!                                  to the discharge value(LSPEED_VAR)
!!     
!-------------------------------------------------------------------------------
!
!*       0.     DECLARATIONS
!               ------------
!
USE MODD_SURF_PAR   ,ONLY: XUNDEF
!
USE MODD_TOPODYN    ,ONLY: XTOPD_STEP, NNCAT, XQTOT, NNMC, XDRIV, XDHIL,    &
                         XSPEEDR, XSPEEDH, XTIME_TOPD_SUB, XQB_RUN_SUB,     &
                         XQB_DR_SUB, XTIME_TOPD_DRAIN_SUB, NNB_TOPD_STEP,   &
                         XTIME_TOPD, NNCAT_SUB, XQTOT_SUB, NCAT_CAT_TO_SUB, &
                         NMASK_OUTLET,XDIST_OUTLET, NNCAT_MAX, XA_SPEED,    &
                         XB_SPEED, XMAX_SPEED, LSPEEDR_VAR
!
USE MODI_GET_LUOUT
USE MODI_OPEN_FILE
USE MODI_CLOSE_FILE
USE MODI_WRITE_FILE_1MAP
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
REAL, DIMENSION(:,:), INTENT(IN) :: PRO     ! Total water for runoff for each pixel (m3/s)
!ludo
REAL, DIMENSION(:,:), INTENT(IN) :: PDR     ! Total water for drainage for each pixel
INTEGER, INTENT(IN)              :: KSTEP   ! current integration step
!
!
!*      0.2    declarations of local variables
!
!
INTEGER                    :: JCAT, JJ,JO,JSCAT ! Loop variables
INTEGER                    :: JSBG,JSBP,JTG,JTP
INTEGER                    :: JSTEP     ! current or future integration steps
REAL                       :: ZTMP
REAL                       :: ZA,ZB
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('ROUTING_SUB',0,ZHOOK_HANDLE)
!
!*       1.0.     Initialization :
!               --------------
!
! 
DO JCAT=1,NNCAT
  !
  IF (LSPEEDR_VAR) THEN 
    ZA=XA_SPEED(JCAT)
    ZB=XB_SPEED(JCAT) !*0.6
    IF (KSTEP>1)THEN
      XSPEEDR(JCAT)=MIN(EXP(ZA*LOG(MAX(XQTOT(JCAT,KSTEP),0.001))+ZB),XMAX_SPEED(JCAT))
    ENDIF
    write(*,*) 'VARIABLE SPEED IS ACTIVATED, XSPEEDR BV,',JCAT,'= ',XSPEEDR(JCAT)
  ENDIF
  !
  DO JJ=1,NNMC(JCAT)
    JO=NMASK_OUTLET(JCAT,JJ)
    !
    !*       2.0    Runoff by geomorpho transfer function
    !               -------------------------------------
    !
    JSCAT=NCAT_CAT_TO_SUB(JCAT,JJ)
    !
    IF (LSPEEDR_VAR) THEN 
      !
      IF ( XSPEEDH(JCAT) > 0.0 .AND. XSPEEDR(JCAT) > 0.0 .AND.&
         XDHIL(JCAT,JJ) < XUNDEF .AND. XDRIV(JCAT,JJ) < XUNDEF) &
        XTIME_TOPD(JCAT,JJ)=XDHIL(JCAT,JJ)/XSPEEDH(JCAT) +XDRIV(JCAT,JJ)/XSPEEDR(JCAT)
    ENDIF
    !
    IF ( PRO(JCAT,JJ) > 0.0 .AND. PRO(JCAT,JJ) < XUNDEF&
         .AND.XTOPD_STEP > 0.0 .AND. XTIME_TOPD(JCAT,JJ) > 0.0 &
         .AND. XTIME_TOPD(JCAT,JJ)< XUNDEF  ) THEN
      !
      JSTEP = NINT(XTIME_TOPD_SUB(JCAT,JJ) / XTOPD_STEP) + KSTEP 
      !
      IF ((JSTEP <= NNB_TOPD_STEP).AND.( JSTEP /= 0))  THEN
        XQTOT_SUB(JCAT,JSCAT,JSTEP) = XQTOT_SUB(JCAT,JSCAT,JSTEP) + PRO(JCAT,JJ)
        XQB_RUN_SUB(JCAT,JSCAT,JSTEP) = XQB_RUN_SUB(JCAT,JSCAT,JSTEP) + PRO(JCAT,JJ)
      ENDIF!JSTEP
      !
    ENDIF!PRO
    !
    !
    !*       3.0    Drainage by geomorpho transfer function
    !               -------------------------------------
    !
    IF ((PDR(JCAT,JJ) > 0.0).AND.(PDR(JCAT,JJ)<XUNDEF)&
         .AND.XTOPD_STEP > 0.0 .AND. XTIME_TOPD(JCAT,JJ) > 0.0 &
         .AND. XTIME_TOPD(JCAT,JJ)< XUNDEF  ) THEN
      JSTEP = NINT(XTIME_TOPD_DRAIN_SUB(JCAT,JJ) / XTOPD_STEP) + KSTEP 
      !
      IF ( JSTEP.LE.NNB_TOPD_STEP.AND.( JSTEP /= 0)) THEN
        XQTOT_SUB(JCAT,JSCAT,JSTEP) = XQTOT_SUB(JCAT,JSCAT,JSTEP) + PDR(JCAT,JJ)
        XQB_DR_SUB(JCAT,JSCAT,JSTEP) = XQB_DR_SUB(JCAT,JSCAT,JSTEP) + PDR(JCAT,JJ)
      ENDIF!JSTEP
      !
    ENDIF!PDR
    !
  ENDDO!JJ
  !
  !*       3.0    Treatment of overlapped sub-catchments
  !               -------------------------------------
  !
  IF (KSTEP==NNB_TOPD_STEP) THEN
    !
    DO JSBP=NNCAT_SUB(JCAT),2,-1 !petit bassin
      !
      DO JSBG=JSBP-1,1,-1 !grand bassin
        !
        DO JTP=1,NNB_TOPD_STEP
          !
          IF (XDIST_OUTLET(JCAT,JSBG,JSBP)/=XUNDEF)THEN
            !
            JTG = NINT(XDIST_OUTLET(JCAT,JSBG,JSBP) / XSPEEDR(JCAT)/ XTOPD_STEP) +  JTP
            !
            IF ( JTG.LE.NNB_TOPD_STEP ) THEN
              !
              XQTOT_SUB(JCAT,JSBG,JTG) = XQTOT_SUB(JCAT,JSBG,JTG) + XQTOT_SUB(JCAT,JSBP,JTP)
              !
            ENDIF!JTG
          ENDIF!XDIST
        ENDDO!JTP
      ENDDO!JSBG
    ENDDO!JSBP
  ENDIF!KSTEP
ENDDO!JCAT

IF (LHOOK) CALL DR_HOOK('ROUTING_SUB',1,ZHOOK_HANDLE)
!
END SUBROUTINE ROUTING_SUB
