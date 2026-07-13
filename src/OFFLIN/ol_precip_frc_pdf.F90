!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!     ######spl
SUBROUTINE OL_PRECIP_FRC_PDF (KNB_ATM,OPRINT,PTSTEP,PPDISTRIB)  
!**************************************************************************
!
!!    PURPOSE
!!    -------
!        Generate disaggregation PDF for precipitation time interpolation
!!
!!**  METHOD
!!    ------
!!
!!    EXTERNAL
!!    --------
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
!!      B. Decharme   *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    03/2020
!
!
USE MODD_SURF_PAR,   ONLY : XUNDEF
!
USE MODI_GET_LUOUT
USE MODI_ABOR1_SFX
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
! global variables
!
INTEGER, INTENT(IN) :: KNB_ATM
LOGICAL, INTENT(IN) :: OPRINT
REAL,    INTENT(IN) :: PTSTEP
!
REAL, DIMENSION(:),  INTENT(OUT) :: PPDISTRIB
!
! local variables
!
REAL, PARAMETER          :: ZEXPO = -3.0
!
REAL, DIMENSION(KNB_ATM) :: ZU, ZV
!
REAL :: ZFACTOR, ZRTSUM, ZRTSUM0, ZNB_ATM
REAL :: ZD0, ZP0, ZP1, ZCOEF
!
INTEGER :: ILUOUT, JFRC
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!========================================================================
!
IF (LHOOK) CALL DR_HOOK('OL_PRECIP_FRC_PDF',0,ZHOOK_HANDLE)
!
CALL GET_LUOUT('OFFLIN',ILUOUT)
!
IF(PTSTEP<=1800.)THEN
   WRITE(ILUOUT,*)'OL_PRECIP_FRC_PDF: PDF Precipitation interpolation should not be used for forcing time step < 30 min'
   CALL ABOR1_SFX('OL_PRECIP_FRC_PDF: PDF Precipitation interpolation should not be used for forcing time step < 30 min')
ELSEIF(PTSTEP<=3600.)THEN
   ZFACTOR = 29.63  ! For 0% dry time steps in a rainy forcing interval
ELSE
   ZFACTOR = 66.6666 ! For 33% dry time steps in a rainy forcing interval
ENDIF
!
ZNB_ATM = REAL(KNB_ATM)
ZRTSUM  = 0.0
ZD0     = ZFACTOR * ZNB_ATM ** ZEXPO
ZCOEF   = (ZEXPO+1.)/ZEXPO
!
DO JFRC = 1, KNB_ATM
!
   ZRTSUM0 = ZRTSUM
   ZP0     = REAL(JFRC-1)/ZNB_ATM
   ZP1     = REAL(JFRC  )/ZNB_ATM
!
   PPDISTRIB(JFRC) = ZD0*ZEXPO/((ZEXPO+1)*(ZP1-ZP0)*100.) &
                   * ((100.*ZP1/ZD0)**ZCOEF - (100.*ZP0/ZD0)**ZCOEF)
!           
   ZRTSUM = ZRTSUM + PPDISTRIB(JFRC)
!
   IF(ZRTSUM>1.0)THEN
     IF(ZRTSUM0<1.0)THEN
       PPDISTRIB(JFRC) = PPDISTRIB(JFRC) + 1.0 - ZRTSUM   ! Stick any remaining weight in last interval
     ELSE
       PPDISTRIB(JFRC) = 0.0
     ENDIF
   ENDIF
!
   WRITE(ILUOUT,FMT='(A23,I2.2,A3,f11.6,A2)')'Precip time disag step ',JFRC,' = ',PPDISTRIB(JFRC)*100.,' %'
   IF(OPRINT) WRITE(*,FMT='(A23,I2.2,A3,f11.6,A2)')'Precip time disag step ',JFRC,' = ',PPDISTRIB(JFRC)*100.,' %'
!
ENDDO
!
! Ensure integral of PDF = 1.0
!
IF (ZRTSUM < 1.0) THEN
   PPDISTRIB(1) = PPDISTRIB(1) + 1.0 - ZRTSUM  
ENDIF
!
IF (LHOOK) CALL DR_HOOK('OL_PRECIP_FRC_PDF',1,ZHOOK_HANDLE)
!
!========================================================================
!
END SUBROUTINE OL_PRECIP_FRC_PDF
