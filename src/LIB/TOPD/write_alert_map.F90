!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!-----------------------------------------------------------------
!     ##########################
      SUBROUTINE WRITE_ALERT_MAP(KCAT,KSTEP,HSTEP)
!     ##########################
!
!!
!!    PURPOSE
!!    -------
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
!!    REFERENCE
!!    ---------
!!     
!!    AUTHOR
!!    ------
!!
!!      B. Vincendon	* Meteo-France *
!!
!!    MODIFICATIONS
!!    -------------
!!
!!      Original   25/12/2014
!-------------------------------------------------------------------------------
!
!*       0.     DECLARATIONS
!               ------------
!
USE MODD_TOPODYN, ONLY : CCAT, NMESHT, NNMC, NNYC, NNXC, XX0, XY0, XDXT, NLINE, &
                         XDHIL, XQTOT_SUB, NCAT_CAT_TO_SUB
!
USE MODD_COUPLING_TOPD, ONLY : XQ2,XQ10,XQ50
!
USE MODI_GET_LUOUT
USE MODI_OPEN_FILE
USE MODI_CLOSE_FILE
!
USE YOMHOOK   ,ONLY : LHOOK, DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
INTEGER,          INTENT(IN) :: KCAT
INTEGER,          INTENT(IN) :: KSTEP
CHARACTER(LEN=*), INTENT(IN) :: HSTEP  ! atmospheric loop index
!
!*      0.2    declarations of local variables
CHARACTER(LEN=50)       :: CNAME
INTEGER                 :: JJ,JI,JSCAT
INTEGER                 :: IINDEX ! reference number of the pixel
INTEGER                 :: IUNIT,ILUOUT
REAL                    :: ZTMP_Q2,ZTMP_Q10,ZTMP_Q50
REAL                    :: ZWORK
REAL, DIMENSION(NMESHT) :: ZALERT
REAL, DIMENSION(NMESHT) :: ZQ
REAL(KIND=JPRB)         :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('WRITE_ALERT_MAP',0,ZHOOK_HANDLE)
!
!*       0.     Initialization:
!               ---------------
!
CALL GET_LUOUT('OFFLIN',ILUOUT)
!
CNAME = TRIM(CCAT(KCAT))//TRIM('.Alert')//HSTEP
write(*,*) 'WRITE_ALERT_MAP : ',CNAME
!
CALL OPEN_FILE('ASCII ',IUNIT,HFILE=CNAME,HFORM='FORMATTED',HACTION='WRITE')
!
WRITE(IUNIT,*) 'Pixel_Ref Alert_level'
!
ZTMP_Q2=300.
ZTMP_Q10=600.
ZTMP_Q50=1200.
ZALERT(:)=0.
ZQ(:)=XQTOT_SUB(KCAT,:,KSTEP)
ZWORK=MINVAL(XDHIL(KCAT,1:NNMC(KCAT)))
!
DO JJ=1,NNMC(KCAT)
  IF(XDHIL(KCAT,JJ)==ZWORK) THEN 
    ZALERT(JJ) = 1.
    JSCAT=NCAT_CAT_TO_SUB(KCAT,JJ)
    IF ((ZQ(JSCAT)>=XQ2(KCAT,JSCAT)).AND.(ZQ(JSCAT)<XQ10(KCAT,JSCAT)))THEN
      ZALERT(JJ) = 2.
    ELSEIF((ZQ(JSCAT)>=XQ10(KCAT,JSCAT)).AND.(ZQ(JSCAT)<XQ50(KCAT,JSCAT)))THEN
      ZALERT(JJ) = 3.
    ELSEIF (ZQ(JSCAT)>=XQ50(KCAT,JSCAT)) THEN
      ZALERT(JJ) = 4.
    ENDIF
  ENDIF
ENDDO
!
DO JJ=1,NNYC(KCAT)
  DO JI=1,NNXC(KCAT)
    IINDEX = (JJ - 1) * NNXC(KCAT) + JI
    IF (NLINE(KCAT,IINDEX)/=0) THEN
      IF (ZALERT(NLINE(KCAT,IINDEX))/=0.) &
        WRITE(IUNIT,*) XX0(KCAT)+(JI-1)*XDXT(KCAT),XY0(KCAT)+(JJ-1)*XDXT(KCAT),&
                     ZALERT(NLINE(KCAT,IINDEX))
    ENDIF
  ENDDO
ENDDO
!
CALL CLOSE_FILE('ASCII ',IUNIT)
!
IF (LHOOK) CALL DR_HOOK('WRITE_ALERT_MAP',1,ZHOOK_HANDLE)
!
END SUBROUTINE WRITE_ALERT_MAP
