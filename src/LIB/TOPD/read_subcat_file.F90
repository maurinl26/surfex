!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!-----------------------------------------------------------------
!     #######################     
      SUBROUTINE READ_SUBCAT_FILE(HPROGRAM,HFILE,HFORM,KCAT)
!     #######################
!
!!****  *READ_SUBCAT_FILE*  
!!
!!    PURPOSE
!!    -------
!     This routine aims at reading subcat files
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
!!
!!      
!!    REFERENCE
!!    ---------
!!
!!    
!!      
!!    AUTHOR
!!    ------
!!
!!      B. Vincendon    * Meteo-France *
!!
!!    MODIFICATIONS
!!    -------------
!!
!!      Original   01/2015
!-------------------------------------------------------------------------------
!
!*       0.     DECLARATIONS
!               ------------
!
USE MODI_GET_LUOUT
USE MODI_OPEN_FILE
USE MODI_CLOSE_FILE
!
USE MODD_COUPLING_TOPD, ONLY : NSUBCAT,XLX,XLY,CSUBCAT,XQ2,XQ10,XQ50
!
USE YOMHOOK   ,ONLY : LHOOK, DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
 CHARACTER(LEN=*),  INTENT(IN)  :: HPROGRAM    !
 CHARACTER(LEN=*),  INTENT(IN)  :: HFILE       ! File to be read
 CHARACTER(LEN=*),  INTENT(IN)  :: HFORM       ! Format of the file to be read
INTEGER,           INTENT(IN)  :: KCAT       ! Number of pixels in the catchment
!
!*      0.2    declarations of local variables
!
!
INTEGER                   :: JJ ! loop control 
INTEGER                   :: IUNIT       ! Unit of the files
INTEGER                   :: ILUOUT      ! Unit of the files
!
REAL                      :: ZWRK        ! work variable
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
CHARACTER(LEN=100)    :: YHEADER    ! Header File to be read
!------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('READ_SUBCAT_FILE',0,ZHOOK_HANDLE)
!
!*       0.2    preparing file openning
!               ----------------------
 CALL GET_LUOUT(HPROGRAM,ILUOUT)
!
 CALL OPEN_FILE(HPROGRAM,IUNIT,HFILE,HFORM,HACTION='READ')
!
READ(IUNIT,*) ZWRK
NSUBCAT(KCAT)=INT(ZWRK)
READ(IUNIT,*)
!
DO JJ=1,NSUBCAT(KCAT)
  READ(IUNIT,*) CSUBCAT(KCAT,JJ),XLX(KCAT,JJ),XLY(KCAT,JJ),&
                        XQ2(KCAT,JJ),XQ10(KCAT,JJ),XQ50(KCAT,JJ)
ENDDO

CALL CLOSE_FILE(HPROGRAM,IUNIT)
!
IF (LHOOK) CALL DR_HOOK('READ_SUBCAT_FILE',1,ZHOOK_HANDLE)
!
END SUBROUTINE READ_SUBCAT_FILE







