!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!-----------------------------------------------------------------
!--------------- special set of characters for SCCS information
!-----------------------------------------------------------------
!
SUBROUTINE WRITE_TIME(ITIME,ISPACE,HSEP,HTDATE)
!
!
!*       0.     DECLARATIONS
!
USE YOMHOOK   ,ONLY : LHOOK, DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
INTEGER, INTENT(IN)             :: ITIME
INTEGER, INTENT(IN)             :: ISPACE
CHARACTER(LEN=*), INTENT(IN)    :: HSEP
CHARACTER(LEN=*), INTENT(INOUT) :: HTDATE
!
!*      0.2    declarations of local variables
CHARACTER(LEN=10)               :: YPAS
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('WRITE_TIME',0,ZHOOK_HANDLE)
!
IF (ITIME.LT.10) THEN
  WRITE(YPAS,'(i1)') ITIME
  IF (ISPACE==1) THEN
    HTDATE=trim(HTDATE)//" 0"//trim(YPAS)//HSEP
  ELSE
    HTDATE=trim(HTDATE)//"0"//trim(YPAS)//HSEP
  ENDIF
ELSE
  IF (ITIME.LT.100) THEN
    WRITE(YPAS,'(i2)') ITIME
  ELSE
    WRITE(YPAS,'(i4)') ITIME
  ENDIF
  IF (ISPACE==1) THEN
    HTDATE=trim(HTDATE)//" "//trim(YPAS)//HSEP
  ELSE
    HTDATE=trim(HTDATE)//trim(YPAS)//HSEP
  ENDIF  
ENDIF
!
IF (LHOOK) CALL DR_HOOK('WRITE_TIME',1,ZHOOK_HANDLE)
!
END SUBROUTINE WRITE_TIME

