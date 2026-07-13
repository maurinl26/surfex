!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!-----------------------------------------------------------------
!     ##########################
      SUBROUTINE WRITE_FILE_1MAP(PVAR,HVAR,KCAT)
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
!!     from WRITE_FILE_MAP
!!
!!    MODIFICATIONS
!!    -------------
!!
!!      Original   march 2013
!-------------------------------------------------------------------------------
!
!*       0.     DECLARATIONS
!               ------------
!
USE MODD_TOPODYN,  ONLY : CCAT, NNCAT, NNYC, NNXC, XX0, XY0, XDXT, XTOPD
!
USE MODD_SURF_PAR, ONLY : XUNDEF,NUNDEF
!
USE MODI_GET_LUOUT
USE MODI_OPEN_FILE
USE MODI_CLOSE_FILE
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
REAL, DIMENSION(:),   INTENT(IN) :: PVAR   ! variable to write in the file
CHARACTER(LEN=30),    INTENT(IN) :: HVAR   ! end name of the file
INTEGER,              INTENT(IN) :: KCAT   ! nb of the "big" catchment
!
!*      0.2    declarations of local variables
!
CHARACTER(LEN=50),DIMENSION(NNCAT) :: CNAME
INTEGER                    :: JJ,JI
INTEGER                    :: IINDEX    ! reference number of the pixel
INTEGER                    :: IUNIT,ILUOUT
REAL                       :: ZOUT      ! pixel not included in the catchment
REAL                       :: ZMIN,ZMAX
REAL                       :: ZX1, ZY1  ! left top and right bottom pixels coordinates
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('WRITE_FILE_1MAP',0,ZHOOK_HANDLE)
!
!*       0.     Initialization:
!               ---------------
!

CALL GET_LUOUT('OFFLIN',ILUOUT)

ZOUT = 0.0
!
CNAME(KCAT) = TRIM(CCAT(KCAT))//TRIM(HVAR)//'.map'
!
WRITE(ILUOUT,*) CNAME(KCAT)
!
CALL OPEN_FILE('ASCII ',IUNIT,HFILE=CNAME(KCAT),HFORM='FORMATTED')
!
!*       1.0    writing header map file
!               --------------------------------------
!
ZX1 = XX0(KCAT)
ZY1 = XY0(KCAT) + ( (NNYC(KCAT)-1) * XDXT(KCAT) )
!
ZMIN = MINVAL(PVAR(:))
ZMAX = MAXVAL(PVAR(:),MASK=(PVAR(:)/=XUNDEF.AND.PVAR(:)/=NUNDEF*1.))
!
write(*,*) 'WRITE_FILE 1MAP=> ZMIN,ZMAX',ZMIN,ZMAX
!
DO JJ=1,5
  WRITE(IUNIT,*)
ENDDO
!
WRITE(IUNIT,*) XX0(KCAT)
WRITE(IUNIT,*) XY0(KCAT)
WRITE(IUNIT,*) NNXC(KCAT) 
WRITE(IUNIT,*) NNYC(KCAT)
WRITE(IUNIT,*) ZOUT
WRITE(IUNIT,*) XDXT(KCAT)
WRITE(IUNIT,*) ZMIN
WRITE(IUNIT,*) ZMAX
!
DO JJ=1,NNYC(KCAT)
  !
  DO JI=1,NNXC(KCAT)
    !
    IINDEX = (JJ - 1) * NNXC(KCAT) + JI
    ZX1 = XX0(KCAT) + ((JI-1) * XDXT(KCAT))
    ZY1 = XY0(KCAT) + ((JJ-1) * XDXT(KCAT))
    !
    IF ( XTOPD(KCAT,IINDEX).EQ.XUNDEF ) THEN
      WRITE(IUNIT,*) ZOUT
    ELSE
      WRITE(IUNIT,*) PVAR(IINDEX)
    ENDIF
    !
  ENDDO
  !
ENDDO
!
CALL CLOSE_FILE('ASCII ',IUNIT)
!
IF (LHOOK) CALL DR_HOOK('WRITE_FILE_1MAP',1,ZHOOK_HANDLE)
!
END SUBROUTINE WRITE_FILE_1MAP
