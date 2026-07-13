!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!-----------------------------------------------------------------
!     #######################
      SUBROUTINE INIT_TOPD_PGD(HPROGRAM)
!     #######################
!
!!****  *INIT_TOPD_PGD*  
!!
!!    PURPOSE
!!    -------
!     This routine aims at initialising the variables 
!     needed of running Topmodel for PGD step.
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
!!      Original   03/2014
!-------------------------------------------------------------------------------
!
!*       0.     DECLARATIONS
!               ------------
!
USE MODD_COUPLING_TOPD, ONLY : LCOUPL_TOPD,LDUMMY_SUBCAT, LSUBCAT
USE MODD_TOPODYN,       ONLY : NNCAT
!
!
USE MODI_GET_LUOUT
USE MODI_INIT_TOPD
USE MODI_DESC_CATCHMENTS
USE MODI_FIND_SUB_CATCH
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
 CHARACTER(LEN=*), INTENT(IN) :: HPROGRAM    !
!
!*      0.2    declarations of local variables
!
!
INTEGER                   :: JCAT     ! loop control
INTEGER                   :: ILUOUT   ! Unit of the files
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('INIT_TOPD_PGD',0,ZHOOK_HANDLE)
!
!*       1    Initialization:
!               ---------------
!
CALL GET_LUOUT(HPROGRAM,ILUOUT)
!
WRITE(ILUOUT,*) 'INITIALISATION INIT_TOPD_PGD'
!
CALL INIT_TOPD('ASCII ')
!
write(*,*) 'INIT_TOPD_PGD ',LCOUPL_TOPD,LDUMMY_SUBCAT,LSUBCAT
DO JCAT=1,NNCAT
  IF (LDUMMY_SUBCAT) THEN
    CALL DESC_CATCHMENTS(JCAT)
  ELSEIF (LSUBCAT) THEN
    CALL FIND_SUB_CATCH(JCAT)
  ENDIF
ENDDO
 !
IF (LHOOK) CALL DR_HOOK('INIT_TOPD_PGD',1,ZHOOK_HANDLE)
!
END SUBROUTINE INIT_TOPD_PGD
