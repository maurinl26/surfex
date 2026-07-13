!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!
!     ##########################
      SUBROUTINE TOPD_TO_DF (IO, NK, NP, NPE, PWG)
!     ##########################
!
!!
!!    PURPOSE
!!    -------
!     This routines updates the soil water content of ISBA DIF afeter TOPODYN
!     lateral distribution  
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
!!       ELYAZIDI/HEYMES/RISTOR * Meteo-France *
!!
!!    MODIFICATIONS
!!    -------------
!!
!!      Original  02/2011
!!      Modif : correction Nature grid considered instead of full grid 02/2017
!-------------------------------------------------------------------------------
!
!*       0.     DECLARATIONS
!               ------------
!
USE MODD_ISBA_OPTIONS_n, ONLY : ISBA_OPTIONS_t
USE MODD_ISBA_n, ONLY : ISBA_NP_t, ISBA_NPE_t, ISBA_NK_t, ISBA_K_t, ISBA_P_t, ISBA_PE_t
!
USE MODD_SURF_PAR,      ONLY : XUNDEF, NUNDEF
USE MODD_COUPLING_TOPD, ONLY : XATOP_NATURE, XFRAC_D3,XWOVSATI_P,XDMAXFC
USE MODD_ISBA_PAR,      ONLY : XWGMIN
!
USE YOMHOOK   ,         ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,         ONLY : JPRB
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
TYPE(ISBA_OPTIONS_t), INTENT(INOUT) :: IO
TYPE(ISBA_NK_t), INTENT(INOUT) :: NK
TYPE(ISBA_NP_t), INTENT(INOUT) :: NP
TYPE(ISBA_NPE_t), INTENT(INOUT) :: NPE
!
 REAL, DIMENSION(:,:), INTENT(IN) :: PWG
!      
!*      0.2    declarations of local variables
!
TYPE(ISBA_K_t), POINTER :: KK
TYPE(ISBA_P_t), POINTER :: PK
TYPE(ISBA_PE_t), POINTER :: PEK
REAL                              :: ZWORK          ! numbers of layers in root and deep zones
INTEGER                           :: ZDEPTH
INTEGER                           :: JI, JL, JP ! loop indexes
REAL(KIND=JPRB)                   :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('TOPD_TO_DF',0,ZHOOK_HANDLE)
!
DO JP=1,IO%NPATCH
  !
  KK => NK%AL(JP)
  PK => NP%AL(JP)
  PEK => NPE%AL(JP)
  !
  IF (PK%NSIZE_P == 0 ) CYCLE
  !
  DO JL = 2,IO%NGROUND_LAYER
    !
    DO JI=1,PK%NSIZE_P
      !
      ZDEPTH=PK%XRUNOFFD(JI)!only on the layers where runoff is authorized
      IF(PK%XDG(JI,JL)<=ZDEPTH.AND.ZDEPTH/=XUNDEF.AND.(XATOP_NATURE(JI)/=0.0)&
                                         .AND.(XATOP_NATURE(JI)/=XUNDEF)) THEN
        ! root layers
        IF (PK%XDZG(JI,JL)/=XUNDEF.AND.PK%XDG2(JI)/=XUNDEF.AND.PK%XDG(JI,JL)/=XUNDEF)&! 
        ZWORK=MIN(PK%XDZG(JI,JL),MAX(0.0,PK%XDG2(JI)-PK%XDG(JI,JL)+PK%XDZG(JI,JL)))
        !
        IF ((PWG(JI,2)/=XUNDEF).AND.(ZWORK>0.).AND.(ZWORK/=XUNDEF)&
                             .AND.(KK%XWSAT(JI,JL)/=XUNDEF) )THEN 
          PEK%XWG(JI,JL)=MAX(PWG(JI,2),XWGMIN) 
          IF (PEK%XWG(JI,JL)>KK%XWSAT(JI,JL))THEN
            XWOVSATI_P(JI) = XWOVSATI_P(JI)+MAX(0.,PEK%XWG(JI,JL) - KK%XWSAT(JI,JL))
            PEK%XWG(JI,JL)=KK%XWSAT(JI,JL) 
          ENDIF
        ENDIF
      ENDIF
      !
    ENDDO
    !
  ENDDO
  !
ENDDO 
!
IF (LHOOK) CALL DR_HOOK('TOPD_TO_DF',1,ZHOOK_HANDLE)

END SUBROUTINE TOPD_TO_DF


