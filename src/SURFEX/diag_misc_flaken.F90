!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!     #########
       SUBROUTINE DIAG_MISC_FLAKE_n (DMF, F, PQ_BOT_FLK, PQ_W_FLK, PQ_SNOW_FLK, PQ_ICE_FLK,    &
                                     PI_ATM_FLK, PI_SNOW_FLK, PI_ICE_FLK, PI_W_FLK, PI_H_FLK,  &
                                     PI_BOT_FLK, PI_INTM_0_H_FLK, PI_INTM_H_D_FLK, PQ_STAR_FLK,&
                                     PRAIN, PSNOW, PEVAP, PTSTEP)
!     ###############################################################################
!
!!****  *DIAG_MISC-FLAKE_n * - additional diagnostics for FLake
!!
!!    PURPOSE
!!    -------
!
!!**  METHOD
!!    ------
!!
!!    REFERENCE
!!    ---------
!!      
!!
!!    AUTHOR
!!    ------
!!     P. Le Moigne 
!!
!!    MODIFICATIONS
!!    -------------
!!      Original     10/2005
!!      P. Le Moigne 05/2023 : temperature profile in sediments
!!      P. Le Moigne 08/2023 : heat storage flux
!!      P. Le Moigne 11/2023 : effective rain P-E
!!------------------------------------------------------------------
!
!
USE MODD_FLAKE_n, ONLY : FLAKE_t
USE MODD_DIAG_MISC_FLAKE_n, ONLY : DIAG_MISC_FLAKE_t
!
USE MODD_SURF_PAR,           ONLY : XUNDEF
!
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
!
TYPE(FLAKE_t), INTENT(INOUT) :: F
TYPE(DIAG_MISC_FLAKE_t), INTENT(INOUT) :: DMF
!
REAL, DIMENSION(SIZE(F%XT_WML)), INTENT(IN) :: PQ_BOT_FLK
REAL, DIMENSION(SIZE(F%XT_WML)), INTENT(IN) :: PQ_W_FLK
REAL, DIMENSION(SIZE(F%XT_WML)), INTENT(IN) :: PQ_SNOW_FLK
REAL, DIMENSION(SIZE(F%XT_WML)), INTENT(IN) :: PQ_ICE_FLK
REAL, DIMENSION(SIZE(F%XT_WML)), INTENT(IN) :: PI_ATM_FLK
REAL, DIMENSION(SIZE(F%XT_WML)), INTENT(IN) :: PI_SNOW_FLK
REAL, DIMENSION(SIZE(F%XT_WML)), INTENT(IN) :: PI_ICE_FLK
REAL, DIMENSION(SIZE(F%XT_WML)), INTENT(IN) :: PI_W_FLK
REAL, DIMENSION(SIZE(F%XT_WML)), INTENT(IN) :: PI_H_FLK
REAL, DIMENSION(SIZE(F%XT_WML)), INTENT(IN) :: PI_BOT_FLK
REAL, DIMENSION(SIZE(F%XT_WML)), INTENT(IN) :: PI_INTM_0_H_FLK
REAL, DIMENSION(SIZE(F%XT_WML)), INTENT(IN) :: PI_INTM_H_D_FLK
REAL, DIMENSION(SIZE(F%XT_WML)), INTENT(IN) :: PQ_STAR_FLK
!
REAL, DIMENSION(SIZE(F%XT_WML)), INTENT(IN) :: PRAIN  ! rainfall (kg/m2/s)
REAL, DIMENSION(SIZE(F%XT_WML)), INTENT(IN) :: PSNOW  ! snowfall (kg/m2/s)
REAL, DIMENSION(SIZE(F%XT_WML)), INTENT(IN) :: PEVAP  ! evaporation (kg/m2/s)
!
REAL, INTENT(IN) :: PTSTEP 
!
!*      0.2    declarations of local variables
!
REAL, DIMENSION(SIZE(DMF%XZWPROF),SIZE(F%XT_WML)) :: ZCSIW      ! Vertical normalized coordinate
REAL, DIMENSION(SIZE(DMF%XZWPROF),SIZE(F%XT_WML)) :: ZSHAPEW    ! Shape function
REAL, DIMENSION(SIZE(DMF%XZSPROF),SIZE(F%XT_WML)) :: ZCSIS      ! Vertical normalized coordinate in the sediment layer
REAL, DIMENSION(SIZE(DMF%XZSPROF),SIZE(F%XT_WML)) :: ZSHAPES    ! Shape function for the sediment layer
!
INTEGER         :: IZW, IZS
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('DIAG_MISC_FLAKE_N',0,ZHOOK_HANDLE)
!
!* Flake temperature profile
!
DMF%XTWPROF(:,:) = XUNDEF
!
IF (DMF%LWATER_PROFILE) THEN
!
   DO IZW=1,SIZE(DMF%XZWPROF)
      WHERE (F%XWATER_DEPTH(:)==F%XH_ML(:))
         ZCSIW(IZW,:) = 0.
      ELSEWHERE
         ZCSIW(IZW,:) = (DMF%XZWPROF(IZW) - F%XH_ML(:))/(F%XWATER_DEPTH(:) - F%XH_ML(:))
      END WHERE
      ZSHAPEW(IZW,:) = (40./3.*F%XCT-20./3.)*ZCSIW(IZW,:)   +     (18.-30.*F%XCT)*ZCSIW(IZW,:)**2 &
                       + (20.*F%XCT-12.)   *ZCSIW(IZW,:)**3+(5./3.-10./3.*F%XCT)*ZCSIW(IZW,:)**4  
   END DO
!
   DO IZW=1,SIZE(DMF%XZWPROF)
      WHERE (F%XH_ML(:) >= DMF%XZWPROF(IZW))
         DMF%XTWPROF(IZW,:) =  F%XT_WML(:) 
      ELSEWHERE (F%XWATER_DEPTH(:) >= DMF%XZWPROF(IZW)) 
         DMF%XTWPROF(IZW,:) = F%XT_WML(:) - (F%XT_WML(:) - F%XT_BOT(:)) * ZSHAPEW(IZW,:)
      END WHERE
   END DO
!
END IF
!
DMF%XTSPROF(:,:) = XUNDEF
!
IF (DMF%LSEDIM_PROFILE) THEN
!
   DO IZS=1,SIZE(DMF%XZSPROF)
      ! upper sediment layer
      WHERE (DMF%XZSPROF(IZS) <= F%XH_B1(:))   
         ZCSIS(IZS,:) = DMF%XZSPROF(IZS)/F%XH_B1(:)
         ZSHAPES(IZS,:) = 2*ZCSIS(IZS,:) - ZCSIS(IZS,:)**2
         DMF%XTSPROF(IZS,:) = F%XT_BOT(:) + (F%XT_B1(:) - F%XT_BOT(:))*ZSHAPES(IZS,:)
      END WHERE
      ! lower sediment layer
      WHERE ((F%XDEPTH_BS(:)>F%XH_B1(:)).AND.(DMF%XZSPROF(IZS) >= F%XH_B1(:) .AND. F%XDEPTH_BS(:) >= DMF%XZSPROF(IZS)))
         ZCSIS(IZS,:) = (DMF%XZSPROF(IZS) - F%XH_B1(:))/(F%XDEPTH_BS(:) - F%XH_B1(:))
         ZSHAPES(IZS,:) = 6*ZCSIS(IZS,:)**2 - 8*ZCSIS(IZS,:)**3 + 3*ZCSIS(IZS,:)**4 
         DMF%XTSPROF(IZS,:) = F%XT_B1(:) + (F%XT_BS(:) - F%XT_B1(:))*ZSHAPES(IZS,:)
      END WHERE
   END DO
END IF
!
DMF%XQBOT_FLK   (:) = XUNDEF
DMF%XQW_FLK     (:) = XUNDEF
DMF%XQSNOW_FLK  (:) = XUNDEF
DMF%XQICE_FLK   (:) = XUNDEF
DMF%XIATM_FLK   (:) = XUNDEF
DMF%XISNOW_FLK  (:) = XUNDEF
DMF%XIICE_FLK   (:) = XUNDEF
DMF%XIW_FLK     (:) = XUNDEF
DMF%XIH_FLK     (:) = XUNDEF
DMF%XIINTM0H_FLK(:) = XUNDEF
DMF%XIINTMHD_FLK(:) = XUNDEF
DMF%XQSTAR_FLK  (:) = XUNDEF
DMF%XQSTO_FLK   (:) = XUNDEF
!
IF (DMF%LFLKFLUX) THEN
   DMF%XQBOT_FLK   (:) = PQ_BOT_FLK(:)
   DMF%XQW_FLK     (:) = PQ_W_FLK(:)
   DMF%XQSNOW_FLK  (:) = PQ_SNOW_FLK(:)
   DMF%XQICE_FLK   (:) = PQ_ICE_FLK(:)
   DMF%XIATM_FLK   (:) = PI_ATM_FLK(:)
   DMF%XISNOW_FLK  (:) = PI_SNOW_FLK(:)
   DMF%XIICE_FLK   (:) = PI_ICE_FLK(:)
   DMF%XIW_FLK     (:) = PI_W_FLK(:)
   DMF%XIH_FLK     (:) = PI_H_FLK(:)
   DMF%XIBOT_FLK   (:) = PI_BOT_FLK(:)
   DMF%XIINTM0H_FLK(:) = PI_INTM_0_H_FLK(:)
   DMF%XIINTMHD_FLK(:) = PI_INTM_H_D_FLK(:)
   DMF%XQSTAR_FLK  (:) = PQ_STAR_FLK(:)
   !
   DMF%XQSTO_FLK   (:) = PQ_W_FLK(:)-PQ_BOT_FLK(:)+PI_W_FLK(:)-PI_BOT_FLK(:)
   !
   DMF%XQSTOC_FLK (:) = DMF%XQSTOC_FLK (:) + DMF%XQSTO_FLK (:)*PTSTEP
   DMF%XQBOTC_FLK (:) = DMF%XQBOTC_FLK (:) + DMF%XQBOT_FLK (:)*PTSTEP
   DMF%XIBOTC_FLK (:) = DMF%XIBOTC_FLK (:) + DMF%XIBOT_FLK (:)*PTSTEP
   DMF%XQWC_FLK (:) = DMF%XQWC_FLK (:) + DMF%XQBOT_FLK (:)*PTSTEP
   DMF%XIWC_FLK (:) = DMF%XIWC_FLK (:) + DMF%XIW_FLK (:)*PTSTEP

ENDIF
!
IF (DMF%LFLKWATER) THEN
   DMF%XPRE  (:) = PRAIN(:)+PSNOW(:)
   DMF%XPME  (:) = DMF%XPRE(:)-PEVAP(:)
   !
   DMF%XPMEC (:) = DMF%XPMEC (:) + DMF%XPME (:)*PTSTEP
   DMF%XPREC (:) = DMF%XPREC (:) + DMF%XPRE (:)*PTSTEP
ENDIF
!
IF (LHOOK) CALL DR_HOOK('DIAG_MISC_FLAKE_N',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------------
!
END SUBROUTINE DIAG_MISC_FLAKE_n
