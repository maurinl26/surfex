!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!     #########
      SUBROUTINE WRITE_DIAG_MISC_FLAKE_n ( DTCO, HSELECT, U, DMF, HPROGRAM)
!     #################################
!
!!****  *WRITE_DIAG_MISC_FLAKE* - writes the FLAKE miscellaneous diagnostic fields
!!
!!    PURPOSE
!!    -------
!!
!!
!!**  METHOD
!!    ------
!!
!!    REFERENCE
!!    ---------
!!
!!
!!    AUTHOR
!!    ------
!!      P. Le Moigne   *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original     10/2004
!!      P. Le Moigne 05/2023 : temperature profile in sediments
!!      P. Le Moigne 08/2023 : heat storage flux
!!      P. Le Moigne 11/2023 : effective P-E
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_DATA_COVER_n, ONLY : DATA_COVER_t
USE MODD_SURF_ATM_n, ONLY : SURF_ATM_t
USE MODD_DIAG_MISC_FLAKE_n, ONLY : DIAG_MISC_FLAKE_t
USE MODD_SURF_PAR, ONLY: LEN_HREC
!
USE MODI_INIT_IO_SURF_n
USE MODI_WRITE_SURF
USE MODI_END_IO_SURF_n
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*       0.1   Declarations of arguments
!              -------------------------
!
TYPE(DATA_COVER_t), INTENT(INOUT) :: DTCO
 CHARACTER(LEN=*), DIMENSION(:), INTENT(IN) :: HSELECT
TYPE(SURF_ATM_t), INTENT(INOUT) :: U
TYPE(DIAG_MISC_FLAKE_t), INTENT(INOUT) :: DMF
!
 CHARACTER(LEN=6),  INTENT(IN)  :: HPROGRAM ! program calling
!
!*       0.2   Declarations of local variables
!              -------------------------------
!
CHARACTER(LEN=LEN_HREC) :: YRECFM   ! Name of the article to be read
CHARACTER(LEN=100) :: YCOMMENT      ! Comment string
INTEGER :: IRESP                    ! IRESP  : return-code if a problem appears
INTEGER :: IZ                       ! loop index
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('WRITE_DIAG_MISC_FLAKE_N',0,ZHOOK_HANDLE)
!
!         Initialisation for IO
!
CALL INIT_IO_SURF_n(DTCO, U, HPROGRAM,'WATER ','FLAKE ','WRITE','FLAKE_DIAGNOSTICS.OUT.nc')
!
!-------------------------------------------------------------------------------
!
!* Flake temperature profile
!
IF (DMF%LWATER_PROFILE) THEN      
   DO IZ=1,SIZE(DMF%XZWPROF)
      WRITE(YRECFM,'(F5.2)') DMF%XZWPROF(IZ)
      YRECFM='TWAT_'//TRIM(ADJUSTL(YRECFM))
      YCOMMENT='X_Y_'//YRECFM//' (K)'
      CALL WRITE_SURF(HSELECT, HPROGRAM,YRECFM,DMF%XTWPROF(IZ,:),IRESP,HCOMMENT=YCOMMENT)
   END DO
END IF
!
!* Flake temperature profile in sediment
!
IF (DMF%LSEDIM_PROFILE) THEN
   DO IZ=1,SIZE(DMF%XZSPROF)
      WRITE(YRECFM,'(F5.2)') DMF%XZSPROF(IZ)
      YRECFM='TSED_'//TRIM(ADJUSTL(YRECFM))
      YCOMMENT='X_Y_'//YRECFM//' (K)'
      CALL WRITE_SURF(HSELECT, HPROGRAM,YRECFM,DMF%XTSPROF(IZ,:),IRESP,HCOMMENT=YCOMMENT)
   END DO
END IF
!
!* Flake heat and radiative fluxes diagnostics
!
IF (DMF%LFLKFLUX) THEN
   YRECFM='QBOT_FLK'
   YCOMMENT='Heat flux through the water-bottom sediment interface [W m^{-2}]'
   CALL WRITE_SURF(HSELECT, HPROGRAM,YRECFM,DMF%XQBOT_FLK(:),IRESP,HCOMMENT=YCOMMENT)
   YRECFM='QW_FLK'
   YCOMMENT='Heat flux through the ice-water or air-water interface [W m^{-2}]'
   CALL WRITE_SURF(HSELECT, HPROGRAM,YRECFM,DMF%XQW_FLK(:),IRESP,HCOMMENT=YCOMMENT)
   YRECFM='QSNOW_FLK'
   YCOMMENT='Heat flux through the air-snow interface [W m^{-2}]'
   CALL WRITE_SURF(HSELECT, HPROGRAM,YRECFM,DMF%XQSNOW_FLK(:),IRESP,HCOMMENT=YCOMMENT)
   YRECFM='QICE_FLK'
   YCOMMENT='Heat flux through the snow-ice or air-ice interface [W m^{-2}]'
   CALL WRITE_SURF(HSELECT, HPROGRAM,YRECFM,DMF%XQICE_FLK(:),IRESP,HCOMMENT=YCOMMENT)
   YRECFM='IATM_FLK'
   YCOMMENT='Radiation flux at the lower boundary of the atmosphere [W m^{-2}]'
   CALL WRITE_SURF(HSELECT, HPROGRAM,YRECFM,DMF%XIATM_FLK(:),IRESP,HCOMMENT=YCOMMENT)
   YRECFM='ISNOW_FLK'
   YCOMMENT='Radiation flux through the air-snow interface [W m^{-2}]'
   CALL WRITE_SURF(HSELECT, HPROGRAM,YRECFM,DMF%XISNOW_FLK(:),IRESP,HCOMMENT=YCOMMENT)
   YRECFM='IICE_FLK'
   YCOMMENT='Radiation flux through the snow-ice or air-ice interface [W m^{-2}]'
   CALL WRITE_SURF(HSELECT, HPROGRAM,YRECFM,DMF%XIICE_FLK(:),IRESP,HCOMMENT=YCOMMENT)
   YRECFM='IW_FLK'
   YCOMMENT='Radiation flux through the ice-water or air-water interface [W m^{-2}]'
   CALL WRITE_SURF(HSELECT, HPROGRAM,YRECFM,DMF%XIW_FLK(:),IRESP,HCOMMENT=YCOMMENT)
   YRECFM='IH_FLK'
   YCOMMENT='Radiation flux through the mixed-layer-thermocline interface [W m^{-2}]'
   CALL WRITE_SURF(HSELECT, HPROGRAM,YRECFM,DMF%XIH_FLK(:),IRESP,HCOMMENT=YCOMMENT)
   YRECFM='IBOT_FLK'
   YCOMMENT='Radiation flux through the water-bottom sediment interface [W m^{-2}]'
   CALL WRITE_SURF(HSELECT, HPROGRAM,YRECFM,DMF%XIBOT_FLK(:),IRESP,HCOMMENT=YCOMMENT)
   YRECFM='IINTM0H_FLK'
   YCOMMENT='Mean radiation flux over the mixed layer [W m^{-1}]'
   CALL WRITE_SURF(HSELECT, HPROGRAM,YRECFM,DMF%XIINTM0H_FLK(:),IRESP,HCOMMENT=YCOMMENT)
   YRECFM='IINTMHD_FLK'
   YCOMMENT='Mean radiation flux over the thermocline [W m^{-1}]'
   CALL WRITE_SURF(HSELECT, HPROGRAM,YRECFM,DMF%XIINTMHD_FLK(:),IRESP,HCOMMENT=YCOMMENT)
   YRECFM='QSTAR_FLK'
   YCOMMENT='A generalized heat flux scale [W m^{-2}]'
   CALL WRITE_SURF(HSELECT, HPROGRAM,YRECFM,DMF%XQSTAR_FLK(:),IRESP,HCOMMENT=YCOMMENT)
   YRECFM='QSTO_FLK'
   YCOMMENT=' Heat storage flux [W m^{-2}]'
   CALL WRITE_SURF(HSELECT, HPROGRAM,YRECFM,DMF%XQSTO_FLK(:),IRESP,HCOMMENT=YCOMMENT)
ENDIF
!
IF (DMF%LFLKWATER) THEN
   YRECFM='PME_FLK'
   YCOMMENT='effective rain P-E over lake [kg m^{-2} s^{-1}]'
   CALL WRITE_SURF(HSELECT, HPROGRAM,YRECFM,DMF%XPME(:),IRESP,HCOMMENT=YCOMMENT)
   YRECFM='PMEC_FLK'
   YCOMMENT='accumulated effective rain P-E over lake [kg m^{-2}]'
   CALL WRITE_SURF(HSELECT, HPROGRAM,YRECFM,DMF%XPMEC(:),IRESP,HCOMMENT=YCOMMENT)
   YRECFM='PRE_FLK'
   YCOMMENT='total precipitation over lake [kg m^{-2} s^{-1}]'
   CALL WRITE_SURF(HSELECT, HPROGRAM,YRECFM,DMF%XPRE(:),IRESP,HCOMMENT=YCOMMENT)
   YRECFM='PREC_FLK'
   YCOMMENT='accumulated total precipitation over lake [kg m^{-2}]'
   CALL WRITE_SURF(HSELECT, HPROGRAM,YRECFM,DMF%XPREC(:),IRESP,HCOMMENT=YCOMMENT)
   YRECFM='QSTOC_FLK'
   YCOMMENT='accumulated heat storage [J m^{-2}]'
   CALL WRITE_SURF(HSELECT, HPROGRAM,YRECFM,DMF%XQSTOC_FLK(:),IRESP,HCOMMENT=YCOMMENT)
   YRECFM='IBOTC_FLK'
   YCOMMENT='accumulated radiation flux through the water-bottom sediment interface [J m^{-2}]'
   CALL WRITE_SURF(HSELECT, HPROGRAM,YRECFM,DMF%XIBOTC_FLK(:),IRESP,HCOMMENT=YCOMMENT)
   YRECFM='QBOTC_FLK'
   YCOMMENT='accumulated heat flux through the water-bottom sediment interface [J m^{-2}]'
   CALL WRITE_SURF(HSELECT, HPROGRAM,YRECFM,DMF%XQBOTC_FLK(:),IRESP,HCOMMENT=YCOMMENT)
   YRECFM='IWC_FLK'
   YCOMMENT='accumulated radiation flux through the ice/air-water sediment interface [J m^{-2}]'
   CALL WRITE_SURF(HSELECT, HPROGRAM,YRECFM,DMF%XIWC_FLK(:),IRESP,HCOMMENT=YCOMMENT)
   YRECFM='QWC_FLK'
   YCOMMENT='accumulated heat flux through the ice/air-water sediment interface [J m^{-2}]'
   CALL WRITE_SURF(HSELECT, HPROGRAM,YRECFM,DMF%XQWC_FLK(:),IRESP,HCOMMENT=YCOMMENT)
ENDIF
!
!-------------------------------------------------------------------------------
!
!         End of IO
!
CALL END_IO_SURF_n(HPROGRAM)
!
IF (LHOOK) CALL DR_HOOK('WRITE_DIAG_MISC_FLAKE_N',1,ZHOOK_HANDLE)
!
END SUBROUTINE WRITE_DIAG_MISC_FLAKE_n
