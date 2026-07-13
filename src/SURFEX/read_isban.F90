!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!SFX_LIC for details. version 1.
!     #########
      SUBROUTINE READ_ISBA_n (DTCO, U, IO, S, K, NK, NP, NPE, NAG, NPGLO, NPEGLO, HPROGRAM, HINIT)
!     ##################################
!
!!****  *READ_ISBA_n* - routine to initialise ISBA variables
!!                         
!!
!!    PURPOSE
!!    -------
!!
!!**  METHOD
!!    ------
!!
!!    EXTERNAL
!!    --------
!!
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
!!      V. Masson   *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    01/2003
!!
!!      READ_SURF for general reading : 08/2003 (S.Malardel)
!!      B. Decharme  2008    : Floodplains
!!      B. Decharme  01/2009 : Optional Arpege deep soil temperature read
!!      A.L. Gibelin   03/09 : modifications for CENTURY model 
!!      A.L. Gibelin    04/2009 : BIOMASS and RESP_BIOMASS arrays 
!!      A.L. Gibelin    06/2009 : Soil carbon variables for CNT option
!!      B. Decharme     09/2012 : suppress NWG_LAYER (parallelization problems)
!!      T. Aspelien     08/2013 : Read diagnostics for assimilation
!!      P. Samuelsson   10/2014 : MEB
!!      A. Druel        02/2019 : streamlines the code, and add TSC and NIRRINUM for irrigation
!!     Séférian/Decharme 08/2016 : fire, carbone leaching and landuse module 
!!      R. Séférian 11/2016 : add cmip6 diagnostics
!!      B. Decharme    02/17 : exact computation of saturation deficit near the leaf surface
!!      B. Decharme 04/2020 : New soil carbon scheme (Morel et al. 2019 JAMES) under CRESPSL = DIF option
!!      B. Decharme 12/2023 : global fields for Land use land cover change (LLULU)
!!      B. Decharme 12/2023 : spliting of previous read_isban.F90 in several routines
!!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_DATA_COVER_n,   ONLY : DATA_COVER_t
USE MODD_SURF_ATM_n,     ONLY : SURF_ATM_t
USE MODD_ISBA_OPTIONS_n, ONLY : ISBA_OPTIONS_t
USE MODD_ISBA_n,         ONLY : ISBA_S_t, ISBA_K_t, ISBA_NK_t, ISBA_NP_t, ISBA_NPE_t, ISBA_P_t, ISBA_PE_t
USE MODD_AGRI_n,         ONLY : AGRI_NP_t
!
USE MODD_SURF_PAR,       ONLY : LEN_HREC
USE MODD_ASSIM,          ONLY : LASSIM,CASSIM_ISBA
!
USE MODI_GET_TYPE_DIM_n
USE MODI_READ_SURF
!
USE MODI_READ_ISBA_PHY_n
USE MODI_READ_ISBA_MEB_n
USE MODI_READ_ISBA_AGS_n
USE MODI_READ_ISBA_CC_n
USE MODI_READ_ISBA_LANDUSE_n
USE MODI_READ_ISBA_NUDGING_n
USE MODI_READ_ISBA_ASSIM_n
!
USE MODI_ABOR1_SFX
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*       0.1   Declarations of arguments
!              -------------------------
!
TYPE(DATA_COVER_t),   INTENT(INOUT) :: DTCO
TYPE(SURF_ATM_t),     INTENT(INOUT) :: U
TYPE(ISBA_OPTIONS_t), INTENT(INOUT) :: IO
TYPE(ISBA_S_t),       INTENT(INOUT) :: S
TYPE(ISBA_K_t),       INTENT(INOUT) :: K
TYPE(ISBA_NK_t),      INTENT(INOUT) :: NK
TYPE(ISBA_NP_t),      INTENT(INOUT) :: NP, NPGLO
TYPE(ISBA_NPE_t),     INTENT(INOUT) :: NPE, NPEGLO
TYPE(AGRI_NP_t),      INTENT(INOUT) :: NAG
!
CHARACTER(LEN=6),     INTENT(IN)    :: HPROGRAM ! calling program
CHARACTER(LEN=3),     INTENT(IN)    :: HINIT    ! choice of fields to initialize
!
!*       0.2   Declarations of local variables
!              -------------------------------
!
CHARACTER(LEN=LEN_HREC) :: YRECFM           ! Name of the article to be read
INTEGER                 :: ILU              ! 1D physical dimension
INTEGER                 :: ISIZE_LMEB_PATCH ! MEB key
INTEGER                 :: IRESP            ! Error code after redding
INTEGER                 :: IVERSION         ! surfex version
INTEGER                 :: IBUGFIX          ! surfex nersion bugfix
LOGICAL                 :: GDIM
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('READ_ISBA_N',0,ZHOOK_HANDLE)
!-------------------------------------------------------------------------------
!
YRECFM='SIZE_NATURE'
CALL GET_TYPE_DIM_n(DTCO, U, 'NATURE',ILU)
!
YRECFM='VERSION'
CALL READ_SURF(HPROGRAM,YRECFM,IVERSION,IRESP)
!
YRECFM='BUG'
CALL READ_SURF(HPROGRAM,YRECFM,IBUGFIX,IRESP)
!
GDIM = (IVERSION>8 .OR. IVERSION==8 .AND. IBUGFIX>0)
IF (GDIM) CALL READ_SURF(HPROGRAM,'SPLIT_PATCH',GDIM,IRESP)
!
!-------------------------------------------------------------------------------
!
!* Read ISBA physical prognostic and semi-prognostic fields
!
CALL READ_ISBA_PHY_n(IO, S, NP, NPE, NAG, NPGLO, NPEGLO,   &
                     HPROGRAM, ILU, IVERSION, IBUGFIX, GDIM)
!
!
IF (HINIT=='PRE'.AND.IO%CISBA=='DIF'.AND.NPE%AL(1)%TSNOW%SCHEME.NE.'3-L' &
                                    .AND.NPE%AL(1)%TSNOW%SCHEME.NE.'CRO' ) THEN
   CALL ABOR1_SFX("READ_ISBA_N: WITH CISBA = DIF, CSNOW MUST BE 3-L OR CRO")
ENDIF
!
!-------------------------------------------------------------------------------
!
!* Read ISBA physical prognostic and semi-prognostic fields
!
ISIZE_LMEB_PATCH=COUNT(IO%LMEB_PATCH(:))
!
IF(ISIZE_LMEB_PATCH>0)THEN
  CALL READ_ISBA_MEB_n(IO, NP, NPE, NPGLO, NPEGLO, &
                       HPROGRAM, ILU, GDIM         )
ENDIF
!
!-------------------------------------------------------------------------------
!
!* Read ISBA Ags prognostic and semi-prognostic fields
!
IF(IO%CPHOTO/='NON')THEN
  CALL READ_ISBA_AGS_n(IO, NP, NPE, NPGLO, NPEGLO,            &
                       HPROGRAM, ILU, IVERSION, IBUGFIX, GDIM )
ENDIF
!
!-------------------------------------------------------------------------------
!
!* Read ISBA Carbon Cycle prognostic and semi-prognostic fields
!
IF(IO%CRESPSL=='CNT'.OR.IO%CRESPSL=='DIF'.OR.IO%LFIRE)THEN
  CALL READ_ISBA_CC_n(IO, NP, NPE, NPGLO, NPEGLO,   &
                      HPROGRAM, ILU, IVERSION, GDIM )
ENDIF
!
!-------------------------------------------------------------------------------
!
!* Read ISBA landuse prognostic and semi-prognostic fields
!
IF(IO%LLULCC)THEN
  CALL READ_ISBA_LANDUSE_n(IO, S, HPROGRAM, ILU, IVERSION, GDIM)
ENDIF
!
!-------------------------------------------------------------------------------
!
!* Assimilation applied to ISBA prognostic and semi-prognostic fields
!
IF(LASSIM.OR.TRIM(CASSIM_ISBA)/="OI")THEN
  CALL READ_ISBA_ASSIM_n(IO, K, NP, NPE, HPROGRAM, ILU, IVERSION, IBUGFIX, GDIM)
ENDIF
!
!-------------------------------------------------------------------------------
!
!* Soil water and/or snow nudging fields
!
IF(HINIT=='ALL'.AND.(IO%CNUDG_WG/='DEF'.OR.IO%LNUDG_SWE))THEN
  CALL READ_ISBA_NUDGING_n(IO, S, K, NK, NP, NPGLO, HPROGRAM, ILU, GDIM)
ENDIF
!
!-------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('READ_ISBA_N',1,ZHOOK_HANDLE)
!-------------------------------------------------------------------------------
!
END SUBROUTINE READ_ISBA_n
