!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!     #########
    SUBROUTINE LANDUSE_SOILGAS(IO, S, NP, NPE, NPGLO, NPEGLO, KI)
!   ###############################################################
!!****  *LANDUSE_SOILGAS*
!!
!!    PURPOSE
!!    -------
!
!     Performs land use land cover change computation at yearly time step
!               
!!**  METHOD
!!    ------
!!
!!    EXTERNAL
!!    --------
!!    none
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      
!!    none
!!
!!    REFERENCE
!!    ---------
!!
!!      
!!    AUTHOR
!!    ------
!!    B. Decharme 08/2023
!!
!!
!-------------------------------------------------------------------------------
!
!*       0.     DECLARATIONS
!               ------------
!
USE MODD_ISBA_OPTIONS_n, ONLY : ISBA_OPTIONS_t
USE MODD_ISBA_n,         ONLY : ISBA_S_t, ISBA_P_t, ISBA_PE_t, &
                                ISBA_NP_t, ISBA_NPE_t
!
USE MODD_SURF_PAR,       ONLY : XUNDEF,NUNDEF
!
USE MODD_CSTS,           ONLY : XTT, XG, XLMTT, XRHOLW
USE MODD_ISBA_PAR,       ONLY : XWGMIN, XWTD_MAXDEPTH
!
USE YOMHOOK   ,          ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,          ONLY : JPRB
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
TYPE(ISBA_OPTIONS_t),  INTENT(INOUT) :: IO
TYPE(ISBA_S_t),        INTENT(INOUT) :: S
TYPE(ISBA_NP_t),       INTENT(INOUT) :: NP, NPGLO
TYPE(ISBA_NPE_t),      INTENT(INOUT) :: NPE, NPEGLO
!
INTEGER,               INTENT(IN)    :: KI
!
!*      0.2    declarations of local arguments
!
TYPE(ISBA_P_t),  POINTER :: PK, PKGLO
TYPE(ISBA_PE_t), POINTER :: PEK, PEKGLO
!
INTEGER, DIMENSION(KI,IO%NPATCH) :: IWG_LAYER ! Number of hydrological layers
!
!
! working table
!
REAL     :: ZWORK1, ZWORK2, ZWORK3, &
            ZLOG, ZWTOT, ZWL,       &
            ZMATPOT, ZMATPOTN
!
INTEGER  :: INL, INP, IDEPTH, IMASK
INTEGER  :: JI, JL, JP
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!-----------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('LANDUSE_SOILGAS',0,ZHOOK_HANDLE)
!-----------------------------------------------------------------
!
!*      1.     Preliminaries
!              -------------
!
INL=IO%NGROUND_LAYER
INP=IO%NPATCH
!
! 1.1 local arguments 
!
IWG_LAYER(:,:) = 0
!
!-----------------------------------------------------------------
!
!*      2.     Compute previous year total soilgas content (kg m-2)
!              ---------------------------------------------------------------
! 
WRITE(*,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
WRITE(*,*)'!!!                                                        !!!'
WRITE(*,*)'!!!                    WARNING    WARNING                  !!!'
WRITE(*,*)'!!!                                                        !!!'
WRITE(*,*)'!!!  Land-use Land Cover Change computation are performed  !!!'
WRITE(*,*)'!!!                                                        !!!'
WRITE(*,*)'!!! SOIL GAS CONSERVATION COMPUTATION NOT YET IMPLEMENTED  !!!'
WRITE(*,*)'!!!                                                        !!!'
WRITE(*,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
!
!
!-----------------------------------------------------------------
!
!*      3.     Compute current year water profile
!              -------------------------
DO JP=1,INP
   PK => NP%AL(JP)
   PEK => NPE%AL(JP)
   DO JL=1,INL
      DO JI=1,PK%NSIZE_P
         IF(PK%NWG_LAYER(JI)/=NUNDEF.AND.PEK%XSGASO2(JI,JL)/=XUNDEF)THEN
           IWG_LAYER(JI,JP) = JL
         ENDIF
      ENDDO
   ENDDO
ENDDO
!
DO JP=1,INP
   !
   PK => NP%AL(JP)
   PEK => NPE%AL(JP)
   !
   DO JL=1,INL
      !
      DO JI=1,PK%NSIZE_P
         !
         ! 3.1 Check consitency
         !
         IF(PK%NWG_LAYER(JI)/=NUNDEF.AND.JL<=PK%NWG_LAYER(JI))THEN
           !
           ! 3.2 Ensure coherent soil gas profile
           ! Method: ???
           !
           IF(PEK%XSGASO2(JI,JL)==XUNDEF) THEN
             !
             IDEPTH = IWG_LAYER(JI,JP)
             !
             PEK%XSGASO2 (JI,JL) = PEK%XSGASO2 (JI,IDEPTH)
             PEK%XSGASCO2(JI,JL) = PEK%XSGASCO2(JI,IDEPTH)
             PEK%XSGASCH4(JI,JL) = PEK%XSGASCH4(JI,IDEPTH)
             !
           ENDIF
           !
         ELSE
           !
           PEK%XSGASO2 (JI,JL) = XUNDEF
           PEK%XSGASCO2(JI,JL) = XUNDEF
           PEK%XSGASCH4(JI,JL) = XUNDEF
           !
         ENDIF
         !
      ENDDO
      !
   ENDDO
   !
ENDDO
!
!-----------------------------------------------------------------
!
!*      4.     Compute total land soil gas storage (kg m-2)
!              -------------------------------------------------
! 
!
!-----------------------------------------------------------------
!
!*      5.    Land-use induced soil gas mass (kg m-2)
!              -----------------------------------
! 
!
!-----------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('LANDUSE_SOILGAS',1,ZHOOK_HANDLE)
!-----------------------------------------------------------------
!
END SUBROUTINE LANDUSE_SOILGAS
