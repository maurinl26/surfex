!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!     #########
      SUBROUTINE WRITESURF_ISBA_LANDUSE_n (HSELECT, IO, S, NP, NPE, KI, HPROGRAM)
!     #####################################
!
!!****  *WRITESURF_ISBA_LANDUSE_n* - writes ISBA prognostic fields
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
!
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
!!      P. LeMoigne 12/2004 : correct dimensionning if more than 10 layers in
!!                            the soil (diffusion version)
!!      B. Decharme  2008    : Floodplains
!!      B. Decharme  01/2009 : Optional Arpege deep soil temperature write
!!      A.L. Gibelin   03/09 : modifications for CENTURY model 
!!      A.L. Gibelin 04/2009 : BIOMASS and RESP_BIOMASS arrays 
!!      A.L. Gibelin 06/2009 : Soil carbon variables for CNT option
!!      B. Decharme  07/2011 : land_use semi-prognostic variables
!!      B. Decharme  09/2012 : suppress NWG_LAYER (parallelization problems)
!!      B. Decharme  09/2012 : write some key for prep_read_external
!!      B. Decharme  04/2013 : Only 2 temperature layer in ISBA-FR
!!      P. Samuelsson 10/2014: MEB
!!      P. Tulet  06/2016 : add XEF et XPFT for MEGAN coupling
!!      M. Leriche 06/2017: comment write XEF & XPFT bug
!!      A. Druel     02/2019 : Add NIRR_TSC and NIRRINUM (with NAG) for irrigation
!!      Séférian/Decharme  08/16  : fire scheme ; change landuse implementation
!!      B. Decharme    02/17 : exact computation of saturation deficit near the leaf surface
!!      B. Decharme    02/21 : explicit soil carbon and gas scheme:browse confirm wa

!!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
!
USE MODN_PREP_SURF_ATM,  ONLY : LWRITE_EXTERN
USE MODD_WRITE_SURF_ATM, ONLY : LSPLIT_PATCH
!
USE MODD_ISBA_OPTIONS_n, ONLY : ISBA_OPTIONS_t
USE MODD_ISBA_n,         ONLY : ISBA_NP_t, ISBA_NPE_t, ISBA_S_t
!
USE MODD_SURF_PAR,       ONLY : LEN_HREC
!
USE MODD_DATA_COVER_PAR, ONLY : NVEGTYPE
!
USE MODI_WRITE_FIELD_1D_PATCH
USE MODI_WRITE_SURF
!
USE YOMHOOK,            ONLY : LHOOK,   DR_HOOK
USE PARKIND1,           ONLY : JPRB
!
IMPLICIT NONE
!
!*       0.1   Declarations of arguments
!              -------------------------
!
CHARACTER(LEN=*), DIMENSION(:), INTENT(IN) :: HSELECT 
!
TYPE(ISBA_OPTIONS_t), INTENT(INOUT) :: IO
TYPE(ISBA_S_t),       INTENT(INOUT) :: S
TYPE(ISBA_NP_t),      INTENT(INOUT) :: NP
TYPE(ISBA_NPE_t),     INTENT(INOUT) :: NPE
INTEGER,              INTENT(IN)    :: KI
!
CHARACTER(LEN=6),    INTENT(IN)    :: HPROGRAM ! program calling
!
!*       0.2   Declarations of local variables
!              -------------------------------
!
!
CHARACTER(LEN=LEN_HREC) :: YRECFM         ! Name of the article to be read
!
CHARACTER(LEN=4 ) :: YLVL
CHARACTER(LEN=3 ) :: YVAR
CHARACTER(LEN=100):: YCOMMENT       ! Comment string
CHARACTER(LEN=25) :: YFORM          ! Writing format
!
INTEGER :: IRESP        ! IRESP  : return-code if a problem appears
INTEGER :: JL, JNC, JP  ! loop counter on levels
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!------------------------------------------------------------------------------
!
!*       2.     Prognostic fields:
!               -----------------
!
IF (LHOOK) CALL DR_HOOK('WRITESURF_ISBA_LANDUSE_N',0,ZHOOK_HANDLE)
!
!
IF(IO%LLULCC.OR.LWRITE_EXTERN)THEN
  DO JL=1,IO%NGROUND_LAYER
     WRITE(YLVL,'(I4)') JL
     YRECFM='DG'//ADJUSTL(YLVL(:LEN_TRIM(YLVL)))
     YFORM='(A6,I1.1,A8)'
     IF (JL >= 10)  YFORM='(A6,I2.2,A8)'
     WRITE(YCOMMENT,FMT=YFORM) 'X_Y_DG',JL,' (m)'
     DO JP = 1,IO%NPATCH
       CALL WRITE_FIELD_1D_PATCH(HSELECT,HPROGRAM,YRECFM,YCOMMENT,JP,&
                      NP%AL(JP)%NR_P,NP%AL(JP)%XDG(:,JL),KI,S%XWORK_WR)    
     ENDDO
  ENDDO
ENDIF
!
!
IF(IO%LLULCC)THEN
  !
  YRECFM='DTLUL'
  YCOMMENT='-'
  CALL WRITE_SURF(HSELECT,HPROGRAM,YRECFM,S%TLULCC,IRESP,HCOMMENT=YCOMMENT)
  !
  DO JL=1,NVEGTYPE
     WRITE(YLVL,'(I4)') JL
     YRECFM='VEGTYPE'//ADJUSTL(YLVL(:LEN_TRIM(YLVL)))
     YCOMMENT='fraction of each vegetation type in the grid cell'
     CALL WRITE_SURF(HSELECT,HPROGRAM,YRECFM,S%XVEGTYPE(:,JL),IRESP,HCOMMENT=YCOMMENT)
  ENDDO
  !
  ! * Water and carbon conservation
  !
  IF(IO%CISBA=='DIF')THEN
    YRECFM = 'WCONSRV'
    YCOMMENT=YRECFM
    CALL WRITE_SURF(HSELECT,HPROGRAM,YRECFM,S%XWCONSRV(:),IRESP,HCOMMENT=YCOMMENT)
  ENDIF
  !
  IF(IO%CPHOTO=='NCB')THEN
    !
    ! * Carbon conservation
    !
    YRECFM = 'CCONSRV'
    YCOMMENT=YRECFM
    CALL WRITE_SURF(HSELECT,HPROGRAM,YRECFM,S%XCCONSRV(:),IRESP,HCOMMENT=YCOMMENT)
    !
    YRECFM = 'FLUATM'
    YCOMMENT=YRECFM
    CALL WRITE_FIELD_GLO2D(S%XFLUATM(:,:))
    !
    YRECFM = 'FLURES'
    YCOMMENT=YRECFM
    CALL WRITE_FIELD_GLO2D(S%XFLURES(:,:))
    !
  ENDIF
  !
  ! * Carbon Managing
  !
  IF(IO%LLULCC_MANAGE)THEN
    !
    YRECFM = 'FLUANT'
    YCOMMENT=YRECFM
    CALL WRITE_FIELD_GLO2D(S%XFLUANT(:,:))
    !
    YRECFM = 'FANTATM'
    YCOMMENT=YRECFM
    CALL WRITE_FIELD_GLO2D(S%XFANTATM(:,:))
    !
    DO JNC=1,IO%NNDECADAL
       !
       WRITE(YLVL,'(I4)') JNC
       YFORM='(A10,I1.1,A10)'
       IF (JNC >= 10)  YFORM='(A10,I2.2,A10)'
       !
       YRECFM='CANTD'//ADJUSTL(YLVL(:LEN_TRIM(YLVL)))
       WRITE(YCOMMENT,FMT=YFORM) 'X_Y_CANTD',JNC,' (kgC/m2)'
       CALL WRITE_FIELD_GLO2D(S%XCSTOCK_DECADAL(:,JNC,:))
       !
       YRECFM='CEXPD'//ADJUSTL(YLVL(:LEN_TRIM(YLVL)))
       WRITE(YCOMMENT,FMT=YFORM) 'X_Y_CEXPD',JNC,' (kgC/m2/yr)'
       CALL WRITE_FIELD_GLO2D(S%XEXPORT_DECADAL(:,JNC,:))
       !
    END DO
    !
    DO JNC=1,IO%NNCENTURY
       !
       WRITE(YLVL,'(I4)') JNC
       YFORM='(A10,I1.1,A10)'
       IF (JNC >= 10)   YFORM='(A10,I2.2,A10)'
       IF (JNC >= 100)  YFORM='(A10,I3.3,A10)'
       !
       YRECFM='CANTC'//ADJUSTL(YLVL(:LEN_TRIM(YLVL)))
       WRITE(YCOMMENT,FMT=YFORM) 'X_Y_CANTC',JNC,' (kgC/m2)'
       CALL WRITE_FIELD_GLO2D(S%XCSTOCK_CENTURY(:,JNC,:))
       !
       YRECFM='CEXPC'//ADJUSTL(YLVL(:LEN_TRIM(YLVL)))
       WRITE(YCOMMENT,FMT=YFORM) 'X_Y_CEXPC',JNC,' (kgC/m2/yr)'
       CALL WRITE_FIELD_GLO2D(S%XEXPORT_CENTURY(:,JNC,:))
       !
    ENDDO
    !
  ENDIF
  !
ENDIF
!
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('WRITESURF_ISBA_LANDUSE_N',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
CONTAINS
!-------------------------------------------------------------------------------
!
SUBROUTINE WRITE_FIELD_GLO2D(PFIELD_IN)
!
IMPLICIT NONE
!
REAL, DIMENSION(:,:), INTENT(IN) :: PFIELD_IN
!
CHARACTER(LEN=LEN_HREC) :: YREC
CHARACTER(LEN=2)        :: YPATCH
!
INTEGER :: IRESP ! IRESP  : return-code if a problem appears
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
IF (LHOOK) CALL DR_HOOK('WRITESURF_ISBA_LANDUSE_N:WRITE_FIELD_GLO2D',0,ZHOOK_HANDLE)
!
IF (LSPLIT_PATCH) THEN
  !
  DO JP = 1,IO%NPATCH
     WRITE(YPATCH,'(I2)') JP
     YREC=ADJUSTL(YRECFM(:LEN_TRIM(YRECFM)))//'P'//ADJUSTL(YPATCH(:LEN_TRIM(YPATCH)))
     CALL WRITE_SURF(HSELECT,HPROGRAM,YREC,PFIELD_IN(:,JP),IRESP,HCOMMENT=YCOMMENT)
  ENDDO
  !
ELSE
  !
  CALL WRITE_SURF(HSELECT,HPROGRAM,YRECFM,PFIELD_IN(:,:),IRESP,HCOMMENT=YCOMMENT)
  !
ENDIF
!
IF (LHOOK) CALL DR_HOOK('WRITESURF_ISBA_LANDUSE_N:WRITE_FIELD_GLO2D',1,ZHOOK_HANDLE)
!
END SUBROUTINE WRITE_FIELD_GLO2D
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE WRITESURF_ISBA_LANDUSE_n
