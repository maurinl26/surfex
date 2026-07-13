!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!#############################################################
SUBROUTINE LANDUSE_NUDGING (IO, S, NK, NP, NPE, SOLD, NPGLO, HPROGRAM, KI)  
!#############################################################
!
!!****  *LANDUSE_NUDGING* - routine to initialize land use for ISBA field
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
!!      B. Decharme   *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    07/2011
!!      Completelly reframed 08/2016 R. Séférian
!!      R. Séférian 10/2016 correct error in landuse computation fields
!!      R. Séférian 11/2016 : add cmip6 diagnostics
!!      J. Colin    12/2017 : add computations in case the water or snow is
!!                            nudged seperately on each patch
!!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_ISBA_OPTIONS_n, ONLY : ISBA_OPTIONS_t
USE MODD_ISBA_n,         ONLY : ISBA_S_t, ISBA_K_t, ISBA_P_t, ISBA_PE_t, ISBA_NK_t, ISBA_NP_t, ISBA_NPE_t
USE MODD_DATA_ISBA_n,    ONLY : DATA_ISBA_t
!
USE MODD_SURF_PAR,       ONLY : XUNDEF, LEN_HREC
!
USE MODD_DATA_COVER_PAR, ONLY : NVEGTYPE
!
USE MODI_GET_LUOUT
USE MODI_LANDUSE_HYDRO_NUDGING 
!
USE MODI_ATTRIBUTE_CLOSEST_VEGTYPE_NUDGING
!
USE MODI_UNPACK_SAME_RANK
!
USE YOMHOOK,  ONLY : LHOOK,   DR_HOOK
USE PARKIND1, ONLY : JPRB
!
IMPLICIT NONE
!
!*       0.1   Declarations of arguments
!              -------------------------
!
TYPE(ISBA_OPTIONS_t),  INTENT(INOUT) :: IO
TYPE(ISBA_S_t),        INTENT(INOUT) :: S, SOLD
TYPE(ISBA_NK_t),       INTENT(INOUT) :: NK
TYPE(ISBA_NP_t),       INTENT(INOUT) :: NP, NPGLO
TYPE(ISBA_NPE_t),      INTENT(INOUT) :: NPE
!
CHARACTER(LEN=6),      INTENT(IN)    :: HPROGRAM          ! program calling surf. schemes
INTEGER,               INTENT(IN)    :: KI
!
!
!*       0.2   Declarations of local arguments on Patch grid
!
!
TYPE(ISBA_K_t),  POINTER :: KK
TYPE(ISBA_P_t),  POINTER :: PK, PNEAR
TYPE(ISBA_PE_t), POINTER :: PEK
!
!*       0.3   Declarations of local arguments on complete ISBA grid
!
!
REAL, DIMENSION(KI,IO%NPATCH) :: ZWORK          ! work array (ISBA grid)
CHARACTER(LEN=LEN_HREC)       :: YRECFM         ! Name of the article to be read
CHARACTER(LEN=4)              :: YLVL
INTEGER                       :: IRESP          ! Error code after redding
!
LOGICAL                       :: GLULU          ! Logical to perform luluccf computation of not
INTEGER                       :: ILUOUT         ! unit of output listing file
INTEGER                       :: JI, JL         ! loop counter
INTEGER                       :: JP, JP_NEAR    ! loop counter
INTEGER                       :: INP, INL       ! dimension
INTEGER                       :: JT, INTIME     ! loop on time (nudging) and size
INTEGER                       :: IGLO           ! Work integer
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('LANDUSE_NUDGING',0,ZHOOK_HANDLE)
!-------------------------------------------------------------------------------
!        
!
!*       0. Global initialisation
!        ------------------------
!
!
CALL GET_LUOUT(HPROGRAM,ILUOUT)
!
INP   =IO%NPATCH
INL   =IO%NGROUND_LAYER
!
INTIME=3
!
IF(IO%CNUDG_WG=='DAY'.OR.IO%LNUDG_SWE)THEN
  SELECT CASE (S%TTIME%TDATE%MONTH)
     CASE(4,6,9,11)
      INTIME=30
     CASE(1,3,5,7:8,10,12)
      INTIME=31
     CASE(2)
      INTIME=29
  END SELECT
ENDIF
!
!
!-------------------------------------------------------------------------------
!        
!
!*       2. Soil moisture nudging
!        ------------------------
!
IF(IO%CNUDG_WG/='DEF')THEN
  !
  ZWORK(:,:)=XUNDEF
  DO JP=1,INP
     PK => NP%AL(JP)
     CALL UNPACK_SAME_RANK(PK%NR_P,PK%XNUDG_WGTOT(:,1,1),ZWORK(:,JP))
  ENDDO
  !
  DO JP=1,INP
     !
     KK => NK%AL(JP)
     PK => NP%AL(JP)
     PEK => NPE%AL(JP)
     !
     GLULU=ANY((PEK%XWG(:,1)/=XUNDEF).AND.(PK%XNUDG_WGTOT(:,1,1)==XUNDEF))
     !
     IF(GLULU)THEN  ! If any new patch has appeared, change the XNUDG_WGTOT climatology
       !
       WRITE(ILUOUT,*)'--------------------------------------------------------------'
       WRITE(ILUOUT,*)'!!!              Patches Distribution has changed          !!!' 
       WRITE(ILUOUT,*)'!!!  Land-use nudging update computation are performed     !!!'
       WRITE(ILUOUT,*)'--------------------------------------------------------------'
       !
       DO JI=1,PK%NSIZE_P
          !
          ! Where there is a change
          IF(PEK%XWG(JI,1)/=XUNDEF.AND.PK%XNUDG_WGTOT(JI,1,1)==XUNDEF)THEN
            !
            IGLO = PK%NR_P(JI)
            CALL ATTRIBUTE_CLOSEST_VEGTYPE_NUDGING(INP,NVEGTYPE,SOLD%XPATCH(IGLO,:),ZWORK(IGLO,:),JP,JP_NEAR) 
            !
            PNEAR => NPGLO%AL(JP_NEAR)
            !
            !* Total water content
            !
            DO JL=1,INL
               DO JT=1,INTIME
                  PK%XNUDG_WGTOT(JI,JL,JT) = PNEAR%XNUDG_WGTOT(IGLO,JL,JT)
               ENDDO
            ENDDO
            !
            DO JL=1,INL
               DO JT=1,INTIME
                  IF(PK%XNUDG_WGTOT(JI,JL,JT)/=XUNDEF.AND.PK%XNUDG_WGTOT(JI,JL,JT)>KK%XWSAT(JI,JL))THEN
                     PK%XNUDG_WGTOT(JI,JL,JT) = KK%XWSAT(JI,JL)
                  ENDIF
               ENDDO
            ENDDO
            !
            ! Correct (complete) water content profile in the new soil layers
            !
            IF(IO%CISBA=='DIF')THEN
              CALL LANDUSE_HYDRO_NUDGING(IO, NK, NP, KI)                             
            ENDIF
            !
          ENDIF ! If a new patch arises
          !
       ENDDO ! JI
       !      
     ENDIF ! If there is any change
     !
  ENDDO ! JP
  !
ENDIF ! If CNUDG_WG/='DEF'
!
!
!-------------------------------------------------------------------------------
!        
!
!*       3. Snow mass nudging
!        --------------------
!
IF(IO%LNUDG_SWE)THEN
  !
  ZWORK(:,:)=XUNDEF
  DO JP=1,INP
     PK => NP%AL(JP)
     CALL UNPACK_SAME_RANK(PK%NR_P,PK%XNUDG_SWE(:,1),ZWORK(:,JP))
  ENDDO
  !
  DO JP=1,INP
     !
     PK => NP%AL(JP)
     PEK => NPE%AL(JP)
     !
     GLULU=ANY((PEK%XWG(:,1)/=XUNDEF).AND.(PK%XNUDG_SWE(:,1)==XUNDEF))
     !
     IF(GLULU)THEN  ! If any new patch has appeared, change the XNUDG_SWE climatology
       !
       WRITE(ILUOUT,*)'--------------------------------------------------------------'
       WRITE(ILUOUT,*)'!!!              Patches Distribution has changed          !!!' 
       WRITE(ILUOUT,*)'!!!  Land-use nudging update computation are performed     !!!'
       WRITE(ILUOUT,*)'--------------------------------------------------------------'
       !
       DO JI=1,KI
          !
          ! Where there is a change
          IF(PEK%XWG(JI,1)/=XUNDEF.AND.PK%XNUDG_SWE(JI,1)==XUNDEF) THEN
            !
            IGLO = PK%NR_P(JI)
            !
            CALL ATTRIBUTE_CLOSEST_VEGTYPE_NUDGING(INP,NVEGTYPE,SOLD%XPATCH(IGLO,:),ZWORK(IGLO,:),JP,JP_NEAR) 
            !
            PNEAR => NPGLO%AL(JP_NEAR)
            !
            !* Total snow water content
            !
            DO JT=1,INTIME
               PK%XNUDG_SWE(JI,JT) = PNEAR%XNUDG_SWE(IGLO,JT)
            ENDDO        
            !
          ENDIF ! If a new patch arises
          !
       ENDDO ! JI
       !
     ENDIF ! If there is any change
     !
  ENDDO ! JP
!
ENDIF ! If LNUDG_SWE_PATCH
!
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('LANDUSE_NUDGING',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
END SUBROUTINE LANDUSE_NUDGING
