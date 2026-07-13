!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!     #########
      SUBROUTINE READ_ISBA_PHY_n (IO, S, NP, NPE, NAG, NPGLO, NPEGLO,  &
                                  HPROGRAM, KI, KVERSION, KBUGFIX, ODIM)
!     ##################################
!
!!****  *READ_ISBA_PHY_n* - routine to initialise ISBA physicals variables
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
!!      B. Decharme   *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original     12/2023 : Split from previous read_isban.F90 routine
!!      P. Le Moigne 12/2023 : Roughness length for heat Zilitinkevich
!!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_ISBA_OPTIONS_n, ONLY : ISBA_OPTIONS_t
USE MODD_ISBA_n,         ONLY : ISBA_S_t, ISBA_NP_t, ISBA_NPE_t, ISBA_P_t, ISBA_PE_t
USE MODD_AGRI_n,         ONLY : AGRI_NP_t
!           
USE MODD_SURF_PAR,       ONLY : XUNDEF, NUNDEF, LEN_HREC
USE MODD_CSTS,           ONLY : XG, XRD, XP00
!
USE MODD_AGRI,           ONLY : LIRRIGMODE
!
USE MODE_READ_SURF_LAYERS
!
USE MODI_READ_SURF
USE MODI_MAKE_CHOICE_ARRAY
USE MODI_PACK_SAME_RANK
!
USE MODI_READ_GR_SNOW
USE MODI_ABOR1_SFX
!
USE MODI_IO_BUFF_CLEAN
USE MODE_THERMOS
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*       0.1   Declarations of arguments
!              -------------------------
!
!
TYPE(ISBA_OPTIONS_t), INTENT(INOUT) :: IO
TYPE(ISBA_S_t),       INTENT(INOUT) :: S
TYPE(ISBA_NP_t),      INTENT(INOUT) :: NP, NPGLO
TYPE(ISBA_NPE_t),     INTENT(INOUT) :: NPE, NPEGLO
TYPE(AGRI_NP_t),      INTENT(INOUT) :: NAG
!
CHARACTER(LEN=6),     INTENT(IN)    :: HPROGRAM ! calling program
INTEGER,              INTENT(IN)    :: KI       ! number of points
INTEGER,              INTENT(IN)    :: KVERSION
INTEGER,              INTENT(IN)    :: KBUGFIX
LOGICAL,              INTENT(IN)    :: ODIM
!
!*       0.2   Declarations of local variables
!              -------------------------------
!
TYPE(ISBA_P_t), POINTER  :: PK, PKGLO
TYPE(ISBA_PE_t), POINTER :: PEK, PEKGLO
!
!
CHARACTER(LEN=LEN_HREC) :: YRECFM         ! Name of the article to be read
CHARACTER(LEN=LEN_HREC) :: YCBIO          ! Name of biomass variable
!
CHARACTER(LEN=4)  :: YLVL
!
REAL, DIMENSION(:,:,:),ALLOCATABLE :: ZWORK3D  ! 3D array to write data in file
REAL, DIMENSION(:,:)  ,ALLOCATABLE :: ZWORK2D  ! 2D array to write data in file
REAL, DIMENSION(:)    ,ALLOCATABLE :: ZWORK1D  ! 1D array to write data in file
!
REAL,DIMENSION(IO%NPATCH) :: ZVLAIMIN
REAL :: ZCOEF, ZPS
!
INTEGER :: IRESP   ! Error code after redding
INTEGER :: IWORK   ! Work integer
!
INTEGER :: JI, JP, JL  ! loop counter
!
INTEGER :: ISIZE_LMEB_PATCH
INTEGER :: ITGL
INTEGER :: IMASK
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('READ_ISBA_PHY_N',0,ZHOOK_HANDLE)
!-------------------------------------------------------------------------------
!
!*       1.     Physical dimension
!               -----------------
!
ISIZE_LMEB_PATCH=COUNT(IO%LMEB_PATCH(:))
!
IF(IO%LTEMP_ARP)THEN
  ITGL=IO%NTEMPLAYER_ARP
ELSEIF(IO%CISBA=='DIF')THEN
  ITGL=IO%NGROUND_LAYER
ELSE
  ITGL=2 !Only 2 temperature layer in ISBA-FR
ENDIF
!
!-------------------------------------------------------------------------------
!
CALL ALLOC_READ_ISBA_PHY(NP,NPE)
IF(IO%LLULU)THEN
  CALL ALLOC_READ_ISBA_PHY(NPGLO,NPEGLO)
ENDIF
!
!-------------------------------------------------------------------------------
ALLOCATE(ZWORK2D(KI,IO%NPATCH))
!-------------------------------------------------------------------------------
!
!*       2.     ISBA prognostic fields:
!               -----------------------
!
!* soil temperatures
!
ALLOCATE(ZWORK3D(KI,ITGL,IO%NPATCH))
CALL READ_SURF_LAYERS(HPROGRAM,'TG',ODIM,ZWORK3D,IRESP)
!
DO JL=1,ITGL
  DO JP = 1,IO%NPATCH
     PEK => NPE%AL(JP)
     PK  => NP%AL(JP)
     CALL PACK_SAME_RANK(PK%NR_P,ZWORK3D(:,JL,JP),PEK%XTG(:,JL))
     IF(IO%LLULU) NPEGLO%AL(JP)%XTG(:,JL)=ZWORK3D(:,JL,JP)
  ENDDO
ENDDO
DEALLOCATE(ZWORK3D)
!
!
!* soil liquid water content
!
ALLOCATE(ZWORK3D(KI,IO%NGROUND_LAYER,IO%NPATCH))
CALL READ_SURF_LAYERS(HPROGRAM,'WG',ODIM,ZWORK3D,IRESP)
!
DO JL=1,IO%NGROUND_LAYER
  DO JP = 1,IO%NPATCH
     PEK => NPE%AL(JP)
     PK  => NP%AL(JP)
     CALL PACK_SAME_RANK(PK%NR_P,ZWORK3D(:,JL,JP),PEK%XWG(:,JL))
     IF(IO%LLULU) NPEGLO%AL(JP)%XWG(:,JL)=ZWORK3D(:,JL,JP)
  ENDDO
ENDDO
DEALLOCATE(ZWORK3D)
!
!
!* soil ice content
!
IF(IO%CISBA=='DIF')THEN
  IWORK=IO%NGROUND_LAYER
ELSE
  IWORK=2 !Only 2 soil ice layer in ISBA-FR
ENDIF
!
ALLOCATE(ZWORK3D(KI,IWORK,IO%NPATCH))
CALL READ_SURF_LAYERS(HPROGRAM,'WGI',ODIM,ZWORK3D,IRESP)
!
DO JL=1,IWORK
  DO JP = 1,IO%NPATCH
     PEK => NPE%AL(JP)
     PK  => NP%AL(JP)
     CALL PACK_SAME_RANK(PK%NR_P,ZWORK3D(:,JL,JP),PEK%XWGI(:,JL))
     IF(IO%LLULU) NPEGLO%AL(JP)%XWGI(:,JL)=ZWORK3D(:,JL,JP)
  ENDDO
ENDDO
DEALLOCATE(ZWORK3D)
!
!-------------------------------------------------------------------------------
!
!*       3.     Vegetation prognostic fields:
!               -----------------------------
!
!* water intercepted on leaves
!
YRECFM = 'WR'
CALL MAKE_CHOICE_ARRAY(HPROGRAM, IO%NPATCH, ODIM, YRECFM, ZWORK2D)
!
DO JP = 1,IO%NPATCH
   PEK => NPE%AL(JP)
   PK  => NP%AL(JP)
   CALL PACK_SAME_RANK(PK%NR_P,ZWORK2D(:,JP),PEK%XWR(:))
   IF(IO%LLULU) NPEGLO%AL(JP)%XWR(:)=ZWORK2D(:,JP)
ENDDO
!
!* vegetation canopy air specific humidity
!
IF(KVERSION>=9.OR.(KVERSION==8.AND.ISIZE_LMEB_PATCH>0))THEN
  !
  YRECFM = 'QC'
  CALL MAKE_CHOICE_ARRAY(HPROGRAM, IO%NPATCH, ODIM, YRECFM, ZWORK2D)
  DO JP = 1,IO%NPATCH
     PEK => NPE%AL(JP)
     PK => NP%AL(JP)
     CALL PACK_SAME_RANK(PK%NR_P,ZWORK2D(:,JP),PEK%XQC(:))
     IF(IO%LLULU) NPEGLO%AL(JP)%XQC(:)=ZWORK2D(:,JP)
  ENDDO     
  !
ELSE
  !
  DO JP = 1,IO%NPATCH
     PEK => NPE%AL(JP)
     PK => NP%AL(JP)
     ALLOCATE(ZWORK1D(PK%NSIZE_P))
     CALL PACK_SAME_RANK(PK%NR_P,S%XZS(:),ZWORK1D(:))      
     DO JI=1,PK%NSIZE_P
        IF(PK%XPATCH(JI)>0.0)THEN
          ZPS=XP00*EXP(-(XG/XRD/PEK%XTG(JI,1))*ZWORK1D(JI))
          PEK%XQC(JI)=QSAT(PEK%XTG(JI,1),ZPS)  
        ENDIF
     ENDDO
     DEALLOCATE(ZWORK1D)
  ENDDO
  !
  IF(IO%LLULU)THEN
    DO JP = 1,IO%NPATCH
       PEKGLO => NPEGLO%AL(JP)
       PKGLO  => NPGLO%AL(JP)
       DO JI=1,PKGLO%NSIZE_P
          IF(PEKGLO%XTG(JI,1)/=XUNDEF)THEN
            ZPS=XP00*EXP(-(XG/XRD/PEKGLO%XTG(JI,1))*S%XZS(JI))
            PEKGLO%XQC(JI)=QSAT(PEKGLO%XTG(JI,1),ZPS)  
          ENDIF
       ENDDO
    ENDDO
  ENDIF
  !
ENDIF
!
!* aerodynamical resistance
!
YRECFM = 'RESA'
CALL MAKE_CHOICE_ARRAY(HPROGRAM, IO%NPATCH, ODIM, YRECFM, ZWORK2D)
!
DO JP = 1,IO%NPATCH
   PK => NP%AL(JP)
   PEK => NPE%AL(JP)
   CALL PACK_SAME_RANK(PK%NR_P,ZWORK2D(:,JP),PEK%XRESA(:))
   IF(IO%LLULU) NPEGLO%AL(JP)%XRESA(:)=ZWORK2D(:,JP)
ENDDO
!
!* friction velocity
!
IF (IO%CZ0HEAT=='Z95'.AND.KVERSION>=9) THEN
  YRECFM = 'USTAR'
  CALL MAKE_CHOICE_ARRAY(HPROGRAM, IO%NPATCH, ODIM, YRECFM, ZWORK2D)
  DO JP = 1,IO%NPATCH
     PK => NP%AL(JP)
     PEK => NPE%AL(JP)
     CALL PACK_SAME_RANK(PK%NR_P,ZWORK2D(:,JP),PEK%XUSTAR(:))
  ENDDO
ELSE
  PEK%XUSTAR(:)=0.0
ENDIF
!
!-------------------------------------------------------------------------------
!
!*       4.     Irrigation :
!               ------------
!
!* Irrigation time step counter (current irrigation + time before another irrigation)
!
IF(LIRRIGMODE)THEN
  !
  YRECFM = 'IRR_TSTEP'
  CALL MAKE_CHOICE_ARRAY(HPROGRAM, IO%NPATCH, ODIM, YRECFM, ZWORK2D)
  DO JP = 1,IO%NPATCH
    ALLOCATE(ZWORK1D(SIZE(NAG%AL(JP)%NIRR_TSC(:),1)))
    CALL PACK_SAME_RANK(NP%AL(JP)%NR_P,ZWORK2D(:,JP),ZWORK1D(:))
    NAG%AL(JP)%NIRR_TSC(:)=NINT(ZWORK1D(:))
    DEALLOCATE(ZWORK1D)
  ENDDO
  !
ENDIF
!
!* Irrigation number (from the beguinning of the season)
!
IF(LIRRIGMODE)THEN
  !
  YRECFM='IRR_NUM'
  CALL MAKE_CHOICE_ARRAY(HPROGRAM, IO%NPATCH, ODIM, YRECFM, ZWORK2D)
  DO JP = 1,IO%NPATCH
    ALLOCATE(ZWORK1D(SIZE(NAG%AL(JP)%NIRRINUM(:),1)))
    CALL PACK_SAME_RANK(NP%AL(JP)%NR_P,ZWORK2D(:,JP),ZWORK1D(:))
    NAG%AL(JP)%NIRRINUM(:)=NINT(ZWORK1D(:))
    DEALLOCATE(ZWORK1D)
  ENDDO
  !
ENDIF
!
!-------------------------------------------------------------------------------
!
!*       5.     Snow mantel prognostic fields:
!               ------------------------------
!
DO JP = 1,IO%NPATCH
   PK  => NP%AL(JP)
   PEK => NPE%AL(JP)
   IF(JP>1)THEN
     PEK%TSNOW%SCHEME = NPE%AL(1)%TSNOW%SCHEME
     PEK%TSNOW%NLAYER = NPE%AL(1)%TSNOW%NLAYER
   ENDIF
   CALL READ_GR_SNOW(HPROGRAM,'VEG','     ',KI,PK%NSIZE_P,PK%NR_P,JP,PEK%TSNOW,KNPATCH=IO%NPATCH)
   CALL IO_BUFF_CLEAN
ENDDO
!
IF(IO%LLULU)THEN
  DO JP = 1,IO%NPATCH
     PKGLO  => NPGLO%AL(JP)
     PEKGLO => NPEGLO%AL(JP)
     IF(JP>1)THEN
       PEKGLO%TSNOW%SCHEME = NPEGLO%AL(1)%TSNOW%SCHEME
       PEKGLO%TSNOW%NLAYER = NPEGLO%AL(1)%TSNOW%NLAYER
     ENDIF
     CALL READ_GR_SNOW(HPROGRAM,'VEG','     ',KI,PKGLO%NSIZE_P,PKGLO%NR_P,JP,PEKGLO%TSNOW,KNPATCH=IO%NPATCH)
     CALL IO_BUFF_CLEAN
  ENDDO
ENDIF
!
!-------------------------------------------------------------------------------
!
!*       6.     Glacier prognostic fields:
!               --------------------------
!
IF(IO%LGLACIER)THEN
  !
  DO JP = 1,IO%NPATCH
     PK  => NP%AL(JP)
     PEK => NPE%AL(JP)
     ALLOCATE(PEK%XICE_STO(PK%NSIZE_P))
  ENDDO
  !
  IF (KVERSION>7 .OR. KVERSION==7 .AND. KBUGFIX>=2) THEN
    YRECFM = 'ICE_STO'
    CALL MAKE_CHOICE_ARRAY(HPROGRAM, IO%NPATCH, ODIM, YRECFM, ZWORK2D)
    DO JP = 1,IO%NPATCH
       PK  => NP%AL(JP)
       PEK => NPE%AL(JP)    
       CALL PACK_SAME_RANK(PK%NR_P,ZWORK2D(:,JP),PEK%XICE_STO(:))
    ENDDO
  ELSE
    DO JP = 1,IO%NPATCH
       PEK => NPE%AL(JP)
       PEK%XICE_STO(:) = 0.0
    ENDDO
  ENDIF
  !
ELSE  
  !
  DO JP = 1,IO%NPATCH
     PEK => NPE%AL(JP)
     ALLOCATE(PEK%XICE_STO(0))
  ENDDO
  !
ENDIF
!
!-------------------------------------------------------------------------------
!
!*       7.  Semi-prognostic variables
!            -------------------------
!
!* patch averaged radiative temperature (K)
!
ALLOCATE(S%XTSRAD_NAT(KI))
IF (KVERSION<6) THEN
  S%XTSRAD_NAT(:)=0.
  DO JP=1,IO%NPATCH
    DO JI = 1,NP%AL(JP)%NSIZE_P
      IMASK = NP%AL(JP)%NR_P(JI)
      S%XTSRAD_NAT(IMASK) = S%XTSRAD_NAT(IMASK)+NPE%AL(JP)%XTG(JI,1)
    ENDDO
  ENDDO
  S%XTSRAD_NAT(:)=S%XTSRAD_NAT(:)/IO%NPATCH
ELSE
  YRECFM='TSRAD_NAT'
  CALL READ_SURF(HPROGRAM,YRECFM,S%XTSRAD_NAT(:),IRESP)
ENDIF
!
!-------------------------------------------------------------------------------
DEALLOCATE(ZWORK2D)
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('READ_ISBA_PHY_N',1,ZHOOK_HANDLE)
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
CONTAINS
!-------------------------------------------------------------------------------
!
SUBROUTINE ALLOC_READ_ISBA_PHY(NA,NAE)
!
IMPLICIT NONE
!
TYPE(ISBA_NP_t),      INTENT(INOUT) :: NA
TYPE(ISBA_NPE_t),     INTENT(INOUT) :: NAE
!
TYPE(ISBA_P_t),  POINTER :: PA
TYPE(ISBA_PE_t), POINTER :: PEA
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
IF (LHOOK) CALL DR_HOOK('READ_ISBA_PHY_N:ALLOC_READ_ISBA_PHY',0,ZHOOK_HANDLE)
!
DO JP = 1,IO%NPATCH
   !
   PEA => NAE%AL(JP)
   PA => NA%AL(JP)
   !
   ALLOCATE(PEA%XTG (PA%NSIZE_P,ITGL            ))
   ALLOCATE(PEA%XWG (PA%NSIZE_P,IO%NGROUND_LAYER))
   ALLOCATE(PEA%XWGI(PA%NSIZE_P,IO%NGROUND_LAYER))
   !
   PEA%XTG (:,:) = XUNDEF
   PEA%XWG (:,:) = XUNDEF
   PEA%XWGI(:,:) = XUNDEF
   !
   ALLOCATE(PEA%XWR   (PA%NSIZE_P))
   ALLOCATE(PEA%XQC   (PA%NSIZE_P))
   ALLOCATE(PEA%XRESA (PA%NSIZE_P))
   ALLOCATE(PEA%XUSTAR(PA%NSIZE_P))
   PEA%XWR   (:) = XUNDEF
   PEA%XQC   (:) = XUNDEF
   PEA%XRESA (:) = 100.
   PEA%XUSTAR(:) = 0.
   !
ENDDO
!
IF (LHOOK) CALL DR_HOOK('READ_ISBA_PHY_N:ALLOC_READ_ISBA_PHY',1,ZHOOK_HANDLE)
!
END SUBROUTINE ALLOC_READ_ISBA_PHY
!
!-------------------------------------------------------------------------------
!
SUBROUTINE ALLOC_READ_MEB(NA,NAE)
!
IMPLICIT NONE
!
TYPE(ISBA_NP_t),      INTENT(INOUT) :: NA
TYPE(ISBA_NPE_t),     INTENT(INOUT) :: NAE
!
TYPE(ISBA_P_t),  POINTER :: PA
TYPE(ISBA_PE_t), POINTER :: PEA
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
IF (LHOOK) CALL DR_HOOK('READ_ISBA_PHY_N:ALLOC_READ_MEB',0,ZHOOK_HANDLE)
!
DO JP = 1,IO%NPATCH
   !
   PEA => NAE%AL(JP)
   PA => NA%AL(JP)
   !
   ALLOCATE(PEA%XWRL (PA%NSIZE_P))
   ALLOCATE(PEA%XWRLI(PA%NSIZE_P))
   ALLOCATE(PEA%XWRVN(PA%NSIZE_P))
   ALLOCATE(PEA%XTV  (PA%NSIZE_P))
   ALLOCATE(PEA%XTL  (PA%NSIZE_P))
   ALLOCATE(PEA%XTC  (PA%NSIZE_P))
   !
   PEA%XWRL (:) = XUNDEF
   PEA%XWRLI(:) = XUNDEF
   PEA%XWRVN(:) = XUNDEF
   PEA%XTV  (:) = XUNDEF
   PEA%XTL  (:) = XUNDEF
   PEA%XTC  (:) = XUNDEF
   !
ENDDO

IF (LHOOK) CALL DR_HOOK('READ_ISBA_PHY_N:ALLOC_READ_MEB',1,ZHOOK_HANDLE)
!
END SUBROUTINE ALLOC_READ_MEB
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE READ_ISBA_PHY_n
