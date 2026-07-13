!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!     #########
      SUBROUTINE READ_ISBA_ASSIM_n (IO, K, NP, NPE, HPROGRAM, KI, KVERSION, KBUGFIX, ODIM)
!     ##################################
!
!!****  *READ_ISBA_ASSIM_n* - routine to initialise ISBA physicals variables
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
!!      Original    12/2023 Split from previous read_isban.F90 routine
!!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_ISBA_OPTIONS_n, ONLY : ISBA_OPTIONS_t
USE MODD_ISBA_n,         ONLY : ISBA_K_t, ISBA_NP_t, ISBA_NPE_t, ISBA_P_t, ISBA_PE_t
!                
USE MODD_ASSIM,          ONLY : LASSIM,CASSIM_ISBA,XAT2M_ISBA,XAHU2M_ISBA,&
                                XAZON10M_ISBA,XAMER10M_ISBA,NIFIC,NVAR, &
                                COBS,NOBSTYPE,CVAR,LPRT,XTPRT,NIVAR,CBIO, &
                                XADDINFL,NENS,XSIGMA,NIE
!                                
USE MODD_SURF_PAR,       ONLY : XUNDEF, LEN_HREC
!
USE MODI_READ_SURF
USE MODI_MAKE_CHOICE_ARRAY
USE MODI_PACK_SAME_RANK
!
USE MODI_ABOR1_SFX
USE MODI_IO_BUFF
!
USE MODE_RANDOM
USE MODE_EKF
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*       0.1   Declarations of arguments
!              -------------------------
!
TYPE(ISBA_OPTIONS_t), INTENT(INOUT) :: IO
TYPE(ISBA_K_t),       INTENT(INOUT) :: K
TYPE(ISBA_NP_t),      INTENT(INOUT) :: NP
TYPE(ISBA_NPE_t),     INTENT(INOUT) :: NPE
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
TYPE(ISBA_P_t), POINTER  :: PK
TYPE(ISBA_PE_t), POINTER :: PEK
!
!
CHARACTER(LEN=LEN_HREC) :: YRECFM         ! Name of the article to be read
CHARACTER(LEN=LEN_HREC) :: YCBIO          ! Name of biomass variable
!
CHARACTER(LEN=2)  :: YNB
!
REAL, DIMENSION(:,:,:),ALLOCATABLE :: ZLAI
REAL, DIMENSION(:)    ,ALLOCATABLE :: ZCOFSWI
!
REAL,DIMENSION(IO%NPATCH) :: ZVLAIMIN
!
REAL    :: ZCOEF
!
INTEGER :: IRESP          ! Error code after redding
!
INTEGER :: IWORK   ! Work integer
!
INTEGER :: JP, JL, JNB, JNL  ! loop counter on layers
INTEGER :: JVAR, JI, JNCANT
!
INTEGER :: IOBS
INTEGER :: IMASK
!
LOGICAL :: GKNOWN
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('READ_ISBA_ASSIM_N',0,ZHOOK_HANDLE)
!-------------------------------------------------------------------------------
!
!*       1.     Initialise
!               ----------
!
IF ( TRIM(CASSIM_ISBA)=="ENKF") THEN
  DO JP = 1,IO%NPATCH
    PK => NP%AL(JP)
    ALLOCATE(PK%XRED_NOISE(PK%NSIZE_P,NVAR))
    PK%XRED_NOISE(:,:) = 0.
  ENDDO
ELSE
  DO JP = 1,IO%NPATCH
    ALLOCATE(NP%AL(JP)%XRED_NOISE(0,0))
  ENDDO  
ENDIF
!
IF  ( TRIM(CASSIM_ISBA)=="ENKF" .OR. (TRIM(CASSIM_ISBA)=="EKF" .AND. LPRT) ) THEN
  ALLOCATE(ZCOFSWI(KI))
  CALL COFSWI(K%XCLAY(:,1),ZCOFSWI)
ELSE
  ALLOCATE(ZCOFSWI(0))
ENDIF
!
!-------------------------------------------------------------------------------
!
!*       2.     ISBA prognostic fields:
!               -----------------------
!
!* soil temperatures
!
IF(IO%LTEMP_ARP)THEN
  IWORK=IO%NTEMPLAYER_ARP
ELSEIF(IO%CISBA=='DIF')THEN
  IWORK=IO%NGROUND_LAYER
ELSE
  IWORK=2 !Only 2 temperature layer in ISBA-FR
ENDIF
!
! Perturb value if requested
!
IF ( TRIM(CASSIM_ISBA)=="EKF" .AND. LPRT ) THEN
  !
  DO JL=1,IWORK
  ! read in control variable
    IF ( (TRIM(CVAR(NIVAR))=="TG1" .AND. JL==1) .OR. (TRIM(CVAR(NIVAR))=="TG2" .AND. JL==2) ) THEN
      DO JP = 1,IO%NPATCH
        PEK => NPE%AL(JP)
        WHERE ( PEK%XTG(:,JL)/=XUNDEF )
          PEK%XTG(:,JL) = PEK%XTG(:,JL) + XTPRT(NIVAR)*PEK%XTG(:,JL)
        ENDWHERE
      ENDDO
    ENDIF
  END DO
  !
ELSEIF ( TRIM(CASSIM_ISBA)=="ENKF" .AND. NIE<NENS+1 ) THEN
  !
  CALL MAKE_ENS_ENKF(IWORK,KI,"TG ",ZCOFSWI,NP)
  !
ENDIF
!
!
!* soil liquid and ice water contents
!
! Perturb value if requested
IF ( TRIM(CASSIM_ISBA)=="EKF" .AND. LPRT ) THEN
   !
   DO JL=1,IO%NGROUND_LAYER
     ! read in control variable
     IF ( (TRIM(CVAR(NIVAR))=="WG1" .AND. JL==1) .OR. & 
          (TRIM(CVAR(NIVAR))=="WG2" .AND. JL==2) .OR. &
          (TRIM(CVAR(NIVAR))=="WG3" .AND. JL==3) .OR. &
          (TRIM(CVAR(NIVAR))=="WG4" .AND. JL==4) .OR. &
          (TRIM(CVAR(NIVAR))=="WG5" .AND. JL==5) .OR. &
          (TRIM(CVAR(NIVAR))=="WG6" .AND. JL==6) .OR. &
          (TRIM(CVAR(NIVAR))=="WG7" .AND. JL==7) .OR. &
          (TRIM(CVAR(NIVAR))=="WG8" .AND. JL==8) ) THEN     
       !
       DO JP = 1,IO%NPATCH
          PEK => NPE%AL(JP)
          PK => NP%AL(JP)
          DO JI = 1,PK%NSIZE_P
             IMASK = PK%NR_P(JI)
             IF (PEK%XWG(JI,JL)/=XUNDEF ) THEN
                PEK%XWG(JI,JL) = PEK%XWG(JI,JL) + XTPRT(NIVAR) * ZCOFSWI(IMASK) 
             ENDIF
          ENDDO
       END DO
       !
     ENDIF
     !
   END DO
   !
ELSEIF ( TRIM(CASSIM_ISBA)=="ENKF" .AND. NIE<NENS+1 ) THEN
  !
  CALL MAKE_ENS_ENKF(IWORK,KI,"WG ",ZCOFSWI,NP)
  !
ENDIF
!
!-------------------------------------------------------------------------------
!
!*       3.     Leaf Area Index:
!               ----------------
!
IF (IO%CPHOTO=='NIT' .OR. IO%CPHOTO=='NCB') THEN
  !
  IF ( TRIM(CASSIM_ISBA)=="EKF" .AND. LPRT ) THEN
    !
    ! read in control variable
    !
    IF ( TRIM(CVAR(NIVAR))=="LAI" ) THEN
      !
      DO JP = 1,IO%NPATCH
        PEK => NPE%AL(JP)
        WHERE ( PEK%XLAI(:)/=XUNDEF ) 
           PEK%XLAI(:) =  PEK%XLAI(:) + XTPRT(NIVAR)* PEK%XLAI(:)
        ENDWHERE
      ENDDO
      !
    ENDIF
    !
  ELSEIF ( TRIM(CASSIM_ISBA)=="ENKF" .AND. NIE<NENS+1 ) THEN
    !
    IF (IO%NPATCH==12) THEN
      ZVLAIMIN = (/0.3,0.3,0.3,0.3,1.0,1.0,0.3,0.3,0.3,0.3,0.3,0.3/)
    ELSE
      ZVLAIMIN = (/0.3/)
    ENDIF
    !
    ALLOCATE(ZLAI(KI,1,IO%NPATCH))
    !
    ZLAI(:,:,:) = 0.
    !
    DO JP = 1,IO%NPATCH
      PEK => NPE%AL(JP)
      PK  => NP%AL(JP)
      ZLAI(1:PK%NSIZE_P,1,JP) = PEK%XLAI(:)
    ENDDO
    !
    CALL MAKE_ENS_ENKF(1,KI,"LAI",ZCOFSWI,NP,ZLAI)
    !
    DO JP = 1,IO%NPATCH
      PEK => NPE%AL(JP)
      PK  => NP%AL(JP)
      PEK%XLAI(:) = MAX(ZVLAIMIN(JP),ZLAI(1:PK%NSIZE_P,1,JP))
    ENDDO
    !
    DEALLOCATE(ZLAI)
    !    
  ENDIF  
END IF
!
!-------------------------------------------------------------------------------
!
!*       4.  ISBA-AGS
!           ---------
!
IF (IO%CPHOTO=='NIT'.OR.IO%CPHOTO=='NCB') THEN
  !
  IF (KVERSION>7 .OR. KVERSION==7 .AND. KBUGFIX>=3) THEN
    YRECFM='BIOMA'
  ELSE
    YRECFM='BIOMASS'
  ENDIF
  !
  DO JNB=1,IO%NNBIOMASS
     !
     WRITE(YNB,'(I1)') JNB
     ! 
     IF ( TRIM(CASSIM_ISBA)=="EKF" .AND. LPRT ) THEN
        !
        YCBIO = YRECFM(:LEN_TRIM(YRECFM))//ADJUSTL(YNB(:LEN_TRIM(YNB)))
        !
        !read in control variable
        !
        IF ( TRIM(CVAR(NIVAR)) == "LAI" .AND. TRIM(CBIO)==TRIM(YCBIO) ) THEN
           DO JP = 1,IO%NPATCH
              PK => NP%AL(JP)
              PEK => NPE%AL(JP)           
              DO JI = 1,PK%NSIZE_P
                 IF(PEK%XBIOMASS(JI,JNB)/=XUNDEF)THEN
                   PEK%XBIOMASS(JI,JNB) = PEK%XBIOMASS(JI,JNB) * ( 1. + XTPRT(NIVAR) )
                 ENDIF
              ENDDO
           ENDDO
        ENDIF
        !
     ELSEIF ( TRIM(CASSIM_ISBA)=="ENKF" .AND. NIE<NENS+1 .AND. .NOT.LASSIM ) THEN
        !
        IF( TRIM(CBIO)==TRIM(YRECFM) ) THEN
          DO JVAR = 1,NVAR
             IF (TRIM(CVAR(JVAR)) == "LAI") THEN
                DO JP = 1,IO%NPATCH
                   PK => NP%AL(JP)
                   PEK => NPE%AL(JP)           
                   DO JI = 1,PK%NSIZE_P
                      PEK%XBIOMASS(JI,JNB) = PEK%XBIOMASS(JI,JNB) + XADDINFL(JVAR)*RANDOM_NORMAL()
                   ENDDO
                ENDDO
                EXIT
             ENDIF
          ENDDO
        ENDIF
        !
    ENDIF     
    !
  ENDDO
  !
  IWORK=0
  IF(IO%CPHOTO=='NCB'.OR.KVERSION<8)THEN
    IWORK=2
  ENDIF
  !
  DO JNB=2,IO%NNBIOMASS-IWORK
     !
     IF ( TRIM(CASSIM_ISBA)=="EKF" .AND. LPRT ) THEN
        !
        WRITE(YNB,'(I1)') JNB
        !
        IF (KVERSION>7 .OR. (KVERSION==7 .AND. KBUGFIX>=3)) THEN
           YRECFM='RESPI'//ADJUSTL(YNB(:LEN_TRIM(YNB)))
        ELSE
           YRECFM='RESP_BIOM'//ADJUSTL(YNB(:LEN_TRIM(YNB)))
        ENDIF    
        !        
        !read in control variable
        IF ( TRIM(CVAR(NIVAR)) == "LAI" .AND. TRIM(CBIO)==TRIM(YRECFM) ) THEN
           !
           DO JP = 1,IO%NPATCH
              PK => NP%AL(JP)
              PEK => NPE%AL(JP)           
              DO JI = 1,PK%NSIZE_P
                 IF(PEK%XRESP_BIOMASS(JI,JNB)/=XUNDEF)THEN
                    PEK%XRESP_BIOMASS(JI,JNB) = PEK%XRESP_BIOMASS(JI,JNB) * XTPRT(NIVAR)
                 ENDIF
              ENDDO
           ENDDO
           !
        ELSEIF ( TRIM(CASSIM_ISBA)=="ENKF" .AND. NIE<NENS+1 .AND. .NOT.LASSIM ) THEN
           !
           IF ( TRIM(CBIO)==TRIM(YRECFM) ) THEN
              DO JVAR = 1,NVAR
                 IF (TRIM(CVAR(JVAR)) == "LAI") THEN
                    DO JP = 1,IO%NPATCH
                       PK => NP%AL(JP)
                       PEK => NPE%AL(JP)           
                       DO JI = 1,PK%NSIZE_P
                          PEK%XRESP_BIOMASS(JI,JNB) = PEK%XRESP_BIOMASS(JI,JNB) + XADDINFL(JVAR)*RANDOM_NORMAL()
                       ENDDO
                    ENDDO
                    EXIT
                 ENDIF
              ENDDO
           ENDIF
           !
        ENDIF
        !
   ENDIF  
   !
  ENDDO
  !
ENDIF
!
DEALLOCATE(ZCOFSWI)
!
!-------------------------------------------------------------------------------
!
!*       5. Near surface assimilation
!           -------------------------
!
IF ( LASSIM ) THEN
  IF ( TRIM(CASSIM_ISBA) == "OI" ) THEN
    IF ( IO%NPATCH /= 1 ) CALL ABOR1_SFX ('Reading of diagnostical values for'&
                       & //'assimilation at the moment only works for one patch for OI')          
    ! Diagnostic fields for assimilation
    IF ( .NOT. ALLOCATED(XAT2M_ISBA)) ALLOCATE(XAT2M_ISBA(KI,1))
    XAT2M_ISBA=XUNDEF
    YRECFM='T2M'
    CALL IO_BUFF(YRECFM,'R',GKNOWN)
    CALL READ_SURF(HPROGRAM,YRECFM,XAT2M_ISBA(:,1),IRESP)

    IF ( .NOT. ALLOCATED(XAHU2M_ISBA)) ALLOCATE(XAHU2M_ISBA(KI,1))
    XAHU2M_ISBA=XUNDEF
    YRECFM='HU2M'
    CALL IO_BUFF(YRECFM,'R',GKNOWN)
    CALL READ_SURF(HPROGRAM,YRECFM,XAHU2M_ISBA(:,1),IRESP)

    IF ( .NOT. ALLOCATED(XAZON10M_ISBA)) ALLOCATE(XAZON10M_ISBA(KI,1))
    XAZON10M_ISBA=XUNDEF
    YRECFM='ZON10M'
    CALL IO_BUFF(YRECFM,'R',GKNOWN)
    CALL READ_SURF(HPROGRAM,YRECFM,XAZON10M_ISBA(:,1),IRESP)

    IF ( .NOT. ALLOCATED(XAMER10M_ISBA)) ALLOCATE(XAMER10M_ISBA(KI,1))
    XAMER10M_ISBA=XUNDEF
    YRECFM='MER10M'
    CALL IO_BUFF(YRECFM,'R',GKNOWN)
    CALL READ_SURF(HPROGRAM,YRECFM,XAMER10M_ISBA(:,1),IRESP)
  ELSEIF ( NIFIC/=NVAR+2 ) THEN
    ! Diagnostic fields for EKF assimilation ("observations")
    DO IOBS = 1,NOBSTYPE
     SELECT CASE (TRIM(COBS(IOBS)))
       CASE("T2M")
         IF ( .NOT. ALLOCATED(XAT2M_ISBA)) ALLOCATE(XAT2M_ISBA(KI,1))
         XAT2M_ISBA=XUNDEF
         YRECFM='T2M'
         CALL IO_BUFF(YRECFM,'R',GKNOWN)
         CALL READ_SURF(HPROGRAM,YRECFM,XAT2M_ISBA(:,1),IRESP)
       CASE("HU2M")
         IF ( .NOT. ALLOCATED(XAHU2M_ISBA)) ALLOCATE(XAHU2M_ISBA(KI,1))
         XAHU2M_ISBA=XUNDEF
         YRECFM='HU2M'
         CALL IO_BUFF(YRECFM,'R',GKNOWN)
         CALL READ_SURF(HPROGRAM,YRECFM,XAHU2M_ISBA(:,1),IRESP)
       CASE("WG1")
         ! This is already read above
       CASE("WG2")
         ! This is already read above
       CASE("LAI")
         ! This is already read above   
       CASE("SWE")
         ! This is handled independently 
       CASE DEFAULT
         CALL ABOR1_SFX("Mapping of "//TRIM(COBS(IOBS))//" is not defined in READ_ISBA_ASSIM_n!")
     END SELECT
    ENDDO
  ENDIF
ENDIF
!
!-------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('READ_ISBA_ASSIM_N',1,ZHOOK_HANDLE)
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
CONTAINS
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
!
SUBROUTINE MAKE_ENS_ENKF(KWORK,KLU,HREC,PCOFSWI,NP,PVAR)
!
USE MODD_ASSIM, ONLY : LENS_GEN, XADDTIMECORR, XADDINFL, XASSIM_WINH
!
USE MODI_ADD_NOISE
USE MODE_RANDOM
!
IMPLICIT NONE
!
INTEGER, INTENT(IN) :: KWORK
INTEGER, INTENT(IN) :: KLU
CHARACTER(LEN=3), INTENT(IN) :: HREC
REAL, DIMENSION(:), INTENT(IN) :: PCOFSWI
TYPE(ISBA_NP_t), INTENT(INOUT) :: NP
REAL, DIMENSION(:,:,:), INTENT(INOUT), OPTIONAL :: PVAR
!
REAL, DIMENSION(:,:) ,ALLOCATABLE :: ZWORK2D  ! 2D array
!
CHARACTER(LEN=LEN_HREC) :: YRECFM         ! Name of the article to be read
CHARACTER(LEN=4) :: YLVL
CHARACTER(LEN=3) :: YVAR
REAL, DIMENSION(KLU) :: ZVAR
REAL :: ZWHITE_NOISE, ZVAR0
INTEGER :: JL, JI, JP, IVAR
LOGICAL :: GPASS
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
IF (LHOOK) CALL DR_HOOK('READ_ISBA_ASSIM_N:MAKE_ENS_ENKF',0,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
ALLOCATE(ZWORK2D(KI,IO%NPATCH))
!-------------------------------------------------------------------------------
!
DO JL=1,KWORK
  !
  IF (KWORK>1) THEN
    WRITE(YLVL,'(I4)') JL
    YRECFM = TRIM(HREC)//ADJUSTL(YLVL(:LEN_TRIM(YLVL)))
  ELSE
    YRECFM = TRIM(HREC)
  ENDIF
  !
  IVAR = 0
  DO JVAR = 1,NVAR
    GPASS = ( TRIM(CVAR(JVAR))==TRIM(YRECFM) )
    IF (GPASS) THEN
      IVAR = JVAR
      EXIT
    ENDIF
  ENDDO
  !
  IF ( GPASS ) THEN
    !
    IF (XADDINFL(IVAR)>0.) THEN
      !
      IF (LASSIM .OR. (.NOT.LENS_GEN .AND. XADDTIMECORR(IVAR)>0.)) THEN
        !
        WRITE(YVAR,'(I3)') IVAR
        YRECFM='RD_NS'//ADJUSTL(YVAR(:LEN_TRIM(YVAR)))
        CALL MAKE_CHOICE_ARRAY(HPROGRAM, IO%NPATCH, ODIM, YRECFM, ZWORK2D)
        DO JP = 1,IO%NPATCH
          CALL PACK_SAME_RANK(NP%AL(JP)%NR_P,ZWORK2D(:,JP),NP%AL(JP)%XRED_NOISE(:,IVAR))
        ENDDO           
        !
        IF (.NOT.LASSIM) THEN
          !
          DO JP = 1,IO%NPATCH
            PK => NP%AL(JP)
            DO JI = 1,NP%AL(JP)%NSIZE_P
              IMASK = PK%NR_P(JI)
              ZWHITE_NOISE = XADDINFL(IVAR)*PCOFSWI(IMASK)*RANDOM_NORMAL()
              CALL ADD_NOISE(XADDTIMECORR(IVAR),XASSIM_WINH,ZWHITE_NOISE,PK%XRED_NOISE(JI,IVAR))
            ENDDO
          ENDDO
          !
          ZCOEF = XASSIM_WINH/24.
          !
        ENDIF
        !
      ELSE
        !
        DO JP = 1,IO%NPATCH
          PK => NP%AL(JP)
          DO JI = 1,NP%AL(JP)%NSIZE_P
            IMASK = PK%NR_P(JI)
            NP%AL(JP)%XRED_NOISE(JI,IVAR) = XADDINFL(IVAR)*PCOFSWI(IMASK)*RANDOM_NORMAL()
          ENDDO
        ENDDO
        !
        ZCOEF = 1. 
        !
      ENDIF
      !
      IF (.NOT.LASSIM) THEN
        !
        DO JP = 1,IO%NPATCH
          !
          ZVAR(:) = 0.
          IF (TRIM(HREC)=='TG') THEN
            ZVAR(1:NP%AL(JP)%NSIZE_P) = NPE%AL(JP)%XTG(:,JL)
          ELSEIF (TRIM(HREC)=='WG') THEN
            ZVAR(1:NP%AL(JP)%NSIZE_P) = NPE%AL(JP)%XWG(:,JL)
          ELSEIF (TRIM(HREC)=='LAI' .AND. PRESENT(PVAR)) THEN
            ZVAR(1:NP%AL(JP)%NSIZE_P) = PVAR(1:NP%AL(JP)%NSIZE_P,JL,JP)
          ELSE
            CALL ABOR1_SFX("READ_ISBA_ASSIM_N: HREC "//HREC//" not permitted")
          ENDIF
          !
          DO JI = 1,NP%AL(JP)%NSIZE_P
            IF ( ZVAR(JI)/=XUNDEF ) THEN
              !
              ZVAR0 = ZVAR(JI)
              !
              ZVAR(JI) = ZVAR(JI) + ZCOEF * NP%AL(JP)%XRED_NOISE(JI,IVAR)
              !
              IF (ZVAR(JI) < 0.) THEN
                IF (LENS_GEN) THEN
                  ZVAR(JI) = ABS(ZVAR(JI))
                ELSE
                  ZVAR(JI) = ZVAR0
                ENDIF
              ENDIF
            ENDIF
          ENDDO
          !
          IF (TRIM(HREC)=='TG') THEN
            NPE%AL(JP)%XTG(:,JL) = ZVAR(1:NP%AL(JP)%NSIZE_P)
          ELSEIF (TRIM(HREC)=='WG') THEN
            NPE%AL(JP)%XWG(:,JL) = ZVAR(1:NP%AL(JP)%NSIZE_P)
          ELSEIF (TRIM(HREC)=='LAI') THEN
            PVAR(1:NP%AL(JP)%NSIZE_P,JL,JP) = ZVAR(1:NP%AL(JP)%NSIZE_P)
          ENDIF
          !
        ENDDO
        !
      ENDIF
      !
    ENDIF
    !
  ENDIF
  !
ENDDO
!
!-------------------------------------------------------------------------------
DEALLOCATE(ZWORK2D)
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('READ_ISBA_ASSIM_N:MAKE_ENS_ENKF',1,ZHOOK_HANDLE)
!
END SUBROUTINE MAKE_ENS_ENKF
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE READ_ISBA_ASSIM_n
