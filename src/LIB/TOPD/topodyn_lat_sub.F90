!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!-----------------------------------------------------------------
!     ###########################
SUBROUTINE TOPODYN_LAT_SUB(PRW,PDEF,PKAPPA,PKAPPAC,OTOPD)
!     ###########################
!
!
!     PURPOSE
!     -------
!     to distribute laterally soil water following topodyn concept
!
!
!     METHOD
!     ------
!
!     EXTERNAL
!     --------
!     none
!
!
!     AUTHOR
!     ------
!
!     G.-M. Saulnier * LTHE * 
!     K. Chancibault * CNRM *
!
!     MODIFICATIONS
!     -------------
!
!     Original    12/2003
!     writing in fortran 90 12/2004
!------------------------------------------------------------------------------------------
!
!*    0.0    DECLARATIONS
!            ------------
USE MODD_TOPODYN,       ONLY : NNCAT,NNCAT_SUB, NMESHT, XDMAXT, XDXT,&
                               XMPARA, NNMC, XCONN, NLINE,&
                               XSLOP,  XDAREA, XLAMBDA, NCAT_CAT_TO_SUB,&
                               NPIX_CAT_TO_SUB,NPIX_SUB_TO_CAT,&
                               NBPIX_SUB, NNCAT_MAX,NNMC
USE MODD_COUPLING_TOPD, ONLY : XWSTOPT, XDTOPT
USE MODD_TOPD_PAR,        ONLY : XSTEPK
!
USE MODD_SURF_PAR,        ONLY : XUNDEF,NUNDEF
!
USE MODI_FLOWDOWN_SUB
USE MODI_ABOR1_SFX
USE MODI_WRITE_FILE_VECMAP
USE MODI_WRITE_FILE_MAP
!
USE YOMHOOK   ,ONLY : LHOOK, DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*    0.1   declarations of arguments
!
REAL, DIMENSION(:,:), INTENT(IN)   :: PRW
REAL, DIMENSION(:,:), INTENT(OUT)  :: PDEF
REAL, DIMENSION(:,:), INTENT(OUT)  :: PKAPPA
REAL, DIMENSION(:), INTENT(OUT)    :: PKAPPAC
LOGICAL, DIMENSION(:), INTENT(OUT) :: OTOPD
!
!*    0.2   declarations of local variables
!
LOGICAL              :: GFOUND  ! logical variable
REAL                 :: ZSOMME
REAL                 :: ZM      ! XMPARA in m
REAL                 :: ZDX     ! XDXT in m
REAL                 :: ZKVAL, ZKVALMIN, ZKVALMAX
REAL                 :: ZDAV    ! Averaged deficit (m)
REAL                 :: ZDAV2   ! Averaged deficit on ZA-ZAS-ZAD (m)
REAL                 :: ZNDMAXAV,ZNKAV ! temporary averaged maximal deficit and averaged similarity index
REAL                 :: ZDMAXAV,ZKAV   ! averaged maximal deficit and averaged similarity index
REAL                 :: ZFUNC
REAL                 :: ZDIF,ZDIFMIN   ! difference calcul
REAL                 :: ZNAS, ZNAD     ! temporary saturated and dry relative catchment area
REAL                 :: ZAS,ZAD        ! saturated and dry relative catchment area 
REAL                 :: ZTMP
!
REAL, DIMENSION(NMESHT) :: ZDMAX     ! XDMAXT in m
REAL, DIMENSION(NMESHT) :: ZRW       ! PRW in m
REAL, DIMENSION(NMESHT) :: ZDINI     ! initial deficit
REAL, DIMENSION(NMESHT) :: ZMASK
REAL, DIMENSION(NMESHT) :: ZKAPPA_PACK, ZDMAX_PACK
REAL, DIMENSION(NNCAT_MAX,NMESHT)  :: ZDEF
REAL, DIMENSION(NNCAT_MAX,NMESHT)  :: ZKAPPA
REAL, DIMENSION(NNCAT_MAX)     :: ZKAPPAC
LOGICAL, DIMENSION(NNCAT_MAX)  :: GTOPD
!
INTEGER              :: J1, J2, JJ,JCAT,JSB !
INTEGER              :: INKAPPA ! number of steps in similarity index distribution
INTEGER              :: INPCON  ! number of connected pixels
INTEGER              :: INAS    ! number of saturated pixels
INTEGER              :: INAD    ! number of dry pixels
INTEGER :: I_DIM
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-----------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('TOPODYN_LAT_SUB',0,ZHOOK_HANDLE)
!
PKAPPA(:,:)= 0.0
ZKAPPA(:,:)= 0.0
!
PKAPPAC(:) = 0.
ZKAPPAC(:) = 0.
!
OTOPD(:) = .TRUE.
GTOPD(:) = .TRUE.
!
WHERE ((XDMAXT(:,:)-PRW(:,:)>=0.).AND.XDMAXT(:,:)-PRW(:,:)<=XDMAXT(:,:))
  PDEF(:,:) = XDMAXT(:,:)-PRW(:,:)
ENDWHERE
WHERE(XDMAXT(:,:)-PRW(:,:)>XDMAXT(:,:))
  PDEF(:,:) = XDMAXT(:,:)
ENDWHERE
WHERE(XDMAXT(:,:)-PRW(:,:)<0.)
  PDEF(:,:) = 0
ENDWHERE
ZDEF(:,:)= XUNDEF
ZKVALMIN = XUNDEF
ZKVALMAX = XUNDEF
!
DO JCAT = 1,NNCAT
  DO JSB=1,NNCAT_SUB(JCAT)
! Changement de variable => passage sur les sbv
    ZRW(:)=XUNDEF
    ZDMAX(:)=XUNDEF
    DO JJ=1,NNMC(JCAT)
      IF (JSB==NCAT_CAT_TO_SUB(JCAT,JJ) .AND. JSB/=NUNDEF .AND. NPIX_CAT_TO_SUB(JCAT,JJ)/=NUNDEF) THEN
        ZRW(NPIX_CAT_TO_SUB(JCAT,JJ))   = PRW(JCAT,JJ)
        ZDMAX(NPIX_CAT_TO_SUB(JCAT,JJ)) = XDMAXT(JCAT,JJ)
        ZDEF(JSB,NPIX_CAT_TO_SUB(JCAT,JJ)) = PDEF(JCAT,JJ)
      ENDIF
    ENDDO
   !*    0.    Initialisation
   !           -------------- 
   ZMASK(:) = 0.0
   INPCON   = 0
   ZDAV     = 0.0
   ZDAV2    = 0.0
   GFOUND   = .FALSE.
   ZDIFMIN  = 99999.
   !
   ZDX = XDXT(JCAT)
   ZM = XMPARA(JCAT)
   ZDINI(:) = ZDMAX(:)
   !
   !*    0.2   definition of the catchment area concerned by lateral distribution
   !           ------------------------------------------------------------------
   DO J1=1,NBPIX_SUB(JCAT,JSB)!pixels du sous bassin
   ! 
     IF ( ZRW(J1)>0.0 .AND. ZRW(J1)/=XUNDEF ) THEN
       ZMASK(J1)=1.0
     ELSE
       ZMASK(J1)=0.0
     ENDIF
     !
   ENDDO
   !
   CALL FLOWDOWN_SUB(NBPIX_SUB(JCAT,JSB),ZMASK,&
                    XCONN(JCAT,:,:),NLINE(JCAT,:),JSB,JCAT)
   WHERE (ZMASK == 0.0) ZMASK = XUNDEF
   I_DIM = COUNT( ZMASK(1:NBPIX_SUB(JCAT,JSB))/=XUNDEF )
   !
   !*    1.    Calcul of hydrological similarity and topographic indexes 
   !           ---------------------------------------------------------
   !*    1.1   Calcul of averaged deficit and initialisation of indexes
   !           --------------------------------------------------------
   ZTMP=0.
   !
   DO J1=1,NBPIX_SUB(JCAT,JSB)!pixels du sous bassin
     !
     IF (ZMASK(J1)/=XUNDEF) THEN
       IF ( ZRW(J1)>0.0 .AND. ZRW(J1)/=XUNDEF ) THEN
         ZKAPPA(JSB,J1) = ZRW(J1) 
         INPCON = INPCON + 1
         ZDINI(J1) = ZDMAX(J1) - ZRW(J1)
       ELSE
         ZKAPPA(JSB,J1) = 0.0 
         INPCON = INPCON + 1
         ZDINI(J1) = ZDMAX(J1) 
       ENDIF
       !
       IF ( ZDINI(J1) <0.0 ) THEN
         ZTMP = ZTMP - ZDINI(J1) !we stock here water above saturation to be
                                 !       distributed among the others pixels
         ZDINI(J1) = 0.
       ENDIF
       !
       ZDAV = ZDAV + ZDINI(J1)
       !
     ELSE
       !
       ZKAPPA(JSB,J1) = XUNDEF
       !
     ENDIF
     !
   ENDDO !J1
   !
   IF (ZTMP>0.) THEN
     write(*,*) COUNT(ZDINI(:)<0.),' pixels avec ZDINI negatif. Volume total :', ZTMP
     WHERE ( ZDINI(:)>0. ) ZDINI(:) = ZDINI(:)-ZTMP/(COUNT(ZDINI(:)>0.))
   ENDIF
   !
   IF (INPCON > NBPIX_SUB(JCAT,JSB)/1000) THEN
     ZDAV = ZDAV / INPCON 
     ZDAV = ZDAV / ZM
     !
     !*    1.2   Propagation of indexes
     !           ----------------------
     CALL FLOWDOWN_SUB(NBPIX_SUB(JCAT,JSB),ZKAPPA(JJ,:),&
                       XCONN(JCAT,:,:),NLINE(JCAT,:),JSB,JCAT)
     !
     !*    1.3   Distribution of indexes
     !           ----------------------
     J2=1
     !
     DO WHILE ( .NOT.GFOUND .AND. J2.LE.NBPIX_SUB(JCAT,JSB) )
       !
       IF (ZMASK(J2)/=XUNDEF.AND.&
           NPIX_SUB_TO_CAT(JCAT,JSB,J2)/=NUNDEF.AND.&
           XSLOP(JCAT,NPIX_SUB_TO_CAT(JCAT,JSB,J2))/=0.) THEN
         !
         GFOUND = .TRUE.
         ZKVAL = ZKAPPA(JSB,J2) * EXP(XLAMBDA(JCAT,NPIX_SUB_TO_CAT(JCAT,JSB,J2)))
         IF (ZKVAL>0.0) THEN
           ZKVAL = LOG(ZKVAL)
           ZKVALMAX = ZKVAL
           ZKVALMIN = ZKVAL
           ZKAPPA(JSB,J2) = ZKVAL
         ENDIF
         !
       ELSE
         ZKAPPA(JSB,J2) = XUNDEF
       ENDIF
       !
       J2 = J2 + 1
       !
     ENDDO
     !     
     DO J1 = J2,NBPIX_SUB(JCAT,JSB)
       !
       IF (ZMASK(J1)/=XUNDEF.AND.&
          NPIX_SUB_TO_CAT(JCAT,JSB,J1)/=NUNDEF.AND.&
          XSLOP(JCAT,NPIX_SUB_TO_CAT(JCAT,JSB,J1))/=0.) THEN
        !
         ZKVAL = ZKAPPA(JSB,J1) * EXP(XLAMBDA(JCAT,NPIX_SUB_TO_CAT(JCAT,JSB,J1)))
         IF (ZKVAL>0.0) THEN
           ZKVAL = LOG(ZKVAL)
           !
           IF (ZKVAL.GT.ZKVALMAX) THEN
             ZKVALMAX = ZKVAL
           ELSEIF (ZKVAL.LT.ZKVALMIN) THEN
             ZKVALMIN = ZKVAL
           ENDIF
           !
           ZKAPPA(JSB,J1) = ZKVAL
         ENDIF
         !
       ELSE
         !
         ZKAPPA(JSB,J1) = XUNDEF
         !
       ENDIF
       !
     ENDDO
     !
     !*    1.4   Calcul of saturation index
     !           --------------------------
     !
     I_DIM = COUNT( ZMASK(1:NBPIX_SUB(JCAT,JSB))/=XUNDEF )
     ZKAPPA_PACK(:) = XUNDEF
     ZDMAX_PACK (:) = XUNDEF
     ZKAPPA_PACK(1:I_DIM) = PACK(ZKAPPA(JSB,1:NBPIX_SUB(JCAT,JSB)),ZMASK(1:NBPIX_SUB(JCAT,JSB))/=XUNDEF)
     ZDMAX_PACK (1:I_DIM) = PACK(ZDMAX    (1:NBPIX_SUB(JCAT,JSB)),ZMASK(1:NBPIX_SUB(JCAT,JSB))/=XUNDEF)
     !
     IF (ZKVALMAX/=XUNDEF .AND. ZKVALMIN/=XUNDEF) THEN
       INKAPPA = INT((ZKVALMAX - ZKVALMIN) / XSTEPK)
       !
       DO J1=1,INKAPPA
         !
         ZKVAL = ZKVALMIN + (XSTEPK * (J1-1))
         INAS = 0
         INAD = 0
         ZNDMAXAV = 0.0
         ZNKAV = 0.0
         !
         DO J2=1,I_DIM      
           !      
           IF ( ZKAPPA_PACK(J2).GE.ZKVAL ) THEN
             ! saturated pixel
             INAS = INAS + 1
           ELSEIF  (ZKAPPA_PACK(J2).LE.( ZKVAL-(ZDMAX_PACK(J2)/ZM)) ) THEN
             ! dry pixel
             INAD = INAD + 1
             ZNDMAXAV = ZNDMAXAV + ZDMAX_PACK(J2)
           ELSE
             ZNKAV = ZNKAV + ZKAPPA_PACK(J2)
           ENDIF
           !
         ENDDO
         ! 
         IF (INAD == 0) THEN
           ZNDMAXAV = 0.0
         ELSE
           ZNDMAXAV = ZNDMAXAV /  REAL(INAD)
         ENDIF
         !
         IF ( INPCON == INAS .OR. INPCON == INAD .OR. INPCON == (INAD+INAS)) THEN
           ZNKAV = 0.0
         ELSE
           ZNKAV = ZNKAV / REAL(INPCON - INAD - INAS)
         ENDIF
         !
         IF (INPCON /= 0) THEN
           ZNAS = REAL(INAS) / REAL(INPCON)
           ZNAD = REAL(INAD) / REAL(INPCON)
         ENDIF
         !
         ZFUNC = (1 - ZNAS - ZNAD) * ( ZKVAL - ZNKAV )
         IF (ZM /= 0.) ZFUNC = ZFUNC + (ZNAD * (ZNDMAXAV / ZM))
         !
         ZDIF = ABS( ZFUNC - ZDAV )
         !
         IF ( ZDIF.LT.ZDIFMIN ) THEN
           ZDIFMIN = ZDIF
           ZKAPPAC(JSB) = ZKVAL
           ZAS = ZNAS
           ZAD = ZNAD
           ZDMAXAV = ZNDMAXAV
           ZKAV = ZNKAV
         ENDIF
         !
       ENDDO   !J1=1,INKAPPA     
       !
       !*    2.     Local deficits calculation
       !            --------------------------
       !
       !*    2.1    New averaged deficit on A-Ad-As
       !            -------------------------------
       ZDAV = ZDAV * ZM
       !
       IF ( ZAS<1. .AND. ZAD<1. .AND. (ZAS + ZAD/=1.) ) THEN
         ZDAV2 = (ZDAV - ZDMAXAV * ZAD) / (1 - ZAS - ZAD)  
       ENDIF
       !
       !*    2.2    Local deficits 
       !            --------------
       ZSOMME=0.0
       DO J1=1,NBPIX_SUB(JCAT,JSB)
         !
         IF ( ZMASK(J1)/=XUNDEF ) THEN
           !
           IF ( (ZKAPPA(JSB,J1).GT.(ZKAPPAC(JSB) - ZDMAX(J1)/ZM)) .AND. (ZKAPPA(JSB,J1).LT.ZKAPPAC(JSB)) ) THEN
             !
             ZDEF(JSB,J1) = ZM * (ZKAV - ZKAPPA(JSB,J1)) + ZDAV2
             IF (ZDEF(JSB,J1) < 0.0) ZDEF(JSB,J1) = 0.0
             !
           ELSEIF ( ZKAPPA(JSB,J1).GE.ZKAPPAC(JSB) ) THEN
             !
             ZDEF(JSB,J1) = 0.0
             !
           ELSEIF ( ZKAPPA(JSB,J1).LE.(ZKAPPAC(JSB) - ZDMAX(J1)/ZM) ) THEN
             !
             ZDEF(JSB,J1) = ZDMAX(J1)
             !
           ENDIF
           !
           ! nouveau contenu en eau total (m)
           ZSOMME = ZSOMME + ( XWSTOPT(JSB,J1)*XDTOPT(JSB,J1) - ZDEF(JSB,J1) )
           !
         ELSE
           !
           ZDEF(JSB,J1) = ZDMAX(J1)
           !
         ENDIF
         !
       ENDDO
       !
       ! variation du contenu en eau total
       DO J1=1,NBPIX_SUB(JCAT,JSB)
         IF (ZDEF(JSB,J1)<0.0) THEN
           WRITE(*,*) 'LAMBDA=',ZKAPPA(JSB,J1),'LAMBDAC=',ZKAPPAC(JSB)
         ENDIF
       ENDDO
       !
       GTOPD(JSB)=.TRUE.
     ENDIF! IF (ZKVALMAX/=XUNDEF .AND. ZKVALMIN/=XUNDEF) THEN
     !
   ELSE
    !
    !  'Pas de redistribution laterale'
    write(*,*) 'Pas de redistribution laterale',JCAT,JSB!
    GTOPD(JSB)=.FALSE.
    !
    ZKAPPA(JSB,:) = XUNDEF
    ZDEF(JSB,:) = ZDINI(:)
    ZKAPPAC(JSB) = XUNDEF
    !
  ENDIF
  ! 
  ! Changement de variable => retour sur les grands bv
  DO J1=1,NBPIX_SUB(JCAT,JSB)
    IF (NPIX_SUB_TO_CAT(JCAT,JSB,J1)/=NUNDEF)THEN
      PDEF(JCAT,NPIX_SUB_TO_CAT(JCAT,JSB,J1))=ZDEF(JSB,J1)
      PKAPPA(JCAT,NPIX_SUB_TO_CAT(JCAT,JSB,J1))=ZKAPPA(JSB,J1)
    ENDIF
  ENDDO
  PKAPPAC(JCAT)=ZKAPPAC(JSB)
  OTOPD(JCAT)=GTOPD(JSB)
 ! 
ENDDO!JSB
ENDDO!JCAT
!
IF (LHOOK) CALL DR_HOOK('TOPODYN_LAT_SUB',1,ZHOOK_HANDLE)
!
END SUBROUTINE TOPODYN_LAT_SUB
