!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!-----------------------------------------------------------------
!--------------- special set of characters for SCCS information
!-----------------------------------------------------------------
!
!     #######################
      SUBROUTINE SELECT_MEMBER(U)
!     #######################
!
!!****  *SELECT_MEMBER*  
!!
!!    PURPOSE
!!    -------
!     This routine aims at selecting members in a pool of scenarii
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
!!    
!!    
!!
!!      
!!    REFERENCE
!!    ---------
!!
!!    
!!      
!!    AUTHOR
!!    ------
!!
!!      B. Vincendon    * Meteo-France *
!!
!!    MODIFICATIONS
!!    -------------
!!
!!      Original   11/2007
!-------------------------------------------------------------------------------
!
!*       0.     DECLARATIONS
!               ------------
USE MODD_SURF_ATM_n, ONLY : SURF_ATM_t
USE MODD_SURF_PAR,   ONLY : XUNDEF
USE MODD_PERT_RAIN
USE MODI_OPEN_FILE
USE MODI_CLOSE_FILE
!
USE YOMHOOK   ,ONLY : LHOOK, DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
TYPE(SURF_ATM_t), INTENT(INOUT) :: U
!
!*      0.2    declarations of local variables
!
!
INTEGER               :: JWRK, JWRK1, JWRK2, JSTP
INTEGER               :: ITMP, ITIMES, INB_POP, INB_LOCMOD
INTEGER               :: ISEED_SIZE,JI,IUNIT
REAL                  :: ZRAND_PROB       ! random probability
REAL,    DIMENSION(:,:), ALLOCATABLE :: ZCUMRAIN_NEW   ! hourly rain on large domain
REAL,    DIMENSION(:),   ALLOCATABLE :: ZCUM   ! hourly rain on large domain
INTEGER, DIMENSION(:),   ALLOCATABLE :: INUM_MEMB
INTEGER, DIMENSION(:),   ALLOCATABLE :: ISEED
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('SELECT_MEMBER',0,ZHOOK_HANDLE)
!
ALLOCATE (ZCUMRAIN_NEW(U%NDIM_FULL,NNB_MEMBERS_LOC))
ALLOCATE (ZCUM(NNB_MEMBERS_LOC))
ZCUMRAIN_NEW(:,:)=0.
ZCUM(:)=0.
DO JWRK1=1,NNB_MEMBERS_LOC
  DO JWRK=1,U%NDIM_FULL
    DO JSTP=1,NNB_STEPS_MODIF
      ZCUMRAIN_NEW(JWRK,JWRK1)=ZCUMRAIN_NEW(JWRK,JWRK1)+&
                           XRAIN_NEW(JSTP,JWRK,JWRK1)
    ENDDO
  ENDDO
  ZCUM(JWRK1)=SUM(ZCUMRAIN_NEW(:,JWRK1),MASK=ZCUMRAIN_NEW(:,JWRK1)/=XUNDEF)/&
  COUNT(ZCUMRAIN_NEW(:,JWRK1)/=XUNDEF)
ENDDO  

INB_LOCMOD=50!nombre de membre de localisation changée
ALLOCATE(NMEMBER_SEL(INB_LOCMOD))
NMEMBER_SEL(:)=0
!
XPROB_DECAL(:)=(XPROB_DECAL(:)-MINVAL(XPROB_DECAL))/(MAXVAL(XPROB_DECAL)-MINVAL(XPROB_DECAL))
ALLOCATE (INUM_MEMB(NNB_MEMBERS_LOC*100))
INUM_MEMB(:)=0
JWRK1=1
DO JWRK=1,NNB_MEMBERS_LOC
  ITIMES=FLOOR(XPROB_DECAL(JWRK)*10)
  DO JWRK2=JWRK1,JWRK1+ITIMES
    INUM_MEMB(JWRK2)=JWRK
  ENDDO
  JWRK1=JWRK1+ITIMES+1
ENDDO
INB_POP=JWRK1

CALL OPEN_FILE('ASCII ',IUNIT,HFILE='random.txt',HFORM='FORMATTED',HACTION='READ')
JWRK=1
DO WHILE(JWRK<=INB_LOCMOD)
!
!!!!! Préparation des nombres aléatoires
CALL RANDOM_SEED(ISEED_SIZE)
ALLOCATE(ISEED(ISEED_SIZE))
READ(IUNIT,*) ITMP
ISEED= ITMP+ 37 * (/ (JI - 1, JI = 1, ISEED_SIZE) /)
CALL RANDOM_SEED(PUT = ISEED)
!
!!!!! Préparation des nombres aléatoires
CALL RANDOM_NUMBER(ZRAND_PROB)
write(*,*) 'SELECT MEMBER',JWRK,ISEED_SIZE,ZRAND_PROB
DEALLOCATE(ISEED)
ITMP=FLOOR(INB_POP*ZRAND_PROB)
IF (JWRK==1) THEN
  NMEMBER_SEL(JWRK)=INUM_MEMB(ITMP)
ELSE
  DO JWRK1=1,JWRK-1
    IF (INUM_MEMB(ITMP)==NMEMBER_SEL(JWRK1)) THEN
      GOTO 10
    ELSE
      NMEMBER_SEL(JWRK)=INUM_MEMB(ITMP)
    ENDIF
  ENDDO
ENDIF
JWRK=JWRK+1
10 CONTINUE
ENDDO
CALL CLOSE_FILE('ASCII ',IUNIT)

DO JWRK=1,INB_LOCMOD
  write(*,*) NMEMBER_SEL(JWRK),ZCUM(NMEMBER_SEL(JWRK))*3600.
ENDDO
NNB_MEMBERS_LOC=INB_LOCMOD
DEALLOCATE (INUM_MEMB)
!
IF (LHOOK) CALL DR_HOOK('SELECT_MEMBER',1,ZHOOK_HANDLE)
!
END SUBROUTINE SELECT_MEMBER
