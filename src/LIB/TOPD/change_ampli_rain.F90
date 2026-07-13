!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!     #######################
      SUBROUTINE CHANGE_AMPLI_RAIN(PRR, KI, KJ, PRR_NEW)
!     #######################
!
!!****  *CHANGE_AMPLI_RAIN*  
!!
!!    PURPOSE
!!    -------
!     This routine aims at computing different modified rainfall fields
!     from an initial one.
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
USE MODD_SURF_PAR,       ONLY : XUNDEF, NUNDEF
USE MODD_OBJ,            ONLY : LINAREA, NNUMREG
USE MODD_PERT_RAIN,      ONLY : NNB_STEPS_MODIF
!
USE MODI_FIND_POINTS
USE MODI_REGIONS_2D
!
USE YOMHOOK   ,ONLY : LHOOK, DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
REAL, DIMENSION(:,:) ,INTENT(IN)    :: PRR     ! Rainfall initial field
INTEGER              ,INTENT(IN)    :: KI
INTEGER              ,INTENT(IN)    :: KJ
REAL, DIMENSION(:,:) ,INTENT(OUT)   :: PRR_NEW     ! Output rainfall fields
!
!*      0.2    declarations of local variables
!
!
INTEGER :: JWRK, JWRK0, JWRK1, JWRK2, JSTP
INTEGER :: IMAX_MEMBERS, IDX, IDY, IREG2, IMODIF1
INTEGER :: IMODIF, IDELTAX, IDELTAY, IND_DECAL, IPOINT
INTEGER :: INDIM
! Objets
! Delineate regions
INTEGER                             :: INBTREG ! Total Number of regions
INTEGER                             :: INBREGMAX ! Total Number of regions
INTEGER                             :: JT, ICNT ! Indexes
INTEGER                             :: ITMP, ITMP2, ITMP9 ! Indexes
INTEGER, DIMENSION(:,:),ALLOCATABLE :: NN_REG_ARO
INTEGER, DIMENSION(:),  ALLOCATABLE :: NTOT_REG_ARO
REAL,    DIMENSION(:),  ALLOCATABLE :: ZDIST_BARY_ARO, ZPENTE_BARY_ARO
REAL,    DIMENSION(:),  ALLOCATABLE :: ZTHRESHOLD
INTEGER, DIMENSION(:,:),ALLOCATABLE :: NX_BARY_ARO, NY_BARY_ARO
REAL,    DIMENSION(:,:),ALLOCATABLE :: ZQ !quantiles
REAL                                :: ZQ_O2, ZQ_O9, ZRAND_PROB
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('CHANEG_AMPLI_RAIN',0,ZHOOK_HANDLE)
!
INDIM=KI*KJ
ALLOCATE (ZTHRESHOLD(2))
ZTHRESHOLD(1) = 2./3600.
ZTHRESHOLD(2) = 9./3600.
!
ALLOCATE(NN_REG_ARO(INDIM,2))
NN_REG_ARO(:,:) = 0
!
ALLOCATE (NTOT_REG_ARO(2))
!
INBREGMAX=30
IF (.NOT.ALLOCATED(LINAREA)) ALLOCATE (LINAREA(KI,KJ))
IF (.NOT.ALLOCATED(NNUMREG)) ALLOCATE (NNUMREG(KI,KJ))
 
ALLOCATE (ZQ(2,100))!
! 1=q5
! 2=q10
! 3=q20
! 4=q25
! 5=q30
! 6=q40
! 7=q50=med
! 8=q60
! 9=q70
!10=q75
!11=q80
!12=q90
!13=q95
ZQ(1,1:5)=0.61
ZQ(1,6:10)=0.68
ZQ(1,11:20)=0.79
ZQ(1,21:25)=0.83
ZQ(1,26:30)=0.88
ZQ(1,31:40)=0.97
ZQ(1,41:50)=1.02
ZQ(1,51:60)=1.10
ZQ(1,61:70)=1.19
ZQ(1,71:75)=1.24
ZQ(1,76:85)=1.31
ZQ(1,86:95)=1.49
ZQ(1,96:100)=1.70
!
ZQ(2,1:5)=0.64
ZQ(2,6:10)=0.71
ZQ(2,11:20)=0.81
ZQ(2,21:25)=0.85
ZQ(2,26:30)=0.89
ZQ(2,31:40)=0.97
ZQ(2,41:50)=1.06
ZQ(2,51:60)=1.12
ZQ(2,61:70)=1.18
ZQ(2,71:75)=1.23
ZQ(2,76:85)=1.27
ZQ(2,86:95)=1.48
ZQ(2,95:100)=1.60
!
CALL RANDOM_NUMBER(ZRAND_PROB)
ITMP=FLOOR(100*ZRAND_PROB)
ZQ_O2=ZQ(1,ITMP)
IF (ZQ_O2<0.01) ZQ_O2=1.0
! 
CALL RANDOM_NUMBER(ZRAND_PROB)
ITMP=FLOOR(100*ZRAND_PROB)
ZQ_O9=ZQ(2,ITMP)
IF (ZQ_O9<0.01) ZQ_O9=1.0
! 
write(*,*) 'coefs,',ZQ_O2,ZQ_O9
!
! Modification of Objects amplitude
!
DO JSTP=1,NNB_STEPS_MODIF
  DO JT=1,2!number of different thresholds (2 et 9)
    NNUMREG(:,:)=0
    CALL FIND_POINTS(PRR(:,JSTP),ZTHRESHOLD(JT),KI,KJ,LINAREA)
    CALL REGIONS_2D(LINAREA,NNUMREG,INBTREG)
    !
    ICNT=0
    DO JWRK=1,INBTREG
      IF (COUNT(NNUMREG==JWRK)<20.) THEN
        WHERE(NNUMREG==JWRK) NNUMREG=0
      ELSE
        ICNT=ICNT+1
        WHERE(NNUMREG==JWRK) NNUMREG=ICNT
      ENDIF   
    ENDDO
    INBTREG=ICNT
    DO JWRK2=1,KJ
      DO JWRK1=1,KI
        NN_REG_ARO((JWRK2-1)*KI+JWRK1,JT)=NNUMREG(JWRK1,JWRK2)
        NTOT_REG_ARO(JT)=INBTREG
      ENDDO !jwrk1
    ENDDO !jwrk2
  ENDDO ! jt
  !
  DO JWRK=1,INDIM
    PRR_NEW(JWRK,JSTP)=PRR(JWRK,JSTP)
    IF ((NN_REG_ARO(JWRK,1)>0).AND.(NN_REG_ARO(JWRK,1)/=NUNDEF)) THEN 
      ! in area where rr>2mm
       PRR_NEW(JWRK,JSTP)=PRR(JWRK,JSTP)*ZQ_O2
    ENDIF
    IF ((NN_REG_ARO(JWRK,2)>0).AND.(NN_REG_ARO(JWRK,2)/=NUNDEF)) THEN 
      ! in area where rr>9mm
      PRR_NEW(JWRK,JSTP)=MAX(PRR(JWRK,JSTP)*ZQ_O9,ZTHRESHOLD(2))
    ENDIF
  ENDDO !JWRK
ENDDO !JSTP
PRR_NEW(:,:)=MIN(PRR_NEW(:,:),110./3600.)

DEALLOCATE(ZQ)
DEALLOCATE(LINAREA)
DEALLOCATE(NNUMREG)
DEALLOCATE(ZTHRESHOLD)
DEALLOCATE(NN_REG_ARO)
DEALLOCATE(NTOT_REG_ARO)
!
IF (LHOOK) CALL DR_HOOK('CHANGE_AMPLI_RAIN',1,ZHOOK_HANDLE)
!
END SUBROUTINE CHANGE_AMPLI_RAIN
