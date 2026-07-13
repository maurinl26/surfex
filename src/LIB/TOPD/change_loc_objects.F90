!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!     #######################
      SUBROUTINE CHANGE_LOC_OBJECTS(PRR,KI,KJ,KSTEP,PRR_NEW)
!     #######################
!
!!****  *CHANGE_LOC_OBJECTS*  
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
USE MODD_SURF_PAR,     ONLY : XUNDEF, NUNDEF
USE MODD_OBJ,          ONLY : LINAREA, NNUMREG, LRADOK
USE MODD_PERT_RAIN,    ONLY : NNB_MAX_MEMBERS, NNB_MEMBERS_LOC, XPROB_DECAL
!
USE MODI_FIND_POINTS
USE MODI_REGIONS_2D
USE MODI_GAUSS_1REAL
!
USE YOMHOOK   ,ONLY : LHOOK, DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
REAL, DIMENSION(:)   ,INTENT(IN)    :: PRR     ! Rainfall initial field
INTEGER              ,INTENT(IN)    :: KI
INTEGER              ,INTENT(IN)    :: KJ
INTEGER              ,INTENT(IN)    :: KSTEP
REAL, DIMENSION(:,:) ,INTENT(INOUT) :: PRR_NEW     ! Output rainfall fields
!
!*      0.2    declarations of local variables
!
!
INTEGER :: JWRK, JWRK0, JWRK1, JWRK2
INTEGER :: IDX, IDY, IREG2, IMODIF1
INTEGER :: IMODIF, IDELTAX, IDELTAY
INTEGER :: IND_DECAL, INDIM, IPOINT, ITMP
REAL    :: ZTMP, ZPROBX, ZPROBY
! Objets
! Delineate regions
INTEGER                             :: INBTREG ! Total Number of regions
INTEGER                             :: INBREGMAX ! Total Number of regions
INTEGER                             :: JT,ICNT ! Indexes
INTEGER, DIMENSION(:,:),ALLOCATABLE :: NN_REG_ARO
INTEGER, DIMENSION(:),  ALLOCATABLE :: NTOT_REG_ARO
REAL,    DIMENSION(:),    ALLOCATABLE :: ZTHRESHOLD
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('CHANGE_LOC_OBJECTS',0,ZHOOK_HANDLE)
!
INDIM=KI*KJ
!
! Isolating individual objects where RR>2mm/h (O2) and 9mm/h(O9)
ALLOCATE (ZTHRESHOLD(1))
ZTHRESHOLD(1)=2./3600.
!
ALLOCATE (NN_REG_ARO(INDIM,1))
NN_REG_ARO(:,:)=0
!
ALLOCATE (NTOT_REG_ARO(1))
!
INBREGMAX=30
IF (.NOT.ALLOCATED(LINAREA)) ALLOCATE (LINAREA(KI,KJ))
IF (.NOT.ALLOCATED(NNUMREG)) ALLOCATE (NNUMREG(KI,KJ))
!

JT=1!only thresholds (2mm/h)
NNUMREG(:,:)=0
CALL FIND_POINTS(PRR(:),ZTHRESHOLD(JT),KI,KJ,LINAREA)
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
!
DO JWRK2=1,KJ
  DO JWRK1=1,KI
    NN_REG_ARO((JWRK2-1)*KI+JWRK1,JT)=NNUMREG(JWRK1,JWRK2)
    NTOT_REG_ARO(JT)=INBTREG
  ENDDO !jwrk1
ENDDO !jwrk2
!
! Finding barycenters of each O2
IF (.NOT.ALLOCATED(LRADOK)) ALLOCATE(LRADOK(INDIM))
!
! Introduction of perturbation
! Perturbation 
IMODIF=0
DO IDELTAX=-20,20,2
  DO IDELTAY=-20,20,2
    IMODIF=IMODIF+1
!Treatment of region with RR1h>2mm
    DO JWRK=1,INDIM
      IF ((NN_REG_ARO(JWRK,1)>0).AND.(NN_REG_ARO(JWRK,1)/=NUNDEF)) THEN 
      ! in area where rr>2mm
        IND_DECAL=MIN(MAX(0,JWRK+(IDELTAY*KI)+IDELTAX),INDIM)
        PRR_NEW(IND_DECAL,IMODIF)=PRR(JWRK)!decal original field
      ENDIF
    ENDDO !jwrk
    !
    CALL GAUSS_1REAL(-2.0,13.0,IDELTAX*1.,ZPROBX)
    CALL GAUSS_1REAL(-1.0,13.0,IDELTAY*1.,ZPROBY)
    XPROB_DECAL(IMODIF)=ZPROBX*ZPROBY
  ENDDO!ideltay
ENDDO !ideltax 

WHERE (PRR_NEW(:,:) ==XUNDEF) 
  PRR_NEW(:,:) =0.
ENDWHERE
NNB_MEMBERS_LOC=IMODIF

DEALLOCATE(ZTHRESHOLD)
DEALLOCATE(NN_REG_ARO)
DEALLOCATE(NTOT_REG_ARO)
DEALLOCATE(LINAREA)
DEALLOCATE(NNUMREG)
DEALLOCATE(LRADOK)
!
IF (LHOOK) CALL DR_HOOK('CHANGE_LOC_OBJECTS',1,ZHOOK_HANDLE)
!
END SUBROUTINE CHANGE_LOC_OBJECTS
