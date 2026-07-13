!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!     ############################
      SUBROUTINE INIT_TYPES_PARAM(U)
!     ############################
!
!!**** *INIT_TYPES_PARAM* initializes cover-field correspondance arrays
!!
!!    PURPOSE
!!    -------
!!
!!    METHOD
!!    ------
!!
!!
!!    EXTERNAL
!!    --------
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!    REFERENCE
!!    ---------
!!
!!    AUTHOR
!!    ------
!!
!!    S.Faroux        Meteo-France
!!
!!    MODIFICATION
!!    ------------
!!
!!    Original    23/03/11
!!
!!    R. Alkama    05/2012 : read 19 vegtypes rather than 12
!!    E. Martin    10/2014 : add status='old' for ecoclimap.bin files
!!    A. Druel     02/2019 : remove part non compatible with TOWN_TO_ROCK - but temp ! #rustine
!!
!----------------------------------------------------------------------------
!
!*    0.     DECLARATION
!            -----------

USE MODD_TYPE_DATE_SURF
USE MODD_SURF_ATM_n, ONLY : SURF_ATM_t
!
USE MODD_DATA_COVER,     ONLY : XDATA_TOWN, XDATA_NATURE, XDATA_SEA, XDATA_WATER,   &
                                XDATA_VEGTYPE, XDATA_GARDEN, XDATA_Z0_TOWN, &
                                XDATA_BLD, XDATA_BLD_HEIGHT, XDATA_WALL_O_HOR, &
                                XDATA_H_TRAFFIC, XDATA_H_INDUSTRY
!
USE MODD_DATA_COVER_PAR, ONLY : NTYPE, NUT_CPHR, NUT_CPMR, NUT_CPLR, NUT_OPHR,     &
                                NUT_OPMR, NUT_OPLR, NUT_LWLR, NUT_LALR, NUT_SPAR,  &
                                NUT_INDU, NVT_NO, NVT_ROCK, NVT_SNOW, NVT_BOBD,    &
                                NVT_TEBD, NVT_TRBD, NVT_TEBE, NVT_TRBE, NVT_BONE,   &  
                                NVT_TENE, NVT_BOND, NVT_SHRB, NVT_BOGR, NVT_GRAS,   &  
                                NVT_TROG, NVT_C3S, NVT_C3W, NVT_C4, NVT_FLTR,    &
                                NVT_FLGR, NVEGTYPE
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*    0.1    Declaration of arguments
!            ------------------------
!
TYPE(SURF_ATM_t),   INTENT(INOUT) :: U

!*    0.2    Declaration of local variables
!            ------------------------------
!
INTEGER                     :: JTYPE,JL
INTEGER, DIMENSION(NVEGTYPE):: LIST_NVT_GRAS, LIST_NVT_NO, LIST_NVT_TEBD, NVEG_LIST 
INTEGER                     :: CHOICE_GRAS, CHOICE_NO, CHOICE_TEBD  
INTEGER                     :: ILUOUT                
!
!*    0.3    Declaration of namelists
!            ------------------------
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('INIT_TYPES_PARAM',0,ZHOOK_HANDLE)
!
! SEA
!
XDATA_SEA(1:NTYPE(1)) = 1.
!
! WATER
!
XDATA_WATER(NTYPE(1)+1:SUM(NTYPE(1:2))) = 1.
!
! NATURE
!
DO JTYPE = 1,NTYPE(3)
  XDATA_NATURE (SUM(NTYPE(1:2))+JTYPE) = 1.
  XDATA_VEGTYPE(SUM(NTYPE(1:2))+JTYPE,JTYPE) = 1.
ENDDO
!
! TOWN
!
XDATA_TOWN(SUM(NTYPE(1:3))+1:SUM(NTYPE(1:4))) = 1.
!
XDATA_GARDEN(NUT_CPHR) = 0.10     
XDATA_GARDEN(NUT_CPMR) = 0.10
XDATA_GARDEN(NUT_CPLR) = 0.10
XDATA_GARDEN(NUT_OPHR) = 0.30     
XDATA_GARDEN(NUT_OPMR) = 0.30
XDATA_GARDEN(NUT_OPLR) = 0.40
XDATA_GARDEN(NUT_LWLR) = 0.10
XDATA_GARDEN(NUT_LALR) = 0.10
XDATA_GARDEN(NUT_SPAR) = 0.70
XDATA_GARDEN(NUT_INDU) = 0.40
!
XDATA_Z0_TOWN(NUT_CPHR) = 3.
XDATA_Z0_TOWN(NUT_CPMR) = 1.
XDATA_Z0_TOWN(NUT_CPLR) = 0.5
XDATA_Z0_TOWN(NUT_OPHR) = 2.0
XDATA_Z0_TOWN(NUT_OPMR) = 0.5
XDATA_Z0_TOWN(NUT_OPLR) = 0.5
XDATA_Z0_TOWN(NUT_LWLR) = 0.25
XDATA_Z0_TOWN(NUT_LALR) = 0.25
XDATA_Z0_TOWN(NUT_SPAR) = 0.5
XDATA_Z0_TOWN(NUT_INDU) = 0.5
!
XDATA_BLD(NUT_CPHR) = 0.5
XDATA_BLD(NUT_CPMR) = 0.55
XDATA_BLD(NUT_CPLR) = 0.55
XDATA_BLD(NUT_OPHR) = 0.30
XDATA_BLD(NUT_OPMR) = 0.30
XDATA_BLD(NUT_OPLR) = 0.30
XDATA_BLD(NUT_LWLR) = 0.75
XDATA_BLD(NUT_LALR) = 0.40
XDATA_BLD(NUT_SPAR) = 0.10
XDATA_BLD(NUT_INDU) = 0.25
!
XDATA_BLD_HEIGHT(NUT_CPHR) = 75
XDATA_BLD_HEIGHT(NUT_CPMR) = 20
XDATA_BLD_HEIGHT(NUT_CPLR) = 5
XDATA_BLD_HEIGHT(NUT_OPHR) = 75  
XDATA_BLD_HEIGHT(NUT_OPMR) = 20
XDATA_BLD_HEIGHT(NUT_OPLR) = 5
XDATA_BLD_HEIGHT(NUT_LWLR) = 3
XDATA_BLD_HEIGHT(NUT_LALR) = 5
XDATA_BLD_HEIGHT(NUT_SPAR) = 5
XDATA_BLD_HEIGHT(NUT_INDU) = 10
!
XDATA_WALL_O_HOR(NUT_CPHR) = 4.0
XDATA_WALL_O_HOR(NUT_CPMR) = 1.3
XDATA_WALL_O_HOR(NUT_CPLR) = 0.9
XDATA_WALL_O_HOR(NUT_OPHR) = 1.4  
XDATA_WALL_O_HOR(NUT_OPMR) = 0.7
XDATA_WALL_O_HOR(NUT_OPLR) = 0.7
XDATA_WALL_O_HOR(NUT_LWLR) = 0.75
XDATA_WALL_O_HOR(NUT_LALR) = 0.24
XDATA_WALL_O_HOR(NUT_SPAR) = 0.36
XDATA_WALL_O_HOR(NUT_INDU) = 0.45
!
XDATA_H_TRAFFIC(NUT_CPHR) = 20
XDATA_H_TRAFFIC(NUT_CPMR) = 10
XDATA_H_TRAFFIC(NUT_CPLR) = 5
XDATA_H_TRAFFIC(NUT_OPHR) =  20 
XDATA_H_TRAFFIC(NUT_OPMR) = 10
XDATA_H_TRAFFIC(NUT_OPLR) = 5
XDATA_H_TRAFFIC(NUT_LWLR) = 5
XDATA_H_TRAFFIC(NUT_LALR) = 5
XDATA_H_TRAFFIC(NUT_SPAR) = 5
XDATA_H_TRAFFIC(NUT_INDU) = 5
!
XDATA_H_INDUSTRY(NUT_CPHR) = 0.
XDATA_H_INDUSTRY(NUT_CPMR) = 0.
XDATA_H_INDUSTRY(NUT_CPLR) = 0.
XDATA_H_INDUSTRY(NUT_OPHR) = 0.
XDATA_H_INDUSTRY(NUT_OPMR) = 0.
XDATA_H_INDUSTRY(NUT_OPLR) = 0.
XDATA_H_INDUSTRY(NUT_LWLR) = 0.
XDATA_H_INDUSTRY(NUT_LALR) = 50.
XDATA_H_INDUSTRY(NUT_SPAR) = 0.
XDATA_H_INDUSTRY(NUT_INDU) = 100.
!
!#rustine Temporary modification to solve ECOSG and urban vegetation problem
!#rustine If NO, GRASS or TEBD absent of the domain, PGD can crash. So, we look for
!#rustine a type of vegetation is present.
IF (ASSOCIATED(U%LCOVER)) THEN
   IF (ASSOCIATED(U%LCOVER).AND.(SIZE(U%LCOVER,1).GT.1).AND.(.NOT. U%LTOWN_TO_ROCK )) THEN
      LIST_NVT_GRAS = (/NVT_GRAS, NVT_BOGR  , NVT_TROG , NVT_C3W , NVT_C3S, NVT_C4, NVT_FLGR, &
                        NVT_SHRB, NVT_FLTR, NVT_ROCK, NVT_NO  , NVT_SNOW, NVT_TEBD, NVT_TENE, &
                        NVT_TEBE, NVT_BOBD, NVT_BONE, NVT_BOND, NVT_TRBD, NVT_TRBE/)
      LIST_NVT_NO   = (/NVT_NO  , NVT_ROCK, NVT_SNOW, NVT_GRAS, NVT_C4  , NVT_C3S , NVT_C3W , &
                        NVT_BOGR, NVT_TROG, NVT_FLGR, NVT_SHRB, NVT_FLTR, NVT_TEBD, NVT_TENE, &
                        NVT_TEBE, NVT_BOBD, NVT_BONE, NVT_BOND, NVT_TRBD, NVT_TRBE/)
      LIST_NVT_TEBD = (/NVT_TEBD, NVT_TENE, NVT_TEBE, NVT_BOBD, NVT_BONE, NVT_BOND, NVT_TRBD, &
                        NVT_TRBE, NVT_SHRB, NVT_FLTR, NVT_FLGR, NVT_GRAS, NVT_C4  , NVT_C3S , &
                        NVT_C3W , NVT_BOGR, NVT_TROG, NVT_NO  , NVT_ROCK, NVT_SNOW/)
      !
      CHOICE_GRAS = 0
      CHOICE_NO   = 0
      CHOICE_TEBD = 0
      DO JL = 1, size(LIST_NVT_GRAS,1)
         IF ( CHOICE_GRAS == 0 .AND. U%LCOVER(SUM(NTYPE(1:2)) + LIST_NVT_GRAS(JL)) ) CHOICE_GRAS = LIST_NVT_GRAS(JL)
         IF ( CHOICE_NO   == 0 .AND. U%LCOVER(SUM(NTYPE(1:2)) + LIST_NVT_NO  (JL)) ) CHOICE_NO   = LIST_NVT_NO  (JL)
         IF ( CHOICE_TEBD == 0 .AND. U%LCOVER(SUM(NTYPE(1:2)) + LIST_NVT_TEBD(JL)) ) CHOICE_TEBD = LIST_NVT_TEBD(JL)
      END DO
      !
      IF ( CHOICE_GRAS == 0 ) THEN
         WRITE(ILUOUT,*) ' '
         WRITE(ILUOUT,*) '**************************************************************'
         WRITE(ILUOUT,*) '**************************************************************'
         WRITE(ILUOUT,*) '*                          BIG CARE                          *'
         WRITE(ILUOUT,*) '**************************************************************'
         WRITE(ILUOUT,*) '* CARE ! YOU ARE IN THE SITUATION WHERE THERE IS NO NATURE ! *'
         WRITE(ILUOUT,*) '* You need no have (with ECOSG) at least one point with one  *'
         WRITE(ILUOUT,*) '* fraction. Else this #rustine do not take into account the  *'
         WRITE(ILUOUT,*) '* fraction of town "forgeted".                               *'
         WRITE(ILUOUT,*) '* An alternative is to run the model without input map file  *'
         WRITE(ILUOUT,*) '* (as LAI, Albedo, ... but constants ! Good luck ;)          *'
         WRITE(ILUOUT,*) '**************************************************************'
         WRITE(ILUOUT,*) '**************************************************************'
         WRITE(ILUOUT,*) ' '
         IF ( CHOICE_NO/=0 .OR. CHOICE_TEBD/=0 ) CALL ABOR1_SFX('PGD_ISBA_PAR: ERROR WITH CONCEPT OF TOWN FRAC TO VEG FRAC')
      !
      ! Il y a dans la région étudiée de la nature, alors:
      ! On transmet les TOWN_FRAC_X aux VEG_FRAC_Y associés selon les régles initialement définies dans 
      ELSE
         XDATA_VEGTYPE(NUT_CPHR,CHOICE_GRAS) = 0.0
         XDATA_VEGTYPE(NUT_CPMR,CHOICE_GRAS) = 0.0
         XDATA_VEGTYPE(NUT_CPLR,CHOICE_GRAS) = 0.0 
         XDATA_VEGTYPE(NUT_OPHR,CHOICE_GRAS) = 0.4  
         XDATA_VEGTYPE(NUT_OPMR,CHOICE_GRAS) = 0.4
         XDATA_VEGTYPE(NUT_OPLR,CHOICE_GRAS) = 0.4
         XDATA_VEGTYPE(NUT_LWLR,CHOICE_GRAS) = 0.5
         XDATA_VEGTYPE(NUT_LALR,CHOICE_GRAS) = 0.5
         XDATA_VEGTYPE(NUT_SPAR,CHOICE_GRAS) = 0.5
         XDATA_VEGTYPE(NUT_INDU,CHOICE_GRAS) = 0.4
         !
         XDATA_VEGTYPE(NUT_CPHR,CHOICE_NO) = 0.0
         XDATA_VEGTYPE(NUT_CPMR,CHOICE_NO) = 0.0
         XDATA_VEGTYPE(NUT_CPLR,CHOICE_NO) = 0.0 
         XDATA_VEGTYPE(NUT_OPHR,CHOICE_NO) = 0.2  
         XDATA_VEGTYPE(NUT_OPMR,CHOICE_NO) = 0.2
         XDATA_VEGTYPE(NUT_OPLR,CHOICE_NO) = 0.2
         XDATA_VEGTYPE(NUT_LWLR,CHOICE_NO) = 0.5
         XDATA_VEGTYPE(NUT_LALR,CHOICE_NO) = 0.5
         XDATA_VEGTYPE(NUT_SPAR,CHOICE_NO) = 0.2
         XDATA_VEGTYPE(NUT_INDU,CHOICE_NO) = 0.6
         !
         XDATA_VEGTYPE(NUT_CPHR,CHOICE_TEBD) = 1.0 
         XDATA_VEGTYPE(NUT_CPMR,CHOICE_TEBD) = 1.0
         XDATA_VEGTYPE(NUT_CPLR,CHOICE_TEBD) = 1.0
         XDATA_VEGTYPE(NUT_OPHR,CHOICE_TEBD) = 0.4  
         XDATA_VEGTYPE(NUT_OPMR,CHOICE_TEBD) = 0.4
         XDATA_VEGTYPE(NUT_OPLR,CHOICE_TEBD) = 0.4
         XDATA_VEGTYPE(NUT_LWLR,CHOICE_TEBD) = 0.0
         XDATA_VEGTYPE(NUT_LALR,CHOICE_TEBD) = 0.0
         XDATA_VEGTYPE(NUT_SPAR,CHOICE_TEBD) = 0.3
         XDATA_VEGTYPE(NUT_INDU,CHOICE_TEBD) = 0.0
         !
      ENDIF
   ENDIF
!
ELSE
   XDATA_VEGTYPE(NUT_CPHR,NVT_GRAS) = 0.0
   XDATA_VEGTYPE(NUT_CPMR,NVT_GRAS) = 0.0
   XDATA_VEGTYPE(NUT_CPLR,NVT_GRAS) = 0.0 
   XDATA_VEGTYPE(NUT_OPHR,NVT_GRAS) = 0.4  
   XDATA_VEGTYPE(NUT_OPMR,NVT_GRAS) = 0.4
   XDATA_VEGTYPE(NUT_OPLR,NVT_GRAS) = 0.4
   XDATA_VEGTYPE(NUT_LWLR,NVT_GRAS) = 0.5
   XDATA_VEGTYPE(NUT_LALR,NVT_GRAS) = 0.5
   XDATA_VEGTYPE(NUT_SPAR,NVT_GRAS) = 0.5
   XDATA_VEGTYPE(NUT_INDU,NVT_GRAS) = 0.4
   !
   XDATA_VEGTYPE(NUT_CPHR,NVT_NO) = 0.0
   XDATA_VEGTYPE(NUT_CPMR,NVT_NO) = 0.0
   XDATA_VEGTYPE(NUT_CPLR,NVT_NO) = 0.0 
   XDATA_VEGTYPE(NUT_OPHR,NVT_NO) = 0.2  
   XDATA_VEGTYPE(NUT_OPMR,NVT_NO) = 0.2
   XDATA_VEGTYPE(NUT_OPLR,NVT_NO) = 0.2
   XDATA_VEGTYPE(NUT_LWLR,NVT_NO) = 0.5
   XDATA_VEGTYPE(NUT_LALR,NVT_NO) = 0.5
   XDATA_VEGTYPE(NUT_SPAR,NVT_NO) = 0.2
   XDATA_VEGTYPE(NUT_INDU,NVT_NO) = 0.6
   !
   XDATA_VEGTYPE(NUT_CPHR,NVT_TEBD) = 1.0 
   XDATA_VEGTYPE(NUT_CPMR,NVT_TEBD) = 1.0
   XDATA_VEGTYPE(NUT_CPLR,NVT_TEBD) = 1.0
   XDATA_VEGTYPE(NUT_OPHR,NVT_TEBD) = 0.4  
   XDATA_VEGTYPE(NUT_OPMR,NVT_TEBD) = 0.4
   XDATA_VEGTYPE(NUT_OPLR,NVT_TEBD) = 0.4
   XDATA_VEGTYPE(NUT_LWLR,NVT_TEBD) = 0.0
   XDATA_VEGTYPE(NUT_LALR,NVT_TEBD) = 0.0
   XDATA_VEGTYPE(NUT_SPAR,NVT_TEBD) = 0.3
   XDATA_VEGTYPE(NUT_INDU,NVT_TEBD) = 0.0
!
END IF
IF (LHOOK) CALL DR_HOOK('INIT_TYPES_PARAM',1,ZHOOK_HANDLE)
!
END SUBROUTINE INIT_TYPES_PARAM
