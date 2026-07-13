!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!-------------------------------------------------------------------------------
!     ####################
      SUBROUTINE DESC_CATCHMENTS(KCAT)
!     ####################
!
!!****  *DESC_CATCHMENTS*  
!!
!!    PURPOSE
!!    -------
!!    This routine aims at :
!!    -finding the sub catchments from the river crossings 
!     The variables indexes correspond to the number of pixel in the square that
!     contains the main catchment. 
!     The .map files are also written with this rule.
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
!!    REFERENCE
!!    ---------
!!
!!    
!!     
!!    AUTHOR
!!    ------
!!
!!      B. Vincendon	* Meteo-France *
!!
!!    MODIFICATIONS
!!    -------------
!!
!!      Original    june 2013
!-------------------------------------------------------------------------------
!
!*       0.     DECLARATIONS
!               ------------
!
USE MODD_TOPODYN        ,ONLY : NNMC, XCONN, CCAT, NLINE, NNPT, XTOPD, NNXC, NNYC, &
                                XDXT, XX0, XY0, NNCAT_MAX, XDIST_OUTLET, XDRIV
USE MODD_COUPLING_TOPD  ,ONLY : NSUBCAT, XLX, XLY
USE MODD_SURF_PAR       ,ONLY : XUNDEF, NUNDEF
!
USE MODI_GET_LUOUT
USE MODI_OPEN_FILE
USE MODI_CLOSE_FILE
USE MODI_WRITE_FILE_1MAP
USE MODI_GET_UPSLOPE
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments

INTEGER, INTENT(IN)                :: KCAT   ! nb of the "big" catchment
!
!*      0.2    declarations of local variables
!
!
INTEGER                        :: IUNIT            !unit of discharge files
INTEGER                        :: INBUPSLOPES      !number of contributive pixels
INTEGER                        :: IJRIV_JUNCTION   !index of river crossings
INTEGER                        :: JI, JJ, JK, JSUB ! Loop variables
INTEGER                        :: IX_RIV, IY_RIV, IZ_RIV
REAL                           :: ZLX_RIV, ZLY_RIV
INTEGER                        :: INB_RIVJUNC
INTEGER                        :: ITMP1
CHARACTER(LEN=30)              :: YVAR             ! variable name
REAL                           :: ZDIST_TMP
!
INTEGER, DIMENSION(NNPT(KCAT)) :: IPIXREC_TO_PIXSUB
INTEGER, DIMENSION(NNPT(KCAT)) :: IRIV_JUNCTION
INTEGER, DIMENSION(NNPT(KCAT)) :: IMASK_UP  ! mask of the sub-catchment
INTEGER, DIMENSION(NNMC(KCAT)) :: ILIST_PIX_SUBCAT!list of pix of the sub-catchment
INTEGER, DIMENSION(NNPT(KCAT)) :: INB_SUBPIX ! number of pixels in the sub-catchment
INTEGER, DIMENSION(NNPT(KCAT)) :: IFILE     ! other maner to sort connected points
!
CHARACTER(LEN=28) :: YFILE
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('DESC_CATCHMENTS',0,ZHOOK_HANDLE)
!
!*       1.0.     Initialization :
!               --------------
IMASK_UP(:) = NUNDEF
IPIXREC_TO_PIXSUB(:) = NUNDEF
IRIV_JUNCTION(:) = NUNDEF
!
!*   1.2. Find pixels contributing to each river pixel
!  -----------------------
DO JI=1,NNMC(KCAT)
  JJ=INT(XCONN(KCAT,JI,1))
  IFILE(JJ)=JI
ENDDO
!
WHERE (XTOPD(KCAT,:)==XUNDEF)
  IRIV_JUNCTION(:)=NUNDEF
ELSEWHERE
  IRIV_JUNCTION(:)=0
ENDWHERE
!
!*   1.2. Find pixels contributing to each river pixel
!  -----------------------
!
IJRIV_JUNCTION=INT(XCONN(KCAT,NNMC(KCAT),1)) !le premier point d'intersection est l'exutoire
INB_RIVJUNC = 1
IRIV_JUNCTION(IJRIV_JUNCTION)=IJRIV_JUNCTION
IMASK_UP(:) = 0
!
CALL GET_UPSLOPE(XCONN(KCAT,:,:),IFILE(:),NNMC(KCAT),IJRIV_JUNCTION,XDXT(KCAT),&
                    NNXC(KCAT),NNYC(KCAT),XTOPD(KCAT,:),&
                    ILIST_PIX_SUBCAT(:),INB_SUBPIX(IJRIV_JUNCTION))
 write(*,*) 'DESC_CATCHMENTS=> JUNCTION num ',INB_RIVJUNC,'Surface of sub catchment (in Km²) :  '&
            ,INB_SUBPIX(IJRIV_JUNCTION)*XDXT(KCAT)/1000.*XDXT(KCAT)/1000.
 
DO JK=1,INB_SUBPIX(IJRIV_JUNCTION)
  JJ=ILIST_PIX_SUBCAT(JK)
  IF (JJ/=0.AND.JJ/=NUNDEF) THEN
    IMASK_UP(JJ)=INB_RIVJUNC
    IPIXREC_TO_PIXSUB(JJ)=JK
  ENDIF
ENDDO
!
WHERE (IMASK_UP(:)==INB_RIVJUNC)
  IRIV_JUNCTION(:)=NLINE(KCAT,IJRIV_JUNCTION)
ENDWHERE
!
DO JSUB=1,NSUBCAT(KCAT)
  !
  ZLX_RIV = XLX(KCAT,JSUB)
  ZLY_RIV = XLY(KCAT,JSUB)
  IX_RIV = FLOOR(ZLX_RIV-XX0(KCAT))/XDXT(KCAT)
  IY_RIV = FLOOR(ZLY_RIV-XY0(KCAT))/XDXT(KCAT)
  IZ_RIV = IX_RIV+(IY_RIV-1)*NNXC(KCAT)
  !
  write(*,*) 'DESC_CATCHMENTS: Lx0= ',XX0(KCAT),' Ly0= ',XY0(KCAT)
  write(*,*) 'DESC_CATCHMENTS: num X: ',NNXC(KCAT),'num Y:',NNYC(KCAT)
  write(*,*) 'DESC_CATCHMENTS: Lx= ',ZLX_RIV,' Ly= ',ZLY_RIV,&
           'index i= ',IX_RIV,'index j= ',IY_RIV,'num pix in dom',IZ_RIV
  !
  INBUPSLOPES=INT(XCONN(KCAT,NLINE(KCAT,IZ_RIV),4))
  IJRIV_JUNCTION=IZ_RIV 
  INB_RIVJUNC= INB_RIVJUNC+1
  !
  CALL GET_UPSLOPE(XCONN(KCAT,:,:),IFILE(:),NNMC(KCAT),IJRIV_JUNCTION,XDXT(KCAT),&
                    NNXC(KCAT),NNYC(KCAT),XTOPD(KCAT,:),&
                    ILIST_PIX_SUBCAT(:),INB_SUBPIX(IJRIV_JUNCTION))
  !
  write(*,*) 'DESC_CATCHMENTS=> JUNCTION num ',INB_RIVJUNC,& !'Nb of pix:',INB_SUBPIX(IJRIV_JUNCTION),&
           'Surface of sub catchment (in Km2) :  ',INB_SUBPIX(IJRIV_JUNCTION)*XDXT(KCAT)/1000.*XDXT(KCAT)/1000.
  !
  DO JK=1,INB_SUBPIX(IJRIV_JUNCTION)
    JJ=ILIST_PIX_SUBCAT(JK)
    !
    IF (JJ/=0.AND.JJ/=NUNDEF) THEN
      !
      IF ((IMASK_UP(JJ)/=0).AND.(IMASK_UP(JJ)/=INB_RIVJUNC))THEN
        ITMP1=IMASK_UP(JJ)
        ZDIST_TMP= XDRIV(KCAT,NLINE(KCAT,IJRIV_JUNCTION))-XDRIV(KCAT,IRIV_JUNCTION(JJ))
      ENDIF
      !
    ENDIF
    !
  ENDDO
  !
  IF (ZDIST_TMP>0.0) XDIST_OUTLET(KCAT,ITMP1,INB_RIVJUNC)=ZDIST_TMP
  
  DO JK=1,INB_SUBPIX(IJRIV_JUNCTION)
    JJ=ILIST_PIX_SUBCAT(JK)
    !
    IF (JJ/=0.AND.JJ/=NUNDEF) THEN
      !
      IMASK_UP(JJ)=INB_RIVJUNC
      IPIXREC_TO_PIXSUB(JJ)=JK
      IRIV_JUNCTION(JJ)=NLINE(KCAT,IJRIV_JUNCTION)
      !
    ENDIF!JJ defined
    !
  ENDDO
  !
ENDDO !JSUB
!
write(*,*) 'IRIV_JUNCTION=0',COUNT(IRIV_JUNCTION(:)==0)
write(*,*) 'IRIV_JUNCTION=XUNDEF',COUNT(IRIV_JUNCTION(:)==NUNDEF)
write(*,*) 'NB rensei IRIV_JUNCTION',INB_RIVJUNC

YVAR='SubCat'
CALL WRITE_FILE_1MAP(IMASK_UP(:)*1.,YVAR,KCAT)
YVAR='Outlet'
CALL WRITE_FILE_1MAP(IRIV_JUNCTION(:)*1.,YVAR,KCAT)
YVAR='SubPix'
CALL WRITE_FILE_1MAP(IPIXREC_TO_PIXSUB(:)*1.,YVAR,KCAT)
!
YFILE = TRIM(CCAT(KCAT))//'_overlap'
CALL OPEN_FILE('ASCII ',IUNIT,YFILE,'FORMATTED',HACTION='WRITE')
!
DO JJ=1,NNCAT_MAX
  WRITE(IUNIT,*) XDIST_OUTLET(KCAT,JJ,:)
ENDDO
!
CALL CLOSE_FILE('ASCII ',IUNIT)
!

!
IF (LHOOK) CALL DR_HOOK('DESC_CATCHMENTS',1,ZHOOK_HANDLE)
!
!-----------------------------------------------------------------!
!
END SUBROUTINE DESC_CATCHMENTS
