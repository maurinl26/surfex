!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!-------------------------------------------------------------------------------
!     ####################
      SUBROUTINE FIND_SUB_CATCH(KCAT)
!     ####################
!
!!****  *FIND_SUB_CATCH*  
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
USE MODD_SURF_PAR,  ONLY : XUNDEF, NUNDEF
!
USE MODD_TOPODYN,   ONLY : NNMC, XCONN, CCAT, NLINE, NNPT, XTOPD, NNXC, NNYC, &
                           XDXT, XX0, XY0, NNCAT_MAX, XDIST_OUTLET, XDHIL, XDRIV
!
USE MODI_GET_LUOUT
USE MODI_OPEN_FILE
USE MODI_CLOSE_FILE
USE MODI_WRITE_FILE_1MAP
USE MODI_GET_UPSLOPE
!
USE YOMHOOK   ,ONLY : LHOOK, DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments

INTEGER, INTENT(IN)              :: KCAT   ! nb of the "big" catchment
!
!
!*      0.2    declarations of local variables
!
!
INTEGER                            :: IUNIT            ! unit of discharge files
INTEGER                            :: INBUPSLOPES      ! number of contributive pixels
INTEGER                            :: IJRIV_JUNCTION   ! index of river crossings
INTEGER                            :: JI,JJ,JK,JL,JSB  ! Loop variables
INTEGER                            :: IX,IY,IZ
INTEGER                            :: IX0,IY0,IZ0
INTEGER                            :: INB_RIVJUNC
INTEGER                            :: ITMP1,ITMP2
INTEGER                            :: JCE              ! Loop variables
CHARACTER(LEN=30)                  :: YVAR             ! variable name
REAL                               :: ZDIST_TMP
!
INTEGER, DIMENSION(NNPT(KCAT))     :: IPIXREC_TO_PIXSUB
INTEGER, DIMENSION(NNPT(KCAT))     :: IRIV_JUNCTION
INTEGER, DIMENSION(NNPT(KCAT))     :: IMASK_UP         ! mask of the sub-catchment
INTEGER, DIMENSION(NNMC(KCAT))     :: ILIST_PIX_SUBCAT ! list of pix of the sub-catchment
INTEGER, DIMENSION(NNPT(KCAT))     :: INB_SUBPIX       ! number of pixels in the sub-catchment
INTEGER, DIMENSION(NNPT(KCAT))     :: IFILE            ! other maner to sort connected points
INTEGER                            :: INB_PIX_TMP
!
CHARACTER(LEN=28) :: YFILE
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('FIND_SUB_CATCH',0,ZHOOK_HANDLE)
!
!*       1.0.     Initialization :
!               --------------
IMASK_UP(:)=NUNDEF
IPIXREC_TO_PIXSUB(:)=NUNDEF
IRIV_JUNCTION(:)=NUNDEF
!
YFILE = TRIM(CCAT(KCAT))//'_outletsXY'
!
CALL OPEN_FILE('ASCII ',IUNIT,YFILE,'FORMATTED',HACTION='WRITE')
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
INB_RIVJUNC=1
IRIV_JUNCTION(IJRIV_JUNCTION)=IJRIV_JUNCTION
IMASK_UP(:)= 0
!
CALL GET_UPSLOPE(XCONN(KCAT,:,:),IFILE(:),NNMC(KCAT),IJRIV_JUNCTION,XDXT(KCAT),&
                    NNXC(KCAT),NNYC(KCAT),XTOPD(KCAT,:),&
                    ILIST_PIX_SUBCAT(:),INB_SUBPIX(IJRIV_JUNCTION))
write(*,*) 'FIND_SUB_CATCH=> JUNCTION num',INB_RIVJUNC,'surface of sub catchment: ',INB_SUBPIX(IJRIV_JUNCTION) 
write(*,*) 'FIND_SUB_CATCH=> IJRIV_JUNCTION ',IJRIV_JUNCTION,SIZE(INB_SUBPIX,1)
!
DO JK=1,INB_SUBPIX(IJRIV_JUNCTION)
  JJ=ILIST_PIX_SUBCAT(JK)
  !
  IF (JJ/=0.AND.JJ/=NUNDEF) THEN
    IMASK_UP(JJ)=INB_RIVJUNC
    IPIXREC_TO_PIXSUB(JJ)=JK
  ENDIF
  !
ENDDO
!
WHERE (IMASK_UP(:)==INB_RIVJUNC)
  IRIV_JUNCTION(:)=NLINE(KCAT,IJRIV_JUNCTION)
ENDWHERE
!!!!
!!! ATTENTION DIFFERENCE AVEC DESC_CATCHMENTS
!!!
DO JI=NNMC(KCAT),1,-1
!!On balaie tous les points du bassin en commencant par l'exutoire
  JCE=0  !on met à 0 les points rivière contributifs pour les compter
  !
  INBUPSLOPES=NINT(XCONN(KCAT,JI,4))
  IZ0=NINT(XCONN(KCAT,JI,1))
  IY0=INT(IZ0*1./NNXC(KCAT))+1
  IX0=IZ0-(IY0-1)*NNXC(KCAT)
  !
  IF(XDHIL(KCAT,JI)==MINVAL(XDHIL(KCAT,:))) THEN 
  !
    DO IX=IX0-1,IX0+1
      DO IY=IY0-1,IY0+1
        IF(IX>0.AND.IY>0.AND.IX<=NNXC(KCAT).AND.IY<=NNYC(KCAT)) THEN
          IZ=IX+(IY-1)*NNXC(KCAT)
          IF (IZ/=IZ0.AND.XTOPD(KCAT,IZ)/=XUNDEF.AND.XTOPD(KCAT,IZ)>=XTOPD(KCAT,IZ0)&
             .AND.XDHIL(KCAT,NLINE(KCAT,IZ))==MINVAL(XDHIL(KCAT,:)))&
             JCE=JCE+1
        ENDIF
      ENDDO
    ENDDO
    !
    IF (JCE>1)THEN
      IJRIV_JUNCTION=IZ0 
      INB_RIVJUNC= INB_RIVJUNC+1
      write(*,*) 'FIND_SUB_CATCH=> passage',JI,'point',IZ0,'x=',IX0,'y=',IY0,'junction',IJRIV_JUNCTION,'JCE>3'
      !
      ILIST_PIX_SUBCAT(:)=NUNDEF
      CALL GET_UPSLOPE(XCONN(KCAT,:,:),IFILE(:),NNMC(KCAT),IJRIV_JUNCTION,XDXT(KCAT),&
                    NNXC(KCAT),NNYC(KCAT),XTOPD(KCAT,:),ILIST_PIX_SUBCAT(:),INB_SUBPIX(IJRIV_JUNCTION))
 
      write(*,*) 'FIND_SUB_CATCH=> JUNCTION num',INB_RIVJUNC,'surface of sub catchment: ',INB_SUBPIX(IJRIV_JUNCTION)
 !
 ! Recherche des recoupements entre le sous bassin qui vient d'être défini et
 ! les sous bassins définis plus tot
 ! On ne sélectionne pas le nouveau sous bassin, si
 ! - il est à cheval sur plusieurs sous bassins précédents
 ! - le nouvel exutoire est trop proche de l'ancien 
 ! - la distance est négative entre les exutoires du nouveau et de l'ancien
 ! Premier balayage pour chercher le premier recoupement et calculer les
 ! distances entre exutoires
 !
      DO JK=1,INB_SUBPIX(IJRIV_JUNCTION)
        JJ=ILIST_PIX_SUBCAT(JK)
        !
        IF (JJ/=0.AND.JJ/=NUNDEF) THEN
          IF ((IMASK_UP(JJ)/=0).AND.(IMASK_UP(JJ)/=INB_RIVJUNC)) THEN
            ITMP1=IMASK_UP(JJ)
            ZDIST_TMP= XDRIV(KCAT,NLINE(KCAT,IJRIV_JUNCTION))-XDRIV(KCAT,IRIV_JUNCTION(JJ))
          ENDIF
        ENDIF
        !
      ENDDO
      !
      IF (1==2) THEN
        20 CONTINUE
        ! Second balayage qui démarre au point où le balayage précedent s'est terminé
        !
        DO JL=JK,INB_SUBPIX(IJRIV_JUNCTION)
          JJ=ILIST_PIX_SUBCAT(JL)
          IF (JJ/=0.AND.JJ/=NUNDEF) THEN
            !
            IF ((IMASK_UP(JJ)/=0).AND.(IMASK_UP(JJ)/=INB_RIVJUNC).AND.&
              XDIST_OUTLET(KCAT,IMASK_UP(JJ),INB_RIVJUNC)==XUNDEF) THEN
              ITMP2=IMASK_UP(JJ)
              IF (ITMP1/=ITMP2)THEN
                write(*,*) 'JJ,ITMP2',JJ,ITMP2
                INB_RIVJUNC= INB_RIVJUNC-1
                GOTO 50 !si sbv sur 2 bassins precedents, on ne le prend pas
              ENDIF
            ENDIF
            !
          ENDIF
        ENDDO
        !
      ENDIF  !(1==2)
      !
      ! Si le second balayage s'est passé sans trouver un autre sous bassin
      ! anterieur, on calcule la superficie du bassin "difference"
      !
      INB_PIX_TMP=COUNT(IMASK_UP(:)==ITMP1)-INB_SUBPIX(IJRIV_JUNCTION)
      IF (INB_PIX_TMP<500) THEN 
        write(*,*) 'INB_PIX_TMP<500'
        INB_RIVJUNC= INB_RIVJUNC-1
        GOTO 50 !si sbv 'difference trop petit, on ne prend pas le nouveau sbv
      ENDIF
      !
      ! Si le second balayage s'est passé sans trouver un autre sous bassin
      ! anterieur, on teste la distance entre exutoires
      IF (ZDIST_TMP>0.0)THEN 
        IF (ZDIST_TMP<1000.)THEN 
          ! si moins de 10 pixels de distance, on ne selectionne pas ce nouveau ss bassin
          write(*,*) 'ZDIST_TMP<1000.'
          INB_RIVJUNC= INB_RIVJUNC-1
          GOTO 50
        ELSE
          XDIST_OUTLET(KCAT,ITMP1,INB_RIVJUNC)=ZDIST_TMP
        ENDIF!0.<XDIST<500.   
      ELSE ! Si la distance est négative,  on ne selectionne pas ce nouveau ss bassin
        write(*,*) 'ZDIST_TMP<0.'
        INB_RIVJUNC= INB_RIVJUNC-1
        GOTO 50
      ENDIF
      !
      ! Si le sous bassin a été sélectionné, NB_RIVJUNCTION a été mis à jour, on
      ! rempli les masques
      DO JK=1,INB_SUBPIX(IJRIV_JUNCTION)
        JJ=ILIST_PIX_SUBCAT(JK)
        !
        IF (JJ/=0.AND.JJ/=NUNDEF) THEN
          IMASK_UP(JJ)=INB_RIVJUNC
          IPIXREC_TO_PIXSUB(JJ)=JK
          IRIV_JUNCTION(JJ)=NLINE(KCAT,IJRIV_JUNCTION)
        ENDIF!JJ defined
        !
      ENDDO
      write(IUNIT,*) 'Sub catchment number ',INB_RIVJUNC,&
         'lambX=',(IX0*XDXT(KCAT))+XX0(KCAT),'lamby=',(IY0*XDXT(KCAT))+XY0(KCAT)
    !
    ENDIF!JCE>3 <=> pt d'intersection
    !
  ENDIF!XDHIL=0 <=> pt riviere
  50 CONTINUE
  !
ENDDO

CALL CLOSE_FILE('ASCII ',IUNIT)

DO JSB=1,INB_RIVJUNC
  write(*,*) 'subcatch',JSB,'nb pix=',COUNT(IMASK_UP(:)==JSB),&
           'area=',COUNT(IMASK_UP(:)==JSB)*XDXT(KCAT)*XDXT(KCAT)/1000000.
ENDDO
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
!
CALL OPEN_FILE('ASCII ',IUNIT,YFILE,'FORMATTED',HACTION='WRITE')
!
DO JJ=1,NNCAT_MAX
  WRITE(IUNIT,*) XDIST_OUTLET(KCAT,JJ,:)
ENDDO
!
CALL CLOSE_FILE('ASCII ',IUNIT)
!

IF (LHOOK) CALL DR_HOOK('FIND_SUB_CATCH',1,ZHOOK_HANDLE)
!-----------------------------------------------------------------!
!
END SUBROUTINE FIND_SUB_CATCH
