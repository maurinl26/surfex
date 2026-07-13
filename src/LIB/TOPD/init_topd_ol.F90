!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!-----------------------------------------------------------------
!     #######################
      SUBROUTINE INIT_TOPD_OL(HPROGRAM)
!     #######################
!
!!****  *INIT_TOPD_OL*  
!!
!!    PURPOSE
!!    -------
!     This routine aims at initialising the variables 
!     needed of running Topmodel for OFFLINE step.
!              
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
!!      Original   03/2014
!-------------------------------------------------------------------------------
!
!*       0.     DECLARATIONS
!               ------------
!
USE MODD_COUPLING_TOPD, ONLY : NNB_STP_RESTART, LSUBCAT
USE MODD_TOPODYN,       ONLY : CCAT, NNCAT, NNB_TOPD_STEP, XTOPD_STEP,&
                               XDXT,NNPT,NX_STEP_ROUT, XSPEEDR,XSPEEDH,&
                               NNMC, NMESHT, NPMAX,NLINE, XDRIV, XDHIL,&
                               XTIME_TOPD, XTIME_TOPD_SUB, XDGRD, &
                               XSPEEDG, XTIME_TOPD_DRAIN, XTIME_TOPD_DRAIN_SUB,&
                               XTANB,XLAMBDA,XQB_DR, XQB_RUN, XQB_DR_SUB,&
                               XQB_RUN_SUB, NRIV_JUNCTION,XDRIV_SUBCAT,&
                               NCAT_CAT_TO_SUB,NPIX_CAT_TO_SUB,NPIX_SUB_TO_CAT,&
                               NNCAT_MAX,NNCAT_SUB,NMASK_OUTLET,&
                               XQTOT_SUB,NBPIX_SUB,XDIST_OUTLET, &
                               XA_SPEED, XB_SPEED, XMAX_SPEED,LSPEEDR_VAR
!
USE MODD_SURF_PAR,       ONLY : XUNDEF, NUNDEF
!
USE MODI_GET_LUOUT
USE MODI_INIT_TOPD
USE MODI_READ_TOPD_HEADER_DTM
USE MODI_READ_TOPD_FILE
USE MODI_READ_TOPD_HEADER_CONNEX
USE MODI_READ_CONNEX_FILE
USE MODI_OPEN_FILE
USE MODI_CLOSE_FILE
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
 CHARACTER(LEN=*), INTENT(IN) :: HPROGRAM    !
!
!*      0.2    declarations of local variables
!
!
INTEGER                   :: JI,JJ,JCAT,JO ! loop control 
INTEGER                   :: IOVER                  ! Unit of the files
INTEGER                   :: ILUOUT                 ! Unit of the files
INTEGER                   :: ISUBCAT, ISUBPIX
!
CHARACTER(LEN=50), DIMENSION(NNCAT) :: YFILE_SUBCAT   !  file names
CHARACTER(LEN=50), DIMENSION(NNCAT) :: YFILE_SUBPIX   !  file names
CHARACTER(LEN=50), DIMENSION(NNCAT) :: YFILE_OUTLET   !  file names
CHARACTER(LEN=28) :: YFILE
!
REAL, DIMENSION(:),ALLOCATABLE    :: ZSUBCAT ! variable read
REAL, DIMENSION(:),ALLOCATABLE    :: ZSUBPIX ! variable read
REAL, DIMENSION(:),ALLOCATABLE    :: ZOUTLET ! variable read
!
REAL :: ZDHIL  ! distance along slope
REAL :: ZDRIV  ! distance along rivers
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('INIT_TOPD_OL',0,ZHOOK_HANDLE)
!


!*       1    Initialization:
!               ---------------
!
 CALL GET_LUOUT(HPROGRAM,ILUOUT)
!
WRITE(ILUOUT,*) 'INITIALISATION INIT_TOPD_OL'
!
CALL INIT_TOPD('ASCII ')
!
  !*      2       Calculations for routing by geomorpho
  !               -------------------------------------
  !         2.1 Computation of A and B Eram parameters
  !         ------------------------------------------
IF (LSPEEDR_VAR) THEN 
  ALLOCATE(XA_SPEED(NNCAT))
  ALLOCATE(XB_SPEED(NNCAT))
  XA_SPEED(:)=0.63
  XB_SPEED(:)=0.082
  !
  DO JCAT=1,NNCAT
    XA_SPEED(JCAT)=0.0377*XLAMBDA(JCAT,NNMC(JCAT))+0.00519*100*XTANB(JCAT,NNMC(JCAT))
    XB_SPEED(JCAT)=0.0877*XLAMBDA(JCAT,NNMC(JCAT))-0.0000805*NNMC(JCAT)*XDXT(JCAT)*XDXT(JCAT)*XDXT(JCAT)/1000000.
  ENDDO
  !
ENDIF
!
ALLOCATE(NX_STEP_ROUT(NNCAT))
ALLOCATE(XTIME_TOPD(NNCAT,NMESHT))
ALLOCATE(XTIME_TOPD_DRAIN(NNCAT,NMESHT))
!
XTIME_TOPD(:,:) = 0.0
XTIME_TOPD_DRAIN(:,:) = 0.0
!
DO JCAT=1,NNCAT
  !
  XMAX_SPEED(JCAT)=XSPEEDR(JCAT)
  !
  IF ( XSPEEDR(JCAT)/=0. .AND. XSPEEDG(JCAT)/=0. ) THEN
    !
    DO JJ=1,NNMC(JCAT)
      IF ( XDHIL(JCAT,JJ)/=XUNDEF .AND. XDRIV(JCAT,JJ)/=XUNDEF ) THEN
        XTIME_TOPD(JCAT,JJ) = XDHIL(JCAT,JJ) / XSPEEDH(JCAT) + &
                                        XDRIV(JCAT,JJ) / XSPEEDR(JCAT)
      ENDIF
      IF ( XDGRD(JCAT,JJ)/=XUNDEF .AND. XDRIV(JCAT,JJ)/=XUNDEF ) THEN
        XTIME_TOPD_DRAIN(JCAT,JJ) = XDGRD(JCAT,JJ) / XSPEEDG(JCAT) + &
                                              XDRIV(JCAT,JJ) / XSPEEDR(JCAT)
      ENDIF
    ENDDO
    !
  ELSE 
    WRITE(ILUOUT,*) 'You have to choose some values for routing velocities'
  ENDIF
  !
  IF (XTOPD_STEP/=0.) &
    NX_STEP_ROUT(JCAT) = INT(MAXVAL(XTIME_TOPD(JCAT,1:NNMC(JCAT))) / XTOPD_STEP) + 1
  !
ENDDO
!
IF ( NNB_STP_RESTART==0 ) NNB_STP_RESTART = MAX(NNB_TOPD_STEP,MAXVAL(NX_STEP_ROUT(:)))
!
ALLOCATE(XQB_DR(NNCAT,NNB_TOPD_STEP))
XQB_DR(:,:)=0.0
ALLOCATE(XQB_RUN(NNCAT,NNB_TOPD_STEP))
XQB_RUN(:,:)=0.0
!
!
!*      3.       Sub-catchments treatement
!               -------------------------------------
!
IF (LSUBCAT) THEN
  ALLOCATE(NPIX_SUB_TO_CAT(NNCAT,NNCAT_MAX,NPMAX))
  NPIX_SUB_TO_CAT(:,:,:) = NUNDEF
  !
  ALLOCATE(NBPIX_SUB(NNCAT,NNCAT_MAX))
  NBPIX_SUB(:,:) = NUNDEF
  !
  ALLOCATE(XQTOT_SUB(NNCAT,NNCAT_MAX,NNB_TOPD_STEP))
  XQTOT_SUB(:,:,:)=0.0
  ALLOCATE(XQB_RUN_SUB(NNCAT,NNCAT_MAX,NNB_TOPD_STEP))
  XQB_RUN_SUB(:,:,:)=0.0
  ALLOCATE(XQB_DR_SUB(NNCAT,NNCAT_MAX,NNB_TOPD_STEP))
  XQB_DR_SUB(:,:,:)=0.0
  !
  ALLOCATE(NNCAT_SUB(NNCAT))
  !
  ALLOCATE(NRIV_JUNCTION(NNCAT,NPMAX))
  ALLOCATE(NCAT_CAT_TO_SUB(NNCAT,NPMAX))
  NCAT_CAT_TO_SUB(:,:) = NUNDEF
  ALLOCATE(NPIX_CAT_TO_SUB(NNCAT,NPMAX))
  NPIX_CAT_TO_SUB(:,:) = NUNDEF
  ALLOCATE(NMASK_OUTLET(NNCAT,NPMAX))
  NMASK_OUTLET(:,:) = NUNDEF
  ALLOCATE(XDRIV_SUBCAT(NNCAT,NPMAX))
  XDRIV_SUBCAT(:,:)=0.0
  ALLOCATE(XTIME_TOPD_SUB(NNCAT,NMESHT))
  ALLOCATE(XTIME_TOPD_DRAIN_SUB(NNCAT,NMESHT))
  !
  !
  !*      4.0       Reading Sub-catchments numbers
  !               -------------------------------------
  DO JCAT=1,NNCAT
    ALLOCATE(ZSUBCAT(NNPT(JCAT)))
    ALLOCATE(ZSUBPIX(NNPT(JCAT)))
    ALLOCATE(ZOUTLET(NNPT(JCAT)))
    !
    YFILE_SUBCAT(JCAT)=TRIM(CCAT(JCAT))//'SubCat.map'
    write(*,*) YFILE_SUBCAT(JCAT)
    CALL READ_TOPD_FILE(HPROGRAM,YFILE_SUBCAT(JCAT),'FORMATTED',NNPT(JCAT),ZSUBCAT)
    DO JJ=1,NNPT(JCAT)
      IF ( NLINE(JCAT,JJ)/=0. .AND. NLINE(JCAT,JJ)/=XUNDEF .AND. ZSUBCAT(JJ)/=XUNDEF) &
        NCAT_CAT_TO_SUB(JCAT,NLINE(JCAT,JJ)) = INT(ZSUBCAT(JJ))
    ENDDO
    NNCAT_SUB(JCAT)=MAXVAL(NCAT_CAT_TO_SUB(JCAT,:),MASK=NCAT_CAT_TO_SUB(JCAT,:)/=NUNDEF)
    DO JJ = 1,NNCAT_SUB(JCAT)
      NBPIX_SUB(JCAT,JJ)=COUNT(NCAT_CAT_TO_SUB(JCAT,:)==JJ)
    ENDDO
   !
   !*      4.1       Reading Pixels of Sub-catchments
   !               -------------------------------------
   !
    YFILE_SUBPIX(JCAT)=TRIM(CCAT(JCAT))//'SubPix.map'
    write(*,*) YFILE_SUBPIX(JCAT)
    CALL READ_TOPD_FILE(HPROGRAM,YFILE_SUBPIX(JCAT),'FORMATTED',NNPT(JCAT),ZSUBPIX)
    DO JJ=1,NNPT(JCAT)
      !
      IF ( NLINE(JCAT,JJ)/=0 .AND.NLINE(JCAT,JJ)/=XUNDEF ) THEN
        ISUBCAT=NCAT_CAT_TO_SUB(JCAT,NLINE(JCAT,JJ))
        ISUBPIX=INT(ZSUBPIX(JJ))
        !
        IF ( ZSUBPIX(JJ)/=XUNDEF.AND.ISUBCAT/=NUNDEF.AND.ISUBPIX/=NUNDEF )THEN
          NPIX_CAT_TO_SUB(JCAT,NLINE(JCAT,JJ)) = ISUBPIX
          NPIX_SUB_TO_CAT(JCAT,ISUBCAT,ISUBPIX)=NLINE(JCAT,JJ)
        ENDIF
        !
      ENDIF
      !
    ENDDO
    !
    !*      4.2       Reading Outlets of Sub-catchments
    !               -------------------------------------
    !
    YFILE_OUTLET(JCAT)=TRIM(CCAT(JCAT))//'Outlet.map'
    write(*,*) YFILE_OUTLET(JCAT)
    CALL READ_TOPD_FILE(HPROGRAM,YFILE_OUTLET(JCAT),'FORMATTED',NNPT(JCAT),ZOUTLET)
    DO JJ=1,NPMAX
      IF ( NLINE(JCAT,JJ)/=0. .AND. NLINE(JCAT,JJ)/=XUNDEF .AND. ZOUTLET(JJ)/=XUNDEF) &
        NMASK_OUTLET(JCAT,NLINE(JCAT,JJ)) = INT(ZOUTLET(JJ))
    ENDDO
    !
    !*      4.3       Reading Overlapped Sub-catchments
    !               -------------------------------------
    YFILE = TRIM(CCAT(JCAT))//'_overlap'
    write(*,*) YFILE
    CALL OPEN_FILE('ASCII ',IOVER,YFILE,'FORMATTED',HACTION='READ')
    DO JJ=1,NNCAT_MAX
      READ(IOVER,*,END=110) XDIST_OUTLET(JCAT,JJ,:)
    ENDDO
110    CALL CLOSE_FILE('ASCII ',IOVER)
   
    !
    DO JJ=1,NNCAT_MAX
      DO JI=1,NNCAT_MAX
        IF (XDIST_OUTLET(JCAT,JJ,JI)/=XUNDEF)&
          write(*,*)JCAT,'XDIST_OUTLET(',JJ,',',JI,')= ',XDIST_OUTLET(JCAT,JJ,JI)
      ENDDO
    ENDDO
    !
    DO JJ=1,NNMC(JCAT)
      JO=NMASK_OUTLET(JCAT,JJ)
      ZDHIL=XDHIL(JCAT,JJ)
      ZDRIV=XDRIV(JCAT,JJ)-XDRIV(JCAT,JO)
      !
      XTIME_TOPD_SUB(JCAT,JJ) = ZDHIL/ XSPEEDH(JCAT) + ZDRIV/ XSPEEDR(JCAT)
      XTIME_TOPD_DRAIN_SUB(JCAT,JJ) = ZDHIL/ XSPEEDG(JCAT) + ZDRIV/ XSPEEDR(JCAT)
    ENDDO!JJ
    !
    DEALLOCATE(ZSUBCAT)
    DEALLOCATE(ZSUBPIX)
    DEALLOCATE(ZOUTLET)
  ENDDO !JCAT
ENDIF !LSUBCAT
!
IF (LHOOK) CALL DR_HOOK('INIT_TOPD_OL',1,ZHOOK_HANDLE)
!
END SUBROUTINE INIT_TOPD_OL
