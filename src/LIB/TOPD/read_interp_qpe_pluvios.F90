!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!     #########################################
      SUBROUTINE READ_INTERP_QPE_PLUVIOS(YSC,HFILE,PFIELD_OUT)
!     #######################
!
!!****  *READ_INTERP_QPE_PLUVIOS*  
!!
!!    PURPOSE
!!    -------
!!     This routine aims at reading forcing variables from ascii files in llv
!      format and interpolates values on SURFEX domain 
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
USE MODD_SURFEX_n, ONLY : SURFEX_t
!
USE MODI_GET_LUOUT
USE MODI_OPEN_FILE
USE MODI_CLOSE_FILE
USE MODI_HOR_INTERPOL
!
USE MODD_PREP,             ONLY : CINGRID_TYPE, CINTERP_TYPE, LINTERP
USE MODD_GRID_GRIB,        ONLY : NNI
USE MODD_TYPE_DATE_SURF
!
USE MODD_SURF_PAR,         ONLY : XUNDEF
USE MODD_PGD_GRID,         ONLY : NL, LLATLONMASK
USE MODI_LATLONMASK
USE MODD_GRID_LATLONREGUL, ONLY : XILATARRAY,XILONARRAY
USE MODD_PGDWORK,          ONLY : XSUMVAL, NSIZE,CATYPE
!
USE MODI_TREAT_FIELD
USE MODI_INTERPOL_FIELD
USE MODI_PACK_SAME_RANK
USE MODI_READ_GRID
USE MODI_READ_SURF
USE MODI_PGD_GRIDTYPE_INIT
!
USE MODE_GRIDTYPE_LONLATVAL
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
TYPE(SURFEX_t), INTENT(INOUT) :: YSC
!
CHARACTER(LEN=*),               INTENT(IN)  :: HFILE       ! radar file to be read
REAL, DIMENSION(:),ALLOCATABLE, INTENT(OUT) :: PFIELD_OUT  ! Forcing parameter read on file
!!
!*      0.2    declarations of local variables
!
INTEGER                   :: ILUOUT      ! Unit of the files
!
REAL, DIMENSION(:,:), ALLOCATABLE      :: ZFLD_INTERPOL !Field interpolated
REAL, DIMENSION(:,:), POINTER          :: ZPNT_INTERPOL
REAL, DIMENSION(:,:), ALLOCATABLE      :: ZFIELD    ! physiographic field on full grid
!
CHARACTER(LEN=20)   :: YFIELD
LOGICAL, DIMENSION(720,360)     :: LLALO
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('READ_INTERP_QPE_PLUVIOS',0,ZHOOK_HANDLE)
!
!*       1.0    preparing file openning
!               ----------------------
CALL GET_LUOUT('ASCII ',ILUOUT)
YSC%UG%G%NGRID_PAR=SIZE(YSC%UG%G%XGRID_PAR,1)
!
CALL PGD_GRIDTYPE_INIT(YSC%UG%G%CGRID,YSC%UG%G%NGRID_PAR,YSC%UG%G%XGRID_PAR)
CINTERP_TYPE='HORIBL'
!
CALL LATLONMASK(YSC%UG%G%CGRID,YSC%UG%G%NGRID_PAR,YSC%UG%G%XGRID_PAR,LLALO)
LLATLONMASK=LLALO
!
NL=SIZE(YSC%UG%G%XLAT,1)
!
ALLOCATE(ZFIELD(NL,1))
ALLOCATE(NSIZE     (NL,1))
ALLOCATE(XSUMVAL   (NL,1))
ZFIELD(:,:) = XUNDEF
NSIZE    (:,:) = 0.
XSUMVAL  (:,:) = 0.
!
YFIELD = '                    '
!
!*       2.0    reading file
!               ----------------------
CALL TREAT_FIELD(YSC%UG, YSC%U, YSC%USS,'ASCII ','SURF  ','ASCLLV','A_MESH',HFILE,YFIELD,ZFIELD)
!
WHERE ((ZFIELD>=999.).OR.(ZFIELD<=-999.))
  ZFIELD(:,:)=0.
  NSIZE(:,:)=0
ENDWHERE
!
CALL INTERPOL_FIELD(YSC%UG, YSC%U,'ASCII ',ILUOUT,NSIZE(:,1),ZFIELD(:,1),YFIELD,PDEF=XUNDEF)
!
write(*,*) 'interp ok read_interp_pluvios'
write(*,*) 'min,max=',MINVAL(ZFIELD),MAXVAL(ZFIELD,ZFIELD/=XUNDEF)
write(*,*) '0.,undef=',COUNT(ZFIELD==0.0),COUNT(ZFIELD==XUNDEF)
!
!*       3.0    interpolation on good grid
!               ----------------------
IF (.NOT.ALLOCATED(PFIELD_OUT)) ALLOCATE (PFIELD_OUT(YSC%U%NDIM_FULL))
PFIELD_OUT(:)=XUNDEF
!
ALLOCATE(ZPNT_INTERPOL(SIZE(ZFIELD,1),1))
ZPNT_INTERPOL(:,1)=ZFIELD(:,1)
! ATTENTION, unité=1/10 de mm.
IF ((INDEX(HFILE,'OHMCV')==0).AND.(INDEX(HFILE,'WRFLLV')==0)) THEN
  WHERE (ZFIELD(:,1)/=XUNDEF)  PFIELD_OUT(:)=ZFIELD(:,1)/10.
ELSE 
  WHERE (ZFIELD(:,1)/=XUNDEF)  PFIELD_OUT(:)=ZFIELD(:,1)
ENDIF
!
DEALLOCATE(ZFIELD)
!
IF (.NOT.ALLOCATED(ZFLD_INTERPOL))ALLOCATE(ZFLD_INTERPOL(YSC%U%NDIM_FULL,1))
XILATARRAY=YSC%UG%G%XLAT
XILONARRAY=YSC%UG%G%XLON
!
DEALLOCATE(NSIZE    )
DEALLOCATE(XSUMVAL  )
!
IF (LHOOK) CALL DR_HOOK('READ_INTERP_QPE_PLUVIOS',1,ZHOOK_HANDLE)
!
END SUBROUTINE READ_INTERP_QPE_PLUVIOS
