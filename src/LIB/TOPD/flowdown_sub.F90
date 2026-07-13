!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!-----------------------------------------------------------------
!     ###################
      SUBROUTINE FLOWDOWN_SUB(KNMC,PVAR,PCONN,KLINE,KSUB,KCAT)
!     ###################
!
!!****  *FLOWDOWN*
!
!!    PURPOSE
!!    -------
! to propagate data between pixels of a catchment in function of its topography
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
!!      K. Chancibault	* CNRM / Meteo-France *
!!      G-M Saulnier    * LTHE *
!!
!!    MODIFICATIONS
!!    -------------
!!
!!      Original   14/01/2005
!-------------------------------------------------------------------------------
!
!*       0.     DECLARATIONS
!               ------------
!
USE MODD_TOPODYN      ,ONLY : NCAT_CAT_TO_SUB,NPIX_CAT_TO_SUB,NPIX_SUB_TO_CAT
USE MODD_SURF_PAR     ,ONLY : NUNDEF
!
USE YOMHOOK   ,ONLY : LHOOK, DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*      0.1     declarations of arguments
!
INTEGER, INTENT(IN) :: KNMC  ! catchment grid points number
REAL, DIMENSION(:), INTENT(INOUT)   :: PVAR  ! variable to propagate
REAL, DIMENSION(:,:), INTENT(IN)    :: PCONN ! catchment grid points connections
INTEGER, DIMENSION(:), INTENT(IN)   :: KLINE ! 
INTEGER, INTENT(IN) :: KSUB  ! Sub-catchment number
INTEGER, INTENT(IN) :: KCAT  ! Catchment  number
!
!*      0.2    declarations of local variables
!
INTEGER                  :: JJ, JI ! work variables
INTEGER                  :: JNUP  ! number of upslope pixels
INTEGER                  :: JCOL  ! third index of the pixel in the array XCONN
INTEGER                  :: JREF  ! index of the upslope pixel in the topo domain
REAL                     :: ZFAC  ! propagation factor between this pixel and the
                                  ! upslope one
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('FLOWDOWN_SUB',0,ZHOOK_HANDLE)
!
DO JJ=1,KNMC
  IF (NPIX_SUB_TO_CAT(KCAT,KSUB,JJ)/=NUNDEF) THEN
    JNUP = INT(PCONN(NPIX_SUB_TO_CAT(KCAT,KSUB,JJ),4))
    DO JI=1,JNUP
      JCOL = ((JI-1)*2) + 5
      JREF = INT(PCONN(NPIX_SUB_TO_CAT(KCAT,KSUB,JJ),JCOL))
      ZFAC = PCONN(NPIX_SUB_TO_CAT(KCAT,KSUB,JJ),JCOL+1)
      IF (NCAT_CAT_TO_SUB(KCAT,KLINE(JREF))==KSUB) &
      PVAR(JJ) = PVAR(JJ) + PVAR(NPIX_CAT_TO_SUB(KCAT,KLINE(JREF))) * ZFAC
    ENDDO
  ENDIF
ENDDO
!
IF (LHOOK) CALL DR_HOOK('FLOWDOWN_SUB',1,ZHOOK_HANDLE)
!
END SUBROUTINE FLOWDOWN_SUB
