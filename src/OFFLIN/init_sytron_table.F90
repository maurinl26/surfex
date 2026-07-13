SUBROUTINE INIT_SYTRON_TABLE(USS,PZS,KI,PLAT,PLON)

!
!!**** *INIT_SYTRON_TABLE
!!                         Initialize table that defines relationship between grid points 
!!                         for Sytron simulation
!!                         Two opposites points belongs to the same massif
!!                         have the same elevation and slope but an opposite aspect 
!!
!!    PURPOSE
!!    -------
!!
!!    METHOD
!!    ------
!!   
!
!!    EXTERNAL
!!    --------
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!    REFERENCE
!!    ---------
!!   
!!
!!    AUTHOR
!!    ------
!!
!!    V. Vionnet        Meteo-France
!!
!!    MODIFICATION
!!    ------------
!!
!!    Original    09/14
!!
!!
!----------------------------------------------------------------------------
!
!*    0.     DECLARATION
!            -----------
!
USE MODD_SSO_n, ONLY : SSO_t
USE MODD_SYTRON_PAR
USE MODD_SURF_PAR,       ONLY : XUNDEF
USE MODD_CSTS,           ONLY : XPI
!
IMPLICIT NONE
!
#ifdef SFX_MPI
INCLUDE "mpif.h"
#endif


!*    0.1    Declaration of arguments
!            ------------------------
!
TYPE(SSO_t), INTENT(IN) :: USS
INTEGER,              INTENT(IN)  :: KI        ! number of points
REAL, DIMENSION(KI),   INTENT(IN) :: PZS      ! orography of this MPI thread (or total domain if Open MP)
REAL,DIMENSION(:),INTENT(IN):: PLAT ! latitudes
REAL,DIMENSION(:),INTENT(IN):: PLON ! longitudes

!
!
!*    0.2    Declaration of local variables
!            ------------------------------
!
INTEGER :: JI, JJ,JS
LOGICAL :: GDISTRIB ! TRUE is the point is concerned by SYTRON
REAL    :: ZASP_OPPOSITE ! Opposite aspect
REAL    :: ZSLOPE_DEG,ZSLOPE_DEG2    ! Slopes in deg
!------------------------------------------------------------------------------------------
!
!*    1.0    Identify opposite points
!            ------------------------------
!
ALLOCATE(NTAB_SYT (KI))
NTAB_SYT(:)=-999 

DO JI=1,KI

GDISTRIB=.FALSE.
ZSLOPE_DEG=180./XPI*ATAN(USS%XSSO_SLOPE(JI))
DO JS=1,SIZE(XSLOPE_SYT,1)
     IF(NINT(ZSLOPE_DEG)==NINT(XSLOPE_SYT(JS))) THEN
          GDISTRIB=.TRUE.
          EXIT 
      ENDIF
ENDDO

!If point JI not concerned by SYTRON computation, exit loop
IF(.NOT. GDISTRIB) THEN 
   NTAB_SYT(JI)=JI
ELSE
DO JJ=1,KI

IF(NTAB_SYT(JJ) ==-999) THEN  ! Point JJ has not been treated yet

ZSLOPE_DEG2=180./XPI*ATAN(USS%XSSO_SLOPE(JJ))
IF(PLAT(JJ)==PLAT(JI) .AND. PLON(JJ)==PLON(JI)) THEN ! Points JJ and JI belong to the same massif
    IF(PZS(JJ) == PZS(JI) .AND. NINT(ZSLOPE_DEG)==NINT(ZSLOPE_DEG2)) THEN   ! Points JJ and JI have the same elevation and the same slope
          ZASP_OPPOSITE = USS%XSSO_DIR(JJ)+180.
          IF(ZASP_OPPOSITE > 360.) ZASP_OPPOSITE=ZASP_OPPOSITE-360.

          IF(ZASP_OPPOSITE==USS%XSSO_DIR(JI)) THEN  ! Point JJ is the opposite point of JI
               NTAB_SYT(JI)=JJ
               NTAB_SYT(JJ)=JI
               EXIT
          ENDIF 
    ENDIF
ENDIF 

ENDIF
ENDDO ! JJ
ENDIF
  
ENDDO  ! JI  

END SUBROUTINE INIT_SYTRON_TABLE

