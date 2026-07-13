!     #########
      SUBROUTINE JULIAN (KYEAR, KMONTH, KDAY, PTIME, &
                         PZDATE)
!     ####################################################################################
!
!!****  *SUNPOS * - routine to compute the julian day of year
!!
!!    
!!     
!!    EXTERNAL
!!    --------
!!      NONE
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!    REFERENCE
!!    ---------
!!      "Radiative Processes in Meteorology and Climatology"  
!!                          (1976)   Paltridge and Platt 
!!
!!    AUTHOR
!!    ------
!!      J.-P. Pinty      * Laboratoire d'Aerologie*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original             16/10/94 
!!      Revised              12/09/95
!!      (J.Stein)            01:04/96  bug correction for ZZEANG     
!!      (K. Suhre)           14/02/97  bug correction for ZLON0     
!!      (V. Masson)          01/03/03  add zenithal angle output
!!      (V. Masson)          14/03/14  avoid discontinuous declination at 00UTC each day
!! 	(M. Dumont)      09/11/15 adapted to return julian day
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_CSTS,          ONLY : XPI, XDAY
USE MODD_SURFEX_OMP, ONLY : NBLOCK, NBLOCKTOT
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
#ifdef AIX64
!$ USE OMP_LIB
#endif
!
IMPLICIT NONE
!
#ifndef AIX64
!$ INCLUDE 'omp_lib.h'
#endif
!
!*       0.1   Declarations of dummy arguments :
!

INTEGER,                      INTENT(IN)   :: KYEAR      ! current year                        
INTEGER,                      INTENT(IN)   :: KMONTH     ! current month                        
INTEGER,                      INTENT(IN)   :: KDAY       ! current day                        
REAL,                         INTENT(IN)   :: PTIME      ! current time                        

!
REAL,            INTENT(OUT)  :: PZDATE    ! Julian day of year 

!
!*       0.2   declarations of local variables
!
!
REAL                                       :: ZUT        ! Universal time
!

INTEGER, DIMENSION(0:11)                   :: IBIS, INOBIS ! Cumulative number of days per month
                                                           ! for bissextile and regular years

!                                            
INTEGER                                    :: JI, JJ

REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------
!
!*       1.    TO COMPUTE THE TRUE SOLAR TIME
!              -------------------------------
!
IF (LHOOK) CALL DR_HOOK('SUNPOS',0,ZHOOK_HANDLE)
ZUT  = MOD( 24.0+MOD(PTIME/3600.,24.0),24.0 )

INOBIS(:) = (/0,31,59,90,120,151,181,212,243,273,304,334/)
IBIS(0:1) = INOBIS(0:1)
DO JI=2,11
  IBIS(JI) = INOBIS(JI)+1
END DO
IF( MOD(KYEAR,4).EQ.0 .AND. (MOD(KYEAR,100).NE.0 .OR. MOD(KYEAR,400).EQ.0)) THEN
  PZDATE = FLOAT(KDAY +   IBIS(KMONTH-1)) - 1 + PTIME/XDAY

ELSE
  PZDATE = FLOAT(KDAY + INOBIS(KMONTH-1)) - 1 + PTIME/XDAY

END IF

!
IF (LHOOK) CALL DR_HOOK('SUNPOS',1,ZHOOK_HANDLE)
!-------------------------------------------------------------------------------
!
END SUBROUTINE JULIAN
