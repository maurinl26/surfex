SUBROUTINE SNOW_DEPLETION_PF(TPSNOW, PWORKPFM2, II, IENS,  IVAR, PSCF)
!##########################
!
!! *SNOW_DEPLETION_PF*
!!
!! Computing Model Equivalent Snow Cover Fraction from model SWE values.
!! In order to assimilate it with ModImLab Snow Cover Fraction Product.
!!
!!    REFERENCE
!!    -------------
!!    Zaitchik et al.,  2009
!!    De Lannoy et al., 2012
!!    Thirel et al.,    2013
!!                                                  
!!    AUTHOR
!!    ------
!!    B. Cluzet       * Meteo France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    18/04/2018

! -----------------------------------------------------------------------------
!
! First version : B. Cluzet, 04/2018
!  
! -----------------------------------------------------------------------------
!
USE MODD_ASSIM,         ONLY : XF, XTAUSCF, XCEILSCF ! global variable to write model SCF values
USE MODD_TYPE_SNOW,     ONLY : SURF_SNOW

IMPLICIT NONE

TYPE(SURF_SNOW), INTENT(INOUT) :: TPSNOW 
REAL, DIMENSION(:,:,:), INTENT(IN) :: PWORKPFM2
INTEGER, INTENT(IN) :: II                    !grid point
INTEGER, INTENT(IN) :: IENS 
INTEGER, INTENT(IN) :: IVAR
REAL, INTENT(INOUT)  :: PSCF                  
!*    0.     Declaration of local variables
!            ------------------------------
!
REAL   :: ZBULKSWE ! bulk swe
REAL   :: ZSCF     ! corresponding snow cover fraction
!REAL, DIMENSION(SIZE(PWSNOW,1), SIZE(PWSNOW,3)) :: ZBULKSWE ! bulk swe
!REAL, DIMENSION(SIZE(PWSNOW,1), SIZE(PWSNOW,3)) :: ZSCF     ! corresponding snow cover fraction

INTEGER :: JJ,JI


!******************************************************************************************

! Computation of bulk SWE
!ZBULKSWE = SUM(TPSNOW%WSNOW, 2)
ZBULKSWE = SUM(TPSNOW%WSNOW(II,:))
!PRINT*, SHAPE(PWSNOW), 'shppwsnow'
PRINT*, SHAPE(ZSCF), 'shpzscf'

ZSCF = MIN( 1. - (EXP(-XTAUSCF * (ZBULKSWE/XCEILSCF) - (ZBULKSWE/XCEILSCF)*EXP(-XTAUSCF))),&
           1. &
           )


PSCF = ZSCF
!PRINT*, XF(:,:,IENS, IVAR), 'xfffifififi'

END SUBROUTINE SNOW_DEPLETION_PF






