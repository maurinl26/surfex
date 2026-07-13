SUBROUTINE COUNT_SNOW_VAR_PF()
!##########################
!
!! *COUNT_SNOW_VAR_PF*
!!
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
! Counting of number of chosen MODIS bands for assimilation and flags for SCA/MODIS
!  - NRETM=0    ! 1 if at least one modis band in CVAR or band ratio
!  - NCOUNTM=0  ! number of modis bands/band ratio
!  - NRETS=0    ! 1 if SWE ('PSB') in CVAR
!  - NRETD=0   ! 1 if assim snow depth
!
! First version : B. Cluzet, 04/2018
!  
! -----------------------------------------------------------------------------
!

USE MODD_ASSIM,         ONLY : CVAR, NVAR, NRETM, NCOUNTM, NRETS, NRETD, NRETR ! CVAR holds list of assimilation variables names (model)
                                       ! good practice : 'PB1' , ... , 'PBN' , 'PSB'
                                       !                 |-modis bands first-|, swe after
                                       
                                       
IMPLICIT NONE                                       

INTEGER :: JVAR
! -----------------------------------------------------------
!


! ******************************************************************************************
! init
NRETM  =0
NRETS  =0
NCOUNTM=0
NRETD  =0
NRETR  =0
!
! search
DO JVAR = 1, NVAR
  NRETM=INDEX(CVAR(JVAR), 'PB')            ! on cherche pb dans cvar(jvar)
  NRETR=INDEX(CVAR(JVAR), 'R')             ! searching for band ratios in CVAR
  IF ((NRETM==1) .OR. (NRETR==1)) NCOUNTM = NCOUNTM+1
  IF (NRETS/=1) THEN
    NRETS=INDEX(CVAR(JVAR), 'PSB')
  ENDIF
  IF(CVAR(JVAR)=="DEP") THEN
    NRETD=1
  ENDIF
END DO
IF (NCOUNTM>0) THEN
  NRETM=1
ELSE
  NRETM=0
ENDIF

! ensuring that the user is not trying to assimilate MODIS and snow depth data at the same time
IF ((NRETM==1) .AND. (NRETD==1)) THEN
  CALL ABOR1_SFX('Attempt to assimilate reflectance AND snow height. Not possible for now.')
ENDIF
! 
END SUBROUTINE COUNT_SNOW_VAR_PF
