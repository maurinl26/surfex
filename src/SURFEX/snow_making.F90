!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!     #########
SUBROUTINE SNOW_MAKING(IO, TPTIME, PTSTEP, PEK, PSR, PTA, PVMOD, PPS, PQA, &
                        PRHOA, PSNOW, PSNOWFALL, PSNOWMAK, PDMKXWBT)



!!    AUTHOR
!!    ------
!!	  P. Sandre           * Meteo-France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original        3/14  P. Spandre
!!      Modified by P.Spandre      (06/2014):
!!      Modified by C. Carmagnola  (11/2018): 
!!        - temperature-dependent snow production rate (according to Hanzer et al., 2014)
!!        - 3 new parameters for snow production in namelist (XPR_A, XPR_B, XPT)
!!        - new diagnostic variable (wet bulb temperature, WBT)
!!        - new prognostic variable (water consumption for snowmaking, MMP)
!!        - when LSELF_PROD=.FALSE., XPROD_SCHEME is the target MMP and does not depend on atmospheric or timing conditions


USE MODD_SNOW_PAR,       ONLY : XRHOSMAX_ES, XPSR_SNOWMAK, XPROD_SCHEME, &
                                XPP_D1, XPP_D2, XPP_D3, XPP_H1, &
                                XPP_H2, XPP_H3, XPP_H4, XPR_A, XPR_B,         &
                                XPTA_SEUIL, XWT, XPT, XPTR, XRHO_SNOWMAK

USE MODD_TYPE_DATE_SURF, ONLY: DATE_TIME

USE MODD_DIAG_MISC_ISBA_n, ONLY : DIAG_MISC_ISBA_t

USE MODD_ISBA_OPTIONS_n, ONLY : ISBA_OPTIONS_t

USE MODD_ISBA_n, ONLY : ISBA_PE_t

IMPLICIT NONE
                            

!*      0.1    declarations of arguments
TYPE(ISBA_OPTIONS_t), INTENT(INOUT) :: IO

TYPE(DATE_TIME), INTENT(IN)         :: TPTIME     ! current date and time
REAL, INTENT(IN)                    :: PTSTEP
!                                      PTSTEP    = time step of the integration

TYPE(ISBA_PE_t), INTENT(INOUT)      :: PEK

REAL, DIMENSION(:), INTENT(IN)      :: PSR, PTA, PVMOD, PPS, PQA
!                                      PSR    = snow rate (SWE) [kg/(m2 s)]
!                                      PTA    = atmospheric temperature at level za (K)
!                                      PVMOD  = modulus of the wind parallel to the orography (m/s)
!                                      PPS    = surface pressure
!                                      PQA    = atmospheric specific humidity
!                                                at level za

REAL, DIMENSION(:), INTENT(IN)      :: PRHOA
!                                      PRHOA     = air density


REAL, DIMENSION(SIZE(PTA)), INTENT(IN)    :: PSNOW
REAL, DIMENSION(SIZE(PTA)), INTENT(INOUT) :: PSNOWFALL
!                                            PSNOW        = snow depth (m) 
!                                            PSNOWFALL    = minimum equivalent snow depth
!                                                           for snow falling during the
!                                                           current time step (m)

REAL, DIMENSION(SIZE(PPS)), INTENT(INOUT) :: PSNOWMAK

REAL, DIMENSION(:), INTENT(OUT)     :: PDMKXWBT
!                                            'M98' older Crocus computation for Martin and Lejeune 1998

!*      0.2    declarations of local variables
INTEGER JJ    ! loop counters

REAL, DIMENSION(SIZE(PTA))          :: ZTC, ZTW, ZEOD, ZTD, ZTAV, ZEOAV, DD,  &
                                       GA
!                                      ZTC	  = Atmospheric temp (°C)						                        p.spandre 2014/03/27
!                                      ZTW	  = Wet bulb temperature (K)						                    p.spandre 2014/03/27
!                                      ZEOD	  = Saturated vapor pressure at dew temp. (kPa)			        p.spandre 2014/06/04
!                                      ZTD	  = Dew Point temp.  (°C)						                        p.spandre 2014/06/04
!                                      ZTAV	  = Average temp. =(ZTD+ZTC)/2  (°C)					              p.spandre 2014/06/04
!                                      ZEOAV	= Saturated vapor pressure at average temp. ZTAV (kPa)		p.spandre 2014/06/04
!                                      DD	    = Slope of saturated vapor pressure curve (kPa/°C)        p.spandre 2014/06/04
!                                      GA	    = Psychrometric constant (kPa/°C)					                p.spandre 2014/06/04

LOGICAL, DIMENSION(SIZE(PTA))       :: LCONDSNOWMAK
LOGICAL, DIMENSION(SIZE(PTA))       :: LTIMESNOWMAK
LOGICAL                             :: PMONTH
LOGICAL                             :: PDAY
REAL, DIMENSION(31,31)              :: PRODTHEO
! 				                             LCONDSNOWMAK = Suitable Atmospheric conditions for snowmaking				      p.spandre 2014/03/28
!				                               LTIMESNOWMAK = Suitable timing conditions for snowmaking				            p.spandre 2014/03/28
! 				                             PMONTH	      = Suitable month for snowmaking = 1. Otherwise 0.				      p.spandre 2014/03/28
! 				                             PDAY	        = Suitable time in the day for snowmaking = 1. Otherwise 0.		p.spandre 2014/03/28
! 				                             PRODTHEO     = Theoretical production for each date (month, day)				    p.spandre 2014/03/28


!*       0.     Initialize variables:

ZTC(:)         = 0.0
ZTW(:)         = 0.0
ZEOD(:)        = 0.0
ZTD(:)         = 0.0
ZTAV(:)        = 0.0
ZEOAV(:)       = 0.0
DD(:)          = 0.0
GA(:)          = 0.0

LCONDSNOWMAK(:)= .FALSE.
LTIMESNOWMAK(:)= .FALSE.
PMONTH         = .FALSE.
PDAY           = .FALSE.
PRODTHEO(:,:)  = 0.0
PDMKXWBT(:)    = 0.0


!
! --- A) WBT computation and atmospheric conditions
!

DO JJ=1,SIZE(PTA)

! --- A1) Calculation of Wet Bulb Temperature according to Jensen

  ZTC(JJ) = PTA(JJ)-273.15		! Calculation of atmospheric temperature (°C)
                                 	! Calculation of dew point temp. TD (°C)
  IF (PQA(JJ) < 0.001) THEN		! Loop to prevent ZEOD from being negative or zero (bug with LOG calculation) 2014/09/04
    ZEOD(JJ)  = 0.001/(0.622+0.001)*PPS(JJ)/1000.                                   ! Vapor pressure at dew point (kPa) [2.9]
  ELSE
    ZEOD(JJ)  = (PQA(JJ)/PRHOA(JJ))/(0.622+0.378*(PQA(JJ)/PRHOA(JJ)))*PPS(JJ)/1000. ! Vapor pressure at dew point (kPa) [2.3]
  ENDIF             ! NB: 	PQA   = air humidity forcing (kg/m3)
	              ! PRHOA = air density 	=> mixing ratio r = m(vapor)/m(air) = m(vapor)/[Volume(air)*Density(air)] = [m(vapor)/Volume(air)]/Density(air) = PQA/PRHOA
  IF (ABS(LOG(ZEOD(JJ))-16.78) < 0.001) THEN      ! Loop to prevent LOG(ZEOD)-16.78 from being zero => divide by zero 2014/09/04
    ZTD(JJ) = (116.9+237.3*LOG(ZEOD(JJ)))/0.001
  ELSE
    ZTD(JJ) = (116.9+237.3*LOG(ZEOD(JJ)))/(16.78-LOG(ZEOD(JJ)))   ! Dew Point temperature (°C)  [7.11] and [7.22]
  ENDIF
								    ! Calculation of the slope of the saturation vapor pressure curve
  ZTAV(JJ) = (ZTD(JJ)+ZTC(JJ))/2                                  ! Average temperature between Dew point and actual conditions (cf. p176,7.19, Jensen)
  ZEOAV(JJ) = EXP((16.78*ZTAV(JJ)-116.9)/(ZTAV(JJ)+237.3))        ! Saturated vapor pressure at average temp. (kPa) [7.11]
  DD(JJ) = 4098.*ZEOAV(JJ)/(ZTAV(JJ)+237.3)**2                    ! Slope of the saturation vapor pressure curve (kPa/°C) [7.13]
								    ! Calculation of psychrometric constant
  GA(JJ) = PPS(JJ)/1000.*0.001013/(0.622*(2.501-2.361/1000.*ZTC(JJ)))    ! Latent heat of vaporization (MJ/kg) [7.1] included into GA formula [7.15]
									   ! Wet bulb temp. [7.19]
  ZTW(JJ) = (GA(JJ)*ZTC(JJ)+DD(JJ)*ZTD(JJ))/(GA(JJ)+DD(JJ))
  ZTW(JJ) = ZTW(JJ)+273.15            				   ! End of Wet Bulb Temperature Calculation

! --- A2) Boolean over atmospheric conditions for snowmaking

  IF (ZTW(JJ) < XPTA_SEUIL .and. PVMOD(JJ) < XWT) THEN
    LCONDSNOWMAK(JJ) = .TRUE.
  ELSE
    LCONDSNOWMAK(JJ) = .FALSE.
  ENDIF

ENDDO

!
! --- B) Timing conditions and production conditions
!
IF (IO%LSNOWMAK_BOOL) THEN

! --- B1) Timing conditions

! --- B1.1) Month conditions

  IF (TPTIME%TDATE%MONTH < 11. .and. TPTIME%TDATE%MONTH > 3.) THEN                                                  ! No production allowed from April to October included
    PMONTH = .FALSE.
  ELSE
    PMONTH = .TRUE.
  ENDIF

! --- B1.2) Day conditions

  IF (TPTIME%TDATE%MONTH*31+TPTIME%TDATE%DAY > XPP_D1 .and. TPTIME%TDATE%MONTH*31+TPTIME%TDATE%DAY < XPP_D2) THEN   ! Production allowed all day from 1st of NOV ... until 15th of DEC
    IF (TPTIME%TIME >= XPP_H1 .and. TPTIME%TIME <= XPP_H2) THEN 
      PDAY = .TRUE.
    ELSE
      PDAY = .FALSE.
    ENDIF
  ELSE
    IF (TPTIME%TIME > XPP_H4 .and. TPTIME%TIME < XPP_H3) THEN                                                       ! No production allowed between 8am and 7pm
      PDAY = .FALSE.
    ELSE
      PDAY = .TRUE.
    ENDIF
  ENDIF

! --- B2) Production conditions

  DO JJ=1,SIZE(PTA)

! --- B2.1) Self-production

    IF (IO%LSELF_PROD .or. (.NOT. IO%LSELF_PROD .and. XPROD_SCHEME(JJ)<0.)) THEN
    
      IF (TPTIME%TDATE%MONTH*31+TPTIME%TDATE%DAY > XPP_D1 .and. TPTIME%TDATE%MONTH*31+TPTIME%TDATE%DAY < XPP_D2) THEN   ! SM possible from 1st of NOV (11*31+1=342), until 15th of DEC (12*31+15=387)
        IF(PEK%TSNOW%MMP(JJ) <= XPT) THEN                                                                               ! Max admissible prod in that period
          XPSR_SNOWMAK = (XPR_A*(ZTW(JJ)-273.15)+XPR_B)*400./3600./3300.
            !WRITE(*,*) '-------------------'                                                          
            !WRITE(*,*) 'kg/m2 -> ',PEK%TSNOW%MMP(JJ)
            !WRITE(*,*) 'XPR_A -> ',XPR_A
            !WRITE(*,*) 'XPR_B -> ',XPR_B
            !WRITE(*,*) 'ZTW(JJ) -> ',ZTW(JJ)
            !WRITE(*,*) 'kg/m2/s -> ',XPSR_SNOWMAK
	    !WRITE(*,*) 's -> ',PEK%TSNOW%MMP(JJ)/XPSR_SNOWMAK
            !WRITE(*,*) '-------------------'
          IO%LPRODSNOWMAK(JJ) = .TRUE.      
        ELSE
          IO%LPRODSNOWMAK(JJ) = .FALSE.
        ENDIF
      ENDIF

      IF (TPTIME%TDATE%MONTH*31+TPTIME%TDATE%DAY > (XPP_D2-1.)) THEN       ! I.e. After December 15
        IF (PSNOW(JJ) < XPTR) THEN                                         ! If HTN < XPTR (m) keep producing
          IO%LPRODSNOWMAK(JJ) = .TRUE.
        ELSE
          IO%LPRODSNOWMAK(JJ) = .FALSE.
        ENDIF
      ENDIF

      IF (TPTIME%TDATE%MONTH*31+TPTIME%TDATE%DAY < XPP_D3) THEN            ! Case Between 1 JAN until 28th of Feb
        IF (PSNOW(JJ) < XPTR) THEN                                         ! If HTN < XPTR (m) keep producing
          IO%LPRODSNOWMAK(JJ) = .TRUE.   
        ELSE
          IO%LPRODSNOWMAK(JJ) = .FALSE.
        ENDIF
      ENDIF
     
! --- B2.2) No self-production

    ELSEIF (.NOT. IO%LSELF_PROD .and. XPROD_SCHEME(JJ)>=0.) THEN

      IF (PEK%TSNOW%MMP(JJ) <= XPROD_SCHEME(JJ)) THEN
        XPSR_SNOWMAK = (XPROD_SCHEME(JJ)-PEK%TSNOW%MMP(JJ))/2./3600.
        IO%LPRODSNOWMAK(JJ) = .TRUE.                                       ! PNPROD = integer : Suitable night for snowmaking : current prod < theo prod at 6pm
      ELSE                                                                 ! Then up to day+1 one can produce
        IO%LPRODSNOWMAK(JJ) = .FALSE.
      ENDIF

    ENDIF

! --- B3) Let's produce...

    IF (IO%LSELF_PROD .or. (.NOT. IO%LSELF_PROD .and. XPROD_SCHEME(JJ)<0.)) THEN

      IF (PDAY .and. PMONTH .and. IO%LPRODSNOWMAK(JJ)) THEN                ! Calendar (month+day timing) + suitable night => Timing conditions = TRUE, let's produce!
        LTIMESNOWMAK(JJ) = .TRUE.
      ELSE
        LTIMESNOWMAK(JJ) = .FALSE.
      ENDIF

    ELSEIF (.NOT. IO%LSELF_PROD .and. XPROD_SCHEME(JJ)>=0.) THEN

      IF (IO%LPRODSNOWMAK(JJ)) THEN 
        LTIMESNOWMAK(JJ) = .TRUE.
      ELSE
        LTIMESNOWMAK(JJ) = .FALSE.
      ENDIF

    ENDIF

  ENDDO

ENDIF
!
! --- C) Production
!
DO JJ=1,SIZE(PTA)

  IF (IO%LSNOWMAK_BOOL) THEN

    IF (IO%LSELF_PROD .or. (.NOT. IO%LSELF_PROD .and. XPROD_SCHEME(JJ)<0.)) THEN

      IF (LCONDSNOWMAK(JJ) .and. LTIMESNOWMAK(JJ)) THEN
        XPSR_SNOWMAK = (XPR_A*(ZTW(JJ)-273.15)+XPR_B)*400./3600./3300.
        PSNOWMAK(JJ)   = XPSR_SNOWMAK*PTSTEP/XRHO_SNOWMAK  
        PEK%TSNOW%MMP(JJ) = PEK%TSNOW%MMP(JJ)+PTSTEP*(XPSR_SNOWMAK*1/0.6)  ! 40% losses in converting water into snow
      ELSE
        PSNOWMAK(JJ)=0.
      ENDIF

    ELSEIF (.NOT. IO%LSELF_PROD .and. XPROD_SCHEME(JJ)>=0.) THEN

      IF (LTIMESNOWMAK(JJ)) THEN
        XPSR_SNOWMAK = (XPROD_SCHEME(JJ)-PEK%TSNOW%MMP(JJ))/2./3600.
        PSNOWMAK(JJ)   = XPSR_SNOWMAK*PTSTEP/XRHO_SNOWMAK  
        PEK%TSNOW%MMP(JJ) = PEK%TSNOW%MMP(JJ)+PTSTEP*XPSR_SNOWMAK    
      ELSE
        PSNOWMAK(JJ)=0.
      ENDIF

    ENDIF

    PDMKXWBT(JJ) = ZTW(JJ)-273.15

  ENDIF

  PSNOWFALL(JJ) = PSR(JJ)*PTSTEP/XRHOSMAX_ES + PSNOWMAK(JJ)                ! Minimum possible snowfall depth (m) + snowmaking depth by P.S 19/11/2013


ENDDO

END SUBROUTINE SNOW_MAKING

!! ======================================================================================================================
!!-----------------------	Snowmaking option by p.spandre	--------------------------------------------------------|
!!															|
!! A.Timing conditions for snowmaking
!!	A.1. Theoretical production
!!
!  IF (IO%LSNOWMAK_BOOL) THEN
!  !
!    DO JJ=1, 30
!      PRODTHEO(11,JJ) = XPROD_SCHEME(1)*30
!    ENDDO
!    DO JJ=1, 31
!      PRODTHEO(12,JJ) = XPROD_SCHEME(2)*31 + PRODTHEO(11,1)
!    ENDDO
!    DO JJ=1, 31
!      PRODTHEO(1,JJ) = XPROD_SCHEME(3)*31 + PRODTHEO(12,1)
!    ENDDO
!    DO JJ=1, 28
!      PRODTHEO(2,JJ) = XPROD_SCHEME(4)*28 + PRODTHEO(1,1)
!    ENDDO
!    DO JJ=1, 31
!      PRODTHEO(3,JJ) = XPROD_SCHEME(5)*31 + PRODTHEO(2,1)
!    ENDDO
!!
!!	A.2. Timing conditions
!!		A.2.1. Month condition
!    IF (TPTIME%TDATE%MONTH < 11. .and. TPTIME%TDATE%MONTH > 3.) THEN      ! No production allowed from april to otober included
!      PMONTH = .FALSE.
!    ELSE
!      PMONTH = .TRUE.
!    ENDIF
!  ! 	  	A.2.2. Daily condition
!    IF (TPTIME%TDATE%MONTH*31+TPTIME%TDATE%DAY > 341. .and. TPTIME%TDATE%MONTH*31+TPTIME%TDATE%DAY < 388.) THEN   ! Production allowed all day from 1st of NOV ... until 15th of DEC
!      PDAY = .TRUE.
!    ELSE
!      IF (TPTIME%TIME > 28800. .and. TPTIME%TIME < 64800.) THEN           ! No production allowed between 8am and 7pm
!        PDAY = .FALSE.
!      ELSE
!        PDAY = .TRUE.
!      ENDIF
!    ENDIF
!
!  !       A.3. Boolean from timing conditions
!    DO JJ=1,SIZE(PTA)
!!-----------------------	SELFPROD option by p.spandre	----------------------------------------|
!!				20150728								|
!!		A.2.4. SELFPROD option
!!
!      IF (IO%LSELF_PROD) THEN
!      !
!      !!!!!!!!!!!!!!!			FORMULATION ECHELLE ALPES			!!!!!!!!!!!!!!!!!!!
!      !
!        IF (TPTIME%TDATE%MONTH*31+TPTIME%TDATE%DAY > 341. .and. TPTIME%TDATE%MONTH*31+TPTIME%TDATE%DAY < 388.) THEN       ! i.e. SM possible from 1st of NOV (11*31+1=342) ... until 15th of DEC (12*31+15=387)
!        !																	! i.e. SM possible even during day time on that period
!          IF(1.0*XPROD_COUNT(JJ)*XPSR_SNOWMAK <= 150.) THEN     ! Max admissible prod in that period 150 kg/m2 
!            IO%LPRODSNOWMAK(JJ) = .TRUE.      ! .and. MOD(TPTIME%TDATE%DAY, 2) == 0.		REMOVED + Installation capacity 50% of snowguns simultaneously => 1 day / 2
!          ELSE
!            IO%LPRODSNOWMAK(JJ) = .FALSE.
!          ENDIF
!        ENDIF
!        !
!        IF (TPTIME%TDATE%MONTH*31+TPTIME%TDATE%DAY > 387.) THEN                         ! i.e. After December 15
!          IF (ZSNOW(JJ) < 0.60) THEN                                                    ! If HTN < 0.6 (m) keep producing
!            IO%LPRODSNOWMAK(JJ) = .TRUE.
!          ELSE
!            IO%LPRODSNOWMAK(JJ) = .FALSE.
!          ENDIF
!        ENDIF
!
!        IF (TPTIME%TDATE%MONTH*31+TPTIME%TDATE%DAY < 92.) THEN      ! Case Between 1 JAN until 28th of Feb.!
!
!          IF (ZSNOW(JJ) < 0.60) THEN      ! If HTN < 0.6 (m) keep producing
!            IO%LPRODSNOWMAK(JJ) = .TRUE.      !  .and. MOD(TPTIME%TDATE%DAY, 2) == 0.		REMOVED + Installation capacity 50% of snowguns simultaneously => 1 day / 2
!          ELSE
!            IO%LPRODSNOWMAK(JJ) = .FALSE.
!          ENDIF
!        ENDIF
!      !!!!!!!!!!!!!!!			FORMULATION ECHELLE ALPES			!!!!!!!!!!!!!!!!!!!
!      ! 
!      ELSE      ! SELF_PROD conditions is FALSE
!        !A.2.3. Suitable night for snowmaking
!        IF (TPTIME%TIME == 64800.) THEN     ! condition at 6pm i.e. for each time step, you compare the total.
!          IF (XPROD_COUNT(JJ) < PRODTHEO(TPTIME%TDATE%MONTH,TPTIME%TDATE%DAY)) THEN
!            IO%LPRODSNOWMAK(JJ) = .TRUE.      ! PNPROD = integer : Suitable night for snowmaking : current prod < theo prod at 6pm		p.spandre 2014/03/28
!          ELSE      ! then up to day+1 one can produce
!            IO%LPRODSNOWMAK(JJ) = .FALSE.
!          ENDIF
!        ENDIF
!      ENDIF
!!													|
!!-----------------------	SELFPROD option by p.spandre	----------------------------------------|
!
!      IF (PDAY .and. PMONTH .and. IO%LPRODSNOWMAK(JJ)) THEN     ! Calendar (month+day timing) + suitable night => Timing conditions = TRUE, let's produce!
!        LTIMESNOWMAK(JJ) = .TRUE.
!      ELSE
!        LTIMESNOWMAK(JJ) = .FALSE.
!      ENDIF
!
!    ENDDO
!  ENDIF
!
!  DO JJ=1,SIZE(PTA)
!    IF (IO%LSNOWMAK_BOOL) THEN
!!
!! B. Atmospheric conditions for snowmaking
!!	B.1. Calculation of Wet Bulb temperature calculation according to Jensen,ASCE, 1990 (added p.spandre 04/06/2014)
!      ZTC(JJ) = PTA(JJ)-273.15       !calculation of atmospheric temperature (°C)
!!		B.1.1 Calculation of dew point temp. TD (°C)
!      IF (PQA(JJ) < 0.001) THEN     ! loop to prevent ZEOD from being negative or zero (bug with LOG calculation) 2014/09/04
!        ZEOD(JJ)  = 0.001/(0.622+0.001)*PPS(JJ)/1000. ! Vapor pressure at dew point (kPa) [2.9]
!      ELSE
!        ZEOD(JJ)  = (PQA(JJ)/PRHOA(JJ))/(0.622+0.378*(PQA(JJ)/PRHOA(JJ)))*PPS(JJ)/1000. ! Vapor pressure at dew point (kPa) [2.3]
!      ENDIF             ! NB: 	PQA   = air humidity forcing (kg/m3)
!!												! 	PRHOA = air density 	=> mixing ratio r = m(vapor)/m(air) = m(vapor)/[Volume(air)*Density(air)] = [m(vapor)/Volume(air)]/Density(air) = PQA/PRHOA
!      IF (ABS(LOG(ZEOD(JJ))-16.78) < 0.001) THEN      ! loop to prevent LOG(ZEOD)-16.78 from being zero => divide by zero 2014/09/04
!        ZTD(JJ) = (116.9+237.3*LOG(ZEOD(JJ)))/0.001
!      ELSE
!        ZTD(JJ) = (116.9+237.3*LOG(ZEOD(JJ)))/(16.78-LOG(ZEOD(JJ)))   ! Dew Point temperature (°C)  [7.11] and [7.22]
!      ENDIF
!!		B.1.2. Calculation of the slope of the saturation vapor pressure curve
!      ZTAV(JJ) = (ZTD(JJ)+ZTC(JJ))/2         ! Average temperature between Dew point and actual conditions (cf. p176,7.19, Jensen)
!      ZEOAV(JJ) = EXP((16.78*ZTAV(JJ)-116.9)/(ZTAV(JJ)+237.3))      ! Saturated vapor pressure at average temp. (kPa) [7.11]
!      DD(JJ) = 4098.*ZEOAV(JJ)/(ZTAV(JJ)+237.3)**2   ! Slope of the saturation vapor pressure curve (kPa/°C) [7.13]
!!		B.1.3. Calculation of psychrometric constant
!      GA(JJ) = PPS(JJ)/1000.*0.001013/(0.622*(2.501-2.361/1000.*ZTC(JJ)))    ! Latent heat of vaporization (MJ/kg) [7.1] included into GA formula [7.15]
!!		B.1.4. Wet bulb temp. [7.19]
!      ZTW(JJ) = (GA(JJ)*ZTC(JJ)+DD(JJ)*ZTD(JJ))/(GA(JJ)+DD(JJ))
!      ZTW(JJ) = ZTW(JJ)+273.15
!!   		End of Wet Bulb Temperature Calculation
!!
!!	B.2. Boolean over atmospheric conditions for snowmaking
!      IF (ZTW(JJ) < XPTA_SEUIL .and. PVMOD(JJ) < 4.2) THEN
!        LCONDSNOWMAK(JJ) = .TRUE.
!      ELSE
!        LCONDSNOWMAK(JJ) = .FALSE.
!      ENDIF
!!
!! C. Boolean over timing + atmospheric conditions									! Production possible even if natural snow falling P.Spandre 2014/03/04
!      IF (LCONDSNOWMAK(JJ) .and. LTIMESNOWMAK(JJ)) THEN
!        ZSNOWMAK(JJ)   = XPSR_SNOWMAK*PTSTEP/XRHO_SNOWMAK  
!        ! snowmaking depth by P.S 19/11/2013
!        XPROD_COUNT(JJ) = XPROD_COUNT(JJ)+PTSTEP
!      ELSE
!        ZSNOWMAK(JJ)=0.
!      ENDIF
!    ENDIF
!    ZSNOWFALL(JJ)      = PSR(JJ)*PTSTEP/XRHOSMAX_ES + ZSNOWMAK(JJ)      ! MINImum possible snowfall depth (m) + snowmaking depth by P.S 19/11/2013
!!
!    IF (IO%LSNOWMAK_BOOL) DMK%XPRODCOUNT(JJ) = XPROD_COUNT(JJ)
!!
!  ENDDO
!!															|
!!-----------------------	Snowmaking option by p.spandre	--------------------------------------------------------|
