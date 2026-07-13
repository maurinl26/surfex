!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!     ######spl
! ###########################################################################
SUBROUTINE OL_TIME_INTERP_ATM (KSURF_STEP,KNB_ATM,PPDISTRIB,                 &
                               PTA,PQA,PWIND,PDIR_SW,PSCA_SW,PLW,            &
                               PSNOW,PRAIN,PPS,PCO2,PDIR,PO3,PAE,            &
                               PIMPWET,PIMPDRY,PSUMZEN,PZEN1,PZEN_INT,PZEN2  )  
! ###########################################################################
!**************************************************************************
!
!!    PURPOSE
!!    -------
!        Time interpolation of the atmospheric forcing
!        So far, it is a simple linear interpolation.
!        More complex interpolation may be added, especially for the atmospheric
!        radiation (Option to use).
!        Output are in the module
!!
!!**  METHOD
!!    ------
!!
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
!!      F. Habets   *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    06/2003
!       B. Decharme 06/2021 new interpolation method (LNEW_TIME_INTERP_ATM) to conserv all fluxes considering :
!                                      (1) current forcing value (if flux : average from t-1 to t)
!                                      (2) next forcing value (if flux : average from t to t+1)
!                                      (3) ueber-next forcing value (if flux : average from t+1 to t+2)
!
! ###########################################################################
!
USE MODN_IO_OFFLINE, ONLY : CINTERP_SW, LLIMIT_QAIR, LNEW_TIME_INTERP_ATM, CTIME_INTERP_PRCP, &
                            LFORCIMP, LFORCATMOTARTES
!
! ###########################################################################
USE MODD_CSTS,       ONLY : XPI, XRD, XRV, XG, XI0, XSURF_EPSILON
! ###########################################################################
USE MODD_SURF_PAR,   ONLY : XUNDEF
USE MODD_FORC_ATM,   ONLY : XTA       ,&! air temperature forcing               (K)
                            XQA       ,&! air specific humidity forcing         (kg/m3)
                            XRHOA     ,&! air density forcing                   (kg/m3)
                            XZS       ,&! orography                             (m)
                            XU        ,&! zonal wind                            (m/s)
                            XV        ,&! meridian wind                         (m/s)
                            XDIR_SW   ,&! direct  solar radiation (on horizontal surf.)
                            XSCA_SW   ,&! diffuse solar radiation (on horizontal surf.)
                            XLW       ,&! longwave radiation (on horizontal surf.)
                            XPS       ,&! pressure at atmospheric model surface (Pa)
                            XPA       ,&! pressure at forcing level             (Pa)
                            XRHOA     ,&! density at forcing level              (kg/m3)
                            XCO2      ,&! CO2 concentration in the air          (kg/kg)
                            XO3       ,&! Ozone
                            XAE       ,&! Aerosol optical depth
                            XIMPWET   ,&! wet deposit coefficient
                            XIMPDRY   ,&! dry deposit coefficient
                            XSNOW     ,&! snow precipitation                    (kg/m2/s)
                            XRAIN     ,&! liquid precipitation                  (kg/m2/s)
                            XZREF       ! height of T,q forcing                 (m)  
!
USE MODI_GET_LUOUT
USE MODI_ABOR1_SFX
!
USE MODE_THERMOS
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
!
#ifndef AIX64
!$ INCLUDE 'omp_lib.h'
#endif
!
! global variables
! ###########################################################################
!
INTEGER,               INTENT(IN)    :: KSURF_STEP, KNB_ATM
!  
REAL, DIMENSION(:),    INTENT(IN)    :: PPDISTRIB
!
REAL, DIMENSION(:,:),  INTENT(IN)    :: PTA,PQA,PWIND
REAL, DIMENSION(:,:),  INTENT(INOUT) :: PDIR_SW,PSCA_SW,PLW
REAL, DIMENSION(:,:),  INTENT(IN)    :: PSNOW,PRAIN,PPS,PCO2,PDIR,PO3,PAE
REAL, DIMENSION(:,:,:),INTENT(IN)    :: PIMPWET,PIMPDRY
REAL, DIMENSION(:),    INTENT(IN)    :: PSUMZEN,PZEN_INT,PZEN1,PZEN2
!
! ###########################################################################
!
! local variables
!
CHARACTER(LEN=1) :: YFLAG ! Type of interpolation
!
REAL, DIMENSION(SIZE(PWIND,1),SIZE(PWIND,2)) :: ZU, ZV
!
REAL :: ZPI, ZCOEF, ZQSAT
!
! ###########################################################################
REAL :: ZINT_SW
REAL :: ZSW1
REAL :: ZSW2
REAL :: ZF_SCA_SW1
REAL :: ZF_SCA_SW2  
REAL :: ZTHEO_SW1
REAL :: ZTHEO_SW2
REAL :: ZF_THEO_SW1
REAL :: ZF_THEO_SW2
REAL :: ZDF_THEO_SW
REAL :: ZDF_SCA_SW
REAL :: ZTHEO_SW
! ###########################################################################
INTEGER :: J,JIMP
INTEGER :: ILUOUT
REAL(KIND=JPRB) :: ZHOOK_HANDLE, ZHOOK_HANDLE_OMP
!========================================================================
!
IF (LHOOK) CALL DR_HOOK('OL_TIME_INTERP_ATM_1',0,ZHOOK_HANDLE)
!
CALL GET_LUOUT('OFFLIN',ILUOUT)
!
ZPI = XPI/180.
!
!------------------------------------------------------
!Compute variation from atmospheric time step J and J+1
!------------------------------------------------------
!
!------------------------------------------------------
! Instantaneous fields (Current value valid at midpoint of interval)
!------------------------------------------------------
!
!Current value valid at midpoint of interval (linear interp as before)
YFLAG = 'I'
!
!Air temperature (K)
CALL INTERPOL_FRC(YFLAG,PTA(:,:),XTA(:))
!
!Air humidity
CALL INTERPOL_FRC(YFLAG,PQA(:,:),XQA(:))
!
!Surface Pressure
CALL INTERPOL_FRC(YFLAG,PPS(:,:),XPS(:))
!
!CO2 Concentration
CALL INTERPOL_FRC(YFLAG,PCO2(:,:),XCO2(:))
!
IF (LFORCATMOTARTES) THEN
!
!  O3
   CALL INTERPOL_FRC(YFLAG,PO3(:,:),XO3(:))
!
!  AE
   CALL INTERPOL_FRC(YFLAG,PAE(:,:),XAE(:))
!
ENDIF
!
IF (LFORCIMP) THEN
   XIMPWET(:,:) = PIMPWET(:,:,2)
   XIMPDRY(:,:) = PIMPDRY(:,:,2)
ENDIF
!
!----------------
! Averaged fluxes
!----------------
!
IF(LNEW_TIME_INTERP_ATM)THEN
! current value is average for period ending at current time
  YFLAG = 'F'
ELSE
! Current value valid at midpoint of interval (linear interp as before)
  YFLAG = 'I'
ENDIF
!
!Zonal wind
WHERE(PWIND(:,:)/=XUNDEF)
     ZU(:,:) = PWIND(:,:) * SIN(PDIR(:,:)*ZPI)
     ZV(:,:) = PWIND(:,:) * COS(PDIR(:,:)*ZPI)
ENDWHERE
CALL INTERPOL_FRC(YFLAG,ZU(:,:),XU(:))
CALL INTERPOL_FRC(YFLAG,ZV(:,:),XV(:))
!
!Longwave radiation
CALL INTERPOL_FRC(YFLAG,PLW(:,:),XLW(:))
!
!------------------------------------------------------
!Shortwave radiation (direct & diffuse)
!------------------------------------------------------
!
IF(LNEW_TIME_INTERP_ATM.OR.CINTERP_SW=='LIN')THEN
  !
  CALL INTERPOL_FRC(YFLAG,PDIR_SW(:,:),XDIR_SW(:,1))
  CALL INTERPOL_FRC(YFLAG,PSCA_SW(:,:),XSCA_SW(:,1))
  !
  IF (LHOOK) CALL DR_HOOK('OL_TIME_INTERP_ATM_1',1,ZHOOK_HANDLE)
  !
ELSE
  !
  IF (LHOOK) CALL DR_HOOK('OL_TIME_INTERP_ATM_1',1,ZHOOK_HANDLE)
  !
  !$OMP PARALLEL PRIVATE(ZHOOK_HANDLE_OMP)
  IF (LHOOK) CALL DR_HOOK('OL_TIME_INTERP_ATM_2',0,ZHOOK_HANDLE_OMP)
  !$OMP DO PRIVATE(J,ZCOEF,&
  !$OMP ZSW1,ZSW2,ZF_SCA_SW1,ZF_SCA_SW2,ZTHEO_SW1,ZTHEO_SW2,&
  !$OMP ZDF_SCA_SW,ZDF_THEO_SW,ZTHEO_SW,ZINT_SW)
  DO J = 1,SIZE(PTA,1)
    !
    IF (PTA(J,1)/=XUNDEF) THEN
       !
      SELECT CASE (CINTERP_SW)
      CASE ('ZEN')
        !
        ! new method of interpolation
        ! ---------------------------
        !
        ZCOEF=REAL(KSURF_STEP-1)/REAL(KNB_ATM)
        !
        ! Negative radiation is set to 0.0
        !
        ZSW1 = PDIR_SW(J,1) + PSCA_SW(J,1)
        !
        IF (ZSW1.LT.-XSURF_EPSILON) THEN
           PDIR_SW(J,1) = 0.0
           PSCA_SW(J,1) = 0.0
        ENDIF
        !
        PDIR_SW(J,1) = MAX(0.0,PDIR_SW(J,1))
        PSCA_SW(J,1) = MAX(0.0,PSCA_SW(J,1))
        !
        ZSW2 = PDIR_SW(J,2) + PSCA_SW(J,2)
        !
        IF (ZSW2.LT.-XSURF_EPSILON) THEN
           PDIR_SW(J,2) = 0.0
           PSCA_SW(J,2) = 0.0
        ENDIF   
        !
        PDIR_SW(J,2) = MAX(0.0,PDIR_SW(J,2))
        PSCA_SW(J,2) = MAX(0.0,PSCA_SW(J,2))
        !
        ! Calculation of total radiation (SW1,2) and the scattered fraction (F_SCA_SW1,2)
        !
        !
        ZSW1 = PDIR_SW(J,1) + PSCA_SW(J,1)
        ZSW2 = PDIR_SW(J,2) + PSCA_SW(J,2)
        !
        !
        IF (ZSW1 .LT. XSURF_EPSILON) THEN
           ZF_SCA_SW1 = 1.
        ELSE
           ZF_SCA_SW1 = PSCA_SW(J,1) / ZSW1
        ENDIF
        !
        IF (ZSW2 .LT. XSURF_EPSILON) THEN
           ZF_SCA_SW2 = 1.
        ELSE
           ZF_SCA_SW2 = PSCA_SW(J,2) / ZSW2
        ENDIF
        !
        !
        IF (ZF_SCA_SW1.LT.-XSURF_EPSILON) CALL ABOR1_SFX("OL_TIME_INTEPOL: Wrong fraction")
        IF (ZF_SCA_SW2.LT.-XSURF_EPSILON) CALL ABOR1_SFX("OL_TIME_INTEPOL: Wrong fraction")
        !
        IF (ZF_SCA_SW1.GT.1.0+XSURF_EPSILON) CALL ABOR1_SFX("OL_TIME_INTEPOL: Wrong fraction")
        IF (ZF_SCA_SW2.GT.1.0+XSURF_EPSILON) CALL ABOR1_SFX("OL_TIME_INTEPOL: Wrong fraction")
        !
        ! Calculation of theoretical radiation (THEO_SW1,2) and ratio between total and theoretical radiation (F_THEO_SW1,2)
        !
        ZTHEO_SW1 = MAX(XI0*COS(PZEN1(J)),0.)
        ZTHEO_SW2 = MAX(XI0*COS(PZEN2(J)),0.)
        !
        !
        ZF_THEO_SW1 = MIN(ZSW1/MAX(1.0E-6,ZTHEO_SW1),1.)
        ZF_THEO_SW2 = MIN(ZSW2/MAX(1.0E-6,ZTHEO_SW2),1.)
        !
        IF ( (ZTHEO_SW1.LT.XSURF_EPSILON) .AND. (ZTHEO_SW2.GT.XSURF_EPSILON) ) THEN
           ZF_THEO_SW1 = ZF_THEO_SW2
        ELSEIF ( (ZTHEO_SW2.LT.XSURF_EPSILON) .AND. (ZTHEO_SW1.GT.XSURF_EPSILON) ) THEN
           ZF_THEO_SW2 = ZF_THEO_SW1  
        ENDIF
        !
        !
        IF (ZF_THEO_SW1.LT.-XSURF_EPSILON) CALL ABOR1_SFX("OL_TIME_INTEPOL: Wrong fraction")
        IF (ZF_THEO_SW2.LT.-XSURF_EPSILON) CALL ABOR1_SFX("OL_TIME_INTEPOL: Wrong fraction")
        !
        IF (ZF_THEO_SW1.GT.1.0+XSURF_EPSILON) CALL ABOR1_SFX("OL_TIME_INTEPOL: Wrong fraction")
        IF (ZF_THEO_SW2.GT.1.0+XSURF_EPSILON) CALL ABOR1_SFX("OL_TIME_INTEPOL: Wrong fraction")
        !
        ! Linear temporal interpolation
        !
        ZDF_SCA_SW  = (1.0 - ZCOEF) * ZF_SCA_SW1  + ZCOEF * ZF_SCA_SW2
        ZDF_THEO_SW = (1.0 - ZCOEF) * ZF_THEO_SW1 + ZCOEF * ZF_THEO_SW2
        !
        !
        ! Calculatin of theoretical radiation based on intermediate zenith angle
        !
        ZTHEO_SW = MAX(XI0*COS(PZEN_INT(J)), 0.)
        !
        !
        ZINT_SW = ZDF_THEO_SW * ZTHEO_SW
        !
        !
        XSCA_SW(J,1) = ZDF_SCA_SW * ZINT_SW
        XDIR_SW(J,1) = (1 - ZDF_SCA_SW) * ZINT_SW
        !
      CASE ('OLD')
        !
        ZCOEF=0.
        IF (PSUMZEN(J)>0.) ZCOEF=MAX((COS(PZEN_INT(J))/PSUMZEN(J)),0.)
        !
        XDIR_SW(J,1) = MIN(PDIR_SW(J,2)*ZCOEF,1300.0*MAX(COS(PZEN_INT(J)),0.))
        !
        XSCA_SW(J,1) = MIN(PSCA_SW(J,2)*ZCOEF,1300.0*MAX(COS(PZEN_INT(J)),0.))
        !
      CASE DEFAULT
        CALL ABOR1_SFX('OL_TIME_INTEPOL: OPTION "'//CINTERP_SW//'" NOT RECOGNIZED FOR CINTERP_SW in NAM_IO_OFFLINE')
      END SELECT
      !
    ENDIF
    !
  ENDDO
  !
  !$OMP END DO
  !
  IF (LHOOK) CALL DR_HOOK('OL_TIME_INTERP_ATM_2',1,ZHOOK_HANDLE_OMP)
  !$OMP END PARALLEL
  !
ENDIF
!
IF (LHOOK) CALL DR_HOOK('OL_TIME_INTERP_ATM_3',0,ZHOOK_HANDLE)
!            
!------------------------------------------------------
! Precipitation
!------------------------------------------------------
!
IF(LNEW_TIME_INTERP_ATM)THEN
  SELECT CASE (CTIME_INTERP_PRCP)
    CASE ('OLD')
!     Current value is applied centered on current time without interpolation  
      YFLAG = 'C'
    CASE ('PDF')
!     Current value is from a PDF that disagregates in time over the forcing interval
      YFLAG = 'P'
    CASE DEFAULT
!     current value is average for period ending at current time
      YFLAG = 'F'
  END SELECT
ELSE
  YFLAG = 'C'
ENDIF
!
CALL INTERPOL_FRC(YFLAG,PRAIN(:,:),XRAIN(:))
CALL INTERPOL_FRC(YFLAG,PSNOW(:,:),XSNOW(:))
!
!
! Check No value data
!---------------------
! Error cases
!
IF ((MINVAL(XTA)  .EQ.XUNDEF).OR.(MINVAL(XQA).EQ.XUNDEF).OR.&
      (MINVAL(XU).EQ.XUNDEF).OR.(MINVAL(XRAIN).EQ.XUNDEF).OR.&
      (MINVAL(XSNOW).EQ.XUNDEF)) THEN  
    WRITE(ILUOUT,*)'MINVAL(XTA),MINVAL(XQA),MINVAL(XU),MINVAL(XRAIN),MINVAL(XSNOW)'
    WRITE(ILUOUT,*)MINVAL(XTA),MINVAL(XQA),MINVAL(XU),MINVAL(XRAIN),MINVAL(XSNOW)
    CALL ABOR1_SFX('OL_TIME_INTERP_ATM: UNDEFINED VALUE IN ATMOSPHERIC FORCING')
ENDIF
!
IF ((MINVAL(XDIR_SW).EQ.XUNDEF).AND.(MINVAL(XSCA_SW).EQ.XUNDEF)) THEN
    WRITE(ILUOUT,*)'MINVAL(XSCA_SW),MINVAL(XDIR_SW)'
    WRITE(ILUOUT,*)MINVAL(XSCA_SW),MINVAL(XDIR_SW)
    CALL ABOR1_SFX('OL_TIME_INTERP_ATM: UNDEFINED VALUE IN ATMOSPHERIC FORCING')
ENDIF
!
IF ((MINVAL(XPS).EQ.XUNDEF).AND.(MINVAL(XZS).EQ.XUNDEF)) THEN
    WRITE(ILUOUT,*)'MINVAL(XPS),MINVAL(XZS)'
    WRITE(ILUOUT,*)MINVAL(XPS),MINVAL(XZS)
    CALL ABOR1_SFX('OL_TIME_INTERP_ATM: UNDEFINED VALUE IN ATMOSPHERIC FORCING')
ENDIF
!
IF (MINVAL(XDIR_SW).EQ.XUNDEF) XDIR_SW(:,:)=0. ! No direct solar radiation
IF (MINVAL(XSCA_SW).EQ.XUNDEF) XSCA_SW(:,:)=0. ! No diffuse solar radiation
IF (MINVAL(XPS)    .EQ.XUNDEF) THEN            ! No surface Pressure 
   WRITE(ILUOUT,*)' OL_TIME_INTERP_ATM: SURFACE PRESSURE COMPUTED FROM ZS'
   XPS(:)  = 101325*(1-0.0065 * XZS(:)/288.15)**5.31
ENDIF
!
!* Forcing level pressure from hydrostatism
!
XPA(:) = XUNDEF
WHERE(XPS(:)/=XUNDEF)
  XPA(:) = XPS(:) - XRHOA(:) * XZREF(:) * XG
ENDWHERE
!
! Limit humidity
!
IF(LLIMIT_QAIR)THEN
  DO J = 1,SIZE(XPS)
    IF(XPS(J)/=XUNDEF)THEN
      ZQSAT  = QSAT(XTA(J),XPS(J))
      XQA(J) = MIN(XQA(J),ZQSAT)
    ENDIF
  ENDDO
ENDIF
!
! Air density and humidity in kg/m3
!
WHERE(XPS(:)/=XUNDEF)
  XRHOA(:) = XPS(:) / ( XTA(:)*XRD * ( 1.+((XRV/XRD)-1.)*XQA(:) ) + XZREF(:)*XG )
  XQA  (:) = XQA(:) * XRHOA(:)
ENDWHERE
!
IF (LHOOK) CALL DR_HOOK('OL_TIME_INTERP_ATM_3',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
!
CONTAINS
!
!-------------------------------------------------------------------------------
!
SUBROUTINE INTERPOL_FRC(HFLAG,PVALIN,PVALOUT)
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
CHARACTER(LEN=1),INTENT(IN) :: HFLAG ! Type of interpolation

!
REAL, DIMENSION(:,:), INTENT(IN) :: PVALIN ! 1 : current forcing value (if flux : average from t-1 to t)
                                           ! 2 : next forcing value (if flux : average from t to t+1)
                                           ! 3 : ueber-next forcing value (if flux : average from t+1 to t+2)
!
REAL, DIMENSION(:), INTENT(OUT)   :: PVALOUT  ! Interpolated forcing data vector
!
!*      0.2    declarations of local variables
!
REAL, DIMENSION(SIZE(PVALIN,1)) :: ZDENOM  ! ZDENOMinator of scaling factor for 
                                           ! conserving interpolation
REAL, DIMENSION(SIZE(PVALIN,1)) :: ZNUMER  ! ZNUMERator of scaling factor for 
                                           ! conserving interpolation
!                                            
REAL            :: ZFAC1, ZFAC2, ZFAC3
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
IF (LHOOK) CALL DR_HOOK('OL_TIME_INTERP_ATM:INTERPOL_FRC',0,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
!
PVALOUT(:) = XUNDEF
!
ZFAC1 = 0.0
ZFAC2 = 0.0
ZFAC3 = 0.0
!
IF (HFLAG == 'I') THEN
!        
!** Current value valid at midpoint of interval (method as before)
!
    ZFAC2 = REAL(KSURF_STEP-1)/REAL(KNB_ATM)
    ZFAC1 = 1.0-ZFAC2
!    
    WHERE(PVALIN(:,1)/=XUNDEF.AND.PVALIN(:,2)/=XUNDEF)
         PVALOUT(:) = PVALIN(:,1)*ZFAC1+PVALIN(:,2)*ZFAC2
    ENDWHERE
!
ELSEIF (HFLAG == 'F') THEN
!    
!** Current value is average for period ending at current time
!
    ZFAC2 = (2.0*KNB_ATM-ABS(REAL(2*KSURF_STEP-KNB_ATM-1)))/(KNB_ATM*2.0)
    ZFAC1 = MAX(1.0-REAL(KSURF_STEP*2+KNB_ATM-1)/(KNB_ATM*2.0),0.0)
    ZFAC3 = MAX(1.0-REAL((KNB_ATM+1-KSURF_STEP)*2+KNB_ATM-1)/(KNB_ATM*2.0),0.0)
!
    ZDENOM(:) = 0.0
    ZNUMER(:) = XUNDEF
    WHERE(PVALIN(:,1)/=XUNDEF.AND.PVALIN(:,2)/=XUNDEF.AND.PVALIN(:,3)/=XUNDEF)
         ZDENOM(:) = 0.5*(PVALIN(:,1)+PVALIN(:,3))+3.0*PVALIN(:,2)
         ZNUMER(:) = 4.0*PVALIN(:,2)
    ENDWHERE
!
    WHERE(ABS(ZDENOM(:)) > EPSILON(ZFAC1))
         PVALOUT(:) = (PVALIN(:,1)*ZFAC1+PVALIN(:,2)*ZFAC2+PVALIN(:,3)*ZFAC3) * ZNUMER(:) / ZDENOM(:)
    ELSEWHERE
         PVALOUT(:) = 0.0
    ENDWHERE
!
ELSEIF (HFLAG == 'C') THEN
!    
!** Current value is applied centered on current time without interpolation
!
    WHERE(PVALIN(:,2)/=XUNDEF)
         PVALOUT(:) = PVALIN(:,2)
    ENDWHERE
!
ELSEIF (HFLAG == 'P') THEN
!
!** Current value is from a PDF that disagregates in time over the forcing interval
!
    PVALOUT(:) = MAX((PVALIN(:,2)*PPDISTRIB(KSURF_STEP)*REAL(KNB_ATM)),0.0)
!
ENDIF
!
!
IF (LHOOK) CALL DR_HOOK('OL_TIME_INTERP_ATM:INTERPOL_FRC',1,ZHOOK_HANDLE)
!
END SUBROUTINE INTERPOL_FRC
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE OL_TIME_INTERP_ATM
