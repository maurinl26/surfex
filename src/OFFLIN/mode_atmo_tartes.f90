!# -*- coding: utf-8 -*-
!!!! M. Dumont 
!!!! v1 21/01/2015

MODULE MODE_ATMO_TARTES
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
CONTAINS
!###################################################################
SUBROUTINE RADIANCE(KNI,KYEAR, KMONTH, KDAY, &
                    KBANDS,OFORCATMOTARTES,PSCA_SW,PDIR_SW)
!
! compute direct and diffuse solar radiance according to irradiance for a given atmospheric profile
!
USE MODD_CONST_ATM, ONLY : JPNBANDS_ATM, JPNLYR_CLEAR, PPZP_CUT, PPHUND, PPMU_THRESHOLD
USE MODI_JULIAN
USE MODD_SNOW_METAMO,  ONLY : XUEPSI
USE MODD_FORC_ATM,  ONLY: &
                          XO3         ,&! Ozone
                          XAE         ,&! Aerosol optical depth
                          XZENITH     ,&! zenithal angle (radian from the vertical)
                          XRHOA       ,&! density at forcing level              (kg/m3)
                          XQA         ,&! air humidity forcing                  (kg/m3)
                          XPS         ,&! pressure at atmospheric model surface (Pa)
                          XTA         ! air temperature forcing               (K)
!
IMPLICIT NONE
!
!REAL, DIMENSION(:), INTENT(IN)::
INTEGER,INTENT(IN)               ::KNI, KYEAR, KMONTH, KDAY, KBANDS
LOGICAL,INTENT(IN)               ::OFORCATMOTARTES
REAL, DIMENSION(:,:), INTENT(INOUT):: PSCA_SW, PDIR_SW
!
! For spectral repartion of direct/diffuse solar irradiance 
REAL, DIMENSION(KNI)     :: ZP_CLOUD, ZTCLOUD55, ZMU, ZD_O3,ZD_AE
REAL, DIMENSION(KNI)     ::  ZINT_SCA_SW, ZINT_DIR_SW, ZINT_TOT_SW
INTEGER, DIMENSION(KNI)  :: KCLOUD_TYPE
!
REAL, DIMENSION(KNI,KBANDS)        :: ZIRR_DIFF, ZIRR_DIR
REAL, DIMENSION(KNI, JPNLYR_CLEAR) :: ZP_CUT
REAL :: ZDATI
INTEGER :: JI
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
IF (LHOOK) CALL DR_HOOK('RADIANCE',0,ZHOOK_HANDLE)
! 
! broadband direct and diffuse irradiance
ZINT_SCA_SW(:)=PSCA_SW(:,1)
ZINT_DIR_SW(:)=PDIR_SW(:,1)
ZINT_TOT_SW(:)=ZINT_DIR_SW(:)+ZINT_SCA_SW(:)
!
PSCA_SW(:,2:KBANDS)=0.
PDIR_SW(:,2:KBANDS)=0.
!
IF (OFORCATMOTARTES) THEN 
  ZD_O3(:)= XO3(:)/100 ! integrated ozone (atm-cm)
  ZD_AE(:)= XAE(:) ! aerosol optical depth
ELSE 
  ZD_O3(:)= 0.3 ! integrated ozone (atm-cm)
  ZD_AE(:)= 0.1 ! aerosol optical depth
ENDIF
!
ZMU(:)=0.5!COS(XZENITH(:))
ZP_CLOUD(:)= 378.*100.! cloud bottom pressure (Pa)
KCLOUD_TYPE(:)= 1! cloud type
!
IF (MAXVAL(ZINT_TOT_SW(:))>XUEPSI) THEN
  !IF (MINVAL(XZENITH(:))>0.) THEN
  IF (MINVAL(ZMU(:))>XUEPSI) THEN
    ! prevent calculation during night
    ! Julian date       
    CALL JULIAN(KYEAR, KMONTH, KDAY, 0., ZDATI)
    !
    IF (MAXVAL(ZINT_TOT_SW(:)-ZINT_DIR_SW(:)).GT.XUEPSI) THEN
      ! case where diffuse and direct radiation are known 
      ! calculation of cloud optical depth
      CALL TAU_CLOUD(ZMU(:),ZINT_SCA_SW(:)/ZINT_TOT_SW(:),ZTCLOUD55(:))
      !
      CALL IRRADIANCE(ZDATI,ZMU(:),XQA(:)/XRHOA(:),ZD_O3(:),ZD_AE(:),   &
                      XPS(:),XTA(:),PPZP_CUT,ZP_CLOUD(:),KCLOUD_TYPE(:),&
                      ZTCLOUD55(:),ZIRR_DIR(:,:),ZIRR_DIFF(:,:))
      !      
      ! normalisation by broadband direct and diffuse irradiance  
      !      
      DO JI=1, KNI
        IF (SUM(ZIRR_DIFF(JI,:))>0.) THEN
        PSCA_SW(JI,:)=ZINT_SCA_SW(JI)*ZIRR_DIFF(JI,:)/SUM(ZIRR_DIFF(JI,:))
        ! XSCA_SW(JI,:)=0.
        END IF 
        IF (SUM(ZIRR_DIR(JI,:))>0.) THEN
        PDIR_SW(JI,:)=ZINT_DIR_SW(JI)*ZIRR_DIR(JI,:)/SUM(ZIRR_DIR(JI,:))
        ! XDIR_SW(JI,:)=(ZINT_DIR_SW(JI)+ZINT_DIR_SW(JI))*ZIRR_DIR(JI,:)/SUM(ZIRR_DIR(JI,:))
        END IF
      ENDDO
      !
    ELSE ! case where only total radiation is know
      !
      ZTCLOUD55(:)= 0. ! cloud optical depth
      CALL IRRADIANCE(ZDATI,ZMU(:),XQA(:)/XRHOA(:),ZD_O3(:),ZD_AE(:),&
                      XPS(:),XTA(:),PPZP_CUT,ZP_CLOUD(:),KCLOUD_TYPE(:),&
                      ZTCLOUD55(:),ZIRR_DIR(:,:),ZIRR_DIFF(:,:))
      ! normalisation by broadband total irradiance  
      DO JI=1, KNI
        IF (SUM(ZIRR_DIFF(JI,:)+ZIRR_DIR(JI,:))>0.) THEN
          PSCA_SW(JI,:)=ZINT_TOT_SW(JI)*ZIRR_DIFF(JI,:)/&
          SUM(ZIRR_DIR(JI,:)+ZIRR_DIR(JI,:))
        END IF 
        IF (SUM(ZIRR_DIR(JI,:)+ZIRR_DIR(JI,:))>0.) THEN
          PDIR_SW(JI,:)=ZINT_TOT_SW(JI)*ZIRR_DIR(JI,:)/&
          SUM(ZIRR_DIR(JI,:)+ZIRR_DIR(JI,:))
        END IF
      ENDDO
      !
    ENDIF ! end of clear cloud/case 
    !!! add 89° threshold for mu to prevent from atmotartes divergence 
    IF (MAXVAL(ZMU(:))<PPMU_THRESHOLD) THEN
      ZINT_SCA_SW(:)=SUM(PSCA_SW(:,:))
      PSCA_SW(:,:)=0.
      PSCA_SW(:,8)=ZINT_SCA_SW(:)
    ENDIF
  ENDIF
ENDIF
!
IF (LHOOK) CALL DR_HOOK('RADIANCE',1,ZHOOK_HANDLE)
END SUBROUTINE 
!###################################################################
SUBROUTINE TAU_CLOUD(PMU,PRATIO,PTCLOUD55)
! diagnose cloud optical depth from diffuse and direct braodband irradiance
!! INPUTS : PMU -> cosine of the solar zenith angle
!!	PRATIO -> ratio of diffuse to total broadband irradiance
!! OUTPUTS : PTCLOUD55 -> cloud optical thickness at 0.55 um
!
USE MODD_CONST_ATM, ONLY :PPCLOUD1, PPCLOUD2, PPCLOUD3, PPCLOUD4, PPCLOUDMAX, PPHUND
IMPLICIT NONE
!
REAL, DIMENSION(:), INTENT(IN):: PMU, PRATIO
REAL, DIMENSION(:), INTENT(OUT):: PTCLOUD55
!
REAL, DIMENSION(SIZE(PMU)) :: ZTEMP
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
IF (LHOOK) CALL DR_HOOK('TAU_CLOUD',0,ZHOOK_HANDLE)
!
ZTEMP(:)= PPCLOUD1*PMU(:)*PMU(:)*PMU(:)+PPCLOUD2*PMU(:)*PMU(:)+PPCLOUD3*PMU(:)+PPCLOUD4
ZTEMP(:)=(1-ZTEMP(:))/(1-MIN(PRATIO(:),1.-1./PPHUND))
ZTEMP(:)=PMU*LOG(MAX(ZTEMP(:),1.))
ZTEMP(:)=MIN(ZTEMP(:),PPCLOUDMAX)
ZTEMP(:)=MAX(ZTEMP(:),0.)
PTCLOUD55(:)=ZTEMP(:)
!
IF (LHOOK) CALL DR_HOOK('TAU_CLOUD',1,ZHOOK_HANDLE)
END SUBROUTINE 
!###################################################################
SUBROUTINE IRRADIANCE(PIDAY,PMU,PHUM,PZO3,PZAE,PZP,PZTEMP,PZP_CUT,PZP_CLOUD,KLCLOUD_TYPE,PTCLOUD55,PIRR_DIR, PIRR_DIFF)
!
! compute direct and diffuse irradiance for a given atmospheric profile
!! INPUTS : PMU -> cosine of the solar zenith angle
!!	PHUM -> specific humidity kg m-3
!! 	PZO3 -> ozone total column atm-cm
!!	PZTEMP -> air temperature, K
!! 	PZAE -> aerosols optical thickness at 550 nm
!! 	PZP -> surface pressure (Pa)
!!	PZP_CUT -> layering (pressure, Pa)
!!	PZP_CLOUD -> cloud bottom pressure, Pa
!!	PTCLOUD55 -> cloud optical thickness at 0.55 um
!!	KLCLOUD_TYPE -> 0 : ice cloud, 1: water cloud
!!	ZIDAY -> day of year 
!! OUTPUTS : PIRR_DIFF -> surface diffuse irradiance 
!!	     PIRR_DIR -> surface direct irradiance 
!
USE MODD_CONST_ATM, ONLY :JPNLYR_CLOUD, JPNLYR_CLEAR, JPNBANDS_ATM,     &
                          PPWAVELENGTHS_ATM, PPWL_REF_CLOUD,PPWRAY_EFF, &
                          PPGRAY, PPALB_SOIL,PPCONV, PPCLOUD_THRES,     &
                          PPHUND, PPHUMA, PPHUMB, PZERO
IMPLICIT NONE 
!
REAL, INTENT(IN) :: PIDAY
REAL,DIMENSION(:), INTENT(IN) :: PMU,PHUM,PZO3, PZAE, PZP,PZP_CLOUD,PZTEMP,PTCLOUD55 ! dimension (npoints)
REAL, DIMENSION(:), INTENT(IN) :: PZP_CUT
INTEGER, DIMENSION(:),INTENT(IN) ::  KLCLOUD_TYPE
REAL, DIMENSION(:,:), INTENT(OUT) :: PIRR_DIR, PIRR_DIFF ! dimension (npoints, nwavelengths)
!
REAL, DIMENSION(SIZE(PMU)) :: ZPHUM, ZMU_TEMP
REAL, DIMENSION(SIZE(PMU)) :: ZQQ55, ZWW55,ZGG55,  ZPALBSOIL, PZSVP
REAL :: ZWL
REAL, DIMENSION(SIZE(PMU), SIZE(PZP_CUT))::ZPZP_CUT
INTEGER :: JI, JICOUNT,ICLOUD, JP,ICLEAR
REAL, DIMENSION(SIZE(PMU)) :: ZTAUGAS_U, ZTAUGAS_O, ZTAUGAS_W,ZTAUAER,ZTAURAY,ZWAER,ZGAER,ZIRR_TOA,ZPSOLFAC
REAL, DIMENSION(SIZE(PMU),JPNLYR_CLEAR) :: ZTAUGAS_TAB,ZTAUAER_TAB,ZTAURAY_TAB
REAL, DIMENSION(SIZE(PMU)) :: ZTCLOUD, ZF_TOT
REAL, DIMENSION(SIZE(PMU),JPNLYR_CLOUD) :: ZTAUPRIME,ZWPRIME,ZGPRIME
REAL, DIMENSION(SIZE(PMU),JPNLYR_CLOUD) :: ZTAU,ZW,ZG
REAL, DIMENSION(SIZE(PMU)) :: ZQQ,ZWW,ZGG, ZTAUMIN
INTEGER, DIMENSION(SIZE(PMU)) :: IATM_LYR
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
IF (LHOOK) CALL DR_HOOK('IRRADIANCE',0,ZHOOK_HANDLE)
!
PIRR_DIFF(:,:)=0.
PIRR_DIR(:,:)=0.
! specific to relative humidity 
call VAPSAT(PZTEMP,PZSVP)
ZPHUM(:)=PPHUND*PHUM(:)*PZP(:)/((PPHUMA + PPHUMB*PHUM(:))*PZSVP(:)*PPHUND)
!ZMU_TEMP(JP)=min(PMU(JP),1.)
!ZMU_TEMP(JP)=max(ZMU_TEMP(JP),1./PPHUND)
ZMU_TEMP(:)=PMU(:)
!
!modify the atmospheric layering according to the surface pressure PZP
DO JI = 1, JPNLYR_CLEAR-1
  DO JP=1,SIZE(ZF_TOT)
    ZPZP_CUT(JP,JI)=PZP(JP)- &
         (((PZERO-PZP_CUT(JI))/(PZERO-PZP_CUT(JPNLYR_CLEAR-1)))&
         *( PZP(JP)-PZP_CUT(JPNLYR_CLEAR-1)))                    
  ENDDO
END DO
!
!
DO JI=1,JPNBANDS_ATM
	ZTAU(:,:)=0.
	ZW(:,:)=0.
	ZG(:,:)=0.
	ZF_TOT(:)=0.
	ZIRR_TOA(:)=0.
	ZWL=PPWAVELENGTHS_ATM(JI)/PPCONV ! wavelentgh in microns
	ZPALBSOIL(:)=PPALB_SOIL(JI)
	! atmospheric properties for clear sky
	CALL OPTICS_CLEAR(PIDAY,ZWL,ZMU_TEMP,ZPHUM,PZO3,PZAE,PZP/PPHUND,      &
	  PZTEMP,ZTAUGAS_U, ZTAUGAS_O, ZTAUGAS_W,ZTAUAER,ZTAURAY,ZWAER,ZGAER, &
	  ZIRR_TOA,ZPSOLFAC)
	CALL SPLIT_NLAYER(PZP/PPHUND,ZPZP_CUT/PPHUND,ZTAUGAS_U, ZTAUGAS_O,    &
	  ZTAUGAS_W,ZTAUAER,ZTAURAY,ZTAUGAS_TAB,ZTAUAER_TAB,ZTAURAY_TAB)
	!
	ZTAU(:,1:JPNLYR_CLEAR)=ZTAUAER_TAB(:,:)+ZTAURAY_TAB(:,:)+ZTAUGAS_TAB(:,:)
  !	
	DO JICOUNT=1,JPNLYR_CLEAR
	  ZW(:,JICOUNT)=(ZTAUAER_TAB(:,JICOUNT)*ZWAER(:)+                     &
	                ZTAURAY_TAB(:,JICOUNT)*PPWRAY_EFF)/ZTAU(:,JICOUNT)
	  ZG(:,JICOUNT)=(ZTAUAER_TAB(:,JICOUNT)*ZWAER(:)*ZGAER(:)+            &
	                ZTAURAY_TAB(:,JICOUNT)*PPWRAY_EFF*PPGRAY)
	  ZG(:,JICOUNT)=ZG(:,JICOUNT)/(ZTAUAER_TAB(:,JICOUNT)*ZWAER(:)+       &
	                ZTAURAY_TAB(:,JICOUNT)*PPWRAY_EFF)
	END DO
	!	
	DO JP=1,SIZE(ZF_TOT)
	  ! 
	  IF (PTCLOUD55(JP)>PPCLOUD_THRES) THEN !cloudy case 	    
	    IATM_LYR(JP)=JPNLYR_CLOUD
	    ZTCLOUD(JP)=0.
	    ZQQ(JP)=0.
	    ZWW(JP)=0.
	    ZGG(JP)=0.
	    ! compute cloud optical properties
	    call CLOUD_PROP(PPWL_REF_CLOUD,KLCLOUD_TYPE(JP),ZQQ55(JP),ZWW55(JP),ZGG55(JP)) !Cloud optical depth at 0.55 microns
	    call CLOUD_PROP(ZWL,KLCLOUD_TYPE(JP),ZQQ(JP),ZWW(JP),ZGG(JP))
	    ZWW(JP)=MIN(PPWRAY_EFF,ZWW(JP)) ! prevent single scattering albedo =1 for ice cloud
	    ZTCLOUD(JP)=PTCLOUD55(JP)*ZQQ(JP)/ZQQ55(JP) ! cloud optical depth
	    	    ! Find index of cloud layer 
	    ICLOUD=MINLOC(ABS(ZPZP_CUT(JP,:)-PZP_CLOUD(JP)),1)
	    ! faire le déplacement des couches + pas vraiment besoin de IATM_LYR en fait
      ZTAU(JP,ICLOUD+1:JPNLYR_CLOUD)=ZTAU(JP,ICLOUD:JPNLYR_CLEAR)
	    ZW(JP,ICLOUD+1:JPNLYR_CLOUD)=ZW(JP,ICLOUD:JPNLYR_CLEAR)
	    ZG(JP,ICLOUD+1:JPNLYR_CLOUD)=ZG(JP,ICLOUD:JPNLYR_CLEAR)
	    ZTAU(JP,ICLOUD)=ZTCLOUD(JP)
	    ZW(JP,ICLOUD)=ZWW(JP)
	    ZG(JP,ICLOUD)=ZGG(JP)
	    !	      
	  ELSE ! no cloudy case 
	    ZTCLOUD(JP)=0.
	    IATM_LYR(JP)=JPNLYR_CLEAR
	  END IF 
	  !
	END DO
	!
  ! delta-eddington approximation
  CALL DELTA(JPNLYR_CLOUD,IATM_LYR,ZTAU,ZW,ZG,ZTAUPRIME,ZWPRIME,ZGPRIME)
  ! surface direct solar irradiance 
  CALL TWO_STR(JPNLYR_CLOUD,IATM_LYR,ZTAUPRIME,ZWPRIME,ZGPRIME,ZPALBSOIL, &
  ZIRR_TOA,ZMU_TEMP,ZF_TOT)
  !
	PIRR_DIR(:,JI)=ZMU_TEMP(:)*ZIRR_TOA(:)*EXP(-(SUM(ZTAU(:,:),2))/ZMU_TEMP(:))
	PIRR_DIFF(:,JI)=MAX(0.,ZF_TOT(:)-PIRR_DIR(:,JI)) !prevent any negative values in the absorption bands 
END DO 
!
IF (LHOOK) CALL DR_HOOK('IRRADIANCE',1,ZHOOK_HANDLE)
END SUBROUTINE
!###################################################################
SUBROUTINE CLOUD_PROP(PWL,KCLOUD_TYPE,PZQQ,PZWW,PZGG)
!
! compute cloud optical properties for each wavelength 
!! INPUTS : PWL -> wavelenght (microns)
!!	KCLOUD_TYPE -> 0-> ice cloud, 1-> water cloud
!! OUTPUTS : PZQQ -> extinction efficiency
!!	     PZWW -> single scattering albedo 
!!	     PZGG -> assymetry factor 
!
USE MODD_CONST_ATM, ONLY :JPNCLOUDWV, PPCLOUDWV1, PPCLOUDWV2, JPNCLOUDWV2,PPQQ, PPQQI, PPWW, PPWWI, PPGG, PPGGI
IMPLICIT NONE 
REAL, INTENT(IN) :: PWL 
INTEGER, INTENT(IN) :: KCLOUD_TYPE
REAL, INTENT(OUT) :: PZQQ,PZWW,PZGG
!
REAL :: ZWMIN, ZWMAX,ZWSTEP,ZFW
INTEGER :: IW 
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
IF (LHOOK) CALL DR_HOOK('CLOUD_PROP',0,ZHOOK_HANDLE)
! intialisation 
ZWMIN=0.
ZWMAX=0.
ZWSTEP=0.
ZFW=0.
IW=0
PZQQ=0.
PZWW=0.
PZGG=0.
! computation of reference wavalenght
ZWMIN=LOG(PPCLOUDWV1)
ZWMAX=LOG(PPCLOUDWV2)
ZWSTEP=(ZWMAX-ZWMIN)/JPNCLOUDWV2
!
! computation of wavelenght increment for spectral interpolation
!
ZFW=1.+(LOG(PWL)-ZWMIN)/ZWSTEP
ZFW=MIN(MAX(ZFW,1.),FLOAT(JPNCLOUDWV2))
IW=INT(ZFW)
ZFW=ZFW-IW
!
IF (KCLOUD_TYPE==1) THEN
  !      
  PZQQ=PPQQI(IW)*(1.-ZFW)+PPQQI(IW+1)*ZFW
  PZWW=PPWWI(IW)*(1.-ZFW)+PPWWI(IW+1)*ZFW
  PZGG=PPGGI(IW)*(1.-ZFW)+PPGGI(IW+1)*ZFW
  !
ELSE
  !  
  PZQQ=PPQQ(IW)*(1.-ZFW)+PPQQ(IW+1)*ZFW
  PZWW=PPWW(IW)*(1.-ZFW)+PPWW(IW+1)*ZFW
  PZGG=PPGG(IW)*(1.-ZFW)+PPGG(IW+1)*ZFW
  !
END IF 
!
IF (LHOOK) CALL DR_HOOK('CLOUD_PROP',1,ZHOOK_HANDLE)
END SUBROUTINE
!###################################################################
SUBROUTINE SPLIT_NLAYER(PZP,PZP_CUT,PTAUGAS_U, PTAUGAS_O, PTAUGAS_W,PTAUAER,PTAURAY,PTAUGAS_TAB,PTAUAER_TAB,PTAURAY_TAB)
!
! split optical depth in N layers at pressure ZP_CUT 
!! INPUTS : PZP -> Pression (mb)
!!	PZP_CUT -> pressure layering (mb)
!!	PTAUGAS_O,W,U -> optical depth of ozone, water vapour and unfiformely mixed gaz
!!	     PTAUER -> aerosols optical depth 
!!	     PTAURAY -> optical depth for rayleigh scattering
!! OUTPUTS : PTAUGAS_TAB -> optical depth of ozone, water vapour and unfiformely mixed gaz two layer
!!	     PTAUER -> aerosols optical depth two layers
!!	     PTAURAY -> optical depth for rayleigh scattering two layers
!
! Note that layer 2 is bottom layer 
!
USE MODD_CONST_ATM, ONLY : JPNPREF,PPPREF,PPVAPREF,PPOZOREF,PPAERREF, JPNLYR_CLEAR ! vertical profile of water vapour, aerosols and ozone
IMPLICIT NONE
REAL,DIMENSION(:), INTENT(IN) :: PZP, PTAUGAS_O,PTAUGAS_U,PTAUGAS_W, PTAUAER, PTAURAY
REAL, DIMENSION(:,:), INTENT(IN) :: PZP_CUT
REAL, DIMENSION (:,:), INTENT(OUT) :: PTAUGAS_TAB, PTAUAER_TAB, PTAURAY_TAB
!
REAL, DIMENSION(SIZE(PZP),JPNLYR_CLEAR) :: ZPOIDS_0, ZPOIDS_W, ZPOIDS_AER
INTEGER :: JICOUNT
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
IF (LHOOK) CALL DR_HOOK('SPLIT_NLAYER',0,ZHOOK_HANDLE)
!
PTAURAY_TAB(:,:)=0.
PTAUAER_TAB(:,:)=0.
PTAUGAS_TAB(:,:)=0.
DO JICOUNT=2,JPNLYR_CLEAR-1
  ! Rayleigh scattering only pressure 
  PTAURAY_TAB(:,JICOUNT)=ABS(PZP_CUT(:,JICOUNT)-PZP_CUT(:,JICOUNT-1))/PZP(:)*PTAURAY(:)
  ! Uniformely mixed gaz only pressure
  PTAUGAS_TAB(:,JICOUNT)=ABS(PZP_CUT(:,JICOUNT)-PZP_CUT(:,JICOUNT-1))/PZP(:)*PTAUGAS_U(:)
END DO 
PTAUGAS_TAB(:,1)=(1.-PZP_CUT(:,1)/PZP(:))*PTAUGAS_U(:)
PTAURAY_TAB(:,1)=(1.-PZP_CUT(:,1)/PZP(:))*PTAURAY(:)
PTAUGAS_TAB(:,JPNLYR_CLEAR)=(PZP_CUT(:,JPNLYR_CLEAR-1)/PZP(:))*PTAUGAS_U(:)
PTAURAY_TAB(:,JPNLYR_CLEAR)=(PZP_CUT(:,JPNLYR_CLEAR-1)/PZP(:))*PTAURAY(:)
!
ZPOIDS_0(:,:)=0.
ZPOIDS_W(:,:)=0.
ZPOIDS_AER(:,:)=0.
!
! Ozone
CALL WEIGHT(PZP,PZP_CUT,PPPREF,PPOZOREF,ZPOIDS_0)
! Water vapour
CALL WEIGHT(PZP,PZP_CUT,PPPREF,PPVAPREF,ZPOIDS_W)
!Aerosols
CALL WEIGHT(PZP,PZP_CUT,PPPREF,PPAERREF,ZPOIDS_AER)
!
DO JICOUNT=1, JPNLYR_CLEAR
	PTAUGAS_TAB(:,JICOUNT)=PTAUGAS_TAB(:,JICOUNT)+ZPOIDS_0(:,JICOUNT)*PTAUGAS_O(:)
	PTAUGAS_TAB(:,JICOUNT)=PTAUGAS_TAB(:,JICOUNT)+ZPOIDS_W(:,JICOUNT)*PTAUGAS_W(:)
	PTAUAER_TAB(:,JICOUNT)=ZPOIDS_AER(:,JICOUNT)*PTAUAER(:)
ENDDO 
IF (LHOOK) CALL DR_HOOK('SPLIT_NLAYER',1,ZHOOK_HANDLE)

END SUBROUTINE
!###################################################################
SUBROUTINE WEIGHT(PZP,PZP_CUT,PP_TAB,PPROFILE_TAB,PPOIDS)
!
! calculate optical depth weigth of each layer 
!! INPUTS : PZP -> Surface Pressure (mb)
!!	PZP_CUT -> pressure layering (mb)
!!	P_TAB -> pressure reference vertical profile
!!	PROFILE_TAB-> content reference vertical profile
!! OUTPUTs : POIDS -> optical weigth of the n layers 

USE MODD_CONST_ATM, ONLY : JPNPREF, JPNLYR_CLEAR
IMPLICIT NONE
REAL, DIMENSION(:),INTENT(IN) :: PZP 
REAL, DIMENSION(:,:), INTENT(IN) :: PZP_CUT
REAL, DIMENSION (JPNPREF), INTENT(IN) :: PP_TAB, PPROFILE_TAB
REAL, DIMENSION(:,:), INTENT(OUT) :: PPOIDS

INTEGER, DIMENSION(SIZE(PZP)) ::I_ZP, I_ZPCUT,I_ZPCUT2
INTEGER :: JICOUNT,  JP 
INTEGER, DIMENSION(SIZE(PZP)) ::INSIZE 
REAL, DIMENSION(SIZE(PZP)) :: ZNORMAL, ZP1, ZP2
REAL, DIMENSION(JPNPREF) :: ZTAB
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
IF (LHOOK) CALL DR_HOOK('WEIGHT',0,ZHOOK_HANDLE)
!
PPOIDS(:,:)=0.
! Indices closer to ZP and Z_CUT 
I_ZP(:)=0
!
DO JP=1, SIZE(I_ZP)
  ZTAB(:)=PP_TAB(:)-PZP(JP)
  I_ZP(JP)=MINLOC(ABS(ZTAB),1)
  I_ZPCUT(JP)=MINLOC(ABS(PP_TAB(:)-PZP_CUT(JP,1)),1)
END DO
!
INSIZE(:)=SIZE(PP_TAB)
! Normalization constante
ZNORMAL=0.
CALL TRAPZ(I_ZP,INSIZE,PP_TAB,PPROFILE_TAB,ZNORMAL)
CALL TRAPZ(I_ZP,I_ZPCUT,PP_TAB,PPROFILE_TAB,PPOIDS(:,1))
PPOIDS(:,1)=PPOIDS(:,1)/ZNORMAL(:)
!
DO JICOUNT=2, JPNLYR_CLEAR-1
  ZP1=0.
  DO JP=1, SIZE(I_ZP)
    I_ZPCUT(JP)=MINLOC(ABS(PP_TAB-PZP_CUT(JP,JICOUNT)),1)
    I_ZPCUT2(JP)=MINLOC(ABS(PP_TAB-PZP_CUT(JP,JICOUNT-1)),1)
  END DO
  CALL TRAPZ(I_ZPCUT2,I_ZPCUT,PP_TAB,PPROFILE_TAB,ZP1)
  PPOIDS(:,JICOUNT)=ZP1(:)/ZNORMAL(:)
END DO
!
PPOIDS(:,JPNLYR_CLEAR)=1.-SUM(PPOIDS(:,:),2)
!
IF (LHOOK) CALL DR_HOOK('WEIGHT',1,ZHOOK_HANDLE)
END SUBROUTINE
!###################################################################
SUBROUTINE TRAPZ(KI_DEP,KI_FIN,PLEVEL,PTAB,PRES)
!
! calculate integral of an array using trapezoidal rule  
!!! INPUTS : KI_DEP, KI_FIN : start and end indices of the integral
!!	    TAB -> array to be integrated
!!          LEVEL -> x-axis
!! OUTPUT : RES : integral 
USE MODD_CONST_ATM, ONLY : JPNPREF
IMPLICIT NONE 
INTEGER, DIMENSION(:),INTENT(IN) :: KI_DEP, KI_FIN
REAL, DIMENSION(JPNPREF), INTENT(IN) :: PTAB, PLEVEL
REAL, DIMENSION(:),INTENT(OUT) :: PRES
!
INTEGER :: ICOUNT, JP
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
IF (LHOOK) CALL DR_HOOK('TRAPZ',0,ZHOOK_HANDLE)
!
PRES=0.
DO JP=1, SIZE(KI_DEP)
  DO ICOUNT=KI_DEP(JP),KI_FIN(JP)-1
	   PRES(JP)=PRES(JP)+(PLEVEL(ICOUNT+1)-PLEVEL(ICOUNT))*0.5*(PTAB(ICOUNT+1)+PTAB(ICOUNT))
  END DO 
END DO
PRES=ABS(PRES)
!
IF (LHOOK) CALL DR_HOOK('TRAPZ',1,ZHOOK_HANDLE)
END SUBROUTINE
!###################################################################
SUBROUTINE OPTICS_CLEAR(PIDAY,PWL,PMU,PHUM,PO3,PAE,PZP,PTEMP,PTAUGAS_U, &
                       PTAUGAS_O, PTAUGAS_W,PTAUAER,PTAURAY,& 
                       PWAER,PGAER,PIRR_TOA,PSOLFAC)
!
! computation of optical parameters for clear sky 
!
!! INPUTS : PMU -> cosine of the solar zenith angle
!!	PWL -> microns nm
!!	PHUM -> relative humidity %%
!! 	PO3 -> ozone total column atm-cm
!!	PTEMP -> air temperature, K
!! 	PAE -> aerosols optical thickness at 550 nm
!! 	PZP -> Pression (mb)
!!	PIDAY -> day of year 
!! OUTPUTS : PTAUGAS -> optical depth of ozone, water vapour and unfiformely mixed gaz
!!	     PTAUER -> aerosols optical depth 
!!	     PTAURAY -> optical depth for rayleigh scattering
!!	     PWAER -> single scattering albedo of aerosols
!!	     PGAER -> assymetry factor of aerosols
!!	     PZIRR_SUN_TOA -> top of atmosphere solar irradiance  W/m²/microns
!!	    
!
USE MODD_CONST_ATM, ONLY : PPPZERO,PPWAVE_SUN_ABS,PPAO,PPAU,PPAW,JPNBANDS_REF_ABS,&
PPWAVE_SUN_REF,JPNBANDS_REF_SUN,PPEXTRA_SUN, PPECCEN,  &
PPDEGPDAY,PPDAYPH, PPPAVO,PPAMFAC, PPMD, PPEARTH_G,    &
PPA1, PPB1, PPC1,PPCLOUD_THRES,PPCONV, PPDEG, PPTEN,   &
PPA2, PPB2,PPC2,PPD2,PPE2,PPF2, PPWAT1, PPHUND, PPWATA,&
PPWATB, PPWATC, PPHO,PPO1,PPO2,PPU1, PPU2, PPU3, PPAER1, &
PPAER2,PPAER3,PPWL_REF_CLOUD, PPPI, ZRAY1, ZRAY2, PPWL_REF_RAY,PPTEMP2
!
IMPLICIT NONE 
!
REAL, INTENT(IN) :: PIDAY,PWL
REAL, DIMENSION(:), INTENT(IN) :: PMU,PHUM,PO3,PAE,PZP,PTEMP
REAL, DIMENSION(:), INTENT(OUT) :: PTAUGAS_U, PTAUGAS_O, PTAUGAS_W,PTAUAER,&
                                    PTAURAY,PWAER,PGAER,PIRR_TOA,PSOLFAC
!
REAL :: ZAO, ZAU,ZAW, ZRSUN,  ZSIG, ZEXTRA
REAL, DIMENSION(SIZE(PZP)) :: ZABS_AER,ZEXT_AER,Z_GAER, ZABS_AER_55,ZEXT_AER_55,Z_GAER_55, ZWAER_55 
REAL, DIMENSION(SIZE(PZP)) :: ZZEN,ZM, ZMPRIME,ZW, ZMO,ZSVP
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
IF (LHOOK) CALL DR_HOOK('OPTICS_CLEAR',0,ZHOOK_HANDLE)
!
ZEXTRA=0.
PTAUGAS_O=0.
PTAUGAS_U=0.
PTAUGAS_W=0.
PTAUAER=0.
PTAURAY=0.
PWAER=0.
PGAER=0.
PSOLFAC=0.
ZRSUN=0.
PIRR_TOA=0.
ZM=0.
ZMPRIME=0.
ZW=0.
ZMO=0.
ZSVP=0.
! Calculate TOA irradiance W/m²/microns
CALL LOCATE(PWL*PPCONV,JPNBANDS_REF_SUN,PPWAVE_SUN_REF,PPEXTRA_SUN,ZEXTRA) 
!! earth-sun distances
ZRSUN=1.-PPECCEN*cos(PPDEGPDAY*(PIDAY-PPDAYPH)*PPPI/PPDEG)    
PSOLFAC(:)=1./ZRSUN**2     
PIRR_TOA(:)=PSOLFAC(:)*ZEXTRA
!PIRR_TOA=ZEXTRA ! disable if you want to account for earth-sun distance variations
!
! select absorption coefficient for ozone, water vapour and uniformely mixed gaz
CALL LOCATE(PWL,JPNBANDS_REF_ABS,PPWAVE_SUN_ABS,PPAO,ZAO)
CALL LOCATE(PWL,JPNBANDS_REF_ABS,PPWAVE_SUN_ABS,PPAU,ZAU)
CALL LOCATE(PWL,JPNBANDS_REF_ABS,PPWAVE_SUN_ABS,PPAW,ZAW)
!
! Relative air mass
ZZEN(:)=ACOS(PMU(:))
ZM(:)=1./(PMU(:)+PPA1*(PPB1-ZZEN(:))**(-PPC1))!Kasten&young1989
ZMPRIME(:)=ZM(:)*PZP(:)/PPPZERO
!
!Rayleigh scattering 
! Bucholtz et al. 1995, rayleigh molecular scattering
ZSIG=0.
IF (PWL<PPWL_REF_RAY) THEN 
    ZSIG=ZRAY1(1)*PWL**(-(ZRAY1(2)+ZRAY1(3)*PWL+ZRAY1(4)/PWL))
ELSE 
    ZSIG=ZRAY2(1)*PWL**(-(ZRAY2(2)+ZRAY2(3)*PWL+ZRAY2(4)/PWL))
END IF 
!
! Rayleigh optical depth calculation from  Bodhaine et al., 1999 
PTAURAY(:)=ZSIG*PPPAVO*PPPZERO*PPTEN/(PPEARTH_G*PPMD)
PTAURAY(:)=PTAURAY(:)*PZP(:)/PPPZERO*PTEMP(:)/PPTEMP2 ! modulation of Rayleigh absorption according to elevation and temperature
!
! link between water vapor column and RH 
call VAPSAT(PTEMP(:),ZSVP(:))
ZW(:) = PPWAT1*(PHUM(:)/PPHUND)*ZSVP(:)/PTEMP(:) !Prata 1996
!
!Water vapor absorption 
PTAUGAS_W(:)=PPWATA*ZAW*ZW(:)*ZM(:)/(1.+PPWATB*ZAW*ZW(:)*ZM(:))**(PPWATC)
!
! Ozone absorption 
! Ozone mass (Iqbal 1983)
ZMO(:)=(1.+PPHO/PPO1)/(PMU(:)**PPO2+PPO2*PPHO/PPO1)**(1./PPO2)
PTAUGAS_O(:)= ZAO*PO3(:)*ZMO(:)*0.5
! Uniformely mixed gaz (leckner)
PTAUGAS_U(:)= (PPU1*ZAU*ZMPRIME(:)/(1.+PPU2*ZAU*ZMPRIME(:))**(PPU3))
!
! aerosols scattering and absportion
! interpolation as a function of zhum
CALL AER_HUMID(PHUM(:),PWL,ZABS_AER,ZEXT_AER,Z_GAER)
!
PGAER=Z_GAER
PWAER(:)=1.-ZABS_AER(:)/ZEXT_AER(:)
! scaling at 0.55 nm (not use)
CALL AER_HUMID(PHUM(:),PPWL_REF_CLOUD,ZABS_AER_55,ZEXT_AER_55,Z_GAER_55)
!
PTAUAER(:)=PAE*ZEXT_AER(:)/ZEXT_AER_55(:)
!
! ! other formulation for aerosols (Shettle and Fenn, 1975)
! IF (PWL<PPAER3) THEN
!       PTAUAER=PAE*PPWL_REF_CLOUD**PPAER1/(PWL)**PPAER2
! ELSE
!       PTAUAER=PAE*(PPWL_REF_CLOUD/PWL)**PPAER1
! ENDIF
! careful no dependance at elevation 
! only taugas and tauray are dependant on elevation
! no stratospheric aerosols
IF (LHOOK) CALL DR_HOOK('OPTICS_CLEAR',1,ZHOOK_HANDLE)
!
END SUBROUTINE
!###################################################################
SUBROUTINE AER_HUMID(PHUM,PWL,PABS_AER,PEXT_AER,PGAER)
!
! computation of single scattering properties of aerosols scaled by humidity
!
! inputs : PHUM -> relative humidity 0-100 %
!	PWL -> wavalenght in microns
! outputs : PABS_AER -> absorption efficiency
!	PEXT_AER -> extinction efficiency
! 	PGAER -> assymetry factor 
USE MODD_CONST_ATM, ONLY :JPNBANDS_REF_AER,PPWAVE_SUN_AER,PPRUR_EXT_1,PPRUR_EXT_2,PPRUR_EXT_3,&
PPRUR_EXT_4,PPRUR_ABS_1,PPRUR_ABS_2,PPRUR_ABS_3,PPRUR_ABS_4,PPRUR_G_1,PPRUR_G_2,PPRUR_G_3,PPRUR_G_4,&
PPSMALL,PPHUM_ZONE,PPHUND, PPHUM1,PPHUM2
IMPLICIT NONE 
!
REAL, INTENT(IN) ::PWl
REAL, DIMENSION(:),INTENT(IN) :: PHUM
REAL, DIMENSION(:),INTENT(OUT) ::PABS_AER,PEXT_AER,PGAER 
!
REAL :: ZEX1, ZEX2, ZW1, ZW2, ZG1, ZG2 
INTEGER :: JJ, JP
REAL, DIMENSION(SIZE(PHUM)) :: ZHUM,ZWT
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
IF (LHOOK) CALL DR_HOOK('AER_HUMID',0,ZHOOK_HANDLE)
!
ZHUM(:)=PHUM(:)/PPHUND
!
DO JP=1, size(ZHUM)
  IF (ZHUM(JP).lt.PPHUM1) THEN
    JJ=1
  ELSE
    IF (ZHUM(JP).lt.PPHUM2) THEN
      JJ=2
    ELSE
      JJ=3
    END IF 
  END IF
  !
  ZWT(JP)= (ZHUM(JP)-PPHUM_ZONE(JJ))/(PPHUM_ZONE(JJ+1)-PPHUM_ZONE(JJ))
  !
  IF (JJ==1) THEN
    CALL LOCATE(PWL,JPNBANDS_REF_AER,PPWAVE_SUN_AER,PPRUR_EXT_1,ZEX1)
    CALL LOCATE(PWL,JPNBANDS_REF_AER,PPWAVE_SUN_AER,PPRUR_EXT_2,ZEX2)
    CALL LOCATE(PWL,JPNBANDS_REF_AER,PPWAVE_SUN_AER,PPRUR_ABS_1,ZW1)
    CALL LOCATE(PWL,JPNBANDS_REF_AER,PPWAVE_SUN_AER,PPRUR_ABS_2,ZW2)
    CALL LOCATE(PWL,JPNBANDS_REF_AER,PPWAVE_SUN_AER,PPRUR_G_1,ZG1)
    CALL LOCATE(PWL,JPNBANDS_REF_AER,PPWAVE_SUN_AER,PPRUR_G_2,ZG2) 
  ELSEIF (JJ==2) THEN 
    CALL LOCATE(PWL,JPNBANDS_REF_AER,PPWAVE_SUN_AER,PPRUR_EXT_2,ZEX1)
    CALL LOCATE(PWL,JPNBANDS_REF_AER,PPWAVE_SUN_AER,PPRUR_EXT_3,ZEX2)
    CALL LOCATE(PWL,JPNBANDS_REF_AER,PPWAVE_SUN_AER,PPRUR_ABS_2,ZW1)
    CALL LOCATE(PWL,JPNBANDS_REF_AER,PPWAVE_SUN_AER,PPRUR_ABS_3,ZW2)
    CALL LOCATE(PWL,JPNBANDS_REF_AER,PPWAVE_SUN_AER,PPRUR_G_2,ZG1)
    CALL LOCATE(PWL,JPNBANDS_REF_AER,PPWAVE_SUN_AER,PPRUR_G_3,ZG2)
  ELSEIF (JJ==3) THEN 
    CALL LOCATE(PWL,JPNBANDS_REF_AER,PPWAVE_SUN_AER,PPRUR_EXT_3,ZEX1)
    CALL LOCATE(PWL,JPNBANDS_REF_AER,PPWAVE_SUN_AER,PPRUR_EXT_4,ZEX2)
    CALL LOCATE(PWL,JPNBANDS_REF_AER,PPWAVE_SUN_AER,PPRUR_ABS_3,ZW1)
    CALL LOCATE(PWL,JPNBANDS_REF_AER,PPWAVE_SUN_AER,PPRUR_ABS_4,ZW2)
    CALL LOCATE(PWL,JPNBANDS_REF_AER,PPWAVE_SUN_AER,PPRUR_G_3,ZG1)
    CALL LOCATE(PWL,JPNBANDS_REF_AER,PPWAVE_SUN_AER,PPRUR_G_4,ZG2)
  END IF 
  !
  ZEX1=MAX(ZEX1,PPSMALL)
  ZEX2=MAX(ZEX2,PPSMALL)
  ZW1=MAX(ZW1,PPSMALL)
  ZW2=MAX(ZW2,PPSMALL)
  ZG1=MAX(ZG1,PPSMALL)
  ZG2=MAX(ZG2,PPSMALL)
  !
  PABS_AER(JP)= ZW1*(ZW2/ZW1)**ZWT(JP)
  PEXT_AER(JP)=ZEX1*(ZEX2/ZEX1)**ZWT(JP)
  PGAER(JP)=ZG1*(ZG2/ZG1)**ZWT(JP)
  !
END DO
IF (LHOOK) CALL DR_HOOK('AER_HUMID',1,ZHOOK_HANDLE)
END SUBROUTINE
!###################################################################
SUBROUTINE  VAPSAT(PZTEMP,PZSVP)
!
! computation of saturation vapor pressure 
!
! input PZTEMP : temperature (K)
!output PZSVP : saturation vapour pressure (mb)
USE MODD_CONST_ATM, ONLY :PPTEMP1, PPTEMP2,PPPA0W,PPPA0
IMPLICIT NONE 
!
REAL, DIMENSION(:), INTENT(IN) :: PZTEMP
REAL, DIMENSION(:),INTENT(OUT) :: PZSVP
!
REAL ::  ZZTEMP
INTEGER :: JP
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
IF (LHOOK) CALL DR_HOOK('VAPSAT',0,ZHOOK_HANDLE)
!
DO JP=1, size(PZTEMP)
  IF (PZTEMP(JP) < PPTEMP1) PZSVP(JP)=0.
  IF (PZTEMP(JP)>PPTEMP2) THEN
	  ZZTEMP=PZTEMP(JP)
	  PZSVP(JP) = PPPA0W(1)+ZZTEMP*(PPPA0W(2)+ZZTEMP*(PPPA0W(3)+ZZTEMP*(PPPA0W(4)&
	        +ZZTEMP*(PPPA0W(5)+ZZTEMP*(PPPA0W(6)+ZZTEMP*PPPA0W(7))))))	
  ELSE
	  ZZTEMP = PZTEMP(JP) - PPTEMP2
	  PZSVP(JP) = PPPA0(1)+ZZTEMP*(PPPA0(2)+ZZTEMP*(PPPA0(3)+ZZTEMP*(PPPA0(4)&
	        +ZZTEMP*(PPPA0(5)+ZZTEMP*(PPPA0(6)+ZZTEMP*PPPA0(7))))))
  END IF
END DO
!
IF (LHOOK) CALL DR_HOOK('VAPSAT',1,ZHOOK_HANDLE)
END SUBROUTINE
!###################################################################
SUBROUTINE PRESSURE(PZALT,PZP)
!
! computation of pressure as a function of elevation  
!
! INPUTS 	PZALT= height a.s.l. (m)
! OUTPUTS	PZP = pressure in hPa
USE MODD_CONST_ATM, ONLY : PPPZERO, PPATMOS_R,PPATMOS_T0,PPEARTH_G,PPSTLAPSE,PPMD,PPREARTH, PPCONV
IMPLICIT NONE 
REAL,DIMENSION(:), INTENT(IN) :: PZALT
REAL,DIMENSION(:), INTENT(OUT) :: PZP
REAL, DIMENSION(SIZE(PZALT)) :: ZH1, ZHB
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
IF (LHOOK) CALL DR_HOOK('PRESSURE',0,ZHOOK_HANDLE)
!
! pressure as a function of altitude US standard atmosphere p. 12 & 3
! stlapse*1000 as in table 4, p. 3  temperature gradient is given in K/km
ZH1(:)=(PPREARTH*PZALT(:))/(PPREARTH+PZALT(:))!Geopotential height	
ZHB(:)=0.
PZP(:) = PPPZERO*( PPATMOS_T0/(PPATMOS_T0+PPSTLAPSE*(ZH1(:)-ZHB(:)))) &
                    **( (PPEARTH_G*PPMD)/(PPATMOS_R*PPSTLAPSE*PPCONV) )
!                    
IF (LHOOK) CALL DR_HOOK('PRESSURE',1,ZHOOK_HANDLE)
END SUBROUTINE
!###################################################################
SUBROUTINE LOCATE(PWL,KNUM,PWAVE,PTAB, PRES)
!
! find values in wavelength table
!
! input 
! PWL : wavelenght in microns
! KNUM : number of element in ZWAVE
! PWAVE : input wavelength reference tab
! PTAB : input value reference tab
! output 
! PRES : value of ZTAB interpolated for ZWL
!
IMPLICIT NONE 
!
REAL, INTENT(IN) :: PWL  
INTEGER, INTENT(IN) :: KNUM
REAL, DIMENSION(KNUM), INTENT(IN) :: PTAB, PWAVE
REAL, INTENT(OUT) :: PRES
!
INTEGER :: JBANDREF
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
IF (LHOOK) CALL DR_HOOK('LOCATE',0,ZHOOK_HANDLE)
!
PRES=0.
IF (PWAVE(1)>PWL) THEN 
  PRES= PTAB(1)
ELSEIF (PWAVE(KNUM)<=PWL) THEN
  PRES= PTAB(KNUM)
ELSE
  DO JBANDREF=2,KNUM
    !
	  IF (PWAVE(JBANDREF)>PWL) THEN  
	      PRES=  ((PWL-PWAVE (JBANDREF-1))* PTAB(JBANDREF)+&
			        (PWAVE(JBANDREF)-PWL)* PTAB(JBANDREF-1))/&
			        (PWAVE(JBANDREF)-PWAVE(JBANDREF-1))
	      EXIT
	  END IF
	  !
  END DO
END IF
!
IF (LHOOK) CALL DR_HOOK('LOCATE',1,ZHOOK_HANDLE)
END SUBROUTINE 
!###################################################################
SUBROUTINE DELTA(KJNLYR,IATM_LYR,PTAU,PW,PG,PTAUPRIME,PWPRIME,PGPRIME)
!
! Perform delta-eddington approximation
!
!! inputs PTAU : optical depth
!!	PW : single scatering albedo
!!	PG : assymetry factor
!! outpus PTAUPRIME : optical depth modified by delta eddington 
!!	  PGPRIME : assymetry factor modified by delta eddington approximation
!!	  PWPRIME : single scattering albedo modified by delta-eddigton approximation
!
IMPLICIT NONE
! 
INTEGER, INTENT(IN) :: KJNLYR
INTEGER, DIMENSION(:),INTENT(IN) :: IATM_LYR
REAL, DIMENSION(:,:), INTENT(IN) :: PG,PW,PTAU
REAL, DIMENSION(:,:), INTENT(OUT) :: PTAUPRIME, PWPRIME, PGPRIME
!
INTEGER :: JLAYER, JP
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
IF (LHOOK) CALL DR_HOOK('DELTA',0,ZHOOK_HANDLE)
!
PTAUPRIME(:,:)=0.
PWPRIME(:,:)=0.
PGPRIME(:,:)=0.
!
DO JLAYER=1,KJNLYR
  DO JP = 1, SIZE(PG,1)
    IF (JLAYER<=IATM_LYR(JP)) THEN
      PGPRIME(JP,JLAYER) = PG(JP,JLAYER)/(1.+PG(JP,JLAYER))!doc equation 12
      PWPRIME(JP,JLAYER)=(1.-PG(JP,JLAYER)**2)*PW(JP,JLAYER)  / &
                   (1.-PW(JP,JLAYER)*PG(JP,JLAYER)**2)!doc equation 16
      PWPRIME(JP,JLAYER)= MIN(PWPRIME(JP,JLAYER), 0.9999999)                    
      PTAUPRIME(JP,JLAYER)=PTAU(JP,JLAYER)*(1.-PW(JP,JLAYER)*PG(JP,JLAYER)**2)
    ENDIF
  ENDDO
END DO
!
IF (LHOOK) CALL DR_HOOK('DELTA',1,ZHOOK_HANDLE)
END SUBROUTINE
!###################################################################
SUBROUTINE TWO_STR(KJNLYR,IATM_LYR,PTAUPRIME,PWPRIME,PGPRIME,PALB_SOIL,PF0,PMU0, PZF_TOT)
!
! Solves the two-stream equation 
!
!! inputs  KJNLYR : number of atmospheric layer 
!! 	  PTAUPRIME : optical depth modified by delta eddington for each layer
!!	  PGPRIME : assymetry factor modified by delta eddington approximation for each layer
!!	  PWPRIME : single scattering albedo modified by delta-eddigton approximation for each layer
!!	PALB_SOIL : soil spectral albedo 
!!	PF0 : indicent TOA irradiance (W/m²/microns)
!!	ZMU0 : cosine of solar zenith angle
!! output PZF_TOT : total irradiance (diffuse + directe) at bottom (W/m²/microns)
USE MODD_CONST_ATM, ONLY : PPPI
USE MODI_TRIDIAG_GROUND_SNOWCRO
IMPLICIT NONE 
!
INTEGER, INTENT(IN) :: KJNLYR 
INTEGER, DIMENSION(:),INTENT(IN) :: IATM_LYR
REAL, DIMENSION(:,:), INTENT(IN) :: PTAUPRIME, PWPRIME, PGPRIME
REAL, DIMENSION(:),INTENT(IN) ::  PF0, PMU0 ,PALB_SOIL
REAl, DIMENSION(:), INTENT(OUT)::PZF_TOT
!
REAL, DIMENSION(size(PTAUPRIME,1),KJNLYR ) :: ZGAMMA1, ZGAMMA2, ZGAMMA3, ZGAMMA4,ZTAU_TOT
REAL :: ZDGP, ZDGM, ZFDIAG2 
REAL, DIMENSION(size(PTAUPRIME,1),KJNLYR ) :: ZKESTAR, ZALBEDO, ZGP, ZGM, ZFDIAG, ZXA, ZXB
REAL, DIMENSION(size(PTAUPRIME,1),2*KJNLYR ) :: ZVECTOR, ZX0_DIR
REAL, DIMENSION(size(PTAUPRIME,1),2*KJNLYR ) :: ZD ! Diagonal of the matrix 
REAL, DIMENSION(size(PTAUPRIME,1),2*KJNLYR) :: ZDM,ZDP ! Diagonal of the matrix 
INTEGER :: JLAYER, JP
REAL :: ZSTAR,ZDEXP 
INTEGER, DIMENSION(size(PTAUPRIME,1)) :: IATM_LYR_TEMP
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
IF (LHOOK) CALL DR_HOOK('TWO_STR',0,ZHOOK_HANDLE)
!
ZGAMMA1(:,:)=0.
ZGAMMA2(:,:)=0.
ZGAMMA3(:,:)=0.
ZGAMMA4(:,:)=0.
ZGP(:,:)=0.
ZGM(:,:)=0.
ZKESTAR(:,:)=0.
ZALBEDO(:,:)=0.
ZGM(:,:)=0.
ZGP(:,:)=0.
ZTAU_TOT(:,:)=0.
!
DO JLAYER=1,KJNLYR
  DO JP=1,SIZE(PMU0)
		IF (JLAYER<=IATM_LYR(JP)) THEN
			ZGAMMA1(JP,JLAYER)= 0.25*(7.-PWPRIME(JP,JLAYER)*(4.+3.*PGPRIME(JP,JLAYER))) !doc equation 38
			ZGAMMA2(JP,JLAYER)=-0.25*(1.-PWPRIME(JP,JLAYER)*(4.-3.*PGPRIME(JP,JLAYER)))!doc equation 39
			ZKESTAR(JP,JLAYER)=SQRT(ZGAMMA1(JP,JLAYER)**2-ZGAMMA2(JP,JLAYER)**2)!doc equation 42
			ZALBEDO(JP,JLAYER)=(ZGAMMA1(JP,JLAYER)-ZKESTAR(JP,JLAYER))/ZGAMMA2(JP,JLAYER)!doc equation 43
      ! GP_GM_VECTORS
			ZGAMMA4(JP,JLAYER)=0.25*(2.+3.*PGPRIME(JP,JLAYER)*PMU0(JP)) !doc equation 27
			ZGAMMA3(JP,JLAYER)=0.25*(2.-3.*PGPRIME(JP,JLAYER)*PMU0(JP)) !doc equation 28
      !
			ZGM(JP,JLAYER)=(PMU0(JP)**2)*PWPRIME(JP,JLAYER)/((ZKESTAR(JP,JLAYER)*PMU0(JP))**2-1.)&
			*((ZGAMMA1(JP,JLAYER)+1./PMU0(JP))*ZGAMMA4(JP,JLAYER)+ZGAMMA2(JP,JLAYER)*ZGAMMA3(JP,JLAYER)) ! G-
			ZGP(JP,JLAYER)=(PMU0(JP)**2)*PWPRIME(JP,JLAYER)/((ZKESTAR(JP,JLAYER)*PMU0(JP))**2-1.)&
			*((ZGAMMA1(JP,JLAYER)-1./PMU0(JP))*ZGAMMA3(JP,JLAYER)+ZGAMMA2(JP,JLAYER)*ZGAMMA4(JP,JLAYER)) !G+
		END IF
	END DO 
END DO
!
ZTAU_TOT(:,1)=PTAUPRIME(:,1)
DO JLAYER=2,KJNLYR
	DO JP=1,SIZE(PMU0)
		IF (JLAYER<=IATM_LYR(JP)) THEN
		  ! cumulated optical depth
			ZTAU_TOT(JP,JLAYER)=ZTAU_TOT(JP,JLAYER-1)+PTAUPRIME(JP,JLAYER) 
		END IF 
	END DO
END DO
!
!! Matric calculation
ZDM(:,:)=0.
ZD(:,:)=0.
ZDP(:,:)=0.
ZFDIAG(:,:)=0.
!
DO JLAYER=1,KJNLYR-1
  DO JP=1,SIZE(PMU0)
		IF (JLAYER<=IATM_LYR(JP)-1) THEN
			!See matrix documentation page 8 and formal expressions page 9
			ZFDIAG(JP,JLAYER)=EXP(-ZKESTAR(JP,JLAYER)*PTAUPRIME(JP,JLAYER))	    
			!Décalage d'un indice vers la droite par rapport au code python
			ZDM(JP,JLAYER*2)=(1.-ZALBEDO(JP,JLAYER)*ZALBEDO(JP,JLAYER+1))*ZFDIAG(JP,JLAYER)
			ZDM(JP,JLAYER*2+1)=(1./ZALBEDO(JP,JLAYER)-ZALBEDO(JP,JLAYER))*1./ZFDIAG(JP,JLAYER)
		      
			ZD(JP,JLAYER*2)=(1.-ZALBEDO(JP,JLAYER+1)/ZALBEDO(JP,JLAYER))*1./ZFDIAG(JP,JLAYER)
			ZD(JP,JLAYER*2+1)=ZALBEDO(JP,JLAYER)-ZALBEDO(JP,JLAYER+1)
			!Décalage d'un indice vers la gauche par rapport au code python
			ZDP(JP,JLAYER*2)=ZALBEDO(JP,JLAYER+1)*ZALBEDO(JP,JLAYER+1)-1.
			ZDP(JP,JLAYER*2+1)=ZALBEDO(JP,JLAYER)-1./ZALBEDO(JP,JLAYER+1)
		END IF 
	END DO
END DO
!  
ZDP(:,1)=1. !Décalage d'un indice vers la gauche par rapport au code python
ZD(:,1)=1.
!
DO JP=1,SIZE(PMU0)
	ZFDIAG2=EXP(-ZKESTAR(JP,IATM_LYR(JP))*PTAUPRIME(JP,IATM_LYR(JP)))
	!Décalage d'un indice vers la droite par rapport au code python
	ZDM(JP,2*IATM_LYR(JP))=ZFDIAG2*(ZALBEDO(JP,IATM_LYR(JP))-PALB_SOIL(JP))
	ZD(JP,2*IATM_LYR(JP))=1./ZFDIAG2*(1./ZALBEDO(JP,IATM_LYR(JP))-PALB_SOIL(JP))
END DO 
!
!! Vector calculation
ZVECTOR(:,1)=-ZGM(:,1)
DO JLAYER=1,KJNLYR-1
	DO JP=1,SIZE(PMU0) 
		IF (JLAYER<=IATM_LYR(JP)-1)	THEN
		ZDGP=ZGP(JP,JLAYER+1)-ZGP(JP,JLAYER) !doc equation 58
		ZDGM=ZGM(JP,JLAYER+1)-ZGM(JP,JLAYER) !doc equation 58
		!see expression doc page 9
		ZVECTOR(JP,2*JLAYER)=(ZDGM-ZALBEDO(JP,JLAYER+1)*ZDGP)&
		                           *EXP(-ZTAU_TOT(JP,JLAYER)/PMU0(JP))
		ZVECTOR(JP,2*JLAYER+1)=(ZDGP-ZALBEDO(JP,JLAYER)*ZDGM)&
		                           *EXP(-ZTAU_TOT(JP,JLAYER)/PMU0(JP))
		END IF                            
	END DO	              
END DO
!
DO JP=1,SIZE(PMU0)
  ZVECTOR(JP,2*IATM_LYR(JP))=(PALB_SOIL(JP)*(ZGM(JP,IATM_LYR(JP))+PMU0(JP))-&
	      ZGP(JP,IATM_LYR(JP)))*EXP(-ZTAU_TOT(JP,IATM_LYR(JP))/PMU0(JP))
END DO
!
!! Solve the linear system
ZX0_DIR(:,:)=0.
CALL TRIDIAG_GROUND_SNOWCRO(ZDM(:,:),ZD(:,:),ZDP(:,:),ZVECTOR(:,:),ZX0_DIR(:,:),2*IATM_LYR(:),2*MAXVAL(IATM_LYR),0)
!
DO JLAYER=1,KJNLYR
	DO JP=1,SIZE(PMU0)
		  IF (JLAYER<=IATM_LYR(JP)) THEN
			  ZXA(JP,JLAYER)=ZX0_DIR(JP,JLAYER*2-1)
			  ZXB(JP,JLAYER)=ZX0_DIR(JP,JLAYER*2)
		END IF 
	END DO
END DO
!
PZF_TOT(:)=0.
DO JP=1,SIZE(PMU0)
  ZSTAR = ZKESTAR(JP,IATM_LYR(JP))*PTAUPRIME(JP,IATM_LYR(JP))
  ZDEXP = EXP(-ZTAU_TOT(JP,IATM_LYR(JP))/PMU0(JP))
  !
  PZF_TOT(JP)=(ZXA(JP,IATM_LYR(JP))*EXP(-ZSTAR) + &
  ZXB(JP,IATM_LYR(JP))*EXP(ZSTAR) + &
  (ZGM(JP,IATM_LYR(JP))+PMU0(JP))*ZDEXP)*PF0(JP)
END DO
!
IF (LHOOK) CALL DR_HOOK('TWO_STR',0,ZHOOK_HANDLE)
END SUBROUTINE
!
END MODULE MODE_ATMO_TARTES
