!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!     #########
SUBROUTINE PUT_ALL_VAR_PF(YSC, YSC2, KJ, KPATCH)
!
! ------------------------------------------------------------------------------------------
!!
!!    PUT_ALL_VAR_PF
!!
!!    PURPOSE
!!    -------
!!    Subroutine to put values of all fields of surfex type YSC2 (at a given point)
!!    into YSC (for distributed Particle Filtering purposes)
!!    WORKS only with LSPLIT_PATCH = .FALSE.
!!
!!    METHOD
!!    ------
!!    Brutal
!!
!!    AUTHOR
!!    ------
!!
!!    B. Cluzet bertrand.cluzetatmeteo.fr
!!
!!    MODIFICATION
!!    ------------
!!
!!    Original         03/2018
!!
!----------------------------------------------------------------------------


USE MODD_SURFEX_n, ONLY : SURFEX_t ! B. Cluzet
USE MODD_ISBA_n, ONLY : ISBA_P_t, ISBA_PE_t



IMPLICIT NONE

TYPE(SURFEX_t), INTENT(INOUT) :: YSC ! analyzed type to update (assimilation)
TYPE(SURFEX_t), INTENT(IN)    :: YSC2 ! initial type to read into




INTEGER, INTENT(IN)           :: KJ, KPATCH ! coordinates at which to put the new value

!
!
!* local variables
TYPE(ISBA_P_t), POINTER :: PK, PK2
TYPE(ISBA_PE_t), POINTER :: PEK, PEK2
INTEGER                       :: ISIZE_LMEB_PATCH

! ******************************************************************************************
!--------------- pointer initialisation
PK => YSC%IM%NP%AL(KPATCH)
PEK => YSC%IM%NPE%AL(KPATCH)   

PK2 => YSC2%IM%NP%AL(KPATCH)
PEK2 => YSC2%IM%NPE%AL(KPATCH)   
!-------------------soil prognostic fields (writesurf_isban)
!* soil temperatures (14)
PEK%XTG(KJ,:) = PEK2%XTG(KJ,:)
!
!* soil liquid water contents (14)
PEK%XWG(KJ,:) = PEK2%XWG(KJ,:)
!
!* soil ice water contents (14)
PEK%XWGI(KJ,:) = PEK2%XWGI(KJ,:)
!
!* water intercepted on leaves (1)
PEK%XWR(KJ) = PEK2%XWR(KJ)
!
!* Glacier ice storage (0?)
!
IF (YSC%IM%O%LGLACIER) THEN
  PEK%XICE_STO(KJ) = PEK2%XICE_STO(KJ)
ENDIF
!
!* Leaf Area Index (0?)
!
IF (YSC%IM%O%CPHOTO/='NON'.AND.YSC%IM%O%CPHOTO/='AGS' .AND. YSC%IM%O%CPHOTO/='AST') THEN
  PEK%XLAI(KJ) = PEK2%XLAI(KJ)
ENDIF


!----------------- MEB FIELDS (0 ?)
!
ISIZE_LMEB_PATCH=COUNT(YSC%IM%O%LMEB_PATCH(:))
!
IF (ISIZE_LMEB_PATCH>0) THEN

  !* water intercepted on canopy vegetation leaves
  PEK%XWRL(KJ) = PEK2%XWRL(KJ)

  !* ice on litter
  PEK%XWRLI(KJ) = PEK2%XWRLI(KJ)  
  
  !* snow intercepted on canopy vegetation leaves
  PEK%XWRVN(KJ) = PEK2%XWRVN(KJ)  
  
  !* canopy vegetation temperature
  PEK%XTV(KJ) = PEK2%XTV(KJ)  
  
  !* litter temperature
  PEK%XTL(KJ) = PEK2%XTL(KJ)  

  !* vegetation canopy air temperature
  PEK%XTC(KJ) = PEK2%XTC(KJ)  

  !* vegetation canopy air specific humidity
  PEK%XQC(KJ) = PEK2%XQC(KJ)  

ENDIF
!!!!!!!!!!!!!!!!!!! Semi-prognostic variables (some!)

!* Fraction for each patch
PK%XPATCH(KJ) = PK2%XPATCH(KJ)  

!* patch averaged radiative temperature (K)
YSC%IM%S%XTSRAD_NAT(KJ) = YSC2%IM%S%XTSRAD_NAT(KJ)  

!* aerodynamical resistance
PEK%XRESA(KJ) = PEK2%XRESA(KJ)  

!--------------------- Land use variables (NONE !!!!!!!!!!!!! (not any acces to OLAND_USE in this rutine so far)

!--------------------- canopy levels (cf. writesurf sso_canopyn). (dépend pas du tile/patch)

!* number
! always the same

!* altitudes (6)
YSC%SB%XZ(KJ,:) = YSC2%SB%XZ(KJ,:)  

!* wind in canopy (6)
YSC%SB%XU(KJ,:) = YSC2%SB%XU(KJ,:)  

!* Tke in canopy (6)
YSC%SB%XTKE(KJ,:) = YSC2%SB%XTKE(KJ,:)  

!----------------------- snow variables (writesurf_gr_sn) (8 + nimpur=10 ici)
!* 1. layered
!-------------
!* snow reservoir
PEK%TSNOW%WSNOW(KJ, :) = PEK2%TSNOW%WSNOW(KJ, :) 

!* density
PEK%TSNOW%RHO(KJ, :) = PEK2%TSNOW%RHO(KJ, :) 

!* heat
PEK%TSNOW%HEAT(KJ, :) = PEK2%TSNOW%HEAT(KJ, :) 

!* age
PEK%TSNOW%AGE(KJ, :) = PEK2%TSNOW%AGE(KJ, :)

!* optical diameter
PEK%TSNOW%DIAMOPT(KJ, :) = PEK2%TSNOW%DIAMOPT(KJ, :) 

!* sphericity
PEK%TSNOW%SPHERI(KJ, :) = PEK2%TSNOW%SPHERI(KJ, :) 

!* hist
PEK%TSNOW%HIST(KJ, :) = PEK2%TSNOW%HIST(KJ, :) 

!* impur 1 AND 2
IF (YSC%IM%O%CSNOWRAD == 'T17' ) THEN
  PEK%TSNOW%IMPUR(KJ, :,:) = PEK2%TSNOW%IMPUR(KJ,:,:) 
ENDIF
!* 2. non-layered
!-------------
!*albedo
PEK%TSNOW%ALB(KJ) = PEK2%TSNOW%ALB(KJ) 

!* MEPRA (NewV8.1)(6) : useless (diags !!)
!PEK%TSNOW%DEP_SUP(KJ) = PEK2%TSNOW%DEP_SUP(KJ)
PEK%TSNOW%DEP_TOT(KJ) = PEK2%TSNOW%DEP_TOT(KJ) ! USELESS and FALSE (need to use the XF), they are diags, hence not read !
!PEK%TSNOW%DEP_HUM(KJ) = PEK2%TSNOW%DEP_HUM(KJ)
!PEK%TSNOW%NAT_LEV(KJ) = PEK2%TSNOW%NAT_LEV(KJ)
!PEK%TSNOW%PRO_SUP_TYP(KJ) = PEK2%TSNOW%PRO_SUP_TYP(KJ)
!PEK%TSNOW%AVA_TYP(KJ) = PEK2%TSNOW%AVA_TYP(KJ)


!--------------------diagnostic variables (write_diag_seb_surf_atmn
! only cumulated values are usually read and others are calculated from it 
!so we only need to update the cumulated ones.
! this is useless if LRESETCUMUL=.TRUE. because it will be set to 0 at next time step.
! but it helps ensureing the consistency between all variables of the PREP.

!tile-averaged
!surface fluxes (9)
YSC%DUC%XRN(KJ)   = YSC2%DUC%XRN(KJ) 
YSC%DUC%XH(KJ)    = YSC2%DUC%XH(KJ) 
YSC%DUC%XLE(KJ)   = YSC2%DUC%XLE(KJ) 
YSC%DUC%XLEI(KJ)  = YSC2%DUC%XLEI(KJ) 
YSC%DUC%XGFLUX(KJ)= YSC2%DUC%XGFLUX(KJ) 
!
YSC%DUC%XSWD(KJ) = YSC2%DUC%XSWD(KJ)
YSC%DUC%XSWU(KJ) = YSC2%DUC%XSWU(KJ)
YSC%DUC%XLWD(KJ) = YSC2%DUC%XLWD(KJ)
YSC%DUC%XLWU(KJ) = YSC2%DUC%XLWU(KJ)

! params at surface
YSC%DU%XTS(KJ) = YSC2%DU%XTS(KJ)
!YSC%DU%XDIAG_TRAD(KJ) =  YSC2%DU%XDIAG_TRAD(KJ)
!YSC%DU%XDIAG_EMIS(KJ) =  YSC2%DU%XDIAG_EMIS(KJ)
!YSC%DU%XSFCO2(KJ) =  YSC2%DU%XSFCO2(KJ)



END SUBROUTINE PUT_ALL_VAR_PF
