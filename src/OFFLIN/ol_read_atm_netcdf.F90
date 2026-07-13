!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!     #########
SUBROUTINE OL_READ_ATM_NETCDF (HSURF_FILETYPE,                            &
                               PTA,PQA,PWIND,PDIR_SW,PSCA_SW,PLW,PSNOW,   &
                               PRAIN,PPS,PCO2,PIMPWET,PIMPDRY,PO3,PAE,PDIR )  
!**************************************************************************
!
!!    PURPOSE
!!    -------
!         Read in the netcdf file the atmospheric forcing for the actual time
!         step KFORC_STEP, and for the next one.
!         The two time step are needed for the time interpolation of the
!         forcing.
!         If the end of the file  is reached, set the two step to the last
!         values.
!         Return undef value if the variable is not present
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
!!      Original     06/2003
!!      P. Le Moigne 10/2004: set INB to 2 because of revised temporal loop in offline.f90:
!!                            time evolution is done at the end of isba time step so first 
!!                            isba computation is done on first forcing time step
!!      P. Le Moigne 10/2005: consistency checking between orographies read from forcing 
!!                            file and from initial file

!          
!
!
!
!
USE MODN_IO_OFFLINE,  ONLY : NIMPUROF,LFORCIMP,LFORCATMOTARTES
USE MODD_IO_SURF_OL, ONLY : XCOUNT
USE MODI_READ_SURF
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!
! global variables
!
!
REAL, DIMENSION(:,:),INTENT(OUT) :: PTA
REAL, DIMENSION(:,:),INTENT(OUT) :: PQA
REAL, DIMENSION(:,:),INTENT(OUT) :: PWIND
REAL, DIMENSION(:,:),INTENT(OUT) :: PDIR_SW
REAL, DIMENSION(:,:),INTENT(OUT) :: PSCA_SW
REAL, DIMENSION(:,:),INTENT(OUT) :: PLW
REAL, DIMENSION(:,:),INTENT(OUT) :: PSNOW
REAL, DIMENSION(:,:),INTENT(OUT) :: PRAIN
REAL, DIMENSION(:,:),INTENT(OUT) :: PPS
REAL, DIMENSION(:,:),INTENT(OUT) :: PCO2
REAL, DIMENSION(:,:),INTENT(OUT) :: PO3
REAL, DIMENSION(:,:),INTENT(OUT) :: PAE
REAL, DIMENSION(:,:,:),INTENT(OUT) :: PIMPWET
REAL, DIMENSION(:,:,:),INTENT(OUT) :: PIMPDRY
REAL, DIMENSION(:,:),INTENT(OUT) :: PDIR
 CHARACTER(LEN=6)    ,INTENT(IN)  :: HSURF_FILETYPE

! local variables
INTEGER                          :: IRET
REAL(KIND=JPRB) :: ZHOOK_HANDLE
INTEGER  :: JIMP
 CHARACTER(LEN=16)  :: NAMEWET
 CHARACTER(LEN=16)  :: NAMEDRY
!
IF (LHOOK) CALL DR_HOOK('OL_READ_ATM_NETCDF',0,ZHOOK_HANDLE)
 CALL READ_SURF(&
                'OFFLIN','Tair',      PTA    (:,1:XCOUNT),IRET)
 CALL READ_SURF(&
                'OFFLIN','Qair',      PQA    (:,1:XCOUNT),IRET)
 CALL READ_SURF(&
                'OFFLIN','Wind',      PWIND  (:,1:XCOUNT),IRET)
 CALL READ_SURF(&
                'OFFLIN','LWdown',    PLW    (:,1:XCOUNT),IRET)
 CALL READ_SURF(&
                'OFFLIN','DIR_SWdown',PDIR_SW(:,1:XCOUNT),IRET)
 CALL READ_SURF(&
                'OFFLIN','SCA_SWdown',PSCA_SW(:,1:XCOUNT),IRET)
 CALL READ_SURF(&
                'OFFLIN','Rainf',     PRAIN  (:,1:XCOUNT),IRET)
 CALL READ_SURF(&
                'OFFLIN','Snowf',     PSNOW  (:,1:XCOUNT),IRET)
 CALL READ_SURF(&
                'OFFLIN','PSurf',     PPS    (:,1:XCOUNT),IRET)
 CALL READ_SURF(&
                'OFFLIN','Wind_DIR',  PDIR   (:,1:XCOUNT),IRET)
 CALL READ_SURF(&
                'OFFLIN','CO2air',    PCO2   (:,1:XCOUNT),IRET)
         
 IF (LFORCATMOTARTES) THEN
 !Read the O3 content and the total optical depth
   CALL READ_SURF(&
                      'OFFLIN','AODTOT',  PAE  (:,1:XCOUNT),IRET)
   CALL READ_SURF(&
                      'OFFLIN','OZONE',  PO3  (:,1:XCOUNT),IRET)
 ENDIF
 IF (LFORCIMP) THEN
 !Read the Impurity content in the forcing File, you have to set the name of netcdf variables to IMPWET1,IMPWET2... for wet deposit 
 !and IMPDRY1,IMPDRY2... for dry deposit. The fluxes in mocages are in g/m²/s.    
   DO JIMP=1,NIMPUROF
    WRITE(NAMEWET,'(A6,I1)') 'IMPWET',JIMP
    WRITE(NAMEDRY,'(A6,I1)') 'IMPDRY',JIMP
     CALL READ_SURF(&
                      'OFFLIN',NAMEWET,  PIMPWET  (:,JIMP,1:XCOUNT),IRET)
     CALL READ_SURF(&
                      'OFFLIN',NAMEDRY,  PIMPDRY   (:,JIMP,1:XCOUNT),IRET) 
   ENDDO                   
 ENDIF  
                    
IF (LHOOK) CALL DR_HOOK('OL_READ_ATM_NETCDF',1,ZHOOK_HANDLE)

END SUBROUTINE OL_READ_ATM_NETCDF
