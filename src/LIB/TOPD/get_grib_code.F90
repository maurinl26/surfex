!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!     #######################
      SUBROUTINE GET_GRIB_CODE(CGRIB_TYPE,CREC_NAME,KNUM_GRIB,KTYPE_GRIB,KLEV1)
!     #######################
!
!!****  *GET_GRIB_CODE*  
!!
!!    PURPOSE
!!    -------
!     This routine aims at reading forcing variables from grib files
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
!!      B. Vincendon    * Meteo-France *
!!
!!    MODIFICATIONS
!!    -------------
!!
!!      Original   11/2007
!-------------------------------------------------------------------------------
!
!*       0.     DECLARATIONS
!               ------------
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
CHARACTER(LEN=*)               , INTENT(IN)  :: CGRIB_TYPE ! Grib type of file to be read
CHARACTER(LEN=5), DIMENSION(10), INTENT(OUT) :: CREC_NAME  ! Grib type of file to be read
INTEGER, DIMENSION(10)         , INTENT(OUT) :: KNUM_GRIB  ! Code of the parameter to get
INTEGER, DIMENSION(10)         , INTENT(OUT) :: KTYPE_GRIB ! Code of type of level
INTEGER, DIMENSION(10)         , INTENT(OUT) :: KLEV1     ! Level
!
!
!*      0.2    declarations of local variables
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('GET_GRIB_CODE',0,ZHOOK_HANDLE)
!
!*       1.0    default values
!               ----------------------
CREC_NAME(1)='T2M'
CREC_NAME(2)='Q2M'
CREC_NAME(3)='RADSW'
CREC_NAME(4)='RADLW'
CREC_NAME(5)='U10M'
CREC_NAME(6)='V10M'
CREC_NAME(7)='PS'
CREC_NAME(8)='RR1'
CREC_NAME(9)='RS1'
CREC_NAME(10)='RG1'

KNUM_GRIB(1)=11
KNUM_GRIB(2)=51
KNUM_GRIB(3)=111
KNUM_GRIB(4)=112
KNUM_GRIB(5)=33
KNUM_GRIB(6)=34
KNUM_GRIB(7)=1
KNUM_GRIB(8)=150
KNUM_GRIB(9)=99
KNUM_GRIB(10)=29

SELECT CASE (CGRIB_TYPE)
  CASE('SAFRAN')
    CREC_NAME(5)='FF10M'
    CREC_NAME(6)=''
    CREC_NAME(7)=''
    !
    KNUM_GRIB(1)=11
    KNUM_GRIB(2)=51
    KNUM_GRIB(3)=116
    KNUM_GRIB(4)=115
    KNUM_GRIB(5)=32
    KNUM_GRIB(6)=0
    KNUM_GRIB(7)=0
    KNUM_GRIB(8)=169
    !
    KTYPE_GRIB(1)=105       
    KTYPE_GRIB(2)=105       
    KTYPE_GRIB(3)=1     
    KTYPE_GRIB(4)=1     
    KTYPE_GRIB(5)=105       
    KTYPE_GRIB(6)=0       
    KTYPE_GRIB(7)=0       
    KTYPE_GRIB(8)=1     
    KTYPE_GRIB(9)=1
    KTYPE_GRIB(10)=1
    !
    KLEV1(1)=2 
    KLEV1(2)=2 
    KLEV1(3)=0 
    KLEV1(4)=0 
    KLEV1(5)=10 
    KLEV1(6)=0 
    KLEV1(7)=0
    KLEV1(8)=0 
    KLEV1(9)=0
    KLEV1(10)=0 
    !
  CASE('AROME ')
    KNUM_GRIB(1)=11
    KNUM_GRIB(2)=51
    KNUM_GRIB(3)=105
    KNUM_GRIB(4)=104
    KNUM_GRIB(5)=33
    KNUM_GRIB(6)=34
    KNUM_GRIB(7)=1
    KNUM_GRIB(8)=150
    KNUM_GRIB(9)=99
    KNUM_GRIB(10)=29
    !
    KTYPE_GRIB(1)=105       
    KTYPE_GRIB(2)=105       
    KTYPE_GRIB(3)=1
    KTYPE_GRIB(4)=1
    KTYPE_GRIB(5)=105       
    KTYPE_GRIB(6)=105       
    KTYPE_GRIB(7)=1       
    KTYPE_GRIB(8)=1   
    KTYPE_GRIB(9)=1
    KTYPE_GRIB(10)=1
    !
    KLEV1(1)=2 
    KLEV1(2)=2 
    KLEV1(3)=0 
    KLEV1(4)=0 
    KLEV1(5)=10 
    KLEV1(6)=10 
    KLEV1(7)=0
    KLEV1(8)=0 
    KLEV1(9)=0
    KLEV1(10)=0 
    !
  CASE('AROMAN')
    KNUM_GRIB(1)=11
    KNUM_GRIB(2)=51
    KNUM_GRIB(3)=105
    KNUM_GRIB(4)=104
    KNUM_GRIB(5)=33
    KNUM_GRIB(6)=34
    KNUM_GRIB(7)=1
    KNUM_GRIB(8)=150
    KNUM_GRIB(9)=99
    KNUM_GRIB(10)=29
    !
    KTYPE_GRIB(1)=105
    KTYPE_GRIB(2)=105
    KTYPE_GRIB(3)=1
    KTYPE_GRIB(4)=1
    KTYPE_GRIB(5)=105
    KTYPE_GRIB(6)=105
    KTYPE_GRIB(7)=1
    KTYPE_GRIB(8)=1
    KTYPE_GRIB(9)=1
    KTYPE_GRIB(10)=1
    !
    KLEV1(1)=2
    KLEV1(2)=2
    KLEV1(3)=0
    KLEV1(4)=0
    KLEV1(5)=10
    KLEV1(6)=10
    KLEV1(7)=0
    KLEV1(8)=0
    KLEV1(9)=0
    KLEV1(10)=0
    !
  CASE('PEAROM')
    KNUM_GRIB(8)=150
    KTYPE_GRIB(1)=105       
    KTYPE_GRIB(2)=105       
    KTYPE_GRIB(3)=1     
    KTYPE_GRIB(4)=1     
    KTYPE_GRIB(5)=105       
    KTYPE_GRIB(6)=105       
    KTYPE_GRIB(7)=1       
    KTYPE_GRIB(8)=1     
    KTYPE_GRIB(9)=1     
    KTYPE_GRIB(10)=1     
    !
    KLEV1(1)=2 
    KLEV1(2)=2 
    KLEV1(3)=0 
    KLEV1(4)=0 
    KLEV1(5)=10 
    KLEV1(6)=10 
    KLEV1(7)=0
    KLEV1(8)=0 
    KLEV1(9)=0
    KLEV1(10)=0 
    !
  CASE('PEAROP')
    KNUM_GRIB(1)=167
    KNUM_GRIB(2)=51 !not read
    KNUM_GRIB(3)=176
    KNUM_GRIB(4)=177
    KNUM_GRIB(5)=165
    KNUM_GRIB(6)=166
    KNUM_GRIB(7)=152
    KNUM_GRIB(8)=84
    KTYPE_GRIB(1)=1       
    KTYPE_GRIB(2)=1
    KTYPE_GRIB(3)=1     
    KTYPE_GRIB(4)=1     
    KTYPE_GRIB(5)=1
    KTYPE_GRIB(6)=1
    KTYPE_GRIB(7)=109       
    KTYPE_GRIB(8)=1     
    KTYPE_GRIB(9)=1     
    KTYPE_GRIB(10)=1     
    KLEV1(1)=0 
    KLEV1(2)=0 
    KLEV1(3)=0 
    KLEV1(4)=0 
    KLEV1(5)=0 
    KLEV1(6)=0 
    KLEV1(7)=1
    KLEV1(8)=0 
    KLEV1(9)=0
    KLEV1(10)=0 
    !
  CASE('ARPEGE')
    KNUM_GRIB(8)=62
    KTYPE_GRIB(1)=111       
    KTYPE_GRIB(2)=111       
    KTYPE_GRIB(3)=111     
    KTYPE_GRIB(4)=111     
    KTYPE_GRIB(5)=111       
    KTYPE_GRIB(6)=111       
    KTYPE_GRIB(7)=111       
    KTYPE_GRIB(8)=111     
    KTYPE_GRIB(9)=111     
    KTYPE_GRIB(10)=111     
    !
    KLEV1(1)=2 
    KLEV1(2)=2 
    KLEV1(3)=0 
    KLEV1(4)=0 
    KLEV1(5)=10 
    KLEV1(6)=10 
    KLEV1(7)=0
    KLEV1(8)=0 
    KLEV1(9)=0
    KLEV1(10)=0 
    !
  CASE('MESONH')
    KNUM_GRIB(8)=62
    KNUM_GRIB(9)=79
    KNUM_GRIB(10)=78
    KTYPE_GRIB(1)=105       
    KTYPE_GRIB(2)=105       
    KTYPE_GRIB(3)=105       
    KTYPE_GRIB(4)=105       
    KTYPE_GRIB(5)=105       
    KTYPE_GRIB(6)=105       
    KTYPE_GRIB(7)=105       
    KTYPE_GRIB(8)=1       
    KTYPE_GRIB(9)=1       
    KTYPE_GRIB(10)=1      
    !
    KLEV1(1)=2 
    KLEV1(2)=2 
    KLEV1(3)=0 
    KLEV1(4)=0 
    KLEV1(5)=10 
    KLEV1(6)=10 
    KLEV1(7)=1  
    KLEV1(8)=0 
    KLEV1(9)=0
    KLEV1(10)=0 
    !
END SELECT
!
IF (LHOOK) CALL DR_HOOK('GET_GRIB_CODE',1,ZHOOK_HANDLE)
!
END SUBROUTINE GET_GRIB_CODE
