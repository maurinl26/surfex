!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source: /home/ducrocq/ANACONV/SOURCES/MYBUG4_4_bug2/DENSITE/regions_2D.f90,v $ $Revision: 1.1 $ $Date: 2003/09/19 11:23:02 $
!-----------------------------------------------------------------
!-----------------------------------------------------------------
! $Source: /home/ducrocq/ANACONV/SOURCES/MYBUG4_4_bug2/DENSITE/regions_2D.f90,v $ $Revision: 1.1 $ $Date: 2003/09/19 11:23:02 $
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!     ###############################
      MODULE MODI_REGIONS_2D
!     ###############################
INTERFACE
!
SUBROUTINE REGIONS_2D(OINAERA,KNUMREG,KNBTREG)
!
LOGICAL, DIMENSION(:,:), INTENT(IN) ::  OINAERA ! LOGICAL to delineate areas (TRUE if in interresting aeras)
INTEGER, DIMENSION(:,:), INTENT(INOUT) :: KNUMREG  ! Numero of  the region for each grid-point
INTEGER, INTENT(INOUT) :: KNBTREG ! Total Number of regions
!
END SUBROUTINE REGIONS_2D
!
END INTERFACE
!
END MODULE MODI_REGIONS_2D
!##############################################
SUBROUTINE REGIONS_2D(OINAERA,KNUMREG,KNBTREG)
!##############################################
!
  !!****  *REGIONS* - routine to determine regions where a logical is true
  !!
  !!    PURPOSE
  !!    -------
  !        The purpose of this routine is to  determine regions of closely  pixels (at least 
  !  one pixel side in common)  where the  logical array  OINAERA is true
  !
  !!**  METHOD
  !!    ------
  !!    based on ini_cart of ASPIC libraries
  !!    REFERENCE
  !!    ---------
  !!
  !!
  !!    AUTHOR
  !!    ------
  !!  	V. Ducrocq   *Meteo France* 
  !!    MODIFICATIONS
  !!    -------------
  !!      Original    3/01/2001
  !-------------------------------------------------------------------------------
  !
  !*       0.    DECLARATIONS
  !              ------------
!
USE YOMHOOK   ,ONLY : LHOOK, DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*       0.1   Declarations of arguments
!
LOGICAL, DIMENSION(:,:), INTENT(IN) ::  OINAERA ! LOGICAL to delineate areas (TRUE if in interresting aeras)
!FO
INTEGER, DIMENSION(:,:), INTENT(INOUT) :: KNUMREG  ! Numero of  the region for each grid-point
INTEGER, INTENT(INOUT) :: KNBTREG ! Total Number of regions
!
!*       0.2   Declarations of local variables
!
INTEGER :: IIMAX, IJMAX ! size of the arrays
INTEGER :: ISUP, IINF
INTEGER :: JILOOP, JJLOOP,JI,JJ ! loop indices
LOGICAL :: GABS 
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-----------------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('REGIONS_2D',0,ZHOOK_HANDLE)
!
IIMAX=SIZE(OINAERA,1)
IJMAX=SIZE(OINAERA,2)
!
KNUMREG(:,:)=0
! Initialization of the first point : I=1, J=1
KNBTREG=0
IF(OINAERA(1,1)) THEN
  KNBTREG=1
  KNUMREG(1,1)=KNBTREG
ENDIF
! 
!* First line  (J=1) 
!------------
DO JILOOP=2,IIMAX
  IF (OINAERA(JILOOP,1) .EQV. OINAERA(JILOOP-1,1)) THEN
    ! We are in the same region as the previous point on the line
    KNUMREG(JILOOP,1) =KNUMREG(JILOOP-1,1) 
  ELSE IF (OINAERA(JILOOP,1)) THEN 
    ! We are not in the same region as the previous point on the line 
    KNBTREG= KNBTREG + 1
    KNUMREG(JILOOP,1) = KNBTREG
  ENDIF   
END DO
!
!* Loop on the following lines
! --------------------
DO JJLOOP= 2, IJMAX
  ! First point of the line 
  ! ----------------
  IF (OINAERA(1,JJLOOP)  .EQV. OINAERA(1,JJLOOP-1) ) THEN
    ! We are in the same region as the previous point in  the column
    KNUMREG( 1,JJLOOP) = KNUMREG(1,JJLOOP-1) 
  ELSE IF (OINAERA(1,JJLOOP) ) THEN 
    ! We are not in the same region as the previous point in the column 
    KNBTREG= KNBTREG + 1
    KNUMREG(1,JJLOOP) = KNBTREG
  ENDIF   
  ! The other points
  ! -------------
  DO JILOOP=2,IIMAX
    ! a. We are not in the same region as the previous point in  the column and as the previous point in the line  : 
    IF ( (OINAERA(JILOOP,JJLOOP)  .NEQV. OINAERA(JILOOP-1,JJLOOP)).AND.   &
           (OINAERA(JILOOP,JJLOOP)  .NEQV. OINAERA(JILOOP,JJLOOP-1)))   THEN
      IF (OINAERA(JILOOP,JJLOOP) ) THEN
        KNBTREG= KNBTREG + 1
        KNUMREG(JILOOP,JJLOOP) = KNBTREG
      ENDIF
      !
      ! b. We are in the same region for  the previous point in  the column and the previous point in the line
      ! but the numbers of the region differ. So we will set the same number :
    ELSE IF ((OINAERA(JILOOP,JJLOOP)  .EQV. OINAERA(JILOOP-1,JJLOOP)).AND.   &
         (OINAERA(JILOOP,JJLOOP)  .EQV. OINAERA(JILOOP,JJLOOP-1)) .AND.           &
         (KNUMREG(JILOOP-1,JJLOOP)  /= KNUMREG(JILOOP,JJLOOP-1) ) ) THEN
 
         ISUP=MAX(KNUMREG(JILOOP,JJLOOP-1), KNUMREG(JILOOP-1,JJLOOP) ) ! the region number to remove
         IINF=MIN(KNUMREG(JILOOP,JJLOOP-1), KNUMREG(JILOOP-1,JJLOOP) ) ! the region number to keep
         KNUMREG(JILOOP,JJLOOP)=IINF
   !  Modify the region number for the previous points in the line 
         DO JI= JILOOP-1,1, -1
          IF (KNUMREG(JI,JJLOOP) == ISUP) THEN
            ! the number region is changed to IINF
            KNUMREG(JI,JJLOOP) =IINF
          ELSE IF  (KNUMREG(JI,JJLOOP) > ISUP)  THEN
            ! we have to lowered by  1, as ISUP region is removed
            KNUMREG(JI,JJLOOP) = KNUMREG(JI,JJLOOP) -1
          END IF
         ENDDO         
    !Modify the region number for  the previous lines 
         DO JJ=JJLOOP-1,1, -1
           GABS=.TRUE.
           DO JI= 1,IIMAX
            IF ( KNUMREG(JI,JJ) == ISUP)  THEN
              GABS=.FALSE.
               ! the number region is changed to IINF
                KNUMREG(JI,JJ) =IINF
            ELSE IF (KNUMREG(JI,JJ) > ISUP)  THEN
               ! we have to lowered by  1, as ISUP region is removed
               KNUMREG(JI,JJ) = KNUMREG(JI,JJ) -1
            ENDIF
           END DO
           ! there is no point in the line associated with region ISUP, So there is  no other point in the previous line  associated with region ISUP
            IF (GABS) EXIT
         ENDDO         
        KNBTREG= KNBTREG -1
   ! c. We are in the same region as  the previous point in  the column :
   ELSE IF  (OINAERA(JILOOP,JJLOOP)  .EQV. OINAERA(JILOOP,JJLOOP-1))   THEN
      KNUMREG( JILOOP,JJLOOP) = KNUMREG(JILOOP,JJLOOP-1) 
   !d. We are in the same region as  the previous point in  the line :
   ELSE 
      KNUMREG( JILOOP,JJLOOP) = KNUMREG(JILOOP-1,JJLOOP) 
  ENDIF
  ENDDO
!
END DO
!
IF (LHOOK) CALL DR_HOOK('REGIONS_2D',1,ZHOOK_HANDLE)
!
END SUBROUTINE REGIONS_2D
