MODULE MODE_QUICK_SORT
! -----------------------------------------------------------------------------
!
! Sorting module
!
!
! First version : B. Cluzet, 2019
! Created copying the sorting routines from get_near _meshes_ign.F90, e.g.:
!     """Modifié par Renaud Lestrigant (02/2016) : changement complet de l'algo
!     de recherche des plus proches voisins.
!     Récupération d'un code sur Internet et adaptation locale.
!     (http://jblevins.org/mirror/amiller/qsort.f90)  """
! -----------------------------------------------------------------------------

CONTAINS
!
RECURSIVE SUBROUTINE QUICK_SORT(PLIST, KORDER)

! Quick sort routine from:
! Brainerd, W.S., Goldberg, C.H. & Adams, J.C. (1990) "Programmer's Guide to
! Fortran 90", McGraw-Hill  ISBN 0-07-000248-7, pages 149-150.
! Modified by Alan Miller to include an associated integer array which gives
! the positions of the elements in the original order.
!
IMPLICIT NONE
!
REAL, DIMENSION (:), INTENT(INOUT)  :: PLIST
INTEGER, DIMENSION (:), INTENT(OUT)  :: KORDER
!
! Local variable
INTEGER :: JI

DO JI = 1, SIZE(PLIST)
  KORDER(JI) = JI
END DO

CALL QUICK_SORT_1(1, SIZE(PLIST), PLIST, KORDER)

END SUBROUTINE QUICK_SORT


RECURSIVE SUBROUTINE QUICK_SORT_1(KLEFT_END, KRIGHT_END, PLIST1, KORDER1)

INTEGER, INTENT(IN) :: KLEFT_END, KRIGHT_END
REAL, DIMENSION (:), INTENT(INOUT)  :: PLIST1
INTEGER, DIMENSION (:), INTENT(INOUT)  :: KORDER1
!     Local variables
INTEGER             :: JI, JJ, ITEMP
REAL                :: ZREF, ZTEMP
INTEGER, PARAMETER  :: IMAX_SIMPLE_SORT_SIZE = 6

IF (KRIGHT_END < KLEFT_END + IMAX_SIMPLE_SORT_SIZE) THEN
  ! Use interchange sort for small PLISTs
  CALL INTERCHANGE_SORT(KLEFT_END, KRIGHT_END, PLIST1, KORDER1)
  !
ELSE
  !
  ! Use partition ("quick") sort
  ! valeur au centre du tableau
  ZREF = PLIST1((KLEFT_END + KRIGHT_END)/2)
  JI = KLEFT_END - 1
  JJ = KRIGHT_END + 1

  DO
    ! Scan PLIST from left end until element >= ZREF is found
    DO
      JI = JI + 1
      IF (PLIST1(JI) >= ZREF) EXIT
    END DO
    ! Scan PLIST from right end until element <= ZREF is found
    DO
      JJ = JJ - 1
      IF (PLIST1(JJ) <= ZREF) EXIT
    END DO


    IF (JI < JJ) THEN
      ! Swap two out-of-order elements
      ZTEMP = PLIST1(JI)
      PLIST1(JI) = PLIST1(JJ)
      PLIST1(JJ) = ZTEMP
      ITEMP = KORDER1(JI)
      KORDER1(JI) = KORDER1(JJ)
      KORDER1(JJ) = ITEMP
    ELSE IF (JI == JJ) THEN
      JI = JI + 1
      EXIT
    ELSE
      EXIT
    END IF
  END DO

  IF (KLEFT_END < JJ) CALL QUICK_SORT_1(KLEFT_END, JJ, PLIST1, KORDER1)
  IF (JI < KRIGHT_END) CALL QUICK_SORT_1(JI, KRIGHT_END,PLIST1,KORDER1)
END IF

END SUBROUTINE QUICK_SORT_1


SUBROUTINE INTERCHANGE_SORT(KLEFT_END, KRIGHT_END, PLIST2, KORDER2)

INTEGER, INTENT(IN) :: KLEFT_END, KRIGHT_END
REAL, DIMENSION (:), INTENT(INOUT)  :: PLIST2
INTEGER, DIMENSION (:), INTENT(INOUT)  :: KORDER2
!     Local variables
INTEGER             :: JI, JJ, ITEMP
REAL                :: ZTEMP

! boucle sur tous les points
DO JI = KLEFT_END, KRIGHT_END - 1
  !
  ! boucle sur les points suivants le point JI
  DO JJ = JI+1, KRIGHT_END
    !
    ! si la distance de JI au point est plus grande que celle de JJ
    IF (PLIST2(JI) > PLIST2(JJ)) THEN
      ! distance de JI au point (la plus grande)
      ZTEMP = PLIST2(JI)
      ! le point JJ est déplacé à l'indice JI dans le tableau 
      PLIST2(JI) = PLIST2(JJ)
      ! le point JI est déplacé à l'indice JJ dans le tableau
      PLIST2(JJ) = ZTEMP
      ! indice du point JI dans le tableau
      ITEMP = KORDER2(JI)
      ! l'indice du point JJ est mis à la place JI
      KORDER2(JI) = KORDER2(JJ)
      ! l'indice du point JI est mis à la place JJ
      KORDER2(JJ) = ITEMP
    END IF
    !
  END DO
  !
END DO

END SUBROUTINE INTERCHANGE_SORT
!
!-------------------------------------------------------------------------------
!
END MODULE MODE_QUICK_SORT
