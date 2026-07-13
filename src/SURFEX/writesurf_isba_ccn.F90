!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!     #########
      SUBROUTINE WRITESURF_ISBA_CC_n (HSELECT, IO, S, NP, NPE, KI, HPROGRAM)
!     #####################################
!
!!****  *WRITESURF_ISBA_CC_n* - writes ISBA Carbon Cycle prognostic fields
!!                        
!!
!!    PURPOSE
!!    -------
!!
!!**  METHOD
!!    ------
!!
!!    EXTERNAL
!!    --------
!!
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!
!!    REFERENCE
!!    ---------
!!
!!
!!    AUTHOR
!!    ------
!!      B. Decharme   *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    12/2023 Split from previous writesurf_isban.F90 routine
!!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_ISBA_OPTIONS_n, ONLY : ISBA_OPTIONS_t
USE MODD_ISBA_n,         ONLY : ISBA_NP_t, ISBA_NPE_t, ISBA_S_t
!
USE MODD_SURF_PAR,       ONLY : LEN_HREC
!
USE MODI_WRITE_FIELD_1D_PATCH
USE MODI_WRITE_SURF
!
USE YOMHOOK,            ONLY : LHOOK,   DR_HOOK
USE PARKIND1,           ONLY : JPRB
!
IMPLICIT NONE
!
!*       0.1   Declarations of arguments
!              -------------------------
!
CHARACTER(LEN=*), DIMENSION(:), INTENT(IN) :: HSELECT 
!
TYPE(ISBA_OPTIONS_t), INTENT(INOUT) :: IO
TYPE(ISBA_S_t),       INTENT(INOUT) :: S
TYPE(ISBA_NP_t),      INTENT(INOUT) :: NP
TYPE(ISBA_NPE_t),     INTENT(INOUT) :: NPE
INTEGER,              INTENT(IN)    :: KI
!
CHARACTER(LEN=6),    INTENT(IN)    :: HPROGRAM ! program calling
!
!*       0.2   Declarations of local variables
!              -------------------------------
!
!
CHARACTER(LEN=LEN_HREC) :: YRECFM         ! Name of the article to be read
!
CHARACTER(LEN=4 ) :: YLVL
CHARACTER(LEN=3 ) :: YVAR
CHARACTER(LEN=100):: YCOMMENT       ! Comment string
CHARACTER(LEN=25) :: YFORM          ! Writing format
!
INTEGER :: IRESP                    ! IRESP  : return-code if a problem appears
INTEGER :: JL, JP, JNL, JNC, JNLV   ! loop counter on levels
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('WRITESURF_ISBA_CC_N',0,ZHOOK_HANDLE)
!
!------------------------------------------------------------------------------
!
YRECFM = 'RESPSL'
YCOMMENT=YRECFM
CALL WRITE_SURF(HSELECT,HPROGRAM,YRECFM,IO%CRESPSL,IRESP,HCOMMENT=YCOMMENT)
!
YRECFM = 'SOILGAS'
YCOMMENT=YRECFM
CALL WRITE_SURF(HSELECT,HPROGRAM,YRECFM,IO%LSOILGAS,IRESP,HCOMMENT=YCOMMENT)
!
IF(IO%LSPINUPCARBS)THEN
  YRECFM='NBYEARSOLD'
  YCOMMENT='yrs'
  CALL WRITE_SURF(HSELECT,HPROGRAM,YRECFM,IO%NNBYEARSOLD,IRESP,HCOMMENT=YCOMMENT)
ENDIF
!
!------------------------------------------------------------------------------
!
IF(IO%CRESPSL=='CNT')THEN
  !
  !* Bulk Soil carbon
  !
  DO JNL=1,IO%NNLITTER
    DO JNLV=1,IO%NNLITTLEVS
      WRITE(YLVL,'(I1,A1,I1)') JNL,'_',JNLV
      YRECFM='LITTER'//ADJUSTL(YLVL(:LEN_TRIM(YLVL)))
      YFORM='(A10,I1.1,A1,I1.1,A8)'
      WRITE(YCOMMENT,FMT=YFORM) 'X_Y_LITTER',JNL,' ',JNLV,' (gC/m2)'
      DO JP = 1,IO%NPATCH
        CALL WRITE_FIELD_1D_PATCH(HSELECT,HPROGRAM,YRECFM,YCOMMENT,JP,&
                        NP%AL(JP)%NR_P,NPE%AL(JP)%XLITTER(:,JNL,JNLV),KI,S%XWORK_WR)    
      ENDDO        
    END DO
  END DO
  !
  DO JNC=1,IO%NNSOILCARB
    WRITE(YLVL,'(I4)') JNC
    YRECFM='SOILCARB'//ADJUSTL(YLVL(:LEN_TRIM(YLVL)))
    YFORM='(A8,I1.1,A8)'
    WRITE(YCOMMENT,FMT=YFORM) 'X_Y_SOILCARB',JNC,' (gC/m2)'
    DO JP = 1,IO%NPATCH
      CALL WRITE_FIELD_1D_PATCH(HSELECT,HPROGRAM,YRECFM,YCOMMENT,JP,&
                      NP%AL(JP)%NR_P,NPE%AL(JP)%XSOILCARB(:,JNC),KI,S%XWORK_WR)    
    ENDDO     
  END DO
  !
  DO JNLV=1,IO%NNLITTLEVS
    WRITE(YLVL,'(I4)') JNLV
    YRECFM='LIGN_STR'//ADJUSTL(YLVL(:LEN_TRIM(YLVL)))
    YFORM='(A12,I1.1,A8)'
    WRITE(YCOMMENT,FMT=YFORM) 'X_Y_LIGNIN_STRUC',JNLV,' (-)'
    DO JP = 1,IO%NPATCH
      CALL WRITE_FIELD_1D_PATCH(HSELECT,HPROGRAM,YRECFM,YCOMMENT,JP,&
                      NP%AL(JP)%NR_P,NPE%AL(JP)%XLIGNIN_STRUC(:,JNLV),KI,S%XWORK_WR)    
    ENDDO       
  END DO
  !
  !----------------------------------------------------------------------------
  !
  ELSEIF(IO%CRESPSL=='DIF')THEN
  !
  !
  !* Multi-layer Soil carbon
  !
  !
  YRECFM='SFLIGN'
  YCOMMENT='X_Y_SURFACE_LIGNIN_STRUC'
  DO JP = 1,IO%NPATCH
      CALL WRITE_FIELD_1D_PATCH(HSELECT,HPROGRAM,YRECFM,YCOMMENT,JP,&
                      NP%AL(JP)%NR_P,NPE%AL(JP)%XSURFACE_LIGNIN_STRUC(:),KI,S%XWORK_WR)    
  ENDDO
  !
  DO JNL=1,IO%NNLITTER
     WRITE(YLVL,'(I1)') JNL
     YRECFM='SFLIT'//ADJUSTL(YLVL(:LEN_TRIM(YLVL)))
     YCOMMENT='X_Y_SURFACE_LITTER'//ADJUSTL(YLVL(:LEN_TRIM(YLVL)))//' (gC/m2)'
     DO JP = 1,IO%NPATCH
        CALL WRITE_FIELD_1D_PATCH(HSELECT,HPROGRAM,YRECFM,YCOMMENT,JP,&
                      NP%AL(JP)%NR_P,NPE%AL(JP)%XSURFACE_LITTER(:,JNL),KI,S%XWORK_WR)    
     ENDDO
  END DO
  !
  DO JL=1,IO%NGROUND_LAYER
     WRITE(YLVL,'(I2)') JL
     YRECFM='DFLIGN'//ADJUSTL(YLVL(:LEN_TRIM(YLVL)))
     YCOMMENT='X_Y_Z_SOIL_LIGNIN_STRUC'//ADJUSTL(YLVL(:LEN_TRIM(YLVL)))
     DO JP = 1,IO%NPATCH
        CALL WRITE_FIELD_1D_PATCH(HSELECT,HPROGRAM,YRECFM,YCOMMENT,JP,&
                      NP%AL(JP)%NR_P,NPE%AL(JP)%XSOILDIF_LIGNIN_STRUC(:,JL),KI,S%XWORK_WR)    
     ENDDO    
  END DO
  !
  DO JL=1,IO%NGROUND_LAYER
     WRITE(YLVL,'(I2)') JL
     YRECFM='DFLIT1L'//ADJUSTL(YLVL(:LEN_TRIM(YLVL)))
     YCOMMENT='X_Y_Z_SOIL_LITTER1'//ADJUSTL(YLVL(:LEN_TRIM(YLVL)))//' (gC/m2)'
     DO JP = 1,IO%NPATCH
        CALL WRITE_FIELD_1D_PATCH(HSELECT,HPROGRAM,YRECFM,YCOMMENT,JP,&
                      NP%AL(JP)%NR_P,NPE%AL(JP)%XSOILDIF_LITTER(:,JL,1),KI,S%XWORK_WR)    
     ENDDO      
  END DO
  DO JL=1,IO%NGROUND_LAYER
     WRITE(YLVL,'(I2)') JL
     YRECFM='DFLIT2L'//ADJUSTL(YLVL(:LEN_TRIM(YLVL)))
     YCOMMENT='X_Y_Z_SOIL_LITTER2'//ADJUSTL(YLVL(:LEN_TRIM(YLVL)))//' (gC/m2)'
     DO JP = 1,IO%NPATCH
        CALL WRITE_FIELD_1D_PATCH(HSELECT,HPROGRAM,YRECFM,YCOMMENT,JP,&
                        NP%AL(JP)%NR_P,NPE%AL(JP)%XSOILDIF_LITTER(:,JL,2),KI,S%XWORK_WR)    
     ENDDO      
  END DO
  !
  DO JL=1,IO%NGROUND_LAYER
     WRITE(YLVL,'(I2)')JL
     YRECFM='DFSOC1L'//ADJUSTL(YLVL(:LEN_TRIM(YLVL)))
     YCOMMENT='X_Y_Z_SOILCARB1'//ADJUSTL(YLVL(:LEN_TRIM(YLVL)))//' (gC/m2)'
     DO JP = 1,IO%NPATCH
        CALL WRITE_FIELD_1D_PATCH(HSELECT,HPROGRAM,YRECFM,YCOMMENT,JP,&
                          NP%AL(JP)%NR_P,NPE%AL(JP)%XSOILDIF_CARB(:,JL,1),KI,S%XWORK_WR)
     ENDDO        
  ENDDO
  DO JL=1,IO%NGROUND_LAYER
         WRITE(YLVL,'(I2)')JL
         YRECFM='DFSOC2L'//ADJUSTL(YLVL(:LEN_TRIM(YLVL)))
         YCOMMENT='X_Y_Z_SOILCARB2'//ADJUSTL(YLVL(:LEN_TRIM(YLVL)))//' (gC/m2)'
         DO JP = 1,IO%NPATCH
            CALL WRITE_FIELD_1D_PATCH(HSELECT,HPROGRAM,YRECFM,YCOMMENT,JP,&
                          NP%AL(JP)%NR_P,NPE%AL(JP)%XSOILDIF_CARB(:,JL,2),KI,S%XWORK_WR)    
         ENDDO        
  ENDDO
  DO JL=1,IO%NGROUND_LAYER
         WRITE(YLVL,'(I2)')JL
         YRECFM='DFSOC3L'//ADJUSTL(YLVL(:LEN_TRIM(YLVL)))
         YCOMMENT='X_Y_Z_SOILCARB3'//ADJUSTL(YLVL(:LEN_TRIM(YLVL)))//' (gC/m2)'
         DO JP = 1,IO%NPATCH
            CALL WRITE_FIELD_1D_PATCH(HSELECT,HPROGRAM,YRECFM,YCOMMENT,JP,&
                          NP%AL(JP)%NR_P,NPE%AL(JP)%XSOILDIF_CARB(:,JL,3),KI,S%XWORK_WR)    
         ENDDO        
  ENDDO
  !
ENDIF
!
!
!* Multi-layer Soil gas
!
IF(IO%CRESPSL=='DIF'.AND.IO%LSOILGAS)THEN
  !
  DO JL=1,IO%NGROUND_LAYER
     WRITE(YLVL,'(I2)') JL
     YRECFM='GASO2L'//ADJUSTL(YLVL(:LEN_TRIM(YLVL)))
     YCOMMENT='X_Y_Z_SGASO2L'//ADJUSTL(YLVL(:LEN_TRIM(YLVL)))//' (g/m3)'
     DO JP = 1,IO%NPATCH
        CALL WRITE_FIELD_1D_PATCH(HSELECT,HPROGRAM,YRECFM,YCOMMENT,JP,&
                      NP%AL(JP)%NR_P,NPE%AL(JP)%XSGASO2(:,JL),KI,S%XWORK_WR)    
     ENDDO    
  END DO
  !
  DO JL=1,IO%NGROUND_LAYER
     WRITE(YLVL,'(I2)') JL
     YRECFM='GASCO2L'//ADJUSTL(YLVL(:LEN_TRIM(YLVL)))
     YCOMMENT='X_Y_Z_SGASCO2L'//ADJUSTL(YLVL(:LEN_TRIM(YLVL)))//' (g/m3)'
     DO JP = 1,IO%NPATCH
        CALL WRITE_FIELD_1D_PATCH(HSELECT,HPROGRAM,YRECFM,YCOMMENT,JP,&
                      NP%AL(JP)%NR_P,NPE%AL(JP)%XSGASCO2(:,JL),KI,S%XWORK_WR)    
     ENDDO
  END DO
  !
  DO JL=1,IO%NGROUND_LAYER
     WRITE(YLVL,'(I2)') JL
     YRECFM='GASCH4L'//ADJUSTL(YLVL(:LEN_TRIM(YLVL)))
     YCOMMENT='X_Y_Z_SGASCH4L'//ADJUSTL(YLVL(:LEN_TRIM(YLVL)))//' (g/m3)'
     DO JP = 1,IO%NPATCH
        CALL WRITE_FIELD_1D_PATCH(HSELECT,HPROGRAM,YRECFM,YCOMMENT,JP,&
                      NP%AL(JP)%NR_P,NPE%AL(JP)%XSGASCH4(:,JL),KI,S%XWORK_WR)    
     ENDDO
  END DO
  !
ENDIF
!
!
!* Fire scheme
!
!
IF(IO%LFIRE)THEN
  !
  YRECFM = 'FIREIND'
  YCOMMENT=YRECFM
  DO JP = 1,IO%NPATCH
     CALL WRITE_FIELD_1D_PATCH(HSELECT,HPROGRAM,YRECFM,YCOMMENT,JP,&
                   NP%AL(JP)%NR_P,NPE%AL(JP)%XFIREIND,KI,S%XWORK_WR)    
  ENDDO
  !
  YRECFM='MOISTLITFIRE'
  YCOMMENT=YRECFM
  DO JP = 1,IO%NPATCH
     CALL WRITE_FIELD_1D_PATCH(HSELECT,HPROGRAM,YRECFM,YCOMMENT,JP,&
                   NP%AL(JP)%NR_P,NPE%AL(JP)%XMOISTLIT_FIRE,KI,S%XWORK_WR)    
  ENDDO
  !
  YRECFM='TEMPLITFIRE'
  YCOMMENT=YRECFM
  DO JP = 1,IO%NPATCH
     CALL WRITE_FIELD_1D_PATCH(HSELECT,HPROGRAM,YRECFM,YCOMMENT,JP,&
                   NP%AL(JP)%NR_P,NPE%AL(JP)%XTEMPLIT_FIRE,KI,S%XWORK_WR)    
  ENDDO
  !
END IF
!
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('WRITESURF_ISBA_CC_N',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE WRITESURF_ISBA_CC_n
