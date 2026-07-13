!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!-------------------------------------------------------------------------------
!     ##################
      MODULE MODD_TOPODYN
!     ##################
!
!!****  *MODD_TOPODYN - declaration of variables used by Topodyn
!!
!!    PURPOSE
!!    -------
!
!!
!!**  IMPLICIT ARGUMENTS
!!    ------------------
!!      None 
!!
!!    REFERENCE
!!    ---------
!!
!!    AUTHOR
!!    ------
!!     F. Habets and K. Chancibault
!!
!!    MODIFICATIONS
!!    -------------
!!      Original       29/09/03
!!      BV: modifications  2006: division in two part (some variables are
!                            now in modd_coupling_topo_n    
!!      BV: modifications  04/2007: addition of XTOPD_STEP and NNB_TOPD_STEP
!
!*       0.   DECLARATIONS
!             ------------
!
USE MODD_TOPD_PAR, ONLY : JPCAT
!
IMPLICIT NONE
!
!-------------------------------------------------------------------------------
! Variables specific to Topodyn
!
 CHARACTER(LEN=15), DIMENSION(JPCAT) :: CCAT     ! base name for topographic files
INTEGER                             :: NNCAT    ! catchments number
!
INTEGER                             :: NNB_TOPD_STEP   ! number of TOPODYN time steps
REAL                                :: XTOPD_STEP      ! TOPODYN time step
!
INTEGER                             :: NMESHT   ! maximal number of catchments meshes

REAL, ALLOCATABLE, DIMENSION(:,:)   :: XDMAXT   ! maximal deficit on TOPODYN grid (m)
REAL, ALLOCATABLE, DIMENSION(:)     :: XDXT     ! catchment grid mesh size (m)
REAL, ALLOCATABLE, DIMENSION(:)     :: XMPARA   ! M parameter on TOPODYN grid (m)
REAL, ALLOCATABLE, DIMENSION(:)     :: XA_SPEED ! A parameter for river speed computation
REAL, ALLOCATABLE, DIMENSION(:)     :: XB_SPEED ! B parameter for river speed computation
REAL, ALLOCATABLE, DIMENSION(:,:)   :: XALERT    ! ALERT LEVEL for each point of the watershed, each time step

INTEGER, ALLOCATABLE, DIMENSION(:)  :: NNMC     ! catchments pixels number
REAL, ALLOCATABLE, DIMENSION(:,:,:) :: XCONN    ! pixels reference number and 
                                                ! connections between
INTEGER, ALLOCATABLE, DIMENSION(:,:):: NLINE    ! second index of the pixel in the array 
                                                ! XCONN
REAL, ALLOCATABLE, DIMENSION(:,:)   :: XTANB    ! pixels topographic slope (Tan(Beta))
REAL, ALLOCATABLE, DIMENSION(:,:)   :: XSLOP    ! pixels topographic slope/length flow

!Variables à priori inutiles
REAL, ALLOCATABLE, DIMENSION(:,:)   :: XDAREA   ! drainage area (aire drainee)

! Variables defining the catchments

INTEGER, ALLOCATABLE, DIMENSION(:)  :: NNXC     ! number of topographic grid points on 
                                                ! abscissa axis
INTEGER, ALLOCATABLE, DIMENSION(:)  :: NNYC     ! number of topographic grid points on ordinate 
                                                ! axis
INTEGER, ALLOCATABLE, DIMENSION(:)  :: NNPT     ! number of pixels in the topographic 
                                                ! domain
INTEGER                             :: NPMAX    ! maximal number of pixels in the 
                                                ! topographic grid

REAL, ALLOCATABLE, DIMENSION(:)     :: XX0,XY0  ! coordinates bottom-left pixel of each 
                                                ! topographic domain

REAL, ALLOCATABLE, DIMENSION(:)     :: XNUL     ! undefined value in topographic files

REAL, ALLOCATABLE, DIMENSION(:,:)   :: XTOPD    ! topographic values in topographic files
REAL, DIMENSION(JPCAT)              :: XRTOP_D2 ! depth used by topodyn for lateral transfers
                                                ! (expressed in ratio of isba d2)
                                                !
! Variables used in routing module 
INTEGER, ALLOCATABLE, DIMENSION(:)  :: NNISO    ! number of time step for the isochrones
REAL, ALLOCATABLE, DIMENSION(:,:)   :: XCISO    ! isochrones routing constants 

REAL, DIMENSION(JPCAT)              :: XQINIT   ! Initial discharge at the outlet of the catchments
REAL, ALLOCATABLE, DIMENSION(:,:)   :: XQTOT    ! Total discharge at the outlet of the catchments

REAL, DIMENSION(JPCAT)              :: XSPEEDR,XSPEEDH ! River and hillslope speed
REAL, DIMENSION(JPCAT)              :: XSPEEDG         ! Ground speed
REAL, DIMENSION(JPCAT)              :: XMAX_SPEED         ! Ground speed
LOGICAL                             :: LSPEEDR_VAR  ! T to modulate river speed according to the discharge
REAL, ALLOCATABLE, DIMENSION(:,:)   :: XDRIV, XDHIL    ! River and hillslope distances
REAL, ALLOCATABLE, DIMENSION(:,:)   :: XDGRD           ! Ground distance
REAL, ALLOCATABLE, DIMENSION(:,:)   :: XTIME_TOPD      ! Time to go to the outlet
                                                       ! at the soil surface
REAL, ALLOCATABLE, DIMENSION(:,:)   :: XTIME_TOPD_DRAIN! Time to go to the outlet in the ground

INTEGER, ALLOCATABLE, DIMENSION(:)  :: NX_STEP_ROUT   ! number of maximal time step to join the outlet of 
                                                ! any catchment

! Variables used in exfiltration module 
REAL, ALLOCATABLE, DIMENSION(:,:)   :: XLAMBDA  ! pure topographic index
REAL, ALLOCATABLE, DIMENSION(:,:)   :: XQB_DR
REAL, ALLOCATABLE, DIMENSION(:,:)   :: XQB_RUN
! for topodyn alone
REAL, ALLOCATABLE, DIMENSION(:)   :: XRI,XRI_PREV! recharge on ISBA grid
REAL, ALLOCATABLE, DIMENSION(:)   :: XSRFULL! reservoir of interception for
!TOPODYN only
REAL, ALLOCATABLE, DIMENSION(:,:) :: XDEFT! pixel deficit
!
! Variables to introduce sub catchments
INTEGER                             :: NNCAT_MAX    ! maximum number of subcatchments
INTEGER, ALLOCATABLE, DIMENSION(:)  :: NNCAT_SUB    ! number of subcatchments for each catchment
                                                    ! maximum value 99
INTEGER, ALLOCATABLE, DIMENSION(:,:):: NBPIX_SUB! For sub catchments, number of pixels in the sub-catchment
INTEGER, ALLOCATABLE, DIMENSION(:,:):: NRIV_JUNCTION! For sub catchments, closest river junction pixel 
INTEGER, ALLOCATABLE, DIMENSION(:,:):: NCAT_CAT_TO_SUB !Mask giving the numero of the sub catchment for each pixel
                                                       !of a Catchment  
INTEGER, ALLOCATABLE, DIMENSION(:,:):: NPIX_CAT_TO_SUB !Mask  giving the numero of the pixel for each pixel
                                                       !of a Catchment 
INTEGER, ALLOCATABLE, DIMENSION(:,:,:):: NPIX_SUB_TO_CAT !Mask  giving the numero of the  pixel of Catchment
                                                       !for each pixel of sub-catchment
INTEGER, ALLOCATABLE, DIMENSION(:,:):: NMASK_OUTLET ! Mask of sub catchments outlet
REAL   , ALLOCATABLE, DIMENSION(:,:,:):: XDIST_OUTLET ! distance between catchment outlet ans sub cat outlets
REAL   , ALLOCATABLE, DIMENSION(:,:):: XDRIV_SUBCAT ! River distances for sub catchments
REAL, ALLOCATABLE, DIMENSION(:,:,:)   :: XQTOT_SUB! Total discharge at the outlet of sub-catchments
REAL, ALLOCATABLE, DIMENSION(:,:,:)   :: XQB_DR_SUB
REAL, ALLOCATABLE, DIMENSION(:,:,:)   :: XQB_RUN_SUB
REAL, ALLOCATABLE, DIMENSION(:,:)   :: XTIME_TOPD_SUB      ! Time to go to any river point
                                                       ! at the soil surface
REAL, ALLOCATABLE, DIMENSION(:,:)   :: XTIME_TOPD_DRAIN_SUB! Time to go to any river point in the ground

!-------------------------------------------------------------------------------------
!
END MODULE MODD_TOPODYN

