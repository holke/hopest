#include "hopest_f.h"
MODULE MODH_Refine_Vars
!===================================================================================================================================
! Contains global variables provided by the mesh routines
!===================================================================================================================================
! MODULES
USE,INTRINSIC :: ISO_C_BINDING
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PUBLIC
SAVE
!-----------------------------------------------------------------------------------------------------------------------------------
INTEGER                          :: Nsuper
INTEGER                          :: refineLevel
INTEGER                          :: refineGeomType
INTEGER                          :: refineBCIndex
REAL,ALLOCATABLE                 :: xi_Nsuper(:) 
INTEGER,ALLOCATABLE              :: RefineList(:)
REAL,ALLOCATABLE                 :: refineBoundary(:)
INTEGER,ALLOCATABLE              :: TreeSidesToRefine(:,:)
REAL                             :: sphereCenter(3),sphereRadius
REAL                             :: shellCenter(3),shellRadius_inner,shellRadius_outer
REAL                             :: boxBoundary(6)
!-----------------------------------------------------------------------------------------------------------------------------------

END MODULE MODH_Refine_Vars
