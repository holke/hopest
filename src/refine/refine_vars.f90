#include "hopest_f.h"
MODULE MOD_Refine_Vars
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
INTEGER                          :: refineLevel
INTEGER                          :: refineType
INTEGER                          :: refineGeomType
INTEGER                          :: refineBCIndex
INTEGER,ALLOCATABLE              :: RefineList(:)
REAL,ALLOCATABLE                 :: refineBoundary(:)
INTEGER,ALLOCATABLE              :: TreeToQuadRefine(:,:)
REAL                             :: sphereCenter(3),sphereRadius
REAL                             :: boxBoundary(6)
!-----------------------------------------------------------------------------------------------------------------------------------

END MODULE MOD_Refine_Vars
