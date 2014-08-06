#include "hopest_f.h"
MODULE MOD_Refine_Vars
!===================================================================================================================================
! Contains global variables provided by the mesh routines
!===================================================================================================================================
! MODULES
USE MOD_p4estBindingTypes
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PUBLIC
SAVE
!-----------------------------------------------------------------------------------------------------------------------------------
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
