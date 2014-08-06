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
INTEGER                          :: refineListType
INTEGER                          :: refineBCIndex
INTEGER,ALLOCATABLE              :: RefineList(:)
REAL,ALLOCATABLE                 :: refineBoundary(:)
INTEGER,ALLOCATABLE              :: TreeToQuadRefine(:,:)
!-----------------------------------------------------------------------------------------------------------------------------------

END MODULE MOD_Refine_Vars
