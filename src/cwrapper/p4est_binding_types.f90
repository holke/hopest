#include "hopest_f.h"

MODULE MOD_P4estBindingTypes
!===================================================================================================================================
! Fortran <-> C wrapper routine for the p4est routines
!===================================================================================================================================
! MODULES
USE, INTRINSIC :: ISO_C_BINDING  
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE

TYPE :: t_p4est_ptr
  TYPE(C_PTR)                :: p4est
  TYPE(C_PTR)                :: mesh
  TYPE(C_PTR)                :: connectivity
END TYPE t_p4est_ptr

END MODULE MOD_P4estBindingTypes
