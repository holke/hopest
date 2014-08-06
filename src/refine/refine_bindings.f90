#include "hopest_f.h"

MODULE MOD_Refine_Binding
!===================================================================================================================================
! Fortran <-> C++ wrapper routine for the P4est Routines
!===================================================================================================================================
! MODULES
!USE MOD_P4estBindingTypes
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE

INTERFACE 

  SUBROUTINE p4_refine_mesh(p4est,refine_function,refine_level,&
                               mesh) BIND(C)
  !=================================================================================================================================
  ! simple refine function, giving the level and if refine_elem < 0 then a conformal refinement is applied.
  !=================================================================================================================================
  ! MODULES
  USE, INTRINSIC :: ISO_C_BINDING  
  ! IMPLICIT VARIABLE HANDLING
  IMPLICIT NONE
  !---------------------------------------------------------------------------------------------------------------------------------
  ! INPUT VARIABLES
  TYPE(C_PTR),INTENT(IN),VALUE          :: p4est
  TYPE(C_FUNPTR),INTENT(IN),VALUE       :: refine_function
  INTEGER(KIND=C_INT),INTENT(IN),VALUE  :: refine_level 
  !---------------------------------------------------------------------------------------------------------------------------------
  ! OUTPUT VARIABLES
  TYPE(C_PTR),INTENT(OUT)               :: mesh
  !=================================================================================================================================
  END SUBROUTINE p4_refine_mesh 

END INTERFACE


END MODULE MOD_Refine_Binding
