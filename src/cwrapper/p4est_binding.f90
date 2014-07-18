#include "hopest_f.h"

MODULE MOD_P4estBinding
!===================================================================================================================================
! Fortran <-> C++ wrapper routine for the P4est Routines
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE


INTERFACE 

  SUBROUTINE test(areal,bint) BIND(C) 
    USE, INTRINSIC :: ISO_C_BINDING  
    IMPLICIT NONE
    REAL( KIND = C_DOUBLE )       :: areal 
    INTEGER( KIND = C_INT)        :: bint 
    !INTEGER( kind = c_int),value  :: bint   ! call by value, wenn kein rückgabewert
  END SUBROUTINE 

  SUBROUTINE p4est_connectivity_treevertex(num_vertices,num_trees,vertices,tree_to_vertex) BIND(C)
    USE, INTRINSIC :: ISO_C_BINDING  
    IMPLICIT NONE
    INTEGER( KIND = C_INT)     :: num_vertices 
    INTEGER( KIND = C_INT)     :: num_trees 
    REAL( KIND = C_DOUBLE )    :: Vertices(3*num_vertices)
    INTEGER( KIND = C_INT)     :: tree_to_vertex(8*num_trees) 
  END SUBROUTINE p4est_connectivity_treevertex 

  !SUBROUTINE ... BIND(C)
  !  USE, INTRINSIC :: ISO_C_BINDING  
  !END SUBROUTINE ...

END INTERFACE

END MODULE MOD_P4estBinding
