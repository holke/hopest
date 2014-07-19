#include "hopest_f.h"

MODULE MOD_P4estBinding
!===================================================================================================================================
! Fortran <-> C++ wrapper routine for the P4est Routines
!===================================================================================================================================
! MODULES
!USE MOD_P4estBindingTypes
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE

INTERFACE 

  SUBROUTINE p4est_connectivity_treevertex(num_vertices,num_trees,vertices,tree_to_vertex,p4est) BIND(C)
    USE, INTRINSIC :: ISO_C_BINDING  
    IMPLICIT NONE
    INTEGER( KIND = C_INT)     :: num_vertices 
    INTEGER( KIND = C_INT)     :: num_trees 
    REAL( KIND = C_DOUBLE )    :: Vertices(3*num_vertices)
    INTEGER( KIND = C_INT)     :: tree_to_vertex(8*num_trees) 
    INTEGER( KIND = C_INT)     :: p4est 
  END SUBROUTINE p4est_connectivity_treevertex 

  SUBROUTINE p4est_refine_mesh(p4est,mesh) BIND(C)
    USE, INTRINSIC :: ISO_C_BINDING  
    IMPLICIT NONE
    INTEGER( KIND = C_INT)     :: p4est 
    INTEGER( KIND = C_INT)     :: mesh 
  END SUBROUTINE p4est_refine_mesh 
  !SUBROUTINE ... BIND(C)
  !  USE, INTRINSIC :: ISO_C_BINDING  
  !END SUBROUTINE ...

END INTERFACE

END MODULE MOD_P4estBinding
