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
    INTEGER( KIND = C_INT),VALUE     :: num_vertices 
    INTEGER( KIND = C_INT),VALUE     :: num_trees 
    REAL( KIND = C_DOUBLE )          :: Vertices(3*num_vertices)
    INTEGER( KIND = C_INT)           :: tree_to_vertex(8*num_trees) 
    TYPE(C_PTR)                      :: p4est
  END SUBROUTINE p4est_connectivity_treevertex 

  SUBROUTINE p4est_refine_mesh(p4est,refine_level,refine_elem,mesh) BIND(C)
    USE, INTRINSIC :: ISO_C_BINDING  
    IMPLICIT NONE
    TYPE(C_PTR),VALUE                :: p4est
    INTEGER( KIND = C_INT),VALUE     :: refine_level 
    INTEGER( KIND = C_INT),VALUE     :: refine_elem 
    TYPE(C_PTR)                      :: mesh
  END SUBROUTINE p4est_refine_mesh 

  SUBROUTINE p4est_save_all(filename, p4est) BIND(C)
    USE, INTRINSIC :: ISO_C_BINDING  
    IMPLICIT NONE
    CHARACTER(KIND=C_CHAR),DIMENSION(*) :: filename
    TYPE(C_PTR),VALUE                   :: p4est
  END SUBROUTINE p4est_save_all

END INTERFACE

!INTERFACE StringToArray
  !MODULE PROCEDURE StringToArray
!END INTERFACE StringToArray

!PUBLIC::StringToArray

!CONTAINS

!FUNCTION StringToArray(str_in)
!!===================================================================================================================================
!! Converts a character of len* into an character array, for C strings
!!===================================================================================================================================
  !IMPLICIT NONE
  !CHARACTER(LEN=*)  :: str_in
  !CHARACTER,TARGET  :: StringToArray(LEN(TRIM(str_in)))
  !INTEGER :: i,lenstr
  !lenstr=LEN(TRIM(str_in))
  !DO i=1,lenstr
    !StringToArray(i)=str_in(i:i)
  !END DO
!END FUNCTION StringToArray


END MODULE MOD_P4estBinding
