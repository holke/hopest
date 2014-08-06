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

  SUBROUTINE p4_initvars() BIND(C)
  !=================================================================================================================================
  ! initialize MPI, SC, p4est
  !=================================================================================================================================
  ! MODULES
  USE, INTRINSIC :: ISO_C_BINDING  
  ! IMPLICIT VARIABLE HANDLING
  IMPLICIT NONE
  !---------------------------------------------------------------------------------------------------------------------------------
  ! INPUT VARIABLES
  !---------------------------------------------------------------------------------------------------------------------------------
  ! OUTPUT VARIABLES
  !=================================================================================================================================
  END SUBROUTINE p4_initvars
  

  SUBROUTINE p4_loadmesh(filename,p4est) BIND(C)
  !=================================================================================================================================
  ! read p4est connectivity from file and build p4est
  !=================================================================================================================================
  ! MODULES
  USE, INTRINSIC :: ISO_C_BINDING  
  ! IMPLICIT VARIABLE HANDLING
  IMPLICIT NONE
  !---------------------------------------------------------------------------------------------------------------------------------
  ! INPUT VARIABLES
  CHARACTER(KIND=C_CHAR),DIMENSION(*) :: filename
  !---------------------------------------------------------------------------------------------------------------------------------
  ! OUTPUT VARIABLES
  TYPE(C_PTR)                         :: p4est
  !=================================================================================================================================
  END SUBROUTINE p4_loadmesh 


  SUBROUTINE p4_connectivity_treevertex(num_vertices,num_trees,vertices,tree_to_vertex,&
                                        num_periodics,JoinFaces,connectivity) BIND(C)
  !=================================================================================================================================
  ! builds up p4est connectivit, using only element connectivity and vertex positions
  !=================================================================================================================================
  ! MODULES
  USE, INTRINSIC :: ISO_C_BINDING  
  ! IMPLICIT VARIABLE HANDLING
  IMPLICIT NONE
  !---------------------------------------------------------------------------------------------------------------------------------
  ! INPUT VARIABLES
  INTEGER( KIND = C_INT),VALUE                :: num_vertices 
  INTEGER( KIND = C_INT),VALUE                :: num_trees 
  REAL( KIND = C_DOUBLE )                     :: Vertices(3,num_vertices)
  INTEGER( KIND = C_INT)                      :: tree_to_vertex(8*num_trees) 
  INTEGER( KIND = C_INT),VALUE                :: num_periodics 
  INTEGER( KIND = C_INT)                      :: JoinFaces(5*num_periodics) 
  !---------------------------------------------------------------------------------------------------------------------------------
  ! OUTPUT VARIABLES
  TYPE(C_PTR)                                 :: connectivity
  !=================================================================================================================================
  END SUBROUTINE p4_connectivity_treevertex 


  SUBROUTINE p4_build_p4est(connectivity,p4est) BIND(C)
  !=================================================================================================================================
  ! builds up p4est connectivit, using only element connectivity and vertex positions
  !=================================================================================================================================
  ! MODULES
  USE, INTRINSIC :: ISO_C_BINDING  
  ! IMPLICIT VARIABLE HANDLING
  IMPLICIT NONE
  !---------------------------------------------------------------------------------------------------------------------------------
  ! INPUT VARIABLES
  TYPE(C_PTR),INTENT(IN),VALUE     :: connectivity
  !---------------------------------------------------------------------------------------------------------------------------------
  ! OUTPUT VARIABLES
  TYPE(C_PTR)                      :: p4est
  !=================================================================================================================================
  END SUBROUTINE p4_build_p4est


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


  SUBROUTINE p4_build_bcs(p4est,num_trees,bcelemmap) BIND(C)
  !=================================================================================================================================
  ! simple refine function, giving the level and if refine_elem < 0 then a conformal refinement is applied.
  !=================================================================================================================================
  ! MODULES
  USE, INTRINSIC :: ISO_C_BINDING  
  ! IMPLICIT VARIABLE HANDLING
  IMPLICIT NONE
  !---------------------------------------------------------------------------------------------------------------------------------
  ! INPUT VARIABLES
  TYPE(C_PTR),VALUE,INTENT(IN)       :: p4est
  INTEGER( KIND = C_INT),VALUE       :: num_trees 
  INTEGER(KIND=C_INT16_T),INTENT(IN) :: bcelemmap(0:5,num_trees)
  !---------------------------------------------------------------------------------------------------------------------------------
  ! OUTPUT VARIABLES
  !=================================================================================================================================
  END SUBROUTINE p4_build_bcs


  SUBROUTINE p4_get_bcs(p4est,tree_to_bc) BIND(C)
  !=================================================================================================================================
  ! simple refine function, giving the level and if refine_elem < 0 then a conformal refinement is applied.
  !=================================================================================================================================
  ! MODULES
  USE, INTRINSIC :: ISO_C_BINDING  
  ! IMPLICIT VARIABLE HANDLING
  IMPLICIT NONE
  !---------------------------------------------------------------------------------------------------------------------------------
  ! INPUT VARIABLES
  TYPE(C_PTR),VALUE,INTENT(IN)        :: p4est
  TYPE(C_PTR),INTENT(OUT)             :: tree_to_bc
  !---------------------------------------------------------------------------------------------------------------------------------
  ! OUTPUT VARIABLES
  !=================================================================================================================================
  END SUBROUTINE p4_get_bcs


  SUBROUTINE p4_get_mesh_info(p4est,mesh,global_num_quadrants,num_half_faces,num_trees) BIND(C)
  !=================================================================================================================================
  ! simple refine function, giving the level and if refine_elem < 0 then a conformal refinement is applied.
  !=================================================================================================================================
  ! MODULES
  USE, INTRINSIC :: ISO_C_BINDING  
  ! IMPLICIT VARIABLE HANDLING
  IMPLICIT NONE
  !---------------------------------------------------------------------------------------------------------------------------------
  ! INPUT VARIABLES
  TYPE(C_PTR),VALUE,INTENT(IN)     :: p4est
  TYPE(C_PTR),VALUE,INTENT(IN)     :: mesh
  !---------------------------------------------------------------------------------------------------------------------------------
  ! OUTPUT VARIABLES
  INTEGER,INTENT(OUT)              :: global_num_quadrants
  INTEGER,INTENT(OUT)              :: num_half_faces
  INTEGER,INTENT(OUT)              :: num_trees
  !=================================================================================================================================
  END SUBROUTINE p4_get_mesh_info


  SUBROUTINE p4_get_quadrants(p4est,mesh,global_num_quadrants,num_half_faces,&
                                 intsize,quad_to_tree,quad_to_quad,quad_to_face,quad_to_half,&
                                 quadcoords,quadlevel) BIND(C)
  !=================================================================================================================================
  ! simple refine function, giving the level and if refine_elem < 0 then a conformal refinement is applied.
  !=================================================================================================================================
  ! MODULES
  USE, INTRINSIC :: ISO_C_BINDING  
  ! IMPLICIT VARIABLE HANDLING
  IMPLICIT NONE
  !---------------------------------------------------------------------------------------------------------------------------------
  ! INPUT VARIABLES
  TYPE(C_PTR),VALUE         :: p4est
  TYPE(C_PTR),VALUE         :: mesh
  INTEGER(KIND=C_INT),VALUE :: global_num_quadrants
  INTEGER(KIND=C_INT),VALUE :: num_half_faces
  !---------------------------------------------------------------------------------------------------------------------------------
  ! OUTPUT VARIABLES
  INTEGER(KIND=4),INTENT(OUT) :: intsize ! P4EST_ROOT_LEN -> int2real transform in parameter space, REAL=1/intsize*INT [0,1]
  TYPE(C_PTR),INTENT(OUT)     :: quad_to_tree
  TYPE(C_PTR),INTENT(OUT)     :: quad_to_quad
  TYPE(C_PTR),INTENT(OUT)     :: quad_to_face
  TYPE(C_PTR),INTENT(OUT)     :: quad_to_half
  INTEGER(KIND=4),INTENT(OUT) :: quadcoords(  3,global_num_quadrants)
  INTEGER(KIND=1),INTENT(OUT) :: quadlevel(     global_num_quadrants)
  !=================================================================================================================================
  END SUBROUTINE p4_get_quadrants


  SUBROUTINE p4_savemesh(filename, p4est) BIND(C)
  !=================================================================================================================================
  ! save the p4est data  to a p4est state file 
  !=================================================================================================================================
  ! MODULES
  USE, INTRINSIC :: ISO_C_BINDING  
  ! IMPLICIT VARIABLE HANDLING
  IMPLICIT NONE
  !---------------------------------------------------------------------------------------------------------------------------------
  ! INPUT VARIABLES
  CHARACTER(KIND=C_CHAR),DIMENSION(*) :: filename
  TYPE(C_PTR),VALUE                   :: p4est
  !---------------------------------------------------------------------------------------------------------------------------------
  ! OUTPUT VARIABLES
  !=================================================================================================================================
  END SUBROUTINE p4_savemesh
  

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
