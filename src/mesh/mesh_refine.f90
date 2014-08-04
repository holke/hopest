#include "hopest_f.h"

MODULE MOD_Mesh_Refine
!===================================================================================================================================
! Add comments please!
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------

! Public Part ----------------------------------------------------------------------------------------------------------------------
INTERFACE RefineMesh
  MODULE PROCEDURE RefineMesh
END INTERFACE

PUBLIC::RefineMesh
!===================================================================================================================================

CONTAINS


SUBROUTINE RefineMesh()
!===================================================================================================================================
! Subroutine to read the mesh from a mesh data file
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Mesh_Vars
USE MOD_p4estBinding
USE, INTRINSIC :: ISO_C_BINDING
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
TYPE(C_FUNPTR)              :: refineFunc
!===================================================================================================================================
IF(MESHInitIsDone) RETURN
SWRITE(UNIT_stdOut,'(A)')'BUILD P4EST MESH AND REFINE ...'
SWRITE(UNIT_StdOut,'(132("-"))')

refineLevel=GETINT('refineLevel','1')
refineType =GETINT('refineType','1') ! default conform refinement

ALLOCATE(RefineList(nElems))
RefineList=0

! Transform input mesh to adapted mesh
! Do refinement and save p4est refine
SELECT CASE(refineType)
CASE(1)
  refineFunc=C_FUNLOC(RefineAll)
CASE(2)
  refineFunc=C_FUNLOC(RefineByList)
CASE(3)
  refineFunc=C_FUNLOC(RefineByGeom)
CASE(4)
  refineFunc=C_FUNLOC(RefineByBC)
  DO iElem=1,nElems
    DO iSide=1,6
      IF(Elems(iElem)%ep%side(iSide)%BCIndex.GT.0)THEN
        RefineList(iElem)=RefineLevel
      END IF
    END DO
  END DO
CASE DEFAULT
  STOP 'refineType is not defined'
END SELECT

CALL p4est_refine_mesh(p4est_ptr%p4est,refineFunc,refineLevel,p4est_ptr%mesh)

SDEALLOCATE(RefineList)

END SUBROUTINE RefineMesh


! These are the refinement functions which are called by p4est

FUNCTION RefineAll(x,y,z,tree,level) BIND(C)
!===================================================================================================================================
! Subroutine to refine the the mesh
!===================================================================================================================================
! MODULES
USE, INTRINSIC :: ISO_C_BINDING
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER(KIND=C_INT32_T),INTENT(IN),VALUE :: x,y,z
INTEGER(KIND=C_INT32_T),INTENT(IN),VALUE :: tree
INTEGER(KIND=C_INT8_T ),INTENT(IN),VALUE :: level
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
INTEGER(KIND=C_INT)                      :: refineAll
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
IF(tree==0)THEN
  refineAll=1
ELSE
  refineAll=0
END IF
END FUNCTION RefineAll


FUNCTION RefineByList(x,y,z,tree,level) BIND(C)
!===================================================================================================================================
! Subroutine to refine the the mesh
!===================================================================================================================================
! MODULES
USE MOD_Mesh_Vars,ONLY: RefineList
USE, INTRINSIC :: ISO_C_BINDING
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER(KIND=C_INT32_T),INTENT(IN),VALUE :: x,y,z
INTEGER(KIND=C_INT32_T),INTENT(IN),VALUE :: tree
INTEGER(KIND=C_INT8_T ),INTENT(IN),VALUE :: level
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
INTEGER(KIND=C_INT)                      :: refineByList
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
refineByList = RefineList(tree)
RefineList(tree)=RefineList(tree)-1
END FUNCTION RefineByList


FUNCTION RefineByGeom(x,y,z,tree,level) BIND(C)
!===================================================================================================================================
! Subroutine to refine the the mesh
!===================================================================================================================================
! MODULES
USE, INTRINSIC :: ISO_C_BINDING
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER(KIND=C_INT32_T),INTENT(IN),VALUE :: x,y,z
INTEGER(KIND=C_INT32_T),INTENT(IN),VALUE :: tree
INTEGER(KIND=C_INT8_T ),INTENT(IN),VALUE :: level
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
INTEGER(KIND=C_INT)                      :: refineByGeom
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
RefineByGeom=1
END FUNCTION RefineByGeom


FUNCTION RefineByGeom(x,y,z,tree,level) BIND(C)
!===================================================================================================================================
! Subroutine to refine the the mesh
!===================================================================================================================================
! MODULES
USE, INTRINSIC :: ISO_C_BINDING
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER(KIND=C_INT32_T),INTENT(IN),VALUE :: x,y,z
INTEGER(KIND=C_INT32_T),INTENT(IN),VALUE :: tree
INTEGER(KIND=C_INT8_T ),INTENT(IN),VALUE :: level
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
INTEGER(KIND=C_INT)                      :: refineByGeom
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
RefineByGeom=1
END FUNCTION RefineByGeom


END MODULE MOD_Mesh_Refine
