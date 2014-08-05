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
USE MOD_Readintools,ONLY:GETINT
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
  refineListType =GETINT('refineListType') ! 
  CALL InitRefineList()
  refineFunc=C_FUNLOC(RefineByList)
CASE(3)
  refineListType =GETINT('refineListType') ! 
  CALL InitRefineGeom()
  refineFunc=C_FUNLOC(RefineByGeom)
CASE(11)
  refineFunc=C_FUNLOC(RefineFirst)
CASE DEFAULT
  STOP 'refineType is not defined'
END SELECT

CALL p4est_refine_mesh(p4est_ptr%p4est,refineFunc,refineLevel,& !IN
                       p4est_ptr%mesh)                              !OUT

SDEALLOCATE(RefineList)

END SUBROUTINE RefineMesh

SUBROUTINE InitRefineList()
!===================================================================================================================================
! init the refinment list
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Mesh_Vars, ONLY: Elems,refineListType,nElems,TreeToQuadRefine
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                     :: iElem,iSide
!===================================================================================================================================
! These are the refinement functions which are called by p4est
SELECT CASE (refineListType)
CASE(1)
  ALLOCATE(TreeToQuadRefine(1:8,1:nElems))
  TreeToQuadRefine=0
  DO iElem=1,nElems
    DO iSide=1,6
      IF (Elems(iElem)%ep%Side(iSide)%sp%BCIndex.EQ.3) THEN
        SELECT CASE (iSide)
          CASE (1) 
            TreeToQuadRefine(1:4,iElem)=1
          CASE (2) 
            TreeToQuadRefine(2,iElem)=1
            TreeToQuadRefine(4,iElem)=1
            TreeToQuadRefine(6,iElem)=1
            TreeToQuadRefine(8,iElem)=1
          CASE (3) 
            TreeToQuadRefine(3,iElem)=1
            TreeToQuadRefine(4,iElem)=1
            TreeToQuadRefine(7,iElem)=1
            TreeToQuadRefine(8,iElem)=1
          CASE (4) 
            TreeToQuadRefine(1,iElem)=1
            TreeToQuadRefine(3,iElem)=1
            TreeToQuadRefine(5,iElem)=1
            TreeToQuadRefine(7,iElem)=1
          CASE (5) 
            TreeToQuadRefine(1,iElem)=1
            TreeToQuadRefine(2,iElem)=1
            TreeToQuadRefine(5,iElem)=1
            TreeToQuadRefine(6,iElem)=1
          CASE (6) 
            TreeToQuadRefine(5:8,iElem)=1
        END SELECT
      END IF
    END DO
  END DO
CASE DEFAULT
  STOP 'refineType is not defined'
END SELECT

END SUBROUTINE InitRefineList


SUBROUTINE InitRefineGeom()
!===================================================================================================================================
! init the geometric refinment
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Mesh_Vars, ONLY: refineBoundary,refineListType
USE MOD_Readintools,ONLY:GETREALARRAY
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                     :: iElem
REAL                     :: XBary(3)
!===================================================================================================================================
! These are the refinement functions which are called by p4est
SELECT CASE (refineListType)
CASE(1)
  ALLOCATE(refineBoundary(6))
  refineBoundary=GETREALARRAY('refineBoundary',6)
CASE(2)
  ALLOCATE(refineBoundary(4))
  refineBoundary=GETREALARRAY('refineBoundary',4)
CASE(3)
  ALLOCATE(refineBoundary(5))
  refineBoundary=GETREALARRAY('refineBoundary',5)
CASE DEFAULT
  STOP 'refineType is not defined'
END SELECT

END SUBROUTINE InitRefineGeom

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
INTEGER(KIND=C_INT) :: refineAll
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
refineAll=1
END FUNCTION RefineAll

FUNCTION RefineByList(x,y,z,tree,level,childID) BIND(C)
!===================================================================================================================================
! Subroutine to refine the the mesh
!===================================================================================================================================
! MODULES
USE, INTRINSIC :: ISO_C_BINDING
USE MOD_Mesh_Vars, ONLY: TreeToQuadRefine
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER(KIND=C_INT32_T),INTENT(IN),VALUE :: x,y,z
INTEGER(KIND=C_INT32_T),INTENT(IN),VALUE :: tree
INTEGER(KIND=C_INT8_T ),INTENT(IN),VALUE :: level
INTEGER(KIND=C_INT ),INTENT(IN),VALUE :: childID
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
INTEGER(KIND=C_INT)                      :: refineByList
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------

IF (level.EQ.0) RefineByList=SUM(TreeToQuadRefine(:,tree))
IF (level.GE.1) RefineByList=TreeToQuadRefine(childID,tree)
END FUNCTION RefineByList


FUNCTION RefineByGeom(x,y,z,tree,level) BIND(C)
!===================================================================================================================================
! Subroutine to refine the the mesh
!===================================================================================================================================
! MODULES
USE, INTRINSIC :: ISO_C_BINDING
USE MOD_Mesh_Vars,ONLY: RefineList,XGeo,Ngeo,refineBoundary
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER(KIND=C_INT32_T),INTENT(IN),VALUE :: x,y,z
INTEGER(KIND=C_INT32_T),INTENT(IN),VALUE :: tree
INTEGER(KIND=C_INT8_T ),INTENT(IN),VALUE :: level
REAL                                     :: XQuadCoord(3),ElemLength(3),ElemFirstCorner(3),VectorToBaryQuad(3)
REAL                                     :: XBaryQuad(3),lengthQuad,test
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
INTEGER(KIND=C_INT) :: refineByGeom
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
XQuadCoord(1)=REAL(x)
XQuadCoord(2)=REAL(y)
XQuadCoord(3)=REAL(z)
XQuadCoord(:)=XQuadCoord(:)/2.**19

ElemFirstCorner(:)=Xgeo(:,0,0,0,tree)
ElemLength(1)=abs(ElemFirstCorner(1)-XGeo(1,Ngeo,0,0,tree))
ElemLength(2)=abs(ElemFirstCorner(2)-XGeo(2,0,Ngeo,0,tree))
ElemLength(3)=abs(ElemFirstCorner(3)-XGeo(3,0,0,Ngeo,tree))
lengthQuad=1./REAL(2**level)

VectorToBaryQuad(:)=XQuadCoord(:)+lengthQuad/2.

XBaryQuad(:)=VectorToBaryQuad(:)*ElemLength(:)+ElemFirstCorner(:)

SELECT CASE (SIZE(refineBoundary))
CASE(4)
  test=SQRT((XBaryQuad(1)-refineBoundary(1))**2+(XBaryQuad(2)-refineBoundary(2))**2+(XBaryQuad(3)-refineBoundary(3))**2)

  IF (test.LE.refineBoundary(4)) THEN
    refineByGeom = 1
  ELSE
    refineByGeom = 0
  END IF
CASE(6)
  ! refineBoundary(xmin,xmax,ymin,ymax,zmin,zmax)
  IF (XBaryQuad(1) .GE. refineBoundary(1) .AND. XBaryQuad(1) .LE. refineBoundary(2) .AND. & 
      XBaryQuad(2) .GE. refineBoundary(3) .AND. XBaryQuad(2) .LE. refineBoundary(4) .AND. &
      XBaryQuad(3) .GE. refineBoundary(5) .AND. XBaryQuad(3) .LE. refineBoundary(6)) THEN
    refineByGeom = 1
  ELSE
    refineByGeom = 0
  END IF
END SELECT
END FUNCTION RefineByGeom


FUNCTION RefineFirst(x,y,z,tree,level) BIND(C)
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
INTEGER(KIND=C_INT)                      :: refineFirst
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
IF(tree.EQ.1)THEN
  refineFirst=1
ELSE
  refineFirst=0
END IF
END FUNCTION RefineFirst

END MODULE MOD_Mesh_Refine
