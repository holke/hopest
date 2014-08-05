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
  refineFunc=C_FUNLOC(RefineByList)
  CALL InitRefineList()
CASE(3)
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
! Fills the list which elemts gets refined
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Mesh_Vars, ONLY: refineListType,Ngeo,Xgeo,RefineList,refineLevel,nElems
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
REAL,ALLOCATABLE          :: refineBoundary(:)
!===================================================================================================================================
! These are the refinement functions which are called by p4est
SELECT CASE (refineListType)
CASE(1)
  ALLOCATE(refineBoundary(6))
  refineBoundary=GETREALARRAY('refineBoundary',6)
  DO iElem=1,nElems
    XBary(1)=SUM(Xgeo(1,:,:,:,iElem))
    XBary(2)=SUM(Xgeo(2,:,:,:,iElem))
    XBary(3)=SUM(Xgeo(3,:,:,:,iElem))
    XBary(:)=XBary(:)*(1./(Ngeo+1)**3)
    IF (XBary(1) .GE. refineBoundary(1) .AND. XBary(1) .LE. refineBoundary(2) .AND. & ! refineBoundary(xmin,xmax,ymin,ymax,zmin,zmax)
        XBary(2) .GE. refineBoundary(3) .AND. XBary(2) .LE. refineBoundary(4) .AND. &
        XBary(3) .GE. refineBoundary(5) .AND. XBary(3) .LE. refineBoundary(6)) THEN
      RefineList(iElem)=refineLevel
    END IF
  END DO
CASE DEFAULT
  STOP 'refineType is not defined'
END SELECT

SDEALLOCATE(refineBoundary)
END SUBROUTINE InitRefineList

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

FUNCTION RefineByList(x,y,z,tree,level) BIND(C)
!===================================================================================================================================
! Subroutine to refine the the mesh
!===================================================================================================================================
! MODULES
USE MOD_Mesh_Vars,ONLY: RefineList,XGeo
USE, INTRINSIC :: ISO_C_BINDING
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER(KIND=C_INT32_T),INTENT(IN),VALUE :: x,y,z
INTEGER(KIND=C_INT32_T),INTENT(IN),VALUE :: tree
INTEGER(KIND=C_INT8_T ),INTENT(IN),VALUE :: level
REAL                                     :: XQuadCoord(3),ElemLength(3),ElemFirstCorner(3)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
INTEGER(KIND=C_INT)                      :: refineByList
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
XQuadCoord(1)=REAL(x)
XQuadCoord(2)=REAL(y)
XQuadCoord(3)=REAL(z)
XQuadCoord(:)=XQuadCoord(:)/2.**19

ElemFirstCorner(:)=Xgeo(:,1,1,1,tree)

IF (RefineList(tree).GT.level) THEN
  refineByList = 1
ELSE
  refineByList = 0
END IF
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
INTEGER(KIND=C_INT) :: refineByGeom
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
RefineByGeom=1
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
