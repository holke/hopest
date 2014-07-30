#include "hopest_f.h"

MODULE MOD_Mesh
!===================================================================================================================================
! Contains subroutines to build (curviilinear) meshes and provide metrics, etc.
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES (PUBLIC)
!-----------------------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------

INTERFACE InitMesh
  MODULE PROCEDURE InitMesh
END INTERFACE

INTERFACE FinalizeMesh
  MODULE PROCEDURE FinalizeMesh
END INTERFACE

PUBLIC::InitMesh
PUBLIC::FinalizeMesh
!===================================================================================================================================

CONTAINS

SUBROUTINE InitMesh()
!===================================================================================================================================
! Read Parameter from inputfile 
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Mesh_Vars
USE MOD_Output_Vars, ONLY:Projectname
USE MOD_p4estBinding
!-----------------------------------------------------------------------------------------------------------------------------------
USE MOD_Mesh_ReadIn,        ONLY:readMesh
USE MOD_ReadInTools,        ONLY:GETLOGICAL,GETINT,GETINTARRAY,CNTSTR,GETSTR
IMPLICIT NONE
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER :: nQuadrants,nHalfFaces
INTEGER(KIND=4),ALLOCATABLE :: QuadToTree(:),QuadToQuad(:,:),QuadToHalf(:,:),QuadCoords(:,:)
INTEGER(KIND=1),ALLOCATABLE :: QuadToFace(:,:),QuadLevel(:)
REAL :: intsize
!===================================================================================================================================
IF(MeshInitIsDone)&
  CALL abort(__STAMP__,&
  'InitMesh not ready to be called or already called.')

SWRITE(UNIT_StdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A)') ' INIT MESH...'

! prepare pointer structure (get nElems, etc.)
MeshFile = GETSTR('MeshFile')
ProjectName=Meshfile(1:INDEX(Meshfile,'.h5')-1)

CALL readMesh(MeshFile) !set nElems

CALL p4est_refine_mesh(p4est_ptr%p4est,2,1,p4est_ptr%mesh,nQuadrants,nHalfFaces)
CALL p4est_save_all(TRIM(ProjectName)//'.p4est'//C_NULL_CHAR,p4est_ptr%p4est)

ALLOCATE(QuadToTree(nQuadrants),QuadToQuad(6,nQuadrants),QuadToFace(6,nQuadrants))
ALLOCATE(QuadCoords(3,nQuadrants),QuadLevel(nQuadrants)) ! big to small flip
QuadToTree=0
QuadToQuad=0
QuadToFace=0
QuadCoords=0
QuadLevel=0
IF(nHalfFaces.GT.0)THEN
  ! do not allocate if mesh is conform
  ALLOCATE(QuadToHalf(4,nHalfFaces))
  QuadToHalf=0
END IF
CALL p4est_get_quadrants(p4est_ptr%p4est,p4est_ptr%mesh,nQuadrants,nHalfFaces,& !IN
                         intsize,QuadToTree,QuadToQuad,QuadToFace,QuadToHalf,QuadCoords,QuadLevel) !OUT

! dealloacte pointers
SWRITE(UNIT_stdOut,'(A)') "NOW CALLING deleteMeshPointer..."
CALL deleteMeshPointer()

!IF(GETLOGICAL('debugmesh','.FALSE.')) CALL WriteDebugMesh()

MeshInitIsDone=.TRUE.
SWRITE(UNIT_stdOut,'(A)')' INIT MESH DONE!'
SWRITE(UNIT_StdOut,'(132("-"))')
END SUBROUTINE InitMesh


SUBROUTINE FinalizeMesh()
!============================================================================================================================
! Deallocate all global interpolation variables.
!============================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Mesh_Vars
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------
!input parameters
!----------------------------------------------------------------------------------------------------------------------------
!output parameters
!----------------------------------------------------------------------------------------------------------------------------
!local variables
!============================================================================================================================
! Deallocate global variables, needs to go somewhere else later
SDEALLOCATE(Xi_NGeo)
SDEALLOCATE(BoundaryName)
SDEALLOCATE(BoundaryType)
MeshInitIsDone = .FALSE.
END SUBROUTINE FinalizeMesh

END MODULE MOD_Mesh
