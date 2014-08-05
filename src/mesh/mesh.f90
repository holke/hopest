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
USE MOD_IO_HDF5
USE MOD_Mesh_Vars
USE MOD_Output_Vars, ONLY:Projectname
USE MOD_p4estBinding
!-----------------------------------------------------------------------------------------------------------------------------------
USE MOD_Mesh_ReadIn,        ONLY:readMeshFromHDF5
USE MOD_Basis,              ONLY:BarycentricWeights
USE MOD_Mesh_Refine,        ONLY:RefineMesh
USE MOD_MeshFromP4EST,      ONLY:BuildMeshFromP4EST,BuildHOMesh,BuildBCs
USE MOD_Output_HDF5,        ONLY:writeMeshToHDF5
USE MOD_ReadInTools,        ONLY:GETINT,GETSTR,GETINTARRAY,CNTSTR
IMPLICIT NONE
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER :: i
!===================================================================================================================================
CALL InitIO()
CALL p4_initvars()

IF(MeshInitIsDone)&
  CALL abort(__STAMP__,&
  'InitMesh not ready to be called or already called.')

SWRITE(UNIT_StdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A)') ' INIT MESH...'

! prepare pointer structure (get nElems, etc.)
MeshFile = GETSTR('MeshFile')
ProjectName=Meshfile(1:INDEX(Meshfile,'_mesh.h5')-1)

! read in boundary conditions, will overwrite BCs from meshfile!
nUserBCs = CNTSTR('BoundaryName',0)
IF(nUserBCs.GT.0)THEN
  ALLOCATE(BoundaryName(1:nUserBCs))
  ALLOCATE(BoundaryType(1:nUserBCs,2))
  DO i=1,nUserBCs
    BoundaryName(i)   = GETSTR('BoundaryName')
    BoundaryType(i,:) = GETINTARRAY('BoundaryType',2) !(/Type,State/)
  END DO
END IF !nUserBCs>0
CALL readMeshFromHDF5(MeshFile) !set nElems

ALLOCATE(Xi_Ngeo(0:Ngeo))
ALLOCATE(wBary_Ngeo(0:Ngeo))

DO i=0,NGeo
  Xi_Ngeo(i)=-1+REAL(i)*2./REAL(NGeo)
END DO
CALL BarycentricWeights(Ngeo,xi_Ngeo,wBary_Ngeo)

CALL RefineMesh()
CALL BuildBCs()
CALL p4_save_all(TRIM(ProjectName)//'.p4est'//C_NULL_CHAR,p4est_ptr%p4est)
CALL BuildMeshFromP4EST()

CALL BuildHOMesh()
!output new mesh
CALL writeMeshToHDF5(TRIM(ProjectName)//'_mesh_p4est.h5')
! dealloacte pointers
SWRITE(UNIT_stdOut,'(A)') "NOW CALLING deleteMeshPointer..."
CALL deleteMeshPointer()

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
