#include "hopest_f.h"

MODULE MOD_HopestMesh
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

INTERFACE HopestMesh
  MODULE PROCEDURE HopestMesh
END INTERFACE

INTERFACE FinalizeHopestMesh
  MODULE PROCEDURE FinalizeHopestMesh
END INTERFACE

PUBLIC::HopestMesh
PUBLIC::FinalizeHopestMesh
!===================================================================================================================================

CONTAINS

SUBROUTINE HopestMesh()
!===================================================================================================================================
! Read Parameter from inputfile 
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_IO_HDF5
USE MOD_P4EST_Vars,         ONLY: p4est,p4estFile
USE MOD_P4EST,              ONLY: InitP4EST,BuildMeshFromP4EST,BuildBCs
USE MOD_P4EST_Binding,      ONLY: p4_savemesh
USE MOD_Mesh_Vars
USE MOD_Mesh,               ONLY: InitMesh,BuildHOMesh
USE MOD_Output_Vars,        ONLY: Projectname
USE MOD_Output_HDF5,        ONLY: writeMeshToHDF5
USE MOD_Mesh_ReadIn,        ONLY: readMeshFromHDF5
USE MOD_Basis,              ONLY: BarycentricWeights
USE MOD_Refine,             ONLY: RefineMesh
USE MOD_ReadInTools,        ONLY: GETINT,GETSTR,GETINTARRAY,CNTSTR
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
SWRITE(UNIT_StdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A)') ' START HOPEST MESH...'

CALL InitMesh()
CALL InitP4EST()
CALL InitIO()

CALL readMeshFromHDF5(MeshFile) !set nElems

CALL RefineMesh()
CALL BuildBCs()
CALL p4_savemesh(TRIM(p4estFile)//C_NULL_CHAR,p4est)
CALL BuildMeshFromP4EST()

CALL BuildHOMesh()
!output new mesh
CALL writeMeshToHDF5(TRIM(ProjectName)//'_mesh_p4est.h5')
! dealloacte pointers
SWRITE(UNIT_stdOut,'(A)') "NOW CALLING deleteMeshPointer..."
CALL deleteMeshPointer()

CALL FinalizeHopestMesh()
SWRITE(UNIT_stdOut,'(A)')' HOPEST MESH DONE!'
SWRITE(UNIT_StdOut,'(132("-"))')
END SUBROUTINE HopestMesh


SUBROUTINE FinalizeHopestMesh()
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
END SUBROUTINE FinalizeHopestMesh

END MODULE MOD_HopestMesh
