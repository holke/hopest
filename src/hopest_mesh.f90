#include "hopest_f.h"

MODULE MODH_HopestMesh
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
USE, INTRINSIC :: ISO_C_BINDING
USE MODH_Globals
USE MODH_IO_HDF5,            ONLY: InitIO
USE MODH_P4EST_Vars,         ONLY: connectivity,p4est,geom,p4estFile
USE MODH_P4EST,              ONLY: InitP4EST,BuildMeshFromP4EST,BuildBCs,FinalizeP4EST
USE MODH_P4EST_Binding,      ONLY: p4_savemesh,p4_build_p4est
USE MODH_Mesh_Vars,          ONLY: MeshFile
USE MODH_Mesh,               ONLY: InitMesh,BuildHOMesh,FinalizeMesh
USE MODH_Analyze,            ONLY: InitAnalyze,Analyze
USE MODH_Output_Vars,        ONLY: Projectname
USE MODH_Output_HDF5,        ONLY: writeMeshToHDF5
USE MODH_Mesh_ReadIn,        ONLY: ReadMeshFromHDF5,ReadMeshHeader
USE MODH_Mesh,               ONLY: DeformMesh
USE MODH_Refine,             ONLY: RefineMesh
IMPLICIT NONE
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
SWRITE(UNIT_StdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A)') ' START HOPEST MESH...'

! If this routine is called then hopest is run in solver mode
hopestMode = 1

CALL InitMesh()
CALL InitP4EST()
CALL InitIO()

CALL ReadMeshHeader(MeshFile)   ! read mesh header file including BCs
CALL ReadMeshFromHDF5(MeshFile) ! read tree connectivity and geometry
CALL p4_build_p4est(connectivity,p4est,geom)


CALL deformMesh()

CALL RefineMesh()               ! perform mesh refinement using p4est
CALL BuildBCs()                 ! store BCs in p4est file
CALL p4_savemesh(TRIM(p4estFile)//C_NULL_CHAR,p4est)
CALL BuildMeshFromP4EST()       ! build p4est mesh in flexi datastructure

CALL BuildHOMesh()

CALL InitAnalyze()
CALL Analyze()
!output new mesh
CALL writeMeshToHDF5(TRIM(ProjectName)//'_mesh_p4est.h5')
! dealloacte pointers
SWRITE(UNIT_stdOut,'(A)') "NOW CALLING deleteMeshPointer..."
!CALL deleteMeshPointer()

CALL FinalizeP4EST()
CALL FinalizeMesh()
CALL FinalizeHopestMesh()

SWRITE(UNIT_stdOut,'(A)')' HOPEST MESH DONE!'
SWRITE(UNIT_StdOut,'(132("-"))')
END SUBROUTINE HopestMesh


SUBROUTINE FinalizeHopestMesh()
!============================================================================================================================
! Deallocate all global interpolation variables.
!============================================================================================================================
! MODULES
USE MODH_Globals
USE MODH_Mesh_Vars
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
END SUBROUTINE FinalizeHopestMesh

END MODULE MODH_HopestMesh
