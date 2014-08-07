#include "hopest_f.h"

MODULE MOD_HopestSolver
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

INTERFACE HopestSolver
  MODULE PROCEDURE HopestSolver
END INTERFACE

INTERFACE PrepareMesh
  MODULE PROCEDURE PrepareMesh
END INTERFACE

INTERFACE FinalizeHopestSolver
  MODULE PROCEDURE FinalizeHopestSolver
END INTERFACE

PUBLIC::HopestSolver
PUBLIC::PrepareMesh
PUBLIC::FinalizeHopestSolver
!===================================================================================================================================

CONTAINS

SUBROUTINE HopestSolver()
!===================================================================================================================================
! Read Parameter from inputfile 
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_IO_HDF5
USE MOD_P4EST_Vars,         ONLY: p4est,p4estFile
USE MOD_P4EST,              ONLY: InitP4EST,BuildMeshFromP4EST
USE MOD_P4EST_Binding,      ONLY: p4_loadmesh
USE MOD_Mesh_Vars
USE MOD_Mesh,               ONLY: InitMesh,BuildHOMesh
USE MOD_Mesh_ReadIn,        ONLY: ReadGeoFromHDF5
USE MOD_ReadInTools,        ONLY: GETINT,GETSTR
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
SWRITE(UNIT_StdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A)') ' INIT MESH...'

CALL InitMesh()
CALL InitP4EST()
CALL InitIO()

CALL p4_loadmesh(TRIM(p4estFile)//C_NULL_CHAR,p4est)
!CALL p4_partition_info !start and end tree and quadrants
CALL BuildMeshFromP4EST()

!STOP 'Weiter gehts nicht'
CALL ReadGeoFromHDF5(MeshFile)
CALL BuildHOMesh()

SWRITE(UNIT_stdOut,'(A)')' INIT MESH DONE!'
SWRITE(UNIT_StdOut,'(132("-"))')
END SUBROUTINE HopestSolver


SUBROUTINE PrepareMesh()
!===================================================================================================================================
! Read Parameter from inputfile 
!===================================================================================================================================
! MODULES
USE MOD_Prepare_Mesh
USE MOD_Mesh_Vars,ONLY: nQuads,nSides,nBCSides
USE MOD_Mesh_Vars,ONLY: ElemToSide,SideToElem,BC,AnalyzeSide
IMPLICIT NONE
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
CALL countSides()
ALLOCATE(ElemToSide(2,6,nQuads))
ALLOCATE(SideToElem(5,nSides))
ALLOCATE(BC(nBCSides))
ALLOCATE(AnalyzeSide(nSides))
CALL setLocalSideIDs()
#ifdef MPI
CALL exchangeFlip()    ! should be already known from p4est
#endif
CALL fillMeshInfo()
END SUBROUTINE PrepareMesh


SUBROUTINE FinalizeHopestSolver()
!============================================================================================================================
! Deallocate all global interpolation variables.
!============================================================================================================================
! MODULES
USE MOD_Mesh_Vars,ONLY: deleteMeshPointer
USE MOD_P4EST,    ONLY: FinalizeP4EST
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------
!input parameters
!----------------------------------------------------------------------------------------------------------------------------
!output parameters
!----------------------------------------------------------------------------------------------------------------------------
!local variables
!============================================================================================================================
!CALL deleteMeshPointer()
CALL FinalizeP4EST()
END SUBROUTINE FinalizeHopestSolver


END MODULE MOD_HopestSolver
