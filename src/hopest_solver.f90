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

INTERFACE FinalizeHopestSolver
  MODULE PROCEDURE FinalizeHopestSolver
END INTERFACE

PUBLIC::HopestSolver
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
USE MOD_P4EST_Vars,         ONLY: p4est
USE MOD_P4EST,              ONLY: InitP4EST,BuildMeshFromP4EST
USE MOD_P4EST_Binding,      ONLY: p4_initvars,p4_loadmesh
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
INTEGER :: i
!===================================================================================================================================
SWRITE(UNIT_StdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A)') ' INIT MESH...'

CALL InitP4EST()
CALL InitIO()
CALL InitMesh()

CALL p4_loadmesh(MeshFile,p4est)
!CALL p4_partition_info !start and end tree and quadrants
CALL BuildMeshFromP4EST()

!STOP 'Weiter gehts nicht'
CALL ReadGeoFromHDF5(MeshFile)

CALL BuildHOMesh()

!CALL FlexiPrepareMesh() ! Suggestion


! dealloacte pointers
SWRITE(UNIT_stdOut,'(A)') "NOW CALLING deleteMeshPointer..."
CALL deleteMeshPointer()

CALL FinalizeHopestSolver()
SWRITE(UNIT_stdOut,'(A)')' INIT MESH DONE!'
SWRITE(UNIT_StdOut,'(132("-"))')
END SUBROUTINE HopestSolver


SUBROUTINE FinalizeHopestSolver()
!============================================================================================================================
! Deallocate all global interpolation variables.
!============================================================================================================================
! MODULES
USE MOD_P4EST,ONLY: FinalizeP4EST
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------
!input parameters
!----------------------------------------------------------------------------------------------------------------------------
!output parameters
!----------------------------------------------------------------------------------------------------------------------------
!local variables
!============================================================================================================================
CALL FinalizeP4EST()
END SUBROUTINE FinalizeHopestSolver


END MODULE MOD_HopestSolver
