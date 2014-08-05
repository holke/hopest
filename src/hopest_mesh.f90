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
USE MOD_Mesh_Vars
USE MOD_Output_Vars, ONLY:Projectname
USE MOD_p4estBinding
!-----------------------------------------------------------------------------------------------------------------------------------
USE MOD_Mesh_ReadIn,        ONLY:readMeshFromHDF5
USE MOD_Mesh_Refine,        ONLY:RefineMesh
USE MOD_MeshFromP4EST,      ONLY:BuildMeshFromP4EST,BuildHOMesh,BuildBCs,GetBCs
USE MOD_Output_HDF5,        ONLY:writeMeshToHDF5
USE MOD_ReadInTools,        ONLY:GETINT,GETSTR
IMPLICIT NONE
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
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
CALL readMeshFromHDF5(MeshFile) !set nElems

CALL RefineMesh()
CALL BuildBCs()
CALL GetBCs()
!STOP 'Bis hier und nicht weiter'
CALL p4_save_all(TRIM(ProjectName)//'.p4est'//C_NULL_CHAR,p4est_ptr%p4est)
CALL BuildMeshFromP4EST()

CALL BuildHOMesh()
!output new mesh
CALL writeMeshToHDF5(TRIM(ProjectName)//'_mesh_p4est.h5')
! dealloacte pointers
SWRITE(UNIT_stdOut,'(A)') "NOW CALLING deleteMeshPointer..."
CALL deleteMeshPointer()

MeshInitIsDone=.TRUE.

CALL FinalizeHopestMesh()
SWRITE(UNIT_stdOut,'(A)')' INIT MESH DONE!'
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
