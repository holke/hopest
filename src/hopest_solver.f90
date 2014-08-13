#include "hopest_f.h"

MODULE MODH_HopestSolver
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
USE MODH_Globals
USE MODH_MPI
USE MODH_IO_HDF5
USE MODH_P4EST_Vars,         ONLY: p4est,p4estFile
USE MODH_P4EST,              ONLY: InitP4EST,BuildMeshFromP4EST
USE MODH_P4EST_Binding,      ONLY: p4_loadmesh
USE MODH_Mesh_Vars
USE MODH_Mesh,               ONLY: InitMesh,BuildHOMesh
USE MODH_Mesh_ReadIn,        ONLY: ReadGeoFromHDF5,ReadMeshHeader
USE MODH_ReadInTools,        ONLY: GETINT,GETSTR
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
LOGICAL :: fileExists
!===================================================================================================================================
SWRITE(UNIT_StdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A)') ' INIT MESH UNO...'

! If this routine is called then hopest is run in solver mode
hopestMode = 2

CALL InitMPI()
CALL InitMesh()
CALL InitP4EST()
CALL InitIO()

INQUIRE(FILE=TRIM(p4estFile), EXIST=fileExists)
IF(.NOT.fileExists) CALL abort(__STAMP__,&
  'p4est Input file does not exists: '//p4estFile)
CALL p4_loadmesh(TRIM(p4estFile)//C_NULL_CHAR,p4est)
!CALL p4_partition_info !start and end tree and quadrants

CALL ReadMeshHeader(MeshFile)   ! get BC information
CALL BuildMeshFromP4EST()       ! build mesh connectivity
CALL ReadGeoFromHDF5(MeshFile)  ! get tree geometry
CALL BuildHOMesh()              ! build quads from trees

SWRITE(UNIT_stdOut,'(A)')' INIT MESH DONE!'
SWRITE(UNIT_StdOut,'(132("-"))')
END SUBROUTINE HopestSolver


SUBROUTINE PrepareMesh()
!===================================================================================================================================
! Read Parameter from inputfile 
!===================================================================================================================================
! MODULES
USE MODH_Prepare_Mesh
USE MODH_Mesh_Vars
IMPLICIT NONE
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
ALLOCATE(ElemToSide(2,6,nElems))
ALLOCATE(SideToElem(5,nSides))
ALLOCATE(BC(nBCSides))
ALLOCATE(AnalyzeSide(nSides))
CALL setLocalSideIDs()
#ifdef MPI
CALL exchangeFlip()    ! should be already known from p4est
#endif
CALL fillMeshInfo()

!MortarStuff
firstMortarSideID = nBCSides+1
lastMortarSideID  = nBCSides+nMortarSides                                                             
ALLOCATE(MortarType(1:nSides))
! FOR PARALLIZATION : ALLOCATE(MortarFlip(1:4,firstMortarSideID:lastMortarSideID))                    
ALLOCATE(Mortar_nbSideID(1:4,firstMortarSideID:lastMortarSideID))                                     
ALLOCATE(Mortar_Flip(1:4,firstMortarSideID:lastMortarSideID))                                         
MortarType=0
Mortar_nbSideID=0                                                                                     
Mortar_Flip=-1                                                                                        

!lower and upper index of U/gradUx/y/z _plus                                                          
!lower and upper index of U/gradUx/y/z _plus                                                          
sideID_minus_lower = 1
sideID_minus_upper = nBCSides+nMortarSides+nInnerSides+nMPISides_MINE                                 
sideID_plus_lower  = nBCSides+nMortarSides+1
sideID_plus_upper  = nBCSides+nMortarSides+nInnerSides+nMPISides
write(*,*) nSides,nBCSides

write(*,*) sideID_minus_lower 
write(*,*) sideID_minus_upper 
write(*,*) sideID_plus_lower  
write(*,*) sideID_plus_upper  

!! dealloacte pointers
!SWRITE(UNIT_stdOut,'(A)') "NOW CALLING deleteMeshPointer..."
!CALL deleteMeshPointer()

END SUBROUTINE PrepareMesh


SUBROUTINE FinalizeHopestSolver()
!============================================================================================================================
! Deallocate all global interpolation variables.
!============================================================================================================================
! MODULES
USE MODH_Mesh_Vars,ONLY: deleteMeshPointer
USE MODH_P4EST,    ONLY: FinalizeP4EST
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


END MODULE MODH_HopestSolver
