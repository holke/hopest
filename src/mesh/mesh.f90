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

INTERFACE SetCurvedInfo
  MODULE PROCEDURE SetCurvedInfo
END INTERFACE

INTERFACE BuildHOMesh
  MODULE PROCEDURE BuildHOMesh
END INTERFACE

INTERFACE FinalizeMesh
  MODULE PROCEDURE FinalizeMesh
END INTERFACE

PUBLIC::InitMesh
PUBLIC::SetCurvedInfo
PUBLIC::BuildHOMesh
PUBLIC::FinalizeMesh
!===================================================================================================================================

CONTAINS

SUBROUTINE InitMesh()
!===================================================================================================================================
! Read Parameter from inputfile 
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Output_Vars, ONLY: Projectname
USE MOD_Mesh_Vars,   ONLY: BoundaryName,BoundaryType,MeshFile,nUserBCs,MeshInitIsDone
USE MOD_ReadInTools, ONLY: GETINT,GETSTR,GETINTARRAY,CNTSTR
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

SWRITE(UNIT_stdOut,'(A)')' INIT MESH DONE!'
SWRITE(UNIT_StdOut,'(132("-"))')
END SUBROUTINE InitMesh


SUBROUTINE SetCurvedInfo()
!===================================================================================================================================
! Set and allocate information related to high order data
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Mesh_Vars,ONLY: NGeo,Xi_NGeo,wBary_NGeo,HexMap,HexMapInv
USE MOD_Mesh_Vars,ONLY: nCurvedNodes 
USE MOD_Basis,    ONLY: BarycentricWeights
!-----------------------------------------------------------------------------------------------------------------------------------
IMPLICIT NONE
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER :: i,j,k,l
!===================================================================================================================================
ALLOCATE(Xi_Ngeo(0:NGeo))
ALLOCATE(wBary_Ngeo(0:Ngeo))
DO i=0,NGeo
  Xi_Ngeo(i)=-1+REAL(i)*2./REAL(NGeo)
END DO
CALL BarycentricWeights(Ngeo,xi_Ngeo,wBary_Ngeo)

! mapping form one-dimensional list [1 ; (Ngeo+1)^3] to tensor-product 0 <= i,j,k <= Ngeo and back
ALLOCATE(HexMap(0:Ngeo,0:Ngeo,0:Ngeo),HexMapInv(3,(Ngeo+1)**3))
l=0
DO k=0,Ngeo ; DO j=0,Ngeo ; DO i=0,Ngeo
  l=l+1
  HexMap(i,j,k)=l
  HexMapInv(:,l)=(/i,j,k/)
END DO ; END DO ; END DO

IF(NGeo.GT.1)THEN
  nCurvedNodes=(NGeo+1)**3
ELSE
  nCurvedNodes=0
END IF

END SUBROUTINE SetCurvedInfo


SUBROUTINE BuildHOMesh()
!===================================================================================================================================
! uses XGeo High order data from trees and interpolates it to the quadrants 
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Mesh_Vars,   ONLY: Ngeo,nElems,nQuadrants,Xgeo,XgeoQuad
USE MOD_Mesh_Vars,   ONLY: wBary_Ngeo,xi_Ngeo
USE MOD_P4EST_Vars,  ONLY: TreeToQuad,QuadCoords,QuadLevel,sIntSize
USE MOD_Basis,       ONLY: LagrangeInterpolationPolys 
USE MOD_ChangeBasis, ONLY: ChangeBasis3D_XYZ 
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                          :: xi0(3)
REAL                          :: dxi,length
REAL,DIMENSION(0:Ngeo,0:Ngeo) :: Vdm_xi,Vdm_eta,Vdm_zeta
INTEGER                       :: StartQuad,EndQuad,nQuads
INTEGER                       :: i,iQuad,iElem 
!===================================================================================================================================
ALLOCATE(XgeoQuad(3,0:Ngeo,0:Ngeo,0:Ngeo,nQuadrants))

DO iElem=1,nElems
  StartQuad = TreeToQuad(1,iElem)+1
  EndQuad   = TreeToQuad(2,iElem)
  nQuads    = TreeToQuad(2,iElem)-TreeToQuad(1,iElem)
  IF(nQuads.EQ.1)THEN !no refinement in this tree
    XgeoQuad(:,:,:,:,StartQuad)=Xgeo(:,:,:,:,iElem)
  ELSE
    DO iQuad=StartQuad,EndQuad
      ! transform p4est first corner coordinates (integer from 0... intsize) to [-1,1] reference element
      xi0(:)=-1.+2.*REAL(QuadCoords(:,iQuad))*sIntSize
      ! length of each quadrant in integers
      length=2./REAL(2**QuadLevel(iQuad))
      ! Build Vandermonde matrices for each parameter range in xi, eta,zeta
      DO i=0,Ngeo
        dxi=0.5*(xi_Ngeo(i)+1.)*Length
        CALL LagrangeInterpolationPolys(xi0(1) + dxi,Ngeo,xi_Ngeo,wBary_Ngeo,Vdm_xi(i,:)) 
        CALL LagrangeInterpolationPolys(xi0(2) + dxi,Ngeo,xi_Ngeo,wBary_Ngeo,Vdm_eta(i,:)) 
        CALL LagrangeInterpolationPolys(xi0(3) + dxi,Ngeo,xi_Ngeo,wBary_Ngeo,Vdm_zeta(i,:)) 
      END DO
      !interpolate tree HO mapping to quadrant HO mapping
      CALL ChangeBasis3D_XYZ(3,Ngeo,Ngeo,Vdm_xi,Vdm_eta,Vdm_zeta,XGeo(:,:,:,:,iElem),XgeoQuad(:,:,:,:,iQuad))
    END DO !iQuad=StartQuad,EndQuad
  END IF !nQuads==1
END DO !iElem=1,nElems
END SUBROUTINE BuildHOMesh


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
