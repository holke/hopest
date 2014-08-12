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

INTERFACE DeformMesh
  MODULE PROCEDURE DeformMesh
END INTERFACE

INTERFACE FinalizeMesh
  MODULE PROCEDURE FinalizeMesh
END INTERFACE

PUBLIC::InitMesh
PUBLIC::SetCurvedInfo
PUBLIC::BuildHOMesh
PUBLIC::DeformMesh
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
USE MOD_Mesh_Vars,   ONLY: BoundaryName,BoundaryType,MeshFile,nUserBCs,Deform
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

IF(CNTSTR('ProjectName',0).EQ.0)THEN
  !default project name frommesh file
  ProjectName=TRIM(Meshfile(1:INDEX(Meshfile,'_mesh.h5')-1))
ELSE
  ProjectName = GETSTR('ProjectName')
END IF

Deform = GETINT('Deform','0')

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
USE MOD_Mesh_Vars,ONLY: NGeo_out,XiCL_NGeo_out,Vdm_CL_EQ_out,HexMap_out
USE MOD_Mesh_Vars,ONLY: nCurvedNodes 
USE MOD_Basis,    ONLY: BarycentricWeights,ChebyGaussLobNodesAndWeights,InitializeVandermonde
USE MOD_ReadInTools, ONLY: GETINT
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
REAL,ALLOCATABLE    :: Xi_Ngeo_out(:),wBaryCL_Ngeo_out(:)
CHARACTER(LEN=5) :: tmpstr
!===================================================================================================================================

ALLOCATE(Xi_Ngeo(0:NGeo))
ALLOCATE(wBary_Ngeo(0:Ngeo))
DO i=0,NGeo
  Xi_Ngeo(i)=-1+REAL(i)*2./REAL(NGeo)
END DO
CALL BarycentricWeights(Ngeo,xi_Ngeo,wBary_Ngeo)

!dirty readin
WRITE(tmpstr,'(I5)')Ngeo
Ngeo_out = GETINT('Ngeo_out',tmpstr)
Ngeo_out = MIN(Ngeo,Ngeo_out) ! should be at maximum Ngeo

ALLOCATE(Xi_Ngeo_out(0:Ngeo_out),wBaryCL_Ngeo_out(0:Ngeo_out))
ALLOCATE(XiCL_Ngeo_out(0:Ngeo_out))
CALL ChebyGaussLobNodesAndWeights(Ngeo_out,XiCL_Ngeo_out)
CALL BarycentricWeights(Ngeo_out,xiCL_Ngeo_out,wBaryCL_Ngeo_out)

!for output from CL to Equidistant points
DO i=0,NGeo_out
  Xi_Ngeo_out(i)=-1+REAL(i)*2./REAL(NGeo_out)
END DO

ALLOCATE(Vdm_CL_EQ_out(0:Ngeo_out,0:Ngeo_out))
CALL InitializeVandermonde(Ngeo_out,Ngeo_out,wBaryCL_Ngeo_out,xiCL_Ngeo_out,xi_Ngeo_Out,Vdm_CL_EQ_out)
! mapping form one-dimensional list [1 ; (Ngeo+1)^3] to tensor-product 0 <= i,j,k <= Ngeo and back
ALLOCATE(HexMap(0:Ngeo,0:Ngeo,0:Ngeo),HexMapInv(3,(Ngeo+1)**3))
l=0
DO k=0,Ngeo ; DO j=0,Ngeo ; DO i=0,Ngeo
  l=l+1
  HexMap(i,j,k)=l
  HexMapInv(:,l)=(/i,j,k/)
END DO ; END DO ; END DO
ALLOCATE(HexMap_out(0:Ngeo_out,0:Ngeo_out,0:Ngeo_out))
l=0
DO k=0,Ngeo_out ; DO j=0,Ngeo_out ; DO i=0,Ngeo_out
  l=l+1
  HexMap_out(i,j,k)=l
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
USE MOD_Mesh_Vars,   ONLY: Ngeo,nElems,nQuads,Xgeo,XgeoQuad
USE MOD_Mesh_Vars,   ONLY: wBary_Ngeo,xi_Ngeo
USE MOD_Mesh_Vars,   ONLY: Ngeo_out,xiCL_Ngeo_out,Vdm_CL_EQ_out
USE MOD_P4EST_Vars,  ONLY: TreeToQuad,QuadCoords,QuadLevel,sIntSize
USE MOD_Basis,       ONLY: LagrangeInterpolationPolys 
USE MOD_ChangeBasis, ONLY: ChangeBasis3D_XYZ ,ChangeBasis3D
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                              :: xi0(3)
REAL                              :: dxi,length
REAL,DIMENSION(0:Ngeo_out,0:Ngeo) :: Vdm_xi,Vdm_eta,Vdm_zeta
REAL                              :: Xgeo_out(3,0:Ngeo_out,0:Ngeo_out,0:Ngeo_out)
INTEGER                           :: StartQuad,EndQuad,nLocalQuads
INTEGER                           :: i,iQuad,iElem 
!===================================================================================================================================
ALLOCATE(XgeoQuad(3,0:Ngeo_out,0:Ngeo_out,0:Ngeo_out,nQuads))

DO iElem=1,nElems
  StartQuad = TreeToQuad(1,iElem)+1
  EndQuad   = TreeToQuad(2,iElem)
  nLocalQuads = TreeToQuad(2,iElem)-TreeToQuad(1,iElem)
  IF((Ngeo.EQ.Ngeo_out).AND.(nLocalQuads.EQ.1))THEN !no refinement in this tree
    XgeoQuad(:,:,:,:,StartQuad)=Xgeo(:,:,:,:,iElem)
  ELSE
    DO iQuad=StartQuad,EndQuad
      ! transform p4est first corner coordinates (integer from 0... intsize) to [-1,1] reference element
      xi0(:)=-1.+2.*REAL(QuadCoords(:,iQuad))*sIntSize
      ! length of each quadrant in integers
      length=2./REAL(2**QuadLevel(iQuad))
      ! Build Vandermonde matrices for each parameter range in xi, eta,zeta
      DO i=0,Ngeo_out
        dxi=0.5*(xiCL_Ngeo_out(i)+1.)*Length
        CALL LagrangeInterpolationPolys(xi0(1) + dxi,Ngeo,xi_Ngeo,wBary_Ngeo,Vdm_xi(i,:)) 
        CALL LagrangeInterpolationPolys(xi0(2) + dxi,Ngeo,xi_Ngeo,wBary_Ngeo,Vdm_eta(i,:)) 
        CALL LagrangeInterpolationPolys(xi0(3) + dxi,Ngeo,xi_Ngeo,wBary_Ngeo,Vdm_zeta(i,:)) 
      END DO
      !interpolate tree HO mapping to quadrant HO mapping (If Ngeo_out < Ngeo: Interpolation error!)
      CALL ChangeBasis3D_XYZ(3,Ngeo,Ngeo_out,Vdm_xi,Vdm_eta,Vdm_zeta,XGeo(:,:,:,:,iElem),Xgeo_out(:,:,:,:))
      ! for output: change points (exactly! Ngeo_out->Ngeo_out) from CL to EQ point distribution
      CALL ChangeBasis3D(3,Ngeo_out,Ngeo_out,Vdm_CL_EQ_out,XGeo_out(:,:,:,:),XgeoQuad(:,:,:,:,iQuad))
    END DO !iQuad=StartQuad,EndQuad
  END IF !nLocalQuads==1
END DO !iElem=1,nElems

IF((Ngeo_out.LT.Ngeo).AND.(Ngeo.GT.1))THEN
  WRITE(*,*)'!!!!!!!!! WARNING: Ngeo_out<Ngeo: non-conforming interfaces are not watertight!!!!'
END IF

END SUBROUTINE BuildHOMesh


SUBROUTINE DeformMesh()
!===================================================================================================================================
! Subroutine to read the mesh from a mesh data file
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Mesh_Vars, ONLY: nElems,XGeo,Ngeo,Deform
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                        :: i,j,k
INTEGER                        :: iElem
REAL                           :: Pi,x(3)
!-----------------------------------------------------------------------------------------------------------------------------------
IF(Deform.EQ.0) RETURN
!deform the mesh
SELECT CASE(Deform)
CASE(1) !sinus -1,1 deformation
  Pi = ACOS(-1.) 
  DO iElem=1,nElems
    DO k=0,NGeo; DO j=0,NGeo; DO i=0,NGeo
      x(:)=Xgeo(:,i,j,k,iElem)
      Xgeo(:,i,j,k,iElem) = x+ 0.1*SIN(Pi*x(1))*SIN(Pi*x(2))*SIN(Pi*x(3))
    END DO; END DO; END DO;
  END DO
CASE DEFAULT
  STOP 'This deform case is not defined'
END SELECT !Deform
END SUBROUTINE DeformMesh


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
INTEGER       :: iElem,iLocSide,iNode
!============================================================================================================================
! Deallocate global variables, needs to go somewhere else later
DO iElem=1,nElems
  DO iLocSide=1,6
    DEALLOCATE(Elems(iElem)%ep%Side(iLocSide)%sp)
  END DO
  DEALLOCATE(Elems(iElem)%ep)
END DO
DEALLOCATE(Elems)
DO iNode=1,nNodes
    DEALLOCATE(Nodes(iNode)%np)
END DO
DEALLOCATE(Nodes)
SDEALLOCATE(XGeo)
SDEALLOCATE(HexMap)
SDEALLOCATE(HexMap_out)
SDEALLOCATE(HexMapInv)
SDEALLOCATE(Xi_NGeo)
SDEALLOCATE(XiCL_NGeo_out)
SDEALLOCATE(wBary_NGeo)
SDEALLOCATE(Vdm_CL_EQ_out)
SDEALLOCATE(BoundaryName)
SDEALLOCATE(BoundaryType)
MeshInitIsDone = .FALSE.
END SUBROUTINE FinalizeMesh

END MODULE MOD_Mesh
