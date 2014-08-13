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
USE MOD_Mesh_Vars,ONLY: NGeo_out,XiCL_NGeo_out,wBaryCL_Ngeo_out,HexMap_out
USE MOD_Mesh_Vars,ONLY: Vdm_01,Vdm_10,Vdm_CL_EQ_out
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
REAL,ALLOCATABLE    :: xi_Ngeo_out(:)
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

!only used in output!!!
ALLOCATE(Vdm_CL_EQ_out(0:Ngeo_out,0:Ngeo_out))
CALL InitializeVandermonde(Ngeo_out,Ngeo_out,wBaryCL_Ngeo_out,xiCL_Ngeo_out,xi_Ngeo_Out,Vdm_CL_EQ_out)


ALLOCATE(Vdm_10(0:Ngeo_out,0:NGeo_out),Vdm_01(0:Ngeo_out,0:NGeo_out))
! change from interval [-1,1] -> [-1,0] Vdm_10 and interval [-1,1]-> [0,1] Vdm__01
CALL InitializeVandermonde(Ngeo_out,Ngeo_out,wBaryCL_Ngeo_out,xiCL_Ngeo_out,-1+0.5*(xiCL_Ngeo_Out+1),Vdm_10)
CALL InitializeVandermonde(Ngeo_out,Ngeo_out,wBaryCL_Ngeo_out,xiCL_Ngeo_out,   0.5*(xiCL_Ngeo_Out+1),Vdm_01)

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
USE MOD_Mesh_Vars,   ONLY: Quads
USE MOD_Mesh_Vars,   ONLY: wBary_Ngeo,xi_Ngeo
USE MOD_Mesh_Vars,   ONLY: Ngeo_out,xiCL_Ngeo_out
USE MOD_P4EST_Vars,  ONLY: TreeToQuad,QuadCoords,QuadLevel,sIntSize
USE MOD_P4EST_Vars,  ONLY: P2H_FaceMap,P_FaceToEdge,P_EdgeToFaces
USE MOD_Basis,       ONLY: LagrangeInterpolationPolys 
USE MOD_ChangeBasis, ONLY: ChangeBasis3D_XYZ
USE MOD_ChangeBasis, ONLY: ChangeBasis2D_XY
USE MOD_Mesh_Vars,   ONLY: Vdm_01,Vdm_10
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
INTEGER                           :: StartQuad,EndQuad,nLocalQuads
INTEGER                           :: i,iQuad,iElem 
!-----------------------------------------------------------------------------------------------------------------------------------
INTEGER                           :: j,k,plus
INTEGER                           :: dir0,dir1,dir2
INTEGER                           :: iEdge,PlocSide1,PlocSide2,dirSide1,dirSide2,pos1,pos2
INTEGER                           :: Pmortar,PlocSide 
INTEGER                           :: MortarType(0:5)
INTEGER                           :: EdgeMarker(0:11)
LOGICAL                           :: MarkForTrans(0:5)
REAL,DIMENSION(0:Ngeo_out,0:Ngeo) :: Vdm_a,Vdm_b
REAL                              :: l_1D(0:Ngeo)
REAL                              :: XgeoSlice(3,0:Ngeo,0:Ngeo)
REAL                              :: XGeoCLFace(3,0:Ngeo_out,0:Ngeo_out,0:5)
REAL                              :: XGeoCLVol(3,0:Ngeo_out,0:Ngeo_out,0:Ngeo_out)
REAL                              :: XgeoCLBigFace(3,0:Ngeo_out,0:Ngeo_out)
REAL                              :: maxdist(0:11)
!===================================================================================================================================
ALLOCATE(XgeoQuad(3,0:Ngeo_out,0:Ngeo_out,0:Ngeo_out,nQuads))

DO iElem=1,nElems
  StartQuad = TreeToQuad(1,iElem)+1
  EndQuad   = TreeToQuad(2,iElem)
  nLocalQuads = TreeToQuad(2,iElem)-TreeToQuad(1,iElem)
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
    CALL ChangeBasis3D_XYZ(3,Ngeo,Ngeo_out,Vdm_xi,Vdm_eta,Vdm_zeta,XGeo(:,:,:,:,iElem),XgeoQuad(:,:,:,:,iQuad))
  END DO !iQuad=StartQuad,EndQuad
END DO !iElem=1,nElems


IF((Ngeo_out.GE.Ngeo).OR.(Ngeo.EQ.1)) THEN
   RETURN
ELSE
  WRITE(*,*)'!!!!!!!!! WARNING: Ngeo_out<Ngeo: non-conforming interfaces are made watertight!!!!'
END IF
! For Ngeo>1 and Ngeo_out < Ngeo, we need to correct the mortar faces!!

!RETURN

DO iElem=1,nElems
  StartQuad = TreeToQuad(1,iElem)+1
  EndQuad   = TreeToQuad(2,iElem)
  nLocalQuads = TreeToQuad(2,iElem)-TreeToQuad(1,iElem)
  DO iQuad=StartQuad,EndQuad
    DO PLocSide=0,5
      MortarType(PlocSide)=Quads(iQuad)%ep%Side(P2H_FaceMap(PLocSide))%sp%MortarType
    END DO !PlocSide
    IF(SUM(ABS(MortarType)).EQ.0) CYCLE !no mortar sides found
    !initialize Face Data (equidistant point distribution )
    XGeoCLVol=XGeoQuad(:,:,:,:,iQuad)
    XGeoCLFace(:,:,:,0)=XGeoCLVol(:,       0,:,:)     
    XGeoCLFace(:,:,:,1)=XGeoCLVol(:,Ngeo_out,:,:)     
    XGeoCLFace(:,:,:,2)=XGeoCLVol(:,:,       0,:)     
    XGeoCLFace(:,:,:,3)=XGeoCLVol(:,:,Ngeo_out,:)     
    XGeoCLFace(:,:,:,4)=XGeoCLVol(:,:,:,       0)     
    XGeoCLFace(:,:,:,5)=XGeoCLVol(:,:,:,Ngeo_out)     
    !Mark already sides which are big mortars
    MarkForTrans=(MortarType.GT.0)  
    !initialize EdgeMarker
    EdgeMarker=0
    DO PLocSide=0,5
      IF(MortarType(PlocSide).LT.0)THEN  !small mortar face:
        ! transform p4est first corner coordinates (integer from 0... intsize) to [-1,1] reference element
        xi0(:)=-1.+2.*REAL(QuadCoords(:,iQuad))*sIntSize
        ! length of each quadrant in integers
        length=2./REAL(2**(QuadLevel(iQuad)))
        PMortar=-MortarType(PlocSide)-1
        !edgemarker<10: 1 edge: PlocSide=EdgeMarker-1, 
        !edgemarker>10: 2edges:PlocSide1=MOD(EdgeMarker,10)-1, PlocSide2=EdgeMarker-10*PlocSide1-1
        EdgeMarker(P_FaceToEdge(:,PlocSide))=EdgeMarker(P_FaceToEdge(:,PlocSide))*10 + PlocSide+1 
        SELECT CASE(PlocSide)
        CASE(0,1) !ximinus,xiplus
          MarkForTrans(2:5)=.TRUE.
          plus=PlocSide ! plus=0: minus side, plus=1: plus Side
          dir0=1
          dir1=2
          dir2=3
        CASE(2,3) !etaminus,etaplus
          MarkForTrans(0:1)=.TRUE.
          MarkForTrans(4:5)=.TRUE.
          plus=PlocSide-2 ! plus=0: minus side, plus=1: plus Side
          dir0=2
          dir1=1
          dir2=3
        CASE(4,5) !zetaminus,zetaplus
          MarkForTrans(0:3)=.TRUE.
          plus=PlocSide-4 !  plus=0: minus side, plus=1: plus Side
          dir0=3
          dir1=1
          dir2=2
        END SELECT !PlocSide
        !position of origin of mortar side in tree
        xi0(dir0)=xi0(dir0)+plus*length
        !position of origin of big (neighbor) mortar side in tree
        SELECT CASE(PMortar)
        CASE(0)! lower left
          !xi0(dir1),xi0(dir2) same
        CASE(1)! lower right
          xi0(dir1)=xi0(dir1)-length
        CASE(2)! upper left 
          xi0(dir2)=xi0(dir2)-length
        CASE(3)! upper right 
          xi0(dir1)=xi0(dir1)-length
          xi0(dir2)=xi0(dir2)-length
        END SELECT !Pmortar
        !extract slice from tree:
        CALL LagrangeInterpolationPolys(xi0(dir0),Ngeo,xi_Ngeo,wBary_Ngeo,l_1D(:)) 
        SELECT CASE(dir0)
        CASE(1)
          XgeoSlice=0.
          DO i=0,Ngeo
            XgeoSlice(:,:,:)=XgeoSlice(:,:,:)+l_1D(i)*Xgeo(:,i,:,:,iElem)
          END DO 
        CASE(2)
          XgeoSlice=0.
          DO j=0,Ngeo
            XgeoSlice(:,:,:)=XgeoSlice(:,:,:)+l_1D(j)*Xgeo(:,:,j,:,iElem)
          END DO 
        CASE(3)
          XgeoSlice=0.
          DO k=0,Ngeo
            XgeoSlice(:,:,:)=XgeoSlice(:,:,:)+l_1D(k)*Xgeo(:,:,:,k,iElem)
          END DO 
        END SELECT !dir0
        !interpolate slice of tree to quadrant big neighbor face (length*2)  EQ (Ngeo) -> CL (Ngeo_out)
        !   build Vdm for interpolation EQ Ngeo -> CL Ngeo_out 
        DO i=0,Ngeo_out
          dxi=(xiCL_Ngeo_out(i)+1.)*Length !large element side (length*2)!!
          CALL LagrangeInterpolationPolys(xi0(dir1) + dxi,Ngeo,xi_Ngeo,wBary_Ngeo,Vdm_a(i,:)) 
          CALL LagrangeInterpolationPolys(xi0(dir2) + dxi,Ngeo,xi_Ngeo,wBary_Ngeo,Vdm_b(i,:)) 
        END DO
        CALL ChangeBasis2D_XY(3,Ngeo,Ngeo_out,Vdm_a,Vdm_b,XGeoSlice(:,:,:),XgeoCLBigFace(:,:,:))
        !transfinite face remap, because big mortar faces will be transfinite too!!
        CALL TransFace(3,Ngeo_out,xiCL_Ngeo_out,XgeoCLBigFace(:,:,:)) 
        ! interplation to small face 
        SELECT CASE(PMortar)
        CASE(0)! lower left [-1,1]^2 -> [-1,0]x[-1,0]
          CALL ChangeBasis2D_XY(3,Ngeo_out,Ngeo_out,Vdm_10,Vdm_10,XgeoCLBigFace(:,:,:),XGeoCLFace(:,:,:,PlocSide))
        CASE(1)! lower right [-1,1]^2 -> [1,0]x[-1,0]
          CALL ChangeBasis2D_XY(3,Ngeo_out,Ngeo_out,Vdm_01,Vdm_10,XgeoCLBigFace(:,:,:),XGeoCLFace(:,:,:,PlocSide))
        CASE(2)! upper left [-1,1]^2 -> [-1,0]x[0,1]
          CALL ChangeBasis2D_XY(3,Ngeo_out,Ngeo_out,Vdm_10,Vdm_01,XgeoCLBigFace(:,:,:),XGeoCLFace(:,:,:,PlocSide))
        CASE(3)! upper right [-1,1]^2 -> [0,1]x[0,1] 
          CALL ChangeBasis2D_XY(3,Ngeo_out,Ngeo_out,Vdm_01,Vdm_01,XgeoCLBigFace(:,:,:),XGeoCLFace(:,:,:,PlocSide))
        END SELECT !Pmortar
        !copy modified faces into volume (produces unique edges)
        SELECT CASE(PlocSide)
        CASE(0)
          XGeoCLVol(:,       0,:,:)=XGeoCLFace(:,:,:,0)     
        CASE(1)
          XGeoCLVol(:,Ngeo_out,:,:)=XGeoCLFace(:,:,:,1)     
        CASE(2)
          XGeoCLVol(:,:,       0,:)=XGeoCLFace(:,:,:,2)     
        CASE(3)
          XGeoCLVol(:,:,Ngeo_out,:)=XGeoCLFace(:,:,:,3)     
        CASE(4)
          XGeoCLVol(:,:,:,       0)=XGeoCLFace(:,:,:,4)     
        CASE(5)
          XGeoCLVol(:,:,:,Ngeo_out)=XGeoCLFace(:,:,:,5)     
        END SELECT !PlocSide
      END IF !smallmortarSide
    END DO !PlocSide=0,5


    !Check if edges are really unique
    maxdist=-1.
    DO iEdge=0,11
      Plocside1=P_EdgeToFaces(1,iEdge)              ! first adjacent local side
      dirside1 =P_EdgeToFaces(2,iEdge)              ! 0: i, 1: j
      pos1     =P_EdgeToFaces(3,iEdge)*NGeo_out     ! 0: 0, 1: N
      Plocside2=P_EdgeToFaces(4,iEdge)              ! second adjacent local side
      dirside2 =P_EdgeToFaces(5,iEdge)              ! 0: i, 1: j
      pos2     =P_EdgeToFaces(6,iEdge)*NGeo_out     ! 0: 0, 1: N
      IF(EdgeMarker(iEdge).GT.10)THEN
        IF(dirside1.EQ.0)THEN !first side in i
          IF(dirside2.EQ.0)THEN !second side in i
             maxdist(iEdge)=MAXVAL( (XGeoCLFace(1,pos1,:,PlocSide1)-XGeoCLFace(1,pos2,:,PlocSide2))**2 &    
                                   +(XGeoCLFace(2,pos1,:,PlocSide1)-XGeoCLFace(2,pos2,:,PlocSide2))**2 &
                                   +(XGeoCLFace(3,pos1,:,PlocSide1)-XGeoCLFace(3,pos2,:,PlocSide2))**2)
          ELSE !second side in j
             maxdist(iEdge)=MAXVAL( (XGeoCLFace(1,pos1,:,PlocSide1)-XGeoCLFace(1,:,pos2,PlocSide2))**2 &    
                                   +(XGeoCLFace(2,pos1,:,PlocSide1)-XGeoCLFace(2,:,pos2,PlocSide2))**2 &
                                   +(XGeoCLFace(3,pos1,:,PlocSide1)-XGeoCLFace(3,:,pos2,PlocSide2))**2)
          END IF
        ELSE !dirside1=1 ->  first side in j
          IF(dirside2.EQ.0)THEN !second side in i
             maxdist(iEdge)=MAXVAL( (XGeoCLFace(1,:,pos1,PlocSide1)-XGeoCLFace(1,pos2,:,PlocSide2))**2 &    
                                   +(XGeoCLFace(2,:,pos1,PlocSide1)-XGeoCLFace(2,pos2,:,PlocSide2))**2 &
                                   +(XGeoCLFace(3,:,pos1,PlocSide1)-XGeoCLFace(3,pos2,:,PlocSide2))**2)
          ELSE !second side in j                       
             maxdist(iEdge)=MAXVAL( (XGeoCLFace(1,:,pos1,PlocSide1)-XGeoCLFace(1,:,pos2,PlocSide2))**2 &    
                                   +(XGeoCLFace(2,:,pos1,PlocSide1)-XGeoCLFace(2,:,pos2,PlocSide2))**2 &
                                   +(XGeoCLFace(3,:,pos1,PlocSide1)-XGeoCLFace(3,:,pos2,PlocSide2))**2)
          END IF
        END IF
      END IF
    END DO !iEdge
    
    IF(ANY(maxDist(:).GT.1.0E-15))THEN
      WRITE(*,*)'WARNING!!! PROBLEMS WITH Ngeo_out'
      WRITE(*,'(A20,12E11.2)')'maxdist',maxdist
      WRITE(*,'(A20,12I11)')'EdgeMarker',EdgeMarker
      STOP
    END IF
   
    DO PLocSide=0,5
      IF(MortarType(PlocSide).GE.0)THEN  !remaining faces: tranfinite mapping
        IF(MarkForTrans(PlocSide)) THEN
          !transfinite face remap!!
          SELECT CASE(PlocSide)
          CASE(0)
            CALL TransFace(3,Ngeo_out,xiCL_Ngeo_out,XGeoCLVol(:,       0,:,:))     
          CASE(1)
            CALL TransFace(3,Ngeo_out,xiCL_Ngeo_out,XGeoCLVol(:,Ngeo_out,:,:))     
          CASE(2)
            CALL TransFace(3,Ngeo_out,xiCL_Ngeo_out,XGeoCLVol(:,:,       0,:))     
          CASE(3)
            CALL TransFace(3,Ngeo_out,xiCL_Ngeo_out,XGeoCLVol(:,:,Ngeo_out,:))     
          CASE(4)
            CALL TransFace(3,Ngeo_out,xiCL_Ngeo_out,XGeoCLVol(:,:,:,       0))     
          CASE(5)
            CALL TransFace(3,Ngeo_out,xiCL_Ngeo_out,XGeoCLVol(:,:,:,Ngeo_out))     
          END SELECT !PlocSide
        END IF
      END IF !mortarSide
    END DO !PlocSide=0,5
    CALL TransVol(3,Ngeo_out,XiCL_Ngeo_out,XGeoCLVol,XgeoQuad(:,:,:,:,iQuad))
  END DO !iQuad=StartQuad,EndQuad
END DO !iElem=1,nElems

END SUBROUTINE BuildHOMesh


SUBROUTINE TransFace(dim1,N_in,xi_in,Face)
!===================================================================================================================================
! Transfinite mapping edge faces -> face, replace inner points only 
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)         :: dim1           ! size of leading dimension
INTEGER,INTENT(IN)         :: N_in           !polynomial degree
REAL,INTENT(IN)            :: xi_in(0:N_in)  !node positions in parameter space [-1,1]
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(INOUT)         :: Face(1:dim1,0:N_in,0:N_in) !Face data, inner points will be overwritten
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                        :: i,j
REAL                           :: xi(0:N_in),xim(0:N_in)
!-----------------------------------------------------------------------------------------------------------------------------------
IF(N_in.EQ.1) RETURN
xi=0.5*(xi_in+1.)  ![-1,1] ->  0...1
xim=1.-xi          ! 1 ... 0
DO j=1,N_in-1
  DO i=1,N_in-1
     Face(:,i,j)=   Face(:,i,0)*xim(j) + Face(:,i,N_in)*xi(j)                  &
                   +Face(:,0,j)*xim(i) + Face(:,N_in,j)*xi(i)                  &
                  -( (Face(:,0,   0)*xim(i) + Face(:,N_in,   0)*xi(i) )*xim(j) &
                    +(Face(:,0,N_in)*xim(j) + Face(:,N_in,N_in)*xi(j) )*xim(i))
  END DO !i
END DO !j
END SUBROUTINE TransFace


SUBROUTINE TransVol(dim1,N_in,xi_in,Vol_in,Vol)
!===================================================================================================================================
! Transfinite mapping 6 faces -> volume, replace face points and move inner points by a transfinite interpolation of face deformation
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)         :: dim1           ! size of leading dimension
INTEGER,INTENT(IN)         :: N_in           !polynomial degree
REAL,INTENT(IN)            :: xi_in(0:N_in)  !node positions in parameter space [-1,1]
REAL,INTENT(IN)            :: Vol_in(1:dim1,0:N_in,0:N_in,0:N_in) ! volume data, but only face data is used
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(INOUT)         :: Vol(1:dim1,0:N_in,0:N_in,0:N_in) !unchanged volume data, face points will be replaced by vol_in, inner points will be moved by transfinite deformation
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                    :: i,j,k
REAL                       :: xi(0:N_in),xim(0:N_in)
REAL                       :: dVol(1:dim1,0:N_in,0:N_in,0:N_in) !Face data, inner points will be overwritten
!-----------------------------------------------------------------------------------------------------------------------------------
IF(N_in.EQ.1) THEN
  Vol=Vol_in
  RETURN ! no inner nodes
END IF
!now N_in>1
!difference between new (Vol_in) and old (Vold)
dVol=Vol_in-Vol !only face data will be used

xi=0.5*(xi_in+1.)  ![-1,1] ->  0...1
xim=1.-xi          ! 1 ... 0
!overwrite inner nodes of dVol by transfinite blending of its faces
DO k=1,N_in-1
  DO j=1,N_in-1
    DO i=1,N_in-1
       dVol(:,i,j,k)=  dVol(:,i,j,0)*xim(k)+dVol(:,i,j,N_in)*xi(k)      & !faces
                      +dVol(:,i,0,k)*xim(j)+dVol(:,i,N_in,k)*xi(j)      &
                      +dVol(:,0,j,k)*xim(i)+dVol(:,N_in,j,k)*xi(i)      &
                      -(  dVol(:,   0,   0,   k)*xim(i)*xim(j)          & !-edges
                        + dVol(:,   0,N_in,   k)*xim(i)* xi(j)          &
                        + dVol(:,   0,   j,   0)*xim(i)       *xim(k)   &
                        + dVol(:,   0,   j,N_in)*xim(i)       * xi(k)   &
                        + dVol(:,N_in,   0,   k)* xi(i)*xim(j)          &
                        + dVol(:,N_in,N_in,   k)* xi(i)* xi(j)          &
                        + dVol(:,N_in,   j,   0)* xi(i)       *xim(k)   &
                        + dVol(:,N_in,   j,N_in)* xi(i)       * xi(k)   &
                        + dVol(:,   i,   0,   0)       *xim(j)*xim(k)   &
                        + dVol(:,   i,   0,N_in)       *xim(j)* xi(k)   &
                        + dVol(:,   i,N_in,   0)       * xi(j)*xim(k)   &
                        + dVol(:,   i,N_in,N_in)       * xi(j)* xi(k) ) &
                      +(  dVol(:,   0,   0,   0)*xim(i)*xim(j)*xim(k)   & ! + corners
                        + dVol(:,   0,   0,N_in)*xim(i)*xim(j)* xi(k)   &
                        + dVol(:,   0,N_in,   0)*xim(i)* xi(j)*xim(k)   &
                        + dVol(:,   0,N_in,N_in)*xim(i)* xi(j)* xi(k)   &
                        + dVol(:,N_in,   0,   0)* xi(i)*xim(j)*xim(k)   &
                        + dVol(:,N_in,   0,N_in)* xi(i)*xim(j)* xi(k)   &
                        + dVol(:,N_in,N_in,   0)* xi(i)* xi(j)*xim(k)   &
                        + dVol(:,N_in,N_in,N_in)* xi(i)* xi(j)* xi(k) )
    END DO !i
  END DO !j
END DO !k
!add deformation to volume mapping
Vol=Vol+dVol  !vol=vol_in is recovered on the faces, inner nodes will be moved by dVol 
END SUBROUTINE TransVol


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
