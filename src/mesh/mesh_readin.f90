#include "hopest_f.h"

MODULE MOD_Mesh_ReadIn
!===================================================================================================================================
! Add comments please!
!===================================================================================================================================
! MODULES
USE MOD_HDF5_Input
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------

! Public Part ----------------------------------------------------------------------------------------------------------------------
INTERFACE ReadMeshFromHDF5
  MODULE PROCEDURE ReadMeshFromHDF5
END INTERFACE

INTERFACE ReadMeshFromHDF5nobuildp4est
  MODULE PROCEDURE ReadMeshFromHDF5nobuildp4est
END INTERFACE

INTERFACE ReadGeoFromHDF5
  MODULE PROCEDURE ReadGeoFromHDF5
END INTERFACE

PUBLIC::ReadMeshFromHDF5
PUBLIC::ReadMeshFromHDF5nobuildp4est
PUBLIC::ReadGeoFromHDF5
!===================================================================================================================================

CONTAINS

SUBROUTINE ReadBCs()
!===================================================================================================================================
! Read boundary conditions from data file
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Mesh_Vars,ONLY:BoundaryName,BoundaryType,nBCs,nUserBCs
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER, ALLOCATABLE           :: BCMapping(:),BCType(:,:)
CHARACTER(LEN=255), ALLOCATABLE:: BCNames(:)
INTEGER                        :: iBC,iUserBC
INTEGER                        :: Offset=0 ! Every process reads all BCs
!===================================================================================================================================
! Read boundary names from data file
CALL GetDataSize(File_ID,'BCNames',nDims,HSize)
nBCs=HSize(1)
DEALLOCATE(HSize)
ALLOCATE(BCNames(nBCs))
ALLOCATE(BCMapping(nBCs))
CALL ReadArray('BCNames',1,(/nBCs/),Offset,1,StrArray=BCNames)  ! Type is a dummy type only
! User may have redefined boundaries in the ini file. So we have to create mappings for the boundaries.
BCMapping=0
IF(nUserBCs .GT. 0)THEN
  DO iBC=1,nBCs
    DO iUserBC=1,nUserBCs
      IF(INDEX(TRIM(BCNames(iBC)),TRIM(BoundaryName(iUserBC))) .NE.0) BCMapping(iBC)=iUserBC
    END DO
  END DO
END IF

! Read boundary types from data file
CALL GetDataSize(File_ID,'BCType',nDims,HSize)
IF(HSize(1).NE.nBCs) STOP 'Problem in readBC'
DEALLOCATE(HSize)
ALLOCATE(BCType(nBCs,4))
offset=0
CALL ReadArray('BCType',2,(/nBCs,4/),Offset,1,IntegerArray=BCType)
! Now apply boundary mappings
IF(nUserBCs .GT. 0)THEN
  DO iBC=1,nBCs
    IF(BCMapping(iBC) .NE. 0)THEN
      SWRITE(Unit_StdOut,'(A,A)')    ' |     Boundary in HDF file found |  ',TRIM(BCNames(iBC))
      SWRITE(Unit_StdOut,'(A,I2,I2)')' |                            was | ',BCType(iBC,1),BCType(iBC,3)
      SWRITE(Unit_StdOut,'(A,I2,I2)')' |                      is set to | ',BoundaryType(BCMapping(iBC),1:2)
      BCType(iBC,1) = BoundaryType(BCMapping(iBC),1)
      BCType(iBC,3) = BoundaryType(BCMapping(iBC),2)
    END IF
  END DO
END IF
IF(ALLOCATED(BoundaryName)) DEALLOCATE(BoundaryName)
IF(ALLOCATED(BoundaryType)) DEALLOCATE(BoundaryType)
ALLOCATE(BoundaryName(nBCs))
ALLOCATE(BoundaryType(nBCs,BC_SIZE))
BoundaryName = BCNames
BoundaryType(:,BC_TYPE)  = BCType(:,BC_TYPE)  
BoundaryType(:,BC_CURVED)= BCType(:,BC_CURVED)
BoundaryType(:,BC_STATE) = BCType(:,BC_STATE) 
BoundaryType(:,BC_ALPHA) = BCType(:,BC_ALPHA) 
SWRITE(UNIT_StdOut,'(132("."))')
SWRITE(Unit_StdOut,'(A,A16,A20,A10,A10,A10,A10)')'BOUNDARY CONDITIONS','|','Name','Type','Curved','State','Alpha'
DO iBC=1,nBCs
  SWRITE(*,'(A,A33,A20,I10,I10,I10,I10)')' |','|',TRIM(BoundaryName(iBC)),BoundaryType(iBC,:)
END DO
SWRITE(UNIT_StdOut,'(132("."))')
DEALLOCATE(BCNames,BCType,BCMapping)
END SUBROUTINE ReadBCs


SUBROUTINE ReadMeshHeader()
!===================================================================================================================================
! Subroutine to read the mesh from a mesh data file
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Mesh_Vars,ONLY: NGeo,useCurveds,nGlobalElems
USE MOD_Mesh,     ONLY: SetCurvedInfo
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER :: BoundaryOrder_mesh
!===================================================================================================================================
CALL GetDataSize(File_ID,'ElemInfo',nDims,HSize)
nGlobalElems=HSize(1) !global number of elements
DEALLOCATE(HSize)

CALL ReadAttribute(File_ID,'BoundaryOrder',1,IntegerScalar=BoundaryOrder_mesh)
NGeo = BoundaryOrder_mesh-1
CALL ReadAttribute(File_ID,'CurvedFound',1,LogicalScalar=useCurveds)

CALL readBCs()
CALL SetCurvedInfo()

END SUBROUTINE ReadMeshHeader


SUBROUTINE ReadMeshFromHDF5nobuildp4est(FileString)
!===================================================================================================================================
! Subroutine to read the mesh from a mesh data file and build p4est_connectivity
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Mesh_Vars
USE MOD_P4EST_Vars,    ONLY: connectivity,p4est,H2P_VertexMap,H2P_FaceMap,geom
USE MOD_P4EST_Binding, ONLY: p4_connectivity_treevertex,p4_build_p4est
USE MOD_P4EST,         ONLY: getHFlip
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
CHARACTER(LEN=*),INTENT(IN)  :: FileString
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                        :: i,j,k,l
INTEGER                        :: BCindex
INTEGER                        :: iElem,ElemID
INTEGER                        :: iNode,jNode,NodeID,SideID
INTEGER                        :: iLocSide,jLocSide
INTEGER                        :: iSide
INTEGER                        :: nCurvedNodes_loc
LOGICAL                        :: oriented
INTEGER                        :: nPeriodicSides 
LOGICAL                        :: fileExists
LOGICAL                        :: doConnection
TYPE(tElem),POINTER            :: aElem
TYPE(tSide),POINTER            :: aSide,bSide
TYPE(tNode),POINTER            :: aNode
TYPE(tNodePtr),POINTER         :: ElemCurvedNode(:,:)
INTEGER,ALLOCATABLE            :: ElemInfo(:,:),SideInfo(:,:),NodeInfo(:)
REAL,ALLOCATABLE               :: NodeCoords(:,:)
                               
INTEGER                        :: nNodeIDs,nSideIDs
! p4est interface
INTEGER                        :: num_vertices
INTEGER                        :: num_trees
INTEGER                        :: num_periodics,iPeriodic,PFlip,HFlip,HFlip_test
INTEGER,ALLOCATABLE            :: tree_to_vertex(:,:)
REAL,ALLOCATABLE               :: vertices(:,:)
INTEGER,ALLOCATABLE            :: JoinFaces(:,:)
!===================================================================================================================================
IF(MESHInitIsDone) RETURN
INQUIRE (FILE=TRIM(FileString), EXIST=fileExists)
IF(.NOT.FileExists)  &
    CALL abort(__STAMP__, &
       'readMesh from data file "'//TRIM(FileString)//'" does not exist')

SWRITE(UNIT_stdOut,'(A)')'READ MESH FROM DATA FILE "'//TRIM(FileString)//'" ...'
SWRITE(UNIT_StdOut,'(132("-"))')
! Open data file
CALL OpenDataFile(FileString,create=.FALSE.,single=.FALSE.)

CALL ReadMeshHeader()

!----------------------------------------------------------------------------------------------------------------------------
!                              ELEMENTS
!----------------------------------------------------------------------------------------------------------------------------

nElems=nGlobalElems   !local number of Elements 
!read local ElemInfo from data file
ALLOCATE(ElemInfo(1:nElems,ELEM_InfoSize))
CALL ReadArray('ElemInfo',2,(/nElems,ELEM_InfoSize/),0,1,IntegerArray=ElemInfo)

ALLOCATE(Elems(1:nElems))

DO iElem=1,nElems
  iSide=ElemInfo(iElem,ELEM_FirstSideInd) !first index -1 in Sideinfo
  iNode=ElemInfo(iElem,ELEM_FirstNodeInd) !first index -1 in NodeInfo
  Elems(iElem)%ep=>GETNEWELEM()
  aElem=>Elems(iElem)%ep
  aElem%Ind    = iElem
  aElem%Type   = ElemInfo(iElem,ELEM_Type)
  aElem%Zone   = ElemInfo(iElem,ELEM_Zone)
END DO

!----------------------------------------------------------------------------------------------------------------------------
!                              NODES
!----------------------------------------------------------------------------------------------------------------------------

!read local Node Info from data file 
nNodeIDs=ElemInfo(nElems,ELEM_LastNodeInd)-ElemInfo(1,ELEM_FirstNodeInd)
ALLOCATE(NodeInfo(1:nNodeIDs))
CALL ReadArray('NodeInfo',1,(/nNodeIDs/),0,1,IntegerArray=NodeInfo)

ALLOCATE(ElemCurvedNode(nCurvedNodes,nElems))

CALL GetDataSize(File_ID,'NodeCoords',nDims,HSize)
nNodes=HSize(1) !global number of unique nodes
DEALLOCATE(HSize)

ALLOCATE(Nodes(1:nNodes)) ! pointer list, entry is known by NodeCoords
DO iNode=1,nNodes
  NULLIFY(Nodes(iNode)%np)
END DO
!assign nodes 
DO iElem=1,nElems
  aElem=>Elems(iElem)%ep
  iNode=ElemInfo(iElem,ELEM_FirstNodeInd) !first index -1 in NodeInfo
  DO jNode=1,8
    iNode=iNode+1
    NodeID=ABS(NodeInfo(iNode))     !global, unique NodeID
    IF(.NOT.ASSOCIATED(Nodes(NodeID)%np))THEN
      ALLOCATE(Nodes(NodeID)%np) 
      Nodes(NodeID)%np%ind=NodeID 
    END IF
    aElem%Node(jNode)%np=>Nodes(NodeID)%np
  END DO
  CALL createSides(aElem)
  IF(NGeo.GT.1)THEN
    nCurvedNodes_loc = ElemInfo(iElem,ELEM_LastNodeInd) - ElemInfo(iElem,ELEM_FirstNodeInd) - 14 ! corner + oriented nodes
    IF(nCurvedNodes.NE.nCurvedNodes_loc) &
      CALL abort(__STAMP__, &
           'Wrong number of curved nodes for hexahedra.')
    DO i=1,nCurvedNodes
      iNode=iNode+1
      NodeID=NodeInfo(iNode) !first oriented corner node
      IF(.NOT.ASSOCIATED(Nodes(NodeID)%np))THEN
        ALLOCATE(Nodes(NodeID)%np)
        Nodes(NodeID)%np%ind=NodeID 
      END IF
      ElemCurvedNode(i,iElem)%np=>Nodes(NodeID)%np
    END DO
  END IF
END DO

!----------------------------------------------------------------------------------------------------------------------------
!                              SIDES
!----------------------------------------------------------------------------------------------------------------------------

nSideIDs=ElemInfo(nElems,ELEM_LastSideInd)-ElemInfo(1,ELEM_FirstSideInd)
!read local SideInfo from data file 
ALLOCATE(SideInfo(1:nSideIDs,SIDE_InfoSize))
CALL ReadArray('SideInfo',2,(/nSideIDs,SIDE_InfoSize/),0,1,IntegerArray=SideInfo)

DO iElem=1,nElems
  aElem=>Elems(iElem)%ep
  iNode=ElemInfo(iElem,ELEM_LastNodeInd) !first index -1 in NodeInfo
  iNode=iNode-6
  iSide=ElemInfo(iElem,ELEM_FirstSideInd) !first index -1 in Sideinfo
  !build up sides of the element using element Nodes and CGNS standard
  ! assign flip
  DO iLocSide=1,6
    aSide=>aElem%Side(iLocSide)%sp
    iSide=iSide+1

    ElemID=SideInfo(iSide,SIDE_nbElemID) !IF nbElemID <0, this marks a mortar master side. 
                                         ! The number (-1,-2,-3) is the Type of mortar
    IF(ElemID.LT.0)THEN ! mortar Sides attached!
      CALL abort(__STAMP__, &
           'Only conforming meshes in readin.')
    END IF
   
    aSide%Elem=>aElem
    oriented=(Sideinfo(iSide,SIDE_ID).GT.0)
    
    aSide%Ind=ABS(SideInfo(iSide,SIDE_ID))
    iNode=iNode+1
    NodeID=NodeInfo(iNode) !first oriented corner node
    IF(oriented)THEN !oriented side
      aSide%flip=0
    ELSE !not oriented
      DO jNode=1,4
        IF(aSide%Node(jNode)%np%ind.EQ.ABS(NodeID)) EXIT
      END DO
      IF(jNode.GT.4) STOP 'NodeID doesnt belong to side'
      aSide%flip=jNode
    END IF

  END DO !i=1,locnSides
END DO !iElem

 
! build up side connection 
DO iElem=1,nElems
  aElem=>Elems(iElem)%ep
  iSide=ElemInfo(iElem,ELEM_FirstSideInd) !first index -1 in Sideinfo
  DO iLocSide=1,6
    aSide=>aElem%Side(iLocSide)%sp
    iSide=iSide+1

    sideID  = ABS(SideInfo(iSide,SIDE_ID))
    elemID  = SideInfo(iSide,SIDE_nbElemID)
    BCindex = SideInfo(iSide,SIDE_BCID)

    doConnection=.TRUE. ! for periodic sides if BC is reassigned as non periodic
    IF(BCindex.NE.0)THEN !BC
      aSide%BCindex = BCindex
      IF(BoundaryType(aSide%BCindex,BC_TYPE).NE.1)THEN ! Reassignement from periodic to non-periodic
        doConnection=.FALSE.
        aSide%flip  =0
        elemID            = 0
      END IF
    ELSE
      aSide%BCindex = 0
    END IF

    IF(.NOT.ASSOCIATED(aSide%connection))THEN
      IF((elemID.NE.0).AND.doConnection)THEN !connection 
        IF((elemID.LE.nElems).AND.(elemID.GE.1))THEN !local connection
          DO jLocSide=1,6
            bSide=>Elems(elemID)%ep%Side(jLocSide)%sp
            IF(bSide%ind.EQ.aSide%ind)THEN
              aSide%connection=>bSide
              bSide%connection=>aSide
              EXIT
            END IF !bSide%ind.EQ.aSide%ind
          END DO !jLocSide
        ELSE !MPI connection
          CALL abort(__STAMP__, &
            ' elemID of neighbor not in global Elem list ')
        END IF
      END IF
    END IF !connection associated
  END DO !iLocSide 
END DO !iElem

DEALLOCATE(ElemInfo,SideInfo,NodeInfo)

! get physical coordinates

ALLOCATE(NodeCoords(nNodes,3))

CALL ReadArray('NodeCoords',2,(/nNodes,3/),0,1,RealArray=NodeCoords)

ALLOCATE(Xgeo(1:3,0:Ngeo,0:Ngeo,0:Ngeo,nElems))
IF(Ngeo.EQ.1)THEN !use the corner nodes
  DO iElem=1,nElems
    aElem=>Elems(iElem)%ep
    Xgeo(:,0,0,0,iElem)=NodeCoords(aElem%Node(1)%np%ind,:)
    Xgeo(:,1,0,0,iElem)=NodeCoords(aElem%Node(2)%np%ind,:)
    Xgeo(:,1,1,0,iElem)=NodeCoords(aElem%Node(3)%np%ind,:)
    Xgeo(:,0,1,0,iElem)=NodeCoords(aElem%Node(4)%np%ind,:)
    Xgeo(:,0,0,1,iElem)=NodeCoords(aElem%Node(5)%np%ind,:)
    Xgeo(:,1,0,1,iElem)=NodeCoords(aElem%Node(6)%np%ind,:)
    Xgeo(:,1,1,1,iElem)=NodeCoords(aElem%Node(7)%np%ind,:)
    Xgeo(:,0,1,1,iElem)=NodeCoords(aElem%Node(8)%np%ind,:)
 END DO !iElem=1,nElems
ELSE
  DO iElem=1,nElems
    aElem=>Elems(iElem)%ep
    l=0
    DO k=0,Ngeo; DO j=0,Ngeo; DO i=0,Ngeo
      l=l+1
      Xgeo(:,i,j,k,iElem)=NodeCoords(ElemCurvedNode(l,iElem)%np%ind,:)
    END DO ; END DO ; END DO 
 END DO !iElem=1,nElems
END IF

CALL CloseDataFile() 

DEALLOCATE(ElemCurvedNode)


! P4est MESH connectivity (should be replaced by connectivity information ?)

! needs unique corner nodes for mesh connectivity
DO iNode=1,nNodes
  Nodes(iNode)%np%tmp=-1
END DO
num_vertices=0
DO iElem=1,nElems
  aElem=>Elems(iElem)%ep
  DO iNode=1,8
    aNode=>aElem%Node(iNode)%np
    IF(aNode%tmp.EQ.-1)THEN
      num_vertices=num_vertices+1
      aElem%Node(iNode)%np%tmp=num_vertices
    END IF
  END DO
END DO !iElem

ALLOCATE(Vertices(3,num_vertices))
DO iNode=1,nNodes
  aNode=>Nodes(iNode)%np
  IF(aNode%tmp.GT.0)THEN
    Vertices(:,aNode%tmp)=NodeCoords(aNode%ind,:)
  END IF
END DO


DEALLOCATE(NodeCoords)

num_trees=nElems
ALLOCATE(tree_to_vertex(8,num_trees))
DO iElem=1,nElems
  aElem=>Elems(iElem)%ep
  DO iNode=1,8
    tree_to_vertex(iNode,iElem)=aElem%Node(H2P_VertexMap(iNode)+1)%np%tmp-1
  END DO
END DO

!periodic Boundaries
num_periodics=0
DO iElem=1,nElems
  aElem=>Elems(iElem)%ep
  DO iLocSide=1,6
    aSide=>aElem%Side(iLocSide)%sp
    IF(aSide%BCIndex.EQ.0) CYCLE ! NO Boundary Condition
    IF((BoundaryType(aSide%BCIndex,BC_TYPE).EQ.1).AND.(aSide%flip.EQ.0))THEN !periodic side 
      num_periodics=num_periodics+1
    END IF
  END DO !iLocSide
END DO !iElem


IF(num_periodics.GT.0) THEN
  ALLOCATE(JoinFaces(5,num_periodics))
  DO iElem=1,nElems
    aElem=>Elems(iElem)%ep
    DO iLocSide=1,6
      aElem%Side(iLocSide)%sp%tmp=H2P_FaceMap(iLocSide)  !local Face ID in p4est
    END DO
  END DO
  
  iperiodic=0
  DO iElem=1,nElems
    aElem=>Elems(iElem)%ep
    DO iLocSide=1,6
      aSide=>aElem%Side(iLocSide)%sp
      IF(aSide%BCIndex.EQ.0) CYCLE ! NO Boundary Condition
      IF((BoundaryType(aSide%BCIndex,BC_TYPE).EQ.1).AND.(aSide%flip.EQ.0))THEN !periodic side 
        HFlip=aSide%connection%flip
        iperiodic=iperiodic+1
        bSide=>aSide%connection
        IF(aSide%tmp.GT.bSide%tmp)THEN
          aSide=>aSide%connection
        END IF
        bSide=>aSide%connection
        !WRITE(*,*)'DEBUG,aSide%tmp,bSide%tmp:',aSide%tmp,bSide%tmp
        JoinFaces(1,iPeriodic)=aSide%Elem%ind-1  !treeID of face with smaller p4est locfaceID
        JoinFaces(2,iPeriodic)=bSide%Elem%ind-1  ! neighbor tree id
        JoinFaces(3,iPeriodic)=aSide%tmp         ! p4est locSideID
        JoinFaces(4,iPeriodic)=bSide%tmp         ! p4est neighbor locSideID
        DO PFlip=0,3
          Hflip_test=getHflip(aSide%tmp,bSide%tmp,PFlip)
          IF(HFlip_test.EQ.HFlip) EXIT
        END DO
        JoinFaces(5,iPeriodic)=PFlip
        !  WRITE(*,*)'DEBUG,JoinFaces,iPeriodic',iPeriodic,JoinFaces(:,iPeriodic)
      END IF
    END DO !iLocSide
  END DO !iElem
END IF !num_periodics>0

CALL p4_connectivity_treevertex(num_vertices,num_trees,vertices,tree_to_vertex, &
                                   num_periodics,JoinFaces,connectivity)

DEALLOCATE(Vertices,tree_to_vertex)
IF(num_periodics.GT.0) DEALLOCATE(JoinFaces) 
 

! COUNT SIDES

nBCSides=0
nSides=0
nPeriodicSides=0
DO iElem=1,nElems
  aElem=>Elems(iElem)%ep
  DO iLocSide=1,6
    aSide=>aElem%Side(iLocSide)%sp
    aSide%tmp=0 
  END DO !iLocSide
END DO !iElem
DO iElem=1,nElems
  aElem=>Elems(iElem)%ep
  DO iLocSide=1,6
    aSide=>aElem%Side(iLocSide)%sp

    IF(aSide%tmp.EQ.0)THEN
      nSides=nSides+1
      aSide%tmp=-1 !used as marker
      IF(ASSOCIATED(aSide%connection)) aSide%connection%tmp=-1
      IF(aSide%BCindex.NE.0)THEN !side is BC or periodic side
        IF(ASSOCIATED(aSide%connection))THEN
          nPeriodicSides=nPeriodicSides+1
        ELSE
          nBCSides=nBCSides+1
        END IF
      END IF
    END IF
  END DO !iLocSide
END DO !iElem


WRITE(*,*)'-------------------------------------------------------'
WRITE(*,'(A22,I8)' )'NGeo:',NGeo
WRITE(*,'(A22,X7L)')'useCurveds:',useCurveds
WRITE(*,'(A22,I8)' )'nElems:',nElems
WRITE(*,'(A22,I8)' )'nNodes:',nNodes
WRITE(*,'(A22,I8)' )'nSides:',nSides
WRITE(*,'(A22,I8)' )'nBCSides:',nBCSides
WRITE(*,'(A22,I8)' )'nPeriodicSides:',nPeriodicSides
WRITE(*,*)'-------------------------------------------------------'

END SUBROUTINE ReadMeshFromHDF5nobuildp4est

SUBROUTINE Buildp4est()
!===================================================================================================================================
! Subroutine to build the p4est and the p4est_geometry from connectivity
!===================================================================================================================================
! MODULES
USE MOD_P4EST_Vars,    ONLY: connectivity,p4est,geom
USE MOD_P4EST_Binding, ONLY: p4_build_p4est
!-----------------------------------------------------------------------------------------------------------------------------------
CALL p4_build_p4est(connectivity,p4est,geom)
END SUBROUTINE Buildp4est

SUBROUTINE ReadMeshFromHDF5(FileString)
!===================================================================================================================================
! Subroutine to read the mesh from a mesh data file, build p4est and p4est_connectivity
!===================================================================================================================================
! INPUT VARIABLES
CHARACTER(LEN=*),INTENT(IN)  :: FileString
!-----------------------------------------------------------------------------------------------------------------------------------
CALL ReadMeshFromHDF5nobuildp4est(FileString)
CALL Buildp4est
END SUBROUTINE ReadMeshFromHDF5

SUBROUTINE ReadGeoFromHDF5(FileString)
!===================================================================================================================================
! Subroutine to read the mesh from a mesh data file
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Mesh_Vars, ONLY: nElems,XGeo,Ngeo,nNodes,HexMap
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
CHARACTER(LEN=*),INTENT(IN)  :: FileString
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                        :: i,j,k
INTEGER                        :: iElem
INTEGER                        :: nNodeIDs
INTEGER                        :: FirstNodeInd,LastNodeInd
LOGICAL                        :: fileExists
INTEGER,ALLOCATABLE            :: ElemInfo(:,:),NodeInfo(:)
REAL,ALLOCATABLE               :: NodeCoords(:,:)
!-----------------------------------------------------------------------------------------------------------------------------------
INQUIRE (FILE=TRIM(FileString), EXIST=fileExists)
IF(.NOT.FileExists)  &
    CALL abort(__STAMP__, &
       'readMesh from data file "'//TRIM(FileString)//'" does not exist')

SWRITE(UNIT_stdOut,'(A)')'READ GEOMETRY DATA FROM DATA FILE "'//TRIM(FileString)//'" ...'
SWRITE(UNIT_StdOut,'(132("-"))')
! Open data file
CALL OpenDataFile(FileString,create=.FALSE.,single=.FALSE.)

CALL ReadMeshHeader()

CALL GetDataSize(File_ID,'NodeCoords',nDims,HSize)
nNodes=HSize(1) !global number of unique nodes
DEALLOCATE(HSize)

!----------------------------------------------------------------------------------------------------------------------------
!                              ELEMENTS
!----------------------------------------------------------------------------------------------------------------------------

!read local ElemInfo from data file
ALLOCATE(ElemInfo(1:nElems,ELEM_InfoSize))
CALL ReadArray('ElemInfo',2,(/nElems,ELEM_InfoSize/),0,1,IntegerArray=ElemInfo)

!----------------------------------------------------------------------------------------------------------------------------
!                              NODES
!----------------------------------------------------------------------------------------------------------------------------

!read local Node Info from data file 
nNodeIDs=ElemInfo(nElems,ELEM_LastNodeInd)-ElemInfo(1,ELEM_FirstNodeInd)
ALLOCATE(NodeInfo(1:nNodeIDs))
CALL ReadArray('NodeInfo',1,(/nNodeIDs/),0,1,IntegerArray=NodeInfo)

! get physical coordinates
ALLOCATE(NodeCoords(nNodes,3))

CALL ReadArray('NodeCoords',2,(/nNodes,3/),0,1,RealArray=NodeCoords)

CALL CloseDataFile() 

ALLOCATE(Xgeo(1:3,0:Ngeo,0:Ngeo,0:Ngeo,nElems))
IF(Ngeo.EQ.1)THEN !use the corner nodes
  DO iElem=1,nElems
    firstNodeInd=ElemInfo(iElem,ELEM_FirstNodeInd) !first index -1 in NodeInfo
    Xgeo(:,0,0,0,iElem)=NodeCoords(NodeInfo(firstNodeInd+1),:)
    Xgeo(:,1,0,0,iElem)=NodeCoords(NodeInfo(firstNodeInd+2),:)
    Xgeo(:,1,1,0,iElem)=NodeCoords(NodeInfo(firstNodeInd+3),:)
    Xgeo(:,0,1,0,iElem)=NodeCoords(NodeInfo(firstNodeInd+4),:)
    Xgeo(:,0,0,1,iElem)=NodeCoords(NodeInfo(firstNodeInd+5),:)
    Xgeo(:,1,0,1,iElem)=NodeCoords(NodeInfo(firstNodeInd+6),:)
    Xgeo(:,1,1,1,iElem)=NodeCoords(NodeInfo(firstNodeInd+7),:)
    Xgeo(:,0,1,1,iElem)=NodeCoords(NodeInfo(firstNodeInd+8),:)
 END DO !iElem=1,nElems
ELSE
  DO iElem=1,nElems
    firstNodeInd=ElemInfo(iElem,ELEM_FirstNodeInd) !first index -1 in NodeInfo
    lastNodeInd=ElemInfo(iElem,ELEM_LastNodeInd)
    IF(LastNodeInd-firstNodeInd-14.NE.(Ngeo+1)**3) STOP 'Problem with curved'
    firstNodeInd=firstNodeInd +8
    DO k=0,Ngeo; DO j=0,Ngeo; DO i=0,Ngeo
      Xgeo(:,i,j,k,iElem)=NodeCoords(NodeInfo(firstNodeInd+HexMap(i,j,k)),:)
    END DO ; END DO ; END DO 
 END DO !iElem=1,nElems
END IF

DEALLOCATE(ElemInfo,NodeInfo,NodeCoords)


END SUBROUTINE ReadGeoFromHDF5


END MODULE MOD_Mesh_ReadIn
