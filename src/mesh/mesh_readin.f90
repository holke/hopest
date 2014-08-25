#include "hopest_f.h"

MODULE MODH_Mesh_ReadIn
!===================================================================================================================================
! Add comments please!
!===================================================================================================================================
! MODULES
USE MODH_HDF5_Input
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------

! Public Part ----------------------------------------------------------------------------------------------------------------------
INTERFACE ReadMeshFromHDF5
  MODULE PROCEDURE ReadMeshFromHDF5
END INTERFACE

INTERFACE ReadGeoFromHDF5
  MODULE PROCEDURE ReadGeoFromHDF5
END INTERFACE

INTERFACE ReadMeshHeader
  MODULE PROCEDURE ReadMeshHeader
END INTERFACE

PUBLIC::ReadMeshFromHDF5
PUBLIC::ReadGeoFromHDF5
PUBLIC::ReadMeshHeader
!===================================================================================================================================

CONTAINS

SUBROUTINE ReadBCs()
!===================================================================================================================================
! Read boundary conditions from data file
!===================================================================================================================================
! MODULES
USE MODH_Globals
USE MODH_Mesh_Vars,ONLY:BoundaryName,BoundaryType,nBCs
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                        :: iBC,Offset
!===================================================================================================================================
! Get number of boundary condtions
CALL GetDataSize(File_ID,'BCNames',nDims,HSize)
nBCs=HSize(1)
DEALLOCATE(HSize)
CALL GetDataSize(File_ID,'BCType',nDims,HSize)
IF(HSize(1).NE.nBCs) STOP 'Problem in readBC'
DEALLOCATE(HSize)

ALLOCATE(BoundaryName(nBCs))
ALLOCATE(BoundaryType(nBCs,BC_SIZE))
offset=0
CALL ReadArray('BCNames',1,(/nBCs/),  Offset,1,StrArray    =BoundaryName)
CALL ReadArray('BCType' ,2,(/nBCs,4/),Offset,1,IntegerArray=BoundaryType)

SWRITE(UNIT_StdOut,'(132("."))')
SWRITE(Unit_StdOut,'(A,A16,A20,A10,A10,A10,A10)')'BOUNDARY CONDITIONS','|','Name','Type','Curved','State','Alpha'
DO iBC=1,nBCs
  SWRITE(*,'(A,A33,A20,I10,I10,I10,I10)')' |','|',TRIM(BoundaryName(iBC)),BoundaryType(iBC,:)
END DO
SWRITE(UNIT_StdOut,'(132("."))')
END SUBROUTINE ReadBCs


SUBROUTINE SetUserBCs()
!===================================================================================================================================
! The user can redefine boundaries in the ini file. We create the mappings for the boundaries.
!===================================================================================================================================
! MODULES
USE MODH_Globals
USE MODH_Mesh_Vars,  ONLY: BoundaryName,BoundaryType,nBCs,nUserBCs
USE MODH_ReadinTools,ONLY: CNTSTR,GETSTR,GETINTARRAY
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                        :: BCMapping(nBCs)
CHARACTER(LEN=255),ALLOCATABLE :: BoundaryNameUser(:)
INTEGER,ALLOCATABLE            :: BoundaryTypeUser(:,:)
INTEGER                        :: iBC,iUserBC
!===================================================================================================================================
! read in boundary conditions, will overwrite BCs from meshfile!
nUserBCs = CNTSTR('BoundaryName',0)
IF(nUserBCs.EQ.0) RETURN

! Read user BC
ALLOCATE(BoundaryNameUser(nUserBCs))
ALLOCATE(BoundaryTypeUser(nUserBCs,2))
DO iBC=1,nUserBCs
  BoundaryNameUser(iBC)   = GETSTR('BoundaryName')
  BoundaryTypeUser(iBC,:) = GETINTARRAY('BoundaryType',2) !(/Type,State/)
END DO

! Override BCs
BCMapping=0
DO iBC=1,nBCs
  DO iUserBC=1,nUserBCs
    IF(INDEX(TRIM(BoundaryNameUser(iUserBC)),TRIM(BoundaryName(iBC))) .NE.0)THEN
      SWRITE(Unit_StdOut,'(A,A)')    ' |     Boundary in HDF file found | ',TRIM(BoundaryName(iBC))
      SWRITE(Unit_StdOut,'(A,I2,I2)')' |                            was | ',BoundaryType(iBC,1),BoundaryType(iBC,3)
      SWRITE(Unit_StdOut,'(A,I2,I2)')' |                      is set to | ',BoundaryTypeUser(iUserBC,1:2)
      BoundaryType(iBC,1) = BoundaryTypeUser(iUserBC,1)
      BoundaryType(iBC,3) = BoundaryTypeUser(iUserBC,2)
    END IF
  END DO
END DO

SWRITE(UNIT_StdOut,'(132("."))')
DEALLOCATE(BoundaryNameUser,BoundaryTypeUser)
END SUBROUTINE SetUserBCs


SUBROUTINE ReadMeshHeader(FileString)
!===================================================================================================================================
! Subroutine to read the mesh from a mesh data file
!===================================================================================================================================
! MODULES
USE MODH_Globals
USE MODH_Mesh_Vars,ONLY: NGeo,useCurveds,nGlobalTrees
USE MODH_Mesh,     ONLY: SetCurvedInfo
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
CHARACTER(LEN=*),INTENT(IN)  :: FileString
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER :: BoundaryOrder_mesh
!===================================================================================================================================
CALL OpenDataFile(FileString,create=.FALSE.,single=.FALSE.)
CALL GetDataSize(File_ID,'ElemInfo',nDims,HSize)
nGlobalTrees=HSize(1) !global number of elements
DEALLOCATE(HSize)

CALL ReadAttribute(File_ID,'BoundaryOrder',1,IntegerScalar=BoundaryOrder_mesh)
NGeo = BoundaryOrder_mesh-1
CALL SetCurvedInfo()

CALL ReadAttribute(File_ID,'CurvedFound',1,LogicalScalar=useCurveds)

CALL readBCs()
IF(hopestMode.EQ.2) CALL setUserBCs()
CALL CloseDataFile() 

END SUBROUTINE ReadMeshHeader


SUBROUTINE ReadMeshFromHDF5(FileString)
!===================================================================================================================================
! Subroutine to read the mesh from a mesh data file and build p4est_connectivity
!===================================================================================================================================
! MODULES
USE MODH_Globals
USE MODH_Mesh_Vars
USE MODH_P4EST_Vars,    ONLY: connectivity,H2P_VertexMap,H2P_FaceMap
USE MODH_P4EST_Binding, ONLY: p4_connectivity_treevertex,p4_build_p4est
USE MODH_P4EST,         ONLY: getHFlip
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
INTEGER                        :: iTree,ElemID
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
INQUIRE (FILE=TRIM(FileString), EXIST=fileExists)
IF(.NOT.FileExists)  &
    CALL abort(__STAMP__, &
       'readMesh from data file "'//TRIM(FileString)//'" does not exist')

SWRITE(UNIT_stdOut,'(A)')'READ MESH FROM DATA FILE "'//TRIM(FileString)//'" ...'
SWRITE(UNIT_StdOut,'(132("-"))')
! Open data file
CALL OpenDataFile(FileString,create=.FALSE.,single=.FALSE.)

!----------------------------------------------------------------------------------------------------------------------------
!                              ELEMENTS
!----------------------------------------------------------------------------------------------------------------------------

nTrees=nGlobalTrees   !local number of Elements 
!read local ElemInfo from data file
ALLOCATE(ElemInfo(1:nTrees,ELEM_InfoSize))
CALL ReadArray('ElemInfo',2,(/nTrees,ELEM_InfoSize/),0,1,IntegerArray=ElemInfo)

ALLOCATE(Trees(1:nTrees))

DO iTree=1,nTrees
  iSide=ElemInfo(iTree,ELEM_FirstSideInd) !first index -1 in Sideinfo
  iNode=ElemInfo(iTree,ELEM_FirstNodeInd) !first index -1 in NodeInfo
  Trees(iTree)%ep=>GETNEWELEM()
  aElem=>Trees(iTree)%ep
  aElem%Ind    = iTree
  aElem%Type   = ElemInfo(iTree,ELEM_Type)
  aElem%Zone   = ElemInfo(iTree,ELEM_Zone)
END DO

!----------------------------------------------------------------------------------------------------------------------------
!                              NODES
!----------------------------------------------------------------------------------------------------------------------------

!read local Node Info from data file 
nNodeIDs=ElemInfo(nTrees,ELEM_LastNodeInd)-ElemInfo(1,ELEM_FirstNodeInd)
ALLOCATE(NodeInfo(1:nNodeIDs))
CALL ReadArray('NodeInfo',1,(/nNodeIDs/),0,1,IntegerArray=NodeInfo)

ALLOCATE(ElemCurvedNode(nCurvedNodes,nTrees))

CALL GetDataSize(File_ID,'NodeCoords',nDims,HSize)
nNodes=HSize(1) !global number of unique nodes
DEALLOCATE(HSize)

ALLOCATE(Nodes(1:nNodes)) ! pointer list, entry is known by NodeCoords
DO iNode=1,nNodes
  NULLIFY(Nodes(iNode)%np)
END DO
!assign nodes 
DO iTree=1,nTrees
  aElem=>Trees(iTree)%ep
  iNode=ElemInfo(iTree,ELEM_FirstNodeInd) !first index -1 in NodeInfo
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
    nCurvedNodes_loc = ElemInfo(iTree,ELEM_LastNodeInd) - ElemInfo(iTree,ELEM_FirstNodeInd) - 14 ! corner + oriented nodes
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
      ElemCurvedNode(i,iTree)%np=>Nodes(NodeID)%np
    END DO
  END IF
END DO

!----------------------------------------------------------------------------------------------------------------------------
!                              SIDES
!----------------------------------------------------------------------------------------------------------------------------

nSideIDs=ElemInfo(nTrees,ELEM_LastSideInd)-ElemInfo(1,ELEM_FirstSideInd)
!read local SideInfo from data file 
ALLOCATE(SideInfo(1:nSideIDs,SIDE_InfoSize))
CALL ReadArray('SideInfo',2,(/nSideIDs,SIDE_InfoSize/),0,1,IntegerArray=SideInfo)

DO iTree=1,nTrees
  aElem=>Trees(iTree)%ep
  iNode=ElemInfo(iTree,ELEM_LastNodeInd) !first index -1 in NodeInfo
  iNode=iNode-6
  iSide=ElemInfo(iTree,ELEM_FirstSideInd) !first index -1 in Sideinfo
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
END DO !iTree

 
! build up side connection 
DO iTree=1,nTrees
  aElem=>Trees(iTree)%ep
  iSide=ElemInfo(iTree,ELEM_FirstSideInd) !first index -1 in Sideinfo
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
        aSide%flip  = 0
        elemID      = 0
      END IF
    ELSE
      aSide%BCindex = 0
    END IF

    IF(.NOT.ASSOCIATED(aSide%connection))THEN
      IF((elemID.NE.0).AND.doConnection)THEN !connection 
        IF((elemID.LE.nTrees).AND.(elemID.GE.1))THEN !local connection
          DO jLocSide=1,6
            bSide=>Trees(elemID)%ep%Side(jLocSide)%sp
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
END DO !iTree

DEALLOCATE(ElemInfo,SideInfo,NodeInfo)

! get physical coordinates

ALLOCATE(NodeCoords(nNodes,3))

CALL ReadArray('NodeCoords',2,(/nNodes,3/),0,1,RealArray=NodeCoords)

ALLOCATE(Xgeo(1:3,0:Ngeo,0:Ngeo,0:Ngeo,nTrees))
IF(Ngeo.EQ.1)THEN !use the corner nodes
  DO iTree=1,nTrees
    aElem=>Trees(iTree)%ep
    Xgeo(:,0,0,0,iTree)=NodeCoords(aElem%Node(1)%np%ind,:)
    Xgeo(:,1,0,0,iTree)=NodeCoords(aElem%Node(2)%np%ind,:)
    Xgeo(:,1,1,0,iTree)=NodeCoords(aElem%Node(3)%np%ind,:)
    Xgeo(:,0,1,0,iTree)=NodeCoords(aElem%Node(4)%np%ind,:)
    Xgeo(:,0,0,1,iTree)=NodeCoords(aElem%Node(5)%np%ind,:)
    Xgeo(:,1,0,1,iTree)=NodeCoords(aElem%Node(6)%np%ind,:)
    Xgeo(:,1,1,1,iTree)=NodeCoords(aElem%Node(7)%np%ind,:)
    Xgeo(:,0,1,1,iTree)=NodeCoords(aElem%Node(8)%np%ind,:)
 END DO !iTree=1,nTrees
ELSE
  DO iTree=1,nTrees
    aElem=>Trees(iTree)%ep
    l=0
    DO k=0,Ngeo; DO j=0,Ngeo; DO i=0,Ngeo
      l=l+1
      Xgeo(:,i,j,k,iTree)=NodeCoords(ElemCurvedNode(l,iTree)%np%ind,:)
    END DO ; END DO ; END DO 
 END DO !iTree=1,nTrees
END IF

CALL CloseDataFile() 

DEALLOCATE(ElemCurvedNode)


! P4est MESH connectivity (should be replaced by connectivity information ?)

! needs unique corner nodes for mesh connectivity
DO iNode=1,nNodes
  Nodes(iNode)%np%tmp=-1
END DO
num_vertices=0
DO iTree=1,nTrees
  aElem=>Trees(iTree)%ep
  DO iNode=1,8
    aNode=>aElem%Node(iNode)%np
    IF(aNode%tmp.EQ.-1)THEN
      num_vertices=num_vertices+1
      aElem%Node(iNode)%np%tmp=num_vertices
    END IF
  END DO
END DO !iTree

ALLOCATE(Vertices(3,num_vertices))
DO iNode=1,nNodes
  aNode=>Nodes(iNode)%np
  IF(aNode%tmp.GT.0)THEN
    Vertices(:,aNode%tmp)=NodeCoords(aNode%ind,:)
  END IF
END DO


DEALLOCATE(NodeCoords)

num_trees=nTrees
ALLOCATE(tree_to_vertex(8,num_trees))
DO iTree=1,nTrees
  aElem=>Trees(iTree)%ep
  DO iNode=1,8
    tree_to_vertex(iNode,iTree)=aElem%Node(H2P_VertexMap(iNode)+1)%np%tmp-1
  END DO
END DO

!periodic Boundaries
num_periodics=0
DO iTree=1,nTrees
  aElem=>Trees(iTree)%ep
  DO iLocSide=1,6
    aSide=>aElem%Side(iLocSide)%sp
    IF(aSide%BCIndex.EQ.0) CYCLE ! NO Boundary Condition
    IF((BoundaryType(aSide%BCIndex,BC_TYPE).EQ.1).AND.(aSide%flip.EQ.0))THEN !periodic side 
      num_periodics=num_periodics+1
    END IF
  END DO !iLocSide
END DO !iTree


IF(num_periodics.GT.0) THEN
  ALLOCATE(JoinFaces(5,num_periodics))
  DO iTree=1,nTrees
    aElem=>Trees(iTree)%ep
    DO iLocSide=1,6
      aElem%Side(iLocSide)%sp%tmp=H2P_FaceMap(iLocSide)  !local Face ID in p4est
    END DO
  END DO
  
  iperiodic=0
  DO iTree=1,nTrees
    aElem=>Trees(iTree)%ep
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
        !WRITE(*,*)'DEBUG,JoinFaces,iPeriodic',iPeriodic,JoinFaces(:,iPeriodic)
      END IF
    END DO !iLocSide
  END DO !iTree
END IF !num_periodics>0

CALL p4_connectivity_treevertex(num_vertices,num_trees,vertices,tree_to_vertex, &
                                   num_periodics,JoinFaces,connectivity)

DEALLOCATE(Vertices,tree_to_vertex)
IF(num_periodics.GT.0) DEALLOCATE(JoinFaces) 
 

! COUNT SIDES

nBCSides=0
nSides=0
nPeriodicSides=0
DO iTree=1,nTrees
  aElem=>Trees(iTree)%ep
  DO iLocSide=1,6
    aSide=>aElem%Side(iLocSide)%sp
    aSide%tmp=0 
  END DO !iLocSide
END DO !iTree
DO iTree=1,nTrees
  aElem=>Trees(iTree)%ep
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
END DO !iTree


WRITE(*,*)'-------------------------------------------------------'
WRITE(*,'(A22,I8)' )'NGeo:',NGeo
WRITE(*,'(A22,X7L)')'useCurveds:',useCurveds
WRITE(*,'(A22,I8)' )'nTrees:',nTrees
WRITE(*,'(A22,I8)' )'nNodes:',nNodes
WRITE(*,'(A22,I8)' )'nSides:',nSides
WRITE(*,'(A22,I8)' )'nBCSides:',nBCSides
WRITE(*,'(A22,I8)' )'nPeriodicSides:',nPeriodicSides
WRITE(*,*)'-------------------------------------------------------'

END SUBROUTINE ReadMeshFromHDF5



SUBROUTINE ReadGeoFromHDF5(FileString)
!===================================================================================================================================
! Subroutine to read the mesh from a mesh data file
!===================================================================================================================================
! MODULES
USE MODH_Globals
USE MODH_Mesh_Vars, ONLY: nTrees,XGeo,Ngeo,nNodes,HexMap
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
INTEGER                        :: iTree
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

CALL GetDataSize(File_ID,'NodeCoords',nDims,HSize)
nNodes=HSize(1) !global number of unique nodes
DEALLOCATE(HSize)

!----------------------------------------------------------------------------------------------------------------------------
!                              ELEMENTS
!----------------------------------------------------------------------------------------------------------------------------

!read local ElemInfo from data file
ALLOCATE(ElemInfo(1:nTrees,ELEM_InfoSize))
CALL ReadArray('ElemInfo',2,(/nTrees,ELEM_InfoSize/),0,1,IntegerArray=ElemInfo)

!----------------------------------------------------------------------------------------------------------------------------
!                              NODES
!----------------------------------------------------------------------------------------------------------------------------

!read local Node Info from data file 
nNodeIDs=ElemInfo(nTrees,ELEM_LastNodeInd)-ElemInfo(1,ELEM_FirstNodeInd)
ALLOCATE(NodeInfo(1:nNodeIDs))
CALL ReadArray('NodeInfo',1,(/nNodeIDs/),0,1,IntegerArray=NodeInfo)

! get physical coordinates
ALLOCATE(NodeCoords(nNodes,3))

CALL ReadArray('NodeCoords',2,(/nNodes,3/),0,1,RealArray=NodeCoords)

CALL CloseDataFile() 

ALLOCATE(Xgeo(1:3,0:Ngeo,0:Ngeo,0:Ngeo,nTrees))
IF(Ngeo.EQ.1)THEN !use the corner nodes
  DO iTree=1,nTrees
    firstNodeInd=ElemInfo(iTree,ELEM_FirstNodeInd) !first index -1 in NodeInfo
    Xgeo(:,0,0,0,iTree)=NodeCoords(NodeInfo(firstNodeInd+1),:)
    Xgeo(:,1,0,0,iTree)=NodeCoords(NodeInfo(firstNodeInd+2),:)
    Xgeo(:,1,1,0,iTree)=NodeCoords(NodeInfo(firstNodeInd+3),:)
    Xgeo(:,0,1,0,iTree)=NodeCoords(NodeInfo(firstNodeInd+4),:)
    Xgeo(:,0,0,1,iTree)=NodeCoords(NodeInfo(firstNodeInd+5),:)
    Xgeo(:,1,0,1,iTree)=NodeCoords(NodeInfo(firstNodeInd+6),:)
    Xgeo(:,1,1,1,iTree)=NodeCoords(NodeInfo(firstNodeInd+7),:)
    Xgeo(:,0,1,1,iTree)=NodeCoords(NodeInfo(firstNodeInd+8),:)
 END DO !iTree=1,nTrees
ELSE
  DO iTree=1,nTrees
    firstNodeInd=ElemInfo(iTree,ELEM_FirstNodeInd) !first index -1 in NodeInfo
    lastNodeInd=ElemInfo(iTree,ELEM_LastNodeInd)
    IF(LastNodeInd-firstNodeInd-14.NE.(Ngeo+1)**3) STOP 'Problem with curved'
    firstNodeInd=firstNodeInd +8
    DO k=0,Ngeo; DO j=0,Ngeo; DO i=0,Ngeo
      Xgeo(:,i,j,k,iTree)=NodeCoords(NodeInfo(firstNodeInd+HexMap(i,j,k)),:)
    END DO ; END DO ; END DO 
 END DO !iTree=1,nTrees
END IF

DEALLOCATE(ElemInfo,NodeInfo,NodeCoords)


END SUBROUTINE ReadGeoFromHDF5


END MODULE MODH_Mesh_ReadIn
