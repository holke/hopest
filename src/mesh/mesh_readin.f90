#include "../hopest_f.h"

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
INTEGER,PARAMETER    :: ElemInfoSize=6        !number of entry in each line of ElemInfo
INTEGER,PARAMETER    :: ELEM_Type=1           !entry position, 
INTEGER,PARAMETER    :: ELEM_Zone=2           
INTEGER,PARAMETER    :: ELEM_FirstSideInd=3
INTEGER,PARAMETER    :: ELEM_LastSideInd=4
INTEGER,PARAMETER    :: ELEM_FirstNodeInd=5
INTEGER,PARAMETER    :: ELEM_LastNodeInd=6

INTEGER,PARAMETER    :: SideInfoSize=4
INTEGER,PARAMETER    :: SIDE_Type=1           !entry position
INTEGER,PARAMETER    :: SIDE_ID=2
INTEGER,PARAMETER    :: SIDE_nbElemID=3
INTEGER,PARAMETER    :: SIDE_BCID=4

INTEGER,ALLOCATABLE  :: ElemInfo(:,:),SideInfo(:,:),NodeInfo(:)
REAL,ALLOCATABLE     :: NodeCoords(:,:)

INTEGER              :: BoundaryOrder_mesh
INTEGER              :: nNodeIDs,nSideIDs
! Public Part ----------------------------------------------------------------------------------------------------------------------
INTERFACE ReadMesh
  MODULE PROCEDURE ReadMesh
END INTERFACE

INTERFACE BuildEdges
  MODULE PROCEDURE BuildEdges
END INTERFACE

PUBLIC::ReadMesh,BuildEdges
!===================================================================================================================================

CONTAINS

SUBROUTINE ReadBCs()
!===================================================================================================================================
! Read boundary conditions from data file
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Mesh_Vars,ONLY:BoundaryName,BoundaryType,nBCs
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
INTEGER                        :: iBC
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

! Read boundary types from data file
CALL GetDataSize(File_ID,'BCType',nDims,HSize)
IF(HSize(1).NE.nBCs) STOP 'Problem in readBC'
DEALLOCATE(HSize)
ALLOCATE(BCType(nBCs,4))
offset=0
CALL ReadArray('BCType',2,(/nBCs,4/),Offset,1,IntegerArray=BCType)

IF(ALLOCATED(BoundaryName)) DEALLOCATE(BoundaryName)
IF(ALLOCATED(BoundaryType)) DEALLOCATE(BoundaryType)
ALLOCATE(BoundaryName(nBCs))
ALLOCATE(BoundaryType(nBCs,3))
BoundaryName = BCNames
BoundaryType(:,BC_TYPE)  = BCType(:,1)
BoundaryType(:,BC_STATE) = BCType(:,3)
BoundaryType(:,BC_ALPHA) = BCType(:,4)
SWRITE(UNIT_StdOut,'(132("."))')
SWRITE(Unit_StdOut,'(A,A16,A20,A10,A10,A10)')'BOUNDARY CONDITIONS','|','Name','Type','State','Alpha'
DO iBC=1,nBCs
  SWRITE(*,'(A,A33,A20,I10,I10,I10)')' |','|',TRIM(BoundaryName(iBC)),BoundaryType(iBC,:)
END DO
SWRITE(UNIT_StdOut,'(132("."))')
DEALLOCATE(BCNames,BCType,BCMapping)
END SUBROUTINE ReadBCs


SUBROUTINE ReadMesh(FileString)
!===================================================================================================================================
! Subroutine to read the mesh from a mesh data file
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Mesh_Vars
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
INTEGER                        :: BCindex
INTEGER                        :: iElem,ElemID
INTEGER                        :: iNode,jNode,NodeID,SideID
INTEGER                        :: iLocSide,jLocSide
INTEGER                        :: iSide
INTEGER                        :: FirstNodeInd,LastNodeInd,FirstSideInd,LastSideInd
INTEGER                        :: nCurvedNodes,nCurvedNodesRef
LOGICAL                        :: oriented
INTEGER                        :: nPeriodicSides 
LOGICAL                        :: fileExists
LOGICAL                        :: doConnection
TYPE(tElem),POINTER            :: aElem
TYPE(tSide),POINTER            :: aSide,bSide
!===================================================================================================================================
IF(MESHInitIsDone) RETURN
IF(MPIRoot)THEN
  INQUIRE (FILE=TRIM(FileString), EXIST=fileExists)
  IF(.NOT.FileExists)  CALL abort(__STAMP__, &
          'readMesh from data file "'//TRIM(FileString)//'" does not exist')
END IF

SWRITE(UNIT_stdOut,'(A)')'READ MESH FROM DATA FILE "'//TRIM(FileString)//'" ...'
SWRITE(UNIT_StdOut,'(132("-"))')
! Open data file
CALL OpenDataFile(FileString,create=.FALSE.,single=.FALSE.)

CALL GetDataSize(File_ID,'ElemInfo',nDims,HSize)
nGlobalElems=HSize(1) !global number of elements
DEALLOCATE(HSize)
nElems=nGlobalElems   !local number of Elements 

CALL GetDataSize(File_ID,'NodeCoords',nDims,HSize)
nNodes=HSize(1) !global number of unique nodes
DEALLOCATE(HSize)

CALL readBCs()

CALL ReadAttribute(File_ID,'BoundaryOrder',1,IntegerScalar=BoundaryOrder_mesh)
NGeo = BoundaryOrder_mesh-1
CALL ReadAttribute(File_ID,'CurvedFound',1,LogicalScalar=useCurveds)
!----------------------------------------------------------------------------------------------------------------------------
!                              ELEMENTS
!----------------------------------------------------------------------------------------------------------------------------

!read local ElemInfo from data file
ALLOCATE(ElemInfo(1:nElems,ElemInfoSize))
CALL ReadArray('ElemInfo',2,(/nElems,ElemInfoSize/),0,1,IntegerArray=ElemInfo)

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


nCurvedNodesRef=(NGeo+1)**3
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
    nCurvedNodes = ElemInfo(iElem,ELEM_LastNodeInd) - ElemInfo(iElem,ELEM_FirstNodeInd) - 14 ! corner + oriented nodes
    IF(nCurvedNodes.NE.nCurvedNodesRef) &
      CALL abort(__STAMP__,'Wrong number of curved nodes.')
    ALLOCATE(aElem%CurvedNode(0:NGeo,0:NGeo,0:NGeo))
    DO k=0,NGeo; DO j=0,NGeo; DO i=0,NGeo
      iNode=iNode+1
      NodeID=NodeInfo(iNode) !first oriented corner node
      IF(.NOT.ASSOCIATED(Nodes(NodeID)%np))THEN
        ALLOCATE(Nodes(NodeID)%np)
        Nodes(NodeID)%np%ind=NodeID 
      END IF
      aElem%CurvedNode(i,j,k)%np=>Nodes(NodeID)%np
    END DO; END DO; END DO
  END IF
END DO

!----------------------------------------------------------------------------------------------------------------------------
!                              SIDES
!----------------------------------------------------------------------------------------------------------------------------

nSideIDs=ElemInfo(nElems,ELEM_LastSideInd)-ElemInfo(1,ELEM_FirstSideInd)
!read local SideInfo from data file 
ALLOCATE(SideInfo(1:nSideIDs,SideInfoSize))
CALL ReadArray('SideInfo',2,(/nSideIDs,SideinfoSize/),0,1,IntegerArray=SideInfo)

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
      CALL abort(__STAMP__,'Only conforming meshes in readin.')
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
!assign to pointer
DO iNode=1,nNodes
  IF(ASSOCIATED(Nodes(iNode)%np))Nodes(iNode)%np%x=NodeCoords(iNode,:)
END DO
DEALLOCATE(NodeCoords)
 
CALL CloseDataFile() 

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

CALL buildEdges()


WRITE(*,*)'-------------------------------------------------------'
WRITE(*,'(A22,I8)')'NGeo:',NGeo
WRITE(*,'(A22,I8)')'useCurveds:',useCurveds
WRITE(*,'(A22,I8)')'nElems:',nElems
WRITE(*,'(A22,I8)')'nEdges:',nEdges
WRITE(*,'(A22,I8)')'nNodes:',nNodes
WRITE(*,'(A22,I8)')'nSides:',nSides
WRITE(*,'(A22,I8)')'nBCSides:',nBCSides
WRITE(*,'(A22,I8)')'nPeriodicSides:',nPeriodicSides
WRITE(*,*)'-------------------------------------------------------'

END SUBROUTINE ReadMesh


SUBROUTINE buildEdges()
!===================================================================================================================================
! Create Edge datastructure, each edge is unique, and has a pointer from each side and from the node with the lower index.
! on the node, a list beginning with node%firstEdge is build up. On the Element sides, a edge pointer array Edge(1:nNodes) is 
! filled, together with their orientation inside the side. Very important: OrientedNodes are used!!!! 
! If the edge is oriented, it goes from orientedNode(i)-> orientedNode(i+1), and 
! If the edge is not  oriented, it goes from orientedNode(i+1)-> orientedNode(i)
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Mesh_Vars,ONLY:tElem,tSide,tEdge,tNode
USE MOD_Mesh_Vars,ONLY:nElems,nEdges,Edges,Elems
USE MOD_Mesh_Vars,ONLY:EdgeToElemNode
USE MOD_Mesh_Vars,ONLY:GETNEWEDGE
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
! OUTPUT VARIABLES
! LOCAL VARIABLES
INTEGER                      :: iNode,iEdge,iElem,EdgeInd
TYPE(tElem),POINTER          :: Elem
TYPE(tEdge),POINTER          :: Edge
TYPE(tNode),POINTER          :: aNode,bNode
INTEGER,ALLOCATABLE          :: Edge_nElems(:)
LOGICAL                      :: edgeFound,isNew,isOriented
!===================================================================================================================================
WRITE(UNIT_stdOut,'(132("~"))')
WRITE(UNIT_stdOut,'(A)')'BUILD EDGES ...'

EdgeInd=0
DO iElem=1,nElems
  Elem=>Elems(iElem)%ep
  DO iEdge=1,12
    aNode=>Elem%Node(EdgeToElemNode(1,iEdge))%np
    bNode=>Elem%Node(EdgeToElemNode(2,iEdge))%np
    Elem%Edge(iEdge)%edp => GETNEWEDGE(aNode,bNode,isNew,isOriented)
    IF(isNew)THEN
      EdgeInd=EdgeInd+1
      Elem%Edge(iEdge)%edp%ind    = EdgeInd
      Elem%Edge(iEdge)%isOriented = isOriented
    END IF
    Elem%Edge(iEdge)%edp%nNbElems=Elem%Edge(iEdge)%edp%nNbElems+1
  END DO
END DO !! ELEMS!!
nEdges=EdgeInd

ALLOCATE(Edges(nEdges))
DO iEdge=1,nEdges
  NULLIFY(Edges(iEdge)%edp)
END DO

DO iElem=1,nElems
  Elem=>Elems(iElem)%ep
  DO iEdge=1,12
    Edge=>Elem%Edge(iEdge)%edp
    EdgeInd=Edge%ind
    IF(.NOT.ASSOCIATED(Edges(EdgeInd)%edp))THEN
      Edges(EdgeInd)%edp=>Edge
    END IF
  END DO
END DO
DO iEdge=1,nEdges
  ALLOCATE(Edges(iEdge)%edp%nbElem(Edges(iEdge)%edp%nNbElems))
  Edges(iEdge)%edp%nNbElems=0
END DO

DO iElem=1,nElems
  Elem=>Elems(iElem)%ep
  DO iEdge=1,12
    Edge=>Elem%Edge(iEdge)%edp
    EdgeInd=Edge%ind
    Edge%nNbElems=Edge%nNbElems+1
    Edge%nbElem(Edge%nNbElems)%ep=>Elem
  END DO
END DO

END SUBROUTINE buildEdges


END MODULE MOD_Mesh_ReadIn
