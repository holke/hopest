#include "hopest_f.h"
MODULE MOD_Output_HDF5
!===================================================================================================================================
! Add comments please!
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE HDF5
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES 
!-----------------------------------------------------------------------------------------------------------------------------------
INTEGER(HID_T)                 :: File_ID
INTEGER(HSIZE_T),POINTER       :: HSize(:)
INTEGER,ALLOCATABLE            :: ElemInfo(:,:),SideInfo(:,:),NodeInfo(:)
REAL,ALLOCATABLE               :: NodeCoords(:,:)
REAL,ALLOCATABLE               :: ElemWeight(:)
INTEGER                        :: ElemCounter(11,2)
INTEGER                        :: nSideIDs,nNodeIDs
INTEGER                        :: nTotalSides,nTotalNodes
! Private Part ---------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------

INTERFACE WriteMeshToHDF5
  MODULE PROCEDURE WriteMeshToHDF5
END INTERFACE

PUBLIC::WriteMeshToHDF5
!===================================================================================================================================

CONTAINS

SUBROUTINE WriteMeshToHDF5(FileString)
!===================================================================================================================================
! Subroutine to write Data to HDF5 format
!===================================================================================================================================
! MODULES
USE MOD_Mesh_Vars
USE MOD_IO_HDF5
USE MOD_HDF5_output
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
CHARACTER(LEN=*),INTENT(IN)    :: FileString
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES

TYPE(tElem),POINTER            :: Elem
TYPE(tSide),POINTER            :: Side
INTEGER                        :: iElem,i
INTEGER                        :: NodeID,iNode
INTEGER                        :: SideID,iLocSide
REAL,ALLOCATABLE               :: ElemBary(:,:)
!===================================================================================================================================
WRITE(*,'(132("~"))')
WRITE(*,'(A)')' WRITE MESH TO HDF5 FILE... ' // TRIM(FileString) 
! Create the file collectively.
CALL OpenHDF5File(FileString,create=.TRUE.,single=.TRUE.)  

!set all node and side indices =0
DO iElem=1,nElems
  Elem=>Elems(iElem)%ep
  DO iNode=1,8
    Elem%Node(iNode)%np%ind=0
  END DO
  DO iNode=1,nCurvedNodes
    Elem%curvedNode(iNode)%np%ind=0
  END DO
  DO iLocSide=1,6
    Side=>Elem%Side(iLocSide)%sp
    Side%ind=0
  END DO !iLocSide=1,6
END DO !iElem=1,nElems

! count Elements , unique sides and nodes are marked with ind=0
nNodeIDs=0 !number of unique nodeIDs
nSideIDs=0 !number of unique side IDs (side and side%connection have the same sideID)
nTotalSides=0   !number of all sides
nTotalNodes=0   !number of all nodes


DO iElem=1,nElems
  Elem=>Elems(iElem)%ep
  ! Count nodes
  DO iNode=1,8
    IF(Elem%Node(iNode)%np%ind.NE.0) CYCLE
    nNodeIDs=nNodeIDs+1
    Elem%Node(iNode)%np%ind=-88888  ! mark no MPI side
  END DO

  DO iNode=1,nCurvedNodes
    IF(Elem%CurvedNode(iNode)%np%ind.NE.0) CYCLE
    nNodeIDs=nNodeIDs+1
    Elem%CurvedNode(iNode)%np%ind=-88888
  END DO

  ! Count sides
  DO iLocSide=1,6
    Side=>Elem%Side(iLocSide)%sp
    IF(Side%ind.EQ.0) THEN
      nSideIDs=nSideIDs+1
      Side%ind=-88888
      IF(ASSOCIATED(Side%connection))THEN      
        IF(Side%connection%ind.EQ.0) nSideIDs=nSideIDs-1 ! count inner and periodic sides only once 
      END IF
    END IF
  END DO
  nTotalNodes = nTotalNodes+8+6+nCurvedNodes ! corner + oriented + curved
  nTotalSides = nTotalSides+6
END DO

!set unique nodes and Side Indices
SideID=0
NodeID=0
DO iElem=1,nElems
  Elem=>Elems(iElem)%ep
  Elem%ind=iElem
  DO iNode=1,8
    IF(Elem%Node(iNode)%np%ind.NE.-88888) CYCLE
    NodeID=NodeID+1
    Elem%Node(iNode)%np%ind=NodeID
  END DO
  DO iNode=1,nCurvedNodes
    IF(Elem%CurvedNode(iNode)%np%ind.NE.-88888) CYCLE
    NodeID=NodeID+1
    Elem%CurvedNode(iNode)%np%ind=NodeID
  END DO

  DO iLocSide=1,6
    Side=>Elem%Side(iLocSide)%sp
    IF(side%ind.EQ.-88888) THEN  ! assign side ID only for non MPI sides and lower MPI sides
      SideID=SideID+1
      Side%ind=SideID
      IF(ASSOCIATED(Side%connection))THEN     
        IF(Side%connection%ind.NE.-88888) Side%connection%ind=SideID !already assigned
      END IF
    END IF
  END DO
END DO !Elem

CALL getMeshInfo() !allocates and fills ElemInfo,SideInfo,NodeInfo,NodeCoords

!WRITE ElemInfo,into (1,nElems)  
CALL WriteArrayToHDF5(File_ID,'ElemInfo',nElems,2,(/nElems,ELEM_InfoSize/),0,IntegerArray=ElemInfo)

!WRITE ElemWeight,into (1,nElems)  
CALL WriteArrayToHDF5(File_ID,'ElemWeight',nElems,1,(/nElems/),0,RealArray=ElemWeight)

ALLOCATE(ElemBary(nElems,3))
DO iElem=1,nElems
  ElemBary(iElem,:)=0.
  DO iNode=1,8
    ElemBary(iElem,:)=ElemBary(iElem,:)+Elems(iElem)%ep%Node(iNode)%np%x
  END DO !iNode
  ElemBary(iElem,:)=ElemBary(iElem,:)*0.125  ! / nNodes=8
END DO !iElem=1,nElem
CALL WriteArrayToHDF5(File_ID,'ElemBarycenters',nElems,2,(/nElems,3/),0,RealArray=ElemBary)
DEALLOCATE(ElemBary)


!WRITE SideInfo,into (1,nTotalSides)   
CALL WriteArrayToHDF5(File_ID,'SideInfo',nTotalSides,2,(/nTotalSides,SIDE_InfoSize/),0,IntegerArray=SideInfo)
DEALLOCATE(SideInfo)

!WRITE NodeInfo,into (1,nTotalNodes) 
CALL WriteArrayToHDF5(File_ID,'NodeInfo',nTotalNodes,1,(/nTotalNodes/),0,IntegerArray=NodeInfo)
DEALLOCATE(NodeInfo)

! WRITE NodeCoords (have to be sorted according to nodemap)
CALL WriteArrayToHDF5(File_ID,'NodeCoords',nNodeIDs,2,(/nNodeIDs,3/),0,RealArray=NodeCoords)
DEALLOCATE(NodeCoords)

!! WRITE NodeCoords,arbitrary ordering by NodeMap
!CALL WriteCoordsToHDF5(File_ID,'NodeCoords',nNodeIDs,(/nNodeIDs,3/),NodeMap,NodeCoords)


CALL WriteArrayToHDF5(File_ID,'ElemCounter',11,2,(/11,2/),0,IntegerArray=ElemCounter)
WRITE(*,*)'Mesh statistics:'
WRITE(*,*)'Element Type | number of elements'
DO i=1,11
  WRITE(*,'(I4,A,I8)') Elemcounter(i,1),'        | ',Elemcounter(i,2)
END DO

!attributes 
CALL WriteAttributeToHDF5(File_ID,'BoundaryOrder',1,IntegerScalar=Ngeo+1)
CALL WriteAttributeToHDF5(File_ID,'CurvedFound',1,LogicalScalar=useCurveds)
! WRITE BC 
CALL WriteArrayToHDF5(File_ID,'BCNames',nBCs,1,(/nBCs/),0,StrArray=BoundaryName)
CALL WriteArrayToHDF5(File_ID,'BCType',nBCs,2,(/nBcs,4/),0,IntegerArray=BoundaryType)

! Close the file.
CALL CloseHDF5File()

DEALLOCATE(ElemInfo)
DEALLOCATE(ElemWeight)

END SUBROUTINE WriteMeshToHDF5



SUBROUTINE getMeshInfo()
!===================================================================================================================================
! Subroutine prepares ElemInfo,Sideinfo,Nodeinfo,NodeCoords arrays
!===================================================================================================================================
! MODULES
USE MOD_Mesh_Vars
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
TYPE(tElem),POINTER            :: Elem
TYPE(tSide),POINTER            :: Side
INTEGER                        :: iNodeP,iNode,NodeID
INTEGER                        :: iElem
INTEGER                        :: iSide,iLocSide 
INTEGER                        :: locnSides,locnNodes
!===================================================================================================================================
!fill ElementInfo. 
ALLOCATE(ElemInfo(1:nElems,ELEM_InfoSize))
ALLOCATE(ElemWeight(1:nElems))
ElemInfo=0
ElemWeight=0.
Elemcounter=0
Elemcounter(:,1)=(/104,204,105,115,205,106,116,206,108,118,208/)
NodeID  = 0 
iSide  = 0
DO iElem=1,nElems
  Elem=>Elems(iElem)%ep
  locnNodes=8+nCurvedNodes
  locnSides=0
  ! for element sides
  DO iLocSide=1,6
    Side=>Elem%Side(iLocSide)%sp
    locnNodes = locnNodes+1 ! wirte only first oriented node of side, if curved all!         
    locnSides = locnSides+1
  END DO !iLocSide=1,6
  SELECT CASE(elem%type)
    CASE(104) !linear tet
      elemcounter(1,2)=elemcounter(1,2)+1
    CASE(204) !spline tet
      elemcounter(2,2)=elemcounter(2,2)+1
    CASE(105) !linear pyr
      elemcounter(3,2)=elemcounter(3,2)+1
    CASE(115) !non-linear pyr
      elemcounter(4,2)=elemcounter(4,2)+1
    CASE(205) !spline pyr
      elemcounter(5,2)=elemcounter(5,2)+1
    CASE(106) !linear prism
      elemcounter(6,2)=elemcounter(6,2)+1
    CASE(116) !non-linear prism
      elemcounter(7,2)=elemcounter(7,2)+1
    CASE(206) !spline prism
      elemcounter(8,2)=elemcounter(8,2)+1
    CASE(108) !linear hex
      elemcounter(9,2)=elemcounter(9,2)+1
    CASE(118) !non-linear hex
      elemcounter(10,2)=elemcounter(10,2)+1
    CASE(208) !spline hex
      elemcounter(11,2)=elemcounter(11,2)+1
  END SELECT
  ElemInfo(iElem,ELEM_Type)         = Elem%Type        ! Element Type
  ElemInfo(iElem,ELEM_Zone)         = Elem%Zone        ! Zone Number
  ElemInfo(iElem,ELEM_FirstSideInd) = iSide            ! first index -1 in SideInfo
  ElemInfo(iElem,ELEM_LastSideInd)  = iSide+locnSides  ! last index in SideInfo
  ElemInfo(iElem,ELEM_FirstNodeInd) = NodeID            ! first index -1 in NodeInfo
  ElemInfo(iElem,ELEM_LastNodeInd)  = NodeID+locnNodes  ! last index in NodeInfo
  NodeID = NodeID + locnNodes
  iSide = iSide + locnSides
  !approximate weight: locnNodes
  ElemWeight(iElem)=REAL(locnNodes)
END DO !iElem


!fill SideInfo
ALLOCATE(SideInfo(1:nTotalSides,SIDE_InfoSize)) 
SideInfo=0 
iSide=0
DO iElem=1,nElems
  Elem=>Elems(iElem)%ep
  DO iLocSide=1,6
    Side=>Elem%Side(iLocSide)%sp
    iSide=iSide+1
    !Side Tpye
    IF(ASSOCIATED(Elem%curvedNode))THEN
      SideInfo(iSide,SIDE_Type)=7            ! Side Type: NL quad
    ELSE
      SideInfo(iSide,SIDE_Type)=5            ! Side Type: bilinear
    END IF
    !Side ID
    SideInfo(iSide,SIDE_ID)=Side%ind
    IF(Side%flip.GT.0) SideInfo(iSide,SIDE_ID)=-SideInfo(iSide,SIDE_ID)           ! side is slave side
    !neighbor Element ID
    IF(ASSOCIATED(Side%Connection))THEN
      SideInfo(iSide,SIDE_nbElemID)=Side%Connection%Elem%ind                   ! Element ID of neighbor Element
    END IF
    !BC ID 
    SideInfo(iSide,SIDE_BCID)=Side%BCIndex                            
  END DO !iLocSide=1,6
END DO !iElem=1,nElems


!fill NodeInfo
ALLOCATE(NodeInfo(1:nTotalNodes))
NodeInfo=-1
NodeID=0
DO iElem=1,nElems
  Elem=>Elems(iElem)%ep
  DO iNode=1,8
    NodeID=NodeID+1
    NodeInfo(NodeID)=Elem%Node(iNode)%np%ind
  END DO
  DO iNode=1,nCurvedNodes
    NodeID=NodeID+1
    NodeInfo(NodeID)=Elem%curvedNode(iNode)%np%ind
  END DO
  DO iLocSide=1,6
    Side=>Elem%Side(iLocSide)%sp
    NodeID=NodeID+1
    IF(side%flip.EQ.0)THEN
      NodeInfo(NodeID)=Side%Node(1)%np%ind
    ELSE 
      NodeInfo(NodeID)=Side%Node(Side%flip)%np%ind
    END IF
  END DO !iLocSide=1,6
END DO !iElem=1,nElems

IF(NodeID.NE.nTotalNodes) &
       CALL abort(__STAMP__,&
          'Sanity check: nNodes not equal to number of nodes in NodeInfo!')

ALLOCATE(NodeCoords(nNodeIDs,3))
NodeCoords=-999.

!fill NodeCoords and Sanity check of curvednodes
DO iElem=1,nElems
  Elem=>Elems(iElem)%ep
  DO iNode=1,8
    Elem%Node(iNode)%np%tmp=1
  END DO
  ! CURVED  
  DO iNode=1,nCurvedNodes
    Elem%CurvedNode(iNode)%np%tmp=1
  END DO
END DO !iElem=1,nElems
NodeID=0
DO iElem=1,nElems
  Elem=>Elems(iElem)%ep
  DO iNode=1,8
    IF(Elem%Node(iNode)%np%tmp.GT.0)THEN
      iNodeP=Elem%Node(iNode)%np%ind
      NodeCoords(iNodeP,:)=Elem%Node(iNode)%np%x
      NodeID=NodeID+1
      Elem%Node(iNode)%np%tmp=0
    END IF                  
  END DO
  ! CURVED  
  DO iNode=1,nCurvedNodes
    IF(Elem%CurvedNode(iNode)%np%tmp.GT.0)THEN
      iNodeP=Elem%CurvedNode(iNode)%np%ind
      NodeCoords(iNodeP,:)=Elem%CurvedNode(iNode)%np%x
      NodeID=NodeID+1
      Elem%CurvedNode(iNode)%np%tmp=0
    END IF                  
  END DO
END DO !iElem=1,nElems

IF(NodeID.NE.nNodeIDs) & 
         CALL abort(__STAMP__,&
                'Sanity check: nNodeIDs not equal to number of nodes in NodeCoords!')

END SUBROUTINE getMeshinfo


END MODULE MOD_Output_HDF5
