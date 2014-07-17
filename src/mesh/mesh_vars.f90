#include "hopest_f.h"
MODULE MOD_Mesh_Vars
!===================================================================================================================================
! Contains global variables provided by the mesh routines
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PUBLIC
SAVE
!-----------------------------------------------------------------------------------------------------------------------------------
! basis
!-----------------------------------------------------------------------------------------------------------------------------------
LOGICAL           :: useCurveds
INTEGER           :: NGeo                        ! polynomial degree of geometric transformation
REAL,ALLOCATABLE  :: Xi_NGeo(:)                  ! 1D equidistant point positions for curved elements (during readin)
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES 
!-----------------------------------------------------------------------------------------------------------------------------------
INTEGER,ALLOCATABLE :: BoundaryType(:,:)
!-----------------------------------------------------------------------------------------------------------------------------------
INTEGER          :: nGlobalElems=0      ! number of elements in mesh
INTEGER          :: nElems=0            ! number of local elements
INTEGER          :: nEdges=0            ! number of unique edges
INTEGER          :: nSides=0            ! =nInnerSides+nBCSides+nMPISides
INTEGER          :: nMortarSides=0      ! 
INTEGER          :: nBCSides=0          ! BCSide index range: sideID \in [1:nBCSides]
INTEGER          :: nNodes=0            ! SIZE of Nodes pointer array, number of unique nodes
INTEGER          :: nBCs=0              ! number of BCs in mesh
!-----------------------------------------------------------------------------------------------------------------------------------
CHARACTER(LEN=255),ALLOCATABLE   :: BoundaryName(:)
CHARACTER(LEN=255)               :: MeshFile        ! name of hdf5 meshfile (write with ending .h5!)
!-----------------------------------------------------------------------------------------------------------------------------------
! USER DEFINED TYPES 
TYPE tNodePtr
  TYPE(tNode),POINTER          :: np                     ! node pointer
END TYPE tNodePtr

TYPE tEdgePtr
  LOGICAL                      :: isOriented
  TYPE(tEdge),POINTER          :: edp                    ! edge pointer
END TYPE tEdgePtr

TYPE tSidePtr
  TYPE(tSide),POINTER          :: sp                     ! side pointer
END TYPE tSidePtr

TYPE tElemPtr
  TYPE(tElem),POINTER          :: ep                     ! Local element pointer
END TYPE tElemPtr

TYPE tElem
  INTEGER                      :: ind             ! global element index
  INTEGER                      :: Type            ! element type (linear/bilinear/curved)
  INTEGER                      :: Zone
  !INTEGER                      :: which_tree
  !INTEGER                      :: quadrant_id
  TYPE(tNodePtr)               :: Node(8)
  TYPE(tSidePtr)               :: Side(6)
  TYPE(tEdgePtr)               :: Edge(12)
  TYPE(tNodePtr),POINTER       :: CurvedNode(:,:,:)
END TYPE tElem

TYPE tSide
  INTEGER                      :: ind             ! global side ID 
  INTEGER                      :: sideID          ! local side ID on Proc 
  INTEGER                      :: tmp 
  INTEGER                      :: NbProc 
  INTEGER                      :: BCindex         ! index in BoundaryType array! 
  INTEGER                      :: flip 
  INTEGER                      :: nMortars        ! number of slave mortar sides associated with master mortar
  INTEGER                      :: MortarType      ! type of mortar: Type1 : 1-4 , Type 2: 1-2 in eta, Type 2: 1-2 in xi
  TYPE(tNodePtr)               :: Node(4)
  TYPE(tEdgePtr)               :: Edge(4)
  TYPE(tElem),POINTER          :: Elem
  TYPE(tSide),POINTER          :: connection
  TYPE(tSidePtr),POINTER       :: MortarSide(:)   ! array of side pointers to slave mortar sides
END TYPE tSide

TYPE tEdge ! provides data structure for local edge
  INTEGER                      :: ind
  INTEGER                      :: tmp
  INTEGER                      :: nNbElems
  TYPE(tNodePtr)               :: Node(2)                ! pointer to node always 2
  TYPE(tEdge),POINTER          :: nextEdge               ! only used to assign edges 
  TYPE(tElemPtr),POINTER       :: nbElem(:)
END TYPE tEdge

TYPE tNode
  TYPE(tEdge),POINTER          :: firstEdge     ! only used to assign edges 
  INTEGER                      :: ind=0         ! global unique node index
  REAL                         :: x(3)=0.
END TYPE tNode
!-----------------------------------------------------------------------------------------------------------------------------------
TYPE(tElemPtr),POINTER         :: Elems(:)
TYPE(tNodePtr),POINTER         :: Nodes(:)
TYPE(tEdgePtr),POINTER         :: Edges(:)
!-----------------------------------------------------------------------------------------------------------------------------------
INTEGER,PARAMETER              :: EdgeToElemNode(1:2,1:12) = RESHAPE((/ 1, 2,&  ! CGNS corner nodes mapped 
                                                                        4, 3,&  ! to p4est edges
                                                                        5, 6,&
                                                                        8, 7,&
                                                                        1, 4,&
                                                                        2, 3,&
                                                                        5, 8,&
                                                                        6, 7,&
                                                                        1, 5,&
                                                                        2, 6,&
                                                                        4, 8,&
                                                                        3, 7 /),(/2,12/))
!-----------------------------------------------------------------------------------------------------------------------------------
#ifdef MPI
#endif /*MPI*/
!-----------------------------------------------------------------------------------------------------------------------------------
LOGICAL          :: MeshInitIsDone =.FALSE.
!===================================================================================================================================
INTERFACE GETNEWEDGE
  MODULE PROCEDURE GETNEWEDGE
END INTERFACE

INTERFACE GETNEWSIDE
  MODULE PROCEDURE GETNEWSIDE
END INTERFACE

INTERFACE GETNEWELEM
  MODULE PROCEDURE GETNEWELEM
END INTERFACE

INTERFACE deleteMeshPointer
  MODULE PROCEDURE deleteMeshPointer
END INTERFACE

CONTAINS


FUNCTION GETNEWEDGE(Node1,Node2,isNew,isOriented)
!===================================================================================================================================
! Create "Edge" with nodes "Node1" and "Node2"
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
TYPE(tNode),POINTER :: Node1, Node2 ! Node pointers
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
LOGICAL             :: isNew
LOGICAL             :: isOriented
TYPE(tEdge),POINTER :: GETNEWEDGE   ! New edge
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
TYPE(tNode),POINTER :: aNode,bNode
TYPE(tEdge),POINTER :: aEdge
LOGICAL             :: edgeFound
!===================================================================================================================================
IF((Node2%ind).GT.(Node1%ind)) THEN
  aNode=>Node1
  bNode=>Node2
  isOriented=.FALSE.
ELSEIF((Node2%ind).LT.(Node1%ind))THEN
  aNode=>Node2
  bNode=>Node1
  isOriented=.TRUE.
ELSE 
  WRITE(*,*) 'Problem with node%ind in GETNEWEDGE'
  WRITE(*,*) 'node IDs',Node1%ind,Node2%ind
  WRITE(*,*) 'node1%x',Node1%x
  WRITE(*,*) 'node2%x',Node2%x
  STOP
END IF

edgeFound=.FALSE.
aEdge=>aNode%firstEdge
DO WHILE (ASSOCIATED(aEdge))    
  IF(aEdge%Node(2)%np%ind .EQ. bNode%ind) THEN
    edgeFound=.TRUE.
    EXIT
  END IF
  aEdge=>aEdge%nextEdge
END DO
IF(edgeFound)THEN
  isNew=.FALSE.
ELSE
  ALLOCATE(aEdge)
  aEdge%Node(1)%np=>aNode
  aEdge%Node(2)%np=>bNode
  NULLIFY(aEdge%nextEdge)
  aEdge%ind=0
  aEdge%nNbElems=0
  ! sort new edge into edge list ( from the from)
  IF(ASSOCIATED(aNode%firstEdge)) THEN
    aEdge%nextEdge=>aNode%firstEdge 
  END IF
  aNode%firstEdge=>aEdge
  isNew=.TRUE.
END IF
GETNEWEDGE=>aEdge

END FUNCTION GETNEWEDGE


FUNCTION GETNEWSIDE()
!===================================================================================================================================
!  
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
TYPE(tSide),POINTER :: getNewSide
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER             :: iNode
!===================================================================================================================================
ALLOCATE(getNewSide)
DO iNode=1,4
  NULLIFY(getNewSide%Node(iNode)%np)
END DO
NULLIFY(getNewSide%Elem)
NULLIFY(getNewSide%MortarSide)
NULLIFY(getNewSide%connection)
getNewSide%sideID=0
getNewSide%ind=0
getNewSide%tmp=0
getNewSide%NbProc=-1
getNewSide%BCindex=0
getNewSide%flip=0
getNewSide%nMortars=0
getNewSide%MortarType=0
END FUNCTION GETNEWSIDE

FUNCTION GETNEWELEM()
!===================================================================================================================================
! 
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
TYPE(tElem),POINTER :: getNewElem
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER             :: iNode,iLocSide
!===================================================================================================================================
ALLOCATE(getNewElem)
DO iNode=1,8
  NULLIFY(getNewElem%Node(iNode)%np)
END DO
DO iLocSide=1,6
  getNewElem%Side(iLocSide)%sp=>getNewSide()
END DO
NULLIFY(getNewElem%CurvedNode)
getNewElem%ind=0
getNewElem%Zone=0
getNewElem%Type=0
END FUNCTION GETNEWELEM


SUBROUTINE createSides(Elem)
!===================================================================================================================================
! if element nodes already assigned, create Sides using CGNS standard
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
TYPE(tElem),POINTER :: Elem
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
!side 1
Elem%Side(1)%sp%Node(1)%np=>Elem%Node(1)%np
Elem%Side(1)%sp%Node(2)%np=>Elem%Node(4)%np
Elem%Side(1)%sp%Node(3)%np=>Elem%Node(3)%np
Elem%Side(1)%sp%Node(4)%np=>Elem%Node(2)%np
!side 2                                    
Elem%Side(2)%sp%Node(1)%np=>Elem%Node(1)%np
Elem%Side(2)%sp%Node(2)%np=>Elem%Node(2)%np
Elem%Side(2)%sp%Node(3)%np=>Elem%Node(6)%np
Elem%Side(2)%sp%Node(4)%np=>Elem%Node(5)%np
!side 3                                    
Elem%Side(3)%sp%Node(1)%np=>Elem%Node(2)%np
Elem%Side(3)%sp%Node(2)%np=>Elem%Node(3)%np
Elem%Side(3)%sp%Node(3)%np=>Elem%Node(7)%np
Elem%Side(3)%sp%Node(4)%np=>Elem%Node(6)%np
!side 4                                    
Elem%Side(4)%sp%Node(1)%np=>Elem%Node(3)%np
Elem%Side(4)%sp%Node(2)%np=>Elem%Node(4)%np
Elem%Side(4)%sp%Node(3)%np=>Elem%Node(8)%np
Elem%Side(4)%sp%Node(4)%np=>Elem%Node(7)%np
!side 5                                    
Elem%Side(5)%sp%Node(1)%np=>Elem%Node(1)%np
Elem%Side(5)%sp%Node(2)%np=>Elem%Node(5)%np
Elem%Side(5)%sp%Node(3)%np=>Elem%Node(8)%np
Elem%Side(5)%sp%Node(4)%np=>Elem%Node(4)%np
!side 6                                                
Elem%Side(6)%sp%Node(1)%np=>Elem%Node(5)%np
Elem%Side(6)%sp%Node(2)%np=>Elem%Node(6)%np
Elem%Side(6)%sp%Node(3)%np=>Elem%Node(7)%np
Elem%Side(6)%sp%Node(4)%np=>Elem%Node(8)%np
END SUBROUTINE createSides


SUBROUTINE deleteMeshPointer()
!===================================================================================================================================
! Deallocates all pointers used for the mesh readin
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER             :: iElem,iLocSide,iNode,nAssocNodes
INTEGER             :: iMortar
TYPE(tElem),POINTER :: aElem
TYPE(tSide),POINTER :: aSide
!===================================================================================================================================
DO iElem=1,nElems
  aElem=>Elems(iElem)%ep
  IF(ASSOCIATED(aElem%CurvedNode)) DEALLOCATE(aElem%curvedNode)
  DO iLocSide=1,6
    aSide=>aElem%Side(iLocSide)%sp
    DEALLOCATE(aSide)
  END DO
  DEALLOCATE(aElem)
END DO
DEALLOCATE(Elems)
nAssocNodes=0
DO iNode=1,nNodes
  IF(ASSOCIATED(Nodes(iNode)%np))THEN
    DEALLOCATE(Nodes(iNode)%np)
    nAssocNodes=nAssocNodes+1
  END IF
END DO
!WRITE(*,*)'DEBUG,nAssocNodes',nAssocNodes
DEALLOCATE(Nodes)
END SUBROUTINE deleteMeshPointer


END MODULE MOD_Mesh_Vars
