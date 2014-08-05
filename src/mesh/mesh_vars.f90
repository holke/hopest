#include "hopest_f.h"
MODULE MOD_Mesh_Vars
!===================================================================================================================================
! Contains global variables provided by the mesh routines
!===================================================================================================================================
! MODULES
USE MOD_p4estBindingTypes
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
REAL,ALLOCATABLE  :: wBary_NGeo(:)               ! barycentric weights from xi_Ngeo
REAL,ALLOCATABLE  :: XGeo(:,:,:,:,:)              ! High order geometry nodes, per element (1:3,0:Ngeo,0:Ngeo,0:Ngeo,nElems)
REAL,ALLOCATABLE  :: XGeoQuad(:,:,:,:,:)              ! High order geometry nodes, per element (1:3,0:Ngeo,0:Ngeo,0:Ngeo,nElems)
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES 
!-----------------------------------------------------------------------------------------------------------------------------------
INTEGER,ALLOCATABLE :: BoundaryType(:,:)
!-----------------------------------------------------------------------------------------------------------------------------------
INTEGER             :: nGlobalElems=0      ! number of elements in mesh
INTEGER             :: nElems=0            ! number of local elements
INTEGER             :: nSides=0            ! =nInnerSides+nBCSides
INTEGER             :: nMortarSides=0      ! 
INTEGER             :: nBCSides=0          ! BCSide index range: sideID \in [1:nBCSides]
INTEGER             :: nNodes=0            ! SIZE of Nodes pointer array, number of unique nodes
INTEGER             :: nBCs=0              ! number of BCs in mesh
INTEGER             :: nUserBCs=0     
INTEGER             :: nCurvedNodes=0      ! number of curved nodes per element = (Ngeo+1)^3
!-----------------------------------------------------------------------------------------------------------------------------------
CHARACTER(LEN=255),ALLOCATABLE   :: BoundaryName(:)
CHARACTER(LEN=255)               :: MeshFile        ! name of hdf5 meshfile (write with ending .h5!)
INTEGER                          :: refineLevel
INTEGER                          :: refineType
INTEGER,ALLOCATABLE              :: RefineList(:)
!-----------------------------------------------------------------------------------------------------------------------------------
! USER DEFINED TYPES 
TYPE tNodePtr
  TYPE(tNode),POINTER          :: np                     ! node pointer
END TYPE tNodePtr

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
  INTEGER                      :: treeID
  !INTEGER                      :: quadrant_id
  TYPE(tNodePtr)               :: Node(8)
  TYPE(tSidePtr)               :: Side(6)
END TYPE tElem

TYPE tSide
  INTEGER                      :: ind             ! global side ID 
  INTEGER                      :: sideID          ! local side ID on Proc 
  INTEGER                      :: tmp 
  INTEGER                      :: NbProc 
  INTEGER                      :: BCindex         ! index in BoundaryType array! 
  INTEGER                      :: flip 
  INTEGER                      :: nMortars        ! number of slave mortar sides associated with master mortar
  INTEGER                      :: MortarType      ! type of mortar: Type1 : 1-4 , Type 2: 1-2 in eta, Type 3: 1-2 in xi
  !TODO: if tSide is small side of a mortar group, mortar type is -1
  TYPE(tNodePtr)               :: Node(4)
  TYPE(tElem),POINTER          :: Elem
  TYPE(tSide),POINTER          :: connection
  TYPE(tSidePtr),POINTER       :: MortarSide(:)   ! array of side pointers to slave mortar sides
END TYPE tSide

TYPE tNode
  INTEGER                      :: ind=0         ! global unique node index
  INTEGER                      :: tmp=0
  !REAL                         :: x(3)=0.
END TYPE tNode
!-----------------------------------------------------------------------------------------------------------------------------------
TYPE(tElemPtr),POINTER         :: Elems(:)
TYPE(tNodePtr),POINTER         :: Nodes(:)
INTEGER,ALLOCATABLE            :: HexMap(:,:,:)
INTEGER,ALLOCATABLE            :: HexMapInv(:,:)
!-----------------------------------------------------------------------------------------------------------------------------------
! P4EST related data structures 
!-----------------------------------------------------------------------------------------------------------------------------------
TYPE(t_p4est_ptr) :: p4est_ptr              ! c pointer derived data type, see MOD_P4estBindingTypes

TYPE(tElemPtr),ALLOCATABLE  :: Quads(:)           ! new element list elements are "quadrants/octants"        
INTEGER                     :: nQuadrants         ! local number of quadrants (here no MPI => all) 
INTEGER                     :: nHalfFaces         ! number of mortar sides
INTEGER(KIND=4)             :: IntSize            ! used to transform INT coords/levels to REAL coords/levels: REAL=1/inssize*INT  [0. ; 1.]
REAL                        :: sIntSize           ! 1./REAL(intsize)
INTEGER(KIND=4),POINTER     :: QuadToTree(:)      ! from quadrant to tree ( ~ new element ID to old element ID) 
INTEGER(KIND=4),ALLOCATABLE :: TreeToQuad(:,:) ! from tree to quad range (2,nElems), entry 1: firstInd-1, entry2:lastInd 
INTEGER(KIND=4),POINTER     :: QuadToQuad(:,:)    ! p4est quadrant connectivity (1:6,1:nQuadrants) => neighbor quadrant
INTEGER(KIND=1),POINTER     :: QuadToFace(:,:)    ! p4est face connectivity (1:6,1:nQuadrants) => neighbor faceId + orientation + non-conform info
INTEGER(KIND=4),POINTER     :: QuadToHalf(:,:)    ! p4est face connectivity for mortars (1:4,1:nHalfFaces), ( ~small sides)
INTEGER(KIND=4),ALLOCATABLE :: QuadCoords(:,:)    ! p4est Integer coordinates of first quadrant node (xyz,nQuadrants)
INTEGER(KIND=1),ALLOCATABLE :: QuadLevel(:)       ! p4est Integer Level of quadrant (use to compute quadrant size

!-----------------------------------------------------------------------------------------------------------------------------------
INTEGER,PARAMETER   :: EdgeToElemNode(1:2,1:12) = RESHAPE((/ 1, 2,&  ! CGNS corner nodes mapped 
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
INTEGER,PARAMETER   :: H2P_FaceMap(1:6)     =  (/4,2,1,3,0,5/)     !mapping from local face order (CGNS) to p4est face
INTEGER,PARAMETER   :: P2H_FaceMap(0:5)     =  (/5,3,2,4,1,6/)     !mapping from local face order (CGNS) to p4est face
INTEGER,PARAMETER   :: H2P_VertexMap(1:8)   =  (/0,1,3,2,4,5,7,6/) !mapping from local node order (CGNS) to p4est node order 
INTEGER,PARAMETER   :: P2H_VertexMap(0:7)   =  (/1,2,4,3,5,6,8,7/) !mapping from local node order (CGNS) to p4est node order 

! mapping from HOPEST node of local sides to P4EST nodes of local sides
INTEGER,PARAMETER   :: H2P_FaceNodeMap(1:4,1:6) = &
                                      RESHAPE((/ 0,2,3,1,&
                                                 0,1,3,2,&
                                                 0,1,3,2,&
                                                 1,0,2,3,&
                                                 0,2,3,1,&
                                                 0,1,3,2 /),(/4,6/))

! mapping from P4EST node of local sides to HOPEST node of local sides
INTEGER,PARAMETER   :: P2H_FaceNodeMap(0:3,0:5) = &
                                      RESHAPE((/ 1,4,2,3,&
                                                 1,2,4,3,&
                                                 1,2,4,3,&
                                                 2,1,3,4,&
                                                 1,4,2,3,&
                                                 1,2,4,3 /),(/4,6/))

! Mapping matrices for computation of same node on adjacent face, see paper Burstedde p4est, 2011
! Node1= P4P(P4Q(P4R(Face0,Face1),orientation),Node0)
INTEGER,PARAMETER   :: P4R(0:5,0:5) = TRANSPOSE(RESHAPE((/ 0,1,1,0,0,1,&
                                                           2,0,0,1,1,0,&
                                                           2,0,0,1,1,0,&
                                                           0,2,2,0,0,1,&
                                                           2,0,0,2,2,0,&
                                                           2,0,0,2,2,0 /),(/6,6/)))

INTEGER,PARAMETER   :: P4Q(0:2,0:3) = TRANSPOSE(RESHAPE((/ 1,2,5,6,&
                                                           0,3,4,7,&
                                                           0,4,3,7 /),(/4,3/)))

INTEGER,PARAMETER   :: P4P(0:7,0:3) = TRANSPOSE(RESHAPE((/ 0,1,2,3,&
                                                           0,2,1,3,&
                                                           1,0,3,2,&
                                                           1,3,0,2,&
                                                           2,0,3,1,&
                                                           2,3,0,1,&
                                                           3,1,2,0,&
                                                           3,2,1,0 /),(/4,8/)))

!-----------------------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------------------
LOGICAL          :: MeshInitIsDone =.FALSE.
!===================================================================================================================================
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
Elem%Side(1)%sp%elem=>Elem
!side 2                                    
Elem%Side(2)%sp%Node(1)%np=>Elem%Node(1)%np
Elem%Side(2)%sp%Node(2)%np=>Elem%Node(2)%np
Elem%Side(2)%sp%Node(3)%np=>Elem%Node(6)%np
Elem%Side(2)%sp%Node(4)%np=>Elem%Node(5)%np
Elem%Side(2)%sp%elem=>Elem
!side 3                                    
Elem%Side(3)%sp%Node(1)%np=>Elem%Node(2)%np
Elem%Side(3)%sp%Node(2)%np=>Elem%Node(3)%np
Elem%Side(3)%sp%Node(3)%np=>Elem%Node(7)%np
Elem%Side(3)%sp%Node(4)%np=>Elem%Node(6)%np
Elem%Side(3)%sp%elem=>Elem
!side 4                                    
Elem%Side(4)%sp%Node(1)%np=>Elem%Node(3)%np
Elem%Side(4)%sp%Node(2)%np=>Elem%Node(4)%np
Elem%Side(4)%sp%Node(3)%np=>Elem%Node(8)%np
Elem%Side(4)%sp%Node(4)%np=>Elem%Node(7)%np
Elem%Side(4)%sp%elem=>Elem
!side 5                                    
Elem%Side(5)%sp%Node(1)%np=>Elem%Node(1)%np
Elem%Side(5)%sp%Node(2)%np=>Elem%Node(5)%np
Elem%Side(5)%sp%Node(3)%np=>Elem%Node(8)%np
Elem%Side(5)%sp%Node(4)%np=>Elem%Node(4)%np
Elem%Side(5)%sp%elem=>Elem
!side 6                                                
Elem%Side(6)%sp%Node(1)%np=>Elem%Node(5)%np
Elem%Side(6)%sp%Node(2)%np=>Elem%Node(6)%np
Elem%Side(6)%sp%Node(3)%np=>Elem%Node(7)%np
Elem%Side(6)%sp%Node(4)%np=>Elem%Node(8)%np
Elem%Side(6)%sp%elem=>Elem
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
TYPE(tElem),POINTER :: aElem
TYPE(tSide),POINTER :: aSide
!===================================================================================================================================
DO iElem=1,nElems
  aElem=>Elems(iElem)%ep
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
