#include "hopest_f.h"
MODULE MODH_Output_HDF5
!===================================================================================================================================
! Add comments please!
!===================================================================================================================================
! MODULES
USE MODH_Globals
USE HDF5
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES 
INTEGER(HID_T)                 :: File_ID
!-----------------------------------------------------------------------------------------------------------------------------------
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
USE MODH_Mesh_Vars
USE MODH_IO_HDF5
USE MODH_HDF5_output
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
TYPE(tElem),POINTER            :: master
TYPE(tSide),POINTER            :: Side
INTEGER,ALLOCATABLE            :: ElemInfo(:,:),SideInfo(:,:),NodeInfo(:)
REAL,ALLOCATABLE               :: ElemBary(:,:)
REAL,ALLOCATABLE               :: NodeCoords(:,:)
REAL,ALLOCATABLE               :: ElemWeight(:)
INTEGER                        :: iQuad,i,j,k
INTEGER                        :: NodeID,iNode
INTEGER                        :: iSide,SideID,iLocSide
INTEGER                        :: ElemCounter(11,2)
INTEGER                        :: nSideIDs,nNodeIDs
INTEGER                        :: nTotalSides,nTotalNodes
INTEGER                        :: locnSides,locnNodes,offsetID
!===================================================================================================================================
WRITE(*,'(132("~"))')
WRITE(*,'(A)')' WRITE MESH TO HDF5 FILE... ' // TRIM(FileString) 

!set all side indices =0
DO iQuad=1,nQuads
  Elem=>Quads(iQuad)%ep
  DO iLocSide=1,6
    Side=>Elem%Side(iLocSide)%sp
    Side%ind=0
  END DO !iLocSide=1,6
END DO !iQuad=1,nQuads

! count Elements , unique sides 
nSideIDs=0 !number of unique side IDs (side and side%connection have the same sideID)


DO iQuad=1,nQuads
  Elem=>Quads(iQuad)%ep

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
END DO

!set unique nodes and Side Indices
SideID=0
NodeID=0
DO iQuad=1,nQuads
  Elem=>Quads(iQuad)%ep
  Elem%ind=iQuad

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

! start output
CALL OpenHDF5File(FileString,create=.TRUE.,single=.TRUE.)  


!-----------------------------------------------------------------
!attributes 
!-----------------------------------------------------------------

CALL WriteAttributeToHDF5(File_ID,'BoundaryOrder',1,IntegerScalar=Ngeo+1)
CALL WriteAttributeToHDF5(File_ID,'CurvedFound',1,LogicalScalar=useCurveds)

!-----------------------------------------------------------------
! WRITE BC 
!-----------------------------------------------------------------
CALL WriteArrayToHDF5(File_ID,'BCNames',nBCs,1,(/nBCs/),0,StrArray=BoundaryName)
CALL WriteArrayToHDF5(File_ID,'BCType',nBCs,2,(/nBcs,4/),0,IntegerArray=BoundaryType)

!-----------------------------------------------------------------
! Barycenters
!-----------------------------------------------------------------
WRITE(*,*)'WRITE Barycenters'

ALLOCATE(ElemBary(nQuads,3))
DO iQuad=1,nQuads
  ElemBary(iQuad,1)=SUM(XgeoQuad(1,:,:,:,iQuad))
  ElemBary(iQuad,2)=SUM(XgeoQuad(2,:,:,:,iQuad))
  ElemBary(iQuad,3)=SUM(XgeoQuad(3,:,:,:,iQuad))
END DO !iQuad=1,nElem
ElemBary(:,:)=ElemBary(:,:)*(1./(Ngeo+1)**3)

CALL WriteArrayToHDF5(File_ID,'ElemBarycenters',nQuads,2,(/nQuads,3/),0,RealArray=ElemBary)
DEALLOCATE(ElemBary)

!-----------------------------------------------------------------
! WRITE NodeCoords  for each element !!!! (multiple nodes!!!)
!-----------------------------------------------------------------
nNodeIDs=(Ngeo+1)**3*nQuads
CALL WriteArrayToHDF5(File_ID,'NodeCoords',nNodeIDs,2,(/nNodeIDs,3/),0,  &
          RealArray=TRANSPOSE(RESHAPE(XgeoQuad,(/3,nNodeIDs/))) )
DEALLOCATE(XgeoQuad)


!-----------------------------------------------------------------
!fill ElementInfo. 
!-----------------------------------------------------------------
WRITE(*,*)'WRITE ElemInfo'

ALLOCATE(ElemInfo(1:nQuads,ELEM_InfoSize))
ElemInfo=0
Elemcounter=0
Elemcounter(:,1)=(/104,204,105,115,205,106,116,206,108,118,208/)
iNode  = 0 
iSide  = 0
DO iQuad=1,nQuads
  Elem=>Quads(iQuad)%ep
  locnNodes=8+nCurvedNodes
  locnSides=0
  ! for element sides
  DO iLocSide=1,6
    Side=>Elem%Side(iLocSide)%sp
    locnNodes = locnNodes+1 ! wirte only first oriented node of side, if curved all!         
    locnSides = locnSides+1

    !+ MORTAR SIDES  !?!?

  END DO !iLocSide=1,6
  ElemInfo(iQuad,ELEM_Type)         = Elem%Type        ! Element Type
  ElemInfo(iQuad,ELEM_Zone)         = Elem%Zone        ! Zone Number
  ElemInfo(iQuad,ELEM_FirstSideInd) = iSide            ! first index -1 in SideInfo
  ElemInfo(iQuad,ELEM_LastSideInd)  = iSide+locnSides  ! last index in SideInfo
  ElemInfo(iQuad,ELEM_FirstNodeInd) = iNode            ! first index -1 in NodeInfo
  ElemInfo(iQuad,ELEM_LastNodeInd)  = iNode+locnNodes  ! last index in NodeInfo
  iNode = iNode + locnNodes
  iSide = iSide + locnSides
  !approximate weight: locnNodes
  CALL AddToElemCounter(Elem%type,ElemCounter)
END DO !iQuad
nTotalNodes = iNode
nTotalSides = iSide

!WRITE ElemInfo,into (1,nQuads)  
CALL WriteArrayToHDF5(File_ID,'ElemInfo',nQuads,2,(/nQuads,ELEM_InfoSize/),0,IntegerArray=ElemInfo)

DEALLOCATE(ElemInfo)

CALL WriteArrayToHDF5(File_ID,'ElemCounter',11,2,(/11,2/),0,IntegerArray=ElemCounter)
WRITE(*,*)'Mesh statistics:'
WRITE(*,*)'Element Type | number of elements'
DO i=1,11
  WRITE(*,'(I4,A,I8)') Elemcounter(i,1),'        | ',Elemcounter(i,2)
END DO

!-----------------------------------------------------------------
! element weights
!-----------------------------------------------------------------
!WRITE ElemWeight,into (1,nQuads)  
ALLOCATE(ElemWeight(1:nQuads))
ElemWeight=1.
CALL WriteArrayToHDF5(File_ID,'ElemWeight',nQuads,1,(/nQuads/),0,RealArray=ElemWeight)
DEALLOCATE(ElemWeight)



!-----------------------------------------------------------------
!fill SideInfo
!-----------------------------------------------------------------
WRITE(*,*)'WRITE SideInfo'

ALLOCATE(SideInfo(1:nTotalSides,SIDE_InfoSize)) 
SideInfo=0 
iSide=0
DO iQuad=1,nQuads
  Elem=>Quads(iQuad)%ep
  DO iLocSide=1,6
    Side=>Elem%Side(iLocSide)%sp
    iSide=iSide+1
    !Side Tpye
    IF(Ngeo.GT.1)THEN
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
END DO !iQuad=1,nQuads

!WRITE SideInfo,into (1,nTotalSides)   
CALL WriteArrayToHDF5(File_ID,'SideInfo',nTotalSides,2,(/nTotalSides,SIDE_InfoSize/),0,IntegerArray=SideInfo)
DEALLOCATE(SideInfo)

!-----------------------------------------------------------------
!fill NodeInfo
!-----------------------------------------------------------------
WRITE(*,*)'WRITE NodeInfo'

! since each element has its own nodes in NodeCoords ( = multiply defined nodes)
! thus each element has its own index range, from [ locNodes*(iQuad-1)+1 : locNodes*iQuad ]

! The node coordinates are written in tensor-product style i,j,k (i inner loop)
! and we deduct the cgns corner and side node mapping

locnNodes=8+6+nCurvedNodes  

master=>GETNEWELEM()
DO iNode=1,8
  ALLOCATE(master%Node(iNode)%np)
END DO
master%Node(1)%np%ind=HexMap(   0,   0,   0)
master%Node(2)%np%ind=HexMap(Ngeo,   0,   0)
master%Node(3)%np%ind=HexMap(Ngeo,Ngeo,   0)
master%Node(4)%np%ind=HexMap(   0,Ngeo,   0)
master%Node(5)%np%ind=HexMap(   0,   0,Ngeo)
master%Node(6)%np%ind=HexMap(Ngeo,   0,Ngeo)
master%Node(7)%np%ind=HexMap(Ngeo,Ngeo,Ngeo)
master%Node(8)%np%ind=HexMap(   0,Ngeo,Ngeo)
CALL createSides(master)

IF(nTotalNodes.NE.locnNodes*nQuads) &
       CALL abort(__STAMP__,&
          'Sanity check: nNodes not equal to locnNodes*nQuads!')



ALLOCATE(NodeInfo(1:nTotalNodes))
NodeInfo=-1
NodeID=0


offsetID=0
locnNodes=(Ngeo+1)**3
DO iQuad=1,nQuads
  Elem=>Quads(iQuad)%ep
  DO iNode=1,8
    NodeID=NodeID+1
    NodeInfo(NodeID)= master%Node(iNode)%np%ind + offsetID !Elem%Node(iNode)%np%ind
  END DO
  DO iNode=1,nCurvedNodes
    NodeID=NodeID+1
    NodeInfo(NodeID)=iNode + offsetID !Elem%curvedNode(iNode)%np%ind
  END DO
  DO iLocSide=1,6
    Side=>Elem%Side(iLocSide)%sp
    NodeID=NodeID+1
    IF(side%flip.EQ.0)THEN
      NodeInfo(NodeID)=master%Side(iLocSide)%sp%Node(1)%np%ind + offsetID !Side%Node(1)%np%ind
    ELSE 
      NodeInfo(NodeID)=master%Side(iLocSide)%sp%Node(Side%flip)%np%ind + offsetID !Side%Node(Side%flip)%np%ind
    END IF
  END DO !iLocSide=1,6
  offsetID=offsetID+locnNodes
END DO !iQuad=1,nQuads

IF(NodeID.NE.nTotalNodes) &
       CALL abort(__STAMP__,&
          'Sanity check: nNodes not equal to number of nodes in NodeInfo!')

!WRITE NodeInfo,into (1,nTotalNodes) 
CALL WriteArrayToHDF5(File_ID,'NodeInfo',nTotalNodes,1,(/nTotalNodes/),0,IntegerArray=NodeInfo)
DEALLOCATE(NodeInfo)



! Close the file.
CALL CloseHDF5File()

WRITE(*,'(A)')' DONE WRITING MESH.'
WRITE(*,'(132("~"))')

END SUBROUTINE WriteMeshToHDF5



SUBROUTINE AddtoElemCounter(elemtype,elemCounter) 
!===================================================================================================================================
!
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)       :: ElemType
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
INTEGER,INTENT(INOUT)    :: ElemCounter(11,2)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
  SELECT CASE(ElemType)
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
    CASE DEFAULT
      STOP 'elem type not defined in elemcounter'
  END SELECT
END SUBROUTINE AddToElemCounter


END MODULE MODH_Output_HDF5
