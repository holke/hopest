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
USE MODH_ChangeBasis, ONLY:ChangeBasis3D
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
REAL,ALLOCATABLE               :: ElemWeight(:)
INTEGER                        :: iElem,i
INTEGER                        :: NodeID,iNode
INTEGER                        :: iSide,SideID,iLocSide,iMortar
INTEGER                        :: ElemCounter(11,2)
INTEGER                        :: nSideIDs,nNodeIDs
INTEGER                        :: nTotalSides,nTotalNodes
INTEGER                        :: locnSides,locnNodes,offsetID
!===================================================================================================================================
WRITE(*,'(132("~"))')
WRITE(*,'(A)')' WRITE MESH TO HDF5 FILE... ' // TRIM(FileString) 

!set all side indices =0
DO iElem=1,nElems
  Elem=>Elems(iElem)%ep
  DO iLocSide=1,6
    Side=>Elem%Side(iLocSide)%sp
    Side%ind=0
  END DO !iLocSide=1,6
END DO !iElem=1,nElems

! count Elements , unique sides 
nSideIDs=0 !number of unique side IDs (side and side%connection have the same sideID)


DO iElem=1,nElems
  Elem=>Elems(iElem)%ep

  ! Count sides
  DO iLocSide=1,6
    Side=>Elem%Side(iLocSide)%sp
    IF(Side%ind.EQ.0) THEN
      IF(Side%MortarType.EQ.0)THEN
        nSideIDs=nSideIDs+1
        Side%ind=-88888
        IF(ASSOCIATED(Side%connection))THEN      
          IF(Side%connection%ind.EQ.0) nSideIDs=nSideIDs-1 ! count inner and periodic sides only once 
        END IF
      ELSEIF(Side%MortarType.GT.0)THEN
        nSideIDs=nSideIDs+1
        Side%ind=-88888
        DO iMortar=1,Side%nMortars
          IF(Side%MortarSide(iMortar)%sp%ind.EQ.0)THEN
            nSideIDs=nSideIDs+1
            Side%MortarSide(iMortar)%sp%ind=-88888
          END IF 
        END DO !iMortar
      ELSE
        nSideIDs=nSideIDs+1
        Side%ind=-88888
      END IF
    END IF
  END DO
END DO

WRITE(*,*)'nSideIDs',nSideIDs

!set unique nodes and Side Indices
SideID=0
NodeID=0
DO iElem=1,nElems
  Elem=>Elems(iElem)%ep
  Elem%ind=iElem
  DO iLocSide=1,6
    Side=>Elem%Side(iLocSide)%sp
    IF(side%ind.EQ.-88888) THEN  ! assign side ID only for non MPI sides and lower MPI sides
      IF(Side%MortarType.EQ.0)THEN
        SideID=SideID+1
        Side%ind=SideID
        IF(ASSOCIATED(Side%connection))THEN      
          IF(Side%connection%ind.EQ.-88888) Side%connection%ind=SideID ! count inner and periodic sides only once 
        END IF
      ELSEIF(Side%MortarType.GT.0)THEN
        SideID=SideID+1
        Side%ind=SideID
        DO iMortar=1,Side%nMortars
          IF(Side%MortarSide(iMortar)%sp%ind.EQ.-88888)THEN
            SideID=SideID+1
            Side%MortarSide(iMortar)%sp%ind=SideID
          END IF 
        END DO !iMortar
      ELSE
        SideID=SideID+1
        Side%ind=SideID
      END IF
    END IF
  END DO !iLocSide
END DO !Elem

DO iElem=1,nElems
  Elem=>Elems(iElem)%ep
  DO iLocSide=1,6
    Side=>Elem%Side(iLocSide)%sp
    IF((Side%ind.LE.0).OR.(Side%ind.GT.nSideIds)) STOP 'Problem with sideID assigment'
  END DO !iLocSide
END DO !Elem

IF(SideID.NE.nSideIDs) STOP' problem: SideID <> nSideIDs'

! start output
CALL OpenHDF5File(FileString,create=.TRUE.,single=.TRUE.)  

!-----------------------------------------------------------------
!attributes 
!-----------------------------------------------------------------
IF(Ngeo_out.EQ.1) THEN
  useCurveds=.FALSE.
  nCurvedNodes=0
ELSE
  nCurvedNodes=(Ngeo_out+1)**3
END IF

CALL WriteAttributeToHDF5(File_ID,'CurvedFound',1,LogicalScalar=useCurveds)
CALL WriteAttributeToHDF5(File_ID,'BoundaryOrder',1,IntegerScalar=Ngeo_out+1)

!-----------------------------------------------------------------
! WRITE BC 
!-----------------------------------------------------------------
CALL WriteArrayToHDF5(File_ID,'BCNames',nBCs,1,(/nBCs/),0,StrArray=BoundaryName)
CALL WriteArrayToHDF5(File_ID,'BCType',nBCs,2,(/nBcs,4/),0,IntegerArray=BoundaryType)

!-----------------------------------------------------------------
! Barycenters
!-----------------------------------------------------------------
WRITE(*,*)'WRITE Barycenters'

ALLOCATE(ElemBary(nElems,3))
DO iElem=1,nElems
  ElemBary(iElem,1)=SUM(XGeoElem(1,:,:,:,iElem))
  ElemBary(iElem,2)=SUM(XGeoElem(2,:,:,:,iElem))
  ElemBary(iElem,3)=SUM(XGeoElem(3,:,:,:,iElem))
END DO !iElem=1,nElem
ElemBary(:,:)=ElemBary(:,:)*(1./(Ngeo_out+1)**3)

CALL WriteArrayToHDF5(File_ID,'ElemBarycenters',nElems,2,(/nElems,3/),0,RealArray=ElemBary)
DEALLOCATE(ElemBary)

!-----------------------------------------------------------------
! WRITE NodeCoords  for each element !!!! (multiple nodes!!!)
!-----------------------------------------------------------------
!transform to equidistant nodes (overwrite!!!):
DO iElem=1,nElems
  CALL ChangeBasis3D(3,Ngeo_out,Ngeo_out,Vdm_CL_EQ_out,XgeoElem(:,:,:,:,iElem),XgeoElem(:,:,:,:,iElem))
END DO
nNodeIDs=(Ngeo_out+1)**3*nElems
CALL WriteArrayToHDF5(File_ID,'NodeCoords',nNodeIDs,2,(/nNodeIDs,3/),0,  &
          RealArray=TRANSPOSE(RESHAPE(XGeoElem,(/3,nNodeIDs/))) )
DEALLOCATE(XGeoElem)


!-----------------------------------------------------------------
!fill ElementInfo. 
!-----------------------------------------------------------------
WRITE(*,*)'WRITE ElemInfo'

ALLOCATE(ElemInfo(1:nElems,ELEM_InfoSize))
ElemInfo=0
Elemcounter=0
Elemcounter(:,1)=(/104,204,105,115,205,106,116,206,108,118,208/)
iNode  = 0 
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

    !Mortar
    SELECT CASE(Side%MortarType)
    CASE(0) !do nothing
    CASE(1)
      locnSides=locnSides+4
    CASE(2,3)
      locnSides=locnSides+2
    END SELECT

  END DO !iLocSide=1,6
  ElemInfo(iElem,ELEM_Type)         = Elem%Type        ! Element Type
  ElemInfo(iElem,ELEM_Zone)         = Elem%Zone        ! Zone Number
  ElemInfo(iElem,ELEM_FirstSideInd) = iSide            ! first index -1 in SideInfo
  ElemInfo(iElem,ELEM_LastSideInd)  = iSide+locnSides  ! last index in SideInfo
  ElemInfo(iElem,ELEM_FirstNodeInd) = iNode            ! first index -1 in NodeInfo
  ElemInfo(iElem,ELEM_LastNodeInd)  = iNode+locnNodes  ! last index in NodeInfo
  iNode = iNode + locnNodes
  iSide = iSide + locnSides
  !approximate weight: locnNodes
  CALL AddToElemCounter(Elem%type,ElemCounter)
END DO !iElem
nTotalNodes = iNode
nTotalSides = iSide

!WRITE ElemInfo,into (1,nElems)  
CALL WriteArrayToHDF5(File_ID,'ElemInfo',nElems,2,(/nElems,ELEM_InfoSize/),0,IntegerArray=ElemInfo)

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
!WRITE ElemWeight,into (1,nElems)  
ALLOCATE(ElemWeight(1:nElems))
ElemWeight=1.
CALL WriteArrayToHDF5(File_ID,'ElemWeight',nElems,1,(/nElems/),0,RealArray=ElemWeight)
DEALLOCATE(ElemWeight)



!-----------------------------------------------------------------
!fill SideInfo
!-----------------------------------------------------------------
WRITE(*,*)'WRITE SideInfo'

ALLOCATE(SideInfo(1:nTotalSides,SIDE_InfoSize)) 
SideInfo=0 
iSide=0
DO iElem=1,nElems
  Elem=>Elems(iElem)%ep
  DO iLocSide=1,6
    Side=>Elem%Side(iLocSide)%sp
    iSide=iSide+1
    !Side Tpye
    IF(Ngeo_out.GT.1)THEN
      SideInfo(iSide,SIDE_Type)=7            ! Side Type: NL quad
    ELSE
      SideInfo(iSide,SIDE_Type)=5            ! Side Type: bilinear
    END IF
    !Side ID
    SideInfo(iSide,SIDE_ID)=Side%ind
    IF(Side%flip.GT.0) SideInfo(iSide,SIDE_ID)=-SideInfo(iSide,SIDE_ID)           ! side is slave side
    !MORTAR
    IF(Side%MortarType.GT.0)THEN
      SideInfo(iSide,SIDE_nbElemID)=-Side%MortarType !marker for attached mortar sides
      SideInfo(iSide,SIDE_BCID)    = Side%BCIndex                            
      IF(Side%flip.NE.0) STOP 'Problem with flip on mortar'
      DO iMortar=1,Side%nMortars
        iSide=iSide+1
        !Side Tpye
        IF(Ngeo_out.GT.1)THEN
          SideInfo(iSide,SIDE_Type)=7            ! Side Type: NL quad
        ELSE
          SideInfo(iSide,SIDE_Type)=5            ! Side Type: bilinear
        END IF
        !Side ID
        SideInfo(iSide,SIDE_ID)      = Side%MortarSide(iMortar)%sp%ind
        SideInfo(iSide,SIDE_nbElemID)= Side%MortarSide(iMortar)%sp%Elem%ind      ! Element ID of neighbor Element
        SideInfo(iSide,SIDE_BCID)    = Side%MortarSide(iMortar)%sp%BCIndex                            
      END DO
    ELSE !no mortar side
      IF(ASSOCIATED(Side%Connection))THEN
        SideInfo(iSide,SIDE_nbElemID)=Side%Connection%Elem%ind                   ! Element ID of neighbor Element
      END IF
      SideInfo(iSide,SIDE_BCID)=Side%BCIndex                            
      IF(Side%MortarType.LT.0) SideInfo(iSide,SIDE_Type)=-SideInfo(iSide,SIDE_Type) !mark sides as belonging to a mortar
    END IF
  END DO !iLocSide=1,6
END DO !iElem=1,nElems

!WRITE SideInfo,into (1,nTotalSides)   
CALL WriteArrayToHDF5(File_ID,'SideInfo',nTotalSides,2,(/nTotalSides,SIDE_InfoSize/),0,IntegerArray=SideInfo)
DEALLOCATE(SideInfo)

!-----------------------------------------------------------------
!fill NodeInfo
!-----------------------------------------------------------------
WRITE(*,*)'WRITE NodeInfo'

! since each element has its own nodes in NodeCoords ( = multiply defined nodes)
! thus each element has its own index range, from [ locNodes*(iElem-1)+1 : locNodes*iElem ]

! The node coordinates are written in tensor-product style i,j,k (i inner loop)
! and we deduct the cgns corner and side node mapping

locnNodes=8+6+nCurvedNodes  

master=>GETNEWELEM()
DO iNode=1,8
  ALLOCATE(master%Node(iNode)%np)
END DO
master%Node(1)%np%ind=HexMap_Out(       0,       0,       0)
master%Node(2)%np%ind=HexMap_Out(Ngeo_out,       0,       0)
master%Node(3)%np%ind=HexMap_Out(Ngeo_out,Ngeo_out,       0)
master%Node(4)%np%ind=HexMap_Out(       0,Ngeo_out,       0)
master%Node(5)%np%ind=HexMap_Out(       0,       0,Ngeo_out)
master%Node(6)%np%ind=HexMap_Out(Ngeo_out,       0,Ngeo_out)
master%Node(7)%np%ind=HexMap_Out(Ngeo_out,Ngeo_out,Ngeo_out)
master%Node(8)%np%ind=HexMap_Out(       0,Ngeo_out,Ngeo_out)
CALL createSides(master)

IF(nTotalNodes.NE.locnNodes*nElems) &
       CALL abort(__STAMP__,&
          'Sanity check: nNodes not equal to locnNodes*nElems!')



ALLOCATE(NodeInfo(1:nTotalNodes))
NodeInfo=-1
NodeID=0


offsetID=0
locnNodes=(Ngeo_out+1)**3
DO iElem=1,nElems
  Elem=>Elems(iElem)%ep
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
END DO !iElem=1,nElems

IF(NodeID.NE.nTotalNodes) &
       CALL abort(__STAMP__,&
          'Sanity check: nNodes not equal to number of nodes in NodeInfo!')

!WRITE NodeInfo,into (1,nTotalNodes) 
CALL WriteArrayToHDF5(File_ID,'NodeInfo',nTotalNodes,1,(/nTotalNodes/),0,IntegerArray=NodeInfo)
DEALLOCATE(NodeInfo)
DO iLocSide=1,6
  DEALLOCATE(master%Side(iLocSide)%sp)
END DO
DO iNode=1,8
  DEALLOCATE(master%Node(iNode)%np)
END DO
DEALLOCATE(master)



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
