#include "defines.f90"
MODULE MOD_Output_HDF5
!===================================================================================================================================
! Add comments please!
!===================================================================================================================================
! MODULES
USE HDF5
USE MOD_IO_HDF5
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES 
!-----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------

INTERFACE WriteMeshToHDF5
  MODULE PROCEDURE WriteMeshToHDF5
END INTERFACE

PUBLIC::WriteMeshToHDF5
!===================================================================================================================================

CONTAINS

SUBROUTINE WriteMeshToHDF5(Mesh,FileString)
!===================================================================================================================================
! Subroutine to write Data to HDF5 format
!===================================================================================================================================
! MODULES
USE MOD_Mesh_Vars,ONLY:tElem,tSide
USE MOD_Mesh_Vars,ONLY:FirstElem
USE MOD_Mesh_Vars,ONLY:N
USE MOD_Output_Vars,ONLY:dosortIJK
USE MOD_Mesh_Vars,ONLY:nUserDefinedBoundaries,BoundaryName,BoundaryType
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
INTEGER                        :: ElemID,SideID,NodeID
INTEGER,ALLOCATABLE            :: IDlist(:)
REAL                           :: sBasisBary(4)
CHARACTER(LEN=255)             :: Meshfile
!===================================================================================================================================
WRITE(UNIT_stdOut,'(132("~"))')
CALL Timer(.TRUE.)
WRITE(UNIT_stdOut,'(A)')' WRITE DATA TO HDF5 FILE...'
! Create the file collectively.
CALL OpenHDF5File(FileString,create=.TRUE.)  

!set all node and side indices =0
DO iElem=1,nElems
  Elem=>Elems(iElem)%ep
  DO i=1,8
    Elem%Node(i)%np%ind=0
  END DO
  IF(useCurveds)THEN
    DO k=0,NGeo; DO j=0,NGeo; DO i=0,NGeo
      Elem%curvedNode(i,j,k)%np%ind=0
    END DO
  END IF
  DO iLocSide=1,6
    Side=>Elem%Side(iLocSide)%sp
    Side%ind=0
  END DO
END DO

! count Elements , unique sides and nodes are marked with ind=0
nNodeIDs=0 !number of unique nodeIDs
nSideIDs=0 !number of unique side IDs (side and side%connection have the same sideID)
nElems=0   !number of elements
nSides=0   !number of all sides
nNodes=0   !number of all nodes

DO iElem=1,nElems
  Elem=>Elems(iElem)%ep
  ! Count nodes
  DO i=1,8
    IF(Elem%Node(i)%np%ind.NE.0) CYCLE
    nNodeIDs=nNodeIDs+1
    Elem%Node(i)%np%ind=-88888  ! mark no MPI side
  END DO

  IF(useCurveds)THEN
    DO k=0,NGeo; DO j=0,NGeo; DO i=0,NGeo
      IF(Elem%CurvedNode(i,j,k)%np%ind.NE.0) CYCLE
      nNodeIDs=nNodeIDs+1
      Elem%CurvedNode(i,j,k)%np%ind=-88888
    END DO
  END IF

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
  nNodes = nNodes+8+6+nCurvedNodes ! corner + oriented + curved
  nSides = nSides+6
END DO

!set unique nodes and Side Indices
ElemID=0
SideID=0
NodeID=0
DO iElem=1,nElems
  Elem=>Elems(iElem)%ep
  ElemID=ElemID+1 
  Elem%ind=ElemID
  DO i=1,8
    IF(Elem%Node(i)%np%ind.NE.-88888) CYCLE
    NodeID=NodeID+1
    Elem%Node(i)%np%ind=NodeID
  END DO
  DO k=0,NGeo; DO j=0,NGeo; DO i=0,NGeo
    IF(Elem%CurvedNode(i,j,k)%np%ind.NE.-88888) CYCLE
    NodeID=NodeID+1
    Elem%CurvedNode(i,j,k)%np%ind=NodeID
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
CALL WriteArrayToHDF5(File_ID,'ElemInfo',nElems,2,(/nElems,ElemInfoSize/),0,IntegerArray=ElemInfo)

!WRITE ElemWeight,into (1,nElems)  
CALL WriteArrayToHDF5(File_ID,'ElemWeight',nElems,1,(/nElems/),0,RealArray=ElemWeight)

CALL WriteArrayToHDF5(File_ID,'ElemBarycenters',nElems,2,(/nElems,3/),0,RealArray=ElemBarycenters)
DEALLOCATE(ElemBarycenters)

IF(dosortIJK)THEN
  ! WRITE element ijk index (for postprocessing of structured/semistructured domains)
  CALL WriteArrayToHDF5(File_ID,'nElems_IJK',3,1,(/3/),0,IntegerArray=nElems_IJK)
  CALL WriteArrayToHDF5(File_ID,'Elem_IJK',nElems,2,(/nElems,3/),0,IntegerArray=Elem_IJK)
  DEALLOCATE(Elem_IJK)
END IF

!WRITE SideInfo,into (1,nSides)   
CALL WriteArrayToHDF5(File_ID,'SideInfo',nSides,2,(/nSides,SideInfoSize/),0,IntegerArray=SideInfo)
DEALLOCATE(SideInfo)

!WRITE NodeInfo,into (1,nNodes) 
CALL WriteArrayToHDF5(File_ID,'NodeInfo',nNodes,1,(/nNodes/),0,IntegerArray=NodeInfo)
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
CALL WriteAttributeToHDF5(File_ID,'BoundaryOrder',1,IntegerScalar=N+1)
CALL WriteAttributeToHDF5(File_ID,'CurvedFound',1,LogicalScalar=CurvedFound)
nBCs=nUserDefinedBoundaries
ALLOCATE(BCNames(nBCs))
ALLOCATE(BCType(nBCs,4))
DO i=1,nBCs
  BCNames(i)=BoundaryName(i) 
  BCType(i,:)=BoundaryType(i,:) 
END DO
! WRITE BC 

CALL WriteArrayToHDF5(File_ID,'BCNames',nBCs,1,(/nBCs/),0,StrArray=BCNames)
CALL WriteArrayToHDF5(File_ID,'BCType',nBCs,2,(/nBcs,4/),0,IntegerArray=BCType)

DEALLOCATE(BCNames)
DEALLOCATE(BCType)
! Close the file.
CALL CloseHDF5File()
CALL Timer(.FALSE.)

DEALLOCATE(ElemInfo)
DEALLOCATE(ElemWeight)

END SUBROUTINE WriteMeshToHDF5



SUBROUTINE getMeshInfo()
!===================================================================================================================================
! Subroutine prepares ElemInfo,Sideinfo,Nodeinfo,NodeCoords arrays
!===================================================================================================================================
! MODULES
USE MOD_Mesh_Vars,ONLY:tElem,tSide
USE MOD_Mesh_Vars,ONLY:FirstElem
USE MOD_Mesh_Basis,ONLY:ElemGeometry,ISORIENTED
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
INTEGER                        :: iNodeP
!===================================================================================================================================
curvedfound=.FALSE.
!fill ElementInfo. 
ALLOCATE(ElemInfo(1:nElems,ElemInfoSize))
ALLOCATE(ElemWeight(1:nElems))
ElemInfo=0
ElemWeight=0.
Elemcounter=0
Elemcounter(:,1)=(/104,204,105,115,205,106,116,206,108,118,208/)
iNode  = 0 
iSide  = 0
iElem  = 0
Elem=>firstElem
DO WHILE(ASSOCIATED(Elem))
  iElem=iElem+1
  locnNodes=Elem%nNodes+Elem%nCurvedNodes
  locnSides=0
  Side=>Elem%firstSide
  DO WHILE(ASSOCIATED(Side))
    locnNodes = locnNodes+1 ! wirte only first oriented node of side, if curved all!         
    locnSides = locnSides+1
    Side=>Side%nextElemSide
  END DO
  !CALL ElemGeometry(Elem)
  IF(Elem%type.GT.200)curvedfound=.TRUE.
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
  ElemInfo(iElem,ELEM_FirstNodeInd) = iNode            ! first index -1 in NodeInfo
  ElemInfo(iElem,ELEM_LastNodeInd)  = iNode+locnNodes  ! last index in NodeInfo
  iNode = iNode + locnNodes
  iSide = iSide + locnSides
  !approximate weight: locnNodes
  ElemWeight(iElem)=REAL(locnNodes)
  Elem=>Elem%nextElem
END DO


!fill SideInfo
ALLOCATE(SideInfo(1:nSides,SideInfoSize)) 
SideInfo=0 
iSide=0
Elem=>firstElem
DO WHILE(ASSOCIATED(Elem))
  Side=>Elem%firstSide
  DO WHILE(ASSOCIATED(Side))
    iSide=iSide+1
    !Side Tpye
    IF(Side%isCurved)THEN
      IF(ASSOCIATED(Side%curvedNode))THEN
        SideInfo(iSide,SIDE_Type)=Side%nNodes+3         ! Side Type: NL tria/quad, 6/7
      ELSE
        SideInfo(iSide,SIDE_Type)=5                     ! Side Type: bilinear quad
      END IF
    ELSE
      SideInfo(iSide,SIDE_Type)=Side%nNodes             ! Side Type: linear 3/4
    END IF
    !Side ID
    SideInfo(iSide,SIDE_ID)=Side%ind
    IF(.NOT.ISORIENTED(Side)) SideInfo(iSide,SIDE_ID)=-SideInfo(iSide,SIDE_ID)           
    !neighbor Element ID
    IF(ASSOCIATED(Side%Connection))THEN
      SideInfo(iSide,SIDE_nbElemID)=Side%Connection%Elem%ind                   ! Element ID of neighbor Element
    END IF
    !BC ID 
    IF(ASSOCIATED(Side%BC))THEN
      SideInfo(iSide,SIDE_BCID)=Side%BC%BCIndex                            
      IF(Side%BC%BCIndex.EQ.0) WRITE(*,*)'DEBUG, Warning, BC ind =0'
    END IF
    Side=>Side%nextElemSide
  END DO
  Elem=>Elem%nextElem
END DO


!fill NodeInfo
ALLOCATE(NodeInfo(1:nNodes))
NodeInfo=-1
iNode=0
Elem=>firstElem
DO WHILE(ASSOCIATED(Elem))
  DO i=1,Elem%nNodes
    iNode=iNode+1
    NodeInfo(iNode)=Elem%Node(i)%np%ind
  END DO
  DO i=1,Elem%nCurvedNodes
    iNode=iNode+1
    NodeInfo(iNode)=Elem%curvedNode(i)%np%ind
  END DO
  Side=>Elem%firstSide
  DO WHILE(ASSOCIATED(Side))
    iNode=iNode+1
    NodeInfo(iNode)=Side%OrientedNode(1)%np%ind
    Side=>Side%nextElemSide
  END DO
  Elem=>Elem%nextElem
END DO

IF(iNode.NE.nNodes) CALL abort(__STAMP__,&
                               'Sanity check: nNodes not equal to number of nodes in NodeInfo!')

! get mapping from node ID to continuous ID in NodeMap 
ALLOCATE(NodeCoords(MAX(1,nNodeIDs),3))
IF(nNodeIDs.GT.0)THEN
  CALL getNodeMap()  !nodes have to be written
ELSE !no node coords have to be written
  ALLOCATE(NodeMap(1))
  NodeMap=-999
END IF
NodeInfo=ABS(NodeInfo) ! negative sign was used for NodeMap 
NodeCoords=-999.
IF(nNodeIDs.LE.1) RETURN

!fill NodeCoords and Sanity check of curvednodes
Elem=>firstElem
DO WHILE(ASSOCIATED(Elem))
  DO i=1,Elem%nNodes
    Elem%Node(i)%np%tmp=1
  END DO
  ! CURVED  
  DO i=1,Elem%nCurvedNodes
    Elem%CurvedNode(i)%np%tmp=1
  END DO
  Elem=>Elem%nextElem
END DO
iNode=0
Elem=>firstElem
DO WHILE(ASSOCIATED(Elem))
  DO i=1,Elem%nNodes
    IF(Elem%Node(i)%np%tmp.GT.0)THEN
      iNodeP=INVMAP(Elem%Node(i)%np%ind,nNodeIDs,NodeMap)  ! index in local Nodes pointer array
      IF(iNodeP.LE.0) STOP 'Problem in INVMAP' 
      NodeCoords(iNodeP,:)=Elem%Node(i)%np%x
      iNode=iNode+1
      Elem%Node(i)%np%tmp=0
    END IF                  
  END DO
  ! CURVED  
  DO i=1,Elem%nCurvedNodes
    IF(Elem%CurvedNode(i)%np%tmp.GT.0)THEN
      iNodeP=INVMAP(Elem%CurvedNode(i)%np%ind,nNodeIDs,NodeMap)  ! index in local Nodes pointer array
      IF(iNodeP.LE.0) STOP 'Problem in INVMAP' 
      NodeCoords(iNodeP,:)=Elem%CurvedNode(i)%np%x
      iNode=iNode+1
      Elem%CurvedNode(i)%np%tmp=0
    END IF                  
  END DO
  Elem=>Elem%nextElem
END DO
DEALLOCATE(NodeMap)

IF(iNode.NE.nNodeIDs) CALL abort(__STAMP__,&
                      'Sanity check: nNodeIDs not equal to number of nodes in NodeCoords!')

END SUBROUTINE getMeshinfo


SUBROUTINE spaceFillingCurve(IDlist)
!===================================================================================================================================
! Subroutine prepares elementlist and barycentric coodrinates for element sorting by space filling curve and maps back to pointer
! srtructure
!===================================================================================================================================
! MODULES
USE MOD_Mesh_Vars,ONLY:tElemPtr
USE MOD_Mesh_Vars,ONLY:FirstElem
USE MOD_Mesh_Vars,ONLY:MeshMode
USE MOD_Output_Vars,ONLY:DebugVisu,dosortijk
USE MOD_SpaceFillingCurve,ONLY:SortElemsBySpaceFillingCurve
USE MOD_sortIJK,ONLY:SortElemsByCoords
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
INTEGER                        :: IDlist(1:nElems)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES

TYPE(tElemPtr)                 :: Elems(1:nElems)
REAL                           :: ElemBary(1:nElems,3)
INTEGER                        :: ElemID
!===================================================================================================================================

Elems(1)%ep=>firstElem
DO ElemID=2,nElems 
  Elems(ElemID)%ep=>Elems(ElemID-1)%ep%nextElem
END DO
DO ElemID=1,nElems
  IDList(ElemID)=ElemID 
  ElemBary(ElemID,:)=0.
  DO iNode=1,Elems(ElemID)%ep%nNodes
    ElemBary(ElemID,:)=ElemBary(ElemID,:)+Elems(ElemID)%ep%Node(iNode)%np%x
  END DO
  ElemBary(ElemID,:)=ElemBary(ElemID,:)/REAL(Elems(ElemID)%ep%nNodes)
END DO

IF(MeshMode.EQ.11)THEN ! for Meshmode=11: structured single block, elem_IJK already defined
  CALL SortElemsBySpaceFillingCurve(nElems,REAL(Elem_IJK),IDList,1) !use IJK for space filling curve
ELSE
  CALL SortElemsBySpaceFillingCurve(nElems,ElemBary,IDList,2)
END IF

NULLIFY(Elems(IDlist(1))%ep%prevElem)
firstElem=>Elems(IDlist(1))%ep
DO ElemID=2,nElems 
  Elems(IDlist(ElemID-1))%ep%nextElem=>Elems(IDList(ElemID))%ep
  Elems(IDlist(ElemID))%ep%prevElem  =>Elems(IDList(ElemID-1))%ep
END DO
IF(DebugVisu)THEN
  WRITE(*,*)'write space filling curve to sfc.dat'
  OPEN(UNIT=200,FILE='sfc.dat',STATUS='REPLACE')
  DO ElemID=1,nElems 
    WRITE(200,'(3E21.6)')ElemBary(IDlist(ElemID),:)
  END DO
  CLOSE(200)
END IF
IF(ALLOCATED(ElemBarycenters)) DEALLOCATE(ElemBarycenters)
ALLOCATE(ElemBarycenters(1:nElems,3))
DO ElemID=1,nElems 
  ElemBarycenters(ElemID,:)=ElemBary(IDlist(ElemID),:)
END DO
NULLIFY(Elems(IDlist(nElems))%ep%nextElem)
IF(MeshMode.EQ.11)THEN ! for Meshmode=11: structured single block, elem_IJK already defined
  !sort by spacefillingcurve
  Elem_IJK(:,1)=Elem_IJK(IDList(:),1)
  Elem_IJK(:,2)=Elem_IJK(IDList(:),2)
  Elem_IJK(:,3)=Elem_IJK(IDList(:),3)
ELSE
  IF(dosortijk) THEN
    ! do also directly the ijk coordinates of the elements
    ALLOCATE(Elem_IJK(1:nElems,3))
    CALL SortElemsByCoords(nElems,ElemBarycenters(:,:),nElems_IJK,Elem_IJK)
  END IF
END IF
END SUBROUTINE spaceFillingCurve



! HFD5 STUFF
SUBROUTINE WriteArrayToHDF5(Loc_ID,ArrayName,nValglobal,Rank,nVal,offset_in,RealArray,IntegerArray,StrArray)
!===================================================================================================================================
! Subroutine to write Data to HDF5 format
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER                        :: Rank,offset_in
INTEGER(HID_T), INTENT(IN)     :: Loc_ID
INTEGER                        :: nValglobal,nVal(Rank)
CHARACTER(LEN=*), INTENT(IN)   :: ArrayName
REAL              ,DIMENSION(Rank),OPTIONAL,INTENT(IN) :: RealArray
INTEGER           ,DIMENSION(Rank),OPTIONAL,INTENT(IN) :: IntegerArray
CHARACTER(LEN=255),DIMENSION(Rank),OPTIONAL,INTENT(IN) :: StrArray
!TYPE(IO_elem),OPTIONAL,INTENT(IN) :: IO_elemArray(nVal(1))
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                        :: iError
!INTEGER(SIZE_T)                :: tmp1,tmp2,typeoffset,typesize
INTEGER(HID_T)                 :: PList_ID,DSet_ID,MemSpace,FileSpace,HDF5DataType
INTEGER(HSIZE_T)               :: Dimsf(Rank),Offset(Rank)
!===================================================================================================================================
LOGWRITE(UNIT_stdOut,'(A,I1.1,A,A,A)')' WRITE ',Rank,'D ARRAY "',TRIM(ArrayName),'" TO HDF5 FILE...'

! Get global array size, always first dimension!!
Dimsf=nVal
Dimsf(1)=nValGlobal 

! Create the data space for the  dataset.
CALL H5SCREATE_SIMPLE_F(Rank, Dimsf, FileSpace, iError)
! Create the dataset with default properties.
IF(PRESENT(RealArray))     HDF5DataType=H5T_NATIVE_DOUBLE
IF(PRESENT(IntegerArray))  HDF5DataType=H5T_NATIVE_INTEGER
IF(PRESENT(StrArray))THEN
  ! Create HDF5 datatype for the character array.
  CALL H5TCOPY_F(H5T_NATIVE_CHARACTER, HDF5DataType, iError)
  SizeSet=255
  CALL H5TSET_SIZE_F(HDF5DataType, SizeSet, iError)
END IF

CALL H5DCREATE_F(Loc_ID, TRIM(ArrayName), HDF5DataType, filespace, DSet_ID, iError)
CALL H5SCLOSE_F(FileSpace, iError)

! Each process defines dataset in memory and writes it to the hyperslab in the file.
Dimsf=nVal  ! Now we need the local array size
Offset(:)  = 0
Offset(1)  = Offset_in
! Create the data space in the memory
IF(Dimsf(1) .NE. 0)THEN
  CALL H5SCREATE_SIMPLE_F(Rank, Dimsf, MemSpace, iError)
ELSE
  CALL H5SCREATE_F(H5S_NULL_F,MemSpace,iError)
END IF
! Select hyperslab in the file.
CALL H5DGET_SPACE_F(DSet_id, FileSpace, iError)
IF(Dimsf(1) .NE. 0)THEN
  CALL H5SSELECT_HYPERSLAB_F(FileSpace, H5S_SELECT_SET_F, Offset, Dimsf, iError)
ELSE
  CALL H5SSELECT_NONE_F(FileSpace,iError)
END IF

! Create property list for collective dataset write
CALL H5PCREATE_F(H5P_DATASET_XFER_F, PList_ID, iError)
!Write the dataset collectively.

IF(PRESENT(IntegerArray))THEN
  CALL H5DWRITE_F(DSet_ID,HDF5DataType,IntegerArray,Dimsf,iError,file_space_id=filespace,mem_space_id=memspace,xfer_prp=PList_ID)
END IF
IF(PRESENT(RealArray))THEN
  CALL H5DWRITE_F(DSet_ID,HDF5DataType,RealArray   ,Dimsf,iError,file_space_id=filespace,mem_space_id=memspace,xfer_prp=PList_ID)
END IF
IF(PRESENT(StrArray))THEN
  CALL H5DWRITE_F(DSet_ID,HDF5DataType,StrArray    ,Dimsf,iError,file_space_id=filespace,mem_space_id=memspace,xfer_prp=PList_ID)
END IF

! Close the property list.
CALL H5PCLOSE_F(PList_ID, iError)
! Close dataspaces.
CALL H5SCLOSE_F(FileSpace, iError)
CALL H5SCLOSE_F(MemSpace, iError)
! Close the dataset.
CALL H5DCLOSE_F(DSet_ID, iError)

LOGWRITE(UNIT_stdOut,*)'...DONE!'

END SUBROUTINE WriteArrayToHDF5


!SUBROUTINE WriteCoordsToHDF5(Loc_ID,ArrayName,nValglobal,nVal,ElementList,CoordArray)
!!==================================================================================================================================
!! Subroutine to write node coords to HDF5 at position specified by elementlist array
!!==================================================================================================================================
!! MODULES
!! IMPLICIT VARIABLE HANDLING
!IMPLICIT NONE
!!----------------------------------------------------------------------------------------------------------------------------------
!! INPUT VARIABLES
!INTEGER,INTENT(IN)             :: ElementList(:)
!INTEGER(HID_T), INTENT(IN)     :: Loc_ID
!INTEGER                        :: nVal(2),nValglobal
!CHARACTER(LEN=*), INTENT(IN)   :: ArrayName
!REAL,DIMENSION(2),INTENT(IN)   :: CoordArray
!!----------------------------------------------------------------------------------------------------------------------------------
!! OUTPUT VARIABLES
!!----------------------------------------------------------------------------------------------------------------------------------
!! LOCAL VARIABLES
!INTEGER                        :: Rank,i,j
!INTEGER(HID_T)                 :: PList_ID,DSet_ID,MemSpace,FileSpace
!INTEGER(HSIZE_T)               :: Dimsf(2)
!INTEGER(SIZE_T)                :: num_elements
!INTEGER(HSIZE_T)               :: Coords(2,nVal(1)*nVal(2))
!!==================================================================================================================================
!rank=2
!LOGWRITE(UNIT_stdOut,'(A,I1.1,A,A,A)')' WRITE ',Rank,'D ARRAY "',TRIM(ArrayName),'" TO HDF5 FILE...'
!
!! Get global array size
!Dimsf=nVal
!Dimsf(1)=nValGlobal 
!! Create the data space for the  dataset.
!CALL H5SCREATE_SIMPLE_F(Rank, Dimsf, FileSpace, iError)
!! Create the dataset with default properties.
!CALL H5DCREATE_F(Loc_ID, TRIM(ArrayName),H5T_NATIVE_DOUBLE, filespace, DSet_ID, iError)
!CALL H5SCLOSE_F(FileSpace, iError)
!
!! Each process defines dataset in memory and writes it to the hyperslab in the file.
!Dimsf=nVal  ! Now we need the local array size
!
!! Create the data space in the memory
!IF(Dimsf(1) .NE. 0)THEN
!  CALL H5SCREATE_SIMPLE_F(Rank, Dimsf, MemSpace, iError)
!ELSE
!  CALL H5SCREATE_F(H5S_NULL_F,MemSpace,iError)
!END IF
!! Select hyperslab in the file.
!CALL H5DGET_SPACE_F(DSet_id, FileSpace, iError)
!
!IF(Dimsf(1) .NE. 0)THEN
!  !select elements from Element List
!  write(*,*) 'dimsf',dimsf
!  num_elements=0
!  DO j=1,Dimsf(2)
!    DO i=1,Dimsf(1)
!      num_elements=num_elements+1
!      Coords(1,num_elements)=ElementList(i)
!      Coords(2,num_elements)=j 
!      IF(ElementList(i).NE.i) write(*,*) i,ElementList(i)
!    END DO
!  END DO
!  CALL H5SSELECT_ELEMENTS_F(FileSpace, H5S_SELECT_SET_F, rank,num_elements,Coords,iError)
!ELSE
!  CALL H5SSELECT_NONE_F(FileSpace,iError)
!END IF
!
!! Create property list for collective dataset write
!CALL H5PCREATE_F(H5P_DATASET_XFER_F, PList_ID, iError)
!!Write the dataset collectively.
!CALL H5DWRITE_F(DSet_ID,H5T_NATIVE_DOUBLE, &
!                      CoordArray,Dimsf,iError,file_space_id=filespace,mem_space_id=memspace,xfer_prp=PList_ID)
!
!! Close the property list.
!CALL H5PCLOSE_F(PList_ID, iError)
!! Close dataspaces.
!CALL H5SCLOSE_F(FileSpace, iError)
!CALL H5SCLOSE_F(MemSpace, iError)
!! Close the dataset.
!CALL H5DCLOSE_F(DSet_ID, iError)
!
!LOGWRITE(UNIT_stdOut,*)'...DONE!'
!
!END SUBROUTINE WriteCoordsToHDF5



SUBROUTINE WriteAttributeToHDF5(Loc_ID_in,AttribName,nVal,DataSetname,RealScalar,IntegerScalar,StrScalar,LogicalScalar, &
                                                                      RealArray,IntegerArray)
!===================================================================================================================================
! Subroutine to write Attributes to HDF5 format of a given Loc_ID, which can be the File_ID,datasetID,goupID. This must be opened
! outside of the routine. If you directly want to write an attribute to a dataset, just provide the name of the dataset
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER(HID_T),INTENT(IN)              :: Loc_ID_in
CHARACTER(LEN=*), INTENT(IN)           :: AttribName
INTEGER,INTENT(IN)                     :: nVal
CHARACTER(LEN=*),OPTIONAL,INTENT(IN)   :: DatasetName
REAL,OPTIONAL,INTENT(IN)               :: RealArray(nVal)
INTEGER,OPTIONAL,INTENT(IN)            :: IntegerArray(nVal)
REAL,OPTIONAL,INTENT(IN)               :: RealScalar
INTEGER,OPTIONAL,INTENT(IN)            :: IntegerScalar
CHARACTER(LEN=255),OPTIONAL,INTENT(IN) :: StrScalar
LOGICAL,OPTIONAL,INTENT(IN)            :: LogicalScalar
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                        :: Rank
INTEGER(HID_T)                 :: DataSpace, Attr_ID,Loc_ID,aType_ID
INTEGER(HSIZE_T), DIMENSION(1) :: Dimsf
INTEGER(SIZE_T)                :: AttrLen
INTEGER                        :: logtoint
!===================================================================================================================================
LOGWRITE(UNIT_stdOut,*)' WRITE ATTRIBUTE "',TRIM(AttribName),'" TO HDF5 FILE...'
IF(PRESENT(DataSetName))THEN
 ! Open dataset
  CALL H5DOPEN_F(File_ID, TRIM(DatasetName),Loc_ID, iError)
ELSE
  Loc_ID=Loc_ID_in
END IF
! Create scalar data space for the attribute.
Rank=1
Dimsf(:)=0
Dimsf(1)=nVal
CALL H5SCREATE_SIMPLE_F(Rank, Dimsf, DataSpace, iError)
! Create the attribute for group Loc_ID.
! Write the attribute data.
IF(PRESENT(RealArray))THEN
  CALL H5ACREATE_F(Loc_ID, TRIM(AttribName), H5T_NATIVE_DOUBLE, DataSpace, Attr_ID, iError)
  CALL H5AWRITE_F(Attr_ID, H5T_NATIVE_DOUBLE, RealArray , Dimsf, iError)
END IF
IF(PRESENT(RealScalar))THEN
  CALL H5ACREATE_F(Loc_ID, TRIM(AttribName), H5T_NATIVE_DOUBLE, DataSpace, Attr_ID, iError)
  CALL H5AWRITE_F(Attr_ID, H5T_NATIVE_DOUBLE, RealScalar , Dimsf, iError)
END IF
IF(PRESENT(IntegerArray))THEN
  CALL H5ACREATE_F(Loc_ID, TRIM(AttribName), H5T_NATIVE_INTEGER, DataSpace, Attr_ID, iError)
  CALL H5AWRITE_F(Attr_ID, H5T_NATIVE_INTEGER, IntegerArray , Dimsf, iError)
END IF
IF(PRESENT(IntegerScalar))THEN
  CALL H5ACREATE_F(Loc_ID, TRIM(AttribName), H5T_NATIVE_INTEGER, DataSpace, Attr_ID, iError)
  CALL H5AWRITE_F(Attr_ID, H5T_NATIVE_INTEGER, IntegerScalar , Dimsf, iError)
END IF
IF(PRESENT(LogicalScalar))THEN
  IF(logicalScalar)THEN
    logtoint=1
  ELSE
    logtoint=0
  END IF
  CALL H5ACREATE_F(Loc_ID, TRIM(AttribName), H5T_NATIVE_INTEGER, DataSpace, Attr_ID, iError)
  CALL H5AWRITE_F(Attr_ID, H5T_NATIVE_INTEGER, logtoint , Dimsf, iError)
END IF
IF(PRESENT(StrScalar))THEN
  ! Create character string datatype for the attribute.
  ! For a attribute character, we have to build our own type with corresponding attribute length
  CALL H5TCOPY_F(H5T_NATIVE_CHARACTER, atype_id, iError)
  AttrLen=255
  CALL H5TSET_SIZE_F(aType_ID, AttrLen, iError)
  CALL H5ACREATE_F(Loc_ID, TRIM(AttribName), aType_ID, DataSpace, Attr_ID, iError)
  CALL H5AWRITE_F(Attr_ID, aType_ID, StrScalar , Dimsf, iError)
END IF
! Close dataspace
CALL H5SCLOSE_F(DataSpace, iError)
! Close the attribute.
CALL H5ACLOSE_F(Attr_ID, iError)
IF(PRESENT(DataSetName))THEN
  ! Close the dataset and property list.
  CALL H5DCLOSE_F(Loc_ID, iError)
END IF
LOGWRITE(UNIT_stdOut,*)'...DONE!'
END SUBROUTINE WriteAttributeToHDF5

END MODULE MOD_Output_HDF5
