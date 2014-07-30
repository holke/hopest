#include "hopest_f.h"

MODULE MOD_HDF5_output
!===================================================================================================================================
! Add comments please!
!===================================================================================================================================
! MODULES
USE MOD_io_HDF5
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------

INTERFACE WriteMeshToHDF5
  MODULE PROCEDURE WriteMeshToHDF5
END INTERFACE

INTERFACE WriteHDF5Header
  MODULE PROCEDURE WriteHDF5Header
END INTERFACE

!INTERFACE WriteArrayToHDF5
!  MODULE PROCEDURE WriteArrayToHDF5
!END INTERFACE

INTERFACE WriteAttribute
  MODULE PROCEDURE WriteAttributeToHDF5
END INTERFACE

PUBLIC :: WriteMeshToHDF5
PUBLIC :: WriteArrayToHDF5
PUBLIC :: WriteHDF5Header
PUBLIC :: WriteAttribute
!===================================================================================================================================

CONTAINS


SUBROUTINE WriteMeshToHDF5(FileString)
!===================================================================================================================================
! Subroutine to write Data to HDF5 format
!===================================================================================================================================
! MODULES
!USE MOD_Mesh_Vars,ONLY:tElem,tSide
!USE MOD_Mesh_Vars,ONLY:FirstElem
!USE MOD_Mesh_Vars,ONLY:readFlowSolution,FlowSol
!USE MOD_Mesh_Vars,ONLY:N
!USE MOD_Output_Vars,ONLY:dosortIJK
!USE MOD_Mesh_Vars,ONLY:nUserDefinedBoundaries,BoundaryName,BoundaryType
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
CHARACTER(LEN=*),INTENT(IN)    :: FileString
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES

!TYPE(tElem),POINTER            :: Elem
!TYPE(tSide),POINTER            :: Side
!INTEGER                        :: ElemID,SideID,NodeID
!INTEGER,ALLOCATABLE            :: IDlist(:)
!REAL                           :: sBasisBary(4)
!CHARACTER(LEN=255)             :: Meshfile
!===================================================================================================================================
!WRITE(UNIT_stdOut,'(132("~"))')
!CALL Timer(.TRUE.)
!WRITE(UNIT_stdOut,'(A)')' WRITE DATA TO HDF5 FILE...'
!! Create the file collectively.
!CALL OpenHDF5File(FileString,create=.TRUE.)  
!
!
!!set all node and side indices =0
!Elem=>firstElem
!DO WHILE(ASSOCIATED(Elem))
!  DO i=1,Elem%nNodes
!    Elem%Node(i)%np%ind=0
!  END DO
!  DO i=1,Elem%nCurvedNodes
!    Elem%curvedNode(i)%np%ind=0
!  END DO
!  Side=>Elem%firstSide
!  DO WHILE(ASSOCIATED(Side))
!    Side%ind=0
!    !CURVED
!    DO i=1,side%nCurvedNodes
!      side%curvedNode(i)%np%ind=0
!    END DO
!    Side=>Side%nextElemSide
!  END DO
!  Elem=>Elem%nextElem
!END DO
!
!! count Elements , unique sides and nodes are marked with ind=0
!nNodeIDs=0 !number of unique nodeIDs
!nSideIDs=0 !number of unique side IDs (side and side%connection have the same sideID)
!nElems=0   !number of elements
!nSides=0   !number of all sides
!nNodes=0   !number of all nodes
!
!Elem=>firstElem
!DO WHILE(ASSOCIATED(Elem))
!  nElems=nElems+1
!  locnNodes=Elem%nNodes+Elem%nCurvedNodes
!  locnSides=0
!  Side=>Elem%firstSide
!
!  ! Count nodes
!  DO i=1,Elem%nNodes
!    IF(Elem%Node(i)%np%ind.NE.0) CYCLE
!    nNodeIDs=nNodeIDs+1
!    Elem%Node(i)%np%ind=-88888  ! mark no MPI side
!  END DO
!  DO i=1,Elem%nCurvedNodes
!    IF(Elem%CurvedNode(i)%np%ind.NE.0) CYCLE
!    nNodeIDs=nNodeIDs+1
!    Elem%CurvedNode(i)%np%ind=-88888
!  END DO
!
!  ! Count sides
!  DO WHILE(ASSOCIATED(Side))
!    locnNodes = locnNodes + 1 ! write only first oriented node of side, if curved all!  
!    locnSides = locnSides + 1
!    IF(side%ind.EQ.0) THEN
!      nSideIDs=nSideIDs+1
!      Side%ind=-88888
!      IF(ASSOCIATED(Side%connection))THEN      
!        IF(Side%connection%ind.EQ.0) nSideIDs=nSideIDs-1 ! count inner and periodic sides only once 
!      END IF
!    END IF
!    Side=>Side%nextElemSide
!  END DO
!  nNodes = nNodes+locnNodes
!  nSides = nSides+locnSides
!  Elem=>Elem%nextElem
!END DO
!
!! prepare sorting by space filling curve
!ALLOCATE(IDlist(1:nElems))
!CALL SpaceFillingCurve(IDList)
!
!!set unique nodes and Side Indices
!ElemID=0
!SideID=0
!NodeID=0
!Elem=>firstElem
!DO WHILE(ASSOCIATED(Elem))
!  ElemID=ElemID+1 
!  Elem%ind=ElemID
!  DO i=1,Elem%nNodes
!    IF(Elem%Node(i)%np%ind.NE.-88888) CYCLE
!    NodeID=NodeID+1
!    Elem%Node(i)%np%ind=NodeID
!  END DO
!  DO i=1,Elem%nCurvedNodes
!    IF(Elem%CurvedNode(i)%np%ind.NE.-88888) CYCLE
!    NodeID=NodeID+1
!    Elem%CurvedNode(i)%np%ind=NodeID
!  END DO
!
!  Side=>Elem%firstSide
!  DO WHILE(ASSOCIATED(Side))
!    IF(side%ind.EQ.-88888) THEN  ! assign side ID only for non MPI sides and lower MPI sides
!      SideID=SideID+1
!      Side%ind=SideID
!      IF(ASSOCIATED(Side%connection))THEN     
!        IF(side%connection%ind.NE.-88888) Side%connection%ind=SideID !already assigned
!      END IF
!    END IF
!    Side=>Side%nextElemSide
!  END DO
!  Elem=>Elem%nextElem
!END DO !Elem
!
!CALL getMeshInfo() !allocates and fills ElemInfo,SideInfo,NodeInfo,NodeCoords
!
!!WRITE ElemInfo,into (1,nElems)  
!CALL WriteArrayToHDF5(File_ID,'ElemInfo',nElems,2,(/nElems,ElemInfoSize/),0,IntegerArray=ElemInfo)
!
!!WRITE ElemWeight,into (1,nElems)  
!CALL WriteArrayToHDF5(File_ID,'ElemWeight',nElems,1,(/nElems/),0,RealArray=ElemWeight)
!
!CALL WriteArrayToHDF5(File_ID,'ElemBarycenters',nElems,2,(/nElems,3/),0,RealArray=ElemBarycenters)
!DEALLOCATE(ElemBarycenters)
!
!IF(dosortIJK)THEN
!  ! WRITE element ijk index (for postprocessing of structured/semistructured domains)
!  CALL WriteArrayToHDF5(File_ID,'nElems_IJK',3,1,(/3/),0,IntegerArray=nElems_IJK)
!  CALL WriteArrayToHDF5(File_ID,'Elem_IJK',nElems,2,(/nElems,3/),0,IntegerArray=Elem_IJK)
!  DEALLOCATE(Elem_IJK)
!END IF
!
!!WRITE SideInfo,into (1,nSides)   
!CALL WriteArrayToHDF5(File_ID,'SideInfo',nSides,2,(/nSides,SideInfoSize/),0,IntegerArray=SideInfo)
!DEALLOCATE(SideInfo)
!
!!WRITE NodeInfo,into (1,nNodes) 
!CALL WriteArrayToHDF5(File_ID,'NodeInfo',nNodes,1,(/nNodes/),0,IntegerArray=NodeInfo)
!DEALLOCATE(NodeInfo)
!
!! WRITE NodeCoords (have to be sorted according to nodemap)
!CALL WriteArrayToHDF5(File_ID,'NodeCoords',nNodeIDs,2,(/nNodeIDs,3/),0,RealArray=NodeCoords)
!DEALLOCATE(NodeCoords)
!
!!! WRITE NodeCoords,arbitrary ordering by NodeMap
!!CALL WriteCoordsToHDF5(File_ID,'NodeCoords',nNodeIDs,(/nNodeIDs,3/),NodeMap,NodeCoords)
!
!
!CALL WriteArrayToHDF5(File_ID,'ElemCounter',11,2,(/11,2/),0,IntegerArray=ElemCounter)
!WRITE(*,*)'Mesh statistics:'
!WRITE(*,*)'Element Type | number of elements'
!DO i=1,11
!  WRITE(*,'(I4,A,I8)') Elemcounter(i,1),'        | ',Elemcounter(i,2)
!END DO
!
!!attributes 
!CALL WriteAttributeToHDF5(File_ID,'BoundaryOrder',1,IntegerScalar=N+1)
!CALL WriteAttributeToHDF5(File_ID,'CurvedFound',1,LogicalScalar=CurvedFound)
!nBCs=nUserDefinedBoundaries
!ALLOCATE(BCNames(nBCs))
!ALLOCATE(BCType(nBCs,4))
!DO i=1,nBCs
!  BCNames(i)=BoundaryName(i) 
!  BCType(i,:)=BoundaryType(i,:) 
!END DO
!! WRITE BC 
!
!CALL WriteArrayToHDF5(File_ID,'BCNames',nBCs,1,(/nBCs/),0,StrArray=BCNames)
!CALL WriteArrayToHDF5(File_ID,'BCType',nBCs,2,(/nBcs,4/),0,IntegerArray=BCType)
!
!DEALLOCATE(BCNames)
!DEALLOCATE(BCType)
!! Close the file.
!CALL CloseHDF5File()
!CALL Timer(.FALSE.)
!
!IF(readFlowSolution)THEN !write flow solution, which was read in from cgns file
!  WRITE(UNIT_stdOut,'(132("~"))')
!  CALL Timer(.TRUE.)
!  WRITE(UNIT_stdOut,'(A)')' WRITE RESTART FILE TO HDF5 FILE...'
!  WRITE(*,*)'DEBUG, nElems= ',nElems 
!  ! Create the file collectively.
!  CALL OpenHDF5File('CGNSrestartfile.h5',create=.TRUE.)  
!  CALL WriteAttributeToHDF5(File_ID,'Time',1,RealScalar=0.0)
!  MeshFile=TRIM(FileString)
!  CALL WriteAttributeToHDF5(File_ID,'MeshFile',1,StrScalar=MeshFile)
!  ElemWeight=1.
!  CALL WriteArrayToHDF5(File_ID,'ElemWeight',nElems,1,(/nElems/),0,RealArray=ElemWeight)
!  ! resort Flow Solution
!  FlowSol(1:nElems,1:5)=FlowSol(IDlist(1:nElems),1:5)
!  ! Halo style mean value: meanvalue/sqrt(refelem%volume)
!  sBasisBary(1)=SQRT(1./6.)
!  sBasisBary(2)=SQRT(1./3.)
!  sBasisBary(3)=SQRT(0.5)
!  sBasisBary(4)=1.
!  iElem=0
!  Elem=>firstElem
!  DO WHILE(ASSOCIATED(Elem))
!    iElem=iElem+1
!    SELECT CASE(Elem%type)
!    CASE(104)
!      FlowSol(iElem,1:5)=FlowSol(iElem,1:5)*sBasisBary(1)
!    CASE(105)
!      FlowSol(iElem,1:5)=FlowSol(iElem,1:5)*sBasisBary(2)
!    CASE(106)
!      FlowSol(iElem,1:5)=FlowSol(iElem,1:5)*sBasisBary(3) 
!    CASE(115)
!      FlowSol(iElem,1:5)=FlowSol(iElem,1:5)*sBasisBary(2) !approximation
!    CASE(116)
!      FlowSol(iElem,1:5)=FlowSol(iElem,1:5)*sBasisBary(3) !approximation
!    CASE(118)
!      FlowSol(iElem,1:5)=FlowSol(iElem,1:5)*sBasisBary(4) !approximation
!  !  CASE (115,116,118)
!  !    WRITE(*,*)'Elem%type',Elem%type
!  !    STOP 'Halo style mean value not jet implemented for trilinear hexa/penta/pyra'
!    END SELECT
!    Elem=>Elem%nextElem
!  END DO
!  CALL WriteArrayToHDF5(File_ID,'DGsolution',nElems,2,(/nElems,5/),0,RealArray=FlowSol(1:nElems,1:5))
!  !use elemInfo as dummy for elemint
!  ElemInfo(:,1)=1 !order
!  DO iElem=1,nElems
!    ElemInfo(iElem,2)=iElem-1
!    ElemInfo(iElem,3)=iElem
!  END DO
!  CALL WriteArrayToHDF5(File_ID,'ElemInt',nElems,2,(/nElems,3/),0,IntegerArray=ElemInfo(1:nElems,1:3))
!  CALL CloseHDF5File()
!  CALL Timer(.FALSE.)
!END IF !readFlowSolution
!DEALLOCATE(ElemInfo)
!DEALLOCATE(ElemWeight)
!DEALLOCATE(IDlist)

END SUBROUTINE WriteMeshToHDF5


SUBROUTINE GenerateFileSkeleton(TypeString,nVar,StrVarNames,MeshFileName)
!===================================================================================================================================
! Subroutine that generates the output file on a single processor and writes all the necessary attributes (better MPI performance)
!===================================================================================================================================
! MODULES
USE MOD_PreProc
USE MOD_Globals
USE MOD_Mesh_Vars  ,ONLY: nGlobalElems
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
CHARACTER(LEN=*),INTENT(IN)    :: TypeString
INTEGER                        :: nVar
CHARACTER(LEN=255)             :: StrVarNames(nVar)
CHARACTER(LEN=*),INTENT(IN)    :: MeshFileName
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!INTEGER(HID_T)                 :: DSet_ID,FileSpace,HDF5DataType
!INTEGER(HSIZE_T)               :: Dimsf(5)
!CHARACTER(LEN=255)             :: FileName,MeshFile255
!===================================================================================================================================
!! Create file
!FileName=TRIM(TIMESTAMP(TRIM(ProjectName)//'_'//TRIM(TypeString),OutputTime))//'.h5'
!CALL OpenDataFile(TRIM(FileName),create=.TRUE.,single=.TRUE.)
!
!! Write file header
!CALL WriteHDF5Header(TRIM(TypeString),File_ID)
!
!! Preallocate the data space for the dataset.
!Dimsf=(/nVar,NOut+1,NOut+1,NOut+1,nGlobalElems/)
!
!CALL H5SCREATE_SIMPLE_F(5, Dimsf, FileSpace, iError)
!! Create the dataset with default properties.
!HDF5DataType=H5T_NATIVE_DOUBLE
!CALL H5DCREATE_F(File_ID,'DG_Solution', HDF5DataType, FileSpace, DSet_ID, iError)
!! Close the filespace and the dataset
!CALL H5DCLOSE_F(Dset_id, iError)
!CALL H5SCLOSE_F(FileSpace, iError)
!
!! Write dataset properties "Time","MeshFile","NextFile","NodeType","VarNames"
!CALL WriteAttributeToHDF5(File_ID,'N',1,IntegerScalar=N)
!CALL WriteAttributeToHDF5(File_ID,'Time',1,RealScalar=OutputTime)
!CALL WriteAttributeToHDF5(File_ID,'MeshFile',1,StrScalar=TRIM(MeshFileName))
!IF(PRESENT(FutureTime))THEN
!  MeshFile255=TRIM(TIMESTAMP(TRIM(ProjectName)//'_'//TRIM(TypeString),FutureTime))//'.h5'
!  CALL WriteAttributeToHDF5(File_ID,'NextFile',1,StrScalar=MeshFile255)
!END IF
!CALL WriteAttributeToHDF5(File_ID,'NodeType',1,StrScalar=NodeType)
!CALL WriteAttributeToHDF5(File_ID,'VarNames',nVar,StrArray=StrVarNames)
!
!IF(NOut.NE.PP_N) CALL WriteAttributeToHDF5(File_ID,'NComputation',1,IntegerScalar=PP_N)
!
!CALL CloseDataFile()
END SUBROUTINE GenerateFileSkeleton


SUBROUTINE WriteHDF5Header(FileType_in,File_ID)
!===================================================================================================================================
! Subroutine to write a distinct file header to each HDF5 file
!===================================================================================================================================
! MODULES
USE MOD_Output_Vars,ONLY:ProgramName,FileVersion,ProjectName
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
CHARACTER(LEN=*),INTENT(IN)              :: FileType_in
INTEGER(HID_T),INTENT(IN)                :: File_ID
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
! Write a small file header to identify a Flexi HDF5 files
! Attributes are program name, file type identifier, project name and version number
CALL WriteAttributeToHDF5(File_ID,'Program'     ,1,StrScalar=TRIM(ProgramName))
CALL WriteAttributeToHDF5(File_ID,'File_Type'   ,1,StrScalar=TRIM(FileType_in))
CALL WriteAttributeToHDF5(File_ID,'Project_Name',1,StrScalar=TRIM(ProjectName))
CALL WriteAttributeToHDF5(File_ID,'File_Version',1,RealScalar=FileVersion)
END SUBROUTINE WriteHDF5Header


SUBROUTINE WriteArrayToHDF5(DataSetName,rank,nValGlobal,nVal,offset,&
                            existing,collective,resizeDim,chunkSize,&
                            RealArray,IntegerArray,StrArray)
!===================================================================================================================================
! Subroutine to write Data to HDF5 format
!===================================================================================================================================
! MODULES
USE MOD_Globals
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
CHARACTER(LEN=*),INTENT(IN)   :: DataSetName
INTEGER,INTENT(IN)            :: rank             ! number of dimensions of the array
INTEGER,INTENT(IN)            :: nValGlobal(rank) ! max size of array in offset dimension
INTEGER,INTENT(IN)            :: nVal(rank)       ! size of complete (local) array to write
INTEGER,INTENT(IN)            :: offset(rank)     ! offset =0, start at beginning of the array
LOGICAL,INTENT(IN)            :: existing
LOGICAL,INTENT(IN)            :: collective       ! use collective writes from all procs
LOGICAL,INTENT(IN),OPTIONAL   :: resizeDim(rank)  ! specify dimensions which can be resized (enlarged)
INTEGER,INTENT(IN),OPTIONAL   :: chunkSize(rank)  ! specify chunksize
REAL   ,INTENT(IN),OPTIONAL            :: RealArray(rank)
INTEGER,INTENT(IN),OPTIONAL            :: IntegerArray(rank)
CHARACTER(LEN=255),INTENT(IN),OPTIONAL :: StrArray(rank)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER(HID_T)                 :: PList_ID,DSet_ID,MemSpace,FileSpace,HDF5DataType,dsetparams
INTEGER(HSIZE_T)               :: Dimsf(Rank),OffsetHDF(Rank),nValMax(Rank)
INTEGER(SIZE_T)                :: SizeSet
LOGICAL                        :: chunky
!===================================================================================================================================
LOGWRITE(*,'(A,I1.1,A,A,A)')' WRITE ',Rank,'D ARRAY "',TRIM(DataSetName),'" TO HDF5 FILE...'

! specify chunk size if desired 
nValMax=nValGlobal
chunky=.FALSE.
CALL H5PCREATE_F(H5P_DATASET_CREATE_F,dsetparams,iError)
IF(PRESENT(chunkSize))THEN
  chunky=.TRUE.
  Dimsf=chunkSize
  CALL H5PSET_CHUNK_F(dsetparams,rank,dimsf,iError)
END IF
! make array extendable in case you want to append something
IF(PRESENT(resizeDim))THEN
  IF(.NOT.PRESENT(chunkSize))&
    CALL abort(__STAMP__,&
    'Chunk size has to be specified when using resizable arrays.')
  nValMax = MERGE(H5S_UNLIMITED_F,nValMax,resizeDim)
END IF

! Create the dataset with default properties.
IF(PRESENT(RealArray))     HDF5DataType=H5T_NATIVE_DOUBLE
IF(PRESENT(IntegerArray))  HDF5DataType=H5T_NATIVE_INTEGER
IF(PRESENT(StrArray))THEN
  ! Create HDF5 datatype for the character array.
  CALL H5TCOPY_F(H5T_NATIVE_CHARACTER, HDF5DataType, iError)
  SizeSet=255
  CALL H5TSET_SIZE_F(HDF5DataType, SizeSet, iError)
END IF

Dimsf = nValGlobal ! we need the global array size
IF(existing)THEN
  CALL H5DOPEN_F(File_ID, TRIM(DatasetName),DSet_ID, iError)
ELSE
  ! Create the data space for the  dataset.
  CALL H5SCREATE_SIMPLE_F(Rank, Dimsf, FileSpace, iError, nValMax)
  CALL H5DCREATE_F(File_ID, TRIM(DataSetName), HDF5DataType, FileSpace, DSet_ID,iError,dsetparams)
  CALL H5SCLOSE_F(FileSpace, iError)
END IF
IF(chunky)THEN
  CALL H5DSET_EXTENT_F(DSet_ID,Dimsf,iError) ! if resizable then dataset may need to be extended
END IF

! Each process defines dataset in memory and writes it to the hyperslab in the file.
Dimsf=nVal  ! Now we need the local array size
OffsetHDF = Offset
! Create the data space in the memory
IF(ANY(Dimsf.EQ.0))THEN
  CALL H5SCREATE_F(H5S_NULL_F,MemSpace,iError)
ELSE
  CALL H5SCREATE_SIMPLE_F(Rank, Dimsf, MemSpace, iError)
END IF
! Select hyperslab in the file.
CALL H5DGET_SPACE_F(DSet_id, FileSpace, iError)
IF(ANY(Dimsf.EQ.0))THEN
  CALL H5SSELECT_NONE_F(FileSpace,iError)
ELSE
  CALL H5SSELECT_HYPERSLAB_F(FileSpace, H5S_SELECT_SET_F, OffsetHDF, Dimsf, iError)
END IF

! Create property list for collective dataset write
CALL H5PCREATE_F(H5P_DATASET_XFER_F, PList_ID, iError)
#ifdef MPI
IF(collective)THEN
  CALL H5PSET_DXPL_MPIO_F(PList_ID, H5FD_MPIO_COLLECTIVE_F,  iError)
ELSE
  CALL H5PSET_DXPL_MPIO_F(PList_ID, H5FD_MPIO_INDEPENDENT_F, iError)
END IF
#endif

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

! Close the property list, dataspaces and dataset.
CALL H5PCLOSE_F(dsetparams, iError)
CALL H5PCLOSE_F(PList_ID, iError)
CALL H5SCLOSE_F(FileSpace, iError)
CALL H5SCLOSE_F(MemSpace, iError)
CALL H5DCLOSE_F(DSet_ID, iError)

LOGWRITE(*,*)'...DONE!'
END SUBROUTINE WriteArrayToHDF5


SUBROUTINE WriteAttributeToHDF5(Loc_ID_in,AttribName,nVal,DataSetname,RealScalar,IntegerScalar,StrScalar,LogicalScalar, &
                                                                      RealArray,IntegerArray,StrArray)
!===================================================================================================================================
! Subroutine to write Attributes to HDF5 format of a given Loc_ID, which can be the File_ID,datasetID,groupID. This must be opened
! outside of the routine. If you directly want to write an attribute to a dataset, just provide the name of the dataset
!===================================================================================================================================
! MODULES
USE MOD_Globals
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
CHARACTER(LEN=255),OPTIONAL,INTENT(IN) :: StrArray(nVal)
REAL,OPTIONAL,INTENT(IN)               :: RealScalar
INTEGER,OPTIONAL,INTENT(IN)            :: IntegerScalar
CHARACTER(LEN=*),OPTIONAL,INTENT(IN)   :: StrScalar
LOGICAL,OPTIONAL,INTENT(IN)            :: LogicalScalar
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                        :: Rank
INTEGER(HID_T)                 :: DataSpace,Attr_ID,Loc_ID,aType_ID
INTEGER(HSIZE_T), DIMENSION(1) :: Dimsf
INTEGER(SIZE_T)                :: AttrLen
INTEGER                        :: logtoint
!===================================================================================================================================
LOGWRITE(*,*)' WRITE ATTRIBUTE "',TRIM(AttribName),'" TO HDF5 FILE...'
IF(PRESENT(DataSetName))THEN
  ! Open dataset
  IF(TRIM(DataSetName).NE.'') CALL H5DOPEN_F(File_ID, TRIM(DatasetName),Loc_ID, iError)
ELSE
  Loc_ID=Loc_ID_in
END IF
! Create scalar data space for the attribute.
Rank=1
Dimsf(:)=0 !???
Dimsf(1)=nVal
CALL H5SCREATE_SIMPLE_F(Rank, Dimsf, DataSpace, iError)
! Create the attribute for group Loc_ID.
! Write the attribute data.
IF(PRESENT(RealArray))THEN
  CALL H5ACREATE_F(Loc_ID, TRIM(AttribName), H5T_NATIVE_DOUBLE, DataSpace, Attr_ID, iError)
  CALL H5AWRITE_F(Attr_ID, H5T_NATIVE_DOUBLE, RealArray, Dimsf, iError)
END IF
IF(PRESENT(RealScalar))THEN
  CALL H5ACREATE_F(Loc_ID, TRIM(AttribName), H5T_NATIVE_DOUBLE, DataSpace, Attr_ID, iError)
  CALL H5AWRITE_F(Attr_ID, H5T_NATIVE_DOUBLE, RealScalar, Dimsf, iError)
END IF
IF(PRESENT(IntegerArray))THEN
  CALL H5ACREATE_F(Loc_ID, TRIM(AttribName), H5T_NATIVE_INTEGER, DataSpace, Attr_ID, iError)
  CALL H5AWRITE_F(Attr_ID, H5T_NATIVE_INTEGER, IntegerArray, Dimsf, iError)
END IF
IF(PRESENT(IntegerScalar))THEN
  CALL H5ACREATE_F(Loc_ID, TRIM(AttribName), H5T_NATIVE_INTEGER, DataSpace, Attr_ID, iError)
  CALL H5AWRITE_F(Attr_ID, H5T_NATIVE_INTEGER, IntegerScalar, Dimsf, iError)
END IF
IF(PRESENT(LogicalScalar))THEN
  IF(logicalScalar)THEN
    logtoint=1
  ELSE
    logtoint=0
  END IF
  CALL H5ACREATE_F(Loc_ID, TRIM(AttribName), H5T_NATIVE_INTEGER, DataSpace, Attr_ID, iError)
  CALL H5AWRITE_F(Attr_ID, H5T_NATIVE_INTEGER, logtoint, Dimsf, iError)
END IF
IF(PRESENT(StrScalar))THEN
  ! Create character string datatype for the attribute.
  ! For a attribute character, we have to build our own type with corresponding attribute length
  CALL H5TCOPY_F(H5T_NATIVE_CHARACTER, atype_id, iError)
  AttrLen=LEN(StrScalar)
  CALL H5TSET_SIZE_F(aType_ID, AttrLen, iError)
  CALL H5ACREATE_F(Loc_ID, TRIM(AttribName), aType_ID, DataSpace, Attr_ID, iError)
  CALL H5AWRITE_F(Attr_ID, aType_ID, StrScalar, Dimsf, iError)
END IF
IF(PRESENT(StrArray))THEN
  ! Create character string array datatype for the attribute.
  ! For a attribute character, we have to build our own type with corresponding attribute length
  CALL H5TCOPY_F(H5T_NATIVE_CHARACTER, atype_id, iError)
  AttrLen=255
  CALL H5TSET_SIZE_F(aType_ID, AttrLen, iError)
  CALL H5ACREATE_F(Loc_ID, TRIM(AttribName), aType_ID, DataSpace, Attr_ID, iError)
  CALL H5AWRITE_F(Attr_ID, aType_ID, StrArray, Dimsf, iError)
END IF
! Close dataspace
CALL H5SCLOSE_F(DataSpace, iError)
! Close the attribute.
CALL H5ACLOSE_F(Attr_ID, iError)
IF(Loc_ID.NE.Loc_ID_in)THEN
  ! Close the dataset and property list.
  CALL H5DCLOSE_F(Loc_ID, iError)
END IF
LOGWRITE(*,*)'...DONE!'
END SUBROUTINE WriteAttributeToHDF5

END MODULE MOD_HDF5_output
