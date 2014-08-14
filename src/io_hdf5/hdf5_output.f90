#include "hopest_f.h"

MODULE MODH_HDF5_output
!===================================================================================================================================
! Add comments please!
!===================================================================================================================================
! MODULES
USE MODH_IO_HDF5
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------

INTERFACE WriteHDF5Header
  MODULE PROCEDURE WriteHDF5Header
END INTERFACE

!INTERFACE WriteArrayToHDF5
!  MODULE PROCEDURE WriteArrayToHDF5
!END INTERFACE

INTERFACE WriteAttributeToHDF5
  MODULE PROCEDURE WriteAttributeToHDF5
END INTERFACE

PUBLIC :: WriteArrayToHDF5
PUBLIC :: WriteHDF5Header
PUBLIC :: WriteAttributeToHDF5
!===================================================================================================================================

CONTAINS


SUBROUTINE WriteHDF5Header(FileType_in,File_ID)
!===================================================================================================================================
! Subroutine to write a distinct file header to each HDF5 file
!===================================================================================================================================
! MODULES
USE MODH_Output_Vars,ONLY:ProgramName,FileVersion,ProjectName
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
CALL WriteAttributeToHDF5(File_ID,'Program'     ,1,StrScalar=ProgramName)
CALL WriteAttributeToHDF5(File_ID,'File_Type'   ,1,StrScalar=FileType_in)
CALL WriteAttributeToHDF5(File_ID,'Project_Name',1,StrScalar=ProjectName)
CALL WriteAttributeToHDF5(File_ID,'File_Version',1,RealScalar=FileVersion)
END SUBROUTINE WriteHDF5Header


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
INTEGER(SIZE_T)                :: SizeSet
!===================================================================================================================================
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

END SUBROUTINE WriteArrayToHDF5




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
END SUBROUTINE WriteAttributeToHDF5


END MODULE MODH_HDF5_output
