#include "hopest_f.h"

MODULE MODH_HDF5_Input
!===================================================================================================================================
! Add comments please!
!===================================================================================================================================
! MODULES
USE MODH_io_hdf5
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------

INTERFACE ISVALIDHDF5FILE
  MODULE PROCEDURE ISVALIDHDF5FILE
END INTERFACE

INTERFACE GetDataSize
  MODULE PROCEDURE GetHDF5DataSize
END INTERFACE

!INTERFACE ReadArray
!  MODULE PROCEDURE ReadArrayFromHDF5
!END INTERFACE

INTERFACE ReadAttribute
  MODULE PROCEDURE ReadAttributeFromHDF5
END INTERFACE

PUBLIC :: ISVALIDHDF5FILE,GetDataSize
PUBLIC :: ReadArray,ReadAttribute
PUBLIC :: File_ID,HSize,nDims        ! Variables that need to be public
PUBLIC :: OpenDataFile,CloseDataFile ! Subroutines that need to be public
!===================================================================================================================================

CONTAINS

FUNCTION ISVALIDHDF5FILE(FileName)
!===================================================================================================================================
! Subroutine to check if a file is a valid Flexi HDF5 file
!===================================================================================================================================
! MODULES
USE MODH_Globals
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
CHARACTER(LEN=*),INTENT(IN)    :: FileName
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
LOGICAL                        :: isValidHDF5File
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                           :: FileVersion
INTEGER(HID_T)                 :: Plist_ID
CHARACTER(LEN=255)             :: ProgramName
!===================================================================================================================================
isValidHDF5File=.TRUE.
iError=0

! Disable error messages
CALL H5ESET_AUTO_F(0, iError)
! Initialize FORTRAN predefined datatypes
CALL H5OPEN_F(iError)
IF(iError.NE.0)THEN
  CALL Abort(__STAMP__,&
       'ERROR: COULD NOT OPEN FILE!')
END IF

! Open HDF5 file
CALL H5FOPEN_F(TRIM(FileName), H5F_ACC_RDONLY_F, File_ID, iError,access_prp = Plist_ID)
CALL H5PCLOSE_F(Plist_ID, iError)
IF(iError.EQ.0) THEN
  isValidHDF5File=.TRUE.
  ! Check program name -------------------------------------------------------------------------------------------------------------
  ! Open the attribute "Program" of root group
  CALL ReadAttributeFromHDF5(File_ID,'Program',1,StrScalar=ProgramName)
  IF(TRIM(ProgramName) .NE. 'Flexi') isValidHDF5File=.FALSE.
 
  ! Check file version -------------------------------------------------------------------------------------------------------------
  ! Open the attribute "File_Version" of root group
  CALL ReadAttributeFromHDF5(File_ID,'File_Version',1,RealScalar=FileVersion)
  IF(FileVersion .LT. 0.1)THEN
    isValidHDF5File=.FALSE.
    SWRITE(UNIT_stdOut,'(A)')' ERROR: FILE VERSION < 0.1, FILE TOO OLD! '
    SWRITE(UNIT_stdOut,'(A)')'        Try performing a restart...'
  END IF
  ! Close the file.
  CALL H5FCLOSE_F(File_ID, iError)
  ! Close FORTRAN predefined datatypes
  CALL H5CLOSE_F(iError)
ELSE
  ! Close FORTRAN predefined datatypes
  isValidHDF5File=.FALSE.
  CALL H5CLOSE_F(iError)
END IF
END FUNCTION ISVALIDHDF5FILE



SUBROUTINE GetHDF5DataSize(Loc_ID,DSetName,nDims,Size)
!===================================================================================================================================
! Subroutine to determine HDF5 datasize
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
CHARACTER(LEN=*)                     :: DSetName
INTEGER(HID_T),INTENT(IN)            :: Loc_ID
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
INTEGER,INTENT(OUT)                  :: nDims
INTEGER(HSIZE_T),POINTER,INTENT(OUT) :: Size(:)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER(HID_T)                       :: DSet_ID,FileSpace
INTEGER(HSIZE_T), POINTER            :: SizeMax(:)
!===================================================================================================================================
! Open the dataset with default properties.
CALL H5DOPEN_F(Loc_ID, TRIM(DSetName) , DSet_ID, iError)
! Get the data space of the dataset.
CALL H5DGET_SPACE_F(DSet_ID, FileSpace, iError)
! Get number of dimensions of data space
CALL H5SGET_SIMPLE_EXTENT_NDIMS_F(FileSpace, nDims, iError)
! Get size and max size of data space
ALLOCATE(Size(nDims),SizeMax(nDims))
CALL H5SGET_SIMPLE_EXTENT_DIMS_F(FileSpace, Size, SizeMax, iError)
CALL H5SCLOSE_F(FileSpace, iError)
CALL H5DCLOSE_F(DSet_ID, iError)

END SUBROUTINE GetHDF5DataSize


SUBROUTINE ReadArray(ArrayName,Rank,nVal,Offset_in,Offset_dim,RealArray,IntegerArray,StrArray)
!===================================================================================================================================
! Subroutine to read arrays of rank "Rank" with dimensions "Dimsf(1:Rank)".
!===================================================================================================================================
! MODULES
USE MODH_Globals
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER                        :: Rank                  ! number of dimensions of the array
INTEGER                        :: offset_in             ! offset =0, start at beginning of the array
INTEGER                        :: offset_dim            ! which dimension is the offset (only one dimension possible here)
INTEGER                        :: nVal(Rank)            ! size of complete (local) array to write
CHARACTER(LEN=*),INTENT(IN)    :: ArrayName
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL              ,DIMENSION(PRODUCT(nVal)),OPTIONAL,INTENT(OUT) :: RealArray
INTEGER           ,DIMENSION(PRODUCT(nVal)),OPTIONAL,INTENT(OUT) :: IntegerArray
CHARACTER(LEN=255),DIMENSION(PRODUCT(nVal)),OPTIONAL,INTENT(OUT) :: StrArray
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER(HID_T)                                                   :: DSet_ID,Type_ID,MemSpace,FileSpace,PList_ID
INTEGER(HSIZE_T)                                                 :: Offset(Rank),Dimsf(Rank)
!===================================================================================================================================
LOGWRITE(*,'(A,I1.1,A,A,A)')'    READ ',Rank,'D ARRAY "',TRIM(ArrayName),'"'
Dimsf=nVal
LOGWRITE(*,*)'Dimsf,Offset=',Dimsf,Offset_in
CALL H5SCREATE_SIMPLE_F(Rank, Dimsf, MemSpace, iError)
CALL H5DOPEN_F(File_ID, TRIM(ArrayName) , DSet_ID, iError)
! Define and select the hyperslab to use for reading.
CALL H5DGET_SPACE_F(DSet_ID, FileSpace, iError)
Offset(:)=0
Offset(offset_dim)=Offset_in
CALL H5SSELECT_HYPERSLAB_F(FileSpace, H5S_SELECT_SET_F, Offset, Dimsf, iError)
! Create property list
CALL H5PCREATE_F(H5P_DATASET_XFER_F, PList_ID, iError)
#ifdef MPI
! Set property list to collective dataset read
CALL H5PSET_DXPL_MPIO_F(PList_ID, H5FD_MPIO_COLLECTIVE_F, iError)
#endif
! Read the data
IF(PRESENT(RealArray))THEN
  CALL H5DREAD_F(DSet_ID,H5T_NATIVE_DOUBLE,&
                     RealArray,Dimsf,iError,mem_space_id=MemSpace,file_space_id=FileSpace,xfer_prp=PList_ID)
END IF
IF(PRESENT(IntegerArray))THEN
  CALL H5DREAD_F(DSet_ID,H5T_NATIVE_INTEGER,&
                  IntegerArray,Dimsf,iError,mem_space_id=MemSpace,file_space_id=FileSpace,xfer_prp=PList_ID)
END IF
IF(PRESENT(StrArray))THEN
  ! Get datatype for the character string array
  CALL H5DGET_TYPE_F(DSet_ID, Type_ID, iError)
  CALL H5DREAD_F(DSet_ID,Type_ID,&
                      StrArray,Dimsf,iError,mem_space_id=MemSpace,file_space_id=FileSpace,xfer_prp=PList_ID)
  CALL H5TCLOSE_F(Type_ID, iError)
END IF

! Close the property list
CALL H5PCLOSE_F(PList_ID,iError)
! Close the file dataspace
CALL H5SCLOSE_F(FileSpace,iError)
! Close the dataset
CALL H5DCLOSE_F(DSet_ID, iError)
! Close the memory dataspace
CALL H5SCLOSE_F(MemSpace,iError)

LOGWRITE(*,*)'...DONE!'
END SUBROUTINE ReadArray



SUBROUTINE ReadAttributeFromHDF5(Loc_ID_in,AttribName,nVal,DatasetName,RealScalar,IntegerScalar,StrScalar,LogicalScalar,&
                                                                       RealArray,IntegerArray,StrArray)
!===================================================================================================================================
! Subroutine to read attributes from HDF5 file.
!===================================================================================================================================
! MODULES
USE MODH_Globals
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER(HID_T), INTENT(IN)           :: Loc_ID_in
INTEGER                              :: nVal
CHARACTER(LEN=*), INTENT(IN)         :: AttribName
CHARACTER(LEN=*),OPTIONAL,INTENT(IN) :: DatasetName
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL              ,OPTIONAL :: RealArray(nVal)
INTEGER           ,OPTIONAL :: IntegerArray(nVal)
REAL              ,OPTIONAL :: RealScalar
INTEGER           ,OPTIONAL :: IntegerScalar
LOGICAL           ,OPTIONAL :: LogicalScalar
CHARACTER(LEN=255),OPTIONAL :: StrScalar
CHARACTER(LEN=255),OPTIONAL :: StrArray(nVal)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER(HID_T)                 :: Attr_ID,Type_ID,Loc_ID
INTEGER(HSIZE_T), DIMENSION(1) :: Dimsf
INTEGER                        :: IntToLog,i
!===================================================================================================================================
LOGWRITE(*,*)' READ ATTRIBUTE "',TRIM(AttribName),'" FROM HDF5 FILE...'
Dimsf(1)=nVal
Loc_ID=Loc_ID_in
IF(PRESENT(DatasetName))THEN
  ! Open dataset
  IF(TRIM(DataSetName).NE.'') CALL H5DOPEN_F(File_ID, TRIM(DatasetName),Loc_ID, iError)
END IF
! Create scalar data space for the attribute.
! Create the attribute for group Loc_ID.
CALL H5AOPEN_F(Loc_ID, TRIM(AttribName), Attr_ID, iError)
! Write the attribute data.
IF(PRESENT(RealArray))THEN
  RealArray=0.
  CALL H5AREAD_F(Attr_ID, H5T_NATIVE_DOUBLE, RealArray, Dimsf, iError)
END IF
IF(PRESENT(RealScalar))THEN
  RealScalar=0.
  CALL H5AREAD_F(Attr_ID, H5T_NATIVE_DOUBLE, RealScalar, Dimsf, iError)
END IF
IF(PRESENT(IntegerArray))THEN
  IntegerArray=0
  CALL H5AREAD_F(Attr_ID, H5T_NATIVE_INTEGER, IntegerArray, Dimsf, iError)
END IF
IF(PRESENT(IntegerScalar))THEN
  IntegerScalar=0
  CALL H5AREAD_F(Attr_ID, H5T_NATIVE_INTEGER, IntegerScalar, Dimsf, iError)
END IF
IF(PRESENT(LogicalScalar))THEN
  LogicalScalar=.FALSE.
  CALL H5AREAD_F(Attr_ID, H5T_NATIVE_INTEGER, IntToLog, Dimsf, iError)
  LogicalScalar=(inttolog.EQ.1)
END IF
IF(PRESENT(StrScalar))THEN
  StrScalar=''
  CALL H5AGET_TYPE_F(Attr_ID, Type_ID, iError)  ! Get HDF5 data type for character string
  CALL H5AREAD_F(Attr_ID, Type_ID, StrScalar, Dimsf, iError)
  CALL H5TCLOSE_F(Type_ID, iError)
  LOGWRITE(*,*)' SCALAR STRING READ "',TRIM(StrScalar)
END IF
IF(PRESENT(StrArray))THEN
  DO i=1,nVal
    StrArray(i)=''
  END DO
  CALL H5AGET_TYPE_F(Attr_ID, Type_ID, iError)  ! Get HDF5 data type for character string
  CALL H5AREAD_F(Attr_ID, Type_ID, StrArray, Dimsf, iError)
  CALL H5TCLOSE_F(Type_ID, iError)
  DO i=1,nVal
    LOGWRITE(*,*)' ARRAY STRING READ "',TRIM(StrArray(i))
  END DO
END IF
! Close the attribute.
CALL H5ACLOSE_F(Attr_ID, iError)
IF(Loc_ID.NE.Loc_ID_in)THEN
  ! Close the dataset and property list.
  CALL H5DCLOSE_F(Loc_ID, iError)
END IF
LOGWRITE(*,*)'...DONE!'
END SUBROUTINE ReadAttributeFromHDF5

END MODULE MODH_HDF5_Input
