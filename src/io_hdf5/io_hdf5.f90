#include "hopest_f.h"

MODULE MODH_IO_HDF5
!===================================================================================================================================
! Add comments please!
!===================================================================================================================================
! MODULES
USE HDF5
USE MODH_Globals,ONLY: iError
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------
INTEGER(HID_T)           :: File_ID
INTEGER(HSIZE_T),POINTER :: HSize(:)
INTEGER                  :: nDims
INTEGER                  :: MPIInfo !for lustre file system

INTERFACE InitIO
  MODULE PROCEDURE InitIO_HDF5
END INTERFACE

INTERFACE OpenDataFile
  MODULE PROCEDURE OpenHDF5File
END INTERFACE

INTERFACE CloseDataFile
  MODULE PROCEDURE CloseHDF5File
END INTERFACE

!===================================================================================================================================

CONTAINS

SUBROUTINE InitIO_HDF5()
!===================================================================================================================================
! Initialize HDF5 IO
!===================================================================================================================================
! MODULES
USE MODH_Globals
USE MODH_ReadInTools,        ONLY:GETLOGICAL,CNTSTR
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
#ifdef MPI
  CALL MPI_Info_Create(MPIInfo, iError)

  !normal case:
  MPIInfo=MPI_INFO_NULL
#endif /* MPI */
END SUBROUTINE InitIO_HDF5



SUBROUTINE OpenHDF5File(FileString,create,single,communicatorOpt)
!===================================================================================================================================
! Open HDF5 file and groups
!===================================================================================================================================
! MODULES
USE MODH_Globals
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
CHARACTER(LEN=*),INTENT(IN)   :: FileString
LOGICAL,INTENT(IN)            :: create
LOGICAL,INTENT(IN)            :: single
INTEGER,INTENT(IN),OPTIONAL   :: communicatorOpt
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER(HID_T)                 :: Plist_ID
!===================================================================================================================================
LOGWRITE(*,'(A)')'  OPEN HDF5 FILE "',TRIM(FileString),'" ...'

! Initialize FORTRAN predefined datatypes
CALL H5OPEN_F(iError)

! Setup file access property list with parallel I/O access (MPI) or with default property list.
CALL H5PCREATE_F(H5P_FILE_ACCESS_F, Plist_ID, iError)
#ifdef MPI
IF(PRESENT(communicatorOpt))THEN
  comm=communicatorOpt
ELSE
  comm=MPI_COMM_WORLD
END IF

IF(.NOT.single)THEN
  CALL H5PSET_FAPL_MPIO_F    (Plist_ID, comm, MPIInfo,       iError)
END IF
#endif /* MPI */

! Open the file collectively.
IF(create)THEN
  CALL H5FCREATE_F(TRIM(FileString), H5F_ACC_TRUNC_F, File_ID, iError, access_prp = Plist_ID)
ELSE
  write(*,*) FileString
  CALL H5FOPEN_F(  TRIM(FileString), H5F_ACC_RDWR_F,  File_ID, iError, access_prp = Plist_ID)
END IF
CALL H5PCLOSE_F(Plist_ID, iError)
LOGWRITE(*,*)'...DONE!'
END SUBROUTINE OpenHDF5File



SUBROUTINE CloseHDF5File()
!===================================================================================================================================
! Close HDF5 file and groups
!===================================================================================================================================
! MODULES
USE MODH_Globals,ONLY:UNIT_stdOut,UNIT_logOut,Logging
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
LOGWRITE(*,'(A)')'  CLOSE HDF5 FILE...'
! Close file
CALL H5FCLOSE_F(File_ID, iError)
! Close FORTRAN predefined datatypes.
CALL H5CLOSE_F(iError)
File_ID=0
LOGWRITE(*,*)'...DONE!'
END SUBROUTINE CloseHDF5File

END MODULE MODH_io_HDF5
