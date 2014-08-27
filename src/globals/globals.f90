#include "hopest_f.h"

MODULE MODH_Globals
!===================================================================================================================================
! Provides parameters, used globally (please use EXTREMLY carefully!) 
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES 
!-----------------------------------------------------------------------------------------------------------------------------------
INTEGER           ::HopestMode ! 1: Hopest run in serial preproc mode with mesh refinement/HDF5 output
                               ! 2: Hopest run in parallel production mode for linking to solver 
                               ! 3: Hopest run in VTK mode for high order visualization

INTEGER,PARAMETER ::UNIT_stdOut=6
INTEGER,PARAMETER ::UNIT_logOut=133
INTEGER           ::UNIT_errOut=999
LOGICAL           ::Logging
CHARACTER(LEN=255)::ErrorFileName='NOT_SET'
INTEGER           ::iError
REAL              ::StartTime
INTEGER           ::myRank
INTEGER           ::nProcessors
LOGICAL           ::MPIRoot
INTEGER           ::comm
#ifdef MPI
#include "mpif.h"
INTEGER           :: MPIStatus(MPI_STATUS_SIZE)
#endif

INTERFACE Abort
  MODULE PROCEDURE Abort
END INTERFACE Abort

INTERFACE INTSTAMP
  MODULE PROCEDURE INTSTAMP
END INTERFACE INTSTAMP

INTERFACE TIMESTAMP
  MODULE PROCEDURE TIMESTAMP
END INTERFACE

INTERFACE RUNTIME
  MODULE PROCEDURE RUNTIME
END INTERFACE

INTERFACE GETFREEUNIT
  MODULE PROCEDURE GETFREEUNIT
END INTERFACE GETFREEUNIT

INTERFACE CreateErrFile
  MODULE PROCEDURE CreateErrFile
END INTERFACE CreateErrFile

!===================================================================================================================================
CONTAINS

SUBROUTINE Abort(SourceFile,SourceLine,CompDate,CompTime,ErrorMessage,IntInfoOpt,RealInfoOpt)
!===================================================================================================================================
! Terminate program correctly if an error has occurred (important in MPI mode!).
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
CHARACTER(LEN=*)                  :: SourceFile      ! Source file where error has occurred
INTEGER                           :: SourceLine      ! Line in source file
CHARACTER(LEN=*)                  :: CompDate        ! Compilation date
CHARACTER(LEN=*)                  :: CompTime        ! Compilation time
CHARACTER(LEN=*)                  :: ErrorMessage    ! Error message
INTEGER,OPTIONAL                  :: IntInfoOpt      ! Error info (integer)
REAL,OPTIONAL                     :: RealInfoOpt     ! Error info (real)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!   There is no way back!
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                           :: IntInfo         ! Error info (integer)
REAL                              :: RealInfo        ! Error info (real)
!===================================================================================================================================
IntInfo  = MERGE(IntInfoOpt ,999 ,PRESENT(IntInfoOpt) )
RealInfo = MERGE(RealInfoOpt,999.,PRESENT(RealInfoOpt))
WRITE(UNIT_stdOut,*)
WRITE(UNIT_stdOut,*)'_____________________________________________________________________________'
WRITE(UNIT_stdOut,*)'Program abort caused on Proc ',myRank,' in File : ',TRIM(SourceFile),' Line ',SourceLine
WRITE(UNIT_stdOut,*)'This file was compiled at ',TRIM(CompDate),'  ',TRIM(CompTime)
WRITE(UNIT_stdOut,'(A10,A)',ADVANCE='NO')'Message: ',TRIM(ErrorMessage)
IF(IntInfo  .NE. 999 ) WRITE(UNIT_stdOut,'(I8)',ADVANCE='NO')IntInfo
IF(RealInfo .NE. 999.) WRITE(UNIT_stdOut,'(E16.8)')RealInfo
WRITE(UNIT_stdOut,*)
WRITE(UNIT_stdOut,'(A,A,A)')'See ',TRIM(ErrorFileName),' for more details'
WRITE(UNIT_stdOut,*)
!CALL delete()
#ifdef MPI
CALL MPI_ABORT(MPI_COMM_WORLD,iError)
#endif
STOP 0001
END SUBROUTINE Abort


SUBROUTINE CreateErrFile()
!===================================================================================================================================
! Open file for error output
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
INTEGER                        :: OpenStat
LOGICAL                        :: isOpen
!===================================================================================================================================
INQUIRE(UNIT=UNIT_errOut,OPENED=isOpen)
IF(.NOT.isOpen)THEN
  OPEN(UNIT=UNIT_errOut,  &
       FILE=ErrorFileName,&
       STATUS='REPLACE',  &
       ACTION='WRITE',    &
       IOSTAT=OpenStat)
END IF
END SUBROUTINE CreateErrFile


FUNCTION INTSTAMP(Nam,Num)
!===================================================================================================================================
! Creates an integer stamp that will afterwards be given to the SOUBRUTINE timestamp
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
CHARACTER(LEN=*)   :: Nam      ! Name
INTEGER            :: Num      ! Number
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
CHARACTER(LEN=200) :: IntStamp ! The stamp
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
WRITE(IntStamp,'(A,A5,I6.6)')TRIM(Nam),'_Proc',Num
END FUNCTION INTSTAMP



FUNCTION TIMESTAMP(Filename,Time)
!===================================================================================================================================
! Creates a timestamp, consistent of a filename (project name + processor) and current time niveau
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
CHARACTER(LEN=*)   :: Filename  ! (file)name
REAL               :: Time      ! time
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
CHARACTER(LEN=255) :: TimeStamp ! the complete timestamp
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER            :: i         ! loop variable
!===================================================================================================================================
WRITE(TimeStamp,'(F15.7)')Time
! Replace spaces with 0's
DO i=1,LEN(TRIM(TimeStamp))
  IF(TimeStamp(i:i).EQ.' ') TimeStamp(i:i)='0'
END DO
TimeStamp=TRIM(Filename)//'_'//TRIM(TimeStamp)
END FUNCTION TIMESTAMP



FUNCTION RUNTIME(Comm)
!===================================================================================================================================
! Calculates current time (own function because of a laterMPI implementation)
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN),OPTIONAL    :: Comm
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL                            :: RunTime
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
!===================================================================================================================================
#ifdef MPI
IF(PRESENT(Comm))THEN
  CALL MPI_BARRIER(Comm,iError)
ELSE
  CALL MPI_BARRIER(MPI_COMM_WORLD,iError)
END IF
RunTime=MPI_WTIME()
#else
CALL CPU_TIME(RunTime)
#endif
END FUNCTION RUNTIME



FUNCTION GETFREEUNIT()
!===================================================================================================================================
! Get unused file unit number
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
INTEGER :: GetFreeUnit ! File unit number
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
LOGICAL :: connected
!===================================================================================================================================
GetFreeUnit=55
INQUIRE(UNIT=GetFreeUnit, OPENED=connected)
IF(connected)THEN
  DO  
    GetFreeUnit=GetFreeUnit+1
    INQUIRE(UNIT=GetFreeUnit, OPENED=connected)
    IF(.NOT.connected)EXIT
  END DO
END IF
END FUNCTION GETFREEUNIT

PURE FUNCTION CROSS(v1,v2)
!===================================================================================================================================
! computes the cross product of to 3 dimensional vectpors: cross=v1 x v2
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
    IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
    REAL,INTENT(IN) :: v1(3)    ! 
    REAL,INTENT(IN) :: v2(3)    ! 
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
    REAL            :: CROSS(3) !
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
!===================================================================================================================================
  CROSS=(/v1(2)*v2(3)-v1(3)*v2(2),v1(3)*v2(1)-v1(1)*v2(3),v1(1)*v2(2)-v1(2)*v2(1)/)
END FUNCTION CROSS

END MODULE MODH_Globals
