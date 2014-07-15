#include "../hopest_f.h"

MODULE MOD_Output
!===================================================================================================================================
! Add comments please!
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES 
!-----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------

INTERFACE InitOutput
  MODULE PROCEDURE InitOutput
END INTERFACE

INTERFACE FinalizeOutput
  MODULE PROCEDURE FinalizeOutput
END INTERFACE

PUBLIC:: InitOutput,FinalizeOutput
!===================================================================================================================================

CONTAINS

SUBROUTINE InitOutput()
!===================================================================================================================================
! Initialize all output (and analyze) variables.
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_Output_Vars
USE MOD_ReadInTools       ,ONLY:GETSTR,GETLOGICAL,GETINT
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                        :: OpenStat
CHARACTER(LEN=8)               :: StrDate
CHARACTER(LEN=10)              :: StrTime
CHARACTER(LEN=255)             :: LogFile
!===================================================================================================================================
IF(OutputInitIsDone)THEN
  CALL abort(__STAMP__,'InitOutput not ready to be called or already called.')
END IF
SWRITE(UNIT_StdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A)') ' INIT OUTPUT...'

! Name for all output files
ProjectName=GETSTR('ProjectName')
Logging    =GETLOGICAL('Logging','.FALSE.')

WRITE(ErrorFileName,'(A,A8,I6.6,A4)')TRIM(ProjectName),'_ERRORS_',myRank,'.out'

OutputFormat = GETINT('OutputFormat','1')
! Open file for logging
IF(Logging)THEN
  WRITE(LogFile,'(A,A1,I6.6,A4)')TRIM(ProjectName),'_',myRank,'.log'
  OPEN(UNIT=UNIT_logOut,  &
       FILE=LogFile,      &
       STATUS='UNKNOWN',  &
       ACTION='WRITE',    &
       POSITION='APPEND', &
       IOSTAT=OpenStat)
  CALL DATE_AND_TIME(StrDate,StrTime)
  WRITE(UNIT_logOut,*)
  WRITE(UNIT_logOut,'(132("#"))')
  WRITE(UNIT_logOut,*)
  WRITE(UNIT_logOut,*)'STARTED LOGGING FOR PROC',myRank,' ON ',StrDate(7:8),'.',StrDate(5:6),'.',StrDate(1:4),' | ',&
                      StrTime(1:2),':',StrTime(3:4),':',StrTime(5:10)
END IF  ! Logging

OutputInitIsDone =.TRUE.
SWRITE(UNIT_stdOut,'(A)')' INIT OUTPUT DONE!'
SWRITE(UNIT_StdOut,'(132("-"))')
END SUBROUTINE InitOutput


SUBROUTINE FinalizeOutput()
!===================================================================================================================================
! Deallocate global variables
!===================================================================================================================================
! MODULES
USE MOD_Output_Vars,ONLY:OutputInitIsDone
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
OutputInitIsDone = .FALSE.
END SUBROUTINE FinalizeOutput

END MODULE MOD_Output
