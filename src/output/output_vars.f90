MODULE MODH_Output_Vars
!===================================================================================================================================
! Contains global variables provided by the output routines
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PUBLIC
SAVE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES 
!-----------------------------------------------------------------------------------------------------------------------------------
REAL,PARAMETER               :: FileVersion=0.1
CHARACTER(LEN=255),PARAMETER :: ProgramName='Hopest'
CHARACTER(LEN=255)           :: ProjectName
INTEGER                      :: outputFormat=0           ! =0: visualization off, >0 visualize
LOGICAL                      :: OutputInitIsDone=.FALSE.
!===================================================================================================================================
END MODULE MODH_Output_Vars
