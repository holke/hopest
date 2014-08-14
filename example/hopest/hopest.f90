#include <hopest_f.h>

PROGRAM Hopest
!===================================================================================================================================
! Control program of the Flexi code. Initialization of the computation
!===================================================================================================================================
! MODULES
USE MODH_Globals
USE MODH_ReadInTools,  ONLY:IgnoredStrings,GETINT
USE MODH_HopestMesh,   ONLY:HopestMesh
USE MODH_HopestSolver, ONLY:HopestSolver,PrepareMesh,FinalizeHopestSolver
USE MODH_MPI,          ONLY:InitMPI
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
REAL    :: Time
!===================================================================================================================================
! Init MPI in here, if used as lib for solver, this is done by the solver
CALL InitMPI()

! Measure init duration
StartTime=RUNTIME()

hopestMode=GETINT('HopestMode')
SELECT CASE(hopestMode)
CASE(1) ! HOPEST reads HDF5 builds P4EST + static refinement
  IF(.NOT.MPIRoot) STOP 'HOPEST Mesh only runs in single mode!' 
  CALL HopestMesh()
CASE(2) ! use HOPEST as INPUT for Flexi
  CALL HopestSolver()
  CALL PrepareMesh()
  CALL FinalizeHopestSolver()
CASE(3) ! HOPEST VTK ?!
END SELECT

CALL IgnoredStrings()

! Measure simulation duration
Time=RUNTIME()
#ifdef MPI
CALL MPI_FINALIZE(iError)
IF(iError .NE. 0) &
  CALL abort(__STAMP__, &
  'MPI finalize error',iError)
#endif
SWRITE(UNIT_stdOut,'(132("="))')
SWRITE(UNIT_stdOut,'(A,F8.2,A)') ' HOPEST FINISHED! [',Time-StartTime,' sec ]'
SWRITE(UNIT_stdOut,'(132("="))')

END PROGRAM Hopest
