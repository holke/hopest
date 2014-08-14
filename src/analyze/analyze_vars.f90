#include "hopest_f.h"
MODULE MODH_Analyze_Vars
!===================================================================================================================================
! Contains global variables provided by the analyze routines
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PUBLIC
SAVE
!-----------------------------------------------------------------------------------------------------------------------------------

LOGICAL          :: checkJacobian
INTEGER          :: Nanalyze
REAL,ALLOCATABLE :: Vdm_analyze(:,:)
REAL,ALLOCATABLE :: D_Ngeo_out(:,:)

END MODULE MODH_Analyze_Vars
