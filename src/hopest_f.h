!===================================================================================================================================
! Here, preprocessor variables for different equation systems and abbreviations for specific expressions are defined
!===================================================================================================================================

! Abbrevations
#define __STAMP__ __FILE__,__LINE__,__DATE__,__TIME__

#ifdef GNU
#  define IEEE_IS_NAN ISNAN
#endif

#ifdef MPI
#  define SWRITE IF(MPIRoot) WRITE
#  define IPWRITE(a,b) WRITE(a,b)myRank,
#else
#  define SWRITE WRITE
#  define IPWRITE WRITE
#endif
#define ERRWRITE(a,b) CALL CreateErrFile(); WRITE(UNIT_errOut,b)
#define LOGWRITE(a,b) IF(Logging) WRITE(UNIT_logOut,b)
#define SDEALLOCATE(A) IF(ALLOCATED(A)) DEALLOCATE(A)

! Predefined "PARAMETER-like" variables

! Entry position in BC
#define BC_SIZE   4
#define BC_TYPE   1
#define BC_CURVED 2
#define BC_STATE  3
#define BC_ALPHA  4

!entry positions in ElemInfo 
#define ELEM_InfoSize     6        /*number of entry in each line of ElemInfo*/
#define ELEM_Type         1        /*entry position in ElemInfo */
#define ELEM_Zone         2           
#define ELEM_FirstSideInd 3
#define ELEM_LastSideInd  4
#define ELEM_FirstNodeInd 5
#define ELEM_LastNodeInd  6

!entry positions in SideInfo 
#define SIDE_InfoSize     4        /*number of entry in each line of SideInfo*/
#define SIDE_Type         1         /*entry position in SideInfo */
#define SIDE_ID           2
#define SIDE_nbElemID     3
#define SIDE_BCID         4
