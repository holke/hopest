#include "hopest_f.h"
MODULE MODH_P4EST_Vars
!===================================================================================================================================
! Contains global variables provided by the mesh routines
!===================================================================================================================================
! MODULES
USE,INTRINSIC :: ISO_C_BINDING
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PUBLIC
SAVE
!-----------------------------------------------------------------------------------------------------------------------------------
! P4EST related data structures 
!-----------------------------------------------------------------------------------------------------------------------------------
CHARACTER(LEN=255)          :: p4estFile          ! name of hdf5 meshfile (write with ending .h5!)

TYPE(C_PTR)                 :: p4est              ! c pointers to p4est structures
TYPE(C_PTR)                 :: mesh               !
TYPE(C_PTR)                 :: connectivity       !
TYPE(C_PTR)                 :: geom               !

INTEGER(KIND=4)             :: nHalfFaces         ! number of mortar sides
INTEGER(KIND=4)             :: IntSize            ! used to transform INT coords/levels to REAL coords/levels: REAL=1/inssize*INT  [0. ; 1.]
REAL                        :: sIntSize           ! 1./REAL(intsize)
INTEGER(KIND=4),POINTER     :: QuadToTree(:)      ! from quadrant to tree ( ~ new element ID to old element ID) 
INTEGER(KIND=4),POINTER     :: QuadToQuad(:,:)    ! p4est quadrant connectivity (1:6,1:nElems) => neighbor quadrant
INTEGER(KIND=1),POINTER     :: QuadToFace(:,:)    ! p4est face connectivity (1:6,1:nElems) => neighbor faceId + orientation + non-conform info
INTEGER(KIND=4),POINTER     :: QuadToHalf(:,:)    ! p4est face connectivity for mortars (1:4,1:nHalfFaces), ( ~small sides)
INTEGER(KIND=4),ALLOCATABLE :: QuadCoords(:,:)    ! p4est Integer coordinates of first quadrant node (xyz,nElems)
INTEGER(KIND=1),ALLOCATABLE :: QuadLevel(:)       ! p4est Integer Level of quadrant (use to compute quadrant size
INTEGER(KIND=C_INT32_T),POINTER :: TreeToBC(:,:)  ! (1...6,1,nTrees) first index is p4est local side ID +1!!

!-----------------------------------------------------------------------------------------------------------------------------------
INTEGER,PARAMETER   :: EdgeToElemNode(1:2,1:12) = RESHAPE((/ 1, 2,&  ! CGNS corner nodes mapped 
                                                             4, 3,&  ! to p4est edges
                                                             5, 6,&
                                                             8, 7,&
                                                             1, 4,&
                                                             2, 3,&
                                                             5, 8,&
                                                             6, 7,&
                                                             1, 5,&
                                                             2, 6,&
                                                             4, 8,&
                                                             3, 7 /),(/2,12/))
INTEGER,PARAMETER   :: H2P_FaceMap(1:6)     =  (/4,2,1,3,0,5/)     !mapping from local face order (CGNS) to p4est face
INTEGER,PARAMETER   :: P2H_FaceMap(0:5)     =  (/5,3,2,4,1,6/)     !mapping from local face order (CGNS) to p4est face
INTEGER,PARAMETER   :: H2P_VertexMap(1:8)   =  (/0,1,3,2,4,5,7,6/) !mapping from local node order (CGNS) to p4est node order 
INTEGER,PARAMETER   :: P2H_VertexMap(0:7)   =  (/1,2,4,3,5,6,8,7/) !mapping from local node order (CGNS) to p4est node order 

! mapping from HOPEST node of local sides to P4EST nodes of local sides
INTEGER,PARAMETER   :: H2P_FaceNodeMap(1:4,1:6) = &
                                      RESHAPE((/ 0,2,3,1,&
                                                 0,1,3,2,&
                                                 0,1,3,2,&
                                                 1,0,2,3,&
                                                 0,2,3,1,&
                                                 0,1,3,2 /),(/4,6/))

! mapping from P4EST node of local sides to HOPEST node of local sides
INTEGER,PARAMETER   :: P2H_FaceNodeMap(0:3,0:5) = &
                                      RESHAPE((/ 1,4,2,3,&
                                                 1,2,4,3,&
                                                 1,2,4,3,&
                                                 2,1,3,4,&
                                                 1,4,2,3,&
                                                 1,2,4,3 /),(/4,6/))

! Mapping matrices for computation of same node on adjacent face, see paper Burstedde p4est, 2011
! Node1= P4P(P4Q(P4R(Face0,Face1),orientation),Node0)
INTEGER,PARAMETER   :: P4R(0:5,0:5) = TRANSPOSE(RESHAPE((/ 0,1,1,0,0,1,&
                                                           2,0,0,1,1,0,&
                                                           2,0,0,1,1,0,&
                                                           0,2,2,0,0,1,&
                                                           0,2,2,0,0,1,&
                                                           2,0,0,2,2,0 /),(/6,6/)))

INTEGER,PARAMETER   :: P4Q(0:2,0:3) = TRANSPOSE(RESHAPE((/ 1,2,5,6,&
                                                           0,3,4,7,&
                                                           0,4,3,7 /),(/4,3/)))

INTEGER,PARAMETER   :: P4P(0:7,0:3) = TRANSPOSE(RESHAPE((/ 0,1,2,3,&
                                                           0,2,1,3,&
                                                           1,0,3,2,&
                                                           1,3,0,2,&
                                                           2,0,3,1,&
                                                           2,3,0,1,&
                                                           3,1,2,0,&
                                                           3,2,1,0 /),(/4,8/)))
INTEGER,PARAMETER   :: H_MortarCase(1:4,1:4) = &                              !  first CGNS node and second CGNS node->Mortar Case [1:8]
                                      TRANSPOSE(RESHAPE((/  0,1,0,2,&                         ! (1,2)->1, (1,4)->2
                                                            3,0,4,0,&                         ! (2,1)->3, (2,3)->4
                                                            0,5,0,6,&                         ! (3,2)->5, (3,4)->6
                                                            7,0,8,0 /),(/4,4/)))              ! (4,1)->7, (4,3)->8
INTEGER,PARAMETER   :: P2H_MortarMap(0:3,1:8) = &                 !p4est mortar ID, MortarCase -> iMortar CGNS
                                      RESHAPE((/ 1,2,3,4,&        ! iMortar = P2H_MortarMap(iPMortar, H_MortarCase( node1, node2) ) 
                                                 1,3,2,4,&
                                                 2,1,4,3,&
                                                 2,4,1,3,&
                                                 4,2,3,1,&
                                                 4,3,2,1,&
                                                 3,1,4,2,&
                                                 3,4,1,2 /),(/4,8/))
INTEGER,PARAMETER   :: P_FaceToEdge(0:3,0:5) = &  !mapping from face edges 0...3 (zordered) for each face 0..5 -> element edges  0..11
                                      RESHAPE((/  4, 6, 8,10,&     
                                                  5, 7, 9,11,&
                                                  0, 2, 8, 9,&
                                                  1, 3,10,11,&
                                                  0, 1, 4, 5,&
                                                  2, 3, 6, 7 /),(/4,6/))
INTEGER,PARAMETER   :: P_EdgeToFaces(1:6,0:11) = & !mapping from element edges  0..11 -> first and second adjacent face 0...6 , i/j direction(0/1) and lower/upper bound(0/1)
                                      RESHAPE((/  2,1,0,4,1,0,&    !edge 0: Face2,j=0-Face4,j=0
                                                  3,1,0,4,1,1,&    !edge 1: Face3,j=0-Face4,j=N
                                                  2,1,1,5,1,0,&    !edge 2: Face2,j=N-Face5,j=0
                                                  3,1,1,5,1,1,&    !edge 3: Face3,j=N-Face5,j=N
                                                  0,1,0,4,0,0,&    !edge 4: Face0,j=0-Face4,i=0
                                                  1,1,0,4,0,1,&    !edge 5: Face1,j=0-Face4,i=N
                                                  0,1,1,5,0,0,&    !edge 6: Face0,j=N-Face5,i=0
                                                  1,1,1,5,0,1,&    !edge 7: Face1,j=N-Face5,i=N
                                                  0,0,0,2,0,0,&    !edge 8: Face0,i=0-Face2,i=0
                                                  1,0,0,2,0,1,&    !edge 9: Face1,i=0-Face2,i=N
                                                  0,0,1,3,0,0,&    !edge10: Face0,i=N-Face3,i=0
                                                  1,0,1,3,0,1 &    !edge11: Face1,i=N-Face3,i=N
                                                  /),(/6,12/))

!===================================================================================================================================


END MODULE MODH_P4EST_Vars
