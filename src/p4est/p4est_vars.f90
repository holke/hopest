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
INTEGER(KIND=4),ALLOCATABLE :: TreeToQuad(:,:)    ! from tree to quad range (2,nElems), entry 1: firstInd-1, entry2:lastInd 
INTEGER(KIND=4),POINTER     :: QuadToQuad(:,:)    ! p4est quadrant connectivity (1:6,1:nQuads) => neighbor quadrant
INTEGER(KIND=1),POINTER     :: QuadToFace(:,:)    ! p4est face connectivity (1:6,1:nQuads) => neighbor faceId + orientation + non-conform info
INTEGER(KIND=4),POINTER     :: QuadToHalf(:,:)    ! p4est face connectivity for mortars (1:4,1:nHalfFaces), ( ~small sides)
INTEGER(KIND=4),ALLOCATABLE :: QuadCoords(:,:)    ! p4est Integer coordinates of first quadrant node (xyz,nQuads)
INTEGER(KIND=1),ALLOCATABLE :: QuadLevel(:)       ! p4est Integer Level of quadrant (use to compute quadrant size

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
                                                           2,0,0,2,2,0,&
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

!===================================================================================================================================


END MODULE MODH_P4EST_Vars
