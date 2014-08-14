#include "hopest_f.h"

MODULE MODH_P4EST
!===================================================================================================================================
! Add comments please!
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------

! Public Part ----------------------------------------------------------------------------------------------------------------------
INTERFACE InitP4EST
  MODULE PROCEDURE InitP4EST
END INTERFACE

INTERFACE BuildMeshFromP4EST
  MODULE PROCEDURE BuildMeshFromP4EST
END INTERFACE

INTERFACE getHFlip
  MODULE PROCEDURE getHFlip
END INTERFACE

INTERFACE BuildBCs
  MODULE PROCEDURE BuildBCs
END INTERFACE

INTERFACE testHOabc
  MODULE PROCEDURE testHOabc
END INTERFACE

INTERFACE FinalizeP4EST
  MODULE PROCEDURE FinalizeP4EST
END INTERFACE

PUBLIC::InitP4EST
PUBLIC::BuildMeshFromP4EST
PUBLIC::getHFlip
PUBLIC::BuildBCs
PUBLIC::testHOabc
PUBLIC::FinalizeP4EST
!===================================================================================================================================

CONTAINS


SUBROUTINE InitP4EST()
!===================================================================================================================================
! Subroutine to translate p4est mesh datastructure to HOPR datastructure
!===================================================================================================================================
! MODULES
USE, INTRINSIC :: ISO_C_BINDING
USE MODH_Globals,       ONLY: hopestMode
USE MODH_P4EST_Vars,    ONLY: p4estFile
USE MODH_P4EST_Binding, ONLY: p4_initvars
USE MODH_Output_Vars,   ONLY: Projectname
USE MODH_ReadInTools,   ONLY: GETSTR
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
CALL p4_initvars()
p4estFile = GETSTR('p4estFile',TRIM(ProjectName)//'.p4est')

END SUBROUTINE InitP4EST



SUBROUTINE BuildMeshFromP4EST()
!===================================================================================================================================
! Subroutine to translate p4est mesh datastructure to HOPR datastructure
!===================================================================================================================================
! MODULES
USE, INTRINSIC :: ISO_C_BINDING
USE MODH_Globals
USE MODH_Mesh_Vars
USE MODH_P4EST_Vars
USE MODH_P4EST_Binding
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
TYPE(C_PTR)                 :: QT,QQ,QF,QH,TB
TYPE(tElem),POINTER         :: aElem,nbElem
TYPE(tSide),POINTER         :: aSide
INTEGER                     :: iElem,iMortar,jMortar,iTree
INTEGER                     :: PSide,PnbSide,nbSide,iSide
INTEGER                     :: nbElemInd
INTEGER                     :: PMortar,PFlip,HFlip,QHInd
INTEGER                     :: iLocSide
INTEGER                     :: StartElem,EndElem
INTEGER                     :: BClocSide,BCindex
INTEGER(KIND=8)             :: offsetElemTmp,nGlobalElemsTmp
LOGICAL                     :: doConnection
!===================================================================================================================================
SWRITE(UNIT_stdOut,'(A)')'GENERATE HOPEST MESH FROM P4EST ...'
SWRITE(UNIT_StdOut,'(132("-"))')

! build p4est mesh
CALL p4_build_mesh(p4est,mesh)
! Get arrays from p4est: use pointers for c arrays (QT,QQ,..), duplicate data for QuadCoords,Level
CALL p4_get_mesh_info(p4est,mesh,nElems,nGlobalElemsTmp,offsetElemTmp,nHalfFaces,nTrees)
offsetElem=offsetElemTmp
nGlobalElems=nGlobalElemsTmp

ALLOCATE(QuadCoords(3,nElems),QuadLevel(nElems)) ! big to small flip
QuadCoords=0
QuadLevel=0

CALL p4_get_quadrants(p4est,mesh,nElems,nHalfFaces,& !IN
                      intsize,QT,QQ,QF,QH,QuadCoords,QuadLevel)              !OUT
sIntSize=1./REAL(Intsize)

CALL C_F_POINTER(QT,QuadToTree,(/nElems/))
CALL C_F_POINTER(QQ,QuadToQuad,(/6,nElems/))
CALL C_F_POINTER(QF,QuadToFace,(/6,nElems/))
IF(nHalfFaces.GT.0) CALL C_F_POINTER(QH,QuadToHalf,(/4,nHalfFaces/))

! Get boundary conditions from p4est
CALL p4_get_bcs(p4est,TB)
CALL C_F_POINTER(TB,TreeToBC,(/6,nTrees/))

!----------------------------------------------------------------------------------------------------------------------------
!             Start to build p4est datastructure in HOPEST
!----------------------------------------------------------------------------------------------------------------------------
!                              ELEMENTS
!----------------------------------------------------------------------------------------------------------------------------

!read local ElemInfo from data file
iSide=0
ALLOCATE(Elems(1:nElems))
DO iElem=1,nElems
  Elems(iElem)%ep=>GETNEWELEM()
  aElem=>Elems(iElem)%ep
  aElem%Ind = iElem
  CALL CreateSides(aElem)
  DO iLocSide=1,6
    iSide=iSide+1
    aSide=>aElem%Side(iLocSide)%sp
    aSide%locSide=iLocSide
    aSide%tmp=0
    aSide%flip=-999
  END DO
END DO

DO iElem=1,nElems
  aElem=>Elems(iElem)%ep
  IF(useCurveds)THEN
    aElem%type=208 ! fully curved
  ELSE
    aElem%type=118 ! trilinear
  END IF
  DO iLocSide=1,6
    aSide=>aElem%Side(iLocSide)%sp
    ! Get P4est local side
    PSide=H2P_FaceMap(iLocSide)
    ! Get P4est neighbour side/flip/morter
    CALL EvalP4ESTConnectivity(QuadToFace(PSide+1,iElem),PnbSide,PFlip,PMortar)
    ! transform p4est orientation to HOPR flip (magic)
    HFlip=GetHFlip(PSide,PnbSide,PFlip)  !Hflip of neighbor side!!!
    IF(PMortar.EQ.4)THEN
      ! Neighbour side is mortar (4 sides), all neighbour element sides have same orientation and local side ind
      QHInd=QuadToQuad(PSide+1,iElem)+1
      aSide%nMortars=4
      aSide%MortarType=1             ! 1->4 case
      aSide%flip=HFlip
      ALLOCATE(aSide%MortarSide(4))
      DO jMortar=0,3
        nbElemInd=QuadToHalf(jMortar+1,QHInd)+1
        nbElem=>Elems(nbElemInd)%ep
        nbSide=P2H_FaceMap(PnbSide)

        iMortar=GetHMortar(jMortar,PSide,PnbSide,PFlip)

        aSide%MortarSide(iMortar)%sp=>nbElem%side(nbSide)%sp
        aSide%MortarSide(iMortar)%sp%flip=HFlip
      END DO ! iMortar
    ELSE
      nbElemInd=QuadToQuad(PSide+1,iElem)+1
      nbElem=>Elems(nbElemInd)%ep
      nbSide=P2H_FaceMap(PnbSide)

      IF((nbElemInd.EQ.iElem).AND.(nbSide.EQ.iLocSide))THEN  !p4est BC side
        ! this is a boundary side (periodic sides are not p4est BCsides!): 
        NULLIFY(aSide%connection)
        aSide%Flip=0
      ELSE
        !this is an inner side (either no mortar or small side mortar)
        aSide%connection=>nbElem%side(nbSide)%sp
        aSide%flip=HFlip
      END IF

      IF(PMortar.NE.-1) aSide%MortarType= - (PMortar+1)  ! Pmortar 0...3, small side belonging to  mortar group -> -1..-4
    END IF ! PMortar
  END DO !iLocSide
END DO !iElem

! set BC from tree to Quad
DO iElem=1,nElems
  aElem=>Elems(iElem)%ep
  iTree=QuadToTree(iElem)+1
  DO iLocSide=1,6
    aSide=>aElem%Side(iLocSide)%sp
    ! Get P4est local side
    PSide=H2P_FaceMap(iLocSide)
    aSide%BCindex=GET_BCINDEX_FROM_TREE(iTree,iElem,Pside)
    IF((aSide%BCindex.GT.0).AND.(ASSOCIATED(aSide%connection)))THEN !check if still periodic side (could be set to another BC!)
      IF(BoundaryType(aSide%BCindex,BC_TYPE).GT.1)THEN !if BCtype=0 & 1 should keep their connections!
        NULLIFY(aSide%connection)
        aSide%flip=0
      END IF
    END IF !check periodic
  END DO !iLocSide
END DO !iElem


nBCSides=0
nMortarSides=0
DO iElem=1,nElems
  aElem=>Elems(iElem)%ep
  DO iLocSide=1,6
    aSide=>aElem%Side(iLocSide)%sp
    IF(aSide%tmp.NE.0) CYCLE
    nSides=nSides+1
    aSide%tmp=-1
    IF(ASSOCIATED(aSide%connection)) aSide%connection%tmp=-1
    IF(aSide%BCindex.NE.0)THEN !side is BC or periodic side
      IF(ASSOCIATED(aSide%connection))THEN
        !nPeriodicSides=nPeriodicSides+1
      ELSE
        IF(aSide%MortarType.EQ.0)THEN !really a BC side
          nBCSides=nBCSides+1
        END IF
      END IF
    END IF
    IF(aSide%MortarType.GT.0) nMortarSides=nMortarSides+1
  END DO
END DO
nInnerSides=nSides-nBCSides-nMortarSides

! set master slave,  element with lower element ID is master (flip=0)
DO iElem=1,nElems
  aElem=>Elems(iElem)%ep
  DO iLocSide=1,6
    aSide=>aElem%Side(iLocSide)%sp
    IF(aSide%MortarType.GT.0)THEN
      aSide%flip=0
      DO iMortar=1,4
        IF(aSide%MortarSide(iMortar)%sp%flip.EQ.0) STOP 'Mortarside flip = 0'
      END DO
    ELSE
      IF(aSide%MortarType.EQ.0)THEN
        IF(ASSOCIATED(aSide%connection))THEN
          IF(aSide%connection%elem%ind.GT.iElem)THEN
            aSide%flip=0
          ELSEIF(aSide%connection%elem%ind.EQ.iElem)THEN
            ! special case: perdiodic within one element
            IF(aSide%locSide.LT.aSide%connection%locSide)&
              aSide%flip=0
          END IF
        END IF
      END IF
    END IF
  END DO
END DO

!sanity check
DO iElem=1,nElems
  aElem=>Elems(iElem)%ep
  DO iLocSide=1,6
    IF(aElem%Side(iLocSide)%sp%flip.LT.0) THEN
      WRITE(*,*) 'flip assignment failed, iElem= ',iElem,', iLocSide= ',iLocSide 
      STOP
    END IF
  END DO
END DO
END SUBROUTINE BuildMeshFromP4EST


SUBROUTINE EvalP4ESTConnectivity(Conn,nbSide,Flip,Mortar)
!===================================================================================================================================
! Subroutine to read the mesh from a mesh data file
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER(KIND=1),INTENT(IN)   :: Conn   ! p4est Side,Flip,Mortar encoding
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
INTEGER,INTENT(OUT)          :: nbSide ! Neighbour side in p4est convention: 0..5
INTEGER,INTENT(OUT)          :: Flip   ! Flip in p4est convention: 0..3
INTEGER,INTENT(OUT)          :: Mortar ! Mortar in p4est convention: 0..3,
                                       ! -1 if conformal, 4 if half-size neighbour
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                      :: tmp
!-----------------------------------------------------------------------------------------------------------------------------------
! The quad_to_quad list stores one value for each local quadrant's face.
! This value is in 0..local_num_quadrants-1 for local quadrants, or in
! local_num_quadrants + (0..ghost_num_quadrants-1) for ghost quadrants.
! The quad_to_face list has equally many entries which are either:
! 1. A value of v = 0..23 indicates one same-size neighbor.
!    This value is decoded as v = r * 6 + nf, where nf = 0..5 is the
!    neigbbor's connecting face number and r = 0..3 is the relative
!    orientation of the neighbor's face, see p8est_connectivity.h.
! 2. A value of v = 24..119 indicates a double-size neighbor.
!    This value is decoded as v = 24 + h * 24 + r * 6 + nf, where
!    r and nf are as above and h = 0..3 is the number of the subface.
! 3. A value of v = -24..-1 indicates four half-size neighbors.
!    In this case the corresponding quad_to_quad index points into the
!    quad_to_half array which stores four quadrant numbers per index,
!    and the orientation of the smaller faces follows from 24 + v.
!    The entries of quad_to_half encode between local and ghost quadrant
!    in the same way as the quad_to_quad values described above.
! A quadrant on the boundary of the forest sees itself and its face number.

SELECT CASE(Conn)
CASE(0:23)   ! 1. conformal neighbour
  nbSide = MOD(Conn,6)       ! 0..5
  Flip   = (Conn-nbSide)/6   ! 0..3
  Mortar = -1
CASE(24:119) ! 2. double-size neighbour
  tmp    = MOD(Conn,24)      ! 0..3
  nbSide = MOD(tmp,6)        ! 0..5  
  Flip   = (tmp-nbSide)/6    ! 0..3
  Mortar = (Conn-tmp-24)/24  ! 0..3 
CASE(-24:-1) ! 3. half-size neighbour
  tmp    = Conn+24
  nbSide = MOD(tmp,6)       ! 0..5
  Flip   = (tmp-nbSide)/6   ! 0..3
  Mortar = 4
CASE DEFAULT
  STOP 'This type of face connectivity does not exist, has to be -24<f<23'
END SELECT

END SUBROUTINE EvalP4ESTConnectivity


FUNCTION GET_BCINDEX_FROM_TREE(iTree,iElem,PSide)
!===================================================================================================================================
! assign either BCindex 0 if side is an inner side, or if the side is on a tree side,
! assign the tree BC index (can be 0 for inner tree sides) 
!===================================================================================================================================
! MODULES
USE MODH_p4est_vars,ONLY:QuadCoords,QuadLevel,TreeToBC,IntSize
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN) :: iTree 
INTEGER,INTENT(IN) :: iElem 
INTEGER,INTENT(IN) :: Pside   !p4est local side ID 0...5
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
INTEGER :: GET_BCINDEX_FROM_TREE
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER :: pos
!-----------------------------------------------------------------------------------------------------------------------------------
GET_BCINDEX_FROM_TREE=0
SELECT CASE(PSide)
CASE(0,2,4) !0: x(1)=0?, 2: x(2)=0?, 4: x(3)=0?
  pos=QuadCoords((PSide+2)/2,iElem)
  IF(pos.EQ.0)THEN
    GET_BCINDEX_FROM_TREE=TreeToBC(Pside+1,iTree)
  END IF
CASE(1,3,5) !1: x(1)+2^l=IntSize?, 3: x(2)+2^l=IntSize?, 5: x(3)+2^l=IntSize?
  pos=QuadCoords((PSide+1)/2,iElem) + IntSize/(2**QuadLevel(iElem))
  IF(pos.EQ.IntSize)THEN
    GET_BCINDEX_FROM_TREE=TreeToBC(Pside+1,iTree)
  END IF
END SELECT
END FUNCTION GET_BCINDEX_FROM_TREE


FUNCTION GetHFlip(PSide0,PSide1,PFlip)
!===================================================================================================================================
! transform an p4est orientation (r = PFlip)  in HOPR Flip, using local Side and neighbor side  
!===================================================================================================================================
! MODULES
USE MODH_P4EST_Vars,ONLY:P2H_FaceMap,H2P_FaceNodeMap,P2H_FaceNodeMap,P4R,P4Q,P4P
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)   :: PSide0  ! local side ID of HOPEST
INTEGER,INTENT(IN)   :: PSide1  ! P4EST neighbour local side id
INTEGER,INTENT(IN)   :: PFlip   ! Neighbour side in p4est convention: 0..5
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
INTEGER              :: GetHFlip   ! Neighbour side in p4est convention: 0..5
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER              :: HSide0,PNode0,PNode1 ! p4est Side,Flip,Mortar encoding
!-----------------------------------------------------------------------------------------------------------------------------------
!1. Get CGNS side from P4 side
HSide0=P2H_FaceMap(PSide0)
!2. First node CGNS -> p4est
PNode0=H2P_FaceNodeMap(1,HSide0)
!3. Get oriented node on neighbour side, formula and matrices see paper Burstedde p4est, 2011
PNode1=P4P(P4Q(P4R(PSide0,PSide1),PFlip),PNode0)
!4. P4EST node -> CGNS
GetHFlip=P2H_FaceNodeMap(PNode1,PSide1)

END FUNCTION GetHFlip


FUNCTION GetHMortar(PMortar,PSide,PnbSide,PFlip)
!===================================================================================================================================
! Transform a p4est mortar ID to a HOPR mortar ID 
!===================================================================================================================================
! MODULES
USE MODH_P4EST_Vars,ONLY:P2H_FaceMap,H2P_FaceNodeMap,P2H_FaceNodeMap,P4R,P4Q,P4P
USE MODH_P4EST_Vars,ONLY:H_MortarCase,P2H_MortarMap
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)   :: PMortar ! Index of mortar side: 0...3
INTEGER,INTENT(IN)   :: PSide   ! local side ID of HOPEST
INTEGER,INTENT(IN)   :: PnbSide ! P4EST neighbour local side id
INTEGER,INTENT(IN)   :: PFlip   ! Neighbour side in p4est convention: 0..5
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
INTEGER              :: GetHMortar   ! Index of mortar side: 1..4
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER              :: PNodeA,PNodeB,HNode1,HNode2
!-----------------------------------------------------------------------------------------------------------------------------------

! Side 1: neighbor
!1. Get node a and b from node 0 and 1 of neighbor side P
PNodeA=P4P(P4Q(P4R(PnbSide,PSide),PFlip),0)
PNodeB=P4P(P4Q(P4R(PnbSide,PSide),PFlip),1)

!2. Get CGNS node 1 and 2 from p4est node a and b 
HNode1=P2H_FaceNodeMap(PNodeA,PSide)
HNode2=P2H_FaceNodeMap(PNodeB,PSide)

!3. 8 possible combinations (MortarCase),using node 1 and 2 and map PMortar index to Hmortar index
GetHMortar = P2H_MortarMap(PMortar, H_MortarCase(HNode1,HNode2)) 

END FUNCTION GetHMortar

SUBROUTINE BuildBCs()
!===================================================================================================================================
! Subroutine to translate p4est mesh datastructure to HOPR datastructure
!===================================================================================================================================
! MODULES
USE, INTRINSIC :: ISO_C_BINDING
USE MODH_Globals
USE MODH_P4EST_Vars,   ONLY: H2P_FaceMap,p4est
USE MODH_Mesh_Vars,    ONLY: tElem,tSide,nTrees,Trees
USE MODH_P4EST_Binding,ONLY: p4_build_bcs
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
TYPE(tElem),POINTER         :: Elem
TYPE(tSide),POINTER         :: Side
INTEGER                     :: iTree,iSide
INTEGER(KIND=C_INT32_T)     :: TreeToBC(0:5,nTrees)
!===================================================================================================================================
TreeToBC=-1
DO iTree=1,nTrees
  Elem=>Trees(iTree)%ep
  DO iSide=1,6
    Side=>Elem%side(iSide)%sp
    TreeToBC(H2P_FaceMap(iSide),iTree)=Side%BCIndex
  END DO
END DO
CALL p4_build_bcs(p4est,nTrees,TreeToBC)
END SUBROUTINE BuildBCs


SUBROUTINE buildHOp4GeometryX(a,b,c,x,y,z,tree)
!===================================================================================================================================
! Subroutine to translate p4est mesh datastructure to HOPR datastructure
!===================================================================================================================================
! MODULES
USE, INTRINSIC :: ISO_C_BINDING
USE MODH_Globals
USE MODH_Basis,        ONLY: LagrangeInterpolationPolys
USE MODH_Mesh_Vars,    ONLY: XGeo,xi_Ngeo,wBary_Ngeo,NGeo
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL( KIND = C_DOUBLE ),INTENT(IN),VALUE    :: a,b,c
P4EST_F90_TOPIDX,INTENT(IN),VALUE    :: tree
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL( KIND = C_DOUBLE ),INTENT(OUT)         :: x,y,z
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL         :: Vdm_xi(0:NGeo),Vdm_eta(0:NGeo),Vdm_zeta(0:NGeo),Vdm_eta_zeta
REAL         :: xi(3),HOabc(3)
INTEGER      :: i,j,k
!-----------------------------------------------------------------------------------------------------------------------------------

xi(1)=-1.+2*a
xi(2)=-1.+2*b
xi(3)=-1.+2*c
CALL LagrangeInterpolationPolys(xi(1),Ngeo,xi_Ngeo,wBary_Ngeo,Vdm_xi(:))
CALL LagrangeInterpolationPolys(xi(2),Ngeo,xi_Ngeo,wBary_Ngeo,Vdm_eta(:))
CALL LagrangeInterpolationPolys(xi(3),Ngeo,xi_Ngeo,wBary_Ngeo,Vdm_zeta(:))
HOabc(:)=0.
DO k=0,NGeo
  DO j=0,NGeo
    Vdm_eta_zeta=Vdm_eta(j)*Vdm_zeta(k)
    DO i=0,NGeo
      HOabc(:)=HOabc(:)+XGeo(:,i,j,k,tree)*Vdm_xi(i)*Vdm_eta_zeta
    END DO
  END DO
END DO
x=HOabc(1)
y=HOabc(2)
z=HOabc(3)

END SUBROUTINE buildHOp4GeometryX

SUBROUTINE testHOabc()
!===================================================================================================================================
! Subroutine to translate p4est mesh datastructure to HOPR datastructure
!===================================================================================================================================
! MODULES
USE, INTRINSIC :: ISO_C_BINDING
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL        :: a,b,c,x,y,z
P4EST_F90_TOPIDX      :: tree
!-----------------------------------------------------------------------------------------------------------------------------------
a=1
b=0
c=0
tree=1

CALL buildHOp4GeometryX(a,b,c,x,y,z,tree)

WRITE(*,*) x,y,z

END SUBROUTINE testHOabc


SUBROUTINE FinalizeP4EST()
!===================================================================================================================================
! Subroutine to translate p4est mesh datastructure to HOPR datastructure
!===================================================================================================================================
! MODULES
USE, INTRINSIC :: ISO_C_BINDING
USE MODH_P4EST_Vars,   ONLY: QuadCoords,QuadLevel,p4est,mesh,connectivity
USE MODH_Mesh_Vars,    ONLY: nElems,Elems
USE MODH_P4EST_Binding
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER           :: iElem,iLocSide
!===================================================================================================================================
! TODO: DEALLOCATE P4EST BC POINTER ARRAY
DO iElem=1,nElems
  DO iLocSide=1,6
    DEALLOCATE(Elems(iElem)%ep%Side(iLocSide)%sp)
  END DO
  DEALLOCATE(Elems(iElem)%ep)
END DO
DEALLOCATE(Elems)
! TODO: DEALLOCATE P4EST / CONNECTIVITY THEMSELVES
CALL p4_destroy_p4est(p4est)
CALL p4_destroy_mesh(mesh)
!CALL p4_destroy_connectivity(connectivity)
SDEALLOCATE(QuadCoords)
SDEALLOCATE(QuadLevel)
! TO NOT TOUCH QuadToTree/Quad/Face/Half -> belongs to p4est

END SUBROUTINE FinalizeP4EST


END MODULE MODH_P4EST
