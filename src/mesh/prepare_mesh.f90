#include "hopest_f.h"

MODULE MODH_Prepare_Mesh
!===================================================================================================================================
! Contains subroutines to build (curviilinear) meshes and provide metrics, etc.
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES (PUBLIC)
!-----------------------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------

INTERFACE setLocalSideIDs
  MODULE PROCEDURE setLocalSideIDs
END INTERFACE

INTERFACE fillMeshInfo
  MODULE PROCEDURE fillMeshInfo
END INTERFACE

PUBLIC::setLocalSideIDs,fillMeshInfo

#ifdef MPI
INTERFACE exchangeFlip
  MODULE PROCEDURE exchangeFlip
END INTERFACE

PUBLIC::exchangeFlip 
#endif
!===================================================================================================================================

CONTAINS


SUBROUTINE setLocalSideIDs()
!===================================================================================================================================
! 
!===================================================================================================================================
! MODULES
USE MODH_Globals
USE MODH_Mesh_Vars,  ONLY: tElem,tSide
USE MODH_Mesh_Vars,  ONLY: nElems,nInnerSides,nSides,nBCSides,offsetElem
USE MODH_Mesh_Vars,  ONLY: Elems,nMPISides_MINE,nMPISides_YOUR,BoundaryType,nBCs
USE MODH_Mesh_Vars,  ONLY: nMortarSides 
#ifdef MPI
USE MODH_ReadInTools,ONLY: GETLOGICAL
USE MODH_MPI_Vars,   ONLY: nNbProcs,NbProc,nMPISides_Proc,nMPISides_MINE_Proc,nMPISides_YOUR_Proc
USE MODH_MPI_Vars,   ONLY: offsetElemMPI,offsetMPISides_MINE,offsetMPISides_YOUR
USE MODH_Mesh_ReadIn,ONLY: Qsort1Int,INVMAP
#endif
IMPLICIT NONE
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER   :: iElem,FirstElemInd,LastElemInd
INTEGER   :: iLocSide,iSide,iInnerSide,iBCSide
INTEGER   :: iMortar,iMortarSide,nMortars
INTEGER   :: i,j
INTEGER   :: PeriodicBCMap(nBCs)       !connected periodic BCs
TYPE(tElem),POINTER :: aElem
TYPE(tSide),POINTER :: aSide
#ifdef MPI
INTEGER               :: iNbProc,ioUnit
INTEGER               :: ProcInfo(4),nNBmax      !for output only
INTEGER,ALLOCATABLE   :: SideIDMap(:)
INTEGER,ALLOCATABLE   :: NBinfo(:,:),NBinfo_glob(:,:,:),nNBProcs_glob(:),Procinfo_glob(:,:),tmparray(:,:)  !for output only
REAL,ALLOCATABLE      :: tmpreal(:,:)
CHARACTER(LEN=10)     :: formatstr
LOGICAL               :: writePartitionInfo
#endif
!===================================================================================================================================
!FirstElemInd= offsetElem+1
!LastElemInd = offsetElem+nElems
FirstElemInd= 1
LastElemInd = nElems
! ----------------------------------------
! Set side IDs to arrange sides:
! 1. BC sides
! 2. inner sides
! 3. MPI sides
! MPI Sides are not included here!
! ----------------------------------------
! Get connection between periodic BCs
PeriodicBCMap=-2
DO i=1,nBCs
  IF((BoundaryType(i,BC_TYPE).NE.1)) PeriodicBCMap(i)=-1 ! not periodic
  IF((BoundaryType(i,BC_TYPE).EQ.1).AND.(BoundaryType(i,BC_ALPHA).GT.0)) PeriodicBCMap(i)=-1 ! slave
  IF((BoundaryType(i,BC_TYPE).EQ.1).AND.(BoundaryType(i,BC_ALPHA).LT.0))THEN
    DO j=1,nBCs
      IF(BoundaryType(j,BC_TYPE).NE.1) CYCLE
      IF(BoundaryType(j,BC_ALPHA).EQ.(-BoundaryType(i,BC_ALPHA))) PeriodicBCMap(i)=j
    END DO
  END IF
END DO
IF(ANY(PeriodicBCMap.EQ.-2))&
  CALL abort(__STAMP__,'Periodic connection not found.')

DO iElem=FirstElemInd,LastElemInd
  aElem=>Elems(iElem)%ep
  DO iLocSide=1,6
    aSide=>aElem%Side(iLocSide)%sp
    nMortars=aSide%nMortars 
    DO iMortar=0,nMortars
      IF(iMortar.GT.0) aSide=>aElem%Side(iLocSide)%sp%mortarSide(iMortar)%sp

      aSide%sideID=-1
      aSide%tmp=0
      ! periodics have two bcs: set to (positive) master bc (e.g. from -1 to 1)
      IF(aSide%BCIndex.GE.1)THEN
        IF(PeriodicBCMap(aSide%BCIndex).NE.-1)&
          aSide%BCIndex=PeriodicBCMap(aSide%BCIndex)
      END IF
    END DO !iMortar
  END DO
END DO


iSide=0
iBCSide=0
iMortarSide=nBCSides
iInnerSide=nBCSides+nMortarSides
DO iElem=FirstElemInd,LastElemInd
  aElem=>Elems(iElem)%ep
  DO iLocSide=1,6
    aSide=>aElem%Side(iLocSide)%sp
    nMortars=aSide%nMortars 
    DO iMortar=0,nMortars
      IF(iMortar.GT.0) aSide=>aElem%Side(iLocSide)%sp%mortarSide(iMortar)%sp

      IF(aSide%sideID.EQ.-1)THEN
        IF(aSide%NbProc.EQ.-1)THEN ! no MPI Sides
          IF(ASSOCIATED(aSide%connection))THEN
            iInnerSide=iInnerSide+1
            iSide=iSide+1
            aSide%SideID=iInnerSide
            IF(aSide%connection%SideID.NE.-1) STOP 'ERROR: SideID of connection already set!'
            aSide%connection%SideID=iInnerSide
          ELSE
            IF(aSide%MortarType.GT.0) THEN
              iMortarSide=iMortarSide+1
              iSide=iSide+1
              aSide%SideID=iMortarSide
            ELSE !this is now a BC side, really!
              iBCSide=iBCSide+1
              iSide=iSide+1
              aSide%SideID=iBCSide
            END IF !mortar
          END IF !associated connection
        END IF ! .NOT. MPISide
      END IF !sideID NE -1
    END DO !iMortar
  END DO ! iLocSide=1,6
END DO !iElem
IF(iSide.NE.nInnerSides+nBCSides+nMortarSides) STOP'not all SideIDs are set!'

nMPISides_MINE=0
nMPISides_YOUR=0
#ifdef MPI
STOP 'no mpi yet'
! SPLITTING MPISides in MINE and YOURS
ALLOCATE(nMPISides_MINE_Proc(1:nNbProcs),nMPISides_YOUR_Proc(1:nNbProcs))
nMPISides_MINE_Proc=0
nMPISides_YOUR_Proc=0
DO iNbProc=1,nNbProcs
  IF(myRank.LT.NbProc(iNbProc)) THEN
    nMPISides_MINE_Proc(iNbProc)=nMPISides_Proc(iNbProc)/2
  ELSE
    nMPISides_MINE_Proc(iNbProc)=nMPISides_Proc(iNbProc)-nMPISides_Proc(iNbProc)/2
  END IF    
  nMPISides_YOUR_Proc(iNbProc)=nMPISides_Proc(iNbProc)-nMPISides_MINE_Proc(iNbProc)
END DO
nMPISides_MINE=SUM(nMPISides_MINE_Proc)
nMPISides_YOUR=SUM(nMPISides_YOUR_Proc)

ALLOCATE(offsetMPISides_YOUR(0:nNbProcs),offsetMPISides_MINE(0:nNbProcs))
offsetMPISides_MINE=0
offsetMPISides_YOUR=0
! compute offset, first all MINE , then all YOUR MPISides
offsetMPISides_MINE(0)=nInnerSides+nBCSides+nMortarSides
DO iNbProc=1,nNbProcs
  offsetMPISides_MINE(iNbProc)=offsetMPISides_MINE(iNbProc-1)+nMPISides_MINE_Proc(iNbProc)
END DO
offsetMPISides_YOUR(0)=offsetMPISides_MINE(nNbProcs)
DO iNbProc=1,nNbProcs
  offsetMPISides_YOUR(iNbProc)=offsetMPISides_YOUR(iNbProc-1)+nMPISides_YOUR_Proc(iNbProc)
END DO
IF(nProcessors.EQ.1) RETURN
DO iNbProc=1,nNbProcs
  ALLOCATE(SideIDMap(nMPISides_Proc(iNbProc)))
  iSide=0
  DO iElem=FirstElemInd,LastElemInd
    aElem=>Elems(iElem)%ep
    DO iLocSide=1,6
      aSide=>aElem%Side(iLocSide)%sp
      nMortars=aSide%nMortars 
      DO iMortar=0,nMortars
        IF(iMortar.GT.0) aSide=>aElem%Side(iLocSide)%sp%mortarSide(iMortar)%sp
        IF(aSide%NbProc.NE.NbProc(iNbProc))CYCLE
        iSide=iSide+1
        SideIDMap(iSide)=aSide%ind !global Side Index 
      END DO !iMortar
    END DO !iLocSide
  END DO !iElem
  CALL Qsort1Int(SideIDMap) !sort by global side index
  DO iElem=FirstElemInd,LastElemInd
    aElem=>Elems(iElem)%ep
    DO iLocSide=1,6
      aSide=>aElem%Side(iLocSide)%sp
      nMortars=aSide%nMortars 
      DO iMortar=0,nMortars
        IF(iMortar.GT.0) aSide=>aElem%Side(iLocSide)%sp%mortarSide(iMortar)%sp
        IF(aSide%NbProc.NE.NbProc(iNbProc))CYCLE
        aSide%SideID=INVMAP(aSide%ind,nMPISides_Proc(iNbProc),SideIDMap) ! get sorted iSide
        IF(myRank.LT.aSide%NbProc)THEN
          IF(aSide%SideID.LE.nMPISides_MINE_Proc(iNbProc))THEN !MINE
            aSide%SideID=aSide%SideID +offsetMPISides_MINE(iNbProc-1)
          ELSE !YOUR
            aSide%SideID=(aSide%SideID-nMPISides_MINE_Proc(iNbProc))+offsetMPISides_YOUR(iNbProc-1)
          END IF
        ELSE
          IF(aSide%SideID.LE.nMPISides_YOUR_Proc(iNbProc))THEN !MINE
            aSide%SideID=aSide%SideID +offsetMPISides_YOUR(iNbProc-1)
          ELSE !YOUR
            aSide%SideID=(aSide%SideID-nMPISides_YOUR_Proc(iNbProc))+offsetMPISides_MINE(iNbProc-1)
          END IF
        END IF !myrank<NbProc
      END DO !iMortar
    END DO !iLocSide
  END DO !iElem
  DEALLOCATE(SideIDMap)
END DO !nbProc(i)

WRITE(formatstr,'(a5,I2,a3)')'(A22,',nNBProcs,'I8)'
LOGWRITE(*,*)'-------------------------------------------------------'
LOGWRITE(*,'(A22,I8)')'nNbProcs:',nNbProcs
LOGWRITE(*,*)'-------------------------------------------------------'
LOGWRITE(*,formatstr)'NbProc:'   ,NbProc
LOGWRITE(*,*)'-------------------------------------------------------'
LOGWRITE(*,formatstr)'nMPISides_Proc:',nMPISides_Proc
LOGWRITE(*,*)'-------------------------------------------------------'
LOGWRITE(*,formatstr)'nMPISides_MINE_Proc:',nMPISides_MINE_Proc
LOGWRITE(*,formatstr)'nMPISides_YOUR_Proc:',nMPISides_YOUR_Proc
WRITE(formatstr,'(a5,I2,a3)')'(A22,',nNBProcs+1,'I8)'
LOGWRITE(*,*)'-------------------------------------------------------'
LOGWRITE(*,formatstr)'offsetMPISides_MINE:',offsetMPISides_MINE
LOGWRITE(*,formatstr)'offsetMPISides_YOUR:',offsetMPISides_YOUR
LOGWRITE(*,*)'-------------------------------------------------------'

writePartitionInfo = GETLOGICAL('writePartitionInfo','.FALSE.')
IF(.NOT.writePartitionInfo) RETURN
!output partitioning info
ProcInfo(1)=nElems
ProcInfo(2)=nSides
ProcInfo(3)=nInnerSides
ProcInfo(4)=nBCSides
IF(MPIroot)THEN
  ALLOCATE(nNBProcs_glob(0:nProcessors-1))
  ALLOCATE(ProcInfo_glob(4,0:nProcessors-1))
  nNBProcs_glob=-99999
  Procinfo_glob=-88888
ELSE
  ALLOCATE(nNBProcs_glob(1)) !dummy for debug
  ALLOCATE(ProcInfo_glob(1,1)) !dummy for debug
END IF !MPIroot 
CALL MPI_GATHER(nNBProcs,1,MPI_INTEGER,nNBProcs_glob,1,MPI_INTEGER,0,MPI_COMM_WORLD,iError)
CALL MPI_GATHER(ProcInfo,4,MPI_INTEGER,ProcInfo_glob,4,MPI_INTEGER,0,MPI_COMM_WORLD,iError)
IF(MPIroot)THEN
  nNBmax=MAXVAL(nNBProcs_glob) !count, total number of columns in table
  ALLOCATE(NBinfo_glob(6,nNBmax,0:nProcessors))
  NBinfo_glob=-77777
ELSE
  ALLOCATE(NBinfo_glob(1,1,1)) !dummy for debug
END IF
CALL MPI_BCAST(nNBmax,1,MPI_INTEGER,0,MPI_COMM_WORLD,iError) 
ALLOCATE(NBinfo(6,nNbmax))
NBinfo=0
NBinfo(1,1:nNBProcs)=NBProc
NBinfo(2,1:nNBProcs)=nMPISides_Proc
NBinfo(3,1:nNBProcs)=nMPISides_MINE_Proc
NBinfo(4,1:nNBProcs)=nMPISides_YOUR_Proc
NBinfo(5,1:nNBProcs)=offsetMPISides_MINE(0:nNBProcs-1)
NBinfo(6,1:nNBProcs)=offsetMPISides_YOUR(0:nNBProcs-1)
CALL MPI_GATHER(NBinfo,6*nNBmax,MPI_INTEGER,NBinfo_glob,6*nNBmax,MPI_INTEGER,0,MPI_COMM_WORLD,iError)
DEALLOCATE(NBinfo)
IF(MPIroot)THEN
  ioUnit=GETFREEUNIT()
  OPEN(UNIT=ioUnit,FILE='partitionInfo.out',STATUS='REPLACE')
  WRITE(ioUnit,*)'Partition Information:'
  WRITE(ioUnit,*)'total number of Procs,',nProcessors
  WRITE(ioUnit,*)'total number of Elems,',SUM(Procinfo_glob(1,:))

  WRITE(ioUnit,'(8(A15))')'Rank','nElems','nSides','nInnerSides','nBCSides','nMPISides','nMPISides_MINE','nNBProcs'
  WRITE(ioUnit,'(A120)')&
      '======================================================================================================================='
  !statistics
  ALLOCATE(tmparray(7,0:3),tmpreal(7,2))
  tmparray(:,0)=0      !tmp
  tmparray(:,1)=0      !mean
  tmparray(:,2)=-1E08  !max
  tmparray(:,3)=1E08   !min
  DO i=0,nProcessors-1
    !actual proc
    tmparray(1,0)=Procinfo_glob(1,i)
    tmparray(2,0)=Procinfo_glob(2,i)
    tmparray(3,0)=Procinfo_glob(3,i)
    tmparray(4,0)=Procinfo_glob(4,i)
    tmparray(5,0)=SUM(NBinfo_glob(2,:,i))
    tmparray(6,0)=SUM(NBinfo_glob(3,:,i))
    tmparray(7,0)=nNBProcs_glob(i)
    DO j=1,7
      !mean
      tmparray(j,1)=tmparray(j,1)+tmparray(j,0)
      !max
      tmparray(j,2)=MAX(tmparray(j,2),tmparray(j,0))
      tmparray(j,3)=MIN(tmparray(j,3),tmparray(j,0))
    END DO
  END DO
  tmpreal(:,1)=REAL(tmparray(:,1))/REAL(nProcessors) !mean in REAL
  tmpreal(:,2)=0.   !RMS
  DO i=0,nProcessors-1
    !actual proc
    tmparray(1,0)=Procinfo_glob(1,i)
    tmparray(2,0)=Procinfo_glob(2,i)
    tmparray(3,0)=Procinfo_glob(3,i)
    tmparray(4,0)=Procinfo_glob(4,i)
    tmparray(5,0)=SUM(NBinfo_glob(2,:,i))
    tmparray(6,0)=SUM(NBinfo_glob(3,:,i))
    tmparray(7,0)=nNBProcs_glob(i)
    DO j=1,7
      tmpreal(j,2)=tmpreal(j,2)+(tmparray(j,0)-tmpreal(j,1))**2 
    END DO
  END DO
  tmpreal(:,2)=SQRT(tmpreal(:,2)/REAL(nProcessors))
  WRITE(ioUnit,'(A15,7(5X,F10.2))')'   MEAN        ',tmpreal(:,1)
  WRITE(ioUnit,'(A120)')&
      '-----------------------------------------------------------------------------------------------------------------------'
  WRITE(ioUnit,'(A15,7(5X,F10.2))')'   RMS         ',tmpreal(:,2)
  WRITE(ioUnit,'(A120)')&
      '-----------------------------------------------------------------------------------------------------------------------'
  WRITE(ioUnit,'(A15,7(5X,I10))')'   MIN         ',tmparray(:,3)
  WRITE(ioUnit,'(A120)')&
      '-----------------------------------------------------------------------------------------------------------------------'
  WRITE(ioUnit,'(A15,7(5X,I10))')'   MAX         ',tmparray(:,2)
  WRITE(ioUnit,'(A120)')&
      '======================================================================================================================='
  DO i=0,nProcessors-1
    WRITE(ioUnit,'(8(5X,I10))')i,Procinfo_glob(:,i),SUM(NBinfo_glob(2,:,i)),SUM(NBinfo_glob(3,:,i)),nNBProcs_glob(i)
    WRITE(ioUnit,'(A120)')&
      '-----------------------------------------------------------------------------------------------------------------------'
  END DO
  WRITE(ioUnit,*)' '
  WRITE(ioUnit,*)'Information per neighbor processor'
  WRITE(ioUnit,*)' '
  WRITE(ioUnit,'(7(A15))')'Rank','NBProc','nMPISides_Proc','nMPISides_MINE','nMPISides_YOUR','offset_MINE','offset_YOUR'
  WRITE(ioUnit,'(A120)')&
      '======================================================================================================================='
  DO i=0,nProcessors-1
    WRITE(ioUnit,'(7(5X,I10))')i,NBinfo_glob(:,1,i)
    DO j=2,nNBProcs_glob(i)
      WRITE(ioUnit,'(A15,6(5X,I10))')' ',NBinfo_glob(:,j,i)
    END DO
    WRITE(ioUnit,'(A120)')&
      '-----------------------------------------------------------------------------------------------------------------------'
  END DO
  DEALLOCATE(tmparray,tmpreal)
  CLOSE(ioUnit) 
END IF !MPIroot
DEALLOCATE(NBinfo_glob,nNBProcs_glob,ProcInfo_glob)
#endif /*MPI*/  
END SUBROUTINE setLocalSideIDs



SUBROUTINE fillMeshInfo()
!===================================================================================================================================
! 
!===================================================================================================================================
! MODULES
USE MODH_Globals
USE MODH_Mesh_Vars,ONLY:tElem,tSide,BoundaryType
USE MODH_Mesh_Vars,ONLY:nElems,offsetElem,nSides,nInnerSides,nBCSides,nMPISides,nMortarSides
USE MODH_Mesh_Vars,ONLY:nMPISides_MINE
USE MODH_Mesh_Vars,ONLY:ElemToSide,SideToElem,BC,AnalyzeSide
USE MODH_Mesh_Vars,ONLY:Elems
USE MODH_Mesh_Vars,ONLY:MortarType,Mortar_nbSideID,Mortar_Flip
#ifdef MPI
USE MODH_MPI_vars
#endif
IMPLICIT NONE
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER             :: iElem,iSide,LocSideID,nSides_flip(0:4)
INTEGER             :: nSides_MortarType(1:3)
INTEGER             :: iMortar,nMortars
TYPE(tElem),POINTER :: aElem
TYPE(tSide),POINTER :: aSide
!===================================================================================================================================
! ELement to Side mapping
nSides_flip=0
DO iElem=1,nElems
  aElem=>Elems(iElem+offsetElem)%ep
  DO LocSideID=1,6
    aSide=>aElem%Side(LocSideID)%sp
    ElemToSide(E2S_SIDE_ID,LocSideID,iElem)=aSide%SideID
    ElemToSide(E2S_FLIP,LocSideID,iElem)   =aSide%Flip
    nSides_flip(aSide%flip)=nSides_flip(aSide%flip)+1
  END DO ! LocSideID
END DO ! iElem

! Side to Element mapping, sorted by SideID
SideToElem=0
DO iElem=1,nElems
  aElem=>Elems(iElem+offsetElem)%ep
  DO LocSideID=1,6
    aSide=>aElem%Side(LocSideID)%sp
    IF(aSide%Flip.EQ.0)THEN !root side
      SideToElem(S2E_ELEM_ID,aSide%SideID)         = iElem !root Element
      SideToElem(S2E_LOC_SIDE_ID,aSide%SideID)     = LocSideID
    ELSE
      SideToElem(S2E_NB_ELEM_ID,aSide%SideID)      = iElem ! element with flipped side
      SideToElem(S2E_NB_LOC_SIDE_ID,aSide%SideID)  = LocSideID
      SideToElem(S2E_FLIP,aSide%SideID)            = aSide%Flip
    END IF
    IF(aSide%sideID .LE. nBCSides)THEN
      BC(aSide%sideID)=aSide%BCIndex
    ENDIF
    AnalyzeSide(aSide%sideID)=aSide%BCIndex
  END DO ! LocSideID
END DO ! iElem


! Mapping of Mortar Master Side to Mortar Slave Side
nSides_MortarType=0
DO iElem=1,nElems
  aElem=>Elems(iElem+offsetElem)%ep
  DO LocSideID=1,6
    aSide=>aElem%Side(LocSideID)%sp
    IF(aSide%nMortars.GT.0)THEN !mortar side
      MortarType(aSide%SideID)=aSide%MortarType
      DO iMortar=1,aSide%nMortars
        Mortar_nbSideID(iMortar,aSide%SideID)=aSide%MortarSide(iMortar)%sp%SideID
        Mortar_flip(iMortar,aSide%SideID)=aSide%MortarSide(iMortar)%sp%Flip
      END DO !iMortar
      nSides_MortarType(aSide%MortarType)=nSides_MortarType(aSide%MortarType)+1
    END IF !mortarSide
  END DO ! LocSideID
END DO ! iElem

#ifdef MPI
IF(MPIroot)THEN
  CALL MPI_REDUCE(MPI_IN_PLACE,nSides_flip,5,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,iError)
ELSE
  CALL MPI_REDUCE(nSides_flip,nSides_flip,5,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,iError)
END IF
IF(MPIroot)THEN
  CALL MPI_REDUCE(MPI_IN_PLACE,nSides_MortarType,3,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,iError)
ELSE
  CALL MPI_REDUCE(nSides_MortarType,nSides_MortarType,3,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,iError)
END IF
#endif /*MPI*/
SWRITE(UNIT_StdOut,'(132("."))')
SWRITE(*,'(A,A34,I0)')' |','nSides with Flip=0     | ',nSides_flip(0)
SWRITE(*,'(A,A34,I0)')' |','nSides with Flip=1     | ',nSides_flip(1)
SWRITE(*,'(A,A34,I0)')' |','nSides with Flip=2     | ',nSides_flip(2)
SWRITE(*,'(A,A34,I0)')' |','nSides with Flip=3     | ',nSides_flip(3)
SWRITE(*,'(A,A34,I0)')' |','nSides with Flip=4     | ',nSides_flip(4)
SWRITE(UNIT_StdOut,'(132("."))')
SWRITE(*,'(A,A34,I0)')' |','nSides of MortarType=1 | ',nSides_MortarType(1)
SWRITE(*,'(A,A34,I0)')' |','nSides of MortarType=2 | ',nSides_MortarType(2)
SWRITE(*,'(A,A34,I0)')' |','nSides of MortarType=3 | ',nSides_MortarType(3)
SWRITE(UNIT_StdOut,'(132("."))')
END SUBROUTINE fillMeshInfo


#ifdef MPI
SUBROUTINE exchangeFlip()
!===================================================================================================================================
! set flip of MINE sides to zero, therefore send flip of MINE to other processor, so that YOUR sides get their corresponding flip>0
!===================================================================================================================================
! MODULES
USE MODH_Globals
USE MODH_Mesh_Vars,  ONLY: tElem,tSide
USE MODH_Mesh_Vars,ONLY:nElems,offsetElem
USE MODH_Mesh_Vars,ONLY:Elems
USE MODH_MPI_vars
IMPLICIT NONE
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER             :: iElem,LocSideID
INTEGER             :: iMortar,nMortars
INTEGER             :: Flip_MINE(offsetMPISides_MINE(0)+1:offsetMPISides_MINE(nNBProcs))
INTEGER             :: Flip_YOUR(offsetMPISides_YOUR(0)+1:offsetMPISides_YOUR(nNBProcs))
INTEGER             :: SendRequest(nNbProcs),RecRequest(nNbProcs)
TYPE(tElem),POINTER :: aElem
TYPE(tSide),POINTER :: aSide
!===================================================================================================================================
IF(nProcessors.EQ.1) RETURN
!fill MINE flip info
DO iElem=1,nElems
  aElem=>Elems(iElem+offsetElem)%ep
  DO LocSideID=1,6
    aSide=>aElem%Side(LocSideID)%sp
    nMortars=aSide%nMortars 
    DO iMortar=0,nMortars
      IF(iMortar.GT.0) aSide=>aElem%Side(LocSideID)%sp%mortarSide(iMortar)%sp
      IF((aSide%SideID.GT.offsetMPISides_MINE(0)       ).AND.&
         (aSide%SideID.LE.offsetMPISides_MINE(nNBProcs)))THEN
        Flip_MINE(aSide%sideID)=aSide%flip
      END IF
    END DO ! iMortar
  END DO ! LocSideID
END DO ! iElem
DO iNbProc=1,nNbProcs
  ! Start send flip from MINE 
  IF(nMPISides_MINE_Proc(iNbProc).GT.0)THEN
    nSendVal    =nMPISides_MINE_Proc(iNbProc)
    SideID_start=OffsetMPISides_MINE(iNbProc-1)+1
    SideID_end  =OffsetMPISides_MINE(iNbProc)
    CALL MPI_ISEND(Flip_MINE(SideID_start:SideID_end),nSendVal,MPI_INTEGER,  &
                    nbProc(iNbProc),0,MPI_COMM_WORLD,SendRequest(iNbProc),iError)
  END IF
  ! Start receive flip to YOUR
  IF(nMPISides_YOUR_Proc(iNbProc).GT.0)THEN
    nRecVal     =nMPISides_YOUR_Proc(iNbProc)
    SideID_start=OffsetMPISides_YOUR(iNbProc-1)+1
    SideID_end  =OffsetMPISides_YOUR(iNbProc)
    CALL MPI_IRECV(Flip_YOUR(SideID_start:SideID_end),nRecVal,MPI_INTEGER,  &
                    nbProc(iNbProc),0,MPI_COMM_WORLD,RecRequest(iNbProc),iError)
  END IF
END DO !iProc=1,nNBProcs
DO iNbProc=1,nNbProcs
  IF(nMPISides_YOUR_Proc(iNbProc).GT.0)CALL MPI_WAIT(RecRequest(iNbProc) ,MPIStatus,iError)
  IF(nMPISides_MINE_Proc(iNBProc).GT.0)CALL MPI_WAIT(SendRequest(iNbProc),MPIStatus,iError)
END DO !iProc=1,nNBProcs
DO iElem=1,nElems
  aElem=>Elems(iElem+offsetElem)%ep
  DO LocSideID=1,6
    aSide=>aElem%Side(LocSideID)%sp
    nMortars=aSide%nMortars 
    DO iMortar=0,nMortars
      IF(iMortar.GT.0) aSide=>aElem%Side(LocSideID)%sp%mortarSide(iMortar)%sp
      IF(aSide%NbProc.EQ.-1) CYCLE !no MPISide
      IF(aSide%SideID.GT.offsetMPISides_YOUR(0))THEN
        IF(aSide%flip.EQ.0)THEN
          IF(Flip_YOUR(aSide%SideID).EQ.0) STOP 'problem in exchangeflip'
          aSide%flip=Flip_YOUR(aSide%sideID)
        END IF
      ELSE
        aSide%flip=0 !MINE MPISides flip=0
      END IF
    END DO ! iMortar
  END DO ! LocSideID
END DO ! iElem
  
END SUBROUTINE exchangeFlip
#endif


END MODULE MODH_Prepare_Mesh
