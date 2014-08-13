#include "hopest_f.h"

MODULE MOD_Analyze
!===================================================================================================================================
! Add comments please!
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES 
!-----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------
INTERFACE InitAnalyze
  MODULE PROCEDURE InitAnalyze
END INTERFACE

INTERFACE Analyze
  MODULE PROCEDURE Analyze
END INTERFACE

PUBLIC::InitAnalyze
PUBLIC::Analyze
!===================================================================================================================================

CONTAINS

SUBROUTINE InitAnalyze()
!===================================================================================================================================
! Basic Analyze initialization. 
!===================================================================================================================================
! MODULES
USE MOD_Mesh_Vars,ONLY:Ngeo_out,Xi_Ngeo_out
USE MOD_Analyze_Vars
USE MOD_ReadInTools,ONLY:GETINT,GETLOGICAL
USE MOD_Basis,    ONLY: BarycentricWeights,InitializeVandermonde,PolynomialDerivativeMatrix
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
CHARACTER(LEN=5)           :: tmpstr
INTEGER                    :: i
REAL,DIMENSION(0:Ngeo_out) :: wBary_Ngeo_out
REAL,ALLOCATABLE           :: xi(:)
!===================================================================================================================================

checkJacobian=GETLOGICAL('checkJacobian','T')

WRITE(tmpstr,'(I5)')2*Ngeo_out
Nanalyze=GETINT('Nanalyze',tmpstr)

ALLOCATE(xi(0:Nanalyze))
DO i=0,Nanalyze
  xi(i)=-1.+2*REAL(i)/REAL(Nanalyze)
END DO
CALL BarycentricWeights(Ngeo_out,xi_Ngeo_out,wBary_Ngeo_out)
ALLOCATE(Vdm_analyze(0:Nanalyze,0:Ngeo_out))
CALL InitializeVandermonde(Ngeo_out,Nanalyze,wBary_Ngeo_out,xi_Ngeo_out,xi,Vdm_analyze)

ALLOCATE(D_Ngeo_out(0:Ngeo_out,0:Ngeo_out))

CALL PolynomialDerivativeMatrix(Ngeo_out,xi_Ngeo_out,D_Ngeo_out)

END SUBROUTINE InitAnalyze


SUBROUTINE Analyze()
!===================================================================================================================================
! Basic Analyze initialization. 
!===================================================================================================================================
! MODULES
USE MOD_Analyze_Vars,ONLY:checkJacobian
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
!===================================================================================================================================
IF(CheckJacobian) CALL CheckJac()
END SUBROUTINE Analyze


SUBROUTINE CheckJac()
!===================================================================================================================================
! Basic Analyze initialization. 
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Mesh_Vars,   ONLY:Ngeo_out,XgeoQuad
USE MOD_Mesh_Vars,   ONLY:nQuads
USE MOD_Analyze_Vars,ONLY:Nanalyze,Vdm_analyze,D_Ngeo_out
USE MOD_ChangeBasis, ONLY:ChangeBasis3D
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
INTEGER                                            :: i,j,k,l,iQuad
REAL                                               :: scaledJacTol,Xbary(3)
REAL                                               :: scaledJac(nQuads),minJac,maxJac
INTEGER                                            :: scaledJacStat(0:10)
REAL,DIMENSION(3,3,0:Ngeo_out,0:Ngeo_out,0:Ngeo_out) :: dX
REAL,DIMENSION(3,3,0:Nanalyze,0:Nanalyze,0:Nanalyze) :: dXa
REAL,DIMENSION(0:Nanalyze,0:Nanalyze,0:Nanalyze)   :: detJac
!===================================================================================================================================
SWRITE(UNIT_stdOut,'   (A)') "Checking Jacobians..."
DO iQuad=1,nQuads
  dX=0.
  DO k=0,Ngeo_out
    DO j=0,Ngeo_out
      DO i=0,Ngeo_out
      ! Matrix-vector multiplication
        DO l=0,Ngeo_out
          dX(1,:,i,j,k)=dX(1,:,i,j,k) + D_Ngeo_out(i,l)*XgeoQuad(:,l,j,k,iQuad)
          dX(2,:,i,j,k)=dX(2,:,i,j,k) + D_Ngeo_out(j,l)*XgeoQuad(:,i,l,k,iQuad)
          dX(3,:,i,j,k)=dX(3,:,i,j,k) + D_Ngeo_out(k,l)*XgeoQuad(:,i,j,l,iQuad)
        END DO !l=0,Ngeo_out
      END DO !i=0,Ngeo_out
    END DO !j=0,Ngeo_out
  END DO !k=0,Ngeo_out
  CALL ChangeBasis3D(3,NGeo_out,Nanalyze,Vdm_analyze,dX(1,:,:,:,:),dXa(1,:,:,:,:))
  CALL ChangeBasis3D(3,NGeo_out,Nanalyze,Vdm_analyze,dX(2,:,:,:,:),dXa(2,:,:,:,:))
  CALL ChangeBasis3D(3,NGeo_out,Nanalyze,Vdm_analyze,dX(3,:,:,:,:),dXa(3,:,:,:,:))
  DO k=0,Nanalyze
    DO j=0,Nanalyze
      DO i=0,Nanalyze
        detJac(i,j,k)= dXa(1,1,i,j,k)*(dXa(2,2,i,j,k)*dXa(3,3,i,j,k) - dXa(3,2,i,j,k)*dXa(2,3,i,j,k)) + &
                       dXa(2,1,i,j,k)*(dXa(3,2,i,j,k)*dXa(1,3,i,j,k) - dXa(1,2,i,j,k)*dXa(3,3,i,j,k)) + &
                       dXa(3,1,i,j,k)*(dXa(1,2,i,j,k)*dXa(2,3,i,j,k) - dXa(2,2,i,j,k)*dXa(1,3,i,j,k))
      END DO !i
    END DO !j
  END DO !k
  maxJac=MAXVAL(ABS(detJac))
  minJac=MINVAL(detJac)
  scaledJac(iQuad)=minJac/maxJac
END DO !iQuad

! Error Section
scaledJacTol=1.0E-6

IF(ANY(scaledJac.LE.scaledJacTol))THEN
  OPEN(UNIT=100,FILE='Jacobian_Error.dat',STATUS='UNKNOWN',ACTION='WRITE')
  WRITE(100,*) 'TITLE="Corrupt Elems Found, i.e. the Jacobian is negative" '
  WRITE(100,*)
  WRITE(100,'(A)')'VARIABLES="xBary","yBary","zBary","scaledJac"'
  WRITE(100,*)
  DO iQuad=1,nQuads
    IF(scaledJac(iQuad).LT.scaledJactol) THEN
      XBary(1)=SUM(XgeoQuad(1,:,:,:,iQuad))
      XBary(2)=SUM(XgeoQuad(2,:,:,:,iQuad))
      XBary(3)=SUM(XgeoQuad(3,:,:,:,iQuad))
      XBary=XBary/REAL((Ngeo_out+1)**3)
      WRITE(100,'(4(4X,E12.5))') Xbary,scaledJac(iQuad)
    END IF
  END DO !iQaud
  CLOSE(UNIT=100)
  WRITE(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
  WRITE(*,*) '!!! WARNING! Elements with negative Jacobian found, see Jacobian_Error.dat for details.'
  WRITE(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
END IF
scaledJacStat(:)=0

DO iQuad=1,nQuads
  i=CEILING(MAX(0.,scaledJac(iQuad)*10))
  scaledJacStat(i)=scaledJacStat(i)+1 
END DO
WRITE(Unit_StdOut,'(A)') ' Number of element with scaled Jacobians ranging between:'
WRITE(Unit_StdOut,'(A)') '   <  0.0  <  0.1  <  0.2  <  0.3  <  0.4  <  0.5  <  0.6  <  0.7  <  0.8  <  0.9  <  1.0 '
DO i=0,10
  WRITE(Unit_StdOut,'(I6,X,A1)',ADVANCE='NO')scaledJacStat(i),'|'
END DO
WRITE(Unit_StdOut,'(A1)')' '
END SUBROUTINE CheckJac

END MODULE MOD_Analyze
