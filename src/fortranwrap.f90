!  This file is part of hopest.
!  hopest is a Fortran/C library and application for high-order mesh
!  preprocessing and interfacing to the p4est apaptive mesh library.
!
!  Copyright (C) 2014 by the developers.
!
!  hopest is free software; you can redistribute it and/or modify
!  it under the terms of the GNU General Public License as published by
!  the Free Software Foundation; either version 2 of the License, or
!  (at your option) any later version.
!
!  hopest is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!  GNU General Public License for more details.
!
!  You should have received a copy of the GNU General Public License
!  along with hopest; if not, write to the Free Software Foundation, Inc.,
!  51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.

! This file wraps some fortran routines to be called by C-Code
! It is important that this is NOT a Fortran Module

#include <hopest_f.h>

SUBROUTINE wrapReadMeshFromHDF5nobuildp4est(hdf5file,hdf5file_len,conn,communicator)
!===================================================================================================================================
! Subroutine to wrap the ReadMeshFromHDF5nobuildp4est routine
!===================================================================================================================================
! MODULES
USE MPI
USE MODH_Globals,       ONLY: comm
USE MODH_IO_HDF5,       ONLY: InitIO
USE MODH_Mesh_ReadIn,   ONLY: ReadMeshHeader,ReadMeshFromHDF5
USE MODH_P4EST_Vars,    ONLY: connectivity
USE, INTRINSIC :: ISO_C_BINDING
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER(C_INT), INTENT(IN), VALUE               :: hdf5file_len
CHARACTER(HDF5FILE_LEN,KIND=C_CHAR), INTENT(IN) :: hdf5file
INTEGER(C_INT), INTENT(IN), VALUE               :: communicator
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
TYPE(C_PTR),INTENT(OUT)                         :: conn
!===================================================================================================================================
! To pass strings from C to Fortran 2003 we pass the string as char* and its
! length as an integer
IF(communicator.EQ.0) THEN
    comm=MPI_COMM_WORLD
ELSE
    comm=MPI_COMM_SELF
END IF
CALL InitIO()
CALL ReadMeshHeader(hdf5file)
CALL ReadMeshFromHDF5(hdf5file)
conn=connectivity
END SUBROUTINE wrapReadMeshFromHDF5nobuildp4est


SUBROUTINE wrapbuildHOp4GeometryX(a,b,c,x,y,z,tree)
!===================================================================================================================================
! Subroutine to wrap the buildHOp4GeometryX routine
!===================================================================================================================================
! MODULES
USE MODH_P4EST,ONLY:buildHOp4GeometryX
USE, INTRINSIC :: ISO_C_BINDING
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL(KIND=C_DOUBLE),INTENT(IN),VALUE    :: a,b,c
P4EST_F90_TOPIDX   ,INTENT(IN),VALUE    :: tree
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL(KIND=C_DOUBLE),INTENT(OUT)         :: x,y,z
!===================================================================================================================================
CALL buildHOp4GeometryX(a,b,c,x,y,z,tree)
END SUBROUTINE wrapbuildHOp4GeometryX

SUBROUTINE wrapInitRefineGeom()
!===================================================================================================================================
! Subroutine to wrap the ReadMeshFromHDF5nobuildp4est routine
!===================================================================================================================================
! MODULES
    USE MODH_Refine,     ONLY: InitRefineGeom
!-----------------------------------------------------------------------------------------------------------------------------------
    call InitRefineGeom()
END SUBROUTINE wrapInitRefineGeom

SUBROUTINE wrapInitRefineBoundaryElems(level)
!===================================================================================================================================
! Subroutine to wrap the InitRefineBoundaryElems routine
!===================================================================================================================================
! MODULES
    USE MODH_Refine,     ONLY: InitRefineBoundaryElems
    USE MODH_Refine_Vars, ONLY: TreeSidesToRefine,refineBCIndex,refineLevel
    USE, INTRINSIC :: ISO_C_BINDING
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
    INTEGER(KIND=C_INT),INTENT(IN),VALUE :: level
!-----------------------------------------------------------------------------------------------------------------------------------
    refineBCIndex=1
    refineLevel=level
    call InitRefineBoundaryElems
END SUBROUTINE wrapInitRefineBoundaryElems

FUNCTION wrapRefineByList(x,y,z,tree,level,childID)
!===================================================================================================================================
! Subroutine to wrap the RefineByList routine
!===================================================================================================================================
! MODULES
    USE MODH_Refine,     ONLY: RefineByList
    USE, INTRINSIC :: ISO_C_BINDING
!-----------------------------------------------------------------------------------------------------------------------------------
    P4EST_F90_QCOORD,INTENT(IN),VALUE :: x,y,z
    P4EST_F90_TOPIDX,INTENT(IN),VALUE :: tree
    P4EST_F90_QLEVEL,INTENT(IN),VALUE :: level
    INTEGER(KIND=C_INT ),INTENT(IN),VALUE    :: childID
    INTEGER :: wrapRefineByList
    wrapRefineByList=RefineByList(x,y,z,tree,level,childID)
END FUNCTION wrapRefineByList

FUNCTION wrapRefineByGeom(x,y,z,tree,level,childID)
!===================================================================================================================================
! Subroutine to wrap the RefineByGeom routine
!===================================================================================================================================
! MODULES
    USE MODH_Refine,     ONLY: RefineByGeom
    USE, INTRINSIC :: ISO_C_BINDING
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
    P4EST_F90_QCOORD,INTENT(IN),VALUE :: x,y,z
    P4EST_F90_TOPIDX,INTENT(IN),VALUE :: tree
    P4EST_F90_QLEVEL,INTENT(IN),VALUE :: level
    INTEGER(KIND=C_INT ),INTENT(IN),VALUE    :: childID
    INTEGER :: wrapRefineByGeom
!-----------------------------------------------------------------------------------------------------------------------------------
    wrapRefineByGeom=RefineByGeom(x,y,z,tree,level,childID)
END FUNCTION wrapRefineByGeom
