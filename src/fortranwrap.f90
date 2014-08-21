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

SUBROUTINE wrapReadMeshFromHDF5nobuildp4est(hdf5file,hdf5file_len,conn)
!===================================================================================================================================
! Subroutine to wrap the ReadMeshFromHDF5nobuildp4est routine
!===================================================================================================================================
! MODULES
USE MODH_IO_HDF5,       ONLY: InitIO
USE MODH_Mesh_ReadIn,   ONLY: ReadMeshHeader,ReadMeshFromHDF5
USE MODH_P4EST_Vars,    ONLY: connectivity
USE, INTRINSIC :: ISO_C_BINDING
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER(C_INT), INTENT(IN), VALUE               :: hdf5file_len
CHARACTER(HDF5FILE_LEN,KIND=C_CHAR), INTENT(IN) :: hdf5file
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
TYPE(C_PTR),INTENT(OUT)                         :: conn
!===================================================================================================================================
! To pass strings from C to Fortran 2003 we pass the string as char* and its
! length as an integer
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

