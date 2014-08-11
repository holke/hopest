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
    USE MODH_Mesh_ReadIn
    USE MODH_IO_HDF5
    USE MODH_Mesh,          ONLY: InitMesh
    USE MODH_P4EST_Vars,    ONLY: connectivity
    USE MODH_P4EST,         ONLY: InitP4est
    USE, intrinsic :: ISO_C_BINDING
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
    integer(C_INT), intent(IN), VALUE               :: hdf5file_len
    character(hdf5file_len,kind=C_CHAR), intent(IN) :: hdf5file
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
    TYPE(C_PTR)                                     :: conn
!-----------------------------------------------------------------------------------------------------------------------------------
    CALL InitIO()
    call ReadMeshFromHDF5nobuildp4est(hdf5file)
    conn=connectivity
END SUBROUTINE wrapReadMeshFromHDF5nobuildp4est

! To pass strings from C to Fortran 2003 we pass the string as char* and its
! length as an integer

SUBROUTINE wrapbuildHOp4GeometryX(a,b,c,x,y,z,tree)
!===================================================================================================================================
! Subroutine to wrap the buildHOp4GeometryX routine
!===================================================================================================================================
! MODULES
    USE MODH_P4EST,ONLY:buildHOp4GeometryX
    USE, intrinsic :: ISO_C_BINDING
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
    REAL( KIND = C_DOUBLE ),INTENT(IN),VALUE    :: a,b,c
    P4EST_F90_TOPIDX,INTENT(IN),VALUE           :: tree
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
    REAL( KIND = C_DOUBLE ),INTENT(OUT)         :: x,y,z
!-----------------------------------------------------------------------------------------------------------------------------------
    call buildHOp4GeometryX(a,b,c,x,y,z,tree)
END SUBROUTINE wrapbuildHOp4GeometryX
