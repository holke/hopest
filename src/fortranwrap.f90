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

!USE MOD_ReadInTools
!USE MOD_Mesh
!USE MOD_MeshVars

! To pass strings from C to Fortran 2003 we pass the string as char* and its
! length as an integer
SUBROUTINE wrapinitmesh(hdf5file,hdf5file_len,conn)
    USE MOD_Mesh
    USE MOD_Mesh_Vars
    USE MOD_Mesh_ReadIn
    USE MOD_IO_HDF5
    USE MOD_MPI,          ONLY:InitMPI

    integer(C_INT), VALUE :: hdf5file_len
    character(hdf5file_len,kind=C_CHAR), intent(IN) :: hdf5file
    TYPE(C_PTR)           :: conn

    CALL InitMPI()
    call InitIO()
    print *,"calling read mesh with file ", hdf5file
    call ReadMesh(hdf5file)
    conn=p4est_ptr%connectivity
    !p4est=p4est_ptr%p4est
    !print *,"Address of p4est ", p4est
    !p4est=p4est_from_initmesh
END SUBROUTINE wrapinitmesh

! To pass strings from C to Fortran 2003 we pass the string as char* and its
! length as an integer
SUBROUTINE wrapfillstrings(inifile,inifile_len)
    USE MOD_ReadInTools
    USE, intrinsic :: ISO_C_BINDING
    implicit none
    integer(C_INT), intent(IN), VALUE :: inifile_len
    character(inifile_len,kind=C_CHAR), intent(IN) :: inifile
    print *,"calling FillString routine with parameter ", inifile
    call fillstrings(inifile)
END SUBROUTINE wrapfillstrings


