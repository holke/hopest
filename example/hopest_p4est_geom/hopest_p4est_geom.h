/*
  This file is part of hopest.
  hopest is a Fortran/C library and application for high-order mesh
  preprocessing and interfacing to the p4est apaptive mesh library.

  Copyright (C) 2014 by the developers.

  hopest is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  hopest is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with hopest; if not, write to the Free Software Foundation, Inc.,
  51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
*/

#ifndef _HOPEST_P4EST_GEOM_H
#define _HOPEST_P4EST_GEOM_H


#define FillStrings_FC \
  HOPEST_FC_FUNC (wrapfillstrings, WRAPFILLSTRINGS)
#define InitMesh_FC \
  HOPEST_FC_FUNC (wrapinitmesh, WRAPINITMESH)

#ifdef __cplusplus
extern              "C"         /* prevent C++ name mangling */
{
#if 0
}
#endif
#endif

void FillStrings_FC(char *inifile,int inifile_len);
void InitMesh_FC(char *hdf5file,int hdf5file_len,p8est_connectivity_t **conn);

#ifdef __cplusplus
#if 0
{
#endif
}
#endif

#endif /* _HOPEST_P4EST_GEOM_H */
