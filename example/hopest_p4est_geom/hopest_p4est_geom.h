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
#define ReadMeshFromHDF5_FC \
  HOPEST_FC_FUNC (wrapreadmeshfromhdf5nobuildp4est,WRAPREADMESHFROMHDF5NOBUILDP4EST)
#define InitRefineBoundaryElems_FC \
  HOPEST_FC_FUNC (wrapinitrefineboundaryelems, WRAPINITREFINEBOUNDARYELEMS)
#define InitRefineGeom_FC \
  HOPEST_FC_FUNC (wrapinitrefinegeom, WRAPINITREFINEGEOM)
#define RefineByList_FC \
  HOPEST_FC_FUNC (wraprefinebylist, WRAPREFINEBYLIST)
#define RefineByGeom_FC \
  HOPEST_FC_FUNC (wraprefinebygeom, WRAPREFINEBYGEOM)

#ifdef __cplusplus
extern              "C"         /* prevent C++ name mangling */
{
#if 0
}
#endif
#endif

void                FillStrings_FC (char *inifile, int inifile_len);
void                ReadMeshFromHDF5_FC (char *hdf5file, int hdf5file_len,
                                         p8est_connectivity_t ** conn);
void                InitRefineGeom_FC (void);
void                InitRefineBoundaryElems_FC(void);
int                 RefineByGeom_FC (p4est_qcoord_t, p4est_qcoord_t,
                                     p4est_qcoord_t, p4est_topidx_t, int8_t,
                                     int);
int                 RefineByList_FC (p4est_qcoord_t, p4est_qcoord_t,
                                     p4est_qcoord_t, p4est_topidx_t, int8_t,
                                     int);

#ifdef __cplusplus
#if 0
{
#endif
}
#endif

#endif /* _HOPEST_P4EST_GEOM_H */
