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

/* This program uses the hopest routine ReadMeshFromHDF5 to read
 * a p4est_connectivity from an hdf5-file
 * using this connectivity we build a p4est and use
 * the hopest routine buildHOp4GeometryX to get the Geometry information from the file */

#include <hopest.h>
#include <string.h>
#include <p8est_extended.h>
#include <p8est_vtk.h>
#include <p8est_bits.h>
#include <p4est/p4est_HO_geometry.h>
#include "hopest_p4est_geom.h"

#include <p4est_to_p8est.h>

int
refine (p4est_t * p4est, p4est_topidx_t which_tree, p4est_quadrant_t * q)
{

  return RefineByList_FC (q->x, q->y, q->z, which_tree + 1, q->level,
                          p8est_quadrant_child_id (q));
}

int
main (int argc, char *argv[])
{
  char               *HDF5File;
  int                 HDF5file_len;
  /* p4est_t *p4est; */
  p4est_connectivity_t *conn;
  p4est_t            *p4est;
  p4est_geometry_t   *geom;
  sc_MPI_Comm         comm;
  char               *vtkfilename, *vtkfilename_temp;
  int                 mpiret, mpirank, level;
  const char         *usage;

  usage = "Arguments: <hdf5-file> \n";

  mpiret = sc_MPI_Init (&argc, &argv);
  SC_CHECK_MPI (mpiret);
  comm = sc_MPI_COMM_WORLD;
  mpiret = sc_MPI_Comm_rank (comm, &mpirank);
  if (argc == 2) {
    HDF5File = argv[1];
    HDF5file_len = strlen (HDF5File);

    ReadMeshFromHDF5_FC (HDF5File, HDF5file_len, &conn);
    P4EST_ASSERT (p4est_connectivity_is_valid (conn));

    p4est = p4est_new_ext (sc_MPI_COMM_WORLD, conn, 0, 0, 1, 0, NULL, NULL);
    geom = P4EST_ALLOC_ZERO (p4est_geometry_t, 1);
    geom->name = "hopest_readfromhdf5";
    geom->X = p4_geometry_X;

    InitRefineBoundaryElems_FC ();
    for (level = 0; level < 5; ++level) {
      p4est_refine (p4est, 0, refine, NULL);
      /* Refinement has lead to up to 8x more elements; redistribute them. */
      p4est_partition (p4est, 0, NULL);
    }

    vtkfilename_temp = P4EST_STRDUP (HDF5File);
    vtkfilename = basename (vtkfilename_temp);
    p4est_vtk_write_file (p4est, geom, vtkfilename);

    P4EST_FREE (vtkfilename_temp);
    p4est_geometry_destroy (geom);
    p4est_destroy (p4est);
    p4est_connectivity_destroy (conn);
  }
  else {
    P4EST_GLOBAL_LERROR (usage);
    sc_abort_collective ("Usage error");
  }
  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);
  return 0;
}
