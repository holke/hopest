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
#include <sc_flops.h>
#include <sc_statistics.h>
#include <sc_options.h>

#include <p4est_to_p8est.h>

enum
{
  TIME_READ_CONNECTIVITY,
  TIME_BUILD_P4EST,
  TIME_REFINE,
  TIME_BALANCE,
  NUM_STATS
};

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
  p4est_geometry_t   *geom = NULL;
  sc_MPI_Comm         comm;
  char               *vtkfilename, *vtkfilename_temp;
  int                 mpiret, mpirank;
  int                 refine_level, level;
  int                 first_argc;
  int                 use_connectivity_bcast, write_vtk;
  sc_flopinfo_t       fi, snapshot;
  sc_statinfo_t       stats[NUM_STATS];
  sc_options_t       *opt;

  mpiret = sc_MPI_Init (&argc, &argv);
  SC_CHECK_MPI (mpiret);
  comm = sc_MPI_COMM_WORLD;
  mpiret = sc_MPI_Comm_rank (comm, &mpirank);
  sc_init (comm, 1, 1, NULL, SC_LP_DEFAULT);

  opt = sc_options_new (argv[0]);
  sc_options_add_string (opt, 'f', "file", (const char **) &HDF5File, NULL,
                         "hdf5-file");
  sc_options_add_int (opt, 'l', "level", &refine_level, 5,
                      "maximum refine level");
  sc_options_add_switch (opt, 'c', "connbcast", &use_connectivity_bcast,
                         "only process zero reads connectivity and broadcasts it to the other processes");
  sc_options_add_switch (opt, 'v', "vtk", &write_vtk, "write vtk file");

  first_argc = sc_options_parse (p4est_package_id, SC_LP_DEFAULT,
                                 opt, argc, argv);
  if (first_argc < 0 || first_argc != argc || HDF5File == NULL) {
    sc_options_print_usage (p4est_package_id, SC_LP_ERROR, opt, NULL);
    return 1;
  }
 /* start overall timing */
  mpiret = sc_MPI_Barrier (comm);
  SC_CHECK_MPI (mpiret);
  sc_flops_start (&fi);

  HDF5file_len = strlen (HDF5File);

  sc_flops_snap (&fi, &snapshot);
  if (use_connectivity_bcast) {
    if (mpirank == 0)
      ReadMeshFromHDF5_FC (HDF5File, HDF5file_len, &conn, 1);
    conn = p4est_connectivity_bcast (conn, 0, comm);
  }
  else
    ReadMeshFromHDF5_FC (HDF5File, HDF5file_len, &conn, 0);
  sc_flops_shot (&fi, &snapshot);
  sc_stats_set1 (&stats[TIME_READ_CONNECTIVITY], snapshot.iwtime,
                 "Read connectivity");

  sc_flops_snap (&fi, &snapshot);
  p4est = p4est_new_ext (comm, conn, 0, 0, 1, 0, NULL, NULL);
  sc_flops_shot (&fi, &snapshot);
  sc_stats_set1 (&stats[TIME_BUILD_P4EST], snapshot.iwtime, "build p4est");

  sc_flops_snap (&fi, &snapshot);
  if (!use_connectivity_bcast) {
    InitRefineBoundaryElems_FC (refine_level);
    for (level = 0; level < refine_level; ++level) {
      p4est_refine (p4est, 0, refine, NULL);
      /* Refinement has lead to up to 8x more elements; redistribute them. */
      p4est_partition (p4est, 0, NULL);
    }
  }
  sc_flops_shot (&fi, &snapshot);
  sc_stats_set1 (&stats[TIME_REFINE], snapshot.iwtime,
                 "refine and partition");

  sc_flops_snap (&fi, &snapshot);
  p4est_balance (p4est, P4EST_CONNECT_FULL, NULL);
  p4est_partition (p4est, 0, NULL);
  sc_flops_shot (&fi, &snapshot);
  sc_stats_set1 (&stats[TIME_BALANCE], snapshot.iwtime, "balance p4est");

  if (write_vtk) {
    if(use_connectivity_bcast) geom=NULL;
    else {
        geom = P4EST_ALLOC_ZERO (p4est_geometry_t, 1);
        geom->name = "hopest_readfromhdf5";
        geom->X = p4_geometry_X;
    }
    vtkfilename_temp = P4EST_STRDUP (HDF5File);
    vtkfilename = basename (vtkfilename_temp);
    p4est_vtk_write_file (p4est, geom, vtkfilename);
    P4EST_FREE (vtkfilename_temp);
    if(geom!=NULL)p4est_geometry_destroy (geom);
  }

  sc_stats_compute (comm, NUM_STATS, stats);
  sc_stats_print (p4est_package_id, SC_LP_STATISTICS, NUM_STATS, stats, 1, 1);
  p4est_destroy (p4est);
  p4est_connectivity_destroy (conn);

  sc_options_destroy (opt);
  sc_finalize ();
  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);
  return 0;
}
