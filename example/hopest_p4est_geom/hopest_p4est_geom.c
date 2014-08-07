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
#include <p4est/p4est_HO_geometry.h>
#include "hopest_p4est_geom.h"

#include <p4est_to_p8est.h>

int main(int argc,char *argv[]){
    char *IniFile;
    int inifile_len;
    /* p4est_t *p4est; */
    p4est_connectivity_t *conn;
    p4est_t              *p4est;
    p4est_geometry_t     *geom;
    int mpiret;

    mpiret = sc_MPI_Init (&argc, &argv);
    SC_CHECK_MPI (mpiret);
    if(argc>1) {
        IniFile=argv[1];
        inifile_len=strlen(IniFile);
        ReadMeshFromHDF5_FC(IniFile,inifile_len,&conn);
        P4EST_ASSERT(p4est_connectivity_is_valid(conn));
        p4est=p4est_new_ext(sc_MPI_COMM_WORLD,conn,0,2,1,0,NULL,NULL);
        geom = P4EST_ALLOC_ZERO (p4est_geometry_t, 1);
        geom->name = "hopest_readfromhdf5";
        geom->X = p4_geometry_X;
        p4est_vtk_write_file (p4est,geom, P4EST_STRING "from_hdf5file_via_hopest");
        p4est_geometry_destroy(geom);
        p4est_destroy(p4est);
        p4est_connectivity_destroy(conn);
    }
    else printf("no input file given.\n");
    mpiret = sc_MPI_Finalize ();
    SC_CHECK_MPI (mpiret);
    return 0;
}
