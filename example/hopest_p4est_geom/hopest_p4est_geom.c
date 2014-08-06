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

#include <hopest.h>
#include <hopest_hdf5.h>
#include <string.h>
#include <p4est.h>
#include "hopest_p4est_geom.h"

int main(int argc,char *argv[]){
    char *IniFile;
    int inifile_len;
    p4est_t *p4est;
    p4est_connectivity_t **conn,*conn2;
    int mpiret;

    mpiret = sc_MPI_Init (&argc, &argv);
    SC_CHECK_MPI (mpiret);
    if(argc>1) {
        IniFile=argv[1];
        inifile_len=strlen(IniFile);
        /*FillStrings_FC(IniFile,inifile_len);*/
        InitMesh_FC(IniFile,inifile_len,conn);
        printf("Got connectivity address %p\n",*conn);
        conn2=*conn;
#if 1
        {
             char *mem=(char*)(conn2);
            printf("[%i]",(int)(mem[0]));
            printf("\n");
        }
#endif
        fflush(stdout);
        /*P4EST_ASSERT(p4est_connectivity_is_valid(*conn));*/

    }
    /* connectivity holen und p4est bauen */
    else printf("no input file given.\n");
    mpiret = sc_MPI_Finalize ();
    SC_CHECK_MPI (mpiret);
    return 0;
}
