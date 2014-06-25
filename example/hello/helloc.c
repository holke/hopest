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
#include <hopest_fortran_test.h>

int
main (int argc, char **argv)
{
  int                 mpiret;
  int                 param;
#ifdef HOPEST_ENABLE_DEBUG
  const int           LP_lib = SC_LP_INFO;
  const int           LP_hopest = SC_LP_DEBUG;
#else
  const int           LP_lib = SC_LP_ESSENTIAL;
  const int           LP_hopest = SC_LP_PRODUCTION;
#endif

  mpiret = sc_MPI_Init (&argc, &argv);
  SC_CHECK_MPI (mpiret);

  sc_init (sc_MPI_COMM_WORLD, 1, 1, NULL, LP_lib);
  p4est_init (NULL, LP_lib);
  hopest_init (NULL, LP_hopest);

  hopest_global_essentialf ("Hopest says %s\n", "hello world");

  param = 2;
  HOPEST_C_AND_FORTRAN_MESSAGE_F77 (&param);

  param = 3;
  HOPEST_FORTRAN_AND_C_MESSAGE_F77 (&param);

  sc_finalize ();
  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);

  return 0;
}
