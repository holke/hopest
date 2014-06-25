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

#ifndef HOPEST_FORTRAN_TEST_H
#define HOPEST_FORTRAN_TEST_H

#include <hopest.h>

#define HOPEST_F_MESSAGE_F77 HOPEST_F77_FUNC_ (hopest_f_message, HOPEST_F_MESSAGE)
#ifdef __cplusplus
extern              "C"         /* prevent C++ name mangling */
#endif
void                HOPEST_F_MESSAGE_F77 (int *i);

#define HOPEST_C_AND_FORTRAN_MESSAGE_F77 HOPEST_F77_FUNC_ (hopest_c_and_fortran_message, HOPEST_C_AND_FORTRAN_MESSAGE)
#ifdef __cplusplus
extern              "C"         /* prevent C++ name mangling */
#endif
void                HOPEST_C_AND_FORTRAN_MESSAGE_F77 (int *i);

#define HOPEST_C_MESSAGE_F77 HOPEST_F77_FUNC_ (hopest_c_message, HOPEST_C_MESSAGE)
#ifdef __cplusplus
extern              "C"         /* prevent C++ name mangling */
#endif
void                HOPEST_C_MESSAGE (int *i);

#define HOPEST_FORTRAN_AND_C_MESSAGE_F77 HOPEST_F77_FUNC_ (hopest_fortran_and_c_message, HOPEST_FORTRAN_AND_C_MESSAGE)
#ifdef __cplusplus
extern              "C"         /* prevent C++ name mangling */
#endif
void                hopest_fortran_and_c_message (int *i);

#endif
