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

/* Functions written in Fortran 77 */
#define hopest_f77_message \
  HOPEST_F77_FUNC_ (hopest_f77_message, HOPEST_F77_MESSAGE)
#define hopest_fortran77_and_c_message \
  HOPEST_F77_FUNC_ (hopest_fortran77_and_c_message, HOPEST_FORTRAN77_AND_C_MESSAGE)

/* Functions written in Fortran 90 */
#define hopest_f90_message \
  HOPEST_FC_FUNC_ (hopest_f90_message, HOPEST_F90_MESSAGE)
#define hopest_fortran90_and_c_message \
  HOPEST_FC_FUNC_ (hopest_fortran90_and_c_message, HOPEST_FORTRAN90_AND_C_MESSAGE)

/* Functions written in C */
/* for Fortran 77 */
#define hopest_c_message_f77 \
  HOPEST_F77_FUNC_ (hopest_c_message_f77, HOPEST_C_MESSAGE_F77)
#define hopest_c_and_fortran_message_f77 \
  HOPEST_F77_FUNC_ (hopest_c_and_fortran_message_f77, HOPEST_C_AND_FORTRAN_MESSAGE_F77)

/* for Fortran 90 */
#define hopest_c_message_f90 \
  HOPEST_FC_FUNC_ (hopest_c_message_f90, HOPEST_C_MESSAGE_F90)
#define hopest_c_and_fortran_message_f90 \
  HOPEST_FC_FUNC_ (hopest_c_and_fortran_message_f90, HOPEST_C_AND_FORTRAN_MESSAGE_F90)

#ifdef __cplusplus
extern              "C"         /* prevent C++ name mangling */
{
#if 0
}
#endif
#endif

/* Functions defined in hopest_c_test.f */
void                hopest_f77_message (int *i);
void                hopest_fortran77_and_c_message (int *i);

/* Functions defined in hopest_c_test.f90 */
void                hopest_f90_message (int *i);
void                hopest_fortran90_and_c_message (int *i);

/* Functions defined in hopest_fortran_test.c */
void                hopest_c_message_f77 (int *i);
void                hopest_c_and_fortran_message_f77 (int *i);
void                hopest_c_message_f90 (int *i);
void                hopest_c_and_fortran_message_f90 (int *i);

#ifdef __cplusplus
#if 0
{
#endif
}
#endif

#endif /* !HOPEST_FORTRAN_TEST_H */
