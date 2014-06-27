c This file is part of hopest.
c hopest is a Fortran/C library and application for high-order mesh
c preprocessing and interfacing to the p4est apaptive mesh library.
c
c Copyright (C) 2014 by the developers.
c
c hopest is free software; you can redistribute it and/or modify
c it under the terms of the GNU General Public License as published by
c the Free Software Foundation; either version 2 of the License, or
c (at your option) any later version.
c
c hopest is distributed in the hope that it will be useful,
c but WITHOUT ANY WARRANTY; without even the implied warranty of
c MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
c GNU General Public License for more details.
c
c You should have received a copy of the GNU General Public License
c along with hopest; if not, write to the Free Software Foundation, Inc.,
c 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.

      program hellof
        print *, "Hopest says hello world"
        call hopest_fortran77_and_c_message (0)
        call hopest_c_and_fortran_message_f77 (1)
      end program hellof
