! This file is part of hopest.
! hopest is a Fortran/C library and application for high-order mesh
! preprocessing and interfacing to the p4est apaptive mesh library.
!
! Copyright (C) 2014 by the developers.
!
! hopest is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 2 of the License, or
! (at your option) any later version.
!
! hopest is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with hopest; if not, write to the Free Software Foundation, Inc.,
! 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.

      subroutine hopest_f90_message (i)
        INTEGER i
        print *, "This line was compiled with the Fortran 90 compiler.", &
     &           " Parameter: ", i
        return
      end

      subroutine hopest_fortran90_and_c_message (i)
        INTEGER i
        call hopest_f90_message (i)
        call hopest_c_message_f90 (i)
        return
      end
