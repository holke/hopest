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

static int          hopest_package_id = -1;

int
hopest_get_package_id (void)
{
  return hopest_package_id;
}

static void
hopest_logv (int category, int priority, const char *fmt, va_list ap)
{
  sc_logv ("unknown", -1, hopest_package_id, category, priority, fmt, ap);
}

void
hopest_logf (int category, int priority, const char *fmt, ...)
{
  va_list             ap;

  va_start (ap, fmt);
  sc_logv ("unknown", -1, hopest_package_id, category, priority, fmt, ap);
  va_end (ap);
}

void
hopest_global_essentialf (const char *fmt, ...)
{
  va_list             ap;

  va_start (ap, fmt);
  hopest_logv (SC_LC_GLOBAL, SC_LP_ESSENTIAL, fmt, ap);
  va_end (ap);
}

void
hopest_global_productionf (const char *fmt, ...)
{
  va_list             ap;

  va_start (ap, fmt);
  hopest_logv (SC_LC_GLOBAL, SC_LP_PRODUCTION, fmt, ap);
  va_end (ap);
}

void
hopest_init (sc_log_handler_t log_handler, int log_threshold)
{
  int                 w;

  hopest_package_id = sc_package_register (log_handler, log_threshold,
                                          "hopest",
                                          "High-order mesh preprocessor");

  w = 24;
  hopest_global_essentialf ("This is %s\n", HOPEST_PACKAGE_STRING);
  hopest_global_productionf ("%-*s %s\n", w, "CPP", HOPEST_CPP);
  hopest_global_productionf ("%-*s %s\n", w, "CPPFLAGS", HOPEST_CPPFLAGS);
  hopest_global_productionf ("%-*s %s\n", w, "F77", HOPEST_F77);
  hopest_global_productionf ("%-*s %s\n", w, "FFLAGS", HOPEST_FFLAGS);
  hopest_global_productionf ("%-*s %s\n", w, "CC", HOPEST_CC);
  hopest_global_productionf ("%-*s %s\n", w, "CFLAGS", HOPEST_CFLAGS);
  hopest_global_productionf ("%-*s %s\n", w, "LDFLAGS", HOPEST_LDFLAGS);
  hopest_global_productionf ("%-*s %s\n", w, "LIBS", HOPEST_LIBS);
}
