/*
  This file is part of con4m.
  con4m is a C library for the parallel adaptive discretization and solution
  of partial differential equations (PDEs).  It relies on p4est and libsc.

  Copyright (C) 2013, 2014 Carsten Burstedde and others.

  con4m is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  con4m is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with con4m; if not, write to the Free Software Foundation, Inc.,
  51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
*/

#include <con4m.h>

static int          con4m_package_id = -1;

int
con4m_get_package_id (void)
{
  return con4m_package_id;
}

static void
con4m_logv (int category, int priority, const char *fmt, va_list ap)
{
  sc_logv ("unknown", -1, con4m_package_id, category, priority, fmt, ap);
}

void
con4m_logf (int category, int priority, const char *fmt, ...)
{
  va_list             ap;

  va_start (ap, fmt);
  sc_logv ("unknown", -1, con4m_package_id, category, priority, fmt, ap);
  va_end (ap);
}

void
con4m_global_essentialf (const char *fmt, ...)
{
  va_list             ap;

  va_start (ap, fmt);
  con4m_logv (SC_LC_GLOBAL, SC_LP_ESSENTIAL, fmt, ap);
  va_end (ap);
}

void
con4m_global_productionf (const char *fmt, ...)
{
  va_list             ap;

  va_start (ap, fmt);
  con4m_logv (SC_LC_GLOBAL, SC_LP_PRODUCTION, fmt, ap);
  va_end (ap);
}

void
con4m_init (sc_log_handler_t log_handler, int log_threshold)
{
  int                 w;

  con4m_package_id = sc_package_register (log_handler, log_threshold,
                                          "con4m",
                                          "Adaptive discretizations");

  w = 24;
  con4m_global_essentialf ("This is %s\n", CON4M_PACKAGE_STRING);
  con4m_global_productionf ("%-*s %s\n", w, "CPP", CON4M_CPP);
  con4m_global_productionf ("%-*s %s\n", w, "CPPFLAGS", CON4M_CPPFLAGS);
  con4m_global_productionf ("%-*s %s\n", w, "CC", CON4M_CC);
  con4m_global_productionf ("%-*s %s\n", w, "CFLAGS", CON4M_CFLAGS);
  con4m_global_productionf ("%-*s %s\n", w, "LDFLAGS", CON4M_LDFLAGS);
  con4m_global_productionf ("%-*s %s\n", w, "LIBS", CON4M_LIBS);
}
