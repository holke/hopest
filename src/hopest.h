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

/*
 * This is a file that every program using con4m should include,
 * directly or indirectly.  It pulls in subpackage definitions.
 */

#ifndef CON4M_H
#define CON4M_H

/* include config headers */
#include <con4m_config.h>
#include <sc_config.h>
#if \
  (defined (CON4M_MPI) && !defined (SC_MPI)) || \
  (!defined (CON4M_MPI) && defined (SC_MPI))
#error "MPI configured differently in con4m and libsc"
#endif
#if \
  (defined (CON4M_MPIIO) && !defined (SC_MPIIO)) || \
  (!defined (CON4M_MPIIO) && defined (SC_MPIIO))
#error "MPI I/O configured differently in con4m and libsc"
#endif

/* indirectly also include p4est_config.h, sc.h, sc_containers.h */
#include <p4est_base.h>
#define _con4m_const _sc_const
#define _con4m_restrict _sc_restrict

/* macros for memory allocation, will abort if out of memory */
#define CON4M_ALLOC(t,n)          (t *) sc_malloc (con4m_get_package_id (), \
                                                   (n) * sizeof(t))
#define CON4M_ALLOC_ZERO(t,n)     (t *) sc_calloc (con4m_get_package_id (), \
                                                   (size_t) (n), sizeof(t))
#define CON4M_REALLOC(p,t,n)      (t *) sc_realloc (con4m_get_package_id (), \
                                                    (p), (n) * sizeof(t))
#define CON4M_STRDUP(s)                 sc_strdup (con4m_get_package_id (), (s))
#define CON4M_FREE(p)                   sc_free (con4m_get_package_id (), (p))

/* some error checking */
#ifdef CON4M_DEBUG
#define CON4M_ASSERT(c) SC_CHECK_ABORT ((c), "Assertion '" #c "'")
#else
#define CON4M_ASSERT(c) SC_NOOP ()
#endif

/* functions for mananging the package identity within libsc */
int                 con4m_get_package_id (void);

/* functions for printing log messages */
void                con4m_logf (int category, int priority,
                                const char *fmt, ...);
void                con4m_global_essentialf (const char *fmt, ...);
void                con4m_global_productionf (const char *fmt, ...);

/* register con4m with libsc and print version and variable information */
void                con4m_init (sc_log_handler_t log_handler,
                                int log_threshold);

#endif /* !CON4M_H */
