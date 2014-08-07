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

/*
 * This is a file that every program using hopest should include,
 * directly or indirectly.  It pulls in subpackage definitions.
 */

#ifndef HOPEST_H
#define HOPEST_H

/* include config headers */
#include <hopest_config.h>
#include <sc_config.h>
#if \
  (defined (HOPEST_MPI) && !defined (SC_MPI)) || \
  (!defined (HOPEST_MPI) && defined (SC_MPI))
#error "MPI configured differently in hopest and libsc"
#endif
#if \
  (defined (HOPEST_MPIIO) && !defined (SC_MPIIO)) || \
  (!defined (HOPEST_MPIIO) && defined (SC_MPIIO))
#error "MPI I/O configured differently in hopest and libsc"
#endif

/* indirectly also include p4est_config.h, sc.h, sc_containers.h */
#include <p4est_base.h>
#define _hopest_const _sc_const
#define _hopest_restrict _sc_restrict

/* macros for memory allocation, will abort if out of memory */
#define HOPEST_ALLOC(t,n)          (t *) sc_malloc (hopest_get_package_id (), \
                                                   (n) * sizeof(t))
#define HOPEST_ALLOC_ZERO(t,n)     (t *) sc_calloc (hopest_get_package_id (), \
                                                   (size_t) (n), sizeof(t))
#define HOPEST_REALLOC(p,t,n)      (t *) sc_realloc (hopest_get_package_id (), \
                                                    (p), (n) * sizeof(t))
#define HOPEST_STRDUP(s)                 sc_strdup (hopest_get_package_id (), (s))
#define HOPEST_FREE(p)                   sc_free (hopest_get_package_id (), (p))

/* some error checking */
#ifdef HOPEST_DEBUG
#define HOPEST_ASSERT(c) SC_CHECK_ABORT ((c), "Assertion '" #c "'")
#else
#define HOPEST_ASSERT(c) SC_NOOP ()
#endif

/* protect against C++ compilation */
#ifdef __cplusplus
extern              "C"         /* prevent C++ name mangling */
{
#if 0
}
#endif
#endif

/* functions for mananging the package identity within libsc */
int                 hopest_get_package_id (void);

/* functions for printing log messages */
void                hopest_logf (int category, int priority,
                                const char *fmt, ...);
void                hopest_global_essentialf (const char *fmt, ...);
void                hopest_global_productionf (const char *fmt, ...);

/* register hopest with libsc and print version and variable information */
void                hopest_init (sc_log_handler_t log_handler,
                                int log_threshold);

#ifdef __cplusplus
#if 0
{
#endif
}
#endif

#endif /* !HOPEST_H */
