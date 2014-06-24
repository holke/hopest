#ifndef HOPEST_FORTRAN_TEST_H
#define HOPEST_FORTRAN_TEST_H

#include <hopest.h>

#define HOPEST_F_MESSAGE_F77 HOPEST_F77_FUNC_ (hopest_f_message, HOPEST_F_MESSAGE)
#ifdef __cplusplus
extern              "C"         /* prevent C++ name mangling */
#endif
void                HOPEST_F_MESSAGE_F77 ();

#endif
