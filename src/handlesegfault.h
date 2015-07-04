#ifndef __handlesegfault__h_
#define __handlesegfault__h_

#include <signal.h>
#include <string.h>


#include <R.h>
#include <Rmath.h>
#include <Rdefines.h>
#include <Rinternals.h>

// Handling Sigfault errors. Prevents the code from exiting R.
#ifdef __cplusplus
extern "C" {
#endif

void handle_segfault_KGML();
void handle_segfault_SBML();

#ifdef __cplusplus
}
#endif


#endif
