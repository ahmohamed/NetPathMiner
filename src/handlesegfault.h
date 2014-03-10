#ifndef __handlesegfault__h_
#define __handlesegfault__h_

#include <signal.h>
#include <string.h>


#include <R.h>
#include <Rmath.h>
#include <Rdefines.h>
#include <Rinternals.h>

// Handling Sigfault errors. Prevents the code from exiting R.
extern void handleKGML();
extern void handle_segfault_SBML();

#endif
