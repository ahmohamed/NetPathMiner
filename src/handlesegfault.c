#include "handlesegfault.h"


// Handling Sigfault errors. Prevents the code from exiting R.
// This is only available for UNIX platforms.
#ifndef WIN_COMPILE
#ifdef HAVE_XML
void segfault_KGML(int signal, siginfo_t *si, void *arg){
	EVAL(lang2(install("registerMemoryErr"), mkString("KGML2igraph")));
	error("Critical memory error in KGML2igraph. Please save your work and restart R.");
}
#endif
#ifdef HAVE_SBML
void segfault_SBML(int signal, siginfo_t *si, void *arg){
	EVAL(lang2(install("registerMemoryErr"), mkString("SBML2igraph")));
	error("Critical memory error in SBML2igraph. Please save your work and restart R.");
}
#endif
#endif

#ifdef HAVE_XML
void handle_segfault_KGML(){
 #ifndef WIN_COMPILE
	struct sigaction sa;

    memset(&sa, 0, sizeof(struct sigaction));
    sigemptyset(&sa.sa_mask);
    sa.sa_sigaction = segfault_KGML;
    sa.sa_flags   = SA_SIGINFO;

    sigaction(SIGSEGV, &sa, NULL);
#endif
}
#endif

#ifdef HAVE_SBML
void handle_segfault_SBML(){
#ifndef WIN_COMPILE    
	struct sigaction sa;

    memset(&sa, 0, sizeof(struct sigaction));
    sigemptyset(&sa.sa_mask);
    sa.sa_sigaction = segfault_SBML;
    sa.sa_flags   = SA_SIGINFO;

    sigaction(SIGSEGV, &sa, NULL);
#endif
}
#endif
