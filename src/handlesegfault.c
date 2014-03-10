#include "handlesegfault.h"

// Handling Sigfault errors. Prevents the code from exiting R.
void segfault_KGML(int signal, siginfo_t *si, void *arg){
	EVAL(lang2(install("registerMemoryErr"), mkString("KGML2igraph")));
	error("Critical memory error in KGML2igraph. Please save your work and restart R.");
}
extern void handleKGML(){
    struct sigaction sa;

    memset(&sa, 0, sizeof(struct sigaction));
    sigemptyset(&sa.sa_mask);
    sa.sa_sigaction = segfault_KGML;
    sa.sa_flags   = SA_SIGINFO;

    sigaction(SIGSEGV, &sa, NULL);
}

void segfault_SBML(int signal, siginfo_t *si, void *arg){
	EVAL(lang2(install("registerMemoryErr"), mkString("SBML2igraph")));
	error("Critical memory error in SBML2igraph. Please save your work and restart R.");
}
extern void handle_segfault_SBML(){
    struct sigaction sa;

    memset(&sa, 0, sizeof(struct sigaction));
    sigemptyset(&sa.sa_mask);
    sa.sa_sigaction = segfault_SBML;
    sa.sa_flags   = SA_SIGINFO;

    sigaction(SIGSEGV, &sa, NULL);

}
