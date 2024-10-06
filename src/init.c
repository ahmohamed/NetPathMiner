#include "init.h"
#include <R_ext/Rdynload.h>

#define ENTRY(name, n)  { #name, (DL_FUNC) &name, n }
static const R_CallMethodDef callMethods[] = {

#ifdef HAVE_SBML
	ENTRY(readsbmlfile, 3),
	ENTRY(readsbml_sign, 3),
#endif
#ifdef HAVE_XML
	ENTRY(readkgmlfile, 2),
	ENTRY(readkgml_sign, 3),
#endif

	ENTRY(expand_complexes, 5),
	{NULL, NULL, 0}
};

static const R_CMethodDef cmethods[] = {
	ENTRY(corEdgeWeights, 7),
	ENTRY(hme3m_R, 17),
	ENTRY(pathMix, 9),
	{NULL, NULL, 0}
};

void
R_init_NetPathMiner(DllInfo *dll)
{
	R_useDynamicSymbols(dll, FALSE);
	R_registerRoutines(dll, cmethods, callMethods, NULL, NULL);
}

