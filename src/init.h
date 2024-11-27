#ifndef __init__h_
#define __init__h_

#include <stdio.h>
#include <stdlib.h>

#ifdef __cplusplus
#include <string.h>
#include <vector>
#include <iostream>
#include <algorithm>
using namespace std;
#endif //__cplusplus

#ifndef R_NO_REMAP
#define R_NO_REMAP
#endif

#include <R.h>
#include <Rmath.h>
#include <Rdefines.h>
#include <Rinternals.h>



void hme3m_R(double *Y, double *X, int *M, double *LAMBDA, double *ALPHA, int *NOBS,
			int *NX, int *HME3MITER, int *PLRITER, double *H, double *PATHPROBS,
			double *PLRPRE, double *THETA, double *BETA, double *PROPORTIONS,
			double *HMEPRE, double *LIKELIHOOD);

void pathMix(int *X, int *M, int *NOBS, int *NX, int *ITER, double *H,
			double *THETA, double *PROPORTIONS, double *LIKELIHOOD);



#ifdef __cplusplus
extern "C" {
#endif

#ifdef HAVE_SBML
	SEXP readsbmlfile(SEXP FILENAME, SEXP ATTR_TERMS, SEXP VERBOSE);
	SEXP readsbml_sign(SEXP FILENAME, SEXP ATTR_TERMS, SEXP VERBOSE);
#endif
#ifdef HAVE_XML
	SEXP readkgmlfile(SEXP FILENAME, SEXP VERBOSE);
	SEXP readkgml_sign(SEXP FILENAME, SEXP EXPAND_COMPLEXES, SEXP VERBOSE);
#endif

SEXP expand_complexes(SEXP ATTR_LS, SEXP EL, SEXP V, SEXP EXPAND, SEXP MISSING);
SEXP pathranker(SEXP node_list, SEXP edge_list, SEXP edge_weights, SEXP rk, SEXP minpathsize);
SEXP scope(SEXP node_list, SEXP edge_list, SEXP edge_weights, SEXP SAMPLEDPATHS, SEXP ALPHA, SEXP ECHO);
SEXP samplepaths(SEXP node_list, SEXP edge_list, SEXP edge_weights, SEXP MAXPATHLENGTH,
				SEXP SAMPLEPATHS, SEXP WARMUPSTEPS);

void corEdgeWeights(double * X, int * EDGELIST, int * SAMEGENE,	double * WEIGHT,
					int *NEDGES, int * NOBS, int * NCOR);

#ifdef __cplusplus
}
#endif



#endif
