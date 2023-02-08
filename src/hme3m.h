#ifndef __hme3m__h_
#define __hme3m__h_
#include "init.h"
#ifndef  USE_FC_LEN_T
# define USE_FC_LEN_T
#endif
#include <R_ext/BLAS.h>
#include <R_ext/Lapack.h>
#ifndef FCONE
# define FCONE
#endif

#define oops(s) { perror((s)); error("Failed to allocate memory"); }
#define MALLOC(s,t) if(((s) = malloc(t)) == NULL) { oops("error: malloc() "); }

void hme3m(double * y,
	double * x,
	int m,
	double lambda,
	double alpha,
	size_t nobs,
	size_t nx,
	int * HME3MITER,
	int * PLRITER,
	double * h,
	double * pathprobs,
	double * plrpre,
	double * theta,
	double * beta ,
	double * PROPORTIONS,
	double * HMEPRE,
	double * LIKELIHOOD);
	
void irls(double *y, 
	double *x,
	int nobs,
	int nx,
	double *w,
	double *beta,
	double *ypre,	
	double lambda,
	double alpha,
	int maxiter);

#endif
