#ifndef __hme3m__h_
#define __hme3m__h_
#include <stdlib.h>
#include <stdio.h>
#include <R.h>
#include <Rmath.h>
#include <R_ext/BLAS.h>
#include <R_ext/Lapack.h>

#include <Rdefines.h>
#include <Rinternals.h>


#define oops(s) { perror((s)); error("Failed to allocate memory"); }
#define MALLOC(s,t) if(((s) = malloc(t)) == NULL) { oops("error: malloc() "); }

extern void hme3m_R(double *Y,
	double *X,
	int *M,
	double *LAMBDA,
	double *ALPHA,
	int *NOBS,
	int *NX,
	int *HME3MITER,
	int *PLRITER,
	double *H,
	double *PATHPROBS,
	double *PLRPRE,
	double *THETA,
	double *BETA,
	double *PROPORTIONS,
	double *HMEPRE,
	double *LIKELIHOOD);

extern void pathMix(int *X,
    int *M,
    int *NOBS,
    int *NX,
    int *ITER,
    double *H,
    double *THETA,
    double *PROPORTIONS,
    double *LIKELIHOOD);


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
