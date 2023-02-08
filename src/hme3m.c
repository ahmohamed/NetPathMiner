#include "hme3m.h"

void hme3m_R(double *Y,
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
	double *LIKELIHOOD) 
{
	size_t nx = (size_t)(*NX);
	size_t nobs = (size_t)(*NOBS);
	int m = (int)(*M);
	double alpha = (double)(*ALPHA);
	double lambda = (double)(*LAMBDA);

	// Run the model
	hme3m(Y,
		X,
		m,
		lambda,
		alpha,
		nobs,
		nx,
		HME3MITER,
		PLRITER,
		H,
		PATHPROBS,
		PLRPRE,
		THETA,
		BETA,
		PROPORTIONS,
		HMEPRE,
		LIKELIHOOD);
}

void hme3m(double * Y,
	double * X,
	int m,
	double lambda,
	double alpha,
	size_t nobs,
	size_t nx,
	int * HME3MITER,
	int * PLRITER,
	double * H,
	double * PATHPROBS,
	double * PLRPRE,
	double * THETA,
	double * BETA,
	double * PROPORTIONS,
	double * HMEPRE,
	double * LIKELIHOOD)
{
	int i,j,k,iter;	
	int CONVERGED = 0;
	double tempval = 0.0;
	double tempval2 = 0.0;

	// temporary allocations
	double * mbeta;
	MALLOC(mbeta,sizeof(double)*nx);
	
	double * mypre; 
	MALLOC(mypre,sizeof(double)*nobs);
	
	double * mw;
	MALLOC(mw,sizeof(double)*nobs);
	
	iter = 0;
	while (CONVERGED == 0) {  
/*----------------------------------------------------
                    E-STEP
------------------------------------------------------*/
		for (i = 0;i < nobs;i = i + 1) {
			tempval = 0.0;
			for (k = 0;k < m;k = k + 1) {
                //tempval += pk*pmx*fits
                tempval = tempval + PROPORTIONS[k]*PATHPROBS[nobs*k + i]*PLRPRE[nobs*k + i]; 
			}

			// update the resposibilities
			for (k = 0;k < m;k = k + 1) {
        		// H = pk*pmx*fix/sum(pk*pmx*fits)
				H[nobs*k + i] = (PROPORTIONS[k]*PATHPROBS[nobs*k + i]*PLRPRE[nobs*k + i]) / tempval; 
            }
        }

/*----------------------------------------------------
                       M-STEP
------------------------------------------------------*/
		tempval2 = 0.0;
		for (k = 0;k < m;k = k + 1) { // for each component
			tempval = 0.0;
			for (i = 0;i < nobs;i = i + 1) tempval = tempval + H[k*nobs + i]; 

            // Estimate new pi_k = sum(pk[k]*pmx[k]*fits[k])
			PROPORTIONS[k] = tempval; 
       		tempval2 = tempval2 + PROPORTIONS[k];
			
			// Estimate a new theta
			for (j = 0;j < nx;j = j + 1) { // for each column
				tempval = 0.0;
				for (i = 0;i < nobs; i = i + 1) { // sum over the rows
					if (X[j*nobs + i] == 1.0) tempval = tempval + H[k*nobs + i];
				}
				THETA[k*nx + j] = tempval / PROPORTIONS[k];
			}       
            
			// Estimate the probability for each path beloning to that component
			for (i = 0;i < nobs;i = i + 1) { // for each row
				tempval = 1.0;
				for (j = 0;j < nx; j = j + 1) { // multiply all thetas along a row to get the probability
					if (X[j*nobs + i] == 1.0) tempval = tempval * THETA[k*nx + j];
				}
				PATHPROBS[k*nobs + i] = tempval;
			}

			// Estimate new PLR models
			for (i = 0;i < nx;i = i + 1) mbeta[i] = 0;
			for (i = 0;i < nobs;i = i + 1) {
				mypre[i] = 0.5;
				mw[i] = H[k*nobs + i];
			}
			
            irls(Y,X,nobs,nx,mw,mbeta,mypre,lambda,alpha,(int)(*PLRITER));
			
            for (i = 0;i < nobs;i = i + 1) PLRPRE[k*nobs + i] = mypre[i];
			for (i = 0;i < nx;i = i + 1) BETA[k*nx + i] = mbeta[i];
		}
 
        // normalize 3M the new mixture proportions
        for (k = 0;k < m;k = k + 1) PROPORTIONS[k] = PROPORTIONS[k]/tempval2;

/*----------------------------------------------------
          Prediction and Convergence Testing
------------------------------------------------------*/
        LIKELIHOOD[iter] = 0.0;
        for (i = 0;i < nobs;i = i + 1) {
            tempval = 0.0;
            tempval2 = 0.0;
            HMEPRE[i] = 0.0;
            for (k = 0;k < m;k = k + 1) {
                //tempval += pk*pmx*fits
                tempval = tempval + H[k*nobs + i]*PROPORTIONS[k]*PATHPROBS[nobs*k + i] * PLRPRE[nobs*k + i]; 

                //tempval2 += pmx
                tempval2 = tempval2 + PATHPROBS[nobs*k + i];
                HMEPRE[i] = HMEPRE[i] + PATHPROBS[nobs*k + i]*PLRPRE[nobs*k + i];
            }

            // update the likelihood
            LIKELIHOOD[iter] = LIKELIHOOD[iter] + log(tempval);  // l = sum(log(pk*pmx*fits))

            // update the predictions
            HMEPRE[i] = HMEPRE[i]/tempval2;
        }

        if (iter > 0) {
            if (fabs(LIKELIHOOD[iter] - LIKELIHOOD[iter-1]) < 0.001 || iter > (int)(*HME3MITER)-1) {
                CONVERGED = 1;
                *HME3MITER = iter + 1;
            }
        }
        iter = iter + 1;
    }

    // update responsibilities for the final time.
    for (i = 0;i < nobs;i = i + 1) {
        tempval = 0.0;
        for (k = 0;k < m;k = k + 1) {
            //tempval += pk*pmx*fits
            tempval = tempval + PROPORTIONS[k]*PATHPROBS[nobs*k + i]*PLRPRE[nobs*k + i]; 
        }
        // update the resposibilities
        for (k = 0;k < m;k = k + 1) {
            // H = pk*pmx*fix/sum(pk*pmx*fits)
            H[nobs*k + i] = (PROPORTIONS[k]*PATHPROBS[nobs*k + i]*PLRPRE[nobs*k + i]) / tempval; 
        }
    }

  	//clean up
	free(mbeta);
	free(mypre);
	free(mw);

	return;
}

void pathMix(int *X,
    int *M,
    int *NOBS,
    int *NX,
    int *ITER,
    double *H,
    double *THETA,
    double *PROPORTIONS,
    double *LIKELIHOOD) 
{
  int CONVERGED = 0;
  int iter = 0;
  int nrow = (int)(*NOBS);
  int ncol = (int)(*NX);
  int m = (int)(*M);
  double tempval = 0.0, tempval2 = 0.0;
  double loglikelihood = 0.0;

  while (CONVERGED == 0) {
    /*-----------------------------------------
         E Step: Compute the responsiblities
    ------------------------------------------*/ 
    for (int i = 0;i < nrow;i = i + 1) {
      
      // for each mixture component   
      tempval = 0.0;
      for (int k = 0;k < m;k = k + 1) {
        // Find the probability of a path
        H[k*nrow + i] = PROPORTIONS[k];
        for (int j = 0;j < ncol;j = j + 1) {
            if (X[j*nrow + i] == 1) H[k*nrow + i] = H[k*nrow + i] * THETA[k*ncol + j];
        }
        // get the sum for that row 
        tempval = tempval + H[k*nrow + i];
      }

      // normalize the responsibilities
      for (int k = 0;k < m;k = k + 1) H[k*nrow + i] = H[k*nrow + i]/tempval;
    } 

    /*-----------------------------------------
       M Step: Update the Path Probabilities
    ------------------------------------------*/ 
    // for each component
    tempval2 = 0.0;
    for (int k = 0;k < m;k = k + 1) {

      // Update the new mixture proportions
      PROPORTIONS[k] = 0.0;
      for (int i = 0;i < nrow;i = i + 1) PROPORTIONS[k] = PROPORTIONS[k] + H[k*nrow + i];
      tempval2 = tempval2 + PROPORTIONS[k];
    
      // for each transition
      for (int j = 0;j < ncol;j = j + 1) {

        // sum the responsibilities where X[i][j] == 1
        tempval = 0.0;
        for (int i = 0;i < nrow;i = i + 1) {
          if (X[j*nrow + i] == 1.0) tempval = tempval + H[k*nrow + i];
        }

        // Update the new transition probability theta[k][j]
        THETA[k*ncol + j] = tempval/PROPORTIONS[k];
      }      
    }
    // Normalize the mixture proportions
    for (int k = 0;k < m;k = k + 1) PROPORTIONS[k] = PROPORTIONS[k]/tempval2;
  
    /*-----------------------------------------
             Check for Convergence 
    ------------------------------------------*/   
    loglikelihood = 0.0;
    for (int i = 0;i < nrow;i = i + 1) {
      tempval = 0.0;
      // for each mixture component    
      for (int k = 0;k < m;k = k + 1) {
        // Find the probability of a path
        tempval2 = H[k*nrow + i] * PROPORTIONS[k];
        for (int j = 0;j < ncol;j = j + 1) {
          if (X[j*nrow + i] == 1)  tempval2 = tempval2 * THETA[k*ncol + j];
        }
        // get the sum of each mixture component for that row 
        tempval = tempval + tempval2;
      }
      // update the log-likelihood
      loglikelihood = loglikelihood + log(tempval);
    } 
    // Store the log-likelihood
    LIKELIHOOD[iter] = loglikelihood;


    if (iter > 0) {
        if (fabs(LIKELIHOOD[iter] - LIKELIHOOD[iter-1]) < 0.001 || iter > (int)(*ITER)-1) {
            CONVERGED = 1;
            *ITER = iter + 1;
        }
    }   
    iter = iter + 1;
  }
}

void irls(double *y, 
	double *x,
	int nobs,
	int nx,
	double *w,
	double *beta,
	double *ypre,	
	double lambda,
	double alpha,
	int maxiter) 
{
	int i,k,iter;
	int NotConverged = 1;
	int fclen = 1;
	double likelihood = 0.0;
	double NEWlikelihood = 0.0;
	double tempval = 0.0;
	double tempval2 = 0.0;
	
	char *transpose = "T", *dont_transpose = "N";
	//char *upper_triangle = "U", *left_side = "L";
	//char *non_unit = "N";
	double double_zero = 0.0, double_one = 1.0;
	int int_one = 1,info = 0;
	
	double * cov;
	MALLOC(cov,sizeof(double)*nx*nx);  // Covariance matrix
	double * weights;  
	MALLOC(weights,sizeof(double)*nobs); // HME3M weights
	double * XW;
	MALLOC(XW,sizeof(double)*nobs*nx); // temporary X*weights
	double * rss;
	MALLOC(rss,sizeof(double)*nobs); // residual sums of squares vectory
	double * bret;
	MALLOC(bret,sizeof(double)*nx); // temporary beta information vector
    double * btemp;
    MALLOC(btemp,sizeof(double)*nx);
    
    int *ipiv;
    MALLOC(ipiv,sizeof(int)*nx*nx);

    double *work;
    MALLOC(work,sizeof(double)*100*nx);  // Following IBMs recommendations
    int lwork = 100*nx;

	//Compute the initial likelihood ypre = XB
	likelihood = 0.0;
	
	// ypre = X %*% B
	F77_NAME(dgemv)(dont_transpose, &nobs, &nx, &double_one, x, &nobs, beta, &int_one, &double_zero, ypre, &int_one, fclen);
    
	for (i = 0;i < nobs;i = i + 1) {
		tempval = ypre[i]; 
		tempval2 = 1/(1+exp(-ypre[i])); 
		ypre[i] = tempval2;

		likelihood = likelihood + y[i]*tempval - log(1 + exp(tempval));

		weights[i] = w[i]*tempval2*(1-tempval2); // HME3M weights = 3M weights * IRLS weights
	}
	
	iter = 0;
	while (NotConverged == 1) {
		for (i = 0;i < nobs;i = i + 1) {
			for (k = 0;k < nx;k = k + 1) {
				XW[k*nobs + i] = x[k*nobs + i] * weights[i];
			}
		}
					
		// t(X) %*% W %*% X
		F77_NAME(dgemm)(transpose, dont_transpose,&nx,&nx,&nobs,&double_one, XW, &nobs, x, &nobs, &double_zero, cov, &nx); 

		//cov = t(X) %*% W %*% X + L
		for (i = 0;i < nx;i = i + 1) cov[i*nx + i] = cov[i*nx + i] + lambda;

		// Compute w*(y - ypre)
		for (i = 0;i < nobs;i = i + 1) rss[i] = w[i]*(y[i] - ypre[i]);
			
		// btemp = t(x)*(w*(y - ypre))
		F77_NAME(dgemv)(transpose, &nobs, &nx, &double_one, x, &nobs, rss, &int_one, &double_zero, btemp, &int_one, fclen);

    // (t(X) %*% W %*% X)^(-1)
    F77_CALL(dgetrf)(&nx,&nx,cov,&nx,ipiv,&info);
    F77_NAME(dgetri)(&nx,cov,&nx,ipiv,work,&lwork,&info);

    // (t(X) %*% W %*% X)^(-1) %*% t(X) %*% W * (y - ypre)
    F77_NAME(dgemv)(dont_transpose, &nx, &nx, &double_one, cov, &nx, btemp, &int_one, &double_zero, bret, &int_one, fclen);

		// beta = beta + bupdate
		for (i = 0;i < nx;i = i + 1) beta[i] = beta[i] + alpha*bret[i];

		// ypre = XB
		F77_NAME(dgemv)(dont_transpose, &nobs, &nx, &double_one, x, &nobs, beta, &int_one, &double_zero, ypre, &int_one, fclen);
		
		//Compute the new likelihood
		NEWlikelihood = 0.0;
		for (i = 0;i < nobs;i = i + 1) {
			tempval = ypre[i];
			tempval2 = 1/(1+exp(-ypre[i]));
			ypre[i] = tempval2;

			NEWlikelihood = NEWlikelihood + y[i]*tempval - log(1 + exp(tempval));
			weights[i] = w[i]*tempval2*(1-tempval2); // HME3M weights = 3M weights * IRLS weights
		}
		if (fabs(likelihood - NEWlikelihood) < 0.01 || iter > maxiter-1) NotConverged = 0;
		likelihood = NEWlikelihood;
		iter = iter+1;
	}

	free(cov);
	free(weights);
	free(XW);
	free(bret);
	free(rss);
    free(ipiv);
    free(work);
}
