#define USE_FC_LEN_T
#include <string>
#include "util.h"

#ifdef _OPENMP
#include <omp.h>
#endif

#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>
#include <R_ext/Linpack.h>
#include <R_ext/Lapack.h>
#include <R_ext/BLAS.h>
#ifndef FCONE
# define FCONE
#endif

//Description: update B and F.
void updateBF(double *B, double *F, double *c, double *C, double *coords, int *nnIndx, int *nnIndxLU, int n, int m, double *theta, int tauSqIndx, int sigmaSqIndx, int phiIndx, int nuIndx, int covModel, double *bk, int nb){
    
  int i, k, l;
  int info = 0;
  int inc = 1;
  double one = 1.0;
  double zero = 0.0;
  char lower = 'L';
  double nu = 0;

  if(getCorName(covModel) == "matern"){
    nu = theta[nuIndx];
  }

  //bk must be nb = 1+(int)floor(alpha) * nthread 
  int threadID = 0;
  double e;
  int mm = m*m;

#ifdef _OPENMP
#pragma omp parallel for private(k, l, info, threadID, e)
#endif
    for(i = 0; i < n; i++){
#ifdef _OPENMP
      threadID = omp_get_thread_num();
#endif
      if(i > 0){
	for(k = 0; k < nnIndxLU[n+i]; k++){
	  e = dist2(coords[i], coords[n+i], coords[nnIndx[nnIndxLU[i]+k]], coords[n+nnIndx[nnIndxLU[i]+k]]);
	  c[m*threadID+k] = theta[sigmaSqIndx]*spCor(e, theta[phiIndx], nu, covModel, &bk[threadID*nb]);
	  for(l = 0; l <= k; l++){
	    e = dist2(coords[nnIndx[nnIndxLU[i]+k]], coords[n+nnIndx[nnIndxLU[i]+k]], coords[nnIndx[nnIndxLU[i]+l]], coords[n+nnIndx[nnIndxLU[i]+l]]);  
	    C[mm*threadID+l*nnIndxLU[n+i]+k] = theta[sigmaSqIndx]*spCor(e, theta[phiIndx], nu, covModel, &bk[threadID*nb]); 
	    if(l == k){
	      C[mm*threadID+l*nnIndxLU[n+i]+k] += theta[tauSqIndx];
	    }
	  }
	}
	F77_NAME(dpotrf)(&lower, &nnIndxLU[n+i], &C[mm*threadID], &nnIndxLU[n+i], &info FCONE); if(info != 0){error("c++ error: dpotrf failed\n");}
	F77_NAME(dpotri)(&lower, &nnIndxLU[n+i], &C[mm*threadID], &nnIndxLU[n+i], &info FCONE); if(info != 0){error("c++ error: dpotri failed\n");}
	F77_NAME(dsymv)(&lower, &nnIndxLU[n+i], &one, &C[mm*threadID], &nnIndxLU[n+i], &c[m*threadID], &inc, &zero, &B[nnIndxLU[i]], &inc FCONE);
	F[i] = theta[sigmaSqIndx] - F77_NAME(ddot)(&nnIndxLU[n+i], &B[nnIndxLU[i]], &inc, &c[m*threadID], &inc) + theta[tauSqIndx];
      }else{
	B[i] = 0;
	F[i] = theta[sigmaSqIndx] + theta[tauSqIndx];
      }
    }

}

extern "C" {
  
  SEXP rNNGPReplicated(SEXP X_r, SEXP p_r, SEXP n_r, SEXP m_r, SEXP coords_r, SEXP covModel_r, SEXP nnIndx_r, SEXP nnIndxLU_r,
		       SEXP beta_r, SEXP theta_r, 
		       SEXP nSamples_r, SEXP nThreads_r, SEXP verbose_r, SEXP nReport_r){
    
    int h, i, j, k, l, s, info, nProtect=0;
    const int inc = 1;
    const double one = 1.0;
    const double negOne = -1.0;
    const double zero = 0.0;
    char const *ntran = "N";
    
    //get args
    double *X = REAL(X_r);
    int p = INTEGER(p_r)[0];
    int n = INTEGER(n_r)[0];
    int m = INTEGER(m_r)[0];
    double *coords = REAL(coords_r);
    int *nnIndx = INTEGER(nnIndx_r);
    int *nnIndxLU = INTEGER(nnIndxLU_r);
    int covModel = INTEGER(covModel_r)[0];
    std::string corName = getCorName(covModel);

    //samples
    double *beta = REAL(beta_r);//n.samples x p
    double *theta = REAL(theta_r);//make sure theta comes in transposed n.theta x n.samples
 
    int nSamples = INTEGER(nSamples_r)[0];
    int nThreads = INTEGER(nThreads_r)[0];
    int verbose = INTEGER(verbose_r)[0];
    int nReport = INTEGER(nReport_r)[0];
    
#ifdef _OPENMP
    omp_set_num_threads(nThreads);
#else
    if(nThreads > 1){
      warning("n.omp.threads > %i, but source not compiled with OpenMP support.", nThreads);
      nThreads = 1;
    }
#endif
    
    if(verbose){
      Rprintf("----------------------------------------\n");
      Rprintf("\tComputing replicates\n");
      Rprintf("----------------------------------------\n");
      Rprintf("NNGP Response model fit with %i observations.\n\n", n);
      Rprintf("Number of covariates %i (including intercept if specified).\n\n", p);
      Rprintf("Using the %s spatial correlation model.\n\n", corName.c_str());
      Rprintf("Using %i nearest neighbors.\n\n", m);
      Rprintf("Number of MCMC samples %i.\n\n", nSamples);
#ifdef _OPENMP
      Rprintf("Source compiled with OpenMP support and model fit using %i thread(s).\n", nThreads);
#else
      Rprintf("Source not compiled with OpenMP support.\n");
#endif
    } 
    
    //parameters
    int nTheta, sigmaSqIndx, tauSqIndx, phiIndx, nuIndx;
    
    if(corName != "matern"){
      nTheta = 3;//sigma^2, tau^2, phi
      sigmaSqIndx = 0; tauSqIndx = 1; phiIndx = 2;
    }else{
      nTheta = 4;//sigma^2, tau^2, phi, nu
      sigmaSqIndx = 0; tauSqIndx = 1; phiIndx = 2; nuIndx = 3;
    }    
    
    //other stuff
    int nIndx = static_cast<int>(static_cast<double>(1+m)/2*m+(n-m-1)*m);
    int mm = m*m;
    double *B = (double *) R_alloc(nIndx, sizeof(double));
    double *F = (double *) R_alloc(n, sizeof(double));
    double *c =(double *) R_alloc(m*nThreads, sizeof(double));
    double *C = (double *) R_alloc(mm*nThreads, sizeof(double));

    //return stuff  
    SEXP repSamples_r;
    PROTECT(repSamples_r = allocMatrix(REALSXP, n, nSamples)); nProtect++;

    //other stuff
    int status = 0;
    double *tmp_m = (double *) R_alloc(m, sizeof(double));
    double z;
    double *v = (double *) R_alloc(n, sizeof(double)); zeros(v, n);
    
    double nuMax = 0;
    if(getCorName(covModel) == "matern"){
      for(i = 0; i < nSamples; i++){
	if(theta[nuIndx] > nuMax){
	  nuMax = theta[nuIndx];
	}
      }
    }

    int nb = 1+static_cast<int>(floor(nuMax));
    double *bk = (double *) R_alloc(nThreads*nb, sizeof(double));
    
    if(verbose){
      Rprintf("------------\n");
      Rprintf("\t\tSampling\n");
      #ifdef Win32
        R_FlushConsole();
      #endif
    }

    GetRNGstate();
    
    for(s = 0; s < nSamples; s++){

      updateBF(B, F, c, C, coords, nnIndx, nnIndxLU, n, m, &theta[s*nTheta], tauSqIndx, sigmaSqIndx, phiIndx, nuIndx, covModel, bk, nb);
      
      for(i = 0; i < n; i++){
	z = rnorm(0.0, 1.0);
      	if(i == 0){
      	  v[i] = sqrt(F[i])*z;
      	}else{
      	  for(j = 0; j < nnIndxLU[n+i]; j++){
      	    tmp_m[j] = v[nnIndx[nnIndxLU[i]+j]];
      	  }
      	  v[i] = F77_NAME(ddot)(&nnIndxLU[n+i], &B[nnIndxLU[i]], &inc, tmp_m, &inc)+sqrt(F[i])*z;
      	}
      }
	
      F77_NAME(dcopy)(&n, v, &inc, &REAL(repSamples_r)[s*n], &inc);
      F77_NAME(dgemv)(ntran, &n, &p, &one, X, &n, &beta[s], &nSamples, &one, &REAL(repSamples_r)[s*n], &inc FCONE);
            
      //report
      if(status == nReport){
      	if(verbose){
      	  Rprintf("Sampled: %i of %i, %3.2f%%\n", s, nSamples, 100.0*s/nSamples);
          #ifdef Win32
      	  R_FlushConsole();
          #endif
      	}
      	status = 0;
      }
      
      status++;
      
      R_CheckUserInterrupt();
    }

    if(verbose){
      Rprintf("Sampled: %i of %i, %3.2f%%\n", nSamples, nSamples, 100.0);
      #ifdef Win32
      R_FlushConsole();
      #endif
    }
    
    PutRNGstate();
    
    //make return object
    SEXP result_r, resultName_r;
    int nResultListObjs = 1;
    
    PROTECT(result_r = allocVector(VECSXP, nResultListObjs)); nProtect++;
    PROTECT(resultName_r = allocVector(VECSXP, nResultListObjs)); nProtect++;

    SET_VECTOR_ELT(result_r, 0, repSamples_r);
    SET_VECTOR_ELT(resultName_r, 0, mkChar("y.rep.samples")); 
        
    namesgets(result_r, resultName_r);
    
    //unprotect
    UNPROTECT(nProtect);
    
    return(result_r);
  }
}
