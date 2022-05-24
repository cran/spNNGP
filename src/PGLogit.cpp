#define USE_FC_LEN_T
#include <string>
#include "util.h"
#include "rpg.h"

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

extern "C" {
  
  SEXP PGLogit(SEXP y_r, SEXP X_r, SEXP p_r, SEXP n_r, SEXP nTrial_r, SEXP betaStarting_r, SEXP nSamples_r, SEXP nThreads_r, SEXP verbose_r){
    
    int i, j, s, info, nProtect=0;
    const int inc = 1;
    const double one = 1.0;
    const double negOne = -1.0;
    const double zero = 0.0;
    char const *lower = "L";
    char const *ntran = "N";
    char const *ytran = "T";

    
    //get args
    double *y = REAL(y_r);
    double *X = REAL(X_r);
    int p = INTEGER(p_r)[0];
    int n = INTEGER(n_r)[0];
    int *nTrial = INTEGER(nTrial_r);
    int nSamples = INTEGER(nSamples_r)[0];
    int nThreads = INTEGER(nThreads_r)[0];
    int verbose = INTEGER(verbose_r)[0];
 
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
      Rprintf("\tModel description\n");
      Rprintf("----------------------------------------\n");
      Rprintf("Logistic regression with Polya-Gamma latent\nvariable fit with %i observations.\n\n", n);
      Rprintf("Number of MCMC samples %i.\n\n", nSamples);
#ifdef _OPENMP
      Rprintf("\nSource compiled with OpenMP support and model fit using %i thread(s).\n\n", nThreads);
#else
      Rprintf("Source not compiled with OpenMP support.\n\n");
#endif
      Rprintf("Sampling ... \n");
    } 
    
    //parameters	
    double *beta = (double *) R_alloc(p, sizeof(double));    
    F77_NAME(dcopy)(&p, REAL(betaStarting_r), &inc, beta, &inc);
    double *omega = (double *) R_alloc(n, sizeof(double));
    double *kappa = (double *) R_alloc(n, sizeof(double));
    double *yStr = (double *) R_alloc(n, sizeof(double));

    //return stuff  
    SEXP betaSamples_r;
    PROTECT(betaSamples_r = allocMatrix(REALSXP, p, nSamples)); nProtect++;
    
    //other stuff
    int pp = p*p;
    int np = n*p;
    double *tmp_pp = (double *) R_alloc(pp, sizeof(double));
    double *tmp_p = (double *) R_alloc(p, sizeof(double));
    double *tmp_p2 = (double *) R_alloc(p, sizeof(double));
    double *tmp_n = (double *) R_alloc(n, sizeof(double));
    double *tmp_np = (double *) R_alloc(np, sizeof(double));
	
    for(i = 0; i < n; i++){
      kappa[i] = y[i] - static_cast<double>(nTrial[i])/2.0;
    }
    
    GetRNGstate();
    
    for(s = 0; s < nSamples; s++){

      //update augs
      for(i = 0; i < n; i++){
	omega[i] = rpg(nTrial[i], F77_NAME(ddot)(&p, &X[i], &n, beta, &inc));
	yStr[i] = kappa[i]/omega[i];
      }
           
      //update beta 
      for(i = 0; i < n; i++){
	tmp_n[i] = yStr[i]*omega[i];
      }
      F77_NAME(dgemv)(ytran, &n, &p, &one, X, &n, tmp_n, &inc, &zero, tmp_p, &inc FCONE); 	  

      for(i = 0; i < n; i++){
	for(j = 0; j < p; j++){
	  tmp_np[j*n+i] = X[j*n+i]*omega[i];
	}
      }

      F77_NAME(dgemm)(ytran, ntran, &p, &p, &n, &one, X, &n, tmp_np, &n, &zero, tmp_pp, &p FCONE FCONE);

      F77_NAME(dpotrf)(lower, &p, tmp_pp, &p, &info FCONE); if(info != 0){error("c++ error: dpotrf here failed\n");}
      F77_NAME(dpotri)(lower, &p, tmp_pp, &p, &info FCONE); if(info != 0){error("c++ error: dpotri here failed\n");}
      F77_NAME(dsymv)(lower, &p, &one, tmp_pp, &p, tmp_p, &inc, &zero, tmp_p2, &inc FCONE);
      F77_NAME(dpotrf)(lower, &p, tmp_pp, &p, &info FCONE); if(info != 0){error("c++ error: dpotrf here failed\n");}
      mvrnorm(beta, tmp_p2, tmp_pp, p);

      //save samples  
      F77_NAME(dcopy)(&p, beta, &inc, &REAL(betaSamples_r)[s*p], &inc);

      R_CheckUserInterrupt();
    }
    
    PutRNGstate();

    //make return object
    SEXP result_r, resultName_r;
    int nResultListObjs = 1;

    PROTECT(result_r = allocVector(VECSXP, nResultListObjs)); nProtect++;
    PROTECT(resultName_r = allocVector(VECSXP, nResultListObjs)); nProtect++;

    SET_VECTOR_ELT(result_r, 0, betaSamples_r);
    SET_VECTOR_ELT(resultName_r, 0, mkChar("p.beta.samples")); 
    
    namesgets(result_r, resultName_r);
    
    //unprotect
    UNPROTECT(nProtect);
    
    return(result_r);
  }
}
