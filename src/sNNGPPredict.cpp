#ifndef R_NO_REMAP
#  define R_NO_REMAP
#endif

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

extern "C" {

  SEXP sNNGPPredict(SEXP X_r, SEXP y_r, SEXP coords_r, SEXP n_r, SEXP p_r, SEXP m_r, SEXP X0_r, SEXP coords0_r, SEXP q_r, SEXP nnIndx0_r, 
		    SEXP betaSamples_r, SEXP thetaSamples_r, SEXP wSamples_r, SEXP nSamples_r, SEXP family_r, SEXP covModel_r, SEXP nThreads_r, SEXP verbose_r, SEXP nReport_r){

    int h, i, j, k, l, s, info, nProtect=0;
    const int inc = 1;
    const double one = 1.0;
    const double negOne = -1.0;
    const double zero = 0.0;
    char const *lower = "L";
    
    //get args
    double *X = REAL(X_r);
    double *y = REAL(y_r);
    double *coords = REAL(coords_r);
    int n = INTEGER(n_r)[0];
    int p = INTEGER(p_r)[0];
    int m = INTEGER(m_r)[0];
    int mm = m*m;

    double *X0 = REAL(X0_r);
    double *coords0 = REAL(coords0_r);
    int q = INTEGER(q_r)[0];

    int *nnIndx0 = INTEGER(nnIndx0_r);        
    double *beta = REAL(betaSamples_r);
    double *theta = REAL(thetaSamples_r);
    double *w = REAL(wSamples_r);
    
    int nSamples = INTEGER(nSamples_r)[0];
    int family = INTEGER(family_r)[0];
    int covModel = INTEGER(covModel_r)[0];
    std::string corName = getCorName(covModel);
    int nThreads = INTEGER(nThreads_r)[0];
    int verbose = INTEGER(verbose_r)[0];
    int nReport = INTEGER(nReport_r)[0];
    
#ifdef _OPENMP
    omp_set_num_threads(nThreads);
#else
    if(nThreads > 1){
      Rf_warning("n.omp.threads > %i, but source not compiled with OpenMP support.", nThreads);
      nThreads = 1;
    }
#endif
    
    if(verbose){
      Rprintf("----------------------------------------\n");
      Rprintf("\tPrediction description\n");
      Rprintf("----------------------------------------\n");
      Rprintf("NNGP Latent model fit with %i observations.\n\n", n);
      Rprintf("Number of covariates %i (including intercept if specified).\n\n", p);
      Rprintf("Using the %s spatial correlation model.\n\n", corName.c_str());
      Rprintf("Using %i nearest neighbors.\n\n", m);
      Rprintf("Number of MCMC samples %i.\n\n", nSamples);
      Rprintf("Predicting at %i locations.\n\n", q);  
#ifdef _OPENMP
      Rprintf("\nSource compiled with OpenMP support and model fit using %i threads.\n", nThreads);
#else
      Rprintf("\n\nSource not compiled with OpenMP support.\n");
#endif
    } 
    
    //parameters
    int nTheta, sigmaSqIndx, tauSqIndx, phiIndx, nuIndx;

    if(family == 1){
      if(corName != "matern"){
	nTheta = 3;//sigma^2, tau^2, phi
	sigmaSqIndx = 0; tauSqIndx = 1; phiIndx = 2;
      }else{
	nTheta = 4;//sigma^2, tau^2, phi, nu
	sigmaSqIndx = 0; tauSqIndx = 1; phiIndx = 2; nuIndx = 3;
      }
    }else{//family is binomial
      if(corName != "matern"){
	nTheta = 2;//sigma^2, phi
	sigmaSqIndx = 0; phiIndx = 1;
      }else{
	nTheta = 3;//sigma^2, phi, nu
	sigmaSqIndx = 0; phiIndx = 1; nuIndx = 2;
      }
    }
    
    //get max nu
    double nuMax = 0;
    int nb = 0;
    
    if(corName == "matern"){
      for(i = 0; i < nSamples; i++){
	if(theta[i*nTheta+nuIndx] > nuMax){
	  nuMax = theta[i*nTheta+nuIndx];
	}
      }

      nb = 1+static_cast<int>(floor(nuMax));
    }

    double *bk = (double *) R_alloc(nThreads*nb, sizeof(double));
    
    double *C = (double *) R_alloc(nThreads*mm, sizeof(double)); zeros(C, nThreads*mm);
    double *c = (double *) R_alloc(nThreads*m, sizeof(double)); zeros(c, nThreads*m);
    double *tmp_m  = (double *) R_alloc(nThreads*m, sizeof(double));
    double phi = 0, nu = 0, sigmaSq = 0, tauSq = 0, d;
    int threadID = 0, status = 0;

    SEXP y0_r, w0_r;
    PROTECT(y0_r = Rf_allocMatrix(REALSXP, q, nSamples)); nProtect++; 
    PROTECT(w0_r = Rf_allocMatrix(REALSXP, q, nSamples)); nProtect++;
    double *y0 = REAL(y0_r);
    double *w0 = REAL(w0_r);
 
    if(verbose){
      Rprintf("-------------------------------------------------\n");
      Rprintf("\t\tPredicting\n");
      Rprintf("-------------------------------------------------\n");
      #ifdef Win32
        R_FlushConsole();
      #endif
    }

    int zIndx = -1;
    double *wZ = (double *) R_alloc(q*nSamples, sizeof(double));

    double *yZ = NULL;
    if(family == 1){
      yZ = (double *) R_alloc(q*nSamples, sizeof(double));
    }
    
    GetRNGstate();
    
    for(i = 0; i < q*nSamples; i++){
      wZ[i] = rnorm(0.0,1.0);
    }
    
    if(family == 1){
      for(i = 0; i < q*nSamples; i++){
	yZ[i] = rnorm(0.0,1.0);
      }
    }
    
    PutRNGstate();

    for(i = 0; i < q; i++){
#ifdef _OPENMP
#pragma omp parallel for private(threadID, phi, nu, sigmaSq, tauSq, k, l, d, info)
#endif     
      for(s = 0; s < nSamples; s++){
#ifdef _OPENMP
	threadID = omp_get_thread_num();
#endif 	
	phi = theta[s*nTheta+phiIndx];
	if(corName == "matern"){
	  nu = theta[s*nTheta+nuIndx];
	}
	sigmaSq = theta[s*nTheta+sigmaSqIndx];

	if(family == 1){
	  tauSq = theta[s*nTheta+tauSqIndx];
	}
	
	for(k = 0; k < m; k++){
	  d = dist2(coords[nnIndx0[i+q*k]], coords[n+nnIndx0[i+q*k]], coords0[i], coords0[q+i]);
	  c[threadID*m+k] = sigmaSq*spCor(d, phi, nu, covModel, &bk[threadID*nb]);
	  for(l = 0; l < m; l++){
	    d = dist2(coords[nnIndx0[i+q*k]], coords[n+nnIndx0[i+q*k]], coords[nnIndx0[i+q*l]], coords[n+nnIndx0[i+q*l]]);
	    C[threadID*mm+l*m+k] = sigmaSq*spCor(d, phi, nu, covModel, &bk[threadID*nb]);
	  }
	}

	F77_NAME(dpotrf)(lower, &m, &C[threadID*mm], &m, &info FCONE); if(info != 0){Rf_error("c++ Rf_error: dpotrf failed\n");}
	F77_NAME(dpotri)(lower, &m, &C[threadID*mm], &m, &info FCONE); if(info != 0){Rf_error("c++ Rf_error: dpotri failed\n");}

	F77_NAME(dsymv)(lower, &m, &one, &C[threadID*mm], &m, &c[threadID*m], &inc, &zero, &tmp_m[threadID*m], &inc FCONE);

	d = 0;
	for(k = 0; k < m; k++){
	  d += tmp_m[threadID*m+k]*w[s*n+nnIndx0[i+q*k]];
	}

	#ifdef _OPENMP
        #pragma omp atomic
        #endif   
	zIndx++;
	
	w0[s*q+i] = sqrt(sigmaSq - F77_NAME(ddot)(&m, &tmp_m[threadID*m], &inc, &c[threadID*m], &inc))*wZ[zIndx] + d;

	if(family == 1){
	  y0[s*q+i] = sqrt(tauSq)*yZ[zIndx] + F77_NAME(ddot)(&p, &X0[i], &q, &beta[s*p], &inc) + w0[s*q+i];
	}else{//binomial
	   y0[s*q+i] = F77_NAME(ddot)(&p, &X0[i], &q, &beta[s*p], &inc) + w0[s*q+i];
	}
	
      }
      
      if(verbose){
	if(status == nReport){
	  Rprintf("Location: %i of %i, %3.2f%%\n", i, q, 100.0*i/q);
          #ifdef Win32
	  R_FlushConsole();
          #endif
	  status = 0;
	}
      }
      status++;
      R_CheckUserInterrupt();
    }
    
    if(verbose){
      Rprintf("Location: %i of %i, %3.2f%%\n", i, q, 100.0*i/q);
      #ifdef Win32
      R_FlushConsole();
      #endif
    }

    //make return object
    SEXP result_r, resultName_r;
    int nResultListObjs = 2;

    PROTECT(result_r = Rf_allocVector(VECSXP, nResultListObjs)); nProtect++;
    PROTECT(resultName_r = Rf_allocVector(VECSXP, nResultListObjs)); nProtect++;

    SET_VECTOR_ELT(result_r, 0, y0_r);
    SET_VECTOR_ELT(resultName_r, 0, Rf_mkChar("p.y.0")); 
    
    SET_VECTOR_ELT(result_r, 1, w0_r);
    SET_VECTOR_ELT(resultName_r, 1, Rf_mkChar("p.w.0"));

    Rf_namesgets(result_r, resultName_r);
    
    //unprotect
    UNPROTECT(nProtect);
    
    return(result_r);
  
  }
}

    
