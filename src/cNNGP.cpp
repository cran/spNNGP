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
void updateConjBF(double *B, double *F, double *c, double *C, double *coords, int *nnIndx, int *nnIndxLU, int n, int m, double phi, double alpha, double nu, int covModel, double *bk, double nuMax){
    
  int i, k, l;
  int info = 0;
  int inc = 1;
  double one = 1.0;
  double zero = 0.0;
  char lower = 'L';
   
  //bk must be 1+(int)floor(alpha) * nthread
  int nb = 1+static_cast<int>(floor(nuMax));
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
	  c[m*threadID+k] = spCor(e, phi, nu, covModel, &bk[threadID*nb]);
	  for(l = 0; l <= k; l++){
	    e = dist2(coords[nnIndx[nnIndxLU[i]+k]], coords[n+nnIndx[nnIndxLU[i]+k]], coords[nnIndx[nnIndxLU[i]+l]], coords[n+nnIndx[nnIndxLU[i]+l]]);
	    C[mm*threadID+l*nnIndxLU[n+i]+k] = spCor(e, phi, nu, covModel, &bk[threadID*nb]); 
	    if(l == k){
	      C[mm*threadID+l*nnIndxLU[n+i]+k] += alpha;
	    }
	  }
	}
	
	F77_NAME(dpotrf)(&lower, &nnIndxLU[n+i], &C[mm*threadID], &nnIndxLU[n+i], &info FCONE); if(info != 0){error("c++ error: dpotrf failed\n");}
	F77_NAME(dpotri)(&lower, &nnIndxLU[n+i], &C[mm*threadID], &nnIndxLU[n+i], &info FCONE); if(info != 0){error("c++ error: dpotri failed\n");}
	F77_NAME(dsymv)(&lower, &nnIndxLU[n+i], &one, &C[mm*threadID], &nnIndxLU[n+i], &c[m*threadID], &inc, &zero, &B[nnIndxLU[i]], &inc FCONE);
	F[i] = 1.0 + alpha - F77_NAME(ddot)(&nnIndxLU[n+i], &B[nnIndxLU[i]], &inc, &c[m*threadID], &inc);
      }else{
	B[i] = 0;
	F[i] = 1.0 + alpha;
      }
    }

}

extern "C" {
  
  SEXP cNNGP(SEXP y_r, SEXP X_r, SEXP p_r, SEXP n_r, SEXP coords_r, SEXP thetaAlpha_r, SEXP nnIndx_r, SEXP nnIndxLU_r,
	     SEXP X0_r, SEXP coords0_r, SEXP n0_r, SEXP nnIndx0_r, SEXP g_r, 
	     SEXP m_r, SEXP sigmaSqIG_r, SEXP covModel_r, SEXP nThreads_r, SEXP verbose_r){
    
    int h, i, j, k, l, s, info, nProtect=0;
    const int inc = 1;
    const double one = 1.0;
    const double negOne = -1.0;
    const double zero = 0.0;
    char const *lower = "L";
    char const *ytran = "T";
    
    //get args
    double *y = REAL(y_r);
    double *X = REAL(X_r);
    int p = INTEGER(p_r)[0];
    int n = INTEGER(n_r)[0];
    double *coords = REAL(coords_r);
    int *nnIndx = INTEGER(nnIndx_r);
    int *nnIndxLU = INTEGER(nnIndxLU_r);
    
    double *X0 = REAL(X0_r);
    int n0 = INTEGER(n0_r)[0];
    double *coords0 = REAL(coords0_r);
    int *nnIndx0 = INTEGER(nnIndx0_r);
    int g = INTEGER(g_r)[0];
    double *thetaAlpha = REAL(thetaAlpha_r); //nTheta x g

    int m = INTEGER(m_r)[0];
    double sigmaSqIGa = REAL(sigmaSqIG_r)[0]; double sigmaSqIGb = REAL(sigmaSqIG_r)[1];
    int covModel = INTEGER(covModel_r)[0];
    std::string corName = getCorName(covModel);
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
      Rprintf("NNGP Conjugate model fit with %i observations.\n\n", n);
      Rprintf("Number of covariates %i (including intercept if specified).\n\n", p);
      Rprintf("Using the %s spatial correlation model.\n\n", corName.c_str());
      Rprintf("Using %i nearest neighbors.\n\n", m);
#ifdef _OPENMP
      Rprintf("Source compiled with OpenMP support and model fit using %i thread(s).\n", nThreads);
#else
      Rprintf("Source not compiled with OpenMP support.\n");
#endif
      Rprintf("------------\n");
      Rprintf("Priors and hyperpriors:\n");
      Rprintf("\tbeta flat.\n");
      Rprintf("\tsigma.sq IG hyperpriors shape=%.5f and scale=%.5f\n", sigmaSqIGa, sigmaSqIGb);
      // Rprintf("------------\n");
      // if(corName == "matern"){
      // 	Rprintf("Considering %i set(s) of phi, nu, and alpha.\n", g);
      // }else{
      // 	Rprintf("Considering %i set(s) of phi and alpha.\n", g);
      // }
      if(n0 > 0){
	Rprintf("------------\n");
	Rprintf("Predicting at %i locations.\n", n0);
      }
    } 
    
    //parameters
    int nTheta = 2; //phi, alpha
    
    if(corName == "matern"){
      nTheta++;//nu
    }
    
    //other stuff
    int nIndx = static_cast<int>(static_cast<double>(1+m)/2*m+(n-m-1)*m);
    int mm = m*m;
    double *B = (double *) R_alloc(nIndx, sizeof(double));
    double *F = (double *) R_alloc(n, sizeof(double));
    double *c =(double *) R_alloc(m*nThreads, sizeof(double));   
    double *C = (double *) R_alloc(mm*nThreads, sizeof(double));
   
    //prediction
    double *tmp_m = (double *) R_alloc(nThreads*m, sizeof(double));
    double *C0 = (double *) R_alloc(nThreads*mm, sizeof(double));
    double *c0 = (double *) R_alloc(nThreads*m, sizeof(double));
    double *w = (double *) R_alloc(nThreads*m, sizeof(double));
    
    //return stuff
    int pp = p*p;
    SEXP beta_r, ab_r, bBInv_r, bb_r, y0Hat_r, y0HatVar_r;
    PROTECT(beta_r = allocMatrix(REALSXP, p, g)); nProtect++;
    PROTECT(ab_r = allocMatrix(REALSXP, nTheta, g)); nProtect++;
    PROTECT(bBInv_r = allocMatrix(REALSXP, pp, g)); nProtect++;
    PROTECT(bb_r = allocMatrix(REALSXP, p, g)); nProtect++;
    PROTECT(y0Hat_r = allocMatrix(REALSXP, n0, g)); nProtect++;
    PROTECT(y0HatVar_r = allocMatrix(REALSXP, n0, g)); nProtect++;

    double *beta = REAL(beta_r);
    double *ab = REAL(ab_r);
    double *y0Hat = REAL(y0Hat_r);
    double *y0HatVar = REAL(y0HatVar_r);
    
    //other stuff
    int mp = m*p;
    double *tmp_pp = (double *) R_alloc(pp, sizeof(double));
    double *tmp_p = (double *) R_alloc(nThreads*p, sizeof(double));
    double *tmp_p2 = (double *) R_alloc(nThreads*p, sizeof(double));
    double *tmp_mp = (double *) R_alloc(nThreads*mp, sizeof(double));
    double phi = 0, nu = 0, alpha = 0, e;
    int threadID = 0;
    
    //get max nu
    double *bk = NULL;
    double nuMax = 0;
    int nb = 0;
     
    if(corName == "matern"){
      for(i = 0; i < g; i++){
	if(thetaAlpha[i*nTheta+2] > nuMax){
	  nuMax = thetaAlpha[i*nTheta+2];
	}
      }
      
      nb = 1+static_cast<int>(floor(nuMax));
      bk = (double *) R_alloc(nThreads*nb, sizeof(double));
    }
    
    if(verbose){
      Rprintf("------------\n");
      Rprintf("\tEstimation for parameter set(s)\n");
      #ifdef Win32
        R_FlushConsole();
      #endif
    }
    
    for(k = 0; k < g; k++){

      phi = thetaAlpha[k*nTheta];
      alpha = thetaAlpha[k*nTheta+1];
      if(covModel == 2){//matern
	nu = thetaAlpha[k*nTheta+2];
      }
            
      //update B and F
      updateConjBF(B, F, c, C, coords, nnIndx, nnIndxLU, n, m, phi, alpha, nu, covModel, bk, nuMax);

      //estimation
      for(i = 0; i < p; i++){
	tmp_p[i] = Q(B, F, &X[n*i], y, n, nnIndx, nnIndxLU); //X'C^{-1}y
	for(j = 0; j <= i; j++){
	  tmp_pp[j*p+i] = Q(B, F, &X[n*j], &X[n*i], n, nnIndx, nnIndxLU);//X'C^{-1}X
	}
      }

      //V = inv(B), V = inv(tmp_pp), V = tmp_pp after dpotri
      F77_NAME(dpotrf)(lower, &p, tmp_pp, &p, &info FCONE); if(info != 0){error("c++ error: dpotrf failed\n");}
      F77_NAME(dpotri)(lower, &p, tmp_pp, &p, &info FCONE); if(info != 0){error("c++ error: dpotri failed\n");}

      //g = solve(B, v), g = solve(tmp_pp, tmp_p), g = beta
      F77_NAME(dsymv)(lower, &p, &one, tmp_pp, &p, tmp_p, &inc, &zero, &beta[k*p], &inc FCONE);

      //a 
      ab[k*2] = sigmaSqIGa + 1.0*n/2.0;
      
      //b
      ab[k*2+1] = sigmaSqIGb + 0.5*(Q(B, F, y, y, n, nnIndx, nnIndxLU) - F77_NAME(ddot)(&p, &beta[k*p], &inc, tmp_p, &inc));

      //save B^{-1} to calculate bB^{−1}/(a − 1) on return also keep \bb for exact sampling
      F77_NAME(dcopy)(&pp, tmp_pp, &inc, &REAL(bBInv_r)[k*pp], &inc);
      F77_NAME(dcopy)(&p, tmp_p, &inc, &REAL(bb_r)[k*p], &inc);
            
      //prediction
#ifdef _OPENMP
#pragma omp parallel for private(j, l, info, threadID, e)
#endif
      for(i = 0; i < n0; i++){
#ifdef _OPENMP   
 	threadID = omp_get_thread_num();
#endif
	//make M
	for(j = 0; j < m; j++){
	  e = dist2(coords0[i], coords0[n0+i], coords[nnIndx0[j*n0+i]], coords[n+nnIndx0[j*n0+i]]);
	  c0[m*threadID+j] = spCor(e, phi, nu, covModel, &bk[threadID*nb]);
	  for(l = 0; l <= j; l++){
	    e = dist2(coords[nnIndx0[j*n0+i]], coords[n+nnIndx0[j*n0+i]], coords[nnIndx0[l*n0+i]], coords[n+nnIndx0[l*n0+i]]);
	    C0[mm*threadID+l*m+j] = spCor(e, phi, nu, covModel, &bk[threadID*nb]);
	    if(l == j){
	      C0[mm*threadID+l*m+j] += alpha;
	    }
	  }
	}

	//make w
	F77_NAME(dpotrf)(lower, &m, &C0[mm*threadID], &m, &info FCONE); if(info != 0){error("c++ error: dpotrf failed\n");}
	F77_NAME(dpotri)(lower, &m, &C0[mm*threadID], &m, &info FCONE); if(info != 0){error("c++ error: dpotri failed\n");}
	F77_NAME(dsymv)(lower, &m, &one, &C0[mm*threadID], &m, &c0[m*threadID], &inc, &zero, &w[m*threadID], &inc FCONE);
	
	//make hat(y)
	for(j = 0; j < m; j++){
	  tmp_m[m*threadID+j] = y[nnIndx0[j*n0+i]] - F77_NAME(ddot)(&p, &X[nnIndx0[j*n0+i]], &n, &beta[k*p], &inc);
	}
	
	y0Hat[k*n0+i] = F77_NAME(ddot)(&p, &X0[i], &n0, &beta[k*p], &inc) + F77_NAME(ddot)(&m, &w[m*threadID], &inc, &tmp_m[m*threadID], &inc);
	
	//make u
	for(j = 0; j < m; j++){
	  F77_NAME(dcopy)(&p, &X[nnIndx0[j*n0+i]], &n, &tmp_mp[mp*threadID+j], &m);
	}
	
	F77_NAME(dgemv)(ytran, &m, &p, &one, &tmp_mp[mp*threadID], &m, &w[m*threadID], &inc, &zero, &tmp_p[p*threadID], &inc FCONE);
	
	for(j = 0; j < p; j++){
	  tmp_p[p*threadID+j] = X0[j*n0+i] - tmp_p[p*threadID+j];
	}
	
	//make v_y and var(y)
	F77_NAME(dsymv)(lower, &p, &one, tmp_pp, &p, &tmp_p[p*threadID], &inc, &zero, &tmp_p2[p*threadID], &inc FCONE);
	
	y0HatVar[k*n0+i] = ab[k*2+1] * (F77_NAME(ddot)(&p, &tmp_p[p*threadID], &inc, &tmp_p2[p*threadID], &inc) + 1.0 + alpha - F77_NAME(ddot)(&m, &w[m*threadID], &inc, &c0[m*threadID], &inc))/(ab[k*2]-1.0);

      }
   
      //report
      if(verbose){
	if(corName == "matern"){
	  Rprintf("Set phi=%.5f, nu=%.5f, and alpha=%.5f\n", thetaAlpha[k*nTheta], thetaAlpha[k*nTheta+2], thetaAlpha[k*nTheta+1]);
	}else{
	  Rprintf("Set phi=%.5f and alpha=%.5f\n", thetaAlpha[k*nTheta], thetaAlpha[k*nTheta+1]);
	}
	
        #ifdef Win32
	R_FlushConsole();
        #endif
      }
         
      R_CheckUserInterrupt();
    }
    
    //make return object
    SEXP result_r, resultName_r;
    int nResultListObjs = 4;
    
    if(n0 > 0){
      nResultListObjs += 2;
    }

    i = 0;
    
    PROTECT(result_r = allocVector(VECSXP, nResultListObjs)); nProtect++;
    PROTECT(resultName_r = allocVector(VECSXP, nResultListObjs)); nProtect++;
    
    SET_VECTOR_ELT(result_r, 0, beta_r);
    SET_VECTOR_ELT(resultName_r, 0, mkChar("beta.hat")); 

    SET_VECTOR_ELT(result_r, 1, ab_r);
    SET_VECTOR_ELT(resultName_r, 1, mkChar("ab")); 
    
    SET_VECTOR_ELT(result_r, 2, bBInv_r);
    SET_VECTOR_ELT(resultName_r, 2, mkChar("bB.inv"));
    
    SET_VECTOR_ELT(result_r, 3, bb_r);
    SET_VECTOR_ELT(resultName_r, 3, mkChar("bb"));
        
    if(n0 > 0){
      SET_VECTOR_ELT(result_r, 4, y0Hat_r);
      SET_VECTOR_ELT(resultName_r, 4, mkChar("y.0.hat"));
      
      SET_VECTOR_ELT(result_r, 5, y0HatVar_r);
      SET_VECTOR_ELT(resultName_r, 5, mkChar("y.0.hat.var")); 
    }
    
    namesgets(result_r, resultName_r);
    
    //unprotect
    UNPROTECT(nProtect);
    
    return(result_r);
  }
}
