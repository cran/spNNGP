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
void updateBF(double *B, double *F, double *c, double *C, double *coords, int *nnIndx, int *nnIndxLU, int n, int m, double sigmaSq, double phi, double nu, int covModel, double *bk, double nuUnifb){
    
  int i, k, l;
  int info = 0;
  int inc = 1;
  double one = 1.0;
  double zero = 0.0;
  char lower = 'L';

  //bk must be 1+(int)floor(alpha) * nthread
  int nb = 1+static_cast<int>(floor(nuUnifb));
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
	  c[m*threadID+k] = sigmaSq*spCor(e, phi, nu, covModel, &bk[threadID*nb]);
	  for(l = 0; l <= k; l++){
	    e = dist2(coords[nnIndx[nnIndxLU[i]+k]], coords[n+nnIndx[nnIndxLU[i]+k]], coords[nnIndx[nnIndxLU[i]+l]], coords[n+nnIndx[nnIndxLU[i]+l]]); 
	    C[mm*threadID+l*nnIndxLU[n+i]+k] = sigmaSq*spCor(e, phi, nu, covModel, &bk[threadID*nb]); 
	  }
	}
	F77_NAME(dpotrf)(&lower, &nnIndxLU[n+i], &C[mm*threadID], &nnIndxLU[n+i], &info FCONE); if(info != 0){error("c++ error: dpotrf failed\n");}
	F77_NAME(dpotri)(&lower, &nnIndxLU[n+i], &C[mm*threadID], &nnIndxLU[n+i], &info FCONE); if(info != 0){error("c++ error: dpotri failed\n");}
	F77_NAME(dsymv)(&lower, &nnIndxLU[n+i], &one, &C[mm*threadID], &nnIndxLU[n+i], &c[m*threadID], &inc, &zero, &B[nnIndxLU[i]], &inc FCONE);
	F[i] = sigmaSq - F77_NAME(ddot)(&nnIndxLU[n+i], &B[nnIndxLU[i]], &inc, &c[m*threadID], &inc);
      }else{
	B[i] = 0;
	F[i] = sigmaSq;
      }
    }

}

extern "C" {
  
  SEXP sNNGP(SEXP y_r, SEXP X_r, SEXP p_r, SEXP n_r, SEXP m_r, SEXP coords_r, SEXP covModel_r, SEXP nnIndx_r, SEXP nnIndxLU_r, SEXP uIndx_r, SEXP uIndxLU_r, SEXP uiIndx_r,
	     SEXP sigmaSqIG_r, SEXP tauSqIG_r, SEXP phiUnif_r, SEXP nuUnif_r, 
	     SEXP betaStarting_r, SEXP sigmaSqStarting_r, SEXP tauSqStarting_r, SEXP phiStarting_r, SEXP nuStarting_r,
	     SEXP sigmaSqTuning_r, SEXP tauSqTuning_r, SEXP phiTuning_r, SEXP nuTuning_r, 
	     SEXP nSamples_r, SEXP nThreads_r, SEXP verbose_r, SEXP nReport_r){
    
    int h, i, j, k, l, s, info, nProtect=0;
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
    int m = INTEGER(m_r)[0];
    double *coords = REAL(coords_r);
    int *nnIndx = INTEGER(nnIndx_r);
    int *nnIndxLU = INTEGER(nnIndxLU_r);
    int *uIndx = INTEGER(uIndx_r);
    int *uIndxLU = INTEGER(uIndxLU_r);
    int *uiIndx = INTEGER(uiIndx_r);
    int covModel = INTEGER(covModel_r)[0];
    std::string corName = getCorName(covModel);
        
    //priors
    double sigmaSqIGa = REAL(sigmaSqIG_r)[0]; double sigmaSqIGb = REAL(sigmaSqIG_r)[1];
    double tauSqIGa = REAL(tauSqIG_r)[0]; double tauSqIGb = REAL(tauSqIG_r)[1]; 
    double phiUnifa = REAL(phiUnif_r)[0]; double phiUnifb = REAL(phiUnif_r)[1];
    
    double nuUnifa = 0, nuUnifb = 0;
    if(corName == "matern"){
      nuUnifa = REAL(nuUnif_r)[0]; nuUnifb = REAL(nuUnif_r)[1]; 
    }
    
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
      Rprintf("\tModel description\n");
      Rprintf("----------------------------------------\n");
      Rprintf("NNGP Latent model fit with %i observations.\n\n", n);
      Rprintf("Number of covariates %i (including intercept if specified).\n\n", p);
      Rprintf("Using the %s spatial correlation model.\n\n", corName.c_str());
      Rprintf("Using %i nearest neighbors.\n\n", m);
      Rprintf("Number of MCMC samples %i.\n\n", nSamples);
      Rprintf("Priors and hyperpriors:\n");
      Rprintf("\tbeta flat.\n");
      Rprintf("\tsigma.sq IG hyperpriors shape=%.5f and scale=%.5f\n", sigmaSqIGa, sigmaSqIGb);
      Rprintf("\ttau.sq IG hyperpriors shape=%.5f and scale=%.5f\n", tauSqIGa, tauSqIGb); 
      Rprintf("\tphi Unif hyperpriors a=%.5f and b=%.5f\n", phiUnifa, phiUnifb);
      if(corName == "matern"){
	Rprintf("\tnu Unif hyperpriors a=%.5f and b=%.5f\n", nuUnifa, nuUnifb);	  
      }
#ifdef _OPENMP
      Rprintf("\nSource compiled with OpenMP support and model fit using %i thread(s).\n", nThreads);
#else
      Rprintf("\n\nSource not compiled with OpenMP support.\n");
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
    
    //starting	
    double *beta = (double *) R_alloc(p, sizeof(double));
    double *theta = (double *) R_alloc(nTheta, sizeof(double));
    
    F77_NAME(dcopy)(&p, REAL(betaStarting_r), &inc, beta, &inc);
    
    theta[sigmaSqIndx] = REAL(sigmaSqStarting_r)[0];
    theta[tauSqIndx] = REAL(tauSqStarting_r)[0];
    theta[phiIndx] = REAL(phiStarting_r)[0];
    if(corName == "matern"){
      theta[nuIndx] = REAL(nuStarting_r)[0];
    }
    
    //tuning and fixed
    double *tuning = (double *) R_alloc(nTheta, sizeof(double));
    
    tuning[sigmaSqIndx] = 0; //not accessed
    tuning[tauSqIndx] = 0; //not accessed  
    tuning[phiIndx] = REAL(phiTuning_r)[0];
    
    if(corName == "matern"){
      tuning[nuIndx] = REAL(nuTuning_r)[0];
    }

    //allocate for the U index vector that keep track of which locations have the i-th location as a neighbor
    int nIndx = static_cast<int>(static_cast<double>(1+m)/2*m+(n-m-1)*m);
    // int *uIndx = (int *) R_alloc(nIndx, sizeof(int)); //U indexes 
    
    // //first column holds the uIndx index for i-th location and second column holds
    // //the number of neighbors for which the i-th location is a neighbor
    // int *uIndxLU = (int *) R_alloc(2*n, sizeof(int));
    
    // //make u index
    // if(verbose){
    //   Rprintf("Building neighbors of neighbor index\n");
    //   #ifdef Win32
    //     R_FlushConsole();
    //   #endif
    // }
    // mkUIndx(n, m, nnIndx, uIndx, uIndxLU);
    
    // //u lists those locations that have the i-th location as a neighbor
    // //then for each of those locations that have i as a neighbor, we need to know the index of i in each of their B vectors (i.e. where does i fall in their neighbor set)
    // int *uiIndx = (int *) R_alloc(nIndx, sizeof(int));
    
    // for(i = 0; i < n; i++){//for each i
    //   for(j = 0; j < uIndxLU[n+i]; j++){//for each location that has i as a neighbor
    // 	k = uIndx[uIndxLU[i]+j];//index of a location that has i as a neighbor
    // 	uiIndx[uIndxLU[i]+j] = which(i, &nnIndx[nnIndxLU[k]], nnIndxLU[n+k]);
    //   }
    // }
 
    //other stuff
    int mm = m*m;
    double *B = (double *) R_alloc(nIndx, sizeof(double));
    double *F = (double *) R_alloc(n, sizeof(double));
    double *BCand = (double *) R_alloc(nIndx, sizeof(double));
    double *FCand = (double *) R_alloc(n, sizeof(double));
    double *c =(double *) R_alloc(m*nThreads, sizeof(double));
    double *C = (double *) R_alloc(mm*nThreads, sizeof(double));

    //return stuff  
    SEXP betaSamples_r, thetaSamples_r, wSamples_r;
    PROTECT(betaSamples_r = allocMatrix(REALSXP, p, nSamples)); nProtect++;
    PROTECT(thetaSamples_r = allocMatrix(REALSXP, nTheta, nSamples)); nProtect++; 
    PROTECT(wSamples_r = allocMatrix(REALSXP, n, nSamples)); nProtect++;
    
    //other stuff
    double logPostCand, logPostCurrent, logDet;
    int accept = 0, batchAccept = 0, status = 0;
    int jj, kk, pp = p*p;
    double *tmp_pp = (double *) R_alloc(pp, sizeof(double));
    double *tmp_p = (double *) R_alloc(p, sizeof(double));
    double *tmp_p2 = (double *) R_alloc(p, sizeof(double));
    double *tmp_n = (double *) R_alloc(n, sizeof(double)); zeros(tmp_n, n);
    double *XtX = (double *) R_alloc(pp, sizeof(double));
    double *w = (double *) R_alloc(n, sizeof(double)); zeros(w, n);
    double a, v, b, e, mu, var, aij, phiCand, nuCand = 0, nu = 0;

    double *bk = (double *) R_alloc(nThreads*(1.0+static_cast<int>(floor(nuUnifb))), sizeof(double));
    
    F77_NAME(dgemm)(ytran, ntran, &p, &p, &n, &one, X, &n, X, &n, &zero, XtX, &p FCONE FCONE);

    if(verbose){
      Rprintf("----------------------------------------\n");
      Rprintf("\t\tSampling\n");
      Rprintf("----------------------------------------\n");
      #ifdef Win32
        R_FlushConsole();
      #endif
    }

    if(corName == "matern"){nu = theta[nuIndx];}
    updateBF(B, F, c, C, coords, nnIndx, nnIndxLU, n, m, theta[sigmaSqIndx], theta[phiIndx], nu, covModel, bk, nuUnifb);
      
    GetRNGstate();
    
    for(s = 0; s < nSamples; s++){

      ///////////////
      //update w 
      ///////////////
      for(i = 0; i < n; i++){
	a = 0;
	v = 0;
	if(uIndxLU[n+i] > 0){//is i a neighbor for anybody
	  for(j = 0; j < uIndxLU[n+i]; j++){//how many location have i as a neighbor
	    b = 0;
	    //now the neighbors for the jth location who has i as a neighbor
	    jj = uIndx[uIndxLU[i]+j]; //jj is the index of the jth location who has i as a neighbor
	    for(k = 0; k < nnIndxLU[n+jj]; k++){// these are the neighbors of the jjth location
	      kk = nnIndx[nnIndxLU[jj]+k];// kk is the index for the jth locations neighbors
	      if(kk != i){//if the neighbor of jj is not i
		b += B[nnIndxLU[jj]+k]*w[kk];//covariance between jj and kk and the random effect of kk
	      }
	    }
	    aij = w[jj] - b;
	    a += B[nnIndxLU[jj]+uiIndx[uIndxLU[i]+j]]*aij/F[jj];
	    v += pow(B[nnIndxLU[jj]+uiIndx[uIndxLU[i]+j]],2)/F[jj];
	  }
	}
	
	e = 0;
	for(j = 0; j < nnIndxLU[n+i]; j++){
	  e += B[nnIndxLU[i]+j]*w[nnIndx[nnIndxLU[i]+j]];
	}
	mu = (y[i] - F77_NAME(ddot)(&p, &X[i], &n, beta, &inc))/theta[tauSqIndx] + e/F[i] + a;
	
	var = 1.0/(1.0/theta[tauSqIndx] + 1.0/F[i] + v);
	w[i] = rnorm(mu*var, sqrt(var));
      }
      
      ///////////////
      //update beta 
      ///////////////
      for(i = 0; i < n; i++){
	tmp_n[i] = (y[i] - w[i])/theta[tauSqIndx];
      }
      F77_NAME(dgemv)(ytran, &n, &p, &one, X, &n, tmp_n, &inc, &zero, tmp_p, &inc FCONE); 	  
      
      for(i = 0; i < pp; i++){
	tmp_pp[i] = XtX[i]/theta[tauSqIndx];
      }
      
      F77_NAME(dpotrf)(lower, &p, tmp_pp, &p, &info FCONE); if(info != 0){error("c++ error: dpotrf failed\n");}
      F77_NAME(dpotri)(lower, &p, tmp_pp, &p, &info FCONE); if(info != 0){error("c++ error: dpotri failed\n");}
      F77_NAME(dsymv)(lower, &p, &one, tmp_pp, &p, tmp_p, &inc, &zero, tmp_p2, &inc FCONE);
      F77_NAME(dpotrf)(lower, &p, tmp_pp, &p, &info FCONE); if(info != 0){error("c++ error: dpotrf failed\n");}
      mvrnorm(beta, tmp_p2, tmp_pp, p);
      
      /////////////////////
      //update tau^2
      /////////////////////
      for(i = 0; i < n; i++){
	tmp_n[i] = y[i] - w[i] - F77_NAME(ddot)(&p, &X[i], &n, beta, &inc);
      }
      
      theta[tauSqIndx] = 1.0/rgamma(tauSqIGa+n/2.0, 1.0/(tauSqIGb+0.5*F77_NAME(ddot)(&n, tmp_n, &inc, tmp_n, &inc)));

      /////////////////////
      //update sigma^2
      /////////////////////
      a = 0;

#ifdef _OPENMP
#pragma omp parallel for private (e, j, b) reduction(+:a, logDet)
#endif
      for(i = 0; i < n; i++){
	if(nnIndxLU[n+i] > 0){
	  e = 0;
	  for(j = 0; j < nnIndxLU[n+i]; j++){
	    e += B[nnIndxLU[i]+j]*w[nnIndx[nnIndxLU[i]+j]];
	  }
	  b = w[i] - e;
	}else{
	  b = w[i];
	}	
	a += b*b/F[i];
      }

      theta[sigmaSqIndx] = 1.0/rgamma(sigmaSqIGa+n/2.0, 1.0/(sigmaSqIGb+0.5*a*theta[sigmaSqIndx]));
    
      ///////////////
      //update theta
      ///////////////
      //current
      if(corName == "matern"){nu = theta[nuIndx];}
      updateBF(B, F, c, C, coords, nnIndx, nnIndxLU, n, m, theta[sigmaSqIndx], theta[phiIndx], nu, covModel, bk, nuUnifb);
      
      a = 0;
      logDet = 0;

#ifdef _OPENMP
#pragma omp parallel for private (e, j, b) reduction(+:a, logDet)
#endif
      for(i = 0; i < n; i++){
	if(nnIndxLU[n+i] > 0){
	  e = 0;
	  for(j = 0; j < nnIndxLU[n+i]; j++){
	    e += B[nnIndxLU[i]+j]*w[nnIndx[nnIndxLU[i]+j]];
	  }
	  b = w[i] - e;
	}else{
	  b = w[i];
	}	
	a += b*b/F[i];
	logDet += log(F[i]);
      }
      
      logPostCurrent = -0.5*logDet - 0.5*a;
      logPostCurrent += log(theta[phiIndx] - phiUnifa) + log(phiUnifb - theta[phiIndx]); 
      if(corName == "matern"){
      	logPostCurrent += log(theta[nuIndx] - nuUnifa) + log(nuUnifb - theta[nuIndx]); 
      }
      
      //candidate
      phiCand = logitInv(rnorm(logit(theta[phiIndx], phiUnifa, phiUnifb), tuning[phiIndx]), phiUnifa, phiUnifb);

      if(corName == "matern"){
      	nuCand = logitInv(rnorm(logit(theta[nuIndx], nuUnifa, nuUnifb), tuning[nuIndx]), nuUnifa, nuUnifb);
      }
      
      updateBF(BCand, FCand, c, C, coords, nnIndx, nnIndxLU, n, m, theta[sigmaSqIndx], phiCand, nuCand, covModel, bk, nuUnifb);
            
      a = 0;
      logDet = 0;
      
#ifdef _OPENMP
#pragma omp parallel for private (e, j, b) reduction(+:a, logDet)
#endif
      for(i = 0; i < n; i++){
	if(nnIndxLU[n+i] > 0){
	  e = 0;
	  for(j = 0; j < nnIndxLU[n+i]; j++){
	    e += BCand[nnIndxLU[i]+j]*w[nnIndx[nnIndxLU[i]+j]];
	  }
	  b = w[i] - e;
	}else{
	  b = w[i];
	  }	
	  a += b*b/FCand[i];
	  logDet += log(FCand[i]);
	}
      
      logPostCand = -0.5*logDet - 0.5*a;      
      logPostCand += log(phiCand - phiUnifa) + log(phiUnifb - phiCand); 
      if(corName == "matern"){
      	logPostCand += log(nuCand - nuUnifa) + log(nuUnifb - nuCand); 
      }

      if(runif(0.0,1.0) <= exp(logPostCand - logPostCurrent)){

	std::swap(BCand, B);
	std::swap(FCand, F);
	
	theta[phiIndx] = phiCand;
	if(corName == "matern"){
	  theta[nuIndx] = nuCand; 
	}
	
	accept++;
	batchAccept++;
      }

      //save samples  
      F77_NAME(dcopy)(&p, beta, &inc, &REAL(betaSamples_r)[s*p], &inc);
      F77_NAME(dcopy)(&nTheta, theta, &inc, &REAL(thetaSamples_r)[s*nTheta], &inc);
      F77_NAME(dcopy)(&n, w, &inc, &REAL(wSamples_r)[s*n], &inc);
      
      //report
      if(status == nReport){
	if(verbose){
	  Rprintf("Sampled: %i of %i, %3.2f%%\n", s, nSamples, 100.0*s/nSamples);
	  Rprintf("Report interval Metrop. Acceptance rate: %3.2f%%\n", 100.0*batchAccept/nReport);
	  Rprintf("Overall Metrop. Acceptance rate: %3.2f%%\n", 100.0*accept/s);
      	  Rprintf("-------------------------------------------------\n");
          #ifdef Win32
      	  R_FlushConsole();
          #endif
      	}
	batchAccept = 0;
	status = 0;
      }
      
      status++;
      
      R_CheckUserInterrupt();
    }

    if(verbose){
      Rprintf("Sampled: %i of %i, %3.2f%%\n", s, nSamples, 100.0);
      Rprintf("Report interval Metrop. Acceptance rate: %3.2f%%\n", 100.0*batchAccept/nReport);
      Rprintf("Overall Metrop. Acceptance rate: %3.2f%%\n", 100.0*accept/nSamples);
      Rprintf("-------------------------------------------------\n");
      #ifdef Win32
      R_FlushConsole();
      #endif
    }
   
    PutRNGstate();

    //make return object
    SEXP result_r, resultName_r;
    int nResultListObjs = 3;
    
    PROTECT(result_r = allocVector(VECSXP, nResultListObjs)); nProtect++;
    PROTECT(resultName_r = allocVector(VECSXP, nResultListObjs)); nProtect++;

    SET_VECTOR_ELT(result_r, 0, betaSamples_r);
    SET_VECTOR_ELT(resultName_r, 0, mkChar("p.beta.samples")); 
    
    SET_VECTOR_ELT(result_r, 1, thetaSamples_r);
    SET_VECTOR_ELT(resultName_r, 1, mkChar("p.theta.samples"));

    SET_VECTOR_ELT(result_r, 2, wSamples_r);
    SET_VECTOR_ELT(resultName_r, 2, mkChar("p.w.samples"));
	
    namesgets(result_r, resultName_r);
    
    //unprotect
    UNPROTECT(nProtect);
    
    return(result_r);
  }
}
