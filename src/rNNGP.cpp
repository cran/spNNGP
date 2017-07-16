#include <string>
#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>
#include <R_ext/Linpack.h>
#include <R_ext/Lapack.h>
#include <R_ext/BLAS.h>
#include "util.h"

#ifdef _OPENMP
#include <omp.h>
#endif

//Description: update B and F.
double updateBF(double *B, double *F, double *c, double *C, double *D, double *d, int *nnIndxLU, int *CIndx, int n, double *theta, int tauSqIndx, int sigmaSqIndx, int phiIndx, int nuIndx, int covModel, double *bk, double nuUnifb){
    
  int i, k, l;
  int info = 0;
  int inc = 1;
  double one = 1.0;
  double zero = 0.0;
  char lower = 'L';
  double logDet = 0;
  double nu = 0;

  if(getCorName(covModel) == "matern"){
    nu = theta[nuIndx];
  }

  //bk must be 1+(int)floor(alpha) * nthread
  int nb = 1+static_cast<int>(floor(nuUnifb));
  int threadID = 0;
  
#ifdef _OPENMP
#pragma omp parallel for private(k, l, info, threadID)
#endif
    for(i = 0; i < n; i++){
#ifdef _OPENMP
      threadID = omp_get_thread_num();
#endif
      if(i > 0){
	for(k = 0; k < nnIndxLU[n+i]; k++){
	  c[nnIndxLU[i]+k] = theta[sigmaSqIndx]*spCor(d[nnIndxLU[i]+k], theta[phiIndx], nu, covModel, &bk[threadID*nb]);
	  for(l = 0; l <= k; l++){
	    C[CIndx[i]+l*nnIndxLU[n+i]+k] = theta[sigmaSqIndx]*spCor(D[CIndx[i]+l*nnIndxLU[n+i]+k], theta[phiIndx], nu, covModel, &bk[threadID*nb]); 
	    if(l == k){
	      C[CIndx[i]+l*nnIndxLU[n+i]+k] += theta[tauSqIndx];
	    }
	  }
	}
	F77_NAME(dpotrf)(&lower, &nnIndxLU[n+i], &C[CIndx[i]], &nnIndxLU[n+i], &info); if(info != 0){error("c++ error: dpotrf failed\n");}
	F77_NAME(dpotri)(&lower, &nnIndxLU[n+i], &C[CIndx[i]], &nnIndxLU[n+i], &info); if(info != 0){error("c++ error: dpotri failed\n");}
	F77_NAME(dsymv)(&lower, &nnIndxLU[n+i], &one, &C[CIndx[i]], &nnIndxLU[n+i], &c[nnIndxLU[i]], &inc, &zero, &B[nnIndxLU[i]], &inc);
	F[i] = theta[sigmaSqIndx] - F77_NAME(ddot)(&nnIndxLU[n+i], &B[nnIndxLU[i]], &inc, &c[nnIndxLU[i]], &inc) + theta[tauSqIndx];
      }else{
	B[i] = 0;
	F[i] = theta[sigmaSqIndx] + theta[tauSqIndx];
      }
    }
    
    for(i = 0; i < n; i++){
      logDet += log(F[i]);
    }

    return(logDet);
}

extern "C" {
  
  SEXP rNNGP(SEXP y_r, SEXP X_r, SEXP p_r, SEXP n_r, SEXP m_r, SEXP coords_r, SEXP covModel_r, 
	     SEXP sigmaSqIG_r, SEXP tauSqIG_r, SEXP phiUnif_r, SEXP nuUnif_r, 
	     SEXP betaStarting_r, SEXP sigmaSqStarting_r, SEXP tauSqStarting_r, SEXP phiStarting_r, SEXP nuStarting_r,
	     SEXP sigmaSqTuning_r, SEXP tauSqTuning_r, SEXP phiTuning_r, SEXP nuTuning_r, 
	     SEXP nSamples_r, SEXP sType_r, SEXP returnNNIndx_r, SEXP nThreads_r, SEXP verbose_r, SEXP nReport_r){
    
    int h, i, j, k, l, s, info, nProtect=0;
    const int inc = 1;
    const double one = 1.0;
    const double negOne = -1.0;
    const double zero = 0.0;
    char const *lower = "L";
    char const *upper = "U";
    char const *ntran = "N";
    char const *ytran = "T";
    char const *rside = "R";
    char const *lside = "L";
    
    //get args
    double *y = REAL(y_r);
    double *X = REAL(X_r);
    int p = INTEGER(p_r)[0];
    int n = INTEGER(n_r)[0];
    int m = INTEGER(m_r)[0];
    double *coords = REAL(coords_r);
    
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
      Rprintf("NNGP Response model fit with %i observations.\n\n", n);
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
    
    tuning[sigmaSqIndx] = REAL(sigmaSqTuning_r)[0];
    tuning[tauSqIndx] = REAL(tauSqTuning_r)[0];  
    tuning[phiIndx] = REAL(phiTuning_r)[0];
    
    if(corName == "matern"){
      tuning[nuIndx] = REAL(nuTuning_r)[0];
    }

    //allocated for the nearest neighbor index vector (note, first location has no neighbors).
    int nIndx = static_cast<int>(static_cast<double>(1+m)/2*m+(n-m-1)*m);
    //int *nnIndx = (int *) R_alloc(nIndx, sizeof(int));  
    SEXP nnIndx_r; PROTECT(nnIndx_r = allocVector(INTSXP, nIndx)); nProtect++; int *nnIndx = INTEGER(nnIndx_r);
    double *d = (double *) R_alloc(nIndx, sizeof(double));
    int *nnIndxLU = (int *) R_alloc(2*n, sizeof(int)); //first column holds the nnIndx index for the i-th location and the second columns holds the number of neighbors the i-th location has (the second column is a bit of a waste but will simplifying some parallelization).

    //make the neighbor index
    if(verbose){
      Rprintf("----------------------------------------\n");
      Rprintf("\tBuilding neighbor index\n");
      #ifdef Win32
        R_FlushConsole();
      #endif
    }

    if(INTEGER(sType_r)[0] == 0){
      mkNNIndx(n, m, coords, nnIndx, d, nnIndxLU);
    }else{
      mkNNIndxTree0(n, m, coords, nnIndx, d, nnIndxLU);
    }

    //other stuff
    double *thetaCand = (double *) R_alloc(nTheta, sizeof(double));
    double *B = (double *) R_alloc(nIndx, sizeof(double));
    double *F = (double *) R_alloc(n, sizeof(double));
    double *c =(double *) R_alloc(nIndx, sizeof(double));
    
    int *CIndx = (int *) R_alloc(2*n, sizeof(int)); //index for D and C.
    for(i = 0, j = 0; i < n; i++){//zero should never be accessed
      j += nnIndxLU[n+i]*nnIndxLU[n+i];
      if(i == 0){
	CIndx[n+i] = 0;
	CIndx[i] = 0;
      }else{
	CIndx[n+i] = nnIndxLU[n+i]*nnIndxLU[n+i]; 
	CIndx[i] = CIndx[n+i-1] + CIndx[i-1];
      }
    }
    
    double *C = (double *) R_alloc(j, sizeof(double)); zeros(C, j);
    double *D = (double *) R_alloc(j, sizeof(double)); zeros(D, j);
    
    for(i = 0; i < n; i++){
      for(k = 0; k < nnIndxLU[n+i]; k++){   
	for(l = 0; l <= k; l++){
	  D[CIndx[i]+l*nnIndxLU[n+i]+k] = dist2(coords[nnIndx[nnIndxLU[i]+k]], coords[n+nnIndx[nnIndxLU[i]+k]], coords[nnIndx[nnIndxLU[i]+l]], coords[n+nnIndx[nnIndxLU[i]+l]]);
	}
      }
    }

    //return stuff  
    SEXP betaSamples_r, thetaSamples_r;
    PROTECT(betaSamples_r = allocMatrix(REALSXP, p, nSamples)); nProtect++;
    PROTECT(thetaSamples_r = allocMatrix(REALSXP, nTheta, nSamples)); nProtect++; 

    //other stuff
    double logPostCand, logPostCurrent, logDetCurrent, logDetCand, QCurrent, QCand, accept = 0, batchAccept = 0, status = 0;
    int pp = p*p;
    double *tmp_pp = (double *) R_alloc(pp, sizeof(double));
    double *tmp_p = (double *) R_alloc(p, sizeof(double));
    double *tmp_p2 = (double *) R_alloc(p, sizeof(double));
    double *tmp_n = (double *) R_alloc(n, sizeof(double));
    double *bk = (double *) R_alloc(nThreads*(static_cast<int>(1.0+nuUnifb)), sizeof(double));
    
    bool thetaUpdate = true;
    
    //update B and F
    logDetCurrent = updateBF(B, F, c, C, D, d, nnIndxLU, CIndx, n, theta, tauSqIndx, sigmaSqIndx, phiIndx, nuIndx, covModel, bk, nuUnifb);
    
    F77_NAME(dgemv)(ntran, &n, &p, &one, X, &n, beta, &inc, &zero, tmp_n, &inc);
    F77_NAME(daxpy)(&n, &negOne, y, &inc, tmp_n, &inc);
    QCurrent = Q(B, F, tmp_n, tmp_n, n, nnIndx, nnIndxLU);

    if(verbose){
      Rprintf("----------------------------------------\n");
      Rprintf("\t\tSampling\n");
      Rprintf("----------------------------------------\n");
      #ifdef Win32
        R_FlushConsole();
      #endif
    }

    GetRNGstate();
    
    for(s = 0; s < nSamples; s++){
      
      if(thetaUpdate){
	
	thetaUpdate = false;
	
	///////////////
	//update beta 
	///////////////
	for(i = 0; i < p; i++){
	  tmp_p[i] = Q(B, F, &X[n*i], y, n, nnIndx, nnIndxLU);
	  for(j = 0; j <= i; j++){
	    tmp_pp[j*p+i] = Q(B, F, &X[n*j], &X[n*i], n, nnIndx, nnIndxLU);
	  }
	}
	
	F77_NAME(dpotrf)(lower, &p, tmp_pp, &p, &info); if(info != 0){error("c++ error: dpotrf failed\n");}
	F77_NAME(dpotri)(lower, &p, tmp_pp, &p, &info); if(info != 0){error("c++ error: dpotri failed\n");}
	F77_NAME(dsymv)(lower, &p, &one, tmp_pp, &p, tmp_p, &inc, &zero, tmp_p2, &inc);
	F77_NAME(dpotrf)(lower, &p, tmp_pp, &p, &info); if(info != 0){error("c++ error: dpotrf failed\n");}
      }
      
      mvrnorm(beta, tmp_p2, tmp_pp, p);
      
      ///////////////
      //update theta
      ///////////////
      F77_NAME(dgemv)(ntran, &n, &p, &one, X, &n, beta, &inc, &zero, tmp_n, &inc);
      F77_NAME(daxpy)(&n, &negOne, y, &inc, tmp_n, &inc);
      
      //current    
      logDetCurrent = updateBF(B, F, c, C, D, d, nnIndxLU, CIndx, n, theta, tauSqIndx, sigmaSqIndx, phiIndx, nuIndx, covModel, bk, nuUnifb);
      
      QCurrent = Q(B, F, tmp_n, tmp_n, n, nnIndx, nnIndxLU);

      logPostCurrent = -0.5*logDetCurrent - 0.5*QCurrent;
      logPostCurrent += log(theta[phiIndx] - phiUnifa) + log(phiUnifb - theta[phiIndx]); 
      logPostCurrent += -1.0*(1.0+sigmaSqIGa)*log(theta[sigmaSqIndx])-sigmaSqIGb/theta[sigmaSqIndx]+log(theta[sigmaSqIndx]);
      logPostCurrent += -1.0*(1.0+tauSqIGa)*log(theta[tauSqIndx])-tauSqIGb/theta[tauSqIndx]+log(theta[tauSqIndx]);

       if(corName == "matern"){
	 logPostCurrent += log(theta[nuIndx] - nuUnifa) + log(nuUnifb - theta[nuIndx]); 
       }
      
      //candidate
      thetaCand[phiIndx] = logitInv(rnorm(logit(theta[phiIndx], phiUnifa, phiUnifb), tuning[phiIndx]), phiUnifa, phiUnifb);
      thetaCand[sigmaSqIndx] = exp(rnorm(log(theta[sigmaSqIndx]), tuning[sigmaSqIndx]));
      thetaCand[tauSqIndx] = exp(rnorm(log(theta[tauSqIndx]), tuning[tauSqIndx]));

      if(corName == "matern"){
	thetaCand[nuIndx] = logitInv(rnorm(logit(theta[nuIndx], nuUnifa, nuUnifb), tuning[nuIndx]), nuUnifa, nuUnifb);
      }
      
      //update B and F
      logDetCand = updateBF(B, F, c, C, D, d, nnIndxLU, CIndx, n, thetaCand, tauSqIndx, sigmaSqIndx, phiIndx, nuIndx, covModel, bk, nuUnifb);
      
      QCand = Q(B, F, tmp_n, tmp_n, n, nnIndx, nnIndxLU);
      
      logPostCand = -0.5*logDetCand - 0.5*QCand;
      logPostCand += log(thetaCand[phiIndx] - phiUnifa) + log(phiUnifb - thetaCand[phiIndx]); 
      logPostCand += -1.0*(1.0+sigmaSqIGa)*log(thetaCand[sigmaSqIndx])-sigmaSqIGb/thetaCand[sigmaSqIndx]+log(thetaCand[sigmaSqIndx]);
      logPostCand += -1.0*(1.0+tauSqIGa)*log(thetaCand[tauSqIndx])-tauSqIGb/thetaCand[tauSqIndx]+log(thetaCand[tauSqIndx]);

       if(corName == "matern"){
	 logPostCand += log(thetaCand[nuIndx] - nuUnifa) + log(nuUnifb - thetaCand[nuIndx]); 
       }
      
      if(runif(0.0,1.0) <= exp(logPostCand - logPostCurrent)){
	thetaUpdate = true;
	dcopy_(&nTheta, thetaCand, &inc, theta, &inc);
	accept++;
	batchAccept++;
      }

      //save samples
      F77_NAME(dcopy)(&p, beta, &inc, &REAL(betaSamples_r)[s*p], &inc);
      F77_NAME(dcopy)(&nTheta, theta, &inc, &REAL(thetaSamples_r)[s*nTheta], &inc);

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
    
    PutRNGstate();

    //make return object
    SEXP result_r, resultName_r;
    int nResultListObjs = 2;

    if(INTEGER(returnNNIndx_r)[0]){
      nResultListObjs++;
    }
    
    PROTECT(result_r = allocVector(VECSXP, nResultListObjs)); nProtect++;
    PROTECT(resultName_r = allocVector(VECSXP, nResultListObjs)); nProtect++;

    SET_VECTOR_ELT(result_r, 0, betaSamples_r);
    SET_VECTOR_ELT(resultName_r, 0, mkChar("p.beta.samples")); 
    
    SET_VECTOR_ELT(result_r, 1, thetaSamples_r);
    SET_VECTOR_ELT(resultName_r, 1, mkChar("p.theta.samples"));

    if(INTEGER(returnNNIndx_r)[0]){
      SET_VECTOR_ELT(result_r, 2, nnIndx_r);
      SET_VECTOR_ELT(resultName_r, 2, mkChar("n.indx")); 
    }
    
    namesgets(result_r, resultName_r);
    
    //unprotect
    UNPROTECT(nProtect);
    
    return(result_r);
  }
}
