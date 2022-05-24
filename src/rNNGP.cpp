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

//Description: update replicated data.
void updateRep(double *B, double *F, int n, double *tmp_m, double *tmp_n, int *nnIndx, int *nnIndxLU){

  int inc = 1;
  int i, j;
  double z;
  
  for(i = 0; i < n; i++){
    z = rnorm(0.0, 1.0);
    if(i == 0){
      tmp_n[i] = sqrt(F[i])*z;
    }else{
      for(j = 0; j < nnIndxLU[n+i]; j++){
	tmp_m[j] = tmp_n[nnIndx[nnIndxLU[i]+j]];
      }
      tmp_n[i] = F77_NAME(ddot)(&nnIndxLU[n+i], &B[nnIndxLU[i]], &inc, tmp_m, &inc)+sqrt(F[i])*z;
    }
  }
  
}


//Description: update B and F.
double updateBF(double *B, double *F, double *c, double *C, double *coords, int *nnIndx, int *nnIndxLU, int n, int m, double *theta, int tauSqIndx, int sigmaSqIndx, int phiIndx, int nuIndx, int covModel, double *bk, double nuUnifb){
    
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
    
    for(i = 0; i < n; i++){
      logDet += log(F[i]);
    }

    return(logDet);
}

extern "C" {
  
  SEXP rNNGP(SEXP y_r, SEXP X_r, SEXP p_r, SEXP n_r, SEXP m_r, SEXP coords_r, SEXP covModel_r, SEXP nnIndx_r, SEXP nnIndxLU_r,
	     SEXP sigmaSqIG_r, SEXP tauSqIG_r, SEXP phiUnif_r, SEXP nuUnif_r, 
	     SEXP betaStarting_r, SEXP sigmaSqStarting_r, SEXP tauSqStarting_r, SEXP phiStarting_r, SEXP nuStarting_r,
	     SEXP sigmaSqTuning_r, SEXP tauSqTuning_r, SEXP phiTuning_r, SEXP nuTuning_r, 
	     SEXP nSamples_r, SEXP nThreads_r, SEXP verbose_r, SEXP nReport_r, SEXP nRep_r, SEXP repIndx_r){
    
    int h, i, j, k, l, s, info, nProtect=0;
    const int inc = 1;
    const double one = 1.0;
    const double negOne = -1.0;
    const double zero = 0.0;
    char const *lower = "L";
    char const *ntran = "N";
    
    //get args
    double *y = REAL(y_r);
    double *X = REAL(X_r);
    int p = INTEGER(p_r)[0];
    int n = INTEGER(n_r)[0];
    int m = INTEGER(m_r)[0];
    double *coords = REAL(coords_r);
    int *nnIndx = INTEGER(nnIndx_r);
    int *nnIndxLU = INTEGER(nnIndxLU_r);
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
    int nRep = INTEGER(nRep_r)[0];
    int *repIndx = INTEGER(repIndx_r);
    
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

    //other stuff
    int nIndx = static_cast<int>(static_cast<double>(1+m)/2*m+(n-m-1)*m);
    int mm = m*m;
    double *thetaCand = (double *) R_alloc(nTheta, sizeof(double));
    double *B = (double *) R_alloc(nIndx, sizeof(double));
    double *F = (double *) R_alloc(n, sizeof(double));
    double *c =(double *) R_alloc(m*nThreads, sizeof(double));
    double *C = (double *) R_alloc(mm*nThreads, sizeof(double));

    //return stuff  
    SEXP betaSamples_r, thetaSamples_r, repSamples_r;
    PROTECT(betaSamples_r = allocMatrix(REALSXP, p, nSamples)); nProtect++;
    PROTECT(thetaSamples_r = allocMatrix(REALSXP, nTheta, nSamples)); nProtect++; 

    if(nRep){
      PROTECT(repSamples_r = allocMatrix(REALSXP, n, nRep)); nProtect++;
      zeros(REAL(repSamples_r), n*nRep);
    }
        
    //other stuff
    double logPostCand, logPostCurrent, logDetCurrent, logDetCand, QCurrent, QCand;
    int accept = 0, batchAccept = 0, status = 0, repCnt = 0;
    int pp = p*p;
    double *tmp_pp = (double *) R_alloc(pp, sizeof(double));
    double *tmp_p = (double *) R_alloc(p, sizeof(double));
    double *tmp_p2 = (double *) R_alloc(p, sizeof(double));
    double *tmp_n = (double *) R_alloc(n, sizeof(double));
    double *tmp_n2 = NULL;
    double *tmp_m = NULL;
    if(nRep){
      tmp_n2 = (double *) R_alloc(n, sizeof(double));
      tmp_m = (double *) R_alloc(m, sizeof(double));
    }
    double *bk = (double *) R_alloc(nThreads*(1.0+static_cast<int>(floor(nuUnifb))), sizeof(double));
    
    bool thetaUpdate = true;
    
    //update B and F
    logDetCurrent = updateBF(B, F, c, C, coords, nnIndx, nnIndxLU, n, m, theta, tauSqIndx, sigmaSqIndx, phiIndx, nuIndx, covModel, bk, nuUnifb);

    F77_NAME(dgemv)(ntran, &n, &p, &one, X, &n, beta, &inc, &zero, tmp_n, &inc FCONE);
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
	
	F77_NAME(dpotrf)(lower, &p, tmp_pp, &p, &info FCONE); if(info != 0){error("c++ error: dpotrf failed\n");}
	F77_NAME(dpotri)(lower, &p, tmp_pp, &p, &info FCONE); if(info != 0){error("c++ error: dpotri failed\n");}
	F77_NAME(dsymv)(lower, &p, &one, tmp_pp, &p, tmp_p, &inc, &zero, tmp_p2, &inc FCONE);
	F77_NAME(dpotrf)(lower, &p, tmp_pp, &p, &info FCONE); if(info != 0){error("c++ error: dpotrf failed\n");}
      }
      
      mvrnorm(beta, tmp_p2, tmp_pp, p);
      
      ///////////////
      //update theta
      ///////////////
      F77_NAME(dgemv)(ntran, &n, &p, &one, X, &n, beta, &inc, &zero, tmp_n, &inc FCONE);
      if(nRep && repIndx[s]){
	F77_NAME(dcopy)(&n, tmp_n, &inc, &REAL(repSamples_r)[repCnt*n], &inc);
      }
      F77_NAME(daxpy)(&n, &negOne, y, &inc, tmp_n, &inc);
      
      //current    
      logDetCurrent = updateBF(B, F, c, C, coords, nnIndx, nnIndxLU, n, m, theta, tauSqIndx, sigmaSqIndx, phiIndx, nuIndx, covModel, bk, nuUnifb);

      //update rep
      if(nRep && repIndx[s]){
	updateRep(B, F, n, tmp_m, tmp_n2, nnIndx, nnIndxLU);//tmp_n is tilde{z}
	F77_NAME(daxpy)(&n, &one, tmp_n2, &inc, &REAL(repSamples_r)[repCnt*n], &inc);
	repCnt++;
      }
      
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
      logDetCand = updateBF(B, F, c, C, coords, nnIndx, nnIndxLU, n, m, thetaCand, tauSqIndx, sigmaSqIndx, phiIndx, nuIndx, covModel, bk, nuUnifb);
      
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
    int nResultListObjs = 2;

    if(nRep){
      nResultListObjs++;
    }
    
    PROTECT(result_r = allocVector(VECSXP, nResultListObjs)); nProtect++;
    PROTECT(resultName_r = allocVector(VECSXP, nResultListObjs)); nProtect++;

    SET_VECTOR_ELT(result_r, 0, betaSamples_r);
    SET_VECTOR_ELT(resultName_r, 0, mkChar("p.beta.samples")); 
    
    SET_VECTOR_ELT(result_r, 1, thetaSamples_r);
    SET_VECTOR_ELT(resultName_r, 1, mkChar("p.theta.samples"));

    if(nRep){
      SET_VECTOR_ELT(result_r, 2, repSamples_r);
      SET_VECTOR_ELT(resultName_r, 2, mkChar("y.rep.samples"));
    }
    
    namesgets(result_r, resultName_r);
    
    //unprotect
    UNPROTECT(nProtect);
    
    return(result_r);
  }
}
