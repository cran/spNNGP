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

void updateConjBF(double *B, double *F, double *R_iS, double *R_NiS, double *R_Ni, double *R_S, double *R_SInv, double *J_i, double *J_Ni, double *Omega_i, double *Omega_iNi, double *tmp_r, double *tmp_m, double *tmp_mr, double *tmp_mm,	  
		  double *coords, double *knots, int *nnIndx, int *nnIndxLU,  int n, int m, int r, double phi, double alpha, double nu, int covModel, double *bk, double nuMax){
    
  int i, k, l;
  int info = 0;
  int inc = 1;
  double one = 1.0;
  double zero = 0.0;
  char lower = 'L';
  char ntran = 'N';
  char ytran = 'T';
  char rside = 'R';
  double e, Omega_ii;
  
  //bk must be 1+(int)floor(alpha) * nthread
  int nb = 1+static_cast<int>(floor(nuMax));
  int threadID = 0;

  int rr = r*r;
  int mm = m*m;
  int mr = m*r;
  
  for(i = 0; i < r; i++){
    for(l = 0; l < r; l++){
      e = dist2(knots[i], knots[r+i], knots[l], knots[r+l]);
      R_S[i*r+l] = spCor(e, phi, nu, covModel, bk); //R(S*)
    }
  }
  
  F77_NAME(dcopy)(&rr, R_S, &inc, R_SInv, &inc);//R(S*)^{-1}
  F77_NAME(dpotrf)(&lower, &r, R_SInv, &r, &info FCONE); if(info != 0){error("c++ error: dpotrf failed 1a\n");}
  F77_NAME(dpotri)(&lower, &r, R_SInv, &r, &info FCONE); if(info != 0){error("c++ error: dpotri failed 2a\n");}
    
#ifdef _OPENMP
#pragma omp parallel for private(k, l, info, threadID, e, Omega_ii)
#endif
  for(i = 0; i < n; i++){
#ifdef _OPENMP   
    threadID = omp_get_thread_num();
#endif

      //R(i, S*) 1xr
      for(l = 0; l < r; l++){
	e = dist2(coords[i], coords[n+i], knots[l], knots[r+l]);
	R_iS[r*threadID+l] = spCor(e, phi, nu, covModel, &bk[threadID*nb]);
      }
      
      //J_i 1xr
      F77_NAME(dsymv)(&lower, &r, &one, R_SInv, &r, &R_iS[r*threadID], &inc, &zero, &J_i[r*threadID], &inc FCONE);
      
      //Omega_ii 1x1
      F77_NAME(dsymv)(&lower, &r, &one, R_S, &r, &J_i[r*threadID], &inc, &zero, &tmp_r[r*threadID], &inc FCONE);
      Omega_ii = 1.0 + alpha - F77_NAME(ddot)(&r, &tmp_r[r*threadID], &inc, &J_i[r*threadID], &inc);
      
      if(i > 0){

      	//R(N_i, S*) m(i)xr
      	for(k = 0; k < nnIndxLU[n+i]; k++){
      	  for(l = 0; l < r; l++){
      	    e = dist2(coords[nnIndx[nnIndxLU[i]+k]], coords[n+nnIndx[nnIndxLU[i]+k]], knots[l], knots[r+l]);
      	    R_NiS[mr*threadID+l*nnIndxLU[n+i]+k] = spCor(e, phi, nu, covModel, &bk[threadID*nb]);
      	  }
      	}
	
      	//J_Ni = R(N_i, S*) R(S*)^{-1} m(i)xr
      	F77_NAME(dsymm)(&rside, &lower, &nnIndxLU[n+i], &r, &one, R_SInv, &r, &R_NiS[mr*threadID], &nnIndxLU[n+i], &zero, &J_Ni[mr*threadID], &nnIndxLU[n+i] FCONE FCONE);

      	//J_Ni R(S*) J_Ni^T m(i)xm(i)
      	F77_NAME(dsymm)(&rside, &lower, &nnIndxLU[n+i], &r, &one, R_S, &r, &J_Ni[mr*threadID], &nnIndxLU[n+i], &zero, &tmp_mr[mr*threadID], &nnIndxLU[n+i] FCONE FCONE);
      	F77_NAME(dgemm)(&ntran, &ytran, &nnIndxLU[n+i], &nnIndxLU[n+i], &r, &one, &tmp_mr[mr*threadID], &nnIndxLU[n+i], &J_Ni[mr*threadID], &nnIndxLU[n+i], &zero, &tmp_mm[mm*threadID], &nnIndxLU[n+i] FCONE FCONE);

      	//R_Ni + alpha I m(i)xm(i)
      	for(k = 0; k < nnIndxLU[n+i]; k++){
      	  for(l = 0; l <= k; l++){
      	    e = dist2(coords[nnIndx[nnIndxLU[i]+k]], coords[n+nnIndx[nnIndxLU[i]+k]], coords[nnIndx[nnIndxLU[i]+l]], coords[n+nnIndx[nnIndxLU[i]+l]]);
      	    R_Ni[mm*threadID+l*nnIndxLU[n+i]+k] = spCor(e, phi, nu, covModel, &bk[threadID*nb]);
      	    if(l == k){
      	      R_Ni[mm*threadID+l*nnIndxLU[n+i]+k] += alpha;
      	    }
      	  }
      	}

      	//Omega_i mxm
      	for(k = 0; k < nnIndxLU[n+i]*nnIndxLU[n+i]; k++){
      	  Omega_i[mm*threadID+k] = R_Ni[mm*threadID+k] - tmp_mm[mm*threadID+k];
      	}

      	//Omega_iNi 1xm 
      	//F77_NAME(dsymv)(&lower, &r, &one, R_S, &r, J_i, &inc, &zero, tmp_r, &inc FCONE); already calculated above
      	F77_NAME(dgemv)(&ntran, &nnIndxLU[n+i], &r, &one, &J_Ni[mr*threadID], &nnIndxLU[n+i], &tmp_r[r*threadID], &inc, &zero, &tmp_m[m*threadID], &inc FCONE);
	
      	for(k = 0; k < nnIndxLU[n+i]; k++){
      	  e = dist2(coords[i], coords[n+i], coords[nnIndx[nnIndxLU[i]+k]], coords[n+nnIndx[nnIndxLU[i]+k]]);
      	  Omega_iNi[m*threadID+k] = spCor(e, phi, nu, covModel, &bk[threadID*nb]) - tmp_m[m*threadID+k]; //R_iNi - J_i R_S J_Ni^T
      	}
	
      	//B_i and F_i
      	F77_NAME(dpotrf)(&lower, &nnIndxLU[n+i], &Omega_i[mm*threadID], &nnIndxLU[n+i], &info FCONE); if(info != 0){error("c++ error: dpotrf failed 3a\n");}
	F77_NAME(dpotri)(&lower, &nnIndxLU[n+i], &Omega_i[mm*threadID], &nnIndxLU[n+i], &info FCONE); if(info != 0){error("c++ error: dpotri failed 4a\n");}
       	F77_NAME(dsymv)(&lower, &nnIndxLU[n+i], &one, &Omega_i[mm*threadID], &nnIndxLU[n+i], &Omega_iNi[m*threadID], &inc, &zero, &B[nnIndxLU[i]], &inc FCONE);
       	F[i] = Omega_ii - F77_NAME(ddot)(&nnIndxLU[n+i], &B[nnIndxLU[i]], &inc, &Omega_iNi[m*threadID], &inc);
      }else{
      	B[i] = 0;
      	F[i] = Omega_ii;
      }

    }
    
}

extern "C" {
  
  SEXP cSLGP(SEXP y_r, SEXP X_r, SEXP p_r, SEXP n_r, SEXP r_r, SEXP coords_r, SEXP knots_r, SEXP thetaAlpha_r, SEXP nnIndx_r, SEXP nnIndxLU_r,
	     SEXP X0_r, SEXP coords0_r, SEXP n0_r, SEXP nnIndx0_r, SEXP g_r, 
	     SEXP m_r, SEXP sigmaSqIG_r, SEXP covModel_r, SEXP nThreads_r, SEXP verbose_r, SEXP getXStr_r){
    
    int h, i, j, k, l, s, info, nProtect=0;
    const int inc = 1;
    const double one = 1.0;
    const double negOne = -1.0;
    const double zero = 0.0;
    char const *lower = "L";
    char const *ntran = "N";
    char const *ytran = "T";
    char const *rside = "R";
     
    //get args
    double *y = REAL(y_r);
    double *X = REAL(X_r);
    int p = INTEGER(p_r)[0];
    int n = INTEGER(n_r)[0];
    int r = INTEGER(r_r)[0];
    double *coords = REAL(coords_r);
    double *knots = REAL(knots_r);

    double *X0 = REAL(X0_r);
    int n0 = INTEGER(n0_r)[0];
    double *coords0 = REAL(coords0_r);
    int *nnIndx0 = INTEGER(nnIndx0_r);
    int g = INTEGER(g_r)[0];
    double *thetaAlpha = REAL(thetaAlpha_r); //nTheta x g
    int *nnIndx = INTEGER(nnIndx_r);
    int *nnIndxLU = INTEGER(nnIndxLU_r);
    
    int m = INTEGER(m_r)[0];
    double sigmaSqIGa = REAL(sigmaSqIG_r)[0]; double sigmaSqIGb = REAL(sigmaSqIG_r)[1];
    int covModel = INTEGER(covModel_r)[0];
    std::string corName = getCorName(covModel);
    int nThreads = INTEGER(nThreads_r)[0];
    int verbose = INTEGER(verbose_r)[0];
    int getXStr = INTEGER(getXStr_r)[0];
  
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
      Rprintf("SLGP Conjugate model fit with %i observations.\n\n", n);
      Rprintf("Number of knots %i.\n\n", r);
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
    	  
    //other stuff
    int nIndx = static_cast<int>(static_cast<double>(1+m)/2*m+(n-m-1)*m);
    double *B = (double *) R_alloc(nIndx, sizeof(double));
    double *F = (double *) R_alloc(n, sizeof(double));

    int mm = m*m;
    int mr = r*m;
    int rr = r*r;
    
    double *R_iS = (double *) R_alloc(nThreads*r, sizeof(double));
    double *R_NiS = (double *) R_alloc(nThreads*mr, sizeof(double));
    double *R_Ni = (double *) R_alloc(nThreads*mm, sizeof(double));
    double *R_S = (double *) R_alloc(rr, sizeof(double));
    double *R_SInv = (double *) R_alloc(rr, sizeof(double));
    double *J_i = (double *) R_alloc(nThreads*r, sizeof(double));
    double *J_Ni = (double *) R_alloc(nThreads*mr, sizeof(double));
    double *Omega_i = (double *) R_alloc(nThreads*mm, sizeof(double));
    double *Omega_iNi = (double *) R_alloc(nThreads*m, sizeof(double));
    double *w = (double *) R_alloc(nThreads*m, sizeof(double));
    double *tmp_r = (double *) R_alloc(nThreads*r, sizeof(double));
    double *tmp_m = (double *) R_alloc(nThreads*m, sizeof(double));
    double *tmp_mr = (double *) R_alloc(nThreads*mr, sizeof(double));
    double *tmp_mm = (double *) R_alloc(nThreads*mm, sizeof(double));
 
    int np = n*p;
    int q = p+r;
    int qq = q*q;
    int mq = m*q;
    int nq = n*q;
        
    double *XStr = (double *) R_alloc(n*q, sizeof(double));
    F77_NAME(dcopy)(&np, X, &inc, XStr, &inc);
    double *tmp_q = (double *) R_alloc(nThreads*q, sizeof(double));
    double *tmp_q2 = (double *) R_alloc(nThreads*q, sizeof(double));
    double *tmp_qq = (double *) R_alloc(nThreads*qq, sizeof(double));
    double *tmp_mq = (double *) R_alloc(nThreads*mq, sizeof(double));
    
    //return stuff
    SEXP beta_r, ab_r, bBInv_r, bb_r, y0Hat_r, y0HatVar_r;
    PROTECT(beta_r = allocMatrix(REALSXP, q, g)); nProtect++;
    PROTECT(ab_r = allocMatrix(REALSXP, nTheta, g)); nProtect++;
    PROTECT(bBInv_r = allocMatrix(REALSXP, qq, g)); nProtect++;
    PROTECT(bb_r = allocMatrix(REALSXP, q, g)); nProtect++;
    PROTECT(y0Hat_r = allocMatrix(REALSXP, n0, g)); nProtect++;
    PROTECT(y0HatVar_r = allocMatrix(REALSXP, n0, g)); nProtect++;

    SEXP XStr_r;
    if(getXStr){//only really want to save and return these when it is the final run, i.e., g = 1, otherwise it could be HUGE if n is large
      PROTECT(XStr_r = allocMatrix(REALSXP, nq, g)); nProtect++; 
    }
	
    double *beta = REAL(beta_r);
    double *ab = REAL(ab_r);
    double *y0Hat = REAL(y0Hat_r);
    double *y0HatVar = REAL(y0HatVar_r);
    double *betaVStrInv = (double *) R_alloc(qq, sizeof(double)); zeros(betaVStrInv, qq);
    for(i = 0; i < p; i++){
      betaVStrInv[i*q+i] = 1.0/100000;
    }
    
    double phi = 0, nu = 0, alpha = 0, e, Omega_ii;
    int threadID = 0;
   
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
#ifdef _OPENMP
      omp_set_num_threads(nThreads);
#endif
      updateConjBF(B, F, R_iS, R_NiS, R_Ni, R_S, R_SInv, J_i, J_Ni, Omega_i, Omega_iNi, tmp_r, tmp_m, tmp_mr, tmp_mm, coords, knots, nnIndx, nnIndxLU, n, m, r, phi, alpha, nu, covModel, bk, nuMax);

      //make V_beta*^{-1} betaVStrInv (add R(S*)^{-1} to the lower diag)
      for(i = 0; i < r; i++){
      	for(l = 0; l < r; l++){
      	  betaVStrInv[p*q+q*i+p+l] = R_SInv[i*r+l];
      	}
      }
      
      //make X* = (X,J)
      //note R_SInv is computed in updateConjBF
      for(i = 0; i < n; i++){
      	for(l = 0; l < r; l++){
      	  e = dist2(coords[i], coords[n+i], knots[l], knots[r+l]);
      	  tmp_r[l] = spCor(e, phi, nu, covModel, bk);//R(S,S*)
      	}

      	F77_NAME(dsymv)(lower, &r, &one, R_SInv, &r, tmp_r, &inc, &zero, &XStr[p*n+i], &n FCONE);
      }

      //estimation
      for(i = 0; i < q; i++){
      	tmp_q[i] = Q(B, F, &XStr[n*i], y, n, nnIndx, nnIndxLU); //X*'Omega^{-1}y
      	for(j = 0; j <= i; j++){
      	  tmp_qq[j*q+i] = Q(B, F, &XStr[n*j], &XStr[n*i], n, nnIndx, nnIndxLU);//X*'Omega^{-1}X*
      	}
      }

      //B = V_beta*^{-1} + X*'Omega^{-1}X*
      for(i = 0; i < qq; i++){
      	tmp_qq[i] = betaVStrInv[i]+tmp_qq[i];
      }

      //b = V_beta*^{-1}mu_beta* + X*'Omega^{-1}y but mu_beta* will be zero

      //V = inv(B), V = inv(tmp_qq), V = tmp_qq after dpotri
      F77_NAME(dpotrf)(lower, &q, tmp_qq, &q, &info FCONE); if(info != 0){error("c++ error: dpotrf failed 5\n");}
      F77_NAME(dpotri)(lower, &q, tmp_qq, &q, &info FCONE); if(info != 0){error("c++ error: dpotri failed 6\n");}

      //g = solve(B, v), g = solve(tmp_pp, tmp_p), g = beta
      F77_NAME(dsymv)(lower, &q, &one, tmp_qq, &q, tmp_q, &inc, &zero, &beta[k*q], &inc FCONE);

      //a 
      ab[k*2] = sigmaSqIGa + 1.0*n/2.0;
      
      //b
      ab[k*2+1] = sigmaSqIGb + 0.5*(Q(B, F, y, y, n, nnIndx, nnIndxLU) - F77_NAME(ddot)(&q, &beta[k*q], &inc, tmp_q, &inc));

      //save B^{-1} to calculate bB^{−1}/(a − 1) on return also keep \bb exact samples
      F77_NAME(dcopy)(&qq, tmp_qq, &inc, &REAL(bBInv_r)[k*qq], &inc);
      F77_NAME(dcopy)(&q, tmp_q, &inc, &REAL(bb_r)[k*q], &inc);
      if(getXStr){
	F77_NAME(dcopy)(&nq, XStr, &inc, &REAL(XStr_r)[k*nq], &inc);
      }
      
      //prediction
#ifdef _OPENMP
#pragma omp parallel for private(j, l, info, threadID, e, Omega_ii)
#endif
      for(i = 0; i < n0; i++){
#ifdef _OPENMP   
 	threadID = omp_get_thread_num();
#endif		
	//R(i, S*) 1xr
	for(l = 0; l < r; l++){
	  e = dist2(coords0[i], coords0[n0+i], knots[l], knots[r+l]);
	  R_iS[r*threadID+l] = spCor(e, phi, nu, covModel, &bk[threadID*nb]);
	}
      
	//J_i 1xr
	F77_NAME(dsymv)(lower, &r, &one, R_SInv, &r, &R_iS[r*threadID], &inc, &zero, &J_i[r*threadID], &inc FCONE);

	//Omega_ii 1x1
	F77_NAME(dsymv)(lower, &r, &one, R_S, &r, &J_i[r*threadID], &inc, &zero, &tmp_r[r*threadID], &inc FCONE);
	Omega_ii = 1.0 + alpha - F77_NAME(ddot)(&r, &tmp_r[r*threadID], &inc, &J_i[r*threadID], &inc);
	
	//R(N_i, S*) mxr
	for(j = 0; j < m; j++){
	  for(l = 0; l < r; l++){
	    e = dist2(coords[nnIndx0[j*n0+i]], coords[n+nnIndx0[j*n0+i]], knots[l], knots[r+l]);
	    R_NiS[mr*threadID+l*m+j] = spCor(e, phi, nu, covModel, &bk[threadID*nb]);
	  }
	}
	
	//J_Ni = R(N_i, S*) R(S*)^{-1} mxr
	F77_NAME(dsymm)(rside, lower, &m, &r, &one, R_SInv, &r, &R_NiS[mr*threadID], &m, &zero, &J_Ni[mr*threadID], &m FCONE FCONE);

	//for Omega_i mxm
	//J_Ni R(S*) J_Ni^T mxm
	F77_NAME(dsymm)(rside, lower, &m, &r, &one, R_S, &r, &J_Ni[mr*threadID], &m, &zero, &tmp_mr[mr*threadID], &m FCONE FCONE);
	F77_NAME(dgemm)(ntran, ytran, &m, &m, &r, &one, &tmp_mr[mr*threadID], &m, &J_Ni[mr*threadID], &m, &zero, &tmp_mm[mm*threadID], &m FCONE FCONE);

	//R_Ni + alpha I mxm
	for(j = 0; j < m; j++){
	  for(l = 0; l <= j; l++){
	    e = dist2(coords[nnIndx0[j*n0+i]], coords[n+nnIndx0[j*n0+i]], coords[nnIndx0[l*n0+i]], coords[n+nnIndx0[l*n0+i]]);
	    R_Ni[mm*threadID+l*m+j] = spCor(e, phi, nu, covModel, &bk[threadID*nb]);
	    if(l == j){
	      R_Ni[mm*threadID+l*m+j] += alpha;
	    }
	  }
	}

	for(j = 0; j < mm; j++){
	  Omega_i[mm*threadID+j] = R_Ni[mm*threadID+j] - tmp_mm[mm*threadID+j];
	}

	//Omega_iNi 1xm (this is z)
	F77_NAME(dsymv)(lower, &r, &one, R_S, &r, &J_i[r*threadID], &inc, &zero, &tmp_r[r*threadID], &inc FCONE); 
	F77_NAME(dgemv)(ntran, &m, &r, &one, &J_Ni[mr*threadID], &m, &tmp_r[r*threadID], &inc, &zero, &tmp_m[m*threadID], &inc FCONE);//(J_i R_S) J_Ni^T
	
	for(j = 0; j < m; j++){
	  e = dist2(coords0[i], coords0[n0+i], coords[nnIndx0[j*n0+i]], coords[n+nnIndx0[j*n0+i]]);//R_iNi
	  Omega_iNi[m*threadID+j] = spCor(e, phi, nu, covModel, &bk[threadID*nb]) - tmp_m[m*threadID+j]; //R_iNi - J_i R_S J_Ni^T
	}
	
	//solve for w
	F77_NAME(dpotrf)(lower, &m, &Omega_i[mm*threadID], &m, &info FCONE); if(info != 0){error("c++ error: dpotrf failed 3\n");}
	F77_NAME(dpotri)(lower, &m, &Omega_i[mm*threadID], &m, &info FCONE); if(info != 0){error("c++ error: dpotri failed 4\n");}
	F77_NAME(dsymv)(lower, &m, &one, &Omega_i[mm*threadID], &m, &Omega_iNi[m*threadID], &inc, &zero, &w[m*threadID], &inc FCONE);
	
      	//make hat(y)
      	for(j = 0; j < m; j++){
	  tmp_m[m*threadID+j] = y[nnIndx0[j*n0+i]] - F77_NAME(ddot)(&q, &XStr[nnIndx0[j*n0+i]], &n, &beta[k*q], &inc);
      	}

	//ith row of X0*
	F77_NAME(dcopy)(&p, &X0[i], &n0, &tmp_q[q*threadID], &inc);
	F77_NAME(dcopy)(&r, &J_i[r*threadID], &inc, &tmp_q[q*threadID+p], &inc);//should be (28)
	
      	y0Hat[k*n0+i] = F77_NAME(ddot)(&q, &tmp_q[q*threadID], &inc, &beta[k*q], &inc) + F77_NAME(ddot)(&m, &w[m*threadID], &inc, &tmp_m[m*threadID], &inc);

	//make u
	for(j = 0; j < m; j++){
	  F77_NAME(dcopy)(&q, &XStr[nnIndx0[j*n0+i]], &n, &tmp_mq[mq*threadID+j], &m);
	}

	F77_NAME(dgemv)(ytran, &m, &q, &one, &tmp_mq[mq*threadID], &m, &w[m*threadID], &inc, &zero, &tmp_q2[q*threadID], &inc FCONE);//X*[N(s0),]w //note typo in FDB et al. 2019, dot should be dgemv

	for(j = 0; j < q; j++){
	  tmp_q[q*threadID+j] = tmp_q[q*threadID+j] - tmp_q2[q*threadID+j];
	}

	//make v_y and var(y)
	F77_NAME(dsymv)(lower, &q, &one, tmp_qq, &q, &tmp_q[q*threadID], &inc, &zero, &tmp_q2[q*threadID], &inc FCONE);

	y0HatVar[k*n0+i] = ab[k*2+1] * (F77_NAME(ddot)(&q, &tmp_q[q*threadID], &inc, &tmp_q2[q*threadID], &inc) + Omega_ii - F77_NAME(ddot)(&m, &w[m*threadID], &inc, &Omega_iNi[m*threadID], &inc))/(ab[k*2]-1.0);
       
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

    if(getXStr){
      nResultListObjs++;
    }
    
    if(n0 > 0){
      nResultListObjs += 2;
    }
    
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

    i = 3;
    if(getXStr){
      i++;
      SET_VECTOR_ELT(result_r, i, XStr_r);
      SET_VECTOR_ELT(resultName_r, i, mkChar("X.str"));
    }

    if(n0 > 0){
      i++;
      SET_VECTOR_ELT(result_r, i, y0Hat_r);
      SET_VECTOR_ELT(resultName_r, i, mkChar("y.0.hat"));

      i++;
      SET_VECTOR_ELT(result_r, i, y0HatVar_r);
      SET_VECTOR_ELT(resultName_r, i, mkChar("y.0.hat.var"));
    }

    namesgets(result_r, resultName_r);
    
    //unprotect
    UNPROTECT(nProtect);
    
    return(result_r);
  }
}
