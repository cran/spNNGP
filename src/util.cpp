#include <string>
#include <limits>
#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>
#include <R_ext/BLAS.h>
#include <R_ext/Utils.h>
#include "util.h"

#ifdef _OPENMP
#include <omp.h>
#endif

void zeros(double *a, int n){
  for(int i = 0; i < n; i++)
    a[i] = 0.0;
}

void mvrnorm(double *des, double *mu, double *cholCov, int dim){
  
  int i;
  int inc = 1;
  double one = 1.0;
  double zero = 0.0;
  
  for(i = 0; i < dim; i++){
    des[i] = rnorm(0, 1);
  }
 
  F77_NAME(dtrmv)("L", "N", "N", &dim, cholCov, &dim, des, &inc);
  F77_NAME(daxpy)(&dim, &one, mu, &inc, des, &inc);
}

double logit(double theta, double a, double b){
  return log((theta-a)/(b-theta));
}

double logitInv(double z, double a, double b){
  return b-(b-a)/(1+exp(z));
}

double dist2(double &a1, double &a2, double &b1, double &b2){
  return(sqrt(pow(a1-b1,2)+pow(a2-b2,2)));
}

void getNNIndx(int i, int m, int &iNNIndx, int &iNN){
  
  if(i == 0){
    iNNIndx = 0;//this should never be accessed
    iNN = 0;
    return;
  }else if(i < m){
    iNNIndx = static_cast<int>(static_cast<double>(i)/2*(i-1));
    iNN = i;
    return;
  }else{
    iNNIndx = static_cast<int>(static_cast<double>(m)/2*(m-1)+(i-m)*m);
    iNN = m;
    return;
  } 
}

void mkUIndx(int n, int m, int* nnIndx, int* uIndx, int* uIndxLU){ 
  
  int iNNIndx, iNN, i, j, k, l, h;
  
  for(i = 0, l = 0; i < n; i++){    
    uIndxLU[i] = l; 
    for(j = 0, h = 0; j < n; j++){   
      getNNIndx(j, m, iNNIndx, iNN);  
      for(k = 0; k < iNN; k++){      	
	if(nnIndx[iNNIndx+k] == i){
	  uIndx[l+h] = j;
	  h++;
	}    
      }
    }
    l += h;
    uIndxLU[n+i] = h; 
  }
}

std::string getCorName(int i){

  if(i == 0){
    return "exponential";
  }else if(i == 1){
    return "spherical";
  }else if(i == 2){
    return "matern";
  }else if(i == 3){
    return "gaussian";
  }else{
    error("c++ error: cov.model is not correctly specified");
  }
  
}

double spCor(double &D, double &phi, double &nu, int &covModel, double *bk){

  //0 exponential
  //1 spherical
  //2 matern
  //3 gaussian
  
  if(covModel == 0){//exponential
    
    return exp(-phi*D);
    
  }else if(covModel == 1){//spherical
    
    if(D > 0 && D <= 1.0/phi){
      return 1.0 - 1.5*phi*D + 0.5*pow(phi*D,3);
    }else if(D >= 1.0/phi){
      return 0.0;
    }else{
      return 1.0;
    }
  }else if(covModel == 2){//matern
    
    //(d*phi)^nu/(2^(nu-1)*gamma(nu))*pi/2*(besselI(d*phi,-nu)-besselI(d*phi, nu))/sin(nu*pi), or
    //(d*phi)^nu/(2^(nu-1)*gamma(nu))*besselK(x=d*phi, nu=nu)
    
    if(D*phi > 0.0){
      return pow(D*phi, nu)/(pow(2, nu-1)*gammafn(nu))*bessel_k_ex(D*phi, nu, 1.0, bk);//thread safe bessel
    }else{
      return 1.0;
    } 
  }else if(covModel == 3){//gaussian
    
    return exp(-1.0*(pow(phi*D,2)));
      
  }else{
    error("c++ error: cov.model is not correctly specified");
  }
}

//which index of b equals a, where b is of length n
int which(int a, int *b, int n){
  int i;
  for(i = 0; i < n; i++){
    if(a == b[i]){
      return(i);
    }
  }

  error("c++ error: which failed");
  return -9999;
}

//Description: computes the quadratic term.
double Q(double *B, double *F, double *u, double *v, int n, int *nnIndx, int *nnIndxLU){
  
  double a, b, q = 0;
  int i, j;
  
#ifdef _OPENMP
#pragma omp parallel for private(a, b, j) reduction(+:q)
#endif  
  for(i = 0; i < n; i++){
    a = 0;
    b = 0;
    for(j = 0; j < nnIndxLU[n+i]; j++){
      a += B[nnIndxLU[i]+j]*u[nnIndx[nnIndxLU[i]+j]];
      b += B[nnIndxLU[i]+j]*v[nnIndx[nnIndxLU[i]+j]];
    }
    q += (u[i] - a)*(v[i] - b)/F[i];
  }
  
  return(q);
}


void printMtrx(double *m, int nRow, int nCol){

  int i, j;

  for(i = 0; i < nRow; i++){
    Rprintf("\t");
    for(j = 0; j < nCol; j++){
      Rprintf("%.10f\t", m[j*nRow+i]);
    }
    Rprintf("\n");
  }
}


void printMtrxInt(int *m, int nRow, int nCol){

  int i, j;

  for(i = 0; i < nRow; i++){
    Rprintf("\t");
    for(j = 0; j < nCol; j++){
      Rprintf("%i\t", m[j*nRow+i]);
    }
    Rprintf("\n");
  }

}
