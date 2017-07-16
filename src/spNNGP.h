#include <R.h>
#include <Rinternals.h>

extern "C" {
  SEXP rNNGP(SEXP y_r, SEXP X_r, SEXP p_r, SEXP n_r, SEXP m_r, SEXP coords_r, SEXP covModel_r, 
	     SEXP sigmaSqIG_r, SEXP tauSqIG_r, SEXP phiUnif_r, SEXP nuUnif_r, 
	     SEXP betaStarting_r, SEXP sigmaSqStarting_r, SEXP tauSqStarting_r, SEXP phiStarting_r, SEXP nuStarting_r,
	     SEXP sigmaSqTuning_r, SEXP tauSqTuning_r, SEXP phiTuning_r, SEXP nuTuning_r, 
	     SEXP nSamples_r, SEXP nThreads_r, SEXP verbose_r, SEXP nReport_r);
  
  SEXP sNNGP(SEXP y_r, SEXP X_r, SEXP p_r, SEXP n_r, SEXP m_r, SEXP coords_r, SEXP covModel_r, 
	     SEXP sigmaSqIG_r, SEXP tauSqIG_r, SEXP phiUnif_r, SEXP nuUnif_r, 
	     SEXP betaStarting_r, SEXP sigmaSqStarting_r, SEXP tauSqStarting_r, SEXP phiStarting_r, SEXP nuStarting_r,
	     SEXP sigmaSqTuning_r, SEXP tauSqTuning_r, SEXP phiTuning_r, SEXP nuTuning_r, 
	     SEXP nSamples_r, SEXP nThreads_r, SEXP verbose_r, SEXP nReport_r);

  SEXP cNNGP(SEXP y_r, SEXP X_r, SEXP p_r, SEXP n_r, SEXP coords_r, SEXP thetaAlpha_r, 
	     SEXP X0_r, SEXP coords0_r, SEXP n0_r, SEXP nnIndx0_r, SEXP g_r, 
	     SEXP m_r, SEXP sigmaSqIG_r, SEXP covModel_r, SEXP nThreads_r, SEXP verbose_r);
  
  SEXP rNNGPPredict(SEXP X_r, SEXP y_r, SEXP coords_r, SEXP n_r, SEXP p_r, SEXP m_r, SEXP X0_r, SEXP coords0_r, SEXP q_r, SEXP nnIndx0_r, 
		    SEXP betaSamples_r, SEXP thetaSamples_r, SEXP nSamples_r, SEXP covModel_r, SEXP nThreads_r, SEXP verbose_r, SEXP nReport_r);
  
  SEXP sNNGPPredict(SEXP X_r, SEXP y_r, SEXP coords_r, SEXP n_r, SEXP p_r, SEXP m_r, SEXP X0_r, SEXP coords0_r, SEXP q_r, SEXP nnIndx0_r, 
		    SEXP betaSamples_r, SEXP thetaSamples_r, SEXP wSamples_r, SEXP nSamples_r, SEXP covModel_r, SEXP nThreads_r, SEXP verbose_r, SEXP nReport_r);
}
