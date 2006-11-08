// gmm.h - Gaussian Mixture Model routines
//
// by Artem Sokolov

#ifndef GMM_H__INCLUDED
#define GMM_H__INCLUDED

int GMM_resp( double** x, double** mu, double*** sigma, double* pi,
	      int N, int p, int K, double** resp, double* LL, int* info );

int GMM_learn( double** x, double** v, int N, int M, int K, int p,
	       double* pi, double** mu, double*** sigma,
	       double cov_reg, double tol, double LLstop, 
	       int min_iter, int max_iter, 
	       int debug_output );
#endif
