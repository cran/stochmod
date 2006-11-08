// hmm.h - headers for HMM routines
//
// by Artem Sokolov
////////////////////////////////////

#ifndef HMM_H__INCLUDED
#define HMM_H__INCLUDED

int HMM_obs( double** x, double* pi, double** A, double** mu, double*** sigma,
	     int N, int K, int p, double** B, int* info );

int HMM_fwbk( double** x, double* pi, double** A, double** mu, double*** sigma, int N, int K, int p,
	       double** alpha, double** beta, double* c, double** gamma, double*** xi );

int HMM_LL( double** x, double* pi, double** A, double** mu, double*** sigma, int N, int K, int p, double* LL );

int HMM_learn( double*** x, double*** v, int nSeqs, int mSeqs, int* N, int* M, int K, int p,
		double* pi, double** A, double** mu, double*** sigma, 
		double tol, double LLstop, int min_iter, int max_iter, int debug_output );

#endif
