// hmmR.c - R interface to C implementation of HMM routines
//
// by Artem Sokolov

#include "matrix.h"
#include "hmm.h"

// Interface to HMM_obs in hmm.c
//
// o_B must be N*K long
// o_status must be 1 long
// o_info must be 1 long
//
// Meaning of o_status:
// 0 - Everything is OK
// 1 - Covariance matrix is singular for state o_info
void HMMobsR( double* i_x, double* i_pi, double* i_A, double* i_mu, double* i_sigma,
	      int* i_N, int* i_K, int* i_p, double* o_B, int* o_status, int* o_info )
{
  int N = *i_N;
  int K = *i_K;
  int p = *i_p;

  // Create local storage
  double** x = make2D( N, p );
  double** A = make2D( K, K );
  double** mu = make2D( K, p );
  double*** sigma = make3D( K, p, p );
  double** B = make2D( N, K );

  // Copy parametesr
  cp_vec_mat( N, p, i_x, x );
  cp_vec_mat( K, K, i_A, A );
  cp_vec_mat( K, p, i_mu, mu );
  cp_vec_3D( K, p, p, i_sigma, sigma );

  // Call the C code
  *o_status = HMM_obs( x, i_pi, A, mu, sigma, N, K, p, B, o_info );

  // Copy the results
  cp_mat_vec( N, K, B, o_B );

  // Free memory
  free2D( x, N );
  free2D( A, K );
  free2D( mu, K );
  free3D( sigma, K, p );
  free2D( B, N );
}

// Interface to HMM_fwbk in hmm.c
//
// o_alpha must be N*K long
// o_beta must be N*K long
// o_c must be N long
// o_gamma must be N*K long
// o_xi must be (N-1)*K*K long
// o_status must be 1 long
//
// Meaning of o_status:
// 0 - Everything is OK
// 1 - Covariance matrix for one of the states is singular
// 2 - Some points are not accounted for by any state
void HMMfwbkR( double* i_x, double* i_pi, double* i_A, double* i_mu, double* i_sigma,
	       int* i_N, int* i_K, int* i_p, double* o_alpha, double* o_beta, double* o_c,
	       double* o_gamma, double* o_xi, int* o_status )
{
  int N = *i_N;
  int K = *i_K;
  int p = *i_p;

  // Create local storage
  double** x = make2D( N, p );
  double** A = make2D( K, K );
  double** mu = make2D( K, p );
  double*** sigma = make3D( K, p, p );
  double** alpha = make2D( N, K );
  double** beta = make2D( N, K );
  double** gamma = make2D( N, K );
  double*** xi = make3D( N-1, K, K );

  // Copy parameters
  cp_vec_mat( N, p, i_x, x );
  cp_vec_mat( K, K, i_A, A );
  cp_vec_mat( K, p, i_mu, mu );
  cp_vec_3D( K, p, p, i_sigma, sigma );

  // Call the C code
  *o_status = HMM_fwbk( x, i_pi, A, mu, sigma, N, K, p, alpha, beta, o_c, gamma, xi );

  // Copy the results
  cp_mat_vec( N, K, alpha, o_alpha );
  cp_mat_vec( N, K, beta, o_beta );
  cp_mat_vec( N, K, gamma, o_gamma );
  cp_3D_vec( N-1, K, K, xi, o_xi );

  // Free memory
  free2D( x, N );
  free2D( A, K );
  free2D( mu, K );
  free3D( sigma, K, p );
  free2D( alpha, N );
  free2D( beta, N );
  free2D( gamma, N );
  free3D( xi, N-1, K );
}

// Interface to HMM_learn in hmm.c
//
// o_pi must be of length K
// o_A must be of length K*K
// o_mu must be of length K*p
// o_sigma must be of length K*p*p
// o_status must be of length 1
// 
// o_status is:
//	0 - Everything is OK
//	1 - Covariance matrix is singular
//	2 - Some samples were not accounted for by any state
void HMMlearnR( double* i_x, double* i_v, int* i_nSeq, int* i_mSeq, 
		int* i_N, int* i_M, int* i_K, int* i_p,
		double* i_pi_init, double* i_A_init, double* i_mu_init, double* i_sigma_init,
		double* i_tol, double* i_LLstop, int* i_min_iter, int* i_max_iter,
		int* i_debug_output, double* o_pi, double* o_A, double* o_mu, double* o_sigma, int* o_status )
{

  int i, j, k;
  int x_offset;
  int v_offset;
  double*** x = NULL;
  double*** v = NULL;

  // Retrieve parameters
  int nSeqs = *i_nSeq; int mSeqs = *i_mSeq; int p = *i_p; int K = *i_K;

  // Create local storage and copy and reshape the arguments
  double** A = make2D( K, K );
  double** mu = make2D( K, p );
  double*** sigma = make3D( K, p, p );

  // x
  x = (double***)malloc( nSeqs * sizeof(double**) );
  int N_total = 0; for( k = 0; k < nSeqs; k++ ) N_total += i_N[k];
  x_offset = 0;
  for( k = 0; k < nSeqs; k++ )
    {
      x[k] = make2D( i_N[k], p );
      for( j = 0; j < p; j++ ) for( i = 0; i < i_N[k]; i++ )
	x[k][i][j] = (i_x+x_offset)[ N_total*j+i ];
      x_offset += i_N[k];
    }
  
  // v
  if( i_v != NULL )
    {
      v = (double***)malloc( mSeqs * sizeof(double**) );
      int M_total = 0; for( k = 0; k < mSeqs; k++ ) M_total += i_M[k];
      v_offset = 0;
      for( k = 0; k < mSeqs; k++ )
	{
	  v[k] = make2D( i_M[k], p );
	  for( j = 0; j < p; j++ ) for( i = 0; i < i_M[k]; i++ )
	    v[k][i][j] = (i_v+v_offset)[ M_total*j+i ];
	  v_offset += i_M[k];
	}
    }

  // Copy initial estimates
  cp_vec_vec( K, i_pi_init, o_pi );
  cp_vec_mat( K, K, i_A_init, A );
  cp_vec_mat( K, p, i_mu_init, mu );
  cp_vec_3D( K, p, p, i_sigma_init, sigma );

  // Call the C function
  *o_status = HMM_learn( x, v, nSeqs, mSeqs, i_N, i_M, K, p, 
			 o_pi, A, mu, sigma,
			 *i_tol, *i_LLstop, *i_min_iter, *i_max_iter, *i_debug_output );

  // Copy the results
  cp_mat_vec( K, K, A, o_A );
  cp_mat_vec( K, p, mu, o_mu );
  cp_3D_vec( K, p, p, sigma, o_sigma );

  // Free up memory
  for( i = 0; i < nSeqs; i++ ) free2D( x[i], i_N[i] ); free( x );
  if( i_v != NULL ) { for( i = 0; i < mSeqs; i++ ) free2D( v[i], i_M[i] ); free( v ); }

  free2D( A, K );
  free2D( mu, K );
  free3D( sigma, K, p );
}
