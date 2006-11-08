// gmmR.c - functions that interface GMM routines between C and R languages
//
// by Artem Sokolov

#include <stdio.h>

#include "matrix.h"
#include "gmm.h"

// Interface to GMM_resp in gmm.c
//
// o_resp must be of length N*K
// o_LL must be of length N
// o_status must be of length 1
// o_info must be of length N
//
// Meaning of o_status:
//	0 - Everything is OK
//	1 - Covariance matrix is singular for a cluster.
//		Cluster index is stored in o_info[0].
//	2 - Samples fall outside of all clusters.
//		o_info will contain 1 for those samples, 0 for samples that were OK
void GMMrespR( 
	      double* i_x, double* i_mu, double* i_sigma, double* i_pi,
	      int* i_N, int* i_p, int* i_K,
	      double* o_resp, double* o_LL,
	      int* o_status, int* o_info
	      )
{
  // Retrieve dimension parameters
  int N = *i_N; int p = *i_p; int K = *i_K;

  // Create local storage
  double** x = make2D( N, p );
  double** mu = make2D( K, p );
  double*** sigma = make3D( K, p, p );
  double** resp = make2D( N, K );

  // Copy inputs to the local storage
  cp_vec_mat( N, p, i_x, x );
  cp_vec_mat( K, p, i_mu, mu );
  cp_vec_3D( K, p, p, i_sigma, sigma );

  // Call the C function
  *o_status = GMM_resp( x, mu, sigma, i_pi, N, p, K, resp, o_LL, o_info );

  // Copy the results from local storage
  cp_mat_vec( N, K, resp, o_resp );

  // Free up memory
  free2D( x, N ); free2D( mu, K );
  free3D( sigma, K, p ); free2D( resp, N );
}


// Interface to GMM_learn in gmm.c
//
// o_pi must be of length K
// o_mu must be of length K*p
// o_sigma must be of length K*p*p
// o_status must be of length 1
//
// Meaning of o_status:
//	0 - Everything is OK
//	1 - Covariance matrix is singular for one of the clusters
//	2 - Some points were not assigned to a cluster during the E-step
//	3 - No points assigned to one of the clusters
void GMMlearnR( double* i_x, double* i_v, int* i_N, int* i_M, int* i_K, int* i_p,
		double* i_mu_init, double* i_sigma_init, double* i_pi_init,
		double* i_cov_reg,
		double* i_tol, double* i_LLstop, 
		int* i_min_iter, int* i_max_iter, int* debug_output,
		double* o_pi, double* o_mu, double* o_sigma, int* o_status )
{
  int i, j;
  double** v = NULL;

  // Retrieve parameters
  int N = *i_N; int M = 0; int p = *i_p; int K = *i_K;

  // Create local storage and copy and reshape the arguments
  double** mu = make2D( K, p );
  double*** sigma = make3D( K, p, p );
  double** x = make2D( N, p );

  cp_vec_mat( K, p, i_mu_init, mu );
  cp_vec_3D( K, p, p, i_sigma_init, sigma );
  cp_vec_vec( K, i_pi_init, o_pi );
  cp_vec_mat( N, p, i_x, x );

  if( i_v != NULL )
    {
      M = *i_M;
      v = make2D( M, p );
      cp_vec_mat( M, p, i_v, v );
    }

  // Call the C function
  *o_status = GMM_learn( x, v, N, M, K, p, 
			 o_pi, mu, sigma,
			 *i_cov_reg, *i_tol, *i_LLstop,
			 *i_min_iter, *i_max_iter, *debug_output );

  // Copy the results
  cp_mat_vec( K, p, mu, o_mu );
  cp_3D_vec( K, p, p, sigma, o_sigma );

  // Free up memory
  free2D( mu, K );
  free3D( sigma, K, p );
  free2D( x, N );
  if( i_v != NULL) free2D( v, M );
}
