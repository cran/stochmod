// gmm.c - Gaussian Mixture Model routines
//
// by Artem Sokolov

#include <Rmath.h>

#include "matrix.h"
#include "stats.h"

// Computes cluster responsibilities and posterior log-likelihood
//
// N - number of samples
// p - observation dimensionality
// K - number of components in the model
//
// Inputs:
//      x - [N x p] matrix of observations
//	mu - [K x p] matrix of component means
//	sigma - [K x p x p] array of covariance matrices for each component
//	pi - [K x 1] vector of mixture coefficients (has to sum up to 1)
//
// Outputs:
//	resp - [N x K] matrix of responsibilities
//	LL   - [N x 1] vector of log likelihoods
//
// Returns:
//	0 - everything is OK
//	1 - covariance matrix for cluster info[0] is singular
//	2 - points which( info == TRUE ) were not assigned to any cluster
int GMM_resp( double** x, double** mu, double*** sigma, double* pi,
	      int N, int p, int K,
	      double** resp, double* LL, int* info )
{
  // Create local storage
  int i, j;
  int status;
  double* p_x = make1D( N );
  double** post_prob = make2D( K, N );

  // Normalize pi
  mult_const_vec( K, 1.0 / vec_sum( K, pi ), pi );

  // Compute posterior probabilities using Bayes rule
  for( i = 0; i < K; i++ )
    {
      status = dmvnorm( N, p, x, mu[i], sigma[i], post_prob[i] );

      // Check for numerical stability
      if( status > 0 )
	{
	  free( p_x );
	  free2D( post_prob, K );
	  info[0] = i+1;
	  return status;
	}

      mult_const_vec( N, pi[i], post_prob[i] );
    }

  // Compute responsibilities by scaling posteriors
  for( i = 0; i < N; i++ ) p_x[i] = 0.0;
  for( i = 0; i < K; i++ )
    for( j = 0; j < N; j++ )
      p_x[j] += post_prob[i][j];

  for( i = 0; i < N; i++ )
    {
      if( p_x[i] == 0 )
	{
	  status = 2;
	  info[i] = 1;
	  for( j = 0; j < K; j++ )
	    resp[i][j] = 0.0;
	}
      else
	{
	  info[i] = 0;
	  for( j = 0; j < K; j++ )
	    resp[i][j] = post_prob[j][i] / p_x[i];
	}
    }

  // Compute log likelihood values
  for( i = 0; i < N; i++ )
    LL[i] = log( p_x[i] );

  // Free up local memory
  free( p_x );
  free2D( post_prob, K );
  return status;
}


// Expectation Maximization for multiple observation sequences
//
// Inputs:
// 	x - [N x p] matrix of training samples
// 	v - [M x p] matrix of validation samples 
// 	N - number of training samples
// 	M - number of validation samples
// 	K - Desired number of components
// 	p - Observation dimensionality
//
// Control variables:
//	tol - relative tolerance for convergence of training log-likelihood
//	LLstop - iterations will stop when training log-likelihood exceeds this value
//	max_iter - maximum number of iterations to perform (use negative number to denote infinity)
//	debug_output - amount of output produced by the algorithm
//		0 - none
//		1 - basic
//		2 - verbose
//
// Outputs:
//	pi - [K x 1] vector of mixture coefficients
// 	mu - [K x p] matrix of component means
//	sigma - [K x p x p] matrix of component covariance matrices
//
// Returns:
//	0 - Everything is OK
//	1 - Covariance matrix is singular for one of the clusters
//	2 - Some points were not assigned to a cluster during the E-step
//	3 - No points assigned to one of the clusters
int GMM_learn( double** x, double** v, int N, int M, int K, int p,
	       double* pi, double** mu, double*** sigma,
	       double cov_reg, double tol, double LLstop, 
	       int min_iter, int max_iter, int debug_output )
{
  if( debug_output > 0 )
    printf( "Currently running C version of the GMM.learn code\n" );

  // Create local storage
  int status = 0;
  int i, j, k, l;
  int iter;
  double LL, LLprev, LLbase;
  double LLval, LLval_prev;

  int* info = malloc( N*sizeof(int) );
  double* LLvec = make1D( N );
  double** resp = make2D( N, K );
  double** x_shifted = make2D( N, p );

  double* v_LLvec = NULL;
  double** v_resp = NULL;
  double* pi_prev = NULL;
  double** mu_prev = NULL;
  double*** sigma_prev = NULL;

  if( M > 0 )
    {
      v_LLvec = make1D( M );
      v_resp = make2D( M, K );
      pi_prev = make1D( K );
      mu_prev = make2D( K, p );
      sigma_prev = make3D( K, p, p );
    }

  // Run EM algorithm
  iter = 0;
  while( 1 )
    {
      // Validation log-likelihood stopping criterion
      if( M > 0 )
	{
	  GMM_resp( v, mu, sigma, pi, M, p, K, v_resp, v_LLvec, info );
	  LLval = vec_sum( M, v_LLvec );
	  if( debug_output > 2 )
	    printf( "Validation log-likelihood: %f\n", LLval_prev );

	  if( iter <= min_iter || LLval >= LLval_prev )
	    {
	      LLval_prev = LLval;
	      cp_vec_vec( K, pi, pi_prev );
	      cp_mat_mat( K, p, mu, mu_prev );
	      cp_3D_3D( K, p, p, sigma, sigma_prev );
	    }
	  else
	    {
	      if( debug_output > 0 )
		{
		  printf( "Final log likelihood of training data after %d iterations is %f\n", iter-1, LLprev );
		  printf( "Validation log likelihood: %f\n", LLval_prev );
		}
	      cp_vec_vec( K, pi_prev, pi );
	      cp_mat_mat( K, p, mu_prev, mu );
	      cp_3D_3D( K, p, p, sigma_prev, sigma );

	      status = 0;
	      break;
	    }
	}
      

      //
      // E Step - compute responsibilities
      //
      status = GMM_resp( x, mu, sigma, pi, N, p, K, resp, LLvec, info );
      if( status != 0 ) break;

      LL = vec_sum( N, LLvec );
      if( isnan(LL) || isinf(LL) )
	{
	  status = 1;
	  break;
	}
	
      if( debug_output > 2 || (debug_output > 1 && (iter%10) == 0) )
	printf( "Log likelihood at iteration %d is %f\n", iter, LL );

      // Verify stopping criteria
      if( K == 1 && iter == 1 ) break;
      if( max_iter > 0 && iter >= max_iter ) break;
      if( LL > LLstop ) break;

      if( iter <= (min_iter-1) ) LLbase = LL;
      else if( (LL-LLbase) < (1+tol)*(LLprev-LLbase) ) break;

      // 
      // M Step - compute the new parameters
      //
      for( k = 0; k < K; k++ )
	{
	  pi[k] = 0.0;
	  for( i = 0; i < N; i++ )
	    pi[k] += resp[i][k];

	  if( pi[k] == 0 )
	    {
	      status = 3;
	      break;
	    }

	  for( j = 0; j < p; j++ )
	    {
	      mu[k][j] = 0.0;
	      for( i = 0; i < N; i++ )
		mu[k][j] += x[i][j] * resp[i][k];
	      mu[k][j] /= pi[k];
	    }

	  for( i = 0; i < N; i++ )
	    for( j = 0; j < p; j++ )
	      x_shifted[i][j] = x[i][j] - mu[k][j];

	  for( i = 0; i < p; i++ )
	    for( j = 0; j < p; j++ )
	      {
		sigma[k][i][j] = 0.0;
		for( l = 0; l < N; l++ )
		  sigma[k][i][j] += resp[l][k] * x_shifted[l][i] * x_shifted[l][j];
		sigma[k][i][j] /= pi[k];
	      }

	  if( cov_reg > 0 )
	    {
	      for( i = 0; i < p; i++ )
		for( j = 0; j < p; j++ )
		  {
		    sigma[k][i][j] *= (1-cov_reg);
		    if( i == j )
		      sigma[k][i][j] += cov_reg;
		  }
	    }
	}
      if( status != 0 ) break;

      mult_const_vec( K, 1.0 / vec_sum( K, pi ), pi );

      LLprev = LL;
      iter++;
    }

  // Display final log likelihood
  if( status == 0 && debug_output > 0 )
    printf( "Final log likelihood after %d iterations is %f\n", iter, LL );

  // Free up memory
  free( info );
  free( LLvec );
  free2D( resp, N );
  free2D( x_shifted, N );

  if( M > 0 )
    {
      free( v_LLvec );
      free2D( v_resp, M );
      free( pi_prev );
      free2D( mu_prev, K );
      free3D( sigma_prev, K, p );
    }

  return status;
}
