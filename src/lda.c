// lda.c - Linear Discriminant Analysis routines
//
// by Artem Sokolov

#include <stdio.h>

#include "matrix.h"

// Trains an LDA classifier
//
// Inputs:
// 	i_x - N x p data matrix of N samples in p dimensions
// 	i_y - N x 1 vector of labels
//	i_N - number of samples
//	i_p - data dimensionality
//	i_K - number of classes
//	i_cov_reg - covariance matrix regularization (towards identity)
//	i_debug_output - level of debug output
//
// Output:
//	o_priors - K x 1 vector of Bayesian priors, computed as fraction of points from each class
//	o_means - K x p matrix of means approximated from the data
//	o_covmat - the common p x p covariance matrix
//	o_weights - K x (p+1) matrix of weights and the bias term for each of the K classes
//
// Memory allocation requirements for outputs:
//	o_priors must be of length K
//	o_means must be of length K*p
//	o_covmat must be of length p*p
//	o_weights must be of length K*(p+1)
//	o_status must be of length 1
//
// Meaning of o_status:
//	0 - Everything is OK
//	1 - The common covariance matrix is singular
void LDAtrain( 
	      double* i_x, int* i_y, int* i_N, int* i_p, int* i_K, double* i_cov_reg,
	      int* i_debug_output,
	      double* o_priors, double* o_means, double* o_covmat, double* o_weights, int* o_status
	      )
{
  if( *i_debug_output > 0 )
    printf( "Currently running C version of the LDA.train code\n" );

  // Create local storage
  int i, j, t, k;
  int* Nk = (int*)malloc( *i_K * sizeof(int) );
  double** x = make2D( *i_N, *i_p );
  double** means = make2D( *i_K, *i_p );
  double** covmat = make2D( *i_p, *i_p );
  double** weights = make2D( *i_K, (*i_p)+1 );

  // Copy inputs to the local storage
  cp_vec_mat( *i_N, *i_p, i_x, x );

  // Init the counts and means to 0
  for( k = 0; k < *i_K; k++ )
    {
      Nk[k] = 0;
      for( i = 0; i < *i_p; i++ )
	means[k][i] = 0.0;
    }

  // Traverse one sample at a time to compute means
  for( t = 0; t < *i_N; t++ )
    {
      k = i_y[t] - 1;

      // Count the sample
      Nk[k]++;

      // Add contribution to the mean
      add_vec_vec( *i_p, 1.0, 1.0, means[k], x[t], means[k] );
    }

  // Compute the means and the priors
  for( k = 0; k < *i_K; k++ )
    {
      o_priors[k] = ((double)Nk[k]) / ((double)*i_N);
      mult_const_vec( *i_p, 1.0 / ((double)Nk[k]), means[k] );
    }

  // Subtract the means from the samples
  for( t = 0; t < *i_N; t++ )
    {
      k = i_y[t] - 1;
      add_vec_vec( *i_p, 1.0, -1.0, x[t], means[k], x[t] );
    }

  // Compute the covariance matrix of the now-normalized samples
  for( i = 0; i < *i_p; i++ )
    for( j = 0; j < *i_p; j++ )
      {
	covmat[i][j] = 0.0;
	for( t = 0; t < *i_N; t++ )
	  covmat[i][j] += x[t][i] * x[t][j];
	covmat[i][j] /= (*i_N) - (*i_K);
      }

  // Regularize the covariance matrix
  if( *i_cov_reg > 0 )
    {
      for( i = 0; i < *i_p; i++ )
	for( j = 0; j < *i_p; j++ )
	  {
	    covmat[i][j] *= (1 - *i_cov_reg);
	    if( i == j )
	      covmat[i][j] += *i_cov_reg;
	  }
    }

  // Copy the covariance matrix from local storage to output
  cp_mat_vec( *i_p, *i_p, covmat, o_covmat );

  // Invert the covariance matrix
  *o_status = covmat_inverse_inplace( *i_p, covmat );
  if( *o_status > 0 )
    {
      free( Nk );
      free2D( x, *i_N );
      free2D( means, *i_K );
      free2D( covmat, *i_p );
      free2D( weights, *i_K );
      return;
    }

  // Compute the weights matrix
  for( k = 0; k < *i_K; k++ )
    {
      mult_mat_vec( *i_p, *i_p, covmat, means[k], weights[k] );
      weights[k][(*i_p)] = -0.5 * dot_vec_vec( *i_p, means[k], weights[k] ) + log( o_priors[k] );
    }

  // Copy the local storage to outputs
  cp_mat_vec( *i_K, *i_p, means, o_means );
  cp_mat_vec( *i_K, (*i_p)+1, weights, o_weights );

  // Free up memory
  free( Nk );
  free2D( x, *i_N );
  free2D( means, *i_K );
  free2D( covmat, *i_p );
  free2D( weights, *i_K );
}


