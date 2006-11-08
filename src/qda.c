// qda.c - Quadratic Discriminant Analysis routines
//
// by Artem Sokolov

#include <stdio.h>

#include "matrix.h"

// Trains a QDA classifier
//
// Inputs:
// 	i_x - N x p data matrix of N samples in p dimensions
// 	i_y - N x 1 vector of labels
//	i_N - number of samples
//	i_p - data dimensionality
//	i_K - number of classes
//	i_debug_output - level of debug output
//
// Output:
//	o_priors - K x 1 vector of Bayesian priors, computed as fraction of points from each class
//	o_means - K x p matrix of means approximated from the data
//	o_covmats - K x p x p array of covariance matrices approximated from the data
//	o_icovmats - K x p x p array of inverse covariance matrices
//	o_bias - K x 1 vector of bias terms for discriminant function computations
//
// Memory allocation requirements for outputs:
//	o_priors must be of length K
//	o_means must be of length K*p
//	o_covmat must be of length K*p*p
//	o_icovmat must be of length K*p*p
//	o_bias must be of lenght K
//	o_status must be of length 1
//	o_info must be of length 1
//
// Meaning of o_status:
//	0 - Everything is OK
//	1 - The covariance matrix for class *o_info is singular
void QDAtrain( 
	      double* i_x, int* i_y, int* i_N, int* i_p, int* i_K, double* i_cov_reg,
	      int* i_debug_output,
	      double* o_priors, double* o_means, double* o_covmats, double* o_icovmats,
	      double* o_bias, int* o_status, int* o_info
	      )
{
  if( *i_debug_output > 0 )
    printf( "Currently running C version of the QDA.train code\n" );

  // Create local storage
  int i, j, t, k;
  int* Nk = (int*)malloc( *i_K * sizeof(int) );
  double** x = make2D( *i_N, *i_p );
  double** means = make2D( *i_K, *i_p );
  double*** covmats = make3D( *i_K, *i_p, *i_p );
  double*** icovmats = make3D( *i_K, *i_p, *i_p );

  // Copy and reshape the inputs
  cp_vec_mat( *i_N, *i_p, i_x, x );

  // Init the counts and means to 0
  for( k = 0; k < *i_K; k++ )
    {
      Nk[k] = 0;
      for( i = 0; i < *i_p; i++ )
	means[k][i] = 0.0;
    }

  // Init the covariance matrices to 0
  for( k = 0; k < *i_K; k++ )
    for( i = 0; i < *i_p; i++ )
      for( j = 0; j < *i_p; j++ )
	covmats[k][i][j] = 0.0;

  // Traverse one sample at a time to compute mean contributions
  for( t = 0; t < *i_N; t++ )
    {
      k = i_y[t] - 1;		// account for 0-based index

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

  // Subtract the means from the samples and compute the covariance matrices
  for( t = 0; t < *i_N; t++ )
    {
      k = i_y[t] - 1;		// account for 0-based index
      add_vec_vec( *i_p, 1.0, -1.0, x[t], means[k], x[t] );

      for( i = 0; i < *i_p; i++ )
	for( j = 0; j < *i_p; j++ )
	  covmats[k][i][j] += x[t][i] * x[t][j];
    }

  for( k = 0; k < *i_K; k++ )
    for( i = 0; i < *i_p; i++ )
      for( j = 0; j < *i_p; j++ )
	covmats[k][i][j] /= (Nk[k] - 1);

  // Regularize the covariance matrices
  if( *i_cov_reg > 0 )
    {
      for( k = 0; k < *i_K; k++ )
	for( i = 0; i < *i_p; i++ )
	  for( j = 0; j < *i_p; j++ )
	    {
	      covmats[k][i][j] *= (1 - *i_cov_reg);
	      if( i == j )
		covmats[k][i][j] += *i_cov_reg;
	    }
    }
  
  // Compute the bias terms
  for( k = 0; k < *i_K; k++ )
    o_bias[k] = log( o_priors[k] ) - 0.5 * covmat_logdet( *i_p, covmats[k] );

  // Compute the covariance matrix inverses
  for( k = 0; k < *i_K; k++ )
    {
      *o_status = covmat_inverse( *i_p, covmats[k], icovmats[k] );
      if( *o_status != 0 )
	{
	  free2D( x, *i_N );
	  free2D( means, *i_K );
	  free3D( covmats, *i_K, *i_p );
	  free3D( icovmats, *i_K, *i_p );
	  *o_info = k;
	  return;
	}
    }

  // Copy the results from local storage to outputs
  cp_mat_vec( *i_K, *i_p, means, o_means );
  cp_3D_vec( *i_K, *i_p, *i_p, covmats, o_covmats );
  cp_3D_vec( *i_K, *i_p, *i_p, icovmats, o_icovmats );

  // Free memory
  free2D( x, *i_N );
  free2D( means, *i_K );
  free3D( covmats, *i_K, *i_p );
  free3D( icovmats, *i_K, *i_p );
}

// Applies previously trained QDA classifier to new data matrix x
//
// Inputs:
// 	i_x - N x p data matrix of N samples in p dimensions
//	i_means - K x p matrix of class centroids
// 	i_icovmats - K x p x p array of inverse covariance matrices
//	i_bias - vector of length containing the bias term for each discriminant function
//	i_N - number of samples in x
//	i_p - dimensionality of x
//	i_K - number of classes
//
// Outputs:
//	o_labels - labels for all samples in i_x
//
// Memory allocation requirements for the outputs:
//	o_labels must be of length N
void QDAtest(
	     double* i_x, double* i_means, double* i_icovmats, double* i_bias,
	     int* i_N, int* i_p, int* i_K, int* o_labels
	     )
{
  // Create local storage
  int t, k, i, j;
  double temp;
  double* x_c = make1D( *i_p );
  double** x = make2D( *i_N, *i_p );
  double** means = make2D( *i_K, *i_p );
  double*** icovmats = make3D( *i_K, *i_p, *i_p );
  double** delta = make2D( *i_N, *i_K );

  // Copy and reshape the inputs to local storage
  cp_vec_mat( *i_N, *i_p, i_x, x );
  cp_vec_mat( *i_K, *i_p, i_means, means );
  cp_vec_3D( *i_K, *i_p, *i_p, i_icovmats, icovmats );

  // Compute the delta function values
  for( t = 0; t < *i_N; t++ )
    for( k = 0; k < *i_K; k++ )
      {
	add_vec_vec( *i_p, 1.0, -1.0, x[t], means[k], x_c );
	delta[t][k] = 0.0;
	for( i = 0; i < *i_p; i++ )
	  for( j = 0; j < *i_p; j++ )
	    delta[t][k] += icovmats[k][i][j] * x_c[i] * x_c[j];
	delta[t][k] *= -0.5;
	delta[t][k] += i_bias[k];
      }

  // Pick the highest discriminant function value for each point
  for( t = 0; t < *i_N; t++ )
    {
      o_labels[t] = 1;
      temp = delta[t][0];
      for( k = 1; k < *i_K; k++ )
	{
	  if( delta[t][k] > temp )
	    {
	      o_labels[t] = k + 1;
	      temp = delta[t][k];
	    }
	}
    }

  // Free up memory
  free( x_c );
  free2D( x, *i_N );
  free2D( means, *i_K );
  free3D( icovmats, *i_K, *i_p );
  free2D( delta, *i_N );
}
