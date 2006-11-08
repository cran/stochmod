// hmm.c - HMM routines
//
// by Artem Sokolov
///////////////////////

#include "matrix.h"
#include "stats.h"

#include "hmm.h"

// Computes observation probabilities given HMM parameters
//
// Inputs;
//	x - [N x p] matrix of data samples
//	pi - [K x 1] vector of priors
//	A - [K x K] transition samples
//	mu - [K x p] matrix of state means
//	sigma - [K x p x p] array of state covariance matrices
//	N - number of samples
//	K - number of states
//	p - observation dimensionality
//
// Output:
//	B - [N x K] matrix of observation probabilities for each state
//	info - a single integer used to communicate additional information
//
// Returns:
//	0 - Everything is OK
//	1 - Covariance matrix for state *info is singular
int HMM_obs( double** x, double* pi, double** A, double** mu, double*** sigma,
	     int N, int K, int p, double** B, int* info )
{
  int status;
  int i, j;
  double* res = make1D( N );

  // Use Normal pdf for each state
  for( i = 0; i < K; i++ )
    {
      status = dmvnorm( N, p, x, mu[i], sigma[i], res );

      // Check for numerical stability
      if( status > 0 )
	{
	  *info = i+1;
	  free( res );
	  return status;
	}

      for( j = 0; j < N; j++ ) B[j][i] = res[j];
    }

  free( res );
  return 0;
}



// Forward-Backward procedure
//
// Inputs:
//	x - [N x p] matrix of data samples
//	pi - [K x 1] vector of priors
//	A - [K x K] transition samples
//	mu - [K x p] matrix of state means
//	sigma - [K x p x p] array of state covariance matrices
//	N - number of samples
//	K - number of states
//	p - observation dimensionality
//
// Outputs:
//	alpha - [N x K[ matrix of forward variable values
//	beta - [N x K] matrix of backward variable values
//	c - [N x 1] vector of scaling coefficients
//	gamma - [N x K] matrix of transitions from state i: gamma[t,i] = P( Qt = i | O,hmmp )
//	xi - [N-1 x K x K] matrix of transitions from state i to state j: xi[t,i,j] = P( Qt = i, Qt+1 = j | O,hmmp )
//
// Returns:
//	0 - Everything is OK
//	1 - Covariance matrix for one of the states is singular
//	2 - Some samples were not accounted for by any state
int HMM_fwbk( double** x, double* pi, double** A, double** mu, double*** sigma, int N, int K, int p,
	      double** alpha, double** beta, double* c, double** gamma, double*** xi )
{
  // Create local storage
  int t, i, j, k;
  int status, info;

  double val;
  double* t1 = make1D( K );
  double** B = make2D( N, K );

  // Compute observation probabilities
  status = HMM_obs( x, pi, A, mu, sigma, N, K, p, B, &info );
  if( status == 1 )
    {
      free2D( B, N ); free( t1 );
      return status;
    }

  // Check for sample assignment to states
  for( i = 0; i < N; i++ )
    {
      if( vec_sum( K, B[i] ) == 0 )
	{
	  free2D( B, N ); free( t1 );
	  return 2;
	}
    }

  // Compute alpha values for the first observation
  mult_vec_vec( K, pi, B[0], alpha[0] );
  c[0] = 1.0 / vec_sum( K, alpha[0] );
  mult_const_vec( K, c[0], alpha[0] );

  // Recursively compute remaining alpha values
  for( i = 1; i < N; i++ )
    {
      mult_mat_mat( 1, K, K, alpha+(i-1), A, &t1 );
      mult_vec_vec( K, t1, B[i], alpha[i] );
      c[i] = 1.0 / vec_sum( K, alpha[i] );
      mult_const_vec( K, c[i], alpha[i] );
    }

  // Compute beta values for the last observation
  cp_const_vec( K, c[N-1], beta[N-1] );		// really c[N-1] * 1
  
  // Recursivey compute remaining beta values
  for( i = (N-2); i >= 0; i-- )
    {
      mult_vec_vec( K, beta[i+1], B[i+1], t1 );
      mult_mat_vec( K, K, A, t1, beta[i] );
      mult_const_vec( K, c[i], beta[i] );
    }

  // Compute and normalize gamma values
  mult_mat_mat_ntr( N, K, alpha, beta, gamma );
  for( i = 0; i < N; i++ )
    {
      val = 1.0 / vec_sum( K, gamma[i] );
      mult_const_vec( K, val, gamma[i] );
    }

  // Compute xi values
  for( t = 0; t < (N-1); t++ )
    for( i = 0; i < K; i++ )
      {
	mult_vec_vec( K, B[t+1], beta[t+1], t1 );
	mult_vec_vec( K, A[i], t1, xi[t][i] );
	mult_const_vec( K, alpha[t][i], xi[t][i] );
      }
  // Normally would divide by xi[t], but scaling takes care of it

  // Free memory
  free2D( B, N );
  free( t1 );
  return 0;
}


// Computes log likelihood of an observation sequence given HMM parametesrs
//
// Inputs:
//	x - [N x p] matrix of the observation samples
//	pi - [K x 1] vector of priors
//	A - [K x K] transition matrix
//	mu - [K x p] matrix of state means
//	sigma - [K x p x p] matrix of state covariance matrices
//	N - number of samples
//	K - number of states
//	p - observation dimensionality
//
// Outputs:
//	return log likelihood
//
// Returns:
//	0 - Everything is OK
//	1 - Covariance matrix for is singular for one of the states
int HMM_LL( double** x, double* pi, double** A, double** mu, double*** sigma, int N, int K, int p, double* LL )
{
  // Need to compute forward variables and the corresponding scaling coefficients

  // Create local storage
  int t, i, j, k;
  int status, info;
  
  double* c = make1D( N );
  double* t1 = make1D( K );
  double** B = make2D( N, K );
  double** alpha = make2D( N, K );

  // Compute observation probabilities
  status = HMM_obs( x, pi, A, mu, sigma, N, K, p, B, &info );

  if( status != 0 )
    {
      free2D( B, N ); free( t1 );
      free( c ); free2D( alpha, N );
      return status;
    }

  // Compute alpha values for the first observation
  mult_vec_vec( K, pi, B[0], alpha[0] );
  c[0] = 1.0 / vec_sum( K, alpha[0] );
  mult_const_vec( K, c[0], alpha[0] );

  // Recursively compute remaining alpha values
  for( i = 1; i < N; i++ )
    {
      mult_mat_mat( 1, K, K, alpha+(i-1), A, &t1 );
      mult_vec_vec( K, t1, B[i], alpha[i] );
      c[i] = 1.0 / vec_sum( K, alpha[i] );
      mult_const_vec( K, c[i], alpha[i] );
    }

  *LL = 0.0;
  for( i = 0; i < N; i++ )
    *LL -= log(c[i]);

  free( c );
  free( t1 );
  free2D( B, N );
  free2D( alpha, N );

  return 0;
}



// Learns maximum likelihood HMM given training and validation data
//
// Inputs:
//	x - [nSeqs x N[i] x p] array of training observation sequences
//	v - [mSeqs x M[i] x p] array of validation observation sequences
//	nSeqs - Number of training sequences
//	mSeqs - Number of validation sequences
//	N - [nSeqs x 1] array of observation sequence lengths for training data
//	M - [mSeqs x 1] array of observation sequence lengths for validation data
//	K - desired number of states
//	p - observation dimensionality
//
// Control variables:
//	tol - relative tolerance for determining convergence
//	LLstop - iterations stop when training log-likelihood exceeds this value
//	max_iter - maximum number of iterations of EM algorithm (use a negative value
//			to denote infinity)
// 
// Amount of debug output:
//	0 - none
//	1 - basic
//	2 - verbose
//
// Outputs:
//	pi - [K x 1] vector of priors
//	A - [K x K] transition matrix
//	mu - [K x p] matrix of state means
//	sigma - [K x p x p] matrix of state covariance matrices
//
// Returns:
//	0 - Everything is OK
//	1 - Covariance matrix is singular
//	2 - Some samples were not accounted for by any state
//	3 - Computational instability
int HMM_learn( double*** x, double*** v, int nSeqs, int mSeqs, int* N, int* M, int K, int p,
	       double* pi, double** A, double** mu, double*** sigma, 
	       double tol, double LLstop, int min_iter, int max_iter, int debug_output )
{
  if( debug_output > 0 )
    printf( "Currently running C version of the HMM.learn code\n" );

  // Create local storage
  int i, j, k, l, t;
  double v1, v2, v3;

  // Concatenate all training data together, some routines require it
  int N_total = 0;
  int x_offset = 0;
  for( i = 0; i < nSeqs; i++ ) N_total += N[i];

  double** xCat = make2D( N_total, p );
  for( i = 0; i < nSeqs; i++ )
    {
      cp_mat_mat( N[i], p, x[i], xCat+x_offset );
      x_offset += N[i];
    }

  int M_total = 0;
  double** vCat = NULL;
  double* pi_prev = NULL;
  double** A_prev = NULL;
  double** mu_prev = NULL;
  double*** sigma_prev = NULL;
  if( v != NULL )
    {
      pi_prev = make1D( K );
      A_prev = make2D( K, K );
      mu_prev = make2D( K, p );
      sigma_prev = make3D( K, p, p );

      for( i = 0; i < mSeqs; i++ ) M_total += M[i];
      vCat = make2D( M_total, p );
      int v_offset = 0;
      for( i = 0; i < mSeqs; i++ )
	{
	  cp_mat_mat( M[i], p, v[i], vCat+v_offset );
	  v_offset += M[i];
	}
    }

  double* c = make1D( N_total );
  double** alpha = make2D( N_total, K );
  double** beta = make2D( N_total, K );
  double** gamma = make2D( N_total, K );
  double** gamma_t = make2D( K, N_total );
  double*** xi = make3D( N_total - nSeqs, K, K );

  double* gamma1 = make1D( K );
  double* gammasum = make1D( K );
  double** xisum = make2D( K, K );

  double* x_c = make1D( p );
  double** cov_mat = make2D( p, p );

  // Run EM
  int status = 0;
  double LL, LLbase, LLprev;
  double LLval, LLval_prev;
  double LLval_temp;
  int iter = 0;
  while( 1 )
    {
      // Validation log-likelihood stopping criterion
      if( v != NULL )
	{
	  LLval = 0.0;
	  for( i = 0; i < mSeqs; i++ )
	    {
	      status = HMM_LL( v[i], pi, A, mu, sigma, M[i], K, p, &LLval_temp );
	      if( status == 0 )
		LLval += LLval_temp;
	      else
		{
		  free2D( xCat, N_total );
		  free( pi_prev );
		  free2D( A_prev, K );
		  free2D( mu_prev, K );
		  free3D( sigma_prev, K, p );
		  free2D( vCat, M_total );
		  free( c );
		  free2D( alpha, N_total );
		  free2D( beta, N_total );
		  free2D( gamma, N_total );
		  free2D( gamma_t, K );
		  free3D( xi, N_total - nSeqs, K );
		  free( gamma1 );
		  free( gammasum );
		  free2D( xisum, K );
		  free( x_c );
		  free2D( cov_mat, p );

		  return status;
		}
	    }
	  if( debug_output > 2 )
	    printf( "Validation log-likelihood: %f\n", LLval );

	  if( iter <= min_iter || LLval >= LLval_prev )
	    {
	      LLval_prev = LLval;
	      cp_vec_vec( K, pi, pi_prev );
	      cp_mat_mat( K, K, A, A_prev );
	      cp_mat_mat( K, p, mu, mu_prev );
	      cp_3D_3D( K, p, p, sigma, sigma_prev );
	    }
	  else
	    {
	      if( debug_output > 0 )
		{
		  printf( "Final log-likelihood of training data after %d iterations is %f\n", iter-1, LLprev );
		  printf( "Validation log-likelihood: %f\n", LLval_prev );
		}
	      cp_vec_vec( K, pi_prev, pi );
	      cp_mat_mat( K, K, A_prev, A );
	      cp_mat_mat( K, p, mu_prev, mu );
	      cp_3D_3D( K, p, p, sigma_prev, sigma );
		  
	      status = 0;
	      break;
	    }
	}

      // Init gamma1, gammasum, and xisum
      for( j = 0; j < K; j++ )
	{
	  gamma1[j] = 0.0;
	  gammasum[j] = 0.0;
	  for( k = 0; k < K; k++ )
	    xisum[j][k] = 0.0;
	}

      // Run each observation sequence through the forward-backward algorithm
      //
      // E Step
      x_offset = 0;
      for( i = 0; i < nSeqs; i++ )
	{
	  status = HMM_fwbk( x[i], pi, A, mu, sigma, N[i], K, p, 
			     alpha+x_offset, beta+x_offset, c+x_offset, gamma+x_offset, xi+x_offset-i );

	  // Assess numerical stability
	  if( status != 0 ) break;

	  // Update gamma1, gammasum, and xisum
	  for( j = 0; j < K; j++ )
	    {
	      gamma1[j] += gamma[x_offset][j];
	      for( k = 0; k < N[i]; k++ )
		gammasum[j] += (gamma+x_offset)[k][j];
	      for( k = 0; k < K; k++ ) for( l = 0; l < (N[i]-1); l++ )
		  xisum[j][k] += (xi+x_offset-i)[l][j][k];
	    }

	  x_offset += N[i];
	}

      // Compute log-likelihood
      LL = 0.0;
      for( i = 0; i < N_total; i++ )
	LL -= log( c[i] );

      if( debug_output > 2 || (debug_output > 1  && (iter % 10) == 0) )
	printf( "Log likelihood at iteration %d is %f\n", iter, LL );

      if( isnan( LL ) || isinf( LL ) )
	{
	  status = 3;
	  break;
	}

      // Verify stopping criteria
      if( max_iter > 0 && iter >= max_iter ) break;
      if( LL > LLstop ) break;
      if( isinf( LL ) ) break;

      if( iter <= (min_iter-1) ) LLbase = LL;
      else if( (LL-LLbase) < (1+tol)*(LLprev-LLbase) ) break;

      // Update parameters
      //
      // M Step

      // Update means 
      mat_transpose( N_total, K, gamma, gamma_t );
      mult_mat_mat( K, N_total, p, gamma_t, xCat, mu );
      for( k = 0; k < K; k++ )
	for( j = 0; j < p; j++ )
	  mu[k][j] /= gammasum[k];

      // Update transition matrix
      cp_mat_mat( K, K, xisum, A );
      for( i = 0; i < K; i++ )
	mult_const_vec( K, 1.0 / vec_sum( K, A[i] ), A[i] );

      // Update priors
      cp_vec_vec( K, gamma1, pi );
      mult_const_vec( K, 1.0 / vec_sum( K, pi ), pi );

      // Update covariance matrices
      for( l = 0; l < K; l++ )
	{
	  for( i = 0; i < p; i++ ) for( j = 0; j < p; j++ ) sigma[l][i][j] = 0.0;
	  for( t = 0; t < N_total; t++ )
	    {
	      add_vec_vec( p, 1.0, -1.0, xCat[t], mu[l], x_c );
	      outer_vec_vec( p, p, x_c, x_c, cov_mat );
	      add_mat_mat( p, p, 1.0, gamma[t][l], sigma[l], cov_mat, sigma[l] );
	    }
	  mult_const_mat( p, p, 1.0 / gammasum[l], sigma[l] );
	}

      LLprev = LL;
      iter++;
    }

  if( status == 0 && debug_output > 0 )
    printf( "Final log-likelihood after %d iterations is %f\n", iter, LL );

  // Free memory
  free2D( xCat, N_total );
  if( v != NULL )
    {
      free( pi_prev );
      free2D( A_prev, K );
      free2D( mu_prev, K );
      free3D( sigma_prev, K, p );
      free2D( vCat, M_total );
    }
  free( c );
  free2D( alpha, N_total );
  free2D( beta, N_total );
  free2D( gamma, N_total );
  free2D( gamma_t, K );
  free3D( xi, N_total - nSeqs, K );
  free( gamma1 );
  free( gammasum );
  free2D( xisum, K );
  free( x_c );
  free2D( cov_mat, p );

  return 0;
}

