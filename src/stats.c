// stats.c - statistical routines
//
// by Artem Sokolov

#include "stats.h"
#include "matrix.h"
#include <math.h>

// Multivariate Normal probability density function
//
// x - [N x p] matrix of N samples in p dimensions
// mu - [p x 1] vector representing the mean
// Sigma - [p x p] covariance matrix
//
// res - [N x 1] will contain results of applying the pdf to each sample
//
// Retruns:
//	0 - Everything is OK
//	1 - Covariance matrix is singular
int dmvnorm( int N, int p, double** x, double* mu, double** Sigma, double* res )
{
  int i;
  double coefficient;
  double exponent;
  double Sigma_det;
  double* diff = make1D( p );
  double* temp = make1D( p );
  double** Sigma_inv = make2D( p, p );

  // Compute sigma inverse and determinant
  cp_mat_mat( p, p, Sigma, Sigma_inv );
  Sigma_det = covmat_det( p, Sigma );

  // Assess numerical stability
  if( isinf( Sigma_det ) || Sigma_det == 0 )
    {
      free( diff );
      free( temp );
      free2D( Sigma_inv, p );
      return 1;
    }

  covmat_inverse_inplace( p, Sigma_inv );

  // Compute the coefficient
  coefficient = pow( (2.0 * M_PI), (p / 2.0) ) * sqrt( Sigma_det );
  coefficient = 1.0 / coefficient;

  // Assess numerical stability
  /*  if( isnan(coefficient) || isinf(coefficient) )
    {
      res[0] = HUGE_VAL;
      free( diff );
      free( temp );
      free2D( Sigma_inv, p );
      return;
      }*/

  // Perform computation one sample at a time
  for( i = 0; i < N; i++ )
    {
      add_vec_vec( p, 1.0, -1.0, x[i], mu, diff );
      mult_mat_mat( 1, p, p, &diff, Sigma_inv, &temp );
      mult_mat_vec( 1, p, &temp, diff, &exponent );
      res[i] = coefficient * exp( -exponent / 2.0 );
    }

  free( diff );
  free( temp );
  free2D( Sigma_inv, p );
  return 0;
}

