// matrix.c - matrix manipulation routines
//
// by Artem Sokolov
//////////////////////////////////////////

#include "matrix.h"

////////////////////////////////////////////////////////
// Memory management routines
////////////////////////////////////////////////////////

// Allocates memory for a 1-D vector of specified length
double* make1D( int n )
{
  double* res = (double*)malloc( n*sizeof(double) );
  return res;
}

// Allocates memory for a 2-D matrix of specified dimensions
//
// m - number of rows
// n - number of columns
double** make2D( int m, int n )
{
  int i;
  double** res = (double**)malloc( m*sizeof(double*) );
  for( i = 0; i < m; i++ )
    res[i] = (double*)malloc( n*sizeof(double) );
  return res;
}

// Allocates memory for a 3-D array of specified dimensions
//
// d1 - length along the first dimension
// d2 - length along the second dimension
// d3 - length along the third dimension
double*** make3D( int d1, int d2, int d3 )
{
  int i, j;
  double*** res = (double***)malloc( d1*sizeof(double**) );
  for( i = 0; i < d1; i++ )
    {
      res[i] = (double**)malloc( d2*sizeof(double*) );
      for( j = 0; j < d2; j++ )
	res[i][j] = (double*)malloc( d3*sizeof(double) );
    }

  return res;
}

// Frees up memory allocated for a 2-D matrix
//
// m - the matrix
// n - number of rows
void free2D( double** m, int n )
{
  int i;
  for( i = 0; i < n; i++ )
    free( m[i] );
  free( m );
}

// Frees up memory allocated for a 3-D array
//
// a - the array
// d1 - length along the first dimension
// d2 - length along the second dimension
void free3D( double*** a, int d1, int d2 )
{
  int i, j;
  for( i = 0; i < d1; i++ )
    {
      for( j = 0; j < d2; j++ )
	free( a[i][j] );
      free( a[i] );
    }
  free( a );
}

////////////////////////////////////////////////////////
// stdout interaction
////////////////////////////////////////////////////////

// Display a vector to stdout
// 
// n - size of the vector
// v - the vector itself
void disp_vec( int n, double* v )
{
  int i;

  for( i = 0; i < n; i++ )
    printf( "%f ", v[i] );
  printf( "\n" );
}

// Displays a matrix to stdout
//
// m - number of rows
// n - number of columns
// a - the matrix itself
void disp_mat( int m, int n, double** a )
{
  int i,j;

  for( i = 0; i < m; i++ )
    {
      for( j = 0; j < n; j++ )
	printf( "%f ", a[i][j] );
      printf( "\n" );
    }
}

// Displays a 3-D array to stdout
//
// d1 - length along the first dimension
// d2 - length along the second dimension
// d3 - length along the third dimension
// a - the array itself
void disp_3D( int d1, int d2, int d3, double*** a )
{
  int i;
  for( i = 0; i < d1; i++ )
    {
      printf( "[%d,,] =\n", i );
      disp_mat( d2, d3, a[i] );
    }
}

////////////////////////////////////////////////////////
// Copy routines
////////////////////////////////////////////////////////

// Replicates a constant across all entries of a vector
// 
// n - vector length
// c - the constant to replicate
// dest - the vector itself
void cp_const_vec( int n, double c, double* dest )
{
  int i;
  for( i = 0; i < n; i++ )
    dest[i] = c;
}

// Copies content of one vector into another
//
// n - size of the vector
// src - source
// dest - destination
void cp_vec_vec( int n, double* src, double* dest )
{
  int i;
  for( i = 0; i < n; i++ )
    dest[i] = src[i];
}

// Copies content of one matrix into another
//
// m - number of rows
// n - number of columns
// src - source
// dest - destination
void cp_mat_mat( int m, int n, double** src, double** dest )
{
  int i,j;
  for( i = 0; i < m; i++ )
    for( j = 0; j < n; j++ )
      dest[i][j] = src[i][j];
}

// Copies content of one 3-D array into another
//
// d1 - length along the first dimension
// d2 - length along the second dimension
// d3 - length along the third dimension
// src - source
// dest - destination
void cp_3D_3D( int d1, int d2, int d3, double*** src, double*** dest )
{
  int i,j,k;
  for( i = 0; i < d1; i++ )
    for( j = 0; j < d2; j++ )
      for( k = 0; k < d3; k++ )
	dest[i][j][k] = src[i][j][k];
}

// Copies and reshape content of a 1-D vector to a 2-D matrix
//   performed column-major to be compatible with R
//
// m - number of rows in dest
// n - number of columns in dest
// src - source
// dest - destination 
void cp_vec_mat( int m, int n, double* src, double** dest )
{
  int i,j;
  for( i = 0; i < m; i++ )
    for( j = 0; j < n; j++ )
      dest[i][j] = src[ j*m+i ];
}

// Copies and reshapes content of a 1-D vector to a 3-D array
//   performed column-major to be compatible with R
//
// d1 - length along the first dimension of dest
// d2 - length along the second dimension of dest
// d3 - length along the third dimension of dest
// src - source
// dest - destination
void cp_vec_3D( int d1, int d2, int d3, double* src, double*** dest )
{
  int i, j, k;
  for( i = 0; i < d1; i++ )
    for( j = 0; j < d2; j++ )
      for( k = 0; k < d3; k++ )
	dest[i][j][k] = src[ k*(d1*d2) + j*d1 + i ];
}

// Copies and reshapes content of a 2-D matrix to a 1-D vector
//   performed column-major to be compatible with R
//
// m - number of rows in src
// n - number of rows in dest
// src - source
// dest - destination
void cp_mat_vec( int m, int n, double** src, double* dest )
{
  int i,j;
  for( i = 0; i < m; i++ )
    for( j = 0; j < n; j++ )
      dest[j*m+i] = src[i][j];
}

// Copies and reshapes content of a 3-D array to a 1-D vector
//   performed column-major to be compatible with R
//
// d1 - length along the first dimension of src
// d2 - length along the second dimension of src
// d3 - length along the third dimension of src
// src - source
// dest - destination
void cp_3D_vec( int d1, int d2, int d3, double*** src, double* dest )
{
  int i,j,k;
  for( i = 0; i < d1; i++ )
    for( j = 0; j < d2; j++ )
      for( k = 0; k < d3; k++ )
	dest[ k*(d1*d2) + j*d1 + i ] = src[i][j][k];
}

// Turns matrix into a diagonal one, by putting c on its diagonal
//   and 0s everywhere else
//
// m - number of rows
// n - number of columns
// c - the diagonal term
// dest - result will be stored here
void cp_diag_mat( int m, int n, int c, double** dest )
{
  int i,j;
  for( i = 0; i < m; i++ )
    for( j = 0; j < m; j++ )
      {
	if( i == j ) dest[i][j] = c;
	else dest[i][j] = 0.0;
      }
}

////////////////////////////////////////////////////////
// Linear algebra routines
////////////////////////////////////////////////////////

// Computes vector sum
//
// n - vector length
// v - the vector itself
//
// output: sum of the entries of v
double vec_sum( int n, double* v )
{
  int i;
  double res = 0.0;
  for( i = 0; i < n; i++ )
    res += v[i];
  return res;
}

// Multiplies vector by a constant factor
//
// n - vector length
// c - scaling coefficients
// v - the vector itself
//
// output: v will contain (c*v)
void mult_const_vec( int n, double c, double* v )
{
  int i;
  for( i = 0; i < n; i++ )
    v[i] *= c;
}

// Multiplies matrix by a constant factor
//
// m - number of rows
// n - numebr of columns
// c - scaling coefficient
// mat - the matrix itself
//
// result will be contained in mat
void mult_const_mat( int m, int n, double c, double** mat )
{
  int i,j;
  for( i = 0; i < m; i++ )
    for( j = 0; j < n; j++ )
      mat[i][j] *= c;
}

// Multiplies two vectors component-wise
//
// n - number of elements
// v1 - first vector
// v2 - second vector
// res - will hold the result
void mult_vec_vec( int n, double* v1, double* v2, double* res )
{
  int i;
  for( i = 0; i < n; i++ )
    res[i] = v1[i] * v2[i];
}

// Matrix-Matrix multiplication
// All memory must have been already allocated
//
// Inputs:
// 	m1 - matrix 1 [m x p]
// 	m2 - matrix 2 [p x n]
//
// Outputs:
// 	res - m1 * m2 [m x n]
void mult_mat_mat( int m, int p, int n, double** m1, double** m2, double** res )
{
  int i, j, k;

  for( i = 0; i < m; i++ )
    for( j = 0; j < n; j++ )
      {
	res[i][j] = 0.0;
	for( k = 0; k < p; k++ )
	    res[i][j] += m1[i][k] * m2[k][j];
      }
}

// Matrix-Matrix multiplication (entry-by-entry)
//	m - number of rows
//	n - number of columns
//	m1 - matrix 1
//	m2 - matrix 2
//  	res - will contain the result
void mult_mat_mat_ntr( int m, int n, double** m1, double** m2, double** res )
{
  int i, j;
  for( i = 0; i < m; i++ )
    for( j = 0; j < n; j++ )
      res[i][j] = m1[i][j] * m2[i][j];
}

// Matrix-Vector multiplication
// All memory must have been already allocated
//
// Inputs:
// 	m1 - matrix [m x p]
//	v - vector [p x 1]
//
// Outputs:
//	res = m1 * v [m x 1]
void mult_mat_vec( int m, int p, double** m1, double* v, double* res )
{
  int i, j;
  
  for( i = 0; i < m; i++ )
    {
      res[i] = 0.0;
      for( j = 0; j < p; j++ )
	res[i] += m1[i][j] * v[j];
    }
}

// Computes the dot product of two vectors
//
// Inputs:
//	n - vector length
//	v1 - vector 1
//	v2 - vector 2
//
// Output:
//	res = t(v1) * (v2)
double dot_vec_vec( int n, const double* v1, const double* v2 )
{
  int i;
  double res = 0.0;
  
  for( i = 0; i < n; i++ )
    res += v1[i] * v2[i];

  return res;
}


// Computes the outer product of two vectors
//
// Inputs:
//	m - number of rows in the result
//	n - number of columns in the result
//	v1 - vector 1 (of length m)
//	v2 - vector 2 (of length n )
//
// Output:
//	res = v1 * t(v2)
void outer_vec_vec( int m, int n, double* v1, double* v2, double** res )
{
  int i,j;

  for( i = 0; i < m; i++ )
    for( j = 0; j < n; j++ )
      res[i][j] = v1[i] * v2[j];
}


// Computes a linear combination of two vectors
// All memory must have been already allocated
//
// Inputs:
//	n - vector length
//	c1 - coefficient for vector 1
//	c2 - coefficient for vector 2
//	v1 - vector 1
//	v2 - vector 2
//
// Outputs:
//	res = c1 * v1 + c2 * v2
void add_vec_vec( int n, double c1, double c2, double* v1, double* v2, double* res )
{
  int i;

  for( i = 0; i < n; i++ )
    res[i] = c1 * v1[i] + c2 * v2[i];
}


// Computes a linear combination of two matrices
// All memory must have been already allocated
//
// Inputs:
//	c1 - coefficient for matrix 1
//	c2 - coefficient for matrix 2
//	m1 - matrix 1 [m x n]
//	m2 - matrix 2 [m x n]
//
// Outputs:
//	res = c1 * m1 + c2 * m2 [m x n]
void add_mat_mat( int m, int n, double c1, double c2, double** m1, double** m2, double** res )
{
  int i, j;
  
  for( i = 0; i < m; i++ )
    for( j = 0; j < n; j++)
      res[i][j] = c1 * m1[i][j] + c2 * m2[i][j];
}


// Computes matrix transpose
// All memory must have already been allocated
//
// m - number of rows
// n - number of columns
// A - the original matrix [m x n]
// B - will hold the result [n x m]
void mat_transpose( int m, int n, double** A, double** B )
{
  int i, j;

  for( i = 0; i < m; i++ )
    for( j = 0; j < n; j++ )
      B[j][i] = A[i][j];
}

// Performs a call to LAPACK to compute a covariance matrix inverse
//
// n - matrix dimension
// A - square matrix
// res - will hold the result
//
// Returns:
// 0 - Everything is OK
// 1 - Matrix is singular
int covmat_inverse( int n, double** A, double** res )
{
  int req = -1;
  double opt = 0;
  int lwork = 0;
  double* work = NULL;
  int* ipvt = (int*)malloc( n * sizeof(int) );
  int info = 0;
  double* a = (double*)malloc( n*n*sizeof(double) );

  // Copy the content of the matrix to a local 1-D array
  cp_mat_vec( n, n, A, a );

  // Perform LU factorization
  F77_CALL(dgetrf)( &n, &n, a, &n, ipvt, &info );

  if( info != 0 )
    {
      free( ipvt );
      free( a );
      return 1;
    }

  // Request for optimal size of work space
  F77_CALL(dgetri)( &n, a, &n, ipvt, &opt, &req, &info );
  lwork = (int)(opt + 0.1);

  work = (double*)malloc( lwork * sizeof(double) );

  // Compute the inverse
  F77_CALL(dgetri)( &n, a, &n, ipvt, work, &lwork, &info );

  if( info != 0 )
    {
      free( ipvt );
      free( a );
      free( work );
      return 1;
    }

  // Copy the local array into res
  cp_vec_mat( n, n, a, res );
  
  // Free the memory and return
  free( ipvt );
  free( a );
  free( work );
  return 0;
}

// Save as above except the result is copied back into A
int covmat_inverse_inplace( int n, double** A )
{
  return covmat_inverse( n, A, A );
}


// Performs a call to LAPACK to compute eigenvalues of a covariance matrix
//
// n - matrix dimension
// A - the matrix
// evl - vector of length n that will contain the eigenvalues
//
// return 0 if succeeded, nonzero otherwise
int covmat_eigenvals( int n, double** A, double* evl )
{
  double* a = (double*)malloc( n*n*sizeof(double) );
  double* wi = (double*)malloc( n*sizeof(double) );

  int req = -1;
  double opt = 0;
  int lwork = 0;
  int ldvl = 1;
  int ldvr = 1;
  double* work = NULL;
  char jobvl = 'N';
  char jobvr = 'N';
  int info;
  
  // Copy the content of the matrix to a local 1-D array
  cp_mat_vec( n, n, A, a );

  // Request the optimal size of work space
  F77_CALL(dgeev)( &jobvl, &jobvr, &n, a, &n, evl, wi, 
		   NULL, &ldvl, NULL, &ldvr, &opt, &req, &info );
  lwork = (int)(opt+0.1);
  work = (double*)malloc( lwork * sizeof(double) );

  // Compute the eigenvalues
  F77_CALL(dgeev)( &jobvl, &jobvr, &n, a, &n, evl, wi, 
		   NULL, &ldvl, NULL, &ldvr, work, &lwork, &info );

  free( a );
  free( wi );
  free( work );
  return info;
}

// Computes determinant of a covariance matrix
//
// n - matrix dimension
// A - the matrix
//
// If failed to compute determinant, returns 0
double covmat_det( int n, double** A )
{
  int i;
  double res = 1.0;
  double* evl = (double*)malloc( n*sizeof(double) );

  int status = covmat_eigenvals( n, A, evl );
  if( status != 0 )
    {
      printf( "Failed to compute determinant\n" );
      free( evl );
      return 0.0;
    }

  // Determinant is product of the eigenvalues
  for( i = 0; i < n; i++ )
    res *= evl[i];

  free( evl );
  return res;
}

// Computes log( determinant ) of a covariance matrix
//
// n - matrix dimension
// A - the matrix
//
// If failed to compute, returns 0
double covmat_logdet( int n, double** A )
{
  int i;
  double res = 0.0;
  double* evl = (double*)malloc( n*sizeof(double) );

  int status = covmat_eigenvals( n, A, evl );
  if( status != 0 )
    {
      printf( "Failed to compute log determinant\n" );
      free( evl );
      return 0.0;
    }

  // Determinant is product of the eigenvalues
  for( i = 0; i < n; i++ )
    res += log(evl[i]);

  free( evl );
  return res;
}

// Computes the covariance matrix of given samples
//
// N - number of samples
// p - dimensionality
// A - [N x p] matrix of samples
// means - (optionally) precomputed means, set to NULL if to be computed
//
// res will contain [p x p] covariance matrix
void cov_mat( int N, int p, double** A, double** res, double* means )
{
  // Create local storage
  int i, j, k;
  double** A_cntr = make2D( N, p );

  // Mean-center the data, computing the means if necessary
  if( means == NULL )
    {
      double* i_means = make1D( p );

      for( i = 0; i < p; i++ )
	{
	  i_means[i] = 0.0;
	  for( j = 0; j < N; j++ )
	    i_means[i] += A[j][i];
	  i_means[i] /= N;
	}

      for( i = 0; i < N; i++) for( j = 0; j < p; j++ )
	A_cntr[i][j] = A[i][j] - i_means[j];

      free( i_means );
    }
  else
    {
      for( i = 0; i < N; i++) for( j = 0; j < p; j++ )
	A_cntr[i][j] = A[i][j] - means[j];
    }

  // Compute the covariance matrix
  for( i = 0; i < p; i++ )
    for( j = 0; j < p; j++ )
      {
	res[i][j] = 0.0;
	for( k = 0; k < N; k++ )
	  res[i][j] += A_cntr[k][i] * A_cntr[k][j];
	res[i][j] /= N-1;
      }

  // Free up memory
  free2D( A_cntr, N );
}
