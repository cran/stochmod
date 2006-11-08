// matrix.h - headers for matrix manipulation routines
//
// by Artem Sokolov
//////////////////////////////////////////////////////

#ifndef MATRIX_H__INCLUDED_C
#define MATRIX_H__INCLUDED_C

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <R_ext/Lapack.h>

// Memory management routines
double* make1D( int n );
double** make2D( int m, int n );
double*** make3D( int d1, int d2, int d3 );
void free2D( double** m, int n );
void free3D( double*** a, int d1, int d2 );

// stdout interaction
void disp_vec( int n, double* v );
void disp_mat( int m, int n, double** a );
void disp_3D( int d1, int d2, int d3, double*** a );

// Copy routines
void cp_const_vec( int n, double c, double* dest );
void cp_vec_vec( int n, double* src, double* dest );
void cp_mat_mat( int m, int n, double** src, double** dest );
void cp_3D_3D( int d1, int d2, int d3, double*** src, double*** dest );
void cp_vec_mat( int m, int n, double* src, double** dest );
void cp_vec_3D( int d1, int d2, int d3, double* src, double*** dest );
void cp_mat_vec( int m, int n, double** src, double* dest );
void cp_3D_vec( int d1, int d2, int d3, double*** src, double* dest );
void cp_diag_mat( int m, int n, int c, double** dest );

// Linear algebra routines
double vec_sum( int n, double* v );
void mult_const_vec( int n, double c, double* v );
void mult_const_mat( int m, int n, double c, double** mat );
void mult_vec_vec( int n, double* v1, double* v2, double* res );
void mult_mat_vec( int m, int p, double** m1, double* v, double* res );
void mult_mat_mat( int m, int p, int n, double** m1, double** m2, double** res );
void mult_mat_mat_ntr( int m, int n, double** m1, double** m2, double** res );
double dot_vec_vec( int n, const double* v1, const double* v2 );
void outer_vec_vec( int m, int n, double* v1, double* v2, double** res );
void add_vec_vec( int n, double c1, double c2, double* v1, double* v2, double* res );
void add_mat_mat( int m, int n, double c1, double c2, double** m1, double** m2, double** res );
void mat_transpose( int m, int n, double** A, double** B );
int covmat_inverse( int n, double** A, double** res );
int covmat_inverse_inplace( int n, double** A );
int covmat_eigenvals( int n, double** A, double* evl );
double covmat_det( int n, double** A );
double covmat_logdet( int n, double** A );
void cov_mat( int N, int p, double** A, double** res, double* means );

#endif
