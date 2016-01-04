#pragma once
#include <math.h>
#include <iostream>
#include "Matrix.h"

/*
Function that implemants singular vector decomposition
of a real matrix. On output M[][] is replaced by the 
Sigma matrix, matrices U and V will contain left and
rifgt singular vectors(normalized) respectively.
*/
void SVD(Matrix<float> &M, Matrix<float> &U, Matrix<float> &V);

/*
 Perform a Housholder reduction of a real 
 symmetric matrix a[][]. On output a[][] is replaced by the 
 orthogonal matrix effecting the transformation. d[] returns 
 the diagonal elements of the tri-diagonal matrix, and e[] the
 off-diagonal elements, with e[0] = 0.
 The function is modified from the version in Numerical recipe.
*/
void reduce(float **a, int n, float *d, float *e);


/*
 Determine eigenvalues and eigenvectors of a real symmetric
 tri-diagonal matrix, or a real, symmetric matrix previously
 reduced by function reduce[] to tri-diagonal form. On input,
 d[] contains the diagonal element and e[] the sub-diagonal
 of the tri-diagonal matrix. On output d[] contains the
 eigenvalues and  e[] is destroyed. If eigenvectors are
 desired z[][] on input contains the identity matrix. If
 eigenvectors of a matrix reduced by reduce() are required,
 then z[][] on input is the matrix output from reduce().
 On output, the k'th column returns the normalized eigenvector
 corresponding to d[k].
 The function is modified from the version in Numerical recipe.
*/
void findEigen(float *d, float *e, int n, float **z);


/* 
Sort eigenvalues in d[0...n-1] and appropriate eigenvectors 
in v[0...n-1][0...n-1]. Method - direct insert sort
*/
void sortEigen(float *d, float **v, int n);


