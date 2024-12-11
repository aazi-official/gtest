#ifndef UTILS_H
#define UTILS_H

#include <math.h>
#include <iostream>
#include <vector>

typedef std::vector<double>               RealVector;
typedef std::vector<std::vector<double> > RealMatrix;

double Det(RealMatrix p, int n);                     // Determinant of a matrix
double Det(double *p, int n); 
double Creat_M(RealMatrix p, int m, int n, int k);  
void matrixInversion(RealMatrix &invMat, RealMatrix Mat, int nRow);
void matrixInversion(RealMatrix &invMat, RealMatrix Mat, int nRow, double det);


double vecMatVec(RealVector&  v1, RealMatrix& mat, RealVector& v2); // Compute v1 * mat *v2, the result is a scalar
RealVector matVec(RealMatrix& mat, RealVector& vec);
void vectorAdd(RealVector& v1, RealVector& v2, RealVector& v3);
#endif