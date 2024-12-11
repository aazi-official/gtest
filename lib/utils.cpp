
#include "utils.h"


double Det(double *p, int n)
{
    int r, c, m;
    int lop = 0;
    double result = 0;
    double mid = 1;
    if (n != 1)
    {
        lop = (n == 2) ? 1 : n;            
        for (m = 0; m < lop; m++)
        {
            mid = 1;           
            for (r = 0, c = m; r < n; r++, c++)
            {
                mid = mid * (*(p+r*n+c%n));
            }
            result += mid;
        }
        for (m = 0; m < lop; m++)
        {
            mid = 1;            
            for (r = 0, c = n-1-m+n; r < n; r++, c--)
            {
                mid = mid * (*(p+r*n+c%n));
            }
            result -= mid;
        }
    }
    else 
        result = *p;
    return result;
}

double Det(RealMatrix p, int n)
{
    int r, c, m;
    int lop = 0;
    double result = 0;
    double mid = 1;
    int rowInd, colInd;
    if (n != 1)
    {
        lop = (n == 2) ? 1 : n;            
        for (m = 0; m < lop; m++)
        {
            mid = 1;           
            for (r = 0, c = m; r < n; r++, c++)
            {
                // mid = mid * (*(p+r*n+c%n));
                colInd = c%n;
                rowInd = (r*n+c%n)/n;
                mid = mid * p[rowInd][colInd];
            }
            result += mid;
        }
        for (m = 0; m < lop; m++)
        {
            mid = 1;           
            for (r = 0, c = n-1-m+n; r < n; r++, c--)
            {
                colInd = c%n;
                rowInd = (r*n+c%n)/n;
                mid = mid * p[rowInd][colInd];
                //mid = mid * p[r*n+c%n];
            }
            result -= mid;
        }
    }
    else 
        result = p[0][0];
    return result;
}

double Creat_M(RealMatrix p, int m, int n, int k)
{
    int len;
    int i, j;
    double mid_result = 0;
    int sign = 1;
    double *p_creat, *p_mid;
    len = (k-1)*(k-1);            
    p_creat = (double*)malloc(len*sizeof(double)); 
    p_mid = p_creat;
    for (i = 0; i < k; i++)
    {
        for (j = 0; j < k; j++)
        {
            if (i != m && j != n) 
            {
                *p_mid++ = p[i][j];
            }
        }
    }
    sign = (m+n)%2 == 0 ? 1 : -1;   
    mid_result = (double)sign*Det(p_creat, k-1);
    free(p_creat);
    return mid_result;
}

void matrixInversion(RealMatrix &invMat, RealMatrix Mat, int nRow){
    double det = 0;
    det = Det(Mat, nRow);
    if(det == 0){
        printf("Matrix determinant equals to zero, cannot compute invers.\n");
    
    }
    for(int i = 0; i < nRow; i++)
        for(int j = 0; j < nRow; j++){
            invMat[j][i] = Creat_M(Mat, i, j, nRow)/det;
    }
}

void matrixInversion(RealMatrix &invMat, RealMatrix Mat, int nRow, double det){
    
    for(int i = 0; i < nRow; i++)
        for(int j = 0; j < nRow; j++){
            invMat[j][i] = Creat_M(Mat, i, j, nRow)/det;
    }
}

double vecMatVec(RealVector& v1, RealMatrix &mat, RealVector& v2){
    int n = v1.size();
    double output = 0.;
    RealVector a;
    a.resize(n);
    for(int i = 0; i < n; i++){
        a[i] = 0;
        for(int j = 0; j < n; j++)
            a[i] += v1[i]* mat[j][i];
        output += a[i]*v2[i];
    }
    return output;
}

RealVector matVec(RealMatrix& mat, RealVector& vec){
    int nRow = mat.size();
    int nCol = mat[0].size();
    RealVector result;
    result.resize(nRow);
    for(int i = 0; i < nRow; i++){
        result[i] = 0;
        for(int j = 0;j < nCol; j++)
            result[i] += mat[i][j]*vec[j];
    }
    return result;
}

void vectorAdd(RealVector& v1, RealVector& v2, RealVector& v3){
    // RealVector results;
    // results.resize(v1.size());
    for(int i = 0; i < v1.size(); i++)
        v3[i] = v1[i]+v2[i];
   
}