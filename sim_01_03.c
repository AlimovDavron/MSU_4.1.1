//
// Created by alimovdavron on 11/12/19.
//

#include "task_01_03.h"
#include <math.h>

#define GET(A, x, y, n) (A+(x)*(n)+(y))

typedef struct {
    int* n,* m;
    double* matrix;
}Matrix;

double* At(Matrix A, int i, int j){
    return A.matrix + i * *(A.m) + j;
}

void printMatrix(Matrix matrix){
    int i, j;
    printf("-----\n");
    for(i = 0; i < *matrix.n; i++){
        for(j = 0; j < *matrix.m; j++){
            printf("%lf ", *GET(matrix.matrix, i, j, *matrix.m));
        }
        printf("\n");
    }
    printf("-----\n");
}

// be sure that target matrix has the same or more memory
void copyMatrix(Matrix original, Matrix target){
    *target.n = *original.n;
    *target.m = *original.m;

    int i, j;
    for(i = 0; i < *original.n; i++){
        for(j = 0; j < *original.m; j++){
            *At(target, i, j) = *At(original, i, j);
        }
    }
}

// requires n*m additional memory
void transpose(Matrix A, Matrix tmp){
    int i, j;
    for(i = 0; i < *A.n; i++){
        for(j = 0; j < *A.m; j++){
            *At(tmp, j, i) = *At(A, i, j);
        }
    }
    *tmp.n = *A.m;
    *tmp.m = *A.n;

    copyMatrix(tmp, A);
}

void multiplyByN(Matrix A, double n){
    int i, j;
    for(i = 0; i < *A.n; i++){
        for(j = 0; j < *A.m; j++){
            *At(A, i, j) *= n;
        }
    }
}

double scalarMultiplication(Matrix a, Matrix b){
    double result = 0; int i;
    for(i = 0; i < *a.n; i++){
        result += *GET(a.matrix, i, 0, 1)
                * *GET(b.matrix, i, 0, 1);
    }
    return result;
}

double calculateNorm(Matrix a){
    int i;
    if(*a.m != 1)
        return -1;
    double result = 0;
    for(i = 0; i < *a.n; i++){
        result += (*GET(a.matrix, i, 0, 1)) *
                  (*GET(a.matrix, i, 0, 1));
    }

    return sqrt(result);
}

void multiply(Matrix A, Matrix B, Matrix C){
    int i, j, k;
    *C.n = *A.n;
    *C.m = *B.m;
    for(i = 0; i < *A.n; i++){
        for(j = 0; j < *B.m; j++){
            for(k = 0; k < *A.m; k++)
                *GET(C.matrix, i, j, *B.m) += *GET(A.matrix, i, k, *A.m) *
                        *GET(B.matrix, k, j, *B.m);
        }
    }
}

void subtract(Matrix A, Matrix B){
    int i, j;
    for(i = 0; i < *A.n; i++){
        for(j = 0; j < *A.m; j++){
            *GET(A.matrix, i, j, *A.m) -=
                    *GET(B.matrix, i, j, *B.m);
        }
    }
}


int sim_memsize_01_03(int n){
    return 1000*n*sizeof(double);
}

void initA1(int iteration, Matrix A, Matrix a1){
    int i, j;
    *a1.m = 1;
    *a1.n = *A.n-iteration-1;
    for(i = iteration+1, j = 0; i < *A.n; i++, j++){
        *GET(a1.matrix, j, 0, 1) = *GET(A.matrix, i, iteration, *A.n);
    }
}

void initX(Matrix a1, Matrix x){
    int i; double norm;
    copyMatrix(a1, x);

    *At(x, 0, 0) -= calculateNorm(a1);

    norm = calculateNorm(x);
    for(i = 0; i < *x.n; i++){
        *At(x, i, 0) /= norm;
    }
}

void initY(Matrix A, Matrix x, Matrix y){
    multiply(A, x, y);
}

void initZ(Matrix x, Matrix y, Matrix z){
    int i; double scalar = scalarMultiplication(x, y);

    copyMatrix(y, z);
    multiplyByN(z, 2);

    for(i = 0; i < *z.n; i++){
        *At(z, i, 0) -= 2 * scalar * (*At(x, i, 0));
    }

}

void calculateSubMatrix(Matrix subMatrix, Matrix x, Matrix z, Matrix tmp){
    transpose(x, tmp);
    multiply(z, x, tmp);
    subtract(subMatrix, tmp);
    transpose(x, tmp);
    transpose(z, tmp);
    multiply(x, z, tmp);
    subtract(subMatrix, tmp);
    transpose(z, tmp);
}

void initSubMatrix(int iteration, Matrix A, Matrix subMatrix){
    int i, j, k, d, n;
    n = *A.n;
    *subMatrix.n = n - iteration - 1;
    *subMatrix.m = n - iteration - 1;
    for(i = iteration+1, k = 0; i < n; i++, k++){
        for(j = iteration+1, d = 0; j < n; j++, d++){
            *At(subMatrix, k, d) = *At(A, i, j);
        }
    }
}


int sim_01_03(int n, double* A, double* tmp, double precision){
    Matrix target, a1, x, y, z, subMatrix, tmpSubMatrix;
    int i, j;

    target.n = (int*)tmp;
    target.m = (int*)tmp+1;
    target.matrix = A;
    *target.n = n;
    *target.m = n;

    a1.matrix = tmp+1;
    a1.n = (int*)(tmp+1+n);
    a1.m = (int*)(tmp+1+n)+1;

    x.matrix=(tmp+2+n);
    x.n = (int*)(tmp+2+2*n);
    x.m = (int*)(tmp+2+2*n)+1;

    y.matrix=(tmp+3+2*n);
    y.n = (int*)(tmp+3+3*n);
    y.m = (int*)(tmp+3+3*n)+1;

    z.matrix=(tmp+4+3*n);
    z.n = (int*)(tmp+4+4*n);
    z.m = (int*)(tmp+4+4*n)+1;

    subMatrix.matrix = (tmp+5*n+4);
    subMatrix.n = (int*)(tmp+5*n+4+n*n);
    subMatrix.m = (int*)(tmp+5*n+4+n*n)+1;

    tmpSubMatrix.matrix = (tmp + 5*n + 5 + n*n);
    tmpSubMatrix.n = (int*)(tmp + 5*n + 5 + 2*n*n);
    tmpSubMatrix.m = (int*)(tmp + 5*n + 5 + 2*n*n) + 1;

    for(i = 0; i < n - 2 -1; i++){
        initA1(i, target, a1);
        initX(a1, x);
        initSubMatrix(i, target, subMatrix);
        initY(subMatrix, x, y);
        initZ(x, y, z);
        calculateSubMatrix(subMatrix, x, z, tmpSubMatrix);
        printMatrix(subMatrix);
    }
}