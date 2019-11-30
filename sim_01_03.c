//
// Created by alimovdavron on 11/12/19.
//

#include "task_01_03.h"
#include <math.h>

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
            printf("%lf ", *At(matrix, i, j));
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
        result += *At(a, i, 0) * *At(b, i, 0);
    }
    return result;
}

double calculateNorm(Matrix a){
    int i;
    if(*a.m != 1)
        return -1;
    double result = 0;
    for(i = 0; i < *a.n; i++){
        result += (*At(a, i, 0)) * (*At(a, i, 0));
    }

    return sqrt(result);
}

void multiply(Matrix A, Matrix B, Matrix C){
    int i, j, k;
    *C.n = *A.n;
    *C.m = *B.m;
    for(i = 0; i < *A.n; i++){
        for(j = 0; j < *B.m; j++){
            *At(C, i, j) = 0;
            for(k = 0; k < *A.m; k++)
                *At(C, i, j) += *At(A, i, k) *
                        *At(B, k, j);
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
        *At(a1, j, 0) = *At(A, i, iteration);
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

void calculateSubMatrix(Matrix subMatrix, Matrix x, Matrix z){
    int i, j;

    for(i = 0; i < *subMatrix.n; i++){
        for(j = 0; j < *subMatrix.m; j++){
            *At(subMatrix, i, j) -= *At(z, i, 0) * *At(x, j, 0);
            *At(subMatrix, i, j) -= *At(x, i, 0) * *At(z, j, 0);
        }
    }
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

void nextStep(int iteration, Matrix target, Matrix a1, Matrix targetSubMatrix){
    int i, j, k, d, n; double norm = calculateNorm(a1);
    n = *target.n;
    for(i = iteration+1, k = 0; i < n; i++, k++){
        for(j = iteration+1, d = 0; j < n; j++, d++){
            *At(target, i, j) = *At(targetSubMatrix, k, d);
        }
    }

    for(i = iteration+1; i < n; i++){
        *At(target, iteration, i) = 0;
        *At(target, i, iteration) = 0;
    }

    *At(target, iteration, iteration+1) = norm;
    *At(target, iteration+1, iteration) = norm;
}


int sim_01_03(int n, double* A, double* tmp, double precision){
    Matrix target, a1, x, y, z, subMatrix, tmpSubMatrix;
    int i;

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

    for(i = 0; i < n - 2; i++){
        initA1(i, target, a1);
        if(calculateNorm(a1) > precision) {
            initX(a1, x);
            initSubMatrix(i, target, subMatrix);
            initY(subMatrix, x, y);
            initZ(x, y, z);
            calculateSubMatrix(subMatrix, x, z);
            nextStep(i, target, a1, subMatrix);
        }
    }
    printMatrix(target);
}