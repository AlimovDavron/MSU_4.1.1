//
// Created by alimovdavron on 11/12/19.
//

#include "task_01_03.h"

int evc_memsize_01_03(int n){
    return 1000*n*sizeof(double);
}

#define GET(A, n,  x, y) (A+(x)*(n)+(y))

void fillMatrix(double* A, int dim, double value){
    int i, j;
    for(i = 0; i < dim; i++){
        for(j = 0; j < dim; j++){
            *GET(A, dim, i, j) = value;
        }
    }
}

void LRTransformation(int n, double* A, double* L, double* R, int precision_iteration){
    int i;

    fillMatrix(L, n, 0);
    fillMatrix(R, n, 0);

    *GET(R, n, 0, 0) = *GET(A, n, 0, 0);
    *GET(R, n, 0, 1) = *GET(A, n, 0, 1);
    *GET(L, n, 0, 0) = 1;

    for(i = 1; i < precision_iteration; i++){
        *GET(L, n, i, i) = 1;
        *GET(L, n, i, i - 1) = *GET(A, n, i, i - 1) / (*GET(R, n, i - 1, i - 1));
        *GET(R, n, i, i) = *GET(A, n, i, i) - *GET(L, n, i, i - 1) * *GET(R, n, i - 1, i);
        if(i != precision_iteration-1)
            *GET(R, n, i, i + 1) = *GET(A, n, i, i + 1) - *GET(L, n, i, i - 1) * *GET(R, n, i - 1, i + 1);
    }
}

void fastMultiplyForLR(int n, double* R, double* L, double* C, int precision_iteration){
    int i;
    for(i = 0; i < precision_iteration; i++){
        if(i != precision_iteration-1) {
            *GET(C, n, i, i) = *GET(R, n, i, i) + *GET(R, n, i, i + 1) * (*GET(L, n, i + 1, i));
            *GET(C, n, i, i + 1) = *GET(R, n, i, i + 1);
        }
        if(i != 0){
            *GET(C, n, i, i - 1) = *GET(R, n, i, i) * (*GET(L, n, i, i - 1));
        }
    }
    *GET(C, n, precision_iteration - 1, precision_iteration - 1) = *GET(R, n, precision_iteration - 1,
                                                                        precision_iteration - 1);
}

int evc_01_03(int n, int max_iterations, double epsilon, double* A, double* E, double* tmp, double precision){
    int i;

    double *L = tmp, *R = tmp + n*n;

    int precision_iteration = n;
    for(i = 0; i < max_iterations; i++) {
        LRTransformation(n, A, L, R, precision_iteration);
        fastMultiplyForLR(n, R, L, A, precision_iteration);
    }
}

