//
// Created by alimovdavron on 11/12/19.
//

#include "task_01_03.h"


int evc_memsize_01_03(int n){
    return 1000*n*sizeof(double);
}

void fillMatrix(Matrix A, double value){
    int i, j;
    for(i = 0; i < *A.n; i++){
        for(j = 0; j < *A.m; j++){
            *At(A, i, j) = value;
        }
    }
}

void LRTransformation(Matrix A, Matrix L, Matrix R){
    int i, n = *A.n;

    fillMatrix(L, 0);
    fillMatrix(R, 0);

    *At(R, 0, 0) = *At(A, 0, 0);
    *At(R, 0, 1) = *At(A, 0, 1);
    *At(L, 0, 0) = 1;

    for(i = 1; i < n; i++){
        *At(L, i, i) = 1;
        *At(L, i, i-1) = *At(A, i, i-1)/(*At(R, i-1, i-1));
        *At(R, i, i) = *At(A, i, i) - *At(L, i, i-1) * *At(R, i-1, i);
        if(i != n-1)
            *At(R, i, i+1) = *At(A, i, i+1) - *At(L, i, i-1) * *At(R, i-1, i+1);
    }
}

void fastMultiplyForLR(Matrix R, Matrix L, Matrix C){
    int i, n = *R.n;
    for(i = 0; i < n; i++){
        if(i != n-1) {
            *At(C, i, i) = *At(R, i, i) + *At(R, i, i + 1) * (*At(L, i + 1, i));
            *At(C, i, i+1) = *At(R, i, i + 1);
        }
        if(i != 0){
            *At(C, i, i-1) = *At(R, i, i)*(*At(L, i, i-1));
        }
    }
    *At(C, n-1, n-1) = *At(R, n-1, n-1);
}

int evc_01_03(int n, int max_iterations, double epsilon, double* A, double* E, double* tmp, double precision){
    Matrix target, L, R;
    int i;

    target.n = (int*)tmp;
    target.m = (int*)tmp+1;
    target.matrix = A;
    *target.n = n;
    *target.m = n;

    L.matrix = (tmp+1);
    L.n = (int*)(tmp+1+n*n);
    L.m = (int*)(tmp+1+n*n)+1;

    R.matrix = (tmp+2+n*n);
    R.n = (int*)(tmp+2+2*n*n);
    R.m = (int*)(tmp+2+2*n*n)+1;

    *L.n = *L.m = *R.m = *R.n = n;


    for(i = 0; i < max_iterations; i++) {
        LRTransformation(target, L, R);
        fastMultiplyForLR(R, L, target);
    }

    printMatrix(target);
}

