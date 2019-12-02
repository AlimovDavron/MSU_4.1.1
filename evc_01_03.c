//
// Created by alimovdavron on 11/12/19.
//

#include "task_01_03.h"

int evc_memsize_01_03(int n){
    return 10*n*n*sizeof(double);
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

int LRTransformation(int n, double* A, double* L, double* R, int precision_iteration, double precision){
    int i;

    fillMatrix(L, n, 0);
    fillMatrix(R, n, 0);

    *GET(R, n, 0, 0) = *GET(A, n, 0, 0);
    *GET(R, n, 0, 1) = *GET(A, n, 0, 1);
    *GET(L, n, 0, 0) = 1;

    for(i = 1; i < precision_iteration; i++){
        *GET(L, n, i, i) = 1;

        if(fabs(*GET(R, n, i - 1, i - 1)) < precision){
            return 0;
        }

        *GET(L, n, i, i - 1) = *GET(A, n, i, i - 1) / (*GET(R, n, i - 1, i - 1));

        if(fabs(*GET(L, n, i, i - 1) ) < precision){
            *GET(L, n, i, i - 1) = 0;
        }

        *GET(R, n, i, i) = *GET(A, n, i, i) - *GET(L, n, i, i - 1) * *GET(R, n, i - 1, i);

        if(fabs(*GET(R, n, i, i) ) < precision){
            *GET(R, n, i, i) = 0;
        }

        if(i != precision_iteration-1) {
            *GET(R, n, i, i + 1) = *GET(A, n, i, i + 1) - *GET(L, n, i, i - 1) * *GET(R, n, i - 1, i + 1);

            if(fabs(*GET(R, n, i, i + 1)) < precision){
                *GET(R, n, i, i + 1) = 0;
            }
        }
    }

    return 1;
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

double getNorm(int n, double* A, int precision_iteration) {
    int i, j;
    double norm = fabs(*GET(A, n, 0, 0)) + fabs(*GET(A, n, 0, 1));
    norm = fmax(norm,
                fabs(*GET(A, n, precision_iteration - 1, precision_iteration - 1)) +
                fabs(*GET(A, n, precision_iteration - 1, precision_iteration - 2)));

    for (i = 1; i < precision_iteration - 1; i++) {
        double tmp = 0;
        for (j = i - 1; j <= i + 1; j++) {
            tmp += fabs(*GET(A, n, i, j));
        }
        norm = fmax(tmp, norm);
    }
    return norm;
}

int checkPrecisionCondition(int n, double* A, double epsilon, int precision_iteration, double norm){
    if(fabs(*GET(A, n, precision_iteration-1, precision_iteration - 2)) > epsilon*norm)
        return 0;
    else return 1;
}

void addValueToDiagonal(int n, double *A, double precision_iteration, double value){
    int i, j;
    for(i = 0; i < precision_iteration; i++){
        *GET(A, n, i, i) += value;
    }
}

int compare(const double* a, const double* b){
    return *a > *b;
}

int evc_01_03(int n, int max_iterations, double epsilon, double* A, double* E, double* tmp, double precision){
    int i, accuracy_achieved = 0;
    double D;
    double *L = tmp, *R = tmp + n*n;

    int precision_iteration = n;
    for(i = 0; i < max_iterations || !max_iterations; i++) {
        while(checkPrecisionCondition(n, A, epsilon, precision_iteration, getNorm(n, A,precision_iteration)) && precision_iteration > 2){
            precision_iteration--;
        }

        if(precision_iteration == 2){
            double a11, a12, a21, a22;
            a11 = *GET(A, n, 0, 0);
            a12 = *GET(A, n, 0, 1);
            a21 = *GET(A, n, 1, 0);
            a22 = *GET(A, n, 1, 1);
            D = sqrt(((a11+a22)*(a11+a22)+4*a12*a21-4*a11*a22));
            *GET(A, n, 0, 0) = (a11+a22+D)/2;
            *GET(A, n, 1, 1) = (a11+a22-D)/2;
            accuracy_achieved = 1;
            break;
        }

        double sk = *GET(A, n, precision_iteration-1, precision_iteration-1);
        addValueToDiagonal(n, A, precision_iteration, -sk);

        if(!LRTransformation(n, A, L, R, precision_iteration, precision))
            return -1;

        fastMultiplyForLR(n, R, L, A, precision_iteration);
        addValueToDiagonal(n, A, precision_iteration, sk);
    }

    for(i = 0; i < n; i++){
        *GET(E, n, 0, i) = *GET(A, n, i, i);
        if(fabs(*GET(E, n, 0, i)) < precision){
            *GET(E, n, 0, i) = 0;
        }
    }

    qsort(E, n, sizeof(double), compare);

    return !accuracy_achieved;
}

