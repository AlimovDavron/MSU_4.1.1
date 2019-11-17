//
// Created by alimovdavron on 11/12/19.
//

#ifndef LRMETHOD_TASK_01_03_H
#define LRMETHOD_TASK_01_03_H

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

int evc_01_03(int n, int max_iterations, double epsilon, double* A, double* E, double* tmp, double precision);
int sim_01_03(int n, double* A, double* tmp, double precision);

int sim_memsize_01_03(int n);
int evc_memsize_01_03(int n);

typedef struct {
    int* n,* m;
    double* matrix;
}Matrix;

double* At(Matrix A, int i, int j);
void copyMatrix(Matrix original, Matrix target); // be sure that target matrix has the same or more memory
void printMatrix(Matrix matrix);
void transpose(Matrix A, Matrix tmp); // requires n*m additional memory
void multiplyByN(Matrix A, double n);
void multiply(Matrix A, Matrix B, Matrix C);

#endif //LRMETHOD_TASK_01_03_H
