//
// Created by alimovdavron on 11/12/19.
//

#ifndef LRMETHOD_TASK_01_03_H
#define LRMETHOD_TASK_01_03_H

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "math.h"

int evc_01_03(int n, int max_iterations, double epsilon, double* A, double* E, double* tmp, double precision);
int sim_01_03(int n, double* A, double* tmp, double precision);

int sim_memsize_01_03(int n);
int evc_memsize_01_03(int n);

#endif //LRMETHOD_TASK_01_03_H
