/*
Random 3-diagonal symmetric matrix spawner with known eigenvalues
Author: Igor Kucherenko <ivkmailbox@gmail.com>

Idea: generating matrix T*D*T1,
where
   D is diagonal matrix,
   T is orthogonal matrix
   T1 = T^-1 = T^t

Now we need to generate 3-diag orthogonal matrix...

  a b
 -b a

devided by sqrt(a^2 + b^2) is orthogonal.

T = O1*O2*...*OM

O2i is odd-twoblock, O2i+1 is even-twoblock.
*/

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

//#define MATLAB_PRINT

#define MAX_N 8000
#define DEFAULT_N 10

#define DEFAULT_M 16

#define DEFAULT_MAX_DIAG_VALUE 10
#define DEFAULT_MAX_ROTATION_TAN 10

#define DEFAULT_MATRIX_FILE_NAME "testgen_s_matrix.txt"
#define DEFAULT_RESULT_FILE_NAME "testgen_s_result.txt"


int usage(FILE* fl)
{
   fprintf(fl, "Usage: testgen_s [options] [dimension [matrix_file [result_file]]]\n");
   fprintf(fl, "Where list of options are:\n");
   fprintf(fl, " -p         make positive definite matrix\n");
   fprintf(fl, " -z num     make singular matrix with <= num zero eigenvalues\n");
   fprintf(fl, " -d num     set max diagonal value (num is positive integer, default %d)\n",DEFAULT_MAX_DIAG_VALUE);
   fprintf(fl, " -t num     set max rotation tangens (num is positive integer, default %d)\n",DEFAULT_MAX_ROTATION_TAN);
   fprintf(fl, " -m num     set number of iterations (num is positive integer, default is %d)\n",DEFAULT_M);
   fprintf(fl, " -s num     set seed for PRNG (unsigned, default time(0))\n");
   fprintf(fl, " -h, -?     print this help and exit\n\n");
   fprintf(fl, "Default dimension is %d, max dimension is %d, ", DEFAULT_N, MAX_N);
   fprintf(fl, "default matrix_file is %s, ", DEFAULT_MATRIX_FILE_NAME);
   fprintf(fl, "default result_file is %s. ", DEFAULT_RESULT_FILE_NAME);
   fprintf(fl, "Both dimension, matrix_file and result_file can not start from \'-\' character. ");
   fprintf(fl, "In the case of both options and parameters omitted program prints this help and exits.\n");
   return 0;
}


/*
Heapsort implementation.
*/
void sieve(double* v, int n, int i)
{
   double p;
   int s;
   int s0 = 2*i + 1;
   int s1 = s0 + 1;
   if ((s1 < n) && (v[s1] > v[s0])) s = s1;
   else s = s0;
   while ((s < n) && (v[i] < v[s])) {
      p = v[i]; v[i] = v[s]; v[s] = p;
      i = s;
      s0 = 2*i + 1;
      s1 = s0 + 1;
      if ((s1 < n) && (v[s1] > v[s0])) s = s1;
      else s = s0;
   }
}

void heapsort(double* v, int N)
{
   int t;
   double p;
   for (t = (N-1)/2; t >= 0; --t) {
      sieve(v,N,t);
   }
   for (t = N-1; t > 0; --t) {
      p = v[t]; v[t] = v[0]; v[0] = p;
      sieve(v,t,0);
   }
}


/*
3d matrix multiplication with special tuning - it accepts NULL values for
diagonals of argument matrixes (NULL main diagonal mean ones, NULL secondary
diagonal mean zeroes, NULL diagonal for result mean do not compute it).

Let A 3d matrix with AD main diagonal, AL lower diagonal, AR upper diagonal;
B is 3d matrix in same notation; C = A*B; then

CD[i] = AL[i-1]*BR[i-1] + AD[i]*BD[i] + AR[i]*BL[i] // i=0 and i = N-1 are exeptions
CR[i] = AD[i]*BR[i] + AR[i]*BD[i+1] // i < N-1
CL[i] = AL[i]*BD[i] + AD[i+1]*BL[i] // i < N-1

WARNING: in general case 3d matrix product is not 3d matrix - so only 3 diagonals
of result will be computed.
*/
void matrix_3d_mult(int N,
double* AL, double* AD, double* AR,
double* BL, double* BD, double* BR,
double* CL, double* CD, double* CR)
{
   int i;
   double AL0, AL1, AD1, AD2, AR1;
   double BR0, BR1, BD1, BD2, BL1;

   if (CL != NULL) {
      for (i = 0; i < N-1; ++ i) {
         // CL[i] = AL[i]*BD[i] + AD[i+1]*BL[i]
         if (AL == NULL) AL1 = 0;
         else AL1 = AL[i];
         if (BD == NULL) BD1 = 1;
         else BD1 = BD[i];
         if (AD == NULL) AD2 = 1;
         else AD2 = AD[i+1];
         if (BL == NULL) BL1 = 0;
         else BL1 = BL[i];
         CL[i] = AL1*BD1 + AD2*BL1;
      }
   }

   if (CD != NULL) {
      for (i = 0; i < N; ++ i) {
         // CD[i] = AL[i-1]*BR[i-1] + AD[i]*BD[i] + AR[i]*BL[i]
         if ((AL == NULL) || (i == 0)) AL0 = 0;
         else AL0 = AL[i-1];
         if ((BR == NULL) || (i == 0)) BR0 = 0;
         else BR0 = BR[i-1];
         if (AD == NULL) AD1 = 1;
         else AD1 = AD[i];
         if (BD == NULL) BD1 = 1;
         else BD1 = BD[i];
         if ((AR == NULL) || (i == N-1)) AR1 = 0;
         else AR1 = AR[i];
         if ((BL == NULL) || (i == N-1)) BL1 = 0;
         else BL1 = BL[i];
         CD[i] = AL0*BR0 + AD1*BD1 + AR1*BL1;
      }
   }

   if (CR != NULL) {
      for (i = 0; i < N-1; ++ i) {
         // CR[i] = AD[i]*BR[i] + AR[i]*BD[i+1]
         if (AD == NULL) AD1 = 1;
         else AD1 = AD[i];
         if (BR == NULL) BR1 = 0;
         else BR1 = BR[i];
         if (AR == NULL) AR1 = 0;
         else AR1 = AR[i];
         if (BD == NULL) BD2 = 1;
         else BD2 = BD[i+1];
         CR[i] = AD1*BR1 + AR1*BD2;
      }
   }
}

void matrix_3d_mult_left(int N,
double* AL, double* AD, double* AR,
double* B,
double* C)
{
   int i, j;
   double AL0, AD1, AR1;

   for (i = 0; i < N; ++ i) {
      if ((AL == NULL) || (i == 0)) AL0 = 0;
      else AL0 = AL[i-1];
      if (AD == NULL) AD1 = 1;
      else AD1 = AD[i];
      if ((AR == NULL) || (i == N-1)) AR1 = 0;
      else AR1 = AR[i];
      if ((i > 0) && (i < N-1)) {
         for (j = 0; j < N; ++ j)
            C[i*N + j] = (AL0*B[(i-1)*N + j] + AD1*B[i*N + j]) + AR1*B[(i+1)*N + j];
      } else if (i == 0) {
         for (j = 0; j < N; ++ j)
            C[j] = AD1*B[j] + AR1*B[N + j];
      } else { // i == N-1
         for (j = 0; j < N; ++ j)
            C[N*(N-1) + j] = AL0*B[(N-2)*N + j] + AD1*B[(N-1)*N + j];
      }
   }
}

void matrix_3d_mult_right(int N,
double* A,
double* BL, double* BD, double* BR,
double* C)
{
   int i, j;
   double BR0, BD1, BL1;

   for (j = 0; j < N; ++ j) {
      if ((BR == NULL) || (j == 0)) BR0 = 0;
      else BR0 = BR[j-1];
      if (BD == NULL) BD1 = 1;
      else BD1 = BD[j];
      if ((BL == NULL) || (j == N-1)) BL1 = 0;
      else BL1 = BL[j];
      if ((j > 0) && (j < N-1)) {
         for (i = 0; i < N; ++ i)
            C[i*N + j] = (A[i*N + (j-1)]*BR0 + A[i*N + j]*BD1) + A[i*N + (j+1)]*BL1;
      } else if (j == 0) {
         for (i = 0; i < N; ++ i)
            C[i*N] = A[i*N]*BD1 + A[i*N + 1]*BL1;
      } else { // j == N-1
         for (i = 0; i < N; ++ i)
            C[i*N + (N-1)] = A[i*N + (N-2)]*BR0 + A[i*N + (N-1)]*BD1;
      }
   }
}

void matrix_symmetrize(int N, double* A)
{
   int i, j;
   double a;
   for (i = 1; i < N; ++ i) {
      for (j = 0; j < i; ++ j) {
         a = (A[N*i + j] + A[N*j + i])/2.;
         A[N*i + j] = a;
         A[N*j + i] = a;
      }
   }
}


void print_matrix_3d(FILE* f, int N, double* AL, double* AD, double* AR)
{
   double AL0, AD1, AR1;
   int i, j;
#ifndef MATLAB_PRINT
   fprintf(f,"%d\n",N);
#else
   fprintf(f,"A = [ ");
#endif
   if (AD == NULL) AD1 = 1;
   else AD1 = AD[0];
   if (AR == NULL) AR1 = 0;
   else AR1 = AR[0];
   fprintf(f,"%20.16lf %20.16lf", AD1, AR1);
   for (j = 2; j < N; ++ j) fprintf(f," 0");
#ifndef MATLAB_PRINT
   fprintf(f,"\n");
#else
   fprintf(f,"; ");
#endif
   for (i = 1; i < N-1; ++ i) {
      if (AL == NULL) AL0 = 0;
      else AL0 = AL[i-1];
      if (AD == NULL) AD1 = 1;
      else AD1 = AD[i];
      if (AR == NULL) AR1 = 0;
      else AR1 = AR[i];
      for (j = 0; j < i-1; ++ j) fprintf(f,"0 ");
      fprintf(f,"%20.16lf %20.16lf %20.16lf", AL0, AD1, AR1);
      for (j = i+2; j < N; ++ j) fprintf(f," 0");
#ifndef MATLAB_PRINT
      fprintf(f,"\n");
#else
      fprintf(f,"; ");
#endif
   }
   for (j = 2; j < N; ++ j) fprintf(f,"0 ");
   if (AL == NULL) AL0 = 0;
   else AL0 = AL[N-2];
   if (AD == NULL) AD1 = 1;
   else AD1 = AD[N-1];
#ifndef MATLAB_PRINT
   fprintf(f,"%20.16lf %20.16lf\n", AL0, AD1);
#else
   fprintf(f,"%20.16lf %20.16lf ]\n", AL0, AD1);
#endif
}


void print_matrix(FILE* f, int N, double* A)
{
#ifndef MATLAB_PRINT
   int i, j;
   fprintf(f,"%d\n",N);
   for (i = 0; i < N; ++ i) {
      fprintf(f,"%14.10lf",A[i*N]);
      for (j = 1; j < N; ++ j)
         fprintf(f," %14.10lf",A[i*N+j]);
      fprintf(f,"\n");
   }
#else
   int i, j;
   fprintf(f,"A = [ ");
   for (i = 0; i < N; ++ i) {
      fprintf(f,"%14.10lf",A[i*N]);
      for (j = 1; j < N; ++ j)
         fprintf(f," %14.10lf",A[i*N+j]);
      if (i < N-1) fprintf(f,"; ");
   }
   fprintf(f," ]\n");
#endif
}



int main(int argc, char** argv)
{
   int N; // system dimension

   FILE* matrix_file; // first output
   FILE* result_file; // second output

   char* matrix_file_name;
   char* result_file_name;

   int max_diag_value;
   int max_rotation_tan;

   int make_positive_matrix;
   int make_singular_matrix;

   int num_diag_zeroes;

   int M;

   int i, j, m;
   double r, a, b;

   double* D; // diagonal of D matrix

   double* OD;
   double* OR;
   double* OL;

   double* A; // first result matrix
   double* B; // second result matix

   unsigned seed;
   int seed_set;


   // setting defaults

   N = DEFAULT_N;

   matrix_file_name = DEFAULT_MATRIX_FILE_NAME;
   result_file_name = DEFAULT_RESULT_FILE_NAME;

   max_diag_value = DEFAULT_MAX_DIAG_VALUE;
   max_rotation_tan = DEFAULT_MAX_ROTATION_TAN;

   make_positive_matrix = 0;
   make_singular_matrix = 0;

   num_diag_zeroes = 0;

   M = DEFAULT_M;

   seed_set = 0;


   // parsing command line

   if (argc == 1) { // no parametrs - no work
      usage(stdout);
      return 0;
   }

   j = 0;
   for (i = 1; i < argc; ++ i) {
      if (argv[i][0] != '-') { // either dimension of file name
         if (j == 0) { // dimension
            if (sscanf(argv[i],"%d",&N) != 1) {
               fprintf(stderr, "Bad dimension parameter (%s)!\n",argv[i]);
               usage(stderr);
               return __LINE__;
            }
            if (N < 2) {
               fprintf(stderr, "Dimension parameter is too small (%d)!\n",N);
               return __LINE__;
            }
            if (N > MAX_N) {
               fprintf(stderr, "Dimension parameter is too big (%d)!\n",N);
               return __LINE__;
            }
            ++ j;
            continue;
         } else if (j == 1) { // matrix file
            matrix_file_name = argv[i];
            ++ j;
            continue;           
         } else if (j == 2) { // result file
            result_file_name = argv[i];
            ++ j;
            continue;            
         } else { // error
            fprintf(stderr, "Bad parameters string - too many non-\'minus\' prefixed parameters!\n");
            usage(stderr);
            return __LINE__;
         }
      } else { // some option
         if (argv[i][1] == 0)  {
            fprintf(stderr, "Bad option %s!\n", argv[i]);
            usage(stderr);
            return __LINE__;
         }
         if (argv[i][2] != 0)  {
            fprintf(stderr, "Bad option %s!\n", argv[i]);
            usage(stderr);
            return __LINE__;
         }
         // non-argumented options
         if ((argv[i][1] == 'h') || (argv[i][1] == '?')) {
            usage(stdout);
            return 0;
         }
         if (argv[i][1] == 'p') {
            make_positive_matrix = 1;
            continue;
         }
         // argumented options
         if (argv[i][1] == 'd') {
            if (i == argc - 1) {
               fprintf(stderr,"Error: missing -d argument!\n");
               usage(stderr);
               return __LINE__;
            }
            if (sscanf(argv[i+1],"%d",&max_diag_value) != 1) {
               fprintf(stderr,"Error: can not read -d argument!\n");
               usage(stderr);
               return __LINE__;
            }
            if (max_diag_value < 2) {
               fprintf(stderr,"Error: bad -d argument (%d)!\n", max_diag_value);
               usage(stderr);
               return __LINE__;
            }
            ++ i;
            continue;
         }
         if (argv[i][1] == 't') {
            if (i == argc - 1) {
               fprintf(stderr,"Error: missing -t argument!\n");
               usage(stderr);
               return __LINE__;
            }
            if (sscanf(argv[i+1],"%d",&max_rotation_tan) != 1) {
               fprintf(stderr,"Error: can not read -t argument!\n");
               usage(stderr);
               return __LINE__;
            }
            if (max_rotation_tan < 1) {
               fprintf(stderr,"Error: bad -t argument (%d)!\n", max_rotation_tan);
               usage(stderr);
               return __LINE__;
            }
            ++ i;
            continue;
         }
         if (argv[i][1] == 'm') {
            if (i == argc - 1) {
               fprintf(stderr,"Error: missing -m argument!\n");
               usage(stderr);
               return __LINE__;
            }
            if (sscanf(argv[i+1],"%d",&M) != 1) {
               fprintf(stderr,"Error: can not read -m argument!\n");
               usage(stderr);
               return __LINE__;
            }
            if (M < 0) {
               fprintf(stderr,"Error: bad -m argument (%d)!\n", M);
               usage(stderr);
               return __LINE__;
            }
            ++ i;
            continue;
         }
         if (argv[i][1] == 'z') {
            if (i == argc - 1) {
               fprintf(stderr,"Error: missing -z argument!\n");
               usage(stderr);
               return __LINE__;
            }
            if (sscanf(argv[i+1],"%d",&num_diag_zeroes) != 1) {
               fprintf(stderr,"Error: can not read -z argument!\n");
               usage(stderr);
               return __LINE__;
            }
            if (num_diag_zeroes < 1) {
               fprintf(stderr,"Error: bad -z argument (%d)!\n", num_diag_zeroes);
               usage(stderr);
               return __LINE__;
            }
            make_singular_matrix = 1;
            ++ i;
            continue;
         }
         if (argv[i][1] == 's') {
            if (i == argc - 1) {
               fprintf(stderr,"Error: missing -s argument!\n");
               usage(stderr);
               return __LINE__;
            }
            if (sscanf(argv[i+1],"%u",&seed) != 1) {
               fprintf(stderr,"Error: can not read -s argument!\n");
               usage(stderr);
               return __LINE__;
            }
            seed_set = 1;
            ++ i;
            continue;
         }
         // unknown option
         fprintf(stderr,"Error: unknown option %s\n", argv[i]);
         usage(stderr);
         return __LINE__;
      }
   }


   // opening files

   matrix_file = fopen(matrix_file_name, "w");

   if (matrix_file == NULL) {
      fprintf(stderr, "Error: could not open matrix file %s for writing!\n",matrix_file_name);
      return __LINE__;
   }

   result_file = fopen(result_file_name, "w");

   if (result_file == NULL) {
      fprintf(stderr, "Error: could not open result file %s for writing!\n",result_file_name);
      fclose(matrix_file);
      return __LINE__;
   }


   // allocating memory

   i = (2*N*N + 4*N - 2)*sizeof(double); // 2*N*N + 4*N - 2
   if ((D = malloc(i)) == NULL) {
      fprintf(stderr, "Error: could not allocate %d bytes of memory!\n",i);
      fclose(matrix_file);
      fclose(result_file);
      return __LINE__;
   }

   A = D + N;
   B = A + (N*N);
   OL= B + (N*N);
   OR=OL + (N-1);
   OD=OR + (N-1);


   // seeding PRNG

   if (seed_set) {
      srand(seed);
   } else {
      srand((unsigned)time(0));
   }


   // spawning D matrix

   // spawning diagonal values
   for (i = 0; i < N; ++ i) {
      D[i] = (rand() % (max_diag_value - 1)) + 1;
   }

   // making non positive
   if (!make_positive_matrix) {
      D[0] *= -1;
      for (i = 1; i < N; ++ i) {
         if (rand() % 2) {
            D[i] *= -1;
         }
      }
   }

   // making singular
   if (make_singular_matrix) {
      D[N-1] = 0;
      for (i = 0; i < num_diag_zeroes-1; ++ i) {
         j = (rand() % (N-2)) + 1;
         D[j] = 0;
      }
   }


   // initing A matrix

   for (i = 0; i < N*N; ++ i) A[i] = 0;

   for (i = 0; i < N; ++ i) A[i*(N+1)] = D[i]; // A = D


   // randomization loop

   for (m = 0; m < M; ++ m) {

      // spawning O matrix
      for (i = 0; i < N-1; ++ i) {
         OD[i] = 1;
         OR[i] = 0;
         OL[i] = 0;
      }
      OD[N-1] = 1;

      for (i = (m % 2); i < N-1; i += 2) {
         a = 1.0 + (rand() % (max_rotation_tan - 1));
         b = 1.0 + (rand() % (max_rotation_tan - 1));
         r = sqrt(a*a + b*b);
         if (rand() % 2) {
            OD[i] = a/r;
            OD[i + 1] = a/r;
         } else {
            OD[i] = -a/r;
            OD[i + 1] = -a/r;
         }
         if (rand() % 2) {
            OR[i] = b/r;
            OL[i] = -b/r;
         } else {
            OR[i] = -b/r;
            OL[i] = b/r;
         }
	   }

      // A = O'*A*O
      matrix_3d_mult_right(N,A,OL,OD,OR,B); // B = A*O
      matrix_3d_mult_left(N,OR,OD,OL,B,A); // A = O'*B
      
      // this to prevent error accumulation due to double non-associativity
      matrix_symmetrize(N,A);
   }


   // printing matrix

   print_matrix(matrix_file,N,A);


   // making result

   heapsort(D, N);


   // printing result

   fprintf(result_file,"%d\n",N);
   for (i = 0; i < N; ++ i) {
      fprintf(result_file,"%1.9lf\n",D[i]);
   }


   // freeing memory

   free(D);


   // closing all

   fclose(matrix_file);
   fclose(result_file);

   return 0;
}

