#include "task_01_03.h"

int fl_d = 0,
        fl_e = 0,
        fl_p = 0,
        fl_t = 0,
        fl_h = 0,
        fl_q = 0,
        max_iter=0;

double precision = 1e-14, epsilon=1e-10;

int _strlen(char const *input){
    int length = 0;
    while(input[length]!='\0')
        length++;
    return length;
}

int validateFile(char* filename)
{
    FILE *f = fopen(filename, "rb");
    return f != NULL;
}

int validateParameters(int argc, char **argv){
    int i, cnt = 0;
    if(argc == 1)
        return 0;

    for(i = 1; i < argc; i++){
        if(argv[i][0] == '-' && _strlen(argv[i]) == 2)
            switch (argv[i][1])
            {
                case 'd':
                    fl_d = 1;
                    break;
                case 'h':
                    fl_h = 1;
                    break;
                case '?':
                    fl_q = 1;
                case 'e':
                    fl_e = 1;
                    break;
                case 'p':
                    fl_p = 1;
                    break;
                case 't':
                    fl_t = 1;
                    break;
                default:
                    return 2;
            }
        else if(argv[i][0] == '-'){
            int max_iter_temp; double precision_temp, epsilon_temp;
            if(sscanf(argv[i], "-max_iter=%d",&max_iter_temp) == 1){
                if(max_iter_temp >= 0){
                    max_iter = max_iter_temp;
                } else {
                    return 3;
                }
            }
            else if(sscanf(argv[i], "-pre=%lf", &precision_temp) == 1){
                precision = precision_temp;
            }
            else if(sscanf(argv[i], "-eps=%lf", &epsilon_temp) == 1){
                epsilon = epsilon_temp;
            } else {
                return 2;
            }
        }
        else{
            cnt++;
        }
    }

    if(cnt > 2)
        return 1;

    return 0;
}

int readInputData(char *inputFile, double** A, int *n)
{
    int i, checkInput;
    FILE *in = fopen(inputFile, "r");

    checkInput = fscanf(in, "%d", n);
    if(checkInput == EOF)
        return 6;

    if(checkInput == 0)
        return 7;

    if(*n <= 0)
        return 7;

    *A = malloc((*n) * (*n) * sizeof(double));

    for(i = 0; i < (*n)*(*n); i++) {
        checkInput = fscanf(in, "%lf", (*A + i));
        if (checkInput == EOF)
            return 8;
        if(checkInput == 0)
            return 9;
    }

    return 0;
}

void writeAnswer(char *outputFile, int n, const double* X, int result){
    int i;
    FILE *out = fopen(outputFile, "w");
    if(result == -1){
        fprintf(out, "%d\n", 0);
    }
    else {
        fprintf(out, "%d\n", n);
        for (i = 0; i < n; i++)
            fprintf(out, "%1.9lf\n", *(X + i));
    }
}

void printHelp(){
    printf("Usage: evc [input_file_name] [output_file_name] [options]\n"
           "Where options include:\n"
           " -d                print debug messages [default OFF]\n"
           " -e                print errors [default OFF]\n"
           " -p                print matrix [default OFF]\n"
           " -t                print execution time [default OFF]\n"
           " -prec=<num>       precision [default - 1e-14]\n"
           " -eps=<num>        epsilon [default - 1e-10]\n"
           " -max_iter=<num>   limit number of iterations [default - 0, i.e. not limit]\n"
           " -h, -?            print this and exit");
}

int main(int argc, char* argv[]) {
    int n, setInput = 0;
    double *A, *tmp, *E;
    char* inputFile = "input.txt";
    char* outputFile = "output.txt";

    switch (validateParameters(argc, argv))
    {
        case 1:
            printf("ValidationError: Wrong syntax of parameters. There are more than two filenames\n");
            return 1;
        case 2:
            printf("ValidationError. Wrong syntax of parameters. There is no such parameter or you haven't"
                   "set value to one of the parameters\n");
            return 2;
        case 3:
            printf("ValidationError. Wrong syntax of parameters. Max_iter must be non negative\n");
            return 3;
        default:
            break;
    }

    if(fl_q || fl_h){
        printHelp();
        return 0;
    }

    for(int i = 1; i < argc; i++){
        if(argv[i][0] != '-'){
            if (!setInput) {
                if(i != 1){
                    if(fl_e) printf("ValidationError: Wrong order of parameters.\n");
                    return 4;
                }
                inputFile = argv[i];
                if (!validateFile(inputFile)) {
                    if (fl_e) printf("ValidationError: There is no such file.\n");
                    return 5;
                }
                setInput = 1;
            } else {
                if(i != 2) {
                    if (fl_e) printf("ValidationError: Wrong order of parameters.\n");
                    return 4;
                }
                outputFile = argv[i];
            }

        }
    }

    switch (readInputData(inputFile, &A, &n))
    {
        case 6:
            if(fl_e) printf("ValidationError. File is empty.\n");
            return 6;
        case 7:
            if(fl_e) printf("ValidationError. n is not a positive integer.\n");
            return 7;
        case 8:
            if(fl_e) printf("ValidationError. Not enough elements in the matrix.\n");
            return 8;
        case 9:
            if(fl_e) printf("ValidationError. One of the elements of the matrix is not a number.\n");
            return 9;
        default:
            break;
    }

    tmp = malloc(sim_memsize_01_03(n));
    E = malloc(n * sizeof(double));

    clock_t begin = clock();
    if(sim_01_03(n, A, tmp, precision) == -1){
        return 10;
    }

    free(tmp);
    tmp = malloc(evc_memsize_01_03(n));

    int result = evc_01_03(n,0, epsilon, A, E ,tmp, precision);
    clock_t end = clock();

    double timeSpent = (double)(end - begin) / CLOCKS_PER_SEC;

    if(fl_t) printf("Execution time: %1.9lf\n", timeSpent);

    writeAnswer(outputFile, n, E, result);
    return 0;
}