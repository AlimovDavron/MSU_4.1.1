#include "task_01_03.h"

#define AT(A, x, y, n) (A+(x)*(n)+(y))

int fl_d = 0,
        fl_e = 0,
        fl_p = 0,
        fl_t = 0,
        fl_h = 0,
        fl_q = 0;

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
                    return 3;
            }
        else if(argv[i][0] == '-')
            return 2;
        else{
            cnt++;
        }
    }

    if(cnt > 2)
        return 1;

    return 0;
}

int readInputData(char *inputFile, double** A, double** B, double** X, int *n)
{
    int i, checkInput;
    FILE *in = fopen(inputFile, "r");


    checkInput = fscanf(in, "%d", n);
    if(checkInput == EOF)
        return 4;
    if(checkInput == 0)
        return 3;

    if(*n <= 0)
        return 1;

    *A = malloc((*n) * (*n) * sizeof(double));
    *B = malloc((*n) * sizeof(double));
    *X = malloc((*n) * sizeof(double));

    for(i = 0; i < (*n)*(*n); i++) {
        checkInput = fscanf(in, "%lf", (*A + i));
        if (checkInput == EOF)
            return 2;
        if(checkInput == 0)
            return 3;
    }

    for(i = 0; i < (*n); i++){
        checkInput = fscanf(in, "%lf", (*B+i));
        if (checkInput == EOF)
            return 5;
        if(checkInput == 0)
            return 3;
    }

    // if(fscanf(in, "%c", &tmp) != EOF)
    //     return 4;

    return 0;
}

void writeAnswer(char *outputFile, int n, const double* X, int notExist){
    int i;
    FILE *out = fopen(outputFile, "w");
    if(notExist){
        fprintf(out, "%d\n", 0);
    }
    else {
        fprintf(out, "%d\n", n);
        for (i = 0; i < n; i++)
            fprintf(out, "%1.9lf\n", *(X + i));
    }
}

void printHelp(){
    printf("Usage: lss [input_file_name] [output_file_name] [options]\n"
           "Where options include:\n"
           "  -d        print debug messages [default OFF]\n"
           "  -e        print errors [default OFF]\n"
           "  -p        print matrix [default OFF]\n"
           "  -t        print execution time [default OFF]\n"
           "  -h, -?    print this and exit\n"
           "Default input_file_name value is lss_00_00_in.txt, default output_file_name value is lss_00_00_out.txt.");
}

void printMatrix(int n, int m, const double* A){
    for(int i = 0; i < n; i++)
    {
        for(int j = 0; j < m; j++)
        {
            printf("%1.9lf\n ", *AT(A, i, j, n));
        } printf("\n");
    }
}

int main(int argc, char* argv[]) {
    int n, setInput = 0;
    double *A, *B, *X, *tmp;
    char* inputFile = "lss_01_03_in.txt";
    char* outputFile = "lss_01_03_out.txt";

    switch (validateParameters(argc, argv))
    {
        case 1:
            if(fl_e) printf("ValidationError: Wrong syntax of parameters.\n");
            return 8;
        case 2:
            if(fl_e) printf("ValidationError. Wrong syntax of parameters.\n");
            return 3;
        case 3:
            if(fl_e) printf("ValidationError. Wrong syntax of parameters. There "
                            "is no such parameter.\n");
            return 4;

        default: break;
    }

    if(fl_q || fl_h){
        printHelp();
        return 0;
    }

    for(int i = 1; i < argc; i++){
        if(argv[i][0] != '-'){
            if(!setInput){
                inputFile = argv[i];
                if(!validateFile(inputFile))
                {
                    if(fl_e) printf("ValidationError: Wrong syntax of parameters.\n");
                    return 2;
                }
                setInput = 1;
            }
            else {
                outputFile = argv[i];
            }
        }
    }

    switch (readInputData(inputFile, &A, &B, &X, &n))
    {
        case 1:
            if(fl_e) printf("ValidationError. Incorrect input.\n");
            return 5;
        case 2:
            if(fl_e) printf("ValidationError. Incorrect number of coefficients. (Not enough) \n");
            return 6;
        case 3:
            if(fl_e) printf("ValidationError. One of the coefficients is not a number \n");
            return 7;
        default: break;
    }



    return 0;
}