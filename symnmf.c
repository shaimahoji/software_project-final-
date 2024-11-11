#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "symnmf.h"


/* Allocating memory for a matrix and initializing all its elements to zero */
Matrix createMatrix(int rows, int cols) {
    Matrix matrix;
    int i, j, k;

    matrix.rows = rows;
    matrix.cols = cols;

    /* Allocate memory for the matrix rows */
    matrix.mat = (double**)malloc(rows * sizeof(double *));
    if (matrix.mat == NULL) {
        fprintf(stderr, "An Error Has Occurred");
        return EXIT_FAILURE;
    }

    /* Allocating memory for each row and initializing its elements to zero */
    for (i = 0; i < rows; i++) {
        matrix.mat[i] = (double *)malloc(cols * sizeof(double));
        if (matrix.mat[i] == NULL) {
            fprintf(stderr,"An Error Has Occurred");

            /* Free previously allocated rows in case of failure */
            for (k = 0; k < i; k++) {
                free(matrix.mat[k]);
            }
            free(matrix.mat);

            return EXIT_FAILURE;
        }

        for (j = 0; j < cols; j++) {
            matrix.mat[i][j] = 0.0;
        }
    }

    return matrix;
}

/* Function to create an empty matrix of size rows x cols */
Matrix readMatrix(char filename[]){
    FILE *file;
    int rows = 0, cols = 0;
    char line[256];
    Matrix mat;
    int i;
    file = fopen(filename, "r"); /* Open the file for reading */
    if (file == NULL) {
        mat.cols = 0;
        mat.rows = 0;
        mat.mat = NULL;
        return mat;}

    /* Determine the number of rows and columns */
    while (fgets(line, sizeof(line), file)) {
        int current_cols = 0;
        char *token = strtok(line, ",");       
        rows++;
        while (token) {
            current_cols++;
            token = strtok(NULL, ",");
        }
        if (current_cols > cols) {
            cols = current_cols;/* Update the maximum columns found */
        }
    }
    mat = createMatrix(rows,cols);
    rewind(file); /* Reset file pointer to the beginning of the file */
    i = 0;
    while (fgets(line, sizeof(line), file) && i < rows) {
        int col = 0;
        char *token = strtok(line, ",");
        while (token!=NULL) {
            mat.mat[i][col++] = (double)atof(token); 
            token = strtok(NULL, ",");}
        i++;
    }
    fclose(file); /* Close the file */
    return mat;
    }

/* Calculating the Euclidean distance between two points */
double euclidDistance(double *point1, double *point2, int n) {
    double sum = 0.0;
    int i;
    for (i = 0 ; i<n ;i++) {
        double diff = point1[i]- point2[i];
        sum += diff * diff;
    }
    return sqrt(sum);
}

/* Function to print a matrix */
void printMatrix(Matrix matrix){
    int i,j;
    for (i = 0; i < matrix.rows; i++) {
        for (j = 0; j < matrix.cols; j++) {
            if (j+1 == matrix.cols){printf("%.4f", matrix.mat[i][j]);}
            else{printf("%.4f,", matrix.mat[i][j]);}
        }
        printf("\n");
    }
}

/* Free allocated memory for a matrix */
void freeMatrix(Matrix matrix){
    int i;
    for (i = 0; i < matrix.rows; i++) {
        free(matrix.mat[i]);
    }
    free(matrix.mat);
}

/* Multiplying two matrices */
Matrix multiplyMatrices(Matrix A, Matrix B) {
    int i,j,k;
    Matrix C;
    
    /* Check if multiplication is possible (A.cols must equal B.rows) */
    if (A.cols != B.rows) {
        fprintf(stderr, "An Error Has Occurred");
        return EXIT_FAILURE;
    }
    
    /* Create the multiplication result matrix, C, with the dimensions A.rows x B.cols */
    C = createMatrix(A.rows, B.cols);
    
    /* Perform matrix multiplication */
    for (i = 0; i < A.rows; i++) {
        for (j = 0; j < B.cols; j++) {
            C.mat[i][j] = 0;  
            for (k = 0; k < A.cols; k++) {
                C.mat[i][j] += A.mat[i][k] * B.mat[k][j];
            }
        }
    }

    return C;
}

/* Helper function to divide each element of matrix A by the corresponding element in matrix B */
Matrix divisionMatrices_ij(Matrix A, Matrix B){

    Matrix C;
    int i,j;

    /* Check if division is possible (B must be square and invertible) */
    if ((B.rows != A.rows)||(A.cols!=B.cols)) {
        fprintf(stderr, "An Error Has Occurred");
        return EXIT_FAILURE;
    }

    C = createMatrix(A.rows,A.cols);

    for (i = 0; i < A.rows; i++){
        for (j = 0; j < A.cols; j++){
                if(B.mat[i][j] == 0){        
                    fprintf(stderr, "An Error Has Occurred");
                    return EXIT_FAILURE;
                }
                C.mat[i][j] = A.mat[i][j]/B.mat[i][j];
            }     
        }
    
    return C;
}

/* Helper function to multiply each element of matrix A by the corresponding element in matrix B */
Matrix multiplyMatrices_ij(Matrix A, Matrix B){

    Matrix C;
    int i, j;

    /* Check if division is possible (B must be square and invertible) */
    if ((B.rows != A.rows)||(A.cols!=B.cols)) {
        fprintf(stderr, "An Error Has Occurred");
        return EXIT_FAILURE;
    }

    C = createMatrix(A.rows,A.cols);
    for (i=0; i<A.rows;i++){
        for (j=0; j<A.cols; j++){
                C.mat[i][j] =A.mat[i][j]*B.mat[i][j];
            }     
        }
    
    return C;
}

/* Helper function to multiply each element of a matrix by a constant */
Matrix constMultMat(Matrix A, double num){

    int i,j;
    for (i = 0; i < A.rows; i++) {
        for (j = 0; j < A.cols; j++) {
            A.mat[i][j] = A.mat[i][j] * num;  
        }
    }
    return A;
}

/* Helper function to add a constant to each element of a matrix */
Matrix constAddMat(Matrix A, double num){

    int i,j;
    for (i = 0; i < A.rows; i++) {
        for (j = 0; j < A.cols; j++) {
            A.mat[i][j] = A.mat[i][j] + num;  
        }
    }
    return A;
}

/* Helper function to compute the transpose of a matrix */
Matrix transposeMatrix(Matrix A) {

    /* Creating the transposed matrix with swithced dimensions */
    Matrix A_T = createMatrix(A.cols, A.rows);

    int i,j;

    /*Filling the values of the transposed matrix*/
    for (i = 0; i < A.rows; i++) {
        for (j = 0; j < A.cols; j++) {
            A_T.mat[j][i] = A.mat[i][j];
        }
    }

    return A_T;
}

/* Update H : 1.4.2*/
Matrix H_next_t(Matrix H_prev_t,Matrix W){
    double beta = 0.5;
    Matrix fract;
    Matrix next_H;
    Matrix H_prev_t_T = transposeMatrix(H_prev_t);
    
    /* The Numerator (top part) of the fraction in the rule (in the project file)*/
    Matrix W_mult_prevH = multiplyMatrices(W,H_prev_t);

    /* Previous H (Hi,j to the power of t) multiplied (matrix multiplication) by its transpose*/
    Matrix pH_pHT_pH = multiplyMatrices(H_prev_t,H_prev_t_T);

    /*The Denominator (bottom part) of the fraction*/
    Matrix pH_pHT_pH2 = multiplyMatrices(pH_pHT_pH,H_prev_t);

    /*The fraction itself*/
    fract = constMultMat(divisionMatrices_ij(W_mult_prevH,pH_pHT_pH2),beta);
    
    fract = constAddMat(fract,(1-beta));

    /*The desrired result (H sub i,j to the power of t+1)*/
    next_H = multiplyMatrices_ij(H_prev_t,fract);
    
    freeMatrix(fract);
    freeMatrix(W_mult_prevH);
    freeMatrix(pH_pHT_pH);
    freeMatrix(pH_pHT_pH2);

    return next_H;
}

/*Convergence : 1.4.3*/
Matrix converge(Matrix init_H,Matrix W){

    Matrix prev_H = init_H;
    int i, j, k; 
    int max_iter = 300;
    double epsilon = 1e-4;
    double frobe_squared = 0.0;

    Matrix current_H = createMatrix(prev_H.rows, prev_H.cols);

    /* Iterating until max_iter or convergence */
    for(i = 0; i < max_iter; i += 1){

        current_H = H_next_t(prev_H, W);

        /* Calculating squared ||.||F */
        for (k = 0; k < current_H.rows; k += 1) {
            for (j = 0; j < current_H.cols; j++) {
                double value = current_H.mat[k][j] - prev_H.mat[k][j];
                frobe_squared += value * value; 
            }
        }

        /* Checking for convergence */
        if(frobe_squared < epsilon){
            break;
        }
        
        frobe_squared = 0.0;
        freeMatrix(prev_H);
        
        prev_H = current_H;
    }
    return current_H;
}

/* Compute the Similarity Matrix */
Matrix sym(Matrix matrix){
    int n  = matrix.rows; 
    Matrix sym = createMatrix(n, n);   /* n is the number of data points */
    int i;
    int j;

    /* Loop through each pair of rows (data points) to compute the similarity */
    for (i = 0; i < n; i += 1) {
        for (j = 0; j < n; j += 1) {

            if(i != j) {
                double dist = euclidDistance(matrix.mat[i],matrix.mat[j],matrix.cols);
                sym.mat[i][j] = exp(-(dist*dist)/2);
            }

            else {
                sym.mat[i][j] = 0; /* Setting diagonal values to 0 */
            }
        }
    }
    
    return sym;  
}

/* Sum of two data points (used for row summation) */
double sum(double *arr, int size) {

    double total = 0.0;
    int i;

    for (i = 0; i < size; i++) {
        total += arr[i];  /* Dereferencing the pointer to get the value */
    }

    return total;
}

/* Compute the diagonal degree Matrix */
Matrix ddg(Matrix matrix){

    /* First, computing the similarity matrix */
    Matrix sym_mat = sym(matrix);

    /* n is number of data points */
    int n  = sym_mat.rows;
    Matrix ddg = createMatrix(n,n);

    int i;
    int j;

    /* Fill in the diagonal of ddg with the sum of corresponding row elements from sym_mat */
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            if(i==j) {
                /* Diagonal elements are the sum of the corresponding row in the similarity matrix */
                ddg.mat[i][j]= sum(sym_mat.mat[i], n);
            }
            else {
                /* Non-diagonal elements are set to 0 */
                ddg.mat[i][j] = 0;
            }
        }
    } 

    freeMatrix(sym_mat);
    return ddg;  
}

/* Compute the normalized similarity matrix */
Matrix norm(Matrix matrix){
    Matrix D = ddg(matrix);
    Matrix inverse_sqrt_D = createMatrix(D.rows,D.cols);
    Matrix A = sym(matrix);
    Matrix invers_mul_A;
    Matrix norm; 
    int i, j;

    for (i = 0; i < D.rows; i += 1){
        for (j = 0; j < D.cols; j += 1){
            if(i==j)
            {
                inverse_sqrt_D.mat[i][j] = (1 / sqrt(D.mat[i][j]));
            } 
             
        }
    }

    invers_mul_A = multiplyMatrices(inverse_sqrt_D,A);
    norm = multiplyMatrices(invers_mul_A,inverse_sqrt_D);
    
    freeMatrix(A);
    freeMatrix(inverse_sqrt_D);
    freeMatrix(D);
    freeMatrix(invers_mul_A);

    return norm;
}

int main(int argc,char *argv[]) {

    char goal[256];
    char filename[256];
    Matrix matrix;

    /* Check if the correct number of arguments is provided */
    if ((argc != 3)){return 1;}

    /* Parse the goal and filename from command-line arguments */
    strcpy(goal, argv[1]);
    strcpy(filename, argv[2]);

    matrix = readMatrix(filename);

    /* Based on the goal,  compute and print the corresponding matrix */
    if (strcmp(goal, "sym") == 0){printMatrix(sym(matrix));}
    if (strcmp(goal, "ddg") == 0){printMatrix(ddg(matrix));}
    if (strcmp(goal, "norm") == 0){printMatrix(norm(matrix));}

    return 0;
}