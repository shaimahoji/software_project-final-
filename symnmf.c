#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "symnmf.h"


/* Allocating space for a matrix and initializing all elements to zero */
Matrix create_matrix(int rows, int cols) {
    Matrix matrix;
    int i, j, k;

    matrix.rows = rows;
    matrix.cols = cols;

    /* Allocate memory for the matrix rows */
    matrix.mat = (float **)malloc(rows * sizeof(float *));
    if (matrix.mat == NULL) {
        fprintf(stderr, "Memory allocation failed for matrix rows.\n");
        matrix.rows = matrix.cols = 0;
        return matrix;
    }

    /* Allocate and initialize each row to zero */
    for (i = 0; i < rows; i++) {
        matrix.mat[i] = (float *)malloc(cols * sizeof(float));
        if (matrix.mat[i] == NULL) {
            fprintf(stderr, "Memory allocation failed for matrix columns.\n");

            /* Free any previously allocated rows in case of failure */
            for (k = 0; k < i; k++) {
                free(matrix.mat[k]);
            }
            free(matrix.mat);

            matrix.rows = matrix.cols = 0;
            return matrix;
        }

        /* Initialize each element in the row to zero */
        for (j = 0; j < cols; j++) {
            matrix.mat[i][j] = 0.0;
        }
    }

    return matrix;
}

/* Function to create an empty matrix of size rows x cols */
Matrix read_matrix(char filename[]){
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

    /* determine the number of rows and columns */
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
    mat = create_matrix(rows,cols);
    rewind(file); /* Reset file pointer to the beginning of the file */
    i = 0;
    while (fgets(line, sizeof(line), file) && i < rows) {
        int col = 0;
        char *token = strtok(line, ",");
        while (token!=NULL) {
            mat.mat[i][col++] = atof(token); /* Convert string to float */
            token = strtok(NULL, ",");}
        i++;
    }
    fclose(file); /* Close the file */
    return mat;
    }

double euclid_distance(float *point1, float *point2, int n) {
    double sum = 0.0;
    int i;
    for (i = 0 ; i<n ;i++) {
        double diff = point1[i]- point2[i];
        sum += diff * diff;
    }
    return sqrt(sum);
}

/* Print the matrix */
void print_matrix(Matrix matrix){
    int i;
    int j;
    for (i = 0; i < matrix.rows; i++) {
        for (j = 0; j < matrix.cols; j++) {
            if (j+1 == matrix.cols){printf("%.4f", matrix.mat[i][j]);}
            else{printf("%.4f,", matrix.mat[i][j]);}
        }
        printf("\n");
    }
}

/* Free allocated memory */
void free_matrix(Matrix matrix){
    int i;
    for (i = 0; i < matrix.rows; i++) {
        free(matrix.mat[i]);
    }
    free(matrix.mat);
}

/* Compute the Similarity Matrix */
Matrix sym(Matrix matrix){
    int n  = matrix.rows; 
    Matrix sym = create_matrix(n, n);   /* the number of data points */
    int i;
    int j;

    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            if(i!=j){
                float dist = euclid_distance(matrix.mat[i],matrix.mat[j],matrix.cols);
                sym.mat[i][j]= exp(-(dist*dist)/2);}
            else{sym.mat[i][j]=0;}
        }
    }
    
    return sym;  
}

/* the sum of two data points (a row in the matrix) */
float sum(float *arr, int size) {
    float total = 0.0;
    int i;
    for (i = 0; i < size; i++) {
        total += arr[i];  /* Dereference the pointer to get the value */
    }
    return total;
}

/* Compute the diagonal degree Matrix */
Matrix ddg(Matrix matrix){
    Matrix sym_mat = sym(matrix);
    int n  = sym_mat.rows;   /* number of data points */
    int i;
    int j;
    Matrix ddg = create_matrix(n,n);

    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            if(i==j){
                ddg.mat[i][j]= sum(sym_mat.mat[i],n);
                }
            else{ddg.mat[i][j]=0;}
        }
    } 
    return ddg;  
}

/* Function to multiply matrices A and B, resulting in matrix C */
Matrix multiply_matrices(Matrix A, Matrix B) {
    int i;
    int j;
    int k;
    Matrix C;
    
    /* Check if multiplication is possible (A.cols must equal B.rows) */
    if (A.cols != B.rows) {
        printf("Matrix multiplication is not possible due to incompatible dimensions.\n");
        exit(EXIT_FAILURE);
    }
    
    /* Create result matrix C with dimensions A.rows x B.cols */
    C = create_matrix(A.rows, B.cols);
    
    /* Perform the multiplication */
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

/* Compute the normalized similarity matrix */
Matrix norm(Matrix matrix){
    Matrix D = ddg(matrix);
    Matrix inverse_sqrt_D = create_matrix(D.rows,D.cols);
    Matrix A = sym(matrix);
    Matrix norm; 
    int i;
    int j;

    for (i=0; i<D.rows;i++){
        for (j=0; j<D.cols; j++){
            if(i==j){
                inverse_sqrt_D.mat[i][j] = 1 / sqrt(D.mat[i][j]);} 
             
        }
    }
    norm = multiply_matrices(multiply_matrices(inverse_sqrt_D,A),inverse_sqrt_D);
    return norm;
}

int main(int argc,char *argv[]) {
    char goal[256];
    char filename[256];
    Matrix matrix;
    
    if ((argc != 3)){return 1;}

    strcpy(goal, argv[1]);
    strcpy(filename, argv[2]);

    matrix = read_matrix(filename);

    if (strcmp(goal, "sym") == 0){print_matrix(sym(matrix));}
    if (strcmp(goal, "ddg") == 0){print_matrix(ddg(matrix));}
    if (strcmp(goal, "norm") == 0){print_matrix(norm(matrix));}

    return 0;
}